"""
v2_study.py
===========

Implements the experiments frozen in V2_PROTOCOL.md (E0 equivalence check, E1 standard cyclic
baseline, E2 monotonic route, E3 reference model). Surrogate only, Chaboche only.

    python v2_study.py --exp e0
    python v2_study.py --exp e1 --workers 8
    python v2_study.py --exp e2 --workers 8 [--cmode own|global|lodo] [--r 0.3 0.4 ...]
    python v2_study.py --exp e3

Results are appended to calibration_sessions/v2/<exp>[_<cmode>].jsonl (one JSON record per task;
finished tasks are skipped on re-run).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)
sys.path.insert(0, HERE)

from calibration_app.core import data_loader, objective          # noqa: E402
from calibration_app.core import search as search_mod            # noqa: E402
from calibration_app.core.pipeline import build_model            # noqa: E402
from calibration_app.utils.physics_gate import DataContext, PhysicsGate  # noqa: E402

E_TRUE = 200000.0
LD = 5
DIAS = (10, 12, 16)
STRAINS = (3, 4, 5, 6)
R_GRID = tuple(round(0.3 + 0.1 * i, 1) for i in range(8))
OUT = os.path.join(PROJ, "calibration_sessions", "v2")
MONO_DIR = os.path.join(PROJ, "synthetic_data", "v2")
SOBOL_N, TPE_N = 4096, 200
SEED = int(os.environ.get("V2_SEED", "42"))      # protocol value 42; others only for seed-noise runs


def bars():
    return [(d, n) for d in DIAS for n in STRAINS]


def cyc_path(d, n):
    return os.path.join(PROJ, "input_data", f"{d}mm dia", f"LD{LD}", f"{n}data.csv")


def area(d):
    return math.pi * d * d / 4.0


def gauge(d):
    return float(d * LD)


# ---------------------------------------------------------------- compliance
def _extract(d, n, C, tag):
    """Write the (optionally compliance-corrected) cycle-1 ramp; return (csv, peak %)."""
    import run_mono_to_loops as m2l
    os.makedirs(MONO_DIR, exist_ok=True)
    csv = os.path.join(MONO_DIR, f"mono_{d}_{n}_{tag}.csv")
    peak = m2l.extract_monotonic(cyc_path(d, n), csv, gauge(d), C)
    return csv, peak


def load_mono(d, n, C, tag):
    csv, _ = _extract(d, n, C, tag)
    return data_loader.load_any_experimental(
        csv, gauge_length_mm=gauge(d), bar_diameter_mm=d, sigma_y0=500.0)


def ramp_modulus_raw(d, n):
    m = load_mono(d, n, 0.0, "raw")
    return data_loader.monotonic_modulus(m["strain"], m["stress"])


def c_own(d, n):
    """Machine compliance (mm/kN) from the bar's own ramp, assuming E = E_TRUE."""
    e_ramp = ramp_modulus_raw(d, n)
    k = 1.0 / e_ramp - 1.0 / E_TRUE                      # 1/MPa
    return k * 1000.0 * gauge(d) / area(d), e_ramp


def compliance_table():
    p = os.path.join(OUT, "compliance.json")
    if os.path.exists(p):
        return json.load(open(p))
    os.makedirs(OUT, exist_ok=True)
    rows = {}
    for d, n in bars():
        c, e = c_own(d, n)
        rows[f"{d}_{n}"] = {"C_own_mm_per_kN": c, "E_ramp_MPa": e}
    json.dump(rows, open(p, "w"), indent=1)
    return rows


def get_C(d, n, mode):
    t = compliance_table()
    if mode == "own":
        return t[f"{d}_{n}"]["C_own_mm_per_kN"]
    if mode == "global":
        return float(np.mean([v["C_own_mm_per_kN"] for v in t.values()]))
    if mode == "lodo":
        return float(np.mean([v["C_own_mm_per_kN"] for k, v in t.items()
                              if int(k.split("_")[0]) != d]))
    if mode == "none":
        return 0.0
    raise ValueError(mode)


# ---------------------------------------------------------------- data / scoring
def load_cyc(d, n, C):
    data = data_loader.load_experimental(
        cyc_path(d, n), gauge_length_mm=gauge(d), bar_diameter_mm=d,
        sigma_y0=500.0, compliance_mm_per_kN=C)
    data["loading_mode"] = "cyclic"
    return data


def window(data, target):
    """Window simulated over cycles 2..max(target); only ``target`` cycles scored."""
    present = [int(c) for c in data["cycles_present"]]
    tgt = [c for c in target if c in present]
    if not tgt:
        return None
    mask = np.asarray(data["cycle"]) <= max(tgt)
    return {"strain": data["strain"][mask], "stress": data["stress"][mask],
            "cycle": data["cycle"][mask], "target_cycles": tgt,
            "burnin_strain": data.get("burnin_strain", np.array([], float))}


def score(model, params, sy0, win, E, weighted=False):
    if win is None:
        return float("nan")
    return float(objective.evaluate_surrogate(
        model, params, win, None, E=E, sy0=sy0, weight_reversals=weighted))


# ---------------------------------------------------------------- calibration
def calibrate(model, win, gate_data, sy0_bounds, bound_override, E, weighted=False):
    """Sobol + TPE, same stages as the app. Returns dict(params, sigma_y0, fit_objective)."""
    ctx = DataContext.from_datasets(gate_data, sigma_y0=sy0_bounds[1])
    gate = PhysicsGate(ctx)
    defaults = model.default_bounds()
    mb = {nme: bound_override.get(nme, defaults[nme]) for nme in model.param_names}
    names = ["sigma_y0"] + list(model.param_names)
    bounds = [tuple(sy0_bounds)] + [mb[nme] for nme in model.param_names]

    def evaluate(x):
        return objective.evaluate_surrogate(
            model, tuple(x[1:]), win, gate, E=E, sy0=float(x[0]),
            weight_reversals=weighted)

    seeds, s1 = search_mod.sobol_seed(evaluate, bounds, n_samples=SOBOL_N,
                                      top_k=10, seed=SEED)
    s2 = search_mod.run_bayesian(evaluate, bounds, names, n_trials=TPE_N,
                                 x0=seeds[0] if seeds else None, seed=SEED)
    best = s1 if s1.best_objective <= s2.best_objective else s2
    x = best.best_params
    sy0 = float(x[0])
    params = tuple(float(v) for v in gate.project(model, x[1:], sy0))
    return {"params": list(params), "sigma_y0": sy0, "fit_objective": float(best.best_objective)}


# ---------------------------------------------------------------- tasks
def task_e1(args):
    d, n, N = args
    C = get_C(d, n, "own")
    cyc = load_cyc(d, n, C)
    model = build_model("chaboche", N)
    cal = calibrate(model, window(cyc, [2, 3]), cyc, (100.0, 600.0), {}, E_TRUE)
    p, s = tuple(cal["params"]), cal["sigma_y0"]
    return {"exp": "e1", "dia": d, "strain": n, "N": N, "C_mm_per_kN": C,
            "sigma_y0": s, "params": cal["params"],
            "calib_2_3": score(model, p, s, window(cyc, [2, 3]), E_TRUE),
            "val_4_6": score(model, p, s, window(cyc, [4, 5, 6]), E_TRUE),
            "all_2_6": score(model, p, s, window(cyc, [2, 3, 4, 5, 6]), E_TRUE),
            "val_4_6_weighted": score(model, p, s, window(cyc, [4, 5, 6]), E_TRUE, True)}


def task_e2(args):
    d, n, N, r, mode = args
    C = get_C(d, n, mode)
    tag = f"{mode}_N{N}_r{r:.2f}"
    mono = load_mono(d, n, C, tag)
    e_corr = data_loader.monotonic_modulus(mono["strain"], mono["stress"])
    cyc = load_cyc(d, n, C)
    model = build_model("chaboche", N)
    sy_hi = r * 500.0
    cal = calibrate(model, data_loader.slice_window(mono, 1), mono,
                    (0.995 * sy_hi, sy_hi), {"Q_inf": (0.0, 1e-6)}, E_TRUE)
    p, s = tuple(cal["params"]), cal["sigma_y0"]
    return {"exp": "e2", "cmode": mode, "dia": d, "strain": n, "N": N, "r": r,
            "C_mm_per_kN": C, "E_ramp_corrected_MPa": e_corr,
            "sigma_y0": s, "params": cal["params"], "ramp_fit": cal["fit_objective"],
            "val_4_6": score(model, p, s, window(cyc, [4, 5, 6]), E_TRUE),
            "all_2_6": score(model, p, s, window(cyc, [2, 3, 4, 5, 6]), E_TRUE)}


def task_e3(args):
    d, n = args
    C = get_C(d, n, "own")
    cyc = load_cyc(d, n, C)
    model = build_model("chaboche", 1)
    p = (1e-3, 1.0, 0.0, 1.0)                              # C1, gamma1, Q_inf, b
    best = None
    for s in range(200, 601, 20):
        v = score(model, p, float(s), window(cyc, [4, 5, 6]), E_TRUE)
        if best is None or v < best[0]:
            best = (v, s)
    s = best[1]
    return {"exp": "e3", "dia": d, "strain": n, "best_sy0": s, "val_4_6": best[0],
            "all_2_6": score(model, p, float(s), window(cyc, [2, 3, 4, 5, 6]), E_TRUE)}


def task_e0(args):
    d, n = args
    old = None
    for line in open(os.path.join(PROJ, "calibration_sessions", "mono_ablation_surrogate.jsonl")):
        r = json.loads(line)
        m = r["cyclic_file"].replace("\\", "/")
        if (m.endswith(f"{d}mm dia/LD5/{n}data.csv") and r["n_backstresses"] == 1
                and r["E_source"] == "ramp" and r.get("bauschinger") == 0.7):
            old = r
    model = build_model("chaboche", 1)
    p, s, E_old = tuple(old["params"]), old["sigma_y0"], old["E_modulus"]
    raw = load_cyc(d, n, 0.0)
    C = get_C(d, n, "own")
    cor = load_cyc(d, n, C)
    tgt = [2, 3, 4, 5, 6]
    out = {"exp": "e0", "dia": d, "strain": n, "C_mm_per_kN": C, "E_old": E_old}
    for wtd in (False, True):
        a = score(model, p, s, window(raw, tgt), E_old, wtd)
        b = score(model, p, s, window(cor, tgt), E_TRUE, wtd)
        out[f"old_weighted{int(wtd)}"] = a
        out[f"new_weighted{int(wtd)}"] = b
        out[f"absdiff_weighted{int(wtd)}"] = abs(a - b)
    return out


# ---------------------------------------------------------------- driver
def _warm():
    m = build_model("chaboche", 1)
    m.simulate(np.linspace(0, 0.01, 50), (1000.0, 10.0, 0.0, 1.0), E=E_TRUE, sy0=300.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True, choices=["e0", "e1", "e2", "e3"])
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--cmode", default="own", choices=["own", "global", "lodo"])
    ap.add_argument("--r", type=float, nargs="+", default=list(R_GRID))
    ap.add_argument("--n", type=int, nargs="+", default=[1, 2, 3])
    args = ap.parse_args()

    compliance_table()
    if args.exp == "e0":
        fn, tasks = task_e0, [(10, 4), (12, 5), (16, 3)]
    elif args.exp == "e1":
        fn, tasks = task_e1, [(d, n, N) for d, n in bars() for N in args.n]
    elif args.exp == "e2":
        fn = task_e2
        tasks = [(d, n, N, r, args.cmode) for d, n in bars() for N in args.n for r in args.r]
    else:
        fn, tasks = task_e3, bars()

    def key_of(t):
        if args.exp == "e1":
            return json.dumps({"exp": "e1", "dia": t[0], "strain": t[1], "N": t[2]}, sort_keys=True)
        if args.exp == "e2":
            return json.dumps({"exp": "e2", "cmode": t[4], "dia": t[0], "strain": t[1],
                               "N": t[2], "r": t[3]}, sort_keys=True)
        return json.dumps({"exp": args.exp, "dia": t[0], "strain": t[1]}, sort_keys=True)

    def key_of_rec(rec):
        if args.exp == "e1":
            return key_of((rec["dia"], rec["strain"], rec["N"]))
        if args.exp == "e2":
            return key_of((rec["dia"], rec["strain"], rec["N"], rec["r"], rec["cmode"]))
        return key_of((rec["dia"], rec["strain"]))

    name = args.exp + (f"_{args.cmode}" if args.exp == "e2" else "")
    if SEED != 42:
        name += f"_seed{SEED}"
    path = os.path.join(OUT, name + ".jsonl")
    done = set()
    if os.path.exists(path):
        for line in open(path):
            if line.strip():
                done.add(key_of_rec(json.loads(line)))

    todo = [t for t in tasks if key_of(t) not in done]
    print(f"{args.exp}: {len(tasks)} tasks, {len(todo)} to run, {args.workers} workers", flush=True)
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    n_done = 0
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_warm) as ex, \
            open(path, "a") as out:
        for rec in ex.map(fn, todo, chunksize=1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            n_done += 1
            if n_done % 10 == 0 or n_done == len(todo):
                print(f"  {n_done}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    print("done ->", path)


if __name__ == "__main__":
    main()
