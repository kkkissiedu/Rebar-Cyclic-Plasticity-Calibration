"""
v3_study.py
===========

Implements V3_PROTOCOL.md: the textbook Chaboche calibration, nothing tuned to the dataset.
Surrogate only, Chaboche only, unweighted objective, no physics gate, E measured from the
calibration cycles, strain used exactly as recorded.

    python v3_study.py --exp a  [--n 1 2 3] [--seeds 1 2 3] [--workers 12]   # V3-A
    python v3_study.py --exp r  [--seeds 1 2 3]                                # V3-R reference
    python v3_study.py --exp b                                                  # V3-B (needs V3-A)

Records go to calibration_sessions/v3/<exp>.jsonl (finished tasks are skipped on re-run).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)
sys.path.insert(0, HERE)

from calibration_app.core import objective                       # noqa: E402
from calibration_app.core import search as search_mod            # noqa: E402
from calibration_app.core.pipeline import build_model            # noqa: E402
import v2_study as v2                                             # noqa: E402  (loaders only)

OUT = os.path.join(PROJ, "calibration_sessions", "v3")
SOBOL_N, TPE_N, NM_N = 8192, 500, 2000
CALIB, VALID = [2, 3], [4, 5, 6]
QUALITY = 0.90
DROP_LO, DROP_HI = 0.05, 0.40
QREL = os.environ.get("V3_QREL") == "1"          # sensitivity run (V3-S), not the protocol run
SUFFIX = "_qrel" if QREL else ""


# ---------------------------------------------------------------- data helpers
def load(d, n):
    return v2.load_cyc(d, n, 0.0)                          # strain exactly as recorded


def good_cycles(data):
    """Cycles 2..6 with at least QUALITY x the median sample count of cycles 2..6."""
    cyc = np.asarray(data["cycle"])
    present = [c for c in data["cycles_present"] if 2 <= c <= 6]
    counts = {c: int(np.count_nonzero(cyc == c)) for c in present}
    if not counts:
        return []
    med = float(np.median(list(counts.values())))
    return [c for c in present if counts[c] >= QUALITY * med]


def window(data, target):
    ok = set(good_cycles(data))
    tgt = [c for c in target if c in ok]
    if not tgt:
        return None
    mask = np.asarray(data["cycle"]) <= max(tgt)
    return {"strain": data["strain"][mask], "stress": data["stress"][mask],
            "cycle": data["cycle"][mask], "target_cycles": tgt,
            "burnin_strain": data.get("burnin_strain", np.array([], float))}


def measure_E(data, cycles=CALIB):
    """Median LSQ slope on the unloading branch after each strain extremum (5-40 % stress drop)."""
    e, s, c = np.asarray(data["strain"]), np.asarray(data["stress"]), np.asarray(data["cycle"])
    slopes = []
    for cyc in cycles:
        m = c == cyc
        ee, ss = e[m], s[m]
        if ee.size < 20:
            continue
        for idx in (int(np.argmax(ee)), int(np.argmin(ee))):
            s0 = ss[idx]
            if abs(s0) < 1.0:
                continue
            seg_e, seg_s = ee[idx:], ss[idx:]
            drop = np.sign(s0) * (s0 - seg_s) / abs(s0)       # >= 0 while unloading
            over = np.nonzero(drop > DROP_HI)[0]
            stop = int(over[0]) if over.size else len(drop)
            sel = np.arange(stop)[(drop[:stop] >= DROP_LO) & (drop[:stop] <= DROP_HI)]
            if sel.size < 4:
                continue
            slope = float(np.polyfit(seg_e[sel], seg_s[sel], 1)[0])
            if slope > 0:
                slopes.append(slope)
    return float(np.median(slopes)) if len(slopes) >= 2 else None


def score(model, params, sy0, win, E):
    if win is None:
        return float("nan")
    return float(objective.evaluate_surrogate(
        model, tuple(params), win, None, E=E, sy0=float(sy0), weight_reversals=False))


# ---------------------------------------------------------------- parameter space
class Space:
    """Generic box (protocol): sigma0 lin, C/gamma/b log10, Q_inf lin; gamma-sorted."""

    def __init__(self, model, S):
        self.model, self.names = model, list(model.param_names)
        self.lo, self.hi, self.kind = [1.0], [float(S)], ["lin"]
        for nme in self.names:
            if nme.startswith("gamma"):
                self.lo.append(-1.0); self.hi.append(4.0); self.kind.append("log")
            elif nme.startswith("C"):
                self.lo.append(0.0); self.hi.append(6.0); self.kind.append("log")
            elif nme == "b":
                self.lo.append(-2.0); self.hi.append(2.0); self.kind.append("log")
            elif nme == "Q_inf" and QREL:                  # sensitivity run only (V3-S)
                self.lo.append(-1.0); self.hi.append(1.0); self.kind.append("qrel")
            elif nme == "Q_inf":
                self.lo.append(-0.5 * S); self.hi.append(float(S)); self.kind.append("lin")
            else:
                raise ValueError(nme)
        self.bounds = list(zip(self.lo, self.hi))
        self.search_names = ["sigma_y0"] + self.names

    def decode(self, x):
        v = [10.0 ** xi if k == "log" else xi for xi, k in zip(x, self.kind)]
        sy0, vals = float(v[0]), dict(zip(self.names, v[1:]))
        if QREL:                                           # Q_inf = q * sigma0 (yield radius >= 0)
            vals["Q_inf"] = vals["Q_inf"] * sy0
        n_bs = sum(1 for nme in self.names if nme.startswith("gamma"))
        pairs = sorted(((vals[f"C{k}"], vals[f"gamma{k}"]) for k in range(1, n_bs + 1)),
                       key=lambda p: -p[1])
        out = []
        for c_, g_ in pairs:
            out += [c_, g_]
        out += [vals["Q_inf"], vals["b"]]
        return sy0, tuple(out)

    def boundary_hits(self, x):
        return [self.search_names[i] for i, (xi, lo, hi) in enumerate(zip(x, self.lo, self.hi))
                if (xi - lo) / (hi - lo) < 0.01 or (xi - lo) / (hi - lo) > 0.99]


def fit(evaluate, space, seed):
    """Sobol -> TPE -> Nelder-Mead polish (protocol). Returns best x and objective."""
    from scipy.optimize import minimize
    seeds, s1 = search_mod.sobol_seed(evaluate, space.bounds, n_samples=SOBOL_N,
                                      top_k=10, seed=seed)
    s2 = search_mod.run_bayesian(evaluate, space.bounds, space.search_names,
                                 n_trials=TPE_N, x0=seeds[0] if seeds else None, seed=seed)
    best = s1 if s1.best_objective <= s2.best_objective else s2
    x0, f0 = np.array(best.best_params, float), float(best.best_objective)
    lo, hi = np.array(space.lo), np.array(space.hi)

    def clipped(x):
        return evaluate(tuple(np.clip(x, lo, hi)))

    res = minimize(clipped, x0, method="Nelder-Mead",
                   options={"maxfev": NM_N, "xatol": 1e-6, "fatol": 1e-6})
    if res.fun < f0:
        x0, f0 = np.clip(res.x, lo, hi), float(res.fun)
    return tuple(float(v) for v in x0), f0


# ---------------------------------------------------------------- tasks
def _calib_stress(wc):
    sel = np.isin(wc["cycle"], wc["target_cycles"])
    return float(np.max(np.abs(wc["stress"][sel])))


def task_a(args):
    d, n, N, seed = args
    data = load(d, n)
    model = build_model("chaboche", N)
    wc, wv = window(data, CALIB), window(data, VALID)
    E = measure_E(data)
    S = _calib_stress(wc)
    sp = Space(model, S)

    def evaluate(x):
        sy0, p = sp.decode(x)
        return objective.evaluate_surrogate(model, p, wc, None, E=E, sy0=sy0,
                                            weight_reversals=False)

    x, f = fit(evaluate, sp, seed)
    sy0, p = sp.decode(x)
    gam = [p[2 * k + 1] for k in range(N)]
    collapse = any(gam[k] / gam[k + 1] < 1.1 for k in range(N - 1))
    return {"exp": "a", "dia": d, "strain": n, "N": N, "seed": seed, "E_meas": E, "S": S,
            "sigma_y0": sy0, "params": list(p), "calib": score(model, p, sy0, wc, E),
            "val": score(model, p, sy0, wv, E),
            "n_val_cycles": 0 if wv is None else len(wv["target_cycles"]),
            "boundary_hits": sp.boundary_hits(x), "gamma_collapse": bool(collapse)}


def task_r(args):
    d, n, seed = args
    data = load(d, n)
    model = build_model("chaboche", 1)
    wc, wv = window(data, CALIB), window(data, VALID)
    E = measure_E(data)
    S = _calib_stress(wc)
    fixed = (1e-6, 1.0, 0.0, 1.0)                           # C1~0, gamma1, Q_inf=0, b

    class OneD:
        bounds = [(1.0, S)]
        lo, hi = [1.0], [S]
        search_names = ["sigma_y0"]

    def evaluate(x):
        return objective.evaluate_surrogate(model, fixed, wc, None, E=E, sy0=float(x[0]),
                                            weight_reversals=False)

    x, f = fit(evaluate, OneD, seed)
    return {"exp": "r", "dia": d, "strain": n, "seed": seed, "E_meas": E, "sigma_y0": x[0],
            "calib": score(model, fixed, x[0], wc, E), "val": score(model, fixed, x[0], wv, E),
            "n_val_cycles": 0 if wv is None else len(wv["target_cycles"])}


def run_b():
    """Cross-amplitude: parameters from bar a applied to bar b (same diameter), cycles 2-6 of b."""
    path = os.path.join(OUT, "a.jsonl")
    recs = [json.loads(l) for l in open(path) if l.strip()]
    data = {(d, n): load(d, n) for d, n in v2.bars()}
    Es = {k: measure_E(v) for k, v in data.items()}
    out = []
    for r in recs:
        model = build_model("chaboche", r["N"])
        for n2 in v2.STRAINS:
            if n2 == r["strain"]:
                continue
            k = (r["dia"], n2)
            w = window(data[k], [2, 3, 4, 5, 6])
            out.append({"exp": "b", "dia": r["dia"], "cal_strain": r["strain"], "test_strain": n2,
                        "N": r["N"], "seed": r["seed"],
                        "err": score(model, r["params"], r["sigma_y0"], w, Es[k])})
    with open(os.path.join(OUT, "b.jsonl"), "w") as f:
        for o in out:
            f.write(json.dumps(o) + "\n")
    print(f"b: {len(out)} records")


# ---------------------------------------------------------------- driver
def _warm():
    m = build_model("chaboche", 1)
    m.simulate(np.linspace(0, 0.01, 50), (1000.0, 10.0, 0.0, 1.0), E=60000.0, sy0=300.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True, choices=["a", "r", "b"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--n", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if args.exp == "b":
        run_b()
        return
    if args.exp == "a":
        fn = task_a
        tasks = [(d, n, N, s) for d, n in v2.bars() for N in args.n for s in args.seeds]
        keyf = lambda t: (t[0], t[1], t[2], t[3])
        keyr = lambda r: (r["dia"], r["strain"], r["N"], r["seed"])
    else:
        fn = task_r
        tasks = [(d, n, s) for d, n in v2.bars() for s in args.seeds]
        keyf = lambda t: (t[0], t[1], t[2])
        keyr = lambda r: (r["dia"], r["strain"], r["seed"])
    path = os.path.join(OUT, args.exp + SUFFIX + ".jsonl")
    done = set()
    if os.path.exists(path):
        done = {keyr(json.loads(l)) for l in open(path) if l.strip()}
    todo = [t for t in tasks if keyf(t) not in done]
    print(f"{args.exp}: {len(tasks)} tasks, {len(todo)} to run, {args.workers} workers", flush=True)
    t0, nd = time.time(), 0
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_warm) as ex, \
            open(path, "a") as out:
        for rec in ex.map(fn, todo, chunksize=1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            nd += 1
            if nd % 6 == 0 or nd == len(todo):
                print(f"  {nd}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    print("done ->", path)


if __name__ == "__main__":
    main()
