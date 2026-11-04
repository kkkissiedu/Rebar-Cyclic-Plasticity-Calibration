"""
v6_study.py
===========

Implements V6_PROTOCOL.md: rebar re-calibration with cycle 1 in the objective, true stress/strain,
fitted E, non-softening (or q*sigma0) isotropic term, refined strain path and the UVC metric.

    python v6_study.py --exp a [--workers 12]     # V6-A calibrate cycles 1-3, validate 4-6
    python v6_study.py --exp r [--workers 12]     # elastic-perfectly-plastic reference
    python v6_study.py --exp b                    # V6-B cross-amplitude (needs V6-A)

Records go to calibration_sessions/v6/<exp>.jsonl (finished tasks are skipped on re-run).
"""

from __future__ import annotations

import argparse
import csv
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

from calibration_app.core.pipeline import build_model            # noqa: E402
import v3_study as v3                                             # noqa: E402  (optimiser stages, E measurement)
import v2_study as v2                                             # noqa: E402  (bar list)

OUT = os.path.join(PROJ, "calibration_sessions", "v6")
E_BAR = 200000.0
MAX_INC = 1e-4
CAL_CYCLES = 3
VAL_CYCLES = [4, 5, 6]
QUALITY = 0.90
VARIANTS = ("nonsoft", "soft")
NS = (1, 2, 3)
SEEDS = (1, 2, 3)


# ---------------------------------------------------------------- data
_CACHE = {}


def load_bar(d, n):
    key = (d, n)
    if key in _CACHE:
        return _CACHE[key]
    path = os.path.join(PROJ, "input_data", f"{d}mm dia", "LD5", f"{n}data.csv")
    cyc, pos, ld = [], [], []
    with open(path, newline="") as f:
        r = csv.reader(f)
        hdr = None
        for raw in r:
            if raw and raw[0].strip().strip('"').lower().startswith("notes"):
                hdr = [c.strip().strip('"').lower() for c in raw]
                break
        ic, ip, il = hdr.index("cycle"), hdr.index("position mm"), hdr.index("load kn")
        for raw in r:
            try:
                cyc.append(int(float(raw[ic]))); pos.append(float(raw[ip])); ld.append(float(raw[il]))
            except (ValueError, IndexError):
                continue
    cyc, pos, ld = np.array(cyc), np.array(pos), np.array(ld)
    L0, A0 = 5.0 * d, math.pi * d * d / 4.0
    eps_eng = (pos - pos[0]) / L0
    sig_eng = (ld - ld[0]) * 1000.0 / A0
    eps_t = np.log1p(eps_eng)
    sig_t = sig_eng * (1.0 + eps_eng)
    counts = {c: int(np.count_nonzero(cyc == c)) for c in np.unique(cyc) if 2 <= c <= 6}
    med = float(np.median(list(counts.values()))) if counts else 0.0
    good = [c for c, k in counts.items() if k >= QUALITY * med]
    val_cycles = [c for c in VAL_CYCLES if c in good]
    out = {"cyc": cyc, "eps_eng": eps_eng, "sig_eng": sig_eng, "eps": eps_t, "sig": sig_t,
           "L0": L0, "A0": A0, "val_cycles": val_cycles}
    _CACHE[key] = out
    return out


def refine(eps):
    """Linear-interpolation refinement to increments <= MAX_INC; returns (fine path, index of original samples)."""
    k = np.maximum(1, np.ceil(np.abs(np.diff(eps)) / MAX_INC).astype(int))
    t = np.concatenate([[0], np.cumsum(k)])
    fine = np.interp(np.arange(t[-1] + 1), t, eps)
    return fine, t


class Prefix:
    """Strain path through the last sample of cycle <= last_cycle, with its refinement."""

    def __init__(self, bar, last_cycle):
        m = int(np.nonzero(bar["cyc"] <= last_cycle)[0][-1]) + 1
        self.m = m
        self.eps = bar["eps"][:m]
        self.sig = bar["sig"][:m]
        self.eps_eng = bar["eps_eng"][:m]
        self.sig_eng = bar["sig_eng"][:m]
        self.cyc = bar["cyc"][:m]
        self.fine, self.idx = refine(self.eps)
        self.dstar = np.abs(np.diff(self.eps))

    def simulate(self, model, E, sy0, params):
        s = model.simulate(self.fine, tuple(params), E=E, sy0=sy0)
        if s is None or not np.all(np.isfinite(s)):
            return None
        return s[self.idx]


def metric(sim, sig, dstar, mask):
    """MSE and S over the samples in mask (accumulated-strain weighted, trapezoid)."""
    e2 = (sim - sig) ** 2
    s2 = sig ** 2
    pair = mask[:-1] & mask[1:]
    w = dstar * pair
    tot = float(w.sum())
    if tot <= 0:
        return float("nan"), float("nan")
    mse = float(np.sum(0.5 * (e2[:-1] + e2[1:]) * w) / tot)
    S = float(np.sum(0.5 * (s2[:-1] + s2[1:]) * w) / tot)
    return mse, S


# ---------------------------------------------------------------- parameter space
class Space:
    def __init__(self, model, smax, variant):
        self.model, self.names, self.variant = model, list(model.param_names), variant
        self.lo, self.hi, self.kind = [20e3, 1.0], [E_BAR, float(smax)], ["lin", "lin"]
        for nme in self.names:
            if nme.startswith("gamma"):
                self.lo.append(-1.0); self.hi.append(4.0); self.kind.append("log")
            elif nme.startswith("C"):
                self.lo.append(0.0); self.hi.append(6.0); self.kind.append("log")
            elif nme == "b":
                self.lo.append(-2.0); self.hi.append(2.0); self.kind.append("log")
            elif nme == "Q_inf" and variant == "soft":
                self.lo.append(-1.0); self.hi.append(1.0); self.kind.append("qrel")
            elif nme == "Q_inf":
                self.lo.append(1e-4); self.hi.append(float(smax)); self.kind.append("lin")
        self.bounds = list(zip(self.lo, self.hi))
        self.search_names = ["E_fit", "sigma_y0"] + self.names

    def decode(self, x):
        v = [10.0 ** xi if k == "log" else xi for xi, k in zip(x, self.kind)]
        E, sy0 = float(v[0]), float(v[1])
        vals = dict(zip(self.names, v[2:]))
        if self.variant == "soft":
            vals["Q_inf"] = vals["Q_inf"] * sy0
        n_bs = sum(1 for nme in self.names if nme.startswith("gamma"))
        pairs = sorted(((vals[f"C{k}"], vals[f"gamma{k}"]) for k in range(1, n_bs + 1)), key=lambda p: -p[1])
        out = []
        for c_, g_ in pairs:
            out += [c_, g_]
        out += [vals["Q_inf"], vals["b"]]
        return E, sy0, tuple(out)

    def boundary_hits(self, x):
        return [self.search_names[i] for i, (xi, lo, hi) in enumerate(zip(x, self.lo, self.hi))
                if (xi - lo) / (hi - lo) < 0.01 or (xi - lo) / (hi - lo) > 0.99]


class EPPSpace:
    """Elastic-perfectly-plastic reference: only E and sigma_y0 are free."""

    def __init__(self, smax):
        self.lo, self.hi = [20e3, 1.0], [E_BAR, float(smax)]
        self.bounds = list(zip(self.lo, self.hi))
        self.search_names = ["E_fit", "sigma_y0"]


def first_tension_peak(eps):
    """Index of the first tension strain maximum (before strain falls 0.2 % below the running maximum)."""
    run = np.maximum.accumulate(eps)
    drop = np.nonzero(run - eps > 0.002)[0]
    stop = int(drop[0]) if drop.size else len(eps)
    return int(np.argmax(eps[:stop]))


def evaluate_bar(model, E, sy0, params, d, n, pf_cal, pf_val):
    bar = load_bar(d, n)
    out = {}
    sim_c = pf_cal.simulate(model, E, sy0, params)
    m_all = np.ones(pf_cal.m, bool)
    mse, S = metric(sim_c, pf_cal.sig, pf_cal.dstar, m_all)
    out["cal_phibar"] = 100.0 * math.sqrt(mse / S)
    m1 = pf_cal.cyc == 1
    mse1, S1 = metric(sim_c, pf_cal.sig, pf_cal.dstar, m1)
    out["cal_phibar_cycle1"] = 100.0 * math.sqrt(mse1 / S1)
    ip = first_tension_peak(pf_cal.eps)
    out["virgin_peak_meas"] = float(pf_cal.sig_eng[ip])
    out["virgin_peak_model"] = float(sim_c[ip] / (1.0 + pf_cal.eps_eng[ip]))
    out["virgin_peak_strain"] = float(pf_cal.eps_eng[ip])
    vc = bar["val_cycles"]
    if vc and pf_val is not None:
        sim_v = pf_val.simulate(model, E, sy0, params)
        mask = np.isin(pf_val.cyc, vc)
        mv, Sv = metric(sim_v, pf_val.sig, pf_val.dstar, mask)
        out["val_phibar"] = 100.0 * math.sqrt(mv / Sv) if np.isfinite(mv) else float("nan")
        eng = sim_v / (1.0 + pf_val.eps_eng)
        r = [math.sqrt(float(np.mean((eng[pf_val.cyc == c] - pf_val.sig_eng[pf_val.cyc == c]) ** 2))) for c in vc]
        out["val_rmse_eng"] = float(np.mean(r))
    else:
        out["val_phibar"] = out["val_rmse_eng"] = float("nan")
    out["n_val_cycles"] = len(vc)
    return out


# ---------------------------------------------------------------- tasks
def task_a(args):
    d, n, N, variant, seed = args
    bar = load_bar(d, n)
    pf_cal = Prefix(bar, CAL_CYCLES)
    pf_val = Prefix(bar, max(bar["val_cycles"])) if bar["val_cycles"] else None
    model = build_model("chaboche", N)
    sp = Space(model, float(np.max(np.abs(pf_cal.sig))), variant)
    all_mask = np.ones(pf_cal.m, bool)

    def f(x):
        E, sy0, p = sp.decode(x)
        s = pf_cal.simulate(model, E, sy0, p)
        if s is None:
            return 1e12
        mse, _ = metric(s, pf_cal.sig, pf_cal.dstar, all_mask)
        return mse if np.isfinite(mse) else 1e12

    x, fobj = v3.fit(f, sp, seed)
    E, sy0, p = sp.decode(x)
    rec = {"exp": "a", "dia": d, "strain": n, "N": N, "variant": variant, "seed": seed,
           "E_fit": E, "sigma_y0": sy0, "params": list(p), "boundary_hits": sp.boundary_hits(x)}
    rec.update(evaluate_bar(model, E, sy0, p, d, n, pf_cal, pf_val))
    gam = [p[2 * k + 1] for k in range(N)]
    rec["gamma_collapse"] = bool(any(gam[k] / gam[k + 1] < 1.1 for k in range(N - 1)))
    return rec


def task_r(args):
    d, n, seed = args
    bar = load_bar(d, n)
    pf_cal = Prefix(bar, CAL_CYCLES)
    pf_val = Prefix(bar, max(bar["val_cycles"])) if bar["val_cycles"] else None
    model = build_model("chaboche", 1)
    sp = EPPSpace(float(np.max(np.abs(pf_cal.sig))))
    fixed = (1e-6, 1.0, 0.0, 1.0)
    all_mask = np.ones(pf_cal.m, bool)

    def f(x):
        s = pf_cal.simulate(model, float(x[0]), float(x[1]), fixed)
        if s is None:
            return 1e12
        mse, _ = metric(s, pf_cal.sig, pf_cal.dstar, all_mask)
        return mse if np.isfinite(mse) else 1e12

    x, fobj = v3.fit(f, sp, seed)
    rec = {"exp": "r", "dia": d, "strain": n, "seed": seed, "E_fit": float(x[0]), "sigma_y0": float(x[1])}
    rec.update(evaluate_bar(model, float(x[0]), float(x[1]), fixed, d, n, pf_cal, pf_val))
    return rec


def run_b():
    """Cross-amplitude: parameters of each diameter's 3 % bar applied to its 4, 5, 6 % bars."""
    recs = [json.loads(l) for l in open(os.path.join(OUT, "a.jsonl")) if l.strip()]
    out = []
    for r in recs:
        if r["strain"] != 3:
            continue
        model = build_model("chaboche", r["N"])
        for n2 in (4, 5, 6):
            bar = load_bar(r["dia"], n2)
            E_test = v3.measure_E(v3.load(r["dia"], n2))
            pf_cal = Prefix(bar, CAL_CYCLES)
            pf_val = Prefix(bar, max(bar["val_cycles"])) if bar["val_cycles"] else None
            ev = evaluate_bar(model, E_test, r["sigma_y0"], r["params"], r["dia"], n2, pf_cal, pf_val)
            ev.update({"exp": "b", "dia": r["dia"], "cal_strain": 3, "test_strain": n2, "N": r["N"],
                       "variant": r["variant"], "seed": r["seed"], "E_test": E_test})
            out.append(ev)
    with open(os.path.join(OUT, "b.jsonl"), "w") as f:
        for o in out:
            f.write(json.dumps(o) + "\n")
    print(f"b: {len(out)} records")


# ---------------------------------------------------------------- driver
def _warm():
    m = build_model("chaboche", 1)
    m.simulate(np.linspace(0, 0.01, 50), (1000.0, 10.0, 10.0, 1.0), E=100000.0, sy0=300.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True, choices=["a", "r", "b"])
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--n", type=int, nargs="+", default=list(NS))
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS))
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    ap.add_argument("--bars", type=int, nargs="+", help="optional: restrict to these diameters")
    ap.add_argument("--time-one", action="store_true")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if args.exp == "b":
        run_b()
        return
    bars = [(d, n) for d, n in v2.bars() if (not args.bars or d in args.bars)]
    if args.exp == "a":
        fn = task_a
        tasks = [(d, n, N, v, s) for d, n in bars for N in args.n for v in args.variants for s in args.seeds]
        keyr = lambda r: (r["dia"], r["strain"], r["N"], r["variant"], r["seed"])
    else:
        fn = task_r
        tasks = [(d, n, s) for d, n in bars for s in args.seeds]
        keyr = lambda r: (r["dia"], r["strain"], r["seed"])
    if args.time_one:
        t0 = time.time()
        r = fn(tasks[0])
        print("one task: %.1f s" % (time.time() - t0), {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != "params"})
        return
    path = os.path.join(OUT, args.exp + ".jsonl")
    done = set()
    if os.path.exists(path):
        done = {keyr(json.loads(l)) for l in open(path) if l.strip()}
    todo = [t for t in tasks if t not in done]
    print(f"{args.exp}: {len(tasks)} tasks, {len(todo)} to run, {args.workers} workers", flush=True)
    t0, nd = time.time(), 0
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_warm) as ex, open(path, "a") as out:
        for rec in ex.map(fn, todo, chunksize=1):
            out.write(json.dumps(rec) + "\n")
            out.flush()
            nd += 1
            if nd % 6 == 0 or nd == len(todo):
                print(f"  {nd}/{len(todo)}  {time.time() - t0:.0f}s", flush=True)
    print("done ->", path)


if __name__ == "__main__":
    main()
