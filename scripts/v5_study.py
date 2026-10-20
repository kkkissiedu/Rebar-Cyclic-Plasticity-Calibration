"""
v5_study.py
===========

Implements V5_PROTOCOL.md: standard Voce-Chaboche (N = 1, 2, 3 backstresses, E and sigma_y0 fitted) on the
UVC authors' extensometer coupon data, with the UVC error metric phi_bar (Eqs. 30-33).

    python v5_study.py --exp a [--workers 12]     # V5-A joint fit of the five protocols
    python v5_study.py --exp b [--workers 12]     # V5-B leave-one-protocol-out

Records go to calibration_sessions/v5/<exp>.jsonl (finished tasks are skipped on re-run).
"""

from __future__ import annotations

import argparse
import csv
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

from calibration_app.core.pipeline import build_model            # noqa: E402
import v3_study as v3                                             # noqa: E402  (optimiser stages)

DATA = os.environ.get("UVC_DATA", os.path.join(PROJ, "data", "Clean_Data"))
OUT = os.path.join(PROJ, "calibration_sessions", "v5")
SETS = {
    "A500": ("A500", "A500_t16"),
    "A992_web": ("A992_Gr50/A992_W14X82_web", "A992_t12"),
    "A992_flange": ("A992_Gr50/A992_W14X82_flange", "A992_t22"),
    "BCP325": ("BCP325", "BCP325_t22"),
    "BCR295": ("BCR295", "BCR295_t22"),
    "HYP400": ("HYP400", "HYP_t27"),
}
PROTOS = ["05inc", "2const", "2inc", "random", "tensile"]
PUBLISHED_VC = {"A992_web": 7.01, "A992_flange": 7.31, "A500": 6.86, "BCP325": 5.05, "BCR295": 8.38, "HYP400": 5.31}
PUBLISHED_UVC = {"A992_web": 5.21, "A992_flange": 4.76, "A500": 4.99, "BCP325": 3.89, "BCR295": 5.21, "HYP400": 4.67}
NS = (1, 2, 3)
SEEDS = (1, 2, 3)


# ---------------------------------------------------------------- data and metric
_CACHE = {}


def load_protocol(steel, proto):
    key = (steel, proto)
    if key in _CACHE:
        return _CACHE[key]
    folder, stem = SETS[steel]
    path = os.path.join(DATA, folder, f"{stem}_{proto}.csv")
    rows = []
    with open(path, newline="") as f:
        r = csv.reader(f)
        next(r)
        for row in r:
            if len(row) >= 2 and row[0].strip() != "" and row[1].strip() != "":
                rows.append((float(row[0]), float(row[1])))
    a = np.array(rows)
    sig, eps = a[:, 0], a[:, 1]
    dstar = np.abs(np.diff(eps))
    tot = float(dstar.sum())
    t = {"sig": sig, "eps": eps, "dstar": dstar, "tot": tot,
         "S": float(np.sum(0.5 * (sig[:-1] ** 2 + sig[1:] ** 2) * dstar) / tot)}
    _CACHE[key] = t
    return t


def mse(sim, t):
    e2 = (sim - t["sig"]) ** 2
    return float(np.sum(0.5 * (e2[:-1] + e2[1:]) * t["dstar"]) / t["tot"])


def evaluate_set(model, E, sy0, params, protos, steel):
    """Return (sum MSE, sum S, per-protocol phi_bar %) or (inf,...) on failure."""
    m_sum = s_sum = 0.0
    per = {}
    for p in protos:
        t = load_protocol(steel, p)
        sim = model.simulate(t["eps"], tuple(params), E=E, sy0=sy0)
        if sim is None or not np.all(np.isfinite(sim)):
            return float("inf"), 1.0, {}
        m = mse(sim, t)
        m_sum += m
        s_sum += t["S"]
        per[p] = 100.0 * float(np.sqrt(m / t["S"]))
    return m_sum, s_sum, per


# ---------------------------------------------------------------- parameter space
class Space:
    """E (linear), sigma_y0 (linear), C/gamma/b (log10), Q_inf (linear, >0); gamma-sorted."""

    def __init__(self, model, smax):
        self.model, self.names = model, list(model.param_names)
        self.lo = [100e3, 1.0]
        self.hi = [300e3, float(smax)]
        self.kind = ["lin", "lin"]
        for nme in self.names:
            if nme.startswith("gamma"):
                self.lo.append(-1.0); self.hi.append(4.0); self.kind.append("log")
            elif nme.startswith("C"):
                self.lo.append(0.0); self.hi.append(6.0); self.kind.append("log")
            elif nme == "b":
                self.lo.append(-2.0); self.hi.append(2.0); self.kind.append("log")
            elif nme == "Q_inf":
                self.lo.append(1e-4); self.hi.append(float(smax)); self.kind.append("lin")
        self.bounds = list(zip(self.lo, self.hi))
        self.search_names = ["E", "sigma_y0"] + self.names

    def decode(self, x):
        v = [10.0 ** xi if k == "log" else xi for xi, k in zip(x, self.kind)]
        E, sy0 = float(v[0]), float(v[1])
        vals = dict(zip(self.names, v[2:]))
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


def calibrate(steel, N, fit_protos, seed):
    model = build_model("chaboche", N)
    smax = max(float(np.max(np.abs(load_protocol(steel, p)["sig"]))) for p in fit_protos)
    sp = Space(model, smax)

    def f(x):
        E, sy0, params = sp.decode(x)
        m, _, _ = evaluate_set(model, E, sy0, params, fit_protos, steel)
        return m if np.isfinite(m) else 1e12

    x, fobj = v3.fit(f, sp, seed)
    E, sy0, params = sp.decode(x)
    return model, sp, x, E, sy0, params, fobj


def phibar(model, E, sy0, params, protos, steel):
    m, s, per = evaluate_set(model, E, sy0, params, protos, steel)
    return 100.0 * float(np.sqrt(m / s)) if np.isfinite(m) else float("nan"), per


# ---------------------------------------------------------------- tasks
def task_a(args):
    steel, N, seed = args
    model, sp, x, E, sy0, params, fobj = calibrate(steel, N, PROTOS, seed)
    pb, per = phibar(model, E, sy0, params, PROTOS, steel)
    gam = [params[2 * k + 1] for k in range(N)]
    return {"exp": "a", "steel": steel, "N": N, "seed": seed, "phibar": pb, "per_protocol": per,
            "E": E, "sigma_y0": sy0, "params": list(params), "boundary_hits": sp.boundary_hits(x),
            "gamma_collapse": bool(any(gam[k] / gam[k + 1] < 1.1 for k in range(N - 1)))}


def task_b(args):
    steel, N, seed, held = args
    fit_p = [p for p in PROTOS if p != held]
    model, sp, x, E, sy0, params, fobj = calibrate(steel, N, fit_p, seed)
    pb_in, _ = phibar(model, E, sy0, params, fit_p, steel)
    pb_out, _ = phibar(model, E, sy0, params, [held], steel)
    return {"exp": "b", "steel": steel, "N": N, "seed": seed, "held": held, "phibar_in": pb_in,
            "phibar_held": pb_out, "E": E, "sigma_y0": sy0, "params": list(params),
            "boundary_hits": sp.boundary_hits(x)}


# ---------------------------------------------------------------- driver
def _warm():
    m = build_model("chaboche", 1)
    m.simulate(np.linspace(0, 0.01, 50), (1000.0, 10.0, 10.0, 1.0), E=200000.0, sy0=300.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True, choices=["a", "b"], help="data folder: set UVC_DATA (default data/Clean_Data)")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--steels", nargs="+", default=list(SETS))
    ap.add_argument("--n", type=int, nargs="+", default=list(NS))
    ap.add_argument("--seeds", type=int, nargs="+", default=list(SEEDS))
    ap.add_argument("--time-one", action="store_true", help="run one task serially and print the time")
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    if args.exp == "a":
        fn = task_a
        tasks = [(s, N, sd) for s in args.steels for N in args.n for sd in args.seeds]
        keyf = lambda t: t
        keyr = lambda r: (r["steel"], r["N"], r["seed"])
    else:
        fn = task_b
        tasks = [(s, N, sd, h) for s in args.steels for N in args.n for sd in args.seeds for h in PROTOS]
        keyf = lambda t: t
        keyr = lambda r: (r["steel"], r["N"], r["seed"], r["held"])
    if args.time_one:
        t0 = time.time()
        r = fn(tasks[0])
        print("one task: %.1f s" % (time.time() - t0), {k: v for k, v in r.items() if k != "params"})
        return
    path = os.path.join(OUT, args.exp + ".jsonl")
    done = set()
    if os.path.exists(path):
        done = {keyr(json.loads(l)) for l in open(path) if l.strip()}
    todo = [t for t in tasks if keyf(t) not in done]
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
