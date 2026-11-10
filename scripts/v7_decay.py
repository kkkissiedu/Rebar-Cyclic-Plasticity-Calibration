"""
v7_decay.py
===========

V7 degradation layer (L2) with the full cross-validation of V7_PROTOCOL.md, development set only.

Peak tension and compression stress of cycle n, normalised by the cycle-2 value, follow a master curve g(x; lambda_p, eps_a) of the life
fraction x = n / Nf. Three forms are compared (inside each training set only, by 5-fold cross-validation over files):
  power_drop  : g = 1 - A x^p - B x^m
  weibull_drop: g = exp(-A x^p) - B x^m
  exp_drop    : g = 1 - A (1 - exp(-a x)) - B (exp(c x) - 1) / (exp(c) - 1)
A and B are linear in (1, L/D / 10, amplitude % / 6); the two shape parameters are global. Tension and compression are fitted separately.
The cycle-2 peaks come from a log-linear model of the amplitude (with or without lambda_p and D, chosen by leave-one-out), and the life from
the L1 layer (v7_life). Reported in two modes: A = decay only (observed cycle-2 peaks and observed Nf) and B = the full chain (predicted
cycle-2 peaks and predicted Nf). Baseline B1 repeats the cycle-2 peak (no degradation). Error = |predicted - observed| / max(observed,
0.1 * observed cycle-2 peak): the floor stops the near-zero fracture cycle from dominating (an implementation detail of "relative error").

    python v7_decay.py            # everything (several minutes)
    python v7_decay.py --quick    # three folds, to time it
"""

from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
from scipy.optimize import least_squares

import v7_data as D
import v7_life as L

MAXPTS = 40
FLOOR = 0.10
VARIANTS = ("power_drop", "weibull_drop", "exp_drop")
P2_VARIANTS = ("eps", "eps_lam", "eps_lam_d")
BOUNDS = {"power_drop": ((0.05, 3.0), (3.0, 60.0)), "weibull_drop": ((0.05, 3.0), (3.0, 60.0)), "exp_drop": ((0.3, 30.0), (2.0, 40.0))}
START = {"power_drop": (0.5, 10.0), "weibull_drop": (0.5, 10.0), "exp_drop": (5.0, 10.0)}


def make_files(feats: list = None) -> list:
    out = []
    for f in (feats if feats is not None else D.build_dev_features()):
        cyc = {c["c"]: c for c in f["cycles"]}
        if 2 not in cyc:
            continue
        cs = [c for c in range(2, f["nf"] + 1) if c in cyc]
        t = np.array([cyc[c]["t"] for c in cs]); k = np.array([cyc[c]["k"] for c in cs])
        c = np.array(cs, float)
        idx = np.unique(np.round(np.linspace(0, len(c) - 1, min(MAXPTS, len(c)))).astype(int))
        out.append({"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": f["nf"], "c": c, "t": t, "k": k,
                    "t2": cyc[2]["t"], "k2": cyc[2]["k"], "cov": np.array([1.0, f["ld"] / 10.0, f["n"] / 6.0]),
                    "xs": c[idx] / f["nf"], "yt": t[idx] / cyc[2]["t"], "yk": k[idx] / cyc[2]["k"]})
    return out


def gcurve(x, a, b, s1, s2, v):
    x = np.clip(x, 0.0, 1.0)
    if v == "power_drop":
        return 1.0 - a * x ** s1 - b * x ** s2
    if v == "weibull_drop":
        return np.exp(-a * x ** s1) - b * x ** s2
    return 1.0 - a * (1.0 - np.exp(-s1 * x)) - b * (np.exp(s2 * x) - 1.0) / (np.exp(s2) - 1.0)


def fit_decay(files: list, key: str, v: str) -> np.ndarray:
    (l1, u1), (l2, u2) = BOUNDS[v]
    th0 = np.array([0.3, 0, 0, 0.3, 0, 0, *START[v]])
    w = [1.0 / np.sqrt(len(f["xs"])) for f in files]

    def res(th):
        return np.concatenate([(gcurve(f["xs"], f["cov"] @ th[:3], f["cov"] @ th[3:6], th[6], th[7], v) - f[key]) * wi
                               for f, wi in zip(files, w)])

    lb = np.array([-5.0] * 6 + [l1, l2]); ub = np.array([5.0] * 6 + [u1, u2])
    return least_squares(res, th0, bounds=(lb, ub), max_nfev=150).x


def pred_g(f: dict, th: np.ndarray, v: str, nf_used: float) -> np.ndarray:
    x = np.minimum(f["c"] / max(nf_used, 2.0), 1.0)
    return np.clip(gcurve(x, f["cov"] @ th[:3], f["cov"] @ th[3:6], th[6], th[7], v), 0.0, 1.5)


def rel_err(pred, obs, ref):
    return np.abs(pred - obs) / np.maximum(obs, FLOOR * ref)


def p2_design(files, variant):
    cols = [np.ones(len(files)), np.log([f["n"] for f in files])]
    if variant in ("eps_lam", "eps_lam_d"):
        cols.append(np.array([f["ld"] for f in files]) * L.LAMBDA)
    if variant == "eps_lam_d":
        cols.append(np.array([f["d"] for f in files], float))
    return np.column_stack(cols)


def fit_p2(files, variant):
    x = p2_design(files, variant)
    return (np.linalg.lstsq(x, np.log([f["t2"] for f in files]), rcond=None)[0],
            np.linalg.lstsq(x, np.log([f["k2"] for f in files]), rcond=None)[0])


def pred_p2(coef, files, variant):
    x = p2_design(files, variant)
    return np.exp(x @ coef[0]), np.exp(x @ coef[1])


def select_p2(train):
    best, best_rms = None, np.inf
    for v in P2_VARIANTS:
        e = []
        for i in range(len(train)):
            rest = train[:i] + train[i + 1:]
            c = fit_p2(rest, v)
            pt, pk = pred_p2(c, [train[i]], v)
            e += [np.log(pt[0] / train[i]["t2"]), np.log(pk[0] / train[i]["k2"])]
        rms = float(np.sqrt(np.mean(np.square(e))))
        if rms < best_rms:
            best, best_rms = v, rms
    return best


def select_decay(train):
    """Variant with the smallest decay-only median error over a 5-fold split of the training files."""
    best, best_err = None, np.inf
    for v in VARIANTS:
        errs = []
        for q in range(5):
            tr = [f for i, f in enumerate(train) if i % 5 != q]
            va = [f for i, f in enumerate(train) if i % 5 == q]
            if not va:
                continue
            tt, tk = fit_decay(tr, "yt", v), fit_decay(tr, "yk", v)
            for f in va:
                errs.append(rel_err(f["t2"] * pred_g(f, tt, v, f["nf"]), f["t"], f["t2"]))
                errs.append(rel_err(f["k2"] * pred_g(f, tk, v, f["nf"]), f["k"], f["k2"]))
        m = float(np.median(np.concatenate(errs)))
        if m < best_err:
            best, best_err = v, m
    return best


def run_fold(train: list, test: list) -> dict:
    v = select_decay(train)
    th_t, th_k = fit_decay(train, "yt", v), fit_decay(train, "yk", v)
    pv = select_p2(train)
    p2c = fit_p2(train, pv)
    rtr = np.array([(f["d"], f["ld"], f["n"], f["nf"]) for f in train], float)
    lv = L.select_variant(rtr)
    lc = L.fit(rtr, lv)
    out = {"variant": v, "A_fw": [], "A_b1": [], "B_fw": [], "B_b1": []}
    for f in test:
        rt = np.array([[f["d"], f["ld"], f["n"], f["nf"]]], float)
        nfp = float(L.predict(lc, rt, lv)[0])
        t2p, k2p = (a[0] for a in pred_p2(p2c, [f], pv))
        for tag, t2, k2, nfu in (("A", f["t2"], f["k2"], f["nf"]), ("B", t2p, k2p, nfp)):
            for obs, ref, base, th in ((f["t"], f["t2"], t2, th_t), (f["k"], f["k2"], k2, th_k)):
                out[tag + "_fw"].append(rel_err(base * pred_g(f, th, v, nfu), obs, ref))
                out[tag + "_b1"].append(rel_err(np.full_like(obs, base), obs, ref))
    return {k: (np.concatenate(x) if isinstance(x, list) else x) for k, x in out.items()}


def folds(files: list) -> dict:
    r = np.array([(f["d"], f["ld"], f["n"], f["nf"]) for f in files], float)
    return L.folds(r)


def main() -> int:
    quick = "--quick" in sys.argv
    files = make_files()
    print(f"development files: {len(files)}")
    groups = folds(files)
    results = {}
    t0 = time.time()
    for group, fl in groups.items():
        if quick and group != "leave one condition out":
            continue
        pool = {k: [] for k in ("A_fw", "A_b1", "B_fw", "B_b1")}
        win = {"A": [], "B": []}
        chosen = {}
        for name, test in (fl[:3] if quick else fl):
            test = [int(i) for i in test]
            train = [f for i, f in enumerate(files) if i not in set(test)]
            res = run_fold(train, [files[i] for i in test])
            chosen[res["variant"]] = chosen.get(res["variant"], 0) + 1
            for k in pool:
                pool[k].append(res[k])
            for tag in ("A", "B"):
                win[tag].append(float(np.median(res[tag + "_fw"]) < np.median(res[tag + "_b1"])))
        results[group] = {"chosen": chosen, "folds": len(win["A"])}
        for tag in ("A", "B"):
            fw, b1 = np.concatenate(pool[tag + "_fw"]), np.concatenate(pool[tag + "_b1"])
            results[group][tag] = {"framework_median_err": float(np.median(fw)), "no_degradation_median_err": float(np.median(b1)),
                                   "share_of_folds_better": float(np.mean(win[tag]))}
        print(f"[{time.time() - t0:6.0f}s] {group:44s} A(decay only) fw {results[group]['A']['framework_median_err']:.3f} vs B1 "
              f"{results[group]['A']['no_degradation_median_err']:.3f} win {results[group]['A']['share_of_folds_better']:.2f} | "
              f"B(full chain) fw {results[group]['B']['framework_median_err']:.3f} vs B1 {results[group]['B']['no_degradation_median_err']:.3f} "
              f"win {results[group]['B']['share_of_folds_better']:.2f} | chosen {chosen}", flush=True)
    if not quick:
        os.makedirs(D.OUT, exist_ok=True)
        with open(os.path.join(D.OUT, "decay_cv.json"), "w") as f:
            json.dump(results, f, indent=1)
        sel = [results[g] for g in ("leave one condition out", "leave one amplitude out")]
        n = [s["folds"] for s in sel]
        share = sum(s["B"]["share_of_folds_better"] * s["folds"] for s in sel) / sum(n)
        print(f"\nACCEPTANCE (decay part, full chain): framework beats B1 on {share:.0%} of leave-one-condition-out and leave-one-amplitude-out "
              f"folds (needs >= 80 %) -> {'PASS' if share >= 0.8 else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
