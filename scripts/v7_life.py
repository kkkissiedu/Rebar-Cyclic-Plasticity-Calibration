"""
v7_life.py
==========

V7 life layer (L1) and its full cross-validation, exactly as frozen in V7_PROTOCOL.md (sections 3 to 5). Development set only.

Model (Koh-Stephens form as in Kashani 2019): eps_a = eps_f * (2 Nf)^alpha with ln eps_f and alpha linear in lambda_p = (L/D) * sqrt(fy/100)
and, in one variant, in D. Fitted as Kashani fits it (ln strain on ln life, linear in the coefficients, interaction terms), then inverted to
predict the life of a held-out condition. Inside every fold the choice between the two variants uses leave-one-out error on the training
rows only. Baselines for life: B2 = lambda_p-only variant (the Kashani form), BG = one global Coffin-Manson line.

    python v7_life.py
"""

from __future__ import annotations

import json
import math
import os

import numpy as np

import v7_data as D

LAMBDA = math.sqrt(D.FY_NOMINAL / 100.0)
ALPHA_MAX = -0.05                 # numerical guard: the fatigue exponent must stay clearly negative
FACTOR2 = math.log10(2.0)         # |log10 ratio| within this = within a factor of 2
VARIANTS = ("lam_d", "lam")       # selectable variants; "global" is a baseline only


def rows() -> np.ndarray:
    """Columns: D, L/D, amplitude (%), cycles to failure."""
    return np.array([(f["d"], f["ld"], f["n"], f["nf"]) for f in D.build_dev_features()], float)


def _x(d, ld, variant):
    lam = np.asarray(ld, float) * LAMBDA
    cols = [np.ones_like(lam)]
    if variant != "global":
        cols.append(lam)
    if variant == "lam_d":
        cols.append(np.asarray(d, float))
    return np.column_stack(cols)


def fit(r: np.ndarray, variant: str) -> np.ndarray:
    x = _x(r[:, 0], r[:, 1], variant)
    ln2n = np.log(2.0 * r[:, 3])[:, None]
    design = np.hstack([x, x * ln2n])
    coef, *_ = np.linalg.lstsq(design, np.log(r[:, 2] / 100.0), rcond=None)
    return coef


def predict(coef: np.ndarray, r: np.ndarray, variant: str) -> np.ndarray:
    x = _x(r[:, 0], r[:, 1], variant)
    k = x.shape[1]
    c, a = x @ coef[:k], np.minimum(x @ coef[k:], ALPHA_MAX)
    return np.exp((np.log(r[:, 2] / 100.0) - c) / a) / 2.0


def log_ratio(r: np.ndarray, pred: np.ndarray) -> np.ndarray:
    return np.log10(pred / r[:, 3])


def select_variant(train: np.ndarray) -> str:
    """Pick the variant with the smaller leave-one-out RMS log-life error on the training rows."""
    best, best_rms = None, np.inf
    for v in VARIANTS:
        err = []
        for i in range(len(train)):
            rest = np.delete(train, i, axis=0)
            err.append(log_ratio(train[i:i + 1], predict(fit(rest, v), train[i:i + 1], v))[0])
        rms = float(np.sqrt(np.mean(np.square(err))))
        if rms < best_rms:
            best, best_rms = v, rms
    return best


def folds(r: np.ndarray) -> dict:
    d, ld, n = r[:, 0], r[:, 1], r[:, 2]
    idx = np.arange(len(r))
    g = {"leave one condition out": [(f"{int(d[i])}/{int(ld[i])}/{int(n[i])}", [i]) for i in idx]}
    g["leave one diameter out"] = [(f"D{int(v)}", idx[d == v]) for v in sorted(set(d))]
    g["leave one L/D out"] = [(f"LD{int(v)}", idx[ld == v]) for v in sorted(set(ld))]
    g["leave one amplitude out"] = [(f"{int(v)}%", idx[n == v]) for v in sorted(set(n))]
    g["leave out D x L/D cell"] = [(f"D{int(a)} LD{int(b)}", idx[(d == a) & (ld == b)]) for a in sorted(set(d)) for b in sorted(set(ld))]
    g["leave out D x amplitude cell"] = [(f"D{int(a)} {int(b)}%", idx[(d == a) & (n == b)]) for a in sorted(set(d)) for b in sorted(set(n))]
    g["leave out L/D x amplitude cell"] = [(f"LD{int(a)} {int(b)}%", idx[(ld == a) & (n == b)]) for a in sorted(set(ld)) for b in sorted(set(n))]
    g["extrapolation: amplitude >= 5 % held out"] = [("n>=5", idx[n >= 5])]
    g["extrapolation: L/D >= 12 held out"] = [("ld>=12", idx[ld >= 12])]
    g["extrapolation: D = 20 mm held out"] = [("D20", idx[d == 20])]
    rng = np.random.default_rng(7)
    k = int(round(0.2 * len(r)))
    g["random 20 % splits (200)"] = [(f"rep{j}", rng.choice(idx, size=k, replace=False)) for j in range(200)]
    return g


def summarise(x: np.ndarray) -> dict:
    return {"n": int(x.size), "median_abs_log10": float(np.median(np.abs(x))), "within_factor2": float(np.mean(np.abs(x) <= FACTOR2)),
            "mean_signed_log10": float(np.mean(x))}


def main() -> int:
    r = rows()
    print(f"development rows: {len(r)}")
    results, chosen = {}, {}
    pooled = {"framework": [], "B2": [], "BG": []}
    for group, fl in folds(r).items():
        acc = {"framework": [], "B2": [], "BG": []}
        for name, test in fl:
            test = np.asarray(test)
            train = np.delete(r, test, axis=0)
            v = select_variant(train)
            chosen.setdefault(group, {}).setdefault(v, 0)
            chosen[group][v] += 1
            rt = r[test]
            acc["framework"].append(log_ratio(rt, predict(fit(train, v), rt, v)))
            acc["B2"].append(log_ratio(rt, predict(fit(train, "lam"), rt, "lam")))
            acc["BG"].append(log_ratio(rt, predict(fit(train, "global"), rt, "global")))
        results[group] = {k: summarise(np.concatenate(v)) for k, v in acc.items()}
        if group in ("leave one condition out", "leave one amplitude out"):
            for k in pooled:
                pooled[k].append(np.concatenate(acc[k]))
    print(f"\n{'group':44s} {'model':10s} {'median |log10|':>14s} {'within x2':>10s} {'bias':>8s}  n")
    for group, res in results.items():
        for k, s in res.items():
            print(f"{group:44s} {k:10s} {s['median_abs_log10']:14.3f} {s['within_factor2']:10.2f} {s['mean_signed_log10']:8.3f}  {s['n']}")
    print("\nvariant chosen inside the folds:", {g: c for g, c in chosen.items()})
    acc_rule = summarise(np.concatenate(pooled["framework"]))
    print(f"\nACCEPTANCE (life part), pooled over leave-one-condition-out and leave-one-amplitude-out: median |log10| = "
          f"{acc_rule['median_abs_log10']:.3f} (limit {FACTOR2:.3f}) -> {'PASS' if acc_rule['median_abs_log10'] <= FACTOR2 else 'FAIL'}")
    for v in VARIANTS + ("global",):
        print(f"fit on all development rows ({v}):", np.round(fit(r, v), 4).tolist())
    os.makedirs(D.OUT, exist_ok=True)
    with open(os.path.join(D.OUT, "life_cv.json"), "w") as f:
        json.dump({"results": results, "chosen": chosen, "acceptance_life": acc_rule}, f, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
