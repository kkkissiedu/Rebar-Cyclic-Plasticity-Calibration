"""
v7_loops_cv.py
==============

Full-chain cross-validation of the V7 framework on the development set (V7_PROTOCOL.md sections 4 and 5, Amendment 1).

Inside every fold, on the training files only: life layer (v7_life), degradation layer (v7_decay), loop core (v7_loop, fitted on the shortest
available L/D) and buckling distortion (v7_buckle). Held-out loops (cycle 2 and the cycles nearest to life fractions 0.1 .. 0.9) are predicted
from PREDICTED life and peaks and compared with the observed loops. Baselines: B1 (cycle-2 loop repeated), B3 (elastic-perfectly-plastic),
core alone scaled to the same predicted peaks. A shape-only variant (observed peaks and observed Nf given) is also reported.

    python v7_loops_cv.py --quick        # one fold, with timing
    python v7_loops_cv.py                # leave-one-condition/diameter/L/D/amplitude out + extrapolation blocks, in parallel
    python v7_loops_cv.py --all          # adds the cell groups and the 200 random splits
"""

from __future__ import annotations

import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import v7_buckle as B
import v7_data as D
import v7_decay as DC
import v7_life as L1
import v7_loop as V

DEFAULT_GROUPS = ("leave one condition out", "leave one diameter out", "leave one L/D out", "leave one amplitude out",
                  "extrapolation: amplitude >= 5 % held out", "extrapolation: L/D >= 12 held out", "extrapolation: D = 20 mm held out")
_FILES = None


def load_all() -> list:
    loops = {(f["d"], f["ld"], f["n"]): f for f in V.build_dev_loops()}
    out = []
    for f in DC.make_files():
        lp = loops.get((f["d"], f["ld"], f["n"]))
        if lp is not None:
            g = dict(f)
            g["loops"] = lp["loops"]
            out.append(g)
    return out


def files() -> list:
    global _FILES
    if _FILES is None:
        _FILES = load_all()
    return _FILES


def epp_unit(eps: float, e_app: float, sy: float) -> np.ndarray:
    path = V.tri_path(eps)
    out, s, prev = np.empty_like(path), 0.0, 0.0
    for i, e in enumerate(path):
        s = max(-sy, min(sy, s + e_app * (e - prev)))
        prev = e
        out[i] = s
    return out[V.PTS:2 * V.PTS + 1]


def fit_epp(train: list, km: float) -> float:
    best, best_err = 300.0, np.inf
    for sy in np.linspace(100.0, 800.0, 36):
        err = 0.0
        for f in train:
            u = epp_unit(f["n"] / 100.0, V.e_app(f["d"], f["ld"], km), sy)
            err += float(np.mean((np.interp(V.S_GRID, V.S_IDEAL, u) - f["loops"][2]["y"]) ** 2))
        if err < best_err:
            best, best_err = sy, err
    return best


def loo_log_rms(rtr: np.ndarray, variant: str) -> float:
    """RMS of the leave-one-out log10 life residuals of the L1 model on the training rows (the scatter used for the life distribution)."""
    err = []
    for i in range(len(rtr)):
        rest = np.delete(rtr, i, axis=0)
        err.append(L1.log_ratio(rtr[i:i + 1], L1.predict(L1.fit(rest, variant), rtr[i:i + 1], variant))[0])
    return float(np.sqrt(np.mean(np.square(err))))


def energy(sig: np.ndarray, e: np.ndarray) -> float:
    return float(abs(np.trapz(sig, e)))


def run_fold(task: tuple) -> list:
    group, name, test_idx = task
    allf = files()
    test_idx = set(int(i) for i in test_idx)
    train = [f for i, f in enumerate(allf) if i not in test_idx]
    test = [allf[i] for i in sorted(test_idx)]
    # L1 life, L2 decay and cycle-2 peaks
    rtr = np.array([(f["d"], f["ld"], f["n"], f["nf"]) for f in train], float)
    lv = L1.select_variant(rtr)
    lc = L1.fit(rtr, lv)
    dv = DC.select_decay(train)
    th_t, th_k = DC.fit_decay(train, "yt", dv), DC.fit_decay(train, "yk", dv)
    pv = DC.select_p2(train)
    p2c = DC.fit_p2(train, pv)
    # L3a core on the shortest available L/D, L3b distortion on all training loops
    short_ld = min(f["ld"] for f in train)
    core = V.fit_joint([f for f in train if f["ld"] == short_ld], max_nfev=150)
    km = math.exp(core[0])
    params = B.fit(train, core, mode="joint")
    sy_epp = fit_epp([f for f in train if f["ld"] == short_ld], km)
    s_log = loo_log_rms(rtr, lv)
    zq, wq = np.polynomial.hermite.hermgauss(9)
    rows = []
    for f in test:
        rt = np.array([[f["d"], f["ld"], f["n"], f["nf"]]], float)
        nfp = float(L1.predict(lc, rt, lv)[0])
        t2p, k2p = (a[0] for a in DC.pred_p2(p2c, [f], pv))
        for c, lp in f["loops"].items():
            one = dict(f, c=np.array([float(c)]))
            gt = float(DC.pred_g(one, th_t, dv, nfp)[0]); gk = float(DC.pred_g(one, th_k, dv, nfp)[0])
            tp, kp = (t2p, k2p) if c == 2 else (t2p * gt, k2p * gk)
            spec = {"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": nfp, "cycle": c, "loop": lp}
            spec2 = dict(spec, cycle=2)
            spec_true = dict(spec, nf=f["nf"])
            preds = {
                "framework": B.predict(spec, core, params, {"t": tp, "k": kp}),
                "core_only": B.core_predict(spec, core, {"t": tp, "k": kp}),
                "B1": B.predict(spec2, core, params, {"t": t2p, "k": k2p}),
                "B3": np.interp(lp["arc"], V.S_IDEAL, epp_unit(f["n"] / 100.0, V.e_app(f["d"], f["ld"], km), sy_epp)),
                "shape_only": B.predict(spec_true, core, params, {"t": lp["t"], "k": lp["k"]}),
            }
            avg = np.zeros(lp["s"].shape, float)                     # expected loop over the log-normal life distribution (Amendment 2)
            for z, w in zip(zq, wq):
                nfi = nfp * 10.0 ** (s_log * math.sqrt(2.0) * z)
                gti = float(DC.pred_g(one, th_t, dv, nfi)[0]); gki = float(DC.pred_g(one, th_k, dv, nfi)[0])
                ti, ki = (t2p, k2p) if c == 2 else (t2p * gti, k2p * gki)
                avg += w / math.sqrt(math.pi) * B.predict(dict(spec, nf=nfi), core, params, {"t": ti, "k": ki})
            preds["framework_avg"] = avg
            eo = energy(lp["s"], lp["e"])
            row = {"group": group, "fold": name, "ld": f["ld"], "d": f["d"], "n": f["n"], "c": int(c), "x": c / f["nf"]}
            for k, p in preds.items():
                row["rmse_" + k] = float(np.sqrt(np.mean((p - lp["s"]) ** 2)))
                row["en_" + k] = float(abs(energy(p, lp["e"]) - eo) / eo * 100.0)
            rows.append(row)
    return rows


def main() -> int:
    quick, run_all = "--quick" in sys.argv, "--all" in sys.argv
    allf = files()
    r = np.array([(f["d"], f["ld"], f["n"], f["nf"]) for f in allf], float)
    groups = L1.folds(r)
    tasks = []
    for g, fl in groups.items():
        if not run_all and g not in DEFAULT_GROUPS:
            continue
        tasks += [(g, name, [int(i) for i in test]) for name, test in fl]
    if quick:
        tasks = tasks[:1]
    print(f"{len(allf)} files, {len(tasks)} folds", flush=True)
    t0 = time.time()
    rows = []
    if quick:
        rows = run_fold(tasks[0])
        print(f"one fold: {time.time() - t0:.0f} s, {len(rows)} loops")
        for k in ("framework", "framework_avg", "core_only", "B1", "B3", "shape_only"):
            print(f"  {k:10s} median RMSE {np.median([x['rmse_' + k] for x in rows]):.1f} MPa, energy error {np.median([x['en_' + k] for x in rows]):.1f} %")
        return 0
    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 4) - 3)) as ex:
        for i, res in enumerate(ex.map(run_fold, tasks, chunksize=1)):
            rows += res
            if (i + 1) % 10 == 0:
                print(f"[{time.time() - t0:6.0f}s] {i + 1}/{len(tasks)} folds", flush=True)
    os.makedirs(D.OUT, exist_ok=True)
    json.dump(rows, open(os.path.join(D.OUT, "loops_cv.json"), "w"))
    models = ("framework", "framework_avg", "core_only", "B1", "B3", "shape_only")
    print(f"\n{'group':44s} " + " ".join(f"{m:>11s}" for m in models) + "   (median loop RMSE, MPa)")
    for g in dict.fromkeys(x["group"] for x in rows):
        sel = [x for x in rows if x["group"] == g]
        print(f"{g:44s} " + " ".join(f"{np.median([x['rmse_' + m] for x in sel]):11.1f}" for m in models))
    print("\nby L/D, all groups pooled (median RMSE):")
    for ld in (5, 8, 10, 12, 15):
        sel = [x for x in rows if x["ld"] == ld]
        print(f"  L/D {ld:2d}: " + " ".join(f"{m} {np.median([x['rmse_' + m] for x in sel]):6.1f}" for m in models))
    print("\nmedian energy error per cycle (%), pooled: " + ", ".join(f"{m} {np.median([x['en_' + m] for x in rows]):.1f}" for m in models))
    sel = [x for x in rows if x["group"] in ("leave one condition out", "leave one amplitude out")]
    folds_ = {}
    for x in sel:
        folds_.setdefault((x["group"], x["fold"]), []).append(x)
    wins_b1 = np.mean([np.median([x["rmse_framework"] for x in v]) < np.median([x["rmse_B1"] for x in v]) for v in folds_.values()])
    wins_core = np.mean([np.median([x["rmse_framework"] for x in v]) < np.median([x["rmse_core_only"] for x in v]) for v in folds_.values()])
    en_ok = np.median([x["en_framework"] for x in sel]) < np.median([x["en_B1"] for x in sel])
    med = lambda rs, k: float(np.median([x[k] for x in rs]))
    c1 = med(sel, "rmse_framework_avg") <= med(sel, "rmse_framework")
    sub = [x for x in rows if x["ld"] == 5 and x["n"] >= 4 and x["c"] > 2]
    c2 = med(sub, "rmse_framework_avg") < med(sub, "rmse_framework")
    print(f"AMENDMENT 2 RULE: (i) leave-one-condition/amplitude-out median {med(sel, 'rmse_framework_avg'):.1f} (averaged) vs {med(sel, 'rmse_framework'):.1f} (point) -> {c1}; "
          f"(ii) L/D 5, 4-6 %, later cycles {med(sub, 'rmse_framework_avg'):.1f} vs {med(sub, 'rmse_framework'):.1f} -> {c2}; averaged version adopted: {c1 and c2}")
    print(f"\nPASS RULE (full chain): beats B1 on {wins_b1:.0%} of folds, beats core-only on {wins_core:.0%} (each needs >= 80 %); "
          f"pooled median energy error below B1: {en_ok} -> {'PASS' if wins_b1 >= 0.8 and wins_core >= 0.8 and en_ok else 'FAIL'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
