"""
v7_transfer.py
==============

New-bar transfer test on the DEVELOPMENT set (V7_PROTOCOL.md section 6). Each development diameter in turn is treated as a new bar: the framework is
fitted on the other two diameters, the two material factors (f_e on the fatigue ductility coefficient, f_s on the stress level) are calibrated only from
that diameter's L/D 5 tests at 2 % and 4 %, and its other 28 files are predicted with and without the factors.

Scale 1.0 is the real data (one steel, so the factors only absorb diameter-specific offsets). Scale 1.2 multiplies every stress of the new diameter by 1.2,
a synthetic 'stronger steel' with the same ductility: it checks that the stress factor recovers a KNOWN change. This tests the machinery of the factors; it
cannot show transfer to a different real steel, because the dataset contains one steel.

    python v7_transfer.py
"""

from __future__ import annotations

import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import v7_data as D
import v7_loops_cv as LC
import v7_predict as P

CALIB = ((5, 2), (5, 4))      # (L/D, amplitude %) of the tests used to calibrate the factors


def perturb(f: dict, s: float) -> dict:
    g = dict(f)
    for k in ("t2", "k2", "t", "k"):
        g[k] = f[k] * s
    g["loops"] = {c: dict(lp, s=lp["s"] * s, y=lp["y"] * s, t=lp["t"] * s, k=lp["k"] * s, w=lp["w"] * s) for c, lp in f["loops"].items()}
    return g


def summarise(outs: list) -> dict:
    loops = [r for o in outs for r in o["loops"]]
    life = np.array([abs(o["log_life_err"]) for o in outs])
    return {"loop_rmse_framework": float(np.median([r["rmse_fw"] for r in loops])), "loop_rmse_B1": float(np.median([r["rmse_b1"] for r in loops])),
            "cycle2_tension_err": float(np.median([o["t2_err"] for o in outs])), "cycle2_compression_err": float(np.median([o["k2_err"] for o in outs])),
            "life_median_abs_log10": float(np.median(life)), "life_within_factor2": float(np.mean(life <= np.log10(2.0)))}


def task(args: tuple) -> dict:
    d_new, scale = args
    allf = LC.files()
    train = [f for f in allf if f["d"] != d_new]
    new = [perturb(f, scale) if scale != 1.0 else f for f in allf if f["d"] == d_new]
    calib = [f for f in new if (f["ld"], f["n"]) in CALIB]
    test = [f for f in new if (f["ld"], f["n"]) not in CALIB]
    m = P.fit_models(train)
    f_e, f_s = P.calib_factors(m, calib)
    res = {"d_new": d_new, "scale": scale, "f_e": f_e, "f_s": f_s, "n_test": len(test)}
    res["no_factors"] = summarise([P.predict_file(m, f) for f in test])
    res["with_factors"] = summarise([P.predict_file(m, f, f_s=f_s, f_e=f_e) for f in test])
    return res


def main() -> int:
    tasks = [(d, s) for d in (10, 12, 20) for s in (1.0, 1.2)]
    with ProcessPoolExecutor(max_workers=len(tasks)) as ex:
        results = list(ex.map(task, tasks))
    os.makedirs(D.OUT, exist_ok=True)
    json.dump(results, open(os.path.join(D.OUT, "transfer_dev.json"), "w"), indent=1)
    print(f"{'new bar':8s} {'stress x':>8s} | {'f_s':>5s} {'f_e':>5s} | loop RMSE MPa: none -> factors (B1 none -> factors) | cycle-2 peak err T/K %: none -> factors | life |log10|: none -> factors")
    for r in results:
        a, b = r["no_factors"], r["with_factors"]
        print(f"D{r['d_new']:<7d} {r['scale']:8.1f} | {r['f_s']:5.2f} {r['f_e']:5.2f} | "
              f"{a['loop_rmse_framework']:5.1f} -> {b['loop_rmse_framework']:5.1f} ({a['loop_rmse_B1']:5.1f} -> {b['loop_rmse_B1']:5.1f}) | "
              f"{100 * a['cycle2_tension_err']:4.1f}/{100 * a['cycle2_compression_err']:4.1f} -> {100 * b['cycle2_tension_err']:4.1f}/{100 * b['cycle2_compression_err']:4.1f} | "
              f"{a['life_median_abs_log10']:.3f} -> {b['life_median_abs_log10']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
