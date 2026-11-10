"""
v7_export.py
============

Refits the frozen V7 framework on ALL 120 files (development set plus the 16 mm bars, after the single sealed evaluation is done) and exports the fitted
parameters for the shipped predictor (RebarCalibrator/calibration_app/predictor/predictor_params.json). The validation numbers in V7_PROTOCOL.md refer to the
90-file fit; the exported parameters are the same procedure refitted on all data, as is usual after cross-validation.

    set V7_UNSEAL=1 & python v7_export.py
"""

from __future__ import annotations

import datetime
import hashlib
import json
import math
import os

import numpy as np

import v7_buckle as B
import v7_data as D
import v7_final as F
import v7_life as L1
import v7_loops_cv as LC
import v7_predict as P

PKG = os.path.join(D.PROJ, "calibration_app", "predictor")


def main() -> int:
    files = LC.files() + F.build_files(D.conditions(sealed=True))
    m = P.fit_models(files)
    final = json.load(open(os.path.join(D.OUT, "final_sealed.json")))
    out = {
        "version": "V7, refit on all 120 files",
        "trained_on": {"files": len(files), "diameters_mm": [10, 12, 16, 20], "L_over_D": [5, 8, 10, 12, 15], "amplitude_pct": [1, 2, 3, 4, 5, 6], "steel": "B500C (Kashani et al. 2019)"},
        "constants": {"fy_nominal": D.FY_NOMINAL, "lambda_ref": B.LAMBDA_REF, "alpha_max": L1.ALPHA_MAX, "delta": B.DELTA, "e_bar": 200000.0},
        "life": {"variant": m["lv"], "coef": np.asarray(m["lc"]).tolist(), "s_log10": m["s_log"]},
        "peak2": {"variant": m["pv"], "coef_t": np.asarray(m["p2c"][0]).tolist(), "coef_k": np.asarray(m["p2c"][1]).tolist()},
        "decay": {"variant": m["dv"], "theta_t": np.asarray(m["th_t"]).tolist(), "theta_k": np.asarray(m["th_k"]).tolist()},
        "core": {"theta": np.asarray(m["core"]).tolist(), "km_N_per_mm": m["km"]},
        "buckling": {"params": np.asarray(m["params"]).tolist()},
        "validation_reference": {"sealed_16mm_headline": final["headline"], "note": "Sealed result refers to the fit on the 90 development files."},
    }
    os.makedirs(PKG, exist_ok=True)
    path = os.path.join(PKG, "predictor_params.json")
    json.dump(out, open(path, "w"), indent=1)
    h = hashlib.sha256(open(path, "rb").read()).hexdigest()
    with open(os.path.join(D.PROJ, "docs", "protocols", "V7_FREEZE.txt"), "a", encoding="utf-8", newline="\n") as f:
        f.write(f"EXPORT {datetime.datetime.now().isoformat(timespec='seconds')}: shipped parameters refit on all 120 files; predictor_params.json sha256 {h}\n")
    print("exported", path, "| life variant", m["lv"], "| decay variant", m["dv"], "| k_m", round(m["km"] / 1000), "kN/mm | life scatter s =", round(m["s_log"], 3))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
