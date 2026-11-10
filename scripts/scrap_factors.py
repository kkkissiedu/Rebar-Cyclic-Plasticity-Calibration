"""
scrap_factors.py
================

Candidate stress-level factors f_s for the scrap-metal bars, from their monotonic tension tests (uts_parse.py output), against the measured B500C tension
values of Moodley, Afshan and De Risi (2024, Structures, Table 2: four tests on 12 and 16 mm B500C, 50 mm clip-on extensometer; mean fy 577, fu 669 MPa).
These are used in place of Kashani's own measured values, which are not available; the Moodley bars are from the same research line but a later study, and
that they match the Kashani 2019 bars is NOT confirmed. Two rules, both HYPOTHESES to be tested against cyclic tests of the same bars:

    f_s(Ry) = Ry_bar / fy_ref          f_s(Rm) = Rm_bar / fu_ref

No factor for life (f_e) is computed: the scrap files give no usable strain at peak force or fracture elongation (Moodley reports 8.5 to 12.6 % at fu).

    python uts_parse.py && python scrap_factors.py   ->  calibration_sessions/uniaxial/scrap_factors.csv
"""

from __future__ import annotations

import csv
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # repository root
SUMMARY = os.path.join(HERE, "calibration_sessions", "uniaxial", "uniaxial_summary.csv")
OUT = os.path.join(HERE, "calibration_sessions", "uniaxial", "scrap_factors.csv")
FY_REF = float(np.mean([567, 553, 591, 597]))    # Moodley et al. (2024) Table 2, B500C 12 and 16 mm, four tests, 50 mm clip-on extensometer
FU_REF = float(np.mean([663, 651, 679, 684]))
sys.path.insert(0, HERE)


def main() -> int:
    rm, ry = defaultdict(list), defaultdict(list)
    for r in csv.DictReader(open(SUMMARY, encoding="utf-8")):
        if "aborted" in r["flags"]:
            continue
        rm[r["source"]].append(float(r["rm_mpa"]))
        if r["ext_ok"] == "True" and r["ry_0.2_mpa"] not in ("nan", ""):
            ry[r["source"]].append(float(r["ry_0.2_mpa"]))
    rows = []
    for s in sorted(rm):
        a = float(np.mean(rm[s]))
        b = float(np.mean(ry[s])) if ry[s] else float("nan")
        rows.append([s, len(rm[s]), round(a), len(ry[s]), round(b) if b == b else "", round(b / FY_REF, 3) if b == b else "", round(a / FU_REF, 3)])
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["source", "n_bars_rm", "rm_mean_mpa", "n_bars_ry", "ry_mean_mpa", "f_s_from_ry", "f_s_from_rm"])
        w.writerows(rows)
    print(f"anchor: fy {FY_REF:.0f}, fu {FU_REF:.0f} MPa (Moodley et al. 2024, B500C 12 and 16 mm; match to the Kashani bars not confirmed)")
    for r in rows:
        print(r)
    try:
        from calibration_app.predictor import predict_bar
        base = predict_bar(12, 5, 3)["cycle2_peaks"]["tension_mpa"]
        print(f"\nExample, 12 mm L/D 5 3 %: fitted-steel cycle-2 tension peak {base:.0f} MPa; with each source's f_s(Rm):")
        for r in rows:
            print(f"  {r[0]:18s} f_s {r[6]:.2f} -> {predict_bar(12, 5, 3, factors=(r[6], 1.0))['cycle2_peaks']['tension_mpa']:.0f} MPa")
    except Exception as exc:                                  # predictor optional here
        print("predictor example skipped:", exc)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
