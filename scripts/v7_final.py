"""
v7_final.py
===========

Single evaluation of the frozen V7 framework on the sealed set (all D = 16 mm files), as set out in V7_PROTOCOL.md (Amendment 3).

    python v7_final.py --dry-run     # proves the code path with the 20 mm development diameter standing in for the sealed set (no sealed data read)
    set V7_UNSEAL=1 & python v7_final.py     # the real, single evaluation; refuses to run twice (the freeze log records its use)

The framework is fitted on all 90 development files, then predicts every sealed file with NO use of its data. A second, separate result applies the
material factors calibrated from that diameter's L/D 5 tests at 2 % and 4 % to the other 28 files (the new-bar claim).
"""

from __future__ import annotations

import datetime
import json
import math
import os
import sys

import numpy as np

import v7_data as D
import v7_decay as DC
import v7_loop as V
import v7_loops_cv as LC
import v7_predict as P

FREEZE = os.path.join(D.PROJ, "docs", "protocols", "V7_FREEZE.txt")
CALIB = ((5, 2), (5, 4))
LIFE_LIMIT = math.log10(2.0)


def build_files(conds: list) -> list:
    feats = [D.extract(*c) for c in conds]            # raises SealedError for 16 mm unless V7_UNSEAL=1
    loops = {(f["d"], f["ld"], f["n"]): V.read_loops(f["d"], f["ld"], f["n"], f["nf"]) for f in feats}
    out = []
    for g in DC.make_files(feats):
        lp = loops.get((g["d"], g["ld"], g["n"]))
        if not lp:
            continue
        ld = {}
        for c, (e, s) in lp.items():
            sa = V.arc_fraction(e)
            ld[int(c)] = {"e": e, "s": s, "arc": sa, "y": np.interp(V.S_GRID, sa, s), "t": float(s.max()), "k": float(-s.min()),
                          "w": float(abs(np.trapz(s, e)))}
        out.append(dict(g, loops=ld))
    return out


def summarise(outs: list) -> dict:
    loops = [r for o in outs for r in o["loops"]]
    fw = [float(np.median([r["rmse_fw"] for r in o["loops"]])) for o in outs]
    b1 = [float(np.median([r["rmse_b1"] for r in o["loops"]])) for o in outs]
    life = np.array([abs(o["log_life_err"]) for o in outs])
    return {"files": len(outs), "loop_rmse_framework": float(np.median([r["rmse_fw"] for r in loops])), "loop_rmse_B1": float(np.median([r["rmse_b1"] for r in loops])),
            "share_files_framework_better_than_B1": float(np.mean(np.array(fw) < np.array(b1))),
            "energy_err_framework": float(np.median([r["en_fw"] for r in loops])), "energy_err_B1": float(np.median([r["en_b1"] for r in loops])),
            "life_median_abs_log10": float(np.median(life)), "life_within_factor2": float(np.mean(life <= LIFE_LIMIT)),
            "cycle2_tension_err": float(np.median([o["t2_err"] for o in outs])), "cycle2_compression_err": float(np.median([o["k2_err"] for o in outs]))}


def breakdown(test: list, outs: list, key: str) -> dict:
    res = {}
    for v in sorted({f[key] for f in test}):
        sel = [o for f, o in zip(test, outs) if f[key] == v]
        loops = [r for o in sel for r in o["loops"]]
        res[int(v)] = {"loop_rmse_framework": float(np.median([r["rmse_fw"] for r in loops])), "loop_rmse_B1": float(np.median([r["rmse_b1"] for r in loops])),
                       "life_median_abs_log10": float(np.median([abs(o["log_life_err"]) for o in sel]))}
    return res


def verdict(s: dict) -> dict:
    return {"life": s["life_median_abs_log10"] <= LIFE_LIMIT, "loops_vs_B1_80pct_of_files": s["share_files_framework_better_than_B1"] >= 0.8,
            "energy_below_B1": s["energy_err_framework"] < s["energy_err_B1"]}


def main() -> int:
    dry, again = "--dry-run" in sys.argv, "--again" in sys.argv
    log = open(FREEZE, encoding="utf-8").read()
    if not dry and "SEALED SET USED" in log and not again:
        print("The sealed set has already been used (see V7_FREEZE.txt). Refusing to run again without --again; any rerun must be reported as post-final.")
        return 2
    dev = LC.files()
    if dry:
        train, test = [f for f in dev if f["d"] != 20], [f for f in dev if f["d"] == 20]
        name = "final_dry_run.json"
    else:
        train, test = dev, build_files(D.conditions(sealed=True))
        name = "final_sealed.json"
    print(f"training files {len(train)}, evaluation files {len(test)}" + (" (DRY RUN: 20 mm development bars stand in)" if dry else " (SEALED 16 mm)"), flush=True)
    m = P.fit_models(train)
    outs = [P.predict_file(m, f) for f in test]
    head = summarise(outs)
    calib = [f for f in test if (f["ld"], f["n"]) in CALIB]
    rest = [f for f in test if (f["ld"], f["n"]) not in CALIB]
    f_e, f_s = P.calib_factors(m, calib)
    transfer = {"f_e": f_e, "f_s": f_s, "no_factors": summarise([P.predict_file(m, f) for f in rest]),
                "with_factors": summarise([P.predict_file(m, f, f_s=f_s, f_e=f_e) for f in rest])}
    result = {"mode": "dry_run" if dry else "sealed", "headline": head, "verdict": verdict(head), "by_L_over_D": breakdown(test, outs, "ld"),
              "by_amplitude": breakdown(test, outs, "n"), "transfer": transfer}
    os.makedirs(D.OUT, exist_ok=True)
    json.dump(result, open(os.path.join(D.OUT, name), "w"), indent=1)
    print("HEADLINE (no factors):", json.dumps(head, indent=1))
    print("VERDICT against the pre-set rules:", result["verdict"])
    print("by L/D:", {k: {a: round(b, 3) for a, b in v.items()} for k, v in result["by_L_over_D"].items()})
    print("by amplitude:", {k: {a: round(b, 3) for a, b in v.items()} for k, v in result["by_amplitude"].items()})
    print(f"TRANSFER: f_e {f_e:.3f}, f_s {f_s:.3f}")
    for k in ("no_factors", "with_factors"):
        t = transfer[k]
        print(f"  {k:13s} loop RMSE {t['loop_rmse_framework']:.1f} MPa (B1 {t['loop_rmse_B1']:.1f}), life |log10| {t['life_median_abs_log10']:.3f}, "
              f"cycle-2 peak error T/K {100 * t['cycle2_tension_err']:.1f}/{100 * t['cycle2_compression_err']:.1f} %")
    if not dry:
        with open(FREEZE, "a", encoding="utf-8", newline="\n") as f:
            f.write(f"SEALED SET USED {datetime.datetime.now().isoformat(timespec='seconds')}: single evaluation of the frozen framework; results in "
                    f"calibration_sessions/v7/final_sealed.json\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
