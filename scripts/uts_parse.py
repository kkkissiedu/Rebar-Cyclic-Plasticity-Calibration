"""
uts_parse.py
============

Reads the universal-testing-machine tension tests of the scrap-metal bars (UNIAXIAL TEST/12mm and /16mm; the Day1 and DAY 2 folders are copies) and writes a
per-bar summary with quality flags. Nothing is excluded silently: every bar is listed with its flags.

Columns of the .GRF text table: F(N) R(MPa) L(mm) A(%) t(s) Lb F/t R/t Lc/t En Lc(mm) Le(mm)[extensometer] Le/t. R is F / A0 with A0 from the .dat header
(float32 at byte 528, from the diameter entered at byte 504). The extensometer gauge length is not recorded in these files, so the elastic slope of
stress against extensometer extension is converted to an IMPLIED gauge length assuming E = 200 GPa (bar steel); a faulty extensometer shows up as an implied
gauge that departs from the consistent value of the other bars, or as a non-linear / noisy elastic fit.

    python uts_parse.py            ->  calibration_sessions/uniaxial/uniaxial_summary.csv
"""

from __future__ import annotations

import csv
import glob
import os
import struct

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# the tension-test files are not part of this repository: point UNIAXIAL_DIR at the folder holding the 12mm and 16mm sub-folders
ROOT = os.environ.get("UNIAXIAL_DIR", os.path.join(os.path.dirname(REPO), "UNIAXIAL TEST"))
OUT = os.path.join(REPO, "calibration_sessions", "uniaxial")
E_STEEL = 200000.0
GAUGE_MM = 50.0                    # assumed extensometer gauge (standard); the clean fits imply 48 to 57 mm, median ~50, but the files do not record it
FRACTURE_FRACTION = 0.5           # a test that never exceeds this share of the group median peak force is flagged aborted


def read_dat(path: str):
    if not os.path.exists(path):
        return float("nan"), float("nan")
    b = open(path, "rb").read()
    return struct.unpack_from("<f", b, 504)[0], struct.unpack_from("<f", b, 528)[0]


def read_grf(path: str) -> dict:
    rows = [l.rstrip("\r\n").split("\t") for l in open(path, encoding="latin-1").read().splitlines()[1:]]

    def col(i):
        return np.array([float(r[i]) if i < len(r) and r[i].strip() else np.nan for r in rows])
    return {"F": col(0), "R": col(1), "t": col(4), "Lc": col(10), "Le": col(11)}


def elastic_fit(R: np.ndarray, Le: np.ndarray, rm: float):
    """Slope of R against extensometer extension on the loading branch, 15 to 45 % of the peak stress."""
    ipk = int(np.nanargmax(R))
    m = np.zeros(R.size, bool)
    m[:ipk] = True
    m &= np.isfinite(Le) & (R > 0.15 * rm) & (R < 0.45 * rm)
    if m.sum() < 8 or np.ptp(Le[m]) < 1e-3:
        return float("nan"), float("nan")
    p = np.polyfit(Le[m], R[m], 1)
    res = R[m] - np.polyval(p, Le[m])
    r2 = 1.0 - res.var() / R[m].var()
    return float(p[0]), float(r2)


def proof_stress(R: np.ndarray, Le: np.ndarray, slope: float, gauge: float, offset: float = 0.002):
    """0.2 % offset proof stress: first crossing of the stress-extension curve with the elastic line shifted by offset x gauge."""
    ipk = int(np.nanargmax(R))
    ok = np.isfinite(Le[:ipk]) & (R[:ipk] > 0)
    x, y = Le[:ipk][ok], R[:ipk][ok]
    mm = np.polyfit(x[(y > 0.15 * y.max()) & (y < 0.45 * R[ipk])], y[(y > 0.15 * y.max()) & (y < 0.45 * R[ipk])], 1)
    line = slope * (x - offset * gauge) + (mm[1])
    d = y - line
    i = np.nonzero((d < 0) & (y > 0.5 * R[ipk] * 0.5))[0]
    if i.size == 0 or i[0] == 0 or x.max() < offset * gauge * 1.5:
        return float("nan")                              # extensometer did not record far enough past yield
    j = i[0]
    return float(y[j - 1] + (y[j] - y[j - 1]) * d[j - 1] / (d[j - 1] - d[j]))


def main() -> int:
    recs = []
    for size in ("12mm", "16mm"):
        for p in sorted(glob.glob(os.path.join(ROOT, size, "*", "*.GRF"))):
            name = os.path.basename(p)[:-4]
            src = os.path.basename(os.path.dirname(p))
            d, a0 = read_dat(p[:-4] + ".dat")
            c = read_grf(p)
            ipk = int(np.nanargmax(c["F"]))
            rm = float(c["R"][ipk])
            sl, r2 = elastic_fit(c["R"], c["Le"], rm)
            recs.append({"source": src, "size": size, "bar": name, "d_mm": d, "area_mm2": a0, "rm_mpa": rm, "fmax_kN": float(c["F"][ipk]) / 1000.0,
                         "slope": sl, "r2": r2, "implied_gauge_mm": E_STEEL / sl if sl == sl and sl > 0 else float("nan"),
                         "le_max_mm": float(np.nanmax(c["Le"])), "t_end_s": float(np.nanmax(c["t"])), "flags": [], "R": c["R"], "Le": c["Le"]})
    for r in recs:
        peers = [x["fmax_kN"] for x in recs if x["size"] == r["size"] and x["source"] == r["source"]]
        if r["fmax_kN"] < FRACTURE_FRACTION * np.median(peers):
            r["flags"].append("aborted_or_partial")
        if r["area_mm2"] != r["area_mm2"]:
            r["flags"].append("no_dat_header_area_unknown")
        if r["slope"] != r["slope"]:
            r["flags"].append("no_usable_extensometer")
        elif r["r2"] < 0.98:
            r["flags"].append("extensometer_elastic_fit_poor")
    good = np.array([r["implied_gauge_mm"] for r in recs if r["slope"] == r["slope"] and r["r2"] >= 0.98 and "aborted_or_partial" not in r["flags"]])
    ref = float(np.median(good)) if good.size else float("nan")
    for r in recs:
        g = r["implied_gauge_mm"]
        if g == g and ref == ref and "extensometer_elastic_fit_poor" not in r["flags"] and abs(g / ref - 1.0) > 0.25:
            r["flags"].append(f"implied_gauge_off_by_{100 * (g / ref - 1):+.0f}pct")
    for r in recs:                                     # E and proof stress only from bars whose extensometer passed every check
        r["ok"] = not any(k in " ".join(r["flags"]) for k in ("extensometer", "implied", "no_usable", "aborted"))
        r["e_gpa"] = r["ry"] = float("nan")
        if r["ok"]:
            r["e_gpa"] = r["slope"] * GAUGE_MM / 1000.0
            r["ry"] = proof_stress(r["R"], r["Le"], r["slope"], GAUGE_MM)
        r["ry_ratio"] = r["rm_mpa"] / r["ry"] if r["ry"] == r["ry"] else float("nan")
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "uniaxial_summary.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["source", "size", "bar", "d_mm", "area_mm2", "rm_mpa", "fmax_kN", "slope_mpa_per_mm", "r2", "implied_gauge_mm", "le_max_mm", "ext_ok",
                    "e_gpa", "ry_0.2_mpa", "rm_over_ry", "flags"])
        for r in recs:
            w.writerow([r["source"], r["size"], r["bar"], f"{r['d_mm']:.2f}", f"{r['area_mm2']:.1f}", f"{r['rm_mpa']:.0f}", f"{r['fmax_kN']:.1f}",
                        f"{r['slope']:.0f}", f"{r['r2']:.3f}", f"{r['implied_gauge_mm']:.1f}", f"{r['le_max_mm']:.2f}", r["ok"],
                        f"{r['e_gpa']:.0f}", f"{r['ry']:.0f}", f"{r['ry_ratio']:.3f}", ";".join(r["flags"])])
    print(f"median implied gauge of the clean fits: {ref:.1f} mm over {good.size} bars; wrote {path}")
    for r in recs:
        print(f"{r['bar']:14s} Rm {r['rm_mpa']:5.0f}  slope {r['slope']:8.0f}  R2 {r['r2']:.3f}  gauge {r['implied_gauge_mm']:7.1f}  E {r['e_gpa']:4.0f}  Ry {r['ry']:4.0f}  Rm/Ry {r['ry_ratio']:.2f}  {','.join(r['flags'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
