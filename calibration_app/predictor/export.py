"""
export.py
=========

Write a predict_bar() result to disk in a form like the Kashani et al. (2019) test files.

CSV (one row per point, all requested cycles): ``Time``, ``Cycle``, ``Position mm``, ``Load kN`` (the same actuator channels the test files carry, with
position = strain x gauge length and load = stress x area), then ``Strain %`` and ``Stress MPa``. ``Time`` runs 0 to 20 s inside each cycle (the 20 s
period of the tests) and continues from the cycle number, so cycle c starts at (c - 1) x 20 s. A few comment lines at the top record the inputs and
the fitted-range warnings. PNG: loops of the requested cycles and the peak-stress history.
"""

from __future__ import annotations

import csv
import math
import os
from typing import Tuple

import numpy as np

PERIOD_S = 20.0


def _geometry(result: dict) -> Tuple[float, float]:
    i = result["inputs"]
    return i["diameter_mm"] * i["L_over_D"], math.pi * i["diameter_mm"] ** 2 / 4.0


def export_csv(result: dict, path: str) -> str:
    gauge, area = _geometry(result)
    i = result["inputs"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"# Predicted loops (V7 framework): D {i['diameter_mm']:g} mm, L/D {i['L_over_D']:g}, amplitude {i['amplitude_pct']:g} %, "
                    f"fy {i['fy_mpa']:g} MPa, strain basis {i['strain_basis']}, f_s {i['f_s']:g}, f_e {i['f_e']:g}"])
        w.writerow([f"# Predicted cycles to failure {result['life']['cycles_to_failure']:.1f} (factor-2 band "
                    f"{result['life']['factor2_band'][0]:.1f} to {result['life']['factor2_band'][1]:.1f})"])
        for msg in result["warnings"]:
            w.writerow([f"# WARNING: {msg}"])
        w.writerow(["Time", "Cycle", "Position mm", "Load kN", "Strain %", "Stress MPa"])
        for c in sorted(result["loops"]):
            lp = result["loops"][c]
            t = np.linspace(0.0, PERIOD_S, len(lp["strain"])) + (c - 1) * PERIOD_S
            for tt, e, s in zip(t, lp["strain"], lp["stress"]):
                w.writerow([f"{tt:.2f}", c, f"{e * gauge:.5f}", f"{s * area / 1000.0:.5f}", f"{e * 100.0:.5f}", f"{s:.3f}"])
    return path


def export_png(result: dict, base_path: str) -> list:
    """Save <base>_loops.png and <base>_peaks.png; returns the two paths."""
    from matplotlib.figure import Figure
    i = result["inputs"]
    title = f"D {i['diameter_mm']:g} mm, L/D {i['L_over_D']:g}, {i['amplitude_pct']:g} % strain"
    base = os.path.splitext(base_path)[0]
    f1 = Figure(figsize=(7, 5.5), dpi=150)
    ax = f1.subplots()
    for c in sorted(result["loops"]):
        lp = result["loops"][c]
        ax.plot(lp["strain"] * 100, lp["stress"], lw=1.3, label=f"cycle {c}")
    ax.set(xlabel="engineering strain (%)", ylabel="stress (MPa)", title=f"Predicted loops, {title}")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    f1.tight_layout(); p1 = base + "_loops.png"; f1.savefig(p1)
    f2 = Figure(figsize=(7, 5.5), dpi=150)
    ax = f2.subplots()
    ax.plot(result["peaks"]["cycle"], result["peaks"]["tension"], label="tension")
    ax.plot(result["peaks"]["cycle"], [-x for x in result["peaks"]["compression"]], label="compression")
    ax.set(xlabel="cycle", ylabel="peak stress (MPa)", title=f"Predicted peak stress, {title}")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    f2.tight_layout(); p2 = base + "_peaks.png"; f2.savefig(p2)
    return [p1, p2]
