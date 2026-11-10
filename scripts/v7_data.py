"""
v7_data.py
==========

Data access for the V7 study (see V7_PROTOCOL.md, frozen 2026-11-10) with the sealed set enforced in code.

Development set: D = 10, 12, 20 mm. Sealed final test: D = 16 mm. Reading any 16 mm file raises SealedError unless the
environment variable V7_UNSEAL=1 is set, which must happen only after the framework and its parameters are frozen.

Per file the features are the cycles to failure (largest cycle number), and for every cycle the peak tension and compression
engineering stress (load / nominal area, initial load subtracted) and the loop energy (area of the stress-strain loop, with strain
= (position - first position) / (L/D * D)).
"""

from __future__ import annotations

import csv
import json
import math
import os

import numpy as np

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # repository root (this file is in scripts/)
import sys as _sys
_sys.path.insert(0, PROJ)                                              # so that calibration_app imports resolve
ROOT = os.path.join(PROJ, "input_data")
OUT = os.path.join(PROJ, "calibration_sessions", "v7")
DIAMETERS = (10, 12, 16, 20)
LDS = (5, 8, 10, 12, 15)
AMPS = (1, 2, 3, 4, 5, 6)
SEALED_D = 16
FY_NOMINAL = 540.0     # MPa; back-calculated from Kashani's tabulated lambda_p, one value for this dataset (lambda_p scales with sqrt(fy))
MIN_POINTS = 20        # a cycle with fewer samples is treated as truncated and skipped


class SealedError(RuntimeError):
    pass


def conditions(sealed: bool = False):
    """(D, L/D, amplitude %) of every condition in the development set, or in the sealed set if sealed=True (needs unseal)."""
    if sealed:
        _require_unseal()
    return [(d, ld, n) for d in DIAMETERS if (d == SEALED_D) == sealed for ld in LDS for n in AMPS]


def _require_unseal() -> None:
    if os.environ.get("V7_UNSEAL") != "1":
        raise SealedError("D = 16 mm is the sealed final test of V7. Set V7_UNSEAL=1 only after the framework is frozen.")


def _path(d: int, ld: int, n: int) -> str:
    return os.path.join(ROOT, f"{d}mm dia", f"LD{ld}", f"{n}data.csv")


def extract(d: int, ld: int, n: int) -> dict:
    """Read one test file and return its per-cycle features."""
    if d == SEALED_D:
        _require_unseal()
    cyc, pos, load = [], [], []
    with open(_path(d, ld, n), newline="") as f:
        r = csv.reader(f)
        idx = None
        for raw in r:
            if raw and raw[0].strip().strip('"').lower().startswith("notes"):
                h = [c.strip().strip('"').lower() for c in raw]
                idx = (h.index("cycle"), h.index("position mm"), h.index("load kn"))
                break
        for raw in r:
            try:
                cyc.append(int(float(raw[idx[0]]))); pos.append(float(raw[idx[1]])); load.append(float(raw[idx[2]]))
            except (ValueError, IndexError):
                continue
    cyc, pos, load = np.array(cyc), np.array(pos), np.array(load)
    a0 = math.pi * d * d / 4.0
    strain = (pos - pos[0]) / (ld * d)
    stress = (load - load[0]) * 1000.0 / a0
    cycles = []
    for c in range(1, int(cyc.max()) + 1):
        m = cyc == c
        if m.sum() < MIN_POINTS:
            continue
        cycles.append({"c": c, "t": float(stress[m].max()), "k": float(-stress[m].min()),
                       "w": float(abs(np.trapz(stress[m], strain[m]))), "np": int(m.sum())})
    return {"d": d, "ld": ld, "n": n, "nf": int(cyc.max()), "cycles": cycles}


def build_dev_features(force: bool = False) -> list:
    """Features of the 90 development files, cached in calibration_sessions/v7/dev_features.json."""
    path = os.path.join(OUT, "dev_features.json")
    if os.path.exists(path) and not force:
        return json.load(open(path))
    os.makedirs(OUT, exist_ok=True)
    feats = [extract(*c) for c in conditions(sealed=False)]
    with open(path, "w") as f:
        json.dump(feats, f)
    return feats


if __name__ == "__main__":
    fs = build_dev_features(force=True)
    print(len(fs), "development files;", sum(len(f["cycles"]) for f in fs), "cycles")
    try:
        extract(16, 5, 3)
    except SealedError as exc:
        print("seal check OK:", exc)
