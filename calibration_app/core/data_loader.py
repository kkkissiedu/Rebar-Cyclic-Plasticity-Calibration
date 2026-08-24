"""
core/data_loader.py
===================

Load experimental low-cycle-fatigue CSVs and prepare calibration windows.

Handles the Kashani/Bristol dataset layout used for pipeline validation and, by
design, any specimen the same rig produces (all diameters x all L/D ratios) plus
the future Ghanaian scrap-rebar data. Nothing here is B500C-specific.

Dataset conventions (Kashani/Bristol; readme in input_data/)
------------------------------------------------------------
* Directory: ``input_data/{D}mm dia/LD{ratio}/{N}data.csv``
* Units are kN and mm ("read 'Volts' as mm").
* File name ``{N}data.csv`` -> N = nominal strain amplitude in percent.
* **Gauge length is set by the folder, not the file name**: gauge = L/D x
  diameter (e.g. 12 mm dia / LD5 -> 60 mm). This is why the strain amplitude in
  mm is ``N/100 * gauge`` and depends on both diameter and L/D.
* Stress = load_kN * 1000 / A, with A = pi (D/2)^2.
* Strain = (Position_mm - ref_position) / gauge, ref = first row of cycle 2.
* Cycle 1 is machine ramp-up: excluded from scoring but kept as a burn-in path
  so the surrogate/FE state is conditioned before cycle 2.

Yield stress is a **per-run input** (companion monotonic CSV or manual entry),
never hardcoded — the material grade is unknown for scrap rebar. This module
only ever reads the data files; it never writes or modifies them.
"""

from __future__ import annotations

import csv
import math
import os
import re
from typing import Optional

import numpy as np

#: Matches the Kashani file-name convention, capturing the integer percent.
_STRAIN_RE = re.compile(r"(\d+)\s*data\.csv$", re.IGNORECASE)
#: Matches a diameter folder like "12mm dia".
_DIA_RE = re.compile(r"(\d+)\s*mm\s*dia", re.IGNORECASE)
#: Matches an L/D folder like "LD5".
_LD_RE = re.compile(r"LD\s*(\d+)", re.IGNORECASE)


def detect_strain_amplitude_pct(path: str) -> Optional[float]:
    """Return the nominal strain amplitude in percent from ``{N}data.csv``.

    Returns ``None`` if the file name does not match the convention (the caller
    should then warn and require a manual amplitude).
    """
    m = _STRAIN_RE.search(os.path.basename(path))
    return float(m.group(1)) if m else None


def detect_geometry(path: str) -> dict:
    """Infer (diameter_mm, ld_ratio, gauge_length_mm) from the folder path.

    Looks for ``{D}mm dia`` and ``LD{ratio}`` components anywhere in the path.
    Returns a dict with whatever could be found; missing keys are ``None`` so the
    caller can fall back to explicit configuration.
    """
    norm = path.replace("\\", "/")
    dia = _DIA_RE.search(norm)
    ld = _LD_RE.search(norm)
    diameter = float(dia.group(1)) if dia else None
    ld_ratio = float(ld.group(1)) if ld else None
    gauge = (ld_ratio * diameter) if (diameter and ld_ratio) else None
    return {"diameter_mm": diameter, "ld_ratio": ld_ratio,
            "gauge_length_mm": gauge}


def bar_area_mm2(diameter_mm: float) -> float:
    """Cross-sectional area of a round bar (mm^2)."""
    return math.pi * (diameter_mm / 2.0) ** 2


def load_monotonic_yield(
    path: str,
    diameter_mm: float,
    offset_strain: float = 0.002,
) -> Optional[float]:
    """Estimate yield stress (MPa) from a companion monotonic tension CSV.

    Accepts either (strain, stress) columns or the same (Position mm, Load kN)
    layout as the cyclic files (converted with ``diameter_mm``). Uses the 0.2%
    (``offset_strain``) offset method if an elastic slope can be estimated;
    otherwise returns ``None`` so the caller falls back to manual entry.

    This is a best-effort helper — for unknown scrap-rebar grades the user can
    always type the yield stress in directly.
    """
    try:
        strain, stress = _read_monotonic(path, diameter_mm)
    except Exception:
        return None
    if strain is None or strain.size < 5:
        return None
    # Elastic modulus from the initial ~0.05% strain portion.
    elastic_mask = strain <= max(0.0005, strain.max() * 0.02)
    if np.count_nonzero(elastic_mask) < 2:
        return None
    E = float(np.polyfit(strain[elastic_mask], stress[elastic_mask], 1)[0])
    if not math.isfinite(E) or E <= 0:
        return None
    offset_line = E * (strain - offset_strain)
    diff = stress - offset_line
    sign_change = np.where(np.diff(np.sign(diff)) != 0)[0]
    if sign_change.size == 0:
        return None
    i = sign_change[0]
    return float(0.5 * (stress[i] + stress[i + 1]))


def _read_monotonic(path: str, diameter_mm: float):
    """Return (strain, stress) arrays from a monotonic CSV, or (None, None)."""
    with open(path, "r", newline="") as f:
        rows = [r for r in csv.reader(f) if r]
    if not rows:
        return None, None
    header = [c.strip().strip('"').lower() for c in rows[0]]

    def col(*names):
        for n in names:
            if n in header:
                return header.index(n)
        return None

    i_strain = col("strain", "strain %", "e")
    i_stress = col("stress", "stress_mpa", "stress mpa")
    i_pos = col("position mm", "displacement", "position")
    i_load = col("load kn", "load", "force kn")

    strain_vals, stress_vals = [], []
    for raw in rows[1:]:
        try:
            if i_strain is not None and i_stress is not None:
                s = float(raw[i_strain]); sig = float(raw[i_stress])
            elif i_pos is not None and i_load is not None:
                s = float(raw[i_pos]); sig = float(raw[i_load])
            else:
                return None, None
        except (ValueError, IndexError):
            continue
        strain_vals.append(s)
        stress_vals.append(sig)
    if not strain_vals:
        return None, None
    strain = np.asarray(strain_vals, float)
    stress = np.asarray(stress_vals, float)
    # If we read Position/Load, convert load to stress (MPa). Position is left
    # as-is; the offset method operates on relative shape.
    if i_strain is None:
        area = bar_area_mm2(diameter_mm)
        stress = stress * 1000.0 / area          # kN -> N -> MPa
    return strain, stress


# --- monotonic-CSV column detection + curve extraction (grade-agnostic) -----

#: Header-token patterns for the force and displacement columns (case-insens.).
_FORCE_PAT = re.compile(r"\b(kn|force|load|newton|n)\b|force|load", re.IGNORECASE)
_DISP_PAT = re.compile(r"\b(mm|disp|displacement|extension|position|stroke|volt)\b|displacement|extension", re.IGNORECASE)


def monotonic_columns(path: str) -> dict:
    """Best-effort auto-detection of the force/displacement columns.

    Returns ``{"header": [...], "force_idx": int|None, "disp_idx": int|None}``.
    When either index is ``None`` the caller should show a column-mapping dialog.
    """
    with open(path, "r", newline="") as f:
        for raw in csv.reader(f):
            if raw and any(cell.strip() for cell in raw):
                header = [c.strip().strip('"') for c in raw]
                break
        else:
            return {"header": [], "force_idx": None, "disp_idx": None}
    force_idx = disp_idx = None
    for i, name in enumerate(header):
        low = name.lower()
        if force_idx is None and _FORCE_PAT.search(low):
            force_idx = i
        elif disp_idx is None and _DISP_PAT.search(low):
            disp_idx = i
    return {"header": header, "force_idx": force_idx, "disp_idx": disp_idx}


def monotonic_curve(
    path: str,
    diameter_mm: float,
    gauge_length_mm: float,
    force_idx: int,
    disp_idx: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Read a monotonic tension CSV into (strain, stress[MPa]) arrays.

    Stress = force_kN * 1000 / area(diameter); strain = displacement / gauge.
    Grade-agnostic: only geometry (from the LCF file) is used, no assumed grade.
    """
    area = bar_area_mm2(diameter_mm)
    disp_vals, force_vals = [], []
    with open(path, "r", newline="") as f:
        rows = list(csv.reader(f))
    started = False
    for raw in rows:
        if not raw or len(raw) <= max(force_idx, disp_idx):
            continue
        try:
            d = float(raw[disp_idx]); fkn = float(raw[force_idx])
        except (ValueError, IndexError):
            if started:
                continue
            continue                              # skip header / preamble
        started = True
        disp_vals.append(d); force_vals.append(fkn)
    disp = np.asarray(disp_vals, float)
    force = np.asarray(force_vals, float)
    strain = disp / gauge_length_mm
    stress = force * 1000.0 / area
    return strain, stress


def yield_from_curve(strain: np.ndarray, stress: np.ndarray,
                     offset_strain: float = 0.002):
    """0.2%-offset yield from a monotonic curve.

    Returns ``(sy0, E, (offset_x, offset_y))`` for plotting, or
    ``(None, None, None)`` if an elastic slope / intercept can't be found.
    """
    if strain is None or strain.size < 5:
        return None, None, None
    elastic_mask = strain <= max(0.0005, float(strain.max()) * 0.02)
    if np.count_nonzero(elastic_mask) < 2:
        return None, None, None
    E = float(np.polyfit(strain[elastic_mask], stress[elastic_mask], 1)[0])
    if not math.isfinite(E) or E <= 0:
        return None, None, None
    offset_y = E * (strain - offset_strain)
    diff = stress - offset_y
    crossings = np.where(np.diff(np.sign(diff)) != 0)[0]
    if crossings.size == 0:
        return None, E, (strain, offset_y)
    i = crossings[0]
    sy0 = float(0.5 * (stress[i] + stress[i + 1]))
    return sy0, E, (strain, offset_y)


def load_experimental(
    path: str,
    *,
    gauge_length_mm: float,
    bar_diameter_mm: float,
    sigma_y0: Optional[float] = None,
    ref_cycle: int = 2,
) -> dict:
    """Load a cyclic LCF CSV and build the calibration data dictionary.

    Parameters
    ----------
    path : str
        Path to a ``{N}data.csv`` file (Kashani layout).
    gauge_length_mm, bar_diameter_mm : float
        Specimen geometry (from :func:`detect_geometry` or explicit config).
    sigma_y0 : float or None
        Per-run yield stress. Stored on the returned dict for the physics gate;
        may be ``None`` here and set later from a monotonic CSV / manual entry.
    ref_cycle : int
        Cycle whose first row defines the zero-strain reference position and the
        first scored cycle (cycle 1 is treated as machine ramp / burn-in).

    Returns
    -------
    dict with keys: ``time, cycle, strain, stress`` (1-D arrays, cycles >=
    ref_cycle only), ``burnin_strain`` (cycle 1), ``ref_position``,
    ``cycles_present``, ``measured_peak_stress``, ``strain_amplitude_pct``,
    ``amp_mm``, ``gauge_length_mm``, ``bar_diameter_mm``, ``sigma_y0``,
    ``data_file``.
    """
    times, cycles, positions, loads = [], [], [], []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = None
        for raw in reader:
            if not raw:
                continue
            joined = ",".join(s.strip() for s in raw).lower()
            if joined.startswith("notes") or joined.startswith('"notes"'):
                header = raw
                break
        if header is None:
            raise RuntimeError(f"Header row (containing 'Notes') not found in {path}")
        clean = [c.strip().strip('"').lower() for c in header]
        try:
            i_t = clean.index("time")
            i_c = clean.index("cycle")
            i_p = clean.index("position mm")
            i_l = clean.index("load kn")
        except ValueError as exc:
            raise RuntimeError(f"Expected columns missing in {path}: {exc}")

        for raw in reader:
            if not raw or len(raw) <= max(i_t, i_c, i_p, i_l):
                continue
            try:
                t = float(raw[i_t]); c = int(float(raw[i_c]))
                pos = float(raw[i_p]); ld = float(raw[i_l])
            except ValueError:
                continue                         # date-stamped first row etc.
            times.append(t); cycles.append(c)
            positions.append(pos); loads.append(ld)

    if not cycles:
        raise RuntimeError(f"No numeric rows parsed from {path}")

    time_arr = np.asarray(times, float)
    cycle_arr = np.asarray(cycles, np.int32)
    pos_arr = np.asarray(positions, float)
    load_arr = np.asarray(loads, float)

    mask_ref = cycle_arr == ref_cycle
    if not np.any(mask_ref):
        raise RuntimeError(f"Reference cycle {ref_cycle} not present in {path}")
    ref_position = float(pos_arr[mask_ref][0])

    area = bar_area_mm2(bar_diameter_mm)

    # Burn-in: cycle 1 strain path (kept for state conditioning, not scored).
    mask_ramp = cycle_arr < ref_cycle
    burnin_strain = ((pos_arr[mask_ramp] - ref_position) / gauge_length_mm
                     if np.any(mask_ramp) else np.array([], float))

    keep = cycle_arr >= ref_cycle
    strain = (pos_arr[keep] - ref_position) / gauge_length_mm
    stress = load_arr[keep] * 1000.0 / area
    cycle_kept = cycle_arr[keep]

    strain_pct = detect_strain_amplitude_pct(path)
    amp_mm = (strain_pct / 100.0 * gauge_length_mm) if strain_pct is not None else None

    # Cycle bookkeeping read entirely from the data (never hardcoded), so the
    # UI can size N_fit, the cycle selector and the .inp cycle count to the file.
    all_cycles = [int(c) for c in np.unique(cycle_arr).tolist()]
    scoreable = [int(c) for c in np.unique(cycle_kept).tolist()
                 if np.count_nonzero(cycle_kept == c) >= 4]
    cycle_summary = {
        "total_cycles": len(all_cycles),
        "first_scoreable": int(ref_cycle),
        "last_cycle": int(max(all_cycles)) if all_cycles else 0,
        "cycles_with_full_data": scoreable,
    }

    return {
        "time": time_arr[keep],
        "cycle": cycle_kept,
        "strain": strain,
        "stress": stress,
        "burnin_strain": burnin_strain,
        "ref_position": ref_position,
        "cycles_present": np.unique(cycle_kept).tolist(),
        "cycle_summary": cycle_summary,
        "measured_peak_stress": float(np.max(np.abs(stress))) if stress.size else 0.0,
        "strain_amplitude_pct": strain_pct,
        "amp_mm": amp_mm,
        "gauge_length_mm": gauge_length_mm,
        "bar_diameter_mm": bar_diameter_mm,
        "sigma_y0": sigma_y0,
        "data_file": path,
    }


def slice_window(data: dict, n_fit: int) -> dict:
    """Slice to the first ``n_fit`` scored cycles for the objective.

    Returns a window dict with ``strain, stress, cycle, target_cycles`` and the
    carried-through ``burnin_strain`` (see :func:`load_experimental`).
    """
    cycles_present = data["cycles_present"]
    if not cycles_present:
        raise RuntimeError("No cycles in data to slice")
    target = cycles_present[:n_fit]
    mask = np.isin(data["cycle"], target)
    return {
        "strain": data["strain"][mask],
        "stress": data["stress"][mask],
        "cycle": data["cycle"][mask],
        "target_cycles": target,
        "burnin_strain": data.get("burnin_strain", np.array([], float)),
    }
