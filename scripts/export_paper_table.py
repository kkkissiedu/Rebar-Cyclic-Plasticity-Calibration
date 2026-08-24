"""
export_paper_table.py
=====================

Build the paper's calibration-results table from session folders.

For each session: model, N, calibrated parameters, sigma_y0 (monotonic +
cyclic + Bauschinger ratio), surrogate objective, Stage-3 FE objective and
peaks, and -- critically for reviewers -- WHICH constraints were active at
the optimum (a parameter sitting on a bound or a physics-gate boundary is
not a free material constant and must be reported as constrained).

Active-constraint detection (relative tolerance 1%):
  * user search bounds (lower/upper) for every parameter and sigma_y0;
  * gamma-separation floor gamma_k/gamma_{k+1} = MIN_GAMMA_RATIO;
  * saturated-yield floor sigma_y0 + Q_inf = yield_fraction * sigma_y0;
  * physics-gate b window.

Usage (from the repo root):

    python scripts/export_paper_table.py                       # all sessions
    python scripts/export_paper_table.py calibration_sessions/<dir> [...]

Writes paper_table.md + paper_table.csv to the project root and prints the
markdown table.
"""

from __future__ import annotations

import csv
import json
import os
import sys

# repo root = parent of scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration_app.utils.physics_gate import (
    MIN_GAMMA_RATIO, DEFAULT_MIN_YIELD_FRACTION, MIN_B, MAX_B)

SESSIONS_ROOT = "calibration_sessions"
REL_TOL = 0.01


def _near(a: float, b: float, scale: float = None) -> bool:
    ref = abs(scale if scale is not None else b) or 1.0
    return abs(a - b) <= REL_TOL * ref


def active_constraints(state: dict) -> list[str]:
    """Human-readable list of constraints active at the session's optimum."""
    flags: list[str] = []
    best = state.get("best_params") or []
    cfg = state.get("config") or {}
    bounds = state.get("bounds") or {}
    n_bs = int(cfg.get("n_backstresses", 3))
    names = [x for k in range(1, n_bs + 1) for x in (f"C{k}", f"gamma{k}")]
    names += ["Q_inf", "b"]
    if len(best) != len(names):
        return ["param-vector length mismatch"]
    vals = dict(zip(names, [float(v) for v in best]))

    for n, v in vals.items():                     # search-bound activity
        if n in bounds:
            lo, hi = bounds[n]
            span = (hi - lo) or 1.0
            if _near(v, lo, span):
                flags.append(f"{n}@lower-bound({lo:g})")
            elif _near(v, hi, span):
                flags.append(f"{n}@upper-bound({hi:g})")

    for k in range(1, n_bs):                      # gamma-separation floor
        g1, g2 = vals[f"gamma{k}"], vals[f"gamma{k + 1}"]
        if g2 > 0 and _near(g1 / g2, MIN_GAMMA_RATIO, MIN_GAMMA_RATIO):
            flags.append(f"gamma{k}/gamma{k + 1}@separation-floor"
                         f"({MIN_GAMMA_RATIO:g})")

    sy0 = state.get("best_sigma_y0") or cfg.get("sigma_y0")
    if sy0:                                       # saturated-yield floor
        floor_q = (DEFAULT_MIN_YIELD_FRACTION - 1.0) * float(sy0)
        if _near(vals["Q_inf"], floor_q, abs(floor_q) or 1.0):
            flags.append(f"Q_inf@yield-floor({DEFAULT_MIN_YIELD_FRACTION:.0%}"
                         f" sigma_y0)")
        lo_sy = float(cfg.get("sy0_lower", 0.0))
        hi_sy = float(cfg.get("sigma_y0", 0.0))
        span = (hi_sy - lo_sy) or 1.0
        if _near(float(sy0), lo_sy, span):
            flags.append(f"sigma_y0@lower-bound({lo_sy:g})")
        elif _near(float(sy0), hi_sy, span):
            flags.append(f"sigma_y0@upper-bound({hi_sy:g})")

    if _near(vals["b"], MIN_B, MAX_B - MIN_B):    # gate b window
        flags.append(f"b@gate-min({MIN_B:g})")
    elif _near(vals["b"], MAX_B, MAX_B - MIN_B):
        flags.append(f"b@gate-max({MAX_B:g})")
    return flags


def row_from_session(path: str) -> dict | None:
    sj = os.path.join(path, "session.json")
    if not os.path.exists(sj):
        return None
    with open(sj, "r") as f:
        state = json.load(f)
    if not state.get("best_params"):
        return None
    cfg = state.get("config") or {}
    fe = state.get("fe_verification") or {}
    n_bs = int(cfg.get("n_backstresses", 3))
    names = [x for k in range(1, n_bs + 1) for x in (f"C{k}", f"gamma{k}")]
    names += ["Q_inf", "b"]
    params = dict(zip(names, state["best_params"]))
    row = {
        "session": os.path.basename(path),
        "data": os.path.basename(state.get("data_file", "?")),
        "material": state.get("material_description", ""),
        "model": state.get("model"),
        "N": n_bs,
        "backend": state.get("backend"),
        "optimiser": state.get("optimiser"),
        "n_fit": cfg.get("n_fit"),
    }
    for k in range(1, 5):
        row[f"C{k}"] = params.get(f"C{k}")
        row[f"gamma{k}"] = params.get(f"gamma{k}")
    row.update({
        "Q_inf": params.get("Q_inf"),
        "b": params.get("b"),
        "sigma_y0_monotonic": state.get("sigma_y0_monotonic"),
        "sigma_y0_cyclic": state.get("best_sigma_y0"),
        "bauschinger_ratio": state.get("bauschinger_ratio"),
        "surrogate_obj_MPa": state.get("best_objective"),
        "fe_obj_MPa": fe.get("objective_MPa"),
        "fe_peak_T_MPa": fe.get("peak_tension_MPa"),
        "fe_peak_C_MPa": fe.get("peak_compression_MPa"),
        "active_constraints": "; ".join(active_constraints(state)) or "none",
    })
    return row


def fmt(v) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def main() -> int:
    dirs = sys.argv[1:]
    if not dirs:
        dirs = [os.path.join(SESSIONS_ROOT, d)
                for d in sorted(os.listdir(SESSIONS_ROOT))
                if os.path.isdir(os.path.join(SESSIONS_ROOT, d))]
    rows = [r for r in (row_from_session(d) for d in dirs) if r]
    if not rows:
        print("No completed sessions found.")
        return 1

    cols = list(rows[0].keys())
    with open("paper_table.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    lines = ["| " + " | ".join(cols) + " |",
             "|" + "|".join("---" for _ in cols) + "|"]
    for r in rows:
        lines.append("| " + " | ".join(fmt(r[c]) for c in cols) + " |")
    md = "\n".join(lines) + "\n"
    with open("paper_table.md", "w") as f:
        f.write("# Calibration results (auto-generated by "
                "export_paper_table.py)\n\n"
                "Parameters flagged in `active_constraints` sit on a search "
                "bound or physics-gate boundary at the optimum and must be "
                "reported as constrained, not free.\n\n" + md)
    print(md)
    print("Wrote paper_table.md + paper_table.csv "
          f"({len(rows)} session(s)).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
