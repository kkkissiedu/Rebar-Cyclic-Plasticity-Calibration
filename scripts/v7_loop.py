"""
v7_loop.py
==========

V7 loop-shape layer (L3), development set only (V7_PROTOCOL.md section 3).

* Machine compliance as one physical parameter: the LVDT strain includes the deflection of a load train of stiffness k_m (N/mm), so the apparent
  modulus is  1 / E_app = 1 / E_bar + (pi D / (4 L/D)) / k_m  (A0/L0 = pi D / (4 L/D) in mm). With an extensometer k_m drops out.
* One joint Chaboche (N = 2) core fitted to the cycle-2 loops of all training files at once (UVC practice: several protocols jointly), with E_app
  from the compliance model. Parameters: k_m, sigma_0, C1, gamma1, C2, gamma2, Q, b.
* Loops at any cycle are the cycle-2 unit loop mapped on the observed strain path by arc-length fraction and scaled affinely so its tension and
  compression peaks equal the predicted peaks of the degradation layer.

    python v7_loop.py            # joint fit on all development files and diagnostics
"""

from __future__ import annotations

import csv
import json
import math
import os
import time

import numpy as np
from scipy.optimize import least_squares

import v7_data as D
from calibration_app.core.pipeline import build_model

E_BAR = 200_000.0
MODEL = build_model("chaboche", 2)          # parameters (C1, gamma1, C2, gamma2, Q_inf, b)
TARGET_X = (0.1, 0.3, 0.5, 0.7, 0.9)        # life fractions of the representative cycles used for loop metrics
NFIT = 80
PTS = 200
S_GRID = np.linspace(0.0, 1.0, NFIT)


def tri_path(eps: float, cycles: int = 2, pts: int = PTS) -> np.ndarray:
    one = np.concatenate([np.linspace(0, eps, pts // 4 + 1)[:-1], np.linspace(eps, -eps, pts // 2 + 1)[:-1],
                          np.linspace(-eps, 0, pts // 4 + 1)[:-1]])
    return np.concatenate([np.tile(one, cycles), [0.0]])


_P = tri_path(1.0)[PTS:2 * PTS + 1]
S_IDEAL = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(_P)))]) / np.abs(np.diff(_P)).sum()


def arc_fraction(e: np.ndarray) -> np.ndarray:
    s = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(e)))])
    return s / s[-1]


def read_loops(d: int, ld: int, n: int, nf: int) -> dict:
    """Cycle 2 and the representative cycles of one file as {cycle: (strain, stress)} (engineering, initial load subtracted)."""
    if d == D.SEALED_D:
        D._require_unseal()
    cyc, pos, load = [], [], []
    with open(D._path(d, ld, n), newline="") as f:
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
    strain = (pos - pos[0]) / (ld * d)
    stress = (load - load[0]) * 1000.0 / (math.pi * d * d / 4.0)
    ok = [c for c in range(2, nf) if (cyc == c).sum() >= D.MIN_POINTS]
    if 2 not in ok:
        return {}
    pick = {2} | {min(ok, key=lambda c: abs(c / nf - x)) for x in TARGET_X}
    return {int(c): (strain[cyc == c], stress[cyc == c]) for c in sorted(pick)}


def build_dev_loops(force: bool = False) -> list:
    path = os.path.join(D.OUT, "dev_loops.json")
    if os.path.exists(path) and not force:
        raw = json.load(open(path))
    else:
        os.makedirs(D.OUT, exist_ok=True)
        raw = []
        for f in D.build_dev_features():
            lp = read_loops(f["d"], f["ld"], f["n"], f["nf"])
            if lp:
                raw.append({"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": f["nf"],
                            "loops": {str(c): {"e": np.round(e, 6).tolist(), "s": np.round(s, 2).tolist()} for c, (e, s) in lp.items()}})
        json.dump(raw, open(path, "w"))
    files = []
    for f in raw:
        loops = {}
        for c, v in f["loops"].items():
            e, s = np.array(v["e"]), np.array(v["s"])
            sa = arc_fraction(e)
            loops[int(c)] = {"e": e, "s": s, "arc": sa, "y": np.interp(S_GRID, sa, s), "t": float(s.max()), "k": float(-s.min()),
                             "w": float(abs(np.trapz(s, e)))}
        files.append({"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": f["nf"], "loops": loops})
    return files


def e_app(d: float, ld: float, km: float) -> float:
    return 1.0 / (1.0 / E_BAR + (math.pi * d / (4.0 * ld)) / km)


NAMES = ("ln k_m", "sigma_0", "ln C1", "ln gamma1", "ln C2", "ln gamma2", "Q", "ln b")
TH0 = np.array([math.log(1.6e5), 250.0, math.log(2e4), math.log(200.0), math.log(2e3), math.log(10.0), 50.0, math.log(5.0)])
LB = np.array([math.log(2e4), 20.0, math.log(500.0), math.log(20.0), math.log(100.0), math.log(1.0), 0.0, math.log(0.05)])
UB = np.array([math.log(1e7), 600.0, math.log(1e5), math.log(2000.0), math.log(2e4), math.log(200.0), 300.0, math.log(50.0)])


def unit_loop(theta: np.ndarray, d: float, ld: float, eps: float) -> np.ndarray:
    """Cycle-2 stress response (on the ideal strain path of amplitude eps) of the joint Chaboche core."""
    p = (math.exp(theta[2]), math.exp(theta[3]), math.exp(theta[4]), math.exp(theta[5]), theta[6], math.exp(theta[7]))
    sig = MODEL.simulate(tri_path(eps), p, E=e_app(d, ld, math.exp(theta[0])), sy0=theta[1])
    return np.asarray(sig, float)[PTS:2 * PTS + 1]


def residuals(theta: np.ndarray, files: list) -> np.ndarray:
    out = []
    for f in files:
        sig = unit_loop(theta, f["d"], f["ld"], f["n"] / 100.0)
        if not np.all(np.isfinite(sig)):
            out.append(np.full(NFIT, 1.0e3))
            continue
        out.append((np.interp(S_GRID, S_IDEAL, sig) - f["loops"][2]["y"]) / math.sqrt(NFIT))
    return np.concatenate(out)


def fit_joint(files: list, th0: np.ndarray = TH0, max_nfev: int = 80) -> np.ndarray:
    return least_squares(residuals, th0, bounds=(LB, UB), args=(files,), max_nfev=max_nfev, x_scale="jac").x


def scale_to_peaks(unit: np.ndarray, t: float, k: float) -> np.ndarray:
    tm, km = float(unit.max()), float(-unit.min())
    a = (t + k) / (tm + km)
    return a * unit + (t - a * tm)


def loop_on_obs(theta: np.ndarray, f: dict, loop: dict, t: float, k: float) -> np.ndarray:
    """Predicted stress at the observed samples of a loop: the unit loop scaled to peaks (t, k), mapped by arc fraction."""
    unit = unit_loop(theta, f["d"], f["ld"], f["n"] / 100.0)
    return np.interp(loop["arc"], S_IDEAL, scale_to_peaks(unit, t, k))


def main() -> int:
    t0 = time.time()
    files = build_dev_loops(force=not os.path.exists(os.path.join(D.OUT, "dev_loops.json")))
    print(f"development files with loops: {len(files)}; representative loops: {sum(len(f['loops']) for f in files)}")
    th = fit_joint(files)
    print(f"joint fit in {time.time() - t0:.0f} s")
    for n_, v in zip(NAMES, th):
        print(f"  {n_:10s} {v:9.3f}   ({math.exp(v):.4g})" if n_.startswith("ln") else f"  {n_:10s} {v:9.3f}")
    km = math.exp(th[0])
    print(f"machine stiffness k_m = {km / 1000:.0f} kN/mm; apparent E (GPa): " +
          ", ".join(f"D{d} LD{ld}: {e_app(d, ld, km) / 1000:.0f}" for d, ld in ((10, 5), (12, 5), (20, 5), (12, 15))))
    r = residuals(th, files) * math.sqrt(NFIT)
    per_file = np.sqrt((r.reshape(len(files), NFIT) ** 2).mean(axis=1))
    print(f"cycle-2 loop RMSE of the joint core, per file: median {np.median(per_file):.1f} MPa, 90th pct {np.percentile(per_file, 90):.1f}, max {per_file.max():.1f}")
    # shape-only error on the representative cycles: unit loop scaled to the OBSERVED peaks
    shape = [float(np.sqrt(np.mean((loop_on_obs(th, f, lp, lp["t"], lp["k"]) - lp["s"]) ** 2)))
             for f in files for c, lp in f["loops"].items()]
    flat = [float(np.sqrt(np.mean((np.full_like(lp["s"], lp["s"].mean()) - lp["s"]) ** 2))) for f in files for lp in f["loops"].values()]
    print(f"shape-only error on {len(shape)} representative loops (observed peaks given): median {np.median(shape):.1f} MPa, "
          f"90th pct {np.percentile(shape, 90):.1f}; stress scale (loop sd) median {np.median(flat):.0f} MPa")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
