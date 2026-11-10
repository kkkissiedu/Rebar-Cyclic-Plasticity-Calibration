"""
framework.py
============

Predictive framework for steel reinforcing bars under constant-amplitude strain cycling (V7 of the research; see docs/V7_PROTOCOL.md and docs/V7_EVIDENCE.md).
Given the diameter D, the slenderness L/D, the strain amplitude and (optionally) the yield strength fy, it predicts:

  1. the cycles to failure Nf                    eps_a = eps_f (2 Nf)^alpha                         (Koh-Stephens form)
        ln eps_f = c0 + c1 lp + c2 D,   alpha = a0 + a1 lp + a2 D,   lp = (L/D) sqrt(fy / 100)       (buckling parameter, Dhakal and Maekawa)
  2. the peak tension and compression stress of cycle 2     ln sigma_2 = b0 + b1 ln(eps_a in %) + b2 lp
  3. the decay of the peaks with the life fraction x = n / Nf   sigma(n) = sigma_2 g(x),
        g = exp(-A x^p) - B x^m   (A, B linear in 1, L/D / 10, eps_a / 6 %; p, m fitted)       (empirical master curve of this work)
  4. the loop of any cycle: a Chaboche (N = 2) material core, integrated numerically, whose elastic slope includes the machine compliance
        1 / E_app = 1 / E_bar + (pi D / (4 L/D)) / k_m,
     distorted by buckling (compression capacity B (1 + u / u0)^(-q), and a tension branch that straightens back to the peak, with parameters
     linear in 1, ln lp, ln eps_a and x), scaled to the predicted peaks, and averaged over the log-normal scatter of life (ASTM E739 treatment).

Fitted on the Kashani et al. (2019) B500C data (D 10 to 20 mm, L/D 5 to 15, amplitude 1 to 6 %, strain from the actuator LVDT). Outside that range, or for another
steel, the results are extrapolations; the material factors f_s (stress level) and f_e (fatigue ductility) exist for rescaling from a few tests, but their use on a
different real steel is untested.
"""

from __future__ import annotations

import json
import math
import os
from typing import Optional, Sequence

import numpy as np

from ..core.pipeline import build_model

PTS = 200
I0, I1, I2 = PTS // 4, 3 * PTS // 4, PTS
_PARAMS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "predictor_params.json")
_ZQ, _WQ = np.polynomial.hermite.hermgauss(9)
_GLB = np.array([-5.0, 0.05, -6.0, math.log(0.005), math.log(0.3)])
_GUB = np.array([5.0, 3.0, 6.0, math.log(2.0), math.log(8.0)])
_MODEL = None
_CACHE: dict = {}


def load_params(path: Optional[str] = None) -> dict:
    with open(path or _PARAMS_FILE, encoding="utf-8") as f:
        return json.load(f)


def _model():
    global _MODEL
    if _MODEL is None:
        _MODEL = build_model("chaboche", 2)
    return _MODEL


def tri_path(eps: float, cycles: int = 2, pts: int = PTS) -> np.ndarray:
    one = np.concatenate([np.linspace(0, eps, pts // 4 + 1)[:-1], np.linspace(eps, -eps, pts // 2 + 1)[:-1],
                          np.linspace(-eps, 0, pts // 4 + 1)[:-1]])
    return np.concatenate([np.tile(one, cycles), [0.0]])


def apparent_modulus(d: float, ld: float, km: float, e_bar: float = 200000.0) -> float:
    return 1.0 / (1.0 / e_bar + (math.pi * d / (4.0 * ld)) / km)


def _lambda_p(ld: float, fy: float) -> float:
    return ld * math.sqrt(fy / 100.0)


def _life(p: dict, d: float, ld: float, eps_pct: float, fy: float, f_e: float = 1.0) -> float:
    lam = _lambda_p(ld, fy)
    x = np.array([1.0, lam] + ([d] if p["life"]["variant"] == "lam_d" else []))
    k = x.size
    c = float(x @ np.array(p["life"]["coef"][:k])) + math.log(f_e)
    a = min(float(x @ np.array(p["life"]["coef"][k:])), p["constants"]["alpha_max"])
    return math.exp((math.log(eps_pct / 100.0) - c) / a) / 2.0


def _peak2(p: dict, d: float, ld: float, eps_pct: float, fy: float) -> tuple:
    cols = [1.0, math.log(eps_pct)]
    if p["peak2"]["variant"] in ("eps_lam", "eps_lam_d"):
        cols.append(_lambda_p(ld, fy))
    if p["peak2"]["variant"] == "eps_lam_d":
        cols.append(d)
    x = np.array(cols)
    return math.exp(float(x @ np.array(p["peak2"]["coef_t"]))), math.exp(float(x @ np.array(p["peak2"]["coef_k"])))


def _gcurve(x, a, b, s1, s2, variant):
    x = np.clip(x, 0.0, 1.0)
    if variant == "power_drop":
        return 1.0 - a * x ** s1 - b * x ** s2
    if variant == "weibull_drop":
        return np.exp(-a * x ** s1) - b * x ** s2
    return 1.0 - a * (1.0 - np.exp(-s1 * x)) - b * (np.exp(s2 * x) - 1.0) / (np.exp(s2) - 1.0)


def _decay(p: dict, ld: float, eps_pct: float, cycles, nf: float, which: str):
    th = np.array(p["decay"]["theta_" + which])
    cov = np.array([1.0, ld / 10.0, eps_pct / 6.0])
    x = np.minimum(np.asarray(cycles, float) / max(nf, 2.0), 1.0)
    return np.clip(_gcurve(x, float(cov @ th[:3]), float(cov @ th[3:6]), th[6], th[7], p["decay"]["variant"]), 0.0, 1.5)


def _unit_loop(p: dict, d: float, ld: float, eps: float, strain_basis: str) -> np.ndarray:
    th = np.array(p["core"]["theta"])
    key = (round(d, 6), round(ld, 6), round(eps, 9), strain_basis)
    if key not in _CACHE:
        e = p["constants"]["e_bar"] if strain_basis == "extensometer" else apparent_modulus(d, ld, math.exp(th[0]), p["constants"]["e_bar"])
        pr = (math.exp(th[2]), math.exp(th[3]), math.exp(th[4]), math.exp(th[5]), th[6], math.exp(th[7]))
        _CACHE[key] = np.asarray(_model().simulate(tri_path(eps), pr, E=e, sy0=th[1]), float)[PTS:2 * PTS + 1]
    return _CACHE[key]


def _scale_to_peaks(unit: np.ndarray, t: float, k: float) -> np.ndarray:
    tm, km = float(unit.max()), float(-unit.min())
    a = (t + k) / (tm + km)
    return a * unit + (t - a * tm)


def _smax(a, b, delta):
    return 0.5 * (a + b + np.sqrt((a - b) ** 2 + delta * delta))


def _distort(unit: np.ndarray, eps: float, g: np.ndarray, t: float, k: float, delta_frac: float) -> np.ndarray:
    """Buckling distortion of the material-core loop (compression capacity envelope and tension straightening branch)."""
    g = np.clip(g, _GLB, _GUB)
    tc, kc = float(unit.max()), float(-unit.min())
    a = t / tc
    if k >= a * kc * 0.999:                                   # no buckling: the core loop, affinely scaled
        return _scale_to_peaks(unit, t, k)
    u0, q = math.exp(g[0]) * eps, g[1]
    w, ze, m = 1.0 / (1.0 + math.exp(-g[2])), math.exp(g[3]), math.exp(g[4])
    sc = a * unit[I0:I1 + 1]
    u = np.linspace(0.0, 2.0 * eps, sc.size)
    neg = np.nonzero(sc <= 0.0)[0]
    j = int(neg[0]) if neg.size else sc.size - 1
    uz = u[j - 1] + (u[j] - u[j - 1]) * sc[j - 1] / (sc[j - 1] - sc[j]) if j > 0 and sc[j - 1] > sc[j] else u[j]
    decay = (1.0 + np.maximum(u - uz, 0.0) / u0) ** (-q)
    delta = delta_frac * k
    lo, hi = 0.2 * k, 1.5 * a * kc + k
    for _ in range(30):                                       # strength of the mechanism solved so that the compression peak equals k
        mid = 0.5 * (lo + hi)
        if -_smax(sc, -mid * decay, delta).min() > k:
            hi = mid
        else:
            lo = mid
    sb = _smax(sc, -0.5 * (lo + hi) * decay, delta)
    sb = sb + (t - sb[0]) * np.linspace(1.0, 0.0, sb.size) ** 8
    sr = float(sb[-1])
    z = np.linspace(0.0, 1.0, 2 * (I2 - I1) + 1)
    h = w * (1.0 - np.exp(-z / ze)) / (1.0 - math.exp(-1.0 / ze)) + (1.0 - w) * z ** m
    nh = I2 - I1
    seg_a = sr + (t - sr) * h[nh:]
    seg_c = sr + (t - sr) * h[:nh + 1]
    return np.concatenate([seg_a[:-1], sb, seg_c[1:]])


def _g_params(p: dict, ld: float, eps_pct: float, x: float, fy: float) -> np.ndarray:
    phi = np.array([1.0, math.log(_lambda_p(ld, fy) / p["constants"]["lambda_ref"]), math.log(eps_pct / 3.0), x])
    return np.array(p["buckling"]["params"]).reshape(5, -1) @ phi


def predict_bar(diameter_mm: float, slenderness: float, amplitude_pct: float, fy_mpa: float = 540.0,
                cycles: Optional[Sequence[int]] = None, strain_basis: str = "lvdt", factors: tuple = (1.0, 1.0),
                params: Optional[dict] = None) -> dict:
    """Predict the life, the peak-stress history and the loops of one bar. factors = (f_s, f_e), both 1 for the fitted steel."""
    p = params or load_params()
    d, ld, ea = float(diameter_mm), float(slenderness), float(amplitude_pct)
    f_s, f_e = factors
    eps = ea / 100.0
    warnings = []
    for name, v, lo, hi in (("diameter", d, 10, 20), ("L/D", ld, 5, 15), ("strain amplitude", ea, 1, 6)):
        if not lo <= v <= hi:
            warnings.append(f"{name} {v:g} is outside the fitted range {lo} to {hi}: extrapolation.")
    if abs(fy_mpa - p["constants"]["fy_nominal"]) > 1e-9:
        warnings.append("fy enters only through the buckling parameter lp; its effect is untested on the data (one steel).")
    if abs(f_s - 1.0) > 1e-9 or abs(f_e - 1.0) > 1e-9:
        warnings.append(f"Material factors f_s = {f_s:g}, f_e = {f_e:g} applied: derived from monotonic tests only and untested against cyclic tests of that steel.")
    if strain_basis == "extensometer":
        warnings.append("Extensometer basis removes the fitted machine compliance (E = 200 GPa); this is untested because the data have no extensometer strain.")
    nf = _life(p, d, ld, ea, fy_mpa, f_e)
    s = p["life"]["s_log10"]
    t2, k2 = (v * f_s for v in _peak2(p, d, ld, ea, fy_mpa))
    n_last = max(3, int(round(nf)))
    n_grid = np.arange(2, n_last + 1)
    peaks = {"cycle": n_grid.tolist(), "tension": (t2 * _decay(p, ld, ea, n_grid, nf, "t")).tolist(),
             "compression": (k2 * _decay(p, ld, ea, n_grid, nf, "k")).tolist()}
    if cycles is None:
        cycles = sorted({2} | {int(min(max(round(fr * nf), 2), max(n_last - 1, 2))) for fr in (0.25, 0.5, 0.75, 0.9)})
    unit = _unit_loop(p, d, ld, eps, strain_basis)
    strain_grid = tri_path(eps)[PTS:2 * PTS + 1]
    loops = {}
    for c in cycles:
        avg = np.zeros(unit.shape, float)
        for z, w in zip(_ZQ, _WQ):
            nfi = nf * 10.0 ** (s * math.sqrt(2.0) * z)
            ti, ki = (t2, k2) if c == 2 else (t2 * float(_decay(p, ld, ea, [c], nfi, "t")[0]), k2 * float(_decay(p, ld, ea, [c], nfi, "k")[0]))
            avg += w / math.sqrt(math.pi) * _distort(unit, eps, _g_params(p, ld, ea, c / nfi, fy_mpa), ti, ki, p["constants"]["delta"])
        loops[int(c)] = {"strain": strain_grid, "stress": avg, "tension_peak": float(avg.max()), "compression_peak": float(-avg.min())}
    th = p["core"]["theta"]
    return {"inputs": {"diameter_mm": d, "L_over_D": ld, "amplitude_pct": ea, "fy_mpa": fy_mpa, "strain_basis": strain_basis, "f_s": f_s, "f_e": f_e},
            "lambda_p": _lambda_p(ld, fy_mpa), "apparent_modulus_gpa": (p["constants"]["e_bar"] if strain_basis == "extensometer"
                                                                     else apparent_modulus(d, ld, math.exp(th[0]), p["constants"]["e_bar"])) / 1000.0,
            "life": {"cycles_to_failure": nf, "factor2_band": (nf / 2.0, nf * 2.0), "log10_scatter": s},
            "cycle2_peaks": {"tension_mpa": t2, "compression_mpa": k2}, "peaks": peaks, "loops": loops, "warnings": warnings}
