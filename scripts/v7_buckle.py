"""
v7_buckle.py  (designed on the development set only; see V7_PROTOCOL.md amendment 1)
==========================================================

Buckling-distortion layer for the V7 loop-shape layer (L3). Maps the no-buckling material-core loop (v7_loop.unit_loop, Chaboche N = 2 with the
machine-compliance apparent modulus) to the loop of a slender bar, given the tension and compression peak stress (t, k) of the cycle.

Model ("plastic-hinge mechanism envelope")
------------------------------------------
Ideal cycle path 0 -> +eps_a -> -eps_a -> 0, three segments A (0 -> +eps_a), B (+eps_a -> -eps_a), C (-eps_a -> 0).

B, compression half, travel u = eps_a - eps in [0, 2 eps_a]:
    s_c(u)   = (t / t_core) * core(u)                                  core loop scaled multiplicatively to the tension peak
    cap(u)   = B * (1 + max(u - u_z, 0) / u0) ** (-q)                  post-buckling capacity, u_z = travel at which s_c crosses zero
    s(u)     = smax(s_c(u), -cap(u))                                   smooth maximum, width delta = DELTA * k
    B is solved so that min s = -k (the given compression peak). If k >= compression peak of the scaled core the bar does not buckle and the
    core loop (affinely scaled to t, k) is returned unchanged (the v7_loop behaviour).
    With q = 1/2 and u0 -> 0 this is the three-hinge mechanism  sigma = c * fy / ((L/D) sqrt(eps_g))  (Gomes-Appleton 1997 form, as given in
    Kashani 2024 review Eq. 3); u0 > 0 regularises the onset (axial-flexure interaction, residual bow).

A and C, tension half, travel v = eps + eps_a in [0, 2 eps_a], z = v / (2 eps_a), from the reversal stress s_r = s(u = 2 eps_a) to the peak t:
    s(z) = s_r + (t - s_r) * H(z),   H(z) = w * (1 - exp(-z / ze)) / (1 - exp(-1 / ze)) + (1 - w) * z ** m
    first term: unloading of the bent bar and reversal of the hinge (fast, saturating); second term: stiffening while the bar straightens
    (convex, reaches the core level only at full straightening = the tension peak).
    Segment C uses z in [0, 1/2] (optionally towards the NEXT cycle's peak t_next), segment A uses z in [1/2, 1] (optionally starting from the
    PREVIOUS cycle's reversal level, scaled by k_prev / k).

Shape parameters g = (ln(u0/eps_a), q, logit w, ln ze, ln m); each is LINEAR in the features
    phi = (1, ln(lambda_p / LAMBDA_REF), ln(amp / 3 %), x = n / Nf)
so params is a (5, 4) matrix (20 numbers). lambda_p = (L/D) sqrt(fy / 100)  (Dhakal-Maekawa 2002 bar buckling parameter).

Interface
---------
    fit(train_files, core_theta, mode="two_stage" | "joint")  -> params (5 x 4 array)
    predict(file_spec, core_theta, params, peaks)             -> stress at the observed strain samples
        file_spec: dict with d, ld, n (amplitude %), nf, and the loop dict (arc, e) under key "loop" plus its cycle number "cycle"
        peaks: dict t, k and optionally t_next, k_prev
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.optimize import least_squares

PROJ = os.path.dirname(os.path.abspath(__file__))
if PROJ not in sys.path:
    sys.path.insert(0, PROJ)
import v7_data as VD      # noqa: E402
import v7_loop as L       # noqa: E402

FY = VD.FY_NOMINAL
LAMBDA_REF = 10.0 * math.sqrt(FY / 100.0)
DELTA = 0.10
COMP_FORM = "power"      # "power" (hinge-mechanism type) | "linear" (Dhakal-Maekawa-type linear softening to a 0.2 floor)
TENS_FORM = "blend"      # "blend" (straightening function H) | "core" (core reloading branch rescaled between s_r and t)
NPAR = 5
G0 = np.array([0.0, 0.5, 0.0, math.log(0.08), math.log(2.0)])
GLB = np.array([-5.0, 0.05, -6.0, math.log(0.005), math.log(0.3)])
GUB = np.array([5.0, 3.0, 6.0, math.log(2.0), math.log(8.0)])
I0, I1, I2 = L.PTS // 4, 3 * L.PTS // 4, L.PTS          # 50, 150, 200
_CORE_CACHE: dict = {}


def lambda_p(ld: float) -> float:
    return ld * math.sqrt(FY / 100.0)


def features(ld: float, amp_pct: float, x: float) -> np.ndarray:
    return np.array([1.0, math.log(lambda_p(ld) / LAMBDA_REF), math.log(amp_pct / 3.0), x])


def core_loop(core_theta, d, ld, n) -> np.ndarray:
    key = (tuple(np.round(core_theta, 9)), d, ld, n)
    if key not in _CORE_CACHE:
        _CORE_CACHE[key] = L.unit_loop(np.asarray(core_theta), d, ld, n / 100.0)
    return _CORE_CACHE[key]


def _smax(a, b, delta):
    return 0.5 * (a + b + np.sqrt((a - b) ** 2 + delta * delta))


def shape_on_ideal(unit: np.ndarray, eps: float, g: np.ndarray, t: float, k: float, t_next=None, k_prev=None) -> np.ndarray:
    """Distorted loop on the ideal path grid (L.S_IDEAL) for shape parameters g and peaks (t, k)."""
    g = np.clip(g, GLB, GUB)
    tc, kc = float(unit.max()), float(-unit.min())
    a = t / tc
    if k >= a * kc * 0.999:                                  # no buckling: the core loop, affinely scaled
        return L.scale_to_peaks(unit, t, k)
    u0 = math.exp(g[0]) * eps
    q = g[1]
    w = 1.0 / (1.0 + math.exp(-g[2]))
    ze = math.exp(g[3])
    m = math.exp(g[4])
    sc = a * unit[I0:I1 + 1]
    u = np.linspace(0.0, 2.0 * eps, sc.size)
    neg = np.nonzero(sc <= 0.0)[0]
    j = int(neg[0]) if neg.size else sc.size - 1
    uz = u[j - 1] + (u[j] - u[j - 1]) * sc[j - 1] / (sc[j - 1] - sc[j]) if j > 0 and sc[j - 1] > sc[j] else u[j]
    if COMP_FORM == "power":
        decay = (1.0 + np.maximum(u - uz, 0.0) / u0) ** (-q)
    else:
        decay = np.maximum(1.0 - np.maximum(u - uz, 0.0) / (u0 * 4.0), 0.2)
    delta = DELTA * k
    lo, hi = 0.2 * k, 1.5 * a * kc + k
    for _ in range(30):                                      # bisection on the mechanism strength B so that the compression peak is k
        mid = 0.5 * (lo + hi)
        if -_smax(sc, -mid * decay, delta).min() > k:
            hi = mid
        else:
            lo = mid
    sb = _smax(sc, -0.5 * (lo + hi) * decay, delta)
    sb = sb + (t - sb[0]) * np.linspace(1.0, 0.0, sb.size) ** 8      # remove the smoothing offset at the tension peak
    sr = float(sb[-1])
    z = np.linspace(0.0, 1.0, 2 * (I2 - I1) + 1)
    if TENS_FORM == "blend":
        h = w * (1.0 - np.exp(-z / ze)) / (1.0 - math.exp(-1.0 / ze)) + (1.0 - w) * z ** m
    else:
        ct = np.concatenate([unit[I1:I2], unit[:I0 + 1]])
        h = (ct - ct[0]) / (ct[-1] - ct[0])
    nh = I2 - I1
    sr_a = sr * (k_prev / k) if k_prev else sr
    seg_a = sr_a + (t - sr_a) * h[nh:]
    seg_c = sr + ((t_next if t_next else t) - sr) * h[:nh + 1]
    return np.concatenate([seg_a[:-1], sb, seg_c[1:]])


def g_of(params: np.ndarray, ld, amp, x) -> np.ndarray:
    P = np.asarray(params).reshape(NPAR, -1)
    return P @ features(ld, amp, x)[: P.shape[1]]


def predict(file_spec: dict, core_theta, params, peaks: dict) -> np.ndarray:
    f = file_spec
    unit = core_loop(core_theta, f["d"], f["ld"], f["n"])
    g = g_of(params, f["ld"], f["n"], f["cycle"] / f["nf"])
    y = shape_on_ideal(unit, f["n"] / 100.0, g, peaks["t"], peaks["k"], peaks.get("t_next"), peaks.get("k_prev"))
    return np.interp(f["loop"]["arc"], L.S_IDEAL, y)


def core_predict(file_spec: dict, core_theta, peaks: dict) -> np.ndarray:
    f = file_spec
    return np.interp(f["loop"]["arc"], L.S_IDEAL, L.scale_to_peaks(core_loop(core_theta, f["d"], f["ld"], f["n"]), peaks["t"], peaks["k"]))


def _loops(files):
    for f in files:
        for c, lp in f["loops"].items():
            yield f, c, lp


def fit_loop(unit, eps, lp, g_start=G0):
    """Per-loop shape parameters (own peaks only)."""
    def res(g):
        return np.interp(L.S_GRID, L.S_IDEAL, shape_on_ideal(unit, eps, g, lp["t"], lp["k"])) - lp["y"]
    best = None
    for g0 in (g_start, np.array([1.5, 0.5, 2.0, math.log(0.05), 0.0]), np.array([-1.0, 0.7, -1.0, math.log(0.1), math.log(3.0)])):
        r = least_squares(res, g0, bounds=(GLB, GUB), x_scale=np.array([1, 0.3, 1, 0.5, 0.5]))
        if best is None or r.cost < best.cost:
            best = r
    return best.x, math.sqrt(2 * best.cost / L.NFIT)


def per_loop_table(files, core_theta) -> list:
    rows = []
    for f, c, lp in _loops(files):
        unit = core_loop(core_theta, f["d"], f["ld"], f["n"])
        buck = lp["k"] < 0.999 * lp["t"] / unit.max() * (-unit.min())
        g, rm = fit_loop(unit, f["n"] / 100.0, lp) if buck else (np.full(NPAR, np.nan), float("nan"))
        rows.append({"d": f["d"], "ld": f["ld"], "n": f["n"], "c": c, "x": c / f["nf"], "g": g, "rmse": rm, "buck": bool(buck)})
    return rows


def fit_two_stage(rows: list, nfeat: int = 4) -> np.ndarray:
    """Stage 2: robust regression of the per-loop parameters on the features (rows = per_loop_table of the TRAINING files)."""
    rr = [r for r in rows if r["buck"]]
    X = np.array([features(r["ld"], r["n"], r["x"])[:nfeat] for r in rr])
    G = np.array([r["g"] for r in rr])
    P = np.zeros((NPAR, nfeat))
    for j in range(NPAR):
        y = G[:, j]
        wts = np.ones(len(y))
        for _ in range(8):                                   # Huber-type reweighting: per-loop parameters are noisy when weakly identified
            beta = np.linalg.lstsq(X * wts[:, None], y * wts, rcond=None)[0]
            r_ = y - X @ beta
            s_ = 1.4826 * np.median(np.abs(r_)) + 1e-9
            wts = np.sqrt(np.minimum(1.0, 1.5 * s_ / (np.abs(r_) + 1e-12)))
        P[j] = beta
    return P


def fit_joint(files, core_theta, p0: np.ndarray, max_nfev: int = 25) -> np.ndarray:
    """Direct least squares of the (5 x nfeat) matrix on all training loops (starting from the two-stage estimate)."""
    items = [(core_loop(core_theta, f["d"], f["ld"], f["n"]), f["n"] / 100.0, features(f["ld"], f["n"], c / f["nf"])[: p0.shape[1]], lp)
             for f, c, lp in _loops(files)]

    def res(p):
        P = p.reshape(p0.shape)
        return np.concatenate([np.interp(L.S_GRID, L.S_IDEAL, shape_on_ideal(u, e, P @ ph, lp["t"], lp["k"])) - lp["y"] for u, e, ph, lp in items])
    return least_squares(res, p0.ravel(), max_nfev=max_nfev, x_scale="jac").x.reshape(p0.shape)


def fit(train_files, core_theta, mode: str = "joint", nfeat: int = 4, rows=None) -> np.ndarray:
    rows = rows if rows is not None else per_loop_table(train_files, core_theta)
    P = fit_two_stage(rows, nfeat)
    return fit_joint(train_files, core_theta, P) if mode == "joint" else P


def rmse_file_loops(files, core_theta, params, neighbours=None):
    """Shape-only RMSE at the observed samples for every loop: (ld, d, n, c, is_cycle2, rmse_model, rmse_core, rmse_flat)."""
    out = []
    for f, c, lp in _loops(files):
        spec = {"d": f["d"], "ld": f["ld"], "n": f["n"], "nf": f["nf"], "cycle": c, "loop": lp}
        pk = {"t": lp["t"], "k": lp["k"]}
        if neighbours is not None:
            pk.update(neighbours.get((f["d"], f["ld"], f["n"], c), {}))
        pm = predict(spec, core_theta, params, pk)
        pc = core_predict(spec, core_theta, pk)
        rm = lambda p: float(np.sqrt(np.mean((p - lp["s"]) ** 2)))
        out.append((f["ld"], f["d"], f["n"], c, c == 2, rm(pm), rm(pc), rm(np.full_like(lp["s"], lp["s"].mean()))))
    return out
