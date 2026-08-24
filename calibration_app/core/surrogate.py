"""
core/surrogate.py
=================

Model-agnostic pure-Python cyclic-plasticity surrogates.

Every constitutive model exposes the same interface (`CyclicPlasticityModel`)
so the optimiser, objective, session system and GUI never need to know which
model they are driving. Adding a new model (e.g. UVC, Ohno-Wang) means adding a
new subclass here and registering it in `MODEL_REGISTRY` — nothing else changes.

Currently implemented
---------------------
* ``ChabocheModel`` — native ABAQUS combined isotropic + nonlinear-kinematic
  hardening with a **configurable number of backstresses** (2, 3 or 4).
  The discrete update is the same fully-implicit radial return ABAQUS uses at a
  uniaxial material point, so identified parameters transfer to ABAQUS exactly.

Planned (added after ChabocheModel is verified — see model_recommendation.md)
-----------------------------------------------------------------------------
* ``UVCModel``      — Updated Voce-Chaboche (Hartloper & Lignos 2021).
* ``OhnoWangModel`` — Ohno-Wang (1993) with critical-state dynamic recovery.

Constitutive equations (Chaboche, small-strain, rate-independent, uniaxial)
---------------------------------------------------------------------------
    Elastic:    sig_dot = E (eps_dot - eps_p_dot)
    Yield:      f = |sig - X| - (sig_y0 + R) <= 0,   X = sum_k Xk
    Flow:       eps_p_dot = dp * sign(sig - X)
    Kinematic:  Xk_dot = Ck eps_p_dot - gamma_k Xk dp     (Armstrong-Frederick)
    Isotropic:  R = Q_inf (1 - exp(-b p))                 (Voce)

Parameter vector layout (Chaboche, N backstresses):
    (C1, gamma1, C2, gamma2, ..., CN, gammaN, Q_inf, b)   -> length 2N + 2
Fixed (not identified): E, sig_y0  (both are per-run inputs, not hardcoded).
"""

from __future__ import annotations

import math
from abc import ABC, abstractmethod
from typing import Sequence

import numpy as np

# Newton solver tolerances for the scalar plastic-multiplier consistency eqn.
_NEWTON_MAX_ITER: int = 40
_NEWTON_TOL: float = 1.0e-9          # relative residual tolerance
_BISECT_MAX_ITER: int = 80


# ---------------------------------------------------------------------------
# Hot-loop integrator kernel.
#
# The per-strain-increment radial return runs a Newton loop with a per-backstress
# inner sum. Written with plain SCALAR arithmetic (no numpy calls inside the
# loop) because the backstress arrays are tiny (2-4 elements): numpy's per-call
# overhead there dominated and made each evaluation ~0.5 s. numba JIT-compiles
# this to machine code for a large speed-up; if numba is unavailable the same
# function runs in pure Python (still far faster than numpy-per-step).
# ---------------------------------------------------------------------------
def _integrate_impl(eps, C, g, Q, b, E, sy0):
    n = eps.shape[0]
    nb = C.shape[0]
    sig_out = np.empty(n, dtype=np.float64)
    X = np.zeros(nb, dtype=np.float64)
    sig = 0.0
    p = 0.0
    if n == 0:
        return sig_out
    sig_out[0] = 0.0
    eps_prev = eps[0]

    for i in range(1, n):
        d_eps = eps[i] - eps_prev
        eps_prev = eps[i]

        sig_tr = sig + E * d_eps
        Xsum = 0.0
        for k in range(nb):
            Xsum += X[k]
        eta_tr = sig_tr - Xsum
        R = Q * (1.0 - np.exp(-b * p))

        if abs(eta_tr) - (sy0 + R) <= 0.0:
            sig = sig_tr                       # elastic step
        else:
            s = 1.0 if eta_tr >= 0.0 else -1.0
            eta_abs = abs(eta_tr)

            Csum = 0.0
            for k in range(nb):
                Csum += C[k]
            qb = Q * b
            denom0 = E + Csum + (qb if qb > -0.1 * E else -0.1 * E)
            if denom0 < E:
                denom0 = E
            dp = eta_abs / denom0
            if dp < 1e-12:
                dp = 1e-12

            converged = False
            for _ in range(40):
                ex = np.exp(-b * (p + dp))
                F = eta_abs - E * dp - sy0 - Q * (1.0 - ex)
                dF = -E - Q * b * ex
                for k in range(nb):
                    one_p = 1.0 + g[k] * dp
                    A = C[k] - g[k] * s * X[k]
                    F -= dp * A / one_p
                    dF -= A / (one_p * one_p)
                if abs(F) < 1e-9 * (sy0 + eta_abs):
                    converged = True
                    break
                if dF == 0.0:
                    break
                dp_new = dp - F / dF
                if dp_new <= 0.0:
                    dp_new = dp * 0.5
                dp = dp_new

            if not converged:                  # bounded bisection fallback
                lo = 0.0
                hi = eta_abs / E
                if hi < 1.0:
                    hi = 1.0
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    Fmid = eta_abs - E * mid - sy0 - Q * (1.0 - np.exp(-b * (p + mid)))
                    for k in range(nb):
                        one_p = 1.0 + g[k] * mid
                        A = C[k] - g[k] * s * X[k]
                        Fmid -= mid * A / one_p
                    if Fmid > 0.0:
                        lo = mid
                    else:
                        hi = mid
                    if hi - lo < 1e-14:
                        break
                dp = 0.5 * (lo + hi)

            for k in range(nb):
                X[k] = (X[k] + C[k] * s * dp) / (1.0 + g[k] * dp)
            sig = sig_tr - E * s * dp
            p = p + dp

        if not np.isfinite(sig):
            for j in range(i, n):
                sig_out[j] = np.nan
            return sig_out
        sig_out[i] = sig

    return sig_out


try:                                            # JIT-compile if numba is present
    from numba import njit as _njit
    _integrate = _njit(cache=True, fastmath=True)(_integrate_impl)
    _HAVE_NUMBA = True
except Exception:                               # pure-Python fallback
    _integrate = _integrate_impl
    _HAVE_NUMBA = False


class CyclicPlasticityModel(ABC):
    """Abstract interface every constitutive-model surrogate implements.

    Concrete subclasses must be *pure functions of their parameters*: no hidden
    state between :meth:`simulate` calls, so the optimiser can evaluate
    candidates in parallel workers safely.
    """

    #: short machine-friendly identifier, e.g. "chaboche" / "uvc" / "ohno_wang"
    key: str = "abstract"
    #: human-readable name for the GUI
    display_name: str = "Abstract model"

    @property
    @abstractmethod
    def param_names(self) -> tuple[str, ...]:
        """Ordered parameter names, matching the vector `simulate` expects."""

    @property
    @abstractmethod
    def param_units(self) -> tuple[str, ...]:
        """Units for each parameter, same order as :attr:`param_names`."""

    def n_params(self) -> int:
        """Number of free (identified) parameters."""
        return len(self.param_names)

    @abstractmethod
    def simulate(
        self,
        strain_path: Sequence[float],
        params: Sequence[float],
        *,
        E: float,
        sy0: float,
    ) -> np.ndarray:
        """Integrate the model along a prescribed total-strain path.

        Parameters
        ----------
        strain_path : 1-D array-like
            Total strain at each time station (dimensionless, NOT percent).
        params : sequence of float
            Model parameters in :attr:`param_names` order.
        E, sy0 : float
            Young's modulus (MPa) and initial yield stress (MPa) — per-run
            inputs, never hardcoded to a steel grade.

        Returns
        -------
        numpy.ndarray of stress (MPa), same length as ``strain_path``. Returns
        an array containing NaNs if a non-finite value is produced, so the
        caller can treat it as a failed evaluation.
        """

    @abstractmethod
    def material_block(
        self,
        params: Sequence[float],
        *,
        E: float,
        nu: float,
        sy0: float,
        name: str = "STEEL",
    ) -> str:
        """Return the ABAQUS material keyword block for these parameters.

        For models with no native ABAQUS support (UVC, Ohno-Wang) this returns
        the ``*MATERIAL`` + ``*DEPVAR`` + ``*USER MATERIAL`` block referencing
        the model's UMAT.
        """

    @abstractmethod
    def default_bounds(self) -> dict[str, tuple[float, float]]:
        """Literature-derived default (lower, upper) bounds per parameter.

        These are *defaults* only; the GUI lets the user widen them (needed for
        batch-variable scrap rebar). See ``utils.physics_gate``.
        """


class ChabocheModel(CyclicPlasticityModel):
    """Combined isotropic/kinematic Chaboche with N Armstrong-Frederick terms.

    Parameters
    ----------
    n_backstresses : int
        Number of kinematic backstresses (2, 3 or 4). Default 3 — the
        minimum-sufficient count for wide-range (1-6%) cyclic strain per the
        literature synthesis (see research_report.md, Q1B).
    """

    key = "chaboche"

    def __init__(self, n_backstresses: int = 3) -> None:
        if n_backstresses not in (2, 3, 4):
            raise ValueError("n_backstresses must be 2, 3 or 4")
        self.n_backstresses = int(n_backstresses)
        self.display_name = f"Chaboche ({self.n_backstresses}-backstress)"

    # -- interface -----------------------------------------------------------
    @property
    def param_names(self) -> tuple[str, ...]:
        names: list[str] = []
        for k in range(1, self.n_backstresses + 1):
            names += [f"C{k}", f"gamma{k}"]
        names += ["Q_inf", "b"]
        return tuple(names)

    @property
    def param_units(self) -> tuple[str, ...]:
        units: list[str] = []
        for _ in range(self.n_backstresses):
            units += ["MPa", "-"]
        units += ["MPa", "-"]
        return tuple(units)

    def _split(self, params: Sequence[float]) -> tuple[np.ndarray, np.ndarray, float, float]:
        """Split a flat parameter vector into (C[], gamma[], Q_inf, b)."""
        p = np.asarray(params, dtype=np.float64)
        n = self.n_backstresses
        if p.size != 2 * n + 2:
            raise ValueError(
                f"expected {2 * n + 2} parameters for {n} backstresses, got {p.size}"
            )
        C = p[0 : 2 * n : 2].copy()
        g = p[1 : 2 * n : 2].copy()
        Q = float(p[2 * n])
        b = float(p[2 * n + 1])
        return C, g, Q, b

    # -- integrator ----------------------------------------------------------
    def simulate(
        self,
        strain_path: Sequence[float],
        params: Sequence[float],
        *,
        E: float,
        sy0: float,
    ) -> np.ndarray:
        C, g, Q, b = self._split(params)
        eps = np.ascontiguousarray(strain_path, dtype=np.float64)
        return _integrate(eps, np.ascontiguousarray(C), np.ascontiguousarray(g),
                          float(Q), float(b), float(E), float(sy0))

    # -- ABAQUS export -------------------------------------------------------
    def material_block(
        self,
        params: Sequence[float],
        *,
        E: float,
        nu: float,
        sy0: float,
        name: str = "STEEL",
    ) -> str:
        """Native ABAQUS combined-Chaboche card.

        ``*PLASTIC`` parameter line is ``sy0, C1, g1, ..., CN, gN`` (1 + 2N
        values). ABAQUS allows at most 8 values per data line, so N=4 (9 values)
        is wrapped onto a continuation line.
        """
        C, g, Q, b = self._split(params)
        vals = [sy0]
        for Ck, gk in zip(C, g):
            vals += [Ck, gk]

        lines = [
            f"*MATERIAL, NAME={name}",
            "*ELASTIC",
            f"{E:.6g}, {nu:.6g}",
            f"*PLASTIC, HARDENING=COMBINED, DATATYPE=PARAMETERS, "
            f"NUMBER BACKSTRESSES={self.n_backstresses}",
        ]
        for i in range(0, len(vals), 8):                 # 8 values per data line
            chunk = vals[i : i + 8]
            lines.append(", ".join(f"{v:.6f}" for v in chunk))
        lines += [
            "*CYCLIC HARDENING, PARAMETERS",
            f"{sy0:.6f}, {Q:.6f}, {b:.6f}",
        ]
        return "\n".join(lines) + "\n"

    # -- default bounds ------------------------------------------------------
    def default_bounds(self) -> dict[str, tuple[float, float]]:
        """Literature envelope (research_report.md §5 / model_recommendation.md
        §3). Structural-steel-class Chaboche calibrations; generous by design
        and meant to be widened by the user for batch-variable scrap rebar."""
        b: dict[str, tuple[float, float]] = {}
        # Per-backstress bounds. Term 1 = fast (near-yield), later terms slower.
        per_term = [
            ((2_000.0, 30_000.0), (15.0, 800.0)),   # C1, gamma1  (Krolo 2016; Hartloper 2021)
            ((500.0, 10_000.0), (2.0, 350.0)),      # C2, gamma2  (Q235; S275; UVC)
            ((200.0, 5_000.0), (1.0, 100.0)),       # C3, gamma3  (S355; SA508)
            ((100.0, 3_000.0), (0.5, 60.0)),        # C4, gamma4  (near-linear tail)
        ]
        for k in range(self.n_backstresses):
            (cb, gb) = per_term[k]
            b[f"C{k+1}"] = cb
            b[f"gamma{k+1}"] = gb
        b["Q_inf"] = (-250.0, 230.0)   # negative allowed (cyclic softening)
        # Isotropic saturation rate. Literature for structural/rebar-class steel
        # clusters at b ~ 3-14 (S355 3.2, S275 4.4, UVC 0.11-14.1); Q235's b=40
        # is an outlier, so the default ceiling is 15. Lower the upper bound in
        # the UI (e.g. 8) for a more classic low-b fit at a small RMSE cost.
        b["b"] = (1.0, 15.0)
        return b


# ---------------------------------------------------------------------------
# Registry — the GUI/optimiser discover available models through this dict.
# UVC and Ohno-Wang are registered here once their subclasses land.
# ---------------------------------------------------------------------------
def make_chaboche(n_backstresses: int = 3) -> ChabocheModel:
    """Factory used by the registry so the GUI can pick the backstress count."""
    return ChabocheModel(n_backstresses=n_backstresses)


MODEL_REGISTRY: dict[str, str] = {
    "chaboche": "Chaboche (combined, N backstresses)",
    "uvc": "Updated Voce-Chaboche (Hartloper & Lignos 2021)",
    "ohno_wang": "Ohno-Wang (1993)",
}


def self_test() -> tuple[bool, str]:
    """+/-1% triangular strain, 3 cycles, default 3-backstress Chaboche.

    Checks the integrator runs and produces finite, physically-plausible peak
    stresses for a mid-range parameter set (sanity only, not a calibration).
    """
    model = ChabocheModel(n_backstresses=3)
    # (C1,g1, C2,g2, C3,g3, Q, b)
    params = (12_000.0, 300.0, 2_000.0, 60.0, 800.0, 12.0, -40.0, 8.0)
    npts = 601
    t = np.linspace(0.0, 3.0, npts)
    tri = 2.0 * np.abs(2.0 * (t - np.floor(t + 0.5))) - 1.0
    eps = 0.01 * tri
    sig = model.simulate(eps, params, E=200_000.0, sy0=500.0)
    if not np.all(np.isfinite(sig)):
        return False, "non-finite stress produced"
    peak = float(np.max(np.abs(sig)))
    ok = 400.0 <= peak <= 900.0
    return ok, f"peak |sigma| = {peak:.1f} MPa ({'PASS' if ok else 'FAIL'})"


if __name__ == "__main__":
    print("surrogate self-test:", self_test()[1])
