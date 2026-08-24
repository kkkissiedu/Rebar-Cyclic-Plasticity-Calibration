"""
core/ohno_wang_model.py
=======================

Ohno-Wang model I (1993) - N. Ohno & J.-D. Wang, "Kinematic hardening rules with
critical state of dynamic recovery, Part I: formulation and basic features for
ratchetting behavior", Int. J. Plasticity 9(3):375-390,
doi:10.1016/0749-6419(93)90042-O.

Constitutive equations (uniaxial)
---------------------------------
Backstress decomposition (Ohno & Wang 1993, Eq. 5):
    alpha = sum_{i=1}^{M} alpha_i

Ohno-Wang I kinematic evolution - critical state of dynamic recovery, nonlinear
form (Ohno & Wang 1993, Eq. 10; the bilinear Eq. 8 is the limit m_i -> inf):
    dalpha_i = C_i*n*dp
               - C_i*(|alpha_i|/r_i)^{m_i} * <n * sign(alpha_i)> * sign(alpha_i)*dp
  where n = sign(sigma - alpha), dp = |d eps_p|, C_i = zeta_i*r_i is the initial
  modulus, r_i is the critical (saturation) value of component i, m_i the
  recovery exponent, and <.> is the Macaulay bracket: the dynamic-recovery term
  is active only while loading toward the critical state (n = sign(alpha_i)).
  m_i = 1 (without the bracket) recovers Armstrong-Frederick/Chaboche;
  m_i -> inf recovers the multilinear (bilinear-per-component) rule.

Isotropic hardening (optional Voce extension; Abdel-Karim 2010,
doi:10.1016/j.ijpvp.2010.02.003 - needed because pure OW-I 1993 is kinematic
only and cannot represent the cyclic softening seen in rebar LCF):
    sigma_y(p) = sigma_y0 + Q*(1 - exp(-b*p))
  Setting **Q = 0 recovers the pure Ohno-Wang I 1993 model** (fixed yield size).

Yield function / flow (Ohno & Wang 1993, Eq. 1-2):
    f = |sigma - alpha| - sigma_y <= 0 ,  d eps_p = dp * sign(sigma - alpha).

Parameter vector (identified; sigma_y0 and E separate/per-run):
    (Q, b, C1, r1, m1, ..., CN, rN, mN)   -> length 2 + 3N

Integration
-----------
Radial return with an explicit backstress update: within each strain increment
the flow direction n is fixed, the plastic multiplier dp is solved from the
isotropic-hardening consistency (Newton), then each alpha_i is advanced one
step of Eq. 10. Each input increment is internally SUBSTEPPED to a maximum of
``_OW_SUBSTEP_EPS`` (1e-4) strain (2026-07-14): the explicit Eq. 10 update
has O(dp) error, so without a fixed substep cap the integrated response
depended on the caller's strain sampling (up to ~12 MPa drift between the
254-pt/cycle experimental path and coarser paths at 3% amplitude). With the
cap, the surrogate and the C++ UMAT (umats/OhnoWang.cpp, same cap) both
converge to the same step-size-independent response. OW sessions calibrated
before 2026-07-14 used the un-substepped integrator - re-run before reporting.

FE / UMAT
---------
``umat_file`` points to ``umats/OhnoWang.cpp`` - a C++ UMAT written for this
project (2026-07-14, no permissive Fortran OW UMAT was available and this
machine has no Intel Fortran; ABAQUS compiles it with MSVC). It mirrors this
integrator exactly (same substep cap); cross-validated to < 0.01 MPa against
this surrogate on identical paths.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .surrogate import CyclicPlasticityModel

_OW_NEWTON_ITER: int = 30
_OW_TOL: float = 1.0e-9
#: Max strain per internal substep (2026-07-14). The explicit Eq. 10
#: backstress update is O(dp)-accurate, so the response must not depend on
#: the caller's sampling; this cap makes surrogate and UMAT (OhnoWang.cpp,
#: identical constant) converge to the same path-independent solution.
_OW_SUBSTEP_EPS: float = 1.0e-4
_OW_MAX_SUBSTEP: int = 1000


def _ow_integrate_impl(eps, C, r, m, Q, b, E, sy0):
    """Uniaxial Ohno-Wang I integrator (explicit backstress, substepped;
    numba-friendly). Mirrors umats/OhnoWang.cpp exactly."""
    n = eps.shape[0]
    nb = C.shape[0]
    sig_out = np.empty(n, dtype=np.float64)
    alpha_k = np.zeros(nb, dtype=np.float64)
    sigma = 0.0
    p = 0.0
    if n == 0:
        return sig_out
    sig_out[0] = 0.0
    eps_prev = eps[0]

    for i in range(1, n):
        d_eps = eps[i] - eps_prev
        eps_prev = eps[i]

        nsub = int(np.ceil(abs(d_eps) / _OW_SUBSTEP_EPS))
        if nsub < 1:
            nsub = 1
        if nsub > _OW_MAX_SUBSTEP:
            nsub = _OW_MAX_SUBSTEP
        de = d_eps / nsub

        for _sub in range(nsub):
            sig_tr = sigma + E * de
            alpha = 0.0
            for k in range(nb):
                alpha += alpha_k[k]
            eta = sig_tr - alpha
            sy = sy0 + Q * (1.0 - np.exp(-b * p))

            if abs(eta) - sy <= 0.0:
                sigma = sig_tr                         # elastic
            else:
                ndir = 1.0 if eta >= 0.0 else -1.0
                eta_abs = abs(eta)
                # Newton for dp with backstress frozen; isotropic consistency:
                # F(dp) = eta_abs - E*dp - sy0 - Q*(1 - exp(-b*(p+dp)))
                dp = (eta_abs - sy) / E
                if dp < 1e-12:
                    dp = 1e-12
                for _ in range(_OW_NEWTON_ITER):
                    ex = np.exp(-b * (p + dp))
                    F = eta_abs - E * dp - sy0 - Q * (1.0 - ex)
                    if abs(F) < _OW_TOL * (sy0 + eta_abs):
                        break
                    dF = -E - Q * b * ex
                    if dF == 0.0:
                        break
                    dp_new = dp - F / dF
                    if dp_new <= 0.0:
                        dp_new = dp * 0.5
                    dp = dp_new

                p = p + dp
                sigma = sig_tr - E * ndir * dp
                # explicit Ohno-Wang I backstress update (Eq. 10)
                for k in range(nb):
                    ak = alpha_k[k]
                    if ak > 0.0:
                        sk = 1.0
                    elif ak < 0.0:
                        sk = -1.0
                    else:
                        sk = 0.0
                    if r[k] > 0.0:
                        recov = (abs(ak) / r[k]) ** m[k]
                    else:
                        recov = 0.0
                    active = 1.0 if (ndir * sk) > 0.0 else 0.0
                    alpha_k[k] = ak + C[k] * ndir * dp - C[k] * recov * active * sk * dp

        if not np.isfinite(sigma):
            for j in range(i, n):
                sig_out[j] = np.nan
            return sig_out
        sig_out[i] = sigma

    return sig_out


try:
    from numba import njit as _njit
    _ow_integrate = _njit(cache=True, fastmath=True)(_ow_integrate_impl)
except Exception:
    _ow_integrate = _ow_integrate_impl


class OhnoWangModel(CyclicPlasticityModel):
    """Ohno-Wang I (1993) with an optional Voce isotropic extension."""

    key = "ohno_wang"

    def __init__(self, n_backstresses: int = 3) -> None:
        if n_backstresses not in (2, 3, 4):
            raise ValueError("n_backstresses must be 2, 3 or 4")
        self.n_backstresses = int(n_backstresses)
        self.display_name = f"Ohno-Wang ({self.n_backstresses}-backstress)"
        # C++ UMAT (MSVC-compiled; no Intel Fortran on this machine)
        self.umat_file = "calibration_app/umats/OhnoWang.cpp"
        #: 3D (C3D8 coupon) UMAT used by "Transfer to CAE" (2026-07-14).
        #: PROPS = E, nu, sy0, Q, b, C_k, r_k, m_k; STATEV = 7 + 6N.
        self.umat_multiaxial_file = \
            "calibration_app/umats/OhnoWangMultiaxial.cpp"

    @property
    def param_names(self) -> tuple[str, ...]:
        names = ["Q_inf", "b"]
        for k in range(1, self.n_backstresses + 1):
            names += [f"C{k}", f"r{k}", f"m{k}"]
        return tuple(names)

    @property
    def param_units(self) -> tuple[str, ...]:
        units = ["MPa", "-"]
        for _ in range(self.n_backstresses):
            units += ["MPa", "MPa", "-"]
        return tuple(units)

    def _split(self, params: Sequence[float]):
        p = np.asarray(params, dtype=np.float64)
        n = self.n_backstresses
        if p.size != 2 + 3 * n:
            raise ValueError(f"expected {2 + 3*n} Ohno-Wang params, got {p.size}")
        Q, b = float(p[0]), float(p[1])
        C = p[2::3][:n].copy()
        r = p[3::3][:n].copy()
        m = p[4::3][:n].copy()
        return Q, b, C, r, m

    def simulate(self, strain_path, params, *, E, sy0):
        Q, b, C, r, m = self._split(params)
        eps = np.ascontiguousarray(strain_path, dtype=np.float64)
        return _ow_integrate(eps, np.ascontiguousarray(C), np.ascontiguousarray(r),
                             np.ascontiguousarray(m), float(Q), float(b),
                             float(E), float(sy0))

    def material_block(self, params, *, E, nu, sy0, name="STEEL"):
        """UMAT-based card. PROPS order: E, sy0, Q, b, C1, r1, m1, ..., CN, rN, mN.
        State vars = 1 (eq. plastic strain) + N (backstresses)."""
        Q, b, C, r, m = self._split(params)
        props = [E, sy0, Q, b]
        for Ck, rk, mk in zip(C, r, m):
            props += [Ck, rk, mk]
        nstatv = 1 + self.n_backstresses
        lines = [
            f"*MATERIAL, NAME={name}",
            f"*DEPVAR\n{nstatv}",
            f"*USER MATERIAL, CONSTANTS={len(props)}",
        ]
        for i in range(0, len(props), 8):
            lines.append(", ".join(f"{v:.6f}" for v in props[i:i + 8]))
        return "\n".join(lines) + "\n"

    def default_bounds(self):
        """Ohno-Wang default search bounds. r_i is the per-component saturation
        backstress; sum r_i approximates the stress above yield at saturation.
        m_i: 1 ~ Armstrong-Frederick, large -> multilinear (Kobayashi & Ohno
        2002 use m1=1, m2=2, m3=5 for 316L). Q negative allowed (cyclic
        softening); Q=0 gives the pure OW-I 1993 model."""
        b: dict = {}
        b["Q_inf"] = (-250.0, 230.0)
        b["b"] = (0.1, 15.0)
        per_term = [
            ((1_000.0, 60_000.0), (10.0, 200.0), (0.5, 20.0)),  # C1, r1, m1 (fast)
            ((500.0, 20_000.0), (5.0, 150.0), (0.5, 20.0)),
            ((200.0, 8_000.0), (2.0, 100.0), (0.5, 20.0)),
            ((100.0, 4_000.0), (1.0, 60.0), (0.5, 20.0)),
        ]
        for k in range(self.n_backstresses):
            cb, rb, mb = per_term[k]
            b[f"C{k+1}"] = cb
            b[f"r{k+1}"] = rb
            b[f"m{k+1}"] = mb
        return b


def self_test() -> tuple[bool, str]:
    """+/-1% triangular strain, 3-backstress Ohno-Wang; check finite,
    physically-plausible peak stress."""
    model = OhnoWangModel(n_backstresses=3)
    # (Q, b, C1,r1,m1, C2,r2,m2, C3,r3,m3)
    params = (-30.0, 6.0,
              40000.0, 40.0, 1.0,
              8000.0, 30.0, 2.0,
              2000.0, 20.0, 5.0)
    npts = 601
    t = np.linspace(0.0, 3.0, npts)
    tri = 2.0 * np.abs(2.0 * (t - np.floor(t + 0.5))) - 1.0
    eps = 0.01 * tri
    sig = model.simulate(eps, params, E=200_000.0, sy0=350.0)
    if not np.all(np.isfinite(sig)):
        return False, "non-finite stress produced"
    peak = float(np.max(np.abs(sig)))
    ok = 300.0 <= peak <= 900.0
    return ok, f"peak |sigma| = {peak:.1f} MPa ({'PASS' if ok else 'FAIL'})"


if __name__ == "__main__":
    print("Ohno-Wang self-test:", self_test()[1])
