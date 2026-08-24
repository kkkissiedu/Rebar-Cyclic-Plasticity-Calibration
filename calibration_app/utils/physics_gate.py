"""
utils/physics_gate.py
======================

Grade-agnostic physical-validity gate for cyclic-plasticity parameters.

Design principle (see project_context memory)
---------------------------------------------
The target material is batch-variable scrap-metal rebar with **unknown** grade,
so the gate must NOT enforce B500C-specific literature values. Instead it checks
each candidate for *physical consistency against the loaded experimental data*:

  * every kinematic saturation rate gamma_k > 0 (Armstrong-Frederick requires it);
  * isotropic rate ``b`` within a sane numeric window (avoids exp overflow /
    degenerate saturation);
  * total kinematic saturation ``sum_k C_k/gamma_k`` (max backstress the model
    can develop above yield) does not exceed a multiple of the stress-above-
    yield actually observed in the data;
  * the saturated yield ``sigma_y0 + Q_inf`` does not collapse below a fraction
    of ``sigma_y0`` (real steel does not lose most of its strength in a handful
    of cycles).

The literature-derived *default search bounds* live on each model
(`model.default_bounds()`); those are the only place grade-informed numbers
appear, and the GUI lets the user widen them. This gate is the hard floor that
even widened bounds cannot bypass, but it is expressed relative to the data, not
to a grade.

Every numeric limit below is a named constant with its rationale, so there are
no magic numbers.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Sequence

# --- grade-agnostic limit constants (ratios/windows, not grade values) ------

#: gamma must be strictly positive for Armstrong-Frederick recall to make sense.
MIN_GAMMA: float = 1.0e-6

#: Isotropic saturation-rate window. b<~0.1 is effectively no saturation over a
#: realistic accumulated-plastic-strain range; b>~50 saturates within the first
#: reversal and risks exp() underflow. Kept wide; the model's default bounds are
#: tighter (2-40) but the user may widen them for scrap rebar.
MIN_B: float = 0.1
MAX_B: float = 60.0

#: The combined saturated backstress (sum C_k/gamma_k) is the maximum kinematic
#: hardening the model can add above initial yield. It should not exceed this
#: multiple of the stress-above-yield actually seen in the data. 1.6x leaves
#: headroom for transient over/undershoot without allowing runaway hardening.
MAX_SATURATION_FACTOR: float = 1.6

#: Saturated yield (sigma_y0 + Q_inf) must stay at least this fraction of the
#: initial yield. B500C softens, but not to nothing; scrap rebar may soften
#: more, so this is a floor, not a grade value. User-adjustable via the ctor.
DEFAULT_MIN_YIELD_FRACTION: float = 0.35

#: Absolute backstop the GUI can never go below, even if the user lowers
#: min_yield_fraction - guards against a degenerate "no yield left" fit.
ABSOLUTE_MIN_YIELD_FRACTION: float = 0.15

#: Minimum ratio gamma_k / gamma_{k+1} between successive backstresses.
#: Chaboche's multi-kinematic decomposition (Chaboche 1986, Int. J.
#: Plasticity 2, Eq. 5: sum of N Armstrong-Frederick terms) is only
#: meaningful when each term covers a DISTINCT strain-range scale - term 1
#: fast/near-yield, later terms progressively slower. Bari & Hassan 2000
#: (Int. J. Plasticity 16, Sec. 4) calibrate with roughly decade-separated
#: gammas. Without a floor on the ratio the optimiser can collapse all
#: gammas to one value (observed: three gammas = 75.16 on 2data), which is
#: an over-parameterised single backstress. Ratio, not absolute - grade-
#: agnostic by construction.
MIN_GAMMA_RATIO: float = 3.0


_C_RE = re.compile(r"^C(\d+)$")
_G_RE = re.compile(r"^gamma(\d+)$")


@dataclass
class DataContext:
    """Facts extracted from the loaded experiment that the gate checks against.

    Attributes
    ----------
    sigma_y0 : float
        Initial yield stress for this run (per-run input; from a companion
        monotonic test or entered manually).
    measured_peak_stress : float
        Maximum |stress| observed across the calibration cycles (MPa).
    """

    sigma_y0: float
    measured_peak_stress: float

    @property
    def stress_above_yield(self) -> float:
        """Observed hardening headroom above initial yield (>= small positive)."""
        return max(self.measured_peak_stress - self.sigma_y0, 1.0)

    @classmethod
    def from_datasets(cls, datasets, sigma_y0: float) -> "DataContext":
        """Build a context from one or more loaded specimens (Option C).

        The saturation cap must be judged against the full range of behaviour
        the model has to reproduce, so when several strain amplitudes are loaded
        together the **cross-amplitude peak stress** (the max over all specimens,
        typically the highest-amplitude test) sets the reference - not a single
        low-amplitude file's small headroom. This keeps the gate consistent with
        the grade-agnostic, data-relative design.

        Parameters
        ----------
        datasets : one data dict or an iterable of them (from
            ``data_loader.load_experimental``); each must expose
            ``measured_peak_stress``.
        sigma_y0 : per-run yield stress.
        """
        if isinstance(datasets, dict):
            datasets = [datasets]
        peaks = [float(d["measured_peak_stress"]) for d in datasets
                 if d.get("measured_peak_stress")]
        peak = max(peaks) if peaks else 0.0
        return cls(sigma_y0=sigma_y0, measured_peak_stress=peak)


@dataclass
class PhysicsGate:
    """Grade-agnostic validity gate, parameterised by the loaded data.

    Parameters
    ----------
    data : DataContext
        The per-run data facts to check candidates against.
    min_yield_fraction : float
        Minimum ``(sigma_y0 + Q_inf) / sigma_y0``. Clamped to at least
        :data:`ABSOLUTE_MIN_YIELD_FRACTION`.
    max_saturation_factor : float
        Cap on ``sum C_k/gamma_k`` as a multiple of observed stress-above-yield.
    """

    data: DataContext
    min_yield_fraction: float = DEFAULT_MIN_YIELD_FRACTION
    max_saturation_factor: float = MAX_SATURATION_FACTOR
    min_b: float = MIN_B
    max_b: float = MAX_B
    _yield_fraction: float = field(init=False)

    def __post_init__(self) -> None:
        self._yield_fraction = max(self.min_yield_fraction,
                                   ABSOLUTE_MIN_YIELD_FRACTION)

    # -- name-convention parsing (model-agnostic for backstress models) ------
    @staticmethod
    def _kinematic_saturation(param_names: Sequence[str],
                              params: Sequence[float]) -> float | None:
        """Sum_k C_k / gamma_k from a (names, values) pair, or None if the
        model doesn't follow the C{k}/gamma{k} naming convention."""
        C: dict[int, float] = {}
        G: dict[int, float] = {}
        for name, val in zip(param_names, params):
            mc = _C_RE.match(name)
            mg = _G_RE.match(name)
            if mc:
                C[int(mc.group(1))] = float(val)
            elif mg:
                G[int(mg.group(1))] = float(val)
        if not C or set(C) != set(G):
            return None
        total = 0.0
        for k in C:
            gk = G[k]
            if gk <= MIN_GAMMA:
                return float("inf")   # forces rejection below
            total += C[k] / gk
        return total

    @staticmethod
    def _get(param_names: Sequence[str], params: Sequence[float],
             key: str) -> float | None:
        for name, val in zip(param_names, params):
            if name == key:
                return float(val)
        return None

    def _headroom(self, sy0: float) -> float:
        """Observed stress-above-yield for a given (possibly calibrated) sy0."""
        return max(self.data.measured_peak_stress - sy0, 1.0)

    # -- projection (repair to the feasible set) -----------------------------
    def project(self, model, params: Sequence[float],
                sy0: float = None) -> tuple:
        """Repair an infeasible candidate onto the feasible boundary.

        This is the mechanism the optimiser uses instead of hard rejection: a
        rejected candidate returning a constant PENALTY gives the optimiser a
        flat landscape (it locks onto the first sampled point). Projecting keeps
        *every* evaluation meaningful - the search explores the whole box while
        only ever simulating physically consistent parameters.

        Repairs, in order (each is idempotent, so projecting a feasible point
        returns it unchanged):
          1. gamma_k -> strictly positive;
          1b. gamma_{k+1} capped at gamma_k / MIN_GAMMA_RATIO (distinct
              timescales; prevents backstress collapse);
          2. b -> clipped into [min_b, max_b];
          3. Q_inf -> raised so sigma_y0 + Q_inf >= yield-fraction * sigma_y0;
          4. all C_k scaled by (cap / saturation) if the kinematic saturation
             sum_k C_k/gamma_k exceeds the data-relative cap (preserves the
             *shape* / relative weighting of the backstresses).
        """
        if sy0 is None:
            sy0 = self.data.sigma_y0
        names = list(model.param_names)
        vals = [float(v) for v in params]
        pos = {n: i for i, n in enumerate(names)}

        for n, i in pos.items():                       # 1. gamma positivity
            if _G_RE.match(n) and vals[i] <= MIN_GAMMA:
                vals[i] = MIN_GAMMA

        # 1b. gamma timescale separation: gamma_{k+1} <= gamma_k / ratio so
        # each backstress keeps a distinct strain-range role (see
        # MIN_GAMMA_RATIO rationale). Walking k upward is idempotent and
        # leaves already-separated candidates untouched.
        g_idx = sorted(
            ((int(_G_RE.match(n).group(1)), i) for n, i in pos.items()
             if _G_RE.match(n)))
        for (_, i_prev), (_, i_next) in zip(g_idx, g_idx[1:]):
            cap_g = vals[i_prev] / MIN_GAMMA_RATIO
            if vals[i_next] > cap_g:
                vals[i_next] = max(cap_g, MIN_GAMMA)

        if "b" in pos:                                 # 2. isotropic-rate window
            vals[pos["b"]] = min(max(vals[pos["b"]], self.min_b), self.max_b)

        if "Q_inf" in pos:                             # 3. saturated-yield floor
            q_min = (self._yield_fraction - 1.0) * sy0
            if vals[pos["Q_inf"]] < q_min:
                vals[pos["Q_inf"]] = q_min

        sat = self._kinematic_saturation(names, vals)  # 4. saturation cap
        if sat is not None and sat > 0.0:
            cap = self.max_saturation_factor * self._headroom(sy0)
            if sat > cap:
                scale = cap / sat
                for n, i in pos.items():
                    if _C_RE.match(n):
                        vals[i] *= scale
        return tuple(vals)

    # -- the gate ------------------------------------------------------------
    def is_valid(self, model, params: Sequence[float]) -> bool:
        """True if ``params`` are physically consistent with the loaded data."""
        return self.reason(model, params) is None

    def reason(self, model, params: Sequence[float]) -> str | None:
        """Return the first failure reason, or None if the candidate is valid.

        ``model`` is any :class:`~calibration_app.core.surrogate.CyclicPlasticityModel`;
        only its ``param_names`` are used, so the gate stays model-agnostic.
        """
        names = tuple(model.param_names)

        # 1. gamma_k > 0
        for name, val in zip(names, params):
            if _G_RE.match(name) and float(val) <= MIN_GAMMA:
                return f"{name} <= 0 (Armstrong-Frederick requires gamma > 0)"

        # 1b. gamma timescale separation (see MIN_GAMMA_RATIO)
        gammas = sorted(
            ((int(_G_RE.match(n).group(1)), float(v))
             for n, v in zip(names, params) if _G_RE.match(n)))
        for (k_prev, g_prev), (k_next, g_next) in zip(gammas, gammas[1:]):
            if g_next > g_prev / MIN_GAMMA_RATIO:
                return (f"gamma{k_next}={g_next:.3g} not separated from "
                        f"gamma{k_prev}={g_prev:.3g} (need ratio >= "
                        f"{MIN_GAMMA_RATIO:g}; backstresses must cover "
                        f"distinct strain scales)")

        # 2. isotropic rate window
        b = self._get(names, params, "b")
        if b is not None and not (self.min_b <= b <= self.max_b):
            return f"b={b:.3g} outside [{self.min_b}, {self.max_b}]"

        # 3. kinematic saturation vs observed hardening headroom
        sat = self._kinematic_saturation(names, params)
        if sat is not None:
            cap = self.max_saturation_factor * self.data.stress_above_yield
            if sat > cap:
                return (f"kinematic saturation {sat:.0f} MPa exceeds "
                        f"{self.max_saturation_factor}x observed stress-above-"
                        f"yield ({cap:.0f} MPa)")

        # 4. saturated yield floor (relative to sigma_y0, not to a grade)
        Q = self._get(names, params, "Q_inf")
        if Q is not None:
            sy0 = self.data.sigma_y0
            if (sy0 + Q) < self._yield_fraction * sy0:
                return (f"saturated yield {sy0 + Q:.0f} MPa below "
                        f"{self._yield_fraction:.0%} of sigma_y0 ({sy0:.0f})")

        return None


def validate_startup(sigma_y0: float, measured_peak_stress: float) -> list[str]:
    """Sanity-check the run configuration against the loaded data at startup.

    Returns a list of human-readable warnings (empty if all consistent). This is
    grade-agnostic: it only flags internal inconsistencies, e.g. a yield stress
    that exceeds the largest stress ever measured (physically impossible for a
    hardening material under these tests).
    """
    warnings: list[str] = []
    if sigma_y0 > measured_peak_stress:
        warnings.append(
            f"sigma_y0 ({sigma_y0:.0f} MPa) exceeds the measured peak stress "
            f"({measured_peak_stress:.0f} MPa) - check the yield input or the "
            f"stress units (load kN / area).")
    if sigma_y0 <= 0:
        warnings.append(f"sigma_y0 ({sigma_y0:.0f}) must be positive.")
    return warnings
