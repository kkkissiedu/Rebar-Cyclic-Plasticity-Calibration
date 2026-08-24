"""
core/objective.py
=================

The single objective function shared by **every** backend and optimiser.

Why this matters
----------------
The study compares four conditions - Surrogate+DE, Surrogate+Bayesian, FE+DE,
FE+Bayesian. For that comparison to be clean, the *scoring* of a stress-strain
response must be byte-for-byte identical regardless of where the response came
from. So the low-level scorer :func:`score_hysteresis` takes only arrays
(experimental vs simulated stress on a common strain grid) and knows nothing
about surrogates or ABAQUS. Both backends resample onto the experimental strain
points and call it.

Objective
---------
For each cycle in the calibration window::

    RMSE_c = sqrt( weighted_mean( (sigma_sim - sigma_exp)^2 ) )
    objective = mean_c( RMSE_c )                     # MPa

Samples in the first/last ``REVERSAL_FRACTION`` of each cycle's strain
arc-length (the Bauschinger knees) are optionally up-weighted.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np

#: Objective value returned for a failed / physically-invalid evaluation.
PENALTY: float = 1.0e9

#: Fraction of each cycle's strain arc-length near the reversals to up-weight.
REVERSAL_FRACTION: float = 0.15
#: Weight applied to reversal-region samples when weighting is enabled.
REVERSAL_WEIGHT: float = 1.5


def _cycle_rmse(
    strain_c: np.ndarray,
    sim_c: np.ndarray,
    exp_c: np.ndarray,
    weight_reversals: bool,
) -> float:
    """Weighted RMSE (MPa) over one cycle's samples."""
    err = sim_c - exp_c
    if not weight_reversals:
        return math.sqrt(float(np.mean(err * err)))
    ds = np.abs(np.diff(strain_c, prepend=strain_c[0]))
    arc = np.cumsum(ds)
    total = arc[-1] if arc[-1] > 0 else 1.0
    frac = arc / total
    w = np.ones_like(err)
    w[(frac <= REVERSAL_FRACTION) | (frac >= 1.0 - REVERSAL_FRACTION)] = REVERSAL_WEIGHT
    den = float(np.sum(w))
    return math.sqrt(float(np.sum(w * err * err)) / den) if den > 0 else float("nan")


def score_hysteresis(
    strain: np.ndarray,
    sim_stress: np.ndarray,
    exp_stress: np.ndarray,
    cycle: np.ndarray,
    target_cycles: Sequence[int],
    weight_reversals: bool = True,
) -> float:
    """Mean per-cycle stress RMSE (MPa) between two responses on a common grid.

    Backend-agnostic: ``sim_stress`` may come from the surrogate or from an
    ABAQUS ODB already resampled onto the experimental ``strain`` points. Both
    arrays must be aligned index-for-index with ``strain``/``cycle``.

    Returns :data:`PENALTY` if no cycle yields a finite RMSE.
    """
    rmses: list[float] = []
    for c in target_cycles:
        idx = np.where(cycle == c)[0]
        if idx.size < 4:
            continue
        r = _cycle_rmse(strain[idx], sim_stress[idx], exp_stress[idx],
                        weight_reversals)
        if math.isfinite(r):
            rmses.append(r)
    if not rmses:
        return PENALTY
    return float(np.mean(rmses))


def evaluate_surrogate(
    model,
    params: Sequence[float],
    window: dict,
    gate=None,
    *,
    E: float,
    sy0: float,
    weight_reversals: bool = True,
) -> float:
    """Objective for the **surrogate** backend.

    Runs the model over the (optional) burn-in cycle-1 path plus the scoring
    window so the state (backstresses, accumulated plastic strain) is conditioned
    before scoring starts at cycle 2 - this removes the artificial cycle-2
    overshoot that a virgin start produces.

    Parameters
    ----------
    model : CyclicPlasticityModel
    params : sequence of float, in ``model.param_names`` order.
    window : dict with keys ``strain``, ``stress``, ``cycle``, ``target_cycles``
        and optionally ``burnin_strain`` (see ``data_loader.slice_window``).
    gate : PhysicsGate or None
        If given, invalid candidates short-circuit to :data:`PENALTY`.
    E, sy0 : float
        Per-run elastic modulus and yield stress.

    Returns
    -------
    float objective (MPa), or :data:`PENALTY` on invalidity / non-finite result.
    """
    # Repair (project) infeasible candidates onto the feasible set rather than
    # returning a constant PENALTY -- a flat penalty landscape makes the
    # optimiser lock onto the first sampled point (all evals tie at PENALTY).
    if gate is not None:
        params = gate.project(model, params, sy0)
    try:
        burnin = window.get("burnin_strain")
        if burnin is not None and np.size(burnin) > 0:
            full = np.concatenate([np.asarray(burnin), window["strain"]])
            sig_full = model.simulate(full, params, E=E, sy0=sy0)
            sim = sig_full[np.size(burnin):]
        else:
            sim = model.simulate(window["strain"], params, E=E, sy0=sy0)
    except Exception:
        return PENALTY
    if sim is None or not np.all(np.isfinite(sim)):
        return PENALTY
    return score_hysteresis(window["strain"], sim, window["stress"],
                            window["cycle"], window["target_cycles"],
                            weight_reversals)


def score_from_abaqus(
    sim_strain: np.ndarray,
    sim_stress: np.ndarray,
    window: dict,
    weight_reversals: bool = True,
) -> float:
    """Objective for the **FE backend**: score an ABAQUS response.

    The ABAQUS run produces its own (strain, stress) trace; here it is resampled
    onto the experimental strain points *per cycle* (matched by arc-length so
    the comparison is reversal-aligned) and then passed to the identical
    :func:`score_hysteresis`. This guarantees the FE and surrogate objectives are
    the same function of the underlying response.
    """
    exp_strain = window["strain"]
    exp_cycle = window["cycle"]
    target = window["target_cycles"]
    n_cyc = len(target)
    if n_cyc == 0 or sim_strain.size < 4:
        return PENALTY

    # Split the ABAQUS trace into n_cyc equal chunks (one per scored cycle).
    chunk = sim_strain.size // n_cyc
    if chunk < 4:
        return PENALTY

    resampled = np.empty_like(exp_strain, dtype=np.float64)
    for k, c in enumerate(target):
        idx = np.where(exp_cycle == c)[0]
        if idx.size < 4:
            resampled[idx] = np.nan
            continue
        e_exp = exp_strain[idx]
        sl = slice(k * chunk, (k + 1) * chunk if k < n_cyc - 1 else sim_strain.size)
        e_sim, s_sim = sim_strain[sl], sim_stress[sl]
        arc_exp = _arc(e_exp)
        arc_sim = _arc(e_sim)
        if arc_exp is None or arc_sim is None:
            resampled[idx] = np.nan
            continue
        resampled[idx] = np.interp(arc_exp, arc_sim, s_sim)

    if not np.all(np.isfinite(resampled)):
        return PENALTY
    return score_hysteresis(exp_strain, resampled, window["stress"],
                            exp_cycle, target, weight_reversals)


def _arc(strain_c: np.ndarray):
    """Normalised cumulative arc-length along a strain segment, or None."""
    ds = np.abs(np.diff(strain_c, prepend=strain_c[0]))
    arc = np.cumsum(ds)
    if arc[-1] <= 0:
        return None
    return arc / arc[-1]
