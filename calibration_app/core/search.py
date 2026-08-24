"""
core/search.py
==============

Backend-agnostic optimisers for parameter identification.

The optimisers never know whether they are driving the pure-Python surrogate or
a real ABAQUS FE evaluation — they are handed a single ``evaluate(params) ->
float`` callable plus the bounds. This is what makes the four study conditions
(Surrogate/FE x DE/Bayesian) share one code path with an identical objective.

Provided
--------
* :func:`sobol_seed`  — space-filling Stage-1 sampling (replaces the old grid,
  which does not scale past ~6 parameters). Returns the best-K seeds.
* :func:`run_differential_evolution` — SciPy DE (Stage-2 global refinement).
* :func:`run_bayesian` — Optuna TPE Bayesian optimisation (default for the
  expensive FE backend; DE remains selectable).

All three accept a ``progress`` callback (called once per evaluation with an
:class:`EvalRecord`) and a ``stop`` :class:`threading.Event` for graceful abort,
so the GUI can show live progress and cancel cleanly.
"""

from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np

EvaluateFn = Callable[[Sequence[float]], float]


@dataclass
class EvalRecord:
    """One objective evaluation, emitted to the ``progress`` callback."""
    params: tuple[float, ...]
    objective: float
    stage: str                      # "sobol" | "de" | "bayesian"
    index: int                      # running evaluation counter within the stage
    best_objective: float           # best seen so far (this run)
    wall_s: float


@dataclass
class SearchResult:
    """Outcome of a search stage."""
    best_params: tuple[float, ...]
    best_objective: float
    n_evals: int
    stage: str
    stopped: bool = False
    history: list[EvalRecord] = field(default_factory=list)


class _StopSearch(Exception):
    """Raised internally to unwind an optimiser when the stop event is set."""


def _bounds_arrays(bounds: Sequence[tuple[float, float]]):
    lo = np.array([b[0] for b in bounds], dtype=np.float64)
    hi = np.array([b[1] for b in bounds], dtype=np.float64)
    return lo, hi


def sobol_seed(
    evaluate: EvaluateFn,
    bounds: Sequence[tuple[float, float]],
    n_samples: int = 4096,
    top_k: int = 10,
    *,
    progress: Optional[Callable[[EvalRecord], None]] = None,
    stop: Optional[threading.Event] = None,
    seed: int = 42,
) -> tuple[list[tuple[float, ...]], SearchResult]:
    """Stage 1: evaluate a Sobol space-filling sample; return the best-K seeds.

    Scales to any parameter count (unlike a grid). ``n_samples`` is rounded up to
    the next power of two internally for Sobol balance.

    Returns
    -------
    (seeds, result) : the top-``top_k`` parameter tuples (best first) and a
        :class:`SearchResult` summarising the stage.
    """
    from scipy.stats import qmc

    dim = len(bounds)
    lo, hi = _bounds_arrays(bounds)
    m = int(np.ceil(np.log2(max(n_samples, 2))))
    sampler = qmc.Sobol(d=dim, scramble=True, seed=seed)
    unit = sampler.random_base2(m=m)
    pts = qmc.scale(unit, lo, hi)

    best_obj = float("inf")
    best_params: tuple[float, ...] = tuple(pts[0])
    evaluated: list[tuple[float, tuple[float, ...]]] = []
    history: list[EvalRecord] = []
    stopped = False

    for i, row in enumerate(pts):
        if stop is not None and stop.is_set():
            stopped = True
            break
        p = tuple(float(x) for x in row)
        t0 = time.perf_counter()
        obj = float(evaluate(p))
        wall = time.perf_counter() - t0
        evaluated.append((obj, p))
        if obj < best_obj:
            best_obj, best_params = obj, p
        rec = EvalRecord(p, obj, "sobol", i + 1, best_obj, wall)
        history.append(rec)
        if progress is not None:
            progress(rec)

    evaluated.sort(key=lambda t: t[0])
    seeds = [p for _, p in evaluated[:top_k]]
    result = SearchResult(best_params, best_obj, len(evaluated), "sobol",
                          stopped, history)
    return seeds, result


def run_differential_evolution(
    evaluate: EvaluateFn,
    bounds: Sequence[tuple[float, float]],
    *,
    x0: Optional[Sequence[float]] = None,
    maxiter: int = 300,
    popsize: int = 20,
    tol: float = 1e-4,
    seed: int = 42,
    progress: Optional[Callable[[EvalRecord], None]] = None,
    stop: Optional[threading.Event] = None,
) -> SearchResult:
    """Stage 2: SciPy differential evolution (reliable global refinement).

    ``workers=1`` deliberately — the objective closes over model/data and is not
    guaranteed picklable on Windows spawn; parallelism is handled at the
    ABAQUS-job level in the FE backend instead.
    """
    from scipy.optimize import differential_evolution

    counter = {"n": 0}
    best = {"obj": float("inf"), "p": tuple(x0) if x0 is not None else None}
    history: list[EvalRecord] = []

    def wrapped(x: np.ndarray) -> float:
        if stop is not None and stop.is_set():
            raise _StopSearch()
        p = tuple(float(v) for v in x)
        t0 = time.perf_counter()
        obj = float(evaluate(p))
        wall = time.perf_counter() - t0
        counter["n"] += 1
        if obj < best["obj"]:
            best["obj"], best["p"] = obj, p
        rec = EvalRecord(p, obj, "de", counter["n"], best["obj"], wall)
        history.append(rec)
        if progress is not None:
            progress(rec)
        return obj

    stopped = False
    try:
        differential_evolution(
            wrapped, list(bounds),
            maxiter=maxiter, popsize=popsize, tol=tol,
            workers=1, polish=True, init="sobol", seed=seed,
            x0=list(x0) if x0 is not None else None,
        )
    except _StopSearch:
        stopped = True

    bp = best["p"] if best["p"] is not None else tuple(b[0] for b in bounds)
    return SearchResult(bp, best["obj"], counter["n"], "de", stopped, history)


def run_bayesian(
    evaluate: EvaluateFn,
    bounds: Sequence[tuple[float, float]],
    param_names: Sequence[str],
    *,
    n_trials: int = 200,
    x0: Optional[Sequence[float]] = None,
    seed: int = 42,
    progress: Optional[Callable[[EvalRecord], None]] = None,
    stop: Optional[threading.Event] = None,
) -> SearchResult:
    """Bayesian optimisation via Optuna (TPE sampler).

    Default optimiser for the FE backend, where each evaluation is an expensive
    ABAQUS run and sample efficiency matters. If ``x0`` is given it is enqueued
    as the first trial so the search starts from the Sobol/DE best.
    """
    import optuna

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    sampler = optuna.samplers.TPESampler(seed=seed)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    if x0 is not None:
        study.enqueue_trial({n: float(v) for n, v in zip(param_names, x0)})

    counter = {"n": 0}
    best = {"obj": float("inf"), "p": None}
    history: list[EvalRecord] = []

    def objective(trial: "optuna.Trial") -> float:
        p = tuple(trial.suggest_float(name, lo, hi)
                  for name, (lo, hi) in zip(param_names, bounds))
        t0 = time.perf_counter()
        obj = float(evaluate(p))
        wall = time.perf_counter() - t0
        counter["n"] += 1
        if obj < best["obj"]:
            best["obj"], best["p"] = obj, p
        rec = EvalRecord(p, obj, "bayesian", counter["n"], best["obj"], wall)
        history.append(rec)
        if progress is not None:
            progress(rec)
        return obj

    def stop_cb(study_, trial_) -> None:
        if stop is not None and stop.is_set():
            study_.stop()

    study.optimize(objective, n_trials=n_trials, callbacks=[stop_cb])

    stopped = stop is not None and stop.is_set()
    bp = best["p"] if best["p"] is not None else tuple(b[0] for b in bounds)
    return SearchResult(bp, best["obj"], counter["n"], "bayesian", stopped, history)


#: Optimiser registry for the GUI selector.
OPTIMISERS: dict[str, str] = {
    "de": "Differential Evolution (SciPy)",
    "bayesian": "Bayesian / TPE (Optuna)",
}
