"""The FE scorer must give ~zero error for an FE trace that equals the test response, whatever the output frame spacing.

    python tests/test_fe_scorer.py

A synthetic hysteretic response is sampled like ABAQUS output (irregular automatic increments, first frame after the
burn-in cutoff not on the cycle boundary) and compared with the same response sampled uniformly. Splitting the trace by
time recovers the match; the old equal point-count split does not.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration_app.core import objective                       # noqa: E402

P, AMP = 20.0, 0.03


def strain_of(t):
    ph = (np.asarray(t, float) % P) / P
    return np.interp(ph, [0, 0.25, 0.5, 0.75, 1.0], [0, AMP, 0, -AMP, 0])


def stress_of(t):
    """Smooth and periodic (no jumps, which no interpolation can reproduce), different on loading and unloading."""
    ph = (np.asarray(t, float) % P) / P
    return 450.0 * np.sin(2 * np.pi * ph) + 120.0 * np.sin(4 * np.pi * ph + 0.5)


def fe_trace(burnin_cycles=1, scored_cycles=3):
    steps, t = [], 0.0
    total = (burnin_cycles + scored_cycles) * P
    while t < total:
        t += 0.13 if len(steps) % 2 == 0 else 0.37
        steps.append(min(t, total))
    tt = np.array([0.0] + steps)
    return tt, strain_of(tt), stress_of(tt)


def window(scored_cycles=3):
    te = np.linspace(0, P, 201)
    return {"strain": np.tile(strain_of(te), scored_cycles), "stress": np.tile(stress_of(te), scored_cycles),
            "cycle": np.repeat(np.arange(2, 2 + scored_cycles), te.size), "target_cycles": list(range(2, 2 + scored_cycles))}


def main() -> int:
    tt, e, s = fe_trace()
    keep = tt >= P - 1e-9
    first = int(np.argmax(keep))
    legacy = (e[keep], s[keep])                       # what the old code scored
    with_prev = slice(first - 1, None)                # what run_candidate(with_time=True) returns now
    w = window()
    old = objective.score_from_abaqus(*legacy, w, False)
    new = objective.score_from_abaqus(e[with_prev], s[with_prev], w, False, sim_time=tt[with_prev], t0=P, period_s=P)
    print(f"equal-chunk split: {old:.2f} MPa   time split: {new:.2f} MPa")
    ok = new < 1.5 and old > 2.0 * new
    # an analysis that stops early must be rejected, not scored on a partial trace
    short = objective.score_from_abaqus(e[with_prev][:-60], s[with_prev][:-60], w, False,
                                        sim_time=tt[with_prev][:-60], t0=P, period_s=P)
    ok = ok and short >= objective.PENALTY
    print("OK" if ok else "FAIL", "(short trace rejected: %s)" % (short >= objective.PENALTY))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
