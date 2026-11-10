"""Regression test: the three constitutive models must give identical stress paths before and after any refactor.

    python tests/test_regression.py            # compare against tests/baseline.json
    python tests/test_regression.py --update   # rewrite the baseline (only when a change in results is intended)

Each case drives a model with a fixed triangular strain path (3 % amplitude, 200 points per cycle, 3 cycles) and the
midpoint of the model's default bounds, then compares summary statistics and a hash of the rounded stress path.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from calibration_app.core.pipeline import build_model            # noqa: E402

BASELINE = os.path.join(HERE, "baseline.json")
CASES = [("chaboche", 1), ("chaboche", 2), ("chaboche", 3), ("uvc", 2), ("uvc", 3), ("ohno_wang", 2)]
E, SY0 = 200_000.0, 300.0


def strain_path(amp: float = 0.03, pts: int = 200, cycles: int = 3) -> np.ndarray:
    one = np.concatenate([np.linspace(0, amp, pts // 4 + 1)[:-1], np.linspace(amp, -amp, pts // 2 + 1)[:-1],
                          np.linspace(-amp, 0, pts // 4 + 1)[:-1]])
    return np.concatenate([np.tile(one, cycles), [0.0]])


def run_case(name: str, n_bs: int) -> dict:
    model = build_model(name, n_bs)
    bounds = model.default_bounds()
    params = tuple(0.5 * (lo + hi) for lo, hi in (bounds[k] for k in model.param_names))
    sig = np.asarray(model.simulate(strain_path(), params, E=E, sy0=SY0), dtype=float)
    return {"peak": float(sig.max()), "trough": float(sig.min()), "rms": float(np.sqrt(np.mean(sig ** 2))),
            "sha": hashlib.sha256(np.round(sig, 6).tobytes()).hexdigest()}


def main() -> int:
    results = {f"{n}_{k}": run_case(n, k) for n, k in CASES}
    if "--update" in sys.argv:
        with open(BASELINE, "w") as f:
            json.dump(results, f, indent=1)
        print("baseline written:", {k: round(v["peak"], 3) for k, v in results.items()})
        return 0
    with open(BASELINE) as f:
        base = json.load(f)
    bad = []
    for key, res in results.items():
        b = base[key]
        same = all(abs(res[s] - b[s]) <= 1e-9 * max(1.0, abs(b[s])) for s in ("peak", "trough", "rms"))
        if not same or res["sha"] != b["sha"]:
            bad.append((key, res["peak"], b["peak"]))
    if bad:
        print("REGRESSION:", bad)
        return 1
    print("OK: all", len(results), "model cases identical to baseline")
    return 0


if __name__ == "__main__":
    sys.exit(main())
