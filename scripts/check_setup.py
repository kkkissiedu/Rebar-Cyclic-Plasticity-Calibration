"""Checks that this machine can run the app and reports which optional tools were found.

    python scripts/check_setup.py
"""

from __future__ import annotations

import importlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

REQUIRED = ("numpy", "scipy", "numba", "matplotlib", "optuna", "tkinter")


def main() -> int:
    ok = True
    print("Python", sys.version.split()[0])
    for name in REQUIRED:
        try:
            mod = importlib.import_module(name)
            print(f"  ok       {name} {getattr(mod, '__version__', '')}".rstrip())
        except Exception as exc:                                    # noqa: BLE001
            ok = False
            print(f"  MISSING  {name}: {exc}")
    from calibration_app import settings
    print("\nExternal tools (only the finite-element parts need them):")
    for key, value in settings.preflight().items():
        print(f"  {key:9s}{value}")
    print("\nSurrogate calibration" + (" is ready." if ok else " cannot run until the missing packages are installed."))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
