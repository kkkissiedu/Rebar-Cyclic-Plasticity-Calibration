"""
main.py
=======

Entry point. Launches the tkinter GUI.

Run (from the project directory, in the cuda_pt conda environment)::

    conda activate cuda_pt
    python -m calibration_app.main
    # or:  python calibration_app/main.py

DPI awareness is enabled on Windows before any Tk widget is created, and the Tk
scaling + named fonts are scaled from the monitor DPI so the UI is crisp on
high-DPI displays.
"""

from __future__ import annotations

import ctypes
import sys


def _enable_dpi_awareness() -> None:
    """Best-effort per-monitor DPI awareness on Windows (no-op elsewhere)."""
    if sys.platform != "win32":
        return
    try:
        ctypes.windll.shcore.SetProcessDpiAwareness(2)   # PER_MONITOR_AWARE
    except Exception:
        try:
            ctypes.windll.user32.SetProcessDPIAware()
        except Exception:
            pass


def main() -> None:
    _enable_dpi_awareness()

    import tkinter as tk
    import tkinter.font as tkfont

    root = tk.Tk()

    # Scale Tk + named fonts from the actual monitor DPI (spec: keep DPI scaling).
    try:
        scale = root.winfo_fpixels("1i") / 72.0
        root.tk.call("tk", "scaling", scale)
        for fname in ("TkDefaultFont", "TkTextFont", "TkFixedFont",
                      "TkMenuFont", "TkHeadingFont"):
            try:
                f = tkfont.nametofont(fname)
                base = abs(f.cget("size")) or 10
                f.configure(size=max(8, int(round(base * min(scale, 1.5)))))
            except Exception:
                pass
    except Exception:
        pass

    from calibration_app.ui.app_window import AppWindow
    AppWindow(root)
    root.mainloop()


if __name__ == "__main__":
    main()
