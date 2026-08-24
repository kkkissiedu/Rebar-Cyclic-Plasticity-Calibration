"""
ui/theme.py
===========

Dark/light theming for the tkinter GUI, with the preference persisted between
sessions.

If the optional ``sv_ttk`` package is installed in the environment it is used for
a modern look; otherwise a manual ttk palette is applied so the app still themes
correctly with only the standard library (the spec's "fall back to ttk styled
manually"). The choice is transparent to the rest of the UI, which only calls
:func:`apply_theme` and reads colours from :data:`PALETTES`.
"""

from __future__ import annotations

import json
import os
from typing import Literal

Mode = Literal["dark", "light"]

#: Where the theme preference is stored (per-user).
_PREF_PATH = os.path.join(os.path.expanduser("~"), ".chaboche_calib_ui.json")

#: Colour palettes for the two modes. Panels read these for matplotlib etc.
PALETTES: dict[str, dict[str, str]] = {
    "dark": {
        "bg": "#1e1e1e", "panel": "#252526", "fg": "#e0e0e0",
        "muted": "#a0a0a0", "accent": "#4ea1ff", "field": "#2d2d30",
        "grid": "#3a3a3a",
        "exp": "#9a9a9a", "best": "#ff5555", "last": "#4ea1ff", "fe": "#2ecc71",
        "ok": "#4caf50", "warn": "#e0a000", "fail": "#e05555",
    },
    "light": {
        "bg": "#f4f4f4", "panel": "#ffffff", "fg": "#1a1a1a",
        "muted": "#606060", "accent": "#0a66c2", "field": "#ffffff",
        "grid": "#cccccc",
        "exp": "#777777", "best": "#d00000", "last": "#0a66c2", "fe": "#1e8e3e",
        "ok": "#2e7d32", "warn": "#a06000", "fail": "#c62828",
    },
}


def load_pref() -> Mode:
    """Return the saved theme mode, defaulting to 'dark'."""
    try:
        with open(_PREF_PATH, "r") as f:
            mode = json.load(f).get("theme", "dark")
            return mode if mode in PALETTES else "dark"
    except (OSError, ValueError):
        return "dark"


def save_pref(mode: Mode) -> None:
    """Persist the chosen theme mode."""
    try:
        with open(_PREF_PATH, "w") as f:
            json.dump({"theme": mode}, f)
    except OSError:
        pass


def apply_theme(root, style, mode: Mode) -> dict[str, str]:
    """Apply ``mode`` to the whole widget tree; return the active palette.

    Parameters
    ----------
    root : tk.Tk
    style : ttk.Style
    mode : "dark" | "light"
    """
    pal = PALETTES[mode]

    used_sv = False
    try:
        import sv_ttk  # optional, modern theme
        sv_ttk.set_theme(mode)
        used_sv = True
    except Exception:
        used_sv = False

    if not used_sv:
        # Manual ttk palette on the 'clam' base (present everywhere).
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure(".", background=pal["panel"], foreground=pal["fg"],
                        fieldbackground=pal["field"], bordercolor=pal["grid"])
        style.configure("TFrame", background=pal["panel"])
        style.configure("TLabelframe", background=pal["panel"],
                        foreground=pal["fg"])
        style.configure("TLabelframe.Label", background=pal["panel"],
                        foreground=pal["muted"])
        style.configure("TLabel", background=pal["panel"], foreground=pal["fg"])
        style.configure("TButton", background=pal["field"], foreground=pal["fg"])
        style.map("TButton",
                  background=[("active", pal["accent"])],
                  foreground=[("active", "#ffffff")])
        style.configure("TCheckbutton", background=pal["panel"],
                        foreground=pal["fg"])
        style.configure("TRadiobutton", background=pal["panel"],
                        foreground=pal["fg"])
        style.configure("Accent.TButton", background=pal["accent"],
                        foreground="#ffffff")
        style.configure("Horizontal.TProgressbar", background=pal["accent"],
                        troughcolor=pal["field"])

    root.configure(bg=pal["bg"])

    # -- Treeview readability + hover/selection legibility (applied in BOTH the
    # sv_ttk and manual-clam cases; later configure/map wins over the base theme).
    # Rows were too short and text too small in the sessions browser; on
    # hover/selection the background lightened but the text stayed light and
    # became unreadable. Force a taller row, a >=10pt font, and BLACK text on any
    # selected/active/pressed state over a mid-blue highlight (readable in both
    # dark and light themes).
    style.configure("Treeview", rowheight=24, font=("Segoe UI", 10),
                    background=pal["field"], fieldbackground=pal["field"],
                    foreground=pal["fg"])
    style.configure("Treeview.Heading", font=("Segoe UI", 10, "bold"))
    style.map("Treeview",
              background=[("selected", "#4a90d9")],
              foreground=[("selected", "black")])
    style.map("TButton",
              foreground=[("active", "black"), ("pressed", "black")])
    # Entry / Combobox / Spinbox active/selected states share the same problem.
    style.map("TEntry", foreground=[("active", "black")])
    style.map("TCombobox",
              foreground=[("active", "black"), ("focus", "black")],
              selectforeground=[("!disabled", "black")],
              selectbackground=[("!disabled", "#4a90d9")])
    style.map("TSpinbox", foreground=[("active", "black")])

    # Status-badge styles (used by progress_panel) in every case.
    for name, col in (("Ok", pal["ok"]), ("Warn", pal["warn"]),
                      ("Fail", pal["fail"]), ("Muted", pal["muted"]),
                      ("Accent", pal["accent"])):
        style.configure(f"{name}.TLabel", background=pal["panel"], foreground=col)
    # Yellow highlight badge for out-of-literature bounds (config_panel bounds
    # editor). Amber background so the row is visibly flagged without blocking.
    style.configure("Badge.TLabel", background=pal["warn"], foreground="#1a1a1a")
    return pal


def toggle(mode: Mode) -> Mode:
    """Return the opposite mode and persist it."""
    new: Mode = "light" if mode == "dark" else "dark"
    save_pref(new)
    return new
