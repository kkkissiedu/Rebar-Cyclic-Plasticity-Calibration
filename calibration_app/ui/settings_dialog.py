"""
settings_dialog.py
==================

A small dialog to point the app at ABAQUS, Visual Studio and the results folder on this machine. Values are saved in
``~/.rebar_calibrator.json`` through :mod:`calibration_app.settings`; environment variables still take precedence.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, ttk

from .. import settings

ROWS = (
    ("abaqus_cmd", "ABAQUS command", "file"),
    ("vcvars", "Visual Studio vcvars64.bat", "file"),
    ("sessions_dir", "Results folder", "dir"),
)


def open_settings(parent: tk.Misc) -> tk.Toplevel:
    win = tk.Toplevel(parent)
    win.title("Settings")
    win.transient(parent)
    values = {}
    for r, (key, label, kind) in enumerate(ROWS):
        ttk.Label(win, text=label).grid(row=r, column=0, sticky="w", padx=8, pady=4)
        var = tk.StringVar(value=settings.get(key) or "")
        values[key] = var
        ttk.Entry(win, textvariable=var, width=62).grid(row=r, column=1, padx=4)

        def browse(k=key, kd=kind):
            path = filedialog.askdirectory() if kd == "dir" else filedialog.askopenfilename()
            if path:
                values[k].set(path)

        ttk.Button(win, text="Browse...", command=browse).grid(row=r, column=2, padx=4)

    status = tk.StringVar()
    ttk.Label(win, textvariable=status, justify="left", wraplength=660).grid(
        row=len(ROWS), column=0, columnspan=3, sticky="w", padx=8, pady=6)

    def refresh() -> None:
        pf = settings.preflight()
        status.set("In use now (environment variables override saved values):\n"
                   f"  ABAQUS: {pf['abaqus']}\n  Visual Studio: {pf['msvc']}\n  Results: {pf['sessions']}\n"
                   "A changed results folder takes effect the next time the app starts. Leave a field empty to "
                   "use automatic discovery.")

    def save() -> None:
        for key, _, _ in ROWS:
            settings.set_value(key, values[key].get().strip() or None)
        refresh()

    buttons = ttk.Frame(win)
    buttons.grid(row=len(ROWS) + 1, column=0, columnspan=3, sticky="e", padx=8, pady=8)
    ttk.Button(buttons, text="Save and check", command=save).pack(side="left", padx=4)
    ttk.Button(buttons, text="Close", command=win.destroy).pack(side="left", padx=4)
    refresh()
    return win
