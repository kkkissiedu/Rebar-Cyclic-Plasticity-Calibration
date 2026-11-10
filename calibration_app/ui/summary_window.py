"""
summary_window.py
=================

One results view shown when a run finishes: the fitted loops against the test, the error of each cycle, the peak stress
of each cycle (tension and compression), the best parameters and, when present, the ABAQUS check. The figure is also
saved as ``summary.png`` in the session folder, so the same view can be reopened or attached to a report.

``build_summary`` does the computation without any GUI (used by the tests); ``show_summary`` draws it.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tkinter as tk
from tkinter import ttk
from typing import Optional

import numpy as np

from ..core import data_loader, pipeline as pipeline_mod

CAVEAT = ("Fit quality is shown on the calibration window only. Fitted values are calibration results, not material "
          "constants, unless they are checked against independent data.")


def build_summary(session_dir: str) -> dict:
    """Load a finished session and compute everything the summary shows."""
    with open(os.path.join(session_dir, "session.json"), encoding="utf-8") as f:
        sj = json.load(f)
    cfg = sj.get("config", {})
    out = {"session": sj.get("name", os.path.basename(session_dir)), "model": sj.get("model"),
           "backend": sj.get("backend"), "objective": sj.get("best_objective"),
           "sigma_y0": sj.get("best_sigma_y0") or cfg.get("sigma_y0"), "params": sj.get("best_params"),
           "param_names": None, "curves": None, "cycles": None, "fe": sj.get("fe_verification"), "note": ""}
    if not out["params"] or not out["model"]:
        out["note"] = "No best parameters were found in this session."
        return out
    model = pipeline_mod.build_model(out["model"], int(cfg.get("n_backstresses", 3)))
    out["param_names"] = list(model.param_names)
    path = cfg.get("data_file", "")
    if not path or not os.path.exists(path):
        out["note"] = f"Data file not found, so no plots: {path}"
        return out
    data = data_loader.load_any_experimental(
        path, gauge_length_mm=float(cfg["gauge_length_mm"]), bar_diameter_mm=float(cfg["diameter_mm"]),
        sigma_y0=float(cfg.get("sigma_y0", 500.0)))
    dw = data_loader.slice_window(data, int(cfg.get("n_fit", 3)))
    E = float(cfg.get("E_modulus") or 200_000.0)
    sy0 = float(out["sigma_y0"])
    strain, burn = dw["strain"], np.asarray(dw.get("burnin_strain", []), float)
    if burn.size:
        sim = model.simulate(np.concatenate([burn, strain]), out["params"], E=E, sy0=sy0)[burn.size:]
    else:
        sim = model.simulate(strain, out["params"], E=E, sy0=sy0)
    sim = np.asarray(sim, float)
    exp = np.asarray(dw["stress"], float)
    cyc = np.asarray(dw["cycle"])
    rows = []
    for c in dw["target_cycles"]:
        m = cyc == c
        if not m.any():
            continue
        rows.append({"cycle": int(c), "rmse": float(np.sqrt(np.mean((sim[m] - exp[m]) ** 2))),
                     "exp_t": float(exp[m].max()), "exp_c": float(-exp[m].min()),
                     "fit_t": float(sim[m].max()), "fit_c": float(-sim[m].min())})
    out["curves"] = {"strain": strain, "exp": exp, "fit": sim, "cycle": cyc}
    out["cycles"] = rows
    return out


def make_figure(s: dict):
    from matplotlib.figure import Figure
    fig = Figure(figsize=(10.5, 6.4), dpi=100)
    if not s.get("curves"):
        ax = fig.add_subplot(111)
        ax.axis("off")
        ax.text(0.5, 0.5, s.get("note") or "Nothing to plot.", ha="center", va="center", wrap=True)
        return fig
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], hspace=0.35, wspace=0.28)
    c = s["curves"]
    a1 = fig.add_subplot(gs[:, 0])
    a1.plot(c["strain"] * 100, c["exp"], color="k", lw=1.0, label="test")
    a1.plot(c["strain"] * 100, c["fit"], color="#D55E00", lw=1.0, label="best fit")
    a1.set_xlabel("strain (%)")
    a1.set_ylabel("stress (MPa)")
    a1.set_title("Loops, calibration window")
    a1.grid(alpha=0.25)
    a1.legend(fontsize=8)
    rows = s["cycles"]
    cs = [r["cycle"] for r in rows]
    a2 = fig.add_subplot(gs[0, 1])
    a2.bar(cs, [r["rmse"] for r in rows], color="#0072B2")
    a2.set_xlabel("cycle")
    a2.set_ylabel("RMSE (MPa)")
    a2.set_title("Error per cycle")
    a2.grid(axis="y", alpha=0.25)
    a3 = fig.add_subplot(gs[1, 1])
    a3.plot(cs, [r["exp_t"] for r in rows], "ko-", ms=3, label="test, tension")
    a3.plot(cs, [r["exp_c"] for r in rows], "ks--", ms=3, mfc="none", label="test, compression")
    a3.plot(cs, [r["fit_t"] for r in rows], "-", color="#D55E00", label="fit, tension")
    a3.plot(cs, [r["fit_c"] for r in rows], "--", color="#D55E00", label="fit, compression")
    a3.set_xlabel("cycle")
    a3.set_ylabel("peak stress (MPa)")
    a3.set_title("Peak stress per cycle")
    a3.grid(alpha=0.25)
    a3.legend(fontsize=7)
    return fig


def _open_folder(path: str) -> None:
    try:
        if sys.platform == "win32":
            os.startfile(path)                                       # noqa: S606
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except OSError:
        pass


def _text(s: dict) -> str:
    lines = [f"Session: {s['session']}", f"Model: {s['model']}   backend: {s['backend']}"]
    if s.get("objective") is not None:
        lines.append(f"Best objective: {float(s['objective']):.2f} MPa")
    if s.get("sigma_y0"):
        lines.append(f"Initial yield: {float(s['sigma_y0']):.1f} MPa")
    for n, v in zip(s.get("param_names") or [], s.get("params") or []):
        lines.append(f"  {n} = {v:.5g}")
    fe = s.get("fe")
    if fe and fe.get("objective_MPa") is not None:
        lines.append(f"ABAQUS check objective: {float(fe['objective_MPa']):.2f} MPa")
    if s.get("note"):
        lines.append(s["note"])
    return "\n".join(lines)


def show_summary(parent: tk.Misc, session_dir: str) -> Optional[tk.Toplevel]:
    """Open the summary window for a finished session and save summary.png beside session.json."""
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    s = build_summary(session_dir)
    fig = make_figure(s)
    png = os.path.join(session_dir, "summary.png")
    try:
        fig.savefig(png, dpi=130)
    except OSError:
        png = ""
    win = tk.Toplevel(parent)
    win.title(f"Results: {s['session']}")
    win.geometry("1180x760")
    canvas = FigureCanvasTkAgg(fig, master=win)
    canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
    canvas.draw()
    box = ttk.Frame(win)
    box.pack(side="bottom", fill="x", padx=8, pady=6)
    ttk.Label(box, text=_text(s), justify="left", font=("Consolas", 9)).pack(side="left", anchor="n")
    right = ttk.Frame(box)
    right.pack(side="right", anchor="n")
    ttk.Label(right, text=CAVEAT, wraplength=380, justify="left").pack(anchor="e", pady=(0, 6))
    ttk.Button(right, text="Open results folder", command=lambda: _open_folder(session_dir)).pack(side="left", padx=4)
    ttk.Button(right, text="Close", command=win.destroy).pack(side="left", padx=4)
    return win
