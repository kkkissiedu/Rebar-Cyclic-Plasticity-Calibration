"""
predictor_window.py
===================

"Predict..." window: enter bar diameter, L/D, strain amplitude (and optionally fy and the two material factors) and get the predicted life,
peak-stress history and loops from calibration_app.predictor (V7 framework). Shows the fitted range and every extrapolation warning.
"""

from __future__ import annotations

import json
import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from ..predictor import predict_bar
from ..predictor.export import export_csv, export_png

NOTE = ("Fitted on one steel (B500C, Kashani et al. 2019), D 10-20 mm, L/D 5-15, 1-6 % strain, LVDT strain. Outside this range, or for another steel, "
        "the output is an extrapolation. Material factors f_s, f_e rescale stress and life from a few tests; their use on a different steel is untested.")


class PredictorWindow:
    def __init__(self, parent: tk.Misc):
        self.win = tk.Toplevel(parent)
        self.win.title("Predict life and loops")
        self.win.geometry("1100x640")
        left = ttk.Frame(self.win)
        left.pack(side="left", fill="y", padx=8, pady=8)
        self.vars = {}
        for label, key, default in (("Diameter D (mm)", "d", "16"), ("Slenderness L/D", "ld", "8"), ("Strain amplitude (%)", "amp", "3"),
                                    ("Yield strength fy (MPa)", "fy", "540"), ("Stress factor f_s", "fs", "1.0"), ("Life factor f_e", "fe", "1.0")):
            row = ttk.Frame(left)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=label, width=22).pack(side="left")
            self.vars[key] = tk.StringVar(value=default)
            ttk.Entry(row, textvariable=self.vars[key], width=9).pack(side="left")
        self._presets = self._load_presets()
        prow = ttk.Frame(left)
        prow.pack(fill="x", pady=(6, 2))
        ttk.Label(prow, text="Steel preset", width=22).pack(side="left")
        self.preset = tk.StringVar(value="Kashani B500C (fitted steel)")
        box = ttk.Combobox(prow, textvariable=self.preset, state="readonly", width=26,
                           values=["Kashani B500C (fitted steel)"] + [p["name"] for p in self._presets.get("presets", [])])
        box.pack(side="left")
        box.bind("<<ComboboxSelected>>", lambda _e: self._apply_preset())
        self.basis = tk.StringVar(value="lvdt")
        ttk.Radiobutton(left, text="Strain basis: LVDT (machine compliance)", variable=self.basis, value="lvdt").pack(anchor="w", pady=(6, 0))
        ttk.Radiobutton(left, text="Strain basis: extensometer", variable=self.basis, value="extensometer").pack(anchor="w")
        row = ttk.Frame(left)
        row.pack(fill="x", pady=2)
        ttk.Label(row, text="Cycles (blank = auto)", width=22).pack(side="left")
        self.vars["cycles"] = tk.StringVar(value="")
        ttk.Entry(row, textvariable=self.vars["cycles"], width=14).pack(side="left")
        self.result = None
        ttk.Button(left, text="Predict", style="Accent.TButton", command=self.run).pack(fill="x", pady=(8, 2))
        brow = ttk.Frame(left)
        brow.pack(fill="x", pady=(0, 8))
        ttk.Button(brow, text="Export CSV...", command=self._export_csv).pack(side="left", expand=True, fill="x", padx=(0, 2))
        ttk.Button(brow, text="Export PNG...", command=self._export_png).pack(side="left", expand=True, fill="x", padx=(2, 0))
        self.text = tk.Text(left, width=44, height=22, wrap="word")
        self.text.pack(fill="both", expand=True)
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.win)
        self.canvas.get_tk_widget().pack(side="right", fill="both", expand=True)
        self._say(NOTE)

    def _say(self, s: str) -> None:
        self.text.configure(state="normal")
        self.text.delete("1.0", "end")
        self.text.insert("1.0", s)
        self.text.configure(state="disabled")

    @staticmethod
    def _load_presets() -> dict:
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "predictor", "scrap_presets.json")
        try:
            with open(path, encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    def _apply_preset(self) -> None:
        """Fill f_s and f_e from the chosen preset (1 and 1 for the fitted steel) and say plainly what they are."""
        name = self.preset.get()
        hit = next((p for p in self._presets.get("presets", []) if p["name"] == name), None)
        self.vars["fs"].set(f"{hit['f_s_rm']:.3f}" if hit else "1.0")
        self.vars["fe"].set(f"{hit['f_e']:g}" if hit else "1.0")
        if not hit:
            self._say(NOTE)
            return
        a = self._presets["anchor"]
        self._say(f"{name}: f_s = {hit['f_s_rm']:.3f} (peak stress {hit['rm_mean_mpa']} MPa from {hit['n_bars_rm']} tension tests against "
                  f"fu {a['fu_mpa']:.0f} MPa); the yield-based value is {hit['f_s_ry']:.3f} ({hit['n_bars_ry']} bars). f_e = 1, no life factor "
                  f"is available.\n\n{self._presets['status']}\n\nAnchor: {a['source']}. {a['note']}\n\nPress Predict.")

    def _default_name(self) -> str:
        i = self.result["inputs"]
        return f"predicted_D{i['diameter_mm']:g}_LD{i['L_over_D']:g}_{i['amplitude_pct']:g}pct"

    def _export_csv(self) -> None:
        if self.result is None:
            self._say("Run Predict first.")
            return
        path = filedialog.asksaveasfilename(parent=self.win, defaultextension=".csv", initialfile=self._default_name() + ".csv",
                                            filetypes=[("CSV", "*.csv")])
        if path:
            export_csv(self.result, path)
            messagebox.showinfo("Export CSV", f"Saved {path}", parent=self.win)

    def _export_png(self) -> None:
        if self.result is None:
            self._say("Run Predict first.")
            return
        path = filedialog.asksaveasfilename(parent=self.win, defaultextension=".png", initialfile=self._default_name() + ".png",
                                            filetypes=[("PNG", "*.png")])
        if path:
            files = export_png(self.result, path)
            messagebox.showinfo("Export PNG", "Saved:\n" + "\n".join(files), parent=self.win)

    def run(self) -> None:
        try:
            v = {k: float(x.get()) for k, x in self.vars.items() if k != "cycles"}
            txt = self.vars["cycles"].get().replace(";", ",").strip()
            cyc = [int(c) for c in txt.split(",") if c.strip()] if txt else None
            r = predict_bar(v["d"], v["ld"], v["amp"], v["fy"], cycles=cyc, strain_basis=self.basis.get(), factors=(v["fs"], v["fe"]))
        except Exception as exc:                                   # bad input: show, do not crash the app
            self._say(f"Cannot predict: {exc}")
            return
        self.result = r
        life = r["life"]
        lines = [f"Cycles to failure: {life['cycles_to_failure']:.0f} (factor-2 band {life['factor2_band'][0]:.0f} to {life['factor2_band'][1]:.0f}; "
                 f"log10 scatter {life['log10_scatter']:.2f})",
                 f"Cycle-2 peaks: tension {r['cycle2_peaks']['tension_mpa']:.0f} MPa, compression {r['cycle2_peaks']['compression_mpa']:.0f} MPa",
                 f"Apparent modulus {r['apparent_modulus_gpa']:.0f} GPa, lambda_p {r['lambda_p']:.1f}", ""]
        lines += [f"WARNING: {w}" for w in r["warnings"]] + ["", NOTE]
        self._say("\n".join(lines))
        self.fig.clear()
        a1, a2 = self.fig.subplots(1, 2)
        for c, lp in r["loops"].items():
            a1.plot(lp["strain"] * 100, lp["stress"], lw=1.2, label=f"cycle {c}")
        a1.set_xlabel("strain (%)"); a1.set_ylabel("stress (MPa)"); a1.set_title("Predicted loops"); a1.grid(alpha=0.25); a1.legend(fontsize=7)
        a2.plot(r["peaks"]["cycle"], r["peaks"]["tension"], label="tension")
        a2.plot(r["peaks"]["cycle"], [-x for x in r["peaks"]["compression"]], label="compression")
        a2.set_xlabel("cycle"); a2.set_ylabel("peak stress (MPa)"); a2.set_title("Peak-stress history"); a2.grid(alpha=0.25); a2.legend(fontsize=7)
        self.fig.tight_layout()
        self.canvas.draw_idle()
