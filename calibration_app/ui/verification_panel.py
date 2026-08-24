"""
ui/verification_panel.py
=========================

Collapsible panel below the plot showing the ABAQUS-vs-surrogate agreement --
the source for the paper's Table R3 (FE calibrated parameters + fit errors).

Two views, selected by the same radio row that drives :meth:`PlotPanel.set_view`:

* "All models"  -- one row per model: surrogate/FE objective + peak T/C.
* single model  -- Metric/Surrogate/ABAQUS FE/Difference detail table for
  peak tension, peak compression, objective and calibrated sigma_y0.

Data arrives incrementally (surrogate summary at ``model_done``, FE summary at
``stage3_result``, which may not run at all e.g. UMAT unavailable) -- missing
fields render as "-" rather than blocking the table.
"""

from __future__ import annotations

import csv
from tkinter import ttk, filedialog
from typing import Optional


class VerificationPanel:
    """Model-comparison + single-model FE-verification table, with CSV export."""

    def __init__(self, parent, on_view_change=None) -> None:
        self.frame = ttk.LabelFrame(parent, text="Verification")
        self._models: dict[str, dict] = {}   # key -> {"surr": {...}|None, "fe": {...}|None}
        self._view = "all"
        self._collapsed = False
        self._on_view_change = on_view_change

        header = ttk.Frame(self.frame)
        header.pack(fill="x", padx=4, pady=(2, 0))
        self.var_view = None
        self._btn_collapse = ttk.Button(header, text="▲ Collapse", width=10,
                                        command=self._toggle_collapse)
        self._btn_collapse.pack(side="right")
        ttk.Button(header, text="Export verification table...",
                  command=self._export).pack(side="right", padx=4)

        self._view_frame = ttk.Frame(header)
        self._view_frame.pack(side="left")

        self._body = ttk.Frame(self.frame)
        self._body.pack(fill="both", expand=True, padx=4, pady=4)
        self.tree = ttk.Treeview(self._body, show="headings", height=5)
        self.tree.pack(fill="both", expand=True)

        self._render_view_selector(["all"])
        self._render()

    # -- wiring ----------------------------------------------------------
    def set_available_models(self, keys: list[str]) -> None:
        """(Re)build the view-selector radios for the models in this run."""
        self._render_view_selector(keys)

    def _render_view_selector(self, keys: list[str]) -> None:
        import tkinter as tk
        for w in self._view_frame.winfo_children():
            w.destroy()
        self.var_view = tk.StringVar(value=self._view)
        ttk.Label(self._view_frame, text="View:").pack(side="left", padx=(0, 4))
        opts = [("all", "All models")] + [(k, self._label(k)) for k in keys]
        for val, text in opts:
            ttk.Radiobutton(self._view_frame, text=text, value=val,
                            variable=self.var_view,
                            command=lambda: self.set_view(self.var_view.get())
                            ).pack(side="left")

    @staticmethod
    def _label(key: str) -> str:
        return {"chaboche": "Chaboche", "uvc": "UVC",
                "ohno_wang": "Ohno-Wang"}.get(key, key)

    # -- data updates ------------------------------------------------------
    def clear(self) -> None:
        self._models = {}
        self._render()

    def update_surrogate(self, key: str, *, objective: Optional[float],
                         peak_tension: Optional[float],
                         peak_compression: Optional[float],
                         sigma_y0: Optional[float]) -> None:
        m = self._models.setdefault(key, {"surr": None, "fe": None})
        m["surr"] = {"objective": objective, "peak_tension": peak_tension,
                     "peak_compression": peak_compression, "sigma_y0": sigma_y0}
        self._render()

    def update_fe(self, key: str, *, objective: Optional[float],
                  peak_tension: Optional[float],
                  peak_compression: Optional[float]) -> None:
        m = self._models.setdefault(key, {"surr": None, "fe": None})
        m["fe"] = {"objective": objective, "peak_tension": peak_tension,
                   "peak_compression": peak_compression}
        self._render()

    # -- view ----------------------------------------------------------------
    def set_view(self, key: str, *, notify: bool = True) -> None:
        self._view = key
        if self.var_view is not None:
            self.var_view.set(key)
        self._render()
        if notify and self._on_view_change is not None:
            self._on_view_change(key)

    def _toggle_collapse(self) -> None:
        self._collapsed = not self._collapsed
        if self._collapsed:
            self._body.pack_forget()
            self._btn_collapse.configure(text="▼ Expand")
        else:
            self._body.pack(fill="both", expand=True, padx=4, pady=4)
            self._btn_collapse.configure(text="▲ Collapse")

    # -- rendering -------------------------------------------------------
    @staticmethod
    def _fmt(v, unit="") -> str:
        return "-" if v is None else f"{v:.1f}{unit}"

    def _render(self) -> None:
        for c in self.tree["columns"]:
            self.tree.heading(c, text="")
        self.tree.delete(*self.tree.get_children())

        if self._view == "all":
            cols = ("model", "surr_obj", "fe_obj", "peak_t", "peak_c")
            heads = ("Model", "Surr obj (MPa)", "FE obj (MPa)", "Peak T (MPa)", "Peak C (MPa)")
            self.tree["columns"] = cols
            for c, h in zip(cols, heads):
                self.tree.heading(c, text=h)
                self.tree.column(c, width=110, anchor="center")
            for key, m in self._models.items():
                surr, fe = m.get("surr") or {}, m.get("fe") or {}
                self.tree.insert("", "end", values=(
                    self._label(key), self._fmt(surr.get("objective")),
                    self._fmt(fe.get("objective")), self._fmt(surr.get("peak_tension")),
                    self._fmt(surr.get("peak_compression"))))
        else:
            cols = ("metric", "surrogate", "fe", "diff")
            heads = ("Metric", "Surrogate", "ABAQUS FE", "Difference")
            self.tree["columns"] = cols
            for c, h in zip(cols, heads):
                self.tree.heading(c, text=h)
                self.tree.column(c, width=130, anchor="center")
            m = self._models.get(self._view, {})
            surr, fe = m.get("surr") or {}, m.get("fe") or {}
            rows = [
                ("Peak tension", surr.get("peak_tension"), fe.get("peak_tension"), "MPa"),
                ("Peak compression", surr.get("peak_compression"), fe.get("peak_compression"), "MPa"),
                ("Objective", surr.get("objective"), fe.get("objective"), ""),
                ("sigma_y0 calibrated", surr.get("sigma_y0"), None, "MPa"),
            ]
            for label, sv, fv, unit in rows:
                diff = None
                if sv is not None and fv is not None:
                    diff = fv - sv
                diff_s = "-" if diff is None else f"{diff:+.1f}{unit}"
                self.tree.insert("", "end", values=(
                    label, self._fmt(sv, unit), self._fmt(fv, unit), diff_s))

    # -- export ------------------------------------------------------------
    def _export(self) -> None:
        path = filedialog.asksaveasfilename(
            defaultextension=".csv", initialfile="verification_table.csv",
            filetypes=[("CSV", "*.csv")])
        if not path:
            return
        cols = self.tree["columns"]
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow([self.tree.heading(c)["text"] for c in cols])
            for iid in self.tree.get_children():
                w.writerow(self.tree.item(iid)["values"])
