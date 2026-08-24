"""
ui/param_panel.py
=================

Side-by-side display of the last-evaluated and best-so-far parameters, with a
"Copy best to clipboard" button that emits a ready-to-paste ABAQUS material
block for the active model.

The rows are built dynamically from ``model.param_names`` so the panel works
unchanged for 2-, 3- or 4-backstress Chaboche (and future UVC / Ohno-Wang). A
``sigma_y0`` row shows the (possibly calibrated) yield stress, and the clipboard
export uses that calibrated value -- not the user's monotonic input.
"""

from __future__ import annotations

from tkinter import ttk
from typing import Optional, Sequence


class ParamPanel:
    """Two dynamic parameter tables (last / best) plus a clipboard export."""

    def __init__(self, parent) -> None:
        self.frame = ttk.LabelFrame(parent, text="Parameters")
        self.frame.columnconfigure(0, weight=1)
        self.frame.columnconfigure(1, weight=1)

        self._sub_last = ttk.LabelFrame(self.frame, text="Last run")
        self._sub_last.grid(row=0, column=0, sticky="nsew", padx=2, pady=2)
        self._sub_best = ttk.LabelFrame(self.frame, text="Best so far")
        self._sub_best.grid(row=0, column=1, sticky="nsew", padx=2, pady=2)

        self._btn_copy = ttk.Button(self.frame, text="Copy best -> ABAQUS block",
                                    command=self._copy_best, state="disabled")
        self._btn_copy.grid(row=1, column=0, columnspan=2, sticky="we",
                            padx=4, pady=4)

        self._model = None
        self._w_last: dict[str, ttk.Label] = {}
        self._w_best: dict[str, ttk.Label] = {}
        self._best_params: Optional[tuple[float, ...]] = None
        self._best_sy0: Optional[float] = None
        self._elastic = {"E": 200000.0, "nu": 0.3, "sy0": 500.0}

    # -- model wiring --------------------------------------------------------
    def set_model(self, model, *, E: float, nu: float, sy0: float) -> None:
        """Rebuild the tables for a model and record elastic constants."""
        self._model = model
        self._elastic = {"E": E, "nu": nu, "sy0": sy0}
        self._w_last = self._build_table(self._sub_last, model)
        self._w_best = self._build_table(self._sub_best, model)
        self._best_params = None
        self._best_sy0 = None
        self._btn_copy.configure(state="disabled")

    @staticmethod
    def _build_table(parent: ttk.LabelFrame, model) -> dict[str, ttk.Label]:
        for child in parent.winfo_children():
            child.destroy()
        widgets: dict[str, ttk.Label] = {}
        rows = list(model.param_names) + ["sigma_y0", "_obj"]
        units = list(model.param_units) + ["MPa", "MPa"]
        labels = list(model.param_names) + ["sigma_y0", "objective"]
        for i, (key, lab, unit) in enumerate(zip(rows, labels, units)):
            ttk.Label(parent, text=lab).grid(row=i, column=0, sticky="w", padx=4)
            val = ttk.Label(parent, text="-", width=12, anchor="e")
            val.grid(row=i, column=1, sticky="e", padx=4)
            ttk.Label(parent, text=unit, style="Muted.TLabel").grid(
                row=i, column=2, sticky="w", padx=4)
            widgets[key] = val
        return widgets

    # -- updates -------------------------------------------------------------
    def set_last(self, params: Sequence[float], objective: float,
                 sigma_y0: Optional[float] = None) -> None:
        self._fill(self._w_last, params, objective, sigma_y0)

    def set_best(self, params: Sequence[float], objective: float,
                 sigma_y0: Optional[float] = None) -> None:
        self._fill(self._w_best, params, objective, sigma_y0)
        self._best_params = tuple(float(p) for p in params)
        if sigma_y0 is not None:
            self._best_sy0 = float(sigma_y0)
        self._btn_copy.configure(state="normal")

    def _fill(self, widgets: dict[str, ttk.Label], params: Sequence[float],
              objective: float, sigma_y0: Optional[float]) -> None:
        if not widgets or self._model is None:
            return
        for name, val in zip(self._model.param_names, params):
            if name in widgets:
                widgets[name].configure(text=f"{val:.4g}")
        if "sigma_y0" in widgets:
            sy = sigma_y0 if sigma_y0 is not None else self._elastic["sy0"]
            widgets["sigma_y0"].configure(text=f"{sy:.4g}")
        widgets["_obj"].configure(text=f"{objective:.3f}")

    # -- clipboard -----------------------------------------------------------
    def _copy_best(self) -> None:
        """Emit the ABAQUS material block for the best params, using the
        CALIBRATED sigma_y0 (falls back to the user input only if sigma_y0 was
        not calibrated)."""
        if self._model is None or self._best_params is None:
            return
        sy0 = self._best_sy0 if self._best_sy0 is not None else self._elastic["sy0"]
        block = self._model.material_block(
            self._best_params, E=self._elastic["E"], nu=self._elastic["nu"],
            sy0=sy0, name="STEEL")
        self.frame.clipboard_clear()
        self.frame.clipboard_append(block)
