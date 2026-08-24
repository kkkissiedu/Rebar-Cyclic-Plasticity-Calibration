"""
ui/plot_panel.py
================

Live stress-strain plot with a **fixed** vertical scale.

The old app let matplotlib autoscale, so a single blown-up candidate stretched
the axis and hid the experimental loops. Here the y-limits are locked to
``+/- Y_SCALE * peak_experimental_stress`` at initialisation and never change
during a run; simulated values that exceed the range are **clipped** rather than
triggering a rescale. A "Reset axes" button restores the locked view after any
manual zoom.

Curves:
* experimental cycles (user-selectable) - grey
* best-so-far - red, bold
* last evaluated - blue, thin
"""

from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

#: Vertical half-range as a multiple of the experimental peak stress.
Y_SCALE: float = 1.5


class PlotPanel:
    """Matplotlib stress-strain panel embedded in a Tk widget."""

    def __init__(self, parent, palette: dict[str, str]) -> None:
        self.pal = palette
        self.fig = Figure(figsize=(6, 4.5), dpi=110)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.widget = self.canvas.get_tk_widget()

        self._data: Optional[dict] = None
        self._shown_cycles: list[int] = []
        self._ylim: tuple[float, float] = (-700.0, 700.0)
        self._xlim: tuple[float, float] = (-1.0, 1.0)
        self._line_last = None
        self._line_best = None
        self._model_lines: dict = {}
        self._fe_lines: dict = {}
        self._fe_annots: dict = {}
        self._view = "all"
        #: Model key currently producing live best/last curves (None = idle).
        #: Only affects labelling/visibility while `_view == "all"` -- in a
        #: single-model view "best"/"last" are unambiguous and untouched.
        self._running_key: Optional[str] = None
        self._locking = False
        self._cids: list = []
        self._style_axes()
        self._connect_lock()

    def _connect_lock(self) -> None:
        """(Re)connect axis-limit callbacks that snap the view back if anything
        tries to rescale it. Reconnected after every ``ax.clear()`` because cla
        may drop the registrations."""
        for cid in self._cids:
            try:
                self.ax.callbacks.disconnect(cid)
            except Exception:
                pass
        self._cids = [self.ax.callbacks.connect("ylim_changed", self._snap),
                      self.ax.callbacks.connect("xlim_changed", self._snap)]

    def _snap(self, _ax) -> None:
        """Force the locked limits back; guarded against recursion."""
        if self._locking or self._data is None:
            return
        self._locking = True
        try:
            self.ax.set_ylim(*self._ylim)
            self.ax.set_xlim(*self._xlim)
        finally:
            self._locking = False

    # -- styling -------------------------------------------------------------
    def _style_axes(self) -> None:
        p = self.pal
        self.fig.patch.set_facecolor(p["panel"])
        self.ax.set_facecolor(p["panel"])
        for spine in self.ax.spines.values():
            spine.set_color(p["grid"])
        self.ax.tick_params(colors=p["fg"])
        self.ax.set_xlabel("Strain (%)", color=p["fg"])
        self.ax.set_ylabel("Stress (MPa)", color=p["fg"])
        self.ax.grid(True, color=p["grid"], alpha=0.4)

    def apply_palette(self, palette: dict[str, str]) -> None:
        """Re-theme the plot (called on dark/light toggle)."""
        self.pal = palette
        self._style_axes()
        if self._data is not None:
            self.set_experimental(self._data, self._shown_cycles)
        self.canvas.draw_idle()

    # -- experimental --------------------------------------------------------
    def set_experimental(self, data: dict, cycles: Optional[Sequence[int]] = None) -> None:
        """Draw the chosen experimental cycles and LOCK the axes.

        The y-scale is fixed from the *full* experimental peak (not just the
        shown cycles) so changing which cycles are displayed never rescales.
        """
        self._data = data
        peak = float(data.get("measured_peak_stress", 0.0)) or 1.0
        self._ylim = (-Y_SCALE * peak, Y_SCALE * peak)
        smax = float(np.max(np.abs(data["strain"]))) * 100.0 if data["strain"].size else 1.0
        self._xlim = (-1.15 * smax, 1.15 * smax)

        avail = data["cycles_present"]
        self._shown_cycles = list(cycles) if cycles else avail[:5]

        self.ax.clear()
        self._style_axes()
        self._connect_lock()          # cla drops callbacks; reconnect the lock
        for c in self._shown_cycles:
            idx = np.where(data["cycle"] == c)[0]
            if idx.size:
                self.ax.plot(data["strain"][idx] * 100.0, data["stress"][idx],
                             color=self.pal["exp"], linewidth=0.8, alpha=0.9)
        self.ax.set_title(f"Experimental cycles {self._shown_cycles} (grey)",
                          color=self.pal["fg"], fontsize=9)
        self._lock_axes()
        self._line_last = None
        self._line_best = None
        self._model_lines = {}
        self._fe_lines = {}
        self._fe_annots = {}
        self._running_key = None
        self.canvas.draw_idle()

    def _lock_axes(self) -> None:
        """Fix the limits and disable autoscaling."""
        self._locking = True
        try:
            self.ax.set_xlim(*self._xlim)
            self.ax.set_ylim(*self._ylim)
            self.ax.set_autoscale_on(False)
        finally:
            self._locking = False

    def reset_axes(self) -> None:
        """Restore the locked view (after a manual zoom/pan)."""
        self._lock_axes()
        self.canvas.draw_idle()

    def set_cycles(self, cycles: Sequence[int]) -> None:
        """Change which experimental cycles are displayed (keeps the y-scale)."""
        if self._data is not None:
            self.set_experimental(self._data, cycles)

    # -- simulated curves ----------------------------------------------------
    def _clip(self, stress: np.ndarray) -> np.ndarray:
        """Clip simulated stress into the locked y-range (no rescale)."""
        return np.clip(stress, self._ylim[0], self._ylim[1])

    def _running_label(self, plain: str, default_color: str) -> tuple[str, str]:
        """In the all-models view, while a model is actively running, label the
        live best/last curves with that model's name instead of the generic
        "best"/"last" (which is ambiguous once several per-model solid curves
        are also on the plot). Single-model views are unaffected."""
        if self._view == "all" and self._running_key:
            lab = self.MODEL_LABELS.get(self._running_key, self._running_key)
            return f"{lab} (running)", self.MODEL_COLORS.get(self._running_key, default_color)
        return plain, default_color

    def _generic_visible(self) -> bool:
        """Whether the generic best/last curves should be shown. They are
        LIVE-PROGRESS indicators: visible only while a model is actively
        running AND the current view includes that model. Once idle, the
        per-model curves from :meth:`add_model_curve` carry the best-fit
        meaning in every view. (2026-07-14 fix: they previously stayed
        visible in single-model views after the run, so switching to any
        per-model graph also showed the LAST run's generic fit on top of
        that model's best + FE curves.)"""
        return (self._running_key is not None
                and self._view in ("all", self._running_key))

    def set_running_model(self, key: Optional[str]) -> None:
        """Record which model is currently producing live best/last curves
        (None = no run active). Relabels the live curves per-model in the
        all-models view and re-derives their visibility (hidden whenever
        idle -- the per-model solid curves from :meth:`add_model_curve`
        already represent each model's best fit once its run finishes)."""
        self._running_key = key
        if self._line_best is not None or self._line_last is not None:
            show = self._generic_visible()
            if self._line_best is not None:
                self._line_best.set_visible(show)
            if self._line_last is not None:
                self._line_last.set_visible(show)
            self._rebuild_legend()
            self.canvas.draw_idle()

    def clear_last(self) -> None:
        """Remove the 'last evaluated' (blue) curve once a run finishes -- it
        is a live progress indicator during calibration; once stages are done
        only the best-fit / FE curves matter on the graph. Independent of
        `_view` (applies in single-model mode too, unlike the all-view-only
        best/last hide in :meth:`set_running_model`)."""
        if self._line_last is not None:
            try:
                self._line_last.remove()
            except Exception:
                pass
            self._line_last = None
            self._rebuild_legend()
            self.canvas.draw_idle()

    def update_last(self, strain: np.ndarray, stress: np.ndarray) -> None:
        """Update the 'last evaluated' curve (blue, thin)."""
        if self._line_last is not None:
            try:
                self._line_last.remove()
            except Exception:
                pass
        label, color = self._running_label("last", self.pal["last"])
        (self._line_last,) = self.ax.plot(
            np.asarray(strain) * 100.0, self._clip(np.asarray(stress)),
            color=color, linewidth=1.0, alpha=0.75, label=label)
        self._line_last.set_visible(self._generic_visible())
        self._rebuild_legend()
        self._lock_axes()
        self.canvas.draw_idle()

    def update_best(self, strain: np.ndarray, stress: np.ndarray) -> None:
        """Update the 'best so far' curve (red, bold)."""
        if self._line_best is not None:
            try:
                self._line_best.remove()
            except Exception:
                pass
        label, color = self._running_label("best", self.pal["best"])
        (self._line_best,) = self.ax.plot(
            np.asarray(strain) * 100.0, self._clip(np.asarray(stress)),
            color=color, linewidth=1.8, label=label)
        self._line_best.set_visible(self._generic_visible())
        self._rebuild_legend()
        self._lock_axes()
        self.canvas.draw_idle()

    # -- multi-model overlay (comparison mode) -------------------------------
    #: Distinct colours per model for the overlaid comparison curves; the same
    #: colour is reused (dashed) for that model's ABAQUS FE verification curve.
    MODEL_COLORS = {"chaboche": "#ff5555", "uvc": "#4ea1ff", "ohno_wang": "#ff9800"}
    MODEL_LABELS = {"chaboche": "Chaboche", "uvc": "UVC", "ohno_wang": "Ohno-Wang"}

    def clear_model_curves(self) -> None:
        """Remove all overlaid per-model best + FE curves (called at run start)."""
        for line in self._model_lines.values():
            try:
                line.remove()
            except Exception:
                pass
        self._model_lines = {}
        self._clear_fe_curves()
        self._running_key = None
        self.canvas.draw_idle()

    def add_model_curve(self, key: str, strain: np.ndarray, stress: np.ndarray) -> None:
        """Overlay one model's best-fit curve in its distinct colour (solid, bold)."""
        if key in self._model_lines:
            try:
                self._model_lines[key].remove()
            except Exception:
                pass
        col = self.MODEL_COLORS.get(key, self.pal["fg"])
        label = self.MODEL_LABELS.get(key, key)
        (self._model_lines[key],) = self.ax.plot(
            np.asarray(strain) * 100.0, self._clip(np.asarray(stress)),
            color=col, linewidth=1.8, label=label)
        visible = self._view in ("all", key)
        self._model_lines[key].set_visible(visible)
        self._rebuild_legend()
        self._lock_axes()
        self.canvas.draw_idle()

    # -- ABAQUS FE verification overlay --------------------------------------
    def _clear_fe_curves(self) -> None:
        for line in self._fe_lines.values():
            try:
                line.remove()
            except Exception:
                pass
        for ann in self._fe_annots.values():
            try:
                ann.remove()
            except Exception:
                pass
        self._fe_lines = {}
        self._fe_annots = {}

    def add_fe_curve(self, key: str, strain: np.ndarray, stress: np.ndarray,
                     objective: float, *, color: Optional[str] = None,
                     label: Optional[str] = None) -> None:
        """Overlay one model's ABAQUS verification curve (dashed).

        ``key`` identifies the curve for later removal/visibility toggling -
        pass the model key ("chaboche"/"uvc"/"ohno_wang") in comparison mode so
        it shares colour with :meth:`add_model_curve`, or any single key (e.g.
        the active model) for a non-comparison run, with ``color`` defaulted to
        the palette's dedicated "fe" green.
        """
        if key in self._fe_lines:
            try:
                self._fe_lines[key].remove()
            except Exception:
                pass
        if key in self._fe_annots:
            try:
                self._fe_annots[key].remove()
            except Exception:
                pass
        col = color or self.MODEL_COLORS.get(key) or self.pal["fe"]
        lab = label or "ABAQUS (FE)"
        strain_pct = np.asarray(strain) * 100.0
        stress_c = self._clip(np.asarray(stress))
        (self._fe_lines[key],) = self.ax.plot(
            strain_pct, stress_c, color=col, linewidth=1.6, linestyle="--",
            label=lab)
        if strain_pct.size:
            i = int(np.argmax(np.abs(stress_c)))
            self._fe_annots[key] = self.ax.annotate(
                f"FE: {objective:.1f} MPa", xy=(strain_pct[i], stress_c[i]),
                xytext=(6, 6), textcoords="offset points", fontsize=7, color=col)
        visible = self._view in ("all", key)
        self._fe_lines[key].set_visible(visible)
        if key in self._fe_annots:
            self._fe_annots[key].set_visible(visible)
        self._rebuild_legend()
        self._lock_axes()
        self.canvas.draw_idle()

    # -- view selector (all models / single model) ---------------------------
    def set_view(self, key: str) -> None:
        """Show only ``key``'s curves ("all" shows every model); curves for
        hidden models stay in memory and reappear when re-selected."""
        self._view = key
        for k, line in self._model_lines.items():
            line.set_visible(self._view in ("all", k))
        for k, line in self._fe_lines.items():
            line.set_visible(self._view in ("all", k))
        for k, ann in self._fe_annots.items():
            ann.set_visible(self._view in ("all", k))
        # Re-derive the generic best/last curves' visibility + label/colour
        # for the new view (live-progress only; hidden whenever idle or the
        # view excludes the running model -- see _generic_visible).
        show_generic = self._generic_visible()
        if self._line_best is not None:
            lab, col = self._running_label("best", self.pal["best"])
            self._line_best.set_label(lab)
            self._line_best.set_color(col)
            self._line_best.set_visible(show_generic)
        if self._line_last is not None:
            lab, col = self._running_label("last", self.pal["last"])
            self._line_last.set_label(lab)
            self._line_last.set_color(col)
            self._line_last.set_visible(show_generic)
        self._rebuild_legend()
        self.canvas.draw_idle()

    def _rebuild_legend(self) -> None:
        handles, labels = self.ax.get_legend_handles_labels()
        pairs = [(h, l) for h, l in zip(handles, labels) if h.get_visible()]
        if pairs:
            self.ax.legend([h for h, _ in pairs], [l for _, l in pairs],
                           loc="lower right", fontsize=8)
        elif self.ax.get_legend() is not None:
            self.ax.get_legend().remove()
