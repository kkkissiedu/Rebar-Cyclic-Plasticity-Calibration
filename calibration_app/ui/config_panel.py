"""
ui/config_panel.py
==================

All run configuration: data + geometry + yield, model / backend / optimiser
selection, search settings (sliders with live labels), a collapsible Advanced
section, and a dynamic parameter-bounds editor.

Grade-agnostic: yield stress is a per-run entry (typed or estimated from a
companion monotonic CSV), and every bound is user-editable so batch-variable
scrap rebar can be calibrated without any B500C assumption.

The panel owns only widgets/state; :meth:`get_config` and :meth:`get_bounds`
hand a plain dict to the pipeline. Data loading is delegated to the ``on_load``
callback supplied by the main window (which owns the data_loader).
"""

from __future__ import annotations

import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import Callable, Optional

from .widgets import Tooltip


def _slider(parent, row, label, lo, hi, init, *, integer=True, fmt="{:.0f}",
            on_change=None):
    """A labelled ttk.Scale with a live value readout. Returns the tk var."""
    ttk.Label(parent, text=label).grid(row=row, column=0, sticky="w", padx=4, pady=2)
    var = tk.IntVar(value=int(init)) if integer else tk.DoubleVar(value=float(init))
    val_lbl = ttk.Label(parent, text=fmt.format(init), width=8, anchor="e",
                        style="Accent.TLabel")
    val_lbl.grid(row=row, column=2, sticky="e", padx=4)

    def _on(_v):
        v = var.get()
        val_lbl.configure(text=fmt.format(int(v) if integer else v))
        if on_change:
            on_change()

    ttk.Scale(parent, from_=lo, to=hi, variable=var, command=_on).grid(
        row=row, column=1, sticky="we", padx=4)
    return var


class ConfigPanel:
    """Configuration controls; queried by the pipeline via get_config/get_bounds."""

    def __init__(self, parent, on_load: Optional[Callable[[str], None]] = None) -> None:
        self.frame = ttk.LabelFrame(parent, text="Configuration")
        self.frame.columnconfigure(1, weight=1)
        self._on_load = on_load
        self._bound_vars: dict[str, tuple[tk.DoubleVar, tk.DoubleVar]] = {}
        r = 0

        # -- data + geometry -------------------------------------------------
        ttk.Label(self.frame, text="Test data CSV:").grid(row=r, column=0, sticky="w", padx=4)
        self.var_data = tk.StringVar()
        ttk.Entry(self.frame, textvariable=self.var_data).grid(row=r, column=1, sticky="we", padx=4)
        ttk.Button(self.frame, text="Browse...", command=self._browse_data).grid(row=r, column=2, padx=4)
        r += 1

        self.lbl_detect = ttk.Label(self.frame, text="(no data loaded)", style="Muted.TLabel")
        self.lbl_detect.grid(row=r, column=0, columnspan=3, sticky="w", padx=4)
        r += 1

        geo = ttk.Frame(self.frame)
        geo.grid(row=r, column=0, columnspan=3, sticky="we", padx=2)
        self.var_diameter = tk.DoubleVar(value=12.0)
        self.var_gauge = tk.DoubleVar(value=60.0)
        self.var_strain = tk.DoubleVar(value=1.0)
        for i, (lab, var, w) in enumerate((
                ("Dia mm", self.var_diameter, 6),
                ("Gauge mm", self.var_gauge, 6),
                ("Strain %", self.var_strain, 6))):
            ttk.Label(geo, text=lab).grid(row=0, column=2 * i, sticky="e", padx=2)
            ttk.Entry(geo, textvariable=var, width=w).grid(row=0, column=2 * i + 1, padx=2)
        r += 1

        # -- FATIGUE_rebuilt.cae target for "Transfer to CAE" ----------------
        # Configurable (never hardcoded elsewhere); FATIGUE.cae is refused by
        # cae_transfer.build_transfer_config regardless of what is typed here.
        ttk.Label(self.frame, text="Rebuilt .cae file:").grid(row=r, column=0, sticky="w", padx=4)
        proj_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        self.var_cae_path = tk.StringVar(
            value=os.path.join(proj_root, "FATIGUE_rebuilt.cae"))
        cae_entry = ttk.Entry(self.frame, textvariable=self.var_cae_path)
        cae_entry.grid(row=r, column=1, sticky="we", padx=4)
        ttk.Button(self.frame, text="Browse...", command=self._browse_cae).grid(row=r, column=2, padx=4)
        Tooltip(cae_entry, "Target for 'Transfer to CAE'. Must be the rebuilt "
                           "model - FATIGUE.cae (the protected original) is "
                           "always refused.")
        r += 1

        # -- material description (free-text, grade-agnostic metadata) -------
        ttk.Label(self.frame, text="Material description:").grid(row=r, column=0, sticky="w", padx=4)
        self.var_material_desc = tk.StringVar(value="")
        md_entry = ttk.Entry(self.frame, textvariable=self.var_material_desc)
        md_entry.grid(row=r, column=1, columnspan=2, sticky="we", padx=4)
        Tooltip(md_entry, "Free-text note on this specimen/batch (e.g. 'Accra "
                          "scrap rebar lot 3'). Stored in session.json and shown "
                          "in the sessions browser, comparison CSV and export.")
        r += 1

        # -- yield stress (per-run, grade-agnostic) --------------------------
        ttk.Label(self.frame, text="Yield sigma_y0 (MPa):").grid(row=r, column=0, sticky="w", padx=4)
        self.var_sy0 = tk.DoubleVar(value=500.0)
        ttk.Entry(self.frame, textvariable=self.var_sy0, width=10).grid(row=r, column=1, sticky="w", padx=4)
        ttk.Button(self.frame, text="From monotonic CSV...",
                   command=self._browse_monotonic).grid(row=r, column=2, padx=4)
        r += 1

        yc = ttk.Frame(self.frame)
        yc.grid(row=r, column=0, columnspan=3, sticky="we", padx=2)
        self.var_calib_sy0 = tk.BooleanVar(value=True)
        chk_sy0 = ttk.Checkbutton(yc, text="Calibrate sigma_y0 (cyclic yield)",
                                  variable=self.var_calib_sy0)
        chk_sy0.grid(row=0, column=0, sticky="w", padx=2)
        Tooltip(chk_sy0, "The cyclic yield-surface size is smaller than the "
                         "monotonic yield (Bauschinger). Fit it in "
                         "[lower, sigma_y0]; sigma_y0 above is the upper limit.")
        ttk.Label(yc, text="lower:").grid(row=0, column=1, sticky="e", padx=2)
        self.var_sy0_lower = tk.DoubleVar(value=150.0)
        ttk.Entry(yc, textvariable=self.var_sy0_lower, width=8).grid(row=0, column=2, padx=2)
        r += 1

        # -- model + backstresses -------------------------------------------
        mb = ttk.Frame(self.frame)
        mb.grid(row=r, column=0, columnspan=3, sticky="we", padx=2, pady=4)
        ttk.Label(mb, text="Model:").grid(row=0, column=0, sticky="w", padx=2)
        self.var_model = tk.StringVar(value="chaboche")
        from ..core.surrogate import MODEL_REGISTRY
        ttk.Combobox(mb, textvariable=self.var_model, width=12, state="readonly",
                     values=list(MODEL_REGISTRY)).grid(row=0, column=1, padx=2)
        ttk.Label(mb, text="Backstresses:").grid(row=0, column=2, sticky="w", padx=2)
        self.var_nbs = tk.IntVar(value=3)
        sp_nbs = ttk.Spinbox(mb, from_=2, to=4, increment=1, width=4,
                             textvariable=self.var_nbs, state="readonly",
                             values=(2, 3, 4))
        sp_nbs.grid(row=0, column=3, padx=2)
        Tooltip(sp_nbs, "2: fast; 3: recommended for 1-6% range; "
                        "4: near-yield precision")
        self.var_run_all = tk.BooleanVar(value=False)
        ttk.Checkbutton(mb, text="Run all 3 models in sequence",
                        variable=self.var_run_all).grid(row=0, column=4, padx=6)
        r += 1

        # -- backend (radio) -------------------------------------------------
        bkf = ttk.Frame(self.frame)
        bkf.grid(row=r, column=0, columnspan=3, sticky="we", padx=2)
        ttk.Label(bkf, text="Backend:").grid(row=0, column=0, sticky="w", padx=2)
        self.var_backend = tk.StringVar(value="surrogate")
        ttk.Radiobutton(bkf, text="Surrogate", variable=self.var_backend,
                        value="surrogate", command=self._on_backend).grid(row=0, column=1, padx=2)
        ttk.Radiobutton(bkf, text="FE (ABAQUS)", variable=self.var_backend,
                        value="fe", command=self._on_backend).grid(row=0, column=2, padx=2)
        self.lbl_fe_warn = ttk.Label(bkf, text="", style="Warn.TLabel")
        self.lbl_fe_warn.grid(row=0, column=3, sticky="w", padx=6)
        r += 1

        # -- optimiser (radio) + per-option eval-count hints -----------------
        opf = ttk.Frame(self.frame)
        opf.grid(row=r, column=0, columnspan=3, sticky="we", padx=2)
        ttk.Label(opf, text="Optimiser:").grid(row=0, column=0, sticky="w", padx=2)
        self.var_optim = tk.StringVar(value="de")
        ttk.Radiobutton(opf, text="Sobol + DE", variable=self.var_optim,
                        value="de", command=self._update_hints).grid(row=0, column=1, sticky="w", padx=2)
        self.lbl_hint_de = ttk.Label(opf, text="", style="Muted.TLabel")
        self.lbl_hint_de.grid(row=0, column=2, sticky="w", padx=2)
        ttk.Radiobutton(opf, text="Sobol + Bayesian", variable=self.var_optim,
                        value="bayesian", command=self._update_hints).grid(row=1, column=1, sticky="w", padx=2)
        self.lbl_hint_bo = ttk.Label(opf, text="", style="Muted.TLabel")
        self.lbl_hint_bo.grid(row=1, column=2, sticky="w", padx=2)
        r += 1

        # -- primary controls (always visible) ------------------------------
        ttk.Label(self.frame, text="Calibration cycles N_fit").grid(row=r, column=0, sticky="w", padx=4, pady=2)
        self.var_nfit = tk.IntVar(value=5)
        self.sp_nfit = ttk.Spinbox(self.frame, from_=1, to=30, increment=1, width=8,
                                   textvariable=self.var_nfit)
        self.sp_nfit.grid(row=r, column=1, sticky="w", padx=4)
        self.lbl_nfit_max = ttk.Label(self.frame, text="(load data)", style="Muted.TLabel")
        self.lbl_nfit_max.grid(row=r, column=2, sticky="w", padx=4)
        r += 1

        self.var_sobol = _slider(self.frame, r, "Sobol samples (2^k)", 8, 15, 12,
                                 fmt="2^{:.0f}", on_change=self._update_hints); r += 1

        # -- collapsible Advanced -------------------------------------------
        self._adv_open = tk.BooleanVar(value=False)
        self._adv_btn = ttk.Button(self.frame, text="+ Advanced settings",
                                   command=self._toggle_adv)
        self._adv_btn.grid(row=r, column=0, columnspan=3, sticky="we", padx=4, pady=(6, 2))
        r += 1
        self._adv = ttk.Frame(self.frame)
        self._adv.grid(row=r, column=0, columnspan=3, sticky="we")
        self._adv.columnconfigure(1, weight=1)
        self._adv.grid_remove()
        r += 1

        a = 0
        self.var_de_maxiter = _slider(self._adv, a, "DE maxiter", 10, 1000, 300,
                                      on_change=self._update_hints); a += 1
        self.var_de_popsize = _slider(self._adv, a, "DE popsize", 5, 60, 20); a += 1
        self.var_trials = _slider(self._adv, a, "Bayesian trials", 20, 1000, 200,
                                  on_change=self._update_hints); a += 1
        self.var_topk = _slider(self._adv, a, "Sobol top-K seeds", 1, 30, 10); a += 1
        self.var_topm = _slider(self._adv, a, "ABAQUS top-M verify", 1, 10, 3); a += 1
        self.var_jobs = _slider(self._adv, a, "Concurrent ABAQUS jobs", 1, 8, 1); a += 1
        self.var_weight = tk.BooleanVar(value=True)
        ttk.Checkbutton(self._adv, text="Weight reversal regions (1.5x)",
                        variable=self.var_weight).grid(row=a, column=0, columnspan=2, sticky="w", padx=4); a += 1
        self.var_stage3 = tk.BooleanVar(value=False)
        ttk.Checkbutton(self._adv, text="Stage 3: ABAQUS verification of top-M",
                        variable=self.var_stage3).grid(row=a, column=0, columnspan=2, sticky="w", padx=4); a += 1

        self._bounds_frame = ttk.LabelFrame(self._adv, text="Parameter bounds")
        self._bounds_frame.grid(row=a, column=0, columnspan=3, sticky="we", padx=4, pady=4)

        self._on_backend()          # initialise hints + FE warning
        self._update_hints()

    # -- dynamic limits from loaded data ------------------------------------
    def set_cycle_limits(self, total_cycles: int) -> None:
        """Cap the N_fit spinbox at (total_cycles - 1); clamp the current value.

        Never hardcodes a ceiling — the max always comes from the loaded file.
        """
        max_fit = max(1, int(total_cycles) - 1)
        self.sp_nfit.configure(to=max_fit)
        if int(self.var_nfit.get()) > max_fit:
            self.var_nfit.set(max_fit)
        self.lbl_nfit_max.configure(text=f"max {max_fit}")

    # -- bounds --------------------------------------------------------------
    def rebuild_bounds(self, model) -> None:
        """Populate the bounds editor from a model's default bounds.

        The model's ``default_bounds`` are the literature envelope (the only
        place grade-informed numbers live). Each row carries a non-blocking
        yellow badge that appears when the user widens a bound past that
        envelope; an "Ack" button clears it (the value is kept either way).
        """
        for child in self._bounds_frame.winfo_children():
            child.destroy()
        self._bound_vars = {}
        self._bound_badges: dict[str, ttk.Label] = {}
        self._bound_ackbtns: dict[str, ttk.Button] = {}
        self._bound_ack: dict[str, bool] = {}
        ttk.Label(self._bounds_frame, text="param").grid(row=0, column=0, padx=2)
        ttk.Label(self._bounds_frame, text="lower").grid(row=0, column=1, padx=2)
        ttk.Label(self._bounds_frame, text="upper").grid(row=0, column=2, padx=2)
        defaults = model.default_bounds()
        self._env_defaults = dict(defaults)
        for i, name in enumerate(model.param_names, start=1):
            lo, hi = defaults[name]
            ttk.Label(self._bounds_frame, text=name).grid(row=i, column=0, sticky="w", padx=2)
            v_lo, v_hi = tk.DoubleVar(value=lo), tk.DoubleVar(value=hi)
            ttk.Entry(self._bounds_frame, textvariable=v_lo, width=10).grid(row=i, column=1, padx=2)
            ttk.Entry(self._bounds_frame, textvariable=v_hi, width=10).grid(row=i, column=2, padx=2)
            badge = ttk.Label(self._bounds_frame, text="[!] outside literature",
                              style="Badge.TLabel")
            badge.grid(row=i, column=3, padx=2)
            badge.grid_remove()
            ackbtn = ttk.Button(self._bounds_frame, text="Ack", width=5,
                                command=lambda n=name: self._ack_bound(n))
            ackbtn.grid(row=i, column=4, padx=2)
            ackbtn.grid_remove()
            self._bound_vars[name] = (v_lo, v_hi)
            self._bound_badges[name] = badge
            self._bound_ackbtns[name] = ackbtn
            # Editing a bound re-checks the envelope live.
            for v in (v_lo, v_hi):
                v.trace_add("write", lambda *_a, n=name: self._check_bound(n))
        self._check_bounds()

    def _check_bound(self, name: str) -> None:
        """Show/hide the yellow badge for one bound row vs the literature envelope."""
        if not hasattr(self, "_env_defaults") or name not in self._env_defaults:
            return
        badge = self._bound_badges.get(name)
        ackbtn = self._bound_ackbtns.get(name)
        if badge is None or ackbtn is None:
            return
        def_lo, def_hi = self._env_defaults[name]
        v_lo, v_hi = self._bound_vars[name]
        try:
            lo, hi = float(v_lo.get()), float(v_hi.get())
        except (tk.TclError, ValueError):
            return                                   # mid-edit partial value
        outside = (hi > def_hi) or (lo < def_lo)
        if outside and not self._bound_ack.get(name, False):
            badge.grid(); ackbtn.grid()
        else:
            badge.grid_remove(); ackbtn.grid_remove()

    def _check_bounds(self) -> None:
        """Re-check every bound row (used after a rebuild)."""
        for name in list(self._bound_vars):
            self._check_bound(name)

    def _ack_bound(self, name: str) -> None:
        """User explicitly acknowledges an out-of-literature bound; clear badge."""
        self._bound_ack[name] = True
        self._check_bound(name)

    def get_bounds(self) -> dict[str, tuple[float, float]]:
        return {n: (float(lo.get()), float(hi.get()))
                for n, (lo, hi) in self._bound_vars.items()}

    # -- config --------------------------------------------------------------
    def get_config(self) -> dict:
        """Return the full run configuration as a plain dict."""
        return {
            "data_file": self.var_data.get(),
            "material_description": self.var_material_desc.get().strip(),
            "diameter_mm": float(self.var_diameter.get()),
            "gauge_length_mm": float(self.var_gauge.get()),
            "strain_pct": float(self.var_strain.get()),
            "sigma_y0": float(self.var_sy0.get()),
            "calibrate_sy0": bool(self.var_calib_sy0.get()),
            "sy0_lower": float(self.var_sy0_lower.get()),
            "model": self.var_model.get(),
            "n_backstresses": int(self.var_nbs.get()),
            "run_all_models": bool(self.var_run_all.get()),
            "backend": self.var_backend.get(),
            "optimiser": self.var_optim.get(),
            "n_fit": int(self.var_nfit.get()),
            "sobol_pow2": int(self.var_sobol.get()),
            "de_maxiter": int(self.var_de_maxiter.get()),
            "de_popsize": int(self.var_de_popsize.get()),
            "n_trials": int(self.var_trials.get()),
            "top_k": int(self.var_topk.get()),
            "abaqus_top_m": int(self.var_topm.get()),
            "concurrent_jobs": int(self.var_jobs.get()),
            "weight_reversals": bool(self.var_weight.get()),
            "run_stage3": bool(self.var_stage3.get()),
            "cae_path": self.var_cae_path.get().strip(),
        }

    def apply_config(self, cfg: dict) -> None:
        """Restore fields from a saved session config (used by 'Load config')."""
        mapping = {
            "model": self.var_model, "n_backstresses": self.var_nbs,
            "backend": self.var_backend, "optimiser": self.var_optim,
            "n_fit": self.var_nfit, "sobol_pow2": self.var_sobol,
            "de_maxiter": self.var_de_maxiter, "de_popsize": self.var_de_popsize,
            "n_trials": self.var_trials, "sigma_y0": self.var_sy0,
            "run_all_models": self.var_run_all,
            "calibrate_sy0": self.var_calib_sy0, "sy0_lower": self.var_sy0_lower,
            "material_description": self.var_material_desc,
        }
        for key, var in mapping.items():
            if key in cfg and cfg[key] is not None:
                try:
                    var.set(cfg[key])
                except Exception:
                    pass
        self._on_backend(); self._update_hints()

    # -- display helpers -----------------------------------------------------
    def show_detected(self, *, strain_pct, gauge, diameter, sy0=None, warn=None) -> None:
        """Update the detected-parameters label and geometry fields."""
        if diameter:
            self.var_diameter.set(round(float(diameter), 3))
        if gauge:
            self.var_gauge.set(round(float(gauge), 3))
        if strain_pct is not None:
            self.var_strain.set(round(float(strain_pct), 3))
        if sy0 is not None:
            self.var_sy0.set(round(float(sy0), 1))
        msg = (f"Detected: strain {strain_pct}%, gauge {gauge} mm, dia {diameter} mm"
               if strain_pct is not None else
               "Filename does not match {N}data.csv - set strain % manually")
        if warn:
            msg += f"  [!] {warn}"
        self.lbl_detect.configure(text=msg,
                                  style="Warn.TLabel" if (warn or strain_pct is None) else "Muted.TLabel")

    # -- internal callbacks --------------------------------------------------
    def _update_hints(self) -> None:
        n = 2 ** int(self.var_sobol.get())
        self.lbl_hint_de.configure(text=f"{n} samples -> {int(self.var_de_maxiter.get())} DE iters")
        self.lbl_hint_bo.configure(text=f"{n} samples -> {int(self.var_trials.get())} acquisitions")

    def _on_backend(self) -> None:
        """Auto-flip optimiser and toggle the FE slow-warning on backend change."""
        if self.var_backend.get() == "fe":
            self.var_optim.set("bayesian")
            self.lbl_fe_warn.configure(text="slow: ~30s/eval")
        else:
            self.var_optim.set("de")
            self.lbl_fe_warn.configure(text="")
        self._update_hints()

    def _toggle_adv(self) -> None:
        if self._adv_open.get():
            self._adv.grid_remove(); self._adv_btn.configure(text="+ Advanced settings")
        else:
            self._adv.grid(); self._adv_btn.configure(text="- Advanced settings")
        self._adv_open.set(not self._adv_open.get())

    def _browse_cae(self) -> None:
        p = filedialog.askopenfilename(title="Select FATIGUE_rebuilt.cae",
                                       filetypes=[("ABAQUS CAE", "*.cae"), ("All", "*.*")])
        if p:
            if os.path.basename(p).lower() == "fatigue.cae":
                messagebox.showwarning(
                    "Protected file",
                    "FATIGUE.cae is the protected original and is never "
                    "modified. Select FATIGUE_rebuilt.cae instead.")
                return
            self.var_cae_path.set(p)

    def _browse_data(self) -> None:
        p = filedialog.askopenfilename(title="Select test data CSV",
                                       filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if p:
            self.var_data.set(p)
            if self._on_load:
                self._on_load(p)

    # -- monotonic yield dialog (auto-detect + mapping + thumbnail) ----------
    def _browse_monotonic(self) -> None:
        p = filedialog.askopenfilename(title="Select monotonic tension CSV",
                                       filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not p:
            return
        from ..core import data_loader
        cols = data_loader.monotonic_columns(p)
        self._monotonic_dialog(p, cols)

    def _monotonic_dialog(self, path: str, cols: dict) -> None:
        """Column-mapping + preview dialog for monotonic yield extraction."""
        from ..core import data_loader
        import matplotlib
        matplotlib.use("TkAgg")
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

        header = cols.get("header") or []
        win = tk.Toplevel(self.frame)
        win.title("Monotonic yield (0.2% offset)")
        win.geometry("560x520")

        top = ttk.Frame(win); top.pack(fill="x", padx=6, pady=6)
        ttk.Label(top, text="Force col:").grid(row=0, column=0, sticky="w")
        v_force = tk.StringVar(value=header[cols["force_idx"]] if cols.get("force_idx") is not None else "")
        ttk.Combobox(top, textvariable=v_force, values=header, width=22,
                     state="readonly").grid(row=0, column=1, padx=4)
        ttk.Label(top, text="Disp col:").grid(row=0, column=2, sticky="w")
        v_disp = tk.StringVar(value=header[cols["disp_idx"]] if cols.get("disp_idx") is not None else "")
        ttk.Combobox(top, textvariable=v_disp, values=header, width=22,
                     state="readonly").grid(row=0, column=3, padx=4)

        fig = Figure(figsize=(5, 3.4), dpi=100)
        ax = fig.add_subplot(111)
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=6)
        info = ttk.Label(win, text="", style="Muted.TLabel"); info.pack(pady=2)

        state = {"sy0": None}

        def recompute():
            try:
                fi = header.index(v_force.get()); di = header.index(v_disp.get())
            except ValueError:
                info.configure(text="Pick both columns.", style="Warn.TLabel"); return
            try:
                strain, stress = data_loader.monotonic_curve(
                    path, float(self.var_diameter.get()), float(self.var_gauge.get()), fi, di)
                sy0, E, offset = data_loader.yield_from_curve(strain, stress)
            except Exception as exc:
                info.configure(text=f"Read failed: {exc}", style="Warn.TLabel"); return
            state["sy0"] = sy0
            ax.clear()
            ax.plot(strain * 100.0, stress, color="tab:blue", lw=1.2, label="monotonic")
            if offset is not None:
                ox, oy = offset
                ax.plot(ox * 100.0, oy, color="0.5", ls="--", lw=0.9, label="0.2% offset")
            if sy0 is not None:
                ax.axhline(sy0, color="tab:red", lw=0.9)
                ax.text(0.02, sy0, f"sigma_y0 = {sy0:.1f} MPa", color="tab:red",
                        va="bottom", fontsize=8)
            ax.set_xlabel("Strain (%)"); ax.set_ylabel("Stress (MPa)")
            ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
            try:
                ax.set_ylim(0, max(stress.max() * 1.1, (sy0 or 1) * 1.3))
            except Exception:
                pass
            canvas.draw_idle()
            info.configure(
                text=(f"Detected sigma_y0 = {sy0:.1f} MPa (E~{E:.0f} MPa)"
                      if sy0 else "Could not find a 0.2% intercept - adjust columns"),
                style="Ok.TLabel" if sy0 else "Warn.TLabel")

        def accept():
            if state["sy0"]:
                self.var_sy0.set(round(state["sy0"], 1))
                self.lbl_detect.configure(
                    text=f"sigma_y0 = {state['sy0']:.1f} MPa from {os.path.basename(path)}",
                    style="Ok.TLabel")
                win.destroy()
            else:
                messagebox.showinfo("No yield", "No yield detected; adjust the columns.")

        btns = ttk.Frame(win); btns.pack(fill="x", pady=6)
        ttk.Button(btns, text="Recompute", command=recompute).pack(side="left", padx=6)
        ttk.Button(btns, text="Accept", command=accept).pack(side="right", padx=6)
        ttk.Button(btns, text="Cancel", command=win.destroy).pack(side="right")
        if v_force.get() and v_disp.get():
            recompute()
