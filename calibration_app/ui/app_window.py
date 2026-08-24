"""
ui/app_window.py
================

Main application window: lays out the plot, progress, parameter and
configuration panels; owns the theme toggle, the run controls, the live GUI
queue poller (plot throttled to <=5 Hz), and the Sessions browser.

The whole right-hand control column lives inside a :class:`ScrollableFrame` so
nothing is clipped on any display; the left plot expands with the window. The
window has a 900x600 minimum size.

Threading model: the calibration runs in a daemon thread
(:class:`~calibration_app.core.pipeline.CalibrationPipeline`) that posts
messages to a ``queue.Queue``; the Tk main thread drains it in :meth:`_poll`.
No Tk call is ever made off the main thread.
"""

from __future__ import annotations

import os
import queue
import threading
import time
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import Optional

import numpy as np

from . import theme as theme_mod
from .plot_panel import PlotPanel
from .progress_panel import ProgressPanel
from .param_panel import ParamPanel
from .config_panel import ConfigPanel
from .verification_panel import VerificationPanel
from .widgets import ScrollableFrame
from ..core import cae_transfer, data_loader, pipeline as pipeline_mod
from ..core.data_loader import slice_window
from ..utils import physics_gate
from ..session.session_manager import SessionManager

E_FIXED = 200_000.0
NU_FIXED = 0.30
PLOT_MIN_INTERVAL = 0.2          # seconds -> <=5 Hz plot refresh
DEFAULT_CYCLES = (2, 3, 4, 5, 6)


class AppWindow:
    """Owns the root window and wires the panels to the pipeline."""

    def __init__(self, root: tk.Tk) -> None:
        self.root = root
        self.mode = theme_mod.load_pref()
        self.style = ttk.Style(root)
        self.pal = theme_mod.apply_theme(root, self.style, self.mode)

        root.title("Cyclic-plasticity calibration - Chaboche / UVC / Ohno-Wang")
        root.geometry("1500x950")
        root.minsize(900, 600)

        self.gui_queue: queue.Queue = queue.Queue()
        self.stop_event = threading.Event()
        self.thread: Optional[threading.Thread] = None
        self.data: Optional[dict] = None
        self.sessions = SessionManager()
        self._last_plot_t = 0.0
        self.total_cycles = 0
        self._last_session_dir: Optional[str] = None   # newest completed run
        self._progress_win: Optional[tk.Toplevel] = None

        self._build_layout()
        self._refresh_model_panels()
        self._poll()

    # -- layout --------------------------------------------------------------
    def _build_layout(self) -> None:
        root = self.root
        root.columnconfigure(0, weight=3)
        root.columnconfigure(1, weight=2)
        root.rowconfigure(1, weight=1)

        # top bar
        top = ttk.Frame(root)
        top.grid(row=0, column=0, columnspan=2, sticky="we", padx=4, pady=4)
        ttk.Button(top, text="Start", style="Accent.TButton",
                   command=self._start).pack(side="left", padx=2)
        self.btn_stop = ttk.Button(top, text="Stop", command=self._stop, state="disabled")
        self.btn_stop.pack(side="left", padx=2)
        ttk.Button(top, text="Sessions...", command=self._open_sessions).pack(side="left", padx=2)
        ttk.Button(top, text="Export best...", command=self._export_best).pack(side="left", padx=2)
        self.btn_transfer = ttk.Button(top, text="Transfer to CAE",
                                       command=self._transfer_to_cae, state="disabled")
        self.btn_transfer.pack(side="left", padx=2)
        self.btn_fecomp = ttk.Button(top, text="Run post-run comparison",
                                     command=self._run_fe_comparison, state="disabled")
        self.btn_fecomp.pack(side="left", padx=2)
        ttk.Button(top, text="Reset axes", command=lambda: self.plot.reset_axes()).pack(side="left", padx=2)
        ttk.Button(top, text="Toggle theme", command=self._toggle_theme).pack(side="right", padx=2)
        ttk.Button(top, text="Apply", command=self._apply_cycles).pack(side="right", padx=2)
        self.var_cycles = tk.StringVar(value=",".join(str(c) for c in DEFAULT_CYCLES))
        self.ent_cycles = ttk.Entry(top, textvariable=self.var_cycles, width=18)
        self.ent_cycles.pack(side="right", padx=2)
        ttk.Label(top, text="Cycles:").pack(side="right")
        self.lbl_avail = ttk.Label(top, text="Available: -", style="Muted.TLabel")
        self.lbl_avail.pack(side="right", padx=6)

        # plot (left, spans both control rows) + verification table below it
        left_col = ttk.Frame(root)
        left_col.grid(row=1, column=0, sticky="nsew", padx=4, pady=4)
        left_col.columnconfigure(0, weight=1)
        left_col.rowconfigure(0, weight=1)

        plot_frame = ttk.LabelFrame(left_col, text="Stress-strain (experimental vs simulation)")
        plot_frame.grid(row=0, column=0, sticky="nsew")
        plot_frame.rowconfigure(0, weight=1); plot_frame.columnconfigure(0, weight=1)
        self.plot = PlotPanel(plot_frame, self.pal)
        self.plot.widget.grid(row=0, column=0, sticky="nsew")

        self.verification = VerificationPanel(
            left_col, on_view_change=lambda k: self.plot.set_view(k))
        self.verification.frame.grid(row=1, column=0, sticky="we", pady=(4, 0))

        # right: scrollable control column
        self.scroll = ScrollableFrame(root, bg=self.pal["panel"])
        self.scroll.outer.grid(row=1, column=1, sticky="nsew", padx=4, pady=4)
        inner = self.scroll.inner
        self.config = ConfigPanel(inner, on_load=self._load_data)
        self.config.frame.pack(fill="x", expand=False, pady=(0, 6))
        self.params = ParamPanel(inner)
        self.params.frame.pack(fill="x", expand=False, pady=(0, 6))
        self.progress = ProgressPanel(inner)
        self.progress.frame.pack(fill="x", expand=False)

        # status bar (A6)
        self.lbl_status = ttk.Label(root, text="No data loaded.", style="Muted.TLabel",
                                    anchor="w")
        self.lbl_status.grid(row=2, column=0, columnspan=2, sticky="we", padx=6)

        # log
        self.log = tk.Text(root, height=6, wrap="word",
                           bg=self.pal["field"], fg=self.pal["fg"], insertbackground=self.pal["fg"])
        self.log.grid(row=3, column=0, columnspan=2, sticky="we", padx=4, pady=4)

    def _log(self, text: str) -> None:
        self.log.insert("end", text + "\n"); self.log.see("end")

    def _current_model(self):
        cfg = self.config.get_config()
        return pipeline_mod.build_model("chaboche", cfg["n_backstresses"])

    def _refresh_model_panels(self) -> None:
        model = self._current_model()
        self.params.set_model(model, E=E_FIXED, nu=NU_FIXED,
                              sy0=float(self.config.var_sy0.get()))
        self.config.rebuild_bounds(model)

    # -- data ----------------------------------------------------------------
    def _load_data(self, path: str) -> None:
        try:
            geo = data_loader.detect_geometry(path)
            strain_pct = data_loader.detect_strain_amplitude_pct(path)
            gauge = geo["gauge_length_mm"] or float(self.config.var_gauge.get())
            dia = geo["diameter_mm"] or float(self.config.var_diameter.get())
            data = data_loader.load_experimental(
                path, gauge_length_mm=gauge, bar_diameter_mm=dia,
                sigma_y0=float(self.config.var_sy0.get()))
            self.data = data
            cs = data["cycle_summary"]
            self.total_cycles = cs["last_cycle"]
            warn = "filename not {N}data.csv" if strain_pct is None else None
            self.config.show_detected(strain_pct=strain_pct, gauge=gauge,
                                      diameter=dia, warn=warn)
            self.config.set_cycle_limits(cs["last_cycle"])
            self.lbl_avail.configure(
                text=f"Available: {cs['first_scoreable']}-{cs['last_cycle']}")

            # default cycles read from data, never hardcoded
            scoreable = cs["cycles_with_full_data"]
            defaults = [c for c in DEFAULT_CYCLES if c in scoreable] or scoreable[:5]
            self.var_cycles.set(",".join(str(c) for c in defaults))
            self.ent_cycles.configure(foreground=self.pal["fg"])

            for w in physics_gate.validate_startup(
                    float(self.config.var_sy0.get()), data["measured_peak_stress"]):
                self._log("[startup] " + w)
            self.plot.set_experimental(data, defaults)
            self._set_status(path, data, strain_pct, gauge, dia)
        except Exception as exc:
            messagebox.showerror("Load failed", str(exc))
            self._log(f"Load failed: {exc}")

    def _restore_session_curves(self, state: dict, cfg: dict) -> None:
        """Replay a loaded session's best-fit + ABAQUS FE curves onto the plot
        and verification panel so a completed run stays visible after 'Load
        this config' (Sessions browser), not just its numeric config."""
        model_key = state.get("model")
        best_params = state.get("best_params")
        if self.data is None or model_key is None or not best_params:
            return
        self.plot.clear_model_curves()
        self.verification.clear()
        self.verification.set_available_models([model_key])
        self.verification.set_view(model_key, notify=False)
        self.plot.set_view(model_key)
        try:
            model = pipeline_mod.build_model(model_key, cfg.get("n_backstresses", 3))
            best_sy0 = state.get("best_sigma_y0") or cfg.get("sigma_y0") or float(self.config.var_sy0.get())
            dw = slice_window(self.data, min(6, cfg.get("n_fit", 6)))
            burn = dw.get("burnin_strain")
            strain = dw["strain"]
            if burn is not None and np.size(burn) > 0:
                full = np.concatenate([np.asarray(burn), strain])
                sig = model.simulate(full, best_params, E=E_FIXED, sy0=best_sy0)[np.size(burn):]
            else:
                sig = model.simulate(strain, best_params, E=E_FIXED, sy0=best_sy0)
            if np.all(np.isfinite(sig)):
                # per-model curve (not generic update_best): generic best/last
                # are live-progress only and stay hidden while idle, so a
                # replayed session must use the view-filterable model curve.
                self.plot.add_model_curve(model_key, strain, sig)
                self.params.set_best(best_params, state.get("best_objective", 0.0), best_sy0)
        except Exception as exc:
            self._log(f"Could not replay session curve: {exc}")
            return
        self.verification.update_surrogate(
            model_key, objective=state.get("best_objective"),
            peak_tension=state.get("best_peak_tension"),
            peak_compression=state.get("best_peak_compression"),
            sigma_y0=state.get("best_sigma_y0"))
        fe = state.get("fe_verification")
        if fe:
            self.plot.add_fe_curve(
                model_key, np.asarray(fe["strain"]), np.asarray(fe["stress"]),
                fe["objective_MPa"], color=self.pal["fe"], label="ABAQUS (FE)")
            self.verification.update_fe(
                model_key, objective=fe.get("objective_MPa"),
                peak_tension=fe.get("peak_tension_MPa"),
                peak_compression=fe.get("peak_compression_MPa"))

    def _set_status(self, path, data, strain_pct, gauge, dia) -> None:
        cs = data["cycle_summary"]
        s = strain_pct if strain_pct is not None else "?"
        self.lbl_status.configure(
            text=(f"Loaded {os.path.basename(path)}: {cs['total_cycles']} cycles, "
                  f"peak {data['measured_peak_stress']:.1f} MPa, gauge {gauge:g} mm, "
                  f"dia {dia:g} mm, strain {s}%"))
        self._log(self.lbl_status.cget("text"))

    def _parse_cycles(self):
        """Return (valid_cycles, invalid_tokens) from the Cycles entry."""
        valid, invalid = [], []
        for tok in self.var_cycles.get().split(","):
            tok = tok.strip()
            if not tok:
                continue
            if tok.isdigit():
                c = int(tok)
                lo, hi = 2, (self.total_cycles or 10 ** 9)
                (valid if lo <= c <= hi else invalid).append(c if lo <= c <= hi else tok)
            else:
                invalid.append(tok)
        return valid, invalid

    def _apply_cycles(self) -> None:
        if self.data is None:
            return
        valid, invalid = self._parse_cycles()
        if invalid:
            self.ent_cycles.configure(foreground=self.pal["fail"])
            self._log(f"Invalid cycles (must be 2-{self.total_cycles}): {invalid}")
            return
        self.ent_cycles.configure(foreground=self.pal["fg"])
        self.plot.set_cycles(valid or list(DEFAULT_CYCLES))

    # -- run control ---------------------------------------------------------
    def _start(self) -> None:
        if self.thread and self.thread.is_alive():
            messagebox.showinfo("Running", "A calibration is already running.")
            return
        if self.data is None:
            path = self.config.var_data.get()
            if path:
                self._load_data(path)
            if self.data is None:
                messagebox.showwarning("No data", "Load a test CSV first.")
                return
        self._refresh_model_panels()
        cfg = self.config.get_config()
        bounds = self.config.get_bounds()
        self.progress.reset()
        self.plot.clear_model_curves()
        run_all = bool(cfg.get("run_all_models"))
        model_keys = ["chaboche", "uvc", "ohno_wang"] if run_all else [cfg["model"]]
        self.verification.clear()
        self.verification.set_available_models(model_keys)
        self.verification.set_view("all" if run_all else model_keys[0], notify=False)
        self.plot.set_view("all" if run_all else model_keys[0])
        self.stop_event.clear()
        pipe = pipeline_mod.CalibrationPipeline(
            self.gui_queue, self.stop_event, self.sessions, self.data, cfg,
            bounds, E=E_FIXED, nu=NU_FIXED)
        self.thread = threading.Thread(target=pipe.run, daemon=True)
        self.thread.start()
        self.btn_stop.configure(state="normal")
        self._log("Pipeline started.")

    def _stop(self) -> None:
        self.stop_event.set()
        self._log("Stop requested; halting after the current evaluation.")

    # -- queue poller --------------------------------------------------------
    def _poll(self) -> None:
        try:
            while True:
                kind, p = self.gui_queue.get_nowait()
                if kind == "status":
                    self._log(p["text"])
                elif kind == "model":
                    self._log(f"Model started: {p['key']}")
                    self.plot.set_running_model(p["key"])
                elif kind == "stage":
                    self.progress.update_stage(
                        p["key"], status=p["status"], done=p.get("done", 0),
                        total=p.get("total", 0), elapsed=p.get("elapsed", 0.0),
                        best=p.get("best"))
                elif kind == "eval":
                    self._on_eval(p)
                elif kind == "model_done":
                    self._last_session_dir = p.get("session")
                    if p.get("best_params"):
                        self.btn_transfer.configure(state="normal")
                    bo = p.get("best_obj"); sy = p.get("best_sigma_y0")
                    extra = f", sigma_y0={sy:.0f} MPa" if sy else ""
                    self._log(f"Model done: {p['key']} best {bo:.3f} MPa{extra}"
                              if bo else f"Model done: {p['key']}")
                    bc = p.get("best_curve")
                    if bc is not None:
                        s, sig = bc
                        self.plot.add_model_curve(p["key"], s / 100.0, sig)
                    self.verification.update_surrogate(
                        p["key"], objective=bo,
                        peak_tension=p.get("peak_tension"),
                        peak_compression=p.get("peak_compression"),
                        sigma_y0=sy)
                    # This model's solid curve now carries the "best fit"
                    # meaning; the generic best/last curves go idle (hidden
                    # in the all-models view) until the next model starts.
                    self.plot.set_running_model(None)
                    # Stages for this model are done -- the live "last
                    # evaluated" (blue) curve was only useful during the
                    # search; drop it so the graph shows only results that
                    # matter (best fit + FE verification).
                    self.plot.clear_last()
                elif kind == "stage3_result":
                    strain = np.asarray(p["strain"]); stress = np.asarray(p["stress"])
                    run_all = bool(self.config.get_config().get("run_all_models"))
                    if run_all:
                        color = self.plot.MODEL_COLORS.get(p["key"])
                        model_label = self.plot.MODEL_LABELS.get(p["key"], p["key"])
                        label = f"{model_label} ABAQUS (FE)"
                    else:
                        color = self.pal["fe"]
                        label = "ABAQUS (FE)"
                    self.plot.add_fe_curve(p["key"], strain, stress, p["objective"],
                                           color=color, label=label)
                    self.verification.update_fe(
                        p["key"], objective=p["objective"],
                        peak_tension=p.get("peak_tension"),
                        peak_compression=p.get("peak_compression"))
                elif kind == "comparison_done":
                    self._log(f"Comparison complete -> "
                              f"{os.path.join(p['parent'], 'model_comparison.csv')}")
                    for r in p.get("summaries", []):
                        bo = r.get("best_objective_MPa")
                        bo_s = f"{bo:.2f} MPa" if bo is not None else "-"
                        self._log(f"  {r['model']}: {bo_s}  "
                                  f"({r.get('n_evals')} evals, {r.get('wall_time_s')}s, "
                                  f"sy0 {r.get('sy0_cyclic')})")
                elif kind == "done":
                    self.btn_stop.configure(state="disabled")
                    self._log("Pipeline finished.")
                    self.plot.set_running_model(None)
                    self.plot.clear_last()
                elif kind == "transfer_done":
                    self._on_transfer_done(p)
                elif kind == "fecomp_done":
                    self._on_fecomp_done(p)
        except queue.Empty:
            pass
        self.root.after(100, self._poll)

    def _on_eval(self, p: dict) -> None:
        if p.get("params") is not None:
            self.params.set_last(p["params"], p["obj"], p.get("sigma_y0"))
        if p.get("best_params") is not None and p.get("best_obj") is not None:
            self.params.set_best(p["best_params"], p["best_obj"],
                                 p.get("best_sigma_y0"))
        # Gate ALL plot work on a single time check so a burst of queued evals
        # (the _poll drain can hand us many in one 100ms tick) triggers at most
        # one redraw per PLOT_MIN_INTERVAL. The params panel above is a cheap
        # label update and stays live; only the matplotlib redraw is throttled.
        now = time.time()
        if (now - self._last_plot_t) < PLOT_MIN_INTERVAL:
            return
        drew = False
        if p.get("best_curve") is not None:
            s, sig = p["best_curve"]
            self.plot.update_best(s / 100.0, sig)   # plot expects strain (frac)
            drew = True
        if p.get("last_curve") is not None:
            s, sig = p["last_curve"]
            self.plot.update_last(s / 100.0, sig)
            drew = True
        if drew:
            self._last_plot_t = now

    # -- Transfer to CAE + post-run FE comparison -----------------------------
    def _progress_dialog(self, title: str, text: str) -> None:
        """Small modal-ish window with an indeterminate bar while ABAQUS runs."""
        self._close_progress()
        win = tk.Toplevel(self.root)
        win.title(title)
        win.geometry("420x110")
        win.transient(self.root)
        win.protocol("WM_DELETE_WINDOW", lambda: None)   # closed by us only
        ttk.Label(win, text=text, wraplength=400).pack(padx=10, pady=(12, 6))
        bar = ttk.Progressbar(win, mode="indeterminate", length=380)
        bar.pack(padx=10, pady=6)
        bar.start(12)
        self._progress_win = win

    def _close_progress(self) -> None:
        if self._progress_win is not None:
            try:
                self._progress_win.destroy()
            except Exception:
                pass
            self._progress_win = None

    def _transfer_to_cae(self) -> None:
        """Write + run transfer_to_cae.py for the last completed session."""
        if not self._last_session_dir:
            messagebox.showinfo("No session", "Complete a calibration first.")
            return
        state = self.sessions.open(self._last_session_dir).state
        cae_path = self.config.get_config().get("cae_path", "")
        try:
            tcfg = cae_transfer.build_transfer_config(
                state, cae_path, session_dir=self._last_session_dir)
        except cae_transfer.TransferError as exc:
            messagebox.showerror("Transfer to CAE", str(exc))
            self._log(f"Transfer to CAE refused: {exc}")
            return
        script = cae_transfer.write_transfer_script(self._last_session_dir, tcfg)
        self._log(f"Transfer script written: {script}")
        self._progress_dialog("Transfer to CAE",
                              f"Running abaqus cae noGUI on\n"
                              f"{os.path.basename(tcfg['cae_path'])} ...")

        def worker():
            ok, out = cae_transfer.run_transfer_script(script)
            if ok and tcfg.get("material_mode") == "user":
                # UMAT models: prebuild standardU.dll in cae/ so the job
                # submits from CAE with no compiler (usub_lib_dir).
                lib_ok, lib_out = cae_transfer.build_umat_library(tcfg)
                ok = ok and lib_ok
                out += "\n--- UMAT library build ---\n" + lib_out
            self.gui_queue.put(("transfer_done",
                                {"ok": ok, "output": out, "tcfg": tcfg,
                                 "session_dir": self._last_session_dir}))

        threading.Thread(target=worker, daemon=True).start()

    def _on_transfer_done(self, p: dict) -> None:
        self._close_progress()
        out, tcfg = p["output"], p["tcfg"]
        if not p["ok"]:
            # Surface the exact CAE error, not just "failed".
            tail = "\n".join(out.strip().splitlines()[-15:])
            messagebox.showerror("Transfer to CAE failed", tail or "no output")
            self._log("Transfer to CAE FAILED:\n" + out)
            return
        log_path = cae_transfer.write_transfer_log(p["session_dir"], tcfg, out)
        self.btn_fecomp.configure(state="normal")
        summary = cae_transfer.transfer_summary(tcfg, out)
        msg = (summary + "\n\n"
               f"sigma_y0 (cyclic) = {tcfg['sigma_y0_cyclic']:.1f} MPa, "
               f"{tcfg['n_cycles']} cycles x {tcfg['period_s']:.0f} s, "
               f"amp {tcfg['amp_mm']:.3f} mm.\n"
               f"Log: {log_path}")
        stats = cae_transfer.parse_mesh_stats(out)
        if stats.get("n_failed"):
            messagebox.showwarning("Transfer to CAE - done (mesh quality)", msg)
        else:
            messagebox.showinfo("Transfer to CAE - done", msg)
        self._log("Transfer to CAE OK.\n" + summary)
        # Copy plot_fe_comparison.py's log next to nothing -- the log lives in
        # the session folder; remind the user of the post-run step instead.
        self._log("After the ABAQUS job completes, click 'Run post-run "
                  "comparison' and select the .odb.")

    def _run_fe_comparison(self) -> None:
        """Run plot_fe_comparison.py on a completed FATIGUE_rebuilt job ODB."""
        initial = None
        if self._last_session_dir:
            cae_dir = os.path.join(self._last_session_dir, "cae")
            if os.path.isdir(cae_dir):
                initial = cae_dir
        odb = filedialog.askopenfilename(
            title="Select completed job ODB", initialdir=initial,
            filetypes=[("ABAQUS output DB", "*.odb"), ("All", "*.*")])
        if not odb:
            return
        exp_csv = None
        if self._last_session_dir:
            exp_csv = self.sessions.open(self._last_session_dir).state.get("data_file")
        if not exp_csv or not os.path.exists(exp_csv):
            exp_csv = self.config.var_data.get()
        if not exp_csv or not os.path.exists(exp_csv):
            messagebox.showwarning("No experiment CSV",
                                   "Load the test data CSV first.")
            return
        # transfer_log.json lives in the session folder; copy it next to the
        # ODB so the standalone script finds the metadata (never overwrites a
        # log the user already placed there).
        if self._last_session_dir:
            src = os.path.join(self._last_session_dir, cae_transfer.TRANSFER_LOG)
            dst = os.path.join(os.path.dirname(odb), cae_transfer.TRANSFER_LOG)
            if os.path.exists(src) and not os.path.exists(dst):
                import shutil
                shutil.copyfile(src, dst)
        proj_root = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        script = os.path.join(proj_root, "plot_fe_comparison.py")
        self._progress_dialog("Post-run comparison",
                              f"Extracting {os.path.basename(odb)} via "
                              f"abaqus python ...")

        def worker():
            import subprocess
            cmd = ["abaqus", "python", script, "--", odb, exp_csv]
            try:
                proc = subprocess.run(cmd, cwd=os.path.dirname(odb),
                                      timeout=900, capture_output=True,
                                      text=True, shell=(os.name == "nt"))
                out = (proc.stdout or "") + "\n" + (proc.stderr or "")
                ok = "COMPARISON_OK" in out
            except Exception as exc:
                out, ok = repr(exc), False
            self.gui_queue.put(("fecomp_done",
                                {"ok": ok, "output": out,
                                 "odb_dir": os.path.dirname(odb)}))

        threading.Thread(target=worker, daemon=True).start()

    def _on_fecomp_done(self, p: dict) -> None:
        self._close_progress()
        if not p["ok"]:
            tail = "\n".join(p["output"].strip().splitlines()[-15:])
            messagebox.showerror("Post-run comparison failed", tail or "no output")
            self._log("Post-run comparison FAILED:\n" + p["output"])
            return
        odb_dir = p["odb_dir"]
        json_path = os.path.join(odb_dir, "fe_comparison.json")
        try:
            import json
            with open(json_path, "r") as f:
                summary = json.load(f)
        except Exception as exc:
            messagebox.showerror("Post-run comparison",
                                 f"Could not read fe_comparison.json: {exc}")
            return
        png = summary.get("png")
        if not png or not os.path.exists(png):
            # ABAQUS python had no matplotlib -- render from the CSV here.
            png = os.path.join(odb_dir, "fe_comparison.png")
            try:
                self._render_png_from_csv(summary.get("csv"), png, summary)
            except Exception as exc:
                self._log(f"Could not render fe_comparison.png: {exc}")
                png = None
        self.verification.update_fe(
            summary.get("model", "chaboche"),
            objective=summary.get("fe_objective_MPa"),
            peak_tension=summary.get("fe_peak_tension_MPa"),
            peak_compression=summary.get("fe_peak_compression_MPa"))
        self._log(f"Post-run FE comparison: objective "
                  f"{summary.get('fe_objective_MPa'):.1f} MPa "
                  f"(surrogate {summary.get('surrogate_objective_MPa')})")
        if png:
            self._show_png_window(png, "FE verification - full model")

    def _render_png_from_csv(self, csv_path: str, png_path: str,
                             summary: dict) -> None:
        """Fallback publication plot when ABAQUS python lacks matplotlib.
        Same content/format as plot_fe_comparison.try_png (DPI=150)."""
        import numpy as _np
        from matplotlib.figure import Figure
        arr = _np.genfromtxt(csv_path, delimiter=",", skip_header=1)
        fe_e, fe_s = arr[:, 1] * 100.0, arr[:, 2]
        exp_e, exp_s = arr[:, 3] * 100.0, arr[:, 4]
        fig = Figure(figsize=(7.0, 5.0), dpi=150)
        ax = fig.add_subplot(111)
        m = _np.isfinite(exp_e) & _np.isfinite(exp_s)
        ax.plot(exp_e[m], exp_s[m], color="0.6", lw=0.9, label="Experiment")
        ax.plot(fe_e, fe_s, color="tab:red", lw=1.2, label="FE simulation")
        dia = summary.get("diameter_mm") or 0
        gauge = summary.get("gauge_length_mm") or 0
        ld = gauge / dia if dia else 0
        strain_pct = summary.get("strain_pct") or 0
        ax.set_title(f"FE verification - {summary.get('model')} - "
                     f"{strain_pct:g}% strain - {dia:g}mm LD{ld:g}")
        sy0 = summary.get("sigma_y0")
        sy0_s = f"{sy0:.0f}" if sy0 else "?"
        ax.text(0.02, 0.02,
                f"Objective: {summary.get('fe_objective_MPa'):.1f} MPa | "
                f"Peak T: {summary.get('fe_peak_tension_MPa'):+.0f} MPa | "
                f"Peak C: {summary.get('fe_peak_compression_MPa'):+.0f} MPa | "
                f"sigma_y0: {sy0_s} MPa",
                transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        ax.set_xlabel("Strain (%)"); ax.set_ylabel("Stress (MPa)")
        ax.grid(True, alpha=0.3); ax.legend(fontsize=8, loc="upper left")
        fig.tight_layout()
        fig.savefig(png_path, dpi=150)

    def _show_png_window(self, png_path: str, title: str) -> None:
        """Display a PNG in a Toplevel inside the app (Tk 8.6 native PNG)."""
        win = tk.Toplevel(self.root)
        win.title(title)
        try:
            img = tk.PhotoImage(file=png_path)
        except Exception as exc:
            ttk.Label(win, text=f"Could not display {png_path}: {exc}").pack(
                padx=10, pady=10)
            return
        lbl = ttk.Label(win, image=img)
        lbl.image = img                      # keep a reference alive
        lbl.pack(padx=4, pady=4)
        ttk.Label(win, text=png_path, style="Muted.TLabel").pack(pady=(0, 4))

    # -- theme ---------------------------------------------------------------
    def _toggle_theme(self) -> None:
        self.mode = theme_mod.toggle(self.mode)
        self.pal = theme_mod.apply_theme(self.root, self.style, self.mode)
        self.plot.apply_palette(self.pal)
        self.scroll.set_bg(self.pal["panel"])
        self.log.configure(bg=self.pal["field"], fg=self.pal["fg"],
                           insertbackground=self.pal["fg"])

    # -- export --------------------------------------------------------------
    def _export_best(self) -> None:
        infos = self.sessions.list_sessions()
        if not infos:
            messagebox.showinfo("No sessions", "Run a calibration first.")
            return
        best = min(infos, key=lambda i: i.best_objective)
        src = os.path.join(best.path, "best_material_block.txt")
        if not os.path.exists(src):
            messagebox.showinfo("Nothing to export", "No best block yet.")
            return
        dst = filedialog.asksaveasfilename(
            defaultextension=".txt", initialfile="best_material_block.txt",
            filetypes=[("Text", "*.txt")])
        if dst:
            with open(src, "r") as f:
                block = f.read()
            with open(dst, "w") as f:
                f.write(block)
            self._log(f"Exported best block from {best.name} -> {dst}")

    # -- sessions browser ----------------------------------------------------
    def _open_sessions(self) -> None:
        win = tk.Toplevel(self.root)
        win.title("Sessions")
        win.geometry("1060x520")
        cols = ("created", "material", "data", "model", "backend", "optimiser", "stage", "best", "evals")
        tree = ttk.Treeview(win, columns=cols, show="headings")
        # Widths sized so headers/content are not truncated by default; ttk
        # Treeview columns remain drag-resizable (stretch keeps the layout filled
        # when the window is maximised).
        for c, w in zip(cols, (155, 150, 160, 95, 95, 95, 90, 100, 75)):
            tree.heading(c, text=c)
            tree.column(c, width=w, anchor="w", stretch=True)
        tree.pack(fill="both", expand=True, side="top")
        for i in self.sessions.list_sessions():
            bo = "-" if i.best_objective == float("inf") else f"{i.best_objective:.3f}"
            tree.insert("", "end", iid=i.path,
                        values=(i.created, i.material_description,
                                os.path.basename(i.data_file), i.model,
                                i.backend, i.optimiser, i.stage, bo, i.n_evals))
        bar = ttk.Frame(win); bar.pack(fill="x", side="bottom")

        def open_folder():
            sel = tree.selection()
            if sel:
                try:
                    os.startfile(sel[0])          # Windows: open the session dir
                except Exception as exc:
                    messagebox.showerror("Open folder", str(exc))

        def load_config():
            sel = tree.selection()
            if not sel:
                return
            st = self.sessions.open(sel[0]).state
            cfg = st.get("config", {})
            self.config.apply_config(cfg)
            if cfg.get("data_file"):
                self.config.var_data.set(cfg["data_file"])
                self._load_data(cfg["data_file"])
            self._refresh_model_panels()
            self._restore_session_curves(st, cfg)
            self._last_session_dir = sel[0]
            if st.get("best_params"):
                self.btn_transfer.configure(state="normal")
            win.destroy()
            self._log(f"Loaded config from {os.path.basename(sel[0])}. "
                      f"Prior sessions are preserved; Start makes a new one.")

        ttk.Button(bar, text="Open folder", command=open_folder).pack(side="left", padx=4, pady=4)
        ttk.Button(bar, text="Load this config", command=load_config).pack(side="left", padx=4, pady=4)
        ttk.Label(bar, text="Restarting always creates a NEW session; old runs are never overwritten.",
                  style="Muted.TLabel").pack(side="right", padx=6)
