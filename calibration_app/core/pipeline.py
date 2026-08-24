"""
core/pipeline.py
================

Run orchestration: ties a model + backend + optimiser + session together and
drives the three stages in a background thread, posting messages to a GUI queue.

Message contract (put on ``gui_queue`` as ``(kind, payload)``)::

    ("status",  {"text": str})
    ("model",   {"key": str, "session": path})            # a model run started
    ("stage",   {"key","status","done","total","elapsed","best"})
    ("eval",    {"params","obj","best_params","best_obj","improved",
                 "last_curve"|None, "best_curve"|None})    # curves = (strain%, stress)
    ("model_done", {"key","session","best_obj","best_params"})
    ("done",    {})

Live-plot curves are always rendered with the (fast) surrogate of the active
model - even when the FE backend is scoring - so plotting stays cheap and both
backends show progress the same way.

Backend selection:
    * surrogate -> objective.evaluate_surrogate
    * fe        -> FEBackend.evaluate, gated first so invalid candidates never
                   launch an ABAQUS job.
"""

from __future__ import annotations

import os
import threading
import time
import traceback

import numpy as np

from . import objective as objective_mod
from . import search as search_mod
from . import surrogate as surrogate_mod
from .abaqus_runner import FEBackend
from .data_loader import slice_window
from ..utils.physics_gate import DataContext, PhysicsGate

#: Cadence (in evaluations) for pushing a "last evaluated" curve to the plot.
LAST_CURVE_EVERY: int = 15


def build_model(model_key: str, n_backstresses: int = 3):
    """Instantiate a model by key (Chaboche / UVC / Ohno-Wang)."""
    if model_key == "chaboche":
        return surrogate_mod.ChabocheModel(n_backstresses=n_backstresses)
    if model_key == "uvc":
        from .uvc_model import UVCModel
        return UVCModel(n_backstresses=n_backstresses)
    if model_key == "ohno_wang":
        from .ohno_wang_model import OhnoWangModel
        return OhnoWangModel(n_backstresses=n_backstresses)
    raise NotImplementedError(f"unknown model '{model_key}'")


class CalibrationPipeline:
    """Executes one or more model calibrations, reporting to a GUI queue."""

    def __init__(self, gui_queue, stop_event: threading.Event, session_manager,
                 data: dict, config: dict, bounds: dict,
                 *, E: float, nu: float) -> None:
        self.q = gui_queue
        self.stop = stop_event
        self.sessions = session_manager
        self.data = data
        self.cfg = config
        self.bounds = bounds
        self.E, self.nu = E, nu
        self.sy0 = float(config["sigma_y0"])
        self.t0 = time.time()
        self._stage_t0 = self.t0        # reset at each stage so ETA is per-stage

    # -- messaging -----------------------------------------------------------
    def _send(self, kind: str, **payload) -> None:
        self.q.put((kind, payload))

    def _status(self, text: str) -> None:
        self._send("status", text=text)

    # -- entry ---------------------------------------------------------------
    def run(self) -> None:
        try:
            run_all = bool(self.cfg.get("run_all_models"))
            model_keys = (["chaboche", "uvc", "ohno_wang"] if run_all
                          else [self.cfg["model"]])
            parent = (self.sessions.new_comparison_dir(self.cfg["data_file"])
                      if run_all else None)
            summaries = []
            for key in model_keys:
                if self.stop.is_set():
                    break
                try:
                    model = build_model(key, self.cfg["n_backstresses"])
                except NotImplementedError as exc:
                    self._status(f"Skipping {key}: {exc}")
                    continue
                summary = self._run_one_model(key, model, subdir_of=parent)
                if summary:
                    summaries.append(summary)
            if run_all and summaries:
                csv_path = self._write_comparison_csv(parent, summaries)
                self._status(f"Model comparison written: {csv_path}")
                self._send("comparison_done", parent=parent, summaries=summaries)
            self._send("done")
        except Exception as exc:                       # never kill the GUI thread
            self._status(f"ERROR: {exc}\n{traceback.format_exc()}")
            self._send("done")

    # -- one model -----------------------------------------------------------
    def _run_one_model(self, key: str, model, subdir_of=None):
        cfg = self.cfg
        self._model_t0 = time.time()
        window = slice_window(self.data, cfg["n_fit"])
        display_window = slice_window(self.data, min(6, cfg["n_fit"]))

        # Grade-agnostic gate; Option C cross-amplitude peak (single file here).
        ctx = DataContext.from_datasets(self.data, sigma_y0=self.sy0)
        gate = PhysicsGate(ctx)

        # Per-model bounds: the user's config bounds where the parameter name
        # matches (the selected model), otherwise the model's own literature
        # defaults -- required so a run-all comparison uses each model's bounds.
        defaults = model.default_bounds()
        model_bounds = {n: self.bounds.get(n, defaults[n]) for n in model.param_names}

        session = self.sessions.new_session(
            cfg["data_file"], key, cfg["backend"], cfg["optimiser"], cfg,
            model_bounds, subdir_of=subdir_of)
        self._send("model", key=key, session=session.dir)
        self._status(f"[{key}] session {session.state['name']}")

        eval_counter = {"n": 0}
        amp_mm = self.data.get("amp_mm") or cfg["strain_pct"] / 100.0 * cfg["gauge_length_mm"]

        # Optionally calibrate sigma_y0 as the leading search variable. The
        # cyclic yield-surface size is smaller than the monotonic yield the user
        # enters (Bauschinger effect), so it is fit in [sy0_lower, sigma_y0].
        calib_sy0 = bool(cfg.get("calibrate_sy0"))
        sy0_fixed = self.sy0
        sy0_lower = float(cfg.get("sy0_lower", 150.0))
        chaboche_bounds = [model_bounds[n] for n in model.param_names]
        if calib_sy0:
            search_names = ["sigma_y0"] + list(model.param_names)
            search_bounds = [(sy0_lower, sy0_fixed)] + chaboche_bounds
        else:
            search_names = list(model.param_names)
            search_bounds = chaboche_bounds

        def split(x):
            """Split a search vector into (sigma_y0, chaboche_params)."""
            if calib_sy0:
                return float(x[0]), tuple(x[1:])
            return sy0_fixed, tuple(x)

        # -- evaluate callable (backend-specific) ---------------------------
        if cfg["backend"] == "fe":
            fe = FEBackend(
                model, window, work_dir=session.abaqus_dir,
                abaqus_cmd="abaqus", E=self.E, nu=self.nu, sy0=self.sy0,
                gauge_length_mm=cfg["gauge_length_mm"], amp_mm=amp_mm,
                n_cycles=len(window["target_cycles"]),
                concurrent_jobs=cfg["concurrent_jobs"],
                weight_reversals=cfg["weight_reversals"])

            def evaluate(x):
                sy0c, p = split(x)
                return fe.evaluate(gate.project(model, p, sy0c), sy0=sy0c)
        else:
            def evaluate(x):
                sy0c, p = split(x)
                return objective_mod.evaluate_surrogate(
                    model, p, window, gate, E=self.E, sy0=sy0c,
                    weight_reversals=cfg["weight_reversals"])

        # -- display curve helper (always surrogate, cheap) -----------------
        def curve(sy0c, params):
            burn = display_window.get("burnin_strain")
            strain = display_window["strain"]
            if burn is not None and np.size(burn) > 0:
                full = np.concatenate([np.asarray(burn), strain])
                sig = model.simulate(full, params, E=self.E, sy0=sy0c)[np.size(burn):]
            else:
                sig = model.simulate(strain, params, E=self.E, sy0=sy0c)
            return strain, sig

        # -- progress handler shared by both stages -------------------------
        def make_progress(stage_total):
            def _p(rec: search_mod.EvalRecord):
                # The objective simulated the PROJECTED params at the candidate's
                # sigma_y0, so report/record those, not the raw box point.
                sy0c, raw_p = split(rec.params)
                proj = tuple(gate.project(model, raw_p, sy0c))
                eval_counter["n"] += 1
                log_params = ([sy0c] + list(proj)) if calib_sy0 else proj
                session.record_eval(eval_counter["n"], rec.stage, rec.objective,
                                    cfg["backend"], rec.wall_s, log_params)
                improved = False
                best_curve = None
                if rec.objective < session.state.get("best_objective", float("inf")):
                    # Report the monotonic (input) vs cyclic (calibrated) yield
                    # and the Bauschinger reduction ratio -- a reportable finding
                    # explaining why monotonic yield != Chaboche sigma_y0.
                    desc = self.cfg.get("material_description")
                    header = ("** Material: %s\n" % desc) if desc else ""
                    if calib_sy0 and sy0_fixed:
                        ratio = sy0c / sy0_fixed
                        header += (
                            "** sigma_y0_monotonic: %.1f MPa  (grade spec / monotonic input)\n"
                            "** sigma_y0_cyclic:    %.1f MPa  (calibrated Chaboche yield-surface radius)\n"
                            "** Bauschinger_ratio:  %.3f  (cyclic / monotonic)\n"
                            % (sy0_fixed, sy0c, ratio))
                    block = header + model.material_block(
                        proj, E=self.E, nu=self.nu, sy0=sy0c, name="STEEL")
                    improved = session.update_best(proj, rec.objective,
                                                   cfg["backend"], block)
                    session.state["best_sigma_y0"] = float(sy0c)
                    session.state["sigma_y0_monotonic"] = float(sy0_fixed)
                    if calib_sy0 and sy0_fixed:
                        session.state["bauschinger_ratio"] = float(sy0c / sy0_fixed)
                    session.save()
                    try:
                        s, sig = curve(sy0c, proj)
                        if np.all(np.isfinite(sig)):
                            best_curve = (s * 100.0, sig)
                    except Exception:
                        best_curve = None
                last_curve = None
                if eval_counter["n"] % LAST_CURVE_EVERY == 0:
                    try:
                        s, sig = curve(sy0c, proj)
                        if np.all(np.isfinite(sig)):
                            last_curve = (s * 100.0, sig)
                    except Exception:
                        last_curve = None
                self._send("stage", key=rec.stage, status="RUNNING",
                           done=rec.index, total=stage_total,
                           elapsed=time.time() - self._stage_t0,
                           best=session.state.get("best_objective"))
                self._send("eval", params=list(proj), obj=rec.objective,
                           sigma_y0=float(sy0c),
                           best_params=session.state.get("best_params"),
                           best_obj=session.state.get("best_objective"),
                           best_sigma_y0=session.state.get("best_sigma_y0"),
                           improved=improved,
                           last_curve=last_curve, best_curve=best_curve)
            return _p

        # -- Stage 1: Sobol seeding -----------------------------------------
        session.set_stage("sobol")
        self._stage_t0 = time.time()
        n_samples = 2 ** int(cfg["sobol_pow2"])
        self._send("stage", key="sobol", status="RUNNING", done=0, total=n_samples,
                   elapsed=0.0)
        seeds, s1 = search_mod.sobol_seed(
            evaluate, search_bounds, n_samples=n_samples, top_k=cfg["top_k"],
            progress=make_progress(n_samples), stop=self.stop)
        self._send("stage", key="sobol",
                   status="DONE" if not s1.stopped else "FAILED",
                   done=s1.n_evals, total=n_samples,
                   elapsed=time.time() - self._stage_t0,
                   best=session.state.get("best_objective"))
        if self.stop.is_set():
            return self._finish_model(key, session, model)

        # -- Stage 2: DE or Bayesian ----------------------------------------
        session.set_stage(cfg["optimiser"])
        self._stage_t0 = time.time()
        x0 = seeds[0] if seeds else None
        if cfg["optimiser"] == "bayesian":
            total2 = int(cfg["n_trials"])
            self._send("stage", key="bayesian", status="RUNNING", done=0, total=total2)
            s2 = search_mod.run_bayesian(
                evaluate, search_bounds, search_names,
                n_trials=total2, x0=x0, progress=make_progress(total2), stop=self.stop)
        else:
            dim = len(search_bounds)
            total2 = int(cfg["de_maxiter"]) * int(cfg["de_popsize"]) * dim
            self._send("stage", key="de", status="RUNNING", done=0, total=total2)
            s2 = search_mod.run_differential_evolution(
                evaluate, search_bounds, x0=x0, maxiter=cfg["de_maxiter"],
                popsize=cfg["de_popsize"], progress=make_progress(total2), stop=self.stop)
        self._send("stage", key=cfg["optimiser"],
                   status="DONE" if not s2.stopped else "FAILED",
                   done=s2.n_evals, total=total2,
                   elapsed=time.time() - self._stage_t0,
                   best=session.state.get("best_objective"))

        # -- Stage 3: ABAQUS verification of the best (either backend) -------
        if self.stop.is_set():
            pass
        elif not cfg.get("run_stage3"):
            self._send("stage", key="abaqus", status="PENDING", done=0, total=0)
            self._status(f"[{key}] Stage 3 not enabled - skipped.")
        else:
            self._verify_top_m(key, model, session, window, amp_mm)

        return self._finish_model(key, session, model)

    # -- Stage 3 -------------------------------------------------------------
    def _verify_top_m(self, key, model, session, window, amp_mm) -> None:
        cfg = self.cfg
        best = session.state.get("best_params")
        if not best:
            self._status(f"[{key}] Stage 3 skipped - no valid best to verify.")
            self._send("stage", key="abaqus", status="FAILED", done=0, total=1)
            return
        self._stage_t0 = time.time()
        self._status(f"[{key}] Stage 3: verifying best params in ABAQUS "
                     f"(single element)...")
        self._send("stage", key="abaqus", status="RUNNING", done=0,
                   total=cfg["abaqus_top_m"])
        fe = FEBackend(
            model, window, work_dir=session.abaqus_dir, abaqus_cmd="abaqus",
            E=self.E, nu=self.nu, sy0=self.sy0,
            gauge_length_mm=cfg["gauge_length_mm"], amp_mm=amp_mm,
            n_cycles=len(window["target_cycles"]),
            concurrent_jobs=cfg["concurrent_jobs"],
            weight_reversals=cfg["weight_reversals"])
        sy0b = session.state.get("best_sigma_y0", self.sy0)
        # Run once and reuse the returned trace both for the objective and for
        # the plot overlay -- fe.evaluate() alone would re-run the same job.
        res = fe.run_candidate(best, sy0=sy0b)
        obj = (objective_mod.score_from_abaqus(res[0], res[1], window,
                                               cfg["weight_reversals"])
               if res is not None else objective_mod.PENALTY)
        if res is None or obj >= objective_mod.PENALTY:
            self._status(f"[{key}] Stage 3 ABAQUS run failed "
                         f"(no ODB / did not converge) - see session abaqus_runs/.")
            self._send("stage", key="abaqus", status="FAILED", done=0,
                       total=cfg["abaqus_top_m"], elapsed=time.time() - self._stage_t0)
            return
        strain, stress = res
        peak_t, peak_c = float(np.max(stress)), float(np.min(stress))
        session.state["fe_verification"] = {
            "objective_MPa": float(obj),
            "strain": strain.tolist(),
            "stress": stress.tolist(),
            "peak_tension_MPa": peak_t,
            "peak_compression_MPa": peak_c,
            "sigma_y0": float(sy0b),
        }
        session.save()
        self._status(f"[{key}] Stage 3 ABAQUS verification objective "
                     f"= {obj:.3f} MPa")
        self._send("stage", key="abaqus", status="DONE", done=1,
                   total=cfg["abaqus_top_m"], elapsed=time.time() - self._stage_t0)
        self._send("stage3_result", key=key, strain=strain.tolist(),
                   stress=stress.tolist(), objective=float(obj),
                   peak_tension=peak_t, peak_compression=peak_c,
                   sigma_y0=float(sy0b))

    def _finish_model(self, key, session, model) -> dict:
        """Finalise a model run: compute its best display curve, emit model_done
        (with the curve so the GUI can overlay all models), and return a summary
        row for the comparison CSV."""
        session.set_stage("done")
        best_params = session.state.get("best_params")
        best_sy0 = session.state.get("best_sigma_y0") or self.sy0
        best_curve = None
        peak_t = peak_c = None
        if best_params:
            try:
                dw = slice_window(self.data, min(6, self.cfg["n_fit"]))
                burn = dw.get("burnin_strain")
                strain = dw["strain"]
                if burn is not None and np.size(burn) > 0:
                    full = np.concatenate([np.asarray(burn), strain])
                    sig = model.simulate(full, best_params, E=self.E, sy0=best_sy0)[np.size(burn):]
                else:
                    sig = model.simulate(strain, best_params, E=self.E, sy0=best_sy0)
                if np.all(np.isfinite(sig)):
                    best_curve = (strain * 100.0, sig)
                    peak_t, peak_c = float(np.max(sig)), float(np.min(sig))
                    session.state["best_peak_tension"] = peak_t
                    session.state["best_peak_compression"] = peak_c
                    session.save()
            except Exception:
                best_curve = None
        self._send("model_done", key=key, session=session.dir,
                   best_obj=session.state.get("best_objective"),
                   best_params=best_params,
                   best_sigma_y0=session.state.get("best_sigma_y0"),
                   best_curve=best_curve,
                   peak_tension=peak_t, peak_compression=peak_c)
        return {
            "model": key,
            "material_description": self.cfg.get("material_description", ""),
            "best_objective_MPa": session.state.get("best_objective"),
            "n_evals": session.logger.count(),
            "wall_time_s": round(time.time() - self._model_t0, 1),
            "backend": self.cfg["backend"],
            "optimiser": self.cfg["optimiser"],
            "sy0_cyclic": session.state.get("best_sigma_y0"),
            "sy0_ratio": session.state.get("bauschinger_ratio"),
        }

    @staticmethod
    def _write_comparison_csv(parent_dir: str, rows: list) -> str:
        """Write the model-comparison summary CSV in the comparison parent dir."""
        import csv
        cols = ["model", "material_description", "best_objective_MPa", "n_evals",
                "wall_time_s", "backend", "optimiser", "sy0_cyclic", "sy0_ratio"]
        path = os.path.join(parent_dir, "model_comparison.csv")
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            for r in rows:
                w.writerow({c: r.get(c) for c in cols})
        return path
