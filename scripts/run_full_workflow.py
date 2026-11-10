"""
run_full_workflow.py
=====================

End-to-end, unattended pipeline: calibrate -> Transfer to CAE -> submit the
full-coupon ABAQUS job -> wait for it -> run the post-run comparison. One
call produces one complete, documentation-ready session folder with no
manual CAE/GUI step in between -- built for generating a batch of results
for the website (see run_batch.py for driving several cases).

Works for both cyclic (Kashani-layout) and monotonic test-data CSVs; the
loading mode is auto-detected from the file (see
core/data_loader.detect_loading_mode) and threaded through calibration,
FE verification and the full-coupon deck automatically.

Usage (cuda_pt python, from the project root):

    python run_full_workflow.py --data input_data/16mm dia/LD5/2data.csv \\
        --model chaboche --sy0 500 [--n-backstresses 3] [--n-fit 5] \\
        [--optimiser bayesian] [--tag "batch 1"] \\
        [--cae-path FATIGUE_rebuilt.cae] [--wait-timeout 14400] \\
        [--skip-transfer] [--manifest calibration_sessions/results.jsonl]

Each stage's status is printed; on success a JSON record (session dir,
model, loading mode, surrogate + FE objectives, peaks, figure paths) is
appended to the manifest file (one JSON object per line) for the batch
driver / website to pick up.
"""

from __future__ import annotations

import argparse
import json
import os
import queue
import subprocess
import sys
import threading
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from calibration_app.core import cae_transfer, data_loader
from calibration_app.core import pipeline as pipeline_mod
from calibration_app.session.session_manager import SessionManager
from run_calibration_headless import drain

E_FIXED = 200_000.0
NU_FIXED = 0.30
DEFAULT_MANIFEST = os.path.join("calibration_sessions", "website_results.jsonl")


def log(msg: str) -> None:
    print(f"[full-workflow] {msg}")


def build_cfg(args, data: dict) -> dict:
    """Config dict matching config_panel.get_config()'s keys/defaults, built
    from scratch (no GUI, no prior session) for a headless run."""
    max_fit = max(1, data["cycle_summary"]["total_cycles"] - 1)
    n_fit = min(int(args.n_fit or 5), max_fit)
    return {
        "data_file": args.data,
        "material_description": args.tag or "",
        "diameter_mm": data["bar_diameter_mm"],
        "gauge_length_mm": data["gauge_length_mm"],
        "strain_pct": data.get("strain_amplitude_pct") or 0.0,
        "sigma_y0": float(args.sy0),
        "calibrate_sy0": not args.no_calibrate_sy0,
        "sy0_lower": float(args.sy0_lower),
        "model": args.model,
        "n_backstresses": int(args.n_backstresses),
        "run_all_models": False,
        "backend": args.backend,
        "optimiser": args.optimiser,
        "n_fit": n_fit,
        "sobol_pow2": int(args.sobol_pow2),
        "de_maxiter": 300,
        "de_popsize": 20,
        "n_trials": int(args.n_trials),
        "top_k": 10,
        "abaqus_top_m": 3,
        "concurrent_jobs": 1,
        "weight_reversals": True,
        "run_stage3": not args.no_stage3,
        "cae_path": args.cae_path,
        "loading_mode": data.get("loading_mode", "cyclic"),
        "E_modulus": float(args.E_resolved),
    }


class _TeeQueue(queue.Queue):
    """Queue that remembers the session dir announced by ``model_done``.

    Picking "the newest session" afterwards is racy when two workflows start
    in the same second (both get ``YYYYMMDD_HHMMSS`` names); the pipeline's
    own message is unambiguous.
    """
    session = None

    def put(self, item, block=True, timeout=None):
        try:
            if item[0] == "model_done":
                self.session = item[1].get("session")
        except Exception:
            pass
        super().put(item, block, timeout)


def run_calibration(cfg: dict, data: dict, bounds: dict = None) -> str:
    """Run the pipeline synchronously; return the new session directory."""
    q = _TeeQueue()
    stop = threading.Event()
    done = threading.Event()
    printer = threading.Thread(target=drain, args=(q, done), daemon=True)
    printer.start()
    pipe = pipeline_mod.CalibrationPipeline(
        q, stop, SessionManager(), data, cfg, bounds=dict(bounds or {}),
        E=float(cfg.get("E_modulus", E_FIXED)), nu=NU_FIXED)
    pipe.run()
    done.set()
    printer.join(timeout=5)
    if q.session:
        return q.session
    return SessionManager().list_sessions()[0].path


def run_transfer_and_solve(session_dir: str, cae_path: str,
                           wait_timeout: int) -> tuple[bool, dict]:
    """Transfer -> (build UMAT lib if needed) -> auto-submit -> wait.

    Returns (ok, transfer_log_dict). Raises nothing -- failures are
    reported via the returned flag + printed diagnostics so a batch driver
    can move on to the next case.
    """
    sm = SessionManager()
    state = sm.open(session_dir).state
    if not state.get("best_params"):
        log("no calibrated best_params -- calibration must have failed")
        return False, {}

    tcfg = cae_transfer.build_transfer_config(
        state, cae_path, session_dir=session_dir, auto_submit=True)
    log(f"transfer target: {os.path.basename(tcfg['cae_path'])} "
        f"(loading_mode={tcfg['loading_mode']}, material_mode={tcfg['material_mode']})")

    # UMAT models: the standardU.dll MUST exist BEFORE the job is submitted.
    # With auto_submit the submit+solve happens INSIDE the transfer script
    # (write_transfer_script/run_transfer_script below), so the library has
    # to be built first -- building it after, as the manual (non-auto)
    # workflow always did, was a real bug: the solve ran with no UMAT and
    # failed at increment 1 ("USER SUBROUTINE UMAT IS MISSING") while still
    # producing an ODB file, which looked like a pass to a naive existence
    # check. build_transfer_config already created cae/ so this is safe to
    # do before the .cae copy itself exists.
    if tcfg.get("material_mode") == "user":
        lib_ok, lib_out = cae_transfer.build_umat_library(tcfg)
        log(f"UMAT library build: {'OK' if lib_ok else 'FAILED'}")
        if not lib_ok:
            log(lib_out[-2000:])
            return False, {}

    script = cae_transfer.write_transfer_script(session_dir, tcfg)
    log(f"submitting + waiting for the ABAQUS solve "
        f"(timeout {wait_timeout/3600:.1f} h) ...")
    t0 = time.time()
    ok, out = cae_transfer.run_transfer_script(script, timeout_s=wait_timeout)
    solve_ok = "SOLVE_OK" in out
    log(f"transfer edit {'OK' if ok else 'FAILED'}, "
        f"solve {'OK' if solve_ok else 'FAILED'} "
        f"({time.time() - t0:.0f} s)")
    if not ok or not solve_ok:
        tail = "\n".join(out.strip().splitlines()[-30:])
        log("--- last 30 lines of transfer/solve log ---\n" + tail)
        return False, {}

    log_path = cae_transfer.write_transfer_log(session_dir, tcfg, out)
    with open(log_path, "r") as f:
        return True, json.load(f)


def run_comparison(cae_dir: str, exp_csv: str, abaqus_cmd: str = "abaqus") -> dict:
    """Run plot_fe_comparison.py against the just-solved ODB; return its
    fe_comparison.json contents, or {} on failure. ``exp_csv`` is the
    original test-data CSV (plot_fe_comparison.py requires it explicitly --
    it only reads transfer_log.json for everything else)."""
    odb = os.path.join(cae_dir, "Calibration_Job.odb")
    if not os.path.exists(odb):
        log(f"no ODB at {odb} -- cannot compare")
        return {}
    proj_root = os.path.dirname(os.path.abspath(__file__))
    script = os.path.join(proj_root, "plot_fe_comparison.py")
    # exp_csv is commonly relative to the project root (CLI --data), but
    # this subprocess runs with cwd=cae_dir -- resolve before that breaks it.
    exp_csv_abs = os.path.abspath(os.path.join(proj_root, exp_csv))
    cmd = [abaqus_cmd, "python", script, "--", odb, exp_csv_abs]
    proc = subprocess.run(cmd, cwd=cae_dir, capture_output=True, text=True,
                          timeout=600, shell=(os.name == "nt"))
    log(f"comparison: {'OK' if 'COMPARISON_OK' in (proc.stdout or '') else 'FAILED'}")
    result_json = os.path.join(cae_dir, "fe_comparison.json")
    if os.path.exists(result_json):
        with open(result_json, "r") as f:
            return json.load(f)
    log("--- plot_fe_comparison.py output ---\n" + (proc.stdout or "")[-2000:]
        + "\n" + (proc.stderr or "")[-1000:])
    return {}


def append_manifest(path: str, record: dict) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(record) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="test-data CSV (cyclic or monotonic)")
    ap.add_argument("--model", required=True, choices=["chaboche", "uvc", "ohno_wang"])
    ap.add_argument("--sy0", required=True, type=float, help="yield stress (MPa), per-run input")
    ap.add_argument("--n-backstresses", type=int, default=3)
    ap.add_argument("--n-fit", type=int, default=None,
                    help="cycles to fit (cyclic only; monotonic always uses 1)")
    ap.add_argument("--sy0-lower", type=float, default=150.0)
    ap.add_argument("--no-calibrate-sy0", action="store_true")
    ap.add_argument("--backend", default="surrogate", choices=["surrogate", "fe"])
    ap.add_argument("--optimiser", default="bayesian", choices=["bayesian", "de"])
    ap.add_argument("--sobol-pow2", type=int, default=12)
    ap.add_argument("--n-trials", type=int, default=200)
    ap.add_argument("--no-stage3", action="store_true",
                    help="skip the single-element FE sanity check before transfer")
    ap.add_argument("--tag", default="")
    ap.add_argument("--cae-path", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "FATIGUE_rebuilt.cae"))
    ap.add_argument("--wait-timeout", type=int, default=4 * 3600,
                    help="seconds to wait for the full-coupon solve (default 4h)")
    ap.add_argument("--skip-transfer", action="store_true",
                    help="calibrate only, no Transfer to CAE / ABAQUS submit")
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST)
    ap.add_argument("--E", type=float, default=E_FIXED,
                    help="elastic modulus (MPa) used by surrogate, FE and the "
                         "CAE transfer (default 200000)")
    ap.add_argument("--E-from-data", action="store_true",
                    help="use the apparent modulus measured from the data's "
                         "own unloading slope instead of --E (cyclic only)")
    args = ap.parse_args()

    data = data_loader.load_any_experimental(
        args.data, gauge_length_mm=data_loader.detect_geometry(args.data)["gauge_length_mm"],
        bar_diameter_mm=data_loader.detect_geometry(args.data)["diameter_mm"],
        sigma_y0=args.sy0)
    log(f"loaded {os.path.basename(args.data)}: loading_mode={data['loading_mode']}, "
        f"peak {data['measured_peak_stress']:.1f} MPa")

    args.E_resolved = args.E
    if args.E_from_data and data.get("loading_mode") != "monotonic":
        e_app = data_loader.apparent_modulus(data)
        if e_app:
            args.E_resolved = e_app
    log(f"elastic modulus: {args.E_resolved:.0f} MPa"
        f"{' (apparent, from data)' if args.E_resolved != args.E else ''}")

    cfg = build_cfg(args, data)
    log(f"calibrating: model={cfg['model']} N={cfg['n_backstresses']} "
        f"backend={cfg['backend']} optimiser={cfg['optimiser']} n_fit={cfg['n_fit']}")
    session_dir = run_calibration(cfg, data)
    log(f"calibration done: {session_dir}")

    record = {
        "session_dir": session_dir,
        "data_file": args.data,
        "model": args.model,
        "loading_mode": data["loading_mode"],
        "tag": args.tag,
    }
    sm = SessionManager()
    st = sm.open(session_dir).state
    record["surrogate_objective_MPa"] = st.get("best_objective")
    record["sigma_y0_cyclic"] = st.get("best_sigma_y0")

    if args.skip_transfer:
        append_manifest(args.manifest, record)
        log(f"skip-transfer requested; manifest updated: {args.manifest}")
        return 0

    ok, tlog = run_transfer_and_solve(session_dir, args.cae_path, args.wait_timeout)
    record["transfer_ok"] = ok
    if not ok:
        append_manifest(args.manifest, record)
        log("FAILED at transfer/solve stage; manifest updated with partial record")
        return 1

    cae_dir = os.path.dirname(tlog["cae_path"])
    fe = run_comparison(cae_dir, tlog.get("data_file") or args.data)
    record["fe_objective_MPa"] = fe.get("fe_objective_MPa")
    record["fe_peak_tension_MPa"] = fe.get("fe_peak_tension_MPa")
    record["fe_peak_compression_MPa"] = fe.get("fe_peak_compression_MPa")
    record["figure_png"] = fe.get("png")
    record["cae_dir"] = cae_dir
    append_manifest(args.manifest, record)
    log(f"DONE. surrogate={record['surrogate_objective_MPa']:.2f} MPa, "
        f"FE={record.get('fe_objective_MPa')}. manifest: {args.manifest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
