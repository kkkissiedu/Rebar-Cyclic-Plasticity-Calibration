"""
run_calibration_headless.py
===========================

Reproducible, GUI-free calibration runs for paper numbers.

Re-runs the pipeline with the exact config + bounds of an existing session
(optionally overriding the model / backstress count), creating a NEW session
directory as always -- nothing is overwritten. Use it to produce the
N=2 vs N=3 justification runs and any reproducibility repeats.

Usage (from the repo root):

    python scripts/run_calibration_headless.py --from-session calibration_sessions/<dir>
        [--n-backstresses 2] [--model chaboche|uvc|ohno_wang]
        [--no-stage3] [--tag "paper N2"]

Prints stage progress and a final summary (params, objective, sigma_y0).
"""

from __future__ import annotations

import argparse
import json
import os
import queue
import sys
import threading

# repo root = parent of scripts/
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration_app.core import data_loader, pipeline as pipeline_mod
from calibration_app.session.session_manager import SessionManager

E_FIXED = 200_000.0
NU_FIXED = 0.30


def drain(q: queue.Queue, stop: threading.Event) -> None:
    """Print pipeline messages until the 'done' message arrives."""
    while True:
        try:
            kind, p = q.get(timeout=0.5)
        except queue.Empty:
            if stop.is_set():
                return
            continue
        if kind == "status":
            print(p["text"])
        elif kind == "stage":
            if p.get("status") in ("DONE", "FAILED"):
                print(f"  stage {p['key']}: {p['status']} "
                      f"({p.get('done', 0)}/{p.get('total', 0)}, "
                      f"best={p.get('best')})")
        elif kind == "model_done":
            print(f"model_done: best={p.get('best_obj')} "
                  f"sigma_y0={p.get('best_sigma_y0')}")
            print(f"session: {p.get('session')}")
        elif kind == "done":
            stop.set()
            return


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from-session", required=True,
                    help="existing session dir whose config/bounds to reuse")
    ap.add_argument("--n-backstresses", type=int, default=None)
    ap.add_argument("--model", default=None,
                    choices=["chaboche", "uvc", "ohno_wang"])
    ap.add_argument("--no-stage3", action="store_true",
                    help="skip the single-element ABAQUS verification")
    ap.add_argument("--tag", default="",
                    help="note stored in material_description for the run")
    args = ap.parse_args()

    with open(os.path.join(args.from_session, "session.json"), "r") as f:
        state = json.load(f)
    cfg = dict(state["config"])
    bounds = {k: tuple(v) for k, v in state.get("bounds", {}).items()}

    if args.n_backstresses:
        cfg["n_backstresses"] = int(args.n_backstresses)
    if args.model:
        cfg["model"] = args.model
    if args.no_stage3:
        cfg["run_stage3"] = False
    if args.tag:
        cfg["material_description"] = (cfg.get("material_description", "")
                                       + " " + args.tag).strip()
    cfg["run_all_models"] = False

    data = data_loader.load_experimental(
        cfg["data_file"], gauge_length_mm=cfg["gauge_length_mm"],
        bar_diameter_mm=cfg["diameter_mm"], sigma_y0=cfg["sigma_y0"])
    print(f"Loaded {os.path.basename(cfg['data_file'])}: "
          f"{data['cycle_summary']['total_cycles']} cycles, "
          f"peak {data['measured_peak_stress']:.1f} MPa")
    print(f"Model {cfg['model']} N={cfg['n_backstresses']}, "
          f"backend {cfg['backend']}, optimiser {cfg['optimiser']}, "
          f"stage3={'on' if cfg.get('run_stage3') else 'off'}")

    q: queue.Queue = queue.Queue()
    stop = threading.Event()
    done = threading.Event()
    printer = threading.Thread(target=drain, args=(q, done), daemon=True)
    printer.start()

    pipe = pipeline_mod.CalibrationPipeline(
        q, stop, SessionManager(), data, cfg, bounds,
        E=E_FIXED, nu=NU_FIXED)
    pipe.run()                      # blocking; messages drained by printer
    done.set()
    printer.join(timeout=5)

    # Summarise the newest session (the one this run just created).
    sm = SessionManager()
    newest = sm.list_sessions()[0]
    sess = sm.open(newest.path)
    st = sess.state
    print("\n=== RESULT ===")
    print(f"session:   {newest.path}")
    print(f"objective: {st.get('best_objective'):.3f} MPa (surrogate)")
    fe = st.get("fe_verification") or {}
    if fe:
        print(f"FE check:  {fe.get('objective_MPa'):.3f} MPa "
              f"(peaks {fe.get('peak_tension_MPa'):+.0f}/"
              f"{fe.get('peak_compression_MPa'):+.0f})")
    print(f"sigma_y0:  {st.get('best_sigma_y0'):.2f} MPa cyclic "
          f"(ratio {st.get('bauschinger_ratio')})")
    model = pipeline_mod.build_model(st["model"], cfg["n_backstresses"])
    for name, val in zip(model.param_names, st.get("best_params", [])):
        print(f"  {name:8s} = {val:.6g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
