"""
run_mono_to_loops.py
====================

Monotonic data in -> full hysteresis loops out.

For one cyclic test file (Kashani layout) this script:

  1. cuts the file's virgin opening ramp (cycle 1, start -> first reversal)
     out as a monotonic tension curve (synthetic_data/.../mono_from_<N>data.csv);
  2. calibrates the model on that monotonic curve ALONE (surrogate, Bayesian);
  3. runs the file's real cyclic strain history (burn-in + cycles 2..) through
     the calibrated model and scores the predicted loops against the
     experimental loops (surrogate objective);
  4. optionally writes the parameters into the full-coupon ABAQUS model under
     the cyclic loading, solves, and runs the standard post-run comparison, so
     the output (fe_comparison.png / .json) is identical in form to a run that
     was calibrated on the cyclic data.

A monotonic curve fixes the TOTAL hardening but not its kinematic / isotropic
split, so ``--prior`` chooses the split:
  kinematic : isotropic hardening off (Q_inf ~ 0); the monotonic curve fixes
              C_k, gamma_k and sigma_y0 and the loops follow from the model.
  free      : Q_inf and b stay free in their default bounds.

The apparent elastic modulus is taken from the file's own reversals by default
(``--no-E-from-data`` restores 200 GPa).

    python run_mono_to_loops.py --cyclic "input_data/12mm dia/LD5/5data.csv" \\
        --prior kinematic [--no-fe] [--manifest calibration_sessions/mono2loops.jsonl]
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import types

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from calibration_app.core import data_loader, objective
from calibration_app.core import pipeline as pipeline_mod
from calibration_app.session.session_manager import SessionManager
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_full_workflow as rfw

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REVERSAL_MM = 0.03          # displacement drop that marks the first reversal


def log(msg: str) -> None:
    print(f"[mono2loops] {msg}", flush=True)


def extract_monotonic(cyc_path: str, out_csv: str, gauge_mm: float) -> float:
    """Write the cycle-1 virgin ramp (zeroed displacement and load) as a
    monotonic CSV; return the peak strain in percent."""
    rows = []
    with open(cyc_path, newline="") as f:
        reader = csv.reader(f)
        hdr = None
        for raw in reader:
            if raw and raw[0].strip().strip('"').lower().startswith("notes"):
                hdr = [c.strip().strip('"').lower() for c in raw]
                break
        it, ic = hdr.index("time"), hdr.index("cycle")
        ip, il = hdr.index("position mm"), hdr.index("load kn")
        for raw in reader:
            try:
                rows.append((float(raw[it]), int(float(raw[ic])),
                             float(raw[ip]), float(raw[il])))
            except (ValueError, IndexError):
                continue
    c1 = [(t, p, l) for t, c, p, l in rows if c == 1]
    t = np.array([r[0] for r in c1])
    disp = np.array([r[1] for r in c1]) - c1[0][1]
    load = np.array([r[2] for r in c1]) - c1[0][2]
    drop = np.nonzero(np.maximum.accumulate(disp) - disp > REVERSAL_MM)[0]
    stop = int(drop[0]) if drop.size else len(disp)
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Time", "Displacement_mm", "Load_kN"])
        for i in range(stop):
            w.writerow([t[i], disp[i], load[i]])
    return float(disp[:stop].max() / gauge_mm * 100.0)


def make_args(data_path, tag, E, cae_path, sy0_lower, sy0=500.0, nbs=3):
    return types.SimpleNamespace(
        data=data_path, tag=tag, sy0=sy0, no_calibrate_sy0=False,
        sy0_lower=sy0_lower, model="chaboche", n_backstresses=nbs,
        backend="surrogate", optimiser="bayesian", n_fit=None,
        sobol_pow2=12, n_trials=200, no_stage3=True, cae_path=cae_path,
        E_resolved=E)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cyclic", required=True,
                    help="Kashani cyclic CSV; its cycle-1 ramp is the monotonic input")
    ap.add_argument("--prior", choices=["kinematic", "free"], default="kinematic")
    ap.add_argument("--no-E-from-data", action="store_true",
                    help="same as --E-source handbook")
    ap.add_argument("--E-source", choices=["ramp", "reversals", "handbook"],
                    default="ramp",
                    help="ramp: slope of the monotonic curve itself (monotonic "
                         "data only, default); reversals: apparent modulus from "
                         "the cyclic file's unloading (uses cyclic data); "
                         "handbook: 200 GPa")
    ap.add_argument("--n-backstresses", type=int, default=3, choices=[1, 2, 3, 4])
    ap.add_argument("--sy0-lower", type=float, default=100.0)
    ap.add_argument("--bauschinger", type=float, default=None,
                    help="Bauschinger prior r: fix the cyclic yield radius at "
                         "r x 500 MPa (500 = the monotonic proof stress used "
                         "as the usual upper bound); the backstresses then "
                         "absorb the rest of the monotonic curve")
    ap.add_argument("--no-fe", action="store_true",
                    help="surrogate-level prediction only (seconds)")
    ap.add_argument("--cae-path", default=os.path.join(PROJ, "FATIGUE_rebuilt.cae"))
    ap.add_argument("--wait-timeout", type=int, default=3600)
    ap.add_argument("--manifest", default=os.path.join(
        "calibration_sessions", "mono2loops_results.jsonl"))
    args = ap.parse_args()

    geo = data_loader.detect_geometry(args.cyclic)
    dia, gauge = geo["diameter_mm"], geo["gauge_length_mm"]
    ld = int(round(gauge / dia))
    cyc = data_loader.load_experimental(args.cyclic, gauge_length_mm=gauge,
                                        bar_diameter_mm=dia, sigma_y0=500.0)
    cyc["loading_mode"] = "cyclic"
    e_src = "handbook" if args.no_E_from_data else args.E_source
    E = data_loader.apparent_modulus(cyc) if e_src == "reversals" else 200000.0
    n = os.path.basename(args.cyclic).split("data")[0]
    sy0_in, sy0_lo = 500.0, args.sy0_lower
    if args.bauschinger:
        sy0_in = args.bauschinger * 500.0
        sy0_lo = sy0_in * 0.995
    case = f"{int(dia)}mm_LD{ld}_{n}pct_{args.prior}_N{args.n_backstresses}" + (
        f"_r{args.bauschinger:.2f}" if args.bauschinger else "")
    log(f"{case}: E = {E:.0f} MPa, sigma_y0 search ({sy0_lo:.0f}, {sy0_in:.0f})")

    mono_csv = os.path.join("synthetic_data", f"{int(dia)}mm dia", f"LD{ld}",
                            f"mono_from_{n}data.csv")
    peak = extract_monotonic(args.cyclic, mono_csv, gauge)
    log(f"monotonic input: {mono_csv} (virgin ramp to {peak:.2f}% strain)")

    mono = data_loader.load_any_experimental(
        mono_csv, gauge_length_mm=gauge, bar_diameter_mm=dia, sigma_y0=500.0)
    if e_src == "ramp":
        E = data_loader.monotonic_modulus(mono["strain"], mono["stress"]) or 200000.0
    log(f"elastic modulus ({e_src}): {E:.0f} MPa")
    tag = f"mono2loops {args.prior}" + (
        f" r={args.bauschinger:.2f}" if args.bauschinger else "")
    cfg = rfw.build_cfg(make_args(mono_csv, tag, E, args.cae_path, sy0_lo, sy0_in,
                                  args.n_backstresses), mono)
    bounds = {"Q_inf": (0.0, 1e-6)} if args.prior == "kinematic" else {}
    mono_sess = rfw.run_calibration(cfg, mono, bounds)
    st = SessionManager().open(mono_sess).state
    params, sy0 = st["best_params"], st["best_sigma_y0"]
    log(f"monotonic fit: objective {st['best_objective']:.2f} MPa, sigma_y0 {sy0:.1f}")

    cfg_cyc = rfw.build_cfg(make_args(args.cyclic, tag, E, args.cae_path, sy0_lo, sy0_in,
                                      args.n_backstresses), cyc)
    win = data_loader.slice_window(cyc, cfg_cyc["n_fit"])
    model = pipeline_mod.build_model("chaboche", args.n_backstresses)
    obj_loops = objective.evaluate_surrogate(model, params, win, None, E=E, sy0=sy0)
    log(f"PREDICTED LOOPS (surrogate), cycles scored {win['target_cycles']}: "
        f"{obj_loops:.2f} MPa")

    record = {
        "kind": "mono_to_loops", "case": case, "cyclic_file": args.cyclic,
        "prior": args.prior, "n_backstresses": args.n_backstresses,
        "E_source": e_src, "bauschinger": args.bauschinger,
        "E_modulus": E, "mono_csv": mono_csv,
        "mono_peak_strain_pct": peak, "mono_session": mono_sess,
        "mono_fit_objective_MPa": st["best_objective"], "sigma_y0": sy0,
        "params": params, "loops_surrogate_objective_MPa": obj_loops,
        "n_fit": cfg_cyc["n_fit"],
    }

    if not args.no_fe:
        name = f"{time.strftime('%Y%m%d_%H%M%S')}_{case}_loops"
        new_dir = os.path.join("calibration_sessions", name)
        os.makedirs(new_dir)
        state2 = dict(st)
        state2.update({"name": name, "data_file": args.cyclic, "config": cfg_cyc,
                       "best_objective": obj_loops,
                       "material_description": f"monotonic-calibrated ({args.prior}) "
                                               f"-> cyclic loops; mono session {mono_sess}"})
        state2.pop("fe_verification", None)
        with open(os.path.join(new_dir, "session.json"), "w") as f:
            json.dump(state2, f, indent=1)
        ok, tlog = rfw.run_transfer_and_solve(new_dir, args.cae_path, args.wait_timeout)
        record["loops_session"] = new_dir
        record["transfer_ok"] = ok
        if ok:
            cae_dir = os.path.dirname(tlog["cae_path"])
            fe = rfw.run_comparison(cae_dir, tlog.get("data_file") or args.cyclic)
            record["fe_objective_MPa"] = fe.get("fe_objective_MPa")
            record["fe_peak_tension_MPa"] = fe.get("fe_peak_tension_MPa")
            record["fe_peak_compression_MPa"] = fe.get("fe_peak_compression_MPa")
            record["figure_png"] = fe.get("png")
            log(f"FULL-COUPON LOOPS: FE objective {record['fe_objective_MPa']} MPa")
    rfw.append_manifest(args.manifest, record)
    log("DONE " + case)
    return 0


if __name__ == "__main__":
    sys.exit(main())
