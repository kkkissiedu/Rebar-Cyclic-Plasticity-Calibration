"""
v3_fe.py
========

V3-F (see V3_PROTOCOL.md): put the frozen V3-A parameter sets (seed 1; N = 1 and 3; 10/12/16 mm at
4 and 5 %) into the full 3-D ABAQUS coupon and compare with the measured loops.

    python v3_fe.py --workers 6 [--only 12 5 1]       # d strain N

Per job: derived session -> Transfer to CAE + solve (run_full_workflow.run_transfer_and_solve) ->
plot_fe_comparison.py -> one record appended to calibration_sessions/v3/f.jsonl.
"""

from __future__ import annotations

import argparse
import copy
import glob
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)
sys.path.insert(0, HERE)

import run_full_workflow as rfw                                   # noqa: E402
import v3_study as v3                                             # noqa: E402
from calibration_app.core.pipeline import build_model            # noqa: E402

OUT = os.path.join(PROJ, "calibration_sessions", "v3", "f.jsonl")
SEED = 1
MAX_INC_DIV = 200.0
CAE = os.path.join(PROJ, "FATIGUE_rebuilt.cae")
JOBS = [(d, n, N) for N in (1, 3) for d in (10, 12, 16) for n in (4, 5)]


def v3_record(d, n, N):
    for line in open(os.path.join(PROJ, "calibration_sessions", "v3", "a.jsonl")):
        r = json.loads(line)
        if (r["dia"], r["strain"], r["N"], r["seed"]) == (d, n, N, SEED):
            return r
    raise KeyError((d, n, N))


def template_state():
    p = sorted(glob.glob(os.path.join(PROJ, "calibration_sessions",
                                      "*10mm_LD5_5pct_kinematic_N1_r0.70_loops", "session.json")))[0]
    return json.load(open(p))


def surrogate_scores(d, n, N, rec):
    data = v3.load(d, n)
    E = v3.measure_E(data)
    model = build_model("chaboche", N)
    out = {}
    for name, tgt in (("sur_4_6", [4, 5, 6]), ("sur_2_6", [2, 3, 4, 5, 6])):
        out[name] = v3.score(model, rec["params"], rec["sigma_y0"], v3.window(data, tgt), E)
    return out, E


def run_job(job):
    d, n, N = job
    rec = v3_record(d, n, N)
    sur, E = surrogate_scores(d, n, N, rec)
    case = f"v3F_{d}mm_LD5_{n}pct_N{N}_seed{SEED}"
    sdir = os.path.join(PROJ, "calibration_sessions", f"{time.strftime('%Y%m%d_%H%M%S')}_{case}")
    os.makedirs(sdir)
    st = copy.deepcopy(template_state())
    cfg = st["config"]
    cfg.update({
        "data_file": f"input_data/{d}mm dia/LD5/{n}data.csv", "diameter_mm": float(d),
        "gauge_length_mm": float(d * 5), "strain_pct": float(n), "n_backstresses": N,
        "n_fit": 5, "E_modulus": float(E), "weight_reversals": False,
        "max_inc_divisor": MAX_INC_DIV, "loading_mode": "cyclic", "cae_path": CAE,
        "material_description": f"V3-F frozen V3-A fit, seed {SEED}, N={N}",
    })
    st.update({"name": os.path.basename(sdir), "data_file": cfg["data_file"], "config": cfg,
               "best_params": list(rec["params"]), "best_sigma_y0": float(rec["sigma_y0"]),
               "best_objective": sur["sur_2_6"],
               "material_description": cfg["material_description"]})
    for k in ("fe_verification", "best_peak_tension", "best_peak_compression",
              "sigma_y0_monotonic", "bauschinger_ratio"):
        st.pop(k, None)
    json.dump(st, open(os.path.join(sdir, "session.json"), "w"), indent=1)

    t0 = time.time()
    record = {"exp": "f", "dia": d, "strain": n, "N": N, "seed": SEED, "E_meas": E,
              "sigma_y0": rec["sigma_y0"], "params": rec["params"], "session": sdir,
              "max_inc_divisor": MAX_INC_DIV, **sur}
    ok, tlog = rfw.run_transfer_and_solve(sdir, CAE, 6 * 3600)
    record["solve_ok"] = bool(ok)
    if ok:
        cae_dir = os.path.dirname(tlog["cae_path"])
        rfw.run_comparison(cae_dir, tlog.get("data_file") or cfg["data_file"])
        fj = os.path.join(cae_dir, "fe_comparison.json")
        if os.path.exists(fj):
            j = json.load(open(fj))
            det = j.get("fe_cycle_details", {})
            record["cycles"] = det
            r_all = [v["rmse"] for k, v in det.items() if v["rmse"] == v["rmse"]]
            r_val = [v["rmse"] for k, v in det.items() if int(k) >= 4 and v["rmse"] == v["rmse"]]
            record["fe_2_6"] = float(np.mean(r_all)) if r_all else float("nan")
            record["fe_4_6"] = float(np.mean(r_val)) if r_val else float("nan")
            record["fe_objective_local_MPa"] = j.get("fe_objective_local_MPa")
            record["fe_peak_tension_MPa"] = j.get("fe_peak_tension_MPa")
            record["fe_peak_compression_MPa"] = j.get("fe_peak_compression_MPa")
            record["png"] = j.get("png")
    record["wall_s"] = time.time() - t0
    with open(OUT, "a") as f:
        f.write(json.dumps(record) + "\n")
    return record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--only", type=int, nargs=3, metavar=("D", "STRAIN", "N"))
    args = ap.parse_args()
    jobs = [tuple(args.only)] if args.only else JOBS
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            if line.strip():
                r = json.loads(line)
                done.add((r["dia"], r["strain"], r["N"]))
    todo = [j for j in jobs if j not in done]
    print(f"V3-F: {len(todo)} jobs to run, {args.workers} concurrent", flush=True)
    # Stagger submissions so the transfer step and licence checkouts do not collide.
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = []
        for j in todo:
            futs.append(ex.submit(run_job, j))
            time.sleep(4)
        for f in futs:
            r = f.result()
            print(f"  done {r['dia']}mm {r['strain']}% N={r['N']} solve_ok={r['solve_ok']} "
                  f"fe_4_6={r.get('fe_4_6')} wall={r['wall_s']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
