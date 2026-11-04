"""
v6_fe.py
========

V6-C (see V6_PROTOCOL.md): the seed-1 "nonsoft" V6-A parameters in the ABAQUS coupon with bar E = 200 GPa and a
series spring (stiffness from the fitted lumped E), same set-up as V4-T. Mesh, geometry, loading and step are copied
from the V3-F input file of the same bar; only the material cards and the spring/E items change.

    python v6_fe.py --workers 8 [--only 10 4 2]        # d strain N

Records go to calibration_sessions/v6/c.jsonl.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PROJ = os.path.dirname(HERE)
sys.path.insert(0, PROJ)
sys.path.insert(0, HERE)

import run_full_workflow as rfw                                   # noqa: E402
import v4_spring_test as v4                                       # noqa: E402

A_JSONL = os.path.join(PROJ, "calibration_sessions", "v6", "a.jsonl")
OUT = os.path.join(PROJ, "calibration_sessions", "v6",
                   "c_soft.jsonl" if os.environ.get("V6_VARIANT") == "soft" else "c.jsonl")
VARIANT, SEED = os.environ.get("V6_VARIANT", "nonsoft"), 1
JOBS = [(d, n, N) for N in (1, 2, 3) for d in (10, 16) for n in (4, 5)]


def v6a_record(d, n, N):
    for line in open(A_JSONL):
        r = json.loads(line)
        if (r["dia"], r["strain"], r["N"], r["variant"], r["seed"]) == (d, n, N, VARIANT, SEED):
            return r
    raise KeyError((d, n, N))


def swap_material(text, N, sy0, params):
    """Replace the *Plastic (combined, N backstresses) and *Cyclic Hardening cards of a V3-F input file."""
    lines = text.splitlines()
    out, i = [], 0
    cg = params[: 2 * N]
    Q, b = params[2 * N], params[2 * N + 1]
    while i < len(lines):
        low = lines[i].strip().lower()
        if low.startswith("*plastic, hardening=combined"):
            out.append(f"*Plastic, hardening=COMBINED, datatype=PARAMETERS, number backstresses={N}")
            vals = [sy0] + list(cg)
            out.append(", ".join(f"{v:.6g}" for v in vals))
            i += 1
            while i < len(lines) and not lines[i].lstrip().startswith("*"):
                i += 1
            continue
        if low.startswith("*cyclic hardening"):
            out.append("*Cyclic Hardening, parameters")
            out.append(f"{sy0:.6g}, {Q:.6g}, {b:.6g}")
            i += 1
            while i < len(lines) and not lines[i].lstrip().startswith("*"):
                i += 1
            continue
        out.append(lines[i])
        i += 1
    return "\n".join(out) + "\n"


def run_job(job):
    d, n, N = job
    rec6 = v6a_record(d, n, N)
    src = v4.v3f_record(d, n, 1)          # same mesh/geometry; any N works as the template
    src_cae = os.path.join(src["session"], "cae")
    k_m = v4.spring_stiffness(d, rec6["E_fit"])
    sdir = os.path.join(PROJ, "calibration_sessions",
                        f"{time.strftime('%Y%m%d_%H%M%S')}_v6C_{VARIANT}_{d}mm_LD5_{n}pct_N{N}_seed{SEED}")
    cae = os.path.join(sdir, "cae")
    os.makedirs(cae)
    text = open(os.path.join(src_cae, "Calibration_Job.inp")).read()
    text = swap_material(text, N, rec6["sigma_y0"], rec6["params"])
    open(os.path.join(cae, "Calibration_Job.inp"), "w").write(v4.edit_inp(text, k_m))
    shutil.copy(os.path.join(src_cae, "transfer_log.json"), cae)
    shutil.copy(os.path.join(src["session"], "session.json"), sdir)
    t0 = time.time()
    proc = subprocess.run("abaqus job=Calibration_Job input=Calibration_Job.inp interactive",
                          cwd=cae, shell=True, capture_output=True, text=True)
    sta = os.path.join(cae, "Calibration_Job.sta")
    ok = os.path.exists(sta) and "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in open(sta, errors="ignore").read()
    rec = {"exp": "c", "dia": d, "strain": n, "N": N, "variant": VARIANT, "seed": SEED, "session": sdir,
           "E_fit": rec6["E_fit"], "K_m_N_per_mm": k_m, "sigma_y0": rec6["sigma_y0"], "params": rec6["params"],
           "point_val_rmse_eng": rec6["val_rmse_eng"], "solve_ok": bool(ok),
           "v3f_fe_4_6": src.get("fe_4_6") if N == 1 else None}
    if not ok:
        rec["abaqus_tail"] = (proc.stdout or "")[-1500:] + (proc.stderr or "")[-500:]
    else:
        rfw.run_comparison(cae, f"input_data/{d}mm dia/LD5/{n}data.csv")
        fj = os.path.join(cae, "fe_comparison.json")
        if os.path.exists(fj):
            j = json.load(open(fj))
            det = j.get("fe_cycle_details", {})
            rec["cycles"] = det
            r_all = [v["rmse"] for v in det.values() if v["rmse"] == v["rmse"]]
            r_val = [v["rmse"] for k, v in det.items() if int(k) >= 4 and v["rmse"] == v["rmse"]]
            rec["fe_2_6"] = float(np.mean(r_all)) if r_all else float("nan")
            rec["fe_4_6"] = float(np.mean(r_val)) if r_val else float("nan")
            rec["fe_objective_local_MPa"] = j.get("fe_objective_local_MPa")
    rec["wall_s"] = time.time() - t0
    with open(OUT, "a") as f:
        f.write(json.dumps(rec) + "\n")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", type=int, nargs=3, metavar=("D", "STRAIN", "N"))
    ap.add_argument("--dry", action="store_true")
    ap.add_argument("--variant", default="nonsoft", choices=("nonsoft", "soft"))
    ap.add_argument("--dias", type=int, nargs="+", default=[10, 16])
    args = ap.parse_args()
    if args.variant != VARIANT:
        sys.exit("set env V6_VARIANT=%s (worker processes re-import this module)" % args.variant)
    jobs = [tuple(args.only)] if args.only else [(d, n, N) for N in (1, 2, 3) for d in args.dias for n in (4, 5)]
    if args.dry:
        d, n, N = jobs[0]
        rec6 = v6a_record(d, n, N)
        src = v4.v3f_record(d, n, 1)
        txt = open(os.path.join(src["session"], "cae", "Calibration_Job.inp")).read()
        new = v4.edit_inp(swap_material(txt, N, rec6["sigma_y0"], rec6["params"]), v4.spring_stiffness(d, rec6["E_fit"]))
        for key in ("*Elastic", "*Plastic", "*Cyclic", "*Spring"):
            k = new.find(key)
            print("----", key); print(new[k:k + 220])
        return
    done = set()
    if os.path.exists(OUT):
        done = {(r["dia"], r["strain"], r["N"]) for r in map(json.loads, open(OUT)) if r.get("solve_ok")}
    todo = [j for j in jobs if j not in done]
    print(f"V6-C: {len(todo)} jobs, {args.workers} concurrent", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = []
        for j in todo:
            futs.append(ex.submit(run_job, j))
            time.sleep(4)
        for f in futs:
            r = f.result()
            print(f"  done {r['dia']}mm {r['strain']}% N={r['N']} ok={r['solve_ok']} fe_4_6={r.get('fe_4_6')} wall={r['wall_s']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
