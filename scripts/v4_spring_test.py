"""
v4_spring_test.py
=================

V4-T (see V4_PROTOCOL.md): re-run frozen V3-F coupon jobs with the bar at E = 200 GPa and the machine
compliance as a spring in series between the actuator and the top face. Material cards, mesh, loading
and step are copied unchanged from the V3-F input files; only the items listed in the protocol change.

    python v4_spring_test.py --workers 8 [--only 10 4 1]        # d strain N

Per job: edit the V3-F input file -> abaqus job (interactive) -> plot_fe_comparison.py -> one record in
calibration_sessions/v3/v4t.jsonl.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
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

V3F = os.path.join(PROJ, "calibration_sessions", "v3", "f.jsonl")
OUT = os.path.join(PROJ, "calibration_sessions", "v3", "v4t.jsonl")
E_BAR = 200000.0
JOBS = [(d, n, N) for N in (1, 3) for d in (10, 16) for n in (4, 5)]
ACT = 900001


def v3f_record(d, n, N):
    for line in open(V3F):
        r = json.loads(line)
        if (r["dia"], r["strain"], r["N"]) == (d, n, N) and r["solve_ok"]:
            return r
    raise KeyError((d, n, N))


def spring_stiffness(d, E_eff):
    """K_m (N/mm) from 1/K_m = (L0/A0)(1/E_eff - 1/E_bar)."""
    L0 = 5.0 * d
    A0 = math.pi * d * d / 4.0
    return 1.0 / ((L0 / A0) * (1.0 / E_eff - 1.0 / E_BAR))


def edit_inp(text: str, k_m: float) -> str:
    lines = text.splitlines()
    out = []
    i = 0
    first = last = None
    n_top_eq = 0
    while i < len(lines):
        ln = lines[i]
        low = ln.strip().lower()
        # 1. bar modulus
        if low == "*elastic":
            out.append(ln)
            out.append(f" {E_BAR:.1f}, 0.3")
            i += 2
            continue
        # 2. replace the top-face history/reaction set by the actuator node + spring
        if low.startswith("*nset, nset=top_auto"):
            j = i + 1
            ids = []
            while j < len(lines) and not lines[j].lstrip().startswith("*"):
                ids += [int(t) for t in lines[j].replace(",", " ").split()]
                j += 1
            if "generate" in low:                     # "first, last, step"
                top_nodes = list(range(ids[0], ids[1] + 1, ids[2] if len(ids) > 2 else 1))
            else:                                     # explicit node list (12 mm deck)
                top_nodes = ids
            first = top_nodes[0]
            out += [
                "*Node",
                f" {ACT}, 0., 0., 0.",
                "*Nset, nset=TOP_AUTO",
                f" {ACT},",
                "*Nset, nset=ACT_NODE",
                f" {ACT},",
                "*Element, type=SPRING2, elset=MACH_SPRING",
                f" {ACT}, Specimen-1.{first}, {ACT}",
                "*Spring, elset=MACH_SPRING", "3, 3", f" {k_m:.6g}",
            ]
            i = j
            continue
        # 3. after the assembly: spring property + plane-grip equations
        if low == "*end assembly":
            for nid in top_nodes[1:]:
                out += ["*Equation", "2", f"Specimen-1.{nid}, 3, 1., Specimen-1.{first}, 3, -1."]
                n_top_eq += 1
            out.append(ln)
            i += 1
            continue
        # 4. drive the actuator node, not the top face
        if low.startswith("*boundary, amplitude=cyclicstrain"):
            out.append(ln)
            lat = lines[i + 1].split(",")[0].strip()      # lateral top constraint set (name differs per deck)
            out += [f"{lat}, 1, 1", f"{lat}, 2, 2", "ACT_NODE, 1, 1", "ACT_NODE, 2, 2",
                    "ACT_NODE, 3, 3, 1."]
            i += 4
            continue
        out.append(ln)
        i += 1
    if first is None or n_top_eq == 0:
        raise RuntimeError("top node set not found in input file")
    return "\n".join(out) + "\n"


def run_job(job):
    d, n, N = job
    src = v3f_record(d, n, N)
    src_cae = os.path.join(src["session"], "cae")
    k_m = spring_stiffness(d, src["E_meas"])
    sdir = os.path.join(PROJ, "calibration_sessions",
                        f"{time.strftime('%Y%m%d_%H%M%S')}_v4T_{d}mm_LD5_{n}pct_N{N}_seed1")
    cae = os.path.join(sdir, "cae")
    os.makedirs(cae)
    text = open(os.path.join(src_cae, "Calibration_Job.inp")).read()
    open(os.path.join(cae, "Calibration_Job.inp"), "w").write(edit_inp(text, k_m))
    shutil.copy(os.path.join(src_cae, "transfer_log.json"), cae)
    shutil.copy(os.path.join(src["session"], "session.json"), sdir)
    t0 = time.time()
    proc = subprocess.run("abaqus job=Calibration_Job input=Calibration_Job.inp interactive",
                          cwd=cae, shell=True, capture_output=True, text=True)
    sta = os.path.join(cae, "Calibration_Job.sta")
    ok = os.path.exists(sta) and "THE ANALYSIS HAS COMPLETED SUCCESSFULLY" in open(sta, errors="ignore").read()
    rec = {"exp": "v4t", "dia": d, "strain": n, "N": N, "session": sdir, "K_m_N_per_mm": k_m,
           "E_bar": E_BAR, "E_eff_V3F": src["E_meas"], "solve_ok": bool(ok),
           "v3f_fe_4_6": src.get("fe_4_6"), "v3f_fe_2_6": src.get("fe_2_6"),
           "sur_4_6": src.get("sur_4_6"), "sur_2_6": src.get("sur_2_6"),
           "v3f_cycles": src.get("cycles")}
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
            rec["png"] = j.get("png")
    rec["wall_s"] = time.time() - t0
    with open(OUT, "a") as f:
        f.write(json.dumps(rec) + "\n")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--only", type=int, nargs=3, metavar=("D", "STRAIN", "N"))
    ap.add_argument("--dry", action="store_true", help="only write/print the edited input, no solve")
    args = ap.parse_args()
    jobs = [tuple(args.only)] if args.only else JOBS
    if args.dry:
        d, n, N = jobs[0]
        src = v3f_record(d, n, N)
        txt = open(os.path.join(src["session"], "cae", "Calibration_Job.inp")).read()
        new = edit_inp(txt, spring_stiffness(d, src["E_meas"]))
        print("K_m =", spring_stiffness(d, src["E_meas"]), "N/mm; lines", len(txt.splitlines()), "->", len(new.splitlines()))
        for key in ("*Elastic", "*Node\n", "*Spring", "*Boundary, amplitude"):
            k = new.find(key)
            print("----", key.strip()); print(new[k:k + 330])
        return
    done = set()
    if os.path.exists(OUT):
        done = {(r["dia"], r["strain"], r["N"]) for r in map(json.loads, open(OUT)) if r.get("solve_ok")}
    todo = [j for j in jobs if j not in done]
    print(f"V4-T: {len(todo)} jobs, {args.workers} concurrent", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = []
        for j in todo:
            futs.append(ex.submit(run_job, j))
            time.sleep(4)
        for f in futs:
            r = f.result()
            print(f"  done {r['dia']}mm {r['strain']}% N={r['N']} ok={r['solve_ok']} "
                  f"fe_4_6={r.get('fe_4_6')} v3f={r['v3f_fe_4_6']} wall={r['wall_s']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
