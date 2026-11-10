"""
run_batch.py
============

Drives run_full_workflow.py through a small, explicitly-listed matrix of
cases, sequentially (one ABAQUS full-coupon solve at a time -- these are
real FE jobs, not something to fan out in parallel on one machine). A
failure in one case is logged and the batch moves on to the next; nothing
stops the whole run.

Usage (cuda_pt python, from the project root):

    python run_batch.py --cases batch_cases.json
        [--wait-timeout 14400] [--manifest calibration_sessions/website_results.jsonl]

``batch_cases.json`` is a JSON list of dicts, each with at least
``data``, ``model``, ``sy0``; any other run_full_workflow.py flag can be
included per-case (e.g. ``"n_backstresses": 3``, ``"tag": "..."``).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

PROJ = os.path.dirname(os.path.abspath(__file__))


def log(msg: str) -> None:
    print(f"[batch] {msg}", flush=True)


def run_case(case: dict, wait_timeout: int, manifest: str, py: str) -> bool:
    cmd = [py, os.path.join(PROJ, "run_full_workflow.py"),
          "--data", case["data"], "--model", case["model"],
          "--sy0", str(case["sy0"]),
          "--n-backstresses", str(case.get("n_backstresses", 3)),
          "--optimiser", case.get("optimiser", "bayesian"),
          "--wait-timeout", str(case.get("wait_timeout", wait_timeout)),
          "--manifest", manifest]
    if case.get("tag"):
        cmd += ["--tag", case["tag"]]
    if case.get("n_fit"):
        cmd += ["--n-fit", str(case["n_fit"])]
    if case.get("no_stage3"):
        cmd += ["--no-stage3"]
    if case.get("E_from_data"):
        cmd += ["--E-from-data"]
    log("START: " + " ".join(cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, cwd=PROJ)
    dt = time.time() - t0
    ok = proc.returncode == 0
    log(f"{'OK' if ok else 'FAILED'} ({dt/60.0:.1f} min): {case['model']} "
        f"{os.path.basename(case['data'])}")
    return ok


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True, help="JSON file: list of case dicts")
    ap.add_argument("--wait-timeout", type=int, default=4 * 3600)
    ap.add_argument("--manifest", default=os.path.join(
        "calibration_sessions", "website_results.jsonl"))
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()

    with open(args.cases, "r") as f:
        cases = json.load(f)
    log(f"{len(cases)} case(s) to run, sequentially")

    results = []
    for i, case in enumerate(cases, 1):
        log(f"--- case {i}/{len(cases)} ---")
        ok = run_case(case, args.wait_timeout, args.manifest, args.python)
        results.append((case, ok))

    n_ok = sum(1 for _, ok in results if ok)
    log(f"BATCH DONE: {n_ok}/{len(results)} succeeded. Manifest: {args.manifest}")
    for case, ok in results:
        log(f"  {'OK  ' if ok else 'FAIL'} {case['model']:10s} {case['data']}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
