"""
run_mono_ablation.py
====================

Surrogate-level ablation of the monotonic -> loops route over many specimens
and backstress counts (N = 3, 2, 1), Chaboche only. For each (specimen, N) it
runs ``run_mono_to_loops.py --no-fe`` (monotonic ramp calibration +
predicted-loop objective) and merges the per-task manifests into one JSONL.

Launches are staggered (>= 2.5 s apart) because session folders are named to
the second and concurrent starts on the same data name would collide.

    python run_mono_ablation.py [--workers 8] [--n 3 2 1] [--r 0.4]
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(PROJ, "calibration_sessions", "mono_ablation_parts")
MERGED = os.path.join(PROJ, "calibration_sessions", "mono_ablation_surrogate.jsonl")


def specimens(lds=(5,), strains=(3, 4, 5, 6)):
    return [(dia, ld, st) for dia in (10, 12, 16) for ld in lds for st in strains]


_lock = threading.Lock()
_last = [0.0]


def run_task(dia, ld, st, nbs, r, py, esrc="ramp"):
    path = f"input_data/{dia}mm dia/LD{ld}/{st}data.csv"
    tag = f"{dia}mm_LD{ld}_{st}pct_N{nbs}_E{esrc}_r{r:.2f}"
    out = os.path.join(PARTS, tag + ".jsonl")
    if os.path.exists(out):
        return tag, 0
    with _lock:                       # stagger launches
        wait = 2.5 - (time.time() - _last[0])
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()
        log = open(os.path.join(PARTS, tag + ".log"), "w")
        proc = subprocess.Popen(
            [py, os.path.join(PROJ, "scripts", "run_mono_to_loops.py"), "--cyclic", path,
             "--prior", "kinematic", "--bauschinger", str(r),
             "--n-backstresses", str(nbs), "--E-source", esrc, "--no-fe",
             "--manifest", out],
            cwd=PROJ, stdout=log, stderr=subprocess.STDOUT)
    rc = proc.wait()
    log.close()
    return tag, rc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--n", type=int, nargs="+", default=[3, 2, 1])
    ap.add_argument("--r", type=float, nargs="+", default=[0.4])
    ap.add_argument("--ld", type=int, nargs="+", default=[5])
    ap.add_argument("--strains", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--E-source", nargs="+", default=["ramp"],
                    choices=["ramp", "reversals", "handbook"])
    ap.add_argument("--python", default=sys.executable)
    args = ap.parse_args()
    os.makedirs(PARTS, exist_ok=True)
    tasks = [(d, l, s, n, r, e)
             for (d, l, s) in specimens(tuple(args.ld), tuple(args.strains))
             for n in args.n for r in args.r for e in args.E_source]
    print(f"{len(tasks)} tasks, {args.workers} workers", flush=True)
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(run_task, d, l, s, n, r, args.python, e)
                for (d, l, s, n, r, e) in tasks]
        for f in futs:
            tag, rc = f.result()
            done += 1
            print(f"[{done}/{len(tasks)}] {tag} rc={rc}", flush=True)
    n_rec = 0
    with open(MERGED, "w") as out:
        for p in sorted(glob.glob(os.path.join(PARTS, "*.jsonl"))):
            with open(p) as f:
                for line in f:
                    if line.strip():
                        out.write(line if line.endswith("\n") else line + "\n")
                        n_rec += 1
    print(f"merged {n_rec} records -> {MERGED}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
