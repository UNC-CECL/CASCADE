#!/usr/bin/env python3
"""
The option A groin grid, M x f x window, run through the unchanged runner.

    python launch_groin_grid.py [--streams 3] [--only 1996:3:0.6] [--kind blocking] [--study TAG]

Sets only the groin and the filing (HAT_IGNORE_SETTINGS=1); every other value
is the code default. Each cell is its own member of the experiment,
raw_runs/experiments/<study>/<member>/; a cell already run is skipped, so a
relaunch resumes. Logs go to grid_logs/ beside this script.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

from __future__ import annotations

import argparse
import os
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
STUDY = "groin/2026-09-29-option-a-grid"   # the first grid; later studies pass --study
KIND = "dipole"                            # 1996->2003 linear ramp then; later studies pass --kind
RAW = REPO / "output" / "raw_runs" / "experiments"
LOGS = HERE / "grid_logs"

M_VALUES = (1, 2, 3, 4, 5, 7, 10)
F_VALUES = (0.2, 0.4, 0.6, 0.8, 1.0)
PERIODS = (1996, 2010)
# -----------------------------------------------------------------------------


# Experiment member for one cell: M<M>_f<f>, or b<b>_f<f> for blocking
def member(M, f):
    if KIND == "blocking":
        return f"b{M:.2f}_f{f:.1f}"
    return f"M{M}_f{f:.1f}"


# Has this cell already run (its run metadata exists)?
def done(period, M, f):
    return any((RAW / STUDY / member(M, f)).glob(
        f"{period}_*/*/*/*_run_metadata.json"))


# Run one cell through the runner, logging to grid_logs/; returns the exit code
def run(period, M, f):
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE",
        "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "true",
        "HAT_GROIN_KIND": KIND,
        **({"HAT_GROIN_BLOCKING_FRACTION": str(float(M))} if KIND == "blocking"
           else {"HAT_GROIN_TRAPPING_RATE_M_YR": str(float(M))}),
        "HAT_GROIN_DETERIORATION_FRACTION": str(f),
        "HAT_RUN_KIND": "experiment",
        "HAT_RUN_TAG": f"{STUDY}/{member(M, f)}",
        "HAT_SHOW_FIGURES": "false",
        "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false",
        "MPLBACKEND": "Agg",
    })
    log = LOGS / f"{STUDY.split('/')[-1]}_{period}_{member(M, f)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(RUNNER)], cwd=REPO,
                               env=env, stdout=fh, stderr=subprocess.STDOUT)
    print(f"{time.strftime('%H:%M:%S')}  {period} {member(M, f)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


# Run: build the cell list, skip finished cells, run the rest in parallel streams
def main():
    global KIND, STUDY
    ap = argparse.ArgumentParser()
    ap.add_argument("--streams", type=int, default=3)
    ap.add_argument("--only", help="one cell, period:M:f, for a smoke test")
    ap.add_argument("--M", type=float, nargs="+",
                    help="override M_VALUES (b values when --kind blocking)")
    ap.add_argument("--kind", choices=("dipole", "blocking"), default="dipole")
    ap.add_argument("--study", default="groin/2026-09-29-option-a-grid",
                    help="experiment tag under raw_runs/experiments/")
    ap.add_argument("--f", type=float, nargs="+", help="override F_VALUES")
    args = ap.parse_args()
    KIND, STUDY = args.kind, args.study
    LOGS.mkdir(exist_ok=True)

    if args.only:
        p, M, f = args.only.split(":")
        M = float(M) if KIND == "blocking" else int(M)
        cells = [(int(p), M, float(f))]
    else:
        # Low M first in both windows, so an early stop still leaves the plausible end
        Ms = [int(m) if float(m).is_integer() else m for m in (args.M or M_VALUES)]
        cells = [(p, M, f) for M in Ms for f in (args.f or F_VALUES)
                 for p in PERIODS]
    todo = [c for c in cells if not done(*c)]
    print(f"{len(cells) - len(todo)} of {len(cells)} cells already done; "
          f"running {len(todo)} in {args.streams} streams", flush=True)

    jobs = queue.Queue()
    for c in todo:
        jobs.put(c)
    failed = []

    # Take cells off the queue until it is empty
    def worker():
        while True:
            try:
                c = jobs.get_nowait()
            except queue.Empty:
                return
            if run(*c) != 0:
                failed.append(c)

    threads = [threading.Thread(target=worker) for _ in range(args.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"\nGRID DONE: {len(todo) - len(failed)}/{len(todo)} clean; "
          f"failed {failed}", flush=True)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
