#!/usr/bin/env python3
"""The option A groin grid, M x f x window, run through the unchanged runner.

WHY THE RUNNER AND NOT HAT_groin_sweep_worker. The worker is a second copy of
the runner, and as of 2026-09-29 it is not the adopted model: it hardcodes the
old waves (2.5 / 8 / 0.7 / 0.1), and its parameter-file repair restores a
snapshot with no per-cell dune ceilings. Driving the runner means every cell
is the adopted model by construction. Runs write their own parameters file,
so they are safe alongside other sessions' runs (Hannah, 2026-09-29).

WHAT IS SET. Only the groin and the filing; every other value is the code
default (HAT_IGNORE_SETTINGS=1, stray HAT_* dropped, as HAT_run_all does):
edgeBE, full_management without historical relocations, option A waves,
v3_trim24 storms. M = 0 is not run: the paired baseline is the adopted matrix
run with the same tokens, which the runner resolves by itself.
    1996  HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin
    2010  HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin

The M = 2-5, f = 0.6 cells are the confirmation runs for the emulator in
groin_stability_option_a.py.

FILING. Run names carry `groin` but not M or f, so each cell is its own
experiment member: raw_runs/experiments/groin/2026-09-29-option-a-grid/M<M>_f<f>/.
A cell whose metadata already exists is skipped, so a relaunch resumes.

    python launch_groin_grid.py [--streams 3] [--only 1996:3:0.6]

Author: Hannah A. Henry, UNC CECL
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

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[4]
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
# The first grid (dipole, 1996->2003 linear ramp: the schedule the runner had
# until 2026-09-29). Later studies pass --study and --kind.
STUDY = "groin/2026-09-29-option-a-grid"
KIND = "dipole"
RAW = REPO / "output" / "raw_runs" / "experiments"
LOGS = HERE / "grid_logs"

M_VALUES = (1, 2, 3, 4, 5, 7, 10)
F_VALUES = (0.2, 0.4, 0.6, 0.8, 1.0)
PERIODS = (1996, 2010)


def member(M, f):
    if KIND == "blocking":
        return f"b{M:.2f}_f{f:.1f}"
    return f"M{M}_f{f:.1f}"


def done(period, M, f):
    return any((RAW / STUDY / member(M, f)).glob(
        f"{period}_*/*/*/*_run_metadata.json"))


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
        # Low M first in both windows, so an early stop still leaves the
        # plausible end of the grid.
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
