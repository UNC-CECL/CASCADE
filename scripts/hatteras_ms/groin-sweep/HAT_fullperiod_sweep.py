#!/usr/bin/env python3
"""Groin sweep over ONE continuous 1984-2024 window, scored on the change profile.

This is the 1967-rig method applied to the production geometry and the full
hindcast span. It exists because the two-period, scalar-fillet approach failed
in a specific and repeatable way, and this one did not.

WHAT WENT WRONG BEFORE, AND WHAT IS DIFFERENT HERE

    scalar target -> profile target
        Matching one number (the fillet) left a RIDGE of equally good (M, f)
        pairs; the ranking then returned whichever cell sat at a grid edge.
        Scoring the per-domain CHANGE PROFILE constrains shape as well as
        magnitude. That is what gave the rig an interior optimum in both knobs.

    two 20-year windows -> one 40-year window
        The fillet builds through 1984-2004 (+52 m) and declines through
        2004-2024 (-76 m). Each 20-year window sees only one of those, so M and
        f trade off within it. One continuous run sees both, and they constrain
        different combinations of the pair.

    coarse then fine
        A wide coarse pass, then a fine pass centred on its best cell. M is
        capped at 80: on the rig every cell at M >= 100 DROWNED the barrier
        partway through, and M >= 70 produced RMSE an order of magnitude above
        its neighbours. Sweeping into that region wastes runs on a model that
        refuses the parameter.

ASSUMPTIONS, STATED

    static background erosion
        `background_erosion` is written once at construction, so a 40-year run
        carries ONE field. The two periods' calibrated edge values differ in
        sign (be1 -41.8 vs +50.3), but over the whole window they largely
        cancel: the surveyed change at D1 is -3.2 m in 40 years. be1 is
        therefore SOLVED against that, not inherited from either period.

    spliced storms
        1984-2003 from the 1984-2004 series, 2004-2024 from the 2004-2024 one,
        with the model-year index remapped. Checked at the join: wave height,
        runup and period agree to within 2%.

    this is a rig, not a hindcast
        One static background field cannot reproduce a field that reverses
        mid-window. The run is fit for calibrating a LOCAL structure, where a
        smooth regional field largely cancels in the profile's shape. It is not
        a substitute for the two-period hindcast.

Usage:
    python HAT_fullperiod_sweep.py [--workers 6] [--be1 -5.0] [--stage coarse|fine|both]

Writes output/groin_sweep/fullperiod_1984_2024/:
    results.csv          one row per cell, ranked by profile RMSE
    <combo>/             shoreline matrix per cell
    figures/             heatmap, best-fit profile, top-N profiles

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from HAT_fullperiod_target import (  # noqa: E402
    END_YEAR,
    FIT_DOMAINS_GIS,
    START_YEAR,
    model_change_profile,
    observed_change_profile,
    profile_rmse,
)

WORKER = _HERE.parent / "HAT_groin_sweep_worker.py"
OUT_ROOT = PROJECT_BASE_DIR / "output" / "groin_sweep" / "fullperiod_1984_2024"
RESULTS_CSV = OUT_ROOT / "results.csv"
FIGURE_DIR = OUT_ROOT / "figures"

STORM_REL = "3-env-forcings/storms/hindcast_storms/1984_2024/1984_2024_storms_spliced.npy"
PRESET = "edgeBE"

# Capped at 80: M >= 100 drowned the barrier on every rig cell, M >= 70 went
# unstable there. The production grid is better buffered so the threshold may
# differ, but sweeping far past a known failure mode buys nothing.
M_COARSE = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
F_COARSE = [0.1, 0.3, 0.5, 0.7, 0.9]
FINE_M_STEP, FINE_F_STEP = 10.0, 0.1

WORKER_TIMEOUT_S = 1800     # a 40-year cell is ~6 min; this is a wide margin


def combo_name(M, fraction):
    return "M0" if M == 0 else f"M{M:g}_f{fraction:.2f}"


def run_cell(M, fraction, be1):
    """Runs one combination in its own process. Returns (name, ok, seconds)."""
    name = combo_name(M, fraction)
    out_dir = OUT_ROOT / name
    if (out_dir / "shoreline_matrix.npy").exists():
        return name, True, 0.0
    out_dir.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env["HAT_SWEEP_END_YEAR"] = str(END_YEAR)
    env["HAT_SWEEP_STORM_FILE"] = STORM_REL
    started = time.perf_counter()
    try:
        done = subprocess.run(
            [sys.executable, str(WORKER), str(START_YEAR), PRESET,
             str(M), str(fraction),
             "none" if be1 is None else str(be1), str(out_dir)],
            capture_output=True, text=True, timeout=WORKER_TIMEOUT_S, env=env)
    except subprocess.TimeoutExpired:
        return name, False, WORKER_TIMEOUT_S
    ok = done.returncode == 0 and (out_dir / "shoreline_matrix.npy").exists()
    if not ok:
        (out_dir / "stderr.txt").write_text(done.stderr[-4000:], encoding="utf-8")
    return name, ok, time.perf_counter() - started


def score_cell(name, observed):
    """Profile RMSE for one finished cell, or None if its matrix is absent."""
    path = OUT_ROOT / name / "shoreline_matrix.npy"
    if not path.exists():
        return None
    model = model_change_profile(np.load(path), GEOMETRY, FIT_DOMAINS_GIS)
    return profile_rmse(model, observed), model


def run_grid(cells, be1, workers, label):
    """Runs a list of (M, f) cells in parallel, printing as each lands."""
    todo = [c for c in cells
            if not (OUT_ROOT / combo_name(*c) / "shoreline_matrix.npy").exists()]
    print(f"\n{'=' * 70}\nSTAGE {label}: {len(cells)} cells, {len(todo)} to run, "
          f"{workers} at a time\n{'=' * 70}")
    if not todo:
        return
    done_count = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_cell, M, f, be1): (M, f) for M, f in todo}
        for future in as_completed(futures):
            name, ok, seconds = future.result()
            done_count += 1
            print(f"  [{done_count:>2}/{len(todo)}] {name:<14} "
                  f"{'ok' if ok else 'FAILED':<7} {seconds / 60:5.1f} min",
                  flush=True)


def collate(observed):
    """Scores every finished cell and writes results.csv, ranked."""
    rows = []
    for cell in sorted(OUT_ROOT.iterdir()):
        if not cell.is_dir() or cell.name == "figures":
            continue
        scored = score_cell(cell.name, observed)
        if scored is None:
            continue
        rmse, model = scored
        name = cell.name
        M = 0.0 if name == "M0" else float(name.split("_")[0][1:])
        fraction = 0.0 if name == "M0" else float(name.split("_f")[1])
        rows.append(dict(combo=name, M=M, fraction=fraction, rmse_m=rmse,
                         **{f"change_D{d}": model[d] for d in FIT_DOMAINS_GIS}))
    # An empty list gives a DataFrame with no columns, and sorting on a column
    # that does not exist raises KeyError -- which is what a first invocation,
    # or a --collate-only before anything has run, would hit.
    if not rows:
        return pd.DataFrame(columns=["combo", "M", "fraction", "rmse_m"])
    frame = pd.DataFrame(rows).sort_values("rmse_m").reset_index(drop=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(RESULTS_CSV, index=False)
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--be1", type=float, default=None,
                        help="solved edge value for the continuous window")
    parser.add_argument("--stage", choices=("coarse", "fine", "both"),
                        default="both")
    parser.add_argument("--collate-only", action="store_true")
    args = parser.parse_args()

    observed = observed_change_profile()
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"continuous groin sweep  {START_YEAR}-{END_YEAR}  {PRESET}")
    print(f"  fit window   D{min(FIT_DOMAINS_GIS)}-D{max(FIT_DOMAINS_GIS)}")
    print(f"  be1          {args.be1}")
    print(f"  observed     {min(observed.values()):+.1f} to "
          f"{max(observed.values()):+.1f} m")

    if not args.collate_only:
        if args.stage in ("coarse", "both"):
            cells = [(0.0, 0.0)] + [(M, f) for M in M_COARSE for f in F_COARSE]
            run_grid(cells, args.be1, args.workers, "coarse")

        frame = collate(observed)
        if args.stage in ("fine", "both") and not frame.empty:
            best = frame[frame.M > 0].iloc[0]
            fine = [(M, f)
                    for M in (best.M - FINE_M_STEP, best.M, best.M + FINE_M_STEP)
                    for f in (best.fraction - FINE_F_STEP, best.fraction,
                              best.fraction + FINE_F_STEP)
                    if M > 0 and 0.0 <= f <= 1.0]
            print(f"\n  coarse best {best.combo} (RMSE {best.rmse_m:.1f} m) "
                  f"-> refining around it")
            run_grid(fine, args.be1, args.workers, "fine")

    frame = collate(observed)
    print(f"\n{'=' * 70}")
    print(f"  {len(frame)} cells scored -> {RESULTS_CSV}")
    if not frame.empty:
        print(frame[["combo", "M", "fraction", "rmse_m"]].head(8).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
