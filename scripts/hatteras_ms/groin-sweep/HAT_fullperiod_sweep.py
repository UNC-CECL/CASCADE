#!/usr/bin/env python3
"""
Groin sweep over one continuous 1984-2024 window, scored on the change profile.

    python scripts/hatteras_ms/groin-sweep/HAT_fullperiod_sweep.py --workers 6 --stage coarse

The 1967-rig method on the production geometry and the full span; one
worker process per cell. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
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
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_groin_sweep_config import GROIN_SWEEP_ROOT  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS as GEOMETRY  # noqa: E402

from HAT_fullperiod_target import (  # noqa: E402
    END_YEAR,
    FIT_DOMAINS_GIS,
    START_YEAR,
    model_change_profile,
    observed_change_profile,
    profile_rmse,
)

WORKER = _HERE.parent / "HAT_groin_sweep_worker.py"
# Scoped by wave climate, so another Hs cannot overwrite the 2.5 m results
def _out_root():
    raw = os.environ.get("HAT_SWEEP_HS", "").strip()
    stem = "fullperiod_1984_2024"
    if raw and float(raw) != 2.5:
        stem += "_Hs" + f"{float(raw):g}".replace(".", "p")
    return GROIN_SWEEP_ROOT / stem


# --- CONFIG ------------------------------------------------------------------
OUT_ROOT = _out_root()
RESULTS_CSV = OUT_ROOT / "results.csv"
FIGURE_DIR = OUT_ROOT / "figures"

import sys as _envsys
from pathlib import Path as _EnvP
_envsys.path.insert(0, str(next(_q for _q in _EnvP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_env_forcings as _env  # noqa: E402
STORM_REL = _env.init_relpath(_env.SPLICED_1984_2024)
PRESET = "edgeBE"

# Capped at 80: M >= 100 drowned every rig cell
M_COARSE = [20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0]
F_COARSE = [0.1, 0.3, 0.5, 0.7, 0.9]
FINE_M_STEP, FINE_F_STEP = 10.0, 0.1

WORKER_TIMEOUT_S = 1800     # a 40-year cell is ~6 min; this is a wide margin
# -----------------------------------------------------------------------------


# A cell's name from (M, f)
def combo_name(M, fraction):
    return "M0" if M == 0 else f"M{M:g}_f{fraction:.2f}"


# Runs one combination in its own process
def run_cell(M, fraction, be1):
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


# Profile RMSE for one finished cell, or None if its matrix is absent
def score_cell(name, observed):
    path = OUT_ROOT / name / "shoreline_matrix.npy"
    if not path.exists():
        return None
    model = model_change_profile(np.load(path), GEOMETRY, FIT_DOMAINS_GIS)
    return profile_rmse(model, observed), model


# Runs a list of (M, f) cells in parallel, printing as each lands
def run_grid(cells, be1, workers, label):
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


# Scores every finished cell and writes results.csv, ranked
def collate(observed):
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
    # An empty list gives a DataFrame with no columns
    if not rows:
        return pd.DataFrame(columns=["combo", "M", "fraction", "rmse_m"])
    frame = pd.DataFrame(rows).sort_values("rmse_m").reset_index(drop=True)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    frame.to_csv(RESULTS_CSV, index=False)
    return frame


# Run: the chosen stage's cells, then the results table
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--be1", type=float, default=None,
                        help="solved edge value for the continuous window")
    parser.add_argument("--stage",
                        choices=("baseline", "coarse", "fine", "both"),
                        default="both",
                        help="baseline runs ONLY the M = 0 cell -- the no-groin "
                             "reference the recorded fillet-decay diagnostic "
                             "is measured on, and the cheapest way to test how "
                             "a forcing change moves it")
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
        if args.stage == "baseline":
            run_grid([(0.0, 0.0)], args.be1, args.workers, "baseline")
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
