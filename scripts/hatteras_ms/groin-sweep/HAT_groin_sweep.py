#!/usr/bin/env python3
"""
Groin / background-erosion grid search for one period and preset.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_sweep.py --period 1996 --preset edgeBE --workers 4

Sweeps M against f (and be1 in period 1 under edgeBE), one worker per
cell, ranked on fillet size. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3]: this file is in hatteras_ms/groin-sweep/; the guard makes a move fail loudly here
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline.run_layout import resolve  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_SWEEP_ROOT,
    END_YEAR,
    GROIN_DOWNDRIFT_GIS,
    GROIN_EXTENT_THRESHOLD_FRAC,
    GROIN_UPDRIFT_GIS,
    OBSERVED_DIFFERENTIAL,
    PERIODS,
    PERIOD_DIFFERENTIAL_IS_REACHABLE,
    PRESETS,
    VALIDATION_REQUIRED_PRESET,
    VALIDATION_TOLERANCE_M_YR,
    be_gis90,
    build_grid,
    combo_dir_name,
    measure_groin_extent,
    measure_fillet,
    OBSERVED_FILLET_M,
    sweep_output_dir,
    validation_run_dir,
)

# --- CONFIG ------------------------------------------------------------------
WORKER = _HERE.parent / "HAT_groin_sweep_worker.py"

# A 20-year run is ~1.5 min; the timeout stops a wedged worker
WORKER_TIMEOUT_S = 900
# -----------------------------------------------------------------------------


# Result log

# Loads every recorded result, or an empty frame if none exist
def load_results(results_jsonl):
    if not results_jsonl.exists():
        return pd.DataFrame(columns=["combo", "differential_err"])
    rows = [json.loads(line) for line in
            results_jsonl.read_text().splitlines() if line.strip()]
    frame = pd.DataFrame(rows)
    # A sweep with only failures has no differential_err column at all, and every caller filters on it
    for column in ("combo", "differential_err"):
        if column not in frame:
            frame[column] = np.nan
    return frame


# Records one result immediately, as a JSON line
def append_result(results_jsonl, row):
    results_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with results_jsonl.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")


# Execution

# Runs one combination in its own subprocess
def run_worker(period, preset, combo, out_root):
    M, be1, fraction = combo
    name = combo_dir_name(M, be1, fraction)
    out_dir = Path(out_root) / name
    failure = dict(M=M, be1=be1, fraction=fraction, combo=name,
                   period=period, preset=preset)

    # One thread per worker, or the pool spends its time context-switching
    env = os.environ.copy()
    env.update({name: "1" for name in (
        "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")})
    env["MPLBACKEND"] = "Agg"

    try:
        proc = subprocess.run(
            [sys.executable, str(WORKER), str(period), preset, str(M),
             str(fraction), "none" if be1 is None else str(be1),
             str(out_dir)],
            capture_output=True, text=True, timeout=WORKER_TIMEOUT_S,
            cwd=str(PROJECT_BASE_DIR), env=env,
        )
    except subprocess.TimeoutExpired:
        return dict(failure, error=f"timeout after {WORKER_TIMEOUT_S}s")

    for line in proc.stdout.splitlines():
        if line.startswith("RESULT_JSON="):
            return json.loads(line[len("RESULT_JSON="):])

    # No result line: the stderr tail is the only diagnostic
    tail = "\n".join((proc.stderr or "").strip().splitlines()[-6:])
    return dict(failure, error=f"exit {proc.returncode}", stderr_tail=tail)


# Caps the pool width at what available RAM will hold
def safe_worker_count(requested):
    try:
        import psutil
    except ImportError:
        print(f"  psutil not installed -- cannot check RAM headroom; "
              f"using {requested} workers as requested")
        return requested

    available_gb = psutil.virtual_memory().available / 1e9
    # 1.8 GB a worker (measured on a 120-domain 20-year run), 2 GB kept for the OS
    fits = max(1, int((available_gb - 2.0) / 1.8))
    if fits < requested:
        print(f"  RAM headroom {available_gb:.1f} GB -> capping pool at "
              f"{fits} workers (requested {requested})")
        return fits
    print(f"  RAM headroom {available_gb:.1f} GB -> {requested} workers fit")
    return requested


# Runs a list of combinations concurrently, recording each as it lands
def run_pool(period, preset, combos, workers, out_root, results_jsonl, label):
    results = []
    if not combos:
        return results

    print(f"\n{label}: {len(combos)} combinations, {workers} at a time")
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_worker, period, preset, c, out_root): c
                   for c in combos}
        printed_stderr = False
        for done, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            append_result(results_jsonl, row)
            results.append(row)
            elapsed = time.perf_counter() - t0
            if "error" in row:
                status = f"FAILED  {row['error']}"
                # The worker writes nothing to stdout when it dies, so its stderr tail is the only record of WHY
                if not printed_stderr and row.get("stderr_tail"):
                    printed_stderr = True
                    print()
                    print("  first failure -- worker stderr:")
                    for line in row["stderr_tail"].splitlines():
                        print(f"    | {line}")
                    print()
            else:
                status = (f"diff {row['differential_m_yr']:+6.3f}  "
                          f"err {row['differential_err']:.3f}  "
                          f"RMSE {row['rmse_window']:.3f}")
            rate = elapsed / done
            eta = (len(combos) - done) * rate / 60.0
            print(f"  [{done:>3}/{len(combos)}] {row['combo']:<22} {status}"
                  f"   ({elapsed / 60:.1f} min, ~{eta:.0f} min left)",
                  flush=True)
    return results


# Drift guard

# Reads the reference matrix run's own M, f and be1 from its metadata
def read_reference_config(period):
    run_dir, problem = validation_run_dir(period)
    if problem:
        return None, problem
    name = run_dir.name
    meta_path = resolve(run_dir, "metadata_json", name)
    rate_csv = resolve(run_dir, "rate_csv", name)
    if not meta_path.exists() or not rate_csv.exists():
        return None, (
            f"reference run incomplete:\n    {run_dir}\n"
            f"  The sweep validates its duplicated model code against this "
            f"run. If it is mid-run, wait for it to finish; otherwise re-run "
            f"it, or pass --skip-validation to sweep without the guard.")

    meta = json.loads(meta_path.read_text())
    source_sink = meta.get("source/sink", {})
    groin = meta.get("groin", {})

    if source_sink.get("preset") != VALIDATION_REQUIRED_PRESET:
        return None, (
            f"reference run used preset {source_sink.get('preset')!r}, not "
            f"{VALIDATION_REQUIRED_PRESET!r}; it does not have the source/sink "
            f"shape the worker builds, so a difference would not mean drift.")
    if not groin.get("enabled"):
        return None, ("reference run has no groin, so it cannot validate the "
                      "worker's groin path.")

    be90_ref = float(source_sink.get("rate_gis90_m_yr"))
    if abs(be90_ref - be_gis90(period)) > 1e-9:
        return None, (
            f"reference run used be90 = {be90_ref:+g} m/yr but the site "
            f"config now holds {be_gis90(period):+g}. Re-run the reference "
            f"against the current table before sweeping -- otherwise every "
            f"combination is fit against a different north end.")

    # The metadata stores the deterioration as prose, so the floor has to be parsed back out
    text = str(groin.get("deterioration", ""))
    match = re.search(r"floor\s+([0-9.]+)", text)
    if "linear_ramp" not in text or match is None:
        return None, (
            f"cannot read the deterioration floor from the reference run's "
            f"metadata ({text!r}); expected 'linear_ramp, floor <f>'.")

    return (float(groin["trapping_rate_m_yr"]),
            float(source_sink["rate_gis1_m_yr"]),
            float(match.group(1))), None


# Checks the worker reproduces the period's reference matrix run exactly
def validate_against_matrix_run(period, workers):
    combo, problem = read_reference_config(period)
    if problem:
        return False, problem

    run_dir, problem = validation_run_dir(period)
    if problem:
        return False, problem
    name = run_dir.name
    out_root = (GROIN_SWEEP_ROOT
                / f"_validation_{period}")
    combo_name = combo_dir_name(*combo)
    print(f"  reference {name}: M = {combo[0]:g}, "
          f"be1 = {combo[1]:g}, f = {combo[2]:g}")

    if not (out_root / combo_name / "result.json").exists():
        print(f"  running validation combination {combo_name} ...", flush=True)
        row = run_worker(period, VALIDATION_REQUIRED_PRESET, combo, out_root)
        append_result(out_root / "validation.jsonl", row)
        if "error" in row:
            return False, f"validation combination failed: {row['error']}"

    swept = pd.read_csv(out_root / combo_name / "shoreline_change_rate.csv")
    published = pd.read_csv(resolve(run_dir, "rate_csv", name))
    merged = swept.merge(published, on="gis_domain",
                         suffixes=("_sweep", "_published"))
    if len(merged) != len(published):
        return False, (f"domain mismatch: sweep has {len(swept)} domains, "
                       f"published run has {len(published)}")

    delta = (merged["change_rate_m_yr_sweep"]
             - merged["change_rate_m_yr_published"]).abs()
    worst = float(delta.max())
    if worst > VALIDATION_TOLERANCE_M_YR:
        offender = merged.loc[delta.idxmax(), "gis_domain"]
        return False, (
            f"DRIFT: the worker no longer reproduces {name}.\n"
            f"  max |difference| {worst:.6g} m/yr at D{offender:.0f} "
            f"(tolerance {VALIDATION_TOLERANCE_M_YR:g})\n"
            f"  The worker builds a DIFFERENT MODEL from the runner. Check, "
            f"in this order:\n"
            f"    1. THE TOPOGRAPHY PRODUCT. HAT_groin_sweep_worker.py must "
            f"pass TOPO_PRODUCT to\n"
            f"       both topo_dirs() and build_domain_file_paths(). Omitting "
            f"it selects\n"
            f"       DEFAULT_PRODUCT (2004-start), so a 1984 sweep silently "
            f"builds on the 2004\n"
            f"       island -- all 90 domains differ, 65 in interior SHAPE. "
            f"This was the cause\n"
            f"       on 2026-08-29, at 0.126 m/yr.\n"
            f"    2. The forcing assembled in the worker's section 2 -- "
            f"setbacks, dunes,\n"
            f"       island offset, background erosion -- against the "
            f"reference run's\n"
            f"       run_metadata.json.\n"
            f"    3. Only then the duplicated code. NOTE that build_cascade is "
            f"IMPORTED from\n"
            f"       cascade_pipeline.hindcast, not copied, and "
            f"run_cascade_simulation differs\n"
            f"       from the shared one only in printing. Verified in sync "
            f"2026-08-29 -- this\n"
            f"       message used to send readers here first and it cost three "
            f"dead ends.")

    return True, (f"reproduces {name} to {worst:.2g} m/yr "
                  f"(tolerance {VALIDATION_TOLERANCE_M_YR:g})")


# Extent, measured against the paired baseline

# Adds the fillet size and its error to every M > 0 row
def attach_fillets(frame, out_root, period):
    geometry = HATTERAS_DOMAINS
    observed = OBSERVED_FILLET_M[period]
    baselines = {}
    for _, row in frame[frame["M"] == 0].iterrows():
        path = Path(out_root) / row["combo"] / "shoreline_matrix.npy"
        if path.exists():
            baselines[row.get("be1")] = np.load(path)

    sizes, errors = [], []
    for _, row in frame.iterrows():
        baseline = baselines.get(row.get("be1"))
        path = Path(out_root) / row["combo"] / "shoreline_matrix.npy"
        if row["M"] == 0 or baseline is None or not path.exists():
            sizes.append(np.nan)
            errors.append(np.nan)
            continue
        fillet = measure_fillet(np.load(path), baseline, geometry)
        sizes.append(fillet)
        errors.append(abs(fillet - observed))

    frame = frame.copy()
    frame["fillet_m"] = sizes
    frame["observed_fillet_m"] = observed
    frame["fillet_err"] = errors
    return frame


# Adds the emergent fillet extent to every M > 0 row
def attach_extents(frame, out_root):
    geometry = HATTERAS_DOMAINS
    baselines = {}
    for _, row in frame[frame["M"] == 0].iterrows():
        path = Path(out_root) / row["combo"] / "shoreline_matrix.npy"
        if path.exists():
            baselines[row.get("be1")] = np.load(path)

    up_m, down_m, peak_m = [], [], []
    for _, row in frame.iterrows():
        baseline = baselines.get(row.get("be1"))
        path = Path(out_root) / row["combo"] / "shoreline_matrix.npy"
        if row["M"] == 0 or baseline is None or not path.exists():
            up_m.append(np.nan); down_m.append(np.nan); peak_m.append(np.nan)
            continue
        extent = measure_groin_extent(
            np.load(path), baseline, geometry,
            GROIN_UPDRIFT_GIS, GROIN_DOWNDRIFT_GIS,
            GROIN_EXTENT_THRESHOLD_FRAC)
        up_m.append(extent["updrift_m"])
        down_m.append(extent["downdrift_m"])
        peak_m.append(extent["peak_m"])

    frame = frame.copy()
    frame["extent_updrift_m"] = up_m
    frame["extent_downdrift_m"] = down_m
    frame["extent_peak_m"] = peak_m
    return frame


# Run: build the grid, run the cells, rank and write
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--period", type=int, required=True, choices=PERIODS)
    parser.add_argument("--preset", required=True, choices=PRESETS)
    parser.add_argument("--workers", type=int, default=8,
                        help="concurrent subprocesses (default 8; capped by "
                             "available RAM)")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the grid and exit without running")
    parser.add_argument("--skip-validation", action="store_true",
                        help="sweep without the code-drift guard")
    args = parser.parse_args()

    period, preset = args.period, args.preset
    out_root = sweep_output_dir(period, preset)
    results_jsonl = out_root / "sweep_results.jsonl"
    results_csv = out_root / "sweep_results.csv"
    grid = build_grid(period, preset)

    print("=" * 72)
    print(f"GROIN SWEEP  {period}-{END_YEAR[period]}  {preset}")
    print("=" * 72)
    print(f"  observed D6-D5   {OBSERVED_DIFFERENTIAL[period]:+.2f} m/yr")
    if not PERIOD_DIFFERENTIAL_IS_REACHABLE[period]:
        print("  NOTE: that target is NEGATIVE. The source/sink pair adds -M "
              "updrift and\n"
              "        +M downdrift, so the modelled differential is "
              "non-negative at any\n"
              "        M >= 0 and cannot reach it. This leg reports a BOUND, "
              "not an optimum:\n"
              "        expect the fit to rail to the smallest cumulative "
              "trapping on the grid.")
    print(f"  grid             {len(grid)} combinations")
    print(f"  output           {out_root}")

    if args.dry_run:
        for combo in grid:
            print("   ", combo_dir_name(*combo))
        return 0

    # Resume
    previous = load_results(results_jsonl)
    scored = set(previous.loc[previous["differential_err"].notna(), "combo"])
    todo = [c for c in grid if combo_dir_name(*c) not in scored]
    print(f"  already scored   {len(scored)}")
    print(f"  to run           {len(todo)}")

    workers = safe_worker_count(args.workers)

    # Drift guard
    if args.skip_validation:
        print("\n  VALIDATION SKIPPED -- results are not guarded against "
              "code drift")
    elif todo:
        print("\nvalidating the worker against the published matrix run")
        ok, message = validate_against_matrix_run(period, workers)
        if not ok:
            print(f"\nSWEEP ABORTED\n  {message}")
            return 2
        print(f"  OK: {message}")

    # Run
    if todo:
        run_pool(period, preset, todo, workers, out_root, results_jsonl,
                 f"{period} {preset}")

    # Collate
    frame = load_results(results_jsonl)
    scored_frame = frame[frame["differential_err"].notna()].copy()
    if scored_frame.empty:
        print("\nno scored combinations; nothing to report")
        return 1

    # A retried combination appears twice in the JSONL
    scored_frame = scored_frame.drop_duplicates(subset="combo", keep="last")
    scored_frame = attach_extents(scored_frame, out_root)
    scored_frame = attach_fillets(scored_frame, out_root, period)
    # Ranked on fillet size, not slope (README)
    if "fillet_err" in scored_frame and scored_frame["fillet_err"].notna().any():
        scored_frame = scored_frame.sort_values("fillet_err")
        _rank_metric = "fillet_err"
    else:
        print("  WARNING: no fillet could be measured (missing M = 0 "
              "baselines?); ranking on the differential instead")
        scored_frame = scored_frame.sort_values("differential_err")
        _rank_metric = "differential_err"
    scored_frame.to_csv(results_csv, index=False)

    failed = frame[frame["differential_err"].isna()]
    failed = failed[~failed["combo"].isin(scored_frame["combo"])]

    print("\n" + "=" * 72)
    print(f"  scored   {len(scored_frame)} / {len(grid)}")
    if len(failed):
        print(f"  FAILED   {len(failed)} (re-run this script to retry them)")
        for _, row in failed.head(5).iterrows():
            print(f"    {row['combo']:<22} {row.get('error')}")
        tail = failed.iloc[0].get("stderr_tail")
        if isinstance(tail, str) and tail.strip():
            print(f"    worker stderr ({failed.iloc[0]['combo']}):")
            for line in tail.splitlines():
                print(f"      | {line}")
    print(f"  written  {results_csv}")

    best = scored_frame.iloc[0]
    print(f"\n  best cell by |differential - observed|:")
    print(f"    {best['combo']}   diff {best['differential_m_yr']:+.3f} "
          f"(target {OBSERVED_DIFFERENTIAL[period]:+.3f}), "
          f"err {best['differential_err']:.3f}")
    print(f"    D1-D12 RMSE {best['rmse_window']:.3f} m/yr, "
          f"bias {best['bias_window']:+.3f} m/yr")
    if not PERIOD_DIFFERENTIAL_IS_REACHABLE[period]:
        print("    ^ this is the grid's LOWER BOUND on trapping, not a fitted "
              "optimum")
    print("\n  M and f are not separable within one period -- run "
          "HAT_groin_joint_fit.py\n  once both periods are swept.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
