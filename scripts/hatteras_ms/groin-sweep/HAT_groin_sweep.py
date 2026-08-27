#!/usr/bin/env python3
"""Groin / background-erosion grid search for ONE period and preset.

Sweeps the groin trapping rate M against the deterioration floor f -- and, in
period 1 under edgeBE, against the GIS 1 background-erosion rate be1 -- scoring
each cell on that period's CoastSat LRR over D1-D12. Every reported number
comes from a real CASCADE run; no surrogate, no interpolation.

This script fits ONE period. M and f are not separable within a single period
(see HAT_groin_sweep_config's JOINT IDENTIFIABILITY note), so the actual
parameter choice is made by `HAT_groin_joint_fit.py`, which intersects the two
periods' surfaces. What this script produces is one period's surface.

WHY M, f AND be1 TOGETHER
    The groin sits at GIS 5.5, four domains from the modelled reach's southern
    boundary. Its downdrift sink does not stay local: at M = 70 it imposes
    roughly -0.5 m/yr at D1 and -0.7 m/yr at D5. The be1 edge source pushes on
    exactly the same domains, so fitting either with the other held fixed
    charges the same erosion to whichever knob happens to be free. f enters
    through the cumulative trapping the schedule actually delivers. All three
    are swept together because none of them is separable from the others.

WHAT IS RANKED, AND WHAT IS ONLY REPORTED
    Ranked:   |differential - observed|, where differential is the modelled
              D6 - D5 rate. It is the only metric that identifies M. Over
              D1-D12 the profile RMSE moves 7% while M moves 4x; the
              differential moves 5x over the same span.
    Reported: RMSE over the D5/D6 pair, RMSE and bias over D1-D12, the
              per-domain modelled rates, and the fillet extent.
    Never fit: the extent. `cascade.groin`'s design argument is that M sets
              amplitude only and the alongshore extent falls out for free, so
              the emergent extent is the module's one independent test. It is
              measured against the paired M = 0 run at the same be1 and
              written to the CSV, and nothing ranks on it.

WHAT THIS SWEEP CANNOT DO, STATED UP FRONT
    Period 1's residual over D1-D12 does not have the shape of a single-domain
    edge source. At the joint least-squares optimum it still runs -1.7 m/yr at
    D6 and +1.7 m/yr at D11: the model under-erodes the Cape Point reach by
    2-3 m/yr, be1's influence decays from 0.121 to 0.013 m/yr per unit across
    the window, and no value of be1 has the reach to flatten that. Expect a
    floor around 1.0 m/yr RMSE. That floor is the finding -- it is what
    `calibBE`'s Cape Point entries were invented to absorb, and holding to
    strict edgeBE leaves it visible instead of hiding it in the interior.

    Period 2's observed D6 - D5 is NEGATIVE (-2.47 m/yr): the updrift domain
    eroded faster than the downdrift one. The source/sink pair cannot produce
    that at any M >= 0, so the period-2 leg reports a BOUND, not an optimum,
    and says so in its summary. Do not read its best cell as a fitted value.

    The model restarts its dipole fresh at the start of each period, from an
    observed surface that already carries 15 years (1984) or 35 years (2004)
    of accumulated fillet. Treating M as one structure-level parameter across
    both windows -- which the joint fit does -- is a modelling decision, not
    something these runs establish.

DESIGN
    One subprocess per combination. An earlier in-process sweep died partway
    through with a Windows access violation (0xC0000005) from state
    accumulating across many Cascade constructions; process isolation makes
    the OS reclaim everything between combinations.

    Resumable. Combinations with a scored result are skipped; rows recorded
    as failed are retried automatically on the next invocation. Each result
    is appended as a JSON line the moment it lands, so an interrupted sweep
    costs at most one run.

    Parallel. NUM_CORES is 1 in the worker (>1 has crashed on this
    configuration) and CASCADE's internal joblib.Parallel therefore does not
    fan out, so concurrent workers do not oversubscribe. The pool is capped if
    the requested width will not fit in available RAM.

    Self-validating. Before reporting anything, the worker's duplicated copy
    of build_cascade / run_cascade_simulation is run against the period's
    published matrix run and the two rate curves differenced. The validation
    cell lives in a per-period `_validation/` directory and is shared by both
    presets, because what it checks is the CODE, not the preset.

Usage:
    python HAT_groin_sweep.py --period 1984 --preset edgeBE [--workers N]
                              [--dry-run] [--skip-validation]

Writes to output/groin_sweep/<start>_<end>_<preset>/:
    sweep_results.jsonl      one JSON line per combination, written as it lands
    sweep_results.csv        the same rows as a table, plus extent, at the end
    <combo>/                 shoreline matrix + rate curve per combination
Nothing is written to output/raw_runs/ or run_index.csv: a sweep combination
is not a run of the scenario matrix and must not be filed as one.

Author: Hannah A. Henry, UNC CECL
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
# parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
# The guard below is what makes a future move fail here, loudly, instead of
# resolving to scripts/scripts and surfacing as a missing data file several
# imports deeper.
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

from HAT_groin_sweep_config import (  # noqa: E402
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

WORKER = _HERE.parent / "HAT_groin_sweep_worker.py"

# A 20-year run is ~1.5 min; 10 minutes is a wide margin for a loaded machine
# and stops one wedged worker from holding the pool open indefinitely.
WORKER_TIMEOUT_S = 900


# =============================================================================
# RESULT LOG
# =============================================================================

def load_results(results_jsonl):
    """Loads every recorded result, or an empty frame if none exist.

    Returns:
        A DataFrame of previous results with guaranteed `combo` and
        `differential_err` columns. Failed rows carry NaN in
        `differential_err` and a message in `error`.
    """
    if not results_jsonl.exists():
        return pd.DataFrame(columns=["combo", "differential_err"])
    rows = [json.loads(line) for line in
            results_jsonl.read_text().splitlines() if line.strip()]
    frame = pd.DataFrame(rows)
    # A sweep with only failures has no differential_err column at all, and
    # every caller filters on it.
    for column in ("combo", "differential_err"):
        if column not in frame:
            frame[column] = np.nan
    return frame


def append_result(results_jsonl, row):
    """Records one result immediately, as a JSON line.

    JSONL rather than appending to the CSV: a failure record has far fewer
    keys than a success record, and `to_csv(mode="a")` writes values
    positionally under the existing header -- so one failure after a run of
    successes silently shifts every later row's columns. JSON lines carry
    their own keys, so a mixed-schema log stays readable and resumable.
    """
    results_jsonl.parent.mkdir(parents=True, exist_ok=True)
    with results_jsonl.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")


# =============================================================================
# EXECUTION
# =============================================================================

def run_worker(period, preset, combo, out_root):
    """Runs one combination in its own subprocess.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".
        combo: An (M, be1, fraction) tuple.
        out_root: Directory the combination's subdirectory goes under.

    Returns:
        A result dict on success, or a failure record carrying `error` and
        the worker's last stderr lines.
    """
    M, be1, fraction = combo
    name = combo_dir_name(M, be1, fraction)
    out_dir = Path(out_root) / name
    failure = dict(M=M, be1=be1, fraction=fraction, combo=name,
                   period=period, preset=preset)

    # One thread per worker. numpy/BLAS defaults to one thread per core, so N
    # concurrent workers each spawn N threads and the pool spends its time
    # context-switching instead of running: measured, four workers took 4-5x
    # longer per cell than one, for almost no net throughput. CASCADE's own
    # NUM_CORES is already 1, so nothing here wants the extra threads.
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

    # No result line: an access violation exits non-zero with nothing on
    # stdout, so the stderr tail is the only diagnostic there is.
    tail = "\n".join((proc.stderr or "").strip().splitlines()[-6:])
    return dict(failure, error=f"exit {proc.returncode}", stderr_tail=tail)


def safe_worker_count(requested):
    """Caps the pool width at what available RAM will hold.

    Each worker builds 120 Barrier3D domains and holds the run's full state.
    Oversubscribing memory does not fail cleanly -- it swaps, and a swapping
    sweep is slower than a serial one.

    Returns:
        The width to actually use. Falls back to `requested` if psutil is not
        installed, since a missing optional dependency should not block a
        sweep the user explicitly sized.
    """
    try:
        import psutil
    except ImportError:
        print(f"  psutil not installed -- cannot check RAM headroom; "
              f"using {requested} workers as requested")
        return requested

    available_gb = psutil.virtual_memory().available / 1e9
    # 1.8 GB/worker, measured from a 120-domain 20-year run, plus 2 GB left
    # for the OS and this process.
    fits = max(1, int((available_gb - 2.0) / 1.8))
    if fits < requested:
        print(f"  RAM headroom {available_gb:.1f} GB -> capping pool at "
              f"{fits} workers (requested {requested})")
        return fits
    print(f"  RAM headroom {available_gb:.1f} GB -> {requested} workers fit")
    return requested


def run_pool(period, preset, combos, workers, out_root, results_jsonl, label):
    """Runs a list of combinations concurrently, recording each as it lands."""
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
                # The worker writes nothing to stdout when it dies, so
                # its stderr tail is the only record of WHY. Printed on
                # the first failure of a sweep rather than only into the
                # results file: a sweep where every cell fails otherwise
                # produces a log of 215 identical "exit 1" lines and no
                # diagnosis, which is how four sweeps came to be broken
                # without anyone being able to say what broke them.
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


# =============================================================================
# DRIFT GUARD
# =============================================================================

def read_reference_config(period):
    """Reads the reference matrix run's own M, f and be1 from its metadata.

    Not pinned in the config module on purpose: the reference run is re-run
    with new values as the fit is refined, and a pinned pair would report
    drift the first time it changed.

    Returns:
        ((M, be1, fraction), None) on success, or (None, message) if the
        reference is missing, mid-run, or not comparable to what the worker
        builds.
    """
    run_dir, problem = validation_run_dir(period)
    if problem:
        return None, problem
    name = run_dir.name
    meta_path = run_dir / f"{name}_run_metadata.json"
    rate_csv = run_dir / f"{name}_shoreline_change_rate.csv"
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

    # The metadata stores the deterioration as prose, so the floor has to be
    # parsed back out. An unreadable string is reported rather than defaulted:
    # guessing the floor would validate the worker against a groin schedule
    # the reference did not run.
    text = str(groin.get("deterioration", ""))
    match = re.search(r"floor\s+([0-9.]+)", text)
    if "linear_ramp" not in text or match is None:
        return None, (
            f"cannot read the deterioration floor from the reference run's "
            f"metadata ({text!r}); expected 'linear_ramp, floor <f>'.")

    return (float(groin["trapping_rate_m_yr"]),
            float(source_sink["rate_gis1_m_yr"]),
            float(match.group(1))), None


def validate_against_matrix_run(period, workers):
    """Checks the worker reproduces the period's reference matrix run exactly.

    The worker holds a third copy of `build_cascade` and
    `run_cascade_simulation` -- the notebook has one, the hindcast runner has
    another. Copies drift, and drift here is silent: the sweep would keep
    producing plausible numbers from a model that no longer matches the
    hindcast. Running the reference's own configuration and differencing the
    rate curves catches that numerically.

    The validation cell lives in a per-period `_validation/` directory shared
    by both presets, and always runs under edgeBE regardless of which preset
    is being swept: what it checks is the duplicated CODE, not the preset.

    Returns:
        (ok, message). ok is False if the sweep must not proceed.
    """
    combo, problem = read_reference_config(period)
    if problem:
        return False, problem

    run_dir, problem = validation_run_dir(period)
    if problem:
        return False, problem
    name = run_dir.name
    out_root = (PROJECT_BASE_DIR / "output" / "groin_sweep"
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
    published = pd.read_csv(run_dir / f"{name}_shoreline_change_rate.csv")
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
            f"DRIFT: the sweep's copy of build_cascade / "
            f"run_cascade_simulation no longer reproduces {name}.\n"
            f"  max |difference| {worst:.6g} m/yr at D{offender:.0f} "
            f"(tolerance {VALIDATION_TOLERANCE_M_YR:g})\n"
            f"  Re-sync HAT_groin_sweep_worker.py section 3 with "
            f"HAT_hindcast_1984_2024.py before trusting any sweep result.")

    return True, (f"reproduces {name} to {worst:.2g} m/yr "
                  f"(tolerance {VALIDATION_TOLERANCE_M_YR:g})")


# =============================================================================
# EXTENT, MEASURED AGAINST THE PAIRED BASELINE
# =============================================================================

def attach_fillets(frame, out_root, period):
    """Adds the fillet size and its error to every M > 0 row.

    THIS IS THE RANKING METRIC, and it is computed here rather than in the
    worker for the same reason the extent is: it needs the paired M = 0 run,
    and pairing is a property of the grid, not of one combination. Each cell
    pairs with the baseline at the SAME be1 -- pairing across be1 would report
    the background-erosion difference as fillet.

    Args:
        frame: Scored results, one row per combination.
        out_root: Directory holding the per-combination subdirectories.
        period: 1984 or 2004, selecting the observed fillet to score against.

    Returns:
        A copy of `frame` with `fillet_m` and `fillet_err` columns. Both are
        NaN on the M = 0 baseline rows, which have no fillet by definition.
    """
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


def attach_extents(frame, out_root):
    """Adds the emergent fillet extent to every M > 0 row.

    Each groin cell is paired with the M = 0 cell at the SAME be1. Pairing
    across be1 would report the background-erosion difference as fillet.
    Nothing here feeds the ranking -- the extent is the pre-registered
    independent check from `cascade.groin`.
    """
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


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
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

    # --- resume ---------------------------------------------------------------
    previous = load_results(results_jsonl)
    scored = set(previous.loc[previous["differential_err"].notna(), "combo"])
    todo = [c for c in grid if combo_dir_name(*c) not in scored]
    print(f"  already scored   {len(scored)}")
    print(f"  to run           {len(todo)}")

    workers = safe_worker_count(args.workers)

    # --- drift guard ----------------------------------------------------------
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

    # --- run ------------------------------------------------------------------
    if todo:
        run_pool(period, preset, todo, workers, out_root, results_jsonl,
                 f"{period} {preset}")

    # --- collate --------------------------------------------------------------
    frame = load_results(results_jsonl)
    scored_frame = frame[frame["differential_err"].notna()].copy()
    if scored_frame.empty:
        print("\nno scored combinations; nothing to report")
        return 1

    # A retried combination appears twice in the JSONL. The last line wins:
    # it is the one whose files are on disk.
    scored_frame = scored_frame.drop_duplicates(subset="combo", keep="last")
    scored_frame = attach_extents(scored_frame, out_root)
    scored_frame = attach_fillets(scored_frame, out_root, period)
    # RANKED ON FILLET SIZE. differential_err is retained and reported,
    # but it scores the fillet's SLOPE, which is near zero once the
    # fillet saturates and therefore nearly uninformative about M. Falls
    # back to the differential only if no fillet could be measured --
    # which means the M = 0 baselines are missing, and the sweep should
    # be re-run rather than ranked on the weaker metric silently.
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
