#!/usr/bin/env python3
"""End-to-end unattended driver: comparison matrix, groin sweeps, joint fit.

Runs the whole thing in one command, in the order the dependencies require,
and can be re-invoked at any point to pick up where it stopped.

THE ORDER, AND WHY IT IS THIS ORDER
    1  archive        The retired run_index.csv is moved aside. Its rows
                      describe runs whose directories no longer exist, and two
                      of them disagree with each other about the background
                      erosion they were fit under, so merging new runs into it
                      would produce an index of mixed provenance.

    2  matrix nogroin 18 runs: the scenarios that are DISTINCT in each period
                      x {zeroBE, edgeBE} x 2 periods, groin off. These are the
                      comparison baselines every groin result is read against.
                      They need no M or f, so they can run before anything is
                      fitted -- and running them first also satisfies the rule
                      that a scenario's no-groin run must exist before its
                      groin run, which is how section 12.3 resolves its paired
                      baseline.

                      18 and not 20: 1984-2004 has no nourishment scheduled
                      (all three projects are 2014 or later), so
                      `full_no_fill` differs from `full_management` in a
                      switch with nothing to act on. The two build the same
                      modules and derive the same RUN_NAME, and the runner's
                      collision guard refuses the second. `scenario_applies`
                      skips those two cells rather than forcing a name onto
                      what would be a duplicate run. The fills axis exists
                      only in period 2.

    3  seed           2 runs: edgeBE / full_management / groin on, one per
                      period, at provisional M and f. These exist only so the
                      sweep has something to validate its duplicated model
                      code against -- the drift guard differences a rate curve
                      against a published matrix run, and there is none on
                      disk. Stage 6 re-runs both at the fitted values.

    4  sweeps         344 cells: 4 sweeps (2 periods x 2 presets). Period 1
                      edgeBE carries the be1 axis and is 215 of them.

    5  joint fit      Intersects the two periods' surfaces per preset. M and f
                      are not separable within one period, so this is where
                      the parameters are actually chosen.

    6  matrix groin   18 runs: the same 18 scenario cells with the groin on,
                      at the fitted (M, f) for their preset.

RESUME
    Every finished job is appended to `driver_manifest.jsonl` the moment it
    lands, keyed by (stage, period, preset, scenario, groin). Re-invoking
    skips anything already recorded complete and retries anything recorded
    failed. The manifest -- not the run directory -- is the source of truth,
    so this file never has to reimplement the runner's RUN_NAME construction,
    which is derived in the runner's section 7.5 from what its modules
    actually built and would drift if copied.

CONCURRENCY
    Matrix runs are SERIAL. Each one appends a row to run_index.csv, and
    concurrent appends to one CSV interleave. The sweep is parallel -- its
    cells write only to their own directories and are collected through a
    single thread -- and it is where nearly all the wall clock goes.

Usage:
    python HAT_run_all.py [--workers N] [--dry-run] [--no-model-state]
                          [--stages 2,4,5]

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline import nourishment  # noqa: E402
from hatteras_site_config import (  # noqa: E402
    HATTERAS_DOMAINS,
    HATTERAS_NOURISHMENT_PROJECTS,
)

from HAT_groin_sweep_config import END_YEAR, PERIODS, PRESETS  # noqa: E402

HINDCAST = _HERE.parent / "HAT_hindcast_1984_2024.py"
SWEEP = _HERE.parent / "HAT_groin_sweep.py"
JOINT_FIT = _HERE.parent / "HAT_groin_joint_fit.py"

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
SWEEP_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep"
JOINT_JSON = SWEEP_DIR / "joint_fit.json"

DRIVER_DIR = PROJECT_BASE_DIR / "output" / "driver"
MANIFEST = DRIVER_DIR / "driver_manifest.jsonl"
LOG_DIR = DRIVER_DIR / "logs"

SCENARIOS = ("natural", "roadway_only", "beachdune_only",
             "full_management", "full_no_fill")


def period_has_fill(period):
    """Whether any nourishment project falls inside a period.

    Resolved with the same `build_schedule` call the runner makes, rather than
    by re-implementing the date filter, so this cannot disagree with what
    section 6 actually builds.
    """
    return bool(nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
        period, END_YEAR[period]).projects)


def scenario_applies(period, scenario):
    """Whether a scenario is a DISTINCT run in this period.

    `full_no_fill` exists to isolate the nourishment fill against
    `full_management`. In a period with no fill scheduled, the two differ in a
    switch that has nothing to act on: they build the same modules, produce
    the same result, and -- correctly -- derive the same RUN_NAME, so the
    runner's collision guard refuses the second one.

    1984-2004 has no fill (all three projects are 2014 or later), so the
    fills axis simply does not exist there. Skipping is the honest handling:
    running it under a forced name would file two identical runs as if they
    were a contrast.

    Returns:
        (applies, reason). reason is None when it applies.
    """
    if scenario == "full_no_fill" and not period_has_fill(period):
        return False, (f"{period}-{END_YEAR[period]} has no nourishment "
                       f"scheduled, so full_no_fill is identical to "
                       f"full_management")
    return True, None

# The seed runs exist only to give the drift guard a reference. Their M and f
# are provisional and are NOT a prior: stage 6 re-runs both cells at the
# fitted values, so nothing downstream ever reads these numbers.
SEED_M, SEED_F = 50.0, 0.9
SEED_PRESET, SEED_SCENARIO = "edgeBE", "full_management"

# A 20-year run is a few minutes; an hour is a wide margin that still stops a
# wedged run from holding an overnight job open until morning.
RUN_TIMEOUT_S = 3600
SWEEP_TIMEOUT_S = 12 * 3600


# =============================================================================
# MANIFEST
# =============================================================================

def job_key(stage, period=None, preset=None, scenario=None, groin=None):
    """Stable identity for one unit of work."""
    return "|".join(str(x) for x in
                    (stage, period, preset, scenario, groin))


def load_manifest():
    """Returns {job_key: row} for everything recorded so far."""
    if not MANIFEST.exists():
        return {}
    done = {}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            done[row["key"]] = row
    return done


def record(row):
    """Appends one job outcome immediately, so an interrupt costs one job."""
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")


def already_done(manifest, key):
    """True only if the job is recorded AND recorded as having succeeded."""
    row = manifest.get(key)
    return bool(row) and row.get("ok") is True


# =============================================================================
# RUNNING ONE HINDCAST
# =============================================================================

def run_hindcast(period, preset, scenario, groin, M, fraction, overwrite,
                 save_state, log_name, dry_run):
    """Runs the hindcast once, driven entirely through the environment.

    The runner reads these through HAT_hindcast_config, so no tracked source
    is edited and the notebook stays byte-identical to what a driven run uses.

    Returns:
        (ok, seconds, detail).
    """
    env = os.environ.copy()
    env.update({
        "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": preset,
        "HAT_SCENARIO": scenario,
        "HAT_GROIN_ENABLED": "true" if groin else "false",
        "HAT_GROIN_TRAPPING_RATE_M_YR": str(M),
        "HAT_GROIN_DETERIORATION_FRACTION": str(fraction),
        "HAT_OVERWRITE": "true" if overwrite else "false",
        "HAT_SAVE_MODEL_STATE": "true" if save_state else "false",
        # Matplotlib must not try to open a window in an unattended run.
        "MPLBACKEND": "Agg",
    })

    if dry_run:
        print(f"      would run: period={period} preset={preset} "
              f"scenario={scenario} groin={groin} M={M} f={fraction}")
        return True, 0.0, "dry-run"

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{log_name}.log"
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(HINDCAST)],
            env=env, cwd=str(PROJECT_BASE_DIR),
            capture_output=True, text=True, timeout=RUN_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return False, time.perf_counter() - t0, f"timeout after {RUN_TIMEOUT_S}s"

    seconds = time.perf_counter() - t0
    log_path.write_text((proc.stdout or "") + "\n--- STDERR ---\n"
                        + (proc.stderr or ""), encoding="utf-8")
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or "").strip().splitlines()[-4:])
        return False, seconds, f"exit {proc.returncode}: {tail}"
    return True, seconds, str(log_path)


def matrix_stage(stage, groin, manifest, args, fits=None):
    """Runs the 20 matrix cells for one groin state, serially.

    Args:
        stage: "matrix_nogroin" or "matrix_groin".
        groin: Whether the groin is attached.
        manifest: Loaded manifest, for resume.
        args: Parsed CLI arguments.
        fits: {preset: {"M": ..., "fraction": ...}} for the groin stage.

    Returns:
        (completed, failed) counts.
    """
    jobs, skipped = [], []
    for period in PERIODS:
        for preset in PRESETS:
            for scenario in SCENARIOS:
                applies, reason = scenario_applies(period, scenario)
                if applies:
                    jobs.append((period, preset, scenario))
                else:
                    skipped.append((period, preset, scenario, reason))

    completed = failed = 0
    print(f"\n{'=' * 72}\nSTAGE {stage}: {len(jobs)} runs (serial)\n{'=' * 72}")
    if skipped:
        print(f"  {len(skipped)} cell(s) skipped as degenerate:")
        for period, preset, scenario, reason in skipped:
            print(f"    {period} {preset:<7} {scenario:<16} -- {reason}")
            key = job_key(stage, period, preset, scenario, groin)
            if not already_done(manifest, key) and not args.dry_run:
                # Recorded as ok so a re-invocation does not retry it, with
                # skipped=True so the summary can tell "did not need running"
                # from "ran and succeeded".
                record(dict(key=key, stage=stage, period=period,
                            preset=preset, scenario=scenario, groin=groin,
                            ok=True, skipped=True, detail=reason,
                            at=datetime.now().isoformat(timespec="seconds")))

    for index, (period, preset, scenario) in enumerate(jobs, start=1):
        key = job_key(stage, period, preset, scenario, groin)
        if already_done(manifest, key):
            print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
                  f"{scenario:<16} groin={groin}  SKIP (done)")
            completed += 1
            continue

        M = fraction = 0.0
        overwrite = False
        if groin:
            fit = (fits or {}).get(preset)
            if not fit:
                print(f"  [{index:>2}/{len(jobs)}] no fitted (M, f) for "
                      f"{preset}; skipping its groin runs")
                failed += 1
                continue
            M, fraction = fit["M"], fit["fraction"]
            # The two seed cells already hold a run at provisional values.
            # Only those are overwritten; the other 18 keep the collision
            # guard that stops a matrix run from silently replacing another.
            overwrite = (preset == SEED_PRESET and scenario == SEED_SCENARIO)

        label = f"{period}_{preset}_{scenario}_{'groin' if groin else 'nogroin'}"
        print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
              f"{scenario:<16} groin={groin}"
              + (f" M={M:g} f={fraction:g}" if groin else ""), flush=True)

        ok, seconds, detail = run_hindcast(
            period, preset, scenario, groin, M, fraction, overwrite,
            not args.no_model_state, label, args.dry_run)
        record(dict(key=key, stage=stage, period=period, preset=preset,
                    scenario=scenario, groin=groin, M=M, fraction=fraction,
                    ok=ok, seconds=round(seconds, 1), detail=detail,
                    at=datetime.now().isoformat(timespec="seconds")))
        if ok:
            completed += 1
            print(f"        done in {seconds / 60:.1f} min")
        else:
            failed += 1
            print(f"        FAILED: {detail}")

    return completed, failed


# =============================================================================
# STAGES
# =============================================================================

def stage_archive(manifest, args):
    """Moves the retired run_index.csv aside, once."""
    key = job_key("archive")
    if already_done(manifest, key):
        print("\nSTAGE archive: SKIP (done)")
        return
    print(f"\n{'=' * 72}\nSTAGE archive\n{'=' * 72}")
    if not RUN_INDEX.exists():
        print("  no run_index.csv to archive")
    elif args.dry_run:
        print(f"  would archive {RUN_INDEX}")
        return
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        target = RAW_RUNS / f"run_index_archive_{stamp}.csv"
        shutil.move(str(RUN_INDEX), str(target))
        print(f"  archived -> {target.name}")
    if not args.dry_run:
        record(dict(key=key, stage="archive", ok=True,
                    at=datetime.now().isoformat(timespec="seconds")))


def stage_seed(manifest, args):
    """Runs the per-period drift-guard reference at provisional M and f."""
    print(f"\n{'=' * 72}\nSTAGE seed: {len(PERIODS)} runs\n{'=' * 72}")
    print(f"  provisional M = {SEED_M:g}, f = {SEED_F:g} -- these exist only "
          f"to give the\n  sweep's drift guard a rate curve to difference "
          f"against. Stage 6 re-runs\n  both cells at the fitted values.")
    completed = failed = 0
    for period in PERIODS:
        key = job_key("seed", period, SEED_PRESET, SEED_SCENARIO, True)
        if already_done(manifest, key):
            print(f"  {period}: SKIP (done)")
            completed += 1
            continue
        print(f"  {period} {SEED_PRESET} {SEED_SCENARIO} groin=True", flush=True)
        ok, seconds, detail = run_hindcast(
            period, SEED_PRESET, SEED_SCENARIO, True, SEED_M, SEED_F,
            False, not args.no_model_state, f"seed_{period}", args.dry_run)
        record(dict(key=key, stage="seed", period=period, preset=SEED_PRESET,
                    scenario=SEED_SCENARIO, groin=True, M=SEED_M,
                    fraction=SEED_F, ok=ok, seconds=round(seconds, 1),
                    detail=detail,
                    at=datetime.now().isoformat(timespec="seconds")))
        if ok:
            completed += 1
            print(f"    done in {seconds / 60:.1f} min")
        else:
            failed += 1
            print(f"    FAILED: {detail}")
    return completed, failed


def stage_sweeps(manifest, args):
    """Runs the four period/preset sweeps."""
    print(f"\n{'=' * 72}\nSTAGE sweeps: 4 sweeps, 344 cells\n{'=' * 72}")
    completed = failed = 0
    # Period 1 edgeBE is 215 of the 344 cells; running it first means the
    # longest leg starts earliest and a morning check finds it either done or
    # visibly progressing, rather than not yet begun.
    order = [(1984, "edgeBE"), (1984, "zeroBE"),
             (2004, "edgeBE"), (2004, "zeroBE")]
    for period, preset in order:
        key = job_key("sweep", period, preset)
        if already_done(manifest, key):
            print(f"  {period} {preset}: SKIP (done)")
            completed += 1
            continue
        command = [sys.executable, str(SWEEP), "--period", str(period),
                   "--preset", preset, "--workers", str(args.workers)]
        if args.dry_run:
            command.append("--dry-run")
        print(f"\n  --> {' '.join(command[1:])}", flush=True)

        LOG_DIR.mkdir(parents=True, exist_ok=True)
        log_path = LOG_DIR / f"sweep_{period}_{preset}.log"
        t0 = time.perf_counter()
        # Streamed to a log rather than captured: a 215-cell sweep prints a
        # progress line per cell, and holding four hours of that in memory to
        # write it at the end would lose all of it if the driver were killed.
        with log_path.open("w", encoding="utf-8") as handle:
            proc = subprocess.run(command, cwd=str(PROJECT_BASE_DIR),
                                  stdout=handle, stderr=subprocess.STDOUT,
                                  timeout=SWEEP_TIMEOUT_S)
        seconds = time.perf_counter() - t0
        ok = proc.returncode == 0
        record(dict(key=key, stage="sweep", period=period, preset=preset,
                    ok=ok, seconds=round(seconds, 1),
                    detail=f"exit {proc.returncode}; log {log_path}",
                    at=datetime.now().isoformat(timespec="seconds")))
        if ok:
            completed += 1
            print(f"      done in {seconds / 60:.1f} min  ({log_path.name})")
        else:
            failed += 1
            print(f"      FAILED exit {proc.returncode} -- see {log_path}")
    return completed, failed


def stage_joint_fit(args):
    """Runs the joint fit and returns the fitted values per preset."""
    print(f"\n{'=' * 72}\nSTAGE joint fit\n{'=' * 72}")
    if args.dry_run:
        print("  would run HAT_groin_joint_fit.py")
        return {}
    proc = subprocess.run([sys.executable, str(JOINT_FIT)],
                          cwd=str(PROJECT_BASE_DIR),
                          capture_output=True, text=True)
    print(proc.stdout)
    if proc.returncode != 0 or not JOINT_JSON.exists():
        print(f"  joint fit produced no result (exit {proc.returncode}); "
              f"stage 6 cannot run")
        return {}
    return json.loads(JOINT_JSON.read_text())


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--workers", type=int, default=8,
                        help="sweep pool width (default 8; RAM-capped)")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan without running anything")
    parser.add_argument("--no-model-state", action="store_true",
                        help="skip the ~160 MB .npz per matrix run "
                             "(saves ~6 GB; figures then need a re-run)")
    parser.add_argument("--stages", default="1,2,3,4,5,6",
                        help="comma-separated stage numbers to run "
                             "(default all)")
    args = parser.parse_args()

    stages = {int(s) for s in args.stages.split(",") if s.strip()}
    started = datetime.now()
    print("=" * 72)
    print(f"HATTERAS FULL RUN  started {started:%Y-%m-%d %H:%M:%S}")
    print("=" * 72)
    print(f"  stages       {sorted(stages)}")
    print(f"  sweep pool   {args.workers}")
    print(f"  model state  {'NOT saved' if args.no_model_state else 'saved'}")
    print(f"  manifest     {MANIFEST}")
    print(f"  logs         {LOG_DIR}")

    manifest = load_manifest()
    totals = {}

    if 1 in stages:
        stage_archive(manifest, args)
    if 2 in stages:
        totals["matrix nogroin"] = matrix_stage(
            "matrix_nogroin", False, manifest, args)
    if 3 in stages:
        totals["seed"] = stage_seed(manifest, args)
    if 4 in stages:
        manifest = load_manifest()
        totals["sweeps"] = stage_sweeps(manifest, args)

    fits = stage_joint_fit(args) if 5 in stages else {}

    if 6 in stages:
        manifest = load_manifest()
        totals["matrix groin"] = matrix_stage(
            "matrix_groin", True, manifest, args, fits=fits)

    elapsed = (datetime.now() - started).total_seconds() / 3600
    print("\n" + "=" * 72)
    print(f"FINISHED  {elapsed:.2f} h")
    print("=" * 72)
    for name, (done, failed) in totals.items():
        flag = "" if not failed else f"   {failed} FAILED"
        print(f"  {name:<16} {done} completed{flag}")
    if fits:
        print("\n  fitted parameters")
        for preset, fit in fits.items():
            bound = (f"   RAILED on {', '.join(fit['at_grid_bound'])}"
                     if fit.get("at_grid_bound") else "")
            print(f"    {preset:<8} M = {fit['M']:g} m/yr, "
                  f"f = {fit['fraction']:g}{bound}")
    print(f"\n  re-invoke this script to retry anything that failed")
    return 0 if not any(f for _, f in totals.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
