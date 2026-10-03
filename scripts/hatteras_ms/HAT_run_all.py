#!/usr/bin/env python3
"""
End-to-end unattended driver: the comparison matrix, the groin sweeps, the joint fit.

    python scripts/hatteras_ms/HAT_run_all.py [--workers N] [--dry-run] [--no-model-state]

Runs every stage in dependency order and resumes where it stopped; each
job recorded in a manifest. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
# The sweep, its config and the joint fit live in groin-sweep/
GROIN_SWEEP_DIR = _HERE.parent / "groin-sweep"
for _path in (SCRIPTS_DIR, _HERE.parent, GROIN_SWEEP_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline import nourishment  # noqa: E402
from cascade_pipeline.roadway import RelocationEvent  # noqa: E402
from cascade_pipeline.run_registry import values_digest  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_BE_PRESETS,
    HATTERAS_DOMAINS,
    HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_ROAD_EVENTS,
)

from HAT_groin_sweep_config import (END_YEAR, GROIN_SWEEP_ROOT, PERIODS,  # noqa: E402
                                    PRESETS)

# --- CONFIG ------------------------------------------------------------------
# The matrix's preset axis is WIDER than the sweep's
ALL_PRESETS = tuple(PRESETS) + tuple(
    name for name in HATTERAS_BE_PRESETS if name not in PRESETS)

HINDCAST = _HERE.parent / "HAT_hindcast_1984_2024.py"
SWEEP = GROIN_SWEEP_DIR / "HAT_groin_sweep.py"
JOINT_FIT = GROIN_SWEEP_DIR / "HAT_groin_joint_fit.py"

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
SWEEP_DIR = GROIN_SWEEP_ROOT
JOINT_JSON = SWEEP_DIR / "joint_fit.json"

PARAMETER_FILE = (PROJECT_BASE_DIR / "data" / "hatteras_init"
                  / "Hatteras-CASCADE-parameters.yaml")

DRIVER_DIR = PROJECT_BASE_DIR / "output" / "logs" / "driver"
MANIFEST = DRIVER_DIR / "driver_manifest.jsonl"
LOG_DIR = DRIVER_DIR / "logs"
LOCK = DRIVER_DIR / "driver.lock"

_TABLE_HINT = (
    "expected a module-level `SCENARIOS = {{...}}` of literal dict(...) calls "
    "in section 3 of {name}; if that changed, this parser changes with it")
# -----------------------------------------------------------------------------


# Reads the runner's SCENARIOS table out of its source, without running it
def _scenario_table():
    tree = ast.parse(HINDCAST.read_text(encoding="utf-8"),
                     filename=str(HINDCAST))
    assign = next(
        (node for node in tree.body
         if isinstance(node, ast.Assign)
         and any(isinstance(target, ast.Name) and target.id == "SCENARIOS"
                 for target in node.targets)),
        None)
    if assign is None or not isinstance(assign.value, ast.Dict):
        raise RuntimeError(_TABLE_HINT.format(name=HINDCAST.name))

    table = {}
    for key, value in zip(assign.value.keys, assign.value.values):
        if not (isinstance(key, ast.Constant) and isinstance(key.value, str)):
            raise RuntimeError(_TABLE_HINT.format(name=HINDCAST.name))
        if not (isinstance(value, ast.Call)
                and isinstance(value.func, ast.Name)
                and value.func.id == "dict"):
            raise RuntimeError(_TABLE_HINT.format(name=HINDCAST.name))
        switches = {}
        for keyword in value.keywords:
            if keyword.arg is None or not isinstance(keyword.value,
                                                     ast.Constant):
                raise RuntimeError(_TABLE_HINT.format(name=HINDCAST.name))
            switches[keyword.arg] = keyword.value.value
        table[key.value] = switches

    if not table:
        raise RuntimeError(_TABLE_HINT.format(name=HINDCAST.name))
    return table


SCENARIO_TABLE = _scenario_table()
SCENARIOS = tuple(SCENARIO_TABLE)


# Whether any nourishment project falls inside a period
def period_has_fill(period):
    return bool(nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
        period, END_YEAR[period]).projects)


# Relocation-event years that fall inside a period
def period_relocation_years(period):
    return tuple(
        event.year for event in HATTERAS_ROAD_EVENTS
        if isinstance(event, RelocationEvent) and event.enabled
        and period <= event.year < END_YEAR[period])


# Whether relocations-on is a DISTINCT run for this cell
def relocation_applies(period, scenario):
    if not SCENARIO_TABLE[scenario].get("roadway", False):
        return False, (f"{scenario} runs with roadway management off, so a "
                       f"relocation has no setback to move; the runner forces "
                       f"the switch off and the run is its non-reloc twin")
    years = period_relocation_years(period)
    if not years:
        return False, (f"no relocation event falls in "
                       f"{period}-{END_YEAR[period]}; the 2022 bridge event "
                       f"is not gated by the switch, so a reloc run is "
                       f"identical to its non-reloc twin")
    return True, None


# Whether a scenario is a DISTINCT run in this period
def scenario_applies(period, scenario):
    if scenario == "full_no_fill" and not period_has_fill(period):
        return False, (f"{period}-{END_YEAR[period]} has no nourishment "
                       f"scheduled, so full_no_fill is identical to "
                       f"full_management")
    return True, None

# The seed runs exist only to give the drift guard a reference
SEED_M, SEED_F = 50.0, 0.9
SEED_PRESET, SEED_SCENARIO = "edgeBE", "full_management"

# A 20-year run is a few minutes
RUN_TIMEOUT_S = 3600
SWEEP_TIMEOUT_S = 12 * 3600


# Exclusion

# Refuses to start while another driver is running
class DriverLock:

    def __init__(self, path=LOCK):
        self.path = path
        self.acquired = False

    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        try:
            handle = os.open(str(self.path),
                             os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            holder = ""
            try:
                holder = self.path.read_text(encoding="utf-8").strip()
            except OSError:
                pass
            raise SystemExit(
                f"another driver holds {self.path}\n"
                f"  {holder}\n"
                f"  Every run rewrites {PARAMETER_FILE.name}; two drivers at "
                f"once corrupt it.\n"
                f"  If no driver is running, delete the lock file and "
                f"re-invoke.")
        with os.fdopen(handle, "w") as stream:
            stream.write(f"pid {os.getpid()} started "
                         f"{datetime.now().isoformat(timespec='seconds')}")
        self.acquired = True
        return self

    def __exit__(self, *_exc):
        if self.acquired:
            try:
                self.path.unlink()
            except OSError:
                pass
        return False


# Manifest

# Fingerprint of the source/sink VALUES a (period, preset) pair imposes
def be_digest_for(period, preset):
    if period is None or preset is None:
        return "empty"
    return values_digest(HATTERAS_BE_PRESETS.get(preset, {}).get(period, {}))


# Stable identity for one unit of work
def job_key(stage, period=None, preset=None, scenario=None, groin=None,
            reloc=None, M=None, fraction=None, hs=None):
    parts = [stage, period, preset, scenario, groin, reloc]
    if groin:
        parts += [M, fraction]
    digest = be_digest_for(period, preset)
    if digest != "empty":
        parts.append(digest)
    # Hs JOINED THE KEY ON 2026-09-01, the fourth instance of the failure this docstring already records three of
    if hs is not None and float(hs) != 2.5:
        parts.append(f"Hs{float(hs):g}")
    return "|".join(str(x) for x in parts)


# Returns {job_key
def load_manifest():
    if not MANIFEST.exists():
        return {}
    done = {}
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if line.strip():
            row = json.loads(line)
            done[row["key"]] = row
    return done


# Appends one job outcome immediately, so an interrupt costs one job
def record(row):
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(row) + "\n")


# True only if the job is recorded AND recorded as having succeeded
def already_done(manifest, key):
    row = manifest.get(key)
    if not (row and row.get("ok") is True):
        return False
    # The key names the start year only; a row from the same start's earlier window is a different job
    name, period = row.get("run_name"), row.get("period")
    if name and period in END_YEAR and not name.startswith(f"HAT_{period}_{END_YEAR[period]}_"):
        return False
    return True


# Running one hindcast

# The line the runner prints once it has derived the run's name
_RUN_NAME_LINE = re.compile(r"^RUN_NAME_BASE\s+'([^']+)'", re.M)


# The run name the runner REPORTED, read back out of its own output
def run_name_from_log(log_path):
    if not log_path:
        return None
    try:
        text = pathlib.Path(log_path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = _RUN_NAME_LINE.search(text)
    return match.group(1) if match else None


# Runs the hindcast once, driven entirely through the environment
def run_hindcast(period, preset, scenario, reloc, groin, M, fraction,
                 overwrite, save_state, log_name, dry_run, hs=None):
    env = {key: value for key, value in os.environ.items()
           if not key.startswith("HAT_")}
    _stray = sorted(k for k in os.environ if k.startswith("HAT_"))
    if _stray:
        print(f"      ignoring {len(_stray)} HAT_* variable(s) from the "
              f"shell: {', '.join(_stray)}")
    env.update({
        # Read by HAT_hindcast_config; see the docstring above.
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": preset,
        "HAT_SCENARIO": scenario,
        # The production matrix, filed under raw_runs/matrix/ (2026-09-16).
        "HAT_RUN_KIND": "matrix",
        "HAT_RELOCATIONS": "true" if reloc else "false",
        "HAT_GROIN_ENABLED": "true" if groin else "false",
        "HAT_GROIN_TRAPPING_RATE_M_YR": str(M),
        "HAT_GROIN_DETERIORATION_FRACTION": str(fraction),
        "HAT_OVERWRITE": "true" if overwrite else "false",
        "HAT_SAVE_MODEL_STATE": "true" if save_state else "false",
        # Matplotlib must not try to open a window in an unattended run.
        "MPLBACKEND": "Agg",
    })
    # Hs is one of the six settings this driver deliberately leaves to the code default
    if hs is not None:
        env["HAT_HS"] = repr(float(hs))

    if dry_run:
        print(f"      would run: period={period} preset={preset} "
              f"scenario={scenario} reloc={reloc} groin={groin} "
              f"M={M} f={fraction}")
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


# Runs the matrix cells for one groin state, serially
def matrix_stage(stage, groin, manifest, args, fits=None):
    jobs, skipped = [], []
    for period in PERIODS:
        for preset in args.presets:
            for scenario in args.scenarios:
                applies, reason = scenario_applies(period, scenario)
                if not applies:
                    skipped.append((period, preset, scenario, False, reason))
                    continue
                # Reloc-off first, for the same reason the no-groin matrix runs before the groin one
                jobs.append((period, preset, scenario, False))
                reloc_ok, reloc_reason = relocation_applies(period, scenario)
                if reloc_ok:
                    jobs.append((period, preset, scenario, True))
                else:
                    skipped.append((period, preset, scenario, True,
                                    reloc_reason))

    completed = failed = 0
    print(f"\n{'=' * 72}\nSTAGE {stage}: {len(jobs)} runs (serial)\n{'=' * 72}")
    # (M, fraction) for this preset's groin runs, or (None, None)
    def _fit_values(preset):
        fit = (fits or {}).get(preset) if groin else None
        if not fit:
            return None, None
        return fit.get("M"), fit.get("fraction")

    if skipped:
        print(f"  {len(skipped)} cell(s) skipped as degenerate:")
        for period, preset, scenario, reloc, reason in skipped:
            arm = "reloc  " if reloc else "       "
            print(f"    {period} {preset:<7} {scenario:<16} {arm}-- {reason}")
            _m, _f = _fit_values(preset)
            key = job_key(stage, period, preset, scenario, groin, reloc,
                          _m, _f, args.hs)
            if not already_done(manifest, key) and not args.dry_run:
                # Recorded ok (no retry) with skipped=True, so the summary tells it from a real success
                record(dict(key=key, stage=stage, period=period,
                            preset=preset, scenario=scenario, groin=groin,
                            reloc=reloc, be_digest=be_digest_for(period, preset),
                            ok=True, skipped=True, detail=reason,
                            at=datetime.now().isoformat(timespec="seconds")))

    for index, (period, preset, scenario, reloc) in enumerate(jobs, start=1):
        _m, _f = _fit_values(preset)
        key = job_key(stage, period, preset, scenario, groin, reloc,
                      _m, _f, args.hs)
        if already_done(manifest, key):
            print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
                  f"{scenario:<16} reloc={reloc!s:<5} groin={groin}  "
                  f"SKIP (done)")
            completed += 1
            continue

        M = fraction = 0.0
        # Off by default: --overwrite is for a directory that exists but is stale
        overwrite = bool(getattr(args, "overwrite", False))
        if groin:
            fit = (fits or {}).get(preset)
            if not fit:
                print(f"  [{index:>2}/{len(jobs)}] no fitted (M, f) for "
                      f"{preset}; skipping its groin runs")
                failed += 1
                continue
            M, fraction = fit["M"], fit["fraction"]
            # The two seed cells already hold a run at provisional values
            overwrite = overwrite or (preset == SEED_PRESET
                                      and scenario == SEED_SCENARIO)

        # A `reloc` token only where the arm carries one, matching how the runner derives RUN_NAME in 7.5
        # The window, not just the start: 2010-2026 must not overwrite the 2010-2024 logs
        label = "_".join(
            [f"{period}_{END_YEAR[period]}", preset, scenario]
            + (["reloc"] if reloc else [])
            + ["groin" if groin else "nogroin"])
        print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
              f"{scenario:<16} reloc={reloc!s:<5} groin={groin}"
              + (f" M={M:g} f={fraction:g}" if groin else ""), flush=True)

        ok, seconds, detail = run_hindcast(
            period, preset, scenario, reloc, groin, M, fraction, overwrite,
            not args.no_model_state, label, args.dry_run)
        # Not recorded on a dry run: a preview must never cancel the job it checks
        if not args.dry_run:
            record(dict(key=key, stage=stage, period=period, preset=preset,
                        scenario=scenario, groin=groin, reloc=reloc, M=M,
                        fraction=fraction, ok=ok, seconds=round(seconds, 1),
                        detail=detail,
                        # Read back from the runner's output, not derived -- see run_name_from_log
                        run_name=run_name_from_log(detail if ok else None),
                        at=datetime.now().isoformat(timespec="seconds")))
        if ok:
            completed += 1
            print(f"        done in {seconds / 60:.1f} min")
        else:
            failed += 1
            print(f"        FAILED: {detail}")

    return completed, failed


# Stages

# Moves the retired run_index.csv aside, once
def stage_archive(manifest, args):
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


# Runs the per-period drift-guard reference at provisional M and f
def stage_seed(manifest, args):
    print(f"\n{'=' * 72}\nSTAGE seed: {len(PERIODS)} runs\n{'=' * 72}")
    print(f"  provisional M = {SEED_M:g}, f = {SEED_F:g} -- these exist only "
          f"to give the\n  sweep's drift guard a rate curve to difference "
          f"against. Stage 6 re-runs\n  both cells at the fitted values.")
    completed = failed = 0
    for period in PERIODS:
        key = job_key("seed", period, SEED_PRESET, SEED_SCENARIO, True,
                      False, hs=args.hs)
        if already_done(manifest, key):
            print(f"  {period}: SKIP (done)")
            completed += 1
            continue
        print(f"  {period} {SEED_PRESET} {SEED_SCENARIO} groin=True", flush=True)
        # Relocations off: the seed only gives the sweep's drift guard a curve
        ok, seconds, detail = run_hindcast(
            period, SEED_PRESET, SEED_SCENARIO, False, True, SEED_M, SEED_F,
            False, not args.no_model_state, f"seed_{period}", args.dry_run,
            hs=args.hs)
        if not args.dry_run:   # see the note in matrix_stage
            record(dict(key=key, stage="seed", period=period,
                        preset=SEED_PRESET, scenario=SEED_SCENARIO,
                        groin=True, reloc=False, M=SEED_M, fraction=SEED_F,
                        ok=ok, seconds=round(seconds, 1), detail=detail,
                        at=datetime.now().isoformat(timespec="seconds")))
        if ok:
            completed += 1
            print(f"    done in {seconds / 60:.1f} min")
        else:
            failed += 1
            print(f"    FAILED: {detail}")
    return completed, failed


# Runs the four period/preset sweeps
def stage_sweeps(manifest, args):
    completed = failed = 0
    # Period 1 edgeBE is 215 of the 344 cells
    order = [(period, preset) for period, preset in
             [(1984, "edgeBE"), (1984, "zeroBE"),
              (2004, "edgeBE"), (2004, "zeroBE")]
             if preset in args.presets]
    print(f"\n{'=' * 72}\nSTAGE sweeps: {len(order)} sweeps\n{'=' * 72}")
    for period, preset in order:
        key = job_key("sweep", period, preset, hs=args.hs)
        if already_done(manifest, key):
            print(f"  {period} {preset}: SKIP (done)")
            completed += 1
            continue
        command = [sys.executable, str(SWEEP), "--period", str(period),
                   "--preset", preset, "--workers", str(args.workers)]
        sweep_env = dict(os.environ)
        if args.hs is not None:
            # Read by HAT_groin_sweep_worker, and by sweep_output_dir, which sends the results to their own directory
            sweep_env["HAT_SWEEP_HS"] = repr(float(args.hs))
        if args.dry_run:
            command.append("--dry-run")
        print(f"\n  --> {' '.join(command[1:])}", flush=True)

        LOG_DIR.mkdir(parents=True, exist_ok=True)
        log_path = LOG_DIR / f"sweep_{period}_{preset}.log"
        t0 = time.perf_counter()
        # Streamed to a log rather than captured
        with log_path.open("w", encoding="utf-8") as handle:
            proc = subprocess.run(command, cwd=str(PROJECT_BASE_DIR),
                                  stdout=handle, stderr=subprocess.STDOUT,
                                  env=sweep_env, timeout=SWEEP_TIMEOUT_S)
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


# Runs the joint fit and returns the fitted values per preset
def stage_joint_fit(args):
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


# Run: every stage asked for, in order
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--workers", type=int, default=8,
                        help="sweep pool width (default 8; RAM-capped)")
    parser.add_argument("--dry-run", action="store_true",
                        help="print the plan without running anything")
    parser.add_argument("--overwrite", action="store_true",
                        help="replace matrix run directories that already "
                             "exist, instead of failing on the runner's "
                             "collision guard. Needed when a run has to be "
                             "REDONE rather than merely completed -- a "
                             "changed BE preset, a changed estimator -- "
                             "where the directory on disk is stale rather "
                             "than missing.")
    parser.add_argument("--no-model-state", action="store_true",
                        help="skip the ~160 MB .npz per matrix run "
                             "(saves ~6 GB; figures then need a re-run)")
    parser.add_argument("--stages", default="1,2,3,4,5,6",
                        help="comma-separated stage numbers to run "
                             "(default all)")
    parser.add_argument("--scenarios", default=",".join(SCENARIOS),
                        help="comma-separated scenarios to cover "
                             "(default all; the reloc arm of each is still "
                             "decided by relocation_applies)")
    parser.add_argument("--hs", type=float, default=None,
                        help="significant wave height for every run and sweep "
                             "cell (default: the code value, 2.5). Anything "
                             "other than 2.5 puts the matrix runs in their own "
                             "directories via the waveHs run-name token and "
                             "the sweeps in their own output directory, so an "
                             "Hs experiment cannot overwrite the 2.5 m one.")
    parser.add_argument("--presets", default=",".join(PRESETS),
                        help="comma-separated source/sink presets to cover "
                             f"(default {','.join(PRESETS)}; also accepts "
                             f"{', '.join(p for p in ALL_PRESETS if p not in PRESETS)}, "
                             f"which the sweep cannot fit)")
    args = parser.parse_args()

    stages = {int(s) for s in args.stages.split(",") if s.strip()}
    # Order follows PRESETS, not the order they were typed
    named = [p.strip() for p in args.presets.split(",") if p.strip()]
    unknown = [p for p in named if p not in ALL_PRESETS]
    if unknown:
        parser.error(f"unknown preset(s) {unknown}; expected some of "
                     f"{list(ALL_PRESETS)}")
    args.presets = tuple(p for p in ALL_PRESETS if p in named)
    if not args.presets:
        parser.error("--presets selected nothing")
    # Said once, here, rather than left for stage 4 to discover
    unsweepable = [p for p in args.presets if p not in PRESETS]
    if unsweepable:
        print(f"  note: {', '.join(unsweepable)} has no sweep grid, so "
              f"stages 4-6 skip it; stages 1-3 run normally")

    # Validated against the runner's table and put in its order
    chosen = [s.strip() for s in args.scenarios.split(",") if s.strip()]
    unknown = [s for s in chosen if s not in SCENARIOS]
    if unknown:
        parser.error(f"unknown scenario(s) {unknown}; expected some of "
                     f"{list(SCENARIOS)}")
    args.scenarios = tuple(s for s in SCENARIOS if s in chosen)
    if not args.scenarios:
        parser.error("--scenarios selected nothing")
    started = datetime.now()
    print("=" * 72)
    print(f"HATTERAS FULL RUN  started {started:%Y-%m-%d %H:%M:%S}")
    print("=" * 72)
    print(f"  stages       {sorted(stages)}")
    print(f"  presets      {', '.join(args.presets)}")
    print(f"  scenarios    {', '.join(args.scenarios)}")
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

    # Stage 6 alone falls back to joint_fit.json for (M, f)
    if 5 in stages:
        fits = stage_joint_fit(args)
    elif JOINT_JSON.exists():
        fits = json.loads(JOINT_JSON.read_text())
        summary = ", ".join(f"{k} M={v.get('M')} f={v.get('fraction')}"
                            for k, v in fits.items())
        print(f"  stage 5 not requested; read fitted values from "
              f"{JOINT_JSON.name}: {summary}")
    else:
        fits = {}

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
    # Held for the whole invocation rather than per stage
    with DriverLock():
        sys.exit(main())
