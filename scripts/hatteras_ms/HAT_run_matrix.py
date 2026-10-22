#!/usr/bin/env python3
"""
End-to-end unattended driver: the comparison matrix, without and with the groin.

    python scripts/hatteras_ms/HAT_run_matrix.py [--dry-run] [--no-model-state] [--stages 1,2,3]

Runs every stage in dependency order and resumes where it stopped; each
job recorded in a manifest. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
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
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline import nourishment  # noqa: E402
from cascade_pipeline.roadway import RelocationEvent  # noqa: E402
from cascade_pipeline.run_registry import values_digest  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_BE_PRESETS,
    HATTERAS_DOMAINS,
    HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_PERIODS,
    HATTERAS_ROAD_EVENTS,
    last_model_year,
)

from HAT_hindcast_config import field_default, groin_strength_for  # noqa: E402  (the groin, in the job key)

# --- CONFIG ------------------------------------------------------------------
# The calibration and test periods, and the presets the matrix covers by default
PERIODS = (1996, 2009)
PRESETS = ("edgeBE", "zeroBE")
ALL_PRESETS = tuple(PRESETS) + tuple(
    name for name in HATTERAS_BE_PRESETS if name not in PRESETS)
END_YEAR = {p: v["end_year"] for p, v in HATTERAS_PERIODS.items()}
LAST_MODEL_YEAR = {p: last_model_year(p) for p in HATTERAS_PERIODS}
# The groin's strength (per window), end year and spread reach the runner as code defaults; they key its runs
GROIN_STRENGTH = {p: groin_strength_for(p) for p in PERIODS}
GROIN_END_YEAR = field_default("groin_end_year")
GROIN_SPREAD_WEIGHTS = field_default("groin_spread_weights")

HINDCAST = _HERE.parent / "HAT_hindcast.py"

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"

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
        period, LAST_MODEL_YEAR[period]).projects)


# Relocation-event years that fall inside a period
def period_relocation_years(period):
    return tuple(
        event.year for event in HATTERAS_ROAD_EVENTS
        if isinstance(event, RelocationEvent) and event.enabled
        and period <= event.year <= LAST_MODEL_YEAR[period])


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

# A 20-year run is a few minutes
RUN_TIMEOUT_S = 3600


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
            reloc=None, hs=None):
    parts = [stage, period, preset, scenario, groin, reloc]
    # A groin run is keyed on the groin's settings, so changing any reruns it
    if groin:
        parts += [f"s{GROIN_STRENGTH[period]:g}", f"off{GROIN_END_YEAR}"]
        # The spread joined the key on 2026-10-21; the plain groin keeps its old key
        if tuple(GROIN_SPREAD_WEIGHTS) != (1.0, 0.0, 0.0):
            parts.append("spread" + "-".join(f"{w:.3g}" for w in GROIN_SPREAD_WEIGHTS))
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
    # A recorded skip is not a result: applicability is re-decided every invocation, and a window can gain a fill
    if row.get("skipped"):
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
def run_hindcast(period, preset, scenario, reloc, groin,
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
              f"scenario={scenario} reloc={reloc} groin={groin}")
        return True, 0.0, "dry-run"

    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log_path = LOG_DIR / f"{log_name}.log"
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, str(HINDCAST)],
            env=env, cwd=str(PROJECT_BASE_DIR),
            capture_output=True, text=True, encoding="utf-8", errors="replace",
            timeout=RUN_TIMEOUT_S)
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
def matrix_stage(stage, groin, manifest, args):
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
    if skipped:
        print(f"  {len(skipped)} cell(s) skipped as degenerate:")
        for period, preset, scenario, reloc, reason in skipped:
            arm = "reloc  " if reloc else "       "
            print(f"    {period} {preset:<7} {scenario:<16} {arm}-- {reason}")
            key = job_key(stage, period, preset, scenario, groin, reloc,
                          args.hs)
            if key not in manifest and not args.dry_run:
                # Recorded ok (no retry) with skipped=True, so the summary tells it from a real success
                record(dict(key=key, stage=stage, period=period,
                            preset=preset, scenario=scenario, groin=groin,
                            reloc=reloc, be_digest=be_digest_for(period, preset),
                            ok=True, skipped=True, detail=reason,
                            at=datetime.now().isoformat(timespec="seconds")))

    for index, (period, preset, scenario, reloc) in enumerate(jobs, start=1):
        key = job_key(stage, period, preset, scenario, groin, reloc,
                      args.hs)
        if already_done(manifest, key):
            print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
                  f"{scenario:<16} reloc={reloc!s:<5} groin={groin}  "
                  f"SKIP (done)")
            completed += 1
            continue

        # Off by default: --overwrite is for a directory that exists but is stale
        overwrite = bool(getattr(args, "overwrite", False))

        # A `reloc` token only where the arm carries one, matching how the runner derives RUN_NAME in 7.5
        # The window, not just the start: 2010-2026 must not overwrite the 2010-2024 logs
        label = "_".join(
            [f"{period}_{END_YEAR[period]}", preset, scenario]
            + (["reloc"] if reloc else [])
            + ["groin" if groin else "nogroin"])
        print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
              f"{scenario:<16} reloc={reloc!s:<5} groin={groin}"
              + (f" strength={GROIN_STRENGTH[period]:g} off={GROIN_END_YEAR}" if groin else ""), flush=True)

        ok, seconds, detail = run_hindcast(
            period, preset, scenario, reloc, groin, overwrite,
            not args.no_model_state, label, args.dry_run)
        # Not recorded on a dry run: a preview must never cancel the job it checks
        if not args.dry_run:
            record(dict(key=key, stage=stage, period=period, preset=preset,
                        scenario=scenario, groin=groin, reloc=reloc,
                        groin_strength=GROIN_STRENGTH[period] if groin else None,
                        groin_end_year=GROIN_END_YEAR if groin else None,
                        groin_spread_weights=list(GROIN_SPREAD_WEIGHTS) if groin else None,
                        ok=ok, seconds=round(seconds, 1),
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


# Run: every stage asked for, in order
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
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
    parser.add_argument("--stages", default="1,2,3",
                        help="comma-separated stage numbers to run: 1 archive the run index, "
                             "2 the no-groin matrix, 3 the groin matrix (default all)")
    parser.add_argument("--scenarios", default=",".join(SCENARIOS),
                        help="comma-separated scenarios to cover "
                             "(default all; the reloc arm of each is still "
                             "decided by relocation_applies)")
    parser.add_argument("--hs", type=float, default=None,
                        help="significant wave height for every run (default: "
                             "the code value). Anything else puts the matrix "
                             "runs in their own directories via the waveHs "
                             "run-name token, so an Hs experiment cannot "
                             "overwrite the default one.")
    parser.add_argument("--presets", default=",".join(PRESETS),
                        help="comma-separated source/sink presets to cover "
                             f"(default {','.join(PRESETS)}; also accepts "
                             f"{', '.join(p for p in ALL_PRESETS if p not in PRESETS)})")
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
    print(f"  groin        strength " + ", ".join(f"{GROIN_STRENGTH[p]:g} ({p})" for p in PERIODS)
          + f", off from {GROIN_END_YEAR}, spread " + ", ".join(f"{w:.3g}" for w in GROIN_SPREAD_WEIGHTS))
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
        manifest = load_manifest()
        totals["matrix groin"] = matrix_stage(
            "matrix_groin", True, manifest, args)

    elapsed = (datetime.now() - started).total_seconds() / 3600
    print("\n" + "=" * 72)
    print(f"FINISHED  {elapsed:.2f} h")
    print("=" * 72)
    for name, (done, failed) in totals.items():
        flag = "" if not failed else f"   {failed} FAILED"
        print(f"  {name:<16} {done} completed{flag}")
    print(f"\n  re-invoke this script to retry anything that failed")
    return 0 if not any(f for _, f in totals.values()) else 1


if __name__ == "__main__":
    # Held for the whole invocation rather than per stage
    with DriverLock():
        sys.exit(main())
