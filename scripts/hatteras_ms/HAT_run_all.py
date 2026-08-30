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

    2  matrix nogroin The scenarios that are DISTINCT in each period x
                      preset x 2 periods x the relocations arm, groin off.
                      These are the comparison baselines every groin result is
                      read against. They need no M or f, so they can run
                      before anything is fitted -- and running them first also
                      satisfies the rule that a scenario's no-groin run must
                      exist before its groin run, which is how section 12.3
                      resolves its paired baseline.

                      Two guards decide which cells exist, and both skip
                      rather than force a name onto what would be a duplicate
                      run:

                      `scenario_applies` drops `full_no_fill` in 1984-2004.
                      All three nourishment projects are 2014 or later, so
                      there it differs from `full_management` in a switch with
                      nothing to act on: the two build the same modules,
                      derive the same RUN_NAME, and the collision guard
                      refuses the second. The fills axis exists only in
                      period 2.

                      `relocation_applies` drops the reloc arm wherever it
                      would be inert. Relocations need roadway management to
                      have a setback to move, so `natural` and
                      `beachdune_only` carry no arm; and both historical
                      events are 1989 and 1999, so 2004-2024 carries none
                      either -- the 2022 bridge event fires whatever the
                      switch says, which is exactly why turning the switch on
                      there would produce an identical run under a name
                      claiming a contrast.

                      Both presets, both periods, groin off: 22 runs.
                      Per preset, 1984-2004 contributes 4 scenarios + 2 reloc
                      arms and 2004-2024 contributes 5 scenarios + 0 arms, so
                      11 each. `--presets zeroBE` covers half of it.

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
    lands, keyed by (stage, period, preset, scenario, groin, reloc) -- plus
    (M, fraction) when the groin is on, and plus the SOURCE/SINK DIGEST when
    the preset imposes any, so a revised fit invalidates the runs made at the
    old values instead of silently reusing them. Re-invoking skips anything
    already recorded complete and retries anything recorded failed. The
    manifest -- not the run directory -- is the source of truth, so this file
    never has to reimplement the runner's RUN_NAME construction, which is
    derived in the runner's section 7.5 from what its modules actually built
    and would drift if copied.

    THE KEY MUST NAME EVERY INPUT THAT CAN CHANGE UNDER A FIXED NAME. That
    rule has now been learned three times: the reloc axis, then (M, fraction)
    for groin runs on 2026-08-24, then the source/sink digest on 2026-08-28.
    Each time the symptom was identical -- a run made at superseded values
    matched an unchanged key and was reported "SKIP (done)" -- and each time
    the skip happens BEFORE the collision guard, so --overwrite cannot reach
    it. If a fourth quantity is ever edited between runs while its preset
    keeps its name, put it in the key rather than remembering to clear the
    manifest by hand.

    A manifest written before the M/f fields existed carries six-field groin
    keys that can never match the eight-field ones. Migrate it rather than
    letting every groin run repeat: the rows already RECORD M and fraction, so
    the new key can be rebuilt from each row exactly.

    The digest is NOT retrofittable the same way: rows written before
    2026-08-28 do not carry `be_digest`, and it cannot be recovered from the
    row because the whole point is that the values have since changed. Those
    rows simply stop matching, which re-runs them -- the safe direction. Only
    zeroBE is unaffected, because an "empty" digest is never appended and
    zeroBE imposes {} in every period.

CONCURRENCY
    Matrix runs are SERIAL. Each one appends a row to run_index.csv, and
    concurrent appends to one CSV interleave. The sweep is parallel -- its
    cells write only to their own directories and are collected through a
    single thread -- and it is where nearly all the wall clock goes.

Usage:
    python HAT_run_all.py [--workers N] [--dry-run] [--no-model-state]
                          [--stages 2,4,5] [--presets zeroBE]

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import ast
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
# The sweep, its config and the joint fit live in groin-sweep/. That is a
# plain directory and not a package -- the hyphen in the name makes it
# unimportable as one -- so its path goes on sys.path and the config is
# imported by module name, exactly as it was when the four files sat beside
# this one. Without this the driver dies on ModuleNotFoundError before it
# reaches stage 1.
GROIN_SWEEP_DIR = _HERE.parent / "groin-sweep"
for _path in (SCRIPTS_DIR, _HERE.parent, GROIN_SWEEP_DIR):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline import nourishment  # noqa: E402
from cascade_pipeline.roadway import RelocationEvent  # noqa: E402
from cascade_pipeline.run_registry import values_digest  # noqa: E402
from hatteras_site_config import (  # noqa: E402
    HATTERAS_BE_PRESETS,
    HATTERAS_DOMAINS,
    HATTERAS_NOURISHMENT_PROJECTS,
    HATTERAS_ROAD_EVENTS,
)

from HAT_groin_sweep_config import END_YEAR, PERIODS, PRESETS  # noqa: E402

# The matrix's preset axis is WIDER than the sweep's. PRESETS is the pair the
# groin sweep fits M and f against, and it stays the default because the
# comparison matrix exists to feed that fit. But a hindcast run is valid under
# any preset the site config defines, and calibBE is the only one that puts a
# source/sink term on the interior domains -- so it is reachable through
# --presets even though no sweep is fitted for it. Order follows PRESETS first
# so the default schedule is unchanged.
ALL_PRESETS = tuple(PRESETS) + tuple(
    name for name in HATTERAS_BE_PRESETS if name not in PRESETS)

HINDCAST = _HERE.parent / "HAT_hindcast_1984_2024.py"
SWEEP = GROIN_SWEEP_DIR / "HAT_groin_sweep.py"
JOINT_FIT = GROIN_SWEEP_DIR / "HAT_groin_joint_fit.py"

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
SWEEP_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep"
JOINT_JSON = SWEEP_DIR / "joint_fit.json"

PARAMETER_FILE = (PROJECT_BASE_DIR / "data" / "hatteras_init"
                  / "Hatteras-CASCADE-parameters.yaml")

DRIVER_DIR = PROJECT_BASE_DIR / "output" / "driver"
MANIFEST = DRIVER_DIR / "driver_manifest.jsonl"
LOG_DIR = DRIVER_DIR / "logs"
LOCK = DRIVER_DIR / "driver.lock"

_TABLE_HINT = (
    "expected a module-level `SCENARIOS = {{...}}` of literal dict(...) calls "
    "in section 3 of {name}; if that changed, this parser changes with it")


def _scenario_table():
    """Reads the runner's SCENARIOS table out of its source, without running it.

    `HAT_hindcast_1984_2024.py` is a script, not a module: importing it runs a
    whole hindcast, so the table cannot simply be imported. It is parsed
    instead, because the literal dict in the runner's section 3 is the one
    thing that decides which management switches a scenario name stands for,
    and a copy typed out here would be free to drift from it.

    Drift in the `roadway` flag would be the damaging kind. `relocation_applies`
    below reads it to decide whether a reloc arm is a distinct run; if it said
    True where the runner says False, the runner would force the switch off,
    derive its NON-reloc name in 7.5, and the arm would collide with the very
    baseline it was meant to be read against.

    Returns:
        {scenario_name: {switch_name: bool}}, in the runner's own order.

    Raises:
        RuntimeError: If the table is absent or is no longer written as
            literal dict(...) calls. Loud on purpose: the alternative is a
            silently empty or partial scenario list, which reads downstream
            as "there was nothing to run".
    """
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


def period_has_fill(period):
    """Whether any nourishment project falls inside a period.

    Resolved with the same `build_schedule` call the runner makes, rather than
    by re-implementing the date filter, so this cannot disagree with what
    section 6 actually builds.
    """
    return bool(nourishment.build_schedule(
        HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS,
        period, END_YEAR[period]).projects)


def period_relocation_years(period):
    """Relocation-event years that fall inside a period.

    Derived from HATTERAS_ROAD_EVENTS rather than listed here, so adding an
    event adds its runs. The window matches the loop's in
    `run_cascade_simulation`: `current_year` runs from `start_year` to
    `start_year + run_years - 1`, so an event dated exactly on a period's end
    year belongs to the NEXT period and never fires in this one.

    BridgeEvents are deliberately not counted. `apply_historical_event`
    applies them whatever `relocations_enabled` says, so they are not part of
    what the switch turns on.
    """
    return tuple(
        event.year for event in HATTERAS_ROAD_EVENTS
        if isinstance(event, RelocationEvent) and event.enabled
        and period <= event.year < END_YEAR[period])


def relocation_applies(period, scenario):
    """Whether relocations-on is a DISTINCT run for this cell.

    Two ways it is not, and either would put two byte-identical runs on disk
    under names claiming to be a contrast:

    ROADWAY MANAGEMENT OFF
        The runner's `_RELOCATIONS_FORCED_OFF` turns the switch off wherever
        nothing manages the road -- a relocation cannot move a setback that no
        RoadwayManager reads. The forced-off run then derives its NON-reloc
        name in 7.5 and the collision guard refuses it as a duplicate.

    NO RELOCATION EVENT IN THE PERIOD
        Both historical relocations are 1989 and 1999, so 2004-2024 has none.
        The 2022 Jug Handle bridge event is not gated by the switch, so a
        period-2 reloc run simulates exactly what its non-reloc twin does.

    Same shape and same motive as `scenario_applies`: skipping is the honest
    handling, because running it under a forced name files two identical runs
    as if they were a contrast.

    Returns:
        (applies, reason). reason is None when it applies.
    """
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
# EXCLUSION
# =============================================================================

class DriverLock:
    """Refuses to start while another driver is running.

    Every run rewrites PARAMETER_FILE and reads it back, so two drivers
    interleave writes to one file and tear each other's reads. The lock is an
    O_EXCL create, which is atomic on Windows and POSIX alike; the PID inside
    is for the human reading the error, not for liveness checking.

    A crashed driver leaves the file behind. That is deliberate -- a stale
    lock reports what to delete rather than being cleared automatically, since
    "the lock looks old" is exactly the reasoning that reintroduces the race.
    """

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


# =============================================================================
# MANIFEST
# =============================================================================

def be_digest_for(period, preset):
    """Fingerprint of the source/sink VALUES a (period, preset) pair imposes.

    Not the preset name -- the numbers behind it. `values_digest` is imported
    from the run registry rather than reimplemented so the manifest key and
    the `be_values_digest` column in run_index.csv can never disagree about
    what a run carried.

    Returns "empty" for a preset that imposes nothing, and for the stages
    (archive, sweep) whose keys carry no period/preset pair.
    """
    if period is None or preset is None:
        return "empty"
    return values_digest(HATTERAS_BE_PRESETS.get(preset, {}).get(period, {}))


def job_key(stage, period=None, preset=None, scenario=None, groin=None,
            reloc=None, M=None, fraction=None):
    """Stable identity for one unit of work.

    The relocations flag joined the key when the reloc axis was added. A key
    written before that has five fields and can never match a six-field one,
    so an older manifest does not mark the new cells done -- which is the
    correct outcome, since those runs were made without the axis and the
    reloc arms among them do not exist.

    M AND f JOINED THE KEY FOR GROIN RUNS ON 2026-08-24, for the same reason
    and after the same failure. The key identified a groin cell by scenario
    alone, so when the fitted values changed the old run still matched and was
    SKIPPED -- and the skip happens before the collision guard, so --overwrite
    could not reach it either. A run made at M = 110, f = 0.4 survived into a
    matrix that was supposed to be M = 50, f = 0.6, reported as "SKIP (done)".
    Including the values means a revised fit invalidates exactly the runs it
    should and leaves the rest alone.

    Only groin runs carry them: a no-groin cell has no M or f, and appending
    "None|None" to its key would invalidate every groin-off run already
    recorded for no reason.

    THE SOURCE/SINK DIGEST JOINED THE KEY ON 2026-08-28, third instance of
    the same failure. The key named the PRESET but not its VALUES, and the
    source/sink table is edited between runs -- that is the whole shape of an
    edge solve, where `edgeBE` keeps its name while GIS 1 and 90 move at every
    Newton step. So a run made at the previous edge values still matched, and
    was reported "SKIP (done)" into a matrix that was supposed to carry the
    new ones. Caught before it corrupted anything only because the edgeBE leg
    was killed and restarted by hand.

    Derived here from (period, preset) rather than passed in, so a call site
    cannot forget it. `values_digest` returns "empty" for a preset that
    imposes nothing, and an "empty" digest is NOT appended -- that keeps every
    zeroBE key byte-identical to what is already on disk, which is correct:
    zeroBE imposes {} in every period and no edit to the calibrated table can
    change what it ran.
    """
    parts = [stage, period, preset, scenario, groin, reloc]
    if groin:
        parts += [M, fraction]
    digest = be_digest_for(period, preset)
    if digest != "empty":
        parts.append(digest)
    return "|".join(str(x) for x in parts)


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

def run_hindcast(period, preset, scenario, reloc, groin, M, fraction,
                 overwrite, save_state, log_name, dry_run):
    """Runs the hindcast once, driven entirely through the environment.

    The runner reads these through HAT_hindcast_config, so no tracked source
    is edited and the notebook stays byte-identical to what a driven run uses.

    HAT_RELOCATIONS is set on EVERY run, including the reloc-off ones, rather
    than only when it departs from the scenario preset. Every named scenario
    sets relocations=False, so an explicit "false" is not a departure and the
    runner reports none -- but `describe()` then records the value as having
    come from the environment, so a run log states that the driver chose it
    instead of leaving the reader to infer a default.

    THE SETTINGS FILE IS SHUT OUT, TWICE OVER
        HAT_IGNORE_SETTINGS=1 stops the runner reading `hat_run.yaml`, so a
        batch run is described entirely by the values set here plus
        HAT_hindcast_config's defaults. The driver sets only nine of the
        fifteen settings; without this, the other six -- offset_mode,
        make_gifs, Hs, sandbags, use_sandbox_cascade, show_figures -- would be
        taken from whatever experiment was last left in the file, and every
        cell of the matrix would silently inherit it.

        Stray HAT_* variables in the calling shell are dropped for the same
        reason, and reported when they are, since `os.environ.copy()` would
        otherwise carry them into every child.

    Returns:
        (ok, seconds, detail).
    """
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
        "HAT_RELOCATIONS": "true" if reloc else "false",
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


def matrix_stage(stage, groin, manifest, args, fits=None):
    """Runs the matrix cells for one groin state, serially.

    Three axes: period x preset x scenario, plus a relocations arm on every
    cell where relocations-on is a distinct run. `scenario_applies` and
    `relocation_applies` between them decide which cells exist; a cell they
    refuse is recorded as skipped WITH its reason rather than dropped
    silently, so the manifest distinguishes "did not need running" from "was
    never scheduled".

    Args:
        stage: "matrix_nogroin" or "matrix_groin".
        groin: Whether the groin is attached.
        manifest: Loaded manifest, for resume.
        args: Parsed CLI arguments. `args.presets` and `args.scenarios`
            select which slice of the matrix this invocation covers.
        fits: {preset: {"M": ..., "fraction": ...}} for the groin stage.

    Returns:
        (completed, failed) counts.
    """
    jobs, skipped = [], []
    for period in PERIODS:
        for preset in args.presets:
            for scenario in args.scenarios:
                applies, reason = scenario_applies(period, scenario)
                if not applies:
                    skipped.append((period, preset, scenario, False, reason))
                    continue
                # reloc-off first, for the same reason the no-groin matrix
                # runs before the groin one: a reloc arm is read against the
                # arm that is identical to it in every other token, and the
                # baseline should be on disk before the arm that cites it.
                jobs.append((period, preset, scenario, False))
                reloc_ok, reloc_reason = relocation_applies(period, scenario)
                if reloc_ok:
                    jobs.append((period, preset, scenario, True))
                else:
                    skipped.append((period, preset, scenario, True,
                                    reloc_reason))

    completed = failed = 0
    print(f"\n{'=' * 72}\nSTAGE {stage}: {len(jobs)} runs (serial)\n{'=' * 72}")
    def _fit_values(preset):
        """(M, fraction) for this preset's groin runs, or (None, None)."""
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
            key = job_key(stage, period, preset, scenario, groin, reloc, _m, _f)
            if not already_done(manifest, key) and not args.dry_run:
                # Recorded as ok so a re-invocation does not retry it, with
                # skipped=True so the summary can tell "did not need running"
                # from "ran and succeeded".
                record(dict(key=key, stage=stage, period=period,
                            preset=preset, scenario=scenario, groin=groin,
                            reloc=reloc, be_digest=be_digest_for(period, preset),
                            ok=True, skipped=True, detail=reason,
                            at=datetime.now().isoformat(timespec="seconds")))

    for index, (period, preset, scenario, reloc) in enumerate(jobs, start=1):
        _m, _f = _fit_values(preset)
        key = job_key(stage, period, preset, scenario, groin, reloc, _m, _f)
        if already_done(manifest, key):
            print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
                  f"{scenario:<16} reloc={reloc!s:<5} groin={groin}  "
                  f"SKIP (done)")
            completed += 1
            continue

        M = fraction = 0.0
        # Off by default: the guard is what stops one matrix run silently
        # replacing another. --overwrite is for the case the guard cannot
        # distinguish -- a directory that exists but is STALE, because the
        # forcing or the reported quantity changed underneath it.
        overwrite = bool(getattr(args, "overwrite", False))
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
            overwrite = overwrite or (preset == SEED_PRESET
                                      and scenario == SEED_SCENARIO)

        # A `reloc` token only where the arm carries one, matching how the
        # runner derives RUN_NAME in 7.5. Spelling it "noreloc" everywhere
        # else would rename every existing log for a switch that is off.
        label = "_".join(
            [str(period), preset, scenario]
            + (["reloc"] if reloc else [])
            + ["groin" if groin else "nogroin"])
        print(f"  [{index:>2}/{len(jobs)}] {period} {preset:<7} "
              f"{scenario:<16} reloc={reloc!s:<5} groin={groin}"
              + (f" M={M:g} f={fraction:g}" if groin else ""), flush=True)

        ok, seconds, detail = run_hindcast(
            period, preset, scenario, reloc, groin, M, fraction, overwrite,
            not args.no_model_state, label, args.dry_run)
        # NOT recorded on a dry run. The manifest means "this job is done",
        # and a dry run did not do it -- a row written here would be read by
        # `already_done` on the next real invocation and skip the very cell
        # the preview was checking. Previewing a plan must never cancel it.
        if not args.dry_run:
            record(dict(key=key, stage=stage, period=period, preset=preset,
                        scenario=scenario, groin=groin, reloc=reloc, M=M,
                        fraction=fraction, ok=ok, seconds=round(seconds, 1),
                        detail=detail,
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
        key = job_key("seed", period, SEED_PRESET, SEED_SCENARIO, True,
                      False)
        if already_done(manifest, key):
            print(f"  {period}: SKIP (done)")
            completed += 1
            continue
        print(f"  {period} {SEED_PRESET} {SEED_SCENARIO} groin=True", flush=True)
        # reloc off: the seed exists only to give the sweep's drift guard a
        # rate curve to difference against, and the sweep's own cells are
        # built without relocations.
        ok, seconds, detail = run_hindcast(
            period, SEED_PRESET, SEED_SCENARIO, False, True, SEED_M, SEED_F,
            False, not args.no_model_state, f"seed_{period}", args.dry_run)
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


def stage_sweeps(manifest, args):
    """Runs the four period/preset sweeps."""
    completed = failed = 0
    # Period 1 edgeBE is 215 of the 344 cells; running it first means the
    # longest leg starts earliest and a morning check finds it either done or
    # visibly progressing, rather than not yet begun.
    order = [(period, preset) for period, preset in
             [(1984, "edgeBE"), (1984, "zeroBE"),
              (2004, "edgeBE"), (2004, "zeroBE")]
             if preset in args.presets]
    print(f"\n{'=' * 72}\nSTAGE sweeps: {len(order)} sweeps\n{'=' * 72}")
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
    parser.add_argument("--presets", default=",".join(PRESETS),
                        help="comma-separated source/sink presets to cover "
                             f"(default {','.join(PRESETS)}; also accepts "
                             f"{', '.join(p for p in ALL_PRESETS if p not in PRESETS)}, "
                             f"which the sweep cannot fit)")
    args = parser.parse_args()

    stages = {int(s) for s in args.stages.split(",") if s.strip()}
    # Order follows PRESETS, not the order they were typed, so two
    # invocations naming the same presets schedule them the same way.
    named = [p.strip() for p in args.presets.split(",") if p.strip()]
    unknown = [p for p in named if p not in ALL_PRESETS]
    if unknown:
        parser.error(f"unknown preset(s) {unknown}; expected some of "
                     f"{list(ALL_PRESETS)}")
    args.presets = tuple(p for p in ALL_PRESETS if p in named)
    if not args.presets:
        parser.error("--presets selected nothing")
    # Said once, here, rather than left for stage 4 to discover: a preset with
    # no sweep grid still runs every matrix cell, and only the stages that
    # need a fitted (M, f) skip it.
    unsweepable = [p for p in args.presets if p not in PRESETS]
    if unsweepable:
        print(f"  note: {', '.join(unsweepable)} has no sweep grid, so "
              f"stages 4-6 skip it; stages 1-3 run normally")

    # Same treatment as --presets: validated against the runner's own table
    # and reordered into its order, so two invocations naming the same
    # scenarios schedule them identically.
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

    # Stage 6 needs fitted (M, f). Stage 5 returns them, but running stage 6
    # ALONE used to hand it {} and it skipped every groin run with "no fitted
    # (M, f)" -- 21 of 22 on 2026-08-24. Falling back to the file makes stage 6
    # runnable on its own, which matters whenever the values in joint_fit.json
    # were NOT produced by the joint fit's own ranking: re-running stage 5 to
    # satisfy stage 6 would overwrite them with the ranking's answer.
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
    # Held for the whole invocation rather than per stage: the parameter file
    # is shared by every stage, not just the matrix.
    with DriverLock():
        sys.exit(main())
