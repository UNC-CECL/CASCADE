#!/usr/bin/env python3
"""Parameter sensitivity for the Hatteras hindcast, either period.

WHY THIS REPLACED THE SCRIPTS THAT WERE HERE
    `HAT_waveSensitivity_1984_2004.py` and
    `HAT_waveheight_Sensitivity_1984_2004.py` both built the model by hand.
    Between them they had: a SyntaxError that made one of them unrunnable,
    input paths that no longer resolved (`RoadSetback_1984.csv` under
    `raw_offset/`, `storms_1984_2004_base.npy` under `hindcast_storms/` --
    neither exists; the only copies live under `old_method_offset/` and
    `old/testing_storms/`), a direct `Cascade(...)` call against the STOCK
    class rather than the sandbox one every hindcast run uses, and a comment
    pinning them to "HAT_hindcast_1984_2024_old version.py". Anything they
    produced would have described a different model on legacy forcing.

    This drives the hindcast instead of reimplementing it. Each cell is one
    ordinary run of `HAT_hindcast_1984_2024.py` with one setting moved through
    the environment, exactly the way `HAT_run_all.py` drives the matrix. There
    is no second copy of the model setup here, so the sweep inherits the period
    table, the measured setbacks, the groin and the LRR estimator
    automatically, and cannot drift from them.

    Nothing in the notebook or the headless .py knows this file exists. Those
    two were touched once, on 2026-09-01, to make the swept values settings
    rather than literals and to give an off-default value a run-name token;
    that is the whole of their involvement, and both changes stand on their own.

BOTH PERIODS COME FREE
    The period is `--start-year`, and everything period-specific is already in
    HATTERAS_PERIODS. Nothing in this file knows what year it is.

CELLS CANNOT COLLIDE WITH THE MATRIX, BUT THE TWO AXES DO IT DIFFERENTLY
    Without some separator every cell would derive the SAME name as the matrix
    run beside it, and the last one to finish would be left wearing the
    production name -- the failure `output/groin_sweep/README.md` documents for
    the rig sweep. Two mechanisms prevent it, and which one applies depends on
    the axis:

      * `relocation_setback` earns a NAME token (`rset40`) from
        `cascade_pipeline.hindcast`, so the cell is a separate directory beside
        the matrix run and a separate row in run_index.csv.

      * The four WAVE axes no longer do. The runner stopped emitting the wave
        token into the run name on 2026-09-01 -- a name describes the SCENARIO
        and a path describes the FORCING -- so a wave cell keeps the matrix
        run's name and is separated by its forcing ARM instead: an `arm` path
        component and the ("run_name", "Hs_m", "arm") index key.

    Both are still one-cell-one-directory. The cells already on disk under
    names like `..._waveHs1p2` predate that change and are left exactly where
    they are; nothing reads the token, and re-running one now files it by arm.

WHAT IS SWEPT
    Five axes, one at a time, each around its calibration value. Only the swept
    parameter moves, so every cell differs from its baseline in exactly one
    thing. Four are the wave climate; the fifth is the relocation target, which
    is not wave physics but is the newest and least settled number in the
    model -- `hat_run.yaml` records that 20 m clears the GIS 11 drowning
    threshold by ONE 10 m cell.

Usage:
    python HAT_hindcast_sensitivity.py --start-year 1984 --param wave_height
    python HAT_hindcast_sensitivity.py --start-year 1984 --param all --dry-run
    python HAT_hindcast_sensitivity.py --start-year 2004 --param wave_height \\
        --values 1.5,2.0,2.5

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts",
              PROJECT_BASE_DIR / "scripts" / "hatteras_ms"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import (  # noqa: E402
    HATTERAS_PERIODS, HATTERAS_ROAD_EVENTS)
from cascade_pipeline.roadway import RelocationEvent  # noqa: E402
from HAT_hindcast_config import field_default  # noqa: E402

HINDCAST = PROJECT_BASE_DIR / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
OUT_ROOT = PROJECT_BASE_DIR / "output" / "sensitivity_analysis"

# The "each domain relocates to its own measured offset" case. The environment
# carries strings, and an EMPTY string reads as unset in HAT_hindcast_config
# (`if raw != ""`), so this case has to be spelled: "" would silently leave the
# cell at the default and the result would be indistinguishable from one.
MEASURED = "measured"


def relocation_arm(start_year, end_year):
    """The arm the relocation-target axis has to be measured in, for a period.

    Period 1 holds the 1989 and 1999 events, and the drowning this axis asks
    about is a prescribed DISPLACEMENT landing on top of the emergent target --
    with the events off, the axis would measure something else and report a
    reassuring flat line. So it forces them on there.

    Period 2 holds no RelocationEvent (only the 2022 Jug Handle BridgeEvent),
    so forcing them on would add a `reloc` token claiming history the period
    does not have, and would name the cells for a baseline the matrix has never
    run. There the axis moves EMERGENT relocations alone, which is a real but
    different question, and the ordinary arm is the right one.

    Read from HATTERAS_ROAD_EVENTS rather than by testing the year, so adding
    an event to that table brings this with it.

    Args:
        start_year: Period start.
        end_year: Period end.

    Returns:
        Environment overrides for the arm, possibly empty.
    """
    events = [event for event in HATTERAS_ROAD_EVENTS
              if isinstance(event, RelocationEvent)
              and start_year <= event.year <= end_year]
    return {"HAT_RELOCATIONS": "true"} if events else {}


# The swept axes. `setting` is the HAT_hindcast_config field, which fixes both
# the environment name (HAT_ + upper case) and the default the name token is
# measured against -- so a value equal to the default produces an untokened
# name, which would be the matrix run. Those cells are skipped, not run.
#
# `arm` is a per-axis callable (start_year, end_year) -> environment overrides,
# for an axis that would be inert or misnamed in the default arm. A callable
# rather than a dict because whether the override is right DEPENDS ON THE
# PERIOD -- see relocation_arm.
SWEEPS = {
    # The Hs range deliberately runs below the physical mean: BRIE's
    # diffusivity scales with Hs, and the sub-1 m cells are what locate the
    # transition between erosion regimes. The USGS hindcast mean for this coast
    # is ~1.2-1.3 m and the morphologically effective height ~1.4-1.5 m, so
    # 2.5 m (the calibration value) is already high, and how much of the fit
    # rests on that is the question this axis asks.
    "wave_height": {
        "setting": "hs",
        "label": "Wave height Hs",
        "units": "m",
        "values": [0.75, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0,
                   4.5, 5.0],
    },
    "wave_period": {
        "setting": "wave_period_s",
        "label": "Wave period",
        "units": "s",
        "values": [7.0, 8.0, 9.0, 10.0],
    },
    "wave_asymmetry": {
        "setting": "wave_asymmetry",
        "label": "Wave asymmetry",
        "units": "",
        "values": [0.6, 0.65, 0.7, 0.75],
    },
    "wave_angle_high_fraction": {
        "setting": "wave_angle_high_fraction",
        "label": "Wave angle high fraction",
        "units": "",
        "values": [0.1, 0.15, 0.2, 0.25, 0.3],
    },

    # NOT wave physics, and swept for a different reason. 20 m was decided on
    # 2026-09-01 because 30 m drowned NC-12 at GIS 11 in all eight 1984-2004
    # reloc arms, and hat_run.yaml records that 20 clears that threshold by ONE
    # 10 m cell. A number chosen with one cell of margin needs its margin
    # measured rather than asserted. 0 is the GIS 85/86 ratchet in its purest
    # form; MEASURED is the pre-2026-08-31 behaviour as it actually was.
    #
    # RUN WITH THE HISTORICAL EVENTS ON. The drowning being asked about is a
    # prescribed 1999 DISPLACEMENT landing on top of the emergent target, so
    # with the events off this axis would measure something else and report a
    # reassuring flat line. Period 1 only: no event falls in 2004-2024, and
    # there the axis moves emergent relocations alone.
    "relocation_setback": {
        "setting": "relocation_setback_m",
        "label": "Relocation target",
        "units": "m behind the dune line",
        "values": [0.0, 10.0, 20.0, 30.0, 40.0, 60.0, MEASURED],
        "arm": relocation_arm,
    },
}


def normalise(value):
    """One swept value as the run will see it: a float, or None for `measured`.

    Comparison against the default has to happen on this, not on what was
    typed. 20, "20" and 20.0 are one cell and MEASURED and None are one cell,
    but `20 == "20"` and `"measured" == None` are both False, and either would
    run a duplicate of the baseline under a name that hides which one it is.
    """
    if value is None or (isinstance(value, str)
                         and value.strip().lower() in ("", "none", MEASURED)):
        return None
    return float(value)


def as_environment_value(value):
    """One swept value, spelled for the environment."""
    return MEASURED if normalise(value) is None else repr(normalise(value))


def build_environment(start_year, sweep, value, args):
    """The environment one sweep cell runs under.

    Mirrors `HAT_run_all.run_once`, including HAT_IGNORE_SETTINGS: without it
    the settings this does not set would come from whatever experiment was last
    left in `hat_run.yaml`, and every cell would silently inherit it. Stray
    HAT_* variables in the calling shell are dropped for the same reason.

    Args:
        start_year: Period start, a key of HATTERAS_PERIODS.
        sweep: One entry of SWEEPS.
        value: The value for this cell.
        args: Parsed CLI arguments.

    Returns:
        The environment dict for subprocess.
    """
    environment = {key: val for key, val in os.environ.items()
                   if not key.startswith("HAT_")}
    environment.update({
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(start_year),
        "HAT_SOURCE_SINK_PRESET": args.preset,
        "HAT_SCENARIO": args.scenario,
        "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "true" if args.groin else "false",
        "HAT_OVERWRITE": "true" if args.overwrite else "false",
        "HAT_MAKE_GIFS": "false",     # 30+ cells x 4 GIFs is files nobody reads
        "HAT_SAVE_MODEL_STATE": "false",
        "MPLBACKEND": "Agg",
    })
    # The axis's own arm first, then the swept value. In this order, so an axis
    # can move the arm it is measured in but cannot overwrite the one value
    # that defines the cell.
    environment.update(arm_for(sweep, start_year, args.end_year))
    environment[env_name(sweep)] = as_environment_value(value)
    return environment


def env_name(sweep):
    """The environment variable HAT_hindcast_config reads for this axis."""
    return "HAT_" + sweep["setting"].upper()


def arm_for(sweep, start_year, end_year):
    """The environment overrides this axis needs for this period, possibly none."""
    arm = sweep.get("arm")
    return {} if arm is None else arm(start_year, end_year)


def cells_for(sweep):
    """The values of a sweep that are worth running.

    A value equal to the calibration default produces no name token, so the run
    would derive -- and collide with -- the matrix run's own directory. That
    cell is not a sensitivity result anyway: it IS the baseline, and reading it
    from there is both free and more honest than re-running it under a name
    that hides which one it is.

    Args:
        sweep: One entry of SWEEPS.

    Returns:
        (runnable_values, skipped_default) where skipped_default is the
        calibration value if it appears in the sweep, else None.
    """
    default = normalise(field_default(sweep["setting"]))
    runnable = [v for v in sweep["values"] if normalise(v) != default]
    present = any(normalise(v) == default for v in sweep["values"])
    return runnable, (default if present else None)


def run_cell(start_year, sweep, value, args):
    """Runs one sweep cell and returns its outcome.

    Args:
        start_year: Period start year.
        sweep: One entry of SWEEPS.
        value: The parameter value for this cell.
        args: Parsed CLI arguments.

    Returns:
        A dict recording the cell, suitable for the manifest.
    """
    environment = build_environment(start_year, sweep, value, args)
    row = dict(setting=sweep["setting"], value=value, ok=True, seconds=0.0,
               detail="dry-run")
    if args.dry_run:
        print(f"      would run {env_name(sweep)}={environment[env_name(sweep)]}")
        return row

    started = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, str(HINDCAST)], env=environment,
        cwd=str(PROJECT_BASE_DIR), capture_output=True, text=True)
    seconds = time.perf_counter() - started
    ok = completed.returncode == 0
    if not ok:
        tail = [line for line in completed.stdout.splitlines()[-40:]
                if line.strip()]
        detail = "\n".join(tail[-6:]) or completed.stderr[-400:]
        print(f"      FAILED exit {completed.returncode}\n{detail}")
    else:
        # The run prints where it landed; echoing it makes the sweep log a map
        # from parameter value to directory without a second convention.
        landed = [line for line in completed.stdout.splitlines()
                  if line.startswith("done ")]
        detail = landed[-1].split(None, 1)[1].strip() if landed else ""
        print(f"      ok  {seconds / 60:.1f} min  {Path(detail).name}")
    row.update(ok=ok, seconds=seconds, detail=detail)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--start-year", type=int, default=1984,
                        choices=sorted(HATTERAS_PERIODS),
                        help="hindcast period to sweep")
    parser.add_argument("--param", default="wave_height",
                        choices=sorted(SWEEPS) + ["all"],
                        help="which axis to sweep (default wave_height)")
    parser.add_argument("--values", default=None,
                        help="comma-separated override for the swept values")
    parser.add_argument("--preset", default="calibBE",
                        help="source/sink preset the cells run under")
    parser.add_argument("--scenario", default="full_management")
    parser.add_argument("--groin", action="store_true", default=True)
    parser.add_argument("--no-groin", dest="groin", action="store_false")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    names = sorted(SWEEPS) if args.param == "all" else [args.param]
    if args.values and args.param == "all":
        parser.error("--values applies to one axis; name it with --param")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = OUT_ROOT / f"sensitivity_{args.start_year}.jsonl"

    args.end_year = HATTERAS_PERIODS[args.start_year].get(
        "end_year", args.start_year + 20)
    print(f"period      {args.start_year}-{args.end_year}")
    print(f"baseline    {args.preset} / {args.scenario} / "
          f"groin {'on' if args.groin else 'off'}")
    print(f"manifest    {manifest}")

    planned = []
    for name in names:
        sweep = dict(SWEEPS[name])
        if args.values:
            sweep["values"] = [v.strip() for v in args.values.split(",")
                               if v.strip()]
        runnable, skipped = cells_for(sweep)
        planned.append((name, sweep, runnable, skipped))
        print(f"\n{sweep['label']:26s} {len(runnable)} cells  "
              f"{runnable}  {sweep['units']}")
        if skipped is not None:
            print(f"{'':26s} {skipped} is the calibration value -- that cell "
                  f"IS the baseline, not re-run here")
        for key, val in arm_for(sweep, args.start_year, args.end_year).items():
            print(f"{'':26s} arm: {key}={val}")

    total = sum(len(r) for _, _, r, _ in planned)
    print(f"\n{total} cells at roughly 2-6 min each\n")

    results = []
    for name, sweep, runnable, _ in planned:
        print(f"==== {sweep['label']}")
        for index, value in enumerate(runnable, 1):
            print(f"  [{index:2d}/{len(runnable)}] {sweep['setting']} = {value}")
            row = run_cell(args.start_year, sweep, value, args)
            row.update(sweep=name, start_year=args.start_year,
                       preset=args.preset, scenario=args.scenario,
                       groin=args.groin,
                       arm=arm_for(sweep, args.start_year, args.end_year))
            results.append(row)
            if not args.dry_run:
                with manifest.open("a", encoding="utf-8") as handle:
                    handle.write(json.dumps(row) + "\n")

    failed = [r for r in results if not r["ok"]]
    print(f"\n{len(results) - len(failed)} ok, {len(failed)} failed")
    for row in failed:
        print(f"  {row['setting']}={row['value']}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
