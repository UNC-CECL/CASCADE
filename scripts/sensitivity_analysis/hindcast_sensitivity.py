#!/usr/bin/env python3
"""
Run the hindcast once per sweep cell, with one setting moved off its calibrated value.

    python scripts/sensitivity_analysis/hindcast_sensitivity.py --start-year 1996 --param wave_height
    python scripts/sensitivity_analysis/hindcast_sensitivity.py --start-year 1996 --param all --dry-run

Each cell is an ordinary run of HAT_hindcast_1984_2024.py with one HAT_ setting
changed through the environment, filed under output/raw_runs/sensitivity/ and
logged to output/calibration/sensitivity/sensitivity_<start>.jsonl for
plot_sensitivity.py. Details: scripts/sensitivity_analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-28
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
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts",
              PROJECT_BASE_DIR / "scripts" / "hatteras_ms"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_PERIODS, HATTERAS_ROAD_EVENTS)
from cascade_pipeline.roadway import RelocationEvent  # noqa: E402
from HAT_hindcast_config import field_default  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_BASE_DIR / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
OUT_ROOT = PROJECT_BASE_DIR / "output" / "calibration" / "sensitivity"

# Spelled out: an empty string reads as unset and would run the default
MEASURED = "measured"
# -----------------------------------------------------------------------------


# The arm the relocation-target axis has to be measured in, for a period
def relocation_arm(start_year, end_year):
    events = [event for event in HATTERAS_ROAD_EVENTS
              if isinstance(event, RelocationEvent)
              and start_year <= event.year <= end_year]
    return {"HAT_RELOCATIONS": "true"} if events else {}


# The swept axes: `setting` is the HAT_hindcast_config field; `arm` adds per-period overrides
SWEEPS = {
    # Wave axes centred on option A, each with one cell past where the model breaks
    "wave_height": {
        "setting": "hs",
        "label": "Wave height Hs",
        "units": "m",
        "values": [0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0],
    },
    "wave_period": {
        "setting": "wave_period_s",
        "label": "Wave period",
        "units": "s",
        "values": [6.0, 7.0, 7.5, 8.0, 9.0, 10.0, 12.0],
    },
    "wave_asymmetry": {
        "setting": "wave_asymmetry",
        "label": "Wave asymmetry",
        "units": "",
        "values": [0.5, 0.55, 0.6, 0.65, 0.7, 0.8],
    },
    "wave_angle_high_fraction": {
        "setting": "wave_angle_high_fraction",
        "label": "Wave angle high fraction",
        "units": "",
        "values": [0.3, 0.4, 0.45, 0.5, 0.55],
    },

    # Relocation target: its margin over the GIS 11 drowning, run with the historical events on
    "relocation_setback": {
        "setting": "relocation_setback_m",
        "label": "Relocation target",
        "units": "m behind the dune line",
        "values": [0.0, 10.0, 20.0, 30.0, 40.0, 60.0, MEASURED],
        "arm": relocation_arm,
    },
}


# One swept value as the run will see it: a float, or None for `measured`
def normalise(value):
    if value is None or (isinstance(value, str)
                         and value.strip().lower() in ("", "none", MEASURED)):
        return None
    return float(value)


# One swept value, spelled for the environment
def as_environment_value(value):
    return MEASURED if normalise(value) is None else repr(normalise(value))


# The environment one sweep cell runs under, as HAT_run_all.run_once builds it
def build_environment(start_year, sweep, value, args):
    environment = {key: val for key, val in os.environ.items()
                   if not key.startswith("HAT_")}
    environment.update({
        "HAT_IGNORE_SETTINGS": "1",
        # A sweep cell, filed under raw_runs/sensitivity/<axis>/ (2026-09-16).
        "HAT_RUN_KIND": "sensitivity",
        "HAT_START_YEAR": str(start_year),
        "HAT_SOURCE_SINK_PRESET": args.preset,
        "HAT_SCENARIO": args.scenario,
        "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "true" if args.groin else "false",
        "HAT_OVERWRITE": "true" if args.overwrite else "false",
        "HAT_MAKE_GIFS": "false",     # 30+ cells x 4 GIFs is files nobody reads
        "HAT_SAVE_MODEL_STATE": "false",
        "MPLBACKEND": "Agg",
        # The child prints non-ASCII; force UTF-8 so the capture decodes
        "PYTHONIOENCODING": "utf-8",
    })
    # The axis's arm first, then the swept value, so the arm cannot overwrite the cell
    environment.update(arm_for(sweep, start_year, args.end_year))
    environment[env_name(sweep)] = as_environment_value(value)
    return environment


# The environment variable HAT_hindcast_config reads for this axis
def env_name(sweep):
    return "HAT_" + sweep["setting"].upper()


# The environment overrides this axis needs for this period, possibly none
def arm_for(sweep, start_year, end_year):
    arm = sweep.get("arm")
    return {} if arm is None else arm(start_year, end_year)


# The values of a sweep worth running: the calibration default is skipped
def cells_for(sweep):
    default = normalise(field_default(sweep["setting"]))
    runnable = [v for v in sweep["values"] if normalise(v) != default]
    present = any(normalise(v) == default for v in sweep["values"])
    return runnable, (default if present else None)


# Run one sweep cell and return its outcome for the manifest
def run_cell(start_year, sweep, value, args):
    environment = build_environment(start_year, sweep, value, args)
    row = dict(setting=sweep["setting"], value=value, ok=True, seconds=0.0,
               detail="dry-run")
    if args.dry_run:
        print(f"      would run {env_name(sweep)}={environment[env_name(sweep)]}")
        return row

    started = time.perf_counter()
    completed = subprocess.run(
        [sys.executable, str(HINDCAST)], env=environment,
        cwd=str(PROJECT_BASE_DIR), capture_output=True, text=True,
        encoding="utf-8", errors="replace")
    seconds = time.perf_counter() - started
    ok = completed.returncode == 0
    if not ok:
        tail = [line for line in completed.stdout.splitlines()[-40:]
                if line.strip()]
        detail = "\n".join(tail[-6:]) or completed.stderr[-400:]
        print(f"      FAILED exit {completed.returncode}\n{detail}")
    else:
        # Echo where the run landed, so the sweep log maps value to directory
        landed = [line for line in completed.stdout.splitlines()
                  if line.startswith("done ")]
        detail = landed[-1].split(None, 1)[1].strip() if landed else ""
        print(f"      ok  {seconds / 60:.1f} min  {Path(detail).name}")
    row.update(ok=ok, seconds=seconds, detail=detail)
    return row


# Run: every cell of the chosen sweeps, appending each outcome to the manifest
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
