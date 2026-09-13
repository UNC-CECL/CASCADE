#!/usr/bin/env python3
r"""
HAT_run_row_insert_set.py
==============================================================================
Runs the 1984-2004 hindcast once per FILL of the 1984 seaward-row insert, so
the fills can be judged on model behaviour, island-wide, under identical
settings. Nothing but the topography differs between the arms.

THE SET (decided with Hannah, 2026-09-04; DISMANTLED 2026-09-07)
    The six-fill comparison ran on 2026-09-04 and answered its question: the
    fill does not matter island-wide (interior RMSE 0.540-0.547 across all
    seven arms) and the insert itself moves every emergent relocation 3-10
    years late. On 2026-09-07 Hannah decided to keep ONLY unmodified
    topography: dune-topo v3-v8 (the layers), every run made on them
    (raw_runs/row-insert/, blocksv4*, islandv5) and
    output/experiments/row_insert_set/ were deleted. Sizes and reasons:
    data/hatteras_init/1-barrier3d-domains/archive_purge_20260907.csv.

    What is still runnable here is the pair of UNMODIFIED arms:
    arm               topography                  the inserted rows are
    original          v1                          absent; the ORIGINAL pick set with its
                                                  v1-era setback CSV (reference)
    none              v2                          absent; the re-pick base, CURRENT;
                                                  GIS 85/86 setback 0

    The retired arms, for reading old metadata (as-built numbers differ, see
    the guide): measured-floor v4, median v5, platform v6, matched-crest v7,
    matched-nocrest v8. v4-v8 shared one footprint (98 rows, 38 domains) and
    one setback CSV. Versions were a plain sequence (v1 original picks, v2
    re-pick base, v3-v8 layers); the guide is
    1984-start/2-domain-reconstruction-1984/DUNE_TOPO_VERSION_GUIDE.md.

SETTINGS
    Exactly the calibration-tree run HAT_1984_2004_calibBE_road_bdm_groin:
    full_management, calibBE, groin on at the code defaults (M 60, f 0.6),
    prescribed relocations inherited from the scenario (off), calibration
    wave climate. hat_run.yaml is ignored (HAT_IGNORE_SETTINGS=1) so a
    half-edited settings file cannot reach the set.

WHERE THE RUNS LAND
    output/raw_runs/row-insert/<arm>/1984_2004/calibBE/<run_name>/

    via HAT_ARM_TAG="row-insert/<arm>" -- a two-level arm, allowed by
    run_registry.arm_component since 2026-09-04 so a set sits under one folder
    named for it. The run name is the same in every arm (the switches are
    identical by design); the arm is what tells them apart, on disk and in the
    `arm` column of run_index.csv.

WHY IT IS A SCRIPT
    Two pieces of GLOBAL state select an arm and both must be put back:
      * the forcing-tree setback CSV, which hatteras_site_config.py hardcodes
        (copied per arm from the version's own folder, restored in `finally`)
      * the topography version, selected through HAT_TOPO_VERSION_1984_START,
        which outranks CURRENT and dies with the subprocess.
    dune-topo/CURRENT is never written.

USAGE
    python HAT_run_row_insert_set.py --dry-run
    python HAT_run_row_insert_set.py                      # both unmodified arms
    python HAT_run_row_insert_set.py --arms median,platform --overwrite
    python HAT_run_row_insert_set.py --relocations 1     # the arm-B partners
==============================================================================
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Anchored by SEARCHING UPWARD for the project root rather than by
# counting parent directories (2026-09-13). A counted depth is correct
# only while the file stays where it was written, and these moved into
# subfolders of hatteras_ms. Six files here already did it this way.
REPO = next(_p for _p in HERE.parents if (_p / 'pyproject.toml').exists())
sys.path.insert(0, str(REPO / "scripts"))
from cascade_pipeline.run_registry import arm_component  # noqa: E402

# The runner stays at the top of hatteras_ms; this driver moved into
# experiments/ on 2026-09-13, so it names the folder rather than its own.
HINDCAST = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
DUNE_TOPO = REPO / "data/hatteras_init/1-barrier3d-domains/1984-start/dune-topo"
LIVE_SETBACK = (REPO / "data/hatteras_init/4-mgmt-forcing/road_offset"
                / "dunestart_offset/1984/RoadSetback_1984_dunestart.csv")
SET = "row-insert"
LOG_DIR = REPO / "output" / "experiments" / "row_insert_set" / "logs"

# arm -> topography version. Ordered: the control first, then the fills in
# the order the fill report discusses them.
ARMS = {
    "original":        "v1",      # the original pick set + its v1-era setbacks (reference)
    "none":            "v2",
    # measured-floor v4, median v5, platform v6, matched-crest v7 and
    # matched-nocrest v8 were RETIRED 2026-09-07: the layers were deleted
    # (see the docstring). A run on them can no longer be built.
}

BASE_ENV = {
    "HAT_IGNORE_SETTINGS": "1",
    "HAT_START_YEAR": "1984",
    "HAT_SOURCE_SINK_PRESET": "calibBE",
    "HAT_SCENARIO": "full_management",
    "HAT_GROIN_ENABLED": "1",
    "HAT_SHOW_FIGURES": "0",
    "HAT_SAVE_MODEL_STATE": "1",     # the comparison reads roadway objects
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--arms", default=",".join(ARMS),
                    help="comma list, default both unmodified arms")
    ap.add_argument("--overwrite", action="store_true",
                    help="empty an arm's existing run directory and re-run it; "
                         "without this an existing result stops that arm")
    ap.add_argument("--no-gifs", action="store_true",
                    help="skip the shoreline animations (minutes per run)")
    ap.add_argument("--relocations", choices=("inherit", "1"), default="inherit",
                    help="inherit: the scenario decides (full_management -> off), the "
                         "set as first run. 1: force the prescribed 1989/1999 events ON; "
                         "the run name gains a `reloc` token and lands in the same arm "
                         "folder, as arm B of HAT_relocation_comparison.py")
    args = ap.parse_args()

    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    for a in arms:
        if a not in ARMS:
            raise SystemExit("unknown arm {!r}; known: {}".format(a, list(ARMS)))
        version = ARMS[a]
        for needed in ("topography", "RoadSetback_1984_dunestart.csv"):
            if not (DUNE_TOPO / version / needed).exists():
                raise SystemExit(
                    "\narm {!r} needs {} in dune-topo/{}, which is missing.\n"
                    "Build it with HAT_insert_seaward_rows.py (the layers) or save "
                    "the road tree's CSV into it (v2).\n".format(a, needed, version))
        arm_component("{}/{}".format(SET, a))     # raises if the tag is malformed

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = LOG_DIR / stamp
    log_dir.mkdir(parents=True, exist_ok=True)
    backup = log_dir / LIVE_SETBACK.name
    shutil.copy2(LIVE_SETBACK, backup)
    print("live setback CSV backed up to {}".format(backup))

    results = {}
    try:
        for arm in arms:
            version = ARMS[arm]
            tag = "{}/{}".format(SET, arm)
            print("\n" + "=" * 78)
            print("ARM {:16s} topo {:4s} -> raw_runs/{}/1984_2004/calibBE/".format(
                arm, version, tag))
            print("=" * 78)
            shutil.copy2(DUNE_TOPO / version / LIVE_SETBACK.name, LIVE_SETBACK)

            env = dict(os.environ)
            env.update(BASE_ENV)
            env["HAT_ARM_TAG"] = tag
            env["HAT_TOPO_VERSION_1984_START"] = version
            env["HAT_OVERWRITE"] = "1" if args.overwrite else "0"
            if args.no_gifs:
                env["HAT_MAKE_GIFS"] = "0"
            if args.relocations == "1":
                env["HAT_RELOCATIONS"] = "1"
            if args.dry_run:
                print("  [dry-run] would run {}\n    with {}".format(
                    HINDCAST.name, {k: env[k] for k in sorted(env)
                                    if k.startswith("HAT_")}))
                continue
            proc = subprocess.run([sys.executable, str(HINDCAST)],
                                  cwd=str(HERE), env=env,
                                  capture_output=True, text=True)
            log = log_dir / "{}.log".format(arm)
            log.write_text((proc.stdout or "") + "\n--- STDERR ---\n"
                           + (proc.stderr or ""), encoding="utf-8")
            results[arm] = proc.returncode
            print("  exit {}   log: {}".format(proc.returncode, log))
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-25:]
                print("\n".join("    " + t for t in tail))
    finally:
        # ALWAYS put the forcing tree back, whatever arm was last.
        shutil.copy2(backup, LIVE_SETBACK)
        print("\nrestored the live setback CSV from {}".format(backup))

    if results:
        print("\n" + "  ".join("{}={}".format(k, v) for k, v in results.items()))


if __name__ == "__main__":
    main()
