#!/usr/bin/env python3
r"""
HAT_run_crest_experiment.py
==============================================================================
Runs the 1984-2004 hindcast three times -- baseline, insert-with-crest-kept,
insert-with-crest-shaved -- so the GIS 84/85/86 row insert can be judged on
model behaviour rather than on cross-sections.

WHY IT IS A SCRIPT AND NOT THREE COMMANDS
    The road-setback CSV is GLOBAL state: `hatteras_site_config.py:142`
    hardcodes its path, so selecting an arm means overwriting a file every other
    reader of this repo also sees. Doing that by hand leaves the tree pointing at
    an experiment arm the moment anything goes wrong. Here the restore is in a
    `finally`.

    THE TOPOGRAPHY IS NOT SELECTED THAT WAY, and the first version of this
    script got it wrong. It wrote `dune-topo/CURRENT` per arm -- and every arm
    still ran on v1, because `resolve_version` reads the EXTRACTOR's VERSION
    literal before it reads CURRENT. Two of three arms were silent duplicates of
    the control. The version now comes from HAT_TOPO_VERSION_1984_START, which
    outranks the extractor, is scoped to one product, and dies with the
    subprocess instead of persisting.

WHY ARMS RATHER THAN RUN NAMES
    All three produce the SAME run name -- the name is derived from the
    management switches, and those are identical by design. `HAT_ARM_TAG` gives
    each its own directory component, and the run index is keyed on
    ("run_name", "Hs_m", "arm"), so nothing overwrites anything. The existing
    calibration-tree run is never touched.

RELOCATIONS ON OR OFF -- TWO DIFFERENT QUESTIONS
    --relocations 1 asks: does an emergent relocation fire BEFORE the prescribed
        1989 Pea Island event and corrupt the base its displacement is added to?
        That is a question about whether the prescribed history is applied to the
        right island.

    --relocations 0 asks: left to itself, WHEN does the module relocate the road,
        and how close is that to 1989? That is a question about model skill, and
        it is the one you cannot ask with the events switched on -- prescribing
        the 1989 relocation and then checking whether the model produces 1989 is
        circular.

    The second needs a road that starts BEHIND the dune. A setback floored to 0
    relocates in year 1 by construction and carries no information about
    anything, which is why the baseline arm is a control here and not a
    prediction.

ARMS
    islandv5       v5  island-wide, every domain with a measured shift
    pea1989base    v1              + the setback CSV as shipped
    (pea1989keep / pea1989lower were removed 2026-09-03 with their
     topography; their run outputs remain under output/raw_runs/)

USAGE
    python HAT_run_crest_experiment.py [--dry-run] [--arms a,b]
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
REPO = HERE.parents[1]
HINDCAST = HERE / "HAT_hindcast_1984_2024.py"

DUNE_TOPO = REPO / "data/hatteras_init/1-barrier3d-domains/1984-start/dune-topo"
CURRENT = DUNE_TOPO / "CURRENT"
LIVE_SETBACK = (REPO / "data/hatteras_init/4-mgmt-forcing/road_offset"
                / "dunestart_offset/1984/RoadSetback_1984_dunestart.csv")

# arm tag -> (topo version, setback CSV source; None = leave the live one)
def _arm(version):
    return (version, DUNE_TOPO / version / "RoadSetback_1984_dunestart.csv")


ARMS = {
    # --- ISLAND-WIDE SCOPE ---------------------------------------------------
    # v5 = v3 + rows at EVERY domain whose measured 1984->1997 shift rounds to
    # >= 1 cell, road or not (38 domains). Same topography treatment as v4 and
    # identical N at all ten block domains -- the only difference is that the
    # 28 other domains with a measured shift are no longer passed over. Under
    # v4 "unchanged" meant two different things along the island: measured-zero
    # inside the blocks, and never-asked outside them.
    "islandv5": _arm("v5"),
    # --- N-ESTIMATE ARMS (scope: both relocation blocks, GIS 9-14 + 84-87) ----
    # Same topography treatment throughout -- translate, crest shaved, measured
    # fill -- so the ONLY thing varying is where N came from. That is the point:
    # the two measurements of N disagree by a factor of ~3 at GIS 85 (65.9 m vs
    # 19.5 m), and the relocation-timing test is what discriminates them.
    # v4 = v3 (2026-09-02 re-pick) + rows. The one to take forward.
    "blocksv4": _arm("v4"),
    # "blocksdate" (v2 = v1 + rows) was REMOVED on 2026-09-03 with v2 itself.
    # v2 was built on the pre-re-pick v1 extraction and superseded by v4. Its
    # completed run output survives at output/raw_runs/blocksdate{,noreloc};
    # only re-running it is no longer possible.
    # --- the earlier GIS 84-86 crest pair, kept so those runs stay reproducible
    # FIVE ARMS REMOVED 2026-09-03 with their topography variants:
    #   blocksduneline / blocksdsas / blocksminimum - the N-SOURCE
    #     comparison, settled by the 1997 dune line. Both v4 and v5 use
    #     `date`, and the dsas source itself is gone.
    #   pea1989keep / pea1989lower - the crest-shaving pair.
    # Their RUN OUTPUTS survive in output/raw_runs/, so
    # HAT_plot_crest_experiment.py, HAT_score_relocation_timing.py and
    # HAT_score_road_position.py all still work on them; only re-running
    # from topography is no longer possible.
    "pea1989base": ("v1", None),
}

BASE_ENV = {
    "HAT_IGNORE_SETTINGS": "1",     # hat_run.yaml must not reach this experiment
    "HAT_START_YEAR": "1984",
    "HAT_SOURCE_SINK_PRESET": "calibBE",
    "HAT_SCENARIO": "full_management",
    "HAT_GROIN_ENABLED": "1",
    "HAT_SHOW_FIGURES": "0",
    "HAT_MAKE_GIFS": "0",           # not needed, and they cost minutes
    "HAT_SAVE_MODEL_STATE": "1",    # the comparison reads roadway objects
    "HAT_OVERWRITE": "1",           # within this arm's own directory only
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--arms", default=",".join(ARMS))
    ap.add_argument("--relocations", choices=("0", "1"), default="1",
                    help="prescribed HATTERAS_ROAD_EVENTS on/off. "
                         "0 appends _noreloc to each arm tag so the "
                         "two experiments cannot overwrite each other.")
    args = ap.parse_args()
    arms = [a.strip() for a in args.arms.split(",") if a.strip()]
    for a in arms:
        if a not in ARMS:
            raise SystemExit("unknown arm {!r}; known: {}".format(a, list(ARMS)))
        # Most of these versions were retired to dune-topo-experiments/ on
        # 2026-09-02. Nothing outside dune-topo/ resolves through topo_dirs or
        # HAT_TOPO_VERSION_1984_START, so the run would otherwise fail deep in
        # the hindcast with a missing-array error that names no cause.
        version = ARMS[a][0]
        if not (DUNE_TOPO / version).is_dir():
            retired = DUNE_TOPO.parent / "dune-topo-experiments" / version
            msg = ("\narm {!r} needs topography version {!r}, which is not in"
                   "\n  {}\n".format(a, version, DUNE_TOPO))
            if retired.is_dir():
                msg += ("It was retired to\n  {}\nMove it back into dune-topo/ "
                        "to re-run this arm; see that folder's README for why "
                        "the friction is deliberate.\n".format(retired))
            else:
                msg += "and is not in dune-topo-experiments/ either.\n"
            raise SystemExit(msg)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = REPO / "output" / "experiments" / "crest_experiment_restore" / stamp
    backup.mkdir(parents=True, exist_ok=True)
    saved_current = CURRENT.read_text(encoding="utf-8") if CURRENT.is_file() else None
    shutil.copy2(LIVE_SETBACK, backup / LIVE_SETBACK.name)
    print("restore copies in {}".format(backup))

    results = {}
    try:
        for arm in arms:
            version, setback_src = ARMS[arm]
            print("\n" + "=" * 78)
            print("ARM {}  |  topo {}  |  setback {}".format(
                arm, version, setback_src.name if setback_src else "as shipped"))
            print("=" * 78)

            # NOTE: CURRENT is deliberately NOT written. See the docstring --
            # it loses to the extractor literal, so writing it here would look
            # like arm selection while doing nothing.
            if setback_src is not None:
                shutil.copy2(setback_src, LIVE_SETBACK)
            else:
                shutil.copy2(backup / LIVE_SETBACK.name, LIVE_SETBACK)

            env = dict(os.environ)
            env.update(BASE_ENV)
            env["HAT_RELOCATIONS"] = args.relocations
            arm_tag = arm if args.relocations == "1" else arm + "noreloc"
            env["HAT_ARM_TAG"] = arm_tag
            env["HAT_TOPO_VERSION_1984_START"] = version
            if args.dry_run:
                print("  [dry-run] would run {}".format(HINDCAST))
                continue
            proc = subprocess.run([sys.executable, str(HINDCAST)],
                                  cwd=str(HERE), env=env,
                                  capture_output=True, text=True)
            log = backup / "{}.log".format(arm_tag)
            log.write_text((proc.stdout or "") + "\n--- STDERR ---\n"
                           + (proc.stderr or ""), encoding="utf-8")
            results[arm] = proc.returncode
            print("  exit {}   log: {}".format(proc.returncode, log))
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-25:]
                print("\n".join("    " + t for t in tail))
    finally:
        # ALWAYS put the tree back. An experiment arm left in CURRENT would make
        # every later road measurement and every later run silently read a
        # fabricated topography.
        if saved_current is not None:
            CURRENT.write_text(saved_current, encoding="utf-8")
        elif CURRENT.is_file():
            CURRENT.unlink()
        shutil.copy2(backup / LIVE_SETBACK.name, LIVE_SETBACK)
        print("\nrestored CURRENT = {!r} and the shipped setback CSV".format(
            (saved_current or "").strip()))

    print("\n" + "  ".join("{}={}".format(k, v) for k, v in results.items()))


if __name__ == "__main__":
    main()
