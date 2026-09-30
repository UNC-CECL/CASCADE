#!/usr/bin/env python3
"""
Run the 1984-2004 hindcast once per crest-experiment arm.

    python scripts/hatteras_ms/experiments/HAT_run_crest_experiment.py [--dry-run] [--arms a,b]

Each arm's setback CSV is swapped in and always put back; the topography comes from
HAT_TOPO_VERSION_1984_START. Only pea1989base remains. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
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
# Repo root, found by searching upward
REPO = next(_p for _p in HERE.parents if (_p / 'pyproject.toml').exists())
# The runner sits at the top of hatteras_ms
HINDCAST = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"

import sys as _b3dsys
from pathlib import Path as _B3DP
_b3dsys.path.insert(0, str(next(_q for _q in _B3DP(__file__).resolve().parents
                                if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _b3d  # noqa: E402
DUNE_TOPO = _b3d.dune_topo_root("1984-start")
CURRENT = DUNE_TOPO / "CURRENT"
import sys as _tvsys
from pathlib import Path as _TVP
_tvsys.path.insert(0, str(next(_q for _q in _TVP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_topo_version as _tv  # noqa: E402
LIVE_SETBACK = _tv.road_setback_file(1984)

# arm tag -> (topo version, setback CSV source; None = leave the live one)
def _arm(version):
    return (version, DUNE_TOPO / version / "RoadSetback_1984_dunestart.csv")


# --- CONFIG ------------------------------------------------------------------
ARMS = {
    # Only the baseline arm remains; None leaves the live setback CSV in place
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
# -----------------------------------------------------------------------------


# Run: every chosen arm, restoring the live setback file after each
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
        # Fail early if an arm's version has left dune-topo/
        version = ARMS[a][0]
        if not (DUNE_TOPO / version).is_dir():
            retired = DUNE_TOPO.parent / "1-extraction" / "dune-topo-experiments" / version
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
    backup = REPO / "output" / "archive" / "2026-09-07_crest-insert-arm-logs" / stamp
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

            # CURRENT is deliberately not written: it loses to the extractor literal
            if setback_src is not None:
                shutil.copy2(setback_src, LIVE_SETBACK)
            else:
                shutil.copy2(backup / LIVE_SETBACK.name, LIVE_SETBACK)

            env = dict(os.environ)
            env.update(BASE_ENV)
            env["HAT_RELOCATIONS"] = args.relocations
            arm_tag = arm if args.relocations == "1" else arm + "noreloc"
            # Filed as an experiment, the member being the arm name without pea1989
            env["HAT_RUN_KIND"] = "experiment"
            env["HAT_RUN_TAG"] = "topography-and-domains/2026-09-02-pea-island-row-insert-control/" + arm_tag.replace("pea1989", "", 1)
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
        # Always put the tree back, or later runs read a fabricated topography
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
