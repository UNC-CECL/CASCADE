#!/usr/bin/env python3
r"""
HAT_run_version_pair.py
==============================================================================
Run ONE hindcast scenario on two dune-topo versions, identically, so the two
runs differ in nothing but the topography and its setback CSV. The pair the
advisor asked for (2026-09-09): v2 (the extraction) against v3 (the 1984
reconstruction), full management, calibBE, groin on, under the modules'
automatic behaviour - and the same pair with the recorded 1989/1999
relocations prescribed, as the control.

WHY A SCRIPT (the same reason HAT_run_row_insert_set.py is one). Two pieces
of global state select a version and both must be put back: the forcing-tree
setback CSV, which hatteras_site_config.py hardcodes, is copied per version
from dune-topo/<version>/ and restored in `finally`; the topography version
goes through HAT_TOPO_VERSION_1984_START, which outranks CURRENT and dies
with the subprocess. hat_run.yaml is ignored (HAT_IGNORE_SETTINGS=1).

WHERE THE RUNS LAND
    output/raw_runs/version-pair/<version>/1984_2004/calibBE/<run_name>/
    via HAT_RUN_KIND=version HAT_RUN_TAG="version-pair/<version>", so the
    runs file under raw_runs/versions/version-pair/<version>/ (2026-09-16).
    The run name is the same for both versions by design; the tag tells them
    apart, on disk and in the `kind`/`tag` columns of run_index.csv.

USAGE
    python HAT_run_version_pair.py --relocations 1          # the prescribed control
    python HAT_run_version_pair.py                          # emergent (the modules decide)
    python HAT_run_version_pair.py --versions v2,v3 --dry-run
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


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))
from cascade_pipeline.run_registry import arm_component  # noqa: E402

HINDCAST = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
DUNE_TOPO = REPO / "data/hatteras_init/1-barrier3d-domains/1984-start/dune-topo"
from site_layer import hat_topo_version as _tv  # noqa: E402
LIVE_SETBACK = _tv.road_setback_file(1984)
SET = "version-pair"
LOG_DIR = REPO / "output" / "experiments" / "version_pair" / "logs"

BASE_ENV = {
    "HAT_IGNORE_SETTINGS": "1",
    "HAT_START_YEAR": "1984",
    "HAT_SOURCE_SINK_PRESET": "calibBE",
    "HAT_SCENARIO": "full_management",
    "HAT_GROIN_ENABLED": "1",
    "HAT_SHOW_FIGURES": "0",
    "HAT_SAVE_MODEL_STATE": "1",     # the comparison reads the saved objects
    "HAT_MAKE_GIFS": "0",
    "MPLBACKEND": "Agg",
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--versions", default="v2,v3")
    ap.add_argument("--relocations", choices=("inherit", "1"), default="inherit",
                    help="inherit: the scenario decides (full_management -> off, the modules act on their own); "
                         "1: the recorded 1989/1999 events prescribed (the run name gains `reloc`)")
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    versions = [v.strip() for v in a.versions.split(",") if v.strip()]
    for v in versions:
        for needed in ("topography", "RoadSetback_1984_dunestart.csv"):
            if not (DUNE_TOPO / v / needed).exists():
                raise SystemExit(f"\n{v} has no {needed} under dune-topo/{v}\n")
        arm_component(f"{SET}/{v}")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_dir = LOG_DIR / stamp
    log_dir.mkdir(parents=True, exist_ok=True)
    backup = log_dir / LIVE_SETBACK.name
    shutil.copy2(LIVE_SETBACK, backup)
    print(f"live setback CSV backed up to {backup}")
    results = {}
    try:
        for v in versions:
            tag = f"{SET}/{v}"
            print("\n" + "=" * 78)
            print(f"VERSION {v}  relocations={a.relocations}  -> raw_runs/versions/{tag}/1984_2004/calibBE/")
            print("=" * 78)
            shutil.copy2(DUNE_TOPO / v / LIVE_SETBACK.name, LIVE_SETBACK)
            env = {k: val for k, val in os.environ.items() if not k.startswith("HAT_")}
            env.update(BASE_ENV)
            env["HAT_RUN_KIND"] = "version"
            env["HAT_RUN_TAG"] = tag
            env["HAT_TOPO_VERSION_1984_START"] = v
            env["HAT_OVERWRITE"] = "1" if a.overwrite else "0"
            if a.relocations == "1":
                env["HAT_RELOCATIONS"] = "1"
            if a.dry_run:
                print("  [dry-run]", {k: env[k] for k in sorted(env) if k.startswith("HAT_")})
                continue
            proc = subprocess.run([sys.executable, str(HINDCAST)], cwd=str(HINDCAST.parent), env=env,
                                  capture_output=True, text=True)
            log = log_dir / f"{v}_reloc{a.relocations}.log"
            log.write_text((proc.stdout or "") + "\n--- STDERR ---\n" + (proc.stderr or ""), encoding="utf-8")
            results[v] = proc.returncode
            print(f"  exit {proc.returncode}   log: {log}")
            if proc.returncode != 0:
                print("\n".join("    " + t for t in (proc.stderr or proc.stdout or "").strip().splitlines()[-25:]))
    finally:
        shutil.copy2(backup, LIVE_SETBACK)
        print(f"\nrestored the live setback CSV from {backup}")
    if results:
        print("  ".join(f"{k}={v}" for k, v in results.items()))


if __name__ == "__main__":
    main()
