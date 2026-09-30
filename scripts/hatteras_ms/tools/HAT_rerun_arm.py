"""
Re-run a set of existing runs under today's code, into an arm, so the two can be compared.

    python scripts/hatteras_ms/tools/HAT_rerun_arm.py --period 1996 --tag <arm>

Reads each run's settings from the index and runs the hindcast with them;
HAT_compare_rerun.py compares the results. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve()
REPO = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
# --- CONFIG ------------------------------------------------------------------
RUNNER = REPO / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
INDEX = REPO / "output" / "raw_runs" / "run_index.csv"

# The switches a run name encodes, recovered from the index columns rather than parsed out of the name
SWITCH_COLUMNS = {
    "roadway_management": "HAT_SCENARIO",         # resolved below
    "beach_dune_management": None,
    "nourishment_fills": None,
    "relocations_enabled": "HAT_RELOCATIONS",
    "groin_enabled": "HAT_GROIN_ENABLED",
}
# -----------------------------------------------------------------------------


# The named scenario matching this run's four management switches
def scenario_for(row):
    road = bool(row["roadway_management"])
    bdm = bool(row["beach_dune_management"])
    fills = bool(row["nourishment_fills"])
    if road and bdm:
        return "full_management" if fills else "full_no_fill"
    if road:
        return "roadway_only"
    if bdm:
        return "beachdune_only"
    return "natural"


# The relocation arm of one period, matrix cells only
def select(index, period):
    d = index[(index["start_year"] == period)
              & (index["relocations_enabled"] == True)      # noqa: E712
              & (index["kind"] == "matrix")]
    # the rset sweep varies the rebuild clearance; it is a separate experiment
    d = d[~d["run_name"].str.contains("_rset")]
    return d.drop_duplicates("run_name")


# The environment that reproduces one indexed run
def env_for(row, tag, topo_version):
    env = dict(os.environ)
    env.update({
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(int(row["start_year"])),
        "HAT_SOURCE_SINK_PRESET": str(row["source_sink_preset"]),
        "HAT_SCENARIO": scenario_for(row),
        "HAT_RELOCATIONS": "true",
        "HAT_GROIN_ENABLED": "true" if row["groin_enabled"] else "false",
        "HAT_HS": str(row["Hs_m"]),
        # An EXPERIMENT, filed under raw_runs/experiments/<tag>/ (2026-09-16).
        "HAT_RUN_KIND": "experiment",
        "HAT_RUN_TAG": tag,
        "HAT_MAKE_GIFS": "false",
        "HAT_SAVE_MODEL_STATE": "false",
        "HAT_OVERWRITE": "true",
        "PYTHONIOENCODING": "utf-8",
        "MPLBACKEND": "Agg",
    })
    if topo_version:
        env["HAT_TOPO_VERSION_1984_START"] = topo_version
    return env


# Run: every selected run, into the arm
def main():
    ap = argparse.ArgumentParser(description="re-run an arm under today's code")
    ap.add_argument("--period", type=int, default=1984)
    ap.add_argument("--tag", "--arm", dest="tag", default="code-checks/2026-09-14-relocation-arm-rerun-new-code",
                    help="experiment tag the re-runs are filed under")
    ap.add_argument("--topo-version", default="v1",
                    help="pin the 1984-start version. Default v1, which is "
                         "what the stored runs used -- so the difference is "
                         "CODE only, not code and island together.")
    ap.add_argument("--limit", type=int)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    index = pd.read_csv(INDEX)
    if "kind" not in index.columns:          # a pre-09-16 index
        index["kind"] = index["arm"].fillna("calibration").map(
            lambda a: "matrix" if a == "calibration" else "other")
    todo = select(index, args.period)
    if args.limit:
        todo = todo.head(args.limit)

    print(f"{len(todo)} run(s) in the {args.period} relocation arm")
    for _, row in todo.iterrows():
        print(f"  {row['run_name']:<50} {scenario_for(row):<16} "
              f"groin={bool(row['groin_enabled'])}")
    if args.list:
        return 0

    log_dir = REPO / "output" / "logs" / "driver" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    for n, (_, row) in enumerate(todo.iterrows(), 1):
        name = row["run_name"]
        log = log_dir / f"rerun_{args.tag.replace('/', '_')}_{name}.log"
        print(f"\n[{n}/{len(todo)}] {name}", flush=True)
        with open(log, "w", encoding="utf-8") as handle:
            result = subprocess.run(
                [sys.executable, str(RUNNER)], cwd=str(RUNNER.parent),
                env=env_for(row, args.tag, args.topo_version),
                stdout=handle, stderr=subprocess.STDOUT)
        print(f"      exit {result.returncode}   log {log.name}", flush=True)
    print(f"\ndone in {(time.time() - started) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
