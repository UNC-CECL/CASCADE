"""
regenerate_all_figures.py
==============================================================================
Redraws every figure published under output/figures/ from its own script, in
the numbered layout's order, then rewrites the index (output/figures/README.md).

    python scripts/figure_making/tools/regenerate_all_figures.py            # everything
    python scripts/figure_making/tools/regenerate_all_figures.py --only 4-model-mechanics
    python scripts/figure_making/tools/regenerate_all_figures.py --list

WHY IT EXISTS
    output/figures/ is gitignored and fed by ~25 scripts, so "are the figures up
    to date?" had no answer short of remembering which script drew what. On
    2026-09-29, 44 of 182 figures were still from 09-17, before option-A waves
    and the 09-28 adoption. This is the answer now: one list of producers, run
    in order, each one's log kept, failures reported at the end rather than
    stopping the rest. Add a producer here when a script starts publishing to
    output/figures/.

    The logs go to output/logs/scratch/figures_<timestamp>/, one per step.
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
S = "scripts/"

# (top folder it fills, label, command). Order = the layout's order.
STEPS = [
    ("1-site", "site maps", [S + "figure_making/island/study_area_figures.py"]),
    ("2-observations", "CoastSat periods", [S + "figure_making/shoreline/plot_coastsat_calibration_periods.py"]),
    ("2-observations", "CoastSat periods, LOWESS 7", [S + "figure_making/shoreline/plot_coastsat_poster.py"]),
    ("2-observations", "shoreline vs dune line", [S + "input_prep/5-scr/4-comparisons/shoreline_vs_duneline/net_change_vs_duneline.py"]),
    ("2-observations", "dune-line positions", [S + "input_prep/5-scr/4-comparisons/duneline_positions/duneline_positions.py"]),
    ("2-observations", "mean shoreline 1996 (DEM-centred)",
     [S + "input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_on_imagery.py",
      "--centred-on", "alace_1996"]),
    ("2-observations", "mean shoreline 2010 (DEM-centred)",
     [S + "input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_on_imagery.py",
      "--centred-on", "usace_2009", "--photo-years", "2008"]),
    ("3-model-inputs", "0 elevation", [S + "figure_making/pipeline/0-elevation/dem_composition_figures.py"]),
    ("3-model-inputs", "1 domain extraction", [S + "figure_making/pipeline/1-barrier3d-domains/domain_extraction_figures.py"]),
    ("3-model-inputs", "1 initial island", [S + "figure_making/island/initialization_figures.py"]),
    ("3-model-inputs", "2 BRIE offset", [S + "figure_making/pipeline/2-brie-offset/offset_build_figures.py"]),
    ("3-model-inputs", "3 storms", [S + "figure_making/pipeline/3-storms/storm_construction_figures.py"]),
    ("3-model-inputs", "4 road setbacks", [S + "figure_making/pipeline/4-mgmt-forcings/road_setback_figures.py"]),
    ("3-model-inputs", "4 management timeline", [S + "figure_making/management/management_timeline_figure.py"]),
    ("3-model-inputs", "4 management rules", [S + "figure_making/management/management_rules_table_figure.py"]),
    ("3-model-inputs", "5 observed target", [S + "figure_making/pipeline/5-scr/observed_target_figures.py"]),
    ("3-model-inputs", "7 source/sink", [S + "figure_making/pipeline/7-source-sink/be_method_figures.py"]),
    ("4-model-mechanics", "model mechanics", [S + "figure_making/model/model_mechanics_figures.py"]),
    ("4-model-mechanics", "storm routing", [S + "figure_making/model/overwash_routing_figures.py"]),
    ("5-results", "hindcast edgeBE", [S + "figure_making/model_output/hindcast_final_figure_lowess.py", "--preset", "edgeBE"]),
    ("5-results", "hindcast zeroBE", [S + "figure_making/model_output/hindcast_final_figure_lowess.py", "--preset", "zeroBE"]),
    ("5-results", "scenario grid", [S + "figure_making/model_output/scenario_grid.py"]),
    ("5-results", "GIS 11 relocation", [S + "figure_making/model_output/gis11_relocation_drown_figure.py"]),
    ("talk", "site maps, projector set", [S + "figure_making/island/study_area_figures.py", "--talk"]),
    ("style", "style sheet", [S + "site_layer/hat_figure_style.py"]),
]
INDEX = [S + "figure_making/tools/figure_index.py"]
FIGURES = REPO / "output" / "figures"
PRODUCERS = FIGURES / "supporting" / "producers.json"   # read by figure_index.py


def png_times() -> dict[str, float]:
    """{path under output/figures: mtime} for every published PNG."""
    return {p.relative_to(FIGURES).as_posix(): p.stat().st_mtime
            for p in FIGURES.rglob("*.png") if "supporting" not in p.parts}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[3])
    ap.add_argument("--only", nargs="+", metavar="FOLDER",
                    help="top folders to redraw, e.g. 4-model-mechanics 5-results")
    ap.add_argument("--list", action="store_true", help="print the producers and stop")
    args = ap.parse_args()
    steps = [s for s in STEPS if not args.only or s[0] in args.only]
    if args.list:
        for folder, label, cmd in steps:
            print(f"{folder:18s} {label:34s} {' '.join(cmd)}")
        return 0
    logs = REPO / "output" / "logs" / "scratch" / f"figures_{datetime.now():%Y%m%d_%H%M%S}"
    logs.mkdir(parents=True, exist_ok=True)
    failed = []
    producers = json.loads(PRODUCERS.read_text(encoding="utf-8")) if PRODUCERS.is_file() else {}
    for i, (folder, label, cmd) in enumerate(steps, 1):
        t0 = time.time()
        before = png_times()
        log = logs / f"{i:02d}_{label.replace(' ', '_').replace('/', '-')}.log"
        with open(log, "w", encoding="utf-8") as fh:
            rc = subprocess.run([sys.executable, *cmd], cwd=REPO, stdout=fh,
                                stderr=subprocess.STDOUT).returncode
        # every PNG this step wrote or rewrote is attributed to its script
        for rel, t in png_times().items():
            if before.get(rel) != t:
                producers[rel] = cmd[0]
        PRODUCERS.parent.mkdir(parents=True, exist_ok=True)
        PRODUCERS.write_text(json.dumps(dict(sorted(producers.items())), indent=1), encoding="utf-8")
        status = "ok" if rc == 0 else f"FAILED (exit {rc})"
        print(f"[{i:2d}/{len(steps)}] {folder:18s} {label:34s} {status}  {time.time() - t0:5.0f}s",
              flush=True)
        if rc:
            failed.append((label, log))
    subprocess.run([sys.executable, *INDEX], cwd=REPO)
    print(f"\nlogs: {logs.relative_to(REPO)}")
    for label, log in failed:
        print(f"FAILED: {label}  ->  {log.relative_to(REPO)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
