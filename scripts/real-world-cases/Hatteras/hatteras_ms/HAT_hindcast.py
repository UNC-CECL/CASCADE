#!/usr/bin/env python3
"""
Hatteras Island CASCADE hindcast: one model run.

HOW TO USE
    1. Choose the run in hat_run.yaml: the period, the source/sink preset, the scenario, groin on or off.
    2. Run this file:
           python scripts/real-world-cases/Hatteras/hatteras_ms/HAT_hindcast.py
    3. The results are in D:\\hatteras_ch1\\model_runs\\output\\raw_runs\\matrix\\<period>\\<preset>\\<run name>\\

WHAT IT DOES
    Step 1  Gather every input the model needs, from the settings    (hindcast/HAT_hindcast_setup.py)
    Step 2  Make the run's output folder
    Step 3  Build the model
    Step 4  Run it, one model year at a time
    Step 5  Check the run and save its results                       (hindcast/HAT_hindcast_outputs.py)
    Step 6  Draw the figures                                          (figure_making/model_output/run_figures.py)

Everything a step works out is kept on `r`, one value each: r.START_YEAR, r.STORM_FILE, r.RUN_DIR, ...
HAT_hindcast.ipynb runs the same steps one at a time, with QC plots in between.
Details: scripts/real-world-cases/Hatteras/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-23
"""

# =============================================================================
# Setup: find the Hatteras code and read the settings (nothing to edit here)
# =============================================================================

import os
import shutil
import sys
import time
from pathlib import Path

# The Hatteras code is not installed, so its folders go on the import path
SCRIPTS_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
for _folder in (SCRIPTS_DIR, SCRIPTS_DIR / "hatteras_ms" / "hindcast",
                SCRIPTS_DIR / "figure_making" / "model_output"):
    if str(_folder) not in sys.path:
        sys.path.insert(0, str(_folder))

# The settings in hat_run.yaml (an HAT_* environment variable overrides any of them)
from HAT_hindcast_config import load_run_config  # noqa: E402
RUN_CONFIG = load_run_config()

# Two settings have to be applied before CASCADE and matplotlib are imported
os.environ["CASCADE_USE_SANDBOX"] = "1" if RUN_CONFIG.use_sandbox_cascade else "0"
SHOW_FIGURES = bool(RUN_CONFIG.show_figures)     # left at null in the yaml: figures are saved, not shown
import matplotlib  # noqa: E402
if not SHOW_FIGURES:
    matplotlib.use("Agg")

# Fixes to CASCADE that the lab has not merged yet, applied before any model is built
import cascade_patch  # noqa: E402
CASCADE_PATCHES = cascade_patch.apply()

from tqdm.auto import tqdm  # noqa: E402
from cascade_pipeline.hindcast import build_cascade, run_cascade_simulation  # noqa: E402
from cascade_pipeline.run_registry import guard_run_dir, run_dir_contents  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_DOMAINS, HATTERAS_RELOCATION_CHECK_2004, HATTERAS_ROAD_EVENTS)
from HAT_hindcast_setup import prepare_run  # noqa: E402
from HAT_hindcast_outputs import pre_run_diagnostics, verify_run, write_run_outputs  # noqa: E402
from run_figures import make_run_figures  # noqa: E402


# =============================================================================
# Step 1. Gather every input the model needs
# =============================================================================
# From the settings: the period, topography, island offsets, storms, sea-level rise, source/sink rates,
# the road, nourishment, the groin, the run's name and the CoastSat target. Each part prints a report.

r = prepare_run(RUN_CONFIG)
r.CASCADE_PATCHES = CASCADE_PATCHES     # recorded in the run's metadata


# =============================================================================
# Step 2. Make the run's output folder
# =============================================================================
# Named after the run. If it already holds a result the run stops, unless hat_run.yaml says overwrite: true
# (which empties the folder first).

r.RUN_NAME = r.RUN_NAME_BASE
r.RUN_DIR = str(r.OUTPUT_BASE_DIR / r.RUN_NAME)
_old_files = run_dir_contents(r.RUN_DIR)        # counted before the guard can empty the folder
guard_run_dir(r.RUN_DIR, overwrite=RUN_CONFIG.overwrite)
print(f"\nRUN_DIR               {r.RUN_DIR}"
      + ("   (OVERWRITE=True)" if RUN_CONFIG.overwrite else ""))
if _old_files:
    print(f"  replaced            {len(_old_files)} file(s) from the previous run")


# =============================================================================
# Step 3. Build the model (year 0)
# =============================================================================

# CASCADE edits its parameter file while it builds, so each run gets its own copy
r.RUN_PARAMETER_FILE = Path(r.RUN_DIR) / f"{r.RUN_NAME}-parameters.yaml"
shutil.copyfile(r.HATTERAS_DATA_BASE / r.PARAMETER_FILE, r.RUN_PARAMETER_FILE)

r.cascade = build_cascade(
    # The run
    name=r.RUN_NAME,
    run_years=r.RUN_YEARS,
    data_base=r.HATTERAS_DATA_BASE,
    parameter_file=str(r.RUN_PARAMETER_FILE),
    num_cores=r.NUM_CORES,
    alongshore_section_count=HATTERAS_DOMAINS.total_domains,

    # The island: topography, dunes and where each domain's shoreline starts
    elevation_file=r.ELEVATION_FILE_PATHS,
    dune_file=r.DUNE_FILE_PATHS,
    enable_shoreline_offset=True,
    shoreline_offset=r.island_offset,
    rmin=r.RMIN, rmax=r.RMAX,                     # dune growth rates
    berm_elevation=r.BERM_ELEVATION,
    MHW=r.MHW_ELEVATION,
    beta=r.BEACH_SLOPE,

    # The forcing: storms, waves, sea level and the source/sink rates
    storm_file=str(r.STORM_FILE),
    wave_height=r.Hs,
    wave_period=r.FIXED_WAVE_PERIOD,
    wave_asymmetry=r.FIXED_WAVE_ASYMMETRY,
    wave_angle_high_fraction=r.FIXED_WAVE_ANGLE_HIGH_FRACTION,
    sea_level_rise_rate=r.SEA_LEVEL_RISE_RATE,
    sea_level_constant=r.SEA_LEVEL_CONSTANT,
    background_erosion=r.BACKGROUND_EROSION_RATES,

    # NC-12: where the road is, how high, and which domains are managed
    roadway_management_on=r.ROADWAY_MANAGEMENT_ON,
    road_setback=r.road_setbacks_full,
    road_ele=r.ROAD_ELEVATION,
    road_width=r.ROAD_WIDTH,
    dune_design_elevation=r.DUNE_DESIGN_ELEVATION,
    dune_minimum_elevation=r.DUNE_MINIMUM_ELEVATION,
    relocation_setback_m=r.RELOCATION_SETBACK_M,
    road_dune_rebuild_trigger_m=r.ROAD_DUNE_REBUILD_TRIGGER_M,

    # The villages: beach/dune management, overwash removal and nourishment
    beach_dune_manager_on=r.BEACH_DUNE_MANAGEMENT_ON,
    overwash_filter=r.OVERWASH_FILTER,
    overwash_to_dune=r.OVERWASH_TO_DUNE,
    nourishment_volume=r.NOURISHMENT_VOLUME_INIT,  # rewritten every year from the fill schedule
    sandbag_management_on=r.SANDBAG_MANAGEMENT_ON,
    sandbag_elevation=r.SANDBAG_ELEVATION,

    # The Buxton groin
    groin_module=r.GROIN_ENABLED,
    groin_updrift_pad=r.GROIN_UPDRIFT_PAD,
    groin_downdrift_pad=r.GROIN_DOWNDRIFT_PAD,
    groin_strength_history=r.GROIN_STRENGTH_HISTORY,
    groin_start_year=r.FIRST_MODEL_YEAR,
    groin_spread_weights=r.GROIN_SPREAD_WEIGHTS,
)

pre_run_diagnostics(r)      # a few read-outs of the built model, printed before it runs


# =============================================================================
# Step 4. Run it, one model year at a time
# =============================================================================
# Each year also applies what really happened that year: the road relocations and the beach fills.

_start = time.perf_counter()
with tqdm(total=r.RUN_YEARS, desc=r.RUN_NAME, unit="yr") as _progress:
    r.cascade = run_cascade_simulation(
        cascade=r.cascade,
        name=r.RUN_NAME,
        run_dir=r.RUN_DIR,
        run_years=r.RUN_YEARS,
        start_year=r.FIRST_MODEL_YEAR,
        geometry=HATTERAS_DOMAINS,
        alongshore_section_count=HATTERAS_DOMAINS.total_domains,
        nourishment_schedule=r.BN_SCHEDULE_APPLIED,
        historical_road_events=HATTERAS_ROAD_EVENTS,
        relocations_enabled=r.ENABLE_HISTORICAL_ROAD_RELOCATIONS,
        setback_check=HATTERAS_RELOCATION_CHECK_2004,
        save_model_state=r.SAVE_MODEL_STATE,
        progress=_progress,
    )
r.RUN_SECONDS = time.perf_counter() - _start
print(f"\nruntime               {r.RUN_SECONDS / 60:.1f} min "
      f"({r.RUN_SECONDS / r.RUN_YEARS:.1f} s per model year)")


# =============================================================================
# Step 5. Check the run and save its results
# =============================================================================

verify_run(r)           # the right number of years, nourishment applied, which roads survived, groin extent
write_run_outputs(r)    # shoreline positions, change rates, skill against CoastSat, and the run's metadata


# =============================================================================
# Step 6. Draw the figures
# =============================================================================
# From what step 5 saved, so a finished run's figures can be redrawn later with run_figures.py

make_run_figures(r.RUN_DIR, gifs=RUN_CONFIG.make_gifs, show=SHOW_FIGURES, cs_series=r.cs_series)
print(f"\ndone                  {r.RUN_DIR}")
