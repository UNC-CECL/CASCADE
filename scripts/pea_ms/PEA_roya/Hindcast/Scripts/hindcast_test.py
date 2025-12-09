# ==============================================================================
# HINDCAST SIMULATION SCRIPT – HATTERAS ISLAND (BUFFER-INCLUSIVE)
# ==============================================================================
# Author: Hannah Henry (UNC-Chapel Hill)
# Last Modified: 2025-10-01
#
# Summary:
# - Runs a CASCADE hindcast with 15 left + 92 real + 15 right domains (TOTAL=122)
# - Mirrors "Ocracoke" logic:
#     * buffers included in files & counts
#     * yearly management checks with "all-or-nothing" gates
#     * shoreline_offset enabled
# - Roads/nourishments OFF by default (natural-only baseline)
# - Robust, project-relative paths & length checks for reproducibility
# ==============================================================================

import os
import logging
import copy
from pathlib import Path
import numpy as np
import pandas as pd
from cascade.cascade import Cascade  # CASCADE core class

# ------------------------------
# Project paths (set once)
# ------------------------------
# Change this to your local CASCADE repo root if different:
BASE = Path(r"/").resolve()
DATA = BASE / "data" / "hatteras_init"
OUT  = BASE / "output"

# Input subfolders
STORMS_DIR   = DATA / "storms" / "hindcast_storms"
DUNES_DIR    = DATA / "dunes"
TOPO_DIR     = DATA / "topography_dunes"
BUFFER_DIR   = DATA / "buffer"
ROADS_DIR    = DATA / "roads"
OFFSET_DIR   = DATA / "island_offset"
PARAMS_FILE  = "Hatteras-CASCADE-parameters.yaml"  # lives under DATA at runtime

# ------------------------------
# Logging
# ------------------------------
OUT.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=str(OUT / "cascade_hindcast.log"),
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
run_name = "Hindcast_1978_1997_buffers"
logging.info(f"Starting simulation: {run_name}")

# ------------------------------
# Geometry & run timing
# ------------------------------
LEFT_BUF = 15
RIGHT_BUF = 15
REAL = 105
TOTAL = LEFT_BUF + REAL + RIGHT_BUF

# Choose ONE absolute start year for shoreline offset selection.
# (This mirrors your labmate’s pattern of selecting a single column by start year.)
START_YEAR = 1978    # <- set to 1997 if you want that start instead
RUN_YEARS  = 23      # yearly timesteps; absolute end depends on model’s internal calendar; this runs 19 years

# ------------------------------
# Storm file
# ------------------------------
# Keep the long horizon file so indexing is always safe
STORM_FILE = STORMS_DIR / "HAT_1978_2022_Final_Hindcast_Storms.npy"

# ------------------------------
# Road & dune offset sources (optional for future management)
# ------------------------------
# Roads (natural-only in this script, but keeping references here)
ROAD_SETBACKS_CSV = ROADS_DIR / "offset" / "1978_Road_Setbacks_CASCADE_Ready.csv"
ROAD_ELEV_CSV     = ROADS_DIR / "elevation" / "1978_road_elevation_transect_summary.csv"

# Dune/shoreline offsets
# Expect columns named as years (e.g., "1978","1988","1997", ...)
DUNE_OFFSET_CSV = OFFSET_DIR / "1978" / "Relative_Dune_1978_Offset.csv"

# ------------------------------
# Helpers
# ------------------------------
def pad_with_buffers(arr, left_val=0, right_val=0):
    """Pad a REAL-length list/array with buffer values on both ends to TOTAL length."""
    arr = np.asarray(arr).tolist()
    assert len(arr) == REAL, f"Expected REAL length {REAL}, got {len(arr)}"
    return [left_val] * LEFT_BUF + arr + [right_val] * RIGHT_BUF

def assert_lengths(**named_lists):
    """Ensure lists/arrays equal TOTAL length."""
    for name, val in named_lists.items():
        if np.isscalar(val):
            # Scalars allowed; Cascade may accept scalars for some params
            continue
        n = len(val)
        assert n == TOTAL, f"{name} length {n} != TOTAL {TOTAL}"

# ------------------------------
# Shoreline/dune offsets (enabled)
# ------------------------------
dune_df = pd.read_csv(DUNE_OFFSET_CSV)
col = str(START_YEAR)
assert col in dune_df.columns, f"Column '{col}' not found in {DUNE_OFFSET_CSV.name}"
# Mirror labmate’s scaling convention (they multiply by 10)
# If your CASCADE branch expects meters instead, remove the *10.
dune_offset_real = dune_df[col].to_numpy() * 10.0  # REAL-length expected
assert len(dune_offset_real) == REAL, "Offset CSV must have 1 value per real domain"

shoreline_offset = pad_with_buffers(dune_offset_real, left_val=0, right_val=0)
shoreline_offset_enabled = False

# ------------------------------
# Background shoreline change (m/yr)
# ------------------------------
# Start with zeros for the 92 real domains; pad to TOTAL with zeros for buffers
background_real = [0.0] * REAL
background_threshold_list = pad_with_buffers(background_real, left_val=0.0, right_val=0.0)

# ------------------------------
# Roads & sandbags (natural-only baseline)
# ------------------------------
road_setbacks = [0] * TOTAL         # units must match your Cascade fork; labmate uses cm after *10
road_ele      = [0.0] * TOTAL       # m MHW
road_cells    = [False] * TOTAL
sandbag_cells = [False] * TOTAL

# ------------------------------
# Elevation & dune files (TOTAL-length with buffers)
# ------------------------------
e_file = []
d_file = []

# Left buffers
buf_dune = BUFFER_DIR / "Sample_1_dune.npy"
buf_topo = BUFFER_DIR / "Sample_1_topography.npy"
for _ in range(LEFT_BUF):
    d_file.append(str(buf_dune))
    e_file.append(str(buf_topo))

# Real domains (1..REAL)
for i in range(1, REAL + 1):
    dune_name = DUNES_DIR / f"domain_{i}_resampled_dune_2009.npy"
    elev_name = TOPO_DIR  / f"domain_{i}_resampled_topography_2009.npy"
    d_file.append(str(dune_name))
    e_file.append(str(elev_name))

# Right buffers
for _ in range(RIGHT_BUF):
    d_file.append(str(buf_dune))
    e_file.append(str(buf_topo))

# Sanity checks
assert_lengths(
    e_file=e_file,
    d_file=d_file,
    road_setbacks=road_setbacks,
    road_ele=road_ele,
    road_cells=road_cells,
    sandbag_cells=sandbag_cells,
    background_threshold_list=background_threshold_list,
    shoreline_offset=shoreline_offset,
)

# ==============================================================================
# Core alongshore-connected wrapper (mirrors labmate logic; no early 'pass')
# ==============================================================================
def alongshore_connected(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    beach_width_threshold,   # TOTAL-length list (m)
    rmin,                    # TOTAL-length list
    rmax,                    # TOTAL-length list
    elevation_file,          # TOTAL-length list of .npy
    dune_file,               # TOTAL-length list of .npy
    dune_design_elevation,   # TOTAL-length list (m)
    dune_minimum_elevation,  # scalar or TOTAL-length list (m MHW)
    road_ele,                # TOTAL-length list (m MHW)
    road_width,              # scalar (m)
    road_setback,            # TOTAL-length list (units per fork; here mirror labmate style)
    overwash_filter,         # scalar (%)
    overwash_to_dune,        # scalar (%)
    nourishment_volume,      # scalar (m^3/m)
    background_erosion,      # TOTAL-length list (m/yr)
    rebuild_dune_threshold,  # scalar (m MHW)
    roadway_management_on,   # TOTAL-length bool list
    beach_dune_manager_on,   # TOTAL-length bool list
    sea_level_rise_rate=0.0056,
    sea_level_constant=True,
    trigger_dune_knockdown=False,
    group_roadway_abandonment=None,
    sandbag_management_on=False,  # TOTAL-length bool list or False
    sandbag_elevation=4,
    enable_shoreline_offset=True,
    shoreline_offset=None,        # TOTAL-length list
):
    # Data directory used internally by Cascade (relative to repo root)
    datadir = str(DATA) + "/"

    cascade = Cascade(
        datadir=datadir,
        name=name,
        storm_file=str(storm_file),
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=PARAMS_FILE,
        wave_height=1.0,
        wave_period=7.0,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,  # TOTAL
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,
        roadway_management_module=roadway_management_on,
        alongshore_transport_module=True,
        beach_nourishment_module=beach_dune_manager_on,
        community_economics_module=False,
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        trigger_dune_knockdown=trigger_dune_knockdown,
        group_roadway_abandonment=group_roadway_abandonment,
        nourishment_interval=None,
        nourishment_volume=nourishment_volume,
        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,
        sandbag_management_on=sandbag_management_on,
        sandbag_elevation=sandbag_elevation,
        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
    )

    # Threshold relative to berm elevation (mirrors labmate's logic; note ×10 conversion)
    dune_rebuild_abs_thresh = rebuild_dune_threshold + (cascade.barrier3d[0].BermEl * 10)

    # ---- Yearly loop ----
    for time_step in range(nt - 1):
        print("\r", "Time Step:", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

        t = cascade.barrier3d[0].time_index
        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_nourish_now  = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):
            if cascade.community_break[iB3D]:
                continue

            if beach_dune_manager_on[iB3D]:
                # Nourishment trigger by beach width
                if cascade.nourishments[iB3D].beach_width[t - 1] < beach_width_threshold[iB3D]:
                    tmp_nourish_now[iB3D] = 1

                # Dune rebuild trigger by min crest elevation (domain-wide)
                dune_crests = cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
                min_dune_crest = (np.min(dune_crests) + cascade.barrier3d[iB3D].BermEl) * 10  # m MHW
                if min_dune_crest < dune_rebuild_abs_thresh:
                    tmp_rebuild_dune[iB3D] = 1

        # All-or-nothing gates (mirror labmate)
        if np.all(tmp_nourish_now[beach_dune_manager_on]) == 1:
            cascade.nourish_now = tmp_nourish_now
        if np.all(tmp_rebuild_dune[beach_dune_manager_on]) == 1:
            cascade.rebuild_dune_now = tmp_rebuild_dune

    # ---- Save outputs ----
    run_dir = OUT / run_name
    run_dir.mkdir(parents=True, exist_ok=True)
    cascade.save(str(run_dir))
    logging.info(f"Simulation completed and saved to {run_dir}")

    return cascade

# ==============================================================================
# Launcher (build TOTAL-length arrays and run)
# ==============================================================================
def alongshore_uniform():
    # Define REAL-length arrays first
    beach_width_threshold_real = [30.0] * REAL       # m
    rmin_real = [0.55] * REAL
    rmax_real = [0.95] * REAL
    dune_design_elev_real = [1.5] * REAL             # m target crest height

    # Pad to TOTAL (buffers get neutral values)
    beach_width_threshold = pad_with_buffers(beach_width_threshold_real)
    rmin = pad_with_buffers(rmin_real)
    rmax = pad_with_buffers(rmax_real)
    dune_design_elevation = pad_with_buffers(dune_design_elev_real)

    # Management toggles (natural-only baseline)
    roads_on_real = [False] * REAL
    nourishments_on_real = [False] * REAL
    roadway_management_on   = pad_with_buffers(roads_on_real, left_val=False, right_val=False)
    beach_dune_manager_on   = pad_with_buffers(nourishments_on_real, left_val=False, right_val=False)

    # Scalars
    dune_minimum_elevation = 1.0   # m MHW (scalar accepted by Cascade)
    road_width = 20                # m
    overwash_filter = 0            # %
    overwash_to_dune = 9           # %
    nourishment_volume = 100       # m^3/m
    rebuild_dune_threshold = 1.0   # m MHW
    sea_level_rise_rate = 0.0056   # m/yr (linear)
    sea_level_constant  = True
    num_cores = 4

    # Final safety checks
    assert_lengths(
        beach_width_threshold=beach_width_threshold,
        rmin=rmin,
        rmax=rmax,
        dune_design_elevation=dune_design_elevation,
        roadway_management_on=roadway_management_on,
        beach_dune_manager_on=beach_dune_manager_on,
    )

    alongshore_connected(
        nt=RUN_YEARS,
        name=run_name,
        storm_file=STORM_FILE,
        alongshore_section_count=TOTAL,
        num_cores=num_cores,
        beach_width_threshold=beach_width_threshold,
        rmin=rmin,
        rmax=rmax,
        elevation_file=e_file,
        dune_file=d_file,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setbacks,
        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,
        nourishment_volume=nourishment_volume,
        background_erosion=background_threshold_list,
        rebuild_dune_threshold=rebuild_dune_threshold,
        roadway_management_on=roadway_management_on,
        beach_dune_manager_on=beach_dune_manager_on,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_constant=sea_level_constant,
        sandbag_management_on=sandbag_cells,
        sandbag_elevation=4,
        shoreline_offset=shoreline_offset,
        enable_shoreline_offset=shoreline_offset_enabled,
    )

# ------------------------------
# Run
# ------------------------------
if __name__ == "__main__":
    alongshore_uniform()
