#!/usr/bin/env python3
"""
HATTERAS ISLAND

Modeled shoreline change rate:
  total_change = shoreline_m[-1, :] - shoreline_m[0, :]
  change_rate  = total_change / (END_YEAR - START_YEAR)   [m/yr]

CoastSat overlays:
  - Reads mean_lrr annual rates (m/yr) from domain_lrr_summary.csv files
  - Maps GIS domain IDs (1–90) to CASCADE padded indices (15–104)
  - Two periods: 1984–2004 and 2004–2024

Overwash filter:
  - Applied per-domain; only set for developed communities, set to 0.4
  - Buxton:              GIS  7– 8 → padded 21–22
  - Avon:                GIS 21–31 → padded 35–45
  - Salvo/Waves/Rodanthe GIS 68–83 → padded 82–97
  - All other domains (undeveloped + buffers) remain at 0
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cascade.cascade import Cascade

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# =============================================================================

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 1    # GIS domain IDs: 1 .. 90
LAST_FILE_NUMBER  = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1  # = 90

TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS          # = 15  (first real domain in padded array)
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS  # = 105

# Roads occupy GIS domains 9–90  →  padded indices 23–104
FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN  = 90
START_ROAD_INDEX  = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS   # = 23
END_ROAD_INDEX    = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS + 1  # = 105

print("=" * 80)
print("HATTERAS ISLAND CASCADE DOMAIN CONFIGURATION")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS} (GIS IDs {FIRST_FILE_NUMBER}..{LAST_FILE_NUMBER})")
print(f"Buffers:        {NUM_BUFFER_DOMAINS} each side")
print(f"TOTAL_DOMAINS:  {TOTAL_DOMAINS}")
print(f"Real index span (padded): [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print(f"Road index span (padded): [{START_ROAD_INDEX}..{END_ROAD_INDEX - 1}]")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 2: DISPLAY OPTION
# =============================================================================

# Set to True  → x-axis shows only the 90 real domains (no buffers)
# Set to False → x-axis shows all 120 padded domains (buffers shaded)
PLOT_REAL_DOMAINS_ONLY = True

# =============================================================================
# SECTION 3: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR   = r"C:\Users\hanna\PycharmProjects\CASCADE"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")
COASTSAT_BASE_DIR  = os.path.join(
    PROJECT_BASE_DIR, "scripts", "input_preperation", "CoastSat"
)

# --- Dune offset files (one per calibration period) ---
# UPDATE these paths when new offset files are generated for 1984/2004 init years
DUNE_OFFSET_FILE_1984 = os.path.join(
    HATTERAS_DATA_BASE,
    "island_offset",
    "hindcast_1984",
    f"Island_Dune_Offsets_1984_PADDED_{TOTAL_DOMAINS}.csv",
)
DUNE_OFFSET_FILE_2004 = os.path.join(
    HATTERAS_DATA_BASE,
    "island_offset",
    "hindcast_1984",
    f"Island_Dune_Offsets_1984_ADDED_{TOTAL_DOMAINS}.csv",
)

# --- Storm files ---
# UPDATE paths once storm files are generated for the new periods
STORM_FILE_1984_2004 = os.path.join(
    HATTERAS_DATA_BASE, "storms", "hindcast_storms", "storms_1984_2004.npy"
)
STORM_FILE_2004_2024 = os.path.join( #NEED TO UPDATE TO 2004-2024
    HATTERAS_DATA_BASE, "storms", "hindcast_storms", "storms_1984_2004.npy"
)

ROAD_SETBACK_FILE = os.path.join( # Updated
    HATTERAS_DATA_BASE, "roads", "offset", "1984", "RoadSetback_1984.csv"
)

# --- Topo/dune init year ---
# 2009 data used for both periods (adjust if you have period-specific init)
TOPO_DUNE_INIT_YEAR = "2009"

# Filename only — CASCADE resolves it relative to datadir (hatteras_init)
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# =============================================================================
# SECTION 4: COASTSAT DATASETS (multi-source)
# =============================================================================
# Column names match domain_lrr_summary.csv output from coastsat_domain_lrr.py:
#   domain_number  – GIS domain IDs (1–90, stored as float)
#   mean_lrr       – domain-mean LRR (m/yr)

COASTSAT_DATASETS = [
    dict(
        label="CoastSat LRR (1984–2004)",
        csv=os.path.join(COASTSAT_BASE_DIR, "1984_2004", "domain_lrr_summary.csv"),
        domain_col="domain_number",
        rate_col="mean_lrr",
        mode="gis",        # domain_number is 1-based GIS ID → maps directly
        flip_alongshore=False,
        flip_sign=False,
    ),
    dict(
        label="CoastSat LRR (2004–2024)",
        csv=os.path.join(
            COASTSAT_BASE_DIR, "2004_2024_specific_dates", "domain_lrr_summary.csv"
        ),
        domain_col="domain_number",
        rate_col="mean_lrr",
        mode="gis",
        flip_alongshore=False,
        flip_sign=False,
    ),
]

# =============================================================================
# SECTION 5: SIMULATION PARAMETERS
# =============================================================================

# -----------------------------------------------------------------------------
# RUN NAME SUFFIX  ← edit this to label each experiment
# The full run name will be:  HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}_Hs{X}
# -----------------------------------------------------------------------------
RUN_NAME_SUFFIX = "SQ_BE"

START_YEAR = 1984
END_YEAR   = 2004
RUN_YEARS  = END_YEAR - START_YEAR  # passed to CASCADE as time_step_count

TO_METERS           = True
SEA_LEVEL_RISE_RATE = 0.003
NUM_CORES           = 4

ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_NOURISHMENT        = False
ENABLE_SANDBAG_PLACEMENT  = False

DUNE_REBUILD_HEIGHT    = 3.0
REBUILD_ELEV_THRESHOLD = 0.01   # dam

FIXED_WAVE_PERIOD              = 8
FIXED_WAVE_ASYMMETRY           = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.2
WAVE_HEIGHTS_TO_TEST           = [2.0]

DOMAIN_TICK_STEP = 5

# -----------------------------------------------------------------------------
# BACKGROUND EROSION RATES — per domain (dam/yr, negative = erosion)
#
# TOGGLE: set USE_BACKGROUND_EROSION = False to zero all rates in one step.
#         Set to True to apply the per-region values defined below.
#
# Padded index to GIS domain ID mapping:
#   Buffers               pad  0-14   (no GIS ID)
#   S. undeveloped        pad 15-20   GIS  1- 6
#   Buxton                pad 21-22   GIS  7- 8
#   Mid undeveloped       pad 23-34   GIS  9-20
#   Avon                  pad 35-45   GIS 21-31
#   Central undeveloped   pad 46-81   GIS 32-67
#   Salvo/Waves/Rodanthe  pad 82-97   GIS 68-83
#   N. undeveloped        pad 98-104  GIS 84-90
#   Buffers               pad 105-119 (no GIS ID)
# -----------------------------------------------------------------------------
USE_BACKGROUND_EROSION = False  # ← flip to False to zero all background rates

# Regional rate values (dam/yr) — only active when USE_BACKGROUND_EROSION = True
_S_RATE     =  0.0   # GIS  1- 6  S. undeveloped
_GIS1_RATE  = -45.0  # GIS  1 override (special; set 0 here to disable)
_BUX_RATE   =  0.0   # GIS  7- 8  Buxton
_M_RATE     =  0.0   # GIS  9-20  Mid undeveloped
_AVN_RATE   =  3.0   # GIS 21-31  Avon
_C_RATE     =  1.0   # GIS 32-67  Central undeveloped
_SWR_RATE   =  0.0   # GIS 68-83  Salvo / Waves / Rodanthe
_N_RATE     =  0.0   # GIS 84-90  N. undeveloped

# Derive working shorthand — DO NOT edit these lines, edit the _*_RATE values above
_scale   = 1.0 if USE_BACKGROUND_EROSION else 0.0
_B       = 0.0                    # Buffer — always 0 regardless of toggle
_GIS1    = _GIS1_RATE  * _scale
_S       = _S_RATE     * _scale
_BUX     = _BUX_RATE   * _scale
_M       = _M_RATE     * _scale
_AVN     = _AVN_RATE   * _scale
_C       = _C_RATE     * _scale
_SWR     = _SWR_RATE   * _scale
_N       = _N_RATE     * _scale

#                         pad:  0      1      2      3      4      5      6      7      8      9
BACKGROUND_EROSION_RATES = [
    _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   # pad   0- 9  | buf (left)
    _B,   _B,   _B,   _B,   _B,   _GIS1,_S,   _S,   _S,   _S,   # pad  10-19  | buf 10-14 | GIS  1- 5
    _S,   _BUX, _BUX, _M,   _M,   _M,   _M,   _M,   _M,   _M,   # pad  20-29  | GIS  6 | 7-8 (Buxton) | 9-15
    _M,   _M,   _M,   _M,   _M,   _AVN, _AVN, _AVN, _AVN, _AVN, # pad  30-39  | GIS 16-20 (Mid) | 21-25 (Avon)
    _AVN, _AVN, _AVN, _AVN, _AVN, _AVN, _C,   _C,   _C,   _C,   # pad  40-49  | GIS 26-31 (Avon) | 32-35 (Central)
    _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   # pad  50-59  | GIS 36-45 (Central)
    _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   # pad  60-69  | GIS 46-55 (Central)
    _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   # pad  70-79  | GIS 56-65 (Central)
    _C,   _C,   _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, # pad  80-89  | GIS 66-67 (Central) | 68-75 (SWR)
    _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _N,   _N,   # pad  90-99  | GIS 76-83 (SWR) | 84-85 (N. undev.)
    _N,   _N,   _N,   _N,   _N,   _B,   _B,   _B,   _B,   _B,   # pad 100-109 | GIS 86-90 (N. undev.) | buf (right)
    _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   # pad 110-119 | buf (right)
]

assert len(BACKGROUND_EROSION_RATES) == TOTAL_DOMAINS, (
    f"BACKGROUND_EROSION_RATES has {len(BACKGROUND_EROSION_RATES)} entries "
    f"but TOTAL_DOMAINS={TOTAL_DOMAINS}. Fix the list before running."
)

FLIP_SIGN_MODEL = True  # flips only the modeled sign (no alongshore reversal)

# -----------------------------------------------------------------------------
# OVERWASH FILTER — per community
# 0 = no filtering (undeveloped / buffer domains)
# 0.4 is based on Laura's article
# -----------------------------------------------------------------------------
OVERWASH_FILTER_DEFAULT         = 0.0   # undeveloped domains + all buffer domains
OVERWASH_FILTER_BUXTON          = 0.4   # GIS  7– 8  (south end of Buxton)
OVERWASH_FILTER_AVON            = 0.4   # GIS 21–31
OVERWASH_FILTER_SALVO_WAVES_ROD = 0.4   # GIS 68–83  (Salvo / Waves / Rodanthe)

# =============================================================================
# SECTION 6: START YEAR CONFIG
# =============================================================================

if START_YEAR == 1984:
    YEAR_COLUMN_INDEX = 0    # first column in dune offset file = 1984 init
    RUN_NAME_BASE     = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}"
    STORM_FILE        = STORM_FILE_1984_2004
    DUNE_OFFSET_FILE  = DUNE_OFFSET_FILE_1984
elif START_YEAR == 2004:
    YEAR_COLUMN_INDEX = 1    # second column in dune offset file = 2004 init
    RUN_NAME_BASE     = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}"
    STORM_FILE        = STORM_FILE_2004_2024
    DUNE_OFFSET_FILE  = DUNE_OFFSET_FILE_2004
else:
    print(f"❌ ERROR: Invalid START_YEAR {START_YEAR}. Must be 1984 or 2004.")
    sys.exit(1)

print("Simulation Configuration:")
print(f"  Start Year: {START_YEAR} | End Year: {END_YEAR} | RUN_YEARS: {RUN_YEARS}")
print(f"  WAVE_HEIGHTS: {WAVE_HEIGHTS_TO_TEST}")
print(f"  USE_BACKGROUND_EROSION: {USE_BACKGROUND_EROSION}  "
      f"(non-zero domains: {sum(1 for v in BACKGROUND_EROSION_RATES if v != 0)})")
print(f"  FLIP_SIGN_MODEL={FLIP_SIGN_MODEL}")
print(f"  PLOT_REAL_DOMAINS_ONLY={PLOT_REAL_DOMAINS_ONLY}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 7: LOAD INPUT DATA
# =============================================================================

print("Loading input data...")

try:
    dune_offset_raw = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    # Single-column CSV returns 1D array; ensure 2D for consistent indexing
    if dune_offset_raw.ndim == 1:
        dune_offset_raw = dune_offset_raw[:, np.newaxis]
    dune_offset_dam = dune_offset_raw[:, YEAR_COLUMN_INDEX] / 10.0  # m → dam
    print(f"✓ Loaded dune offsets: {dune_offset_dam.size} domains (dam)")

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=",")
    print(f"✓ Loaded road setbacks: {road_setbacks_raw.size} values")

except FileNotFoundError as e:
    print(f"❌ CRITICAL ERROR: Missing data file: {e.filename}")
    sys.exit(1)
except Exception as e:
    print(f"❌ CRITICAL ERROR loading data: {e}")
    sys.exit(1)

print(f"✓ Background erosion array: {len(BACKGROUND_EROSION_RATES)} domains "
      f"| non-zero: {sum(1 for v in BACKGROUND_EROSION_RATES if v != 0)}")

# Road setbacks: padded array (zeros outside road span)
road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX : START_ROAD_INDEX + num_road_values] = \
    road_setbacks_raw[:num_road_values]
print(f"✓ Road setback array prepared ({TOTAL_DOMAINS} domains)")

ROADWAY_MANAGEMENT_ON    = [ENABLE_ROADWAY_MANAGEMENT] * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON    = [ENABLE_SANDBAG_PLACEMENT]  * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [ENABLE_NOURISHMENT]        * TOTAL_DOMAINS

# =============================================================================
# SECTION 8: ELEVATION + DUNE FILE LISTS
# =============================================================================

print("Generating elevation + dune profile file paths...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS      = []

for _ in range(START_REAL_INDEX):  # left buffers
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_dune.npy"))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy"))

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):  # real domains
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    DUNE_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, "dunes", TOPO_DUNE_INIT_YEAR,
                     f"domain_{file_num}_dune_{TOPO_DUNE_INIT_YEAR}.npy"))
    ELEVATION_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, "topography", TOPO_DUNE_INIT_YEAR,
                     f"domain_{file_num}_topography_{TOPO_DUNE_INIT_YEAR}.npy"))

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):  # right buffers
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_dune.npy"))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy"))

print(f"✓ Generated {len(ELEVATION_FILE_PATHS)} elevation file paths")
print(f"✓ Generated {len(DUNE_FILE_PATHS)} dune file paths")
print("=" * 80 + "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_x_s_TS(b3d):
    """Extract shoreline time series from a Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError(
        "No shoreline time series found (x_s_TS / _x_s_TS) on Barrier3D object."
    )


def build_shoreline_matrix(cascade, to_meters=True):
    """Build [time × domain] shoreline matrix from a CASCADE object."""
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt   = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, ndom), dtype=float)
    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])
    if to_meters:
        shoreline *= 10.0  # dam → m
    return shoreline


def load_and_map_coastsat_to_cascade(
    coastsat_csv, domain_col, rate_col,
    start_real_index, end_real_index,
    first_file_number, last_file_number,
    mode="gis", flip_alongshore=False,
):
    """
    Load a CoastSat domain_lrr_summary CSV and map GIS domain IDs
    to CASCADE padded indices.

    Parameters
    ----------
    coastsat_csv : str
        Path to domain_lrr_summary.csv (columns: domain_number, mean_lrr, ...)
    domain_col : str
        Column with GIS domain IDs (1–90), typically 'domain_number' (stored as float)
    rate_col : str
        Column with LRR rates in m/yr, typically 'mean_lrr'
    mode : str
        'gis'     → domain IDs are 1-based GIS IDs (default for CoastSat)
        'real_id' → domain IDs are the raw file numbers (same as first_file_number..last_file_number)
        'auto'    → infer from data range

    Returns
    -------
    obs_x    : int array of CASCADE padded domain indices
    obs_rate : float array of annual rates (m/yr)
    """
    if coastsat_csv is None or not os.path.exists(coastsat_csv):
        print(f"⚠️  CoastSat file not found: {coastsat_csv}")
        return None, None

    cs = pd.read_csv(coastsat_csv)
    if domain_col not in cs.columns or rate_col not in cs.columns:
        raise ValueError(
            f"CoastSat CSV missing required columns. "
            f"Need '{domain_col}' and '{rate_col}'. "
            f"Found: {list(cs.columns)}"
        )

    # domain_number is stored as float (1.0, 2.0, ...) — convert safely
    raw_dom  = pd.to_numeric(cs[domain_col], errors="coerce").to_numpy()
    raw_rate = pd.to_numeric(cs[rate_col],   errors="coerce").to_numpy(dtype=float)

    m = np.isfinite(raw_dom) & np.isfinite(raw_rate)
    raw_dom  = raw_dom[m].astype(int)
    raw_rate = raw_rate[m]

    if len(raw_dom) == 0:
        print("⚠️  CoastSat file loaded but no valid rows found.")
        return None, None

    # Auto-detect mode if requested
    if mode == "auto":
        if raw_dom.min() >= first_file_number and raw_dom.max() <= last_file_number:
            mode_use = "real_id"
        else:
            mode_use = "gis"
        print(f"  COASTSAT_DOMAIN_MODE='auto' → inferred '{mode_use}' "
              f"for {os.path.basename(coastsat_csv)}")
    else:
        mode_use = mode

    if mode_use == "real_id":
        obs_x = start_real_index + (raw_dom - first_file_number)
    elif mode_use == "gis":
        obs_x = start_real_index + (raw_dom - 1)
    else:
        raise ValueError("mode must be 'auto', 'gis', or 'real_id'")

    keep     = (obs_x >= start_real_index) & (obs_x < end_real_index)
    obs_x    = obs_x[keep].astype(int)
    obs_rate = raw_rate[keep]

    if len(obs_x) == 0:
        print("⚠️  CoastSat mapped but nothing fell within the real-domain range.")
        return None, None

    if flip_alongshore:
        obs_x = start_real_index + (end_real_index - 1 - obs_x)
        print(f"  *** CoastSat alongshore direction flipped for "
              f"{os.path.basename(coastsat_csv)} ***")

    order    = np.argsort(obs_x)
    obs_x    = obs_x[order]
    obs_rate = obs_rate[order]

    print(f"  ✓ Loaded {len(obs_rate)} CoastSat points  |  "
          f"rate range: {obs_rate.min():.2f}..{obs_rate.max():.2f} m/yr  |  "
          f"CASCADE idx: {obs_x.min()}..{obs_x.max()}")

    return obs_x, obs_rate


# =============================================================================
# CASCADE RUNNER
# =============================================================================

def run_cascade_simulation(
    nt, name, storm_file, alongshore_section_count, num_cores,
    rmin, rmax, elevation_file, dune_file,
    dune_design_elevation, dune_minimum_elevation,
    road_ele, road_width, road_setback,
    overwash_filter, overwash_to_dune,
    nourishment_volume, background_erosion,
    roadway_management_on, beach_dune_manager_on,
    sea_level_rise_rate, sea_level_constant,
    sandbag_management_on, sandbag_elevation,
    enable_shoreline_offset, shoreline_offset,
    wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction,
):
    datadir = HATTERAS_DATA_BASE

    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=PARAMETER_FILE,

        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=wave_asymmetry,
        wave_angle_high_fraction=wave_angle_high_fraction,

        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,

        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,

        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,

        roadway_management_module=roadway_management_on,
        beach_nourishment_module=beach_dune_manager_on,
        sandbag_management_on=sandbag_management_on,
        alongshore_transport_module=True,
        community_economics_module=False,

        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        sandbag_elevation=sandbag_elevation,

        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,

        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,

        nourishment_volume=nourishment_volume,
        nourishment_interval=None,
    )

    for time_step in range(nt - 1):
        print(f"\rYear {time_step + 1}/{nt}", end="", flush=True)
        cascade.update()
        if getattr(cascade, "b3d_break", False):
            print(f"\n⚠️  Model stopped at year {time_step + 1} (b3d_break)")
            break

    save_path = os.path.join(OUTPUT_BASE_DIR, name)
    os.makedirs(save_path, exist_ok=True)
    cascade.save(save_path)
    print(f"\n✓ Saved: {save_path}")

    return cascade


# =============================================================================
# MAIN
# =============================================================================

def main():
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS
    DUNE_DESIGN_ELEVATION  = [DUNE_REBUILD_HEIGHT]    * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    ROAD_ELEVATION = 1.45
    ROAD_WIDTH     = 20.0

    # ── Per-domain overwash filter ───────────────────────────────────────────
    def _gis_to_pad(gis_id):
        return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)

    OVERWASH_FILTER = np.full(TOTAL_DOMAINS, OVERWASH_FILTER_DEFAULT)

    # Buxton: GIS 7–8 → padded 21–22
    OVERWASH_FILTER[_gis_to_pad(7)  : _gis_to_pad(8)  + 1] = OVERWASH_FILTER_BUXTON
    # Avon: GIS 21–31 → padded 35–45
    OVERWASH_FILTER[_gis_to_pad(21) : _gis_to_pad(31) + 1] = OVERWASH_FILTER_AVON
    # Salvo/Waves/Rodanthe: GIS 68–83 → padded 82–97
    OVERWASH_FILTER[_gis_to_pad(68) : _gis_to_pad(83) + 1] = OVERWASH_FILTER_SALVO_WAVES_ROD

    OVERWASH_FILTER = list(OVERWASH_FILTER)

    print("Overwash filter configuration:")
    print(f"  Domains with active filter (> 0): "
          f"{sum(1 for v in OVERWASH_FILTER if v > 0)} of {TOTAL_DOMAINS}")
    for label, gis_start, gis_end, val in [
        ("Buxton",              7,  8,  OVERWASH_FILTER_BUXTON),
        ("Avon",               21, 31,  OVERWASH_FILTER_AVON),
        ("Salvo/Waves/Rod.",   68, 83,  OVERWASH_FILTER_SALVO_WAVES_ROD),
    ]:
        p0, p1 = _gis_to_pad(gis_start), _gis_to_pad(gis_end)
        print(f"  {label:<22s} GIS {gis_start:2d}–{gis_end:2d} "
              f"→ padded {p0:3d}–{p1:3d}  filter={val}")
    print()

    OVERWASH_TO_DUNE   = 9
    NOURISHMENT_VOLUME = 0
    SANDBAG_ELEV       = 0

    SEA_LEVEL_CONSTANT = True

    time_span_years = END_YEAR - START_YEAR
    if time_span_years == 0:
        time_span_years = None

    # ── Load CoastSat (multi-source), mapped to PADDED CASCADE indices ───────
    cs_series = []
    for ds in COASTSAT_DATASETS:
        print(f"Loading CoastSat: {ds['label']}")
        x, r = load_and_map_coastsat_to_cascade(
            coastsat_csv=ds["csv"],
            domain_col=ds["domain_col"],
            rate_col=ds["rate_col"],
            start_real_index=START_REAL_INDEX,
            end_real_index=END_REAL_INDEX,
            first_file_number=FIRST_FILE_NUMBER,
            last_file_number=LAST_FILE_NUMBER,
            mode=ds.get("mode", "gis"),
            flip_alongshore=ds.get("flip_alongshore", False),
        )
        if x is None or r is None or len(r) == 0:
            continue
        if ds.get("flip_sign", False):
            r = -1.0 * r
        cs_series.append(
            dict(label=ds.get("label", os.path.basename(ds["csv"])), x=x, rate=r)
        )

    # ── Run CASCADE + compute modeled change rates ───────────────────────────
    rate_profiles = {}

    for Hs in WAVE_HEIGHTS_TO_TEST:
        run_name_hs = f"{RUN_NAME_BASE}_Hs{Hs:.1f}".replace(".", "p")

        cascade = run_cascade_simulation(
            nt=RUN_YEARS,
            name=run_name_hs,
            storm_file=STORM_FILE,
            alongshore_section_count=TOTAL_DOMAINS,
            num_cores=NUM_CORES,

            rmin=RMIN, rmax=RMAX,

            elevation_file=ELEVATION_FILE_PATHS,
            dune_file=DUNE_FILE_PATHS,

            dune_design_elevation=DUNE_DESIGN_ELEVATION,
            dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,

            road_ele=ROAD_ELEVATION,
            road_width=ROAD_WIDTH,
            road_setback=road_setbacks_full,

            overwash_filter=OVERWASH_FILTER,
            overwash_to_dune=OVERWASH_TO_DUNE,

            nourishment_volume=NOURISHMENT_VOLUME,
            background_erosion=BACKGROUND_EROSION_RATES,

            roadway_management_on=ROADWAY_MANAGEMENT_ON,
            beach_dune_manager_on=NOURISHMENT_MANAGEMENT_ON,

            sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
            sea_level_constant=SEA_LEVEL_CONSTANT,

            sandbag_management_on=SANDBAG_MANAGEMENT_ON,
            sandbag_elevation=SANDBAG_ELEV,

            enable_shoreline_offset=True,
            shoreline_offset=dune_offset_dam,

            wave_height=Hs,
            wave_period=FIXED_WAVE_PERIOD,
            wave_asymmetry=FIXED_WAVE_ASYMMETRY,
            wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,
        )

        shoreline_m = build_shoreline_matrix(cascade, to_meters=TO_METERS)
        nt_actual, ndom = shoreline_m.shape

        denom = time_span_years if time_span_years is not None else max(nt_actual - 1, 1)
        total_change = shoreline_m[-1, :] - shoreline_m[0, :]
        change_rate  = total_change / float(denom)

        if FLIP_SIGN_MODEL:
            change_rate *= -1.0

        rate_profiles[Hs] = change_rate  # length = TOTAL_DOMAINS (padded)

        # ── Per-run output folder ────────────────────────────────────────────
        run_dir = os.path.join(OUTPUT_BASE_DIR, run_name_hs)
        os.makedirs(run_dir, exist_ok=True)

        # ── Save CSV of modeled rates ────────────────────────────────────────
        csv_path = os.path.join(run_dir, f"{run_name_hs}_shoreline_change_rate.csv")
        pd.DataFrame({
            "cascade_padded_index": np.arange(TOTAL_DOMAINS),
            "gis_domain_id": (
                [None] * START_REAL_INDEX
                + list(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))
                + [None] * (TOTAL_DOMAINS - END_REAL_INDEX)
            ),
            "model_rate_m_per_yr": change_rate,
        }).to_csv(csv_path, index=False)
        print(f"✓ Saved rate CSV: {csv_path}")

        # ====================================================================
        # PLOT
        # ====================================================================

        if PLOT_REAL_DOMAINS_ONLY:
            # ── REAL DOMAINS ONLY ────────────────────────────────────────────
            gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)  # 1..90

            fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)

            real_rate = change_rate[START_REAL_INDEX:END_REAL_INDEX]
            ax.plot(gis_ids, real_rate, linewidth=2, label=f"Model Hs={Hs} m")

            for s in cs_series:
                gis_x = s["x"] - START_REAL_INDEX + FIRST_FILE_NUMBER
                ax.plot(
                    gis_x, s["rate"],
                    linestyle="--", marker="x", linewidth=1.5, markersize=6,
                    label=s["label"],
                )

            community_spans = [
                (7,  8,  "Buxton"),
                (21, 31, "Avon"),
                (68, 83, "Salvo/Waves/Rod."),
            ]
            for gis_lo, gis_hi, comm_label in community_spans:
                ax.axvspan(gis_lo - 0.5, gis_hi + 0.5,
                           alpha=0.10, color="steelblue", zorder=0)
                ax.text((gis_lo + gis_hi) / 2, ax.get_ylim()[0],
                        comm_label, ha="center", va="bottom",
                        fontsize=7, style="italic", color="steelblue")

            ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

            be_label = "on" if USE_BACKGROUND_EROSION else "off"
            xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

            ax.set_xlabel("GIS Domain ID (1–90)")
            ax.set_ylabel("Shoreline change rate (m/yr)")
            ax.set_title(
                f"Modeled vs CoastSat Shoreline Change Rate – Hatteras Island | "
                f"Real domains only | SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                f"{START_YEAR}–{END_YEAR} | Hs={Hs} m | BE={be_label}"
            )
            ax.grid(alpha=0.3)
            ax.legend()

            fig_suffix = "REAL_DOMAINS_ONLY"

        else:
            # ── ALL DOMAINS (BUFFERS INCLUDED) ───────────────────────────────
            domain_numbers = np.arange(TOTAL_DOMAINS)

            fig, ax = plt.subplots(figsize=(22, 6), constrained_layout=True)

            ax.axvspan(0, START_REAL_INDEX - 0.5, alpha=0.12, color="red", label="Buffer")
            ax.axvspan(END_REAL_INDEX - 0.5, TOTAL_DOMAINS - 1, alpha=0.12, color="red")

            community_spans_pad = [
                (_gis_to_pad(7),  _gis_to_pad(8),  "Buxton"),
                (_gis_to_pad(21), _gis_to_pad(31), "Avon"),
                (_gis_to_pad(68), _gis_to_pad(83), "Salvo/Waves/Rod."),
            ]
            for pad_lo, pad_hi, comm_label in community_spans_pad:
                ax.axvspan(pad_lo - 0.5, pad_hi + 0.5,
                           alpha=0.10, color="steelblue", zorder=0)

            ax.plot(domain_numbers, change_rate, linewidth=2, label=f"Model Hs={Hs} m")

            for s in cs_series:
                ax.plot(
                    s["x"], s["rate"],
                    linestyle="--", marker="x", linewidth=1.5, markersize=6,
                    label=s["label"],
                )

            ax.axvline(START_REAL_INDEX, linestyle="--", linewidth=1, color="k", alpha=0.5)
            ax.axvline(END_REAL_INDEX - 1, linestyle="--", linewidth=1, color="k", alpha=0.5)
            ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

            y_top = ax.get_ylim()[1] if ax.get_ylim()[1] != 1.0 else 5.0
            ax.text((0 + START_REAL_INDEX) / 2, 0, "Left\nbuffer",
                    ha="center", va="center", fontsize=8, style="italic", alpha=0.6)
            ax.text((START_REAL_INDEX + END_REAL_INDEX) / 2, 0, "Real island (GIS 1–90)",
                    ha="center", va="center", fontsize=9, fontweight="bold", alpha=0.5)
            ax.text((END_REAL_INDEX + TOTAL_DOMAINS) / 2, 0, "Right\nbuffer",
                    ha="center", va="center", fontsize=8, style="italic", alpha=0.6)

            xticks = np.arange(0, TOTAL_DOMAINS, DOMAIN_TICK_STEP)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)
            ax.set_xlabel("CASCADE domain index (including buffers, 0–119)")

            top_ax = ax.secondary_xaxis("top")
            top_positions = []
            top_labels    = []
            for gis_id in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP):
                j = START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)
                top_positions.append(j)
                top_labels.append(str(gis_id))
            top_ax.set_xticks(top_positions)
            top_ax.set_xticklabels(top_labels, fontsize=9)
            top_ax.set_xlabel(f"GIS Domain ID ({FIRST_FILE_NUMBER}–{LAST_FILE_NUMBER})")

            be_label = "on" if USE_BACKGROUND_EROSION else "off"
            ax.set_ylabel("Shoreline change rate (m/yr)")
            ax.set_title(
                f"Modeled vs CoastSat Shoreline Change Rate – Hatteras Island | "
                f"All domains (buffers included) | SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                f"{START_YEAR}–{END_YEAR} | Hs={Hs} m | BE={be_label}"
            )
            ax.grid(alpha=0.3)
            ax.legend()

            fig_suffix = "ALL_DOMAINS_WITH_BUFFERS"

        fig_out = os.path.join(run_dir, f"{run_name_hs}_shoreline_change_rate_{fig_suffix}.png")
        fig.savefig(fig_out, dpi=300, bbox_inches="tight")
        print(f"✓ Saved plot: {fig_out}")
        plt.show()

    if not rate_profiles:
        print("❌ No successful runs; nothing to plot.")
        sys.exit(1)


if __name__ == "__main__":
    main()
