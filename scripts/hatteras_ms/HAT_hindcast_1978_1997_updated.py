#!/usr/bin/env python3
"""
HATTERAS ISLAND

Modeled shoreline change rate:
  total_change = shoreline_m[-1, :] - shoreline_m[0, :]
  change_rate  = total_change / (END_YEAR - START_YEAR)   [m/yr]

DSAS overlays:
  - Reads LRR annual rates (m/yr) from one or more CSVs
  - Maps GIS domain IDs (1–90) to CASCADE padded indices (15–104)
  - DSAS is filtered to real-domain span only

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

DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE,
    "island_offset",
    "hindcast_1978_1997",
    f"Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv",
)

STORM_FILE_1978_1997 = os.path.join(HATTERAS_DATA_BASE, "storms", "hindcast_storms", "storms_1978_1997.npy")
STORM_FILE_1997_2019 = os.path.join(HATTERAS_DATA_BASE, "storms", "hindcast_storms", "storms_1997_2019.npy")

ROAD_SETBACK_FILE = os.path.join(HATTERAS_DATA_BASE, "roads", "offset", "1978", "RoadSetback_1978.csv")

# Filename only — CASCADE resolves it relative to datadir (hatteras_init)
PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# =============================================================================
# SECTION 4: DSAS DATASETS (multi-source)
# =============================================================================

DSAS_DATASETS = [
    dict(
        label="DSAS (LRR, 1978–1997)",
        csv=os.path.join(
            HATTERAS_DATA_BASE, "shoreline_change",
            "dsas_1978_1997_domain_means_SIMPLE.csv"
        ),
        domain_col="domain_id",       # GIS domain IDs (1–90)
        rate_col="annual_rate_m_per_yr",
        mode="auto",                  # "auto" | "gis" | "real_id"
        flip_alongshore=False,
        flip_sign=False,
    ),
    # Add a second period's DSAS file here when ready, e.g.:
    # dict(
    #     label="DSAS (LRR, 1997–2019)",
    #     csv=os.path.join(HATTERAS_DATA_BASE, "shoreline_change",
    #                      "dsas_1997_2019_domain_means_SIMPLE.csv"),
    #     domain_col="domain_id",
    #     rate_col="annual_rate_m_per_yr",
    #     mode="auto",
    #     flip_alongshore=False,
    #     flip_sign=False,
    # ),
]

# =============================================================================
# SECTION 5: SIMULATION PARAMETERS
# =============================================================================

# -----------------------------------------------------------------------------
# RUN NAME SUFFIX  ← edit this to label each experiment
# The full run name will be:  HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}_Hs{X}
# Example suffixes: "SQ_withBE", "filter_test", "Hs_sweep_v2", "groin_test"
# All outputs (model files, CSV, plot) are saved under output/raw_runs/{run_name}/
# -----------------------------------------------------------------------------
RUN_NAME_SUFFIX = "SQ_BE"

START_YEAR = 1978
END_YEAR   = 1997
RUN_YEARS  = END_YEAR - START_YEAR  # passed to CASCADE as time_step_count

TO_METERS         = True
SEA_LEVEL_RISE_RATE = 0.003
NUM_CORES         = 4

ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_NOURISHMENT        = False
ENABLE_SANDBAG_PLACEMENT  = False

DUNE_REBUILD_HEIGHT   = 3.0
REBUILD_ELEV_THRESHOLD = 0.01   # dam

FIXED_WAVE_PERIOD             = 8
FIXED_WAVE_ASYMMETRY          = 0.7
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.2
WAVE_HEIGHTS_TO_TEST          = [2.0]

DOMAIN_TICK_STEP = 5

# -----------------------------------------------------------------------------
# BACKGROUND EROSION RATES — per domain (dam/yr, negative = erosion)
# Edit individual domain values directly in the list below.
# Named labels are provided as convenient shorthand — use them as-is or
# replace any entry with a literal value to override a specific domain.
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
_B   = 0.0   # Buffer  (always 0 - do not change)
_S   = 0.0   # GIS  1- 6  N. undeveloped
_BUX = 0.0   # GIS  7- 8  Buxton
_M   = 0.0   # GIS  9-20  Mid undeveloped
_AVN = 3.0   # GIS 21-31  Avon
_C   = 1.0   # GIS 32-67  Central undeveloped
_SWR = 0.0   # GIS 68-83  Salvo / Waves / Rodanthe
_N   = 0.0   # GIS 84-90  S. undeveloped

#                         pad:  0     1     2     3     4     5     6     7     8     9
BACKGROUND_EROSION_RATES = [
    _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   # buf   0- 9
    _B,   _B,   _B,   _B,   _B,   -45,   _S,   _S,   _S,   _S,   # buf  10-14 | GIS  1- 5
    _S,   _BUX, _BUX, _M,   _M,   _M,   _M,   _M,   _M,   _M,   # GIS  6 | 7-8 (Buxton) | 9-15
    _M,   _M,   _M,   _M,   _M,   _AVN, _AVN, _AVN, _AVN, _AVN, # GIS 16-20 | 21-25 (Avon)
    _AVN, _AVN, _AVN, _AVN, _AVN, _AVN, 3,  3,   _C,   _C,   # GIS 26-31 (Avon) | 32-35 (Central)
    _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   # GIS 36-45  (Central)
    _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   _C,   # GIS 46-55  (Central)
    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   # GIS 56-65  (Central)
    1,   1,   1, 1, 1, 1, 1, 1, 1, 1, # GIS 66-67 (Central) | 68-75 (SWR)
    _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _SWR, _N,   _N,   # GIS 76-83 (SWR) | 84-85 (S. undeveloped)
    _N,   _N,   _N,   _N,   _N,   _B,   _B,   _B,   _B,   _B,   # GIS 86-90 | buf 105-109
    _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   _B,   # buf 110-119
]

FLIP_SIGN_MODEL = True  # flips only the modeled sign (no alongshore reversal)

# -----------------------------------------------------------------------------
# OVERWASH FILTER — per community
# 0 = no filtering (undeveloped / buffer domains)
# 0.4 is based on Lauras article
# -----------------------------------------------------------------------------
OVERWASH_FILTER_DEFAULT         = 0.0   # undeveloped domains + all buffer domains
OVERWASH_FILTER_BUXTON          = 0.4   # GIS  7– 8  (south end of Buxton)
OVERWASH_FILTER_AVON            = 0.4   # GIS 21–31
OVERWASH_FILTER_SALVO_WAVES_ROD = 0.4   # GIS 68–83  (Salvo / Waves / Rodanthe)

# =============================================================================
# SECTION 6: START YEAR CONFIG
# =============================================================================

if START_YEAR == 1978:
    YEAR_COLUMN_INDEX = 0
    RUN_NAME_BASE = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}"
    STORM_FILE    = STORM_FILE_1978_1997
elif START_YEAR == 1997:
    YEAR_COLUMN_INDEX = 1
    RUN_NAME_BASE = f"HAT_{START_YEAR}_{END_YEAR}_{RUN_NAME_SUFFIX}"
    STORM_FILE    = STORM_FILE_1997_2019
else:
    print(f"❌ ERROR: Invalid START_YEAR {START_YEAR}. Must be 1978 or 1997.")
    sys.exit(1)

print("Simulation Configuration:")
print(f"  Start Year: {START_YEAR} | End Year: {END_YEAR} | RUN_YEARS: {RUN_YEARS}")
print(f"  WAVE_HEIGHTS: {WAVE_HEIGHTS_TO_TEST}")
print(f"  BACKGROUND_EROSION: see BACKGROUND_EROSION_RATES list in Section 5")
print(f"  FLIP_SIGN_MODEL={FLIP_SIGN_MODEL}")
print(f"  PLOT_REAL_DOMAINS_ONLY={PLOT_REAL_DOMAINS_ONLY}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 7: LOAD INPUT DATA
# =============================================================================

print("Loading input data...")

try:
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10.0  # m → dam
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

ROADWAY_MANAGEMENT_ON   = [ENABLE_ROADWAY_MANAGEMENT]  * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON   = [ENABLE_SANDBAG_PLACEMENT]   * TOTAL_DOMAINS
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
        os.path.join(HATTERAS_DATA_BASE, "dunes", "2009", f"domain_{file_num}_dune_2009.npy"))
    ELEVATION_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, "topography", "2009", f"domain_{file_num}_topography_2009.npy"))

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


def load_and_map_dsas_to_cascade(
    dsas_csv, domain_col, rate_col,
    start_real_index, end_real_index,
    first_file_number, last_file_number,
    mode="auto", flip_alongshore=False,
):
    """
    Load a DSAS CSV and map domain IDs to CASCADE padded indices.

    Returns
    -------
    obs_x    : int array of CASCADE padded domain indices
    obs_rate : float array of annual rates (m/yr)
    """
    if dsas_csv is None or not os.path.exists(dsas_csv):
        print(f"⚠️  DSAS file not found: {dsas_csv}")
        return None, None

    dsas = pd.read_csv(dsas_csv)
    if domain_col not in dsas.columns or rate_col not in dsas.columns:
        raise ValueError(
            f"DSAS CSV missing required columns. Need '{domain_col}' and '{rate_col}'. "
            f"Found: {list(dsas.columns)}"
        )

    raw_dom  = dsas[domain_col].to_numpy()
    raw_rate = dsas[rate_col].to_numpy(dtype=float)

    m = np.isfinite(raw_dom.astype(float)) & np.isfinite(raw_rate)
    raw_dom  = raw_dom[m].astype(int)
    raw_rate = raw_rate[m]

    if len(raw_dom) == 0:
        print("⚠️  DSAS file loaded but no valid rows found.")
        return None, None

    # Auto-detect whether domain IDs are real file numbers or 1-based GIS IDs
    if mode == "auto":
        if raw_dom.min() >= first_file_number and raw_dom.max() <= last_file_number:
            mode_use = "real_id"
        else:
            mode_use = "gis"
        print(f"  DSAS_DOMAIN_MODE='auto' → inferred '{mode_use}' "
              f"for {os.path.basename(dsas_csv)}")
    else:
        mode_use = mode

    if mode_use == "real_id":
        obs_x = start_real_index + (raw_dom - first_file_number)
    elif mode_use == "gis":
        obs_x = start_real_index + (raw_dom - 1)
    else:
        raise ValueError("mode must be 'auto', 'gis', or 'real_id'")

    keep    = (obs_x >= start_real_index) & (obs_x < end_real_index)
    obs_x   = obs_x[keep].astype(int)
    obs_rate = raw_rate[keep]

    if len(obs_x) == 0:
        print("⚠️  DSAS mapped but nothing fell within the real-domain range.")
        return None, None

    if flip_alongshore:
        obs_x = start_real_index + (end_real_index - 1 - obs_x)
        print(f"  *** DSAS alongshore direction flipped for "
              f"{os.path.basename(dsas_csv)} ***")

    order    = np.argsort(obs_x)
    obs_x    = obs_x[order]
    obs_rate = obs_rate[order]

    print(f"  ✓ Loaded {len(obs_rate)} DSAS points  |  "
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
    # Converts a GIS domain ID to the corresponding padded CASCADE index.
    # padded_index = START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)
    def _gis_to_pad(gis_id):
        return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)

    OVERWASH_FILTER = np.full(TOTAL_DOMAINS, OVERWASH_FILTER_DEFAULT)

    # Buxton: GIS 7–8 → padded 21–22
    OVERWASH_FILTER[_gis_to_pad(7)  : _gis_to_pad(8)  + 1] = OVERWASH_FILTER_BUXTON
    # Avon: GIS 21–31 → padded 35–45
    OVERWASH_FILTER[_gis_to_pad(21) : _gis_to_pad(31) + 1] = OVERWASH_FILTER_AVON
    # Salvo/Waves/Rodanthe: GIS 68–83 → padded 82–97
    OVERWASH_FILTER[_gis_to_pad(68) : _gis_to_pad(83) + 1] = OVERWASH_FILTER_SALVO_WAVES_ROD

    OVERWASH_FILTER = list(OVERWASH_FILTER)  # CASCADE expects a list

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

    # ── Load DSAS (multi-source), mapped to PADDED CASCADE indices ──────────
    dsas_series = []
    for ds in DSAS_DATASETS:
        print(f"Loading DSAS: {ds['label']}")
        x, r = load_and_map_dsas_to_cascade(
            dsas_csv=ds["csv"],
            domain_col=ds["domain_col"],
            rate_col=ds["rate_col"],
            start_real_index=START_REAL_INDEX,
            end_real_index=END_REAL_INDEX,
            first_file_number=FIRST_FILE_NUMBER,
            last_file_number=LAST_FILE_NUMBER,
            mode=ds.get("mode", "auto"),
            flip_alongshore=ds.get("flip_alongshore", False),
        )
        if x is None or r is None or len(r) == 0:
            continue
        if ds.get("flip_sign", False):
            r = -1.0 * r
        dsas_series.append(
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

        # ── Per-run output folder (same folder CASCADE already saved to) ─────
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

            for s in dsas_series:
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

            xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)
            ax.set_xticks(xticks)
            ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

            ax.set_xlabel("GIS Domain ID (1–90)")
            ax.set_ylabel("Shoreline change rate (m/yr)")
            ax.set_title(
                f"Modeled vs DSAS Shoreline Change Rate – Hatteras Island | "
                f"Real domains only | SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                f"{START_YEAR}–{END_YEAR} | Hs={Hs} m"
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

            for s in dsas_series:
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

            ax.set_ylabel("Shoreline change rate (m/yr)")
            ax.set_title(
                f"Modeled vs DSAS Shoreline Change Rate – Hatteras Island | "
                f"All domains (buffers included) | SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
                f"{START_YEAR}–{END_YEAR} | Hs={Hs} m"
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
