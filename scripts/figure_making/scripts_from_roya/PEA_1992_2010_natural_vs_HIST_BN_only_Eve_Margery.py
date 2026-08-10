#!/usr/bin/env python3
"""
PEA ISLAND: 1992–2010
Natural vs Historical Beach Nourishment Only

This version adds diagnostics for threshold-based beach nourishment:

1. Historical BN is still applied from manual year/domain/volume arrays.
2. Automatic beach-width threshold BN can be turned ON/OFF with:
       APPLY_AUTOMATIC_BN_THRESHOLD
3. The code logs:
       - where/when beach width falls below threshold
       - whether threshold BN was actually applied to CASCADE
       - how many managed domains triggered each year
4. The code also plots:
       - yearly relative shoreline change + BN
       - yearly absolute shoreline position + BN
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cascade.cascade import Cascade
import yaml


# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# =============================================================================

BERM_ELEVATION_M = 1.7
MHW_M = 0.36

NUM_REAL_DOMAINS = 41
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 80
LAST_FILE_NUMBER = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1

TOTAL_DOMAINS = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS

START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS

print("=" * 80)
print("PEA ISLAND CASCADE DOMAIN CONFIGURATION")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS} IDs {FIRST_FILE_NUMBER}..{LAST_FILE_NUMBER}")
print(f"Buffers:        {NUM_BUFFER_DOMAINS} each side")
print(f"TOTAL_DOMAINS:  {TOTAL_DOMAINS}")
print(f"Real index span: [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print("=" * 80 + "\n")


# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE"
PEA_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "PeaIsland_init")
OUTPUT_BASE_DIR = os.path.join(PROJECT_BASE_DIR, "comparison", "raw_runs", "Eve_Margery_BN_Sink","nosink_eve_margery")

DUNE_OFFSET_FILE = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/Island_offset/Island_Dune_Offsets_1992_2010_PADDED_71.csv"

STORM_FILE = os.path.join(
    PEA_DATA_BASE,
    "Storms",
    "/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/Storms/storms_1992_2010.npy",
)

ROAD_SETBACK_FILE = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/Roads/raw_offset/1978_road_setback_CLEAN.csv"

PARAMETER_FILE = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/data/PeaIsland_init/barrier3d-NoMngment-parameters.yaml"

COASTSAT_CSV = r"/Users/rsahrae/PycharmProjects/PeaIsland_Hindcast/CASCADE/scripts/input_prep/CoastSat/1992_2010/coastsat_smoothed_table.csv"

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)


# =============================================================================
# SECTION 3: SIMULATION PARAMETERS
# =============================================================================

START_YEAR = 1992
END_YEAR = 2010
RUN_YEARS = END_YEAR - START_YEAR

TO_METERS = True
SEA_LEVEL_RISE_RATE = 0.0048
SEA_LEVEL_CONSTANT = True

NUM_CORES = 4

DUNE_REBUILD_HEIGHT = 3.0
REBUILD_ELEV_THRESHOLD = 0.01
REBUILD_DUNE_THRESHOLD = 1.0

BEACH_WIDTH_THRESHOLD_VALUE = 30.0

# =============================================================================
# IMPORTANT SWITCH
# =============================================================================
# False = historical BN only.
# True  = historical BN + threshold-based BN diagnostic/application.
#
# To see whether threshold BN was effective, set this to True.
# To keep historical BN only, set this to False.
# =============================================================================

APPLY_AUTOMATIC_BN_THRESHOLD = False

FIXED_WAVE_PERIOD = 7
FIXED_WAVE_ASYMMETRY = 0.8
FIXED_WAVE_ANGLE_HIGH_FRACTION = 0.1
WAVE_HEIGHTS_TO_TEST = [1, 2]

DOMAIN_TICK_STEP = 5
FLIP_SIGN_MODEL = True

YEAR_COLUMN_INDEX = 0

DOMAIN_LENGTH_M = 500.0

MAKE_YEARLY_BN_GIF = True
GIF_DURATION_SECONDS = 4

print("Simulation Configuration:")
print(f"  Start Year: {START_YEAR} | End Year: {END_YEAR} | RUN_YEARS: {RUN_YEARS}")
print(f"  Storm file: {STORM_FILE}")
print(f"  WAVE_HEIGHTS: {WAVE_HEIGHTS_TO_TEST}")
print(f"  DUNE_REBUILD_HEIGHT: {DUNE_REBUILD_HEIGHT} m")
print(f"  REBUILD_ELEV_THRESHOLD: {REBUILD_ELEV_THRESHOLD} m")
print(f"  REBUILD_DUNE_THRESHOLD: {REBUILD_DUNE_THRESHOLD} m")
print(f"  BEACH_WIDTH_THRESHOLD_VALUE: {BEACH_WIDTH_THRESHOLD_VALUE} m")
print(f"  APPLY_AUTOMATIC_BN_THRESHOLD: {APPLY_AUTOMATIC_BN_THRESHOLD}")
print(f"  FLIP_SIGN_MODEL={FLIP_SIGN_MODEL}")
print("=" * 80 + "\n")


# =============================================================================
# SECTION 4: MANUAL HISTORICAL NOURISHMENT ARRAYS
# =============================================================================

# BN_YEARS_BY_DOMAIN = {
#     111: [2003],
#     112: [2003],
#     113: [2002, 2003, 2009],
#     114: [1992, 2002, 2003, 2009],
#     115: [1992, 2002, 2003, 2004],
#     116: [1992, 1993, 2001, 2002, 2003, 2004, 2008, 2009],
#     117: [1992, 1993, 1995, 2001, 2002, 2003, 2004, 2008, 2009],
#     118: [1992, 2001, 2003, 2004, 2008, 2009],
#     119: [2001, 2004, 2008, 2009],
#     120: [],
# }
# 
# BN_VOLUME_M3_BY_DOMAIN = {
#     111: [102954.25],
#     112: [102954.25],
# 
#     113: [96763.6, 102954.25, 236628.75],
# 
#     114: [431200.0, 96763.6, 102954.25, 236628.75],
# 
#     115: [431200.0, 96763.6, 102954.25, 74972.4],
# 
#     116: [
#         49146.6667,
#         173294.0,
#         102741.25,
#         96763.6,
#         102954.25,
#         74972.4,
#         316731.5,
#         236628.75,
#     ],
# 
#     117: [
#         49146.6667,
#         173294.0,
#         52184.8,
#         102741.25,
#         96763.6,
#         102954.25,
#         74972.4,
#         316731.5,
#         236628.75,
#     ],
# 
#     118: [
#         49146.6667,
#         102741.25,
#         102954.25,
#         74972.4,
#         187431.0,
#         316731.5,
#     ],
# 
#     119: [
#         102741.25,
#         74972.4,
#         187431.0,
#         316731.5,
#     ],
# }

# Combined Beach Nourishment volumes
# Units: m^3 per domain
# Sources added together:
#   1) screenshot/picture source
#   2) manual-array source from PeaIsland_BN_investigation.py comparison

bn_years = [1992, 1993, 1995, 2001, 2002, 2003, 2004, 2008, 2009]

bn_volume_by_domain_combined = {
    111: [0, 0, 0, 0, 0, 102954.250, 0, 0, 0],

    112: [0, 0, 0, 0, 0, 201347.013, 0, 0, 0],

    113: [0, 0, 0, 0, 96763.600, 201347.013, 0, 0, 263329.667],

    114: [201968.000, 0, 0, 0, 189240.027, 201347.013, 0, 0, 392555.166],

    115: [374105.997, 0, 0, 0, 189240.027, 201347.013, 74972.400, 0, 129225.499],

    116: [374105.997, 173294.000, 0, 102741.250, 189240.027, 201347.013, 146622.964, 373160.592, 392555.166],

    117: [374105.997, 338929.076, 52185.000, 200930.355, 189240.027, 201347.013, 146622.964, 373160.592, 392555.166],

    118: [374105.997, 165635.076, 49872.678, 200930.355, 92476.427, 201347.013, 146622.964, 373160.592, 392555.166],

    119: [0, 0, 0, 200930.355, 0, 98392.763, 146622.964, 373160.592, 392555.166],

    120: [0, 0, 0, 98189.105, 0, 0, 71650.564, 121079.342, 129225.499],
}
# =============================================================================
# SECTION 5: LOAD INPUT DATA
# =============================================================================

print("Loading input data...")

try:
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10.0

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=",")

    print(f"Loaded dune offsets: {dune_offset_dam.size} domains")
    print(f"Loaded road setbacks: {road_setbacks_raw.size} values")

except FileNotFoundError as e:
    print(f"CRITICAL ERROR: Missing data file: {e.filename}")
    sys.exit(1)
except Exception as e:
    print(f"CRITICAL ERROR loading data: {e}")
    sys.exit(1)

BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS
#
# BACKGROUND_EROSION_RATES = [ #buffers
#                             0,0,0,0,0,
#                             0,0,0,0,0,
#                             0,0,0,0,0,
#
#                             0,0,0,0,0,
#                             -0,0,0,0,0,
#                             -0,-0,0,0,0,
#                             0,0,-0,0,-0,
#                             0,0,0,0,0,
#                             0,0,0,0,0,
#                             0,0,0,0,0,
#                             0,0,0,0,-9,-18,
#
#                             #buffers:
#                             0,0,0,0,0,
#                             0,0,0,0,0,
#                             0,0,0,0,0]
# # =============================================================================
# ROAD SETBACK ARRAY
# =============================================================================

road_setbacks_full = np.zeros(TOTAL_DOMAINS)

num_road_values = min(len(road_setbacks_raw), NUM_REAL_DOMAINS - 1)

road_setbacks_full[
    START_REAL_INDEX:START_REAL_INDEX + num_road_values
] = road_setbacks_raw[:num_road_values]


# =============================================================================
# ROADWAY MANAGEMENT DOMAINS
# =============================================================================

ROADWAY_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
ROADWAY_MANAGEMENT_ON[START_REAL_INDEX:END_REAL_INDEX] = [False] * NUM_REAL_DOMAINS

NO_ROAD_REAL_IDS = [120]

for real_id in NO_ROAD_REAL_IDS:
    padded_index = START_REAL_INDEX + (real_id - FIRST_FILE_NUMBER)
    ROADWAY_MANAGEMENT_ON[padded_index] = False
    road_setbacks_full[padded_index] = 0.0

print("\nRoadway management summary:")
for real_id in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1):
    padded_index = START_REAL_INDEX + (real_id - FIRST_FILE_NUMBER)
    print(
        f"Domain {real_id}: "
        f"padded index={padded_index}, "
        f"roadway_management={ROADWAY_MANAGEMENT_ON[padded_index]}, "
        f"road_setback={road_setbacks_full[padded_index]:.3f}"
    )

print("=" * 80 + "\n")


# =============================================================================
# SECTION 6: ELEVATION + DUNE FILE LISTS
# =============================================================================

print("Generating elevation + dune file paths...")
#
# ELEVATION_FILE_PATHS = []
# DUNE_FILE_PATHS = []
#
# for _ in range(START_REAL_INDEX):
#     DUNE_FILE_PATHS.append(
#         os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_dune_Pea.npy")
#     )
#     ELEVATION_FILE_PATHS.append(
#         os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_topography_Pea.npy")
#     )
#
# for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
#     file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
#
#     DUNE_FILE_PATHS.append(
#         os.path.join(
#             PEA_DATA_BASE,
#             "Dunes",
#             "2011_modified",
#             f"resampled_2011_domains_{file_num}_dune.npy",
#         )
#     )
#
#     ELEVATION_FILE_PATHS.append(
#         os.path.join(
#             PEA_DATA_BASE,
#             "Topo",
#             "2011_modified",
#             f"resampled_2011_domains_{file_num}_topography.npy",
#         )
#     )
#
# for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
#     DUNE_FILE_PATHS.append(
#         os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_dune_Pea.npy")
#     )
#     ELEVATION_FILE_PATHS.append(
#         os.path.join(PEA_DATA_BASE, "Buffer", "sample_1_topography_Pea.npy")
#     )
#
# print(f"Generated {len(ELEVATION_FILE_PATHS)} elevation files")
# print(f"Generated {len(DUNE_FILE_PATHS)} dune files")
# print("=" * 80 + "\n")


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def real_domain_to_padded_index(real_domain):
    return START_REAL_INDEX + (real_domain - FIRST_FILE_NUMBER)


def padded_index_to_real_domain(padded_index):
    if START_REAL_INDEX <= padded_index < END_REAL_INDEX:
        return FIRST_FILE_NUMBER + (padded_index - START_REAL_INDEX)
    return None


def real_domain_to_x(real_domain_id):
    return START_REAL_INDEX + (real_domain_id - FIRST_FILE_NUMBER)


def get_x_s_TS(b3d):
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt = len(get_x_s_TS(b3d_list[0]))

    shoreline = np.zeros((nt, ndom), dtype=float)

    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])

    if to_meters:
        shoreline *= 10.0

    return shoreline


# def build_relative_shoreline_change_matrix(cascade, to_meters=True, flip_sign=True):
#     shoreline_m = build_shoreline_matrix(cascade, to_meters=to_meters)
#     shoreline_change_m = shoreline_m - shoreline_m[0, :]
#
#     if flip_sign:
#         shoreline_change_m *= -1.0
#
#     return shoreline_change_m


def load_and_map_coastsat_to_cascade(
    csv_path,
    domain_col,
    rate_col,
    start_real_index,
    end_real_index,
    first_file_number,
    last_file_number,
    mode="real_id",
    flip_alongshore=False,
):
    if csv_path is None or not os.path.exists(csv_path):
        print(f"CoastSat file not found: {csv_path}")
        return None, None

    df = pd.read_csv(csv_path)

    if domain_col not in df.columns or rate_col not in df.columns:
        raise ValueError(
            f"CoastSat CSV missing columns. Need '{domain_col}' and '{rate_col}'. "
            f"Found columns: {list(df.columns)}"
        )

    raw_dom = pd.to_numeric(df[domain_col], errors="coerce").to_numpy()
    raw_rate = pd.to_numeric(df[rate_col], errors="coerce").to_numpy()

    valid = np.isfinite(raw_dom) & np.isfinite(raw_rate)
    raw_dom = raw_dom[valid].astype(int)
    raw_rate = raw_rate[valid]

    if mode == "real_id":
        obs_x = start_real_index + (raw_dom - first_file_number)
    elif mode == "gis":
        obs_x = start_real_index + (raw_dom - 1)
    else:
        raise ValueError("mode must be 'gis' or 'real_id'")

    keep = (obs_x >= start_real_index) & (obs_x < end_real_index)
    obs_x = obs_x[keep].astype(int)
    obs_rate = raw_rate[keep]

    if len(obs_x) == 0:
        print("CoastSat mapped, but no points fell within real-domain range.")
        return None, None

    if flip_alongshore:
        obs_x = start_real_index + (end_real_index - 1 - obs_x)

    order = np.argsort(obs_x)
    obs_x = obs_x[order]
    obs_rate = obs_rate[order]

    print(f"Loaded CoastSat obs from {os.path.basename(csv_path)}: {len(obs_rate)} points")
    print(f"CoastSat rate range: {obs_rate.min():.2f} to {obs_rate.max():.2f} m/yr")
    print(f"CoastSat mapped CASCADE index range: {obs_x.min()} to {obs_x.max()}")

    return obs_x, obs_rate


def build_nourishment_arrays_from_manual_inputs():
    nourishment_on_by_year = {}
    nourishment_volume_by_year = {}

    for year in range(START_YEAR, END_YEAR + 1):
        nourishment_on_by_year[year] = np.zeros(TOTAL_DOMAINS)
        nourishment_volume_by_year[year] = [0.0] * TOTAL_DOMAINS

    for real_domain, volumes in bn_volume_by_domain_combined.items():

        years = bn_years

        if len(years) != len(volumes):
            raise ValueError(
                f"Domain {real_domain}: bn_years and volume list must have same length."
            )

        padded_index = real_domain_to_padded_index(real_domain)

        if padded_index < 0 or padded_index >= TOTAL_DOMAINS:
            print(f"Skipping domain {real_domain}: outside CASCADE range")
            continue

        for year, total_volume_m3 in zip(years, volumes):

            if year < START_YEAR or year > END_YEAR:
                continue

            volume_m3_per_m = float(total_volume_m3) / DOMAIN_LENGTH_M

            nourishment_on_by_year[year][padded_index] = 1
            nourishment_volume_by_year[year][padded_index] = volume_m3_per_m

    print("\nHistorical nourishment schedule from manual arrays:")
    for year in range(START_YEAR, END_YEAR + 1):
        active = np.where(nourishment_on_by_year[year] == 1)[0]
        if len(active) > 0:
            real_domains = [
                FIRST_FILE_NUMBER + (idx - START_REAL_INDEX)
                for idx in active
                if START_REAL_INDEX <= idx < END_REAL_INDEX
            ]

            total_volume = (
                np.sum(np.asarray(nourishment_volume_by_year[year], dtype=float))
                * DOMAIN_LENGTH_M
            )

            print(
                f"  {year}: domains {real_domains}, "
                f"total volume = {total_volume:,.0f} m3"
            )

    return nourishment_on_by_year, nourishment_volume_by_year


def add_site_annotations(ax):
    ymin, ymax = ax.get_ylim()
    xmin, xmax = ax.get_xlim()

    yrange = ymax - ymin
    xrange = xmax - xmin

    x_dom_80 = real_domain_to_x(80)
    x_dom_81 = real_domain_to_x(81)
    x_dom_110 = real_domain_to_x(110)
    x_dom_120 = real_domain_to_x(120)

    x_rodanthe = 0.5 * (x_dom_80 + x_dom_81)

    site_line_kwargs = dict(
        color="0.55",
        linestyle="--",
        linewidth=2.0,
        alpha=0.60,
        zorder=0,
    )

    ax.axvline(x_dom_120, **site_line_kwargs)
    ax.axvline(x_rodanthe, **site_line_kwargs)
    ax.axvline(x_dom_110, **site_line_kwargs)

    ax.text(
        x_rodanthe,
        ymax - 0.035 * yrange,
        "Rodanthe",
        fontsize=10,
        color="0.35",
        ha="center",
        va="top",
    )

    ax.text(
        x_dom_110,
        ymax - 0.035 * yrange,
        "Pea Island\nVisitor Center",
        fontsize=10,
        color="0.35",
        ha="center",
        va="top",
    )

    ax.text(
        x_dom_120,
        ymax - 0.035 * yrange,
        "Oregon Inlet\nGroin",
        fontsize=10,
        color="0.35",
        ha="center",
        va="top",
    )

    # ax.annotate(
    #     "Oregon Inlet\nTerminal Groin\nDomain 120",
    #     xy=(x_dom_120, ymax - 0.20 * yrange),
    #     xytext=(x_dom_120 - 10, ymax - 0.02 * yrange),
    #     fontsize=9,
    #     color="0.22",
    #     ha="center",
    #     va="top",
    #     arrowprops=dict(
    #         arrowstyle="-|>",
    #         color="0.35",
    #         lw=1.2,
    #         shrinkA=4,
    #         shrinkB=4,
    #         connectionstyle="arc3,rad=-0.18",
    #     ),
    #     bbox=dict(
    #         boxstyle="round,pad=0.25",
    #         fc="white",
    #         ec="0.70",
    #         alpha=0.88,
    #     ),
    # )

    # ax.annotate(
    #     "Rodanthe\nDomains 80–81",
    #     xy=(x_rodanthe, ymin + 0.35 * yrange),
    #     xytext=(x_rodanthe + 8, ymin + 0.15 * yrange),
    #     fontsize=9,
    #     color="0.22",
    #     ha="center",
    #     va="center",
    #     arrowprops=dict(
    #         arrowstyle="-|>",
    #         color="0.35",
    #         lw=1.2,
    #         shrinkA=4,
    #         shrinkB=4,
    #         connectionstyle="arc3,rad=0.22",
    #     ),
    #     bbox=dict(
    #         boxstyle="round,pad=0.25",
    #         fc="white",
    #         ec="0.70",
    #         alpha=0.88,
    #     ),
    # )
    #
    # ax.annotate(
    #     "Pea Island\nVisitor Center\nDomain 110",
    #     xy=(x_dom_110, ymax - 0.35 * yrange),
    #     xytext=(x_dom_110 - 8, ymax - 0.13 * yrange),
    #     fontsize=9,
    #     color="0.22",
    #     ha="center",
    #     va="top",
    #     arrowprops=dict(
    #         arrowstyle="-|>",
    #         color="0.35",
    #         lw=1.2,
    #         shrinkA=4,
    #         shrinkB=4,
    #         connectionstyle="arc3,rad=0.18",
    #     ),
    #     bbox=dict(
    #         boxstyle="round,pad=0.25",
    #         fc="white",
    #         ec="0.70",
    #         alpha=0.88,
    #     ),
    # )

    x_arrow = xmax - 0.035 * xrange

    ax.annotate(
        "Accretion  +",
        xy=(x_arrow, ymax - 0.20 * yrange),
        xytext=(x_arrow, ymax - 0.34 * yrange),
        fontsize=9,
        color="0.25",
        ha="center",
        va="center",
        arrowprops=dict(
            arrowstyle="-|>",
            color="0.35",
            lw=1.2,
        ),
    )

    ax.annotate(
        "Erosion  −",
        xy=(x_arrow, ymin + 0.20 * yrange),
        xytext=(x_arrow, ymin + 0.34 * yrange),
        fontsize=9,
        color="0.25",
        ha="center",
        va="center",
        arrowprops=dict(
            arrowstyle="-|>",
            color="0.35",
            lw=1.2,
        ),
    )


def add_real_domain_top_axis(ax):
    top_ax = ax.secondary_xaxis("top")
    top_ax.set_xlabel(
        "Real Pea Island domain ID",
        fontsize=10,
        fontweight="bold",
        labelpad=8,
    )

    top_positions = []
    top_labels = []

    for dom_id in range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP):
        j = START_REAL_INDEX + (dom_id - FIRST_FILE_NUMBER)
        top_positions.append(j)
        top_labels.append(str(dom_id))

    top_ax.set_xticks(top_positions)
    top_ax.set_xticklabels(top_labels, fontsize=9)

    return top_ax


# =============================================================================
# YEARLY PLOTTING FUNCTIONS
# =============================================================================

# def plot_yearly_relative_shoreline_and_bn(
#     cascades_by_label,
#     hist_nourish_volume_by_year,
#     output_dir,
#     filename_prefix,
#     flip_sign=True,
#     make_gif=True,
#     gif_duration_seconds=4,
# ):
#     yearly_dir = os.path.join(output_dir, filename_prefix)
#     os.makedirs(yearly_dir, exist_ok=True)
#
#     domain_numbers = np.arange(TOTAL_DOMAINS)
#
#     shoreline_change_by_label = {}
#     max_nt = None
#
#     for label, cascade in cascades_by_label.items():
#         shoreline_change_m = build_relative_shoreline_change_matrix(
#             cascade,
#             to_meters=True,
#             flip_sign=flip_sign,
#         )
#
#         shoreline_change_by_label[label] = shoreline_change_m
#
#         if max_nt is None:
#             max_nt = shoreline_change_m.shape[0]
#         else:
#             max_nt = min(max_nt, shoreline_change_m.shape[0])
#
#     if max_nt is None:
#         print("No cascades available for yearly plotting.")
#         return []
#
#     all_vals = []
#
#     for shoreline_change_m in shoreline_change_by_label.values():
#         vals = shoreline_change_m[:max_nt, :]
#         vals = vals[np.isfinite(vals)]
#         all_vals.extend(vals)
#
#     if len(all_vals) > 0:
#         y_min = np.nanmin(all_vals)
#         y_max = np.nanmax(all_vals)
#
#         if y_max > y_min:
#             y_pad = 0.15 * (y_max - y_min)
#         else:
#             y_pad = 5.0
#
#         y_min -= y_pad
#         y_max += y_pad
#     else:
#         y_min, y_max = -10.0, 10.0
#
#     all_bn_total_m3 = []
#
#     for year in range(START_YEAR, END_YEAR + 1):
#         if year in hist_nourish_volume_by_year:
#             vols_total_m3 = (
#                 np.asarray(hist_nourish_volume_by_year[year], dtype=float)
#                 * DOMAIN_LENGTH_M
#             )
#             all_bn_total_m3.extend(vols_total_m3[vols_total_m3 > 0])
#
#     if len(all_bn_total_m3) > 0:
#         max_bn_volume = max(all_bn_total_m3)
#     else:
#         max_bn_volume = 1.0
#
#     preferred_bn_label = "Historical BN only"
#
#     if preferred_bn_label in shoreline_change_by_label:
#         bn_curve_label = preferred_bn_label
#     else:
#         bn_curve_label = list(shoreline_change_by_label.keys())[-1]
#
#     png_files = []
#
#     for t in range(max_nt):
#         model_year = START_YEAR + t
#
#         if model_year in hist_nourish_volume_by_year:
#             bn_total_m3 = (
#                 np.asarray(hist_nourish_volume_by_year[model_year], dtype=float)
#                 * DOMAIN_LENGTH_M
#             )
#         else:
#             bn_total_m3 = np.zeros(TOTAL_DOMAINS)
#
#         active_bn = np.where(bn_total_m3 > 0)[0]
#
#         fig, (ax, ax_bn) = plt.subplots(
#             2,
#             1,
#             figsize=(14, 8),
#             sharex=True,
#             gridspec_kw={"height_ratios": [3.2, 1.15]},
#             constrained_layout=True,
#         )
#
#         for label, shoreline_change_m in shoreline_change_by_label.items():
#             ax.plot(
#                 domain_numbers,
#                 shoreline_change_m[t, :],
#                 linewidth=2.2,
#                 label=label,
#             )
#
#         ax.axvline(START_REAL_INDEX, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
#         ax.axvline(END_REAL_INDEX - 1, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
#         ax.axhline(0.0, linestyle="--", linewidth=1.2, color="0.25", alpha=0.75)
#
#         if len(active_bn) > 0:
#             x_bn = []
#             y_bn = []
#
#             for idx in active_bn:
#                 idx = int(idx)
#
#                 x_bn.append(idx)
#                 y_bn.append(shoreline_change_by_label[bn_curve_label][t, idx])
#
#                 ax.axvline(
#                     idx,
#                     linestyle=":",
#                     linewidth=1.0,
#                     color="red",
#                     alpha=0.45,
#                     zorder=2,
#                 )
#
#             ax.scatter(
#                 x_bn,
#                 y_bn,
#                 s=35,
#                 marker="o",
#                 color="red",
#                 edgecolor="black",
#                 linewidth=0.6,
#                 zorder=10,
#                 label="Historical BN domain",
#             )
#
#         ax.set_xlim(0, TOTAL_DOMAINS - 1)
#         ax.set_ylim(y_min, y_max)
#
#         if flip_sign:
#             ax.set_ylabel(
#                 f"Relative shoreline change from {START_YEAR} (m)\n"
#                 "Accretion = positive, erosion = negative",
#                 fontsize=11,
#                 fontweight="bold",
#             )
#         else:
#             ax.set_ylabel(
#                 f"Raw shoreline-position change from {START_YEAR} (m)",
#                 fontsize=11,
#                 fontweight="bold",
#             )
#
#         ax.set_title(
#             f"Pea Island relative shoreline change with historical BN | {model_year}",
#             fontsize=14,
#             fontweight="bold",
#             pad=12,
#         )
#
#         ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
#         ax.legend(loc="best", fontsize=9, frameon=True, framealpha=0.92)
#
#         add_site_annotations(ax)
#         add_real_domain_top_axis(ax)
#
#         ax_bn.bar(
#             domain_numbers,
#             bn_total_m3,
#             width=0.8,
#             color="tab:blue",
#             alpha=0.85,
#         )
#
#         total_bn_this_year = np.sum(bn_total_m3)
#
#         if len(active_bn) > 0:
#             bn_title = (
#                 f"Historical beach nourishment in {model_year}: YES | "
#                 f"Total volume = {total_bn_this_year:,.0f} m³"
#             )
#         else:
#             bn_title = f"Historical beach nourishment in {model_year}: NO"
#
#         ax_bn.set_title(bn_title, fontsize=12, fontweight="bold")
#         ax_bn.set_ylabel("Historical BN\nvolume\n(m³/domain)", fontsize=10, fontweight="bold")
#         ax_bn.set_xlabel("CASCADE domain index including buffers", fontsize=11, fontweight="bold")
#         ax_bn.set_ylim(0, max_bn_volume * 1.25)
#         ax_bn.grid(True, axis="y", linestyle=":", linewidth=0.8, alpha=0.35)
#
#         xticks = np.arange(0, TOTAL_DOMAINS, DOMAIN_TICK_STEP)
#         ax_bn.set_xticks(xticks)
#         ax_bn.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)
#
#         for idx in active_bn:
#             idx = int(idx)
#
#             real_domain = padded_index_to_real_domain(idx)
#             volume = bn_total_m3[idx]
#
#             ax_bn.text(
#                 idx,
#                 volume + max_bn_volume * 0.025,
#                 f"D{real_domain}\n{volume/1000:.0f}k",
#                 ha="center",
#                 va="bottom",
#                 fontsize=8,
#                 color="0.25",
#             )
#
#         fig_out = os.path.join(yearly_dir, f"{filename_prefix}_{model_year}.png")
#         fig.savefig(fig_out, dpi=220, bbox_inches="tight")
#         plt.close(fig)
#
#         png_files.append(fig_out)
#
#     print("\nSaved yearly relative shoreline + historical BN plots to:")
#     print(yearly_dir)
#
#     if make_gif:
#         try:
#             import imageio.v2 as imageio
#
#             images = [imageio.imread(f) for f in png_files]
#             gif_out = os.path.join(output_dir, f"{filename_prefix}.gif")
#
#             imageio.mimsave(
#                 gif_out,
#                 images,
#                 duration=gif_duration_seconds,
#                 loop=0,
#             )
#
#             print("\nSaved yearly relative shoreline + historical BN GIF:")
#             print(gif_out)
#
#         except ImportError:
#             print("\nimageio is not installed. GIF was not created.")
#             print("Install it with:")
#             print("pip install imageio")
#
#     return png_files
#
#
# def plot_yearly_absolute_shoreline_and_bn(
#     cascades_by_label,
#     hist_nourish_volume_by_year,
#     output_dir,
#     filename_prefix,
#     make_gif=True,
#     gif_duration_seconds=4,
# ):
#     yearly_dir = os.path.join(output_dir, filename_prefix)
#     os.makedirs(yearly_dir, exist_ok=True)
#
#     domain_numbers = np.arange(TOTAL_DOMAINS)
#
#     shoreline_position_by_label = {}
#     max_nt = None
#
#     for label, cascade in cascades_by_label.items():
#         shoreline_m = build_shoreline_matrix(
#             cascade,
#             to_meters=True,
#         )
#
#         shoreline_position_by_label[label] = shoreline_m
#
#         if max_nt is None:
#             max_nt = shoreline_m.shape[0]
#         else:
#             max_nt = min(max_nt, shoreline_m.shape[0])
#
#     if max_nt is None:
#         print("No cascades available for absolute shoreline plotting.")
#         return []
#
#     all_vals = []
#
#     for shoreline_m in shoreline_position_by_label.values():
#         vals = shoreline_m[:max_nt, :]
#         vals = vals[np.isfinite(vals)]
#         all_vals.extend(vals)
#
#     if len(all_vals) > 0:
#         y_min = np.nanmin(all_vals)
#         y_max = np.nanmax(all_vals)
#
#         if y_max > y_min:
#             y_pad = 0.08 * (y_max - y_min)
#         else:
#             y_pad = 50.0
#
#         y_min -= y_pad
#         y_max += y_pad
#     else:
#         y_min, y_max = -6000.0, -5000.0
#
#     all_bn_total_m3 = []
#
#     for year in range(START_YEAR, END_YEAR + 1):
#         if year in hist_nourish_volume_by_year:
#             vols_total_m3 = (
#                 np.asarray(hist_nourish_volume_by_year[year], dtype=float)
#                 * DOMAIN_LENGTH_M
#             )
#             all_bn_total_m3.extend(vols_total_m3[vols_total_m3 > 0])
#
#     if len(all_bn_total_m3) > 0:
#         max_bn_volume = max(all_bn_total_m3)
#     else:
#         max_bn_volume = 1.0
#
#     preferred_bn_label = "Historical BN only"
#
#     if preferred_bn_label in shoreline_position_by_label:
#         bn_curve_label = preferred_bn_label
#     else:
#         bn_curve_label = list(shoreline_position_by_label.keys())[-1]
#
#     png_files = []
#
#     for t in range(max_nt):
#         model_year = START_YEAR + t
#
#         if model_year in hist_nourish_volume_by_year:
#             bn_total_m3 = (
#                 np.asarray(hist_nourish_volume_by_year[model_year], dtype=float)
#                 * DOMAIN_LENGTH_M
#             )
#         else:
#             bn_total_m3 = np.zeros(TOTAL_DOMAINS)
#
#         active_bn = np.where(bn_total_m3 > 0)[0]
#
#         fig, (ax, ax_bn) = plt.subplots(
#             2,
#             1,
#             figsize=(14, 6),
#             sharex=True,
#             gridspec_kw={"height_ratios": [3.2, 1.15]},
#             constrained_layout=True,
#         )
#
#         for label, shoreline_m in shoreline_position_by_label.items():
#             ax.plot(
#                 domain_numbers,
#                 shoreline_m[t, :],
#                 linewidth=2.2,
#                 label=label,
#             )
#
#         ax.axvline(START_REAL_INDEX, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
#         ax.axvline(END_REAL_INDEX - 1, linestyle="--", linewidth=1.0, color="0.45", alpha=0.70)
#
#         if len(active_bn) > 0:
#             x_bn = []
#             y_bn = []
#
#             for idx in active_bn:
#                 idx = int(idx)
#
#                 x_bn.append(idx)
#                 y_bn.append(shoreline_position_by_label[bn_curve_label][t, idx])
#
#                 ax.axvline(
#                     idx,
#                     linestyle=":",
#                     linewidth=1.0,
#                     color="red",
#                     alpha=0.45,
#                     zorder=2,
#                 )
#
#             ax.scatter(
#                 x_bn,
#                 y_bn,
#                 s=35,
#                 marker="o",
#                 color="red",
#                 edgecolor="black",
#                 linewidth=0.6,
#                 zorder=10,
#                 label="Historical BN domain",
#             )
#
#         ax.set_xlim(0, TOTAL_DOMAINS - 1)
#         ax.set_ylim(y_min, y_max)
#
#         ax.set_ylabel(
#             "Absolute shoreline position, x_s (m)",
#             fontsize=11,
#             fontweight="bold",
#         )
#
#         ax.set_title(
#             f"Pea Island absolute shoreline position with historical BN | {model_year}",
#             fontsize=14,
#             fontweight="bold",
#             pad=12,
#         )
#
#         ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
#         ax.legend(loc="best", fontsize=9, frameon=True, framealpha=0.92)
#
#         add_site_annotations(ax)
#         add_real_domain_top_axis(ax)
#
#         ax_bn.bar(
#             domain_numbers,
#             bn_total_m3,
#             width=0.8,
#             color="tab:blue",
#             alpha=0.85,
#         )
#
#         total_bn_this_year = np.sum(bn_total_m3)
#
#         if len(active_bn) > 0:
#             bn_title = (
#                 f"Historical beach nourishment in {model_year}: YES | "
#                 f"Total volume = {total_bn_this_year:,.0f} m³"
#             )
#         else:
#             bn_title = f"Historical beach nourishment in {model_year}: NO"
#
#         ax_bn.set_title(bn_title, fontsize=12, fontweight="bold")
#         ax_bn.set_ylabel("Historical BN\nvolume\n(m³/domain)", fontsize=10, fontweight="bold")
#         ax_bn.set_xlabel("CASCADE domain index including buffers", fontsize=11, fontweight="bold")
#         ax_bn.set_ylim(0, max_bn_volume * 1.25)
#         ax_bn.grid(True, axis="y", linestyle=":", linewidth=0.8, alpha=0.35)
#
#         xticks = np.arange(0, TOTAL_DOMAINS, DOMAIN_TICK_STEP)
#         ax_bn.set_xticks(xticks)
#         ax_bn.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)
#
#         for idx in active_bn:
#             idx = int(idx)
#
#             real_domain = padded_index_to_real_domain(idx)
#             volume = bn_total_m3[idx]
#
#             ax_bn.text(
#                 idx,
#                 volume + max_bn_volume * 0.025,
#                 f"D{real_domain}\n{volume/1000:.0f}k",
#                 ha="center",
#                 va="bottom",
#                 fontsize=8,
#                 color="0.25",
#             )
#
#         fig_out = os.path.join(yearly_dir, f"{filename_prefix}_{model_year}.png")
#         fig.savefig(fig_out, dpi=220, bbox_inches="tight")
#         plt.close(fig)
#
#         png_files.append(fig_out)
#
#     print("\nSaved yearly absolute shoreline position + historical BN plots to:")
#     print(yearly_dir)
#
#     if make_gif:
#         try:
#             import imageio.v2 as imageio
#
#             images = [imageio.imread(f) for f in png_files]
#             gif_out = os.path.join(output_dir, f"{filename_prefix}.gif")
#
#             imageio.mimsave(
#                 gif_out,
#                 images,
#                 duration=gif_duration_seconds,
#                 loop=0,
#             )
#
#             print("\nSaved yearly absolute shoreline position + historical BN GIF:")
#             print(gif_out)
#
#         except ImportError:
#             print("\nimageio is not installed. Absolute shoreline GIF was not created.")
#             print("Install it with:")
#             print("pip install imageio")
#
#     return png_files


# =============================================================================
# CASCADE RUNNER
# =============================================================================

def run_cascade_simulation(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    beach_width_threshold,
    rmin,
    rmax,
    elevation_file,
    dune_file,
    dune_design_elevation,
    dune_minimum_elevation,
    road_ele,
    road_width,
    road_setback,
    overwash_filter,
    overwash_to_dune,
    nourishment_volume,
    background_erosion,
    rebuild_dune_threshold,
    roadway_management_on,
    beach_dune_manager_on,
    sea_level_rise_rate,
    sea_level_constant,
    enable_shoreline_offset,
    shoreline_offset,
    wave_height,
    wave_period,
    wave_asymmetry,
    wave_angle_high_fraction,
    historical_nourishment_on_by_year=None,
    historical_nourishment_volume_by_year=None,
):
    datadir = os.path.join(PROJECT_BASE_DIR, "data", "PeaIsland_init")

    dune_rebuild_log = []
    historical_nourishment_log = []
    threshold_bn_crossing_log = []
    threshold_bn_effective_log = []

    with open(PARAMETER_FILE, "r") as f:
        yaml_params = yaml.safe_load(f)

    print("\nPARAMETER CHECK:")
    print(f"  YAML file: {PARAMETER_FILE}")
    print(f"  YAML MHW: {yaml_params.get('MHW')}")
    print(f"  YAML BermEl: {yaml_params.get('BermEl')}")
    print(f"  CASCADE MHW argument: {MHW_M}")
    print(f"  CASCADE berm_elevation argument: {BERM_ELEVATION_M}")
    print(f"  APPLY_AUTOMATIC_BN_THRESHOLD: {APPLY_AUTOMATIC_BN_THRESHOLD}")

    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=PARAMETER_FILE,

        berm_elevation=BERM_ELEVATION_M,
        MHW=MHW_M,

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
        alongshore_transport_module=True,

        beach_nourishment_module=beach_dune_manager_on,

        community_economics_module=False,

        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,

        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,

        nourishment_interval=None,
        nourishment_volume=nourishment_volume,

        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,

        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
    )

    print("\nBarrier3D berm elevation check:")
    print("cascade.barrier3d[0].BermEl raw =", cascade.barrier3d[0].BermEl)
    print("cascade.barrier3d[0].BermEl meters =", cascade.barrier3d[0].BermEl * 10)

    for i in range(START_REAL_INDEX, END_REAL_INDEX):
        real_id = FIRST_FILE_NUMBER + (i - START_REAL_INDEX)
        print(
            f"Domain {real_id}: "
            f"BermEl_raw={cascade.barrier3d[i].BermEl}, "
            f"BermEl_m={cascade.barrier3d[i].BermEl * 10:.3f}"
        )

    cascade._sandbag_management_on = [False] * alongshore_section_count
    cascade._sandbag_elevation = [0.0] * alongshore_section_count

    dune_rebuild_threshold_abs = rebuild_dune_threshold + (
        cascade.barrier3d[0].BermEl * 10
    )

    print(f"\n{name}: absolute dune rebuild threshold = {dune_rebuild_threshold_abs:.3f} m MHW")

    managed_mask = np.asarray(beach_dune_manager_on, dtype=bool)
    n_managed = int(np.sum(managed_mask))

    for time_step in range(nt - 1):
        current_year = START_YEAR + time_step

        cascade.nourish_now = np.zeros(alongshore_section_count)
        cascade.rebuild_dune_now = np.zeros(alongshore_section_count)

        # ---------------------------------------------------------------------
        # HISTORICAL BN
        # ---------------------------------------------------------------------
        if historical_nourishment_on_by_year is not None:
            if current_year in historical_nourishment_on_by_year:

                nourish_now = np.asarray(
                    historical_nourishment_on_by_year[current_year],
                    dtype=float,
                )

                nourish_volume = np.asarray(
                    historical_nourishment_volume_by_year[current_year],
                    dtype=float,
                )

                if np.any(nourish_now == 1):
                    print(f"\n{name}: applying HISTORICAL BN in {current_year}")

                    cascade.nourish_now = nourish_now

                    for iB3D in range(alongshore_section_count):
                        if nourish_now[iB3D] == 1:

                            if hasattr(cascade.nourishments[iB3D], "_nourishment_volume"):
                                cascade.nourishments[iB3D]._nourishment_volume = nourish_volume[iB3D]
                            elif hasattr(cascade.nourishments[iB3D], "nourishment_volume"):
                                cascade.nourishments[iB3D].nourishment_volume = nourish_volume[iB3D]
                            else:
                                raise AttributeError(
                                    "Could not find nourishment volume attribute on "
                                    f"cascade.nourishments[{iB3D}]"
                                )

                            real_domain = padded_index_to_real_domain(iB3D)

                            historical_nourishment_log.append(
                                dict(
                                    run_name=name,
                                    model_year=current_year,
                                    time_step=time_step,
                                    padded_index=iB3D,
                                    real_domain=real_domain,
                                    nourishment_type="historical",
                                    nourishment_volume_m3_per_m=float(nourish_volume[iB3D]),
                                    nourishment_volume_total_m3=float(nourish_volume[iB3D]) * DOMAIN_LENGTH_M,
                                )
                            )

        print(f"\r{name}: Time Step {time_step}", end="", flush=True)

        cascade.update()

        if cascade.b3d_break:
            break

        t = cascade.barrier3d[0].time_index

        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_threshold_nourish_now = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):

            if cascade.community_break[iB3D]:
                continue

            if not beach_dune_manager_on[iB3D]:
                continue

            # -----------------------------------------------------------------
            # THRESHOLD BN DIAGNOSTIC
            # -----------------------------------------------------------------
            if APPLY_AUTOMATIC_BN_THRESHOLD:

                beach_width_now = float(cascade.nourishments[iB3D].beach_width[t - 1])
                threshold_now = float(beach_width_threshold[iB3D])

                if beach_width_now < threshold_now:

                    tmp_threshold_nourish_now[iB3D] = 1

                    real_domain = padded_index_to_real_domain(iB3D)

                    threshold_bn_crossing_log.append(
                        dict(
                            run_name=name,
                            model_year=START_YEAR + time_step + 1,
                            time_step=time_step + 1,
                            padded_index=iB3D,
                            real_domain=real_domain,
                            beach_width_m=beach_width_now,
                            beach_width_threshold_m=threshold_now,
                            threshold_crossed=True,
                        )
                    )

            # -----------------------------------------------------------------
            # DUNE REBUILD LOGIC
            # -----------------------------------------------------------------
            DuneDomainCrest = (
                cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            )

            DuneCrestMin = (
                np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
            ) * 10

            if DuneCrestMin < dune_rebuild_threshold_abs:
                tmp_rebuild_dune[iB3D] = 1

                real_domain = padded_index_to_real_domain(iB3D)

                dune_rebuild_log.append(
                    dict(
                        run_name=name,
                        model_year=START_YEAR + time_step + 1,
                        time_step=time_step + 1,
                        padded_index=iB3D,
                        real_domain=real_domain,
                        dune_crest_min_m=float(DuneCrestMin),
                        rebuild_threshold_m=float(dune_rebuild_threshold_abs),
                        roadway_management_on=bool(roadway_management_on[iB3D]),
                        beach_dune_manager_on=bool(beach_dune_manager_on[iB3D]),
                        road_setback=float(road_setback[iB3D]),
                    )
                )

        # ---------------------------------------------------------------------
        # APPLY THRESHOLD BN USING CURRENT OCRACOKE-STYLE np.all LOGIC
        # ---------------------------------------------------------------------
        if APPLY_AUTOMATIC_BN_THRESHOLD and n_managed > 0:

            n_triggered = int(np.sum(tmp_threshold_nourish_now[managed_mask] == 1))

            threshold_bn_applied = False

            if np.all(tmp_threshold_nourish_now[managed_mask] == 1):
                cascade.nourish_now = tmp_threshold_nourish_now
                threshold_bn_applied = True

            active_threshold_bn = np.where(tmp_threshold_nourish_now == 1)[0]

            for idx in active_threshold_bn:
                idx = int(idx)

                real_domain = padded_index_to_real_domain(idx)

                if hasattr(cascade.nourishments[idx], "beach_width"):
                    beach_width_logged = float(cascade.nourishments[idx].beach_width[t - 1])
                else:
                    beach_width_logged = np.nan

                threshold_bn_effective_log.append(
                    dict(
                        run_name=name,
                        model_year=START_YEAR + time_step + 1,
                        time_step=time_step + 1,
                        padded_index=idx,
                        real_domain=real_domain,
                        beach_width_m=beach_width_logged,
                        beach_width_threshold_m=float(beach_width_threshold[idx]),
                        threshold_crossed=True,
                        n_managed_domains=n_managed,
                        n_triggered_domains=n_triggered,
                        threshold_bn_applied_to_model=threshold_bn_applied,
                        applied_rule="np.all managed domains below threshold",
                    )
                )

            if n_triggered > 0:
                print(
                    f"\n{name}: Threshold BN diagnostic year {START_YEAR + time_step + 1}: "
                    f"{n_triggered}/{n_managed} managed domains below threshold | "
                    f"applied_to_model={threshold_bn_applied}"
                )

        # ---------------------------------------------------------------------
        # APPLY DUNE REBUILD USING OCRACOKE-STYLE np.all LOGIC
        # ---------------------------------------------------------------------
        if n_managed > 0:
            if np.all(tmp_rebuild_dune[managed_mask] == 1):
                cascade.rebuild_dune_now = tmp_rebuild_dune

    save_path = os.path.join(OUTPUT_BASE_DIR, name)
    os.makedirs(save_path, exist_ok=True)

    cascade.save(save_path)
    print(f"\nSaved: {save_path}")

    cascade.dune_rebuild_log = dune_rebuild_log
    cascade.historical_nourishment_log = historical_nourishment_log
    cascade.threshold_bn_crossing_log = threshold_bn_crossing_log
    cascade.threshold_bn_effective_log = threshold_bn_effective_log

    # -------------------------------------------------------------------------
    # SAVE DUNE REBUILD LOG
    # -------------------------------------------------------------------------
    if len(dune_rebuild_log) > 0:
        rebuild_df = pd.DataFrame(dune_rebuild_log)
        rebuild_csv = os.path.join(
            OUTPUT_BASE_DIR,
            f"{name}_dune_rebuild_log.csv",
        )
        rebuild_df.to_csv(rebuild_csv, index=False)

        print("\n" + "=" * 80)
        print(f"DUNE REBUILD CONDITION SUMMARY FOR RUN: {name}")
        print("=" * 80)
        print(f"Number of domain-year records below threshold: {len(rebuild_df)}")
        print("\nBy year:")
        print(rebuild_df.groupby("model_year").size())
        print("\nBy real domain:")
        print(rebuild_df.groupby("real_domain").size())
        print(f"\nSaved dune rebuild log: {rebuild_csv}")
    else:
        print(f"\n{name}: No dune rebuild condition was detected.")

    # -------------------------------------------------------------------------
    # SAVE HISTORICAL BN LOG
    # -------------------------------------------------------------------------
    if len(historical_nourishment_log) > 0:
        hist_bn_df = pd.DataFrame(historical_nourishment_log)
        hist_bn_csv = os.path.join(
            OUTPUT_BASE_DIR,
            f"{name}_historical_BN_log.csv",
        )
        hist_bn_df.to_csv(hist_bn_csv, index=False)

        print("\n" + "=" * 80)
        print(f"HISTORICAL BN SUMMARY FOR RUN: {name}")
        print("=" * 80)
        print(f"Number of domain-year historical BN records: {len(hist_bn_df)}")
        print("\nBy year:")
        print(hist_bn_df.groupby("model_year").size())
        print("\nBy real domain:")
        print(hist_bn_df.groupby("real_domain").size())
        print(f"\nSaved historical BN log: {hist_bn_csv}")
    else:
        print(f"\n{name}: No historical BN was applied.")

    # -------------------------------------------------------------------------
    # SAVE THRESHOLD CROSSING LOG
    # -------------------------------------------------------------------------
    if len(threshold_bn_crossing_log) > 0:
        threshold_cross_df = pd.DataFrame(threshold_bn_crossing_log)
        threshold_cross_csv = os.path.join(
            OUTPUT_BASE_DIR,
            f"{name}_threshold_BN_crossing_log.csv",
        )
        threshold_cross_df.to_csv(threshold_cross_csv, index=False)

        print("\n" + "=" * 80)
        print(f"THRESHOLD BN CROSSING SUMMARY FOR RUN: {name}")
        print("=" * 80)
        print(f"Number of threshold-crossing domain-year records: {len(threshold_cross_df)}")
        print("\nBy year:")
        print(threshold_cross_df.groupby("model_year").size())
        print("\nBy real domain:")
        print(threshold_cross_df.groupby("real_domain").size())
        print(f"\nSaved threshold BN crossing log: {threshold_cross_csv}")
    else:
        print(f"\n{name}: No beach-width threshold crossings were detected.")

    # -------------------------------------------------------------------------
    # SAVE THRESHOLD EFFECTIVE LOG
    # -------------------------------------------------------------------------
    if len(threshold_bn_effective_log) > 0:
        threshold_effective_df = pd.DataFrame(threshold_bn_effective_log)

        threshold_effective_csv = os.path.join(
            OUTPUT_BASE_DIR,
            f"{name}_threshold_BN_effective_log.csv",
        )

        threshold_effective_df.to_csv(threshold_effective_csv, index=False)

        print("\n" + "=" * 80)
        print(f"THRESHOLD BN EFFECTIVENESS SUMMARY FOR RUN: {name}")
        print("=" * 80)

        print(f"Number of threshold-crossing domain-year records: {len(threshold_effective_df)}")

        print("\nBy year:")
        print(
            threshold_effective_df
            .groupby("model_year")
            .agg(
                triggered_domains=("padded_index", "count"),
                managed_domains=("n_managed_domains", "max"),
                applied_to_model=("threshold_bn_applied_to_model", "max"),
            )
        )

        print("\nBy real domain:")
        print(
            threshold_effective_df
            .groupby("real_domain")
            .agg(
                times_triggered=("model_year", "count"),
                ever_applied_to_model=("threshold_bn_applied_to_model", "max"),
            )
        )

        print(f"\nSaved threshold BN effectiveness log: {threshold_effective_csv}")

        # plot_threshold_bn_effectiveness(
        #     threshold_bn_effective_log=threshold_bn_effective_log,
        #     output_dir=OUTPUT_BASE_DIR,
        #     run_name=name,
        # )

    else:
        print(f"\n{name}: No threshold BN effectiveness records were created.")
        print(f"APPLY_AUTOMATIC_BN_THRESHOLD={APPLY_AUTOMATIC_BN_THRESHOLD}")

    return cascade


def build_elevation_dune_file_paths(dune_year):
    elevation_paths = []
    dune_paths = []

    # left buffers
    for _ in range(START_REAL_INDEX):
        dune_paths.append(
            os.path.join(
                PEA_DATA_BASE,
                "Buffer",
                "sample_1_dune_Pea.npy",
            )
        )

        elevation_paths.append(
            os.path.join(
                PEA_DATA_BASE,
                "Buffer",
                "sample_1_topography_Pea.npy",
            )
        )

    # real domains
    for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
        file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)

        # DUNES → 1996 OR 2011
        dune_paths.append(
            os.path.join(
                PEA_DATA_BASE,
                "Dunes",
                f"{dune_year}_dunes",
                f"resampled_{dune_year}_domains_{file_num}_dune.npy",
            )
        )

        # TOPO → ALWAYS 2011
        elevation_paths.append(
            os.path.join(
                PEA_DATA_BASE,
                "Topo",
                "2011_modified",
                f"resampled_2011_domains_{file_num}_topography.npy",
            )
        )

    # right buffers
    for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
        dune_paths.append(
            os.path.join(
                PEA_DATA_BASE,
                "Buffer",
                "sample_1_dune_Pea.npy",
            )
        )

        elevation_paths.append(
            os.path.join(
                PEA_DATA_BASE,
                "Buffer",
                "sample_1_topography_Pea.npy",
            )
        )

    return elevation_paths, dune_paths


# =============================================================================
# MAIN
# =============================================================================

def main():
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS

    DUNE_DESIGN_ELEVATION = [DUNE_REBUILD_HEIGHT] * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    BEACH_WIDTH_THRESHOLD = [BEACH_WIDTH_THRESHOLD_VALUE] * TOTAL_DOMAINS

    ROAD_ELEVATION_REAL = [
        1.0776, 1.2634, 1.4387, 1.1686, 1.3991,
        1.2738, 1.2564, 1.3914, 1.8516, 1.3563,
        0.8953, 0.7497, 0.9807, 0.7779, 0.8456,
        0.8496, 0.9289, 0.9280, 1.0394, 1.3601,
        1.2931, 1.4435, 1.5030, 1.4476, 1.3255,
        1.3430, 1.1604, 1.6557, 1.0107, 1.1600,
        1.2845, 0.9283, 1.8228, 1.4769, 1.2812,
        1.2199, 1.1808, 0.9350, 1.3043, 1.6772,
        1.0,
    ]

    ROAD_ELEVATION = [0.0] * TOTAL_DOMAINS
    ROAD_ELEVATION[START_REAL_INDEX:END_REAL_INDEX] = ROAD_ELEVATION_REAL

    ROAD_WIDTH = [10.0] * TOTAL_DOMAINS

    OVERWASH_FILTER = 0
    OVERWASH_TO_DUNE = 9

    NOURISHMENT_VOLUME = 100

    time_span_years = END_YEAR - START_YEAR

    HIST_NOURISH_ON, HIST_NOURISH_VOLUME = build_nourishment_arrays_from_manual_inputs()

    x, r = load_and_map_coastsat_to_cascade(
        csv_path=COASTSAT_CSV,
        domain_col="domain",
        rate_col="cs_lrr_smooth_1992_2010",
        start_real_index=START_REAL_INDEX,
        end_real_index=END_REAL_INDEX,
        first_file_number=FIRST_FILE_NUMBER,
        last_file_number=LAST_FILE_NUMBER,
        mode="real_id",
        flip_alongshore=False,
    )

    coastsat_series = []
    if x is not None and r is not None:
        coastsat_series.append(
            dict(
                label="CoastSat LRR 1992–2010",
                x=x,
                rate=r,
            )
        )

    # -------------------------------------------------------------------------
    # SCENARIOS TO RUN
    # -------------------------------------------------------------------------
    # This script now runs only:
    #   1) Natural
    #   2) Historical beach nourishment only
    #
    # Roadway management is OFF in both scenarios.
    # For the BN-only scenario, the beach/dune manager is ON only for the
    # real domains that have historical nourishment entries. This allows
    # historical BN to be applied without turning on roadway management.
    # -------------------------------------------------------------------------

    NATURAL_ROADWAY_ON = [False] * TOTAL_DOMAINS
    NATURAL_BEACH_DUNE_ON = [False] * TOTAL_DOMAINS

    BN_ONLY_ROADWAY_ON = [False] * TOTAL_DOMAINS
    BN_ONLY_BEACH_DUNE_ON = [False] * TOTAL_DOMAINS

    for real_domain in bn_volume_by_domain_combined.keys():
        padded_index = real_domain_to_padded_index(real_domain)
        if 0 <= padded_index < TOTAL_DOMAINS:
            BN_ONLY_BEACH_DUNE_ON[padded_index] = True

    scenarios = [
        dict(
            scenario_name="natural",
            label="Natural",
            roadway_on=NATURAL_ROADWAY_ON,
            beach_dune_on=NATURAL_BEACH_DUNE_ON,
            use_historical_nourishment=False,
        ),
        dict(
            scenario_name="HIST_BN_only",
            label="Historical BN only",
            roadway_on=BN_ONLY_ROADWAY_ON,
            beach_dune_on=BN_ONLY_BEACH_DUNE_ON,
            use_historical_nourishment=True,
        ),
    ]

    rate_profiles = {}
    cascades_by_label = {}

    DUNE_INIT_YEARS = [2011]

    for dune_init_year in DUNE_INIT_YEARS:

        elevation_file_paths, dune_file_paths = (
            build_elevation_dune_file_paths(dune_init_year)
        )

        for scenario in scenarios:
            for Hs in WAVE_HEIGHTS_TO_TEST:

                run_name = (
                    f"PEA_{START_YEAR}_{END_YEAR}_"
                    f"{scenario['scenario_name']}_"
                    f"{dune_init_year}_duneline_Hs{Hs:.1f}"
                ).replace(".", "p")

                if APPLY_AUTOMATIC_BN_THRESHOLD:
                    run_name += "_thresholdBN_ON"
                else:
                    run_name += "_thresholdBN_OFF"

                plot_label = (
                    f"{scenario['label']} | {dune_init_year} dunes"
                )

                cascade = run_cascade_simulation(
                    nt=RUN_YEARS,
                    name=run_name,
                    storm_file=STORM_FILE,
                    alongshore_section_count=TOTAL_DOMAINS,
                    num_cores=NUM_CORES,

                    beach_width_threshold=BEACH_WIDTH_THRESHOLD,

                    rmin=RMIN,
                    rmax=RMAX,

                    elevation_file=elevation_file_paths,
                    dune_file=dune_file_paths,

                    dune_design_elevation=DUNE_DESIGN_ELEVATION,
                    dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,

                    road_ele=ROAD_ELEVATION,
                    road_width=ROAD_WIDTH,
                    road_setback=road_setbacks_full,

                    overwash_filter=OVERWASH_FILTER,
                    overwash_to_dune=OVERWASH_TO_DUNE,

                    nourishment_volume=NOURISHMENT_VOLUME,
                    background_erosion=BACKGROUND_EROSION_RATES,

                    rebuild_dune_threshold=REBUILD_DUNE_THRESHOLD,

                    roadway_management_on=scenario["roadway_on"],
                    beach_dune_manager_on=scenario["beach_dune_on"],

                    sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
                    sea_level_constant=SEA_LEVEL_CONSTANT,

                    enable_shoreline_offset=True,
                    shoreline_offset=dune_offset_dam,

                    wave_height=Hs,
                    wave_period=FIXED_WAVE_PERIOD,
                    wave_asymmetry=FIXED_WAVE_ASYMMETRY,
                    wave_angle_high_fraction=FIXED_WAVE_ANGLE_HIGH_FRACTION,

                    historical_nourishment_on_by_year=(
                        HIST_NOURISH_ON
                        if scenario["use_historical_nourishment"]
                        else None
                    ),

                    historical_nourishment_volume_by_year=(
                        HIST_NOURISH_VOLUME
                        if scenario["use_historical_nourishment"]
                        else None
                    ),
                )

                shoreline_m = build_shoreline_matrix(
                    cascade,
                    to_meters=TO_METERS,
                )

                total_change = shoreline_m[-1] - shoreline_m[0]

                change_rate = (
                        total_change / float(time_span_years)
                )

                if FLIP_SIGN_MODEL:
                    change_rate *= -1

                rate_profiles[(plot_label, Hs)] = change_rate

    # =============================================================================
    # YEARLY RELATIVE SHORELINE CHANGE + HISTORICAL BN PLOTS
    # =============================================================================

    threshold_tag = "thresholdBN_ON" if APPLY_AUTOMATIC_BN_THRESHOLD else "thresholdBN_OFF"

    yearly_prefix = (
        f"PEA_{START_YEAR}_{END_YEAR}_yearly_relative_shoreline_Natural_vs_HIST_BN_only_{threshold_tag}"
    )

    # plot_yearly_relative_shoreline_and_bn(
    #     cascades_by_label=cascades_by_label,
    #     hist_nourish_volume_by_year=HIST_NOURISH_VOLUME,
    #     output_dir=OUTPUT_BASE_DIR,
    #     filename_prefix=yearly_prefix,
    #     flip_sign=True,
    #     make_gif=MAKE_YEARLY_BN_GIF,
    #     gif_duration_seconds=GIF_DURATION_SECONDS,
    # )

    # =============================================================================
    # YEARLY ABSOLUTE SHORELINE POSITION + HISTORICAL BN PLOTS
    # =============================================================================

    absolute_yearly_prefix = (
        f"PEA_{START_YEAR}_{END_YEAR}_yearly_absolute_shoreline_Natural_vs_HIST_BN_only_{threshold_tag}"
    )

    # plot_yearly_absolute_shoreline_and_bn(
    #     cascades_by_label=cascades_by_label,
    #     hist_nourish_volume_by_year=HIST_NOURISH_VOLUME,
    #     output_dir=OUTPUT_BASE_DIR,
    #     filename_prefix=absolute_yearly_prefix,
    #     make_gif=MAKE_YEARLY_BN_GIF,
    #     gif_duration_seconds=GIF_DURATION_SECONDS,
    # )

    # =============================================================================
    # FINAL SHORELINE RATE PLOT
    # =============================================================================
    # =============================================================================
    # FINAL SHORELINE RATE PLOT
    # =============================================================================

    # Plot ONLY real Pea Island domains 80–120
    real_domain_numbers = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    fig, ax = plt.subplots(figsize=(14, 6), constrained_layout=True)

    PLOT_DUNE_YEARS = [2011]
    PLOT_SCENARIOS = ["natural", "Historical BN only"]
    # "Historical BN only"
    PLOT_HS_VALUES = [1, 2]

    for (plot_label, Hs), rate in rate_profiles.items():

        scenario_name, dune_text = plot_label.split(" | ")
        dune_year = int(dune_text.replace(" dunes", ""))

        if dune_year not in PLOT_DUNE_YEARS:
            continue

        if scenario_name not in PLOT_SCENARIOS:
            continue

        if Hs not in PLOT_HS_VALUES:
            continue

        rate_real = np.asarray(rate)[START_REAL_INDEX:END_REAL_INDEX]

        ax.plot(
            real_domain_numbers,
            rate_real,
            linewidth=2.2,
            label=f"{plot_label}, Hs={Hs} m",
        )
    # Plot CoastSat only on real domain IDs
    for s in coastsat_series:
        # CoastSat x is currently mapped to padded CASCADE index.
        # Convert padded index back to real Pea Island domain ID.
        coastsat_real_domains = FIRST_FILE_NUMBER + (np.asarray(s["x"]) - START_REAL_INDEX)

        keep_real = (
                (coastsat_real_domains >= FIRST_FILE_NUMBER)
                & (coastsat_real_domains <= LAST_FILE_NUMBER)
        )

        ax.plot(
            coastsat_real_domains[keep_real],
            np.asarray(s["rate"])[keep_real],
            linestyle="--",
            marker="o",
            linewidth=2,
            markersize=5,
            label=s["label"],
        )

    ax.axhline(0.0, linestyle="--", linewidth=2, color="0.25", alpha=0.75)

    # x-axis ticks are now real domains, not padded CASCADE indices
    xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, DOMAIN_TICK_STEP)

    ax.set_xticks(xticks)
    ax.set_xticklabels(
        [str(i) for i in xticks],
        rotation=45,
        ha="right",
        fontsize=9,
    )

    ax.set_xlim(FIRST_FILE_NUMBER, LAST_FILE_NUMBER)

    ax.set_xlabel(
        "Real Pea Island domain ID",
        fontsize=12,
        fontweight="bold",
    )

    ax.set_ylabel(
        "Shoreline change rate (m/yr)",
        fontsize=12,
        fontweight="bold",
    )

    ax.set_title(
        f"Pea Island shoreline change rates | "
        f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | "
        f"{START_YEAR}–{END_YEAR}",
        fontsize=14,
        fontweight="bold",
        pad=14,
    )

    ax.grid(True, which="major", linestyle=":", linewidth=0.8, alpha=0.35)
    ax.minorticks_on()
    ax.grid(True, which="minor", linestyle=":", linewidth=0.5, alpha=0.18)

    # -------------------------------------------------------------------------
    # Site annotations using REAL domain IDs
    # -------------------------------------------------------------------------
    ymin, ymax = ax.get_ylim()
    yrange = ymax - ymin

    site_line_kwargs = dict(
        color="0.55",
        linestyle="--",
        linewidth=2.0,
        alpha=0.60,
        zorder=0,
    )

    x_rodanthe = 80.5
    x_visitor = 110
    x_groin = 120

    ax.axvline(x_rodanthe, **site_line_kwargs)
    ax.axvline(x_visitor, **site_line_kwargs)
    ax.axvline(x_groin, **site_line_kwargs)

    ax.text(
        x_rodanthe,
        ymax - 0.035 * yrange,
        "Rodanthe",
        fontsize=8.5,
        color="0.35",
        ha="center",
        va="top",
    )

    ax.text(
        x_visitor,
        ymax - 0.035 * yrange,
        "Pea Island\nVisitor Center",
        fontsize=8.5,
        color="0.35",
        ha="center",
        va="top",
    )

    ax.text(
        x_groin,
        ymax - 0.035 * yrange,
        "Oregon Inlet\nGroin",
        fontsize=8.5,
        color="0.35",
        ha="center",
        va="top",
    )

    # -------------------------------------------------------------------------
    # y-axis limits based only on real-domain values
    # -------------------------------------------------------------------------
    all_y_values = []

    for rate in rate_profiles.values():
        rate_real = np.asarray(rate)[START_REAL_INDEX:END_REAL_INDEX]
        finite_rate = rate_real[np.isfinite(rate_real)]
        all_y_values.extend(finite_rate)

    for s in coastsat_series:
        coastsat_real_domains = FIRST_FILE_NUMBER + (np.asarray(s["x"]) - START_REAL_INDEX)
        obs_rate = np.asarray(s["rate"])

        keep_real = (
                (coastsat_real_domains >= FIRST_FILE_NUMBER)
                & (coastsat_real_domains <= LAST_FILE_NUMBER)
                & np.isfinite(obs_rate)
        )

        all_y_values.extend(obs_rate[keep_real])

    if len(all_y_values) > 0:
        y_min = np.nanmin(all_y_values)
        y_max = np.nanmax(all_y_values)

        if y_max > y_min:
            y_pad = 0.22 * (y_max - y_min)
        else:
            y_pad = 1.0

        ax.set_ylim(y_min - y_pad, y_max + y_pad)

    ax.legend()

    out_png = os.path.join(
        OUTPUT_BASE_DIR,
        f"PEA_{START_YEAR}_{END_YEAR}_shoreline_change_rates_real_domains_only_hannahstorm.png",
    )

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()

    print(f"\nSaved real-domain-only shoreline rate plot:")
    print(out_png)


if __name__ == "__main__":
    main()