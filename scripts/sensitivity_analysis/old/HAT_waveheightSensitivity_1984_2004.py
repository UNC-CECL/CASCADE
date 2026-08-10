#!/usr/bin/env python3
"""
CASCADE Hatteras Island: Wave Climate Sensitivity Analysis
==========================================================
Tests all four wave climate parameters across physically realistic ranges
for Hatteras Island, NC.  Hit run — the script iterates automatically.

Each parameter is swept independently (one-at-a-time) while all others
hold their calibration baseline.  For each parameter value the script:
  - Runs a full CASCADE simulation
  - Saves the NPZ comparison and a run metadata file
  - Generates an annotated publication figure with CoastSat overlay

After all sweeps complete, a 2x2 session overview figure is written.

Wave parameter ranges (Hatteras-specific justification)
--------------------------------------------------------
Wave height (Hs)
  1.0, 1.2, 1.5, 1.8, 2.0 m
  USGS hindcast mean ~1.2-1.3 m; morphologically effective Heff ~1.4-1.5 m
  under Rayleigh (1.19x mean).  Upper end explores storm-biased forcing.

Wave period (Tp)
  7, 8, 9, 10 s
  WIS station ST63228 modal Tp range for Hatteras; 8 s is the calibrated default.

Wave asymmetry
  0.6, 0.65, 0.7, 0.75
  US East Coast barrier island range.  0.7 is the calibrated default.

Wave angle high fraction
  0.1, 0.15, 0.2, 0.25, 0.3
  0.2 is the recommended Hatteras default; sweep tests sensitivity to +/-0.1.

Output folder structure
-----------------------
comparison/sensitivity_analysis/
  HAT_1984_2004_waveSensitivity_YYYYMMDD_HHMMSS/
    wave_height/
      HAT_1984_2004_wvSens_wave_height_1p0/    <- NPZ + per-run metadata
      HAT_1984_2004_wvSens_wave_height_1p2/
      ...
      wave_height_sensitivity.png
      wave_height_results.csv
    wave_period/ ...
    wave_asymmetry/ ...
    wave_angle_high_fraction/ ...
    session_overview.png
    session_metadata.txt

Author: Hannah Henry (UNC Chapel Hill)
"""

import os
import sys
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.transforms import blended_transform_factory
from statsmodels.nonparametric.smoothers_lowess import lowess
from cascade.cascade import Cascade

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION
# =============================================================================

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 1
LAST_FILE_NUMBER  = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1   # 90

TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS                                          # 15
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS                        # 105

FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN  = 90
START_ROAD_INDEX  = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS   # 23
END_ROAD_INDEX    = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS + 1  # 105

print("=" * 80)
print("CASCADE WAVE CLIMATE SENSITIVITY ANALYSIS - HATTERAS ISLAND")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS} (GIS IDs {FIRST_FILE_NUMBER}..{LAST_FILE_NUMBER})")
print(f"Buffer Domains: {NUM_BUFFER_DOMAINS} each side  |  Total: {TOTAL_DOMAINS}")
print(f"Real index span (padded): [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR   = r"/"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "comparison", "sensitivity_analysis")
COASTSAT_BASE_DIR  = os.path.join(
    PROJECT_BASE_DIR, "scripts", "input_prep", "CoastSat"
)

# Dune raw_offset file — 1984 initialisation
DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE,
    "2-brie-offset",
    "hindcast_1984",
    f"Island_Dune_Offsets_1984_PADDED_{TOTAL_DOMAINS}.csv",
)

# Storm file — 1984_2004 period
STORM_FILE = os.path.join(
    HATTERAS_DATA_BASE,
    "storms", "hindcast_storms", "base_storms",
    "storms_1984_2004_base.npy",
)

# Road setbacks — 1984 initialisation
ROAD_SETBACK_FILE = os.path.join(
    HATTERAS_DATA_BASE, "road_offset", "raw_offset", "1984", "RoadSetback_1984.csv"
)

# Topo/dune groin_init files
TOPO_DUNE_INIT_YEAR = "2009"
TOPO_DUNE_SUBFOLDER = r"2009\2009_v2"

PARAMETER_FILE = "Hatteras-CASCADE-parameters.yaml"

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# =============================================================================
# SECTION 3: COASTSAT OBSERVED REFERENCE
# =============================================================================
# Both periods are loaded; the active period (matching START_YEAR) is plotted
# as a solid reference line.  The secondary period is shown dashed/lighter.

COASTSAT_DATASETS = [
    dict(
        label        = "CoastSat LRR (1984_2004)",
        period_start = 1984,
        csv          = os.path.join(COASTSAT_BASE_DIR, "1984_2004", "domain_lrr_summary.csv"),
        domain_col   = "domain_number",
        rate_col     = "mean_lrr",
    ),
    dict(
        label        = "CoastSat LRR (2004-2024)",
        period_start = 2004,
        csv          = os.path.join(
            COASTSAT_BASE_DIR, "2004_2024_specific_dates", "domain_lrr_summary.csv"
        ),
        domain_col   = "domain_number",
        rate_col     = "mean_lrr",
    ),
]

# =============================================================================
# SECTION 4: GEOGRAPHIC ANNOTATION CONSTANTS
# =============================================================================

ANN_TOWN_SPANS = {
    "Buxton":      (7,  8),
    "Avon":        (21, 31),
    "Tri-Village": (68, 83),
}
ANN_VILLAGE_LINES = {"Salvo": 69, "Waves": 74, "Rodanthe": 80}
ANN_PIERS         = {"Avon Pier": 26, "Rodanthe Pier": 79}
ANN_GROINS        = {"Buxton Groin": 6}
ANN_WIMBLE_SHOALS = (60, 74)

ANN_C_TOWN_SPAN    = "#90AFC5"
ANN_C_WIMBLE       = "#E0A800"
ANN_C_VILLAGE_LINE = "0.40"
ANN_C_PIER         = "#1565C0"
ANN_C_GROIN        = "#B71C1C"

# CoastSat line colours in sensitivity plots
CS_ACTIVE_COLOR = "#333333"   # dark gray — active period (solid)
CS_REF_COLOR    = "#AAAAAA"   # light gray — reference period (dashed)

# =============================================================================
# LOESS SMOOTHING REFERENCE CURVES
# =============================================================================
# Window sizes match coastsat_smoothed_methods_comparison_line.py exactly:
#   3.5 km = 7 domains  |  5.0 km = 10 domains
# These are computed from the active CoastSat period and overlaid on every
# sensitivity figure so you can judge which Hs best reproduces each target.
DOMAIN_SPACING_M = 500   # metres per CASCADE domain
LOESS_WINDOWS    = {7: 3.5, 10: 5.0}   # {n_domains: window_km}
LOESS_COLORS     = {7: "#0072B2", 10: "#CC79A7"}   # blue=7-dom, pink=10-dom

# =============================================================================
# SECTION 5: FIXED SIMULATION PARAMETERS
# =============================================================================
# These match HAT_hindcast_1984_2024_old version.py exactly and are held constant
# while wave parameters are varied one at a time.

START_YEAR = 1984
END_YEAR   = 2004
RUN_YEARS  = END_YEAR - START_YEAR   # 20

TO_METERS = True

# Sea level rise
SEA_LEVEL_RISE_RATE = 0.004   # m/yr — from duck_rslr_analysis.py
SEA_LEVEL_CONSTANT  = True

# Vertical datum
BERM_ELEVATION = 1.7    # m NAVD88
MHW_ELEVATION  = 0.36   # m NAVD88

NUM_CORES = 1   # >1 causes crashes

# Morphodynamics
DUNE_REBUILD_HEIGHT    = 3.0
REBUILD_ELEV_THRESHOLD = 0.01   # dam
OVERWASH_TO_DUNE_VAL   = 9
SANDBAG_ELEV           = 0

# Management (all off — unmanaged baseline for calibration)
ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_NOURISHMENT        = False
ENABLE_SANDBAG_PLACEMENT  = False
NOURISHMENT_VOLUME        = 0

# Sign convention
FLIP_SIGN_MODEL = True

# Background erosion — all zero for wave sensitivity sweeps
# (BE rates are calibrated separately after wave climate is fixed)
BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS

# Overwash filter — community values matching hindcast run script
OVERWASH_FILTER_DEFAULT         = 0.0
OVERWASH_FILTER_BUXTON          = 0.4   # GIS  7-8
OVERWASH_FILTER_AVON            = 0.4   # GIS 21-31
OVERWASH_FILTER_SALVO_WAVES_ROD = 0.4   # GIS 68-83

DOMAIN_TICK_STEP = 5

# =============================================================================
# SECTION 6: WAVE SENSITIVITY CONFIGURATION  <- primary edit point
# =============================================================================
#
# BASELINE VALUES — held fixed when a different parameter is being swept.
# These match the current calibrated defaults in HAT_hindcast_1984_2024_old version.py.
#
BASELINE_WAVE_HEIGHT              = 2.0   # m
BASELINE_WAVE_PERIOD              = 8     # s
BASELINE_WAVE_ASYMMETRY           = 0.7
BASELINE_WAVE_ANGLE_HIGH_FRACTION = 0.2

#
# SENSITIVITY RANGES — physically realistic for Hatteras Island.
#
# To disable a parameter sweep, set  enabled: False.
# To add or remove values, edit the  values  list.
#
SENSITIVITY_ANALYSES = {
    "wave_height": {
        "label":   "Wave Height",
        "units":   "m",
        "enabled": True,
        # Broad sweep from well below physical mean up to storm-dominated forcing.
        # USGS hindcast mean ~1.2-1.3 m; morphologically effective ~1.4-1.5 m.
        # Values below 1.0 m are physically unrealistic for Hatteras but are
        # included because BRIE diffusivity scales with Hs — testing the full
        # range identifies where the model transitions between erosion regimes.
        "values":  [0.1, 0.25, 0.5, 0.8, 1.0, 1.2, 1.5, 1.8, 2.0, 2.5, 3.0],
    },
    "wave_period": {
        "label":   "Wave Period",
        "units":   "s",
        "enabled": False,   # disabled — wave height sweep only
        "values":  [7, 8, 9, 10],
    },
    "wave_asymmetry": {
        "label":   "Wave Asymmetry",
        "units":   "",
        "enabled": False,   # disabled — wave height sweep only
        "values":  [0.6, 0.65, 0.7, 0.75],
    },
    "wave_angle_high_fraction": {
        "label":   "Wave Angle High Fraction",
        "units":   "",
        "enabled": False,   # disabled — wave height sweep only
        "values":  [0.1, 0.15, 0.2, 0.25, 0.3],
    },
}

# Whether to save CASCADE NPZ comparison for each run.
# True  = full reproducibility but ~190 MB per run (~3.5 GB for full sweep).
# False = saves disk space; rate CSVs and metadata are always written.
SAVE_CASCADE_NPZ = True

# =============================================================================
# PRINT CONFIGURATION SUMMARY
# =============================================================================

enabled_params = [k for k, v in SENSITIVITY_ANALYSES.items() if v["enabled"]]
total_runs     = sum(len(v["values"]) for v in SENSITIVITY_ANALYSES.values() if v["enabled"])

print("Wave Sensitivity Configuration:")
print(f"  Period:           {START_YEAR}-{END_YEAR}  ({RUN_YEARS} years)")
print(f"  Parameters:       {len(enabled_params)} enabled")
print(f"  Total runs:       {total_runs}")
print(f"  SLR:              {SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr")
print(f"  Save NPZ:         {SAVE_CASCADE_NPZ}")
print()
print("  Baseline (held fixed when not varied):")
print(f"    Hs              = {BASELINE_WAVE_HEIGHT} m")
print(f"    Tp              = {BASELINE_WAVE_PERIOD} s")
print(f"    Asymmetry       = {BASELINE_WAVE_ASYMMETRY}")
print(f"    Angle high frac = {BASELINE_WAVE_ANGLE_HIGH_FRACTION}")
print()
for pname, pcfg in SENSITIVITY_ANALYSES.items():
    flag = "ENABLED " if pcfg["enabled"] else "disabled"
    print(f"  [{flag}] {pcfg['label']:30s}  {pcfg['values']}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 7: LOAD FIXED INPUT DATA
# =============================================================================

print("Loading input data...")

try:
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=",")
    dune_offset_dam = dune_offset_all / 10.0   # m -> dam (single-column file)
    print(f"  Dune offsets:  {dune_offset_dam.size} domains (dam)")

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=",")
    print(f"  Road setbacks: {road_setbacks_raw.size} values")

except FileNotFoundError as e:
    print(f"CRITICAL ERROR: Missing data file: {e.filename}")
    sys.exit(1)
except Exception as e:
    print(f"CRITICAL ERROR loading data: {e}")
    sys.exit(1)

# Road setback padded array
road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX:START_ROAD_INDEX + num_road_values] = \
    road_setbacks_raw[:num_road_values]

ROADWAY_MANAGEMENT_ON     = [ENABLE_ROADWAY_MANAGEMENT] * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON     = [ENABLE_SANDBAG_PLACEMENT]  * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [ENABLE_NOURISHMENT]        * TOTAL_DOMAINS

print(f"  Data loaded OK.\n")

# =============================================================================
# SECTION 8: ELEVATION + DUNE FILE LISTS
# =============================================================================

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS      = []

for _ in range(START_REAL_INDEX):
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_dune.npy"))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy"))

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    DUNE_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, "dunes", TOPO_DUNE_SUBFOLDER,
                     f"domain_{file_num}_dune_{TOPO_DUNE_INIT_YEAR}.npy"))
    ELEVATION_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, "topography", TOPO_DUNE_SUBFOLDER,
                     f"domain_{file_num}_topography_{TOPO_DUNE_INIT_YEAR}.npy"))

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_dune.npy"))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, "buffer", "sample_1_topography.npy"))

print(f"  Elevation file paths: {len(ELEVATION_FILE_PATHS)}")
print(f"  Dune file paths:      {len(DUNE_FILE_PATHS)}")
print("=" * 80 + "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def _gis_to_pad(gis_id):
    """Convert a 1-based GIS domain ID to CASCADE padded index."""
    return START_REAL_INDEX + (gis_id - FIRST_FILE_NUMBER)


def get_x_s_TS(b3d):
    """Extract shoreline time series from a Barrier3D object."""
    for attr in ("x_s_TS", "_x_s_TS"):
        if hasattr(b3d, attr):
            return np.asarray(getattr(b3d, attr), dtype=float)
    raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    """Build [time x domain] shoreline matrix from CASCADE object."""
    b3d_list = cascade.barrier3d
    nt = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, len(b3d_list)), dtype=float)
    for j, b3d in enumerate(b3d_list):
        shoreline[:, j] = get_x_s_TS(b3d)
    if to_meters:
        shoreline *= 10.0   # dam -> m
    return shoreline


def load_coastsat(ds):
    """
    Load one CoastSat domain_lrr_summary CSV and map to CASCADE padded indices.
    Returns (padded_x, rates_m_per_yr) or (None, None) if file is missing.
    """
    path = ds["csv"]
    if not os.path.exists(path):
        print(f"    CoastSat file not found (skipped): {os.path.basename(path)}")
        return None, None

    df = pd.read_csv(path)
    # Strip filename prefix that sometimes gets prepended to header
    df.columns = [c.split("csv")[-1] if "csv" in c else c for c in df.columns]

    dom_col  = ds["domain_col"]
    rate_col = ds["rate_col"]
    if dom_col not in df.columns or rate_col not in df.columns:
        return None, None

    raw_dom  = pd.to_numeric(df[dom_col],  errors="coerce")
    raw_rate = pd.to_numeric(df[rate_col], errors="coerce")
    m        = raw_dom.notna() & raw_rate.notna()
    raw_dom  = raw_dom[m].astype(int).values
    raw_rate = raw_rate[m].values

    # GIS domain IDs (1-based) -> padded indices
    obs_x    = START_REAL_INDEX + (raw_dom - 1)
    keep     = (obs_x >= START_REAL_INDEX) & (obs_x < END_REAL_INDEX)
    obs_x    = obs_x[keep].astype(int)
    obs_rate = raw_rate[keep]

    order    = np.argsort(obs_x)
    return obs_x[order], obs_rate[order]


# =============================================================================
# LOESS SMOOTHING HELPERS
# =============================================================================
# Logic mirrors coastsat_smoothed_methods_comparison_line.py so window sizes
# are directly comparable.

def _loess_frac(window_km, n_points):
    """Convert a physical window width (km) to a LOESS frac for n_points."""
    k = (window_km * 1000.0) / DOMAIN_SPACING_M
    return float(np.clip(k / n_points, 0.02, 1.0))


def _apply_loess(x, y, frac):
    """Apply LOESS smoother. Returns smoothed y at same x positions."""
    valid = ~np.isnan(y)
    if valid.sum() < 5:
        return y.copy()
    smoothed = np.full_like(y, np.nan, dtype=float)
    result   = lowess(y[valid], x[valid], frac=frac, return_sorted=True)
    smoothed[valid] = np.interp(x[valid], result[:, 0], result[:, 1])
    return smoothed


def compute_loess_curves(cs_series):
    """
    Compute 7-domain and 10-domain LOESS smoothed reference curves from the
    active CoastSat period.  Window sizes match LOESS_WINDOWS (Section 4).

    Returns
    -------
    dict : {n_domains -> {'x': gis_id_array, 'rate': smoothed_rate_array}}
    """
    loess_curves = {}
    for cs in cs_series:
        if not cs["active"]:
            continue
        gis_x = (cs["x"] - START_REAL_INDEX + FIRST_FILE_NUMBER).astype(float)
        for n_dom, km in LOESS_WINDOWS.items():
            frac     = _loess_frac(km, len(gis_x))
            smoothed = _apply_loess(gis_x, cs["rate"], frac)
            loess_curves[n_dom] = {"x": gis_x, "rate": smoothed}
            print(f"  LOESS {n_dom}-domain ({km} km, frac={frac:.3f}): "
                  f"range {smoothed[~np.isnan(smoothed)].min():+.2f}–"
                  f"{smoothed[~np.isnan(smoothed)].max():+.2f} m/yr")
    return loess_curves


# =============================================================================
# GEOGRAPHIC ANNOTATION HELPERS
# =============================================================================

def add_geographic_annotations(ax):
    """Add community spans, village lines, piers, groin, and Wimble Shoals."""
    trans = blended_transform_factory(ax.transData, ax.transAxes)

    # 1. Wimble Shoals
    wlo, whi = ANN_WIMBLE_SHOALS
    ax.axvspan(wlo - 0.5, whi + 0.5,
               color=ANN_C_WIMBLE, alpha=0.10, zorder=0,
               hatch="///", edgecolor=ANN_C_WIMBLE, linewidth=0)
    ax.text((wlo + whi) / 2.0, 0.04,
            "Wimble Shoals\ninfluence", transform=trans,
            ha="center", va="bottom", fontsize=7, color="#7A5800", style="italic",
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.80))

    # 2. Community spans
    for span_label, (d_lo, d_hi) in ANN_TOWN_SPANS.items():
        ax.axvspan(d_lo - 0.5, d_hi + 0.5,
                   color=ANN_C_TOWN_SPAN, alpha=0.14, zorder=0)
        ax.text((d_lo + d_hi) / 2.0, 0.90,
                span_label, transform=trans,
                ha="center", va="top", fontsize=8, color="0.25", fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.85))

    # 3. Village center lines
    for vname, dom in ANN_VILLAGE_LINES.items():
        ax.axvline(dom, color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--",
                   alpha=0.65, zorder=1)
        ax.text(dom, 0.84, vname, transform=trans,
                ha="center", va="top", fontsize=7.5, color="0.30",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 4. Piers
    for pname, dom in ANN_PIERS.items():
        ax.axvline(dom, color=ANN_C_PIER, lw=1.0, ls="-.", alpha=0.80, zorder=2)
        ax.text(dom, 0.76, pname, transform=trans,
                ha="center", va="top", fontsize=7, color=ANN_C_PIER, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))

    # 5. Groins
    for gname, dom in ANN_GROINS.items():
        ax.axvline(dom, color=ANN_C_GROIN, lw=1.1, ls=":", alpha=0.85, zorder=2)
        ax.text(dom, 0.76, gname, transform=trans,
                ha="center", va="top", fontsize=7, color=ANN_C_GROIN, rotation=90,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.80))


def annotation_legend_handles():
    """Proxy legend artists for the geographic annotation layers."""
    return [
        Patch(fc=ANN_C_TOWN_SPAN, alpha=0.30, label="Community"),
        Patch(fc=ANN_C_WIMBLE, alpha=0.25, hatch="///",
              edgecolor=ANN_C_WIMBLE, linewidth=0, label="Wimble Shoals"),
        Line2D([0], [0], color=ANN_C_VILLAGE_LINE, lw=0.9, ls="--",
               label="Village center"),
        Line2D([0], [0], color=ANN_C_PIER, lw=1.0, ls="-.", label="Pier"),
        Line2D([0], [0], color=ANN_C_GROIN, lw=1.1, ls=":", label="Groin"),
    ]


# =============================================================================
# RUN METADATA
# =============================================================================

def write_run_metadata(run_dir, run_name, wh, wp, wa, wahf, param_name, param_value):
    """
    Write a complete record of all parameters to {run_name}_run_metadata.txt.
    Format matches HAT_hindcast_1984_2024_old version.py for consistency.
    """
    lines = [
        "# CASCADE Run Metadata — Wave Sensitivity Analysis",
        "# Generated automatically — do not edit by hand.",
        "",
        "[Run Identity]",
        f"run_name              = {run_name}",
        f"timestamp             = {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"script                = HAT_waveSensitivity_1984_2004.py",
        f"sensitivity_param     = {param_name}",
        f"sensitivity_value     = {param_value}",
        "",
        "[Period]",
        f"start_year            = {START_YEAR}",
        f"end_year              = {END_YEAR}",
        f"run_years             = {RUN_YEARS}",
        "",
        "[Wave Climate]",
        f"wave_height_m         = {wh}",
        f"wave_period_s         = {wp}",
        f"wave_asymmetry        = {wa}",
        f"wave_angle_high_frac  = {wahf}",
        "",
        "[Baseline (held fixed for this sweep)]",
        f"baseline_Hs           = {BASELINE_WAVE_HEIGHT} m",
        f"baseline_Tp           = {BASELINE_WAVE_PERIOD} s",
        f"baseline_asymmetry    = {BASELINE_WAVE_ASYMMETRY}",
        f"baseline_angle_high   = {BASELINE_WAVE_ANGLE_HIGH_FRACTION}",
        "",
        "[Sea Level Rise]",
        f"slr_rate_mm_per_yr    = {SEA_LEVEL_RISE_RATE * 1000:.2f}",
        f"slr_constant          = {SEA_LEVEL_CONSTANT}",
        "",
        "[Vertical Datum]",
        f"berm_elevation_m      = {BERM_ELEVATION}",
        f"mhw_elevation_m       = {MHW_ELEVATION}",
        "",
        "[Storm]",
        f"storm_file            = {os.path.basename(STORM_FILE)}",
        "",
        "[Topo / Dune Initialisation]",
        f"topo_dune_init_year   = {TOPO_DUNE_INIT_YEAR}",
        f"topo_dune_subfolder   = {TOPO_DUNE_SUBFOLDER}",
        "",
        "[Morphodynamics]",
        f"dune_rebuild_height_m = {DUNE_REBUILD_HEIGHT}",
        f"rebuild_elev_thresh   = {REBUILD_ELEV_THRESHOLD}  # dam",
        f"overwash_to_dune      = {OVERWASH_TO_DUNE_VAL}",
        f"flip_sign_model       = {FLIP_SIGN_MODEL}",
        f"background_erosion    = 0.0 (all domains — BE not applied in sensitivity sweeps)",
        "",
        "[Overwash Filter]",
        f"default               = {OVERWASH_FILTER_DEFAULT}",
        f"buxton                = {OVERWASH_FILTER_BUXTON}   (GIS 7-8)",
        f"avon                  = {OVERWASH_FILTER_AVON}   (GIS 21-31)",
        f"tri_village           = {OVERWASH_FILTER_SALVO_WAVES_ROD}   (GIS 68-83)",
        "",
        "[Management]",
        "roadway / nourishment / sandbag = False (unmanaged baseline)",
        "",
        "[Domain Configuration]",
        f"num_real_domains      = {NUM_REAL_DOMAINS}",
        f"num_buffer_domains    = {NUM_BUFFER_DOMAINS}",
        f"total_domains         = {TOTAL_DOMAINS}",
    ]

    out_path = os.path.join(run_dir, f"{run_name}_run_metadata.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out_path


# =============================================================================
# CASCADE RUNNER
# =============================================================================

def run_cascade_simulation(
    name, save_dir,
    wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction,
):
    """
    Run a single CASCADE simulation.  All fixed parameters are taken from
    module-level constants; only wave parameters vary between calls.
    """
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS
    DUNE_DESIGN_ELEVATION  = [DUNE_REBUILD_HEIGHT]    * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    ROAD_ELEVATION = 1.45
    ROAD_WIDTH     = 20.0

    # Per-domain overwash filter — matches hindcast run script
    ow_filter = np.full(TOTAL_DOMAINS, OVERWASH_FILTER_DEFAULT)
    ow_filter[_gis_to_pad(7)  : _gis_to_pad(8)  + 1] = OVERWASH_FILTER_BUXTON
    ow_filter[_gis_to_pad(21) : _gis_to_pad(31) + 1] = OVERWASH_FILTER_AVON
    ow_filter[_gis_to_pad(68) : _gis_to_pad(83) + 1] = OVERWASH_FILTER_SALVO_WAVES_ROD
    ow_filter = list(ow_filter)

    cascade = Cascade(
        HATTERAS_DATA_BASE,
        name,
        storm_file=STORM_FILE,
        elevation_file=ELEVATION_FILE_PATHS,
        dune_file=DUNE_FILE_PATHS,
        parameter_file=PARAMETER_FILE,

        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=wave_asymmetry,
        wave_angle_high_fraction=wave_angle_high_fraction,

        sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
        sea_level_rise_constant=SEA_LEVEL_CONSTANT,

        background_erosion=BACKGROUND_EROSION_RATES,
        alongshore_section_count=TOTAL_DOMAINS,
        time_step_count=RUN_YEARS,

        min_dune_growth_rate=RMIN,
        max_dune_growth_rate=RMAX,
        num_cores=NUM_CORES,

        roadway_management_module=ROADWAY_MANAGEMENT_ON,
        beach_nourishment_module=NOURISHMENT_MANAGEMENT_ON,
        sandbag_management_on=SANDBAG_MANAGEMENT_ON,
        alongshore_transport_module=True,
        community_economics_module=False,

        road_ele=ROAD_ELEVATION,
        road_width=ROAD_WIDTH,
        road_setback=road_setbacks_full,

        dune_design_elevation=DUNE_DESIGN_ELEVATION,
        dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,
        sandbag_elevation=SANDBAG_ELEV,

        overwash_filter=ow_filter,
        overwash_to_dune=OVERWASH_TO_DUNE_VAL,

        enable_shoreline_offset=True,
        shoreline_offset=dune_offset_dam,

        nourishment_volume=NOURISHMENT_VOLUME,
        nourishment_interval=None,

        berm_elevation=BERM_ELEVATION,
        MHW=MHW_ELEVATION,
    )

    for time_step in range(RUN_YEARS - 1):
        print(f"\r    Year {time_step + 1}/{RUN_YEARS}", end="", flush=True)
        cascade.update()
        if getattr(cascade, "b3d_break", False):
            print(f"\n    Model stopped at year {time_step + 1} (b3d_break)")
            break

    if SAVE_CASCADE_NPZ:
        os.makedirs(save_dir, exist_ok=True)
        cascade.save(save_dir)
        print(f"\n    Saved NPZ: {save_dir}")
    else:
        print()

    return cascade


# =============================================================================
# PARAMETER SWEEP FUNCTION
# =============================================================================

def run_sweep(param_name, param_cfg, param_output_dir, cs_series, loess_curves):
    """
    Run the full sweep for one wave parameter.

    Parameters
    ----------
    param_name       : str  — key in SENSITIVITY_ANALYSES
    param_cfg        : dict — {label, units, values, enabled}
    param_output_dir : str  — comparison folder for this parameter
    cs_series        : list — loaded CoastSat series dicts

    Returns
    -------
    dict : {param_value -> change_rate array (TOTAL_DOMAINS,)}
    """
    print(f"\n{'='*80}")
    print(f"SWEEP: {param_cfg['label']}")
    print(f"  Values: {param_cfg['values']}")
    print(f"  Output: {param_output_dir}")
    print(f"{'='*80}\n")

    os.makedirs(param_output_dir, exist_ok=True)
    rate_profiles = {}   # param_value -> change_rate (TOTAL_DOMAINS,)
    results_rows  = []

    for param_value in param_cfg["values"]:
        val_str  = str(param_value).replace(".", "p")
        run_name = f"HAT_{START_YEAR}_{END_YEAR}_wvSens_{param_name}_{val_str}"
        run_dir  = os.path.join(param_output_dir, run_name)

        print(f"  -- {param_cfg['label']} = {param_value}{param_cfg['units']}  "
              f"({param_cfg['values'].index(param_value) + 1}/{len(param_cfg['values'])})")

        # Build wave parameters for this run
        wh   = param_value if param_name == "wave_height"              else BASELINE_WAVE_HEIGHT
        wp   = param_value if param_name == "wave_period"              else BASELINE_WAVE_PERIOD
        wa   = param_value if param_name == "wave_asymmetry"           else BASELINE_WAVE_ASYMMETRY
        wahf = param_value if param_name == "wave_angle_high_fraction" else BASELINE_WAVE_ANGLE_HIGH_FRACTION

        try:
            cascade = run_cascade_simulation(
                name=run_name,
                save_dir=run_dir,
                wave_height=wh,
                wave_period=wp,
                wave_asymmetry=wa,
                wave_angle_high_fraction=wahf,
            )
        except Exception as e:
            print(f"\n    RUN FAILED: {e}")
            import traceback; traceback.print_exc()
            continue

        # Compute shoreline change rate
        shoreline_m  = build_shoreline_matrix(cascade, to_meters=TO_METERS)
        nt_actual, _ = shoreline_m.shape
        denom        = RUN_YEARS if RUN_YEARS > 0 else max(nt_actual - 1, 1)
        change_rate  = (shoreline_m[-1, :] - shoreline_m[0, :]) / float(denom)
        if FLIP_SIGN_MODEL:
            change_rate *= -1.0

        rate_profiles[param_value] = change_rate

        # ── RMSE against each LOESS reference curve ──────────────────────────
        rmse_by_window = {}
        for n_dom, lc in loess_curves.items():
            # lc["x"] is in GIS domain IDs (1–90); map to padded indices for indexing
            obs_pad  = (lc["x"] - FIRST_FILE_NUMBER + START_REAL_INDEX).astype(int)
            keep     = (obs_pad >= START_REAL_INDEX) & (obs_pad < END_REAL_INDEX)
            obs_pad  = obs_pad[keep]
            obs_rate = lc["rate"][keep]
            mod_rate = change_rate[obs_pad]
            valid    = ~np.isnan(obs_rate)
            rmse     = float(np.sqrt(np.mean((mod_rate[valid] - obs_rate[valid]) ** 2)))
            rmse_by_window[n_dom] = rmse
            print(f"    RMSE vs {n_dom}-domain LOESS: {rmse:.3f} m/yr")

        # Write per-run metadata and rate CSV
        os.makedirs(run_dir, exist_ok=True)
        meta_path = write_run_metadata(run_dir, run_name, wh, wp, wa, wahf,
                                       param_name, param_value)
        print(f"    Metadata: {os.path.basename(meta_path)}")

        csv_path = os.path.join(run_dir, f"{run_name}_shoreline_change_rate.csv")
        pd.DataFrame({
            "cascade_padded_index": np.arange(TOTAL_DOMAINS),
            "gis_domain_id": (
                [None] * START_REAL_INDEX
                + list(range(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1))
                + [None] * (TOTAL_DOMAINS - END_REAL_INDEX)
            ),
            "model_rate_m_per_yr": change_rate,
        }).to_csv(csv_path, index=False)
        print(f"    Rate CSV: {os.path.basename(csv_path)}")

        real_rates = change_rate[START_REAL_INDEX:END_REAL_INDEX]
        for j, dom_idx in enumerate(range(START_REAL_INDEX, END_REAL_INDEX)):
            results_rows.append({
                "param_name":            param_name,
                "param_value":           param_value,
                "wave_height":           wh,
                "wave_period":           wp,
                "wave_asymmetry":        wa,
                "wave_angle_high_frac":  wahf,
                "gis_domain_id":         j + 1,
                "cascade_padded_index":  dom_idx,
                "model_rate_m_per_yr":   real_rates[j],
                **{f"rmse_loess_{n}dom": rmse_by_window.get(n, np.nan)
                   for n in LOESS_WINDOWS},
            })

    if not rate_profiles:
        print(f"  No successful runs for {param_cfg['label']} — skipping plot.")
        return rate_profiles

    # ── Save results CSV ──────────────────────────────────────────────────────
    results_df = pd.DataFrame(results_rows)
    results_csv = os.path.join(param_output_dir, f"{param_name}_results.csv")
    results_df.to_csv(results_csv, index=False)
    print(f"  Results CSV: {results_csv}")

    # ── RMSE ranked summary ───────────────────────────────────────────────────
    # One RMSE per (param_value, loess_window) — collapse the per-domain rows
    rmse_summary = (
        results_df
        .groupby("param_value")
        [[f"rmse_loess_{n}dom" for n in LOESS_WINDOWS]]
        .first()   # RMSE is the same for all domain rows within a run
        .reset_index()
    )
    print(f"\n  {'─'*60}")
    print(f"  RMSE SUMMARY — {param_cfg['label']} sweep")
    print(f"  {'─'*60}")
    header = f"  {'Value':>8s}"
    for n in LOESS_WINDOWS:
        header += f"   RMSE vs {n}-dom LOESS"
    print(header)
    for _, row in rmse_summary.iterrows():
        line = f"  {row['param_value']:>8}"
        for n in LOESS_WINDOWS:
            line += f"   {row[f'rmse_loess_{n}dom']:>21.3f} m/yr"
        print(line)
    for n in LOESS_WINDOWS:
        col     = f"rmse_loess_{n}dom"
        best_v  = rmse_summary.loc[rmse_summary[col].idxmin(), "param_value"]
        best_r  = rmse_summary[col].min()
        print(f"  → Best fit ({n}-domain LOESS): "
              f"{param_cfg['label']} = {best_v}{param_cfg['units']}  "
              f"(RMSE={best_r:.3f} m/yr)")
    print(f"  {'─'*60}\n")

    # ── Generate annotated sensitivity figure ─────────────────────────────────
    _plot_parameter_sweep(
        param_name    = param_name,
        param_cfg     = param_cfg,
        rate_profiles = rate_profiles,
        cs_series     = cs_series,
        loess_curves  = loess_curves,
        output_dir    = param_output_dir,
    )

    return rate_profiles


# =============================================================================
# PLOTTING
# =============================================================================

def _plot_parameter_sweep(param_name, param_cfg, rate_profiles, cs_series,
                          loess_curves, output_dir):
    """
    Annotated figure for one parameter sweep:
      - CoastSat active period (gray reference)
      - Model lines coloured by parameter value (viridis)
      - Full geographic annotation layer
    """
    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)
    n_model = len(rate_profiles)
    cmap    = plt.get_cmap("viridis", n_model)
    values  = param_cfg["values"]

    fig, ax = plt.subplots(figsize=(14, 5.8), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Geographic annotations
    add_geographic_annotations(ax)

    # CoastSat raw reference lines (faint — contextual background)
    data_handles = []
    for cs in cs_series:
        gis_x = cs["x"] - START_REAL_INDEX + FIRST_FILE_NUMBER
        if cs["active"]:
            ax.plot(gis_x, cs["rate"], color=CS_ACTIVE_COLOR,
                    lw=1.0, ls="-", alpha=0.30, zorder=3)   # faint — LOESS takes focus
        else:
            ax.plot(gis_x, cs["rate"], color=CS_REF_COLOR,
                    lw=0.9, ls="--", alpha=0.25, zorder=2)

    # LOESS smoothed reference curves (7-domain and 10-domain)
    for n_dom, lc in loess_curves.items():
        col = LOESS_COLORS[n_dom]
        km  = LOESS_WINDOWS[n_dom]
        ax.plot(lc["x"], lc["rate"],
                color=col, lw=2.6, ls="-", zorder=6,
                label=f"CoastSat LOESS {n_dom}-domain ({km} km)")
        data_handles.append(
            Line2D([0], [0], color=col, lw=2.6, ls="-",
                   label=f"CoastSat LOESS {n_dom}-domain ({km} km)")
        )

    # Faint raw CoastSat legend proxy (so readers know it's there)
    data_handles.append(
        Line2D([0], [0], color=CS_ACTIVE_COLOR, lw=1.0, alpha=0.40,
               label="CoastSat raw LRR (active period)")
    )

    # Model lines
    for k, pv in enumerate(values):
        if pv not in rate_profiles:
            continue
        col       = cmap(k / max(n_model - 1, 1))
        real_rate = rate_profiles[pv][START_REAL_INDEX:END_REAL_INDEX]
        lbl       = f"{param_cfg['label']} = {pv}{param_cfg['units']}"
        ax.plot(gis_ids, real_rate, color=col, lw=1.8, alpha=0.85, zorder=5,
                label=lbl)
        data_handles.append(Line2D([0], [0], color=col, lw=1.8, label=lbl))

    # Zero reference
    ax.axhline(0, color="#2c2c2c", lw=1.1, ls="--", alpha=0.55, zorder=3)

    # Axes formatting
    ax.set_xlim(FIRST_FILE_NUMBER - 0.5, LAST_FILE_NUMBER + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="major", direction="in", length=5, labelsize=10)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.1)

    # Accretion / erosion labels
    all_vals = np.concatenate(
        [rate_profiles[pv][START_REAL_INDEX:END_REAL_INDEX] for pv in rate_profiles]
        + ([cs["rate"] for cs in cs_series] if cs_series else [])
    )
    ypad = (all_vals.max() - all_vals.min()) * 0.07
    ax.set_ylim(all_vals.min() - ypad, all_vals.max() + ypad)

    ybot, ytop     = ax.get_ylim()
    zero_frac      = (0 - ybot) / (ytop - ybot)
    ax.text(1.0, zero_frac + (1 - zero_frac) / 2, "Accretion \u25b2",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")
    ax.text(1.0, zero_frac / 2, "Erosion \u25bc",
            transform=ax.transAxes, fontsize=9, color="#555555",
            ha="right", va="center", style="italic")

    # Labels
    ax.set_xlabel("CASCADE Model Domain (500 m alongshore)",
                  fontsize=11, fontweight="bold", labelpad=6)
    ax.set_ylabel("Shoreline Change Rate (m/yr)",
                  fontsize=11, fontweight="bold", labelpad=6)
    ax.text(0.0, 1.01, "\u2190 S  |  Cape Hatteras",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="left", va="bottom", style="italic", clip_on=False)
    ax.text(1.0, 1.01, "Pea Island  |  N \u2192",
            transform=ax.transAxes, fontsize=9, color="#444444",
            ha="right", va="bottom", style="italic", clip_on=False)
    ax.set_title(
        f"Wave Climate Sensitivity: {param_cfg['label']}  |  Hatteras Island, NC  |  "
        f"{START_YEAR}\u2013{END_YEAR}  |  Baseline: "
        f"Hs={BASELINE_WAVE_HEIGHT}m  Tp={BASELINE_WAVE_PERIOD}s  "
        f"asym={BASELINE_WAVE_ASYMMETRY}  \u03b1={BASELINE_WAVE_ANGLE_HIGH_FRACTION}",
        fontsize=11, fontweight="bold", pad=10, color="#1a2a3a",
    )

    # Legend
    ax.legend(
        handles=data_handles + annotation_legend_handles(),
        loc="lower center", bbox_to_anchor=(0.5, 0.02),
        fontsize=8.5, framealpha=0.95, edgecolor="#cccccc", ncol=4,
    )

    out = os.path.join(output_dir, f"{param_name}_sensitivity.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"  Plot: {out}")
    plt.close(fig)


def _plot_session_overview(all_rate_profiles, cs_series, session_dir):
    """
    2x2 figure: one panel per wave parameter, all runs overlaid.
    CoastSat active period shown as reference on every panel.
    """
    params = [k for k in SENSITIVITY_ANALYSES if SENSITIVITY_ANALYSES[k]["enabled"]
              and k in all_rate_profiles and all_rate_profiles[k]]
    if not params:
        return

    ncols = 2
    nrows = (len(params) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(14 * ncols / 2, 5.5 * nrows),
                             constrained_layout=True)
    fig.patch.set_facecolor("white")
    axes_flat = np.array(axes).flatten()

    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    for ax_idx, param_name in enumerate(params):
        ax         = axes_flat[ax_idx]
        param_cfg  = SENSITIVITY_ANALYSES[param_name]
        profiles   = all_rate_profiles[param_name]
        n_model    = len(profiles)
        cmap       = plt.get_cmap("viridis", n_model)

        ax.set_facecolor("white")
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(1.1)

        # CoastSat active period reference
        for cs in cs_series:
            if cs["active"]:
                gis_x = cs["x"] - START_REAL_INDEX + FIRST_FILE_NUMBER
                ax.plot(gis_x, cs["rate"], color=CS_ACTIVE_COLOR,
                        lw=1.6, ls="-", zorder=4, label="CoastSat (active)")

        # Model lines
        for k, pv in enumerate(param_cfg["values"]):
            if pv not in profiles:
                continue
            col       = cmap(k / max(n_model - 1, 1))
            real_rate = profiles[pv][START_REAL_INDEX:END_REAL_INDEX]
            ax.plot(gis_ids, real_rate, color=col, lw=1.5, alpha=0.82, zorder=5,
                    label=f"{pv}{param_cfg['units']}")

        ax.axhline(0, color="#2c2c2c", lw=0.9, ls="--", alpha=0.5)

        ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
        ax.tick_params(axis="both", which="major", direction="in", length=4, labelsize=9)
        ax.tick_params(axis="both", which="minor", direction="in", length=2)
        ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.35, color="gray")

        ax.set_xlabel("GIS Domain ID", fontsize=9, fontweight="bold")
        ax.set_ylabel("Rate (m/yr)", fontsize=9, fontweight="bold")
        ax.set_title(f"Sensitivity to {param_cfg['label']}", fontsize=10,
                     fontweight="bold", color="#1a2a3a")
        ax.legend(loc="lower right", fontsize=7.5, ncol=2, framealpha=0.9)

    # Hide any unused panels
    for ax in axes_flat[len(params):]:
        ax.set_visible(False)

    fig.suptitle(
        f"Wave Climate Sensitivity Overview  |  Hatteras Island, NC  |  "
        f"{START_YEAR}\u2013{END_YEAR}",
        fontsize=13, fontweight="bold", color="#1a2a3a",
    )

    out = os.path.join(session_dir, "session_overview.png")
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"\nSession overview: {out}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    session_ts  = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    session_dir = os.path.join(
        OUTPUT_BASE_DIR,
        f"HAT_{START_YEAR}_{END_YEAR}_waveSensitivity_{session_ts}"
    )
    os.makedirs(session_dir, exist_ok=True)

    print(f"\nSession folder: {session_dir}\n")

    # ── Load CoastSat once ────────────────────────────────────────────────────
    print("Loading CoastSat reference data...")
    cs_series = []
    for ds in COASTSAT_DATASETS:
        x, r = load_coastsat(ds)
        if x is None:
            continue
        cs_series.append(dict(
            label  = ds["label"],
            x      = x,
            rate   = r,
            active = (ds["period_start"] == START_YEAR),
        ))
        print(f"  {ds['label']}: {len(r)} domains  "
              f"{'(active)' if ds['period_start'] == START_YEAR else '(ref)'}")
    print()

    # ── Compute LOESS smoothed reference curves (7-domain and 10-domain) ──────
    print("Computing LOESS smoothed reference curves...")
    loess_curves = compute_loess_curves(cs_series)
    print()

    # ── Run each enabled sweep ────────────────────────────────────────────────
    start_time        = datetime.datetime.now()
    all_rate_profiles = {}
    completed         = []

    for param_name, param_cfg in SENSITIVITY_ANALYSES.items():
        if not param_cfg["enabled"]:
            print(f"  Skipping {param_cfg['label']} (disabled)")
            continue

        param_output_dir = os.path.join(session_dir, param_name)

        try:
            profiles = run_sweep(
                param_name       = param_name,
                param_cfg        = param_cfg,
                param_output_dir = param_output_dir,
                cs_series        = cs_series,
                loess_curves     = loess_curves,
            )
            all_rate_profiles[param_name] = profiles
            completed.append(param_name)
        except Exception as e:
            print(f"\nFAILED: {param_cfg['label']}: {e}")
            import traceback; traceback.print_exc()

    # ── Session overview figure ───────────────────────────────────────────────
    if completed:
        _plot_session_overview(all_rate_profiles, cs_series, session_dir)

    # ── Session metadata ──────────────────────────────────────────────────────
    duration = datetime.datetime.now() - start_time
    meta_lines = [
        "# Wave Climate Sensitivity Session Metadata",
        f"timestamp             = {session_ts}",
        f"session_dir           = {session_dir}",
        f"period                = {START_YEAR}-{END_YEAR}",
        f"total_runs            = {total_runs}",
        f"completed_parameters  = {len(completed)}",
        f"duration              = {str(duration).split('.')[0]}",
        "",
        "[Baseline]",
        f"Hs                    = {BASELINE_WAVE_HEIGHT} m",
        f"Tp                    = {BASELINE_WAVE_PERIOD} s",
        f"asymmetry             = {BASELINE_WAVE_ASYMMETRY}",
        f"angle_high_fraction   = {BASELINE_WAVE_ANGLE_HIGH_FRACTION}",
        "",
        "[Fixed Parameters]",
        f"SLR                   = {SEA_LEVEL_RISE_RATE * 1000:.2f} mm/yr",
        f"berm_elevation        = {BERM_ELEVATION} m",
        f"mhw_elevation         = {MHW_ELEVATION} m",
        f"storm_file            = {os.path.basename(STORM_FILE)}",
        f"background_erosion    = 0.0 (all domains)",
        f"save_npz              = {SAVE_CASCADE_NPZ}",
        "",
        "[Sweeps]",
    ]
    for pname in SENSITIVITY_ANALYSES:
        pcfg = SENSITIVITY_ANALYSES[pname]
        done = pname in completed
        meta_lines.append(
            f"  {'OK' if done else 'FAILED':6s}  {pcfg['label']:35s}  {pcfg['values']}"
        )

    meta_path = os.path.join(session_dir, "session_metadata.txt")
    with open(meta_path, "w") as f:
        f.write("\n".join(meta_lines) + "\n")
    print(f"Session metadata: {meta_path}")

    # ── Final summary ─────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("WAVE SENSITIVITY ANALYSIS COMPLETE")
    print("=" * 80)
    print(f"Duration:  {str(duration).split('.')[0]}")
    print(f"Completed: {len(completed)} / {len(enabled_params)} parameters")
    for pname in completed:
        n_ok = len(all_rate_profiles.get(pname, {}))
        n_total = len(SENSITIVITY_ANALYSES[pname]["values"])
        print(f"  {SENSITIVITY_ANALYSES[pname]['label']:35s}  {n_ok}/{n_total} runs")
    print(f"\nAll outputs: {session_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
