#!/usr/bin/env python3
"""
CASCADE Hatteras Island: Complete Multi-parameter Sensitivity Analysis
======================================================================

Updated to match HAT_hindcast_1978_1997.py (corrected version):
  - FLIP_SIGN_MODEL = True  → change_rate *= -1.0 applied after rate computation
  - BACKGROUND_EROSION_VALUE = 0  (all domains, matching hindcast)
  - SEA_LEVEL_RISE_RATE = 0.003 m/yr (3 mm/yr, matching hindcast)
  - Baseline wave params updated: Hs=2.5, asymmetry=0.8, angle_high_fraction=0.1
  - Simplified run loop (no inactive dune-rebuild management logic)
  - Background erosion array built programmatically (not hardcoded)

Sensitivity analyses for ALL wave parameters:
  - Wave Height
  - Wave Period
  - Wave Asymmetry
  - Wave Angle High Fraction

Each parameter gets its own timestamped output folder with:
  - CASCADE NPZ outputs
  - Sensitivity plot with DSAS overlay
  - CSV data files
  - Run metadata

Author: Hannah Henry (UNC Chapel Hill)
Date: February 2026

Results for the analysis as of 2/20/2026

======================================================================
"""

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from cascade.cascade import Cascade
from datetime import datetime

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION - HATTERAS ISLAND
# =============================================================================

NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15

FIRST_FILE_NUMBER = 1
LAST_FILE_NUMBER  = 90

FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN  = 90

# Automatic calculations
TOTAL_DOMAINS    = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX = NUM_BUFFER_DOMAINS                                           # 15
END_REAL_INDEX   = START_REAL_INDEX + NUM_REAL_DOMAINS                         # 105

START_ROAD_INDEX = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS   # 23
END_ROAD_INDEX   = (LAST_ROAD_DOMAIN  - 1) + NUM_BUFFER_DOMAINS + 1  # 105

print("=" * 80)
print("CASCADE COMPLETE SENSITIVITY ANALYSIS - HATTERAS ISLAND")
print("=" * 80)
print(f"Real Domains:   {NUM_REAL_DOMAINS} (GIS IDs {FIRST_FILE_NUMBER}..{LAST_FILE_NUMBER})")
print(f"Buffer Domains: {NUM_BUFFER_DOMAINS} each side")
print(f"Total Domains:  {TOTAL_DOMAINS} (including buffers)")
print(f"Real index span (padded): [{START_REAL_INDEX}..{END_REAL_INDEX - 1}]")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR   = r'C:\Users\hanna\PycharmProjects\CASCADE'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, 'output', 'sensitivity_analysis')

DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE,
    'island_offset',
    'hindcast_1978_1997',
    f'Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv',
)

STORM_FILE_1978_1997 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1978_1997.npy')
STORM_FILE_1997_2019 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1997_2019.npy')
ROAD_SETBACK_FILE    = os.path.join(HATTERAS_DATA_BASE, 'roads', 'offset', '1978', 'RoadSetback_1978.csv')

PARAMETER_FILE = 'Hatteras-CASCADE-parameters.yaml'

DSAS_DATA_FILE = os.path.join(
    HATTERAS_DATA_BASE, 'shoreline_change', 'dsas_1978_1997_domain_means_SIMPLE.csv'
)

# =============================================================================
# SECTION 3: SIMULATION PARAMETERS
# =============================================================================

START_YEAR = 1978
END_YEAR   = 1997
RUN_YEARS  = END_YEAR - START_YEAR   # 19

TO_METERS = True

# Management features (all disabled for unmanaged baseline)
ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_NOURISHMENT        = False
ENABLE_SANDBAG_PLACEMENT  = False

# --- UPDATED to match hindcast ---
SEA_LEVEL_RISE_RATE = 0.003          # m/yr (3 mm/yr) — was 0.004 in previous sensitivity
BACKGROUND_EROSION_VALUE = 0         # dam/yr — was -1.091 in previous sensitivity
FLIP_SIGN_MODEL = True               # flips only modeled sign (no alongshore reversal)

DUNE_REBUILD_HEIGHT    = 3.0
REBUILD_ELEV_THRESHOLD = 0.01        # dam
NUM_CORES = 4

# =============================================================================
# SECTION 4: SENSITIVITY ANALYSIS CONFIGURATION
# =============================================================================

# Baseline values (used when a parameter is NOT being varied)
# --- UPDATED to match hindcast calibrated values ---
BASELINE_WAVE_HEIGHT              = 2.5   # was 2.0
BASELINE_WAVE_PERIOD              = 7
BASELINE_WAVE_ASYMMETRY           = 0.8   # was 0.7
BASELINE_WAVE_ANGLE_HIGH_FRACTION = 0.1   # was 0.2

SENSITIVITY_ANALYSES = {
    'wave_height': {
        'values': [1.5, 2.0, 2.5],
        'label': 'Wave Height',
        'units': 'm',
        'enabled': True,
    },
    'wave_period': {
        'values': [6, 7, 8, 10],
        'label': 'Wave Period',
        'units': 's',
        'enabled': True,
    },
    'wave_asymmetry': {
        'values': [0.6, 0.7, 0.8],
        'label': 'Wave Asymmetry',
        'units': '',
        'enabled': True,
    },
    'wave_angle_high_fraction': {
        'values': [0.1, 0.2, 0.4],
        'label': 'Wave Angle High Fraction',
        'units': '',
        'enabled': True,
    },
}

enabled_count = sum(1 for c in SENSITIVITY_ANALYSES.values() if c['enabled'])
total_runs    = sum(len(c['values']) for c in SENSITIVITY_ANALYSES.values() if c['enabled'])

print("Sensitivity Analysis Configuration:")
print(f"  Parameters to analyze: {enabled_count} / {len(SENSITIVITY_ANALYSES)}")
print(f"  Total CASCADE runs: {total_runs}")
for param, config in SENSITIVITY_ANALYSES.items():
    if config['enabled']:
        print(f"    ✓ ENABLED: {config['label']} — {len(config['values'])} values: {config['values']}")
print(f"\nBaseline Parameters (when not varied):")
print(f"  Wave Height:              {BASELINE_WAVE_HEIGHT} m")
print(f"  Wave Period:              {BASELINE_WAVE_PERIOD} s")
print(f"  Wave Asymmetry:           {BASELINE_WAVE_ASYMMETRY}")
print(f"  Wave Angle High Fraction: {BASELINE_WAVE_ANGLE_HIGH_FRACTION}")
print(f"\nPhysical Parameters:")
print(f"  Sea Level Rise Rate:      {SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr")
print(f"  Background Erosion:       {BACKGROUND_EROSION_VALUE} dam/yr (all domains)")
print(f"  FLIP_SIGN_MODEL:          {FLIP_SIGN_MODEL}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 5: DYNAMIC CONFIGURATION
# =============================================================================

if START_YEAR == 1978:
    YEAR_COLUMN_INDEX = 0
    RUN_NAME_BASE = f'HAT_{START_YEAR}_{END_YEAR}_sensitivity'
    STORM_FILE    = STORM_FILE_1978_1997
elif START_YEAR == 1997:
    YEAR_COLUMN_INDEX = 1
    RUN_NAME_BASE = f'HAT_{START_YEAR}_{END_YEAR}_sensitivity'
    STORM_FILE    = STORM_FILE_1997_2019
else:
    print(f"❌ ERROR: Invalid START_YEAR {START_YEAR}. Must be 1978 or 1997.")
    sys.exit(1)

os.chdir(PROJECT_BASE_DIR)
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# =============================================================================
# SECTION 6: DATA LOADING
# =============================================================================

print("Loading input data...")

try:
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10.0   # m → dam
    print(f"✓ Loaded dune offsets: {dune_offset_dam.size} domains (dam)")

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=',')
    print(f"✓ Loaded road setbacks: {road_setbacks_raw.size} values")

except FileNotFoundError as e:
    print(f"❌ CRITICAL ERROR: Missing data file: {e.filename}")
    sys.exit(1)
except Exception as e:
    print(f"❌ CRITICAL ERROR loading data: {e}")
    sys.exit(1)

# Background erosion: 0 for buffer domains, BACKGROUND_EROSION_VALUE for real island
# Built programmatically to match hindcast exactly
BACKGROUND_EROSION_RATES = (
    [0.0] * START_REAL_INDEX
    + [BACKGROUND_EROSION_VALUE] * NUM_REAL_DOMAINS
    + [0.0] * (TOTAL_DOMAINS - END_REAL_INDEX)
)

# Road setback padded array
road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX:START_ROAD_INDEX + num_road_values] = \
    road_setbacks_raw[:num_road_values]
print(f"✓ Road setback array prepared ({TOTAL_DOMAINS} domains)")

ROADWAY_MANAGEMENT_ON    = [ENABLE_ROADWAY_MANAGEMENT] * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON    = [ENABLE_SANDBAG_PLACEMENT]  * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [ENABLE_NOURISHMENT]        * TOTAL_DOMAINS

# =============================================================================
# SECTION 7: ELEVATION AND DUNE FILES
# =============================================================================

print("Generating elevation + dune profile file paths...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS      = []

for _ in range(START_REAL_INDEX):  # left buffers
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy'))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))

for i_list in range(START_REAL_INDEX, END_REAL_INDEX):  # real domains
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)
    DUNE_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, 'dunes', '2009', f'domain_{file_num}_dune_2009.npy'))
    ELEVATION_FILE_PATHS.append(
        os.path.join(HATTERAS_DATA_BASE, 'topography', '2009', f'domain_{file_num}_topography_2009.npy'))

for _ in range(END_REAL_INDEX, TOTAL_DOMAINS):  # right buffers
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy'))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))

print(f"✓ Generated {len(ELEVATION_FILE_PATHS)} elevation file paths")
print(f"✓ Generated {len(DUNE_FILE_PATHS)} dune file paths")
print("=" * 80 + "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_x_s_TS(b3d):
    """Extract shoreline time series from Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.asarray(b3d.x_s_TS, dtype=float)
    if hasattr(b3d, "_x_s_TS"):
        return np.asarray(b3d._x_s_TS, dtype=float)
    raise AttributeError("No shoreline time series found (x_s_TS / _x_s_TS) on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    """Build [time × domain] shoreline matrix from CASCADE object."""
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt   = len(get_x_s_TS(b3d_list[0]))
    shoreline = np.zeros((nt, ndom), dtype=float)
    for j in range(ndom):
        shoreline[:, j] = get_x_s_TS(b3d_list[j])
    if to_meters:
        shoreline *= 10.0   # dam → m
    return shoreline

# =============================================================================
# CASCADE SIMULATION FUNCTION
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
    save_dir,
):
    """Run a single CASCADE simulation with specified parameters."""

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

    print(f"✓ CASCADE initialized")
    print(f"  BermEl: {cascade.barrier3d[0].BermEl:.3f} dam ({cascade.barrier3d[0].BermEl * 10:.1f} m)")

    # Simplified run loop (matches hindcast — no management active)
    for time_step in range(nt - 1):
        print(f"\rYear {time_step + 1}/{nt}", end="", flush=True)
        cascade.update()
        if getattr(cascade, "b3d_break", False):
            print(f"\n⚠️  Model stopped at year {time_step + 1} (b3d_break)")
            break

    os.makedirs(save_dir, exist_ok=True)
    cascade.save(save_dir)
    print(f"\n✓ Saved: {save_dir}")

    return cascade

# =============================================================================
# RUN SINGLE SENSITIVITY ANALYSIS
# =============================================================================

def run_sensitivity_for_parameter(param_name, param_config, observed_rates_full, has_observed_data):
    """Run sensitivity analysis for a single parameter."""

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_folder_name = f"{param_name}_{timestamp}"

    RUN_OUTPUT_DIR  = os.path.join(OUTPUT_BASE_DIR, run_folder_name)
    CASCADE_RUNS_DIR = os.path.join(RUN_OUTPUT_DIR, 'cascade_runs')
    PLOTS_DIR        = os.path.join(RUN_OUTPUT_DIR, 'plots')
    DATA_DIR         = os.path.join(RUN_OUTPUT_DIR, 'data')

    os.makedirs(CASCADE_RUNS_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR,        exist_ok=True)
    os.makedirs(DATA_DIR,         exist_ok=True)

    print("\n" + "=" * 80)
    print(f"SENSITIVITY ANALYSIS: {param_config['label']}")
    print("=" * 80)
    print(f"Output:         {RUN_OUTPUT_DIR}")
    print(f"Values to test: {param_config['values']}")
    print("=" * 80 + "\n")

    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS
    DUNE_DESIGN_ELEVATION  = [DUNE_REBUILD_HEIGHT]    * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS

    ROAD_ELEVATION   = 1.45
    ROAD_WIDTH       = 20.0
    OVERWASH_FILTER  = 0
    OVERWASH_TO_DUNE = 9
    NOURISHMENT_VOLUME = 0
    SANDBAG_ELEVATION  = 0

    time_span_years = END_YEAR - START_YEAR or None

    rate_profiles = {}
    run_metadata  = []

    for param_value in param_config['values']:
        param_str = str(param_value).replace(".", "p")
        run_name  = f"{RUN_NAME_BASE}_{param_name}_{param_str}"
        save_dir  = os.path.join(CASCADE_RUNS_DIR, run_name)

        print(f"\n{'='*80}")
        print(f"RUNNING {param_config['label']} = {param_value}")
        print(f"{'='*80}\n")

        # Set wave parameters — one varies, rest hold baseline from hindcast
        if param_name == 'wave_height':
            wh, wp, wa, wahf = param_value, BASELINE_WAVE_PERIOD, BASELINE_WAVE_ASYMMETRY, BASELINE_WAVE_ANGLE_HIGH_FRACTION
        elif param_name == 'wave_period':
            wh, wp, wa, wahf = BASELINE_WAVE_HEIGHT, param_value, BASELINE_WAVE_ASYMMETRY, BASELINE_WAVE_ANGLE_HIGH_FRACTION
        elif param_name == 'wave_asymmetry':
            wh, wp, wa, wahf = BASELINE_WAVE_HEIGHT, BASELINE_WAVE_PERIOD, param_value, BASELINE_WAVE_ANGLE_HIGH_FRACTION
        elif param_name == 'wave_angle_high_fraction':
            wh, wp, wa, wahf = BASELINE_WAVE_HEIGHT, BASELINE_WAVE_PERIOD, BASELINE_WAVE_ASYMMETRY, param_value

        try:
            cascade = run_cascade_simulation(
                nt=RUN_YEARS,
                name=run_name,
                storm_file=STORM_FILE,
                alongshore_section_count=TOTAL_DOMAINS,
                num_cores=NUM_CORES,
                rmin=RMIN,
                rmax=RMAX,
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
                sea_level_constant=True,
                sandbag_management_on=SANDBAG_MANAGEMENT_ON,
                sandbag_elevation=SANDBAG_ELEVATION,
                enable_shoreline_offset=True,
                shoreline_offset=dune_offset_dam,
                wave_height=wh,
                wave_period=wp,
                wave_asymmetry=wa,
                wave_angle_high_fraction=wahf,
                save_dir=save_dir,
            )

            shoreline_m = build_shoreline_matrix(cascade, to_meters=TO_METERS)
            nt_actual, ndom = shoreline_m.shape

            denom = time_span_years if time_span_years is not None else max(nt_actual - 1, 1)
            total_change = shoreline_m[-1, :] - shoreline_m[0, :]
            change_rate  = total_change / float(denom)

            # Apply sign flip — matches hindcast FLIP_SIGN_MODEL = True
            if FLIP_SIGN_MODEL:
                change_rate *= -1.0

            rate_profiles[param_value] = change_rate

            run_metadata.append({
                'parameter':              param_name,
                'value':                  param_value,
                'wave_height':            wh,
                'wave_period':            wp,
                'wave_asymmetry':         wa,
                'wave_angle_high_fraction': wahf,
                'mean_rate_real_domains': np.mean(change_rate[START_REAL_INDEX:END_REAL_INDEX]),
                'std_rate_real_domains':  np.std(change_rate[START_REAL_INDEX:END_REAL_INDEX]),
                'save_dir':               save_dir,
            })

        except Exception as e:
            print(f"\n{'='*80}")
            print(f"❌ RUN FAILED for {param_config['label']} = {param_value}: {e}")
            print(f"{'='*80}")
            import traceback
            traceback.print_exc()
            continue

    if len(rate_profiles) == 0:
        print(f"❌ No successful runs for {param_config['label']}; skipping.")
        return

    # ── Save results ─────────────────────────────────────────────────────────
    print(f"\n{'='*80}")
    print("SAVING RESULTS")
    print(f"{'='*80}")

    results_records = []
    for param_value in param_config['values']:
        if param_value in rate_profiles:
            rates = rate_profiles[param_value]
            for domain_idx in range(START_REAL_INDEX, END_REAL_INDEX):
                domain_id = (domain_idx - START_REAL_INDEX) + 1
                results_records.append({
                    'parameter':              param_name,
                    'parameter_value':        param_value,
                    'domain_id':              domain_id,
                    'cascade_index':          domain_idx,
                    'modeled_rate_m_per_yr':  rates[domain_idx],
                    'observed_rate_m_per_yr': (observed_rates_full[domain_idx]
                                               if has_observed_data else np.nan),
                })

    results_df = pd.DataFrame(results_records)
    results_csv_path = os.path.join(DATA_DIR, 'detailed_results.csv')
    results_df.to_csv(results_csv_path, index=False)
    print(f"✓ Saved: {results_csv_path}")

    summary_df = pd.DataFrame(run_metadata)
    summary_csv_path = os.path.join(DATA_DIR, 'run_summary.csv')
    summary_df.to_csv(summary_csv_path, index=False)
    print(f"✓ Saved: {summary_csv_path}")

    # ── Generate plot ─────────────────────────────────────────────────────────
    print(f"\n{'='*80}")
    print("GENERATING PLOT")
    print(f"{'='*80}")

    # x-axis: GIS domain IDs 1–90 (real domains only, matching hindcast plot style)
    gis_ids = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1)

    fig, ax = plt.subplots(figsize=(14, 5), constrained_layout=True)

    OBSERVED_COLOR = '#1f4788'
    if has_observed_data:
        obs_real = observed_rates_full[START_REAL_INDEX:END_REAL_INDEX]
        ax.plot(gis_ids, obs_real, 'o-',
                color=OBSERVED_COLOR, linewidth=2.5, markersize=5,
                label='Observed (DSAS)', zorder=10, alpha=0.9)

    colors = plt.cm.viridis(np.linspace(0, 0.9, len(param_config['values'])))
    for i, param_value in enumerate(param_config['values']):
        if param_value in rate_profiles:
            real_rate = rate_profiles[param_value][START_REAL_INDEX:END_REAL_INDEX]
            label = f"{param_config['label']}={param_value}{param_config['units']}"
            ax.plot(gis_ids, real_rate, linewidth=1.5, alpha=0.7,
                    color=colors[i], label=label)

    ax.axhline(0.0, linestyle="--", linewidth=1, color="gray", alpha=0.7)

    xticks = np.arange(FIRST_FILE_NUMBER, LAST_FILE_NUMBER + 1, 5)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

    ax.set_xlabel("GIS Domain ID (1–90)", fontsize=11)
    ax.set_ylabel("Shoreline change rate (m/yr)", fontsize=11)
    ax.set_title(
        f"Hatteras Island: Modeled vs. Observed Shoreline Change Rate\n"
        f"Sensitivity to {param_config['label']} | {START_YEAR}–{END_YEAR} | "
        f"SLR={SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr | FLIP_SIGN={FLIP_SIGN_MODEL}",
        fontsize=12, fontweight='bold',
    )
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=9, ncol=2)

    fig_path = os.path.join(PLOTS_DIR, f'{param_name}_sensitivity_plot.png')
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    print(f"✓ Saved plot: {fig_path}")
    plt.close()

    # ── Save run info ─────────────────────────────────────────────────────────
    info_path = os.path.join(RUN_OUTPUT_DIR, 'run_info.txt')
    with open(info_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CASCADE SENSITIVITY ANALYSIS RUN INFORMATION\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Run Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Study Site: Hatteras Island, North Carolina\n")
        f.write(f"Time Period: {START_YEAR}-{END_YEAR} ({RUN_YEARS} years)\n\n")
        f.write(f"Parameter Varied: {param_config['label']} ({param_name})\n")
        f.write(f"Values Tested: {param_config['values']}\n\n")
        f.write(f"Baseline Parameters (matching HAT_hindcast_1978_1997.py):\n")
        f.write(f"  Wave Height:              {BASELINE_WAVE_HEIGHT} m\n")
        f.write(f"  Wave Period:              {BASELINE_WAVE_PERIOD} s\n")
        f.write(f"  Wave Asymmetry:           {BASELINE_WAVE_ASYMMETRY}\n")
        f.write(f"  Wave Angle High Fraction: {BASELINE_WAVE_ANGLE_HIGH_FRACTION}\n")
        f.write(f"  Sea Level Rise Rate:      {SEA_LEVEL_RISE_RATE * 1000:.1f} mm/yr\n")
        f.write(f"  Background Erosion:       {BACKGROUND_EROSION_VALUE} dam/yr\n")
        f.write(f"  FLIP_SIGN_MODEL:          {FLIP_SIGN_MODEL}\n\n")
        f.write(f"Successful Runs: {len(rate_profiles)} / {len(param_config['values'])}\n")
        f.write("=" * 80 + "\n")

    print(f"✓ Saved run info: {info_path}")
    print(f"\n✅ COMPLETED: {param_config['label']}")

# =============================================================================
# MAIN
# =============================================================================

def main():
    """Run all enabled sensitivity analyses."""

    # Load DSAS observed data once
    print("\nLoading DSAS observed data...")
    try:
        dsas_data = pd.read_csv(DSAS_DATA_FILE)
        observed_rates_full = np.zeros(TOTAL_DOMAINS)
        for _, row in dsas_data.iterrows():
            domain_id  = int(row['domain_id'])
            rate       = row['annual_rate_m_per_yr']
            array_idx  = START_REAL_INDEX + (domain_id - 1)
            observed_rates_full[array_idx] = rate
        print(f"✓ Loaded DSAS data: {len(dsas_data)} domains")
        has_observed_data = True
    except Exception as e:
        print(f"⚠️  Could not load DSAS data: {e}")
        observed_rates_full = None
        has_observed_data   = False

    # Run each enabled sensitivity analysis
    start_time = datetime.now()
    completed_analyses = []

    for param_name, param_config in SENSITIVITY_ANALYSES.items():
        if not param_config['enabled']:
            print(f"\n⏭️  Skipping {param_config['label']} (disabled)")
            continue

        try:
            run_sensitivity_for_parameter(
                param_name, param_config, observed_rates_full, has_observed_data
            )
            completed_analyses.append(param_name)
        except Exception as e:
            print(f"\n❌ FAILED: {param_config['label']}")
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()

    duration = datetime.now() - start_time

    print("\n" + "=" * 80)
    print("🎉 ALL SENSITIVITY ANALYSES COMPLETE")
    print("=" * 80)
    print(f"Total time: {duration}")
    print(f"Completed:  {len(completed_analyses)} / {enabled_count}")
    for param_name in completed_analyses:
        print(f"  ✓ {SENSITIVITY_ANALYSES[param_name]['label']}")
    print(f"\nAll outputs in: {OUTPUT_BASE_DIR}")
    print("=" * 80)


if __name__ == '__main__':
    main()
