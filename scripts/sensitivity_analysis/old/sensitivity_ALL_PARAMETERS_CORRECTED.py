#!/usr/bin/env python3
"""
CASCADE Hatteras Island: Complete Multi-parameter Sensitivity Analysis
======================================================================

This script automatically runs CASCADE sensitivity analyses for ALL wave parameters:
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

NUM_REAL_DOMAINS = 90
NUM_BUFFER_DOMAINS = 15

FIRST_DOMAIN_NUMBER = 1
LAST_DOMAIN_NUMBER = 90

# Road Configuration
FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN = 90

# =============================================================================
# AUTOMATIC CALCULATIONS
# =============================================================================

TOTAL_DOMAINS = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS

START_ROAD_INDEX = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
END_ROAD_INDEX = (LAST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS + 1

print("=" * 80)
print("CASCADE COMPLETE SENSITIVITY ANALYSIS - HATTERAS ISLAND")
print("=" * 80)
print(f"Real Domains: {NUM_REAL_DOMAINS} domains (domain {FIRST_DOMAIN_NUMBER} to {LAST_DOMAIN_NUMBER})")
print(f"Buffer Domains: {NUM_BUFFER_DOMAINS} on each side")
print(f"Total Domains: {TOTAL_DOMAINS} (including buffers)")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 2: FILE PATHS
# =============================================================================

PROJECT_BASE_DIR = r'/'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')

DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE,
    'island_offset',
    'hindcast_1978_1997',
    f'Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv'
)

STORM_FILE_1978_1997 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1978_1997.npy')
STORM_FILE_1997_2019 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1997_2019.npy')
ROAD_SETBACK_FILE = os.path.join(HATTERAS_DATA_BASE, 'roads', 'offset', '1978', 'RoadSetback_1978.csv')

PARAMETER_FILE = 'Hatteras-CASCADE-parameters.yaml'  # Must match filename in your data/hatteras_init/

DSAS_DATA_FILE = os.path.join(HATTERAS_DATA_BASE, 'shoreline_change', 'dsas_1978_1997_domain_means_SIMPLE.csv')

# =============================================================================
# SECTION 3: SIMULATION PARAMETERS
# =============================================================================

START_YEAR = 1978
END_YEAR = 1997
RUN_YEARS = 19
TO_METERS = True

# Management features (disabled for natural state)
ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_NOURISHMENT = False
ENABLE_SANDBAG_PLACEMENT = False
ENABLE_DUNE_REBUILDING = False

SEA_LEVEL_RISE_RATE = 0.004  # m/yr (4 mm/yr)
DUNE_REBUILD_HEIGHT = 3.0
REBUILD_ELEV_THRESHOLD = 0.01  # dam
NUM_CORES = 4

# =============================================================================
# SECTION 4: SENSITIVITY ANALYSIS CONFIGURATION
# =============================================================================

# Baseline values (used when parameter is NOT being varied)
BASELINE_WAVE_HEIGHT = 2.0  # ← Changed from 1.0 to Murray's suggested value
BASELINE_WAVE_PERIOD = 7
BASELINE_WAVE_ASYMMETRY = 0.7  # ← Changed from 0.8 to mid-range
BASELINE_WAVE_ANGLE_HIGH_FRACTION = 0.2

SENSITIVITY_ANALYSES = {
    'wave_height': {
        'values': [1.5, 2.0, 2.5],  # ← Focus on Murray's range
        'label': 'Wave Height',
        'units': 'm',
        'enabled': True,
    },
    'wave_period': {
        'values': [6, 7, 8, 10],  # ← Simplified
        'label': 'Wave Period',
        'units': 's',
        'enabled': True,
    },
    'wave_asymmetry': {
        'values': [0.6, 0.7, 0.8],  # ← Simplified
        'label': 'Wave Asymmetry',
        'units': '',
        'enabled': True,
    },
    'wave_angle_high_fraction': {
        'values': [0.1, 0.2, 0.4],  # ← Simplified
        'label': 'Wave Angle High Fraction',
        'units': '',
        'enabled': True,
    },
}

# Count how many analyses will run
enabled_count = sum(1 for config in SENSITIVITY_ANALYSES.values() if config['enabled'])
total_runs = sum(len(config['values']) for config in SENSITIVITY_ANALYSES.values() if config['enabled'])

print("Sensitivity Analysis Configuration:")
print(f"  Parameters to analyze: {enabled_count} / {len(SENSITIVITY_ANALYSES)}")
print(f"  Total CASCADE runs: {total_runs}")
for param, config in SENSITIVITY_ANALYSES.items():
    status = "✓ ENABLED" if config['enabled'] else "✗ DISABLED"
    if config['enabled']:
        print(f"    {status}: {config['label']} - {len(config['values'])} values")
print(f"\nBaseline Parameters (when not varied):")
print(f"  Wave Height: {BASELINE_WAVE_HEIGHT} m")
print(f"  Wave Period: {BASELINE_WAVE_PERIOD} s")
print(f"  Wave Asymmetry: {BASELINE_WAVE_ASYMMETRY}")
print(f"  Wave Angle High Fraction: {BASELINE_WAVE_ANGLE_HIGH_FRACTION}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 5: DYNAMIC CONFIGURATION
# =============================================================================

if START_YEAR == 1978:
    YEAR_COLUMN_INDEX = 0
    RUN_NAME_BASE = f'HAT_{START_YEAR}_{END_YEAR}_natural'
    STORM_FILE = STORM_FILE_1978_1997
elif START_YEAR == 1997:
    YEAR_COLUMN_INDEX = 1
    RUN_NAME_BASE = f'HAT_{START_YEAR}_{END_YEAR}_natural'
    STORM_FILE = STORM_FILE_1997_2019
else:
    print(f"❌ ERROR: Invalid START_YEAR {START_YEAR}. Must be 1978 or 1997.")
    sys.exit(1)

os.chdir(PROJECT_BASE_DIR)

# =============================================================================
# SECTION 6: DATA LOADING
# =============================================================================

print("Loading input data...")

try:
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10.0
    print(f"✓ Loaded dune offsets: {dune_offset_dam.size} domains")

    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=',')
    print(f"✓ Loaded road setbacks: {road_setbacks_raw.size} values")

except FileNotFoundError as e:
    print(f"❌ CRITICAL ERROR: Missing data file: {e.filename}")
    sys.exit(1)
except Exception as e:
    print(f"❌ CRITICAL ERROR loading data: {e}")
    sys.exit(1)

BACKGROUND_EROSION_RATES = [
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  # Domains 0-9
     0,  0,  0,  0,  0,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 10-19
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 20-29
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 30-39
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 40-49
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 50-59
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 60-69
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 70-79
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 80-89
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  # Domains 90-99
     -1.091,  -1.091,  -1.091,  -1.091,  -1.091,  0,  0,  0,  0,  0,  # Domains 100-109
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  # Domains 110-119
]

road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX:START_ROAD_INDEX + num_road_values] = road_setbacks_raw[:num_road_values]
print(f"✓ Road setback array prepared ({TOTAL_DOMAINS} domains)")

ROADWAY_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [False] * TOTAL_DOMAINS

# =============================================================================
# SECTION 7: ELEVATION AND DUNE FILES
# =============================================================================

print("Generating elevation profile file paths...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS = []

# Domain file numbering
FIRST_FILE_NUMBER = 1

# --- South Buffer Domains ---
for i in range(START_REAL_INDEX):
    dune_path = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy')
    elev_path = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy')
    DUNE_FILE_PATHS.append(dune_path)
    ELEVATION_FILE_PATHS.append(elev_path)

# --- Real Domains ---
for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    # Calculate domain file number
    file_num = FIRST_FILE_NUMBER + (i_list - START_REAL_INDEX)

    # File paths with correct naming: domain_{file_num}_dune_2009.npy
    dune_path = os.path.join(
        HATTERAS_DATA_BASE,
        'dunes',
        '2009',
        f'domain_{file_num}_dune_2009.npy'
    )
    elev_path = os.path.join(
        HATTERAS_DATA_BASE,
        'topography',
        '2009',
        f'domain_{file_num}_topography_2009.npy'
    )

    DUNE_FILE_PATHS.append(dune_path)
    ELEVATION_FILE_PATHS.append(elev_path)

    # Verify first file exists (one-time check)
    if i_list == START_REAL_INDEX:
        if not os.path.exists(dune_path):
            print(f"⚠️  WARNING: First dune file not found: {dune_path}")
        if not os.path.exists(elev_path):
            print(f"⚠️  WARNING: First elevation file not found: {elev_path}")

# --- North Buffer Domains ---
for i in range(END_REAL_INDEX, TOTAL_DOMAINS):
    dune_path = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy')
    elev_path = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy')
    DUNE_FILE_PATHS.append(dune_path)
    ELEVATION_FILE_PATHS.append(elev_path)

print(f"✓ Generated {len(ELEVATION_FILE_PATHS)} elevation file paths")
print(f"✓ Generated {len(DUNE_FILE_PATHS)} dune file paths")
print("=" * 80 + "\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_x_s_TS(b3d):
    """Extract shoreline time series from Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No shoreline time series found on Barrier3D object.")

def build_shoreline_matrix(cascade, to_meters=True):
    """Build absolute shoreline[t, domain] matrix from Barrier3D x_s_TS."""
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    
    x_s_first = get_x_s_TS(b3d_list[0])
    nt = len(x_s_first)
    
    shoreline_abs = np.zeros((nt, ndom))
    
    for i_dom in range(ndom):
        x_s_ts = get_x_s_TS(b3d_list[i_dom])
        if to_meters:
            shoreline_abs[:, i_dom] = np.array(x_s_ts) * 10.0
        else:
            shoreline_abs[:, i_dom] = np.array(x_s_ts)
    
    return shoreline_abs

# =============================================================================
# CASCADE SIMULATION FUNCTION
# =============================================================================

def run_cascade_simulation(
    nt, name, storm_file, alongshore_section_count, num_cores,
    beach_width_threshold, rmin, rmax, elevation_file, dune_file,
    dune_design_elevation, dune_minimum_elevation, rebuild_dune_threshold,
    road_ele, road_width, road_setback, overwash_filter, overwash_to_dune,
    nourishment_volume, background_erosion, roadway_management_on,
    beach_dune_manager_on, sandbag_management_on, sandbag_elevation,
    sea_level_rise_rate, sea_level_constant, enable_shoreline_offset,
    shoreline_offset, wave_height, wave_period, wave_asymmetry,
    wave_angle_high_fraction, save_dir,
):
    """Run a single CASCADE simulation with specified parameters."""
    
    print("\n" + "=" * 80)
    print("INITIALIZING CASCADE")
    print("=" * 80)
    
    # CASCADE requires datadir as first positional argument
    datadir = HATTERAS_DATA_BASE
    
    # CASCADE initialization matches HAT_1978_1997_natural.py exactly
    cascade = Cascade(
        datadir,  # POSITIONAL argument 1
        name,     # POSITIONAL argument 2
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file=PARAMETER_FILE,  # Just filename, not full path
        
        # Wave parameters
        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=wave_asymmetry,
        wave_angle_high_fraction=wave_angle_high_fraction,
        
        # Sea level rise
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
        
        # Erosion and time
        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,
        
        # Dune growth rates
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        
        # Computational
        num_cores=num_cores,
        
        # Management modules
        roadway_management_module=roadway_management_on,
        beach_nourishment_module=beach_dune_manager_on,
        sandbag_management_on=sandbag_management_on,
        alongshore_transport_module=True,
        community_economics_module=False,
        
        # Infrastructure parameters
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        
        # Dune parameters
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        sandbag_elevation=sandbag_elevation,
        
        # Overwash parameters
        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,
        
        # Shoreline offset
        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
        
        # Nourishment parameters
        nourishment_volume=nourishment_volume,
        nourishment_interval=None,
    )
    
    print(f"✓ CASCADE initialized successfully")
    print(f"  BermEl: {cascade.barrier3d[0].BermEl:.3f} dam ({cascade.barrier3d[0].BermEl * 10:.1f} m)")
    
    print("\n" + "=" * 80)
    print("RUNNING SIMULATION")
    print("=" * 80)
    
    dune_rebuild_threshold_val = rebuild_dune_threshold + (cascade.barrier3d[0].BermEl * 10)
    
    for time_step in range(nt - 1):
        print(f"\rYear {time_step + 1}/{nt}", end="", flush=True)
        cascade.update()
        
        if cascade.b3d_break:
            print(f"\n⚠️  Simulation stopped at year {time_step + 1}: Barrier3D break condition")
            break
        
        t = cascade.barrier3d[0].time_index
        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_nourish_now = np.zeros(alongshore_section_count)
        
        for iB3D in range(alongshore_section_count):
            if cascade.community_break[iB3D] or not beach_dune_manager_on[iB3D]:
                continue
            
            if cascade.nourishments[iB3D].beach_width[t - 1] < beach_width_threshold[iB3D]:
                tmp_nourish_now[iB3D] = 1
            
            DuneDomainCrest = cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            DuneCrestMin = (np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl) * 10
            
            if DuneCrestMin < dune_rebuild_threshold_val:
                tmp_rebuild_dune[iB3D] = 1
        
        if np.any(tmp_nourish_now):
            cascade.nourish_now = tmp_nourish_now
        if np.any(tmp_rebuild_dune):
            cascade.rebuild_dune_now = tmp_rebuild_dune
    
    print("\n✓ Simulation completed")
    
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)
    
    try:
        os.makedirs(save_dir, exist_ok=True)
        cascade.save(save_dir)
        print(f"✓ Results saved to: {save_dir}")
    except Exception as e:
        print(f"❌ Error saving results: {e}")
        raise
    
    return cascade

# =============================================================================
# RUN SINGLE SENSITIVITY ANALYSIS
# =============================================================================

def run_sensitivity_for_parameter(param_name, param_config, observed_rates_full, has_observed_data):
    """Run sensitivity analysis for a single parameter."""
    
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    run_folder_name = f"{param_name}_{timestamp}"
    
    RUN_OUTPUT_DIR = os.path.join(PROJECT_BASE_DIR, 'output', 'sensitivity_analysis', run_folder_name)
    CASCADE_RUNS_DIR = os.path.join(RUN_OUTPUT_DIR, 'cascade_runs')
    PLOTS_DIR = os.path.join(RUN_OUTPUT_DIR, 'plots')
    DATA_DIR = os.path.join(RUN_OUTPUT_DIR, 'data')
    
    os.makedirs(CASCADE_RUNS_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    
    print("\n" + "=" * 80)
    print(f"SENSITIVITY ANALYSIS: {param_config['label']}")
    print("=" * 80)
    print(f"Output: {RUN_OUTPUT_DIR}")
    print(f"Values to test: {param_config['values']}")
    print("=" * 80 + "\n")
    
    # Barrier3D parameters
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS
    BEACH_WIDTH_THRESHOLD = [30] * TOTAL_DOMAINS
    DUNE_DESIGN_ELEVATION = [DUNE_REBUILD_HEIGHT] * TOTAL_DOMAINS
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS
    REBUILD_DUNE_THRESHOLD = np.full(TOTAL_DOMAINS, 1.0)
    
    ROAD_ELEVATION = 1.45
    ROAD_WIDTH = 20.0
    OVERWASH_FILTER = 0
    OVERWASH_TO_DUNE = 9
    NOURISHMENT_VOLUME = 0
    SANDBAG_ELEVATION = 0
    
    time_span_years = END_YEAR - START_YEAR
    if time_span_years == 0:
        time_span_years = None
    
    rate_profiles = {}
    run_metadata = []
    
    for param_value in param_config['values']:
        param_str = str(param_value).replace(".", "p")
        run_name = f"{param_name}_{param_str}"
        save_dir = os.path.join(CASCADE_RUNS_DIR, run_name)
        
        print(f"\n{'='*80}")
        print(f"RUNNING {param_config['label']} = {param_value}")
        print(f"{'='*80}\n")
        
        # Set wave parameters
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
                beach_width_threshold=BEACH_WIDTH_THRESHOLD,
                rmin=RMIN,
                rmax=RMAX,
                elevation_file=ELEVATION_FILE_PATHS,
                dune_file=DUNE_FILE_PATHS,
                dune_design_elevation=DUNE_DESIGN_ELEVATION,
                dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,
                rebuild_dune_threshold=REBUILD_DUNE_THRESHOLD,
                road_ele=ROAD_ELEVATION,
                road_width=ROAD_WIDTH,
                road_setback=road_setbacks_full,
                overwash_filter=OVERWASH_FILTER,
                overwash_to_dune=OVERWASH_TO_DUNE,
                nourishment_volume=NOURISHMENT_VOLUME,
                background_erosion=BACKGROUND_EROSION_RATES,
                roadway_management_on=ROADWAY_MANAGEMENT_ON,
                beach_dune_manager_on=NOURISHMENT_MANAGEMENT_ON,
                sandbag_management_on=SANDBAG_MANAGEMENT_ON,
                sandbag_elevation=SANDBAG_ELEVATION,
                sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
                sea_level_constant=True,
                enable_shoreline_offset=True,
                shoreline_offset=dune_offset_dam,
                wave_height=wh,
                wave_period=wp,
                wave_asymmetry=wa,
                wave_angle_high_fraction=wahf,
                save_dir=save_dir,
            )
            
            shoreline_abs = build_shoreline_matrix(cascade, to_meters=TO_METERS)
            nt, ndom = shoreline_abs.shape
            
            denom = time_span_years if time_span_years is not None else (nt - 1)
            
            total_change = shoreline_abs[-1, :] - shoreline_abs[0, :]
            change_rate = total_change / float(denom)
            
            rate_profiles[param_value] = change_rate
            
            run_metadata.append({
                'parameter': param_name,
                'value': param_value,
                'wave_height': wh,
                'wave_period': wp,
                'wave_asymmetry': wa,
                'wave_angle_high_fraction': wahf,
                'mean_rate': np.mean(change_rate[START_REAL_INDEX:END_REAL_INDEX]),
                'std_rate': np.std(change_rate[START_REAL_INDEX:END_REAL_INDEX]),
                'save_dir': save_dir,
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
    
    # Save results
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
                    'parameter': param_name,
                    'parameter_value': param_value,
                    'domain_id': domain_id,
                    'cascade_index': domain_idx,
                    'modeled_rate_m_per_yr': rates[domain_idx],
                    'observed_rate_m_per_yr': observed_rates_full[domain_idx] if has_observed_data else np.nan,
                })
    
    results_df = pd.DataFrame(results_records)
    results_csv_path = os.path.join(DATA_DIR, 'detailed_results.csv')
    results_df.to_csv(results_csv_path, index=False)
    print(f"✓ Saved: {results_csv_path}")
    
    summary_df = pd.DataFrame(run_metadata)
    summary_csv_path = os.path.join(DATA_DIR, 'run_summary.csv')
    summary_df.to_csv(summary_csv_path, index=False)
    print(f"✓ Saved: {summary_csv_path}")
    
    # Generate plot
    print(f"\n{'='*80}")
    print("GENERATING PLOT")
    print(f"{'='*80}")
    
    domain_indices = np.arange(TOTAL_DOMAINS)
    fig, ax = plt.subplots(figsize=(16, 6))
    
    # Color options in comments
    # - Navy Blue: '#1f4788' ← CURRENT
    # - Dark Red: '#c1272d'
    # - Dark Green: '#2d7f5e'
    # - Dark Purple: '#5c3c92'
    # - Charcoal: '#3d3d3d'
    # - Black: 'black'
    OBSERVED_COLOR = '#1f4788'
    
    if has_observed_data:
        ax.plot(domain_indices, observed_rates_full, 'o-', 
                color=OBSERVED_COLOR, linewidth=2.5, markersize=5, 
                label='Observed (DSAS)', zorder=10, alpha=0.9)
    
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(param_config['values'])))
    for i, param_value in enumerate(param_config['values']):
        if param_value in rate_profiles:
            label = f"{param_config['label']}={param_value}{param_config['units']}"
            ax.plot(domain_indices, rate_profiles[param_value], 
                   label=label, linewidth=1.5, alpha=0.7, color=colors[i])
    
    ax.axvline(START_REAL_INDEX, linestyle="--", color='gray', linewidth=1, alpha=0.5)
    ax.axvline(END_REAL_INDEX - 1, linestyle="--", color='gray', linewidth=1, alpha=0.5)
    ax.axhline(0.0, linestyle="--", color='gray', linewidth=1, alpha=0.7)
    
    tick_step = 5
    xticks = np.arange(0, TOTAL_DOMAINS, tick_step)
    ax.set_xticks(xticks)
    ax.set_xticklabels([str(i) for i in xticks], rotation=0)
    
    ax.set_xlabel("CASCADE domain index", fontsize=11)
    ax.set_ylabel("Shoreline change rate (m/yr)", fontsize=11)
    ax.set_title(
        f"Hatteras Island: Modeled vs. Observed Shoreline Change Rate\nSensitivity to {param_config['label']} | {START_YEAR}–{END_YEAR}",
        fontsize=13, fontweight='bold', pad=18
    )
    ax.grid(alpha=0.3)
    ax.legend(loc='best', fontsize=9, ncol=2)
    
    top_ax = ax.secondary_xaxis('top')
    top_ax.set_xlabel("Real domain ID (1–90)", fontsize=10)
    top_tick_positions = []
    top_tick_labels = []
    for dom_id in range(1, 91, 10):
        j = START_REAL_INDEX + (dom_id - 1)
        top_tick_positions.append(j)
        top_tick_labels.append(str(dom_id))
    top_ax.set_xticks(top_tick_positions)
    top_ax.set_xticklabels(top_tick_labels, fontsize=9)
    
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig_path = os.path.join(PLOTS_DIR, f'{param_name}_sensitivity_plot.png')
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    print(f"✓ Saved plot: {fig_path}")
    plt.close()
    
    # Save run info
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
        f.write(f"Baseline Parameters:\n")
        f.write(f"  Wave Height: {BASELINE_WAVE_HEIGHT} m\n")
        f.write(f"  Wave Period: {BASELINE_WAVE_PERIOD} s\n")
        f.write(f"  Wave Asymmetry: {BASELINE_WAVE_ASYMMETRY}\n")
        f.write(f"  Wave Angle High Fraction: {BASELINE_WAVE_ANGLE_HIGH_FRACTION}\n")
        f.write(f"  Sea Level Rise Rate: {SEA_LEVEL_RISE_RATE} m/yr\n\n")
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
            domain_id = int(row['domain_id'])
            rate = row['annual_rate_m_per_yr']
            array_idx = START_REAL_INDEX + (domain_id - 1)
            observed_rates_full[array_idx] = rate
        print(f"✓ Loaded DSAS data: {len(dsas_data)} domains")
        has_observed_data = True
    except Exception as e:
        print(f"⚠️  Could not load DSAS data: {e}")
        observed_rates_full = None
        has_observed_data = False
    
    # Run each enabled sensitivity analysis
    start_time = datetime.now()
    completed_analyses = []
    
    for param_name, param_config in SENSITIVITY_ANALYSES.items():
        if not param_config['enabled']:
            print(f"\n⏭️  Skipping {param_config['label']} (disabled)")
            continue
        
        try:
            run_sensitivity_for_parameter(param_name, param_config, observed_rates_full, has_observed_data)
            completed_analyses.append(param_name)
        except Exception as e:
            print(f"\n❌ FAILED: {param_config['label']}")
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
    
    end_time = datetime.now()
    duration = end_time - start_time
    
    # Final summary
    print("\n" + "=" * 80)
    print("🎉 ALL SENSITIVITY ANALYSES COMPLETE")
    print("=" * 80)
    print(f"Total time: {duration}")
    print(f"Completed analyses: {len(completed_analyses)} / {enabled_count}")
    for param_name in completed_analyses:
        print(f"  ✓ {SENSITIVITY_ANALYSES[param_name]['label']}")
    print(f"\nAll outputs in: {os.path.join(PROJECT_BASE_DIR, 'output', 'sensitivity_analysis')}")
    print("=" * 80)

if __name__ == '__main__':
    main()
