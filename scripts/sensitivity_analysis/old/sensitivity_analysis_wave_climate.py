"""
CASCADE Sensitivity Analysis: Wave Climate Parameters
======================================================

This script automatically runs CASCADE multiple times with different wave climate
parameter values to understand their influence on model output.

Purpose:
--------
- Test wave_height, wave_period, wave_asymmetry, wave_angle_high_fraction
- Keep all other parameters constant (storms, background erosion, etc.)
- Calculate metrics for each run (mean, std, R², RMSE)
- Generate comparison plots

Author: Hannah Henry (UNC Chapel Hill)
Date: January 2026
Based on: HAT_1978_1997_natural.py
======================================================
"""

import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys
import glob
from pathlib import Path
from cascade.cascade import Cascade
from datetime import datetime

# =============================================================================
# CONFIGURATION
# =============================================================================

# --- Base Configuration (from your main script) ---
PROJECT_BASE_DIR = r'/'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_BASE_DIR = os.path.join(PROJECT_BASE_DIR, 'output', 'sensitivity_analysis')

# Domain setup
NUM_REAL_DOMAINS = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS

# Simulation setup
START_YEAR = 1978
RUN_YEARS = 19
STORM_FILE = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1978_1997.npy')

# --- Baseline Wave Parameters ---
# These should match YOUR CURRENT CASCADE setup
# When testing one parameter, all others stay at these baseline values
BASELINE_PARAMS = {
    'wave_height': 1.0,
    'wave_period': 7,              # Your current value
    'wave_asymmetry': 0.8,         # Your current value
    'wave_angle_high_fraction': 0.2,  # Your current value
}

# --- Sensitivity Parameter Ranges ---
# These ranges are physically realistic for Hatteras Island
# Ranges span from current values to suggested improvements based on regional wave climate
SENSITIVITY_RANGES = {
    'wave_height': [0.5, 0.75, 1.0, 1.25, 1.5],  # meters
    
    'wave_period': [6, 7, 8, 9, 10],  # seconds
    # Current: 7s, Suggested: 8-9s (Atlantic swell-dominated)
    
    'wave_asymmetry': [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8],  # 0.5 = symmetric
    # Current: 0.8 (strongly directional), Suggested: 0.55-0.6 (more bidirectional)
    
    'wave_angle_high_fraction': [0.1, 0.2, 0.3, 0.4, 0.5],  # fraction
    # Current: 0.2, Suggested: 0.3-0.4 (more oblique wave approaches)
}

# --- DSAS Observed Data ---
# Path to your DSAS CSV file with observed shoreline change rates
# Expected columns: 'domain_id' (GIS domain 1-90), 'annual_rate_m_per_yr' (LRR rate)
DSAS_DATA_FILE = os.path.join(
    PROJECT_BASE_DIR, 
    'data', 
    'hatteras_init', 
    'shoreline_change', 
    'dsas_1978_1997_domain_means_SIMPLE.csv'
)
# Set to None if you don't have DSAS data yet - script will still work without it
# DSAS_DATA_FILE = None

# --- Computational Settings ---
NUM_CORES = 4

# --- Control Which Parameters to Test ---
# Set to False to skip a parameter (saves time if you only want to test specific ones)
TEST_WAVE_HEIGHT = True
TEST_WAVE_PERIOD = True
TEST_WAVE_ASYMMETRY = True
TEST_WAVE_ANGLE_HIGH_FRACTION = True

# =============================================================================
# DATA LOADING (copied from main script)
# =============================================================================

print("=" * 80)
print("SENSITIVITY ANALYSIS SETUP")
print("=" * 80)
print(f"Project directory: {PROJECT_BASE_DIR}")
print(f"Output directory: {OUTPUT_BASE_DIR}")
print(f"\nBaseline Parameters:")
for key, val in BASELINE_PARAMS.items():
    print(f"  {key}: {val}")
print("=" * 80 + "\n")

# Load dune offsets
DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE, 'island_offset', 'hindcast_1978_1997',
    f'Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv'
)
dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')
dune_offset_dam = dune_offset_all[:, 0] / 10  # Convert to decameters

# Load road setbacks
ROAD_SETBACK_FILE = os.path.join(HATTERAS_DATA_BASE, 'roads', 'offset', '1978', 'RoadSetback_1978.csv')
road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=',')

# Setup road setback array
FIRST_ROAD_DOMAIN = 9
LAST_ROAD_DOMAIN = 90
START_ROAD_INDEX = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
END_ROAD_INDEX = (LAST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS + 1
road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX:START_ROAD_INDEX + num_road_values] = road_setbacks_raw[:num_road_values]

# Background erosion (all zeros as per advisor requirements)
BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS

# Management flags (all False for natural state)
ROADWAY_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [False] * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON = [False] * TOTAL_DOMAINS

print("✓ Input data loaded successfully\n")

# =============================================================================
# ELEVATION AND DUNE FILE PATHS (copied from main script)
# =============================================================================

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS = []

# South buffer
for i in range(START_REAL_INDEX):
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy'))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))

# Real domains
for i_list in range(START_REAL_INDEX, END_REAL_INDEX):
    file_num = 1 + (i_list - START_REAL_INDEX)
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'dunes', '2009', f'domain_{file_num}_dune_2009.npy'))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'topography', '2009', f'domain_{file_num}_topography_2009.npy'))

# North buffer
for i in range(END_REAL_INDEX, TOTAL_DOMAINS):
    DUNE_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy'))
    ELEVATION_FILE_PATHS.append(os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy'))

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def load_dsas_observed_data():
    """
    Load DSAS observed shoreline change rates.
    Uses Hannah's DSAS file format from AllShoreline_with_DSAS_HAH_V2_LRR.py
    
    Returns:
        np.array: Observed change rates for real domains (length = NUM_REAL_DOMAINS)
        None if file doesn't exist
    """
    if DSAS_DATA_FILE is None or not os.path.exists(DSAS_DATA_FILE):
        print("⚠️  Warning: DSAS data file not found. R² and RMSE will not be calculated.")
        print(f"   Expected file: {DSAS_DATA_FILE}")
        print("   Set DSAS_DATA_FILE = None in config if you don't have this data yet.")
        return None
    
    try:
        # Load DSAS CSV
        dsas = pd.read_csv(DSAS_DATA_FILE)
        
        # Column names from Hannah's file
        # domain_id: GIS domain numbers (1-90)
        # annual_rate_m_per_yr: LRR rate in m/yr
        gis_domains = dsas['domain_id'].to_numpy()
        obs_rate_all = dsas['annual_rate_m_per_yr'].to_numpy()
        
        # Map GIS domains (1-90) → CASCADE domains (15-104)
        cascade_domains = gis_domains + NUM_BUFFER_DOMAINS  # Add 15 for left buffer
        
        # Keep only observations that fall on real island CASCADE domains
        mask = (cascade_domains > NUM_BUFFER_DOMAINS) & (cascade_domains < (NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS))
        cascade_domains_valid = cascade_domains[mask]
        obs_rate_valid = obs_rate_all[mask]
        
        # Create array with NaN for all real domains
        observed_rates = np.full(NUM_REAL_DOMAINS, np.nan)
        
        # Fill in observed values at their corresponding real domain indices
        for cascade_idx, rate in zip(cascade_domains_valid, obs_rate_valid):
            real_domain_idx = cascade_idx - NUM_BUFFER_DOMAINS  # Convert to 0-indexed real domain
            if 0 <= real_domain_idx < NUM_REAL_DOMAINS:
                observed_rates[real_domain_idx] = rate
        
        # Check how many observations we have
        n_obs = np.sum(~np.isnan(observed_rates))
        print(f"✓ Loaded DSAS data: {n_obs} observations out of {NUM_REAL_DOMAINS} domains")
        
        if n_obs < NUM_REAL_DOMAINS:
            print(f"   Note: {NUM_REAL_DOMAINS - n_obs} domains have no DSAS observations")
        
        return observed_rates
        
    except Exception as e:
        print(f"⚠️  Warning: Could not load DSAS data: {e}")
        print(f"   File: {DSAS_DATA_FILE}")
        return None


def get_x_s_TS(b3d):
    """
    Get shoreline time series from a Barrier3D object.
    Handles both public and private attribute naming.
    """
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No shoreline time series found on Barrier3D object.")


def extract_shoreline_change_rates(cascade_output_dir):
    """
    Extract average shoreline change rates from CASCADE output.
    Uses Hannah's proven method from AllShoreline_with_DSAS_HAH_V2_LRR.py
    
    Args:
        cascade_output_dir (str): Path to CASCADE output directory
        
    Returns:
        np.array: Shoreline change rates for real domains only (m/yr)
    """
    # The NPZ file is named after the run, not "CASCADE.npz"
    # Look for any .npz file in the directory
    import glob
    npz_files = glob.glob(os.path.join(cascade_output_dir, '*.npz'))
    
    if len(npz_files) == 0:
        raise FileNotFoundError(f"No .npz file found in: {cascade_output_dir}")
    
    # Use the first (should be only) npz file found
    cascade_file = npz_files[0]
    
    # Load CASCADE object
    data = np.load(cascade_file, allow_pickle=True)
    cascade = data["cascade"][0]
    
    # Build shoreline matrix using Hannah's method
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt = len(get_x_s_TS(b3d_list[0]))
    
    shoreline = np.zeros((nt, ndom))
    
    for j in range(ndom):
        xs = get_x_s_TS(b3d_list[j])
        shoreline[:, j] = xs
    
    # Convert dam → m
    shoreline = shoreline * 10.0
    
    # Calculate change rates for real domains only (exclude buffers)
    initial_position = shoreline[0, START_REAL_INDEX:END_REAL_INDEX]
    final_position = shoreline[-1, START_REAL_INDEX:END_REAL_INDEX]
    
    # Change rate = (final - initial) / years
    change_rates = (final_position - initial_position) / RUN_YEARS  # m/yr
    
    return change_rates


def calculate_metrics(modeled, observed=None):
    """
    Calculate performance metrics for model output.
    Handles NaN values in observed data (domains without DSAS observations).
    
    Args:
        modeled (np.array): Modeled change rates (m/yr)
        observed (np.array): Observed change rates (m/yr), optional, may contain NaN
        
    Returns:
        dict: Dictionary of metrics
    """
    metrics = {
        'mean_change': np.mean(modeled),
        'median_change': np.median(modeled),
        'std_dev': np.std(modeled),
        'min_change': np.min(modeled),
        'max_change': np.max(modeled),
        'range': np.max(modeled) - np.min(modeled),
    }
    
    # Add comparison metrics if observed data available
    if observed is not None and len(observed) == len(modeled):
        # Remove NaN values for comparison
        valid_mask = ~np.isnan(observed)
        
        if np.sum(valid_mask) > 0:
            obs_valid = observed[valid_mask]
            mod_valid = modeled[valid_mask]
            
            # R-squared
            ss_res = np.sum((obs_valid - mod_valid) ** 2)
            ss_tot = np.sum((obs_valid - np.mean(obs_valid)) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan
            
            # RMSE
            rmse = np.sqrt(np.mean((obs_valid - mod_valid) ** 2))
            
            # Mean absolute error
            mae = np.mean(np.abs(obs_valid - mod_valid))
            
            # Bias
            bias = np.mean(mod_valid - obs_valid)
            
            metrics.update({
                'r_squared': r_squared,
                'rmse': rmse,
                'mae': mae,
                'bias': bias,
                'n_comparisons': len(obs_valid),
            })
        else:
            print("   Warning: No valid observations for comparison")
    
    return metrics


def run_single_cascade_simulation(wave_params, run_name):
    """
    Run CASCADE with specified wave parameters.
    
    Args:
        wave_params (dict): Dictionary of wave parameters
        run_name (str): Name for this simulation run
        
    Returns:
        str: Path to output directory
    """
    print(f"\n  Running: {run_name}")
    print(f"    wave_height={wave_params['wave_height']}")
    print(f"    wave_period={wave_params['wave_period']}")
    print(f"    wave_asymmetry={wave_params['wave_asymmetry']}")
    print(f"    wave_angle_high_fraction={wave_params['wave_angle_high_fraction']}")
    
    save_path = os.path.join(OUTPUT_BASE_DIR, run_name)
    
    # Check if already completed - look for any .npz file
    import glob
    existing_npz = glob.glob(os.path.join(save_path, '*.npz'))
    if len(existing_npz) > 0:
        print(f"    ✓ Already completed (using existing results)")
        return save_path
    
    # Create output directory
    os.makedirs(save_path, exist_ok=True)
    
    # Change to project directory
    os.chdir(PROJECT_BASE_DIR)
    
    # Initialize CASCADE
    datadir = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
    
    cascade = Cascade(
        datadir,
        run_name,
        storm_file=STORM_FILE,
        elevation_file=ELEVATION_FILE_PATHS,
        dune_file=DUNE_FILE_PATHS,
        parameter_file="Hatteras-CASCADE-parameters.yaml",
        
        # Wave parameters (from sensitivity test)
        wave_height=wave_params['wave_height'],
        wave_period=wave_params['wave_period'],
        wave_asymmetry=wave_params['wave_asymmetry'],
        wave_angle_high_fraction=wave_params['wave_angle_high_fraction'],
        
        # Other parameters (constant)
        sea_level_rise_rate=0.005,
        sea_level_rise_constant=True,
        background_erosion=BACKGROUND_EROSION_RATES,
        alongshore_section_count=TOTAL_DOMAINS,
        time_step_count=RUN_YEARS,
        min_dune_growth_rate=[0.55] * TOTAL_DOMAINS,
        max_dune_growth_rate=[0.95] * TOTAL_DOMAINS,
        num_cores=NUM_CORES,
        
        # Management modules (all False)
        roadway_management_module=ROADWAY_MANAGEMENT_ON,
        beach_nourishment_module=NOURISHMENT_MANAGEMENT_ON,
        sandbag_management_on=SANDBAG_MANAGEMENT_ON,
        alongshore_transport_module=True,
        community_economics_module=False,
        
        # Infrastructure
        road_ele=1.45,
        road_width=20.0,
        road_setback=road_setbacks_full,
        
        # Dune parameters
        dune_design_elevation=[3.0] * TOTAL_DOMAINS,
        dune_minimum_elevation=[0.01] * TOTAL_DOMAINS,
        sandbag_elevation=[0] * TOTAL_DOMAINS,
        
        # Overwash
        overwash_filter=0,
        overwash_to_dune=9,
        
        # Shoreline offset
        enable_shoreline_offset=True,
        shoreline_offset=dune_offset_dam,
        
        # Nourishment
        nourishment_volume=0,
        nourishment_interval=None,
    )
    
    # Run simulation
    for time_step in range(RUN_YEARS - 1):
        cascade.update()
        if cascade.b3d_break:
            print(f"    ⚠️  Model break at year {time_step + 1}")
            break
    
    # Save results - try multiple methods
    try:
        cascade.save(save_path)
        
        # Verify the file was created (it will be named after run_name, not CASCADE.npz)
        import glob
        npz_files = glob.glob(os.path.join(save_path, '*.npz'))
        
        if len(npz_files) == 0:
            # If CASCADE.save didn't work, try manual save with the run name
            print(f"    ⚠️  CASCADE.save() didn't create file, trying manual save...")
            npz_filename = f"{run_name}.npz"
            np.savez_compressed(
                os.path.join(save_path, npz_filename),
                cascade=np.array([cascade], dtype=object)
            )
        
        print(f"    ✓ Completed and saved")
    except Exception as e:
        print(f"    ❌ Error saving CASCADE: {e}")
        # Still return the path so we can record the error
    
    return save_path


# =============================================================================
# MAIN SENSITIVITY ANALYSIS
# =============================================================================

def run_sensitivity_analysis():
    """
    Main function to run sensitivity analysis on all wave parameters.
    """
    print("\n" + "=" * 80)
    print("STARTING SENSITIVITY ANALYSIS")
    print("=" * 80)
    
    # Load observed data if available
    observed_rates = load_dsas_observed_data()
    
    # Storage for results
    all_results = []
    
    # Create timestamp for this analysis run
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Test each parameter
    for param_name, param_values in SENSITIVITY_RANGES.items():
        
        # Check if we should test this parameter
        skip_param = False
        if param_name == 'wave_height' and not TEST_WAVE_HEIGHT:
            skip_param = True
        elif param_name == 'wave_period' and not TEST_WAVE_PERIOD:
            skip_param = True
        elif param_name == 'wave_asymmetry' and not TEST_WAVE_ASYMMETRY:
            skip_param = True
        elif param_name == 'wave_angle_high_fraction' and not TEST_WAVE_ANGLE_HIGH_FRACTION:
            skip_param = True
            
        if skip_param:
            print(f"\n{'='*80}")
            print(f"SKIPPING: {param_name}")
            print(f"{'='*80}")
            continue
        
        print(f"\n{'='*80}")
        print(f"TESTING PARAMETER: {param_name}")
        print(f"{'='*80}")
        print(f"Values to test: {param_values}")
        
        for value in param_values:
            # Create parameter set for this run
            current_params = BASELINE_PARAMS.copy()
            current_params[param_name] = value
            
            # Generate run name
            run_name = f"sens_{param_name}_{value}_{timestamp}"
            
            try:
                # Run CASCADE
                output_dir = run_single_cascade_simulation(current_params, run_name)
                
                # Extract results
                modeled_rates = extract_shoreline_change_rates(output_dir)
                
                # Calculate metrics
                metrics = calculate_metrics(modeled_rates, observed_rates)
                
                # Store results
                result_entry = {
                    'parameter': param_name,
                    'value': value,
                    'run_name': run_name,
                    'output_dir': output_dir,
                    **metrics
                }
                all_results.append(result_entry)
                
                # Print summary
                print(f"    Metrics:")
                print(f"      Mean: {metrics['mean_change']:.3f} m/yr")
                print(f"      Std: {metrics['std_dev']:.3f} m/yr")
                print(f"      Range: {metrics['range']:.3f} m/yr")
                if 'r_squared' in metrics:
                    print(f"      R²: {metrics['r_squared']:.3f}")
                    print(f"      RMSE: {metrics['rmse']:.3f} m/yr")
                
            except FileNotFoundError as e:
                print(f"    ❌ OUTPUT FILE ERROR: {e}")
                # Store error entry with just the basic info
                all_results.append({
                    'parameter': param_name,
                    'value': value,
                    'run_name': run_name,
                    'error': 'Output file not found'
                })
            except Exception as e:
                print(f"    ❌ ERROR: {e}")
                import traceback
                traceback.print_exc()
                # Store error entry
                all_results.append({
                    'parameter': param_name,
                    'value': value,
                    'run_name': run_name,
                    'error': str(e)
                })
    
    # Save results to CSV
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)
    
    results_df = pd.DataFrame(all_results)
    results_file = os.path.join(OUTPUT_BASE_DIR, f'sensitivity_results_{timestamp}.csv')
    results_df.to_csv(results_file, index=False)
    print(f"✓ Results saved to: {results_file}")
    
    # Print summary of results
    successful_runs = results_df[~results_df.get('error', pd.Series([None]*len(results_df))).notna()]
    failed_runs = results_df[results_df.get('error', pd.Series([None]*len(results_df))).notna()]
    
    print(f"\n📊 Run Summary:")
    print(f"   Successful: {len(successful_runs)} / {len(results_df)}")
    print(f"   Failed: {len(failed_runs)} / {len(results_df)}")
    
    if len(failed_runs) > 0:
        print(f"\n⚠️  Failed runs:")
        for _, row in failed_runs.iterrows():
            print(f"   - {row['parameter']}={row['value']}: {row.get('error', 'Unknown error')}")
    
    return results_df, timestamp


# =============================================================================
# VISUALIZATION
# =============================================================================

def create_sensitivity_plots(results_df, timestamp, observed_rates=None):
    """
    Create visualization plots for sensitivity analysis results.
    
    Args:
        results_df (pd.DataFrame): Results from sensitivity analysis
        timestamp (str): Timestamp for file naming
        observed_rates (np.array): Observed DSAS rates (optional)
    """
    print("\n" + "=" * 80)
    print("GENERATING PLOTS")
    print("=" * 80)
    
    # Create plots directory
    plots_dir = os.path.join(OUTPUT_BASE_DIR, f'plots_{timestamp}')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Remove any error entries
    if 'mean_change' not in results_df.columns:
        print("❌ ERROR: No successful runs to plot. All runs failed.")
        print("   Check the output directories to see what files were created.")
        print("   The CASCADE.save() method may not be working correctly.")
        return
    
    results_df = results_df[~results_df['mean_change'].isna()].copy()
    
    if len(results_df) == 0:
        print("❌ ERROR: No successful runs to plot. All runs failed.")
        print("   Check the output directories to see what files were created.")
        return
    
    print(f"Successfully completed runs: {len(results_df)} out of 22 total")
    
    # Get unique parameters tested
    parameters = results_df['parameter'].unique()
    
    # Set style
    sns.set_style("whitegrid")
    plt.rcParams['figure.figsize'] = (12, 8)
    
    # --- Plot 1: Parameter vs Mean Change ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    axes = axes.flatten()
    
    # Calculate observed mean for reference line
    obs_mean = np.nanmean(observed_rates) if observed_rates is not None else None
    
    for idx, param in enumerate(parameters):
        if idx >= 4:
            break
        ax = axes[idx]
        param_data = results_df[results_df['parameter'] == param].sort_values('value')
        
        # Plot model mean
        ax.plot(param_data['value'], param_data['mean_change'], 'o-', linewidth=2.5, markersize=10, 
                color='steelblue', label='Model Mean', zorder=3)
        
        # Highlight current baseline value
        baseline_val = BASELINE_PARAMS[param]
        baseline_row = param_data[param_data['value'] == baseline_val]
        if len(baseline_row) > 0:
            ax.plot(baseline_val, baseline_row['mean_change'].values[0], 'r*', 
                   markersize=20, label='Current Baseline', zorder=4)
        
        # Add zero line
        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5, linewidth=1.5, label='Zero Change')
        
        # Add observed mean line if available
        if obs_mean is not None:
            ax.axhline(y=obs_mean, color='green', linestyle=':', alpha=0.7, linewidth=2, 
                      label=f'Observed Mean ({obs_mean:.2f} m/yr)')
        
        ax.set_xlabel(param.replace('_', ' ').title(), fontsize=13, fontweight='bold')
        ax.set_ylabel('Mean Shoreline Change (m/yr)', fontsize=13, fontweight='bold')
        ax.set_title(f'Sensitivity to {param.replace("_", " ").title()}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=9)
    
    plt.suptitle('Wave Parameter Sensitivity: Mean Shoreline Change', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'sensitivity_mean_change.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Saved: sensitivity_mean_change.png")
    plt.close()
    
    # --- Plot 2: Parameter vs Spatial Variability (Std Dev) ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))
    axes = axes.flatten()
    
    # Calculate observed std dev for reference line
    obs_std = np.nanstd(observed_rates) if observed_rates is not None else None
    
    for idx, param in enumerate(parameters):
        if idx >= 4:
            break
        ax = axes[idx]
        param_data = results_df[results_df['parameter'] == param].sort_values('value')
        
        # Plot model std dev
        ax.plot(param_data['value'], param_data['std_dev'], 'o-', linewidth=2.5, markersize=10, 
                color='darkorange', label='Model Std Dev', zorder=3)
        
        # Highlight current baseline value
        baseline_val = BASELINE_PARAMS[param]
        baseline_row = param_data[param_data['value'] == baseline_val]
        if len(baseline_row) > 0:
            ax.plot(baseline_val, baseline_row['std_dev'].values[0], 'r*', 
                   markersize=20, label='Current Baseline', zorder=4)
        
        # Add observed std dev line if available
        if obs_std is not None:
            ax.axhline(y=obs_std, color='green', linestyle=':', alpha=0.7, linewidth=2, 
                      label=f'Observed Std Dev ({obs_std:.2f} m/yr)')
        
        ax.set_xlabel(param.replace('_', ' ').title(), fontsize=13, fontweight='bold')
        ax.set_ylabel('Spatial Variability (Std Dev, m/yr)', fontsize=13, fontweight='bold')
        ax.set_title(f'Spatial Variability vs {param.replace("_", " ").title()}', fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=9)
    
    plt.suptitle('Wave Parameter Sensitivity: Spatial Variability', 
                 fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'sensitivity_spatial_variability.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Saved: sensitivity_spatial_variability.png")
    plt.close()
    
    # --- Plot 3: Parameter vs R² (if available) ---
    if 'r_squared' in results_df.columns and not results_df['r_squared'].isna().all():
        fig, axes = plt.subplots(2, 2, figsize=(16, 11))
        axes = axes.flatten()
        
        for idx, param in enumerate(parameters):
            if idx >= 4:
                break
            ax = axes[idx]
            param_data = results_df[results_df['parameter'] == param].sort_values('value')
            
            # Plot R²
            ax.plot(param_data['value'], param_data['r_squared'], 'o-', linewidth=2.5, markersize=10, 
                    color='forestgreen', label='R² (Model vs Observed)', zorder=3)
            
            # Highlight current baseline value
            baseline_val = BASELINE_PARAMS[param]
            baseline_row = param_data[param_data['value'] == baseline_val]
            if len(baseline_row) > 0:
                ax.plot(baseline_val, baseline_row['r_squared'].values[0], 'r*', 
                       markersize=20, label='Current Baseline', zorder=4)
            
            # Add reference lines
            ax.axhline(y=0, color='red', linestyle='--', alpha=0.7, linewidth=1.5, 
                      label='R²=0 (No better than mean)')
            ax.axhline(y=0.5, color='orange', linestyle=':', alpha=0.5, linewidth=1, 
                      label='R²=0.5 (Moderate fit)')
            ax.axhline(y=0.8, color='green', linestyle=':', alpha=0.5, linewidth=1, 
                      label='R²=0.8 (Good fit)')
            
            ax.set_xlabel(param.replace('_', ' ').title(), fontsize=13, fontweight='bold')
            ax.set_ylabel('R² with Observations', fontsize=13, fontweight='bold')
            ax.set_title(f'Model Fit vs {param.replace("_", " ").title()}', fontsize=14, fontweight='bold')
            
            # Auto-scale y-axis to show negative values
            y_min = min(param_data['r_squared'].min(), -0.5)
            y_max = max(param_data['r_squared'].max(), 0.5)
            ax.set_ylim([y_min - 0.1, y_max + 0.1])
            
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=8)
        
        plt.suptitle('Wave Parameter Sensitivity: Model Fit (R²)', 
                     fontsize=16, fontweight='bold', y=0.995)
        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, 'sensitivity_r_squared.png'), dpi=300, bbox_inches='tight')
        print(f"✓ Saved: sensitivity_r_squared.png")
        plt.close()
    
    # --- Plot 4: Comprehensive Metrics Comparison ---
    # Show mean, std dev, and R² together for easy comparison
    if observed_rates is not None:
        fig, axes = plt.subplots(3, 4, figsize=(20, 14))
        
        obs_mean = np.nanmean(observed_rates)
        obs_std = np.nanstd(observed_rates)
        
        for col_idx, param in enumerate(parameters):
            if col_idx >= 4:
                break
            param_data = results_df[results_df['parameter'] == param].sort_values('value')
            baseline_val = BASELINE_PARAMS[param]
            
            # Row 1: Mean Change
            ax1 = axes[0, col_idx]
            ax1.plot(param_data['value'], param_data['mean_change'], 'o-', linewidth=2, markersize=8, color='steelblue')
            baseline_row = param_data[param_data['value'] == baseline_val]
            if len(baseline_row) > 0:
                ax1.plot(baseline_val, baseline_row['mean_change'].values[0], 'r*', markersize=15)
            ax1.axhline(y=obs_mean, color='green', linestyle=':', linewidth=2, alpha=0.7)
            ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
            ax1.set_ylabel('Mean Change\n(m/yr)', fontsize=11, fontweight='bold')
            if col_idx == 0:
                ax1.text(-0.3, 0.5, 'Model Output\nMean', transform=ax1.transAxes, 
                        fontsize=12, fontweight='bold', va='center', rotation=90)
            ax1.grid(alpha=0.3)
            ax1.set_title(param.replace('_', ' ').title(), fontsize=12, fontweight='bold')
            
            # Row 2: Std Dev
            ax2 = axes[1, col_idx]
            ax2.plot(param_data['value'], param_data['std_dev'], 'o-', linewidth=2, markersize=8, color='darkorange')
            baseline_row = param_data[param_data['value'] == baseline_val]
            if len(baseline_row) > 0:
                ax2.plot(baseline_val, baseline_row['std_dev'].values[0], 'r*', markersize=15)
            ax2.axhline(y=obs_std, color='green', linestyle=':', linewidth=2, alpha=0.7)
            ax2.set_ylabel('Spatial Variability\n(Std Dev, m/yr)', fontsize=11, fontweight='bold')
            if col_idx == 0:
                ax2.text(-0.3, 0.5, 'Model Variability', transform=ax2.transAxes, 
                        fontsize=12, fontweight='bold', va='center', rotation=90)
            ax2.grid(alpha=0.3)
            
            # Row 3: R²
            ax3 = axes[2, col_idx]
            if 'r_squared' in param_data.columns:
                ax3.plot(param_data['value'], param_data['r_squared'], 'o-', linewidth=2, markersize=8, color='forestgreen')
                baseline_row = param_data[param_data['value'] == baseline_val]
                if len(baseline_row) > 0:
                    ax3.plot(baseline_val, baseline_row['r_squared'].values[0], 'r*', markersize=15)
                ax3.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
                y_min = min(param_data['r_squared'].min(), -0.5)
                y_max = max(param_data['r_squared'].max(), 0.5)
                ax3.set_ylim([y_min - 0.1, y_max + 0.1])
            ax3.set_ylabel('Model Fit\n(R²)', fontsize=11, fontweight='bold')
            ax3.set_xlabel('Parameter Value', fontsize=11, fontweight='bold')
            if col_idx == 0:
                ax3.text(-0.3, 0.5, 'vs Observations', transform=ax3.transAxes, 
                        fontsize=12, fontweight='bold', va='center', rotation=90)
            ax3.grid(alpha=0.3)
        
        # Add legend
        from matplotlib.lines import Line2D
        legend_elements = [
            Line2D([0], [0], marker='o', color='steelblue', linewidth=2, markersize=8, label='Model Output'),
            Line2D([0], [0], marker='*', color='red', linewidth=0, markersize=15, label='Current Baseline'),
            Line2D([0], [0], color='green', linewidth=2, linestyle=':', label='Observed Value'),
            Line2D([0], [0], color='red', linewidth=1.5, linestyle='--', label='R²=0 (Threshold)')
        ]
        fig.legend(handles=legend_elements, loc='upper center', ncol=4, fontsize=11, 
                  bbox_to_anchor=(0.5, 0.98), frameon=True, fancybox=True, shadow=True)
        
        plt.suptitle('Comprehensive Sensitivity Analysis: All Metrics', 
                     fontsize=18, fontweight='bold', y=0.995)
        plt.tight_layout(rect=[0, 0, 1, 0.97])
        plt.savefig(os.path.join(plots_dir, 'sensitivity_comprehensive.png'), dpi=300, bbox_inches='tight')
        print(f"✓ Saved: sensitivity_comprehensive.png")
        plt.close()
    
    # --- Plot 5: Tornado Diagram (Relative Sensitivity) ---
    # Calculate range of outputs for each parameter
    sensitivity_summary = []
    for param in parameters:
        param_data = results_df[results_df['parameter'] == param]
        output_range = param_data['mean_change'].max() - param_data['mean_change'].min()
        sensitivity_summary.append({
            'parameter': param,
            'output_range': output_range
        })
    
    sensitivity_df = pd.DataFrame(sensitivity_summary).sort_values('output_range', ascending=True)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(sensitivity_df))
    ax.barh(y_pos, sensitivity_df['output_range'], color='steelblue')
    ax.set_yticks(y_pos)
    ax.set_yticklabels([p.replace('_', ' ').title() for p in sensitivity_df['parameter']])
    ax.set_xlabel('Range of Mean Change (m/yr)', fontsize=12)
    ax.set_title('Relative Parameter Sensitivity (Tornado Diagram)', fontsize=14, fontweight='bold')
    ax.grid(True, axis='x', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'tornado_diagram.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Saved: tornado_diagram.png")
    plt.close()
    
    # --- Plot 5: Summary Table ---
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.axis('tight')
    ax.axis('off')
    
    # Create summary table
    summary_data = []
    for param in parameters:
        param_data = results_df[results_df['parameter'] == param]
        row = [
            param.replace('_', ' ').title(),
            f"{param_data['mean_change'].min():.3f} to {param_data['mean_change'].max():.3f}",
            f"{param_data['std_dev'].mean():.3f}",
            f"{param_data['range'].mean():.3f}",
        ]
        if 'r_squared' in param_data.columns:
            row.append(f"{param_data['r_squared'].mean():.3f}")
        summary_data.append(row)
    
    columns = ['Parameter', 'Mean Change Range', 'Avg Std Dev', 'Avg Range']
    if 'r_squared' in results_df.columns:
        columns.append('Avg R²')
    
    table = ax.table(cellText=summary_data, colLabels=columns, cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header
    for i in range(len(columns)):
        table[(0, i)].set_facecolor('#4472C4')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    plt.title('Sensitivity Analysis Summary', fontsize=14, fontweight='bold', pad=20)
    plt.savefig(os.path.join(plots_dir, 'summary_table.png'), dpi=300, bbox_inches='tight')
    print(f"✓ Saved: summary_table.png")
    plt.close()
    
    print(f"\n✓ All plots saved to: {plots_dir}")


# =============================================================================
# RUN SCRIPT
# =============================================================================

if __name__ == '__main__':
    try:
        # Run sensitivity analysis
        results_df, timestamp = run_sensitivity_analysis()
        
        # Load observed data for plotting
        observed_rates = load_dsas_observed_data()
        
        # Create plots
        create_sensitivity_plots(results_df, timestamp, observed_rates)
        
        print("\n" + "=" * 80)
        print("✓ SENSITIVITY ANALYSIS COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"\nResults saved to: {OUTPUT_BASE_DIR}")
        print(f"Check the 'plots_{timestamp}' folder for visualizations")
        
    except Exception as e:
        print("\n" + "=" * 80)
        print(f"❌ SENSITIVITY ANALYSIS FAILED: {e}")
        print("=" * 80)
        import traceback
        traceback.print_exc()
        sys.exit(1)
