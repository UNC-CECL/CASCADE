"""
CASCADE Model: Hatteras Barrier Island Natural State Simulation

This script runs a hindcast simulation of barrier island evolution WITHOUT
human management interventions (roads, nourishment, etc.).

============================================================================
MODIFY THESE NUMBERS FOR YOUR STUDY AREA
============================================================================

To change your domain setup, update SECTION 1 below:
1. Set NUM_REAL_DOMAINS to your total number of real domains
2. Set FIRST_ROAD_DOMAIN and LAST_ROAD_DOMAIN for where roads exist
3. NUM_BUFFER_DOMAINS is always 15 on each side (don't change unless needed)

The script automatically calculates all indices and array sizes!

============================================================================
Author: Hannah Henry (UNC Chapel Hill)
Date: December 2025
CASCADE Version: 1
============================================================================
"""

import copy
import numpy as np
import os
import sys
from cascade.cascade import Cascade

# =============================================================================
# SECTION 1: DOMAIN CONFIGURATION - **MODIFY THIS FOR YOUR STUDY AREA**
# =============================================================================

# --- Core Domain Setup ---
NUM_REAL_DOMAINS = 90  # Total number of real barrier island domains (south to north)
NUM_BUFFER_DOMAINS = 15  # Number of buffer domains on EACH side (typically 15)

# --- Road Configuration ---
FIRST_ROAD_DOMAIN = 9  # First real domain number with road infrastructure
LAST_ROAD_DOMAIN = 90  # Last real domain number with road infrastructure

# --- Domain File Numbering ---
FIRST_FILE_NUMBER = 1  # Starting file number for domain files (usually 1)

# =============================================================================
# AUTOMATIC CALCULATIONS - **DO NOT MODIFY**
# =============================================================================

# Total domains including buffers
TOTAL_DOMAINS = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS

# Python array indices for real domains (buffers are 0 to NUM_BUFFER_DOMAINS-1)
START_REAL_INDEX = NUM_BUFFER_DOMAINS
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS

# Python array indices for road domains
# Convert domain numbers to array indices: domain N → index (N - 1 + NUM_BUFFER_DOMAINS)
START_ROAD_INDEX = (FIRST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS
END_ROAD_INDEX = (LAST_ROAD_DOMAIN - 1) + NUM_BUFFER_DOMAINS + 1  # +1 for exclusive end

# Last file number
LAST_FILE_NUMBER = FIRST_FILE_NUMBER + NUM_REAL_DOMAINS - 1

# Display configuration for verification
print("=" * 80)
print("CASCADE DOMAIN CONFIGURATION")
print("=" * 80)
print(f"Real Domains: {NUM_REAL_DOMAINS} domains (domain {FIRST_FILE_NUMBER} to {LAST_FILE_NUMBER})")
print(f"Buffer Domains: {NUM_BUFFER_DOMAINS} on each side")
print(f"Total Domains: {TOTAL_DOMAINS} (including buffers)")
print(f"\nRoad Configuration:")
print(f"  Road exists in domains: {FIRST_ROAD_DOMAIN} to {LAST_ROAD_DOMAIN}")
print(f"  Road array indices: {START_ROAD_INDEX} to {END_ROAD_INDEX - 1}")
print(f"\nArray Structure:")
print(f"  Indices 0-{NUM_BUFFER_DOMAINS - 1}: South buffer")
print(f"  Indices {START_REAL_INDEX}-{END_REAL_INDEX - 1}: Real domains")
print(f"  Indices {END_REAL_INDEX}-{TOTAL_DOMAINS - 1}: North buffer")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 2: FILE PATHS - Update these for your system
# =============================================================================

PROJECT_BASE_DIR = r'C:\Users\hanna\PycharmProjects\CASCADE'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_BASE_DIR = os.path.join(PROJECT_BASE_DIR, 'output', 'raw_runs')

# Input data files - update filenames if needed
DUNE_OFFSET_FILE = os.path.join(
    HATTERAS_DATA_BASE,
    'island_offset',
    'hindcast_1978_1997',
    f'Island_Dune_Offsets_1978_1997_PADDED_{TOTAL_DOMAINS}.csv'
)
STORM_FILE_1978_1997 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'modified', 'storms_1978_1997_1.5x_intensity.npy')
STORM_FILE_1997_2022 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'modified', 'HAT_1978_2022_Final_Hindcast_Storms.npy')
ROAD_SETBACK_FILE = os.path.join(HATTERAS_DATA_BASE, 'roads', 'offset', '1978', 'RoadSetback_1978.csv')

# =============================================================================
# SECTION 3: SIMULATION PARAMETERS - Adjust as needed
# =============================================================================

START_YEAR = 1978  # Simulation start year (1978 or 1997)
RUN_YEARS = 19  # Number of years to simulate

# --- MANAGEMENT MODULE TOGGLES ---
# Set these to True/False to enable/disable management interventions
# For "Natural State" simulation: keep all as False
# For management scenarios: set relevant modules to True

ENABLE_ROADWAY_MANAGEMENT = False  # Road relocation, abandonment decisions
ENABLE_NOURISHMENT = False  # Beach nourishment interventions
ENABLE_SANDBAG_PLACEMENT = False  # Sandbag management
ENABLE_DUNE_REBUILDING = False  # Artificial dune reconstruction

# Note: Even with management disabled, you must still provide road setback data
# (used for initial configuration even if roads aren't managed)

# --- Physical Parameters ---
SEA_LEVEL_RISE_RATE = 0.0056  # m/yr - Adjust for different SLR scenarios
DUNE_REBUILD_HEIGHT = 3.0  # m - Target height for dune rebuilding
REBUILD_ELEV_THRESHOLD = 0.01  # dam - Minimum dune elevation threshold

# --- Computational Parameters ---
NUM_CORES = 4  # Number of CPU cores for parallel processing

# =============================================================================
# SECTION 4: DYNAMIC CONFIGURATION BASED ON START YEAR
# =============================================================================

if START_YEAR == 1978:
    YEAR_COLUMN_INDEX = 0  # Column in dune offset CSV for 1978
    RUN_NAME = f'HAT_{START_YEAR}_{START_YEAR + RUN_YEARS}_mod_storms_be1'
    STORM_FILE = STORM_FILE_1978_1997

elif START_YEAR == 1997:
    YEAR_COLUMN_INDEX = 1  # Column in dune offset CSV for 1997
    RUN_NAME = f'HAT_{START_YEAR}_{START_YEAR + RUN_YEARS}_Natural_State'
    STORM_FILE = STORM_FILE_1997_2022

else:
    print(f"❌ ERROR: Invalid START_YEAR {START_YEAR}. Must be 1978 or 1997.")
    sys.exit(1)

# Set working directory
os.chdir(PROJECT_BASE_DIR)

print(f"Simulation Configuration:")
print(f"  Start Year: {START_YEAR}")
print(f"  Duration: {RUN_YEARS} years")
print(f"  Run Name: {RUN_NAME}")
print(f"  Output Directory: {os.path.join(OUTPUT_BASE_DIR, RUN_NAME)}")
print(f"  Sea Level Rise Rate: {SEA_LEVEL_RISE_RATE} m/yr")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 5: DATA LOADING
# =============================================================================

print("Loading input data...")

try:
    # --- Load Dune Offsets (Shoreline Positions) ---
    dune_offset_all = np.loadtxt(DUNE_OFFSET_FILE, skiprows=1, delimiter=',')

    # CRITICAL: Convert meters to decameters (CASCADE internal unit)
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10

    print(f"✓ Loaded dune offsets: {dune_offset_dam.size} domains")
    print(f"  Range: {np.min(dune_offset_dam):.1f} to {np.max(dune_offset_dam):.1f} dam")

    # Verify correct size
    if dune_offset_dam.size != TOTAL_DOMAINS:
        print(f"⚠️  WARNING: Expected {TOTAL_DOMAINS} values, got {dune_offset_dam.size}")

    # --- Load Road Setbacks ---
    road_setbacks_raw = np.loadtxt(ROAD_SETBACK_FILE, skiprows=1, delimiter=',')

    # IMPORTANT: Keep road setbacks in METERS (no conversion)
    # CASCADE expects road setbacks in meters
    print(f"✓ Loaded road setbacks: {road_setbacks_raw.size} values")
    print(f"  Range: {np.min(road_setbacks_raw):.1f} to {np.max(road_setbacks_raw):.1f} m")

    # Calculate expected number of road domains
    expected_road_domains = LAST_ROAD_DOMAIN - FIRST_ROAD_DOMAIN + 1
    if road_setbacks_raw.size != expected_road_domains:
        print(f"⚠️  WARNING: Expected {expected_road_domains} road setback values, got {road_setbacks_raw.size}")

except FileNotFoundError as e:
    print(f"❌ CRITICAL ERROR: Missing data file: {e.filename}")
    print(f"   Check that all file paths in SECTION 2 are correct.")
    sys.exit(1)
except Exception as e:
    print(f"❌ CRITICAL ERROR loading data: {e}")
    sys.exit(1)


# Load in background erosion
# Background erosion rates for CASCADE (m/yr)
# Generated from: background_erosion_rates.csv
# Domains 0-14: Left buffer (0.0)
# Domains 15-104: Real island (calculated rates)
# Domains 105-119: Right buffer (0.0)

# Using 90% and smoothing at 1
BACKGROUND_EROSION_RATES = [
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  # Domains 0-9
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000, -4.9427, -3.0023, -3.7637, -3.1265, -2.4542,  # Domains 10-19
    -0.4817,  0.0523, -1.3607, -1.2815, -1.8179, -1.3211, -1.2617, -0.9251, -0.1232, -0.5447,  # Domains 20-29
    -0.6880, -0.1997,  0.2485,  0.6391,  1.8145,  1.8181,  1.7047,  1.8919,  1.0153,  0.9793,  # Domains 30-39
    -0.7001, -0.8657, -0.0197,  0.5635,  1.6381,  2.6551,  3.5029,  2.7127,  1.3078,  0.7795,  # Domains 40-49
     0.5203,  0.4483,  0.8473, -0.5705,  0.6673,  1.3321,  0.8083,  0.6373, -0.7109, -0.4847,  # Domains 50-59
     0.3103,  0.3007, -0.1889, -0.6065, -1.4267, -1.1897,  0.4393, -0.4427, -1.0691, -1.1573,  # Domains 60-69
    -1.0853,  0.8803,  1.5661,  0.9505,  1.7911,  2.1637,  3.2563,  3.9151,  4.0789,  4.0645,  # Domains 70-79
     2.7199,  1.2673,  1.0715,  1.7857,  1.7821,  0.4550, -2.1365, -1.6451, -0.8252, -1.1015,  # Domains 80-89
     0.5383,  0.3943, -1.8647, -2.6333, -1.6645, -1.0975, -2.1383, -2.6342, -1.5677, -1.9133,  # Domains 90-99
    -1.1327, -0.9737,  0.3313, -0.7577, -0.7457,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  # Domains 100-109
     0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  # Domains 110-119
]

# FIX: CASCADE uses opposite sign convention
BACKGROUND_EROSION_RATES = [-x for x in BACKGROUND_EROSION_RATES]

# --- Initialize Road Setback Array ---
# Create array of correct length for all domains
road_setbacks_full = np.zeros(TOTAL_DOMAINS)

# Insert loaded road setbacks into the correct positions
num_road_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX:START_ROAD_INDEX + num_road_values] = road_setbacks_raw[:num_road_values]

print(f"✓ Road setback array prepared ({TOTAL_DOMAINS} domains)")

# --- Initialize Management Flag Arrays ---
# Use the toggles from SECTION 3 to set management flags
# Note: These are lists where each element corresponds to one domain
# For natural state, all values are False
# For managed scenarios, you could set specific domains to True

ROADWAY_MANAGEMENT_ON = [ENABLE_ROADWAY_MANAGEMENT] * TOTAL_DOMAINS
SANDBAG_MANAGEMENT_ON = [ENABLE_SANDBAG_PLACEMENT] * TOTAL_DOMAINS
NOURISHMENT_MANAGEMENT_ON = [ENABLE_NOURISHMENT] * TOTAL_DOMAINS

# Display management status
management_status = "ENABLED" if any([ENABLE_ROADWAY_MANAGEMENT, ENABLE_NOURISHMENT,
                                      ENABLE_SANDBAG_PLACEMENT, ENABLE_DUNE_REBUILDING]) else "DISABLED"
print(f"✓ Management arrays initialized (management: {management_status})")
if management_status == "ENABLED":
    print(f"  - Roadway management: {ENABLE_ROADWAY_MANAGEMENT}")
    print(f"  - Beach nourishment: {ENABLE_NOURISHMENT}")
    print(f"  - Sandbag placement: {ENABLE_SANDBAG_PLACEMENT}")
    print(f"  - Dune rebuilding: {ENABLE_DUNE_REBUILDING}")
print("=" * 80 + "\n")

# =============================================================================
# SECTION 6: ELEVATION AND DUNE PROFILE FILE PATHS
# =============================================================================

print("Generating elevation profile file paths...")

ELEVATION_FILE_PATHS = []
DUNE_FILE_PATHS = []

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

    # File paths (update naming pattern if your files are named differently)
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
# SECTION 7: CASCADE SIMULATION FUNCTION
# =============================================================================

def run_cascade_simulation(
        nt, name, storm_file, alongshore_section_count, num_cores,
        beach_width_threshold, rmin, rmax, elevation_file, dune_file,
        dune_design_elevation, dune_minimum_elevation, road_ele, road_width,
        road_setback, overwash_filter, overwash_to_dune, nourishment_volume,
        background_erosion, rebuild_dune_threshold, roadway_management_on,
        beach_dune_manager_on, sea_level_rise_rate, sea_level_constant,
        sandbag_management_on, sandbag_elevation,
        enable_shoreline_offset, shoreline_offset,
):
    """
    Initialize and run the CASCADE model simulation.

    This function handles:
    1. CASCADE model initialization with all parameters
    2. Time-stepping loop (yearly updates)
    3. Management decision logic (disabled for Natural State)
    4. Result saving

    Args:
        nt (int): Number of time steps (years) to simulate
        name (str): Run name for output directory
        storm_file (str): Path to storm time series file
        alongshore_section_count (int): Total number of domains
        ... (see function signature for all parameters)

    Returns:
        cascade: Initialized and run CASCADE model object
    """

    print("=" * 80)
    print("INITIALIZING CASCADE MODEL")
    print("=" * 80)
    print(f"  Time steps: {nt}")
    print(f"  Domains: {alongshore_section_count}")
    print(f"  CPU cores: {num_cores}")
    print(f"  Sea level rise: {sea_level_rise_rate} m/yr")

    # --- Initialize CASCADE ---
    datadir = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")

    try:
        cascade = Cascade(
            datadir,
            name,
            storm_file=storm_file,
            elevation_file=elevation_file,
            dune_file=dune_file,
            parameter_file="Hatteras-CASCADE-parameters.yaml",

            # Wave parameters
            wave_height=1,
            wave_period=7,
            wave_asymmetry=0.8,
            wave_angle_high_fraction=0.2,

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

            # Management modules (all False for Natural State)
            roadway_management_module=roadway_management_on,
            beach_nourishment_module=beach_dune_manager_on,
            sandbag_management_on=sandbag_management_on,
            alongshore_transport_module=True,  # Alongshore transport always enabled
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

            # Shoreline offset (initial shoreline position)
            enable_shoreline_offset=enable_shoreline_offset,
            shoreline_offset=shoreline_offset,

            # Nourishment parameters
            nourishment_volume=nourishment_volume,
            nourishment_interval=None,
        )

        print(f"✓ CASCADE initialized successfully")
        print(f"  BermEl: {cascade.barrier3d[0].BermEl:.3f} dam ({cascade.barrier3d[0].BermEl * 10:.1f} m)")

    except Exception as e:
        print(f"❌ CASCADE initialization failed: {e}")
        raise

    # --- Time-Stepping Loop ---
    print("\n" + "=" * 80)
    print("RUNNING SIMULATION")
    print("=" * 80)

    dune_rebuild_threshold_val = rebuild_dune_threshold + (cascade.barrier3d[0].BermEl * 10)

    for time_step in range(nt - 1):
        print(f"\rYear {time_step + 1}/{nt}", end="", flush=True)

        cascade.update()

        # Check for model break condition
        if cascade.b3d_break:
            print(f"\n⚠️  Simulation stopped at year {time_step + 1}: Barrier3D break condition")
            break

        # --- Management Logic (Disabled for Natural State) ---
        # This code runs but has no effect since all management flags are False
        t = cascade.barrier3d[0].time_index
        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_nourish_now = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):
            # Skip if domain is broken or management is off
            if cascade.community_break[iB3D] or not beach_dune_manager_on[iB3D]:
                continue

            # Check beach width threshold for nourishment
            if cascade.nourishments[iB3D].beach_width[t - 1] < beach_width_threshold[iB3D]:
                tmp_nourish_now[iB3D] = 1

            # Check dune elevation threshold for rebuilding
            DuneDomainCrest = cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            DuneCrestMin = (np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl) * 10

            if DuneCrestMin < dune_rebuild_threshold_val:
                tmp_rebuild_dune[iB3D] = 1

        # Trigger management actions if thresholds met
        if np.any(tmp_nourish_now):
            cascade.nourish_now = tmp_nourish_now
        if np.any(tmp_rebuild_dune):
            cascade.rebuild_dune_now = tmp_rebuild_dune

    print("\n✓ Simulation completed")

    # --- Save Results ---
    print("\n" + "=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)

    save_path = os.path.join(OUTPUT_BASE_DIR, RUN_NAME)

    try:
        os.makedirs(save_path, exist_ok=True)
        cascade.save(save_path)
        print(f"✓ Results saved to: {save_path}")
    except Exception as e:
        print(f"❌ Error saving results: {e}")
        raise

    return cascade


# =============================================================================
# SECTION 8: MAIN EXECUTION
# =============================================================================

def main():
    """
    Main execution function - sets up all parameters and runs simulation.
    """

    # --- Define Parameter Arrays (all length = TOTAL_DOMAINS) ---

    # Dune growth parameters (parabolic curve coefficients)
    RMIN = [0.55] * TOTAL_DOMAINS
    RMAX = [0.95] * TOTAL_DOMAINS

    # Beach and dune thresholds
    BEACH_WIDTH_THRESHOLD = [30] * TOTAL_DOMAINS  # meters
    DUNE_DESIGN_ELEVATION = [DUNE_REBUILD_HEIGHT] * TOTAL_DOMAINS  # meters
    DUNE_MINIMUM_ELEVATION = [REBUILD_ELEV_THRESHOLD] * TOTAL_DOMAINS  # dam
    REBUILD_DUNE_THRESHOLD = np.full(TOTAL_DOMAINS, 1.0)  # meters

    # Road infrastructure parameters (ignored for Natural State)
    ROAD_ELEVATION = 1.45  # m above MHW
    ROAD_WIDTH = 20.0  # meters

    # Overwash parameters
    OVERWASH_FILTER = 0  # Percentage of overwash removed (0 = none)
    OVERWASH_TO_DUNE = 9  # Overwash distribution parameter

    # Nourishment parameters (ignored for Natural State)
    NOURISHMENT_VOLUME = 0  # m^3/m
    SANDBAG_ELEVATION = 0  # meters

    # Sea level rise
    SEA_LEVEL_CONSTANT = True  # Linear SLR (True) or variable (False)

    # --- Run Simulation ---
    print("STARTING CASCADE SIMULATION")
    print("=" * 80 + "\n")

    try:
        cascade = run_cascade_simulation(
            nt=RUN_YEARS,
            name=RUN_NAME,
            storm_file=STORM_FILE,
            alongshore_section_count=TOTAL_DOMAINS,
            num_cores=NUM_CORES,

            # Growth and threshold parameters
            beach_width_threshold=BEACH_WIDTH_THRESHOLD,
            rmin=RMIN,
            rmax=RMAX,

            # File paths
            elevation_file=ELEVATION_FILE_PATHS,
            dune_file=DUNE_FILE_PATHS,

            # Dune parameters
            dune_design_elevation=DUNE_DESIGN_ELEVATION,
            dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,
            rebuild_dune_threshold=REBUILD_DUNE_THRESHOLD,

            # Road parameters
            road_ele=ROAD_ELEVATION,
            road_width=ROAD_WIDTH,
            road_setback=road_setbacks_full,

            # Overwash parameters
            overwash_filter=OVERWASH_FILTER,
            overwash_to_dune=OVERWASH_TO_DUNE,

            # Nourishment parameters
            nourishment_volume=NOURISHMENT_VOLUME,

            # Erosion parameters
            background_erosion=BACKGROUND_EROSION_RATES,

            # Management flags (all False for Natural State)
            roadway_management_on=ROADWAY_MANAGEMENT_ON,
            beach_dune_manager_on=NOURISHMENT_MANAGEMENT_ON,
            sandbag_management_on=SANDBAG_MANAGEMENT_ON,
            sandbag_elevation=SANDBAG_ELEVATION,

            # Sea level rise
            sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
            sea_level_constant=SEA_LEVEL_CONSTANT,

            # Shoreline offset (initial position)
            enable_shoreline_offset=True,
            shoreline_offset=dune_offset_dam,
        )

        print("\n" + "=" * 80)
        print("✓ SIMULATION COMPLETED SUCCESSFULLY")
        print("=" * 80)
        print(f"Output saved to: {os.path.join(OUTPUT_BASE_DIR, RUN_NAME)}")

        return cascade

    except Exception as e:
        print("\n" + "=" * 80)
        print(f"❌ SIMULATION FAILED: {e}")
        print("=" * 80)
        import traceback
        traceback.print_exc()
        sys.exit(1)


# =============================================================================
# RUN SCRIPT
# =============================================================================

if __name__ == '__main__':
    main()