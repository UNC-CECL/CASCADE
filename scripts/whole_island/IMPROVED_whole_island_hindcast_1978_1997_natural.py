"""
CASCADE Model - Hatteras Barrier Island Simulation
===================================================

Purpose:
    Simulate barrier island evolution for Cape Hatteras, NC using the CASCADE
    (Coastal Barrier And Shoreface Dynamics Evolution) model.

Configuration:
    - Natural State (no human management)
    - 118 real domains + 30 buffer domains = 148 total
    - Domains numbered 1-118 (south to north)
    - Time period: 1978-1997 (20 years)

Author: Hannah Henry
Date: 11/24/2025
"""

# ==============================================================================
# IMPORTS
# ==============================================================================

import numpy as np
import os
import sys
from cascade.cascade import Cascade

# ==============================================================================
# SECTION 1: USER CONFIGURATION
# ==============================================================================
# This is the main section you'll modify for different simulation scenarios

# --- Simulation Settings ---
START_YEAR = 1978  # Options: 1978 or 1997
RUN_YEARS = 20  # Number of years to simulate
NUM_CORES = 4  # Parallel processing cores (adjust for your machine)

# --- Management Scenarios ---
# Set to True/False to enable/disable different management types
ENABLE_ROADWAY_MANAGEMENT = False  # Road maintenance and relocation
ENABLE_BEACH_NOURISHMENT = False  # Beach sand addition
ENABLE_SANDBAG_MANAGEMENT = False  # Emergency sandbag placement

# --- Physical Parameters ---
SEA_LEVEL_RISE_RATE = 0.0056  # m/yr (5.6 mm/yr)
SEA_LEVEL_RISE_CONSTANT = True  # True = linear, False = accelerating
BACKGROUND_EROSION_RATE = 0.0  # m/yr uniform background erosion

# --- Dune Management Parameters (used if beach nourishment enabled) ---
DUNE_REBUILD_HEIGHT = 3.0  # m above MHW - target dune height
DUNE_REBUILD_THRESHOLD = 1.0  # m - trigger rebuilding if dune drops below
BEACH_WIDTH_THRESHOLD = 30.0  # m - trigger nourishment if beach narrows

# ==============================================================================
# SECTION 2: FILE PATHS
# ==============================================================================
# Update these paths if your data is in a different location

# --- Base Directories ---
PROJECT_BASE = r'C:\Users\hanna\PycharmProjects\CASCADE'
DATA_BASE = os.path.join(PROJECT_BASE, 'data', 'hatteras_init')
OUTPUT_BASE = os.path.join(PROJECT_BASE, 'output', 'raw_runs')

# --- Input Data Files ---
DUNE_OFFSETS_CSV = os.path.join(
    DATA_BASE, 'island_offset', '1978_whole_island',
    'Whole_Island_Dune_Offsets_1978_1997_PADDED_148.csv'
)
ROAD_SETBACKS_CSV = os.path.join(
    DATA_BASE, 'roads', 'offset', '1978_whole_island',
    'CASCADE_Whole_Island_RoadSetback_1978.csv'
)

# --- Storm Files ---
STORMS_1978_1997 = os.path.join(
    DATA_BASE, 'storms', 'hindcast_storms', 'storms_1978_1997.npy'
)
STORMS_1997_2019 = os.path.join(
    DATA_BASE, 'storms', 'hindcast_storms', 'storms_1997_2019.npy'
)

# --- Elevation Profile Directories ---
DUNE_PROFILES_DIR = os.path.join(DATA_BASE, 'dunes', '2009')
TOPO_PROFILES_DIR = os.path.join(DATA_BASE, 'topography', '2009')
BUFFER_PROFILES_DIR = os.path.join(DATA_BASE, 'buffer')

# ==============================================================================
# SECTION 3: DOMAIN CONFIGURATION
# ==============================================================================
# These define the spatial structure of your model domain

# --- Domain Counts ---
NUM_BUFFER_SOUTH = 15  # Buffer domains at south end
NUM_REAL_DOMAINS = 118  # Real island domains
NUM_BUFFER_NORTH = 15  # Buffer domains at north end
TOTAL_DOMAINS = NUM_BUFFER_SOUTH + NUM_REAL_DOMAINS + NUM_BUFFER_NORTH  # 148

# --- Domain Indices (Python 0-based indexing) ---
START_REAL_INDEX = NUM_BUFFER_SOUTH  # 15
END_REAL_INDEX = START_REAL_INDEX + NUM_REAL_DOMAINS  # 133

# --- File Numbering (domains numbered 1-118) ---
START_FILE_NUM = 1
END_FILE_NUM = NUM_REAL_DOMAINS  # 118

# --- Road Configuration ---
# Domains where roads exist (adjust based on actual road extent)
START_ROAD_INDEX = START_REAL_INDEX  # Roads start at first real domain
END_ROAD_INDEX = END_REAL_INDEX  # Roads end at last real domain

# ==============================================================================
# SECTION 4: AUTOMATIC CONFIGURATION
# ==============================================================================
# This section sets up parameters based on your choices above
# You shouldn't need to modify this section

# --- Select Data Files Based on Start Year ---
if START_YEAR == 1978:
    YEAR_COLUMN = 0  # Column in CSV for 1978 data
    RUN_NAME = 'Whole_HAT_1978_1997_Natural_State'
    STORM_FILE = STORMS_1978_1997
elif START_YEAR == 1997:
    YEAR_COLUMN = 1  # Column in CSV for 1997 data
    RUN_NAME = 'Whole_HAT_1997_2019_Natural_State'
    STORM_FILE = STORMS_1997_2019
else:
    sys.exit(f"ERROR: START_YEAR must be 1978 or 1997, not {START_YEAR}")

# --- Create Management Flag Arrays ---
ROADWAY_MGMT_ARRAY = [ENABLE_ROADWAY_MANAGEMENT] * TOTAL_DOMAINS
NOURISHMENT_MGMT_ARRAY = [ENABLE_BEACH_NOURISHMENT] * TOTAL_DOMAINS
SANDBAG_MGMT_ARRAY = [ENABLE_SANDBAG_MANAGEMENT] * TOTAL_DOMAINS

# --- Set Working Directory ---
os.chdir(PROJECT_BASE)

print("\n" + "=" * 80)
print("CASCADE SIMULATION CONFIGURATION")
print("=" * 80)
print(f"Run Name: {RUN_NAME}")
print(f"Time Period: {START_YEAR} - {START_YEAR + RUN_YEARS}")
print(f"Total Domains: {TOTAL_DOMAINS} ({NUM_REAL_DOMAINS} real + {NUM_BUFFER_SOUTH + NUM_BUFFER_NORTH} buffer)")
print(
    f"Management: {'ENABLED' if any([ENABLE_ROADWAY_MANAGEMENT, ENABLE_BEACH_NOURISHMENT, ENABLE_SANDBAG_MANAGEMENT]) else 'DISABLED (Natural State)'}")
print("=" * 80 + "\n")

# ==============================================================================
# SECTION 5: DATA LOADING
# ==============================================================================

print("Loading input data...")

try:
    # --- Load Dune Offsets (Initial Shoreline Positions) ---
    # Units: CSV is in meters, convert to decameters for CASCADE
    dune_offsets_raw = np.loadtxt(DUNE_OFFSETS_CSV, skiprows=1, delimiter=',')
    dune_offsets = dune_offsets_raw[:, YEAR_COLUMN] / 10  # Convert m → dam
    print(f"✓ Loaded dune offsets: {len(dune_offsets)} domains")
    print(f"  Range: {np.min(dune_offsets):.1f} to {np.max(dune_offsets):.1f} dam")

    # --- Load Road Setbacks (Distance from Dune to Road) ---
    # Units: CSV is in meters, convert to decameters for CASCADE
    road_setbacks_raw = np.loadtxt(ROAD_SETBACKS_CSV, skiprows=1, delimiter=',') / 10
    print(f"✓ Loaded road setbacks: {len(road_setbacks_raw)} values")
    print(f"  Range: {np.min(road_setbacks_raw):.1f} to {np.max(road_setbacks_raw):.1f} dam")

except FileNotFoundError as e:
    sys.exit(f"ERROR: Could not find data file: {e.filename}")

# --- Create Full Road Setback Array ---
# Extend to all domains (148), inserting zeros for buffer domains
road_setbacks_full = np.zeros(TOTAL_DOMAINS)
num_values = min(len(road_setbacks_raw), END_ROAD_INDEX - START_ROAD_INDEX)
road_setbacks_full[START_ROAD_INDEX:START_ROAD_INDEX + num_values] = road_setbacks_raw[:num_values]

# --- Create Background Erosion Array ---
background_erosion = np.full(TOTAL_DOMAINS, BACKGROUND_EROSION_RATE)

print("✓ Data arrays prepared\n")

# ==============================================================================
# SECTION 6: ELEVATION PROFILE FILE PATHS
# ==============================================================================

print("Generating elevation profile file paths...")

elevation_files = []
dune_files = []
files_missing = 0

# --- South Buffer Domains (use generic profile) ---
for i in range(NUM_BUFFER_SOUTH):
    elevation_files.append(os.path.join(BUFFER_PROFILES_DIR, 'sample_1_topography.npy'))
    dune_files.append(os.path.join(BUFFER_PROFILES_DIR, 'sample_1_dune.npy'))

# --- Real Domains (use domain-specific profiles) ---
for i in range(START_REAL_INDEX, END_REAL_INDEX):
    domain_num = (i - START_REAL_INDEX) + START_FILE_NUM  # Maps to 1-118

    elev_path = os.path.join(TOPO_PROFILES_DIR, f'domain_{domain_num}_topography_2009.npy')
    dune_path = os.path.join(DUNE_PROFILES_DIR, f'domain_{domain_num}_dune_2009.npy')

    # Check first domain to verify files exist
    if i == START_REAL_INDEX:
        if not os.path.exists(elev_path):
            print(f"  ⚠️  WARNING: Missing elevation file: {elev_path}")
            files_missing += 1
        if not os.path.exists(dune_path):
            print(f"  ⚠️  WARNING: Missing dune file: {dune_path}")
            files_missing += 1

    elevation_files.append(elev_path)
    dune_files.append(dune_path)

# --- North Buffer Domains (use generic profile) ---
for i in range(NUM_BUFFER_NORTH):
    elevation_files.append(os.path.join(BUFFER_PROFILES_DIR, 'sample_1_topography.npy'))
    dune_files.append(os.path.join(BUFFER_PROFILES_DIR, 'sample_1_dune.npy'))

print(f"✓ Generated {len(elevation_files)} elevation file paths")
print(f"✓ Generated {len(dune_files)} dune file paths")
if files_missing > 0:
    print(f"  ⚠️  WARNING: {files_missing} files not found!")
print()


# ==============================================================================
# SECTION 7: CASCADE MODEL EXECUTION
# ==============================================================================

def run_cascade_simulation():
    """
    Initialize and run the CASCADE model with configured parameters.

    Returns:
        cascade: Completed CASCADE model object with simulation results
    """

    # --- Create Parameter Arrays (all must be length TOTAL_DOMAINS) ---
    params = {
        # Dune growth rates (parabolic growth curve parameters)
        'rmin': np.full(TOTAL_DOMAINS, 0.55),
        'rmax': np.full(TOTAL_DOMAINS, 0.95),

        # Beach/dune management thresholds
        'beach_width_threshold': np.full(TOTAL_DOMAINS, BEACH_WIDTH_THRESHOLD),
        'dune_design_elevation': np.full(TOTAL_DOMAINS, DUNE_REBUILD_HEIGHT),
        'dune_minimum_elevation': np.full(TOTAL_DOMAINS, 0.01),  # dam
        'rebuild_dune_threshold': np.full(TOTAL_DOMAINS, DUNE_REBUILD_THRESHOLD),

        # Road parameters (ignored in natural state, but required)
        'road_elevation': 1.45,  # m MHW
        'road_width': 20,  # m
        'road_setback': road_setbacks_full,

        # Overwash parameters
        'overwash_filter': 0,  # Fraction of overwash removed (0 = none)
        'overwash_to_dune': 9,  # Overwash distribution parameter

        # Nourishment parameters
        'nourishment_volume': 0,  # m³/m (0 = none)
        'sandbag_elevation': 0,  # m (0 = none)
    }

    print("=" * 80)
    print("INITIALIZING CASCADE MODEL")
    print("=" * 80)
    print(f"  Time steps: {RUN_YEARS}")
    print(f"  Domains: {TOTAL_DOMAINS}")
    print(f"  CPU cores: {NUM_CORES}")
    print(f"  Sea level rise: {SEA_LEVEL_RISE_RATE} m/yr")

    # --- Initialize CASCADE ---
    cascade = Cascade(
        datadir=os.path.join(PROJECT_BASE, 'data', 'hatteras_init'),
        name=RUN_NAME,

        # Input files
        storm_file=STORM_FILE,
        elevation_file=elevation_files,
        dune_file=dune_files,
        parameter_file="Hatteras-CASCADE-parameters.yaml",

        # Wave climate
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,

        # Sea level
        sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
        sea_level_rise_constant=SEA_LEVEL_RISE_CONSTANT,

        # Background processes
        background_erosion=background_erosion,

        # Domain configuration
        alongshore_section_count=TOTAL_DOMAINS,
        time_step_count=RUN_YEARS,
        num_cores=NUM_CORES,

        # Dune dynamics
        min_dune_growth_rate=params['rmin'],
        max_dune_growth_rate=params['rmax'],

        # Management modules
        roadway_management_module=ROADWAY_MGMT_ARRAY,
        beach_nourishment_module=NOURISHMENT_MGMT_ARRAY,
        sandbag_management_on=SANDBAG_MGMT_ARRAY,
        alongshore_transport_module=True,
        community_economics_module=False,

        # Road configuration
        road_ele=params['road_elevation'],
        road_width=params['road_width'],
        road_setback=params['road_setback'],

        # Dune configuration
        dune_design_elevation=params['dune_design_elevation'],
        dune_minimum_elevation=params['dune_minimum_elevation'],
        sandbag_elevation=params['sandbag_elevation'],

        # Overwash configuration
        overwash_filter=params['overwash_filter'],
        overwash_to_dune=params['overwash_to_dune'],

        # Shoreline offset (initial position adjustment)
        enable_shoreline_offset=True,
        shoreline_offset=dune_offsets,

        # Nourishment configuration
        nourishment_volume=params['nourishment_volume'],
        nourishment_interval=None,
        trigger_dune_knockdown=False,
        group_roadway_abandonment=None,
    )

    print(f"✓ CASCADE initialized successfully")
    print(f"  Initial berm elevation: {cascade.barrier3d[0].BermEl:.3f} dam")

    # --- Run Time Loop ---
    print("\n" + "=" * 80)
    print("RUNNING SIMULATION")
    print("=" * 80)

    for t in range(RUN_YEARS - 1):
        print(f"\r  Time step: {t + 1:3d} / {RUN_YEARS}", end="")

        cascade.update()

        if cascade.b3d_break:
            print(f"\n⚠️  Simulation stopped at year {t + 1} - barrier height drowned")
            break

        # --- Management Decision Logic ---
        # (Only active if management is enabled)
        if ENABLE_BEACH_NOURISHMENT or ENABLE_ROADWAY_MANAGEMENT:
            time_index = cascade.barrier3d[0].time_index
            trigger_nourishment = np.zeros(TOTAL_DOMAINS)
            trigger_rebuild = np.zeros(TOTAL_DOMAINS)

            dune_threshold = params['rebuild_dune_threshold'][0] + (cascade.barrier3d[0].BermEl * 10)

            for i in range(TOTAL_DOMAINS):
                if cascade.community_break[i]:
                    continue

                # Check beach width
                if NOURISHMENT_MGMT_ARRAY[i]:
                    beach_width = cascade.nourishments[i].beach_width[time_index - 1]
                    if beach_width < params['beach_width_threshold'][i]:
                        trigger_nourishment[i] = 1

                # Check dune height
                if NOURISHMENT_MGMT_ARRAY[i]:
                    dune_crest = cascade.barrier3d[i].DuneDomain[time_index - 1, :, :].max(axis=1)
                    min_crest = (np.min(dune_crest) + cascade.barrier3d[i].BermEl) * 10
                    if min_crest < dune_threshold:
                        trigger_rebuild[i] = 1

            if np.any(trigger_nourishment):
                cascade.nourish_now = trigger_nourishment
            if np.any(trigger_rebuild):
                cascade.rebuild_dune_now = trigger_rebuild

    print("\n✓ Simulation completed successfully\n")
    return cascade


def save_results(cascade):
    """
    Save CASCADE simulation results to disk.

    Args:
        cascade: Completed CASCADE model object
    """
    print("=" * 80)
    print("SAVING RESULTS")
    print("=" * 80)

    save_path = os.path.join(OUTPUT_BASE, RUN_NAME)

    try:
        os.makedirs(save_path, exist_ok=True)
        print(f"✓ Output directory: {save_path}")

        cascade.save(save_path)
        print(f"✓ Results saved successfully")

    except Exception as e:
        print(f"❌ ERROR saving results: {e}")
        raise


# ==============================================================================
# SECTION 8: MAIN EXECUTION
# ==============================================================================

if __name__ == '__main__':
    """
    Main execution block - runs when script is executed directly.
    """

    try:
        # Run the simulation
        cascade = run_cascade_simulation()

        # Save results
        save_results(cascade)

        print("\n" + "=" * 80)
        print("SIMULATION COMPLETE")
        print("=" * 80)
        print(f"Results saved to: {os.path.join(OUTPUT_BASE, RUN_NAME)}")
        print("=" * 80 + "\n")

    except KeyboardInterrupt:
        print("\n\n⚠️  Simulation interrupted by user")
        sys.exit(1)

    except Exception as e:
        print(f"\n\n❌ SIMULATION FAILED: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)