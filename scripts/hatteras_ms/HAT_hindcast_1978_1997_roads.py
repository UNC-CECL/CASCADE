import copy
import numpy as np
import os
import sys
from cascade.cascade import Cascade

# This script is modified from 'Natural State' to test Roadway Management.
# The ROADWAY_MANAGEMENT_ON flag is set to True, and the actual road setback data
# (road_setbacks_full) is passed to the simulation, enabling road maintenance logic.

# =============================================================================
# --- 1. PHYSICAL & DOMAIN CONSTANTS ---
# =============================================================================

# Total alongshore domains (15 Left Buffer + 105 Real + 15 Right Buffer)
TOTAL_DOMAINS = 135
NUMBER_BARRIER3D_MODELS = TOTAL_DOMAINS

# Real domain indices for data loading (105 cells)
START_REAL_DOMAIN_INDEX = 15
END_REAL_DOMAIN_INDEX = 120    # Exclusive index
START_FILE_NUMBER = 30         # Corresponds to the first domain file ID (e.g., domain_30...)

# Road Affected Indices (Used for defining the extent of the road setback array)
START_ROAD_INDEX = 29
END_ROAD_INDEX = 128
# The Road setback array must be TOTAL_DOMAINS long (135), but only the
# section START_ROAD_INDEX to END_ROAD_INDEX is typically loaded with non-zero values.

# =============================================================================
# --- 2. FILE PATH & DIRECTORY CONSTANTS ---
# =============================================================================

# Base directory for the CASCADE project and data
PROJECT_BASE_DIR = r'C:\Users\hanna\PycharmProjects\CASCADE'
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, 'data', 'hatteras_init')
OUTPUT_BASE_DIR = os.path.join(PROJECT_BASE_DIR, 'output', 'raw_runs')

# Input data file paths
DUNE_LOAD_ALL = os.path.join(HATTERAS_DATA_BASE, 'island_offset', 'Dune_Offsets_1978_1997_PADDED_135.csv')
STORM_FILE_1978_1997 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1978_1997.npy')
STORM_FILE_1997_2022 = os.path.join(HATTERAS_DATA_BASE, 'storms', 'hindcast_storms', 'storms_1997_2019.npy')
ROAD_LOAD_FILE = os.path.join(HATTERAS_DATA_BASE, 'roads', 'offset', '1978_road_setback_2row_FINAL.csv')

# =============================================================================
# --- 3. SIMULATION CONFIGURATION ---
# =============================================================================

START_YEAR = 1978  # Simulation start year (1978 or 1997)
RUN_YEARS = 20     # Total years to simulate (nt)

# Dune/Beach Management Parameters
DUNE_REBUILD_ELEVATION_M = 3.0  # (m) Design height for rebuild
REBUILD_ELEV_THRESHOLD_DAM = 0.01  # (dam) Minimum dune elevation threshold

# --- DYNAMIC CONFIGURATION (DEPENDENT ON START_YEAR) ---
if START_YEAR == 1978:
    YEAR_COLUMN_INDEX = 0  # Column in Dune_Offsets file for 1978 data
    RUN_NAME = 'HAT_1978_1997_Roadway_Management' # <<< CHANGED: Management run name
    STORM_FILE = STORM_FILE_1978_1997
    ROAD_LOAD_NAME = ROAD_LOAD_FILE

elif START_YEAR == 1997:
    YEAR_COLUMN_INDEX = 1  # Column in Dune_Offsets file for 1997 data
    RUN_NAME = 'HAT_1997_2020_Roadway_Management' # <<< CHANGED: Management run name
    STORM_FILE = STORM_FILE_1997_2022
    ROAD_LOAD_NAME = ROAD_LOAD_FILE

else:
    print(f"Error: Invalid START_YEAR {START_YEAR}. Must be 1978 or 1997.")
    sys.exit(1)

# Set the working directory to the project base
os.chdir(PROJECT_BASE_DIR)

print(f"Configuration loaded for Start Year: {START_YEAR} | Run Name: {RUN_NAME}")
print(f"Output will be saved to: {os.path.join(OUTPUT_BASE_DIR, RUN_NAME)}")

# =============================================================================
# --- 4. DATA LOADING AND PREPARATION ---
# =============================================================================

# --- A. Island Offset Data (Shoreline Position) ---
try:
    dune_offset_all = np.loadtxt(DUNE_LOAD_ALL, skiprows=1, delimiter=',')
    # Extract the column corresponding to the start year and convert m -> dam
    dune_offset_dam = dune_offset_all[:, YEAR_COLUMN_INDEX] / 10
    print(f"Shoreline offsets loaded ({dune_offset_dam.size} domains), converted to dam.")

    # Load road setbacks (required input structure, m -> dam)
    road_setbacks_raw = np.loadtxt(ROAD_LOAD_NAME, skiprows=1, delimiter=',') / 10
    print(f"Road setbacks loaded ({road_setbacks_raw.size} domains), converted to dam.")

except FileNotFoundError as e:
    print(f"CRITICAL FILE ERROR: Missing required data file: {e.filename}")
    sys.exit(1)

# --- B. Road Setback Array Initialization (135 Domains) ---
# Create the full-length array required by the CASCADE alongshore initialization.
# This array holds the initial alongshore position of the road in decameters (dam).
road_setbacks_full = [0.0] * TOTAL_DOMAINS
# Insert the loaded road setbacks into the active domain indices
road_setbacks_full[START_ROAD_INDEX:END_ROAD_INDEX] = copy.deepcopy(road_setbacks_raw)
road_setbacks_full = np.array(road_setbacks_full)
print(f"Padded road setbacks array created (length {road_setbacks_full.size}).")


# --- C. Management and Erosion Flag Arrays (135 Domains) ---
# Roadway management is enabled for the entire domain.
ROADWAY_MANAGEMENT_ON = [True] * NUMBER_BARRIER3D_MODELS # <<< CHANGED: Set to True
SANDBAG_MANAGEMENT_ON = [False] * NUMBER_BARRIER3D_MODELS
NOURISHMENT_MANAGEMENT_ON = [False] * NUMBER_BARRIER3D_MODELS
BACKGROUND_EROSION_RATES = [0.0] * TOTAL_DOMAINS  # (m/yr) Background erosion rate

# =============================================================================
# --- 5. TOPOGRAPHY AND DUNE PROFILE PATH GENERATION ---
# =============================================================================

# Arrays to store file paths for the 135 domains
ELEVATION_FILE_PATHS = []      # Topography/Backbarrier elevation paths (.npy)
DUNE_FILE_PATHS = []           # Dune profile paths (.npy)

# 1. Left Buffer Domains (15 domains: Index 0 to 14)
# Use a static sample profile for buffer zones.
for _ in range(0, START_REAL_DOMAIN_INDEX):
    dune_name = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy')
    elev_name = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy')
    DUNE_FILE_PATHS.append(dune_name)
    ELEVATION_FILE_PATHS.append(elev_name)

# 2. Real Domains (105 domains: Index 15 to 119)
# Use the unique, domain-specific resampled elevation profiles.
for i_list in range(START_REAL_DOMAIN_INDEX, END_REAL_DOMAIN_INDEX):
    # Calculate the file number ID (e.g., 30, 31, 32...)
    file_num = START_FILE_NUMBER + (i_list - START_REAL_DOMAIN_INDEX)

    dune_path_template = os.path.join(HATTERAS_DATA_BASE, 'dunes', '2009', f'domain_{file_num}_dune_2009.npy')
    DUNE_FILE_PATHS.append(dune_path_template)

    elev_path_template = os.path.join(HATTERAS_DATA_BASE, 'topography_dunes', '2009',
                                     f'domain_{file_num}_topography_2009.npy')
    ELEVATION_FILE_PATHS.append(elev_path_template)

# 3. Right Buffer Domains (15 domains: Index 120 to 134)
# Use a static sample profile for buffer zones.
for _ in range(END_REAL_DOMAIN_INDEX, TOTAL_DOMAINS):
    dune_name = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_dune.npy')
    elev_name = os.path.join(HATTERAS_DATA_BASE, 'buffer', 'sample_1_topography.npy')
    DUNE_FILE_PATHS.append(dune_name)
    ELEVATION_FILE_PATHS.append(elev_name)


# =============================================================================
# --- 6. CORE SIMULATION FUNCTION (alongshore_connected remains unchanged) ---
# =============================================================================

def alongshore_connected(
        nt, name, storm_file, alongshore_section_count, num_cores,
        beach_width_threshold, rmin, rmax, elevation_file, dune_file,
        dune_design_elevation, dune_minimum_elevation, road_ele, road_width,
        road_setback, overwash_filter, overwash_to_dune, nourishment_volume,
        background_erosion, rebuild_dune_threshold, roadway_management_on,
        beach_dune_manager_on, sea_level_rise_rate, sea_level_constant,
        trigger_dune_knockdown=False, group_roadway_abandonment=None,
        sandbag_management_on=False, sandbag_elevation=0,
        enable_shoreline_offset=True, shoreline_offset=None,
):
    """
    Initializes and runs the CASCADE alongshore-connected simulation.
    (This function remains identical to the original script)
    """
    # ------------------ INITIALIZE CASCADE OBJECT ------------------
    datadir = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="Hatteras-CASCADE-parameters.yaml", # Default parameter file
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
        background_erosion=background_erosion,
        alongshore_section_count=alongshore_section_count,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=num_cores,

        # Management Modules: Controlled by input flags
        roadway_management_module=roadway_management_on,
        beach_nourishment_module=beach_dune_manager_on,
        sandbag_management_on=sandbag_management_on,
        alongshore_transport_module=True, # Always True for alongshore simulation
        community_economics_module=False,

        # Road and Dune Parameters
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        sandbag_elevation=sandbag_elevation,

        # Overwash/Shoreline Parameters
        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,
        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
        nourishment_volume=nourishment_volume, # m^3/m
        nourishment_interval=None,
        trigger_dune_knockdown=trigger_dune_knockdown,
        group_roadway_abandonment=group_roadway_abandonment,
    )

    # ------------------ TIME-STEP LOOP ------------------
    # The management check logic remains, but only the modules that are ON
    # (i.e., roadway management) will execute their logic within Cascade.update().
    # The original management check for nourishment/dune rebuild is retained for
    # completeness if that module were to be turned on later.

    # Dune Rebuild Threshold: Calculated relative to the initial berm elevation
    dune_rebuild_elevation_threshold = rebuild_dune_threshold + (
            cascade.barrier3d[0].BermEl * 10
    )

    for time_step in range(nt - 1):
        # Progress indicator
        print("\r", "Time Step: ", time_step + 1, " / ", nt, end="")
        cascade.update() # Core physics update, includes management logic

        if cascade.b3d_break:
            print("\nSimulation halted due to internal Barrier3D break condition.")
            break

        # Management Check Logic (Only for Beach/Dune Management, which is OFF)
        t = cascade.barrier3d[0].time_index
        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_nourish_now = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):
            if cascade.community_break[iB3D] or not beach_dune_manager_on[iB3D]:
                continue # Skip if domain broken or management is explicitly off for this domain

            # Check for beach nourishment threshold
            if (
                        cascade.nourishments[iB3D].beach_width[t - 1]
                        < beach_width_threshold[iB3D]
            ):
                tmp_nourish_now[iB3D] = 1

            # Check for dune rebuilding threshold
            DuneDomainCrest = (
                cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            )
            DuneCrestMin = (
                                np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
                           ) * 10

            if DuneCrestMin < dune_rebuild_elevation_threshold:
                tmp_rebuild_dune[iB3D] = 1

        if np.any(tmp_nourish_now) == 1:
            cascade.nourish_now = tmp_nourish_now
        if np.any(tmp_rebuild_dune) == 1:
            cascade.rebuild_dune_now = tmp_rebuild_dune

    # ------------------ SAVE RESULTS ------------------
    SAVE_PATH = os.path.join(OUTPUT_BASE_DIR, RUN_NAME)
    try:
        os.makedirs(SAVE_PATH, exist_ok=True)
        print(f"\n\nOutput directory created or verified: {SAVE_PATH}")
    except OSError as e:
        print(f"\nWarning: Could not create output directory {SAVE_PATH}. Error: {e}")

    try:
        cascade.save(SAVE_PATH)
        print("CASCADE simulation results successfully saved.")
    except Exception as e:
        print(f"CRITICAL ERROR: Failed to save CASCADE results. Error: {e}")

    return cascade


# =============================================================================
# --- 7. ROADWAY MANAGEMENT CONFIGURATION ---
# =============================================================================

def alongshore_roadway_management(
        run_name, s_file, rebuild_elev_threshold_val, roadway_on_flags, road_setbacks_dam
):
    """
    Configures and initiates the 'Roadway Management' simulation.

    This function sets all necessary input arrays, enables the roadway management
    module, and passes the actual initial road position.
    """

    # --- 1. INPUT PARAMETER ARRAYS (Length 135) ---

    # Dune Growth Parameters
    RMIN = [0.55] * NUMBER_BARRIER3D_MODELS
    RMAX = [0.95] * NUMBER_BARRIER3D_MODELS

    # Dune/Beach Configuration
    BEACH_WIDTH_THRESHOLD = [30] * NUMBER_BARRIER3D_MODELS
    DUNE_DESIGN_ELEVATION = [DUNE_REBUILD_ELEVATION_M] * NUMBER_BARRIER3D_MODELS
    DUNE_MINIMUM_ELEVATION = [rebuild_elev_threshold_val] * NUMBER_BARRIER3D_MODELS
    REBUILD_DUNE_THRESHOLD = np.full(NUMBER_BARRIER3D_MODELS, 1.0) # 1.0 m (arbitrary)

    # Road Parameters (These are now active and use the loaded data)
    ROAD_ELE = 1.45  # m MHW
    ROAD_WIDTH = 20  # m
    ROAD_SETBACK = road_setbacks_dam # <<< CRITICAL CHANGE: Use the loaded array (in dam)

    # --- 2. GLOBAL PHYSICAL CONSTANTS ---
    NUM_CORES = 4
    SEA_LEVEL_RISE_RATE = 0.0056 # (m/yr)
    SEA_LEVEL_CONSTANT = True
    OVERWASH_FILTER = 0      # Overwash filtering factor
    OVERWASH_TO_DUNE = 9     # Overwash distribution parameter
    NOURISHMENT_VOLUME = 0   # m^3/m (ignored)
    SANDBAG_ELEVATION = 0    # m (ignored)

    # --- 3. SHORELINE OFFSET CONFIGURATION ---
    SHORELINE_OFFSET_ENABLED = True
    # Initial shoreline offset relative to the alongshore grid (in dam)
    SHORELINE_OFFSET_DAM = dune_offset_dam

    # Run the alongshore simulation function
    alongshore_connected(
        nt=RUN_YEARS,
        name=run_name,
        storm_file=s_file,
        alongshore_section_count=NUMBER_BARRIER3D_MODELS,
        num_cores=NUM_CORES,
        beach_width_threshold=BEACH_WIDTH_THRESHOLD,
        rmin=RMIN,
        rmax=RMAX,
        elevation_file=ELEVATION_FILE_PATHS,
        dune_file=DUNE_FILE_PATHS,
        dune_design_elevation=DUNE_DESIGN_ELEVATION,
        dune_minimum_elevation=DUNE_MINIMUM_ELEVATION,
        road_ele=ROAD_ELE,
        road_width=ROAD_WIDTH,
        road_setback=ROAD_SETBACK, # <<< Now uses the actual road position
        overwash_filter=OVERWASH_FILTER,
        overwash_to_dune=OVERWASH_TO_DUNE,
        nourishment_volume=NOURISHMENT_VOLUME,
        background_erosion=BACKGROUND_EROSION_RATES,
        rebuild_dune_threshold=REBUILD_DUNE_THRESHOLD,
        roadway_management_on=roadway_on_flags, # <<< Now uses [True] flags
        beach_dune_manager_on=NOURISHMENT_MANAGEMENT_ON,
        sea_level_rise_rate=SEA_LEVEL_RISE_RATE,
        sea_level_constant=SEA_LEVEL_CONSTANT,
        sandbag_management_on=SANDBAG_MANAGEMENT_ON,
        sandbag_elevation=SANDBAG_ELEVATION,
        shoreline_offset=SHORELINE_OFFSET_DAM,
        enable_shoreline_offset=SHORELINE_OFFSET_ENABLED,
    )


# =============================================================================
# --- 8. MAIN EXECUTION BLOCK ---
# =============================================================================

if __name__ == '__main__':
    # Execution is gated by 'if __name__ == '__main__': ' to allow clean importing
    # of functions or variables without immediate execution.
    print("\n--- Starting CASCADE Roadway Management Simulation ---")
    alongshore_roadway_management( # Renamed function call
        run_name=RUN_NAME,
        s_file=STORM_FILE,
        rebuild_elev_threshold_val=REBUILD_ELEV_THRESHOLD_DAM,
        roadway_on_flags=ROADWAY_MANAGEMENT_ON,
        road_setbacks_dam=road_setbacks_full # Pass the actual loaded road position
    )
    print("--- Simulation Complete ---")