# Define model simulations to run for modeling Hatteras
import copy
import pandas as pd
import numpy as np
import os
import sys  # Added for error exiting
from cascade.cascade import Cascade

# --- DEFINITIVE MODEL CONSTANTS ---

# Total number of domains: 15 Left Buffer + 105 Real + 15 Right Buffer
TOTAL_DOMAINS = 135

# --- REAL DOMAIN INDICES (105 Cells, Model ID 16 to 120) ---
# The first Real Domain (Model ID 16) is at Python Index 15.
START_REAL_DOMAIN_INDEX = 15
# The slice is exclusive, so index 120 is needed to include Model ID 120 (Python Index 119).
END_REAL_DOMAIN_INDEX = 120
# Total Real Domains = END_REAL_DOMAIN_INDEX - START_REAL_DOMAIN_INDEX = 105

# --- ROAD AFFECTED INDICES (99 Cells, Model ID 30 to 128) ---
# The road's effect starts at Model ID 30 (Python Index 29).
START_ROAD_INDEX = 29
# The road's effect ends at Model ID 128 (Python Index 127).
END_ROAD_INDEX = 128
# Total Road Affected Cells = END_ROAD_INDEX - START_ROAD_INDEX = 99
# ----------------------------------

# Define the starting file number for the real domains (domain_30.npy, domain_31.npy, etc.)
START_FILE_NUMBER = 30


# Set data paths
# NOTE: Keeping the os.chdir here as it was in your original script
os.chdir(r'C:\Users\hanna\PycharmProjects\CASCADE')
# Specify variables to use in calling function

buffer_enabled = True
island_grid_number = 105  # 105 real domains from GIS

# --- FILE PATH CONSTANTS ---
# Dune Offsets (Combined 1978 and 1997 offsets in one file)
DUNE_LOAD_ALL = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\Dune_Offsets_1978_1997_PADDED_135.csv'

# Storm Files (Separate for each period)
STORM_FILE_1978_1997 = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\storms_1978_1997.npy'
STORM_FILE_1997_2022 = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\storms_1997_2019.npy'

# Road Setback Files (Separate for each period)
ROAD_LOAD_1978 = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_2row_FINAL.csv'
ROAD_LOAD_1997 = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\roads\offset\1978_road_setback_2row_FINAL.csv' #Need to update with 1997 values

# --- CONFIGURATION (Set hindcast start year here) ---
start_year = 1978

# --- DYNAMIC VARIABLE ASSIGNMENT ---

# dune_load_name is now always the combined file
dune_load_name = DUNE_LOAD_ALL # This is NOT currently in DAM...

if start_year == 1978:
    # 1978 is the first column (index 0) in Dune_Offsets_1978_1997_CASCADE_Input.csv
    year_column_index = 0
    run_name = 'HAT_1978_1997_Hindcast'
    s_file = STORM_FILE_1978_1997
    road_load_name = ROAD_LOAD_1978

elif start_year == 1997:
    # 1997 is the second column (index 1) in Dune_Offsets_1978_1997_CASCADE_Input.csv
    year_column_index = 1
    # NOTE: Updated run_name to reflect Hatteras, as this script uses Hatteras data paths
    run_name = 'HAT_1997_2020_Hindcast_Final_3'
    s_file = STORM_FILE_1997_2022
    road_load_name = ROAD_LOAD_1997

else:
    print(f"Error: Invalid start_year {start_year}. Must be 1978 or 1997.")
    sys.exit(1)


# --- DATA LOADING ---

# Load the combined dune offsets file
# Assumes dune offset file has shape (N_domains, 2)
dune_offset_all = np.loadtxt(dune_load_name, skiprows=1, delimiter=',')

# Select the correct column for the initial dune offset based on the start year
# This is the single-column array used by CASCADE
dune_offset = dune_offset_all[:, year_column_index] / 10  # Convert m → dam
dune_offset_c = copy.deepcopy(dune_offset)

print(f"Dune offsets converted: {np.min(dune_offset):.1f} to {np.max(dune_offset):.1f} dam")

# Load the road setbacks for the chosen year
road_setbacks = np.loadtxt(road_load_name, skiprows=1, delimiter=',') / 10  # Convert m → dam; testing Benton originally had *10
print(f"Road offsets converted: {np.min(road_setbacks):.1f} to {np.max(road_setbacks):.1f} dam")

print(f"Configuration loaded for start_year: {start_year}")
print(f"Road Setbacks loaded from: {road_load_name}")
print(f"Storm File loaded: {s_file}")

rebuild_elev_threshold = 0.01
Dune_Rebuilding_Height = 3


# --- ROAD SETBACKS: EXTENDED TO 135 DOMAINS, APPLIED TO ROAD RANGE (30-128) ---
# WARNING: Check your units! If road_setbacks is loaded in meters (m), multiplying by 10
# will make them decameters (dam) or 10x too large if CASCADE expects meters (m).
# If the CSV is in meters, remove the '*10'.
# We will keep it as you wrote it, assuming it's intended for a unit correction.
road_setbacks = road_setbacks/10 # *10 ---> This is from Bentons code, may need to convert to DAM instead of multiply?

# Initialize an array of size 135 (TOTAL_DOMAINS) with 0s
r_s = [0]*TOTAL_DOMAINS
# Apply the loaded road setbacks (99 values) to the active road range [29:128]
# Python indices 29 through 127 correspond to Model IDs 30 through 128 (99 domains)
r_s[START_ROAD_INDEX:END_ROAD_INDEX] = copy.deepcopy(road_setbacks)
road_setbacks = r_s
# -------------------------------------------------------------------------------


# Define which B3D island models have management features enabled
# --- ROAD AND SANDBAG CELLS: EXTENDED TO 135 DOMAINS, APPLIED TO ROAD RANGE (30-128) ---
road_cells = [False] * TOTAL_DOMAINS
# 99 cells (index 29 to 127)
road_cells[START_ROAD_INDEX:END_ROAD_INDEX] = [True] * (END_ROAD_INDEX - START_ROAD_INDEX)

sandbag_cells = [False] * TOTAL_DOMAINS
# 99 cells (index 29 to 127)
sandbag_cells[START_ROAD_INDEX:END_ROAD_INDEX] = [True] * (END_ROAD_INDEX - START_ROAD_INDEX)
# ---------------------------------------------------------------------------------------


run_years = 20

# --- BACKGROUND erosion rates: 135 DOMAINS ---
background_threhold_list = [0.0] * TOTAL_DOMAINS

# -------------------------------------------------------------------

# --- 4. INPUT DATA: VERTICAL ELEVATION PROFILE FILE PATHS (.npy files) ---

# Lists to store the file paths for dune and topography profiles for all 135 domains
e_file = []  # Topography/Backbarrier elevation paths
d_file = []  # Dune profile paths

# Base directory for the Hatteras data
HATTERAS_BASE = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init'

# 1. Left Buffer Domains (15 domains: Index 0 to 14)
# These domains repeat the profile from Sample_1
for i_list in range(0, START_REAL_DOMAIN_INDEX):
    dune_name = os.path.join(HATTERAS_BASE, 'buffer', 'sample_1_dune.npy')
    elev_name = os.path.join(HATTERAS_BASE, 'buffer', 'sample_1_topography.npy')
    d_file.append(dune_name)
    e_file.append(elev_name)

# 2. Real Domains (105 domains: Index 15 to 119)
# CRITICAL FIX 1: The 'file_num' calculation now uses the defined START_FILE_NUMBER (30).
for i_list in range(START_REAL_DOMAIN_INDEX, END_REAL_DOMAIN_INDEX):
    # Calculate file_num: starts at 30 (for list index 15) and ends at 134 (for list index 119)
    file_num = START_FILE_NUMBER + (i_list - START_REAL_DOMAIN_INDEX)

    # Dune Profile (2009 profiles for the dune structure)
    dune_path_template = os.path.join(HATTERAS_BASE, 'dunes', '2009', f'domain_{file_num}_resampled_dune_2009.npy')
    d_file.append(dune_path_template)

    # Topography Profile (2009 profiles for the backbarrier)
    elev_path_template = os.path.join(HATTERAS_BASE, 'topography', '2009',
                                      f'domain_{file_num}_resampled_topography_2009.npy')
    e_file.append(elev_path_template)

# 3. Right Buffer Domains (15 domains: Index 120 to 134)
# These domains also repeat the profile from Sample_1
for i_list in range(END_REAL_DOMAIN_INDEX, TOTAL_DOMAINS):
    dune_name = os.path.join(HATTERAS_BASE, 'buffer', 'sample_1_dune.npy')
    elev_name = os.path.join(HATTERAS_BASE, 'buffer', 'sample_1_topography.npy')
    d_file.append(dune_name)
    e_file.append(elev_name)


def alongshore_connected(
    nt,
    name,
    storm_file,
    alongshore_section_count,
    num_cores,
    beach_width_threshold,  # not a parameter in cascade, for triggering: must be list
    rmin,  # the remaining variables are arrays
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
    sea_level_rise_rate=0.0037,  # not an array
    sea_level_constant=True,  # not an array
    trigger_dune_knockdown=False,
    group_roadway_abandonment=None,
    sandbag_management_on=False,
    sandbag_elevation=5,
    enable_shoreline_offset=True,
    shoreline_offset = [0] ,
):
    # ###############################################################################
    # 9 - connect cascade domains (human management) with AST
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "./data/hatteras_init/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="Hatteras-CASCADE-parameters.yaml",
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
        roadway_management_module=roadway_management_on,
        alongshore_transport_module=True,  # couple brie
        beach_nourishment_module=beach_dune_manager_on,
        community_economics_module=False,  # no community dynamics
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        trigger_dune_knockdown=trigger_dune_knockdown,
        group_roadway_abandonment=group_roadway_abandonment,
        nourishment_interval=None,  # yrs
        nourishment_volume=nourishment_volume,  # m^3/m
        overwash_filter=overwash_filter,  # % overwash removed
        overwash_to_dune=overwash_to_dune,
        sandbag_management_on = sandbag_management_on,
        sandbag_elevation = sandbag_elevation,
        enable_shoreline_offset=enable_shoreline_offset,
        shoreline_offset=shoreline_offset,
    )

    # --------- LOOP ---------

    # after each year, check the beach width and dune elevation and decide if you want to nourish or rebuild the dune
    # next year with nourish_now parameter; just use first B3D domain, since all berm elevations are equivalent
    dune_rebuild_threshold = rebuild_dune_threshold + (
        cascade.barrier3d[0].BermEl * 10
    )  # if rebuild_dune_threshold=0.3, this is the same threshold for abs. min elevation as in RoadwayManager (m MHW)

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        if cascade.b3d_break:
            break

        t = cascade.barrier3d[0].time_index
        tmp_rebuild_dune = np.zeros(alongshore_section_count)
        tmp_nourish_now = np.zeros(alongshore_section_count)

        for iB3D in range(alongshore_section_count):
            # don't do any beach/dune management activities if the barrier has become too narrow to sustain a community
            if cascade.community_break[iB3D]:
                pass
            # and only manage beach/dune if it is turned on
            elif beach_dune_manager_on[iB3D]:
                if (
                    cascade.nourishments[iB3D].beach_width[t - 1]
                    < beach_width_threshold[iB3D]
                ):
                    # cascade.nourish_now[iB3D] = 1
                    tmp_nourish_now[iB3D] = 1

                DuneDomainCrest = (
                    cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
                )  # Maximum height of each row in dune domain [dam]
                DuneCrestMin = (
                    np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
                ) * 10  # m MHW

                if DuneCrestMin < dune_rebuild_threshold:
                    # cascade.rebuild_dune_now[iB3D] = 1
                    tmp_rebuild_dune[iB3D] = 1

        # only nourish or rebuild dune if all segments fall below threshold (more realistic)
        # NOTE: This logic ensures ALL managed segments must fail before any management occurs.
        if np.all(tmp_nourish_now[beach_dune_manager_on]) == 1:
            cascade.nourish_now = tmp_nourish_now
        if np.all(tmp_rebuild_dune[beach_dune_manager_on]) == 1:
            cascade.rebuild_dune_now = tmp_rebuild_dune

    # --------- SAVE ---------
    # Define the output directory relative to the current script's location
    # Using os.path.join ensures paths work correctly on Windows, Mac, and Linux.
    SAVE_DIRECTORY = r'C:\Users\hanna\PycharmProjects\CASCADE\output\Hindcast_1978_1997_buffers'

    # --- 1. Create the Directory (Robustly) ---
    # Ensure the directory exists before saving.
    # exist_ok=True prevents an error if the directory is already there.
    try:
        os.makedirs(SAVE_DIRECTORY, exist_ok=True)
        print(f"Output directory created or verified: {SAVE_DIRECTORY}")
    except OSError as e:
        # This is a fallback in case permissions fail or other weird OS errors occur
        print(f"Warning: Could not create output directory {SAVE_DIRECTORY}. Error: {e}")

    # --- 2. Save the Results ---
    # If cascade.save() expects a directory path:
    try:
        cascade.save(SAVE_DIRECTORY)
        print("CASCADE simulation results successfully saved.")
    except Exception as e:
        print(f"CRITICAL ERROR: Failed to save CASCADE results. Error: {e}")

    # --- 3. Return
    return cascade


def alongshore_uniform(run_name, s_file, rebuild_elev_threshold):
    # variables that DO NOT change among runs
    number_barrier3d_models = TOTAL_DOMAINS # Corrected to 135
    beach_width_threshold = [30] * number_barrier3d_models
    rmin = [0.55] * number_barrier3d_models
    rmax = [0.95] * number_barrier3d_models
    elevation_file = e_file
    dune_file = d_file
    storm_file = s_file
    dune_design_elevation = [Dune_Rebuilding_Height] * number_barrier3d_models
    num_cores = 4  # for my laptop, max is 15
    dune_minimum_elevation = rebuild_elev_threshold
    road_ele = 1.45  # m MHW (Need to update for Hatteras, make it per domain)
    road_width = 20  # m (2 lane road on Hatteras)
    road_setback = road_setbacks  # m (now 135 elements)
    overwash_filter = 0
    overwash_to_dune = 9
    nourishment_volume = 100  # m^3/m
    background_erosion = background_threhold_list # m/yr, background shoreline erosion (135 elements)
    rebuild_dune_threshold = 1  # m
    sandbag_management_on = sandbag_cells # (now 135 elements, True for 99 road domains)
    sandbag_elevation = 4 # m

    # baseline models for comparison -- all roadways ----------------------------------------
    roads_on = road_cells # (now 135 elements, True for 99 road domains)
    nourishments_on = [False] * number_barrier3d_models
    sea_level_rise_rate = 0.0056
    sea_level_constant = True  # linear

    # Island offsets
    shoreline_offset_enabled = True
    shoreline_offset = dune_offset

    alongshore_connected(
        nt=run_years,
        name=run_name,
        storm_file=storm_file,
        alongshore_section_count=number_barrier3d_models,
        num_cores=num_cores,
        beach_width_threshold=beach_width_threshold,
        rmin=rmin,
        rmax=rmax,  # rave = 0.75
        elevation_file=elevation_file,
        dune_file=dune_file,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,
        nourishment_volume=nourishment_volume,
        background_erosion=background_erosion,
        rebuild_dune_threshold=rebuild_dune_threshold,
        roadway_management_on=roads_on,
        beach_dune_manager_on=nourishments_on,
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_constant=sea_level_constant,
        sandbag_management_on=sandbag_management_on,
        sandbag_elevation=sandbag_elevation,
        shoreline_offset=shoreline_offset,
        enable_shoreline_offset=shoreline_offset_enabled,
    )


alongshore_uniform(run_name=run_name,s_file=s_file,rebuild_elev_threshold = rebuild_elev_threshold)