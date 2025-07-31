# ==============================================================================
# HINDCAST SIMULATION SCRIPT – HATTERAS ISLAND (1978–1997)
# ==============================================================================

"""
This script runs a historical hindcast simulation of barrier island evolution
on Hatteras Island, North Carolina, using the CASCADE (Coastal cAStal Community-lAnDscape Evolution) model.

The simulation spans from 1978 to 1997 and incorporates site-specific storm,
topography, and management data to evaluate geomorphic change over time.

This setup is adapted from an Ocracoke Island configuration by Benton Franklin,
and modified for Hatteras Island using updated input layers and domain parameters.

Key outputs include shoreline position, dune evolution, road impacts, and
erosion under historical forcing. This script is intended to support reproducible
barrier island modeling and scenario testing.

Author: Hannah Henry
Affiliation: University of North Carolina at Chapel Hill
Date Created: 2024
Last Modified: 2025
"""

# ------------------------------
# Import required libraries
# ------------------------------
import os
import copy
import numpy as np
from cascade.cascade import Cascade  # CASCADE core model class

# ------------------------------
# Set working directory for local execution (adjust if needed)
# ------------------------------
os.chdir('C:\\Users\\hanna\\PycharmProjects\\CASCADE')

# ------------------------------
# Define simulation name and core configuration files
# ------------------------------
run_name = 'Hindcast_1978_1997'

# Path to storm event input file (NumPy array of historical storm data)
s_file = 'C:\\Users\\hanna\\PycharmProjects\\CASCADE\\data\\hatteras_init\\storms\\hindcast_storms\\HAT_1978_2022_Final_Hindcast_Storms.npy'

# ------------------------------
# Define island configuration and input file paths
# ------------------------------

buffer_enabled = True  # Enables domain buffers at both ends for alongshore transport

island_grid_number = 92        # Number of actual model domains for Hatteras (not including buffers)
Total_B3D_Number = 70          # Total number of Barrier3D segments including buffers [⚠️ TODO: verify this value]

# Paths to road and dune offset CSVs
road_load_name = 'C:\\Users\\frank\\OneDrive - University of North Carolina at Chapel Hill\\Chapter 3\\Revised_Offshore_Datum\\Corrected_Road_Offsets.csv'  # ⚠️ TODO: Update path
dune_load_name = 'C:\\Users\\hanna\\PycharmProjects\\CASCADE\\data\\hatteras_init\\dunes\\1978_relative_offsets.csv'  # ✅ Already updated

# Load road and dune offset arrays
road_setbacks = np.loadtxt(road_load_name, skiprows=1, delimiter=',')  # ⚠️ Path needs to point to Hatteras road data
dune_offset = np.loadtxt(dune_load_name, skiprows=1, delimiter=',')    # Each column corresponds to a year

# Create a deep copy to preserve original dune_offset
dune_offset_c = copy.deepcopy(dune_offset)

# ------------------------------
# Set starting year and select corresponding dune offset column
# ------------------------------

start_year = 1997  # Year corresponding to the dune offset column you want to use

# Determine which column of dune_offset to use based on the start year
if start_year == 1974:
    year = 0
elif start_year == 1988:
    year = 1
else:
    year = 2  # Default fallback for 1997 or other year

# ------------------------------
# Set number of simulation years
# ------------------------------
run_years = 23  # Number of time steps (years) for this hindcast run


# ------------------------------
# Final background shoreline change (m/yr) per real domain (length = 69)
# Used to simulate long-term erosion/accretion outside storm events
# ------------------------------
background_threhold_list = [
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    60, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    -14, 0, 0, 0, 0,
    -10, 0, 0, 0, 0,
    0, 0, 0, -30,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0
]

# Alternate versions of background shoreline erosion (m/yr) per domain
# Each list has 69 values (one per real B3D domain), describing localized shoreline change

# background_threhold_list = [15.0,55.0,20,5,-2, ...]
# background_threhold_list = [100,85.0,15,5,0, ...]
# background_threhold_list = [5,60,30,-5,-5, ...]


# ------------------------------
# Process Road and Dune Offset Inputs
# ------------------------------

# Convert road and dune offset values from meters to centimeters (CASCADE expects cm)
# Multiply by 10 to convert from meters to centimeters
road_setbacks = road_setbacks[:, 0] * 10                 # Use first column of road setbacks
dune_offset = dune_offset[:, year] * 10                  # Use selected year column for dune offset

# ----------------------------------------
# Configure Road Setbacks for Model Domains
# ----------------------------------------

# Initialize road setback list with 0s for all domains (including buffers)
# Only real domains should have actual road data assigned
r_s = [0] * Total_B3D_Number

# Assign road setback values to central domains (15–54), excluding 15 buffers on each end
# ✅ These indices assume 15 buffer domains on each side and 40 real domains with roads.
# ⚠️ TODO: If buffer size or domain count changes, update this range accordingly.
r_s[15:55] = copy.deepcopy(road_setbacks)  # Assign road data to middle 40 domains

# Replace road_setbacks list with final version for model input
road_setbacks = r_s

# ----------------------------------------
# Identify Which Domains Contain a Road
# ----------------------------------------

# Create Boolean list of road presence per domain (True = road present)
road_cells = [False] * Total_B3D_Number

# ✅ These indices assume 15 buffers on each side and 40 road-containing domains
# ⚠️ TODO: Confirm this matches your data. Adjust if domain count changes.
road_cells[15:55] = [True] * 40  # Mark True where road exists

# ✅ TODO REMINDERS:
# - The indices 15:55 assume 15 buffer domains on each side and 40 real domains with roads.
# - Total_B3D_Number must equal 15 + 92 + 15 = 122 if you’re using 92 real domains and 15 buffers each side.
# - Make sure the length of road_setbacks[:, 0] is 40 to avoid index mismatch.

# ------------------------------
# Load Dune and Elevation Arrays for Each Domain
# ------------------------------

# Initialize empty lists to store file paths for dune and elevation arrays
e_file = []  # Elevation files
d_file = []  # Dune files

# ----------------------------------------
# Add Left Buffer Domains (15 domains)
# ----------------------------------------
# These are placeholder domains at the updrift end of the island to allow
# alongshore sediment transport to stabilize before entering real domains.

# ⚠️ TODO: Replace with actual Hatteras buffer dune/elevation file paths if available
for i in range(0, 15):
    dune_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\dunes\\Sample_1_dune.npy'
    elev_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\elevations\\Sample_1_topography.npy'
    d_file.append(dune_name)
    e_file.append(elev_name)

# ----------------------------------------
# Add Real Model Domains (92 domains)
# ----------------------------------------
# These are the actual model domains representing Hatteras Island

# ⚠️ TODO: Replace range(11, 50) with correct index range for your 92 real domains
# ⚠️ TODO: Update file paths to point to your Hatteras elevation/dune .npy files
for i in range(11, 50):  # ⚠️ Placeholder range – must be updated
    dune_name = f'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\dunes\\Sample_{i}_dune.npy'
    elev_name = f'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\elevations\\Sample_{i}_topography.npy'
    d_file.append(dune_name)
    e_file.append(elev_name)

# ----------------------------------------
# Add Right Buffer Domains (15 domains)
# ----------------------------------------
# These are placeholder domains at the downdrift end of the island to support
# sediment transport beyond the modeled region.

# ⚠️ TODO: Replace with actual Hatteras buffer dune/elevation file paths if available
for i in range(0, 15):
    dune_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\dunes\\Sample_1_dune.npy'
    elev_name = 'C:\\Users\\frank\\PycharmProjects\\CASCADE\\data\\Ocracoke_init_data\\elevations\\Sample_1_topography.npy'
    d_file.append(dune_name)
    e_file.append(elev_name)

# ✅ TODO REMINDERS:
# - Replace hardcoded file paths with your own Hatteras file structure (use `data/hatteras_init/...`)
# - Replace range(11, 50) with the correct loop to load your 92 real domain files
# - Consider using relative paths or path-building logic with `os.path.join()` for portability

# ------------------------------
# Define CASCADE Simulation Function for Alongshore-Connected Domains
# ------------------------------

def alongshore_connected(
    nt,                                 # Number of timesteps (e.g., years)
    name,                               # Run name (used for output folders/naming)
    storm_file,                         # Path to .npy file containing storm time series
    alongshore_section_count,          # Number of alongshore sections (e.g., 92 real domains)
    num_cores,                          # Number of cores to use (for parallelism, if implemented)

    # --- Site-specific coastal management thresholds ---
    beach_width_threshold,             # List of beach width thresholds triggering management actions
    rmin,                              # Array of minimum road elevations (e.g., from DEM)
    rmax,                              # Array of maximum road elevations (if used)

    # --- Elevation and dune topography ---
    elevation_file,                    # List of .npy files (one per domain) with topography
    dune_file,                         # List of .npy files with dune positions
    dune_design_elevation,             # Array of target dune crest elevations
    dune_minimum_elevation,           # Array of minimum allowed dune elevations

    # --- Roadway features and management ---
    road_ele,                          # Array of mean road elevations per domain
    road_width,                        # Constant or array: width of road (in meters)
    road_setback,                      # Array of road setback distances
    overwash_filter,                   # Array for overwash event filtering (if implemented)
    overwash_to_dune,                  # Array of overwash redistribution to dunes (boolean/float)
    nourishment_volume,                # Array: volume of nourishment added (m³/m)

    # --- Background erosion/accretion ---
    background_erosion,               # Array of background shoreline change (m/yr)

    # --- Management triggers and toggles ---
    rebuild_dune_threshold,           # Threshold for rebuilding dunes (e.g., % loss)
    roadway_management_on,            # Boolean: turn on/off road management logic
    beach_dune_manager_on,            # Boolean: enable beach/dune maintenance module

    # --- Sea level rise settings ---
    sea_level_rise_rate=0.0037,       # Sea level rise per year (m)
    sea_level_constant=True,          # Toggle for constant vs. dynamic sea level rise

    # --- Optional/experimental toggles ---
    trigger_dune_knockdown=False,     # Simulate dune breaching or flattening (if used)
    group_roadway_abandonment=None,   # Used for spatially grouped road retreat (optional)
    sandbag_management_on=False,      # Toggle for simulating sandbag policy
    sandbag_elevation=5,              # Threshold elevation for sandbag deployment
    enable_shoreline_offset=False,    # Toggle: use shoreline offset feature (if used)
    shoreline_offset=[0],             # Array of shoreline offsets per domain (optional)
):
    """
    Core simulation wrapper to run a CASCADE hindcast or forecast across multiple
    alongshore-connected model segments. This allows for detailed spatial and temporal
    simulation of barrier island evolution, integrating management and storm responses.

    ⚠️ You do not need to edit this function definition unless you're modifying
    CASCADE’s architecture. Inputs are passed in from your script below.

    """

    # (Function implementation goes here — in your case, likely a call to `Cascade.run()`)

    pass  # Placeholder if implementation is omitted here

    # ==============================================================================
    # Initialize CASCADE model with site-specific parameters and inputs
    # ==============================================================================

    # ------------------------------
    # Verify and set the model data directory
    # ------------------------------
    datadir = "./data/hatteras_init/"  # ✅ Make sure this points to your main input folder

    # ------------------------------
    # Create the CASCADE model instance
    # ------------------------------
    cascade = Cascade(
        datadir=datadir,  # Directory containing site-specific input data
        name=name,  # Simulation run name (used for output labeling)

        # --- Core storm and terrain inputs ---
        storm_file=storm_file,  # Path to .npy file with storm time series
        elevation_file=elevation_file,  # List of .npy files with domain elevation arrays
        dune_file=dune_file,  # List of .npy files with domain dune crest arrays

        # --- Simulation parameters ---
        parameter_file="Hatteras-CASCADE-parameters.yaml",  # YAML config file for physical settings
        wave_height=1,  # Mean offshore wave height (m)
        wave_period=7,  # Mean offshore wave period (s)
        wave_asymmetry=0.8,  # Asymmetry of wave direction distribution
        wave_angle_high_fraction=0.2,  # Fraction of waves from oblique angles

        # --- Sea level rise settings ---
        sea_level_rise_rate=sea_level_rise_rate,  # Rate of SLR (m/yr)
        sea_level_rise_constant=sea_level_constant,  # Toggle for constant vs. changing SLR

        # --- Background shoreline change ---
        background_erosion=background_erosion,  # Array of background shoreline change (m/yr)

        # --- Grid and time settings ---
        alongshore_section_count=alongshore_section_count,  # Total number of domains (real + buffer)
        time_step_count=nt,  # Number of simulation time steps (years)
        min_dune_growth_rate=rmin,  # Min dune crest height per domain
        max_dune_growth_rate=rmax,  # Max dune crest height per domain
        num_cores=num_cores,  # Number of cores to use (for parallelism)

        # --- Module toggles (physical and management) ---
        roadway_management_module=roadway_management_on,  # Enables road elevation/setback rules
        alongshore_transport_module=True,  # Enables sediment transport between domains (BRIE)
        beach_nourishment_module=beach_dune_manager_on,  # Enables beach/dune nourishment actions
        community_economics_module=False,  # Not used in this hindcast (no economics)

        # --- Road configuration and thresholds ---
        road_ele=road_ele,  # Array of mean road elevations
        road_width=road_width,  # Width of the road (m)
        road_setback=road_setback,  # Array of road setback distances

        # --- Dune design thresholds ---
        dune_design_elevation=dune_design_elevation,  # Target dune height (m)
        dune_minimum_elevation=dune_minimum_elevation,  # Minimum allowed dune height (m)

        # --- Optional management settings ---
        trigger_dune_knockdown=trigger_dune_knockdown,  # Enable dune removal in certain conditions
        group_roadway_abandonment=group_roadway_abandonment,  # Grouped road removal logic (if used)
        nourishment_interval=None,  # Not used here (set to None)
        nourishment_volume=nourishment_volume,  # Array of nourishment volumes (m³/m)
        overwash_filter=overwash_filter,  # Overwash filter setting (percent)
        overwash_to_dune=overwash_to_dune,  # Whether to route overwash back to dunes

        # --- Sandbag settings (optional) ---
        sandbag_management_on=sandbag_management_on,  # Toggle for sandbag module
        sandbag_elevation=sandbag_elevation,  # Elevation threshold for sandbag deployment

        # --- Shoreline offset settings (optional) ---
        enable_shoreline_offset=enable_shoreline_offset,  # Toggle for shoreline offset tool
        shoreline_offset=shoreline_offset  # Array of shoreline offsets per domain
    )

    # ==============================================================================
    # Time Loop: Run simulation year by year and apply management actions
    # ==============================================================================

    # ------------------------------
    # Set absolute dune rebuild threshold
    # ------------------------------
    # Threshold is relative to the berm elevation of the first domain (assumes uniform berm)
    # If rebuild_dune_threshold = 0.3, then dune rebuilding will trigger when dune elevation < BermEl + 0.3 (m MHW)
    dune_rebuild_threshold = rebuild_dune_threshold + (cascade.barrier3d[0].BermEl * 10)

    # ------------------------------
    # Loop over time steps (years)
    # ------------------------------
    for time_step in range(nt - 1):

        # Display progress
        print("\r", "Time Step: ", time_step, end="")

        # Advance model one year
        cascade.update()

        # Stop simulation early if barrier failure flag is raised
        if cascade.b3d_break:
            break

        # Get the current time index (starts at 1)
        t = cascade.barrier3d[0].time_index

        # Initialize management decision arrays
        tmp_rebuild_dune = np.zeros(alongshore_section_count)  # Placeholder for dune rebuilding decisions
        tmp_nourish_now = np.zeros(alongshore_section_count)  # Placeholder for nourishment decisions

        # ------------------------------
        # Check conditions for each domain
        # ------------------------------
        for iB3D in range(alongshore_section_count):

            # Skip if community is no longer viable (e.g., due to repeated overwash)
            if cascade.community_break[iB3D]:
                continue

            # Check beach width and dune crest height for management triggers
            elif beach_dune_manager_on[iB3D]:

                # --- Beach width check: Nourishment trigger ---
                if cascade.nourishments[iB3D].beach_width[t - 1] < beach_width_threshold[iB3D]:
                    tmp_nourish_now[iB3D] = 1

                # --- Dune elevation check: Rebuild trigger ---
                # Get dune crest height (max elevation in each cross-shore row, then min across domain)
                dune_crests = cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
                min_dune_crest = (np.min(dune_crests) + cascade.barrier3d[iB3D].BermEl) * 10  # Convert to m MHW

                if min_dune_crest < dune_rebuild_threshold:
                    tmp_rebuild_dune[iB3D] = 1

        # ------------------------------
        # Apply management decisions
        # ------------------------------
        # Only act if all managed domains fall below thresholds (more realistic assumption)
        if np.all(tmp_nourish_now[beach_dune_manager_on]) == 1:
            cascade.nourish_now = tmp_nourish_now

        if np.all(tmp_rebuild_dune[beach_dune_manager_on]) == 1:
            cascade.rebuild_dune_now = tmp_rebuild_dune

    # ==============================================================================
    # Save Simulation Output
    # ==============================================================================

    # ------------------------------
    # Create structured output directory for this run
    # ------------------------------
    base_output_dir = "output"
    run_output_dir = os.path.join(base_output_dir, run_name)

    # Create directories if they don't exist
    os.makedirs(run_output_dir, exist_ok=True)

    # Save all output files into the run-specific folder
    cascade.save(run_output_dir)

    # Return to previous directory level (optional)
    os.chdir("..")

    # Return cascade object for post-processing or inspection
    return cascade

# ==============================================================================
# Barrier3D Configuration and Model Execution – Hatteras (92 segments + buffers)
# ==============================================================================

def alongshore_uniform():
    """
    Configure and execute a baseline simulation across 92 Barrier3D domains with 15
    buffer segments on each side (total = 122 CASCADE domains).
    """

    # ------------------------------
    # Core Configuration Parameters
    # ------------------------------

    # Number of real Barrier3D models (not including buffer segments)
    number_barrier3d_models = 92

    # ⚠️ Make sure the following lists have length = 92 (NOT 69):
    beach_width_threshold = [30] * number_barrier3d_models  # Beach nourishment trigger
    rmin = [0.55] * number_barrier3d_models                 # Minimum dune growth rate
    rmax = [0.95] * number_barrier3d_models                 # Maximum dune growth rate
    dune_design_elevation = [1.5] * number_barrier3d_models # Design dune height
    background_erosion = background_threhold_list           # Domain-specific erosion rates

    # Input file lists from previous section (length = 122: 15 + 92 + 15)
    elevation_file = e_file
    dune_file = d_file

    storm_file = s_file
    shoreline_offset = dune_offset  # Offset for each domain in the current year
    shoreline_offset_enabled = True

    # ------------------------------
    # Road & Management Parameters
    # ------------------------------
    num_cores = 4  # Adjust based on your hardware
    road_ele = 1.45     # m MHW – typical elevation from 1997 LiDAR
    road_width = 20     # Width in meters
    road_setback = road_setbacks      # Full-length list (122 entries)
    sandbag_management_on = road_cells  # Boolean mask for sandbag zones (122 entries)
    sandbag_elevation = 4  # m MHW

    # Thresholds for triggering management actions
    overwash_filter = 0           # 0 = residential (no filtering)
    overwash_to_dune = 9         # % of overwash used to rebuild dunes
    nourishment_volume = 100     # m³/m added during nourishment
    rebuild_dune_threshold = 1   # m MHW for rebuilding dune

    # ------------------------------
    # Management Activation Flags
    # ------------------------------
    roads_on = road_cells                                # Where roads are managed
    nourishments_on = [False] * number_barrier3d_models  # Toggle nourishment per domain

    # Sea-level rise configuration
    sea_level_rise_rate = 0.0056  # Linear rate (m/yr)
    sea_level_constant = True     # Set to False for exponential rise

    # ------------------------------
    # Run CASCADE Model
    # ------------------------------
    alongshore_connected(
        nt=run_years,
        name=run_name,
        storm_file=storm_file,
        alongshore_section_count=number_barrier3d_models,  # <-- 92 domains
        num_cores=num_cores,
        beach_width_threshold=beach_width_threshold,
        rmin=rmin,
        rmax=rmax,
        elevation_file=elevation_file,
        dune_file=dune_file,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=1.0,
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

# ==============================================================================
# Execute Simulation
# ==============================================================================

alongshore_uniform()
