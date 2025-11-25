# ~******* CASCADE RUN SCRIPT ********~
# This script is configured to run a single simulation with 44 connected domains.

import os
import numpy as np
from cascade.cascade import Cascade


# ###############################################################################
# # Main function for running a connected multi-domain simulation
# ###############################################################################

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
        sea_level_rise_rate=0.004,  # not an array
        sea_level_constant=True,  # not an array
        trigger_dune_knockdown=False,
        group_roadway_abandonment=None,
):
    """
    Initializes and runs a CASCADE simulation with multiple, interconnected
    barrier domains coupled with an alongshore sediment transport model.
    """
    # --------- INITIALIZE ---------
    datadir = "/Users/rsahrae/PycharmProjects/CASCADE_Coupling/CASCADE/data/PeaIsland"  # Assumes data is in a 'data' subdirectory
    cascade = Cascade(
        datadir,
        name,
        storm_file='/Users/rsahrae/PycharmProjects/CASCADE_Coupling/CASCADE/data/cascade-default-storms.npy',
        elevation_file=elevation_file,
        dune_file=dune_file,
        parameter_file="/Users/rsahrae/PycharmProjects/CASCADE_Coupling/CASCADE/data/notes-barrier3d-default-parameters.yaml",
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
        community_economics_module=False,
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
        trigger_dune_knockdown=trigger_dune_knockdown,
        group_roadway_abandonment=group_roadway_abandonment,
        nourishment_interval=None,
        nourishment_volume=nourishment_volume,
        overwash_filter=overwash_filter,
        overwash_to_dune=overwash_to_dune,
        inlet_module= True,
    )

    # --------- LOOP ---------
    for time_step in range(nt - 1):
        print(f"\rTime Step: {time_step}", end="")
        cascade.update()
        # current_width = cascade.barrier3d[2].InteriorWidth_AvgTS[-1] * 10
        # print(f"Time: {time_step + 1} yrs, Domain 3 Width: {current_width:.2f} m")
        if cascade.b3d_break:
            continue

    # --------- SAVE ---------
    save_directory = "/Users/rsahrae/PycharmProjects/CASCADE_Coupling/CASCADE/data/PeaIsland/Run_Output"
    cascade.save(save_directory)

    return cascade


# ###############################################################################
# # MAIN EXECUTION BLOCK: CONFIGURE AND RUN YOUR SIMULATION
# ###############################################################################

if __name__ == "__main__":
    # ------------------- ⚙️ USER CONFIGURATION -------------------

    # Set the total number of domains you are connecting
    number_of_domains = 15  # For your domains 78 through 100

    # Set the parameters for the simulation run
    simulation_time = 100  # Number of years to run
    output_name = "Connected_PeaIsland_15_Domain_Natural_Run"
    number_of_cores = 4  # Adjust based on your computer's capability

    # Set the physical parameters to be applied UNIFORMLY to all 23 domains
    # For a natural run (no human intervention), set management flags to False.

    # Dune growth rates [m/yr]
    uniform_rmin = 0.25
    uniform_rmax = 0.65

    # Background shoreline change [m/yr, negative for erosion]
    uniform_background_erosion = -1.0

    # Sea Level Rise [m/yr]
    slr_rate = 0.004

    # Management settings (set to False for a natural run)
    manage_roads = False
    manage_nourishments = False

    # ------------------- END OF CONFIGURATION --------------------

    # Programmatically generate the lists of your filenames
    elevation_file = [f"raster_domain_{i}_topography_2019.npy" for i in range(78, 93)]
    dune_file = [f"raster_domain_{i}_dune_2019.npy" for i in range(78, 93)]

    # Create the uniform parameter lists based on your settings above
    rmin_list = [uniform_rmin] * number_of_domains
    rmax_list = [uniform_rmax] * number_of_domains
    background_erosion_list = [uniform_background_erosion] * number_of_domains
    roads_on_list = [manage_roads] * number_of_domains
    nourishments_on_list = [manage_nourishments] * number_of_domains

    # Create placeholder lists for disabled modules (these values won't be used)
    beach_width_threshold_list = [30] * number_of_domains
    dune_design_elevation_list = [3.0] * number_of_domains
    dune_minimum_elevation_list = [2.0] * number_of_domains
    road_ele_list = [1] * number_of_domains
    overwash_filter_list = [0] * number_of_domains
    nourishment_volume_list = [0] * number_of_domains

    print(f"Starting CASCADE simulation: {output_name}")
    print(f"Connecting {number_of_domains} domains for {simulation_time} years...")

    # Call the function to run the simulation
    alongshore_connected(
        nt=simulation_time,
        name=output_name,
        storm_file='/Users/rsahrae/PycharmProjects/CASCADE_Coupling/CASCADE/data/cascade-default-storms.npy',
        alongshore_section_count=number_of_domains,
        num_cores=number_of_cores,
        elevation_file=elevation_file,
        dune_file=dune_file,
        rmin=rmin_list,
        rmax=rmax_list,
        background_erosion=background_erosion_list,
        roadway_management_on=roads_on_list,
        beach_dune_manager_on=nourishments_on_list,
        sea_level_rise_rate=slr_rate,
        # Pass in placeholder lists for unused parameters
        beach_width_threshold=beach_width_threshold_list,
        dune_design_elevation=dune_design_elevation_list,
        dune_minimum_elevation=dune_minimum_elevation_list,
        road_ele=road_ele_list,
        road_width=20,
        road_setback=20,
        overwash_filter=overwash_filter_list,
        overwash_to_dune=0,
        nourishment_volume=nourishment_volume_list,
        rebuild_dune_threshold=1,
    )

    print(f"\n\nSimulation '{output_name}' complete! Check the 'Run_Output' directory.")