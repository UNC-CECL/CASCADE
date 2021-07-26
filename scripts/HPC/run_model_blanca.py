# Model runfile for CASCADE simulations on the HPC
# Written by K.Anarde

# remember if I move to a different computer to $ pip install -e . in the brie and B3D directories for the BMI

import numpy as np

from cascade.cascade import Cascade  # the new class


def RUN_6_B3D_Rave_SLR_pt004_Humans(
    nt,
    rmin,
    rmax,
    name,
    run_road_mgmt,
    storm_file,
    elevation_file,
    dune_file,
    road_ele=1.7,  # dummy values for the no management runs
    road_width=30,
    road_setback=30,
    dune_design_elevation=3.7,
    dune_minimum_elevation=2.2,
    percent_water_cells_sensitivity=None,
):

    # ###############################################################################
    # 6 - B3D with dune management (just the roadways module)
    # ###############################################################################
    # GOAL: Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 1000 years
    # or until the barrier drowns. All other modules (beach nourishment, community dynamics) turned off.

    # --------- INITIALIZE ---------
    datadir = "B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        alongshore_section_count=1,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=run_road_mgmt,
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        dune_design_elevation=dune_design_elevation,
        dune_minimum_elevation=dune_minimum_elevation,
    )

    if percent_water_cells_sensitivity is not None:
        cascade.roadways[
            0
        ]._percent_water_cells_touching_road = percent_water_cells_sensitivity

    # --------- LOOP ---------

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        # if cascade.road_break or cascade.b3d_break:
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "Run_Output/"
    cascade.save(save_directory)  # for now, this is a list

    return cascade
