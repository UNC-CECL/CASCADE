# Lexi Van Blunk
# 03/06/2026  updated for outwasher release

# code for running CASCADE for the 0%, 50%, and 100% outwash to shoreface scenario, 100 storm series
# change the dune growth rate
# specify datadir folders
# names of files will update automatically with the changes to the specified variables

import os
import time

from cascade.cascade import Cascade

# input datadir where the 100 storms are located
datadir = "C:/Users/Lexi/PycharmProjects/CASCADE/data/init_outwasher/"
storm_dir = "C:/Users/Lexi/PycharmProjects/CASCADE/data/init_outwasher/overwash_storms/"  # overwash storms

# ---------------------------------- set model parameters that change per run ------------------------------------------
r_dune_growth = 0.35  # 0.25 or 0.35
washout_array = [100, 50, 0]

# automatically set min and max r values based on dune growth rate selection
if r_dune_growth == 0.25:
    rname = "r025"
    min_dune_r = 0.05
    max_dune_r = 0.45
elif r_dune_growth == 0.35:
    rname = "r035"
    min_dune_r = 0.25
    max_dune_r = 0.45

# default beach slope
beach_slope = 0.006

for percent_washout_to_shoreface in washout_array:
    # save to pycharm folder
    save_dir = (
        "C:/Users/Lexi/PycharmProjects/CASCADE/data/results/{}/outwash{}/".format(
            rname, percent_washout_to_shoreface
        )
    )

    # --------------------------------- running overwash scenario with all overwash storms --------------------------------------
    list_storm_files = os.listdir(storm_dir)
    total_storms = len(list_storm_files)
    for storm_num in range(1, total_storms + 1):
        print("\r", "Storm Number: ", storm_num, end="")
        overwash_storm = storm_dir + (f"storm_{storm_num}.npy")

        ### 0%, 50%, or 100% washout to shoreface ------------------------------------------------------------------------------------
        # initialize class
        cascade_outwash = Cascade(
            datadir,
            name="outwash{}_storm_num{}".format(
                percent_washout_to_shoreface, storm_num
            ),
            elevation_file="outwash-default-elevation.npy",
            dune_file="outwash-default-dunes.npy",
            parameter_file="outwash-parameters.yaml",
            storm_file=overwash_storm,
            num_cores=1,  # cascade can run in parallel, can never specify more cores than that
            roadway_management_module=False,
            alongshore_transport_module=False,
            beach_nourishment_module=False,
            community_economics_module=False,
            outwash_module=True,
            alongshore_section_count=1,
            time_step_count=101,
            wave_height=1,  # ---------- for BRIE and Barrier3D --------------- #
            wave_period=7,
            wave_asymmetry=0.8,
            wave_angle_high_fraction=0.2,
            bay_depth=3.0,
            s_background=0.001,
            berm_elevation=1.46,
            MHW=0.36,
            beta=beach_slope,
            sea_level_rise_rate=0.004,
            sea_level_rise_constant=True,
            background_erosion=0.0,
            min_dune_growth_rate=min_dune_r,
            max_dune_growth_rate=max_dune_r,
            road_ele=1.7,  # ---------- roadway management --------------- #
            road_width=30,
            road_setback=30,
            dune_design_elevation=3.7,
            dune_minimum_elevation=2.2,
            trigger_dune_knockdown=False,
            group_roadway_abandonment=None,
            nourishment_interval=None,  # ---------- beach and dune ("community") management --------------- #
            nourishment_volume=300.0,
            overwash_filter=40,
            overwash_to_dune=10,
            number_of_communities=1,  # ---------- coastal real estate markets (in development) --------------- #
            sand_cost=10,
            taxratio_oceanfront=1,
            external_housing_market_value_oceanfront=6e5,
            external_housing_market_value_nonoceanfront=4e5,
            fixed_cost_beach_nourishment=2e6,
            fixed_cost_dune_nourishment=2e5,
            nourishment_cost_subsidy=10e6,
            house_footprint_x=15,
            house_footprint_y=20,
            beach_full_cross_shore=70,
            outwash_storms_file="outwash-default-storms.npy",  # --------- outwasher (in development) ------------ #
            percent_washout_to_shoreface=percent_washout_to_shoreface,
            outwash_beach_file="outwash-default-beach.npy",
        )

        # run the time loop/update function

        for time_step in range(cascade_outwash._nt - 1):
            print("\r", "Time Step: ", time_step + 1, end="")
            cascade_outwash.update()
            if cascade_outwash.b3d_break:
                break

        # save variables
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        cascade_outwash.save(save_dir)
