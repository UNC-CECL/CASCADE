# Lexi Van Blunk
# 7/13/2023

# ### Katherine will need to change the datadir on line 14 and save_dir on line 43

# code for running CASCADE for the overwash only scenario, 100 storm series
# change the storm interval, dune growth rate, and configuration below
# names of files will update automatically with the changes to the specified variables

import time
from cascade.cascade import Cascade

# input datadir where the 100 storms are located
datadir = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/"

# ---------------------------------- set model parameters that change per run ------------------------------------------
storm_interval = 20   # 20 or 10 years
r_dune_growth = 0.35  # 0.25 or 0.35
config = 4            # 1, 2, 3, or 4

# automatically set min and max r values based on dune growth rate selection
if r_dune_growth == 0.25:
    rname = "r025"
    min_dune_r = 0.05
    max_dune_r = 0.45
elif r_dune_growth == 0.35:
    rname = "r035"
    min_dune_r = 0.25
    max_dune_r = 0.45

# automatically set beach slope based on configuration
if config == 1:
    beach_slope = 0.03
elif config == 2:
    beach_slope = 0.013
elif config == 3:
    beach_slope = 0.002
elif config == 4:
    beach_slope = 0.006

# save to pycharm folder
save_dir_b3d = "C:/Users/Lexi/PycharmProjects/CASCADE/cascade/data/outwash_data/storms/slope0pt03/run_output/{0}/overwash_only/".format(rname)

# -------------------- model parameters that are constant throughout the runs ------------------------------------------
ki = 8.75E-3
C = 0.0134

# --------------------------------- running overwash scenario with all 100 storms --------------------------------------
for storm_num in range(1, 101):
    overwash_storm = "StormSeries_100yrs_inclusive_NCB_Berm1pt46m_Slope0pt03_{0}.npy".format(storm_num)

    # ### barrier3D only, outwash module set to false ------------------------------------------------------------------
    # initialize class
    cascade_b3d_only = Cascade(
        datadir,
        name="config{0}_b3d_startyr1_interval{1}yrs_Slope0pt03_{2}".format(config, storm_interval, storm_num),
        elevation_file="NCB-default-elevation-config{0}-damMHW.npy".format(config),
        dune_file="NCB-default-dunes-config{0}-dam.npy".format(config),
        parameter_file="outwash-parameters.yaml",
        storm_file=overwash_storm,
        num_cores=1,  # cascade can run in parallel, can never specify more cores than that
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,
        outwash_module=False,
        alongshore_section_count=1,
        time_step_count=101,
        wave_height=1, # ---------- for BRIE and Barrier3D --------------- #
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
        nourishment_interval=None, # ---------- beach and dune ("community") management --------------- #
        nourishment_volume=300.0,
        overwash_filter=40,
        overwash_to_dune=10,
        number_of_communities=1, # ---------- coastal real estate markets (in development) --------------- #
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
        outwash_storms_file="outwash_storms_startyr_1_interval_{0}yrs.npy".format(storm_interval),  # --------- outwasher (in development) ------------ #
        percent_washout_to_shoreface=100,
        outwash_beach_file="NCB-default-beach-config{0}-damMHW.npy".format(config),
        dune_flow_dynamics="full",
        ki_value=ki,
        c=C,
    )

    # run the time loop/update function
    t0 = time.time()

    for time_step in range(cascade_b3d_only._nt - 1):
        print("\r", "Time Step: ", time_step + 1, end="")
        cascade_b3d_only.update()
        if cascade_b3d_only.b3d_break:
            break

    t1 = time.time()
    t_total_seconds = t1 - t0
    t_total_minutes = t_total_seconds / 60
    t_total_hours = t_total_seconds / 3600

    print("\n", round(t_total_minutes))
    print(round(t_total_hours, 1))

    # save variables
    cascade_b3d_only.save(save_dir_b3d)