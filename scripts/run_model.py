# run file for

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
----------------------------------------------------"""

# remember if I move to a different computer to $ pip install -e . in the brie and B3D directories for the BMI

import numpy as np
import os
import time

from scripts import CASCADE_plotters as CASCADEplt

from cascade.cascade import Cascade  # the new class

from barrier3d.tools.input_files import (
    yearly_storms,
    gen_dune_height_start,
    gen_alongshore_variable_rmin_rmax,
)

# for laptop and desktop, use all but one core; on supercomputer, use all cores
# num_cores = multiprocessing.cpu_count() - 1

# # ###############################################################################
# # runs
# # ###############################################################################


def RUN_1_CASCADE_LTA_COMPARISON(ny, nt, name):
    # ###############################################################################
    # 1 - CASCADE_LTA_COMPARISON
    # ###############################################################################
    # GOAL: highlight different processes in models with alongshore homogenous dune line, 3000 year simulation
    #
    # --------- INITIAL CONDITIONS ---------
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.002  # m/yr
    # ny = 6  #12 # number of alongshore sections (6=3 km for 3000 yr run, 12=6 km for 1500 yr run)
    # nt = 3000  # 3000  #1500 # timesteps for 3000 morphologic years
    rmin = 0.35  # minimum growth rate for logistic dune growth (can be a list)
    rmax = 0.85  # maximum growth rate for logistic dune growth (can be a list)

    # --------- INITIALIZE ---------
    # datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/" # iMAC
    datadir = (
        "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
    )
    brie, barrier3d = CASCADE.initialize(
        datadir,
        name,
        wave_height,
        wave_period,
        asym_frac,
        high_ang_frac,
        slr,
        ny,
        nt,
        rmin,
        rmax,
    )

    # --------- LOOP ---------
    brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)

    # --------- RUN LTA COMPARISON ---------
    # We need these parameters to be as similar as possible to "storm conditions" for B3D. NOTE: the LTA model does not work
    # well when you set the width and height directly from B3D so we set the critical barrier width to the width of the
    # (initial) B3D Interior Domain
    w_b_crit = 450  # critical barrier width [m]
    h_b_crit = 1.9  # (should equal B3D original BermEl in the yaml file, not what is presented in B3D (minus the MHW)
    Qow_max = 20  # max overwash flux [m3/m/yr]
    # Another NOTE: the LTA14 overwash model does best with a smaller cell size (dt) and time step (dt), so the values used
    # to define the grid size and time loop in the B3D run (i.e., ny and nt) are modified within CASCADE for BRIE (LTA14)
    brieLTA = CASCADE.LTA(
        name,
        wave_height,
        wave_period,
        asym_frac,
        high_ang_frac,
        slr,
        ny,
        nt,
        w_b_crit,
        h_b_crit,
        Qow_max,
    )

    # --------- SAVE ---------
    # #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    CASCADE.save(
        brie, barrier3d, save_directory, name
    )  # this returns the barrier3d model without the BMI

    return brieLTA


def RUN_2_AlongshoreVarGrowthParam_Alternating(name):
    # ###############################################################################
    # 2 - variable alongshore dune growth parameters
    # ###############################################################################
    # GOAL: what is the effect of the alongshore variability of dunes (15-30 km)?
    #   - vary the growth parameter by varying rmin and rmax, but keep difference (range) constant
    #        - [rmin = 0.35, raverage = 0.6, and rmax = 0.85 everywhere as control case] with diffusive wave parameters
    #        (look at Brie paper to see what conditions are considered diffusive, or high angle)
    #        - THIS RUN: 2 B3Ds at raverage = 0.45 (or 0.3) and 2 B3Ds at raverage=0.75 (or 0.9), all along the barrier, check that
    #        raverage is 0.6 across the barrier; np.mean([0.25, 0.65]) = 0.45 and np.mean([0.55, 0.95]) = 0.75
    #   - hypothesis is that it will prevent punctuated retreat

    # --------- INITIAL CONDITIONS ---------
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.3  # fraction of waves approaching from higher than 45 degrees
    slr = 0.002  # m/yr
    ny = 32  # number of alongshore sections (30=15 km, 60=30 km, 32=16 km)
    nt = 1000  # timesteps for 1000 morphologic years
    rmin = [
        0.25,
        0.25,
        0.55,
        0.55,
    ]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmin = rmin * int(ny / len(rmin))
    rmax = [
        0.65,
        0.65,
        0.95,
        0.95,
    ]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    rmax = rmax * int(ny / len(rmax))

    # --------- INITIALIZE ---------
    # #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    brie, barrier3d = CASCADE.initialize(
        datadir,
        name,
        wave_height,
        wave_period,
        asym_frac,
        high_ang_frac,
        slr,
        ny,
        nt,
        rmin,
        rmax,
    )

    # --------- LOOP ---------
    brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)

    # --------- SAVE ---------
    # #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    b3d = CASCADE.save(
        brie, barrier3d, save_directory, name
    )  # this returns the barrier3d model without the BMI


def RUN_3_AlongshoreVarGrowthParam_Gradient(slr, nt, name):
    # ###############################################################################
    # 3 - variable alongshore dune growth parameters (gradient)
    # ###############################################################################
    # GOAL: what is the effect of the alongshore variability of dunes?
    #        - THIS RUN: make gradient in raverage across the barrier and reduce the grid size to 6 km
    #        - Increased SLR to 0.004

    # --------- INITIAL CONDITIONS ---------
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    ny = 12  # number of alongshore sections (12 = 6 km)
    rmin = [
        0.25,
        0.25,
        0.25,
        0.35,
        0.35,
        0.35,
        0.45,
        0.45,
        0.45,
        0.55,
        0.55,
        0.55,
    ]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmax = [
        0.65,
        0.65,
        0.65,
        0.75,
        0.75,
        0.75,
        0.85,
        0.85,
        0.95,
        0.95,
        0.95,
        0.95,
    ]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    # rave = [0.45, 0.45, 0.45, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.75, 0.75, 0.75]  # to help me remember the average

    # --------- INITIALIZE ---------
    # datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/" # iMAC
    datadir = (
        "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
    )
    brie, barrier3d = CASCADE.initialize(
        datadir,
        name,
        wave_height,
        wave_period,
        asym_frac,
        high_ang_frac,
        slr,
        ny,
        nt,
        rmin,
        rmax,
    )

    # --------- LOOP ---------
    brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)

    # --------- SAVE ---------
    # save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE/Run_Output"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    CASCADE.save(
        brie, barrier3d, save_directory, name
    )  # this returns the barrier3d model without the BMI


def RUN_4_B3D_Rave_SLR_pt004(rmin, rmax, name):

    # ###############################################################################
    # 4 - check B3D dune growth parameters
    # ###############################################################################
    # GOAL: check which growth rate parameters in B3D undergo punctuated retreat (increased SLR to 0.004)

    # --------- INITIAL CONDITIONS ---------
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.004  # m/yr
    ny = 2  # number of alongshore sections (NOTE: this is just a dummy variable for this run)
    nt = 10000  # timesteps for 1000 morphologic years
    # rave = [0.45, 0.55., 0.65, 0.75]  # to help me remember the average

    # --------- INITIALIZE ---------
    # datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    brie, barrier3d = CASCADE.initialize(
        datadir,
        name,
        wave_height,
        wave_period,
        asym_frac,
        high_ang_frac,
        slr,
        ny,
        nt,
        rmin,
        rmax,
    )

    # --------- LOOP ---------
    # just use the first B3D grid and update B3D without brie coupling
    Time = time.time()

    for time_step in range(brie._nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        barrier3d[0].update()
        barrier3d[0].update_dune_domain()

    SimDuration = time.time() - Time
    print()
    print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

    # --------- SAVE ---------
    # save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE/Run_Output"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    CASCADE.save(
        brie, barrier3d, save_directory, name
    )  # this returns the barrier3d model without the BMI

    # os.chdir(save_directory)
    # b3d = barrier3d[0]
    # filename = name + ".npz"
    # np.savez(filename, barrier3d=b3d)

    # ===================================================
    # 7: Calculate shoreline change periodicity
    Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(
        b3d._x_s_TS
    )
    print("Barrier Punc = " + str(Punc) + " , Periodicity = " + str(Periodicity))

    # 2: Shoreline positions over time
    TMAX = b3d.time_index - 1
    CASCADEplt.plot_ShorelinePositions(b3d._x_s_TS[0:TMAX], b3d._x_b_TS[0:TMAX])

    return Periodicity, AvgFastDur, AvgSlowDur, Punc


def RUN_4_CASCADE_noAST_Rave_SLR_pt004(
    nt,
    rmin,
    rmax,
    name,
    storm_file,
    elevation_file,
    dune_file,
):

    # ###############################################################################
    # 4 - check B3D dune growth parameters
    # ###############################################################################
    # GOAL: check which growth rate parameters in B3D undergo punctuated retreat (increased SLR to 0.004)
    # later GOAL: look at the range of barrier behavioer for different dune growth rates -- 10k runs
    # AST turned off, and beach and dune management modules turned off

    # --------- INITIALIZE ---------
    datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        storm_file=storm_file,
        elevation_file=elevation_file,
        dune_file=dune_file,
        wave_height=1.0,
        wave_period=7,  # s (lowered from 10 s to reduce k_sf)
        wave_asymmetry=0.8,  # fraction approaching from left
        wave_angle_high_fraction=0.2,  # fraction of waves approaching from higher than 45 degrees
        sea_level_rise_rate=0.004,  # m/yr
        alongshore_section_count=1,  # only one B3D domain
        time_step_count=nt,
        min_dune_growth_rate=rmin,  # rave = [0.45, 0.55., 0.65, 0.75]  # to help me remember the average
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=False,  # no beach nourishment
        community_dynamics_module=False,  # no community dynamics
    )

    # --------- LOOP ---------
    Time = time.time()

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        # if cascade.road_break or cascade.b3d_break:
        if cascade.b3d_break:
            break

    SimDuration = time.time() - Time
    print()
    print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

    # --------- SAVE ---------
    save_directory = "/Users/katherineanarde/RESEARCH/RUN_OUTPUT"
    cascade.save(save_directory)  # for now, this is a list

    # # ===================================================
    # # 7: Calculate shoreline change periodicity
    # Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(
    #     b3d._x_s_TS
    # )
    # print("Barrier Punc = " + str(Punc) + " , Periodicity = " + str(Periodicity))
    #
    # # 2: Shoreline positions over time
    # TMAX = b3d.time_index - 1
    # CASCADEplt.plot_ShorelinePositions(b3d._x_s_TS[0:TMAX], b3d._x_b_TS[0:TMAX])
    #
    # return Periodicity, AvgFastDur, AvgSlowDur, Punc
    return cascade


def RUN_5_AlongshoreVarGrowthParam_half():

    # ###############################################################################
    # 5 - variable alongshore dune growth parameters (half/half)
    # ###############################################################################
    # GOAL: what is the effect of the alongshore variability of dunes?
    #        - THIS RUN: make half the barrier have different raverage
    #        - Increased SLR to 0.004

    # --------- INITIAL CONDITIONS ---------
    name = "5-VarGrowthParam_half_pt4SLR_1500yrs"
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.004  # m/yr
    ny = 12  # number of alongshore sections (12 = 6 km)
    nt = 1500  # timesteps for 1000 morphologic years
    rmin = [
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.55,
        0.55,
        0.55,
        0.55,
        0.55,
        0.55,
    ]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmax = [
        0.65,
        0.65,
        0.65,
        0.65,
        0.65,
        0.65,
        0.95,
        0.95,
        0.95,
        0.95,
        0.95,
        0.95,
    ]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    # rave = [0.45, 0.45, 0.45, 0.75, 0.75, 0.75]  # to help me remember the average

    # --------- INITIALIZE ---------
    # datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    brie, barrier3d = CASCADE.initialize(
        datadir,
        name,
        wave_height,
        wave_period,
        asym_frac,
        high_ang_frac,
        slr,
        ny,
        nt,
        rmin,
        rmax,
    )

    # --------- LOOP ---------
    # just use the first B3D grid and update B3D without brie coupling
    Time = time.time()

    for time_step in range(brie._nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        barrier3d[0].update()
        barrier3d[0].update_dune_domain()

    SimDuration = time.time() - Time
    print()
    print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

    # --------- SAVE ---------
    # save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE/Run_Output"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    CASCADE.save(
        brie, barrier3d, save_directory, name
    )  # this returns the barrier3d model without the BMI

    os.chdir(save_directory)
    b3d = barrier3d[0]
    filename = name + ".npz"
    np.savez(filename, barrier3d=b3d)

    # NOTE: the following code is from the iMac - not sure which one I used
    # def RUN_5_AlongshoreVarGrowthParam_half(name, ny, rmin, rmax):
    #     # name = '5-VarGrowthParam_half_pt4SLR_1500yrs'
    #     wave_height = 1.0  # m
    #     wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    #     asym_frac = 0.8  # fraction approaching from left
    #     high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    #     slr = 0.004  # m/yr
    #     nt = 1500  # timesteps for 1000 morphologic years
    #     # ny = 12  # number of alongshore sections (12 = 6 km)
    #     # rmin = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.55, 0.55, 0.55, 0.55, 0.55,
    #     #         0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    #     # rmax = [0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.95, 0.95, 0.95, 0.95, 0.95,
    #     #         0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    #     # # rave = [0.45, 0.45, 0.45, 0.75, 0.75, 0.75]  # to help me remember the average
    #
    #     # --------- INITIALIZE ---------
    #     datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"  # iMAC
    #     # datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
    #     brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt,
    #                                          rmin,
    #                                          rmax, datadir)
    #
    #     # --------- LOOP ---------
    #     brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)
    #
    #     # --------- SAVE ---------
    #     save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE"
    #     # save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    #     CASCADE.save(brie, barrier3d, save_directory, name)  # this returns the barrier3d model without the BMI


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
    # or until the barrier drowns. All other modules (beach nourishment, community dyanmics) turned off.

    # --------- INITIALIZE ---------
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
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
    Time = time.time()

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        # if cascade.road_break or cascade.b3d_break:
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    cascade.save(save_directory)  # for now, this is a list

    return cascade


def RUN_7_B3D_Rave_variableSLR_Humans(
    nt,
    rmin,
    rmax,
    name,
    road_ele,
    road_width,
    road_setback,
    dune_design_elevation,
    dune_minimum_elevation,
    run_road_mgmt,
    storm_file,
    elevation_file,
    dune_file,
    sea_level_rise_rate,
    sea_level_constant,
):

    # ###############################################################################
    # 7 - same as RUN 6 but with variable rates of SLR (i.e., no AST, no nourishment)
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
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
        sea_level_rise_rate=sea_level_rise_rate,
        sea_level_rise_constant=sea_level_constant,
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

    # --------- LOOP ---------
    Time = time.time()

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        # if cascade.road_break or cascade.b3d_break:
        if cascade.b3d_break:
            break

    # --------- SAVE ---------
    save_directory = "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    cascade.save(save_directory)  # for now, this is a list

    return cascade


def RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
    nt,
    rmin,
    rmax,
    name,
    dune_design_elevation,
    storm_file,
    elevation_file,
    dune_file,
    overwash_filter,
    nourishment_volume,
    beach_width_threshold,
    background_erosion,
):

    # ###############################################################################
    # 8 - nourish beach, rebuild dunes, and remove overwash from barrier interior
    # ###############################################################################

    # --------- INITIALIZE ---------
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
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
        sea_level_rise_constant=True,
        background_erosion=background_erosion,
        alongshore_section_count=1,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=1,
        roadway_management_module=False,  # no roadway management
        alongshore_transport_module=False,  # no brie coupling
        beach_nourishment_module=True,
        community_dynamics_module=False,  # no community dynamics
        dune_design_elevation=dune_design_elevation,
        nourishment_interval=None,  # yrs
        nourishment_volume=nourishment_volume,  # m^3/m
        overwash_filter=overwash_filter,  # % overwash removed
    )

    # --------- LOOP ---------
    Time = time.time()

    iB3D = 0  # we only have one Barrier3D domain here

    # after each year, check the beach width and dune elevation and decide if you want to nourish or rebuild the dune
    # next year with nourish_now parameter
    dune_rebuild_threshold = 0.3 + (
        cascade.barrier3d[iB3D].BermEl * 10
    )  # same threshold for absolute minimum elevation as in RoadwayManager (m MHW)

    for time_step in range(nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()
        # if cascade.community_break or cascade.b3d_break:
        if cascade.b3d_break:
            break

        # stop managing if the barrier becomes too narrow to sustain a community
        if cascade.community_break:
            pass
        else:
            t = cascade.barrier3d[iB3D].time_index

            if cascade.nourishments[iB3D].beach_width[t - 1] < beach_width_threshold:
                cascade.nourish_now[iB3D] = 1

            DuneDomainCrest = (
                cascade.barrier3d[iB3D].DuneDomain[t - 1, :, :].max(axis=1)
            )  # Maximum height of each row in dune domain [dam]
            # DuneRestart = cascade.barrier3d[iB3D].DuneRestart
            # DuneDomainCrest[DuneDomainCrest < DuneRestart] = DuneRestart
            DuneCrestMin = (
                np.min(DuneDomainCrest) + cascade.barrier3d[iB3D].BermEl
            ) * 10  # m MHW

            if DuneCrestMin < dune_rebuild_threshold:
                cascade.rebuild_dune_now[iB3D] = 1

    # --------- SAVE ---------
    save_directory = "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    cascade.save(save_directory)  # for now, this is a list

    return cascade


# # ###############################################################################
# # plotters
# # ###############################################################################


def PLOT_1_CASCADE_LTA_COMPARISON(brieLTA, name, save_directory):
    # --------- plot ---------
    # load the simulation if previously saved
    os.chdir("//")
    filename = name + ".npz"
    output = np.load(filename, allow_pickle=True)
    b3d = output["barrier3d"]
    # brie = output['brie']
    # brie = brie[0]

    # 1: Animation Frames of Barrier and Dune Elevation
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    ny = len(b3d)
    CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX, name)

    # ===================================================

    # 4: Cross-shore transects for both brieLTA and B3d
    iB3D = 0
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    time_step = [0, int(TMAX / 2), TMAX - 2]
    CASCADEplt.plot_ModelTransects(b3d, brieLTA, time_step, iB3D)

    # ===================================================

    # 5: Statistics from B3d
    TMAX = b3d[0].time_index - 1
    iB3D = 1
    CASCADEplt.plot_statistics(b3d, iB3D, TMAX)

    # 6: Statistics from BrieLTA
    TMAX = int((b3d[0].time_index - 1) / brieLTA._dt)
    iB3D = 1
    iB3D_BRIE = iB3D * int((b3d[0]._BarrierLength * 10) / brieLTA._dy)
    CASCADEplt.plot_statistics_BRIE(brieLTA, iB3D, TMAX)

    # 6: Statistics from both models for AGU presentation
    iB3D = 1
    TMAX = b3d[0].time_index - 1
    iBRIE = iB3D * int((b3d[0]._BarrierLength * 10) / brieLTA._dy)
    TMAX_BRIE = int((b3d[0].time_index - 1) / brieLTA._dt)
    CASCADEplt.plot_statisticsAGU(b3d, brieLTA, iB3D, iBRIE, TMAX, TMAX_BRIE)

    # ===================================================

    # 2: Shoreline positions over time (#6 in Barrier3D_Functions)
    TMAX = 1001
    CASCADEplt.plot_ShorelinePositions(b3d[0]._x_s_TS[0:TMAX], b3d[0]._x_b_TS[0:TMAX])

    # 3: Cross-Shore Transect for one subgrid every 100 m for last time step
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    CASCADEplt.plot_XShoreTransects(b3d[0], TMAX)


def PLOT_3_AlongshoreVarGrowthParam_gradient(name, save_directory):
    # --------- plot ---------
    filename = name + ".npz"
    os.chdir("/Run_Output")
    output = np.load(filename, allow_pickle=True)
    b3d = output["barrier3d"]
    CASCADE_b3d = (
        b3d  # CASCADE run: with AST (just rename to not get confused for #10 below)
    )

    # # 1: Animation Frames of Barrier and Dune Elevation
    # TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    # ny = len(b3d)
    # CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX, name)

    # ===================================================
    # 10: Differences in punctuated retreat for CASCADE vs B3D (AST model vs no AST)

    # individual B3D models (one 500-m domain) corresponding to the growth rates used in the run above
    output = np.load("4-B3D_Rave_pt45_SLR_pt004.npz", allow_pickle=True)
    b3d_pt45 = output["barrier3d"]
    output = np.load("4-B3D_Rave_pt55_SLR_pt004.npz", allow_pickle=True)
    b3d_pt55 = output["barrier3d"]
    output = np.load("4-B3D_Rave_pt65_SLR_pt004.npz", allow_pickle=True)
    b3d_pt65 = output["barrier3d"]
    output = np.load("4-B3D_Rave_pt75_SLR_pt004.npz", allow_pickle=True)
    b3d_pt75 = output["barrier3d"]

    b3d_only = []

    # the first 3 cells in the CASCADE run are rave = 0.45
    ny = 3
    for iB3D in range(ny):
        b3d_only.append(b3d_pt45[0])

    # the next 3 cells in the CASCADE run are rave = 0.55
    for iB3D in range(ny):
        b3d_only.append(b3d_pt55[0])

    # the next 3 cells in the CASCADE run are rave = 0.65
    for iB3D in range(ny):
        b3d_only.append(b3d_pt65[0])

    # the next 3 cells in the CASCADE run are rave = 0.55
    for iB3D in range(ny):
        b3d_only.append(b3d_pt75[0])

    CASCADEplt.plot_punctuated_difference(CASCADE_b3d, b3d_only, ny=12)

    # # 5: Statistics from B3d
    # TMAX = b3d[0].time_index - 1
    # iB3D = 0
    # CASCADEplt.plot_statistics(b3d, iB3D, TMAX)
    #
    # # 2: Shoreline positions over time
    # iB3D = 0
    # CASCADEplt.plot_ShorelinePositions(
    #     b3d[iB3D]._x_s_TS[0:TMAX], b3d[iB3D]._x_b_TS[0:TMAX]
    # )
    #
    # # 2: Shoreline positions over time
    # iB3D = 9
    # CASCADEplt.plot_ShorelinePositions(
    #     b3d[iB3D]._x_s_TS[0:TMAX], b3d[iB3D]._x_b_TS[0:TMAX]
    # )
    #
    # # 2: Shoreline change rate (AGU version)
    # # CASCADEplt.plot_ShorelineChangeRate_AGU(b3d1, b3d2)


def PLOT_5_AlongshoreVarGrowthParam_half(name, save_directory, ny):

    # --------- plot ---------
    filename = name + ".npz"
    os.chdir("/Run_Output")
    output = np.load(filename, allow_pickle=True)
    CASCADE_b3d = output["barrier3d"]  # CASCADE run: with AST

    # ===================================================
    # 10: Differences in punctuated retreat for CASCADE vs B3D (AST model vs no AST)

    # individual B3D models (one 500-m domain) corresponding to the growth rates used in the run above
    output = np.load("4-B3D_Rave_pt45_SLR_pt004.npz", allow_pickle=True)
    b3d_pt45 = output["barrier3d"]
    output = np.load("4-B3D_Rave_pt75_SLR_pt004.npz", allow_pickle=True)
    b3d_pt75 = output["barrier3d"]

    b3d_only = []

    # the first 6 cells in the CASCADE run are rave = 0.45
    nHalf = int(ny / 2)
    for iB3D in range(nHalf):
        b3d_only.append(b3d_pt45[0])

    # the next 6 cells in the CASCADE run are rave = 0.75
    for iB3D in range(nHalf):
        b3d_only.append(b3d_pt75[0])

    CASCADEplt.plot_punctuated_difference(CASCADE_b3d, b3d_only, ny=ny)

    # 1: Animation Frames of Barrier and Dune Elevation
    TMAX = CASCADE_b3d[0].time_index - 1  # just in case the barrier drowned
    ny = len(CASCADE_b3d)
    CASCADEplt.plot_ElevAnimation(CASCADE_b3d, ny, save_directory, TMAX, name)


def PLOT_5_Nonlinear_Dynamics_B3D(name, save_directory):

    # --------- plot ---------
    # filename = name + ".npz"
    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )
    # output = np.load(filename, allow_pickle=True)
    # CASCADE_b3d = output["barrier3d"]  # CASCADE run: with AST
    output = np.load(
        "4-B3D_Rave_pt45_SLR_pt004_10k-yrs.npz", allow_pickle=True
    )  # B3D low growth rate run, no AST
    b3d_pt45 = output["barrier3d"]
    output = np.load(
        "4-B3D_Rave_pt55_SLR_pt004_10k-yrs.npz", allow_pickle=True
    )  # B3D high growth rate run, no AST
    b3d_pt55 = output["barrier3d"]
    output = np.load(
        "4-B3D_Rave_pt65_SLR_pt004_10k-yrs.npz", allow_pickle=True
    )  # B3D high growth rate run, no AST
    b3d_pt65 = output["barrier3d"]
    output = np.load(
        "4-B3D_Rave_pt75_SLR_pt004_10k-yrs.npz", allow_pickle=True
    )  # B3D high growth rate run, no AST
    b3d_pt75 = output["barrier3d"]
    output = np.load(
        "4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs.npz", allow_pickle=True
    )  # B3D high growth rate run, no AST
    cascade_pt75_v2 = output["cascade"][0]
    b3d_pt75_v2 = cascade_pt75_v2.barrier3d
    # output = np.load(
    #     "4-B3D_Rave_pt75_SLR_pt004.npz", allow_pickle=True
    # )  # B3D high growth rate run, no AST
    # b3d_pt75 = output["barrier3d"]

    tmin = 0  # 500
    tmax = 10000  # 1000

    # individual dune growth rates
    ib3d = 0
    (
        BarrierWidth_45,
        DuneCrestMean_45,
        BarrierHeight_45,
    ) = CASCADEplt.plot_nonlinear_stats(b3d_pt45, ib3d, tmin, tmax)
    (
        BarrierWidth_55,
        DuneCrestMean_55,
        BarrierHeight_55,
    ) = CASCADEplt.plot_nonlinear_stats(b3d_pt55, ib3d, tmin, tmax)
    (
        BarrierWidth_65,
        DuneCrestMean_65,
        BarrierHeight_65,
    ) = CASCADEplt.plot_nonlinear_stats(b3d_pt65, ib3d, tmin, tmax)
    (
        BarrierWidth_75,
        DuneCrestMean_75,
        BarrierHeight_75,
        bw_rate_75,
        bh_rate_75,
        sc_rate_75,
        DuneCrestMin_75,
        DuneCrestMax_75,
    ) = CASCADEplt.plot_nonlinear_stats(b3d_pt75, ib3d, tmin, tmax)
    (
        BarrierWidth_75_v2,
        DuneCrestMean_75_v2,
        BarrierHeight_75_v2,
        bw_rate_75_v2,
        bh_rate_75_v2,
        sc_rate_75_v2,
        DuneCrestMin_75_v2,
        DuneCrestMax_75_v2,
    ) = CASCADEplt.plot_nonlinear_stats(b3d_pt75_v2, ib3d, tmin, 5724)

    CASCADEplt.nonlinear_comparison(
        DuneCrestMean_45,
        BarrierWidth_45,
        DuneCrestMean_55,
        BarrierWidth_55,
        DuneCrestMean_65,
        BarrierWidth_65,
        DuneCrestMean_75,
        BarrierWidth_75,
    )

    # save to text file for use in Matlab
    np.savetxt(
        "Outputs_pt45.txt",
        (
            np.array(BarrierWidth_45),
            np.array(BarrierHeight_45),
            np.array(DuneCrestMean_45),
        ),
    )

    # save final interior and dune domain from last time step to csv (note, dunes are low, maybe set road at 1.7 m)
    with open("b3d_pt45_10kyrs-elevations.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(b3d_pt45[0].InteriorDomain)  # save in decameters
    with open("b3d_pt45_10kyrs-dunes.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            b3d_pt45[0].DuneDomain[-1, :, 0]
        )  # save in decameters, just first row

    t = 6741  # low
    t = 6039  # high
    with open("b3d_pt45_6039yrs_high-elevations.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(b3d_pt45[0].DomainTS[6039])  # save in decameters
    with open("b3d_pt45_6039yrs_high-dunes.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            b3d_pt45[0].DuneDomain[6039, :, 0]
        )  # save in decameters, just first row

    t = 4929  # low
    t = 8793  # high
    with open("b3d_pt75_8793yrs_high-elevations.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(b3d_pt75[0].DomainTS[8793])  # save in decameters
    with open("b3d_pt75_8793yrs_high-dunes.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            b3d_pt75[0].DuneDomain[8793, :, 0]
        )  # save in decameters, just first row

    # now this is for the phase plots
    np.savetxt(
        "Outputs_pt75.txt",
        (
            np.array(BarrierWidth_75),
            np.array(BarrierHeight_75),
            np.array(DuneCrestMean_75),
        ),
    )
    # save final interior and dune domain from last time step to csv (note, dunes are low, maybe set road at 1.5 m)
    with open("b3d_pt75_10kyrs-elevations.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerows(b3d_pt75[0]._InteriorDomain)  # save in decameters
    with open("b3d_pt75_10kyrs-dunes.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(
            b3d_pt75[0].DuneDomain[-1, :, 0]
        )  # save in decameters, just first row

    # half/half run
    ib3d = 10
    CASCADEplt.plot_nonlinear_stats(CASCADE_b3d, ib3d, tmin, tmax)


def PLOT_5_Nonlinear_Dynamics_B3D_CNH(name_prefix, tmin, tmax, plot_name):

    os.chdir("/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output")

    # --------- plot ---------
    output = np.load(name_prefix + ".npz", allow_pickle=True)
    b3d = output["barrier3d"]

    # individual dune growth rates
    ib3d = 0  # this is just B3D, so no alongshore grid cells
    (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
    ) = CASCADEplt.plot_nonlinear_stats(b3d, ib3d, tmin, tmax)

    # save to text file for use in Matlab
    np.savetxt(
        name_prefix + ".txt",
        (
            np.array(BarrierWidth),
            np.array(BarrierHeight),
            np.array(DuneCrestMean),
        ),
    )

    # also make the gif
    directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    CASCADEplt.plot_ElevAnimation(b3d, 1, directory, TMAX=tmax, name=plot_name)


def PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
    name_prefix,
    tmin,
    tmax_roadways,
    tmax_sim,
    plot_name,
    run_road_mgmt,
):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    # --------- plot ---------
    output = np.load(name_prefix + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d

    # individual dune growth rates
    ib3d = 0  # this plotter is just for one B3D domain, no alongshore grid cells

    if run_road_mgmt:
        post_storm_dunes = cascade.roadways[ib3d]._post_storm_dunes
        post_storm_interior = cascade.roadways[ib3d]._post_storm_interior
        design_height = cascade.roadways[ib3d]._dune_design_elevation_TS
        rebuild_threshold = cascade.roadways[ib3d]._dune_minimum_elevation_TS
        road_elevation = cascade.roadways[ib3d]._road_ele_TS
        dunes_rebuilt = cascade.roadways[ib3d]._dunes_rebuilt_TS
        road_relocated = cascade.roadways[ib3d]._road_relocated_TS
    else:
        post_storm_dunes = None
        post_storm_interior = None
        design_height = None
        rebuild_threshold = None
        road_elevation = None
        dunes_rebuilt = None
        road_relocated = None

    (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bw_rate,
        bh_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
    ) = CASCADEplt.plot_nonlinear_stats_RoadwayManager(
        b3d,
        ib3d,
        tmax_roadways=tmax_roadways,
        tmax_sim=tmax_sim,
        post_storm_dunes=post_storm_dunes,
        post_storm_interior=post_storm_interior,
        design_height=design_height,
        rebuild_threshold=rebuild_threshold,
        road_elevation=road_elevation,
        dunes_rebuilt=dunes_rebuilt,
        road_relocated=road_relocated,
    )

    # # save to text file for use in Matlab
    # np.savetxt(
    #     name_prefix + ".txt",
    #     (
    #         np.array(BarrierWidth),
    #         np.array(BarrierHeight),
    #         np.array(DuneCrestMean),
    #     ),
    # )

    # # also make the gif
    # directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    # # CASCADEplt.plot_ElevAnimation(b3d, 1, directory, TMAX=tmax, name=plot_name)
    # if cascade.roadways is not None:  # added the roadways class
    #     CASCADEplt.plot_ElevAnimation_Humans_Roadways(
    #         cascade,
    #         1,
    #         directory,
    #         TMAX=tmax_sim,
    #         TMAX_roadways=tmax_roadways,
    #         name=plot_name,
    #     )
    # else:
    #     CASCADEplt.plot_ElevAnimation_Humans(
    #         cascade, 1, directory, TMAX=tmax, name=plot_name
    #     )

    return (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bw_rate,
        bh_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        cascade,
    )


def PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
    name_prefix, tmax_management, tmax_sim, plot_name
):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    # --------- plot ---------
    output = np.load(name_prefix + ".npz", allow_pickle=True)
    cascade = output["cascade"]
    cascade = cascade[0]
    b3d = cascade.barrier3d

    ib3d = 0  # this is just B3D, so no alongshore grid cells
    rebuild_threshold = 0.3 + (
        b3d[ib3d].BermEl * 10
    )  # min dune height above the berm [m MHW], same as in RoadwayManager

    (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,
    ) = CASCADEplt.plot_nonlinear_stats_BeachDuneManager(
        b3d,
        ib3d,
        tmax_management=tmax_management,
        tmax_sim=tmax_sim,
        nourishments=cascade.nourishments,
        post_storm_dunes=cascade.nourishments[ib3d]._post_storm_dunes,
        post_storm_x_s=cascade.nourishments[ib3d]._post_storm_x_s,
        post_storm_s_sf=cascade.nourishments[ib3d]._post_storm_s_sf,
        post_storm_ave_interior_width=cascade.nourishments[
            ib3d
        ]._post_storm_ave_interior_width,
        post_storm_ave_interior_height=cascade.nourishments[
            ib3d
        ]._post_storm_ave_interior_height,
        post_storm_beach_width=cascade.nourishments[ib3d]._post_storm_beach_width,
        post_storm_Qow=cascade.nourishments[ib3d]._post_storm_Qow,
        design_elevation=cascade.nourishments[ib3d]._dune_design_elevation,  # m MHW,
        rebuild_threshold=rebuild_threshold,
        dunes_rebuilt=cascade.nourishments[ib3d]._dunes_rebuilt_TS,
    )

    # # also make the gif
    # directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    # CASCADEplt.plot_ElevAnimation_Humans_BeachDuneManager(
    #     cascade, 1, directory, TMAX=tmax_management, name=plot_name, TMAX_SIM=tmax_sim
    # )

    return (
        BarrierWidth,
        DuneCrestMean,
        BarrierHeight,
        bh_rate,
        bw_rate,
        sc_rate,
        DuneCrestMin,
        DuneCrestMax,
        shoreline_position,
        shoreface_slope,
        beach_width,
        overwash,
        dune_toe,
        cascade,
    )


def PLOT_7_Initial_CNH_Topographies(name_prefix_list):

    os.chdir(
        "/Users/KatherineAnardeWheels/Research/BARis/UNC/CNH/CASCADE_save_dir/Run_Output"
    )

    cascade = []

    for i in range(0, len(name_prefix_list)):
        output = np.load(name_prefix_list[i] + ".npz", allow_pickle=True)
        csc8d = output["cascade"]
        cascade.append(csc8d[0])

    CASCADEplt.fig2_initialCNH_topo(cascade)

    return


# # ###############################################################################
# # runs
# # ###############################################################################

# record of non-human runs -------------------------------------------------------------------------------------
def early_runs():
    RUN_1_CASCADE_LTA_COMPARISON(
        ny=6, nt=3000, name="1-CASCADE_LTA_COMPARISON_3km_3000yr"
    )
    RUN_1_CASCADE_LTA_COMPARISON(
        ny=12, nt=1500, name="1-CASCADE_LTA_COMPARISON_6km_1500yr"
    )

    RUN_2_AlongshoreVarGrowthParam_Alternating(name="2-VarGrowthParam_Alternating")

    RUN_3_AlongshoreVarGrowthParam_Gradient(
        slr=0.002, nt=500, name="3-VarGrowthParam_grad_pt2HAF_pt2SLR_500yrs"
    )
    RUN_3_AlongshoreVarGrowthParam_Gradient(
        slr=0.002, nt=1500, name="3-VarGrowthParam_grad_pt2HAF_pt2SLR_1500yrs"
    )
    RUN_3_AlongshoreVarGrowthParam_Gradient(
        slr=0.004, nt=1500, name="3-VarGrowthParam_grad_pt2HAF_pt4SLR_1500yrs"
    )

    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.25, rmax=0.65, name="4-B3D_Rave_pt45_SLR_pt004"
    )  # rave = 0.45
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.25, rmax=0.65, name="4-B3D_Rave_pt45_SLR_pt004_v2"
    )  # rave = 0.45, used initial topography from pt45 year 6741 (10k run) and 3000 year storm
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.35, rmax=0.75, name="4-B3D_Rave_pt55_SLR_pt004"
    )  # rave = 0.55
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.45, rmax=0.85, name="4-B3D_Rave_pt65_SLR_pt004"
    )  # rave = 0.65
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.55, rmax=0.95, name="4-B3D_Rave_pt75_SLR_pt004"
    )  # rave = 0.75
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.55, rmax=0.95, name="4-B3D_Rave_pt75_SLR_pt004_v2"
    )  # rave = 0.75, used initial topography from pt45 year 6741 (10k run) and 3000 year storm
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.55, rmax=0.95, name="4-B3D_Rave_pt75_SLR_pt004_10k-yrs"
    )  # rave = 0.75
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.55, rmax=0.95, name="4-B3D_Rave_pt75_SLR_pt004_10k-yrs"
    )  # rave = 0.75
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.45, rmax=0.85, name="4-B3D_Rave_pt65_SLR_pt004_10k-yrs"
    )  # rave = 0.65
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.35, rmax=0.75, name="4-B3D_Rave_pt55_SLR_pt004_10k-yrs"
    )  # rave = 0.55
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.25, rmax=0.65, name="4-B3D_Rave_pt45_SLR_pt004_10k-yrs"
    )  # rave = 0.45

    # original RUN5, 6 km
    rmin = [
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.55,
        0.55,
        0.55,
        0.55,
        0.55,
        0.55,
    ]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmax = [
        0.65,
        0.65,
        0.65,
        0.65,
        0.65,
        0.65,
        0.95,
        0.95,
        0.95,
        0.95,
        0.95,
        0.95,
    ]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    # rave = [0.45, 0.45, 0.45, 0.75, 0.75, 0.75]  # to help me remember the average
    RUN_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs", ny=12, rmin=rmin, rmax=rmax
    )

    # 1 km
    rmin = [0.25, 0.55]
    rmax = [0.65, 0.95]
    # rave = [0.45, 0.75]  # to help me remember the average
    RUN_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_1km", ny=2, rmin=rmin, rmax=rmax
    )
    RUN_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_1km_v2", ny=2, rmin=rmin, rmax=rmax
    )  # used initial topography from pt45 year 6741 (10k run) and 3000 year storm

    # 3 km
    rmin = [0.25, 0.25, 0.25, 0.55, 0.55, 0.55]
    rmax = [0.65, 0.65, 0.65, 0.95, 0.95, 0.95]
    # rave = [0.45, 0.45, 0.45, 0.75, 0.75, 0.75]  # to help me remember the average
    RUN_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_3km", ny=6, rmin=rmin, rmax=rmax
    )
    RUN_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_3km_v2", ny=6, rmin=rmin, rmax=rmax
    )  # used initial topography from pt45 year 6741 (10k run) and 3000 year storm

    # 12 km
    ny = 24
    rmin = [0.25] * int(ny / 2) + [0.55] * int(ny / 2)
    rmax = [0.65] * int(ny / 2) + [0.95] * int(ny / 2)
    # rave = [0.45, 0.45, 0.45, 0.75, 0.75, 0.75]  # to help me remember the average
    RUN_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_12km", ny=ny, rmin=rmin, rmax=rmax
    )


def cascade_10kyr_sensitivity():

    cascade_10kyr_pt45_01 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.25,  # rave = 0.45 (but not 0.5 spaced like in Reeves et al., 2021 -- arbitrary)
        rmax=0.65,
        name="4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_01",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_02 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_02",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_02.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_03 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_03",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_03.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_04 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_04",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_04.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt45_05 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.25,  # rave = 0.45
        rmax=0.65,
        name="4-CASCADE_noAST_Rave_pt45_SLR_pt004_10k-yrs_05",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_05.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_01 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_01",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_02 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_02",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_02.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_03 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_03",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_03.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_04 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_04",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_04.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_05 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_05",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_05.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_Cbbr0pt5 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_Cbb0pt5",
        storm_file="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_old_storms = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_OLD_STORMS",
        storm_file="Default_StormTimeSeries_10k-yr.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    cascade_10kyr_pt75_old_storms_Cbbr0pt5 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_OLD_STORMS_Cbb0pt5",
        storm_file="Default_StormTimeSeries_10k-yr.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )

    # manually changed the berm elevation to 2.0 in the yaml
    cascade_10kyr_pt75_old_storms_BermEl2 = RUN_4_CASCADE_noAST_Rave_SLR_pt004(
        nt=10000,
        rmin=0.55,  # rave = 0.75
        rmax=0.95,
        name="4-CASCADE_noAST_Rave_pt75_SLR_pt004_10k-yrs_OLD_STORMS_BermEl2",
        storm_file="Default_StormTimeSeries_10k-yr.npy",
        elevation_file="InitElevHog.npy",
        dune_file="DuneStart_1000dam.npy",
    )


# record of non-human plots -------------------------------------------------------------------------------------
def early_plots():

    PLOT_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_6km",
        save_directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output",
        ny=12,
    )

    PLOT_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_1km",
        save_directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output",
        ny=2,
    )

    PLOT_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_3km",
        save_directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output",
        ny=6,
    )

    PLOT_5_AlongshoreVarGrowthParam_half(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs_12km",
        save_directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output",
        ny=24,
    )

    PLOT_3_AlongshoreVarGrowthParam_gradient(
        name="3-VarGrowthParam_grad_pt2HAF_pt4SLR_1500yrs",
        save_directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output",
    )

    # these plots actually don't use run 5, just the individual B3D runs
    PLOT_5_Nonlinear_Dynamics(
        name="5-VarGrowthParam_half_pt4SLR_1500yrs",
        save_directory="/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output",
    )


def cascade_10kyr_plots():

    # didn't finish
    CASCADEplt.supp_10kyr_timeseries()


# record of B3D time series -------------------------------------------------------------------------------------
def time_series():

    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"

    StormSeries_NormDist_10kyrs_01 = yearly_storms_from_MSSM(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_01",
    )

    StormSeries_NormDist_10kyrs_02 = yearly_storms_from_MSSM(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_02",
    )

    StormSeries_NormDist_10kyrs_03 = yearly_storms_from_MSSM(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_03",
    )

    StormSeries_NormDist_10kyrs_04 = yearly_storms_from_MSSM(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_04",
    )

    StormSeries_NormDist_10kyrs_05 = yearly_storms_from_MSSM(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=10000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_10kyrs_VCR_Berm1pt9m_Slope0pt04_05",
    )

    StormSeries_NormDist_1kyrs_01 = yearly_storms_from_MSSM(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01",
    )

    StormSeries_NormDist_1kyrs_02 = yearly_storms(
        datadir=datadir,
        storm_list_name="StormList_20k_VCR_Berm1pt9m_Slope0pt04.csv",  # can by .py or .csv
        mean_yearly_storms=8.3,
        SD_yearly_storms=5.9,
        MHW=0.46,  # m NAVD88
        StormStart=2,
        BermEl=1.9,  # m NAVD88, just used for plotting
        model_years=1000,
        bPlot=True,
        bSave=True,
        output_filename="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_02",
    )

    def BermEl_2m():
        name = "StormTimeSeries_10k-yr.npy"
        yearly_storms_from_MSSM(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=10000,  # note, this is the number of storms contained in the MSSM model. probably should make more.
        )

        name = "StormTimeSeries_3000yr.npy"
        yearly_storms_from_MSSM(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=3000,
        )

        name = "StormTimeSeries_1000yr.npy"
        yearly_storms_from_MSSM(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=1000,
        )

        name = "StormTimeSeries_200yr.npy"
        yearly_storms_from_MSSM(
            datadir=datadir,
            name=name,
            storm_list_name="VCRStormList.npy",
            mean_storm=8.3,
            SD_storm=5.9,
            MHW=0.46,
            StormStart=2,
            BermEl=1.9,
            model_years=200,
        )

    name = "DuneStart_1000dam.npy"
    gen_dune_height_start(datadir, name, Dstart=0.5, ny=1000)

    name = "growthparam_1000dam.npy"
    gen_alongshore_variable_rmin_rmax(datadir, name, rmin=0.35, rmax=0.85, ny=1000)


# record of human runs -------------------------------------------------------------------------------------
def human_runs():
    def SLR_sensitivity():

        # misnomer here: only running the natural scenario
        cascade_pt75_low_SLR0pt008 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="6-B3D_Rave_pt75_Natural_low_0pt008SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt75_low_SLR0pt012 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="6-B3D_Rave_pt75_Natural_low_0pt012SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt75_high_SLR0pt008 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="6-B3D_Rave_pt75_Natural_high_0pt008SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt75_high_SLR0pt012 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="6-B3D_Rave_pt75_Natural_high_0pt012SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_low_SLR0pt008 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="6-B3D_Rave_pt45_Natural_low_0pt008SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_low_SLR0pt012 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="6-B3D_Rave_pt45_Natural_low_0pt012SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_high_SLR0pt008 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="6-B3D_Rave_pt45_Natural_high_0pt008SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=True,
        )

        cascade_pt45_high_SLR0pt012 = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="6-B3D_Rave_pt45_Natural_high_0pt012SLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=True,
        )

        # start of accelerated SLR scenarios
        # for the accelerated SLR scenario, I had to hard code the parameters that correspond to the
        # Rohling et al. (2013) 68% upper bound for AD2000-2200. SLRR starts at 0.003 m/yr and ends at 0.022 m/yr;
        # matches with the bounds of RCP8.5 SLR by 2100 and 2200
        cascade_pt75_low_SLRacc = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="6-B3D_Rave_pt75_Natural_low_AccSLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # dummy
            sea_level_constant=False,  # accelerated
        )

        cascade_pt75_high_SLRacc = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="6-B3D_Rave_pt75_Natural_high_AccSLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_829yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=False,
        )

        cascade_pt45_low_SLRacc = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="6-B3D_Rave_pt45_Natural_low_AccSLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.012,  # m/yr
            sea_level_constant=False,
        )

        cascade_pt45_high_SLRacc = RUN_7_B3D_Rave_variableSLR_Humans(
            nt=200,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="7-B3D_Rave_pt45_Natural_high_AccSLR",
            road_ele=None,
            road_width=None,
            road_setback=None,
            dune_design_elevation=None,
            dune_minimum_elevation=None,
            run_road_mgmt=False,
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_802yrs_high-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            sea_level_rise_rate=0.008,  # m/yr
            sea_level_constant=False,
        )

    def roadways():
        def pt75():

            cascade_pt75_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Natural_low",
                run_road_mgmt=False,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 162 years, 20.0% of road borders water
            # lower initial roadway for "low" case vs "high" case
            cascade_pt75_h2m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low",
                road_ele=1.2,  # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 78 years, 20.0% of road borders water
            cascade_pt75_h3m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_low",
                road_ele=1.2,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=4.2,  # m MHW, rebuild to 3 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, so really COMPLETELY erode before rebuilding (this is essentially the dune toe)
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 301 years
            cascade_pt75_h1m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_low",
                road_ele=1.2,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=2.2,  # m MHW, rebuild to 1 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            cascade_pt75_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Natural_high",
                run_road_mgmt=False,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # able to set the highway slightly higher
            # Roadway width drowned at 531 years, 20.0% of road borders water
            # Barrier has HEIGHT DROWNED at t = 582 years
            cascade_pt75_h2m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_high",
                road_ele=2.0,  # 1.7,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=4.0,  # 3.7,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.5,  # 2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 417 years
            cascade_pt75_h3m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_high",
                road_ele=2.0,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=5.0,  # m MHW, rebuild to 3 m dune above the roadway
                dune_minimum_elevation=2.5,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Island is to narrow for roadway to be relocated. Roadway eaten up by dunes at 537 years
            cascade_pt75_h1m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_high",
                road_ele=2.0,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.0,  # m MHW, rebuild to 1 m dune above the roadway
                dune_minimum_elevation=2.5,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

        def pt45():

            cascade_pt45_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Natural_low",
                run_road_mgmt=False,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # note that this roadway is higher than the pt75 "low" scenario (1.2 m)
            # Roadway width drowned at 580 years, 20.0% of road borders water
            cascade_pt45_h2m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low",
                road_ele=1.7,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 430 years, 20.0% of road borders water
            # Barrier has HEIGHT DROWNED at t = 442 years
            cascade_pt45_h3m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_low",
                road_ele=1.7,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 595 years, 20.0% of road borders water
            cascade_pt45_h1m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_low",
                road_ele=1.7,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            cascade_pt45_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Natural_high",
                run_road_mgmt=False,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 631 years, 20.0% of road borders water
            cascade_pt45_h2m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_high",
                road_ele=2.0,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=4.0,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=2.5,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 522 years, 20.0% of road borders water
            # Barrier has HEIGHT DROWNED at t = 532 years
            cascade_pt45_h3m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_high",
                road_ele=2.0,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=5.0,  # m MHW, rebuild to 3 m dune above the roadway
                dune_minimum_elevation=2.5,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

            # Roadway width drowned at 731 years, 20.0% of road borders water
            cascade_pt45_h1m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=1000,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_high",
                road_ele=2.0,
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.0,  # m MHW, rebuild to 1 m dune above the roadway
                dune_minimum_elevation=2.5,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
            )

        def roadway_sensitivity_abandonment_criteria():

            # test the sensitivity of varying the number of water cells that border the roadway as a metric to stop
            # managing the road for the most extreme barrier trajectory (high dune growth rate, low barrier)

            # Roadway width drowned at 159 years, 10.0% of road borders water
            cascade_pt75_h2m_low_10percent = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=500,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_10percent",
                road_ele=1.2,  # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                percent_water_cells_sensitivity=0.1,
            )

            # Roadway width drowned at 162 years, 20.0% of road borders water
            cascade_pt75_h2m_low_20percent = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=500,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_20percent",
                road_ele=1.2,  # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                percent_water_cells_sensitivity=0.2,
            )

            # Roadway width drowned at 165 years, 30.0% of road borders water
            cascade_pt75_h2m_low_30percent = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=500,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_30percent",
                road_ele=1.2,  # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                percent_water_cells_sensitivity=0.3,
            )

            # Roadway width drowned at 173 years, 40.0% of road borders water
            cascade_pt75_h2m_low_40percent = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=500,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_40percent",
                road_ele=1.2,  # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                percent_water_cells_sensitivity=0.4,
            )

            # Roadway width drowned at 182 years, 50.0% of road borders water
            cascade_pt75_h2m_low_50percent = RUN_6_B3D_Rave_SLR_pt004_Humans(
                nt=500,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_50percent",
                road_ele=1.2,  # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
                road_width=20,  # m
                road_setback=20,  # m
                dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
                dune_minimum_elevation=1.7,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                run_road_mgmt=True,
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                percent_water_cells_sensitivity=0.5,
            )

        def old_overwash_model():
            def pt45():
                # used 3000 year storm time series for the _low and _high runs
                cascade_pt45_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=550,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_550yrs_Natural_low",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                )

                cascade_pt45_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=750,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_750yrs_Natural_high",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                )

                # note that this roadway is higher than the pt75 "low" scenario, and equal to the "high" scenario (just higher
                # starting topoagraphy); drowned at 388 years
                cascade_pt45_h2m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=387,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_387yrs_Roadways_2mDune_20mSetback_20mWidth_low",
                    road_ele=1.7,  # average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (I don't want to set it lower than that)
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                # note that this roadway is higher than the pt75 "high" scenario (just higher starting topoagraphy); drowned at 515 years
                cascade_pt45_h2m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=515,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_515yrs_Roadways_2mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                # drowned at 370 years
                cascade_pt45_h3m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=369,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_369yrs_Roadways_3mDune_20mSetback_20mWidth_low",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                cascade_pt45_h3m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=473,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_473yrs_Roadways_3mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                cascade_pt45_h1m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=503,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_503yrs_Roadways_1mDune_20mSetback_20mWidth_low",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                cascade_pt45_h1m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=700,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_700yrs_Roadways_1mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

            def pt75():
                # REMEMBER TO SWITCH TOPOGRAHPHY FILES ------------------------
                # used 1000 year storm time series for the _low and _high runs

                # ACCIDENTALLY reran and saved over this one when checking the new roadways class...
                # will want to check with new plots that they are the same
                cascade_pt75_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Natural_low",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                    storm_file="Default_StormTimeSeries_1000yr.npy",
                    elevation_file="b3d_pt75_4929yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                cascade_pt75_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=450,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    # name="6-B3D_Rave_pt75_450yrs_Natural_high",
                    name="test2",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                    storm_file="Default_StormTimeSeries_1000yr.npy",
                    elevation_file="b3d_pt75_8793yrs_high-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                # lowered roadway and decreased setback to accommodate low-lying barrier, roadway drowned at 98 from back bay
                cascade_pt75_h2m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=97,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    # name="6-B3D_Rave_pt75_97yrs_Roadways_2mDune_20mSetback_20mWidth_low",
                    name="test",
                    road_ele=1.4,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.4,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=1.9,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                    storm_file="Default_StormTimeSeries_1000yr.npy",
                    elevation_file="b3d_pt75_4929yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                # able to set the highway slightly higher (0.3 m), kept setback the same; drowned at 428 years
                cascade_pt75_h2m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=427,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_425yrs_Roadways_2mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

                # lowered roadway and decreased setback to accommodate low-lying barrier, drowned at 71 years
                cascade_pt75_h3m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=70,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    # name="6-B3D_Rave_pt75_70yrs_Roadways_3mDune_20mSetback_20mWidth_low",
                    name="test",
                    road_ele=1.4,  # average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (I don't want to set it lower than that)
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.4,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=1.9,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                    storm_file="Default_StormTimeSeries_3000yr.npy",
                    elevation_file="b3d_pt75_4929yrs_low-elevations.csv",
                    dune_file="barrier3d-default-dunes.npy",
                )

                # able to set the highway slightly higher (0.3 m), kept setback the same; drowned at 358 years
                cascade_pt75_h3m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=357,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_357yrs_Roadways_3mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

                # lowered roadway and decreased setback to accommodate low-lying barrier, roadway drowned at 117 years
                cascade_pt75_h1m_low = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=116,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_116yrs_Roadways_1mDune_20mSetback_20mWidth_low",
                    road_ele=1.4,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.4,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=1.9,  # m MHW, allow dune to erode down to 0.5 m above the roadway
                    run_road_mgmt=True,
                )

                # able to set the highway slightly higher (0.3 m), kept setback the same; drowned at 439 years
                cascade_pt75_h1m_high = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=438,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_438yrs_Roadways_1mDune_20mSetback_20mWidth_high",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=20,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

            # OLD versions using the final output from 10,000 year run -------------------------------
            def old_10k_initial_elevation():
                b3d_pt45 = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Natural_v2",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                )

                # v1 - didn't drown roadway, v2 - roadway drowned at 160, v3 - new class
                # b3d_pt45_h2m, dunes_rebuilt, road_overwash_volume = RUN_6_B3D_Rave_SLR_pt004_Humans(
                cascade_pt45_h2m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=159,  # 200
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_20mWidth_v3",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                    run_road_mgmt=True,
                )

                b3d_pt45_h3m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 3.7 m
                )

                # v1 drowned at 91 years, v2 - roadway drowned at 101 years
                # b3d_pt45_h1m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                cascade_pt45_h1m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=100,  # 90
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2_classtest2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.4 m
                    run_road_mgmt=True,
                )

                # increase road width
                b3d_pt45_h4 = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_30mWidth_v2",
                    road_ele=1.7,
                    road_width=30,  # m
                    road_setback=40,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                )

                # change setback distance
                # v1 drowned at 183 years
                b3d_pt45_h5 = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,  # 182
                    rmin=0.25,
                    rmax=0.65,  # rave = 0.45
                    name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_30mSetback_20mWidth",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=30,  # m
                    dune_design_elevation=3.7,  # m MHW, 2 m dune above the roadway
                    dune_minimum_elevation=2.7,  # m MHW, 1 m dune above the roadway
                )

                # REMEMBER TO SWITCH TOPOGRAHPHY FILES ------------------------
                b3d_pt75 = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Natural",
                    road_ele=None,
                    road_width=None,
                    road_setback=None,
                    dune_design_elevation=None,
                    dune_minimum_elevation=None,
                    run_road_mgmt=False,
                )

                # v2, roadway drowned at 157 years
                b3d_pt75_h2m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=156,  # 200
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Roadways_2mDune_40mSetback_20mWidth_v2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
                )

                b3d_pt75_h3m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=200,
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=4.7,  # m MHW, rebuild to 3 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 3.7 m
                )

                # v1 drowned at 88 years, v2 drowned at 105 years, v3 roadway couldn't be relocated at 408 years
                b3d_pt75_h1m = RUN_6_B3D_Rave_SLR_pt004_Humans(
                    nt=407,  # 87, 104
                    rmin=0.55,
                    rmax=0.95,  # rave = 0.75
                    name="6-B3D_Rave_pt75_407yrs_Roadways_1mDune_40mSetback_20mWidth",
                    road_ele=1.7,
                    road_width=20,  # m
                    road_setback=40,  # m
                    dune_design_elevation=2.7,  # m MHW, rebuild to 1 m dune above the roadway
                    dune_minimum_elevation=2.2,  # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.4 m
                    run_road_mgmt=True,
                )

    def nourishments():

        # note, here we keep all other variables the same for comparison to the roadways scenarios; note that the dune
        # is only rebuilt in the BeachDuneManager when it is totally wiped out (we specify 0.3 m above the berm)

        # Community reached minimum width, drowned at 178 years; roadway scenario drowned at 162 years
        cascade_pt75_h2m_low_nourishment_residential = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_residential",
            dune_design_elevation=3.2,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=40,  # corresponds with residential
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=0.0,
        )

        # Community reached minimum width, drowned at 90 years
        cascade_pt75_h2m_low_nourishment_commercial = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial",
            dune_design_elevation=3.2,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=90,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=0.0,
        )

        # Community reached minimum width, drowned at 90 years
        cascade_pt75_h2m_low_nourishment_commercial_background_erosion = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
            dune_design_elevation=3.2,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=90,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=-0.25,  # m/yr, background shoreline erosion
        )

        # Community reached minimum width, drowned at 90 years
        cascade_pt75_h2m_low_nourishment_commercial_background_erosion_1m = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.55,
            rmax=0.95,  # rave = 0.75
            name="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
            dune_design_elevation=3.2,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=90,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=-1.0,  # m/yr, background shoreline erosion
        )

        # roadway scenario drowned at 404 years
        # Community reached minimum width, drowned at 496 years
        cascade_pt45_h2m_low_nourishment_residential = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential",
            dune_design_elevation=3.7,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=40,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=0.0,
        )

        # Community reached minimum width, drowned at 421 years
        # Barrier has HEIGHT DROWNED at t = 458 years
        cascade_pt45_h2m_low_nourishment_commercial = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial",
            dune_design_elevation=3.7,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=90,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=0.0,
        )

        # Community reached minimum width, drowned at 421 years
        # Barrier has HEIGHT DROWNED at t = 454 years
        cascade_pt45_h2m_low_nourishment_commercial_background_erosion_pt25m = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
            dune_design_elevation=3.7,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=90,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=-0.25,  # m/yr, background shoreline erosion
        )

        # Community reached minimum width, drowned at 421 years
        cascade_pt45_h2m_low_nourishment_commercial_background_erosion_1m = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
            nt=500,
            rmin=0.25,
            rmax=0.65,  # rave = 0.45
            name="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
            dune_design_elevation=3.7,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
            storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
            elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
            dune_file="barrier3d-default-dunes.npy",
            overwash_filter=90,  # corresponds with commercial
            nourishment_volume=100,  # m^3/m
            beach_width_threshold=30,  # m
            background_erosion=-1.0,  # m/yr, background shoreline erosion
        )

        def topo_only():
            # we only run 10 years of the following runs because we use them for plotting the initial topo figure for
            # the CNH simulations
            cascade_pt45_h2m_high_nourishment_commercial = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
                nt=10,
                rmin=0.25,
                rmax=0.65,  # rave = 0.45
                name="8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial",
                dune_design_elevation=3.7,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt45_802yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                overwash_filter=90,  # corresponds with commercial
                nourishment_volume=100,  # m^3/m
                beach_width_threshold=30,  # m
                background_erosion=0.0,
            )

            cascade_pt75_h2m_high_nourishment_commercial = RUN_8_CASCADE_Rave_SLR_pt004_Nourishment(
                nt=10,
                rmin=0.55,
                rmax=0.95,  # rave = 0.75
                name="8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial",
                dune_design_elevation=3.2,  # m MHW, keep dune design height the same as 2 m dune above the "roadway"
                storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
                elevation_file="b3d_pt75_829yrs_high-elevations.csv",
                dune_file="barrier3d-default-dunes.npy",
                overwash_filter=90,  # corresponds with commercial
                nourishment_volume=100,  # m^3/m
                beach_width_threshold=30,  # m
                background_erosion=0.0,
            )


# record of human plots -------------------------------------------------------------------------------------
def human_plots():
    def old_overwash_model():

        # rave = 0.45 runs, low
        def old_10k_initial_elevation():
            # save text files for phase plots and basic statistics
            PLOT_5_Nonlinear_Dynamics_B3D_CNH(
                name_prefix="6-B3D_Rave_pt45_200yrs_Natural_v2",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_plots",
            )
            PLOT_5_Nonlinear_Dynamics_B3D_CNH(
                name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_20mWidth_v2",
                tmin=0,
                tmax=159,  # 200
                plot_name="b3d_pt45_h2m_plots_v2",
            )
            PLOT_5_Nonlinear_Dynamics_B3D_CNH(
                name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_3mDune_40mSetback_20mWidth",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_h3m_plots",
            )
            PLOT_5_Nonlinear_Dynamics_B3D_CNH(
                name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2",
                tmin=0,
                tmax=100,  # 90
                plot_name="b3d_pt45_h1m_plots_v2",
            )
            PLOT_5_Nonlinear_Dynamics_B3D_CNH(
                name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_30mWidth",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_h4_plots",
            )
            PLOT_5_Nonlinear_Dynamics_B3D_CNH(
                name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_30mSetback_20mWidth",
                tmin=0,
                tmax=182,
                plot_name="b3d_pt45_h5_plots",
            )

        # rave = 0.45 runs, low
        def pt45_low():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_550yrs_Natural_low",
                tmin=0,
                tmax=550,
                plot_name="b3d_pt45_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_387yrs_Roadways_2mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax=387,
                plot_name="b3d_pt45_h2m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_369yrs_Roadways_3mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax=369,
                plot_name="b3d_pt45_h3m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_503yrs_Roadways_1mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax=503,
                plot_name="b3d_pt45_h1m_plots_low",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[515, 503, 387, 369],
            )

        # rave = 0.45 runs, high
        def pt45_high():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_750yrs_Natural_high",
                tmin=0,
                tmax=750,
                plot_name="b3d_pt45_plots_high",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_515yrs_Roadways_2mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax=515,
                plot_name="b3d_pt45_h2m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_473yrs_Roadways_3mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax=473,
                plot_name="b3d_pt45_h3m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_700yrs_Roadways_1mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax=700,
                plot_name="b3d_pt45_h1m_plots_high",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[715, 700, 515, 473],
            )

        # rave = 0.75 runs, low
        def pt75_low():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_200yrs_Natural_low",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_97yrs_Roadways_2mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax=97,
                plot_name="b3d_pt75_h2m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_70yrs_Roadways_3mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax=70,
                plot_name="b3d_pt75_h3m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_116yrs_Roadways_1mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax=116,  # 87, 104
                plot_name="b3d_pt75_h1m_plots_low",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=150,
            )

        # rave = 0.75 runs, high
        def pt75_high():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_450yrs_Natural_high",
                tmin=0,
                tmax=450,
                plot_name="b3d_pt75_plots_high",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_425yrs_Roadways_2mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax=425,
                plot_name="b3d_pt75_h2m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_357yrs_Roadways_3mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax=357,
                plot_name="b3d_pt75_h3m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_438yrs_Roadways_1mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax=438,  # 87, 104
                plot_name="b3d_pt75_h1m_plots_high",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=450,
            )

    def roadways():

        # rave = 0.75 runs, low
        def pt75_low():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low",
                tmin=0,
                tmax_roadways=1000,  # dummy
                tmax_sim=1000,
                plot_name="b3d_pt75_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax_roadways=161,
                tmax_sim=1000,
                plot_name="b3d_pt75_h2m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax_roadways=77,  # drowned at 78
                tmax_sim=1000,
                plot_name="b3d_pt75_h3m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax_roadways=300,
                tmax_sim=1000,
                plot_name="b3d_pt75_h1m_plots_low",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[750, 750, 750, 750],
                tmax_management=[
                    0,
                    300,
                    161,
                    77,
                ],  # h3, h1, h2 - 300, 161, 77 roadways drowned
            )

        # rave = 0.75 runs, high
        def pt75_high():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_high",
                tmin=0,
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt75_plots_high",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax_roadways=530,
                tmax_sim=582,  # this may break, want to try 581?
                plot_name="b3d_pt75_h2m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_3mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax_roadways=416,
                tmax_sim=1000,
                plot_name="b3d_pt75_h3m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_1mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax_roadways=536,
                tmax_sim=1000,
                plot_name="b3d_pt75_h1m_plots_high",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[
                    750,
                    750,
                    582,
                    750,
                ],  # # h3, h1, h2 - 536, 530, 416 roadways drowned
                tmax_management=[0, 536, 530, 416],
            )

        # rave = 0.45 runs, low
        def pt45_low():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_low",
                tmin=0,
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt45_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax_roadways=579,
                tmax_sim=1000,
                plot_name="b3d_pt45_h2m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax_roadways=429,
                tmax_sim=442,
                plot_name="b3d_pt45_h3m_plots_low",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_low",
                tmin=0,
                tmax_roadways=594,
                tmax_sim=1000,
                plot_name="b3d_pt45_h1m_plots_low",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[750, 750, 750, 442],
                tmax_management=[0, 594, 579, 429],
            )

        # rave = 0.45 runs, high
        def pt45_high():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_high",
                tmin=0,
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt45_plots_high",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_h2m,
                DuneCrestMean_h2m,
                BarrierHeight_h2m,
                bw_rate_h2m,
                bh_rate_h2m,
                sc_rate_h2m,
                DuneCrestMin_h2m,
                DuneCrestMax_h2m,
                cascade_h2m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_2mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax_roadways=630,
                tmax_sim=1000,
                plot_name="b3d_pt45_h2m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h3m,
                DuneCrestMean_h3m,
                BarrierHeight_h3m,
                bw_rate_h3m,
                bh_rate_h3m,
                sc_rate_h3m,
                DuneCrestMin_h3m,
                DuneCrestMax_h3m,
                cascade_h3m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_3mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax_roadways=521,
                tmax_sim=532,
                plot_name="b3d_pt45_h3m_plots_high",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_h1m,
                DuneCrestMean_h1m,
                BarrierHeight_h1m,
                bw_rate_h1m,
                bh_rate_h1m,
                sc_rate_h1m,
                DuneCrestMin_h1m,
                DuneCrestMax_h1m,
                cascade_h1m,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Roadways_1mDune_20mSetback_20mWidth_high",
                tmin=0,
                tmax_roadways=730,
                tmax_sim=1000,
                plot_name="b3d_pt45_h1m_plots_high",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                cascade=[cascade_nat, cascade_h1m, cascade_h2m, cascade_h3m],
                DuneCrestMin=[
                    DuneCrestMin_nat,
                    DuneCrestMin_h1m,
                    DuneCrestMin_h2m,
                    DuneCrestMin_h3m,
                ],
                DuneCrestMax=[
                    DuneCrestMax_nat,
                    DuneCrestMax_h1m,
                    DuneCrestMax_h2m,
                    DuneCrestMax_h3m,
                ],
                BarrierHeight=[
                    BarrierHeight_nat,
                    BarrierHeight_h1m,
                    BarrierHeight_h2m,
                    BarrierHeight_h3m,
                ],
                BarrierWidth=[
                    BarrierWidth_nat,
                    BarrierWidth_h1m,
                    BarrierWidth_h2m,
                    BarrierWidth_h3m,
                ],
                DuneCrestMean=[
                    DuneCrestMean_nat,
                    DuneCrestMean_h1m,
                    DuneCrestMean_h2m,
                    DuneCrestMean_h3m,
                ],
                TMAX=[750, 750, 750, 532],
                tmax_management=[0, 730, 630, 521],
            )

        # supplementary material
        def sensitivity_abandonment_criteria():
            (
                BarrierWidth_nat,
                DuneCrestMean_nat,
                BarrierHeight_nat,
                bw_rate_nat,
                bh_rate_nat,
                sc_rate_nat,
                DuneCrestMin_nat,
                DuneCrestMax_nat,
                cascade_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low",
                tmin=0,
                tmax_roadways=500,  # dummy
                tmax_sim=500,
                plot_name="b3d_pt75_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_10,
                DuneCrestMean_10,
                BarrierHeight_10,
                bw_rate_10,
                bh_rate_10,
                sc_rate_10,
                DuneCrestMin_10,
                DuneCrestMax_10,
                cascade_10,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_10percent",
                tmin=0,
                tmax_roadways=158,
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_10percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_20,
                DuneCrestMean_20,
                BarrierHeight_20,
                bw_rate_20,
                bh_rate_20,
                sc_rate_20,
                DuneCrestMin_20,
                DuneCrestMax_20,
                cascade_20,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_20percent",
                tmin=0,
                tmax_roadways=161,  # drowned at 162
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_20percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_30,
                DuneCrestMean_30,
                BarrierHeight_30,
                bw_rate_30,
                bh_rate_30,
                sc_rate_30,
                DuneCrestMin_30,
                DuneCrestMax_30,
                cascade_30,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_30percent",
                tmin=0,
                tmax_roadways=164,  # drowned at 165
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_30percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_40,
                DuneCrestMean_40,
                BarrierHeight_40,
                bw_rate_40,
                bh_rate_40,
                sc_rate_40,
                DuneCrestMin_40,
                DuneCrestMax_40,
                cascade_40,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_40percent",
                tmin=0,
                tmax_roadways=172,  # drowned at 173
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_40percent",
                run_road_mgmt=True,
            )

            (
                BarrierWidth_50,
                DuneCrestMean_50,
                BarrierHeight_50,
                bw_rate_50,
                bh_rate_50,
                sc_rate_50,
                DuneCrestMin_50,
                DuneCrestMax_50,
                cascade_50,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Roadways_2mDune_20mSetback_20mWidth_low_50percent",
                tmin=0,
                tmax_roadways=181,  # drowned at 182
                tmax_sim=500,
                plot_name="b3d_pt75_h2m_plots_low_50percent",
                run_road_mgmt=True,
            )

            CASCADEplt.plot_nonlinear_stats_low_high_sensitivity(
                cascade=[
                    cascade_10,
                    cascade_20,
                    cascade_30,
                    cascade_40,
                    cascade_50,
                ],
                BarrierHeight=[
                    BarrierHeight_10,
                    BarrierHeight_20,
                    BarrierHeight_30,
                    BarrierHeight_40,
                    BarrierHeight_50,
                ],
                BarrierWidth=[
                    BarrierWidth_10,
                    BarrierWidth_20,
                    BarrierWidth_30,
                    BarrierWidth_40,
                    BarrierWidth_50,
                ],
                TMAX=[500, 500, 500, 500, 500],
                tmax_roadways=[158, 161, 164, 172, 181],
            )

    def slr_sensitivity():
        def pt45_high():
            (
                BarrierWidth_pt45_high_SLRacc,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_high_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_high_AccSLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_high_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_high_0pt012SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_high_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_high_0pt012SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_high_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_high_0pt008SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_high_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_high_0pt008SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_high_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_high_0pt004SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_high_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_high",  # this is the 0.004 case
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_high_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt45_high_0pt004SLR,
                cascade_pt45_high_0pt008SLR,
                cascade_pt45_high_0pt012SLR,
                cascade_pt45_high_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

        def pt45_low():
            (
                BarrierWidth_pt45_low_SLRacc,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_low_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_low_AccSLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_low_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_0pt012SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_low_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_low_0pt012SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_low_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_0pt008SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_low_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_low_0pt008SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_low_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_0pt004SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt45_low_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_low",  # this is the 0.004 case
                tmin=0,
                tmax=200,
                plot_name="b3d_pt45_Natural_low_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt45_low_0pt004SLR,
                cascade_pt45_low_0pt008SLR,
                cascade_pt45_low_0pt012SLR,
                cascade_pt45_low_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

        def pt75_low():
            (
                BarrierWidth_pt75_low_SLRacc,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_low_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low_AccSLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_low_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_0pt012SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_low_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low_0pt012SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_low_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_0pt008SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_low_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low_0pt008SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_low_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_0pt004SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_low_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low",  # this is the 0.004 case
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_low_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt75_low_0pt004SLR,
                cascade_pt75_low_0pt008SLR,
                cascade_pt75_low_0pt012SLR,
                cascade_pt75_low_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

        def pt75_high():
            (
                BarrierWidth_pt75_high_SLRacc,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_high_SLRacc,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_high_AccSLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_high_AccSLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_high_0pt012SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_high_0pt012SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_high_0pt012SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_high_0pt012SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_high_0pt008SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_high_0pt008SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_high_0pt008SLR",
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_high_0pt008SLR",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_high_0pt004SLR,
                _,
                _,
                _,
                _,
                _,
                _,
                _,
                cascade_pt75_high_0pt004SLR,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_high",  # this is the 0.004 case
                tmin=0,
                tmax=200,
                plot_name="b3d_pt75_Natural_high_0pt004SLR",
                run_road_mgmt=False,
            )

            cascade = [
                cascade_pt75_high_0pt004SLR,
                cascade_pt75_high_0pt008SLR,
                cascade_pt75_high_0pt012SLR,
                cascade_pt75_high_SLRacc,
            ]
            TMAX = [200, 200, 200, 200]
            CASCADEplt.fig3_slr_sensitivity(
                cascade,  # lists
                TMAX,
            )

    def nourishments():
        def pt75_low():
            (
                BarrierWidth_pt75_nat,
                DuneCrestMean_pt75_nat,
                BarrierHeight_pt75_nat,
                bw_rate_pt75_nat,
                bh_rate_pt75_nat,
                sc_rate_pt75_nat,
                DuneCrestMin_pt75_nat,
                DuneCrestMax_pt75_nat,
                cascade_pt75_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt75_Natural_low",
                tmin=0,
                tmax_roadways=500,  # dummy
                tmax_sim=500,
                plot_name="b3d_pt75_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt75_low_40pc,
                DuneCrestMean_pt75_low_40pc,
                BarrierHeight_pt75_low_40pc,
                bh_rate_pt75_low_40pc,
                bw_rate_pt75_low_40pc,
                sc_rate_pt75_low_40pc,
                DuneCrestMin_pt75_low_40pc,
                DuneCrestMax_pt75_low_40pc,
                shoreline_position_pt75_low_40pc,
                shoreface_slope_pt75_low_40pc,
                beach_width_pt75_low_40pc,
                overwash_pt75_low_40pc,
                dune_toe_pt75_low_40pc,
                cascade_pt75_low_40pc,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_residential",
                tmax_management=178,  # drowned at 178
                tmax_sim=500,
                plot_name="b3d_pt75_Nourishment_2mDune_lowEle_residential",
            )

            (
                BarrierWidth_pt75_low_90pc,
                DuneCrestMean_pt75_low_90pc,
                BarrierHeight_pt75_low_90pc,
                bh_rate_pt75_low_90pc,
                bw_rate_pt75_low_90pc,
                sc_rate_pt75_low_90pc,
                DuneCrestMin_pt75_low_90pc,
                DuneCrestMax_pt75_low_90pc,
                shoreline_position_pt75_low_90pc,
                shoreface_slope_pt75_low_90pc,
                beach_width_pt75_low_90pc,
                overwash_pt75_low_90pc,
                dune_toe_pt75_low_90pc,
                cascade_pt75_low_90pc,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial",
                tmax_management=90,  # community drowned at 90 years
                tmax_sim=500,
                plot_name="b3d_pt75_Nourishment_2mDune_lowEle_commercial",
            )

            (
                BarrierWidth_pt75_low_90pc_backerosion_pt25m,
                DuneCrestMean_pt75_low_90pc_backerosion_pt25m,
                BarrierHeight_pt75_low_90pc_backerosion_pt25m,
                bh_rate_pt75_low_90pc_backerosion_pt25m,
                bw_rate_pt75_low_90pc_backerosion_pt25m,
                sc_rate_pt75_low_90pc_backerosion_pt25m,
                DuneCrestMin_pt75_low_90pc_backerosion_pt25m,
                DuneCrestMax_pt75_low_90pc_backerosion_pt25m,
                shoreline_position_pt75_low_90pc_backerosion_pt25m,
                shoreface_slope_pt75_low_90pc_backerosion_pt25m,
                beach_width_pt75_low_90pc_backerosion_pt25m,
                overwash_pt75_low_90pc_backerosion_pt25m,
                dune_toe_pt75_low_90pc_backerosion_pt25m,
                cascade_pt75_low_90pc_backerosion_pt25m,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
                tmax_management=90,  # community drowned at 90 years
                tmax_sim=500,
                plot_name="b3d_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
            )

            (
                BarrierWidth_pt75_low_90pc_backerosion_1m,
                DuneCrestMean_pt75_low_90pc_backerosion_1m,
                BarrierHeight_pt75_low_90pc_backerosion_1m,
                bh_rate_pt75_low_90pc_backerosion_1m,
                bw_rate_pt75_low_90pc_backerosion_1m,
                sc_rate_pt75_low_90pc_backerosion_1m,
                DuneCrestMin_pt75_low_90pc_backerosion_1m,
                DuneCrestMax_pt75_low_90pc_backerosion_1m,
                shoreline_position_pt75_low_90pc_backerosion_1m,
                shoreface_slope_pt75_low_90pc_backerosion_1m,
                beach_width_pt75_low_90pc_backerosion_1m,
                overwash_pt75_low_90pc_backerosion_1m,
                dune_toe_pt75_low_90pc_backerosion_1m,
                cascade_pt75_low_90pc_backerosion_1m,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
                tmax_management=90,  # community drowned at 90 years
                tmax_sim=500,
                plot_name="b3d_pt75_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
            )

            rebuild_threshold = 0.3 + (cascade_pt75_low_40pc.barrier3d[0].BermEl * 10)

            def version2_background_erosion():
                CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                    cascade=[
                        cascade_pt75_nat,
                        # cascade_pt75_low_40pc,  # only residential here
                        cascade_pt75_low_90pc,
                        cascade_pt75_low_90pc_backerosion_pt25m,
                        cascade_pt75_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMin=[
                        DuneCrestMin_pt75_nat,
                        # DuneCrestMin_pt75_low_40pc,
                        DuneCrestMin_pt75_low_90pc,
                        DuneCrestMin_pt75_low_90pc_backerosion_pt25m,
                        DuneCrestMin_pt75_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMax=[
                        DuneCrestMax_pt75_nat,
                        # DuneCrestMax_pt75_low_40pc,
                        DuneCrestMax_pt75_low_90pc,
                        DuneCrestMax_pt75_low_90pc_backerosion_pt25m,
                        DuneCrestMax_pt75_low_90pc_backerosion_1m,
                    ],
                    BarrierHeight=[
                        BarrierHeight_pt75_nat,
                        # BarrierHeight_pt75_low_40pc,
                        BarrierHeight_pt75_low_90pc,
                        BarrierHeight_pt75_low_90pc_backerosion_pt25m,
                        BarrierHeight_pt75_low_90pc_backerosion_1m,
                    ],
                    BarrierWidth=[
                        BarrierWidth_pt75_nat,
                        # BarrierWidth_pt75_low_40pc,
                        BarrierWidth_pt75_low_90pc,
                        BarrierWidth_pt75_low_90pc_backerosion_pt25m,
                        BarrierWidth_pt75_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMean=[
                        DuneCrestMean_pt75_nat,
                        # DuneCrestMean_pt75_low_40pc,
                        DuneCrestMean_pt75_low_90pc,
                        DuneCrestMean_pt75_low_90pc_backerosion_pt25m,
                        DuneCrestMean_pt75_low_90pc_backerosion_1m,
                    ],
                    shoreline_position=[
                        cascade_pt75_nat.barrier3d[
                            0
                        ].x_s_TS,  # wasn't saved yet in roadways, update later
                        # shoreline_position_pt75_low_40pc,
                        shoreline_position_pt75_low_90pc,
                        shoreline_position_pt75_low_90pc_backerosion_pt25m,
                        shoreline_position_pt75_low_90pc_backerosion_1m,
                    ],
                    overwash=[
                        cascade_pt75_nat.barrier3d[
                            0
                        ].QowTS,  # wasn't saved yet in roadways, update later
                        # overwash_pt75_low_40pc,
                        overwash_pt75_low_90pc,
                        overwash_pt75_low_90pc_backerosion_pt25m,
                        overwash_pt75_low_90pc_backerosion_1m,
                    ],
                    dune_toe=[
                        [0],  # dummy
                        # dune_toe_pt75_low_40pc,
                        dune_toe_pt75_low_90pc,
                        dune_toe_pt75_low_90pc_backerosion_pt25m,
                        dune_toe_pt75_low_90pc_backerosion_1m,
                    ],
                    TMAX=[
                        500,
                        # 500,
                        500,
                        500,
                        500,
                    ],
                    tmax_management=[
                        500,  # dummy
                        # 178,
                        90,
                        90,
                        90,
                    ],
                    roadways_on=False,
                    nourishment_on=True,
                    rebuild_threshold=rebuild_threshold,  # min dune height above the berm [m MHW], same as in RoadwayManager
                    scenarios=[
                        "natural",
                        "commercial",
                        "comm, 0.25 m/yr",
                        "comm, 1 m/yr",
                    ],
                )

        def pt45_low():
            (
                BarrierWidth_pt45_nat,
                DuneCrestMean_pt45_nat,
                BarrierHeight_pt45_nat,
                bw_rate_pt45_nat,
                bh_rate_pt45_nat,
                sc_rate_pt45_nat,
                DuneCrestMin_pt45_nat,
                DuneCrestMax_pt45_nat,
                cascade_pt45_nat,
            ) = PLOT_5_Nonlinear_Dynamics_CASCADE_B3Donly_RoadwayManager(
                name_prefix="6-B3D_Rave_pt45_Natural_low",
                tmin=0,
                tmax_roadways=1000,
                tmax_sim=1000,
                plot_name="b3d_pt45_plots_low",
                run_road_mgmt=False,
            )

            (
                BarrierWidth_pt45_low_40pc,
                DuneCrestMean_pt45_low_40pc,
                BarrierHeight_pt45_low_40pc,
                bh_rate_pt45_low_40pc,
                bw_rate_pt45_low_40pc,
                sc_rate_pt45_low_40pc,
                DuneCrestMin_pt45_low_40pc,
                DuneCrestMax_pt45_low_40pc,
                shoreline_position_pt45_low_40pc,
                shoreface_slope_pt45_low_40pc,
                beach_width_pt45_low_40pc,
                overwash_pt45_low_40pc,
                dune_toe_pt45_low_40pc,
                cascade_pt45_low_40pc,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_residential",
                tmax_management=496,  # drowned at 496
                tmax_sim=500,
                plot_name="b3d_pt45_Nourishment_2mDune_lowEle_residential",
            )

            (
                BarrierWidth_pt45_low_90pc,
                DuneCrestMean_pt45_low_90pc,
                BarrierHeight_pt45_low_90pc,
                bh_rate_pt45_low_90pc,
                bw_rate_pt45_low_90pc,
                sc_rate_pt45_low_90pc,
                DuneCrestMin_pt45_low_90pc,
                DuneCrestMax_pt45_low_90pc,
                shoreline_position_pt45_low_90pc,
                shoreface_slope_pt45_low_90pc,
                beach_width_pt45_low_90pc,
                overwash_pt45_low_90pc,
                dune_toe_pt45_low_90pc,
                cascade_pt45_low_90pc,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial",
                tmax_management=421,  # community drowned at 421
                tmax_sim=458,  # barrier height drowned at 458
                plot_name="b3d_pt45_Nourishment_2mDune_lowEle_commercial",
            )

            (
                BarrierWidth_pt45_low_90pc_backerosion_pt25m,
                DuneCrestMean_pt45_low_90pc_backerosion_pt25m,
                BarrierHeight_pt45_low_90pc_backerosion_pt25m,
                bh_rate_pt45_low_90pc_backerosion_pt25m,
                bw_rate_pt45_low_90pc_backerosion_pt25m,
                sc_rate_pt45_low_90pc_backerosion_pt25m,
                DuneCrestMin_pt45_low_90pc_backerosion_pt25m,
                DuneCrestMax_pt45_low_90pc_backerosion_pt25m,
                shoreline_position_pt45_low_90pc_backerosion_pt25m,
                shoreface_slope_pt45_low_90pc_backerosion_pt25m,
                beach_width_pt45_low_90pc_backerosion_pt25m,
                overwash_pt45_low_90pc_backerosion_pt25m,
                dune_toe_pt45_low_90pc_backerosion_pt25m,
                cascade_pt45_low_90pc_backerosion_pt25m,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
                tmax_management=421,  # community drowned at 421
                tmax_sim=454,  # barrier height drowned at 454
                plot_name="b3d_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_pt25m",
            )

            (
                BarrierWidth_pt45_low_90pc_backerosion_1m,
                DuneCrestMean_pt45_low_90pc_backerosion_1m,
                BarrierHeight_pt45_low_90pc_backerosion_1m,
                bh_rate_pt45_low_90pc_backerosion_1m,
                bw_rate_pt45_low_90pc_backerosion_1m,
                sc_rate_pt45_low_90pc_backerosion_1m,
                DuneCrestMin_pt45_low_90pc_backerosion_1m,
                DuneCrestMax_pt45_low_90pc_backerosion_1m,
                shoreline_position_pt45_low_90pc_backerosion_1m,
                shoreface_slope_pt45_low_90pc_backerosion_1m,
                beach_width_pt45_low_90pc_backerosion_1m,
                overwash_pt45_low_90pc_backerosion_1m,
                dune_toe_pt45_low_90pc_backerosion_1m,
                cascade_pt45_low_90pc_backerosion_1m,
            ) = PLOT_6_Nonlinear_Dynamics_CASCADE_B3Donly_Nourishments(
                name_prefix="8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
                tmax_management=421,  # community drowned at 421
                tmax_sim=500,
                plot_name="b3d_pt45_Nourishment_2mDune_lowEle_commercial_backerosion_1m",
            )

            rebuild_threshold = 0.3 + (cascade_pt45_low_40pc.barrier3d[0].BermEl * 10)

            def version1_residential():
                CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                    cascade=[
                        cascade_pt45_nat,
                        cascade_pt45_low_40pc,  # only residential here
                        cascade_pt45_low_90pc,
                        # cascade_pt45_low_90pc_backerosion_pt25m,
                        cascade_pt45_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMin=[
                        DuneCrestMin_pt45_nat,
                        DuneCrestMin_pt45_low_40pc,
                        DuneCrestMin_pt45_low_90pc,
                        # DuneCrestMin_pt45_low_90pc_backerosion_pt25m,
                        DuneCrestMin_pt45_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMax=[
                        DuneCrestMax_pt45_nat,
                        DuneCrestMax_pt45_low_40pc,
                        DuneCrestMax_pt45_low_90pc,
                        # DuneCrestMax_pt45_low_90pc_backerosion_pt25m,
                        DuneCrestMax_pt45_low_90pc_backerosion_1m,
                    ],
                    BarrierHeight=[
                        BarrierHeight_pt45_nat,
                        BarrierHeight_pt45_low_40pc,
                        BarrierHeight_pt45_low_90pc,
                        # BarrierHeight_pt45_low_90pc_backerosion_pt25m,
                        BarrierHeight_pt45_low_90pc_backerosion_1m,
                    ],
                    BarrierWidth=[
                        BarrierWidth_pt45_nat,
                        BarrierWidth_pt45_low_40pc,
                        BarrierWidth_pt45_low_90pc,
                        # BarrierWidth_pt45_low_90pc_backerosion_pt25m,
                        BarrierWidth_pt45_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMean=[
                        DuneCrestMean_pt45_nat,
                        DuneCrestMean_pt45_low_40pc,
                        DuneCrestMean_pt45_low_90pc,
                        # DuneCrestMean_pt45_low_90pc_backerosion_pt25m,
                        DuneCrestMean_pt45_low_90pc_backerosion_1m,
                    ],
                    shoreline_position=[
                        cascade_pt45_nat.barrier3d[
                            0
                        ].x_s_TS,  # wasn't saved yet in roadways, update later
                        shoreline_position_pt45_low_40pc,
                        shoreline_position_pt45_low_90pc,
                        # shoreline_position_pt45_low_90pc_backerosion_pt25m,
                        shoreline_position_pt45_low_90pc_backerosion_1m,
                    ],
                    overwash=[
                        cascade_pt45_nat.barrier3d[
                            0
                        ].QowTS,  # wasn't saved yet in roadways, update later
                        overwash_pt45_low_40pc,
                        overwash_pt45_low_90pc,
                        # overwash_pt45_low_90pc_backerosion_pt25m,
                        overwash_pt45_low_90pc_backerosion_1m,
                    ],
                    dune_toe=[
                        [0],  # dummy
                        dune_toe_pt45_low_40pc,
                        dune_toe_pt45_low_90pc,
                        # dune_toe_pt45_low_90pc_backerosion_pt25m,
                        dune_toe_pt45_low_90pc_backerosion_1m,
                    ],
                    TMAX=[
                        500,
                        500,
                        458,
                        # 454,
                        500,
                    ],
                    tmax_management=[
                        500,  # dummy
                        496,
                        421,
                        # 421,
                        421,
                    ],
                    roadways_on=False,
                    nourishment_on=True,
                    rebuild_threshold=rebuild_threshold,
                    # min dune height above the berm [m MHW], same as in RoadwayManager
                    scenarios=["natural", "residential", "commercial", "comm, 1 m/yr"],
                )

            def version2_background_erosion():
                CASCADEplt.plot_nonlinear_stats_mgmt_array4(
                    cascade=[
                        cascade_pt45_nat,
                        # cascade_pt45_low_40pc,  # only residential here
                        cascade_pt45_low_90pc,
                        cascade_pt45_low_90pc_backerosion_pt25m,
                        cascade_pt45_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMin=[
                        DuneCrestMin_pt45_nat,
                        # DuneCrestMin_pt45_low_40pc,
                        DuneCrestMin_pt45_low_90pc,
                        DuneCrestMin_pt45_low_90pc_backerosion_pt25m,
                        DuneCrestMin_pt45_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMax=[
                        DuneCrestMax_pt45_nat,
                        # DuneCrestMax_pt45_low_40pc,
                        DuneCrestMax_pt45_low_90pc,
                        DuneCrestMax_pt45_low_90pc_backerosion_pt25m,
                        DuneCrestMax_pt45_low_90pc_backerosion_1m,
                    ],
                    BarrierHeight=[
                        BarrierHeight_pt45_nat,
                        # BarrierHeight_pt45_low_40pc,
                        BarrierHeight_pt45_low_90pc,
                        BarrierHeight_pt45_low_90pc_backerosion_pt25m,
                        BarrierHeight_pt45_low_90pc_backerosion_1m,
                    ],
                    BarrierWidth=[
                        BarrierWidth_pt45_nat,
                        # BarrierWidth_pt45_low_40pc,
                        BarrierWidth_pt45_low_90pc,
                        BarrierWidth_pt45_low_90pc_backerosion_pt25m,
                        BarrierWidth_pt45_low_90pc_backerosion_1m,
                    ],
                    DuneCrestMean=[
                        DuneCrestMean_pt45_nat,
                        # DuneCrestMean_pt45_low_40pc,
                        DuneCrestMean_pt45_low_90pc,
                        DuneCrestMean_pt45_low_90pc_backerosion_pt25m,
                        DuneCrestMean_pt45_low_90pc_backerosion_1m,
                    ],
                    shoreline_position=[
                        cascade_pt45_nat.barrier3d[
                            0
                        ].x_s_TS,  # wasn't saved yet in roadways, update later
                        # shoreline_position_pt45_low_40pc,
                        shoreline_position_pt45_low_90pc,
                        shoreline_position_pt45_low_90pc_backerosion_pt25m,
                        shoreline_position_pt45_low_90pc_backerosion_1m,
                    ],
                    overwash=[
                        cascade_pt45_nat.barrier3d[
                            0
                        ].QowTS,  # wasn't saved yet in roadways, update later
                        # overwash_pt45_low_40pc,
                        overwash_pt45_low_90pc,
                        overwash_pt45_low_90pc_backerosion_pt25m,
                        overwash_pt45_low_90pc_backerosion_1m,
                    ],
                    dune_toe=[
                        [0],  # dummy
                        # dune_toe_pt45_low_40pc,
                        dune_toe_pt45_low_90pc,
                        dune_toe_pt45_low_90pc_backerosion_pt25m,
                        dune_toe_pt45_low_90pc_backerosion_1m,
                    ],
                    TMAX=[
                        500,
                        # 500,
                        458,
                        454,
                        500,
                    ],
                    tmax_management=[
                        500,  # dummy
                        # 496,
                        421,
                        421,
                        421,
                    ],
                    roadways_on=False,
                    nourishment_on=True,
                    rebuild_threshold=rebuild_threshold,  # min dune height above the berm [m MHW], same as in RoadwayManager
                    scenarios=[
                        "natural",
                        "commercial",
                        "comm, 0.25 m/yr",
                        "comm, 1 m/yr",
                    ],
                )

    def initial_topos():

        PLOT_7_Initial_CNH_Topographies(
            [
                "8-B3D_Rave_pt45_Nourishment_2mDune_lowEle_commercial",
                "8-B3D_Rave_pt45_Nourishment_2mDune_highEle_commercial",
                "8-B3D_Rave_pt75_Nourishment_2mDune_lowEle_commercial",
                "8-B3D_Rave_pt75_Nourishment_2mDune_highEle_commercial",
            ]
        )
