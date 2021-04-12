# run file for

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
----------------------------------------------------"""

# remember if I move to a different computer to $ pip install -e . in the brie and B3D directories for the BMI

import numpy as np
import os
import multiprocessing
import time

from scripts import CASCADE_plotters as CASCADEplt

from CASCADE import Cascade  # the new class

# import CASCADE as CASCADE

from Tools.Barrier3D_MakeTimeSeries import (
    storms_per_year_from_MSSM_output,
    gen_dune_height_start,
    gen_alongshore_variable_rmin_rmax,
)

# for laptop and desktop, use all but one core; on supercomputer, use all cores
num_cores = multiprocessing.cpu_count() - 1

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
    ny = 12  # number of alongshore sections (12 = 6 km)
    nt = 1000  # timesteps for 1000 morphologic years
    # rave = [0.45]  # to help me remember the average

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


def RUN_6_B3D_Rave_SLR_pt004_Humans(
    nt,
    rmin,
    rmax,
    name,
    road_ele,
    road_width,
    road_setback,
    artificial_max_dune_ele,
    artificial_min_dune_ele,
    run_road_mgmt,
):

    # ###############################################################################
    # 6 - B3D with human modifications
    # ###############################################################################
    # GOAL: Use the starting interior domain from the 10,000 yr runs for each dune growth rate and run for 200 years
    # with and without human modifications. NOTE, currently have to hard-code elevations and storm filenames, but just
    # so I remember, I decided to use the 3000 yr storm time series.

    # --------- INITIALIZE ---------
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    cascade = Cascade(
        datadir,
        name,
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        alongshore_section_count=1,
        time_step_count=nt,
        min_dune_growth_rate=rmin,
        max_dune_growth_rate=rmax,
        num_cores=5,
        roadway_management_module=run_road_mgmt,
        alongshore_transport_module=False,  # no brie coupling
        road_ele=road_ele,
        road_width=road_width,
        road_setback=road_setback,
        artificial_max_dune_ele=artificial_max_dune_ele,
        artificial_min_dune_ele=artificial_min_dune_ele,
    )

    # --------- LOOP ---------
    Time = time.time()

    for time_step in range(nt - 1):
        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", "Time Step: ", time_step, end="")
        cascade.update()

    # --------- SAVE ---------
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
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
    os.chdir("/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output")
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
    # output = np.load(
    #     "4-B3D_Rave_pt75_SLR_pt004.npz", allow_pickle=True
    # )  # B3D high growth rate run, no AST
    # b3d_pt75 = output["barrier3d"]

    tmin = 9000  # 500
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
    ) = CASCADEplt.plot_nonlinear_stats(b3d_pt75, ib3d, tmin, tmax)

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
        rmin=0.35, rmax=0.75, name="4-B3D_Rave_pt55_SLR_pt004"
    )  # rave = 0.55
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.45, rmax=0.85, name="4-B3D_Rave_pt65_SLR_pt004"
    )  # rave = 0.65
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.55, rmax=0.95, name="4-B3D_Rave_pt75_SLR_pt004"
    )  # rave = 0.75
    Periodicity, AvgFastDur, AvgSlowDur, Punc = RUN_4_B3D_Rave_SLR_pt004(
        rmin=0.55, rmax=0.95, name="4-B3D_Rave_pt75_SLR_pt004_10k-yrs"
    )  # rave = 0.75


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


# record of B3D time series -------------------------------------------------------------------------------------
def time_series():
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"

    name = "StormTimeSeries_10k-yr.npy"
    storms_per_year_from_MSSM_output(
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
    storms_per_year_from_MSSM_output(
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
    storms_per_year_from_MSSM_output(
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
    storms_per_year_from_MSSM_output(
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

    b3d_pt45 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=200,
        rmin=0.25,
        rmax=0.65,  # rave = 0.45
        name="6-B3D_Rave_pt45_200yrs_Natural_v2",
        road_ele=None,
        road_width=None,
        road_setback=None,
        artificial_max_dune_ele=None,
        artificial_min_dune_ele=None,
    )

    # v1 - didn't drown roadway, v2 - roadway drowned at 160, v3 - new class
    # b3d_pt45_h1, dunes_rebuilt, road_overwash_volume = RUN_6_B3D_Rave_SLR_pt004_Humans(
    cascade_pt45_h1 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=159,  # 200
        rmin=0.25,
        rmax=0.65,  # rave = 0.45
        name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_20mWidth_v3",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=20,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=3.7,  # m NAVD88, rebuild to 2 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
        run_road_mgmt=True,
    )

    b3d_pt45_h2 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=200,
        rmin=0.25,
        rmax=0.65,  # rave = 0.45
        name="6-B3D_Rave_pt45_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=20,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=4.7,  # m NAVD88, rebuild to 3 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 3.7 m
    )

    # v1 drowned at 91 years, v2 - roadway drowned at 101 years
    # b3d_pt45_h3 = RUN_6_B3D_Rave_SLR_pt004_Humans(
    cascade_pt45_h3 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=100,  # 90
        rmin=0.25,
        rmax=0.65,  # rave = 0.45
        name="6-B3D_Rave_pt45_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2_classtest2",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=20,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=2.7,  # m NAVD88, rebuild to 1 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 2.4 m
        run_road_mgmt=True,
    )

    b3d_pt45_h4 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=200,
        rmin=0.25,
        rmax=0.65,  # rave = 0.45
        name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_40mSetback_30mWidth_v2",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=30,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=3.7,  # m NAVD88, rebuild to 2 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
    )

    # # v1 drowned at 183 years
    # b3d_pt45_h5 = RUN_6_B3D_Rave_SLR_pt004_Humans(
    #     nt=200,  # 182
    #     rmin=0.25,
    #     rmax=0.65,  # rave = 0.45
    #     name="6-B3D_Rave_pt45_200yrs_Roadways_2mDune_30mSetback_20mWidth",
    #     road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
    #     road_width=20,  # m
    #     road_setback=30,  # m
    #     artificial_max_dune_ele=3.7,  # m NAVD88, 2 m dune above the roadway
    #     artificial_min_dune_ele=2.7,  # m NAVD88, 1 m dune above the roadway
    # )

    # REMEMBER TO SWITCH TOPOGRAHPHY FILES
    b3d_pt75 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=200,
        rmin=0.55,
        rmax=0.95,  # rave = 0.75
        name="6-B3D_Rave_pt75_200yrs_Natural",
        road_ele=None,
        road_width=None,
        road_setback=None,
        artificial_max_dune_ele=None,
        artificial_min_dune_ele=None,
    )

    # v2, roadway drowned at 157 years
    b3d_pt75_h1, dunes_rebuilt, road_overwash_volume = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=156,  # 200
        rmin=0.55,
        rmax=0.95,  # rave = 0.75
        name="6-B3D_Rave_pt75_200yrs_Roadways_2mDune_40mSetback_20mWidth_v2",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=20,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=3.7,  # m NAVD88, rebuild to 2 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
    )

    b3d_pt75_h2 = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=200,
        rmin=0.55,
        rmax=0.95,  # rave = 0.75
        name="6-B3D_Rave_pt75_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=20,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=4.7,  # m NAVD88, rebuild to 3 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 3.7 m
    )

    # v1 drowned at 88 years, v2 drowned at 105 years
    b3d_pt75_h3, dunes_rebuilt, road_overwash_volume = RUN_6_B3D_Rave_SLR_pt004_Humans(
        nt=104,  # 87
        rmin=0.55,
        rmax=0.95,  # rave = 0.75
        name="6-B3D_Rave_pt75_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2",
        road_ele=1.7,  # 1.7 m NAVD88 for pt45, 1.5 m NAVD88 for pt75 (average of NC-12 is 1.3 m NAVD88), berm ele is 1.4 m NAV
        road_width=20,  # m
        road_setback=40,  # m
        artificial_max_dune_ele=2.7,  # m NAVD88, rebuild to 1 m dune above the roadway
        artificial_min_dune_ele=2.2,  # m NAVD88, allow dune to erode down to 0.5 m above the roadway, v1 = 2.4 m
    )


# record of human plots -------------------------------------------------------------------------------------
def human_plots():
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
        plot_name="b3d_pt45_h1_plots_v2",
    )
    PLOT_5_Nonlinear_Dynamics_B3D_CNH(
        name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_3mDune_40mSetback_20mWidth",
        tmin=0,
        tmax=200,
        plot_name="b3d_pt45_h2_plots",
    )
    PLOT_5_Nonlinear_Dynamics_B3D_CNH(
        name_prefix="6-B3D_Rave_pt45_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2",
        tmin=0,
        tmax=100,  # 90
        plot_name="b3d_pt45_h3_plots_v2",
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
    PLOT_5_Nonlinear_Dynamics_B3D_CNH(
        name_prefix="6-B3D_Rave_pt75_200yrs_Natural",
        tmin=0,
        tmax=200,
        plot_name="b3d_pt75_plots",
    )
    PLOT_5_Nonlinear_Dynamics_B3D_CNH(
        name_prefix="6-B3D_Rave_pt75_200yrs_Roadways_2mDune_40mSetback_20mWidth_v2",
        tmin=0,
        tmax=156,  # 200
        plot_name="b3d_pt75_h1_plots_v2",
    )
    PLOT_5_Nonlinear_Dynamics_B3D_CNH(
        name_prefix="6-B3D_Rave_pt75_200yrs_Roadways_3mDune_40mSetback_20mWidth_v2",
        tmin=0,
        tmax=200,
        plot_name="b3d_pt75_h2_plots_v2",
    )
    PLOT_5_Nonlinear_Dynamics_B3D_CNH(
        name_prefix="6-B3D_Rave_pt75_200yrs_Roadways_1mDune_40mSetback_20mWidth_v2",
        tmin=0,
        tmax=104,  # 87
        plot_name="b3d_pt75_h3_plots_v2",
    )
