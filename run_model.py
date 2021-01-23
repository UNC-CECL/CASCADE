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

import CASCADE_plotters as CASCADEplt

import CASCADE as CASCADE

# for laptop and desktop, use all but one core; on supercomputer, use all cores
num_cores = multiprocessing.cpu_count() - 1

def RUN_1_CASCADE_LTA_COMPARISON():
    # ###############################################################################
    # 1 - CASCADE_LTA_COMPARISON
    # ###############################################################################
    # GOAL: highlight different processes in models with alongshore homogenous dune line, 3000 year simulation
    #
    # --------- INITIAL CONDITIONS ---------
    name = '1-CASCADE_LTA_COMPARISON_3km_3000yr'
    #name = '1-CASCADE_LTA_COMPARISON_6km_1500yr'
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.004  # m/yr
    ny = 6  #12 # number of alongshore sections (6=3 km for 3000 yr run, 12=6 km for 1500 yr run)
    nt = 3000  # 3000  #1500 # timesteps for 3000 morphologic years
    rmin = 0.35  # minimum growth rate for logistic dune growth (can be a list)
    rmax = 0.85  # maximum growth rate for logistic dune growth (can be a list)

    # --------- INITIALIZE ---------
    # datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/" # iMAC
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
    brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)

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
    brieLTA = CASCADE.LTA(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, w_b_crit, h_b_crit, Qow_max)

    # --------- SAVE ---------
    # #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    CASCADE.save(brie, barrier3d, save_directory, name) # this returns the barrier3d model without the BMI

    return brieLTA

def PLOT_1_CASCADE_LTA_COMPARISON(brieLTA, name, save_directory):
    # --------- plot ---------
    # load the simulation if previously saved
    os.chdir('/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/')
    filename = name + '.npz'
    output = np.load(filename, allow_pickle=True)
    b3d = output['barrier3d']
    # brie = output['brie']
    # brie = brie[0]

    # 1: Animation Frames of Barrier and Dune Elevation
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    ny = len(b3d)
    CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX, name)

    #===================================================

    # 4: Cross-shore transects for both brieLTA and B3d
    iB3D = 0
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    time_step = [0, int(TMAX/2), TMAX-2]
    CASCADEplt.plot_ModelTransects(b3d, brieLTA, time_step, iB3D)

    #===================================================

    # 5: Statistics from B3d
    TMAX = b3d[0].time_index - 1
    iB3D = 1
    CASCADEplt.plot_statistics(b3d, iB3D, TMAX)

    # 6: Statistics from BrieLTA
    TMAX = int((b3d[0].time_index - 1) / brieLTA._dt)
    iB3D = 1
    iB3D_BRIE = iB3D * int( (b3d[0]._BarrierLength * 10) / brieLTA._dy)
    CASCADEplt.plot_statistics_BRIE(brieLTA, iB3D, TMAX)

    # 6: Statistics from both models for AGU presentation
    iB3D = 1
    TMAX = b3d[0].time_index - 1
    iBRIE = iB3D * int( (b3d[0]._BarrierLength * 10) / brieLTA._dy)
    TMAX_BRIE = int((b3d[0].time_index - 1) / brieLTA._dt)
    CASCADEplt.plot_statisticsAGU(b3d, brieLTA, iB3D, iBRIE, TMAX, TMAX_BRIE)

    #===================================================

    # 2: Shoreline positions over time (#6 in Barrier3D_Functions)
    TMAX = 1001
    CASCADEplt.plot_ShorelinePositions(b3d[0]._x_s_TS[0:TMAX], b3d[0]._x_b_TS[0:TMAX])

    # 3: Cross-Shore Transect for one subgrid every 100 m for last time step
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    CASCADEplt.plot_XShoreTransects(b3d[0], TMAX)

def RUN_2_AlongshoreVarGrowthParam_Alternating():
    # ###############################################################################
    # 2 - variable alongshore dune growth parameters
    # ###############################################################################
    # GOAL: what is the effect of the alongshore variability of dunes (15-30 km)?
    #   - vary the growth parameter by varying rmin and rmax, but keep difference (range) constant
    #        - [rmin = 0.35, raverage = 0.6, and rmax = 0.85 everywhere as control case] with diffusive wave parameters
    #        (look at Brie paper to see what conditions are considered diffusive, or high angle)
    #        - THIS RUN: 2 B3Ds at raverage = 0.45 (or 0.3) and 2 B3Ds at raverage=0.75 (or 0.9), all along the barrier, check that
    #        raverage is 0.6 across the barrier; np.mean([0.25, 0.65]) = 0.45 and np.mean([0.55, 0.95]) = 0.75
    #   - show an outcome, but any future conclusions need to come from averaging 20+ different storm scenarios
    #   - hypothesis is that it will prevent punctuated retreat

    # --------- INITIAL CONDITIONS ---------
    name = '2-AlongshoreVarGrowthParam_pt3HAF'
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.3  # fraction of waves approaching from higher than 45 degrees
    slr = 0.002  # m/yr
    ny = 32  # number of alongshore sections (30=15 km, 60=30 km, 32=16 km)
    nt = 1000  # timesteps for 1000 morphologic years
    rmin = [0.25, 0.25, 0.55, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmin = rmin * int(ny/len(rmin))
    rmax = [0.65, 0.65, 0.95, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    rmax = rmax * int(ny/len(rmax))

    # --------- INITIALIZE ---------
    # #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)

    # --------- LOOP ---------
    brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)

    # --------- SAVE ---------
    # #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
    b3d = CASCADE.save(brie, barrier3d, save_directory, name) # this returns the barrier3d model without the BMI

def RUN_3_AlongshoreVarGrowthParam_Gradient():
    # ###############################################################################
    # 3 - variable alongshore dune growth parameters (gradient)
    # ###############################################################################
    # GOAL: what is the effect of the alongshore variability of dunes?
    #        - THIS RUN: make gradient in raverage across the barrier and reduce the grid size to 6 km
    #        - Increased SLR to 0.004

    # --------- INITIAL CONDITIONS ---------
    name = '3-VarGrowthParam_grad_pt2HAF_pt4SLR_1500yrs'
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.004  # m/yr
    ny = 12  # number of alongshore sections (12 = 6 km)
    nt = 1500  # timesteps for 1000 morphologic years
    rmin = [0.25, 0.25, 0.25, 0.35, 0.35, 0.35, 0.45, 0.45, 0.45, 0.55, 0.55, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmax = [0.65, 0.65, 0.65, 0.75, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95, 0.95, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
    # rave = [0.45, 0.45, 0.45, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.75, 0.75, 0.75]  # to help me remember the average

    # --------- INITIALIZE ---------
    # datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/" # iMAC
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
    brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)

    # --------- LOOP ---------
    brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)

    # --------- SAVE ---------
    #save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE/Run_Output"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    CASCADE.save(brie, barrier3d, save_directory, name) # this returns the barrier3d model without the BMI

def PLOT_3_AlongshoreVarGrowthParam_pt2HAF_gradient(name, save_directory):
    # --------- plot ---------
    filename = name + '.npz'
    os.chdir('/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output')
    output = np.load(filename, allow_pickle=True)
    b3d = output['barrier3d']
    brie = output['brie']
    brie = brie[0]

    # 1: Animation Frames of Barrier and Dune Elevation
    TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
    ny = len(b3d)
    CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX, name)

    # ===================================================
    # 7: Calculate shoreline change periodicity
    Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(b3d[0]._x_s_TS)

    # 5: Statistics from B3d
    TMAX = b3d[0].time_index - 1
    iB3D = 0
    CASCADEplt.plot_statistics(b3d, iB3D, TMAX)

    # 2: Shoreline positions over time
    iB3D = 0
    CASCADEplt.plot_ShorelinePositions(b3d[iB3D]._x_s_TS[0:TMAX], b3d[iB3D]._x_b_TS[0:TMAX])

    # 2: Shoreline positions over time
    iB3D = 9
    CASCADEplt.plot_ShorelinePositions(b3d[iB3D]._x_s_TS[0:TMAX], b3d[iB3D]._x_b_TS[0:TMAX])

    # 2: Shoreline change rate (AGU version)
    # CASCADEplt.plot_ShorelineChangeRate_AGU(b3d1, b3d2)

# # ###############################################################################
# # 4 - variable alongshore dune growth parameters (gradient, coalesce in middle high)
# # ###############################################################################
# # GOAL: what is the effect of the alongshore variability of dunes?
# #        - THIS RUN: make gradient in raverage such that they coalesce in the middle
#
# # --------- INITIAL CONDITIONS ---------
# name = '4-AlongshoreVarGrowthParam_pt2HAF_LHL_1500yrs'
# wave_height = 1.0  # m
# wave_period = 7  # s (lowered from 10 s to reduce k_sf)
# asym_frac = 0.8  # fraction approaching from left
# high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
# slr = 0.002  # m/yr
# ny = 12  # number of alongshore sections (12 = 6 km)
# nt = 1500  # timesteps for 1000 morphologic years
# rmin = [0.25, 0.35, 0.35, 0.45, 0.45, 0.55, 0.55, 0.45, 0.45, 0.35, 0.35, 0.25]  # minimum growth rate for logistic dune growth (list for alongshore variability)
# rmax = [0.65, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95, 0.85, 0.85, 0.75, 0.75, 0.65]  # maximum growth rate for logistic dune growth (list for alongshore variability)
#
# # THE FOLLOWING RUNS ARE ON SUPERCOMPUTER
#
# # ###############################################################################
# # 5 - variable alongshore dune growth parameters (gradient, coalesce in middle low)
# # ###############################################################################
# # GOAL: what is the effect of the alongshore variability of dunes?
# #        - THIS RUN: make gradient in raverage such that they coalesce in the middle
#
# # --------- INITIAL CONDITIONS ---------
# name = '5-AlongshoreVarGrowthParam_pt2HAF_HLH_1500yrs'
# wave_height = 1.0  # m
# wave_period = 7  # s (lowered from 10 s to reduce k_sf)
# asym_frac = 0.8  # fraction approaching from left
# high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
# slr = 0.002  # m/yr
# ny = 12  # number of alongshore sections (12 = 6 km)
# nt = 1500  # timesteps for 1000 morphologic years
# rmin = [0.55, 0.45, 0.45, 0.35, 0.35, 0.25, 0.25, 0.35, 0.35, 0.45, 0.45, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
# rmax = [0.95, 0.85, 0.85, 0.75, 0.75, 0.65, 0.65, 0.75, 0.75, 0.85, 0.85, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
#
#
# # ###############################################################################
# # 6 - variable alongshore dune growth parameters (gradient, opposite of #3)
# # ###############################################################################
# # GOAL: does directionality matter for the "dunes that are most important may not be by you" finding? trying to
# # understand the alongshore transport better in BRIE
# #        - THIS RUN: reverse the growth conditions of #3 (i.e., high to low)
#
# # --------- INITIAL CONDITIONS ---------
# name = '3-AlongshoreVarGrowthParam_pt2HAF_gradient_1500yrs'
# wave_height = 1.0  # m
# wave_period = 7  # s (lowered from 10 s to reduce k_sf)
# asym_frac = 0.8  # fraction approaching from left
# high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
# slr = 0.002  # m/yr
# ny = 12  # number of alongshore sections (12 = 6 km)
# nt = 1500  # timesteps for 1000 morphologic years
# rmin = [0.25, 0.25, 0.25, 0.35, 0.35, 0.35, 0.45, 0.45, 0.45, 0.55, 0.55, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
# rmax = [0.65, 0.65, 0.65, 0.75, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95, 0.95, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
# # rave = [0.45, 0.45, 0.45, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.75, 0.75, 0.75]  # to help me remember the average


def RUN_4_B3D_Rave_pt45_SLR_pt004():
    # ###############################################################################
    # 4 - check B3D low-end dune growth parameters
    # ###############################################################################
    # GOAL: check that the low-growth rate parameter in B3D doesn't undergo punctuated retreat (increased SLR to 0.004)

    # --------- INITIAL CONDITIONS ---------
    name = '4-B3D_Rave_pt45_SLR_pt004'
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.004  # m/yr
    ny = 12  # number of alongshore sections (12 = 6 km)
    nt = 1000  # timesteps for 1000 morphologic years
    rmin = 0.25  # 0.2  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmax = 0.65  # 0.4  # maximum growth rate for logistic dune growth (list for alongshore variability)
    # rave = [0.45]  # to help me remember the average

    # --------- INITIALIZE ---------
    #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin,
                                         rmax, datadir)

    # --------- LOOP ---------
    # just use the first B3D grid and update with B3D
    Time = time.time()

    for time_step in range(brie._nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", 'Time Step: ', time_step, end="")
        barrier3d[0].update()

    SimDuration = time.time() - Time
    print()
    print('Elapsed Time: ', SimDuration, 'sec')  # Print elapsed time of simulation

    # --------- SAVE ---------
    #save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE/Run_Output"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    os.chdir(save_directory)
    b3d = barrier3d[0]._model  # just save the model, not the BMI
    filename = name + '.npz'
    np.savez(filename, barrier3d=b3d)

    # ===================================================
    # 7: Calculate shoreline change periodicity
    Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(b3d._x_s_TS)

    # HOORAY, NO PUNCTUATED RETREAT

    # 2: Shoreline positions over time
    TMAX = b3d.time_index - 1
    CASCADEplt.plot_ShorelinePositions(b3d._x_s_TS[0:TMAX], b3d._x_b_TS[0:TMAX])

def RUN_5_B3D_Rave_pt75_SLR_pt004():
    # ###############################################################################
    # 5 - check B3D high-end dune growth parameters
    # ###############################################################################
    # GOAL: check that the high-growth rate parameter in B3D DOES undergo punctuated retreat (increased SLR to 0.004)

    # --------- INITIAL CONDITIONS ---------
    name = '5-B3D_Rave_pt75_SLR_pt004'
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.004  # m/yr
    ny = 12  # number of alongshore sections (12 = 6 km)
    nt = 1000  # timesteps for 1000 morphologic years
    rmin = 0.55  # minimum growth rate for logistic dune growth (list for alongshore variability)
    rmax = 0.95  # maximum growth rate for logistic dune growth (list for alongshore variability)
    # rave = [0.75]  # to help me remember the average

    # --------- INITIALIZE ---------
    #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/"
    datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
    brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin,
                                         rmax, datadir)

    # --------- LOOP ---------
    # just use the first B3D grid and update with B3D
    Time = time.time()

    for time_step in range(brie._nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", 'Time Step: ', time_step, end="")
        barrier3d[0].update()

    SimDuration = time.time() - Time
    print()
    print('Elapsed Time: ', SimDuration, 'sec')  # Print elapsed time of simulation

    # --------- SAVE ---------
    #save_directory = "/Users/katherineanarde/PycharmProjects/CASCADE/Run_Output"
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
    os.chdir(save_directory)
    b3d = barrier3d[0]._model  # just save the model, not the BMI
    filename = name + '.npz'
    np.savez(filename, barrier3d=b3d)

    # ===================================================
    # 7: Calculate shoreline change periodicity
    Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(b3d._x_s_TS)

    # 2: Shoreline positions over time
    TMAX = b3d.time_index - 1
    CASCADEplt.plot_ShorelinePositions(b3d._x_s_TS[0:TMAX], b3d._x_b_TS[0:TMAX])

