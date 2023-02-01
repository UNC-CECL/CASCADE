"""

test that new version of BMI against older versions
v1 - AGU version
v2 - Jan 26, updates to CASCADE, B3D, and Brie

"""

import numpy as np
import os
import multiprocessing

import matplotlib.pyplot as plt

from CASCADE import Cascade

# NOTE: I saved the 100 year (15 core) simulation for #1 and 500 year simulation for #3 for debugging comparison for V2
v1_500yr = "3-AlongshoreVarGrowthParam_pt2HAF_gradient_500yrs_15cores_v1.npz"
v2_500yr = "3-AlongshoreVarGrowthParam_pt2HAF_gradient_500yrs_15cores_v2.npz"
v2_100yr = "1-CASCADE_LTA_COMPARISON_3km_100yr_v2.npz"

new_100_version = "1-CASCADE_LTA_COMPARISON_3km_100yr_v3"
new_500_version = "3-AlongshoreVarGrowthParam_pt2HAF_gradient_500yrs_15cores_v3"  # 100 is a sufficient test


def run_1_100years(name=new_100_version):

    # --------- INITIAL CONDITIONS ---------
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.002  # m/yr
    ny = 6  # 12 # number of alongshore sections (6=3 km for 3000 yr run, 12=6 km for 1500 yr run)
    nt = 100  # 3000  #1500 # timesteps for 3000 morphologic years
    rmin = 0.35  # minimum growth rate for logistic dune growth (can be a list)
    rmax = 0.85  # maximum growth rate for logistic dune growth (can be a list)

    # all but 1 core
    num_cores = multiprocessing.cpu_count() - 1
    datadir = "cascade/data/pathways_data/barrier3d-default-parameters.yaml"  # laptop
    brie, barrier3d_15cores = CASCADE.initialize(
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
        datadir,
    )
    # --------- LOOP ---------
    brie, barrier3d_15cores = CASCADE.time_loop(brie, barrier3d_15cores, num_cores)
    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output/BMI_Version_Tests"
    CASCADE.save(brie, barrier3d_15cores, save_directory, name)


def run_3_500years(name=new_500_version):

    # --------- INITIAL CONDITIONS ---------
    wave_height = 1.0  # m
    wave_period = 7  # s (lowered from 10 s to reduce k_sf)
    asym_frac = 0.8  # fraction approaching from left
    high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
    slr = 0.002  # m/yr
    ny = 12  # number of alongshore sections (12 = 6 km)
    nt = 500  # timesteps
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

    # all but 1 core
    num_cores = multiprocessing.cpu_count() - 1
    datadir = "cascade/data/pathways_data/"
    brie, barrier3d_3_15cores = CASCADE.initialize(
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
        datadir,
    )
    # --------- LOOP ---------
    brie, barrier3d_3_15cores = CASCADE.time_loop(brie, barrier3d_3_15cores, num_cores)

    save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output/BMI_Version_Tests"
    CASCADE.save(brie, barrier3d_3_15cores, save_directory, name)


def plot_500_diff(filename_old=v1_500yr, filename_new=new_500_version):

    # load new version
    os.chdir("/Run_Output/BMI_Version_Tests")
    output = np.load(filename_new, allow_pickle=True)
    b3d_new = output["barrier3d"]

    # load old 500 year version
    output = np.load(filename_old, allow_pickle=True)
    b3d_original = output["barrier3d"]

    # check that the outputs are all the same
    i = 11
    plt.figure()
    plt.plot(b3d_original[i].x_s_TS, "b")
    plt.plot(b3d_new[i].x_s_TS, "g")

    plt.figure()
    plt.plot(b3d_original[i].x_b_TS, "b")
    plt.plot(b3d_new[i].x_b_TS, "g")

    # for v1 vs v2, the outputs diverge different late in the simulation:
    # 1) there is slightly less overwash in this new version
    # 2) there is slightly more shoreline retreat
    # See notes: I'm not too concerned - these differences probably come from my changes (when debugging) the shoreline code
    plt.figure()
    plt.plot(b3d_original[i].QowTS, "b")
    plt.plot(b3d_new[i].QowTS, "g")

    plt.figure()
    plt.plot(b3d_original[i].QsfTS, "b")
    plt.plot(b3d_new[i].QsfTS, "g")


def plot_100_diff(filename_old=v2_100yr, filename_new=new_100_version):

    # load new version
    os.chdir("/Run_Output/BMI_Version_Tests")
    output = np.load(filename_new, allow_pickle=True)
    b3d_new = output["barrier3d"]

    # load old 500 year version
    output = np.load(filename_old, allow_pickle=True)
    b3d_original = output["barrier3d"]

    # check that the outputs are all the same
    i = 0
    plt.figure()
    plt.plot(b3d_original[i].x_s_TS, "b")
    plt.plot(b3d_new[i].x_s_TS, "g")

    plt.figure()
    plt.plot(b3d_original[i].x_b_TS, "b")
    plt.plot(b3d_new[i].x_b_TS, "g")

    # for v1 vs v2, the outputs diverge different late in the simulation:
    # 1) there is slightly less overwash in this new version
    # 2) there is slightly more shoreline retreat
    # See notes: I'm not too concerned - these differences probably come from my changes (when debugging) the shoreline code
    plt.figure()
    plt.plot(b3d_original[i].QowTS, "b")
    plt.plot(b3d_new[i].QowTS, "g")

    plt.figure()
    plt.plot(b3d_original[i].QsfTS, "b")
    plt.plot(b3d_new[i].QsfTS, "g")
