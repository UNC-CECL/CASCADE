import multiprocessing

import matplotlib.pyplot as plt

import CASCADE as CASCADE

# use the default run to test the parallel code for different core options for a 100 year simulation and
# 6 alongshore sections

# ###############################################################################
# 1 - Original CASCADE_LTA_COMPARISON
# ###############################################################################
# GOAL: highlight different processes in models with alongshore homogenous dune line, 3000 year simulation
#
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

# 6 cores
# --------- INITIALIZE ---------
name = "1-CASCADE_Parallel_6cores"
num_cores = 6
datadir = "cascade/data/pathways_data/barrier3d-parameters.yaml"  # laptop
brie, barrier3d_6cores = CASCADE.initialize(
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
brie, barrier3d_6cores = CASCADE.time_loop(brie, barrier3d_6cores, num_cores)

# all but 1 core
name = "1-CASCADE_LTA_COMPARISON_3km_100yr_Parallel_15cores"
num_cores = multiprocessing.cpu_count() - 1
datadir = "cascade/data/pathways_data/barrier3d-parameters.yaml"  # laptop
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
save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
b3d_1_15cores = CASCADE.save(
    brie, barrier3d_15cores, save_directory, name
)  # this returns the barrier3d model without the BMI

# only 1 core (not in parallel)
name = "1-CASCADE_NoParallel"
num_cores = 1
datadir = "cascade/data/pathways_data/barrier3d-parameters.yaml"  # laptop
brie, barrier3d_1cores = CASCADE.initialize(
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
brie, barrier3d_1cores = CASCADE.time_loop(brie, barrier3d_1cores, num_cores)

# make sure the outputs are the same...they are!
plt.figure()
plt.plot(barrier3d_1cores[0]._model.x_s_TS, "b")
plt.plot(barrier3d_6cores[0]._model.x_s_TS, "g")
plt.plot(barrier3d_15cores[0]._model.x_s_TS, "r")

plt.figure()
plt.plot(barrier3d_1cores[0]._model.x_b_TS, "b")
plt.plot(barrier3d_6cores[0]._model.x_b_TS, "g")
plt.plot(barrier3d_15cores[0]._model.x_b_TS, "r")
