import numpy as np
import os
import multiprocessing

import matplotlib.pyplot as plt

import CASCADE as CASCADE
import CASCADE_plotters as CASCADEplt

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
ny = 6  #12 # number of alongshore sections (6=3 km for 3000 yr run, 12=6 km for 1500 yr run)
nt = 100  # 3000  #1500 # timesteps for 3000 morphologic years
rmin = 0.35  # minimum growth rate for logistic dune growth (can be a list)
rmax = 0.85  # maximum growth rate for logistic dune growth (can be a list)

# 6 cores
# --------- INITIALIZE ---------
name = '1-CASCADE_Parallel_6cores'
num_cores = 6
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"  # laptop
brie, barrier3d_6cores = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)
# --------- LOOP ---------
brie, barrier3d_6cores = CASCADE.time_loop(brie, barrier3d_6cores, num_cores)

# all but 1 core
name = '1-CASCADE_LTA_COMPARISON_3km_100yr_Parallel_15cores'
num_cores = multiprocessing.cpu_count() - 1
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"  # laptop
brie, barrier3d_15cores = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)
# --------- LOOP ---------
brie, barrier3d_15cores = CASCADE.time_loop(brie, barrier3d_15cores, num_cores)
save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
b3d_1_15cores = CASCADE.save(brie, barrier3d_15cores, save_directory, name) # this returns the barrier3d model without the BMI

# only 1 core (not in parallel)
name = '1-CASCADE_NoParallel'
num_cores = 1
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"  # laptop
brie, barrier3d_1cores = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)
# --------- LOOP ---------
brie, barrier3d_1cores = CASCADE.time_loop(brie, barrier3d_1cores, num_cores)

# make sure the outputs are the same...they are!
plt.figure()
plt.plot(barrier3d_1cores[0]._model.x_s_TS, 'b')
plt.plot(barrier3d_6cores[0]._model.x_s_TS, 'g')
plt.plot(barrier3d_15cores[0]._model.x_s_TS, 'r')

plt.figure()
plt.plot(barrier3d_1cores[0]._model.x_b_TS, 'b')
plt.plot(barrier3d_6cores[0]._model.x_b_TS, 'g')
plt.plot(barrier3d_15cores[0]._model.x_b_TS, 'r')

# --------- Now, make sure the outputs agree with earlier versions of the code ---------
# I only have one 500 year run saved from the earlier version of the code, and it is #3
# Note, a 6 km run for 500 years took 2 hours with the new parallelization

# ###############################################################################
# 3 - variable alongshore dune growth parameters (gradient)
# ###############################################################################
# GOAL: what is the effect of the alongshore variability of dunes?
#        - THIS RUN: make gradient in raverage across the barrier and reduce the grid size to 6 km

# --------- INITIAL CONDITIONS ---------
name = '3-AlongshoreVarGrowthParam_pt2HAF_gradient_500yrs_15cores'
wave_height = 1.0  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 12  # number of alongshore sections (12 = 6 km)
nt = 500  # timesteps
rmin = [0.25, 0.25, 0.25, 0.35, 0.35, 0.35, 0.45, 0.45, 0.45, 0.55, 0.55, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
rmax = [0.65, 0.65, 0.65, 0.75, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95, 0.95, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
# rave = [0.45, 0.45, 0.45, 0.55, 0.55, 0.55, 0.65, 0.65, 0.65, 0.75, 0.75, 0.75]  # to help me remember the average

# all but 1 core
num_cores = multiprocessing.cpu_count() - 1
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
brie, barrier3d_3_15cores = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)
# --------- LOOP ---------
brie, barrier3d_3_15cores = CASCADE.time_loop(brie, barrier3d_3_15cores, num_cores)

save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output"
b3d_3_15cores = CASCADE.save(brie, barrier3d_3_15cores, save_directory, name) # this returns the barrier3d model without the BMI

# load earlier version
filename = '3-VarGrowthParam_grad_pt2HAF_pt2SLR_500yrs.npz'
os.chdir('/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/Run_Output')
output = np.load(filename, allow_pickle=True)
b3d_original = output['barrier3d']

# check that the outputs are all the same
i = 11
plt.figure()
plt.plot(b3d_original[i].x_s_TS, 'b')
plt.plot(b3d_3_15cores[i].x_s_TS, 'g')

plt.figure()
plt.plot(b3d_original[i].x_b_TS, 'b')
plt.plot(b3d_3_15cores[i].x_b_TS, 'g')

# the outputs diverge different late in the simulation:
# 1) there is slightly less overwash in this new version
# 2) there is slightly more shoreline retreat
# See notes: I'm not too concerned - these differences probably come from my changes (when debugging) the shoreline code
plt.figure()
plt.plot(b3d_original[i].QowTS, 'b')
plt.plot(b3d_3_15cores[i].QowTS, 'g')

plt.figure()
plt.plot(b3d_original[i].QsfTS, 'b')
plt.plot(b3d_3_15cores[i].QsfTS, 'g')

# 1: Animation Frames of Barrier and Dune Elevation
TMAX = b3d_3_15cores[0].time_index - 1  # just in case the barrier drowned
ny = len(b3d_3_15cores)
name = '3-AlongshoreVarGrowthParam_pt2HAF_gradient_500yrs_15cores'
CASCADEplt.plot_ElevAnimation(b3d_3_15cores, ny, save_directory, TMAX, name)

# NOTE: I saved the 100 year (15 core) simulation for #1 and 500 year simulation for #3 for future debugging
# comparisions