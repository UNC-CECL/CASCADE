# Model runfile for CASCADE simulations
# Written by K.Anarde

import numpy as np
import os

import CASCADE_plotters as CASCADEplt

import CASCADE as CASCADE

# ###############################################################################
# 1 - CASCADE_LTA_COMPARISON
# ###############################################################################
# GOAL: highlight different processes in models with alongshore homogenous dune line, 3000 year simulation
#
# --------- INITIAL CONDITIONS ---------
name = '1-CASCADE_LTA_COMPARISON_3000yr'
wave_height = 1  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 6  # number of alongshore sections (6=3 km for testing AST, make 30 for inlets=15 km)
nt = 3000  # timesteps for 3000 morphologic years
rmin = 0.35 # minimum growth rate for logistic dune growth (can be a list)
rmax = 0.85 # maximum growth rate for logistic dune growth (can be a list)

# --------- INITIALIZE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)

# --------- LOOP ---------
brie, barrier3d = CASCADE.time_loop(brie, barrier3d)

# --------- RUN LTA ---------
w_b_crit = 450  # critical barrier width [m]
h_b_crit = 1.9  # (should equal B3D original BermEl in the yaml file, not what is presented in B3D (minus the MHW)
Qow_max = 20  # max overwash flux [m3/m/yr]
brieLTA = CASCADE.LTA(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, w_b_crit, h_b_crit, Qow_max)

# --------- SAVE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
b3d = CASCADE.save(brie, barrier3d, save_directory, name) # this returns the barrier3d model without the BMI

# --------- plot ---------
# load the 3,000 year run from the other computer
filename = name + '.npz'
output = np.load(filename, allow_pickle=True)
b3d = output['barrier3d']
brie = output['brie']
brie = brie[0]

# 1: Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for CASCADE)
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
ny = len(b3d)
CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX, name)

#===================================================

# 2: Shoreline positions over time (#6 in Barrier3D_Functions)
CASCADEplt.plot_ShorelinePositions(b3d[0]._x_s_TS, b3d[0]._x_b_TS)

#===================================================

# 3: Cross-Shore Transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
CASCADEplt.plot_XShoreTransects(b3d[0], TMAX)

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

# 6: Statistics from BrieLTA [LEFT OFF HERE]
TMAX = int((b3d[0].time_index - 1) / brieLTA._dt)
iB3D = 1
iB3D = iB3D * int( (b3d[0]._BarrierLength * 10) / brieLTA._dy)
CASCADEplt.plot_statistics_BRIE(brieLTA, iB3D, TMAX)

# ===================================================
# 7: Calculate shoreline change periodicity
Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(b3d[0]._x_s_TS)

# ###############################################################################
# 2 - variable alongshore dune growth parameters
# ###############################################################################
# GOAL: what is the effect of the alongshore variability of dunes (15-30 km)?
#   - vary the growth parameter by varying rmin and rmax, but keep difference (range) constant
#        - [rmin = 0.35, raverage = 0.6, and rmax = 0.85 everywhere as control case] with diffusive wave parameters
#        (look at Brie paper to see what conditions are considered diffusive, or high angle)
#        - 2 B3Ds at raverage = 0.45 (or 0.3) and 2 B3Ds at raverage=0.75 (or 0.9), all along the barrier, check that
#        raverage is 0.6 across the barrier; np.mean([0.25, 0.65]) = 0.45 and np.mean([0.55, 0.95]) = 0.75
#   - show A outcome, but any conclusion need to come from averaging of different storm scenarios
#   - hypothesis is that it will prevent punctuated retreat

# --------- INITIAL CONDITIONS ---------
name = '2-AlongshoreVarGrowthParam_pt3HAF'
wave_height = 1  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.3  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 32  # number of alongshore sections (30=15 km, 60=30 km, 32 = 16 km)
nt = 1000  # timesteps for 1000 morphologic years
rmin = [0.25, 0.25, 0.55, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
rmin = rmin * int(ny/len(rmin))
rmax = [0.65, 0.65, 0.95, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
rmax = rmax * int(ny/len(rmax))

# --------- INITIALIZE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)

# --------- LOOP ---------
brie, barrier3d = CASCADE.time_loop(brie, barrier3d)

# --------- SAVE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
b3d = CASCADE.save(brie, barrier3d, save_directory, name) # this returns the barrier3d model without the BMI

# --------- plot ---------
filename = name + '.npz'
os.chdir('/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/')
output = np.load(filename, allow_pickle=True)
b3d = output['barrier3d']
brie = output['brie']
brie = brie[0]

# 1: Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for CASCADE)
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
ny = len(b3d)
CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX, name)

#===================================================

# 2: Shoreline positions over time (#6 in Barrier3D_Functions)
CASCADEplt.plot_ShorelinePositions(b3d[0]._x_s_TS, b3d[0]._x_b_TS)

#===================================================

# 3: Cross-Shore Transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
CASCADEplt.plot_XShoreTransects(b3d[0], TMAX)

# ===================================================
# 7: Calculate shoreline change periodicity
Periodicity, AvgFastDur, AvgSlowDur, Punc = CASCADEplt.calc_ShorelinePeriodicity(b3d[0]._x_s_TS)

# ###############################################################################
# 2 - variable alongshore dune growth parameters
# ###############################################################################
# GOAL: what is the effect of the alongshore variability of dunes (15-30 km)?
#   - vary the growth parameter by varying rmin and rmax, but keep difference (range) constant
#        - [rmin = 0.35, raverage = 0.6, and rmax = 0.85 everywhere as control case] with diffusive wave parameters
#        (look at Brie paper to see what conditions are considered diffusive, or high angle)
#        - 2 B3Ds at raverage = 0.45 (or 0.3) and 2 B3Ds at raverage=0.75 (or 0.9), all along the barrier, check that
#        raverage is 0.6 across the barrier; np.mean([0.25, 0.65]) = 0.45 and np.mean([0.55, 0.95]) = 0.75
#   - show A outcome, but any conclusion need to come from averaging of different storm scenarios
#   - hypothesis is that it will prevent punctuated retreat

# --------- INITIAL CONDITIONS ---------
name = '3-AlongshoreVarGrowthParam_p23HAF_gradient'
wave_height = 1  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.7  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 12  # number of alongshore sections (30=15 km, 60=30 km, 32 = 16 km)
nt = 500  # timesteps for 1000 morphologic years
rmin = [0.25, 0.25, 0.25, 0.35, 0.35, 0.35, 0.45, 0.45, 0.45, 0.55, 0.55, 0.55]  # minimum growth rate for logistic dune growth (list for alongshore variability)
#rmin = rmin * int(ny/len(rmin))
rmax = [0.65, 0.65, 0.65, 0.75, 0.75, 0.75, 0.85, 0.85, 0.95, 0.95, 0.95, 0.95]  # maximum growth rate for logistic dune growth (list for alongshore variability)
#rmax = rmax * int(ny/len(rmax))

# --------- INITIALIZE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
brie, barrier3d = CASCADE.initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir)

# --------- LOOP ---------
brie, barrier3d = CASCADE.time_loop(brie, barrier3d)

# --------- SAVE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
save_directory = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/"
b3d = CASCADE.save(brie, barrier3d, save_directory, name) # this returns the barrier3d model without the BMI
