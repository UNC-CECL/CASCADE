# Model runfile for CASCADE simulations
# Written by K.Anarde

import numpy as np

import CASCADE_plotters as CASCADEplt

import CASCADE as CASCADE

# ###############################################################################
# 1 - CASCADE_LTA_COMPARISON
# ###############################################################################
#
# --------- INITIAL CONDITIONS ---------
name = 'CASCADE_LTA_COMPARISON'
wave_height = 1  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 6  # number of alongshore sections (6=3 km for testing AST, make 30 for inlets=15 km)
nt = 10  # timesteps for 200 morphologic years
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
filename = 'CASCADE_BRIE_COMPARISON.npz'
output = np.load(filename, allow_pickle=True)
b3d = output['barrier3d']

# 1: Animation Frames of Barrier and Dune Elevation (#4 in Barrier3D_Functions, modified here for CASCADE)
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
ny = len(b3d)
CASCADEplt.plot_ElevAnimation(b3d, ny, save_directory, TMAX)

#===================================================

# 2: Shoreline positions over time (#6 in Barrier3D_Functions)
CASCADEplt.plot_ShorelinePositions(b3d[0]._x_s_TS, b3d[0]._x_b_TS)

#===================================================

# 3: Cross-Shore Transect for one subgrid every 100 m for last time step (#5 in Barrier3D_Functions)
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
CASCADEplt.plot_XShoreTransects(b3d[0], TMAX)

#===================================================

# 4: Cross-shore transects for both brie and B3d
iB3D = 0
TMAX = b3d[0].time_index - 1  # just in case the barrier drowned
time_step = [0, int(TMAX/2), TMAX-2]
CASCADEplt.plot_ModelTransects(b3d, brieLTA, time_step, iB3D)

#===================================================

