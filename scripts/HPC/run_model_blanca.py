# Model runfile for CASCADE simulations
# Written by K.Anarde

# remember if I move to a different computer to $ pip install -e . in the brie and B3D directories for the BMI

# import os

# import CASCADE_plotters as CASCADEplt

import CASCADE as CASCADE

# ###############################################################################
# 4 - variable alongshore dune growth parameters (gradient, coalesce in middle)
# ###############################################################################
# GOAL: what is the effect of the alongshore variability of dunes?
#        - THIS RUN: make gradient in raverage such that they coalesce in the middle

# --------- INITIAL CONDITIONS ---------
name = "4-AlongshoreVarGrowthParam_pt2HAF_LHL_1500yrs"
wave_height = 1  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 12  # number of alongshore sections (12 = 6 km)
nt = 1500  # timesteps for 1000 morphologic years
rmin = [
    0.25,
    0.35,
    0.35,
    0.45,
    0.45,
    0.55,
    0.55,
    0.45,
    0.45,
    0.35,
    0.35,
    0.25,
]  # minimum growth rate for logistic dune growth (list for alongshore variability)
rmax = [
    0.65,
    0.75,
    0.75,
    0.85,
    0.85,
    0.95,
    0.95,
    0.85,
    0.85,
    0.75,
    0.75,
    0.65,
]  # maximum growth rate for logistic dune growth (list for alongshore variability)

# --------- INITIALIZE ---------
datadir = "barrier3d-parameters.yaml"
# datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/barrier3d-parameters.yaml"
brie, barrier3d = CASCADE.initialize(
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
brie, barrier3d = CASCADE.time_loop(brie, barrier3d)

# --------- SAVE ---------
# #datadir = "/Users/katherineanarde/PycharmProjects/CASCADE/"
b3d = CASCADE.save(
    brie, barrier3d, name
)  # this returns the barrier3d model without the BMI
