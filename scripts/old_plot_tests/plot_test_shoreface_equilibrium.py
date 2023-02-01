import multiprocessing

import matplotlib.pyplot as plt

import CASCADE as CASCADE
from barrier3d import Barrier3dBmi
import time

# Someone asked a question during AGU about whether or not I start the shoreface out of equilibrium.
# When I went back and checked, B3D (both V1 and the BMI) start in equilibrium, but CASCADE starts
# slightly out of equilibrium. Check that this is a precision issue.

num_cores = multiprocessing.cpu_count() - 1

# ###############################################################################
# 1 - CASCADE_LTA_COMPARISON
# ###############################################################################
# GOAL: highlight different processes in models with alongshore homogenous dune line, 3000 year simulation
#
# --------- INITIAL CONDITIONS ---------
name = "1-CASCADE_LTA_COMPARISON_3km_3000yr"
# name = '1-CASCADE_LTA_COMPARISON_6km_1500yr'
wave_height = 1.0  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
slr = 0.002  # m/yr
ny = 6  # 12 # number of alongshore sections (6=3 km for 3000 yr run, 12=6 km for 1500 yr run)
nt = 100  # 3000  #1500 # timesteps for 3000 morphologic years
rmin = 0.35  # minimum growth rate for logistic dune growth (can be a list)
rmax = 0.85  # maximum growth rate for logistic dune growth (can be a list)

# --------- INITIALIZE ---------
datadir = "cascade/data/pathways_data/barrier3d-default-parameters.yaml"  # laptop
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

# initial brie conditions
print(brie._s_sf_eq)
print(
    brie._s_sf_save[0, 0]
)  # updated on 1/14 to represent the actual starting s_sf and not s_sf_eq as Jaap wrote it

# initial barrier3d conditions
print(barrier3d[0]._model._s_sf_eq)
print(
    barrier3d[0]._model._s_sf_TS
)  # ok, so this is slightly below the equilibrium slope...why?

# ok, I understand what is happening: the equilibrium slope in brie is calculated using an equation that takes into
# account the settling velocity and wave period; in Barrier3D, s_sf_eq is an input and Ian just sets it to the initial
# s_sf calculated for the other input parameters BUT you could start the model out of equilibrium by starting with
# an initial shoreface length that is larger or smaller than that needed to achieve equilibrium. Regardless, brie also
# starts slightly out of equilibrium for these initial conditions, so we are good!

# one last test on shoreline change
brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)

# now run the same input barrier3d file without the AST model
# create an instance of the new BMI class, which is the model
barrier3d_noAST = Barrier3dBmi()
barrier3d_noAST.initialize(datadir)

# increase time step
Time = time.time()
for time_step in range(1, barrier3d_noAST._model._TMAX):
    barrier3d_noAST.update()  # update the model by a time step

    # Print time step to screen
    print("\r", "Time Step: ", time_step, end="")

SimDuration = time.time() - Time
print()
print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

# compare the shoreline change for the first B3D domain from the B3D+BRIE coupling and the B3D only run
plt.figure()
plt.plot(barrier3d_noAST._model.x_s_TS, "b")
plt.plot(barrier3d[0]._model.x_s_TS, "g")
plt.xlabel("Time")
plt.ylabel("Shoreline Position")
