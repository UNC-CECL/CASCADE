import multiprocessing

import CASCADE as CASCADE

# for laptop and desktop, use all but one core; on supercomputer, use all cores
num_cores = multiprocessing.cpu_count() - 1

# --------- INITIAL CONDITIONS ---------
name = "test_drowning"
wave_height = 1.0  # m
wave_period = 7  # s (lowered from 10 s to reduce k_sf)
asym_frac = 0.8  # fraction approaching from left
high_ang_frac = 0.2  # fraction of waves approaching from higher than 45 degrees
slr = 0.1  # m/yr
ny = 6  # number of alongshore sections (12 = 6 km)
nt = 1000  # timesteps for 1000 morphologic years
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
datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"  # laptop
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
brie, barrier3d = CASCADE.time_loop(brie, barrier3d, num_cores)
