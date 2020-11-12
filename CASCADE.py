# Model guts to couple barrier3d and brie
# Written by K.Anarde

from yaml import full_load, dump
import numpy as np
import time
import os

from barrier3d import Barrier3dBmi
from brie import Brie

# temporary placement - I want this added to the BMI but keep getting an object error
def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)

    doc[var_name] = new_vals

    with open(file_name, 'w') as f:
        dump(doc, f)


###############################################################################
# initial conditions for Brie
###############################################################################

def initialize(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, rmin, rmax, datadir):

    # start by initializing brie because it has physical parameters related to wave climate that we will use to calculate
    # params for B3D
    brie = Brie()  # initialize class

    # update the initial conditions
    brie._name = name
    brie._plot_on = False
    brie._make_gif = False
    brie._ast_model_on = True  # shoreface formulations on
    brie._inlet_model_on = False  # inlet model off
    brie._barrier_model_on = False  # overwash model off
    brie._b3d_barrier_model_on = True  # B3d overwash model on

    # wave climate parameters
    brie._wave_height = wave_height  # m
    brie._wave_period = wave_period  # s
    brie._wave_asym = asym_frac  # fraction approaching from left
    brie._wave_high = high_ang_frac  # fraction of waves approaching from higher than 45 degrees

    # barrier model parameters (the following are needed for other calculations even if the barrier model is off)
    brie._slr = slr  # m/yr
    brie._s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
    brie._z = 10  # initial sea level (for tracking SL, Eulerian reference frame)
    brie._bb_depth = 3  # back-barrier depth

    # model setup
    brie._dy = 500  # m, length of alongshore section (same as B3D)
    brie._ny = ny  # number of alongshore sections (6=3 km for testing AST, make 30 for inlets=15 km)
    brie._dt = 1  # yr, timestep (same as B3D)
    brie._nt = nt  # timesteps for 200 morphologic years
    brie._dtsave = 1  # save spacing (every year)

    # inlet parameters (use default)
    brie._Jmin = 10000  # minimum inlet spacing [m]
    brie._a0 = 0.5  # amplitude of tide [m]
    brie._marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

    # get dependent variables
    Brie.dependent(brie)

    ###############################################################################
    # initial conditions for Barrier3D
    ###############################################################################

    # for each B3D subgrid, set the initial shoreface geometry equal to what is set in brie (some random perturbations);
    # all other B3D variables are set equal
    barrier3d = []

    for iB3D in range(brie._ny):

        barrier3d.append(Barrier3dBmi())  # initialize each class of the BMI

        # update yaml file (these are the only variables that I'm like to change from default)
        set_yaml('TMAX', brie._nt, datadir)  # [yrs] Duration of simulation (if brie._dt = 1 yr, set to ._nt)
        set_yaml('BarrierLength', brie._dy, datadir)  # [m] Static length of island segment (comprised of 10x10 cells)
        set_yaml('DShoreface', brie._d_sf, datadir)  # [m] Depth of shoreface (set to brie depth, function of wave height)
        set_yaml('LShoreface', float(brie._x_s[iB3D] - brie._x_t[iB3D]),
                 datadir)  # [m] Length of shoreface (calculate from brie variables, shoreline - shoreface toe)
        set_yaml('ShorefaceToe', float(brie._x_t[iB3D]), datadir)  # [m] Start location of shoreface toe
        #set_yaml('BermEl', 1.9 , datadir) # [m] Static elevation of berm (NOTE: if this is changed, the MSSM storm list and storm time series needs to be remade)
        set_yaml('BayDepth', brie._bb_depth, datadir)  # [m] Depth of bay behind island segment (set to brie bay depth)
        #set_yaml('MHW', 0.46, datadir)  # [m] Elevation of Mean High Water (NOTE: if this is changed, the storm time series needs to be remade)
        set_yaml('DuneParamStart', True, datadir)  # Dune height will come from external file
        set_yaml('Rat', 0.0,
                 datadir)  # [m / y] corresponds to Qat in LTA formulations (!!! must set to 0 because Qs is calculated in brie !!!)
        set_yaml('RSLR_Constant', True,
                 datadir)  # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
        set_yaml('RSLR_const', brie._slr, datadir)  # [m / y] Relative sea-level rise rate
        #set_yaml('beta', 0.04, datadir)  # Beach slope for runup calculations
        set_yaml('k_sf', float(brie._k_sf),
                 datadir)  # [m^3 / m / y] Shoreface flux rate constant (function of wave parameters from brie)
        set_yaml('s_sf_eq', float(brie._s_sf_eq),
                 datadir)  # Equilibrium shoreface slope (function of wave and sediment parameters from brie)
        if np.size(rmin) > 1:
            set_yaml('GrowthParamStart', False, datadir)  # Dune growth parameter will come from external file
            set_yaml('rmin', rmin[iB3D], datadir)  # Minimum growth rate for logistic dune growth
            set_yaml('rmax', rmin[iB3D], datadir)  # Maximum growth rate for logistic dune growth
        else:
            set_yaml('GrowthParamStart', True, datadir)  # Dune growth parameter will come from external file
            set_yaml('rmin', rmin, datadir)  # Minimum growth rate for logistic dune growth
            set_yaml('rmax', rmax, datadir)  # Maximum growth rate for logistic dune growth

        barrier3d[iB3D].initialize(datadir)

        # debugging: check that the shoreface toe, shoreline, back-barrier, and h_b are correct between the two models
        # brie._x_t[iB3D] - (barrier3d[iB3D]._model._x_t_TS[0] * 10)  # this isn't zero; rounding error
        # brie._x_s[iB3D] - (barrier3d[iB3D]._model._x_s_TS[0] * 10)

        # now update brie x_b, x_b_save[:,0], h_b, h_b_save[:,0] from B3D so all the initial conditions are the same
        brie._x_b[iB3D] = (barrier3d[iB3D]._model._x_b_TS[0] * 10)
        brie._h_b[iB3D] = (barrier3d[iB3D]._model._h_b_TS[0] * 10)
        brie._x_b_save[iB3D, 0] = brie._x_b[iB3D]
        brie._h_b_save[iB3D, 0] = brie._h_b[iB3D]

    return brie, barrier3d

###############################################################################
# time loop
###############################################################################

def time_loop(brie, barrier3d):
    # preallocate arrays for shoreline and barrier height change used in time loops
    x_t_dt, x_s_dt, x_b_dt, h_b_dt = [
        np.zeros(np.size(barrier3d)).astype(float) for _ in range(4)
    ]

    Time = time.time()

    for time_step in range(brie._nt-1):

        # Print time step to screen
        print("\r", 'Time Step: ', time_step, end="")

        for iB3D in range(brie._ny):
            # advance B3D by one time step (this is time_index = 2 at start of loop)
            barrier3d[iB3D].update()

            # get values for passing to brie (all in dam)
            x_t_TS = barrier3d[iB3D]._values["shoreface_toe_position"]()
            x_s_TS = barrier3d[iB3D]._values["shoreline_position"]()
            x_b_TS = barrier3d[iB3D]._values["back_barrier_shoreline_position"]()
            h_b_TS = barrier3d[iB3D]._values["height_of_barrier"]()

            # calculate the diff in shoreface toe, shorelines (dam), height of barrier and convert to m
            x_t_dt[iB3D] = (x_t_TS[-1] - x_t_TS[-2]) * 10
            x_s_dt[iB3D] = (x_s_TS[-1] - x_s_TS[-2]) * 10
            x_b_dt[iB3D] = (x_b_TS[-1] - x_b_TS[-2]) * 10
            h_b_dt[iB3D] = (h_b_TS[-1] - h_b_TS[-2]) * 10

        # pass values from B3D subdomains to brie for use in second timestep
        # (there has to be a better way to do this with the BMI)
        brie._x_t_dt = x_t_dt
        brie._x_s_dt = x_s_dt
        brie._x_b_dt = x_b_dt
        brie._h_b_dt = h_b_dt

        # update brie one time step (this is time index = 2 at start of loop)
        brie.update()

        # loop to pass x_s and x_b (maybe this will be important for inlets with the x_b_fldt) back to B3D from Brie
        # (convert back to dam)
        for iB3D in range(brie._ny):
            barrier3d[iB3D]._model._x_s = brie._x_s[iB3D] / 10
            barrier3d[iB3D]._model._x_s_TS[-1] = brie._x_s[iB3D] / 10
            # barrier3d[iB3D]._model._x_b = brie._x_b[iB3D] / 10   # maybe will need this if tidal inlets on?
            # barrier3d[iB3D]._model._x_b_TS[-1] = brie._x_b[iB3D] / 10

    SimDuration = time.time() - Time
    print()
    print('Elapsed Time: ', SimDuration, 'sec')  # Print elapsed time of simulation

    return brie, barrier3d

###############################################################################
# save data
###############################################################################

def save(brie, barrier3d, directory, name):

    filename = name + '.npz'

    b3d = []

    # save object without BMI
    for iB3D in range(brie._ny):
        b3d.append(barrier3d[iB3D]._model)

    os.chdir(directory)
    np.savez(filename, barrier3d=b3d, brie=brie)

    return b3d

###############################################################################
# run brie with LTA for comparison
###############################################################################

def LTA(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, w_b_crit, h_b_crit, Qow_max):

    brieLTA = Brie()  # initialize class

    # update the initial conditions
    brieLTA._name = name
    brieLTA._plot_on = False
    brieLTA._make_gif = False
    brieLTA._ast_model_on = True  # shoreface formulations on
    brieLTA._inlet_model_on = False  # inlet model off
    brieLTA._barrier_model_on = True  # overwash model on
    brieLTA._b3d_barrier_model_on = False  # B3d overwash model off

    # wave climate parameters
    brieLTA._wave_height = wave_height  # m (lowered from 1 m to reduce k_sf)
    brieLTA._wave_period = wave_period  # s (lowered from 10 s to reduce k_sf)
    brieLTA._wave_asym = asym_frac  # fraction approaching from left
    brieLTA._wave_high = high_ang_frac  # fraction of waves approaching from higher than 45 degrees

    # barrier model parameters (the following are needed for other calculations even if the barrier model is off)
    brieLTA._slr = slr  # m/yr
    brieLTA._s_background = 0.001  # background slope (for shoreface toe position & inlet calculations)
    brieLTA._z = 10  # initial sea level (for tracking SL, Eulerian reference frame)
    brieLTA._bb_depth = 3  # back-barrier depth

    # also need these to be as similar as possible to "storm conditions" for B3D
    # NOTE: the LTA model does not work well when you set the width and height from B3D so we set the critical barrier
    # width to the width of the B3D Interior Domain
    brieLTA._w_b_crit = w_b_crit  # critical barrier width [m]
    brieLTA._h_b_crit = h_b_crit  # (should equal B3D original BermEl above)
    brieLTA._Qow_max = Qow_max  # max overwash flux [m3/m/yr]

    # model setup
    brieLTA._dy = 100  # m, length of alongshore section (same as B3D)
    brieLTA._ny = ny * 10  # number of alongshore sections (10=6 km for testing AST, make 30 for inlets=15 km)
    brieLTA._dt = 0.05  # yr, timestep (same as B3D)
    brieLTA._nt = int(nt/brieLTA._dt)  # timesteps for 200 morphologic years
    brieLTA._dtsave = 20  # save spacing (every year for 0.05 time step)

    # inlet parameters (use default)
    brieLTA._Jmin = 10000  # minimum inlet spacing [m]
    brieLTA._a0 = 0.5  # amplitude of tide [m]
    brieLTA._marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

    # get dependent variables
    Brie.dependent(brieLTA)

    for time_step in range(int(brieLTA._nt) - 1):
        brieLTA.update()

    return brieLTA