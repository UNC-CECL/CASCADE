# Model coupling of Barrier3D and BRIE

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
----------------------------------------------------"""

from yaml import full_load, dump
from joblib import Parallel, delayed
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

    """initialize both BRIE and Barrier3D"""

    # start by initializing BRIE b/c it has parameters related to wave climate that we use to initialize B3D

    # parameters that we need to initialize in BRIE for coupling (not necessarily default values), but I won't be
    # modifying often (or ever) for CASCADE
    ast_model = True  # shoreface formulations on
    barrier_model = False  # LTA14 overwash model off
    inlet_model = False  # inlet model off
    b3d = True  # B3d overwash model on

    # barrier model parameters (the following are needed for other calculations even if the barrier model is off)
    s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
    z = 10.0  # initial sea level (for tracking SL, Eulerian reference frame)
    bb_depth = 3.0  # back-barrier depth
    h_b_crit = 1.9  # critical barrier height for overwash, used also to calculate shoreline diffusivity;
    # we set equal to the static elevation of berm (NOTE: if the berm elevation is changed, the MSSM storm list and
    # storm time series needs to be remade)

    # inlet parameters (use default; these are here to remind me later that they are important and I can change)
    Jmin = 10000  # minimum inlet spacing [m]
    a0 = 0.5  # amplitude of tide [m]
    marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

    # model setup
    dy = 500  # m, length of alongshore section (same as B3D)
    dt = 1  # yr, timestep (same as B3D)
    dtsave = 1  # save spacing (every year)

    brie = Brie(
        name=name,
        ast_model=ast_model,
        barrier_model=barrier_model,
        inlet_model=inlet_model,
        b3d=b3d,
        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=asym_frac,
        wave_angle_high_fraction=high_ang_frac,
        sea_level_rise_rate=slr,
        sea_level_initial=z,
        barrier_height_critical=h_b_crit,
        tide_amplitude=a0,
        back_barrier_marsh_fraction=marsh_cover,
        back_barrier_depth=bb_depth,
        xshore_slope=s_background,
        inlet_min_spacing=Jmin,
        alongshore_section_length=dy,
        alongshore_section_count=ny,
        time_step=dt,
        time_step_count=nt,
        save_spacing=dtsave
    )  # initialize class

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
        set_yaml('DShoreface', brie._d_sf,
                 datadir)  # [m] Depth of shoreface (set to brie depth, function of wave height)
        set_yaml('LShoreface', float(brie._x_s[iB3D] - brie._x_t[iB3D]),
                 datadir)  # [m] Length of shoreface (calculate from brie variables, shoreline - shoreface toe)
        set_yaml('ShorefaceToe', float(brie._x_t[iB3D]), datadir)  # [m] Start location of shoreface toe
        # set_yaml('BermEl', 1.9 , datadir) # [m] Static elevation of berm (NOTE: if this is changed, the MSSM storm list and storm time series needs to be remade)
        set_yaml('BayDepth', brie._bb_depth, datadir)  # [m] Depth of bay behind island segment (set to brie bay depth)
        # set_yaml('MHW', 0.46, datadir)  # [m] Elevation of Mean High Water (NOTE: if this is changed, the storm time series needs to be remade)
        set_yaml('DuneParamStart', True, datadir)  # Dune height will come from external file
        set_yaml('Rat', 0.0,
                 datadir)  # [m / y] corresponds to Qat in LTA formulations (!!! must set to 0 because Qs is calculated in brie !!!)
        set_yaml('RSLR_Constant', True,
                 datadir)  # Relative sea-level rise rate will be constant, otherwise logistic growth function used for time series
        set_yaml('RSLR_const', brie._slr, datadir)  # [m / y] Relative sea-level rise rate
        # set_yaml('beta', 0.04, datadir)  # Beach slope for runup calculations
        set_yaml('k_sf', float(brie._k_sf),
                 datadir)  # [m^3 / m / y] Shoreface flux rate constant (function of wave parameters from brie)
        set_yaml('s_sf_eq', float(brie._s_sf_eq),
                 datadir)  # Equilibrium shoreface slope (function of wave and sediment parameters from brie)
        if np.size(rmin) > 1:
            set_yaml('GrowthParamStart', False, datadir)  # Dune growth parameter will come from external file
            set_yaml('rmin', rmin[iB3D], datadir)  # Minimum growth rate for logistic dune growth
            set_yaml('rmax', rmax[iB3D], datadir)  # Maximum growth rate for logistic dune growth
        else:
            set_yaml('GrowthParamStart', True, datadir)  # Dune growth parameter will come from external file
            set_yaml('rmin', rmin, datadir)  # Minimum growth rate for logistic dune growth
            set_yaml('rmax', rmax, datadir)  # Maximum growth rate for logistic dune growth

        barrier3d[iB3D].initialize(datadir)

        # debugging: check that the shoreface toe and shoreline are correct between the two models
        #print(brie._x_t[iB3D] - (barrier3d[iB3D]._model._x_t_TS[0] * 10))  # this isn't always zero; rounding error
        #print(brie._x_s[iB3D] - (barrier3d[iB3D]._model._x_s_TS[0] * 10))

        # now update brie x_b, x_b_save[:,0], h_b, h_b_save[:,0] from B3D so all the initial conditions are the same
        brie._x_b[iB3D] = (barrier3d[iB3D]._model._x_b_TS[0] * 10)  # the shoreline position + average interior width
        brie._h_b[iB3D] = (barrier3d[iB3D]._model._h_b_TS[0] * 10)  # average height of the interior domain
        brie._x_b_save[iB3D, 0] = brie._x_b[iB3D]
        brie._h_b_save[iB3D, 0] = brie._h_b[iB3D]

    return brie, barrier3d


###############################################################################
# time loop
###############################################################################

def time_loop(brie, barrier3d, num_cores):
    Time = time.time()

    # parallelize update function for each B3D sub-grid (spread overwash routing algorithms to different cores)
    def batchB3D(subB3D):

        subB3D.update()

        # get values for passing to brie (all in dam) [ UPDATE THIS WHEN BMI IS WORKING ]
        x_t_TS = subB3D._values["shoreface_toe_position"]()
        x_s_TS = subB3D._values["shoreline_position"]()
        x_b_TS = subB3D._values["back_barrier_shoreline_position"]()
        h_b_TS = subB3D._values["height_of_barrier"]()

        # calculate the diff in shoreface toe, shorelines (dam), height of barrier and convert to m
        sub_x_t_dt = (x_t_TS[-1] - x_t_TS[-2]) * 10
        sub_x_s_dt = (x_s_TS[-1] - x_s_TS[-2]) * 10
        sub_x_b_dt = (x_b_TS[-1] - x_b_TS[-2]) * 10
        sub_h_b_dt = (h_b_TS[-1] - h_b_TS[-2]) * 10

        return sub_x_t_dt, sub_x_s_dt, sub_x_b_dt, sub_h_b_dt

    for time_step in range(brie._nt - 1):

        # Print time step to screen (NOTE: time_index in each model is time_step+1)
        print("\r", 'Time Step: ', time_step, end="")

        # Advance B3D by one time step
        # --- (NOTE: B3D initializes at time_index = 1 and then updates the time_index at
        # the end of the update function. The dune and shoreline changes are saved at the beginning of the next time
        # step because they come from BRIE) ---

        batch_output = Parallel(n_jobs=num_cores)(
            delayed(batchB3D)(barrier3d[iB3D]) for iB3D in range(brie._ny)
        )  # set n_jobs=1 for no parallel processing (debugging) and -2 for all but 1 CPU

        # reshape output from parallel processing and convert from tuple to list
        x_t_dt, x_s_dt, x_b_dt, h_b_dt = zip(*batch_output)
        x_t_dt = list(x_t_dt)
        x_s_dt = list(x_s_dt)
        x_b_dt = list(x_b_dt)
        h_b_dt = list(h_b_dt)

        # pass values from B3D subdomains to brie for use in second time step
        # (there has to be a better way to do this with the BMI, but for now, access protected variables)
        brie._x_t_dt = x_t_dt
        brie._x_s_dt = x_s_dt
        brie._x_b_dt = x_b_dt
        brie._h_b_dt = h_b_dt

        # update brie one time step (this is time_index = 2 at start of loop)
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
    brie_save = []

    # save object without BMI
    for iB3D in range(brie._ny):
        b3d.append(barrier3d[iB3D]._model)
    brie_save = brie_save.append(brie)

    os.chdir(directory)
    np.savez(filename, barrier3d=b3d, brie=brie_save)

    return b3d


###############################################################################
# run brie with LTA for comparison
###############################################################################

def LTA(name, wave_height, wave_period, asym_frac, high_ang_frac, slr, ny, nt, w_b_crit, h_b_crit, Qow_max):

    # update the initial conditions
    ast_model = True  # shoreface formulations on
    barrier_model = True  # LTA14 overwash model on
    inlet_model = False  # inlet model off
    b3d = False  # B3d overwash model on

    # barrier model parameters
    s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
    z = 10.0  # initial sea level (for tracking SL, Eulerian reference frame)
    bb_depth = 3.0  # back-barrier depth

    # inlet parameters (use default; these are here to remind me later that they are important and I can change)
    Jmin = 10000  # minimum inlet spacing [m]
    a0 = 0.5  # amplitude of tide [m]
    marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

    # model setup
    dy = 100  # m, length of alongshore section (NOT the same as B3D, but overwash model performs better with dy=100 m)
    ny = ny * int(
        500 / dy)  # number of alongshore sections (NOTE, currently hard-coded for B3D dy = 500 m)
    dt = 0.05  # yr, timestep (NOT the same as B3D, but again, LTA14 performs better with dt = 0.05 yr)
    nt = int(nt / dt)  # equivalent timesteps to B3D
    dtsave = 20  # save spacing (equivalent of yearly for 0.05 time step)

    brieLTA = Brie(
        name=name,
        ast_model=ast_model,
        barrier_model=barrier_model,
        inlet_model=inlet_model,
        b3d=b3d,
        wave_height=wave_height,
        wave_period=wave_period,
        wave_asymmetry=asym_frac,
        wave_angle_high_fraction=high_ang_frac,
        sea_level_rise_rate=slr,
        sea_level_initial=z,
        barrier_height_critical=h_b_crit,
        barrier_width_critical=w_b_crit,
        max_overwash_flux=Qow_max,
        tide_amplitude=a0,
        back_barrier_marsh_fraction=marsh_cover,
        back_barrier_depth=bb_depth,
        xshore_slope=s_background,
        inlet_min_spacing=Jmin,
        alongshore_section_length=dy,
        alongshore_section_count=ny,
        time_step=dt,
        time_step_count=nt,
        save_spacing=dtsave
    )  # initialize class

    for time_step in range(int(brieLTA._nt) - 1):
        brieLTA.update()

    return brieLTA




