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

from barrier3d import Barrier3d
from brie import Brie  # need to update this for the new BMI (i.e., BrieBMI)
from bulldozer import bulldoze, rebuild_dunes

# temporary placement - I want this added to the BMI but keep getting an object error
def set_yaml(var_name, new_vals, file_name):
    with open(file_name) as f:
        doc = full_load(f)

    doc[var_name] = new_vals

    with open(file_name, "w") as f:
        dump(doc, f)


###############################################################################
# initial conditions for Brie
###############################################################################


def initialize(
    datadir,
    name,
    wave_height=1,
    wave_period=7,
    asym_frac=0.8,
    high_ang_frac=0.2,
    slr=0.004,
    ny=6,
    nt=500,
    rmin=0.35,
    rmax=0.85,
):
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
        save_spacing=dtsave,
    )  # initialize class

    ###############################################################################
    # initial conditions for Barrier3D
    ###############################################################################

    # for each B3D subgrid, set the initial shoreface geometry equal to what is set in brie (some random perturbations);
    # all other B3D variables are set equal
    barrier3d = []

    for iB3D in range(brie.ny):

        fid = datadir + "barrier3d-parameters.yaml"

        # update yaml file (these are the only variables that I'm like to change from default)
        set_yaml(
            "TMAX", nt, fid
        )  # [yrs] Duration of simulation (if brie._dt = 1 yr, set to ._nt)
        set_yaml(
            "BarrierLength", dy, fid
        )  # [m] Static length of island segment (comprised of 10x10 cells)
        set_yaml(
            "DShoreface", brie.d_sf, fid
        )  # [m] Depth of shoreface (set to brie depth, function of wave height)
        set_yaml(
            "LShoreface", float(brie.x_s[iB3D] - brie.x_t[iB3D]), fid
        )  # [m] Length of shoreface (calculate from brie variables, shoreline - shoreface toe)
        set_yaml(
            "ShorefaceToe", float(brie.x_t[iB3D]), fid
        )  # [m] Start location of shoreface toe
        # set_yaml('BermEl', 1.9 , datadir) # [m] Static elevation of berm
        # (NOTE: if BermEl is changed, the MSSM storm list and storm time series needs to be remade)
        set_yaml(
            "BayDepth", bb_depth, fid
        )  # [m] Depth of bay behind island segment (set to brie bay depth)
        # set_yaml('MHW', 0.46, datadir)  # [m] Elevation of Mean High Water
        # (NOTE: if MHW is changed, the storm time series needs to be remade)
        set_yaml(
            "DuneParamStart", True, fid
        )  # Dune height will come from external file
        set_yaml(
            "Rat", 0.0, fid
        )  # [m / y] corresponds to Qat in LTA (!!! must set to 0 because Qs is calculated in brie !!!)
        set_yaml(
            "RSLR_Constant", True, fid
        )  # Relative sea-level rise rate will be constant, otherwise logistic growth function used
        set_yaml("RSLR_const", slr, fid)  # [m / y] Relative sea-level rise rate
        # set_yaml('beta', 0.04, datadir)  # Beach slope for runup calculations
        set_yaml(
            "k_sf", float(brie.k_sf), fid
        )  # [m^3 / m / y] Shoreface flux rate constant (function of wave parameters from brie)
        set_yaml(
            "s_sf_eq", float(brie.s_sf_eq), fid
        )  # Equilibrium shoreface slope (function of wave and sediment parameters from brie)
        set_yaml(
            "GrowthParamStart", False, fid
        )  # Dune growth parameter WILL NOT come from external file
        if np.size(rmin) > 1:
            set_yaml(
                "rmin", rmin[iB3D], fid
            )  # Minimum growth rate for logistic dune growth
            set_yaml(
                "rmax", rmax[iB3D], fid
            )  # Maximum growth rate for logistic dune growth
        else:
            set_yaml("rmin", rmin, fid)  # Minimum growth rate for logistic dune growth
            set_yaml("rmax", rmax, fid)  # Maximum growth rate for logistic dune growth

        barrier3d.append(Barrier3d.from_yaml(datadir))

        # now update brie x_b, x_b_save[:,0], h_b, h_b_save[:,0] from B3D so all the initial conditions are the same
        # NOTE: interestingly here we don't need to have a "setter" in the property class for x_b, h_b, etc. because we
        #  are only replacing certain indices but added for completeness
        brie.x_b[iB3D] = (
            barrier3d[iB3D].x_b_TS[0] * 10
        )  # the shoreline position + average interior width
        brie.h_b[iB3D] = (
            barrier3d[iB3D].h_b_TS[0] * 10
        )  # average height of the interior domain
        brie.x_b_save[iB3D, 0] = brie.x_b[iB3D]
        brie.h_b_save[iB3D, 0] = brie.h_b[iB3D]

    return brie, barrier3d


###############################################################################
# time loop
###############################################################################


def time_loop(
    brie,
    barrier3d,
    num_cores,
    nt,
    road_ele=None,
    road_width=None,
    road_setback=None,
    artificial_max_dune_ele=None,
    artificial_min_dune_ele=None,
):
    Time = time.time()
    road_setback = [
        road_setback
    ] * brie.ny  # create an array for variable road setback distances for the human dynamics module

    # parallelize update function for each B3D sub-grid (spread overwash routing algorithms to different cores)
    def batchB3D(subB3D):

        subB3D.update()

        # calculate the diff in shoreface toe, shoreline, and height of barrier (dam)
        sub_x_t_dt = (subB3D.x_t_TS[-1] - subB3D.x_t_TS[-2]) * 10
        sub_x_s_dt = (subB3D.x_s_TS[-1] - subB3D.x_s_TS[-2]) * 10
        sub_h_b_dt = (subB3D.h_b_TS[-1] - subB3D.h_b_TS[-2]) * 10

        return sub_x_t_dt, sub_x_s_dt, sub_h_b_dt, subB3D

    for time_step in range(nt - 1):

        # because of the parallelization of b3d, I don't want to throw an error and exit if a barrier has drowned (will
        # lose the b3d model), so instead I just added a boolean to check for drowning (might need to do this in B3D too)
        if brie.drown == False:

            # Print time step to screen (NOTE: time_index in each model starts at 1)
            print("\r", "Time Step: ", brie.time_index + 1, end="")

            # Advance B3D by one time step
            # NOTE: B3D initializes at time_index = 1 and then updates the time_index after update_dune_domain

            batch_output = Parallel(n_jobs=num_cores, max_nbytes="10M")(
                delayed(batchB3D)(barrier3d[iB3D]) for iB3D in range(brie.ny)
            )  # set n_jobs=1 for no parallel processing (debugging) and -2 for all but 1 CPU; note that joblib uses a
            # threshold on the size of arrays passed to the workers; we use 'None' to disable memory mapping of large arrays

            # reshape output from parallel processing and convert from tuple to list
            x_t_dt, x_s_dt, h_b_dt, b3d = zip(*batch_output)
            x_t_dt = list(x_t_dt)
            x_s_dt = list(x_s_dt)
            h_b_dt = list(h_b_dt)
            barrier3d = list(b3d)

            # pass shoreline and shoreface values from B3D subdomains to brie for use in second time step
            brie.x_t_dt = x_t_dt
            brie.x_s_dt = x_s_dt
            brie.x_b_dt = 0  # dummy variable, will set x_b later
            brie.h_b_dt = h_b_dt

            # update brie one time step (this is time_index = 2 at start of loop)
            brie.update()

            for iB3D in range(brie.ny):

                # pass shoreline position back to B3D from Brie (convert from m to dam)
                barrier3d[iB3D].x_s = brie.x_s[iB3D] / 10
                barrier3d[iB3D].x_s_TS[-1] = brie.x_s[iB3D] / 10

                # update dune domain in B3D (erode/prograde) based on shoreline change from Brie
                barrier3d[iB3D].update_dune_domain()

                # update back-barrier shoreline location in BRIE based on new shoreline + average interior width in B3D
                brie.x_b[iB3D] = barrier3d[iB3D].x_b_TS[-1] * 10
                brie.x_b_save[iB3D, brie.time_index - 1] = brie.x_b[iB3D]

        else:
            # if the barrier drowns in brie, exit the time loop (although, I expect this to happen in B3D first)
            break

        # human dynamics: remove overwash from 16 m-wide section at the start of the island interior -- corresponding to
        # the roadway -- after each model year and place that material on the dune line (only affects B3D)
        if road_ele is None or road_width is None:
            pass
        else:
            # call human dynamics module
            for iB3D in range(brie.ny):

                barrier3d[iB3D], road_setback[iB3D] = roadway_management(
                    barrier3d[iB3D],
                    road_ele,
                    road_width,
                    road_setback[iB3D],
                    artificial_max_dune_ele,
                    artificial_min_dune_ele,
                )

    SimDuration = time.time() - Time
    print()
    print("Elapsed Time: ", SimDuration, "sec")  # Print elapsed time of simulation

    return brie, barrier3d, road_setback


###############################################################################
# save data
###############################################################################


def save(brie, barrier3d, directory, name):

    filename = name + ".npz"

    b3d = []
    brie_save = []

    # save object
    for iB3D in range(brie.ny):
        b3d.append(barrier3d[iB3D])

    brie_save.append(brie)

    os.chdir(directory)
    np.savez(filename, barrier3d=b3d, brie=brie_save)

    return b3d


###############################################################################
# run brie with LTA for comparison
###############################################################################


def LTA(
    name,
    wave_height,
    wave_period,
    asym_frac,
    high_ang_frac,
    slr,
    ny,
    nt,
    w_b_crit,
    h_b_crit,
    Qow_max,
):
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
        500 / dy
    )  # number of alongshore sections (NOTE, currently hard-coded for B3D dy = 500 m)
    dt = 0.05  # yr, timestep (NOT the same as B3D, but again, LTA14 performs better with dt = 0.05 yr)
    nt = int(nt / dt)  # equivalent timesteps to B3D
    dtsave = int(1 / dt)  # save spacing (equivalent of yearly for 0.05 time step)

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
        save_spacing=dtsave,
    )  # initialize class

    for time_step in range(int(brieLTA.nt) - 1):
        brieLTA.update()

    return brieLTA


###############################################################################
# run human dynamics modules
###############################################################################


def roadway_management(
    barrier3d,
    road_ele,
    road_width,
    road_setback,
    artificial_max_dune_ele,
    artificial_min_dune_ele,
):

    # if dune line moved by one grid cell, substract that amount from the setback distance to keep road in
    # the same place
    if barrier3d.dune_migration:  # check if the dune has been forced to migrate in B3D
        road_setback = road_setback - 10

    # bulldoze that road!
    new_dune_domain, new_xyz_interior_domain, road_overwash_removal = bulldoze(
        road_ele=road_ele,  # m NAVD88
        road_width=road_width,
        road_setback=road_setback,
        xyz_interior_grid=barrier3d.InteriorDomain,  # interior domain from this last time step
        yxz_dune_grid=barrier3d.DuneDomain[
            barrier3d.time_index - 1, :, :
        ],  # dune domain from this last time step
        dx=10,
        dy=10,
        dz=10,
    )

    # rebuild artifical dunes
    if artificial_max_dune_ele is None or artificial_min_dune_ele is None:
        pass
    else:
        # in B3D, dune height is the height above the berm crest (keep this in m)
        artificial_max_dune_height = artificial_max_dune_ele - (barrier3d.BermEl * 10)
        artificial_min_dune_height = artificial_min_dune_ele - (barrier3d.BermEl * 10)

        # if front row of dunes is less than the height that dunes are typically rebuilt
        # to -- as measured above the berm crest -- then rebuild artificial dunes back to this height
        if np.max(new_dune_domain[:, 0]) < (artificial_max_dune_height / 10):
            new_dune_domain = rebuild_dunes(
                new_dune_domain,
                max_dune_height=artificial_max_dune_height,
                min_dune_height=artificial_min_dune_height,
                dz=10,
                rng=True,
            )

    # update class variables
    barrier3d.DuneDomain[barrier3d.time_index - 1, :, :] = new_dune_domain
    barrier3d.InteriorDomain = new_xyz_interior_domain
    barrier3d.DomainTS[barrier3d.time_index - 1] = new_xyz_interior_domain

    # set dune growth rate to zero for next time step (turn off natural dynamics)
    barrier3d.growthparam = np.zeros([1, np.size(barrier3d.growthparam)])

    return barrier3d, road_setback
