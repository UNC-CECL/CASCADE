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
import copy

from barrier3d import Barrier3d
from brie import Brie
from bulldozer import bulldoze, rebuild_dunes, set_growth_parameters


def batchB3D(subB3D):

    """parallelize update function for each B3D sub-grid (spread overwash routing algorithms to different cores)"""

    subB3D.update()

    # calculate the diff in shoreface toe, shoreline, and height of barrier (dam)
    sub_x_t_dt = (subB3D.x_t_TS[-1] - subB3D.x_t_TS[-2]) * 10
    sub_x_s_dt = (subB3D.x_s_TS[-1] - subB3D.x_s_TS[-2]) * 10
    sub_h_b_dt = (subB3D.h_b_TS[-1] - subB3D.h_b_TS[-2]) * 10

    return sub_x_t_dt, sub_x_s_dt, sub_h_b_dt, subB3D


def roadway_management(
    barrier3d,
    road_ele,
    road_width,
    road_setback,
    artificial_max_dune_ele,
    artificial_min_dune_ele,
    original_growth_param,
):

    dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters

    # if dune line moved by one grid cell, substract that amount from the setback distance to keep road in
    # the same place
    if (
        barrier3d.dune_migration
    ):  # check if the dune has been forced to migrate in B3D (NOTE: I think I may have been able to just use the B3D variable shoreline_change here!)
        road_setback = road_setback - 10

        # raise error if setback becomes negative (i.e., no longer a part of the island interior)
        if road_setback < 0:
            raise CascadeError(
                "Roadway drowned in the beach at t = {time} years".format(
                    time=barrier3d.time_index
                )
            )

    # bulldoze that road and put bulldozed sand back on the dunes
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
    road_overwash_removal = (
        road_overwash_removal * dm3_to_m3
    )  # convert from dm^3 to m^3

    # dune management: rebuild artifical dunes!
    if artificial_max_dune_ele is None or artificial_min_dune_ele is None:
        pass
    else:
        # in B3D, dune height is the height above the berm crest (keep this in m)
        artificial_max_dune_height = artificial_max_dune_ele - (barrier3d.BermEl * 10)
        artificial_min_dune_height = artificial_min_dune_ele - (barrier3d.BermEl * 10)

        # if any dune cell in the front row of dunes is less than a minimum threshold -- as measured above the berm
        # crest -- then rebuild artificial dunes (first row up to the max dune height, second row to the minimum?
        # -- nah, lets make it the max as well)
        if np.min(new_dune_domain[:, 0]) < (
            artificial_min_dune_height / 10
        ):  # dune domain is in decameters; previously had np.max() and got a lot of small dune breaches that didn't get
            # repaired -- not realistic
            new_dune_domain = rebuild_dunes(  # decameters
                new_dune_domain,
                max_dune_height=artificial_max_dune_height,  # in m
                min_dune_height=artificial_max_dune_height,  # in m
                dz=10,  # specifies dune domain is in decameters
                rng=True,  # adds stochasticity to dune height (seeded)
            )
            dunes_rebuilt = 1  # track when dunes are rebuilt
        else:
            dunes_rebuilt = 0

    # set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the natural eq.
    # dune height (Dmax)
    new_growth_parameters = set_growth_parameters(
        new_dune_domain,
        barrier3d.Dmax,
        barrier3d.growthparam,
        original_growth_param=original_growth_param,  # use original growth rates for resetting values
    )

    # update class variables
    barrier3d.DuneDomain[barrier3d.time_index - 1, :, :] = new_dune_domain
    barrier3d.InteriorDomain = new_xyz_interior_domain
    barrier3d.DomainTS[barrier3d.time_index - 1] = new_xyz_interior_domain
    barrier3d.growthparam = new_growth_parameters

    return (
        barrier3d,
        road_setback,
        dunes_rebuilt,
        road_overwash_removal,
    )


class CascadeError(Exception):
    pass


class Cascade:
    def __init__(
        self,
        datadir,
        name="default",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        alongshore_section_count=6,
        time_step_count=500,
        min_dune_growth_rate=0.35,
        max_dune_growth_rate=0.85,
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=True,
        road_ele=1.7,
        road_width=30,
        road_setback=30,
        artificial_max_dune_ele=3.7,
        artificial_min_dune_ele=2.2,
    ):
        """initialize both BRIE and Barrier3D

        Parameters
        ----------
        datadir: string
            Name of directory where Barrier3D input file is located
        name: string, optional
            Name of simulation
        wave_height: float, optional
            Mean offshore significant wave height [m].
        wave_period: float, optional
            Mean wave period [s].
        wave_asymmetry: float, optional
            Fraction of waves approaching from left (looking onshore).
        wave_angle_high_fraction: float, optional
            Fraction of waves approaching from angles higher than 45 degrees.
        sea_level_rise_rate: float, optional
            Rate of sea_level rise [m/yr].
        alongshore_section_count: int, optional
            Number of alongshore sections.
        time_step_count: int, optional
            Number of time steps.
        min_dune_growth_rate: float, optional
            Minimum dune growth rate [unitless, for Houser growth rate formulation]
        max_dune_growth_rate: float, optional
            Maximum dune growth rate [unitless, for Houser growth rate formulation]
        num_cores: int, optional
            Number of (parallel) processing cores to be used
        roadway_management: boolean, optional
            If True, use roadway management module (includes dune management)
        road_ele: float, optional
            Elevation of the roadway [m NAVD88]
        road_width: int, optional
            Width of roadway [m]
        road_setback: int, optional
            Setback of roadway from the inital dune line [m]
        artificial_max_dune_ele: float, optional
            Elevation to which dune is rebuilt to [m NAVD88]
        artificial_min_dune_ele: float, optional
            Elevation threshold which triggers rebuilding of dune [m NAVD88]

        Examples
        --------
        >>> from CASCADE import Cascade
        >>> datadir = "/Users/KatherineAnardeWheels/PycharmProjects/CASCADE/B3D_Inputs/"
        >>> cascade = Cascade(datadir)
        """

        self._ny = alongshore_section_count
        self._nt = time_step_count
        self._rmin = min_dune_growth_rate
        self._rmax = max_dune_growth_rate
        self._wave_height = wave_height
        self._wave_period = wave_period
        self._wave_assymetry = wave_asymmetry
        self._wave_angle_high_fraction = wave_angle_high_fraction
        self._sea_level_rise_rate = sea_level_rise_rate
        self._slr_constant = True
        self._num_cores = num_cores
        self._roadway_management_module = roadway_management_module
        self._alongshore_transport_module = alongshore_transport_module
        self._filename = name

        ###############################################################################
        # initial conditions for BRIE
        ###############################################################################

        # parameters that we need to initialize in BRIE for coupling (not necessarily default values), but I won't be
        # modifying often (or ever) for CASCADE
        self._brie_ast_model = True  # shoreface formulations on
        self._brie_barrier_model = False  # LTA14 overwash model off
        self._brie_inlet_model = False  # inlet model off
        self._b3d_barrier_model = True  # B3d overwash model on

        # barrier model parameters (the following are needed for other calculations even if the barrier model is off)
        self._s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
        self._z = 10.0  # initial sea level (for tracking SL, Eulerian reference frame)
        self._bb_depth = 3.0  # back-barrier depth
        self._h_b_crit = 1.9  # critical barrier height for overwash, used also to calculate shoreline diffusivity;
        # we set equal to the static elevation of berm (NOTE: if the berm elevation is changed, the MSSM storm list and
        # storm time series needs to be remade)

        # inlet parameters (use default; these are here to remind me later that they are important and I can change)
        self._Jmin = 10000  # minimum inlet spacing [m]
        self._a0 = 0.5  # amplitude of tide [m]
        self._marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism

        # grid and time step params
        self._dy = 500  # m, length of alongshore section (same as B3D)
        self._dt = 1  # yr, timestep (same as B3D)
        self._dtsave = 1  # save spacing (every year)

        # start by initializing BRIE b/c it has parameters related to wave climate that we use to initialize B3D
        self._brie = Brie(
            name=name,
            ast_model=self._brie_ast_model,
            barrier_model=self._brie_barrier_model,
            inlet_model=self._brie_inlet_model,
            b3d=self._b3d_barrier_model,
            wave_height=self._wave_height,
            wave_period=self._wave_period,
            wave_asymmetry=self._wave_assymetry,
            wave_angle_high_fraction=self._wave_angle_high_fraction,
            sea_level_rise_rate=self._sea_level_rise_rate,
            sea_level_initial=self._z,
            barrier_height_critical=self._h_b_crit,
            tide_amplitude=self._a0,
            back_barrier_marsh_fraction=self._marsh_cover,
            back_barrier_depth=self._bb_depth,
            xshore_slope=self._s_background,
            inlet_min_spacing=self._Jmin,
            alongshore_section_length=self._dy,
            alongshore_section_count=self._ny,
            time_step=self._dt,
            time_step_count=self._nt,
            save_spacing=self._dtsave,
        )  # initialize class

        ###############################################################################
        # initial conditions for Barrier3D
        ###############################################################################

        # for each B3D subgrid, set the initial shoreface geometry equal to what is set in brie (some random
        # perturbations); all other B3D variables are set equal
        barrier3d = []

        for iB3D in range(self._ny):

            fid = datadir + "barrier3d-parameters.yaml"

            # update yaml file (these are the only variables that I'm like to change from default)
            self.set_yaml(
                "TMAX", self._nt, fid
            )  # [yrs] Duration of simulation (if brie._dt = 1 yr, set to ._nt)
            self.set_yaml(
                "BarrierLength", self._dy, fid
            )  # [m] Static length of island segment (comprised of 10x10 cells)
            self.set_yaml(
                "DShoreface", self._brie.d_sf, fid
            )  # [m] Depth of shoreface (set to brie depth, function of wave height)
            self.set_yaml(
                "LShoreface", float(self._brie.x_s[iB3D] - self._brie.x_t[iB3D]), fid
            )  # [m] Length of shoreface (calculate from brie variables, shoreline - shoreface toe)
            self.set_yaml(
                "ShorefaceToe", float(self._brie.x_t[iB3D]), fid
            )  # [m] Start location of shoreface toe
            # set_yaml('BermEl', 1.9 , datadir) # [m] Static elevation of berm
            # (NOTE: if BermEl is changed, the MSSM storm list and storm time series needs to be remade)
            self.set_yaml(
                "BayDepth", self._bb_depth, fid
            )  # [m] Depth of bay behind island segment (set to brie bay depth)
            # set_yaml('MHW', 0.46, datadir)  # [m] Elevation of Mean High Water
            # (NOTE: if MHW is changed, the storm time series needs to be remade)
            self.set_yaml(
                "DuneParamStart", True, fid
            )  # Dune height will come from external file
            self.set_yaml(
                "Rat", 0.0, fid
            )  # [m / y] corresponds to Qat in LTA (!!! must set to 0 because Qs is calculated in brie !!!)
            self.set_yaml(
                "RSLR_Constant", self._slr_constant, fid
            )  # Relative sea-level rise rate will be constant, otherwise logistic growth function used
            self.set_yaml(
                "RSLR_const", self._sea_level_rise_rate, fid
            )  # [m / y] Relative sea-level rise rate
            # set_yaml('beta', 0.04, datadir)  # Beach slope for runup calculations
            self.set_yaml(
                "k_sf", float(self._brie.k_sf), fid
            )  # [m^3 / m / y] Shoreface flux rate constant (function of wave parameters from brie)
            self.set_yaml(
                "s_sf_eq", float(self._brie.s_sf_eq), fid
            )  # Equilibrium shoreface slope (function of wave and sediment parameters from brie)
            self.set_yaml(
                "GrowthParamStart", False, fid
            )  # Dune growth parameter WILL NOT come from external file
            if np.size(self._rmin) > 1:
                self.set_yaml(
                    "rmin", self._rmin[iB3D], fid
                )  # Minimum growth rate for logistic dune growth
                self.set_yaml(
                    "rmax", self._rmax[iB3D], fid
                )  # Maximum growth rate for logistic dune growth
            else:
                self.set_yaml(
                    "rmin", self._rmin, fid
                )  # Minimum growth rate for logistic dune growth
                self.set_yaml(
                    "rmax", self._rmax, fid
                )  # Maximum growth rate for logistic dune growth

            barrier3d.append(Barrier3d.from_yaml(datadir))

            # now update brie x_b, x_b_save[:,0], h_b, h_b_save[:,0] from B3D so all the initial conditions are the same
            # NOTE: interestingly here we don't need to have a "setter" in the property class for x_b, h_b, etc. because we
            #  are only replacing certain indices but added for completeness
            self._brie.x_b[iB3D] = (
                barrier3d[iB3D].x_b_TS[0] * 10
            )  # the shoreline position + average interior width
            self._brie.h_b[iB3D] = (
                barrier3d[iB3D].h_b_TS[0] * 10
            )  # average height of the interior domain
            self._brie.x_b_save[iB3D, 0] = self._brie.x_b[iB3D]
            self._brie.h_b_save[iB3D, 0] = self._brie.h_b[iB3D]

        # save array of Barrier3d models
        self._barrier3d = barrier3d

        ###############################################################################
        # initial conditions for human dynamics module
        ###############################################################################

        self._artificial_max_dune_ele = artificial_max_dune_ele
        self._artificial_min_dune_ele = artificial_min_dune_ele
        self._road_width = road_width
        self._road_ele = road_ele
        self._road_setback = [
            road_setback
        ] * self._ny  # needs to be an array due to variable shoreline change
        self._dunes_rebuilt = np.zeros(
            (self._ny, self._nt)
        )  # - 9999  # keep track of when dunes are rebuilt (boolean)
        self._road_overwash_volume = np.zeros(
            (self._ny, self._nt)
        )  # - 9999  # keep track of total overwash removed from roadway [m^3]
        self._post_storm_dunes = [
            None
        ] * self._nt  # keep track of post-storm dune impacts before humans
        self._post_storm_interior = [
            None
        ] * self._nt  # keep track of post-storm interior impacts before humans
        self._growth_params = [
            None
        ] * self._nt  # keep track of when dune growth parameters set to zero b/c of rebuild height
        self._growth_params[0] = [
            self._barrier3d[iB3D].growthparam for iB3D in range(self._ny)
        ]

    @classmethod
    def set_yaml(self, var_name, new_vals, file_name):
        with open(file_name) as f:
            doc = full_load(f)
        doc[var_name] = new_vals
        with open(file_name, "w") as f:
            dump(doc, f)

    @property
    def brie(self):
        return self._brie

    @brie.setter
    def brie(self, value):
        self._brie = value

    @property
    def barrier3d(self):
        return self._barrier3d

    @barrier3d.setter
    def barrier3d(self, value):
        self._barrier3d = value

    @property
    def roadway_management_module(self):
        return self._roadway_management_module

    @roadway_management_module.setter
    def roadway_management_module(self, value):
        self._roadway_management_module = value

    @property
    def post_storm_dunes(self):
        return self._post_storm_dunes

    @property
    def post_storm_interior(self):
        return self._post_storm_interior

    ###############################################################################
    # time loop
    ###############################################################################

    def update(self):
        """Update Cascade by a single time step"""

        # Because of the parallelization of b3d, I don't want to throw an error and exit if a barrier has drowned (will
        # lose the b3d model), so instead I just added a boolean to check for drowning here from the last time step in
        # brie. Note that this will stay false if brie is not used for AST (i.e., a B3D only run)
        if self._brie.drown == False:

            # Advance B3D by one time step; NOTE: B3D initializes at time_index = 1 and then updates the time_index
            # after update_dune_domain
            batch_output = Parallel(n_jobs=self._num_cores, max_nbytes="10M")(
                delayed(batchB3D)(self._barrier3d[iB3D]) for iB3D in range(self._ny)
            )  # set n_jobs=1 for no parallel processing (debugging) and -2 for all but 1 CPU; note that joblib uses a
            # threshold on the size of arrays passed to the workers; we use 'None' to disable memory mapping of large arrays

            # reshape output from parallel processing and convert from tuple to list
            x_t_dt, x_s_dt, h_b_dt, b3d = zip(*batch_output)
            x_t_dt = list(x_t_dt)
            x_s_dt = list(x_s_dt)
            h_b_dt = list(h_b_dt)
            self._barrier3d = list(b3d)

            # use brie to connect B3D subgrids with alongshore sediment transport
            if self._alongshore_transport_module:

                # pass shoreline and shoreface values from B3D subdomains to brie for use in second time step
                self._brie.x_t_dt = x_t_dt
                self._brie.x_s_dt = x_s_dt
                self._brie.x_b_dt = 0  # dummy variable, will set x_b later
                self._brie.h_b_dt = h_b_dt

                # update brie one time step (this is time_index = 2 at start of loop)
                self._brie.update()

                for iB3D in range(self._ny):
                    # pass shoreline position back to B3D from Brie (convert from m to dam)
                    self._barrier3d[iB3D].x_s = self._brie.x_s[iB3D] / 10
                    self._barrier3d[iB3D].x_s_TS[-1] = self._brie.x_s[iB3D] / 10

                    # update dune domain in B3D (erode/prograde) based on shoreline change from Brie
                    self._barrier3d[iB3D].update_dune_domain()

                    # update back-barrier shoreline location in BRIE based on new shoreline + average interior width in B3D
                    self._brie.x_b[iB3D] = self._barrier3d[iB3D].x_b_TS[-1] * 10
                    self._brie.x_b_save[
                        iB3D, self._brie.time_index - 1
                    ] = self._brie.x_b[iB3D]
            else:

                # no alongshore sediment transport; update dune domain
                for iB3D in range(self._ny):
                    self._barrier3d[iB3D].update_dune_domain()

            # human dynamics: remove overwash from roadway after each model year, place on the dune, rebuild dunes if
            # necessary, and check if dunes should grow naturally
            if self._roadway_management_module:

                post_storm_interior = []
                post_storm_dunes = []
                new_growth_params = []

                for iB3D in range(self._ny):

                    # save post-storm dune and interior domain before human modificiations
                    post_storm_interior.append(
                        copy.deepcopy(self._barrier3d[iB3D].InteriorDomain)
                    )  # hoping this makes a deep copy
                    post_storm_dunes.append(
                        copy.deepcopy(
                            self._barrier3d[iB3D].DuneDomain[
                                self._barrier3d[iB3D].time_index - 1, :, :
                            ]
                        )
                    )

                    # call roadway manamgement module
                    (
                        self._barrier3d[iB3D],
                        self._road_setback[iB3D],
                        track_dune_rebuild,
                        road_overwash_removal,
                    ) = roadway_management(
                        self._barrier3d[iB3D],
                        self._road_ele,
                        self._road_width,
                        self._road_setback[iB3D],
                        self._artificial_max_dune_ele,
                        self._artificial_min_dune_ele,
                        self._growth_params[0][
                            iB3D
                        ],  # original dune growth rates for resetting values
                    )

                    # keep track of when dunes rebuilt, road overwash volume, and new growth rates
                    self._dunes_rebuilt[
                        iB3D, self._barrier3d[iB3D].time_index - 1
                    ] = track_dune_rebuild
                    self._road_overwash_volume[
                        iB3D, self._barrier3d[iB3D].time_index - 1
                    ] = road_overwash_removal  # total overwash removed from roadway [m^3]
                    new_growth_params.append(
                        copy.deepcopy(self._barrier3d[iB3D].growthparam)
                    )

                # save lists for each time step to preallocated variables
                self._post_storm_interior[
                    self._barrier3d[0].time_index - 1
                ] = post_storm_interior
                self._post_storm_dunes[
                    self._barrier3d[0].time_index - 1
                ] = post_storm_dunes
                self._growth_params[
                    self._barrier3d[0].time_index - 1
                ] = new_growth_params

    ###############################################################################
    # save data
    ###############################################################################

    def save(self, directory):

        filename = self._filename + ".npz"

        csc8d = []
        csc8d.append(self)

        # b3d = []
        # brie_save = []
        #
        # # save object
        # for iB3D in range(self._ny):
        #     b3d.append(self._barrier3d[iB3D])
        #
        # brie_save.append(self._brie)

        os.chdir(directory)
        np.savez(filename, cascade=csc8d)
        # np.savez(filename, barrier3d=b3d, brie=brie_save)

        # return b3d


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
