# Model coupling of Barrier3D and BRIE

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
----------------------------------------------------"""

from yaml import full_load, dump
from joblib import Parallel, delayed
import numpy as np
import os

from barrier3d import Barrier3d
from brie import Brie

from roadway_manager import RoadwayManager
from alongshore_coupler import AlongshoreCoupler


def batchB3D(subB3D):

    """parallelize update function for each B3D sub-grid (spread overwash routing algorithms to different cores)"""

    subB3D.update()

    # calculate the diff in shoreface toe, shoreline, and height of barrier (dam)
    sub_x_t_dt = (subB3D.x_t_TS[-1] - subB3D.x_t_TS[-2]) * 10
    sub_x_s_dt = (subB3D.x_s_TS[-1] - subB3D.x_s_TS[-2]) * 10
    sub_h_b_dt = (subB3D.h_b_TS[-1] - subB3D.h_b_TS[-2]) * 10

    return sub_x_t_dt, sub_x_s_dt, sub_h_b_dt, subB3D


class CascadeError(Exception):
    pass


class Cascade:
    def __init__(
        self,
        datadir,
        name="default",
        storm_file="barrier3d-default-storms.npy",
        elevation_file="barrier3d-default-elevation.npy",
        dune_file="barrier3d-default-dunes.npy",
        parameter_file="barrier3d-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        alongshore_section_count=6,
        time_step_count=200,
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
        self._storm_file = storm_file
        self._elevation_file = elevation_file
        self._dune_file = dune_file
        self._parameter_file = parameter_file

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

            fid = datadir + self._parameter_file

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

            # external file names used for initialization
            self.set_yaml("storm_file", self._storm_file, fid)
            self.set_yaml("dune_file", self._dune_file, fid)
            self.set_yaml("elevation_file", self._elevation_file, fid)

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
        # initial conditions for human dynamics module, ast module
        ###############################################################################

        self._artificial_max_dune_ele = artificial_max_dune_ele
        self._artificial_min_dune_ele = artificial_min_dune_ele
        self._road_width = road_width
        self._road_ele = road_ele
        self._orig_road_setback = road_setback
        self._road_relocation_setback = road_setback
        self._road_break = 0  # will = 1 if roadway drowns

        # initialize RoadwayManager (always, just in case we want to add a road during the simulation)
        self._roadways = []

        for iB3D in range(self._ny):
            self._roadways.append(
                RoadwayManager(
                    road_elevation=self._road_ele,
                    road_width=self._road_width,
                    road_setback=self._orig_road_setback,
                    road_relocation_setback=self._road_relocation_setback,
                    dune_design_height=self._artificial_max_dune_ele,
                    dune_minimum_height=self._artificial_min_dune_ele,
                    time_step_count=self._nt,
                    original_growth_param=self._barrier3d[iB3D].growthparam,
                )
            )

        if self._alongshore_transport_module:
            self._ast_coupler = AlongshoreCoupler(self._ny)

    @classmethod
    def set_yaml(self, var_name, new_vals, file_name):
        with open(file_name) as f:
            doc = full_load(f)
        doc[var_name] = new_vals
        with open(file_name, "w") as f:
            dump(doc, f)

    @property
    def road_break(self):
        return self._road_break

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
    def roadways(self):
        return self._roadways

    @roadways.setter
    def roadways(self, value):
        self._roadways = value

    @property
    def roadway_management_module(self):
        return self._roadway_management_module

    @roadway_management_module.setter
    def roadway_management_module(self, value):
        self._roadway_management_module = value

    ###############################################################################
    # time loop
    ###############################################################################

    def update(self):
        """Update Cascade by a single time step"""

        # Check for drowning here from the last time step in brie. Note that this will stay false if brie is not used
        # for AST (i.e., a B3D only run).
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

            # use brie to connect B3D subgrids with alongshore sediment transport; otherwise, just update dune domain
            if self._alongshore_transport_module:
                self._ast_coupler.update(
                    self._brie, self._barrier3d, x_t_dt, x_s_dt, h_b_dt
                )
            else:
                for iB3D in range(self._ny):
                    self._barrier3d[iB3D].update_dune_domain()

            # human dynamics: remove overwash from roadway after each model year, place on the dune, rebuild dunes if
            # necessary, and check if dunes should grow naturally
            if self._roadway_management_module:
                for iB3D in range(self._ny):
                    # if the user changed any of the roadway or dune parameters, update them in the module before advancing
                    self._roadways[iB3D].road_ele = self._road_ele
                    self._roadways[iB3D].road_width = self._road_width
                    self._roadways[
                        iB3D
                    ].road_relocation_setback = self._road_relocation_setback
                    self._roadways[
                        iB3D
                    ].artificial_max_dune_ele = self._artificial_max_dune_ele
                    self._roadways[
                        iB3D
                    ].artificial_min_dune_ele = self._artificial_min_dune_ele

                    # call roadway manamgement module
                    self._roadways[iB3D].update(self._barrier3d[iB3D])

                    # if the road drowned or barrier was too narrow to be relocated, break
                    if (
                        self._roadways[iB3D].drown_break
                        or self._roadways[iB3D].narrow_break
                    ):
                        self._road_break = 1
                        return

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


# def LTA(
#     name,
#     wave_height,
#     wave_period,
#     asym_frac,
#     high_ang_frac,
#     slr,
#     ny,
#     nt,
#     w_b_crit,
#     h_b_crit,
#     Qow_max,
# ):
#     # update the initial conditions
#     ast_model = True  # shoreface formulations on
#     barrier_model = True  # LTA14 overwash model on
#     inlet_model = False  # inlet model off
#     b3d = False  # B3d overwash model on
#
#     # barrier model parameters
#     s_background = 0.001  # background slope (for shoreface toe position, back-barrier & inlet calculations)
#     z = 10.0  # initial sea level (for tracking SL, Eulerian reference frame)
#     bb_depth = 3.0  # back-barrier depth
#
#     # inlet parameters (use default; these are here to remind me later that they are important and I can change)
#     Jmin = 10000  # minimum inlet spacing [m]
#     a0 = 0.5  # amplitude of tide [m]
#     marsh_cover = 0.5  # % of backbarrier covered by marsh and therefore does not contribute to tidal prism
#
#     # model setup
#     dy = 100  # m, length of alongshore section (NOT the same as B3D, but overwash model performs better with dy=100 m)
#     ny = ny * int(
#         500 / dy
#     )  # number of alongshore sections (NOTE, currently hard-coded for B3D dy = 500 m)
#     dt = 0.05  # yr, timestep (NOT the same as B3D, but again, LTA14 performs better with dt = 0.05 yr)
#     nt = int(nt / dt)  # equivalent timesteps to B3D
#     dtsave = int(1 / dt)  # save spacing (equivalent of yearly for 0.05 time step)
#
#     brieLTA = Brie(
#         name=name,
#         ast_model=ast_model,
#         barrier_model=barrier_model,
#         inlet_model=inlet_model,
#         b3d=b3d,
#         wave_height=wave_height,
#         wave_period=wave_period,
#         wave_asymmetry=asym_frac,
#         wave_angle_high_fraction=high_ang_frac,
#         sea_level_rise_rate=slr,
#         sea_level_initial=z,
#         barrier_height_critical=h_b_crit,
#         barrier_width_critical=w_b_crit,
#         max_overwash_flux=Qow_max,
#         tide_amplitude=a0,
#         back_barrier_marsh_fraction=marsh_cover,
#         back_barrier_depth=bb_depth,
#         xshore_slope=s_background,
#         inlet_min_spacing=Jmin,
#         alongshore_section_length=dy,
#         alongshore_section_count=ny,
#         time_step=dt,
#         time_step_count=nt,
#         save_spacing=dtsave,
#     )  # initialize class
#
#     for time_step in range(int(brieLTA.nt) - 1):
#         brieLTA.update()
#
#     return brieLTA
