# Model coupling of Barrier3D and BRIE

# ~******* CASCADE ********~

"""----------------------------------------------------
Copyright (C) 2020 Katherine Anarde
----------------------------------------------------"""

from pathlib import Path
from joblib import Parallel, delayed
import numpy as np
import os

from roadway_manager import RoadwayManager
from beach_nourisher import BeachNourisher
from brie_coupler import BrieCoupler, initialize_equal, batchB3D
from chome_coupler import ChomeCoupler


class CascadeError(Exception):
    pass


class Cascade:
    def module_lists(
        self,
        artificial_max_dune_ele,
        artificial_min_dune_ele,
        road_width,
        road_ele,
        road_setback,
        nourishment_interval,
        nourishment_volume,
    ):
        """Configures lists to account for multiple barrier3d domains from single input variables; used in modules"""

        if np.size(artificial_max_dune_ele) > 1:
            self._artificial_max_dune_ele = (
                artificial_max_dune_ele  # list of floats option
            )
        else:
            self._artificial_max_dune_ele = [artificial_max_dune_ele] * self._ny
        if np.size(artificial_min_dune_ele) > 1:
            self._artificial_min_dune_ele = artificial_min_dune_ele
        else:
            self._artificial_min_dune_ele = [artificial_min_dune_ele] * self._ny
        if np.size(road_width) > 1:
            self._road_width = road_width
        else:
            self._road_width = [road_width] * self._ny
        if np.size(road_ele) > 1:
            self._road_ele = road_ele
        else:
            self._road_ele = [road_ele] * self._ny
        if np.size(road_setback) > 1:
            self._orig_road_setback = road_setback
            self._road_relocation_setback = road_setback
        else:
            self._orig_road_setback = [road_setback] * self._ny
            self._road_relocation_setback = [road_setback] * self._ny
        if np.size(nourishment_interval) > 1:
            self._nourishment_interval = nourishment_interval
        else:
            self._nourishment_interval = [nourishment_interval] * self._ny
        if np.size(nourishment_volume) > 1:
            self._nourishment_volume = nourishment_volume
        else:
            self._nourishment_volume = [nourishment_volume] * self._ny

        return

    def __init__(
        self,
        datadir,
        name="default",
        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="b3d_pt45_8750yrs_low-elevations.csv",  # associated with average dune growth rate of 0.45
        dune_file="barrier3d-default-dunes.npy",
        parameter_file="barrier3d-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        alongshore_section_count=6,
        time_step_count=200,
        min_dune_growth_rate=0.25,  # average is 0.45, low dune growth rate
        max_dune_growth_rate=0.65,
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=True,
        beach_nourishment_module=True,
        community_dynamics_module=False,
        road_ele=1.7,
        road_width=30,
        road_setback=30,
        artificial_max_dune_ele=3.7,
        artificial_min_dune_ele=2.2,
        nourishment_interval=None,
        nourishment_volume=300.0,
        number_of_communities=1,
        sand_cost=10,
        taxratio_oceanfront=1,
        external_housing_market_value_oceanfront=6e5,
        external_housing_market_value_nonoceanfront=4e5,
        fixed_cost_beach_nourishment=2e6,
        fixed_cost_dune_nourishment=2e5,
        nourishment_cost_subsidy=10e6,
        house_footprint=15,
    ):
        """initialize models (Barrier3D, BRIE, CHOME) and human dynamics modules

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
        sea_level_rise_constant: boolean, optional
            Rate of sea_level rise is constant if True, otherwise a logistic growth function is used.
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
        roadway_management_module: boolean, optional
            If True, use roadway management module (overwash removal, road relocation, dune management)
        alongshore_transport_module: boolean, optional
            If True, couple Barrier3D with BRIE to use diffusive model for AST
        community_dynamics_module: boolean, optional
            If True, couple with CHOME, a community decision making model; requires nourishment module
        beach_nourishment_module: boolean, optional
            If True, use nourishment module (nourish shoreface, rebuild dunes)
        road_ele: float or list of floats, optional
            Elevation of the roadway [m NAVD88]
        road_width: int or list of int, optional
            Width of roadway [m]
        road_setback: int or list of int, optional
            Setback of roadway from the inital dune line [m]
        artificial_max_dune_ele: float or list of floats, optional
            Elevation to which dune is rebuilt to [m NAVD88]
        artificial_min_dune_ele: float or list of floats, optional
            Elevation threshold which triggers rebuilding of dune [m NAVD88]
        nourishment_interval: int or list of ints, optional
             Interval that nourishment occurs [yrs]
        nourishment_volume: float or list of float, optional
             Volume of nourished sand along cross-shore transect [m^3/m]
        number_of_communities: int, optional
            Number of communities (CHOME model instances) described by the alongshore section count (Barrier3D grids)
        sand_cost: int, optional
            Unit cost of sand $/m^3
        taxratio_oceanfront: float, optional
            The proportion of tax burden placed on oceanfront row for beach management
        external_housing_market_value_oceanfront:  , optional
            Value of comparable housing option outside the coastal system
        external_housing_market_value_nonoceanfront:  , optional
            Value of comparable housing option outside the coastal system
        fixed_cost_beach_nourishment: int, optional
            Fixed cost of 1 nourishment project
        fixed_cost_dune_nourishment: int, optional
            Fixed cost of building dunes once
        nourishment_cost_subsidy: int, optional
            Subsidy on cost of entire nourishment plan

        Examples
        --------
        >>> from CASCADE import Cascade
        >>> datadir = "./B3D_Inputs/"
        >>> cascade = Cascade(datadir)
        """

        self._ny = alongshore_section_count
        self._nt = time_step_count
        self._rmin = min_dune_growth_rate
        self._rmax = max_dune_growth_rate
        self._wave_height = wave_height
        self._wave_period = wave_period
        self._wave_asymmetry = wave_asymmetry
        self._wave_angle_high_fraction = wave_angle_high_fraction
        self._sea_level_rise_rate = sea_level_rise_rate
        self._slr_constant = sea_level_rise_constant
        self._num_cores = num_cores
        self._roadway_management_module = roadway_management_module
        self._alongshore_transport_module = alongshore_transport_module
        self._beach_nourishment_module = beach_nourishment_module
        self._community_dynamics_module = community_dynamics_module
        self._filename = name
        self._storm_file = storm_file
        self._elevation_file = elevation_file
        self._dune_file = dune_file
        self._parameter_file = parameter_file

        ###############################################################################
        # initialize brie and barrier3d classes
        ###############################################################################

        self._brie_coupler = BrieCoupler(
            name,
            self._wave_height,
            self._wave_period,
            self._wave_asymmetry,
            self._wave_angle_high_fraction,
            self._sea_level_rise_rate,
            self._ny,
            self._nt,
        )

        self._barrier3d, self._brie_coupler._brie = initialize_equal(
            datadir,
            self._brie_coupler._brie,
            self._slr_constant,
            self._rmin,
            self._rmax,
            self._parameter_file,
            self._storm_file,
            self._dune_file,
            self._elevation_file,
        )

        ###############################################################################
        # initialize human dynamics modules
        ###############################################################################

        # configure self to create lists of these variables
        self.module_lists(
            artificial_max_dune_ele,
            artificial_min_dune_ele,
            road_width,
            road_ele,
            road_setback,
            nourishment_interval,
            nourishment_volume,
        )

        self._number_of_communities = number_of_communities
        self._b3d_break = 0  # true if barrier in Barrier3D height or width drowns
        self._road_break = 0  # true if roadway drowns
        self._nourish_now = [0] * self._ny  # boolean for triggering nourishment
        self._rebuild_dune_now = [
            0
        ] * self._ny  # boolean for triggering dune rebuilding

        if self._community_dynamics_module:
            if self._beach_nourishment_module == False:
                CascadeError(
                    "Beach nourishment module must be set to `TRUE` to couple with CHOME"
                )
            else:
                self._chome_coupler = ChomeCoupler(
                    barrier3d=self._barrier3d,
                    total_time=self._nt,
                    alongshore_length_b3d=self._brie_coupler._brie._dy,  # this is the barrier3d default, 500 m
                    dune_design_elevation=self._artificial_max_dune_ele,
                    number_of_communities=self._number_of_communities,
                    name=self._filename,
                    sand_cost=sand_cost,
                    taxratio_oceanfront=taxratio_oceanfront,
                    external_housing_market_value_oceanfront=external_housing_market_value_oceanfront,
                    external_housing_market_value_nonoceanfront=external_housing_market_value_nonoceanfront,
                    fixed_cost_beach_nourishment=fixed_cost_beach_nourishment,
                    fixed_cost_dune_nourishment=fixed_cost_dune_nourishment,
                    nourishment_cost_subsidy=nourishment_cost_subsidy,
                    house_footprint=house_footprint,
                )  # contains the CHOME model instances, one per community

        # initialize RoadwayManager and BeachNourisher modules (always, just in case we want to add a road or start
        # nourishing during the simulation)
        self._roadways = []
        self._nourishments = []

        for iB3D in range(self._ny):
            self._roadways.append(
                RoadwayManager(
                    road_elevation=self._road_ele[iB3D],
                    road_width=self._road_width[iB3D],
                    road_setback=self._orig_road_setback[iB3D],
                    road_relocation_setback=self._road_relocation_setback[iB3D],
                    dune_design_elevation=self._artificial_max_dune_ele[iB3D],
                    dune_minimum_elevation=self._artificial_min_dune_ele[iB3D],
                    time_step_count=self._nt,
                    original_growth_param=self._barrier3d[iB3D].growthparam,
                )
            )
            self._nourishments.append(
                BeachNourisher(
                    nourishment_interval=self._nourishment_interval[iB3D],
                    nourishment_volume=self._nourishment_volume[iB3D],
                    initial_beach_width=int(
                        self._barrier3d[iB3D].BermEl / self._barrier3d[iB3D]._beta
                    )
                    * 10,  # m
                    dune_design_elevation=self._artificial_max_dune_ele[iB3D],
                    time_step_count=self._nt,
                    original_growth_param=self._barrier3d[iB3D].growthparam,
                )
            )

    @property
    def road_break(self):
        return self._road_break

    @property
    def b3d_break(self):
        return self._b3d_break

    @property
    def barrier3d(self):
        return self._barrier3d

    @barrier3d.setter
    def barrier3d(self, value):
        self._barrier3d = value

    @property
    def roadways(self):
        return self._roadways

    @property
    def chome(self):
        return self._chome_coupler.chome

    @property
    def brie(self):
        return self._brie_coupler.brie

    @property
    def roadway_management_module(self):
        return self._roadway_management_module

    @roadway_management_module.setter
    def roadway_management_module(self, value):
        self._roadway_management_module = value

    @property
    def beach_nourishment_module(self):
        return self._beach_nourishment_module

    @beach_nourishment_module.setter
    def beach_nourishment_module(self, value):
        self._beach_nourishment_module = value

    @property
    def nourishments(self):
        return self._nourishments

    @property
    def nourish_now(self):
        return self._nourish_now

    @nourish_now.setter
    def nourish_now(self, value):
        self._nourish_now = value

    @property
    def rebuild_dune_now(self):
        return self._rebuild_dune_now

    @rebuild_dune_now.setter
    def rebuild_dune_now(self, value):
        self._rebuild_dune_now = value

    @property
    def nourishment_interval(self):
        return self._nourishment_interval

    @nourishment_interval.setter
    def nourishment_interval(self, value):
        self._nourishment_interval = value

    @property
    def nourishment_volume(self):
        return self._nourishment_volume

    @nourishment_volume.setter
    def nourishment_volume(self, value):
        self._nourishment_volume = value

    ###############################################################################
    # time loop
    ###############################################################################

    def update(self):
        """Update Cascade by a single time step"""

        # Check for drowning here from the last time step in brie. Note that this will stay false if brie is not used
        # for AST (i.e., a B3D only run).
        if self._brie_coupler._brie.drown == True:
            return

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
            self._brie_coupler.update_ast(self._barrier3d, x_t_dt, x_s_dt, h_b_dt)
        else:
            for iB3D in range(self._ny):
                self._barrier3d[iB3D].update_dune_domain()

        # check also for width/height drowning in B3D (would occur in update_dune_domain); important for B3D only runs
        for iB3D in range(self._ny):
            if self._barrier3d[iB3D].drown_break == 1:
                self._b3d_break = 1
                return

        # human dynamics modules
        # RoadwayManager: remove overwash from roadway after each model year, place on the dune, rebuild dunes if
        # fall below height threshold, and check if dunes should grow naturally
        if self._roadway_management_module:
            for iB3D in range(self._ny):
                # call roadway management module; if the user changed any of the roadway or dune parameters,
                # they are updated in the update function
                self._roadways[iB3D].update(
                    barrier3d=self._barrier3d[iB3D],
                    road_ele=self._road_ele[iB3D],
                    road_width=self._road_width[iB3D],
                    road_relocation_setback=self._road_relocation_setback[iB3D],
                    dune_design_elevation=self._artificial_max_dune_ele[iB3D],
                    dune_minimum_elevation=self._artificial_min_dune_ele[iB3D],
                )

                # if the road drowned or barrier was too narrow to be relocated, break
                if (
                    self._roadways[iB3D].drown_break
                    or self._roadways[iB3D].narrow_break
                ):
                    self._road_break = 1
                    return

        if self._community_dynamics_module:
            [
                self._nourish_now,
                self._rebuild_dune_now,
                self._nourishment_volume,
            ] = self._chome_coupler.update(
                barrier3d=self._barrier3d,
                nourishments=self._nourishments,
                dune_design_elevation=self._artificial_max_dune_ele,
            )

        # BeachNourisher: if interval specified, nourish at that interval, otherwise wait until told with nourish_now to
        # nourish or rebuild_dunes_now to rebuild dunes. Resets ..._now parameters to zero (false) after nourishment.
        if self._beach_nourishment_module:
            for iB3D in range(self._ny):
                [
                    self._nourish_now[iB3D],
                    self._rebuild_dune_now[iB3D],
                ] = self._nourishments[iB3D].update(
                    barrier3d=self._barrier3d[iB3D],
                    dune_design_elevation=self._artificial_max_dune_ele[iB3D],
                    nourish_now=self._nourish_now[iB3D],
                    rebuild_dune_now=self._rebuild_dune_now[iB3D],
                    nourishment_interval=self._nourishment_interval[iB3D],
                    nourishment_volume=self._nourishment_volume[iB3D],
                )

    ###############################################################################
    # save data
    ###############################################################################

    def save(self, directory):

        filename = self._filename + ".npz"

        csc8d = []
        csc8d.append(self)

        os.chdir(directory)
        np.savez(filename, cascade=csc8d)


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
