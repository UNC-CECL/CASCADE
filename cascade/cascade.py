from pathlib import Path

import numpy
from joblib import Parallel, delayed
import numpy as np
import os
import math

from .roadway_manager import RoadwayManager, set_growth_parameters
from .beach_dune_manager import BeachDuneManager
from .brie_coupler import BrieCoupler, initialize_equal, batchB3D
from .chom_coupler import ChomCoupler
from bmftc import Bmftc


class CascadeError(Exception):
    pass


class Cascade:
    def module_lists(
        self,
        dune_design_elevation,
        dune_minimum_elevation,
        road_width,
        road_ele,
        road_setback,
        nourishment_interval,
        nourishment_volume,
        overwash_filter,
        overwash_to_dune,
        roadway_management_module,
        beach_nourishment_module,
    ):
        """Configures lists to account for multiple barrier3d domains from single input variables; used in modules"""

        if np.size(dune_design_elevation) > 1:
            self._dune_design_elevation = dune_design_elevation  # list of floats option
        else:
            self._dune_design_elevation = [dune_design_elevation] * self._ny
        if np.size(dune_minimum_elevation) > 1:
            self._dune_minimum_elevation = dune_minimum_elevation
        else:
            self._dune_minimum_elevation = [dune_minimum_elevation] * self._ny
        if np.size(road_width) > 1:
            self._road_width = road_width
        else:
            self._road_width = [road_width] * self._ny
        if np.size(road_ele) > 1:
            self._road_ele = road_ele
        else:
            self._road_ele = [road_ele] * self._ny
        if np.size(road_setback) > 1:
            self._road_setback = road_setback
        else:
            self._road_setback = [road_setback] * self._ny
        if np.size(nourishment_interval) > 1:
            self._nourishment_interval = nourishment_interval
        else:
            self._nourishment_interval = [nourishment_interval] * self._ny
        if np.size(nourishment_volume) > 1:
            self._nourishment_volume = nourishment_volume
        else:
            self._nourishment_volume = [nourishment_volume] * self._ny
        if np.size(overwash_filter) > 1:
            self._overwash_filter = overwash_filter
        else:
            self._overwash_filter = [overwash_filter] * self._ny
        if np.size(overwash_to_dune) > 1:
            self._overwash_to_dune = overwash_to_dune
        else:
            self._overwash_to_dune = [overwash_to_dune] * self._ny
        if np.size(roadway_management_module) > 1:
            self._roadway_management_module = roadway_management_module
        else:
            self._roadway_management_module = [roadway_management_module] * self._ny
        if np.size(beach_nourishment_module) > 1:
            self._beach_nourishment_module = beach_nourishment_module
        else:
            self._beach_nourishment_module = [beach_nourishment_module] * self._ny

        return

    def reset_dune_growth_rates(
        self,
        original_growth_param,
        iB3D,
    ):
        """Reset dune growth paramerters to original values after human abandonment; dune heights must drop below Dmax
        before reset"""

        time_index = self._barrier3d[iB3D].time_index

        new_growth_parameters = set_growth_parameters(
            self._barrier3d[iB3D].DuneDomain[time_index - 1, :, :],
            self._barrier3d[iB3D].Dmax,
            self._barrier3d[iB3D].growthparam,
            original_growth_param=original_growth_param,  # use original growth rates for resetting values
        )

        return new_growth_parameters

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
        bay_depth=3.0,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=6,
        time_step_count=200,
        min_dune_growth_rate=0.25,  # average is 0.45, a low dune growth rate
        max_dune_growth_rate=0.65,
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=True,
        beach_nourishment_module=True,
        community_dynamics_module=False,
        marsh_dynamics = False,
        enable_shoreline_offset=False,
        shoreline_offset=[],
        road_ele=1.7,  # ---------- the rest of these variables are for the human dynamics modules --------------- #
        road_width=30,
        road_setback=30,
        dune_design_elevation=3.7,
        dune_minimum_elevation=2.2,
        nourishment_interval=None,
        nourishment_volume=300.0,
        overwash_filter=40,
        overwash_to_dune=10,
        number_of_communities=1,
        sand_cost=10,
        taxratio_oceanfront=1,
        external_housing_market_value_oceanfront=6e5,
        external_housing_market_value_nonoceanfront=4e5,
        fixed_cost_beach_nourishment=2e6,
        fixed_cost_dune_nourishment=2e5,
        nourishment_cost_subsidy=10e6,
        house_footprint_x=15,
        house_footprint_y=20,
        beach_full_cross_shore=70,
    ):
        """

        CASCADE: The CoAStal Community-lAnDscape Evolution model

        Couples Barrier3D (Reeves et al., 2019), the Barrier Inlet Model (BRIE; Nienhuis and Lorenzo Trueba, 2019), &
        C-HOM (Williams et al., in prep)

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
        bay_depth: float, optional
            Depth of back-barrier bay [m]
        sea_level_rise_rate: float, optional
            Rate of sea_level rise [m/yr].
        sea_level_rise_constant: boolean, optional
            Rate of sea_level rise is constant if True, otherwise a logistic growth function is used.
        background_erosion: float,
            Rate of shoreline retreat attributed to gradients in alongshore transport; (-) = erosion, (+) = acc [m / y]
        alongshore_section_count: int, optional
            Number of alongshore sections.
        time_step_count: int, optional
            Number of time steps.
        min_dune_growth_rate: float or list of floats, optional
            Minimum dune growth rate [unitless, for Houser growth rate formulation]
        max_dune_growth_rate: float or list of floats, optional
            Maximum dune growth rate [unitless, for Houser growth rate formulation]
        num_cores: int, optional
            Number of (parallel) processing cores to be used
        roadway_management_module: boolean or list of booleans, optional
            If True, use roadway management module (overwash removal, road relocation, dune management)
        alongshore_transport_module: boolean, optional
            If True, couple Barrier3D with BRIE to use diffusive model for AST
        community_dynamics_module: boolean, optional
            If True, couple with CHOM, a community decision making model; requires nourishment module
        beach_nourishment_module: boolean or list of booleans, optional
            If True, use nourishment module (nourish shoreface, rebuild dunes)
        marsh_dynamics: boolean, optional
            If True, use the PyBMFT module to represent marsh and bay dynamics
        enable_shoreline_offset: bool, optional
            State whether you want a shoreline offset [True / False]
        shoreline_offset: list, optional
            The alongshore offset between different Barrier3d sections [m]
        road_ele: float or list of floats, optional
            Elevation of the initial roadway [m MHW] and after road relocations
        road_width: int or list of int, optional
            Width of initial roadway [m] and road relocations
        road_setback: int or list of int, optional
            Setback of initial roadway from the inital dune line [m] and after road relocations
        dune_design_elevation: float or list of floats, optional
            Elevation to which dune is initially rebuilt to [m MHW] and after road relocations
        dune_minimum_elevation: float or list of floats, optional
            Elevation threshold which triggers rebuilding of dune [m MHW]
        nourishment_interval: int or list of ints, optional
             Interval that nourishment occurs [yrs]
        nourishment_volume: float or list of float, optional
             Volume of nourished sand along cross-shore transect [m^3/m]
        overwash_filter: float or list of floats,
            Percent overwash removed from barrier interior [40-90% (residential-->commercial) from Rogers et al., 2015]
        overwash_to_dune: float or list of floats,
            Percent overwash removed from barrier interior to dunes [%, overwash_filter+overwash_to_dune <=100]
        number_of_communities: int, optional
            Number of communities (CHOM model instances) described by the alongshore section count (Barrier3D grids)
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
        beach_full_cross_shore: int, optional
            The cross-shore extent (meters) of fully nourished beach (i.e., the community desired beach width) [m]


        Examples
        --------
        >>> from cascade.cascade import Cascade
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
        self._background_erosion = background_erosion
        self._num_cores = num_cores
        self._alongshore_transport_module = alongshore_transport_module
        self._community_dynamics_module = community_dynamics_module
        self._filename = name
        self._storm_file = storm_file
        self._elevation_file = elevation_file
        self._dune_file = dune_file
        self._parameter_file = parameter_file
        self._number_of_communities = number_of_communities
        self._b3d_break = 0  # true if barrier in barrier3d height or width drowns -- if this happens, entire sim stops
        self._road_break = [
            0
        ] * self._ny  # true if roadway drowns from bay reaching roadway
        self._community_break = [
            0
        ] * self._ny  # true if community breaks due to minimum barrier width
        self._nourish_now = [0] * self._ny  # triggers nourishment
        self._rebuild_dune_now = [0] * self._ny  # triggers dune rebuilding
        self._initial_beach_width = [0] * self._ny
        # New offset shoreline additions
        self._enable_shoreline_offset = enable_shoreline_offset
        self._shoreline_offset = shoreline_offset
        self._marsh_dynamics = marsh_dynamics

        ###############################################################################
        # initialize brie and barrier3d model classes
        ###############################################################################

        # initialize brie: used for initial shoreface calculations, AST (optional), tidal inlets (eventually)
        self._brie_coupler = BrieCoupler(
            name=name,
            wave_height=self._wave_height,
            wave_period=self._wave_period,
            wave_asymmetry=self._wave_asymmetry,
            wave_angle_high_fraction=self._wave_angle_high_fraction,
            sea_level_rise_rate=self._sea_level_rise_rate,
            back_barrier_depth=bay_depth,
            ny=self._ny,
            nt=self._nt,
        )

        # Create offset shorelines in BRIE
        self._brie_coupler.offset_shoreline(
            enable_shoreline_offset=self._enable_shoreline_offset,
            offset_values=self._shoreline_offset,
            ny=self._ny,
        )

        # initialize barrier3d models (number set by brie ny above) and make both brie and barrier3d classes equivalent
        self._barrier3d = initialize_equal(
            datadir=datadir,
            brie=self._brie_coupler._brie,
            slr_constant=self._slr_constant,
            rmin=self._rmin,  # can be array
            rmax=self._rmax,  # can be array
            background_erosion=self._background_erosion,  # can be array
            parameter_file=self._parameter_file,
            storm_file=self._storm_file,
            dune_file=self._dune_file,  # can be array
            elevation_file=self._elevation_file,  # can be array
        )

        ###############################################################################
        # initialize marsh dynamic modules
        ###############################################################################
        if self._marsh_dynamics:
            self._bmftc = []
            self._BMFTC_Break = False # Initialize BMFTC break variable
            # Initialize blank PyBMFT variables as lists
            self._name = []
            self._x_b_TS = []
            self._LandscapeTypeWidth_TS = []
            self._bay_overwash_carryover = []  # [m^3] Volume of overwash deposition into back-barrier bay from previous year that did not fill new cell up to sea level; is added to overwash bay dep in following year
            self._x_s_offset = []  # Initial location of B in PyBMFT-C relative to x_s_initial in Barrier3D
            self._cumul_len_change = []
            self._delta_fetch_TS = []
            self._OWspread = 0  # [%] Percentage of overwash past marsh edge that is spread across bay

            # initialize PyBMFT models (number set by brie ny above)
            for iB3D in range(self._ny):
                self._bmftc.append(
                    Bmftc(
                        name="back-barrier",
                        time_step_count=time_step_count,
                        relative_sea_level_rise=self._barrier3d[iB3D]._RSLR[1]*1000,
                        reference_concentration=60,
                        slope_upland=0.005,
                        bay_fetch_initial=5000,
                        forest_width_initial_fixed=False,
                        forest_width_initial=5000,  # 5000 accomodates 250 yrs at R=15 and S=0.001
                        wind_speed=6,
                        forest_on=False,
                        filename_equilbaydepth="/Users/ceclmac/PycharmProjects/PyBMFT-C/Input/PyBMFT-C/Equilibrium Bay Depth.mat",
                        filename_marshspinup="/Users/ceclmac/PycharmProjects/PyBMFT-C/Input/PyBMFT-C/MarshStrat_all_RSLR1_CO50_width500.mat",
                        marsh_width_initial=500,
                    )
                )

        # Equalize Barrier3D/PyBMFT-C Values of Identical Parameters
            #for iB3D in range(self._ny):
                #self._bmftc[iB3D]._dur =(self._barrier3d[iB3D]._TMAX - 1)
                #self._bmftc[iB3D]._RSLRi = (self._barrier3d[iB3D]._RSLR[1]*1000)
                #self._barrier3d[iB3D]._BayDepth = self._bmftc[iB3D].Bay_depth[self._bmftc[iB3D].startyear - 1] / 10
                #BRIE._back_barrier_depth = _barrier3d._model._BayDepth

            for iB3D in range(self._ny):
                # ===========================================
                # Add initial barrier topography from Barrier3D to initial "forest" (i.e., subaerial) portion of PyBMFT-C transect
                b3d_transect = np.mean(self._barrier3d[iB3D].InteriorDomain,axis=1) * 10  # Take average across alongshore dimension, convert to m (vertical dimension)
                x = np.linspace(1, len(b3d_transect) * 10, num=len(b3d_transect) * 10)
                xp = np.linspace(1, len(b3d_transect), num=len(b3d_transect)) * 10
                xp = xp - 5
                b3d_transect = np.interp(x, xp, b3d_transect)  # Interpolate from dam to m (horizontal dimension)
                x_f = np.where(b3d_transect < (self._barrier3d[iB3D].SL * 10))[0][0] - 1  # [m] Distance of first interior (subaerial) cell from B3D ocean shoreline (excluding dunes/beach)
                b3d_transect = b3d_transect[:x_f]
                b3d_transect = np.flip(b3d_transect)
                b3d_transect = b3d_transect + self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear - 1] + self._bmftc[iB3D].amp + (self._bmftc[iB3D].RSLRi / 1000)  # Convert vertical datums

                # Adjust size of Barrier3D topo to fit PyBMFT-C "forest" section
                BB_forest_len = len(self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear, self._bmftc[iB3D].x_f:])
                if len(b3d_transect) > BB_forest_len:
                    subtract = len(b3d_transect) - BB_forest_len
                    b3d_transect = b3d_transect[:-subtract]
                elif len(b3d_transect) < BB_forest_len:
                    add = np.ones([BB_forest_len - len(b3d_transect)]) * (
                                self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear] + self._bmftc[iB3D].amp)
                    b3d_transect = np.append(b3d_transect, add)

                # Replace initial subaerial elevation in PyBMFT-C with Barrier3D initial barrier elevation
                self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear - 1,
                self._bmftc[iB3D].x_f:] = b3d_transect  # Replace!


            # ===========================================
            # Populate blank PyBMFT list variables
                self._name.append(name)
                self._x_b_TS.append(np.zeros([self._bmftc[iB3D].dur]))
                self._LandscapeTypeWidth_TS.append(np.zeros([self._bmftc[iB3D].dur, 4]))
                self._bay_overwash_carryover.append(0)  # [m^3] Volume of overwash deposition into back-barrier bay from previous year that did not fill new cell up to sea level; is added to overwash bay dep in following year
                initial_subaerial_width =self._bmftc[iB3D].B - self._bmftc[iB3D].x_f
                self._x_s_offset.append(initial_subaerial_width - (self._barrier3d[iB3D].InteriorWidth_AvgTS[-1] * 10))  # Initial location of B in PyBMFT-C relative to x_s_initial in Barrier3D
                self._cumul_len_change.append([0])
                #self._OWspread.append(0)  # [%] Percentage of overwash past marsh edge that is spread across bay
                self._delta_fetch_TS.append([])
        #
        ###############################################################################
        # initialize human dynamics modules
        ###############################################################################

        # configure `self` to create lists of these variables; time series of these variables are saved in modules
        self.module_lists(
            dune_design_elevation=dune_design_elevation,
            dune_minimum_elevation=dune_minimum_elevation,
            road_width=road_width,
            road_ele=road_ele,
            road_setback=road_setback,
            nourishment_interval=nourishment_interval,
            nourishment_volume=nourishment_volume,
            overwash_filter=overwash_filter,
            overwash_to_dune=overwash_to_dune,
            roadway_management_module=roadway_management_module,
            beach_nourishment_module=beach_nourishment_module,
        )

        if self._community_dynamics_module:
            if not any(self._beach_nourishment_module):
                CascadeError(
                    "Beach nourishment module must be set to `TRUE` to couple with CHOM"
                )
            else:
                self._chom_coupler = ChomCoupler(
                    barrier3d=self._barrier3d,
                    total_time=self._nt,
                    alongshore_length_b3d=self._brie_coupler._brie._dy,  # this is the barrier3d default, 500 m
                    dune_design_elevation=self._dune_design_elevation,
                    number_of_communities=self._number_of_communities,
                    name=self._filename,
                    sand_cost=sand_cost,
                    taxratio_oceanfront=taxratio_oceanfront,
                    external_housing_market_value_oceanfront=external_housing_market_value_oceanfront,
                    external_housing_market_value_nonoceanfront=external_housing_market_value_nonoceanfront,
                    fixed_cost_beach_nourishment=fixed_cost_beach_nourishment,
                    fixed_cost_dune_nourishment=fixed_cost_dune_nourishment,
                    nourishment_cost_subsidy=nourishment_cost_subsidy,
                    house_footprint_x=house_footprint_x,
                    house_footprint_y=house_footprint_y,
                    beach_full_cross_shore=beach_full_cross_shore,
                )  # contains the CHOM model instances, one per community

        # initialize RoadwayManager and BeachDuneManager modules
        # (always, just in case we want to add a road or start nourishing during the simulation)
        self._roadways = []
        self._nourishments = []
        self._bay_side_breaches = []

        for iB3D in range(self._ny):
            self._roadways.append(
                RoadwayManager(
                    initial_road_elevation=self._road_ele[iB3D],
                    road_width=self._road_width[iB3D],
                    road_setback=self._road_setback[iB3D],
                    initial_dune_design_elevation=self._dune_design_elevation[iB3D],
                    initial_dune_minimum_elevation=self._dune_minimum_elevation[iB3D],
                    time_step_count=self._nt,
                    original_growth_param=self._barrier3d[iB3D].growthparam,
                )
            )

            self._initial_beach_width[iB3D] = (
                int(self._barrier3d[iB3D].BermEl / self._barrier3d[iB3D]._beta) * 10
            )
            self._nourishments.append(
                BeachDuneManager(
                    nourishment_interval=self._nourishment_interval[iB3D],
                    nourishment_volume=self._nourishment_volume[iB3D],
                    initial_beach_width=self._initial_beach_width[iB3D],
                    dune_design_elevation=self._dune_design_elevation[iB3D],
                    time_step_count=self._nt,
                    original_growth_param=self._barrier3d[iB3D].growthparam,
                    overwash_filter=self._overwash_filter[iB3D],
                    overwash_to_dune=self._overwash_to_dune[iB3D],
                )
            )

        # use the initial beach width as a check on the barrier3d user input for mulitple domains; the beach width
        # must be the same for all domains because there is only one storm file, which is made for a set berm
        # elevation and beach slope
        if all(
            elem == self._initial_beach_width[0] for elem in self._initial_beach_width
        ):
            pass
        else:
            CascadeError(
                "Berm elevation and beach slope must be equivalent for all Barrier3D domains"
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
    def chom(self):
        return self._chom_coupler.chom

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

    @property
    def community_break(self):
        return self._community_break

    ###############################################################################
    # time loop
    ###############################################################################

    def update(self, Time_step):
        self._time_step = Time_step
        if self._marsh_dynamics ==False:
            """Update Cascade by a single time step"""
            self._time_step = Time_step

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

            # use brie to connect B3D subgrids with alongshore sediment transport; otherwise, just update (erode/prograde)
            # dune domain
            if self._alongshore_transport_module:
                self._brie_coupler.update_ast(
                    self._barrier3d, x_t_dt, x_s_dt, h_b_dt
                )  # also updates dune domain
            else:
                for iB3D in range(self._ny):
                    self._barrier3d[iB3D].update_dune_domain()

            # check also for width/height drowning in B3D (would occur in update_dune_domain)
            for iB3D in range(self._ny):
                if self._barrier3d[iB3D].drown_break == 1:
                    self._b3d_break = 1
                    return

        ###############################################################################
        # Backbarrier marsh module
        ###############################################################################
        # ~~~~~~~~~~~~~~ PyBMFT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Run PyBMFT module to represent marsh growth and erosion from the back-bay
        if self._marsh_dynamics:
            for iB3D in range(self._ny):
                print('Run PyBMFT ' +str(iB3D))
                """Update BarrierBMFT by one time step"""

                # ===================================================================================================================================================================================================================================
                # ===================================================================================================================================================================================================================================
                # Advance PyBMFT-C back-barrier marshes
                self._bmftc[iB3D].update()

                # Check if marsh has completely drowned or basin is completely full
                if self._bmftc[iB3D].drown_break == 1:
                    self._bmftc[iB3D]._dur = self._time_step
                    self._bmftc[iB3D]._endyear = self._bmftc[iB3D].startyear + self._time_step
                    self._BMFTC_Break = True
                    print("PyBMFT-C Simulation Break: marsh has completely drowned or basin is completely full")
                    return  # If so, end simulation

                # ===================================================================================================================================================================================================================================
                # ===================================================================================================================================================================================================================================
                # Update fetch and marsh point locations from PyBMFT-C bay erosion/deposition processes

                # Calculate change in fetch from erosion of both marshes
                delta_fetch_BB = self._bmftc[iB3D].bfo - self._bmftc[iB3D].fetch[
                    self._bmftc[iB3D].startyear + self._time_step - 1]  # [m] Back-barrier marsh

                self._delta_fetch_TS[iB3D].append(delta_fetch_BB)

                # Determine change in x_b location
                self._x_b_TS[iB3D][self._time_step] = self._bmftc[iB3D].x_b  # Save to array

                # ===================================================================================================================================================================================================================================
                # ===================================================================================================================================================================================================================================
                # Adjust bay depth in Barrier3D according to depth calculated in PyBMFT-C
                self._barrier3d[iB3D]._BayDepth = np.mean([self._bmftc[iB3D].db]) / 10

                # ===================================================================================================================================================================================================================================
                # ===================================================================================================================================================================================================================================
                # Add marsh from PyBMFT-C to Barrier3D

                # Extract and convert marsh elevation from PyBMFT-C
                marsh_transect = self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                                 self._bmftc[iB3D].x_m: self._bmftc[iB3D].x_f + 1]  # Marsh elevation from PyBMFT-C
                if len(marsh_transect) >= 1:
                    len_marsh_transect = 10 * (
                                (len(marsh_transect) + 5) // 10)  # Cross-shore length of marsh rounded to nearest dam
                    self._cumul_len_change[iB3D].append(
                        self._cumul_len_change[iB3D][-1] + (len_marsh_transect - len(marsh_transect)))
                    x = np.linspace(1, len(marsh_transect) / 10, num=int((len_marsh_transect / 10)))
                    xp = np.linspace(1, len(marsh_transect) / 10, num=int(len(marsh_transect)))
                    marsh_transect = np.interp(x, xp,marsh_transect)  # Interpolate marsh elevation from m to dam in the horizontal dimension
                    marsh_transect = marsh_transect - (self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + self._time_step - 1] + self._bmftc[iB3D].amp)  # Make marsh elevation relative to MHW datum
                    marsh_transect = marsh_transect / 10  # Convert from m to dam in the vertial dimension
                    marsh_transect = np.flip(marsh_transect)
                StartDomainWidth = np.shape(self._barrier3d[iB3D].InteriorDomain)[0]  # Width of interior domain from last time step

                # Find barrier interior widths for each dam alongshore
                InteriorWidth = [0] * self._barrier3d[iB3D].BarrierLength
                for bl in range(self._barrier3d[iB3D].BarrierLength):
                    width = next((index for index, value in enumerate(self._barrier3d[iB3D].InteriorDomain[:, bl]) if
                                  value <= self._barrier3d[iB3D].SL), StartDomainWidth)
                    width = width - 1
                    if width < 0:
                        width = 0
                    InteriorWidth[bl] = width

                # Update Barrier3D Domain Sizes
                Target_width_barriermarsh = self._bmftc[iB3D].B - self._bmftc[iB3D].x_m - self._x_s_offset[iB3D]  # [m] Target width of barrier-marsh
                Target_width_barriermarsh = math.ceil(Target_width_barriermarsh / 10)  # [dam]
                addRows = Target_width_barriermarsh - StartDomainWidth + 1  # Number of rows to add (if positive) or subtract (if negative) from Barrier3D domain

                if addRows > 0:
                    # Update interior domain size
                    Marsh_Addition = np.ones([addRows, self._barrier3d[iB3D].BarrierLength]) * -self._barrier3d[
                        iB3D]._BayDepth
                    Zero_Addition = np.zeros([addRows, self._barrier3d[iB3D].BarrierLength])
                    NewDomain = np.vstack([self._barrier3d[iB3D].InteriorDomain, Marsh_Addition])
                    # Update size of shrub domains, too
                    self._barrier3d[iB3D]._ShrubDomainFemale = np.vstack(
                        [self._barrier3d[iB3D]._ShrubDomainFemale, Zero_Addition])
                    self._barrier3d[iB3D]._ShrubDomainMale = np.vstack(
                        [self._barrier3d[iB3D]._ShrubDomainMale, Zero_Addition])
                    self._barrier3d[iB3D]._ShrubDomainDead = np.vstack(
                        [self._barrier3d[iB3D]._ShrubDomainDead, Zero_Addition])
                    self._barrier3d[iB3D]._ShrubPercentCover = np.vstack(
                        [self._barrier3d[iB3D]._ShrubPercentCover, Zero_Addition])
                    self._barrier3d[iB3D]._DeadPercentCover = np.vstack(
                        [self._barrier3d[iB3D]._DeadPercentCover, Zero_Addition])
                    self._barrier3d[iB3D]._BurialDomain = np.vstack(
                        [self._barrier3d[iB3D]._BurialDomain, Zero_Addition])
                    self._barrier3d[iB3D]._ShrubDomainAll = self._barrier3d[iB3D]._ShrubDomainFemale + self._barrier3d[
                        iB3D]._ShrubDomainMale
                elif addRows < 0:
                    # Update interior domain size
                    NewDomain = self._barrier3d[iB3D].InteriorDomain[:addRows, :]
                    # Update size of shrub domains, too
                    self._barrier3d[iB3D]._ShrubDomainFemale = self._barrier3d[iB3D]._ShrubDomainFemale[:addRows, :]
                    self._barrier3d[iB3D]._ShrubDomainMale = self._barrier3d[iB3D]._ShrubDomainMale[:addRows, :]
                    self._barrier3d[iB3D]._ShrubDomainDead = self._barrier3d[iB3D]._ShrubDomainDead[:addRows, :]
                    self._barrier3d[iB3D]._ShrubPercentCover = self._barrier3d[iB3D]._ShrubPercentCover[:addRows, :]
                    self._barrier3d[iB3D]._DeadPercentCover = self._barrier3d[iB3D]._DeadPercentCover[:addRows, :]
                    self._barrier3d[iB3D]._BurialDomain = self._barrier3d[iB3D]._BurialDomain[:addRows, :]
                    self._barrier3d[iB3D]._ShrubDomainAll = self._barrier3d[iB3D]._ShrubDomainFemale + self._barrier3d[
                        iB3D]._ShrubDomainMale
                else:
                    NewDomain = self._barrier3d[iB3D].InteriorDomain  # Domains stay same size

                if len(marsh_transect) >= 1:
                    # Update Marsh In Barrier3D
                    x_marsh = Target_width_barriermarsh + 1  # [dam] Cross-shore location of marsh edge relative to interior domain
                    for w in range(self._barrier3d[iB3D].BarrierLength):
                        width_diff = x_marsh - (InteriorWidth[w] + len(marsh_transect))
                        if width_diff < 0:
                            MarshTransect = marsh_transect[:-int(abs(width_diff))]  # [dam]
                        elif width_diff > 0:
                            add = np.ones([int(abs(width_diff))]) * marsh_transect[
                                -1]  # Set additional marsh cells to elevation of last marsh
                            MarshTransect = np.append(marsh_transect, add)  # [dam]
                        else:
                            MarshTransect = marsh_transect  # [dam]

                        InteriorTransect = NewDomain[:InteriorWidth[w], w]  # [dam]
                        BarrierMarshTransect = np.append(InteriorTransect, MarshTransect)  # Combine interior and marsh

                        NewDomain[:len(BarrierMarshTransect), w] = BarrierMarshTransect
                        NewDomain[len(BarrierMarshTransect):, w] = (self._barrier3d[iB3D].SL - np.mean(
                            [self._bmftc[iB3D].db])) / 10

                self._barrier3d[iB3D].InteriorDomain = NewDomain
                # ===================================================================================================================================================================================================================================
                # ===================================================================================================================================================================================================================================
                # Advance Barrier3D
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

                # use brie to connect B3D subgrids with alongshore sediment transport; otherwise, just update (erode/prograde)
                # dune domain
                if self._alongshore_transport_module:
                    self._brie_coupler.update_ast(
                        self._barrier3d, x_t_dt, x_s_dt, h_b_dt
                    )  # also updates dune domain
                else:
                    for iB3D in range(self._ny):
                        self._barrier3d[iB3D].update_dune_domain()

                # check also for width/height drowning in B3D (would occur in update_dune_domain)
                for iB3D in range(self._ny):
                    if self._barrier3d[iB3D].drown_break == 1:
                        self._b3d_break = 1
                        return

                # ===================================================================================================================================================================================================================================
                # ===================================================================================================================================================================================================================================
                # Update PyBMFT-C transect elevation based on Barrier3D elevation change

                shoreline_change = self._barrier3d[iB3D].x_s_TS[-1] - self._barrier3d[iB3D].x_s_TS[-2]
                self._x_s_offset[iB3D] = self._x_s_offset[iB3D] + (shoreline_change * 10)
                start_b3d = np.mean(NewDomain,axis=1) * 10  # Barrier3D domain before update, averaged across alongshore dimension, converted to m (vertical dimension)
                end_b3d = np.mean(self._barrier3d[iB3D].InteriorDomain,axis=1) * 10  # Barrier3D domain after update, averaged across alongshore dimension, converted to m (vertical dimension)

                # Update start domain size to match end domain
                sc_b3d = self._barrier3d[iB3D].ShorelineChangeTS[1]  # Shoreline change [dam] from Barrier3D model update (this timestep)
                if sc_b3d < 0:  # Shoreline erosion
                    start_b3d = start_b3d[abs(sc_b3d):]  # Trim off front
                elif sc_b3d > 0:  # Shoreline progradation
                    add = np.zeros([sc_b3d])
                    start_b3d = np.append(add, start_b3d)  # Add zeros to front

                if len(start_b3d) < len(end_b3d):
                    add = np.ones([len(end_b3d) - len(start_b3d)]) * np.mean([self._bmftc[iB3D].db]) * -1  # Add bay cells
                    start_b3d = np.append(start_b3d, add)
                elif len(start_b3d) > len(end_b3d):
                    subtract = len(end_b3d) - len(start_b3d)
                    start_b3d = start_b3d[:subtract]

                # Calculate change in elevation from Barrier3D update
                end_b3d = end_b3d + (self._barrier3d[iB3D].RSLR[self._time_step] * 10)  # Offset sea-level rise from Barrier3D so that it isn't counted twice (i.e. RSLR already taken into account in PyBMFT-C)
                elevation_change_b3d = end_b3d - start_b3d  # Change in elevation across transect after Barrier3d update; [dam] horizontal dimension, [m] vertical dimentsion

                # Interpolate from dam to m (horizontal dimension)
                x = np.linspace(1, len(elevation_change_b3d) * 10, num=len(elevation_change_b3d) * 10)
                xp = np.linspace(1, len(elevation_change_b3d), num=len(elevation_change_b3d)) * 10
                xp = xp - 5
                elevation_change_b3d = np.interp(x, xp, elevation_change_b3d)
                off = int(abs(math.floor(self._x_s_offset[iB3D])))  # [m] Offset of barrier shoreline and B
                # Incorporate elevation change from Barrier3D into back-barrier instance of PyBMFT-C
                if int(math.floor(self._x_s_offset[iB3D])) < 0:
                    elevation_change_b3d = np.flip(elevation_change_b3d[off:])  # Flip orientation
                    marsh_barrier_width = (self._bmftc[iB3D].B - self._bmftc[iB3D].x_m)
                    x_m_change = abs(math.floor(len(elevation_change_b3d) - marsh_barrier_width))  # Location of marsh edge within elevation_change_b3d

                    # Add subaerial elevation change
                    self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,-len(elevation_change_b3d[x_m_change:]):] += elevation_change_b3d[x_m_change:] # Store mass of overwash mineral sediment deposited across transect
                    self._bmftc[iB3D].mineral_dep[self._bmftc[iB3D].startyear + self._time_step,-len(elevation_change_b3d[x_m_change:]):] += (elevation_change_b3d[x_m_change:] * self._bmftc[iB3D].rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

                    # Determine volume of sed deposited past initial marsh edge and into bay
                    sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (1 - self._OWspread)  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
                    sum_bay_dep = np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
                    self._bmftc[iB3D]._Fow_min = max(0, sum_bay_dep * self._bmftc[iB3D].rhos)  # [kg/yr] Overwash deposition into bay, volume converted to mass

                    # Add volume of carryover from last time step
                    sum_marsh_dep += self._bay_overwash_carryover  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

                    # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
                    new_marsh_height = self._bmftc[iB3D].db

                    # Determine distance of marsh progradation from overwash deposition
                    progradation_actual = sum_marsh_dep / new_marsh_height  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
                    progradation = int(max(math.floor(progradation_actual), 0))  # Round to nearest FULL meter
                    self._bay_overwash_carryover = (progradation_actual - progradation) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

                    if progradation > 0:
                        # Add subaqueous elevation change
                        self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m - progradation: self._bmftc[iB3D].x_m] += new_marsh_height
                        # Store mass of overwash mineral sediment deposited across transect
                        self._bmftc[iB3D].mineral_dep[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m - progradation: self._bmftc[iB3D].x_m] += (
                                new_marsh_height * self._bmftc[
                            iB3D].rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

                    # Spread overwash bay flux evenly across bay bottom
                    bay_accrete = sum_bay_dep / (self._bmftc[
                                                     iB3D].bfo - progradation)  # [m] Vertical accretion of bay bottom from overwash deposition in bay
                    self._bmftc[iB3D]._db = self._bmftc[iB3D].db + bay_accrete  # Update bay depth

                elif int(math.floor(self._x_s_offset[iB3D])) > 0:
                    elevation_change_b3d = np.flip(elevation_change_b3d)
                    marsh_barrier_width = (self._bmftc[iB3D].B - self._bmftc[iB3D].x_m)
                    x_m_change = abs(math.floor(len(elevation_change_b3d) - (marsh_barrier_width - off)))  # Location of marsh edge within elevation_change_b3d

                    # Add subaerial elevation change
                    self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                    self._bmftc[iB3D].B - off - len(elevation_change_b3d[x_m_change:]): self._bmftc[iB3D].B - off] += elevation_change_b3d[x_m_change:]
                    # Store mass of overwash mineral sediment deposited across transect
                    self._bmftc[iB3D].mineral_dep[self._bmftc[iB3D].startyear + self._time_step,
                    self._bmftc[iB3D].B - off - len(elevation_change_b3d[x_m_change:]): self._bmftc[iB3D].B - off] += (elevation_change_b3d[x_m_change:] * self._bmftc[iB3D].rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

                    # Determine volume of sed deposited past initial marsh edge and into bay
                    sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (1 - self._OWspread)  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
                    sum_bay_dep = np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
                    self._bmftc[iB3D]._Fow_min = max(0, sum_bay_dep * self._bmftc[iB3D].rhos)  # [kg/yr] Overwash deposition into bay, volume converted to mass

                    # Add volume of carryover from last time step
                    sum_marsh_dep += self._bay_overwash_carryover  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

                    # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
                    new_marsh_height = self._bmftc[iB3D].db

                    # Determine distance of marsh progradation from overwash deposition
                    progradation_actual = sum_marsh_dep / new_marsh_height  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
                    progradation = int(max(math.floor(progradation_actual), 0))  # Round to nearest FULL meter
                    self._bay_overwash_carryover = ( progradation_actual - progradation) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

                    if progradation > 0:
                        # Add subaqueous elevation change
                        self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m - progradation: self._bmftc[iB3D].x_m] += new_marsh_height
                        # Store mass of overwash mineral sediment deposited across transect
                        self._bmftc[iB3D].mineral_dep[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m - progradation: self._bmftc[iB3D].x_m] += (new_marsh_height * self._bmftc[
                            iB3D].rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

                    # Spread 50% of overwash bay flux evenly across bay bottom
                    bay_accrete = sum_bay_dep / (self._bmftc[iB3D].bfo - progradation)  # [m] Vertical accretion of bay bottom from overwash deposition in bay
                    self._bmftc[iB3D]._db = self._bmftc[iB3D].db + bay_accrete  # Update bay depth

                    # Remove barrier at front and set to msl to account for shoreline change
                    self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step, -off:] = \
                    self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + self._time_step] + self._bmftc[iB3D].amp

                else:
                    elevation_change_b3d = np.flip(elevation_change_b3d)
                    marsh_barrier_width = (self._bmftc[iB3D].B - self._bmftc[iB3D].x_m)
                    x_m_change = abs(math.floor(len(elevation_change_b3d) - marsh_barrier_width))  # Location of marsh edge within elevation_change_b3d

                    # Add subaerial elevation change
                    self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,-len(elevation_change_b3d[x_m_change:]):] += elevation_change_b3d[x_m_change:]
                    # Store mass of overwash mineral sediment deposited across transect
                    self._bmftc[iB3D].mineral_dep[self._bmftc[iB3D].startyear + self._time_step,-len(elevation_change_b3d[x_m_change:]):] += (elevation_change_b3d[x_m_change:] * self._bmftc[iB3D].rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

                    # Determine volume of sed deposited past initial marsh edge and into bay
                    sum_marsh_dep = np.sum(elevation_change_b3d[:x_m_change]) * (1 - self._OWspread)  # [m^3] Volume of overwash deposition landward of marsh edge deposited as marsh
                    sum_bay_dep = np.sum(elevation_change_b3d[:x_m_change]) * self._OWspread  # [m^3] Volume of overwash deposition landward of marsh edge deposited across bay bottom
                    self._bmftc[iB3D]._Fow_min = max(0, sum_bay_dep * self._bmftc[iB3D].rhos)  # [kg/yr] Overwash deposition into bay, volume converted to mass

                    # Add volume of carryover from last time step
                    sum_marsh_dep += self._bay_overwash_carryover  # [m^3] Bay deposition from previous time step that wasn't enough to fully fill bay cell up to sea level

                    # Calculate height of deposition needed to bring bay bottom up to avg marsh elevation
                    new_marsh_height = self._bmftc[iB3D].db

                    # Determine distance of marsh progradation from overwash deposition
                    progradation_actual = sum_marsh_dep / new_marsh_height  # [m] Amount of marsh progradation, in which all overwash dep in bay fills first bay cell, then second, and so on until no more sediment. Assumes overwash is not spread out over bay.
                    progradation = int(max(math.floor(progradation_actual), 0))  # Round to nearest FULL meter
                    self._bay_overwash_carryover = (progradation_actual - progradation) * new_marsh_height  # Save leftover volume of sediment to be added to sum_bay_dep in following time step

                    if progradation > 0:
                        # Add subaqueous elevation change
                        self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m - progradation: self._bmftc[iB3D].x_m] += new_marsh_height
                        # Store mass of overwash mineral sediment deposited across transect
                        self._bmftc[iB3D].mineral_dep[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m - progradation: self._bmftc[iB3D].x_m] += (new_marsh_height * self._bmftc[iB3D].rhos * 1000)  # [g] Mass of pure mineral sediment deposited by overwash

                    # Spread 50% of overwash bay flux evenly across bay bottom
                    bay_accrete = sum_bay_dep / (self._bmftc[iB3D].bfo - progradation)  # [m] Vertical accretion of bay bottom from overwash deposition in bay
                    self._bmftc[iB3D]._db = self._bmftc[iB3D].db + bay_accrete  # Update bay depth

                # Calculate new marsh and "forest" edge positions after overwash
                self._bmftc[iB3D]._x_m = self._bmftc[iB3D].x_m - progradation
                try:
                    self._bmftc[iB3D]._x_f = max(self._bmftc[iB3D].x_m + 1, np.where(
                        self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step, :] >
                        self._bmftc[iB3D].msl[
                            self._bmftc[iB3D].startyear + self._time_step] + self._bmftc[iB3D].amp - self._bmftc[
                            iB3D].Dmin + 0.03)[0][0])
                except IndexError:
                    self._bmftc[iB3D]._x_f = self._bmftc[iB3D].B
                    # If x_f can't be found, barrier has drowned
                    self._bmftc[iB3D]._dur = self._time_step
                    self._bmftc[iB3D]._endyear = self._bmftc[iB3D].startyear + self._time_step
                    self._BMFTC_Break = True
                    print("PyBMFT-C Simulation Break: marsh has completely drowned or basin is completely full")
                    return  # End simulation

                # Store new positions
                self._bmftc[iB3D].Marsh_edge[self._bmftc[iB3D].startyear + self._time_step] = self._bmftc[ iB3D].x_m  # Save to array
                self._bmftc[iB3D].Forest_edge[self._bmftc[iB3D].startyear + self._time_step] = self._bmftc[iB3D].x_f  # Save to array

                # Determine change in x_b location
                self._x_b_TS[iB3D][self._time_step] = self._bmftc[iB3D].x_b  # Save to array

                # Determine new fetch based on change in opposite marsh - both fetches should be exactly the same!
                self._bmftc[iB3D]._bfo = self._bmftc[iB3D].bfo - progradation
                self._bmftc[iB3D].fetch[self._bmftc[iB3D].startyear + self._time_step] = self._bmftc[
                    iB3D].bfo  # Save to array

                # Update marsh scarp height parameter
                self._bmftc[iB3D]._dmo = self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + self._time_step] + \
                                         self._bmftc[iB3D].amp - \
                                         self._bmftc[iB3D].elevation[
                                             self._bmftc[iB3D].startyear + self._time_step, self._bmftc[iB3D].x_m]

                # Store landscape type widths for this time step
                if int(math.floor(self._x_s_offset[iB3D])) < 0:
                    barrier_width = len(self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                                        self._bmftc[iB3D].x_f:]) + off
                elif int(math.floor(self._x_s_offset[iB3D])) > 0:
                    barrier_width = len(self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                                        self._bmftc[iB3D].x_f:]) - off
                else:
                    barrier_width = len(self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                                        self._bmftc[iB3D].x_f:])
                BB_marsh_width = (
                        self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m: self._bmftc[iB3D].x_f] >
                        self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + self._time_step] + self._bmftc[iB3D].amp -
                        self._bmftc[iB3D].Dmax).sum()
                BB_marsh_pond_width = (
                        self._bmftc[iB3D].elevation[self._bmftc[iB3D].startyear + self._time_step,
                        self._bmftc[iB3D].x_m: self._bmftc[iB3D].x_f] <
                        self._bmftc[iB3D].msl[self._bmftc[iB3D].startyear + self._time_step] + self._bmftc[iB3D].amp -
                        self._bmftc[iB3D].Dmax).sum()
                self._LandscapeTypeWidth_TS[iB3D][self._time_step, :] = [barrier_width, BB_marsh_width, self._bmftc[iB3D].bfo, BB_marsh_pond_width]

        ###############################################################################
        # human dynamics modules
        ###############################################################################

        # ~~~~~~~~~~~~~~ RoadwayManager ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Remove overwash from roadway after each model year, place on the dune, rebuild dunes if
        # fall below height threshold, and check if dunes should grow naturally
        for iB3D in range(self._ny):

            if self._roadway_management_module[iB3D]:

                # if the roadway drowned or was too narrow for the road to be relocated, stop managing the road!
                # NOTE: dune heights must drop below Dmax before reset, so while calling reset_dune_growth rates seems
                # redundant, it doesn't slow us down computationally, so just do it
                if (
                    self._roadways[iB3D].drown_break
                    or self._roadways[iB3D].relocation_break
                ):
                    self._road_break[iB3D] = 1

                    # set dune growth rates back to original only when dune elevation is less than equilibrium
                    self._barrier3d[iB3D].growthparam = self.reset_dune_growth_rates(
                        original_growth_param=self._roadways[
                            iB3D
                        ]._original_growth_param,
                        iB3D=iB3D,
                    )
                else:
                    # manage that road!
                    self._roadways[iB3D].road_relocation_width = self._road_width[
                        iB3D
                    ]  # type: float
                    self._roadways[iB3D].road_relocation_setback = self._road_setback[
                        iB3D
                    ]
                    self._roadways[iB3D].update(self._barrier3d[iB3D])

                # update x_b to include a fake beach width and the dune line; we add a fake beach width for coupling
                # with the beach nourishment module below (i.e., if half the domain is initialized with roadways and
                # the other half with a community, I want them to start with the same beach back barrier position)
                self._barrier3d[iB3D].x_b_TS[-1] = (
                    self._barrier3d[iB3D].x_s
                    + self._barrier3d[iB3D].InteriorWidth_AvgTS[-1]
                    + np.size(
                        self._barrier3d[iB3D].DuneDomain, 2
                    )  # dune domain width in dam
                    + (self._initial_beach_width[iB3D] / 10)  # dam
                )

        # ~~~~~~~~~~~~~~ CHOM coupler ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # provide agents in the Coastal Home Ownership Model (CHOM) with variables describing the physical environment
        # -- including barrier elevation, beach width, dune height, shoreline erosion rate -- who then decide if it is
        # a nourishment year, the corresponding nourishment volume, and whether or not the dune should be rebuilt
        if self._community_dynamics_module:

            for iB3D in range(self._ny):

                # if barrier was too narrow to sustain a community in the last time step (from the BeachDuneManager),
                # stop the coupling with CHOM (i.e., end human mangement); dune growth rates are reset below in the
                # BeachDuneManager loop
                if self._nourishments[iB3D].narrow_break:
                    self._community_break[iB3D] = 1

            # update chom using all barrier3d grids, even if some have stopped being managed
            self._chom_coupler.dune_design_elevation = self._dune_design_elevation
            [
                self._nourish_now,
                self._rebuild_dune_now,
                self._nourishment_volume,
            ] = self._chom_coupler.update(
                barrier3d=self._barrier3d,
                nourishments=self._nourishments,
                community_break=self._community_break,
            )

        # ~~~~~~~~~~~~~~ BeachDuneManager ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # If interval specified, nourish at that interval, otherwise wait until told with nourish_now to nourish
        # or rebuild_dunes_now to rebuild dunes. Resets any "now" parameters to false after nourishment. Module
        # also filters overwash deposition for residential or commercial communities (user specified) and bulldozes
        # some of remaining overwash to dunes.
        for iB3D in range(self._ny):

            if self._beach_nourishment_module[iB3D]:
                # if barrier was too narrow to sustain a community in the last time step, stop managing beach and dunes!
                # NOTE: dune heights must drop below Dmax before reset, so while calling reset_dune_growth rates seems
                # redundant, it doesn't slow us down computationally, so just do it
                if self._nourishments[iB3D].narrow_break:
                    self._community_break[iB3D] = 1

                    # set dune growth rates back to original only when dune elevation is less than equilibrium
                    self._barrier3d[iB3D].growthparam = self.reset_dune_growth_rates(
                        original_growth_param=self._nourishments[
                            iB3D
                        ]._original_growth_param,
                        iB3D=iB3D,
                    )

                # else manage that community!
                else:
                    self._nourishments[
                        iB3D
                    ].dune_design_elevation = self._dune_design_elevation[
                        iB3D
                    ]  # m MHW
                    self._nourishments[
                        iB3D
                    ].nourishment_volume = self._nourishment_volume[iB3D]
                    [
                        self._nourish_now[iB3D],
                        self._rebuild_dune_now[iB3D],
                    ] = self._nourishments[iB3D].update(
                        barrier3d=self._barrier3d[iB3D],
                        nourish_now=self._nourish_now[iB3D],
                        rebuild_dune_now=self._rebuild_dune_now[iB3D],
                        nourishment_interval=self._nourishment_interval[iB3D],
                    )

                # update x_b to include a beach width and the dune line; after the community is abandoned, we set the
                # beach width for the remaining time steps to the last managed beach width in order to not have a huge
                # jump in the back-barrier position in Barrier3D
                self._barrier3d[iB3D].x_b_TS[-1] = (
                    self._barrier3d[iB3D].x_s
                    + self._barrier3d[iB3D].InteriorWidth_AvgTS[-1]
                    + np.size(
                        self._barrier3d[iB3D].DuneDomain, 2
                    )  # dune domain width in dam
                    + (
                        self._nourishments[iB3D].beach_width[
                            self._barrier3d[iB3D].time_index - 1
                        ]
                        / 10
                    )  # dam
                )

        ###############################################################################
        # update brie for any human modifications to the barrier
        ###############################################################################
        if self._alongshore_transport_module:
            [x_t, x_s, x_b, h_b, s_sf] = [np.zeros(self._ny) for _ in range(5)]

            for iB3D in range(self._ny):
                # make lists of the barrier geometry variables that have been changed (and needed to calculate shoreline
                # diffusivity in brie)
                x_t[iB3D] = self._barrier3d[iB3D].x_t_TS[-1]
                x_s[iB3D] = self._barrier3d[iB3D].x_s_TS[-1]
                x_b[iB3D] = self._barrier3d[iB3D].x_b_TS[-1]
                h_b[iB3D] = self._barrier3d[iB3D].h_b_TS[-1]
                s_sf[iB3D] = self._barrier3d[iB3D].s_sf_TS[-1]

            self._brie_coupler.update_brie_for_human_modifications(
                x_t, x_s, x_b, h_b, s_sf
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
