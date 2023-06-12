from joblib import Parallel, delayed
import numpy as np
import os

from .roadway_manager import RoadwayManager, set_growth_parameters
from .beach_dune_manager import BeachDuneManager
from .outwasher_lateral_transport import Outwasher
from .brie_coupler import BrieCoupler, initialize_equal, batchB3D
from .chom_coupler import ChomCoupler

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
        outwash_module
    ):
        """Configures lists to account for multiple Barrier3D domains from single input variables; used in modules"""

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
        if np.size(outwash_module) > 1:
            self._outwash_module = outwash_module
        else:
            self._outwash_module = [outwash_module] * self._ny

        return

    def reset_dune_growth_rates(
        self,
        original_growth_param,
        iB3D,
    ):
        """Reset dune growth parameters to original values after human abandonment; dune heights must drop below Dmax
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
        elevation_file="barrier3d-default-elevation.npy",
        dune_file="barrier3d-default-dunes.npy",
        parameter_file="barrier3d-default-parameters.yaml",
        storm_file="cascade-default-storms.npy",  # same as "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy"
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=True,
        beach_nourishment_module=True,
        community_economics_module=False,
        outwash_module=True,
        alongshore_section_count=6,
        time_step_count=200,
        wave_height=1, # ---------- for BRIE and Barrier3D --------------- #
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        bay_depth=3.0,
        s_background=0.001,
        berm_elevation=1.9,
        MHW=0.46,
        beta=0.04,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        min_dune_growth_rate=0.25,
        max_dune_growth_rate=0.65,
        road_ele=1.7,  # ---------- roadway management --------------- #
        road_width=30,
        road_setback=30,
        dune_design_elevation=3.7,
        dune_minimum_elevation=2.2,
        trigger_dune_knockdown=False,
        group_roadway_abandonment=None,
        nourishment_interval=None, # ---------- beach and dune ("community") management --------------- #
        nourishment_volume=300.0,
        overwash_filter=40,
        overwash_to_dune=10,
        number_of_communities=1, # ---------- coastal real estate markets (in development) --------------- #
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
        outwash_storms_file="outwash_storms10.npy",  # --------- outwasher (in development) ------------ #
        outwash_beach_file="NCB-default_beach.npy",
        percent_washout_to_shoreface=100,
        dune_flow_dynamics="full",
        outwasher_substep=20,
        ki_value=7.5E-3,
        cx=10,
        kL=1,
    ):
        """

        CASCADE: The CoAStal Community-lAnDscape Evolution model

        Couples Barrier3D (Reeves et al., 2019) with the Barrier Inlet Environment Model (BRIE; Nienhuis and Lorenzo
        Trueba, 2019) & the Coastal Home Ownership Model (CHOM), an agent-based model for coastal real estate markets
        (Williams et al., in prep)

        Parameters
        ----------
        datadir: string
            Name of directory where Barrier3D and Outwasher input file is located
        name: string, optional
            Name of simulation
        wave_height: float, optional
            Deepwater significant wave height [m]
        wave_period: float, optional
            Peak wave period [s]
        wave_asymmetry: float, optional
            Fraction of waves approaching from left (looking onshore)
        wave_angle_high_fraction: float, optional
            Fraction of waves approaching from angles higher than 45 degrees
        bay_depth: float, optional
            Depth of back-barrier bay [m]
        s_background: float, optional
            Background slope (for shoreface toe position, back-barrier & inlet calculations)
        berm_elevation: float, optional
            Static elevation of berm [m NAVD88]; needs to be 1.9 m if using the default storm list, time series
        MHW: float, optional
            Elevation of mean high water [m NAVD88]; needs to be 0.46 m NAVD88 if using default storm list, time series
        beta: float, optional
            Beach slope for runup calculations, needs to be 0.04 if using the default storm list, time series
        sea_level_rise_rate: float, optional
            Rate of sea_level rise (SLR) [m/yr]
        sea_level_rise_constant: boolean, optional
            If True, linear SLR; otherwise a logistic growth function is used for acc SLR (max 200 year simulations)
        background_erosion: float,
            Rate of shoreline retreat attributed to gradients in alongshore transport; (-) = erosion, (+) = acc [m / y]
        alongshore_section_count: int, optional
            Number of alongshore sections
        time_step_count: int, optional
            Number of time steps
        min_dune_growth_rate: float or list of floats, optional
            Minimum dune growth rate [unitless]; for Houser et al., (2015) growth rate formulation
        max_dune_growth_rate: float or list of floats, optional
            Maximum dune growth rate [unitless]; for Houser et al., (2015) growth rate formulation
        num_cores: int, optional
            Number of (parallel) processing cores to be used; helpful to have >1 for multiple Barrier3D segments
        roadway_management_module: boolean or list of booleans, optional
            If True, use roadway management module (overwash removal, road relocation, dune management)
        alongshore_transport_module: boolean or list of booleans, optional
            If True, couple Barrier3D with BRIE to use diffusive alongshore sediment transport module
        community_economics_module: boolean or list of booleans, optional
            If True, couple with CHOM, a community decision making model; requires nourishment module (in development)
        beach_nourishment_module: boolean or list of booleans, optional
            If True, use nourishment module (nourish shoreface, rebuild dunes)
        road_ele: float or list of floats, optional
            Elevation of the initial roadway [m MHW]
        road_width: int or list of int, optional
            Width of roadway [m]
        road_setback: int or list of int, optional
            Setback of roadway from the inital dune line and after road relocations [m]
        dune_design_elevation: float or list of floats, optional
            Elevation to which dune is initially rebuilt [m MHW]
        dune_minimum_elevation: float or list of floats, optional
            Elevation threshold which triggers rebuilding of dune for roadway management [m MHW]
        trigger_dune_knockdown: boolean, optional
            Resets the dune elevation to the initial condition (time zero) after roadway abandonment
        group_roadway_abandonment: list of ints greater than zero, optional
            Groups roadways together into segments for abandonment (i.e., if one roadway is abandoned they all are)
        nourishment_interval: int or list of ints, optional
             Interval that nourishment occurs [yrs]
        nourishment_volume: float or list of float, optional
             Volume of nourished sand along cross-shore transect [m^3/m]
        overwash_filter: float or list of floats,
            Percent overwash removed from barrier interior [40-90% (residential-->commercial) from Rogers et al., 2015]
        overwash_to_dune: float or list of floats,
            Percent overwash removed from barrier interior to dunes [%]; overwash_filter+overwash_to_dune <=100
        number_of_communities: int, optional
            Number of communities (CHOM model instances) described by the alongshore section count (Barrier3D models)
        sand_cost: int, optional
            Unit cost of sand $/m^3
        taxratio_oceanfront: float, optional
            The proportion of tax burden placed on oceanfront row for beach management
        external_housing_market_value_oceanfront: float, optional
            Value of comparable housing option outside the coastal system
        external_housing_market_value_nonoceanfront: float, optional
            Value of comparable housing option outside the coastal system
        fixed_cost_beach_nourishment: int, optional
            Fixed cost of 1 nourishment project
        fixed_cost_dune_nourishment: int, optional
            Fixed cost of building dunes once
        nourishment_cost_subsidy: int, optional
            Subsidy on cost of entire nourishment plan
        beach_full_cross_shore: int, optional
            The cross-shore extent (meters) of fully nourished beach (i.e., the community desired beach width) [m]
        outwash_storms: string, optional
            Filename of outwash storm series (npy file)
        washout_to_shoreface: bool
            if True, washout is used to nourish the shoreface
        outwash_module: boolean or list of booleans, optional
            If True, use outwash module (force a bay-side surge event)

        Examples
        --------
        >>> from cascade.cascade import Cascade
        >>> datadir = "./cascade/data/"
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
        self._community_economics_module = community_economics_module
        self._filename = name
        self._storm_file = storm_file
        self._elevation_file = elevation_file
        self._dune_file = dune_file
        self._parameter_file = parameter_file
        self._number_of_communities = number_of_communities
        self._b3d_break = 0
        self._road_break = [
            0
        ] * self._ny
        self._community_break = [
            0
        ] * self._ny
        self._nourish_now = [0] * self._ny  # user can trigger nourishment in time loop
        self._rebuild_dune_now = [0] * self._ny  # user can trigger the dune to be rebuilt in time loop
        self._trigger_dune_knockdown = (
            trigger_dune_knockdown  # user can force the dunes to be knocked down in time loop
        )
        self._initial_beach_width = [0] * self._ny
        self._group_roadway_abandonment = group_roadway_abandonment

        # initialization errors
        if (berm_elevation != 1.9 or MHW !=0.46 or beta !=0.04) and storm_file=="barrier3d-default-storms.npy":
            raise CascadeError(
                "The default storms only apply for a berm elevation=1.9 m NAVD88, MHW=0.46 m NAVD88 & beach slope=0.04."
            )
        if (sea_level_rise_constant is False) and (time_step_count > 200):
            raise CascadeError(
                "The sigmoidal accelerated SLR formulation used in this model by Rohling et al., (2013) should not be"
                "extended beyond 200 years"
            )

        ###############################################################################
        # initialize BRIE and Barrier3D model classes
        ###############################################################################

        # initialize BRIE: used for initial shoreface calculations, AST (optional), tidal inlets (eventually)
        self._brie_coupler = BrieCoupler(
            name=name,
            wave_height=self._wave_height,
            wave_period=self._wave_period,
            wave_asymmetry=self._wave_asymmetry,
            wave_angle_high_fraction=self._wave_angle_high_fraction,
            sea_level_rise_rate=self._sea_level_rise_rate,
            back_barrier_depth=bay_depth,
            s_background=s_background,
            h_b_crit=(berm_elevation - MHW),
            ny=self._ny,
            nt=self._nt,
        )

        # initialize Barrier3D models (number set by brie_ny) and make both "brie" and "barrier3d" classes equivalent
        self._barrier3d = initialize_equal(
            datadir=datadir,
            brie=self._brie_coupler._brie,
            slr_constant=self._slr_constant,
            rmin=self._rmin,  # can be array
            rmax=self._rmax,  # can be array
            background_erosion=self._background_erosion,  # can be array
            MHW=MHW,
            berm_elevation=berm_elevation,
            beta=beta,
            parameter_file=self._parameter_file,
            storm_file=self._storm_file,
            dune_file=self._dune_file,  # can be array
            elevation_file=self._elevation_file,  # can be array
        )

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
            outwash_module=outwash_module,
        )

        if self._community_economics_module:
            if not any(self._beach_nourishment_module):
                raise CascadeError(
                    "Beach nourishment module must be set to `TRUE` to couple with CHOM"
                )
            else:
                self._chom_coupler = ChomCoupler(
                    barrier3d=self._barrier3d,
                    total_time=self._nt,
                    alongshore_length_b3d=self._brie_coupler._brie._dy,  # this is the Barrier3D default, 500 m
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
        # (always initialize just in case we want to add a road or start nourishing during the simulation)
        self._roadways = []
        self._nourishments = []

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

        # use the initial beach width as a check on the Barrier3D user input for mulitple domains; the beach width
        # must be the same for all domains because there is only one storm file, which is made for a set berm
        # elevation and beach slope
        if all(
            elem == self._initial_beach_width[0] for elem in self._initial_beach_width
        ):
            pass
        else:
            raise CascadeError(
                "Berm elevation and beach slope must be equivalent for all Barrier3D domains"
            )

        ###############################################################################
        # initialize outwasher
        ###############################################################################
        # (always initialize just in case...)
        self._outwash = []

        for iB3D in range(self._ny):
            self._outwash.append(
                Outwasher(datadir=datadir,
                          outwash_storms_file=outwash_storms_file,
                          time_step_count=self._nt,
                          berm_elev=self._barrier3d[iB3D].BermEl,
                          barrier_length=self._barrier3d[iB3D].BarrierLength,
                          sea_level=self._barrier3d[iB3D].SL,
                          bay_depth=self._barrier3d[iB3D].BayDepth,
                          interior_domain=self._barrier3d[iB3D].InteriorDomain,
                          dune_domain=self._barrier3d[iB3D].DuneDomain[self._barrier3d[iB3D].time_index - 1, :, :],
                          percent_washout_to_shoreface=percent_washout_to_shoreface,
                          outwash_beach_file=outwash_beach_file,
                          dune_flow_dynamics=dune_flow_dynamics,
                          substep=outwasher_substep,
                          sediment_flux_coefficient_Ki=ki_value,
                          lateral_trans_coeff_kL=kL,
                          cx=cx)
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

    @property
    def outwash(self):
        return self._outwash

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

    def update(self):
        """Update cascade by a single time step"""

        # check for drowning from the last time step in brie. Note that this will stay false if brie is not used for AST
        if self._brie_coupler._brie.drown == True:
            return

        # advance B3D by one time step (B3D initializes at time_index = 1 and then updates the time_index after
        # update_dune_domain). Set n_jobs=1 for no parallel processing (debugging) and -2 for all but 1 CPU;
        # note that joblib uses a threshold on the size of arrays passed to the workers
        batch_output = Parallel(n_jobs=self._num_cores, max_nbytes="10M")(
            delayed(batchB3D)(self._barrier3d[iB3D]) for iB3D in range(self._ny)
        )

        # reshape output from parallel processing and convert from tuple to list
        x_t_dt, x_s_dt, h_b_dt, b3d = zip(*batch_output)
        x_t_dt = list(x_t_dt)
        x_s_dt = list(x_s_dt)
        h_b_dt = list(h_b_dt)
        self._barrier3d = list(b3d)

        # use brie to connect B3D models with AST; otherwise, just update (erode/prograde) dune domain
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
        # human dynamics modules
        ###############################################################################

        # ~~~~~~~~~~~~~~ RoadwayManager ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Remove overwash from roadway after each model year, place on the dune, rebuild dunes if
        # fall below height threshold, and check if dunes should grow naturally.
        for iB3D in range(self._ny):

            if self._roadway_management_module[iB3D]:

                # if the roadway drowned or was too narrow for the road to be relocated, stop managing the road!
                # NOTE: dune heights must drop below Dmax before reset, so while calling reset_dune_growth rates seems
                # redundant, it doesn't slow us down computationally, so just do it
                if (
                    self._roadways[iB3D].drown_break
                    or self._roadways[iB3D].relocation_break
                ):
                    # if it was specified that roadways are abandonded in groups, abandon the group
                    if self._group_roadway_abandonment is not None:

                        # find the group indices
                        group_roadways = np.where(
                            np.array(self._group_roadway_abandonment)
                            == self._group_roadway_abandonment[iB3D]
                        )[0]
                        self._road_break[group_roadways[0] : group_roadways[-1] + 1] = [
                            1
                        ] * len(group_roadways)

                        # label the other group indices as broken
                        for iRoad in group_roadways:
                            if self._roadways[iB3D].drown_break:
                                self._roadways[iRoad].drown_break = 1
                            else:
                                self._roadways[iRoad].relocation_break = 1

                            # set dune growth rates back to original only when dune elevation is less than equilibrium
                            self._barrier3d[
                                iRoad
                            ].growthparam = self.reset_dune_growth_rates(
                                original_growth_param=self._roadways[
                                    iRoad
                                ]._original_growth_param,
                                iB3D=iRoad,
                            )

                    else:
                        self._road_break[iB3D] = 1

                        # set dune growth rates back to original only when dune elevation is less than equilibrium
                        self._barrier3d[
                            iB3D
                        ].growthparam = self.reset_dune_growth_rates(
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
                    self._roadways[iB3D].update(
                        self._barrier3d[iB3D], self._trigger_dune_knockdown
                    )

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

        # ~~~~~~~~~~~~~~ CHOM coupler (in development) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Provide agents in the Coastal Home Ownership Model (CHOM) with variables describing the physical environment
        # -- including barrier elevation, beach width, dune height, shoreline erosion rate -- who then decide if it is
        # a nourishment year, the corresponding nourishment volume, and whether or not the dune should be rebuilt
        if self._community_economics_module:

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
        # outwash module
        ###############################################################################
        # Simulates bay-to-ocean flow for one storm. This modifies the B3D interior domain and adjusts the shoreline
        # position.
        for iB3D in range(self._ny):

            if self._outwash_module[iB3D]:
                self._outwash[iB3D]._interior_domain = self._barrier3d[iB3D].InteriorDomain
                self._outwash[iB3D]._dune_domain = self._barrier3d[iB3D].DuneDomain[self._barrier3d[iB3D].time_index - 1]
                self._outwash[iB3D].update(b3d=self._barrier3d[iB3D])

        ###############################################################################
        # update BRIE for any human modifications to the barrier
        ###############################################################################
        if self._alongshore_transport_module:
            [x_t, x_s, x_b, h_b, s_sf] = [np.zeros(self._ny) for _ in range(5)]

            for iB3D in range(self._ny):
                # make lists of the barrier geometry variables that have been changed (and needed to calculate shoreline
                # diffusivity in BRIE)
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
