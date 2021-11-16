"""Beach and Dune Manager

This module provides functions for modifying a barrier segment from Barrier3D -- consisting of 1+ rows of dune cells
and a separate interior grid -- for beach and dune management decisions to protect a community, including:
    1) beach and shoreface nourishment at intervals or when user specified,
    2) filtering of overwash deposition for development and placement of overwash on dunes,
    3) dune rebuilding when beach nourished or user specified,
    4) limiting dune migration into the community.

References
----------

.. [1] Ashton, A. D., & Lorenzo-Trueba, J. (2018). Morphodynamics of barrier response to sea-level rise. In
        Barrier dynamics and response to changing climate (pp. 277-304). Springer, Cham.
.. [2] Rogers, L. J., Moore, L. J., Goldstein, E. B., Hein, C. J., Lorenzo‐Trueba, J., & Ashton, A. D. (2015).
        Anthropogenic controls on overwash deposition: Evidence and consequences. Journal of Geophysical Research:
        Earth Surface, 120(12), 2609-2624.

Notes
---------
The alongshore length of the barrier segment in Barrier3D is time-invariant, whereas the barrier interior
width -- and number of cross-shore cells -- varies dynamically due to storm impacts and RSLR.

In this module, we do not reduce the dune design elevation by SLR -- as in the RoadwayManager -- because here the
dune design height is always relative to mean sea level, and not the elevation of the roadway (which decreases with
RSLR). Note that dune elevations are decreased by SLR.

Barrier3D does not resolve the beach. To simulate a changing beach width with nourishment and subsequent shoreline
change, we establish a beach width based on the initial configuration of the beach (based on its slope) from Barrier3D
or user input. This beach width is then modified dynamically by nourishment and shoreface dynamics. We also turn dune
migration off in Barrier3D to keep the barrier location fixed, only allowing shoreline change to affect the beach width
and not seaward progradation of the dune line with each nourishment and or landward migration of the dune line with
shoreline erosion.

"""
import numpy as np
import copy
from .roadway_manager import rebuild_dunes, set_growth_parameters

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


def shoreface_nourishment(
    x_s, x_t, nourishment_volume, average_barrier_height, shoreface_depth, beach_width
):
    r"""
    Following the formulation used in Ashton and Lorenzo-Trueba (2018), we apply a nourishment volume (in m^3/m) along
    the entire shoreface, represented in Barrier3D by a single cross-shore transect. This results in a new shoreline
    position, beach width, and shoreface slope.

    Note, all parameters below must be in the same units. Show meters below, but for use with Barrier3D, decameters.

    Parameters
    ----------
    x_s: float
        Shoreline position [m]
    x_t: float
        Shoreface toe position [m]
    nourishment_volume: float
        Volume of sand per nourishment [m^3/m]
    average_barrier_height: float
        Average barrier height [m]
    shoreface_depth: float
        Depth of shoreface [m]
    beach_width: float
        beach_width [m]

    Returns
    -------
    new_x_s: float
        New shoreline position after nourishment [m]
    s_sf: float
        Shoreface slope after nourishment [unitless]
    beach_width: float
        Beach width after nourishment [m]
    """

    # move shoreline back
    new_x_s = x_s - (2 * nourishment_volume) / (
        2 * average_barrier_height + shoreface_depth
    )
    # nourishment_volume = beach_width * (2 * average_barrier_height + shoreface_depth) / 2

    # calculate new shoreface slope
    s_sf = shoreface_depth / (new_x_s - x_t)

    # beach width variable
    beach_width = beach_width + (x_s - new_x_s)

    return new_x_s, s_sf, beach_width


def resize_interior_domain(
    pre_storm_interior, post_storm_interior, bay_depth, dune_migration
):
    r"""
    Resize the pre- or post-storm interior domains if they are not the same size by padding with bay cells so the two
    domains can be easily differenced

    Parameters
    ----------
    pre_storm_interior: grid
        Pre-storm Barrier3D interior domain [dam * dam * dam]
    post_storm_interior: grid
        Post-storm Barrier3D interior domain [dam * dam * dam]
    bay_depth: float
        Bay depth [dam]
    dune_migration: float
        The number of grid cells [1 dam each] that the dunes migrated

    Returns
    -------
    grid
        pre_storm_interior: resized
        post_storm_interior: resized

    """

    # if pre-storm domain is larger than post-storm, remove all rows of bay without any deposition from the domain
    if np.size(pre_storm_interior, 0) > np.size(post_storm_interior, 0):

        # check first if the dunes migrated this past time step, this will make the interior domain smaller
        if dune_migration != 0:
            # remove the number of rows corresponding to the number of cells the dunes migrated from the pre-storm
            # domain (really the last time step); this happens if the user allows the beach width to fall below a min
            # threshold (10 m), which turns dune migration back on and allows for dune erosion in the post-storm domain
            for i in range(0, abs(int(dune_migration))):
                pre_storm_interior = np.delete(pre_storm_interior, 0, axis=0)
                if dune_migration > 0:
                    break  # break if the dune line aggrades, not set up for this

        if all(x <= -bay_depth for x in pre_storm_interior[-1, :]):
            pre_storm_interior = np.delete(pre_storm_interior, (-1), axis=0)

    # if post-storm domain larger than pre-storm, add rows to the bay of the pre-storm domain
    if np.size(post_storm_interior, 0) > np.size(pre_storm_interior, 0):
        number_rows = np.size(post_storm_interior, 0) - np.size(pre_storm_interior, 0)
        pre_storm_interior = np.concatenate(
            (
                pre_storm_interior,
                (np.zeros([number_rows, np.size(post_storm_interior, 1)]) - bay_depth),
            ),
            axis=0,
        )

    return pre_storm_interior, post_storm_interior


def filter_overwash(
    overwash_filter,
    post_storm_xyz_interior_grid,
    pre_storm_xyz_interior_grid,
    post_storm_yxz_dune_grid,
    artificial_maximum_dune_height,
):
    r"""
    Remove a percentage of overwash from the barrier interior, representative of the effect of development filtering
    overwash (Rogers et al., 2015).

    Following Rogers et al., 2015, we suggest that `overwash_filter` be set between 40 and 90%, where the lower
    (upper) end is representative of the effect of residential (commercial) buildings in reducing overwash deposition.

    Overwash that is removed from the barrier interior is placed evenly across the dunes. Note that this is different
    from the RoadwayManager, where overwash is only placed along the adjacent dunes in `bulldoze`.

    Note, all parameters below must be in the same units. For use with Barrier3D, decameters.

    Parameters
    ----------
    overwash_filter: int,
        Percent overwash removed from barrier interior [40-90% (residential-->commercial) from Rogers et al., 2015]
    post_storm_xyz_interior_grid: grid,
        Interior barrier island topography [for Barrier3d, decameters MHW]
    pre_storm_xyz_interior_grid: grid,
        Interior barrier island topography [for Barrier3d, decameters MHW]
    post_storm_yxz_dune_grid: grid,
        Dune topography [for Barrier3d, decameters above the berm elevation]
    artificial_maximum_dune_height: float,
        The maximum dune height than can be created by bulldozer after a storm [for Barrier3d, dam above berm elevation]

    Returns
    -------
    array of float
        new_interior_domain: in units of xyz_interior_grid
        new_dune_domain: in units of yxz_dune_grid
        total_overwash_removal: in units of dx*dy*dz (for Barrier3D, dam^3)

    """

    # remove sand from island interior (only account for positive values)
    overwash_deposition = post_storm_xyz_interior_grid - pre_storm_xyz_interior_grid
    overwash_deposition[overwash_deposition < 0] = 0

    # filter overwash deposition
    overwash_removal = overwash_deposition * (overwash_filter / 100)
    new_interior_domain = post_storm_xyz_interior_grid - overwash_removal

    # spread overwash removed from interior equally over all dune cells
    number_dune_cells = np.size(post_storm_yxz_dune_grid)
    total_overwash_removal = np.sum(overwash_removal)
    overwash_to_dune = total_overwash_removal / number_dune_cells
    new_dune_domain = post_storm_yxz_dune_grid + overwash_to_dune

    # don't allow dunes to exceed a maximum height (limits 10-m dunes after big storms...yikes!)
    new_dune_domain[
        new_dune_domain > artificial_maximum_dune_height
    ] = artificial_maximum_dune_height

    return (
        new_dune_domain,
        new_interior_domain,
        total_overwash_removal,
    )


def width_drown_checks(
    time_index,
    average_barrier_width,
    minimum_community_width=50,
):
    r"""Checks on whether a community -- specified by a minimum width -- is eaten up by the back-barrier. If it is,
    cease management.

    Parameters
    ----------
    time_index: int
        Time index for drowning error message
    average_barrier_width: float
        The average barrier width from the last time step [m]
    minimum_community_width: float
        The minimum width to sustain a community (i.e., a road plus a house footprint) [m]

    Returns
    -------
    boolean
        narrow_break: no room for the road plus a home...what else makes a community if not a road+home? Make room for
        an ark!

    """

    # initialize the break booleans as False
    narrow_break = 0

    # if the community gets eaten up by the back-barrier, then it is lost...and should no longer be managed
    if average_barrier_width <= minimum_community_width:
        narrow_break = 1
        print(
            "Community reached minimum width, drowned at {time} years".format(
                time=time_index - 1
            )
        )

    return narrow_break


class BeachDuneManager:
    """Manage those beaches and dunes to sustain a community!

    Examples
    --------
    # >>> from cascade.beach_dune_manager import BeachDuneManager
    # >>> nourish = BeachDuneManager()
    # >>> nourish.update(barrier3d, nourish_now, nourishment_interval)
    """

    def __init__(
        self,
        nourishment_interval=None,
        nourishment_volume=100,
        initial_beach_width=30,
        dune_design_elevation=3.7,
        time_step_count=500,
        original_growth_param=None,
        overwash_filter=40,
    ):
        """The BeachDuneManager module.

         Parameters
         ----------
        nourishment_interval: optional
             Interval that nourishment occurs [yrs]
        nourishment_volume: optional
             Volume of nourished sand along cross-shore transect [m^3/m]
        initial_beach_width: int, optional
             Initial beach width [m]
        dune_design_elevation: float, optional
            Elevation to which dune is rebuilt to [m MHW]; does not change with RSLR.
        time_step_count: int, optional
            Number of time steps.
        original_growth_param: optional
            Dune growth parameters from first time step of barrier3d, before human modifications [unitless]
        overwash_filter: float,
            Percent overwash removed from barrier interior [40-90% (residential-->commercial) from Rogers et al., 2015]
        """

        self._nourishment_volume = nourishment_volume
        self._nourishment_interval = nourishment_interval
        self._nourishment_counter = nourishment_interval
        self._beach_width_threshold = 0  # m, triggers dune migration to turn back on
        self._dune_design_elevation = dune_design_elevation
        self._original_growth_param = original_growth_param
        self._nt = time_step_count
        self._narrow_break = 0  # boolean for tracking drowning
        self._time_index = 1
        self._overwash_removal = True  # boolean for turning overwash removal on and off; for sensitivity testing
        self._overwash_filter = overwash_filter

        # variables that we don't change, but maybe someone else will want to in the future
        self._minimum_community_width = 50  # m, parameterized for Nags Head, NC
        self._artificial_maximum_dune_height = 4  # m, parameterized for Nags Head, NC
        # (maximum dune height above the berm crest that can be built...period

        # time series
        self._beach_width = [np.nan] * self._nt
        self._beach_width[0] = initial_beach_width  # m
        self._nourishment_TS = np.zeros(self._nt)
        self._dunes_rebuilt_TS = np.zeros(self._nt)
        self._nourishment_volume_TS = np.zeros(self._nt)  # m^3/m
        self._rebuild_dune_volume_TS = np.zeros(self._nt)  # m^3
        self._growth_params = [
            np.nan
        ] * self._nt  # track when dune growth parameters set to zero b/c of rebuild height
        self._growth_params[0] = original_growth_param
        self._overwash_volume_removed = np.zeros(
            self._nt
        )  # total overwash removed from the barrier [m^3], following the filtering effect of Rogers et al., 2015

        # also keep track of post-storm dune and interior impacts before human modifications, as well as pre-nourishment
        # shoreface configuration and beach width
        self._post_storm_dunes = [None] * self._nt
        self._post_storm_interior = [None] * self._nt
        self._post_storm_ave_interior_width = [None] * self._nt
        self._post_storm_ave_interior_height = [None] * self._nt
        self._post_storm_x_s = [None] * self._nt
        self._post_storm_s_sf = [None] * self._nt
        self._post_storm_beach_width = [None] * self._nt
        self._post_storm_Qow = [None] * self._nt
        self._post_storm_x_b = [None] * self._nt
        self._dune_migration_on = [np.nan] * self._nt
        self._dune_migration_on[0] = False

    def update(
        self,
        barrier3d,
        nourish_now,
        rebuild_dune_now,
        nourishment_interval,
    ):

        self._time_index = barrier3d.time_index

        if self._original_growth_param is None:
            self._original_growth_param = barrier3d.growthparam

        # if nourishment interval was updated in CASCADE, update here; otherwise just update the counter if it exists
        if self._nourishment_interval != nourishment_interval:
            self._nourishment_interval = nourishment_interval
            self._nourishment_counter = (
                self._nourishment_interval
            )  # reset the counter for the desired interval
        if self._nourishment_counter is not None:
            self._nourishment_counter -= 1

        # reduce beach width by the amount of shoreline change from last time step; if the beach width reaches zero,
        # turn dune migration in B3D back on -- otherwise keep it off (we don't want the dune line to prograde
        # because we have fake houses there!)
        self._beach_width[self._time_index - 1] = (
            self._beach_width[self._time_index - 2]
            - (barrier3d.x_s_TS[-1] - barrier3d.x_s_TS[-2]) * 10  # m
        )
        if (
            self._beach_width[self._time_index - 1] <= self._beach_width_threshold
        ):  # set to zero
            barrier3d.dune_migration_on = True
            barrier3d.SCRagg[self._time_index - 1] = (
                barrier3d.x_s % 1
            ) * -1  # set the shoreline change aggregate to the dam (cell) fraction
            if self._beach_width[self._time_index - 1] < 0:
                self._beach_width[self._time_index - 1] = 0
        else:
            barrier3d.dune_migration_on = False

        # save post-storm dune and interior impacts before human modifications, as well as pre-nourishment
        # shoreface configuration (essentially a 0.5 yr time step)
        self._post_storm_interior[self._time_index - 1] = copy.deepcopy(
            barrier3d.InteriorDomain
        )
        self._post_storm_dunes[self._time_index - 1] = copy.deepcopy(
            barrier3d.DuneDomain[self._time_index - 1, :, :]
        )
        self._post_storm_x_s[self._time_index - 1] = copy.deepcopy(barrier3d.x_s)
        self._post_storm_s_sf[self._time_index - 1] = copy.deepcopy(
            barrier3d.s_sf_TS[-1]
        )
        # in b3d, x_b = x_s + InteriorWidth_Avg and does not include dune domain or beach width, so we add that here
        # and save as "post storm" -- prior to any human changes in beach width or barrier interior width
        self._post_storm_x_b[self._time_index - 1] = (
            barrier3d.x_s
            + barrier3d.InteriorWidth_AvgTS[-1]
            + np.size(barrier3d.DuneDomain, 2)  # dune domain width in dam
            + (self._beach_width[self._time_index - 1] / 10)  # in dam
        )
        self._post_storm_Qow[self._time_index - 1] = copy.deepcopy(
            barrier3d.QowTS[-1]
        )  # m^3/m
        self._post_storm_ave_interior_width[self._time_index - 1] = copy.deepcopy(
            barrier3d.InteriorWidth_AvgTS[-1]
        )
        self._post_storm_ave_interior_height[self._time_index - 1] = copy.deepcopy(
            barrier3d.h_b_TS[-1]
        )
        self._post_storm_beach_width[self._time_index - 1] = self._beach_width[
            self._time_index - 1
        ]  # save post-storm beach width just in case it is modified later by beach nourishment

        # check for community width drowning prior to any management actions; if community cannot be sustained,
        # don't manage and exit; reset dune migration; dune growth reset in CASCADE
        average_barrier_width = barrier3d.InteriorWidth_AvgTS[-1] * 10
        self._narrow_break = width_drown_checks(
            self._time_index,
            average_barrier_width,  # m
            self._minimum_community_width,  # m
        )
        if self._narrow_break == 1:
            barrier3d.dune_migration_on = True
            self._dune_migration_on[
                self._time_index - 1
            ] = barrier3d.dune_migration_on  # keep track!
            barrier3d.SCRagg[self._time_index - 1] = (
                barrier3d.x_s % 1
            ) * -1  # set the shoreline change aggregate to the dam (cell) fraction
            return nourish_now, rebuild_dune_now

        # ------------------------------------- mgmt -------------------------------------

        # remove a percentage of overwash from the interior, representative of a community removing overwash from
        # roadways, driveways, and development filtering overwash; place this overwash back on the dunes
        if self._overwash_removal:

            # barrier3d doesn't save the pre-storm interior for each time step: therefore, use the output from the
            # previous time step [-2] and subtract SLR to obtain the pre-storm interior. Bay can't be deeper than
            # BayDepth (roughly equivalent to constant back-barrier slope), so update.
            pre_storm_interior = (
                barrier3d.DomainTS[self._time_index - 2]
                - barrier3d.RSLR[self._time_index - 1]
            )
            pre_storm_interior[
                pre_storm_interior < -barrier3d.BayDepth
            ] = -barrier3d.BayDepth

            pre_storm_interior, post_storm_interior = resize_interior_domain(
                pre_storm_interior=pre_storm_interior,
                post_storm_interior=barrier3d.DomainTS[self._time_index - 1],
                bay_depth=barrier3d.BayDepth,
                dune_migration=barrier3d.ShorelineChangeTS[
                    self._time_index - 1
                ],  # if -, dune erodes into interior [dam]
            )

            (
                new_yxz_dune_domain,  # [dam]
                new_xyz_interior_domain,  # [dam]
                barrier_overwash_removed,  # [dam^3]
            ) = filter_overwash(
                overwash_filter=self._overwash_filter,
                post_storm_xyz_interior_grid=post_storm_interior,  # dam
                pre_storm_xyz_interior_grid=pre_storm_interior,  # dam
                post_storm_yxz_dune_grid=barrier3d.DuneDomain[
                    self._time_index - 1, :, :
                ],  # dune domain from this last time step [dam]
                artificial_maximum_dune_height=self._artificial_maximum_dune_height
                / 10,  # dam
            )

            self._overwash_volume_removed[self._time_index - 1] = (
                barrier_overwash_removed * dm3_to_m3
            )  # convert from dm^3 to m^3
            net_overwash = barrier3d.QowTS[-1] - (
                self._overwash_volume_removed[self._time_index - 1]
                / (barrier3d.BarrierLength * 10)
            )  # m^3/m, post-storm - human removal
            new_domain_width, _, new_ave_interior_width = barrier3d.FindWidths(
                new_xyz_interior_domain, barrier3d.SL
            )
            new_ave_interior_height = np.average(
                new_xyz_interior_domain[new_xyz_interior_domain >= barrier3d.SL]
            )

            # set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the
            # natural eq. dune height (Dmax) -- do this because overwash volumes can be very large
            new_growth_parameters = set_growth_parameters(
                new_yxz_dune_domain,  # in dam
                barrier3d.Dmax,  # in dam
                barrier3d.growthparam,
                original_growth_param=self._original_growth_param,  # use original growth rates for resetting values
            )
            self._growth_params[self._time_index - 1] = copy.deepcopy(
                new_growth_parameters
            )

            # update Barrier3D class variables
            barrier3d.DuneDomain[self._time_index - 1, :, :] = copy.deepcopy(
                new_yxz_dune_domain
            )
            barrier3d.InteriorDomain = new_xyz_interior_domain
            barrier3d.DomainTS[self._time_index - 1] = new_xyz_interior_domain
            barrier3d._DomainWidth = new_domain_width
            barrier3d.growthparam = new_growth_parameters
            barrier3d.QowTS[-1] = net_overwash  # m^3/m
            barrier3d.InteriorWidth_AvgTS[-1] = new_ave_interior_width  # dam
            barrier3d.h_b_TS[-1] = new_ave_interior_height  # dam
            # note x_b_TS is modified at the end

        # if specified, rebuild dune (if using nourishment counter option, nourishes dune automatically)
        if rebuild_dune_now or self._nourishment_counter == 0:

            # in B3D, dune height is the height above the berm crest
            dune_design_height = self._dune_design_elevation - (barrier3d.BermEl * 10)
            new_dune_domain, rebuild_dune_volume = rebuild_dunes(
                barrier3d.DuneDomain[self._time_index - 1, :, :],  # dam
                max_dune_height=dune_design_height,  # in m
                min_dune_height=dune_design_height,  # in m
                dz=10,  # specifies dune domain is in dam
                rng=True,  # adds stochasticity to dune height (seeded)
            )
            self._dunes_rebuilt_TS[self._time_index - 1] = 1
            self._rebuild_dune_volume_TS[self._time_index - 1] = (
                rebuild_dune_volume * dm3_to_m3
            )  # m^3

            # [COMING SOON] check beach width and find associated DMAX, change DMAX and feed into next function

            # set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the
            # natural eq. dune height (Dmax)
            new_growth_parameters = set_growth_parameters(
                new_dune_domain,  # in dam
                barrier3d.Dmax,  # in dam
                barrier3d.growthparam,
                original_growth_param=self._original_growth_param,  # use original growth rates for resetting values
            )
            self._growth_params[self._time_index - 1] = copy.deepcopy(
                new_growth_parameters
            )

            # update Barrier3d class dune variables
            barrier3d.DuneDomain[self._time_index - 1, :, :] = copy.deepcopy(
                new_dune_domain
            )
            barrier3d.growthparam = new_growth_parameters

            # reset rebuild_dune_now parameter; if nourishment_counter is what triggered the rebuild, it is reset in the
            # nourishment section below
            rebuild_dune_now = 0

        # finally, nourish the shoreface
        if nourish_now or self._nourishment_counter == 0:
            (
                barrier3d.x_s,  # save over class variables
                barrier3d.s_sf_TS[-1],
                self._beach_width[self._time_index - 1],
            ) = shoreface_nourishment(
                barrier3d.x_s,  # in dam
                barrier3d.x_t,  # in dam
                self._nourishment_volume / 100,  # convert m^3/m to dam^3/dam
                barrier3d.h_b_TS[-1],  # in dam
                barrier3d.DShoreface,  # in dam
                self._beach_width[self._time_index - 1] / 10,  # convert m to dam
            )
            self._beach_width[self._time_index - 1] *= 10  # convert dam back to m
            self._nourishment_TS[self._time_index - 1] = 1
            self._nourishment_volume_TS[self._time_index - 1] = self._nourishment_volume

            # reset counter if its what triggered nourishment and nourish_now parameter
            if self._nourishment_counter is not None:
                self._nourishment_counter = self._nourishment_interval
            nourish_now = 0

            # set dune migration off after nourishment (we don't want the dune line to prograde if
            # the beach width was previously less than threshold)
            barrier3d.dune_migration_on = False

        # update shoreline time series (just in case nourishment happened) and back-barrier location
        barrier3d.x_s_TS[-1] = barrier3d.x_s
        barrier3d.x_b_TS[-1] = (
            barrier3d.x_s
            + barrier3d.InteriorWidth_AvgTS[-1]
            + np.size(barrier3d.DuneDomain, 2)  # dune domain width in dam
            + (self._beach_width[self._time_index - 1] / 10)  # in dam
        )
        self._dune_migration_on[
            self._time_index - 1
        ] = barrier3d.dune_migration_on  # keep track!

        return nourish_now, rebuild_dune_now

    @property
    def beach_width(self):
        return self._beach_width

    @beach_width.setter
    def beach_width(self, value):
        self._beach_width = value

    @property
    def dune_design_elevation(self):
        return self._dune_design_elevation

    @dune_design_elevation.setter
    def dune_design_elevation(self, value):
        self._dune_design_elevation = value

    @property
    def nourishment_volume(self):
        return self._nourishment_volume

    @nourishment_volume.setter
    def nourishment_volume(self, value):
        self._nourishment_volume = value

    @property
    def overwash_removal(self):
        return self._overwash_removal

    @overwash_removal.setter
    def overwash_removal(self, value):
        self._overwash_removal = value

    @property
    def nourishment_volume_TS(self):
        return self._nourishment_volume_TS

    @property
    def rebuild_dune_volume_TS(self):
        return self._rebuild_dune_volume_TS

    @property
    def overwash_volume_removed(self):
        return self._overwash_volume_removed

    @property
    def narrow_break(self):
        return self._narrow_break
