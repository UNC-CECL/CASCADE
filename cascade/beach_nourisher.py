"""Beach Nourisher

This module provides functions ....  ADD DESCRIPTION HERE ONCE FINALIZED

References
----------

.. [1] Ashton, A. D., & Lorenzo-Trueba, J. (2018). Morphodynamics of barrier response to sea-level rise. In
        Barrier dynamics and response to changing climate (pp. 277-304). Springer, Cham.
.. [2] Rogers, L. J., Moore, L. J., Goldstein, E. B., Hein, C. J., Lorenzo‐Trueba, J., & Ashton, A. D. (2015).
        Anthropogenic controls on overwash deposition: Evidence and consequences. Journal of Geophysical Research:
        Earth Surface, 120(12), 2609-2624.


Notes
---------

"""
import numpy as np
import copy
from .roadway_manager import rebuild_dunes, set_growth_parameters

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


def shoreface_nourishment(
    x_s, x_t, nourishment_volume, average_barrier_height, shoreface_depth, beach_width
):
    r"""
    Following the formulation used in Ashton, A. D., & Lorenzo-Trueba, J. (2018), we apply a nourishment volume
    (in m^3/m) along the entire shoreface.

    Note, all parameters below must be in the same units. Show meters below, but for use with barrier3D, decameters.

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
        New shoreline position after nourishment
    s_sf: float
        Shoreface slope after nourishment
    beach_width: float
        Beach width after nourishment
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


def resize_interior_domain(pre_storm_interior, post_storm_interior, bay_depth):
    """If the pre-storm interior domains are not the same size, resize appropriately so we can easily diff"""

    # if pre-storm domain is larger than post-storm, remove all rows of bay without any deposition from the domain
    if np.size(pre_storm_interior, 0) > np.size(post_storm_interior, 0):
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
):
    r"""
    Remove a percentage of overwash from the interior, representative of the effect of development filtering
    overwash (Rogers et al., 2015).

    Following Rogers et al., 2015, we suggest that `overwash_filter` be set between 40 and 90%, where the lower
    (upper) end is representative of the effect of residential (commercial) buildings in reducing overwash deposition.
    Overwash that is removed from the barrier interior is placed evenly across the dunes.

    Note, all parameters below must be in the same units. For use with barrier3D, decameters.

    See also `bulldoze` in roadway_manager.py

    Parameters
    ----------
    overwash_filter: float,
        Percent overwash removed from barrier interior [40-90% (residential-->commercial) from Rogers et al., 2015]
    post_storm_xyz_interior_grid: grid,
        Interior barrier island topography [for barrier3d, decameters MHW (NAVD88)]
    pre_storm_xyz_interior_grid: grid,
        Interior barrier island topography [for barrier3d, decameters MHW (NAVD88)]
    post_storm_yxz_dune_grid: grid,
        Dune topography [for barrier3d, decameters]

    Returns
    -------
    array of float
        new_interior_domain: in units of xyz_interior_grid
        new_dune_domain: in units of yxz_dune_grid
        overwash removal: in units of dz^3

    """

    # remove sand from island interior (only account for positive values)
    overwash_deposition = post_storm_xyz_interior_grid - pre_storm_xyz_interior_grid
    overwash_deposition[overwash_deposition < 0] = 0

    # filter overwash deposition
    overwash_removal = overwash_deposition * (overwash_filter / 100)
    new_interior_domain = post_storm_xyz_interior_grid - overwash_removal

    # spread overwash removed from interior equally over all dune cells
    number_dune_cells = np.size(post_storm_yxz_dune_grid, 1)
    overwash_removal = sum(overwash_removal)
    overwash_to_dune = np.transpose(
        [overwash_removal / number_dune_cells] * number_dune_cells
    )
    new_dune_domain = post_storm_yxz_dune_grid + overwash_to_dune

    return (
        new_dune_domain,
        new_interior_domain,
        np.sum(overwash_removal),
    )


def width_drown_checks(
    time_index,
    average_barrier_width,
    minimum_community_width,
):
    # initialize the break booleans as False
    narrow_break = 0

    # if the community gets eaten up by the back-barrier, then it is lost
    if average_barrier_width <= minimum_community_width:
        narrow_break = 1
        print(
            "Community reached minimum width, drowned at {time} years".format(
                time=time_index - 1
            )
        )

    return narrow_break


class BeachNourisher:
    """Nourish the beach

    Examples
    --------
    # >>> from cascade.beach_nourisher import BeachNourisher
    # >>> nourish = BeachNourisher()
    # >>> nourish.update(barrier3d, nourish_now, nourishment_interval)
    """

    def __init__(
        self,
        nourishment_interval=None,
        nourishment_volume=100,
        initial_beach_width=10,
        dune_design_elevation=3.7,
        time_step_count=500,
        original_growth_param=None,
        overwash_filter=40,
    ):
        """The BeachNourisher module.

         Parameters
         ----------
        nourishment_interval: optional
             Interval that nourishment occurs [yrs]
        nourishment_volume: optional
             Volume of nourished sand along cross-shore transect [m^3/m]
        initial_beach_width: int, optional
             initial beach width [m]
        dune_design_elevation: float, optional
            Elevation to which dune is rebuilt to [m NAVD88]
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
        self._beach_width_threshold = 10  # m
        self._dune_design_elevation = dune_design_elevation
        self._original_growth_param = original_growth_param
        self._nt = time_step_count
        self._narrow_break = 0  # boolean for tracking drowning
        self._time_index = 1
        self._overwash_removal = True  # boolean for turning overwash removal on and off; for sensitivity testing
        self._overwash_filter = overwash_filter
        self._minimum_community_width = 80  # m

        # time series
        self._beach_width = [None] * self._nt
        self._beach_width[0] = initial_beach_width  # m
        self._nourishment_TS = np.zeros(self._nt)
        self._dunes_rebuilt_TS = np.zeros(self._nt)
        self._nourishment_volume_TS = np.zeros(self._nt)  # m^3/m
        self._rebuild_dune_volume_TS = np.zeros(self._nt)  # m^3
        self._growth_params = [
            None
        ] * self._nt  # track when dune growth parameters set to zero b/c of rebuild height
        self._growth_params[0] = original_growth_param
        self._post_storm_dunes = [
            None
        ] * self._nt  # keep track of post-storm dune impacts before human modifications
        self._post_storm_interior = [
            None
        ] * self._nt  # keep track of post-storm interior impacts before human modifications
        self._overwash_volume_removed = np.zeros(
            self._nt
        )  # total overwash removed from the barrier [m^3], following the filtering effect of Rogers et al., 2015

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

        # if nourishment interval was updated in CASCADE, update here
        if self._nourishment_interval != nourishment_interval:
            self._nourishment_interval = nourishment_interval
            self._nourishment_counter = (
                self._nourishment_interval
            )  # reset the counter for the desired interval

        # save post-storm dune and interior domain before human modifications
        self._post_storm_interior[self._time_index - 1] = copy.deepcopy(
            barrier3d.InteriorDomain
        )  # hoping this makes a deep copy
        self._post_storm_dunes[self._time_index - 1] = copy.deepcopy(
            barrier3d.DuneDomain[self._time_index - 1, :, :]
        )

        # reduce beach width by the amount of shoreline change from last time step;  if the beach width is less than
        # some threshold (we use 10 m), turn dune migration in Barrier3D back on
        self._beach_width[self._time_index - 1] = (
            self._beach_width[self._time_index - 2]
            - (barrier3d.x_s_TS[-1] - barrier3d.x_s_TS[-2]) * 10  # m
        )
        if self._beach_width[self._time_index - 1] <= self._beach_width_threshold:
            barrier3d.dune_migration_on = True

        # update nourishment counter
        if self._nourishment_counter is not None:
            self._nourishment_counter -= 1

        # check for barrier width drowning; if barrier drowns, exit program
        average_barrier_width = barrier3d.InteriorWidth_AvgTS[-1] * 10
        self._narrow_break = width_drown_checks(
            self._time_index,
            average_barrier_width,  # m
            self._minimum_community_width,  # m
        )
        if self._narrow_break == 1:
            return

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
            )

            self._overwash_volume_removed[self._time_index - 1] = (
                barrier_overwash_removed * dm3_to_m3
            )  # convert from dm^3 to m^3

            # update class variables
            barrier3d.DuneDomain[self._time_index - 1, :, :] = new_yxz_dune_domain
            barrier3d.InteriorDomain = new_xyz_interior_domain
            barrier3d.DomainTS[self._time_index - 1] = new_xyz_interior_domain

        # if specified, rebuild dune and nourish shoreface
        if rebuild_dune_now or self._nourishment_counter == 0:
            artificial_max_dune_height = self._dune_design_elevation - (
                barrier3d.BermEl * 10
            )
            new_dune_domain, rebuild_dune_volume = rebuild_dunes(
                barrier3d.DuneDomain[self._time_index - 1, :, :],  # dam
                max_dune_height=artificial_max_dune_height,  # in m
                min_dune_height=artificial_max_dune_height,  # in m
                dz=10,  # specifies dune domain is in dam
                rng=True,  # adds stochasticity to dune height (seeded)
            )
            self._dunes_rebuilt_TS[self._time_index - 1] = 1
            self._rebuild_dune_volume_TS[self._time_index - 1] = (
                rebuild_dune_volume * dm3_to_m3
            )

            # [COMING SOON] check beach width and find associated DMAX, change DMAX and feed into next function

            # set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the
            # natural eq. dune height (Dmax)
            new_growth_parameters = set_growth_parameters(
                new_dune_domain,
                barrier3d.Dmax,
                barrier3d.growthparam,
                original_growth_param=self._original_growth_param,  # use original growth rates for resetting values
            )
            self._growth_params[self._time_index - 1] = copy.deepcopy(
                new_growth_parameters
            )

            # update barrier3d dune variables
            barrier3d.DuneDomain[self._time_index - 1, :, :] = new_dune_domain
            barrier3d.growthparam = new_growth_parameters

            # reset rebuild_dune_now parameter; if nourishment_counter is what triggered the rebuild, it is reset in the
            # nourishment section below
            rebuild_dune_now = 0

        if nourish_now or self._nourishment_counter == 0:
            (
                barrier3d.x_s,
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

            # set dune migration off after nourishment (we don't want the dune line to prograde)
            barrier3d.dune_migration_on = False

        # update shoreline time series (just in case nourishment happened) and back-barrier location
        # (this needs to be done every year because b3d doesn't include dune domain width when calculating back barrier)
        barrier3d.x_s_TS[-1] = barrier3d.x_s
        barrier3d.x_b_TS[-1] = (
            barrier3d.x_s
            + barrier3d.InteriorWidth_AvgTS[-1]
            + np.size(barrier3d.DuneDomain, 2)  # dune domain width in dam
            + self._beach_width[self._time_index - 1] / 10  # in dam
        )

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
