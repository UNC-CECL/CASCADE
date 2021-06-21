"""Beach Nourisher

This module provides functions ....  ADD DESCRIPTION HERE ONCE FINALIZED

References
----------

.. [1]


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
    r"""ADD DESCRIPTION HERE
    Note, all parameters below must be in the same units. For use with Barrier3D, decameters.

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
        New shoreline position after nourishment (rounded to the nearest decameter)
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


class BeachNourisher:
    """Nourish the beach

    Examples
    --------
    >>> from cascade.beach_nourisher import BeachNourisher
    >>> nourish = BeachNourisher()
    >>> nourish.update(barrier3d=, nourish_now=, nourishment_interval=)
    """

    def __init__(
        self,
        nourishment_interval=None,
        nourishment_volume=100,
        initial_beach_width=10,
        dune_design_elevation=3.7,
        time_step_count=500,
        original_growth_param=None,
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
        """

        self._nourishment_volume = nourishment_volume
        self._nourishment_interval = nourishment_interval
        self._nourishment_counter = nourishment_interval
        self._beach_width_threshold = 10  # m
        self._dune_design_elevation = dune_design_elevation
        self._original_growth_param = original_growth_param
        self._nt = time_step_count
        self._drown_break = 0
        self._narrow_break = 0
        self._time_index = 1

        # time series
        self._beach_width = [None] * self._nt
        self._beach_width[0] = initial_beach_width  # m
        self._nourishment_TS = np.zeros(self._nt)
        self._dunes_rebuilt_TS = np.zeros(self._nt)
        self._nourishment_volume_TS = np.zeros(self._nt)  # m^3/m
        self._rebuild_dune_volume_TS = np.zeros(self._nt)   # m^3
        self._growth_params = [
            None
        ] * self._nt  # track when dune growth parameters set to zero b/c of rebuild height
        self._growth_params[0] = original_growth_param

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

        # reduce beach width by the amount of shoreline change from last time step
        self._beach_width[self._time_index - 1] = (
            self._beach_width[self._time_index - 2]
            - (barrier3d.x_s_TS[-1] - barrier3d.x_s_TS[-2]) * 10  # m
        )

        # if the beach width is less than 10 m, turn dune migration in Barrier3D back on
        if self._beach_width[self._time_index - 1] <= self._beach_width_threshold:
            barrier3d.dune_migration_on = True

        # update counter and check for nourishment and dune rebuild flags
        if self._nourishment_counter is not None:
            self._nourishment_counter -= 1

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
            self._rebuild_dune_volume_TS[self._time_index - 1] = rebuild_dune_volume * dm3_to_m3

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