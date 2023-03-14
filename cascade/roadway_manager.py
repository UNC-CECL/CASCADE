"""Manage roadways

This module provides functions for modifying a barrier segment from Barrier3D -- consisting of 1+ rows of dune cells, a
separate interior grid, and an idealized shoreface -- for roadway management decisions, including:
    1) overwash removal from the roadway after storms and placement on the dune line,
    2) road relocation landward when the dunes migrate over the roadway,
    3) dune rebuilding when the dunes fall below a minimum height.

References
----------

.. [1] Velasquez-Montoya, L., Sciaudone, E. J., Smyre, E., & Overton, M. F. (2021). Vulnerability Indicators for
    Coastal Roadways Based on Barrier Island Morphology and Shoreline Change Predictions. Natural Hazards Review, 22(2),
    04021003.
.. [2] Sciaudone, E. J., Velasquez-Montoya, L., Smyre, E. A., & Overton, M. F. (2016). Pea Island, North Carolina. Shore
    & Beach, 84(2), 10.
.. [3] Vinent, O. D., & Moore, L. J. (2015). Barrier island bistability induced by biophysical interactions.
    Nature Climate Change, 5(2), 158-162.

Notes
---------
The alongshore length of a barrier segment in Barrier3D is time-invariant, whereas the barrier interior
width -- and number of cross-shore cells -- varies dynamically due to storm impacts and SLR.

Because SLR is simulated using a Lagrangian reference frame in Barrier3D, the roadway and dune elevations are reduced
by SLR for each time step.

"""
import copy

import numpy as np
from scipy.interpolate import interp2d

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


def bulldoze(
    time_index,
    xyz_interior_grid,
    yxz_dune_grid,
    road_ele=2.0,
    road_width=20,
    road_setback=30,
    dx=10,
    dy=10,
    dz=10,
    drown_threshold=0,
    percent_water_cells_touching_road=0.2,
):
    r"""
    Remove overwash from roadway and put it back on the adjacent dune. Spreads sand evenly across adjacent dune cells.

    Check for width drowning of roadway (i.e., when a water cell touches the roadway).

    Parameters
    ----------
    time_index: int,
        Time index for drowning error message
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; for Barrier3d, dz=10, decameters MHW]
    yxz_dune_grid:
        Dune topography [z units specified by dz; for Barrier3d, dz=10, decameters above the berm elevation]
    road_ele: float
        Road elevation [m; needs to be in same reference frame as xyz; for Barrier3d, decameters MHW]
    road_width: int
        Width of roadway [m]
    road_setback: int
        Setback distance of roadway from edge of interior domain [m]
    dx: int
        Cross-shore discretization of x [default is dx=10, dam]
    dy: int
        Alongshore discretization of y [default is dy=10, dam]
    dz: int
        Vertical discretization of z [default is dz=10, dam]
    drown_threshold: float
        Elevation threshold for roadway drowning [m; needs to be in same reference frame as xyz]
    percent_water_cells_touching_road: float
        Fraction of cells below drown_threshold

    Returns
    -------
    array of float
        new_road_domain: in units of dx, dy, dz
        new_dune_domain: in units of dx, dy, dz
        overwash removal: in units of dx*dy*dz (for Barrier3D, dam^3)
    int
        roadway_drown: flag for if water cells border the road on either side

    """

    # road parameters
    road_start = int(
        road_setback / dy
    )  # grid index for start of roadway in interior domain (for B3D, convert to dam)
    road_width = int(road_width / dx)
    road_end = road_start + road_width
    road_ele = (
        road_ele / dz
    )  # convert to units of grid (NOTE: in B3D default simulation, berm is 1.44 m MHW)

    # remove sand from roadway (only account for positive values)
    old_road_domain = xyz_interior_grid[road_start:road_end, :]
    new_road_domain = np.zeros((road_width, np.size(old_road_domain, 1))) + road_ele
    road_overwash_removal = sum(old_road_domain - new_road_domain)
    road_overwash_removal[road_overwash_removal < 0] = 0

    # spread overwash removed from roadway equally over the adjacent dune cells
    number_dune_cells = np.size(yxz_dune_grid, 1)
    overwash_to_dune = np.transpose(
        [road_overwash_removal / number_dune_cells] * number_dune_cells
    )
    new_dune_domain = yxz_dune_grid + overwash_to_dune

    xyz_interior_grid[
        road_start:road_end, :
    ] = new_road_domain  # update interior domain

    # check if any water cells border the road on either side
    number_border_cells = np.size(xyz_interior_grid[road_end, :])
    bayside_water_cells = (
        np.count_nonzero((xyz_interior_grid[road_end + 1, :] * dz) <= drown_threshold)
        / number_border_cells
    )

    if road_start > 0:
        seaside_water_cells = (
            np.count_nonzero(
                (xyz_interior_grid[road_start - 1, :] * dz) <= drown_threshold
            )
            / number_border_cells
        )
    else:
        seaside_water_cells = 0

    # for debugging
    # if seaside_water_cells > 0:
    #     print(
    #         " ALERT: {water_cells}% of seaside road borders water".format(
    #             water_cells=np.array(seaside_water_cells) * 100
    #         )
    #     )
    # if bayside_water_cells > 0:
    #     print(
    #         " ALERT: {water_cells}% of bayside road borders water".format(
    #             water_cells=np.array(bayside_water_cells) * 100
    #         )
    #     )

    if (seaside_water_cells > percent_water_cells_touching_road) or (
        bayside_water_cells > percent_water_cells_touching_road
    ):
        roadway_drown = True
        print(
            "Roadway width drowned at {time} years, {water}% of road borders water".format(
                time=time_index - 1, water=percent_water_cells_touching_road * 100
            )
        )
    else:
        roadway_drown = False

    return (
        new_dune_domain,
        xyz_interior_grid,
        np.sum(road_overwash_removal),
        roadway_drown,
    )


def rebuild_dunes(
    yxz_dune_grid, max_dune_height=3.0, min_dune_height=2.4, dz=10, rng=True
):
    r"""Rebuild dunes if they fall below a minimum elevation after a storm. Note that while artificial dune geometry
    of a given height is constrained by the angle of repose of beach, here we just assume the dunes are built to a width
    capable of maintaining dunes of the stated height (default is two grid cells in Barrier3D = 20 m).

    If the min and max dune heights differ, a linear gradient is applied from the first to last dune row, with small
    random perturbations.

    From Valasquez et al., (2020), the average elevation of the road along NC-12 is 1.3 m (NAVD88); they find that in
    order for the road to not be vulnerable to overwash, the dune crest must be higher than 4.3 m (NAVD88), so here the
    default max_dune_height is 3 m. Note that dune height in Barrier3D is measured as the height above the dune toe
    (berm elevation).

    Parameters
    ----------
    yxz_dune_grid:
        Dune topography [z units specified by dz; for Barrier3D, dz=10, decameters above the berm elevation]
    max_dune_height: float
        Maximum dune height for dune rebuilding [m]
    min_dune_height: float
        Minimum dune height for dune rebuilding [m]
    dz: int
        Vertical discretization of z [default is dz=10, dam]
    rng: boolean
        If True, add random perturbations alongshore to dune height

    Returns
    -------
    new_dune_domain: float or array of float
        New yxz dune domain in new_road_domain: in units of dx, dy, dz
    rebuild_dune_volume: float
        Volume of sand for dune rebuild, in units of dx*dy*dz

    """

    # old_dune_domain = copy.deepcopy(yxz_dune_grid)
    old_dune_domain = yxz_dune_grid

    # convert from m to grid z discretization
    max_height = max_dune_height / dz
    min_height = min_dune_height / dz
    ny = np.size(yxz_dune_grid, 0)
    nx = np.size(yxz_dune_grid, 1)

    if rng:
        # add some random perturbations to dune heights
        RNG = np.random.default_rng(seed=1973)

        dune_start_max = np.ones([ny]) * (
            max_height + (-0.01 + (0.01 - (-0.01)) * RNG.random(ny))
        )
        dune_start_min = np.ones([ny]) * (
            min_height + (-0.01 + (0.01 - (-0.01)) * RNG.random(ny))
        )
    else:
        dune_start_max = np.ones([ny]) * max_height
        dune_start_min = np.ones([ny]) * min_height

    # linearly interpolate from front row (max height) to back row (min height)
    x = [0, nx - 1]
    y = [np.arange(0, ny, 1)]
    z = np.transpose([dune_start_max, dune_start_min])
    f = interp2d(x, y, z)
    new_dune_domain = f(np.arange(0, nx, 1), np.arange(0, ny, 1))
    rebuild_dune_volume = np.sum(new_dune_domain - old_dune_domain)

    return new_dune_domain, rebuild_dune_volume


def set_growth_parameters(
    yxz_dune_grid,
    Dmax,
    growthparam,
    original_growth_param=None,
    rmin=0.35,
    rmax=0.85,
):
    r"""Set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the natural eq.
    dune height (Dmax).

    We understand from modeling work (Duran and Moore, 2013) and from empirical evidence (Houser et al., 2015)
    that dunes can reach a maximum height due to negative feedacks between the (cross-shore) wind field and the dune.
    For non-normal wind incidence, it has been suggested that dunes may continue to grow in height, albeit
    at a very slow rate (Davidson Arnott et al., 2018).

    Parameters
    ----------
    yxz_dune_grid:
        Dune topography [units must be the same as Dmax]
    Dmax: float
        Maximum natural equilibrium dune height [default in Barrier3D is decameters]
    growthparam:
        growth parameters from last time step, used in Houser formulation [unitless]
    original_growth_param: float, optional
        growth parameters from first time step, before humans interfered [unitless]
    rmin: float, optional
        Minimum growth rate - used if original_growth_parm not provided [unitless]
    rmax: float, optional
        Maximum growth rate - used if original_growth_parm not provided [unitless]

    Returns
    -------
    new_growth_param: array of float
        New growth parameter array that accounts for human modifications to dune height above/below the equilibrium

    """

    ny = np.size(yxz_dune_grid, 0)
    new_growth_param = np.copy(growthparam)
    rng = np.random.default_rng(seed=1973)

    for idune in range(ny):
        if yxz_dune_grid[idune, 0] > Dmax:  # if dune height above dmax, don't grow
            new_growth_param[0, idune] = 0

        else:
            # if dune is now below the Dmax (was formerly above), make sure it has a growth rate either the same as
            # before (if original growth rate provided) or a random number between rmin and rmax
            if growthparam[0, idune] == 0:
                if original_growth_param is not None:
                    new_growth_param[0, idune] = original_growth_param[0, idune]
                else:
                    new_growth_param[0, idune] = rmin + (rmax - rmin) * rng.random()

    return new_growth_param


def get_road_relocation_elevation(
    time_index,
    xyz_interior_grid,
    road_setback,
    road_width,
    dx=10,
    dy=10,
    dz=10,
):
    r"""
    For a given setback distance, check what the road elevation would be if relocated at grade.

    If zero MSL, the roadway can't be relocated (height drowned).

    Parameters
    ----------
    time_index: int,
        Time index for drowning error message
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; for Barrier3d, dz=10, decameters MHW]
    road_width: int
        Width of roadway [m]
    road_setback: int
        Setback distance of roadway from edge of interior domain [m]
    dx: int
        Cross-shore discretization of x [default is dx=10, dam]
    dy: int
        Alongshore discretization of y [default is dy=10, dam]
    dz: int
        Vertical discretization of z [default is dz=10, dam]

    Returns
    -------
    array of float
        road_ele: in units of dz
    bool
        roadway_drown: if the roadway elevation would be zero MSL, the roadway can't be constructed or relocated

    """
    roadway_drown = 0

    # road parameters
    road_start = int(
        road_setback / dy
    )  # grid index for start of roadway in interior domain (for B3D, convert to dam)
    road_width = int(road_width / dx)
    road_end = road_start + road_width

    # calculate average elevation
    road_domain = xyz_interior_grid[road_start:road_end, :]
    road_ele = np.mean(road_domain) * dz  # dam to meters

    # if the roadway elevation would be zero MSL, the roadway can't be constructed or relocated
    if road_ele <= 0:
        roadway_drown = 1
        print(
            "Roadway cannot be relocated at {time} years b/c the road would be at or below MSL".format(
                time=time_index - 1
            )
        )

    return road_ele, roadway_drown


def road_relocation_checks(
    time_index,
    dune_migrated,
    road_setback,
    road_relocation_setback,
    road_relocation_width,
    average_barrier_width,
):
    r"""Check if the roadway needs to be relocated due to dune migration, and if there is room for roadway relocation.

    Parameters
    ----------
    time_index: int
        Time index for drowning error message
    dune_migrated: int
        Number of meters dune migrated in Barrier3D; + if progrades, - if erodes [m]
    road_setback: float
        The setback distance of roadway from edge of interior domain from the last time step [m]
    road_relocation_setback: float
        The setback distance specified for roadway relocation [m]
    road_relocation_width: float
        The road width specified for roadway relocation [m]
    average_barrier_width: float
        The average barrier width from the last time step [m]

    Returns
    -------
    bool
        road_relocated: roadway was relocated due to dune migration
        relocation_break: no room for relocation of the roadway
    float
        road_setback: updated setback distance for dune migration and road relocation

    """

    # initialize the break booleans as False
    relocation_break = 0
    road_relocated = 0

    # if dune line eroded or prograded, subtract (erode) or add (prograde) to the setback to keep road in the same place
    if dune_migrated != 0:
        road_setback = road_setback + dune_migrated

        # with this shoreline change and dune migration, check if the roadway needs to be relocated
        if road_setback < 0:
            road_relocated = 1

            # relocate the road only if the width of the island allows it
            if (
                road_relocation_setback
                + (
                    2 * road_relocation_width
                )  # bay shoreline buffer of one roadway width
                > average_barrier_width
            ):
                relocation_break = 1
                print(
                    "Island is too narrow for roadway to be relocated. Roadway eaten up by dunes at {time} years".format(
                        time=time_index - 1
                    )  # -1 because B3D advances time step at end of dune_update
                )
            else:
                road_setback = road_relocation_setback

    return road_relocated, road_setback, relocation_break


class RoadwayManager:
    """Manage the road!

    Examples
    --------
    # >>> from cascade.roadway_manager import RoadwayManager
    # >>> roadways = RoadwayManager()
    # >>> roadways.update(barrier3d, trigger_dune_knockdown)
    """

    def __init__(
        self,
        initial_road_elevation=1.7,
        road_width=30,
        road_setback=30,
        initial_dune_design_elevation=3.7,
        initial_dune_minimum_elevation=2.2,
        time_step_count=500,
        original_growth_param=None,
    ):
        """The RoadwayManager module

        Parameters
        ----------
        initial_road_elevation: float, optional
            Initial elevation of the roadway [m MHW]; road relocation elevations are built to grade based on setback.
        road_width: int, optional
            Width of roadway [m]. Also used for relocation width.
        road_setback: int, optional
            Setback of roadway from the inital dune line [m]. Also used for relocation setback.
        initial_dune_design_elevation: float, optional
            Elevation which dune is originally rebuilt to when road established [m MHW]. Used for rebuild dune height.
        initial_dune_minimum_elevation: float, optional
            Elevation threshold which originally triggers rebuilding of dune [m MHW]. Used for min dune height.
        time_step_count: int, optional
            Number of time steps.
        original_growth_param: optional
            Dune growth parameters from first time step of Barrier3d, before human modifications [unitless]

        """

        self._road_width = road_width
        self._road_setback = road_setback
        self._road_ele = initial_road_elevation
        self._dune_design_elevation = initial_dune_design_elevation  # can be `None` if user doesn't want to rebuild
        self._dune_minimum_elevation = initial_dune_minimum_elevation
        self._original_growth_param = original_growth_param
        self._nt = time_step_count
        self._drown_break = 0
        self._relocation_break = 0
        self._time_index = 1
        self._absolute_minimum_dune_height = 0.3
        self._percent_water_cells_touching_road = 0.2

        # set relocation parameters to original values
        self._road_relocation_width = (
            road_width  # can be updated outside `update` within cascade
        )
        self._road_relocation_setback = (
            road_setback  # can be updated outside `update` within cascade
        )

        # user can specify that dune rebuilding is off with `None`: mostly for debugging and sensitivity testing
        if (
            self._dune_design_elevation is not None
            and self._dune_minimum_elevation is not None
        ):
            self._relocation_dune_design_height_above_road = (
                self._dune_design_elevation - self._road_ele
            )  # m
            self._relocation_dune_minimum_height_above_road = (
                self._dune_minimum_elevation - self._road_ele
            )  # m

        # time series
        self._road_setback_TS = np.zeros(
            self._nt
        )  # changes with time w/dune migration and user input
        self._road_setback_TS[0] = road_setback
        self._road_width_TS = np.zeros(
            self._nt
        )  # could change with time due to user input with road relocation
        self._road_width_TS[0] = road_width
        self._road_ele_TS = np.zeros(
            self._nt
        )  # changes with time due to lagrangian SLR and user input
        self._road_ele_TS[0] = self._road_ele
        self._dune_design_elevation_TS = np.zeros(
            self._nt
        )  # changes with time due to lagrangian SLR and user input
        self._dune_design_elevation_TS[0] = self._dune_design_elevation
        self._dune_minimum_elevation_TS = np.zeros(
            self._nt
        )  # changes with time due to lagrangian SLR and user input
        self._dune_minimum_elevation_TS[0] = self._dune_minimum_elevation
        self._dunes_rebuilt_TS = np.zeros(self._nt)  # when dunes are rebuilt (boolean)
        self._road_relocated_TS = np.zeros(self._nt)  # when road is relocated (boolean)
        self._rebuild_dune_volume_TS = np.zeros(
            self._nt
        )  # sand for rebuilding dunes [m^3]
        self._road_overwash_volume = np.zeros(
            self._nt
        )  # total overwash removed from roadway [m^3]
        self._percent_below_min = [
            None
        ] * self._nt  # keep track of what percent of the dune elevations fall below minimum threshold
        self._growth_params = [
            np.nan
        ] * self._nt  # when dune growth parameters set to zero b/c of rebuild height
        self._growth_params[0] = original_growth_param

        # also keep track of post-storm dune and interior impacts before human modifications
        self._post_storm_dunes = [None] * self._nt
        self._post_storm_interior = [None] * self._nt
        self._post_storm_ave_interior_height = [None] * self._nt

    def update(self, barrier3d, trigger_dune_knockdown):
        self._time_index = barrier3d.time_index

        if self._original_growth_param is None:
            self._original_growth_param = barrier3d.growthparam

        # save post-storm dune and interior domain before human modifications (essentially a 0.5 yr time step)
        self._post_storm_interior[self._time_index - 1] = copy.deepcopy(
            barrier3d.InteriorDomain
        )
        self._post_storm_dunes[self._time_index - 1] = copy.deepcopy(
            barrier3d.DuneDomain[self._time_index - 1, :, :]
        )
        self._post_storm_ave_interior_height[self._time_index - 1] = copy.deepcopy(
            barrier3d.h_b_TS[-1]
        )

        ###############################################################################
        # roadway checks for relocation, drowning; update for SLR
        ###############################################################################

        # check if road relocation is needed
        average_barrier_width = barrier3d.InteriorWidth_AvgTS[-1] * 10  # m
        dune_migration = (
            barrier3d.ShorelineChangeTS[self._time_index - 1] * 10
        )  # if +, dune progrades; -, dune erodes into interior [m]
        [
            road_relocated,
            self._road_setback,
            self._relocation_break,
        ] = road_relocation_checks(
            self._time_index,
            dune_migration,
            self._road_setback,  # current road setback, m
            self._road_relocation_setback,  # setback specified for relocation, m
            self._road_relocation_width,  # width specified for relocation, m
            average_barrier_width,  # current width, m
        )

        # if road can't be relocated, no longer manage and exit; dune growth parameters reset to original in CASCADE
        if self._relocation_break == 1:
            # an adaptation solution may be to knock down the dunes so that they are small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return

        # if road is relocated, get the new road elevation (built at grade) and update dune elevations which are
        # dependent on the road elevation; otherwise, decrease all elevations (m MHW) this year by the SLR increment
        if road_relocated:
            self._road_width = self._road_relocation_width
            self._road_ele, self._drown_break = get_road_relocation_elevation(
                self._time_index,
                xyz_interior_grid=barrier3d.InteriorDomain,  # interior domain from this last time step, dam
                road_setback=self._road_setback,  # m
                road_width=self._road_width,  # m
                dx=10,
                dy=10,
                dz=10,  # specifies interior is in dam
            )

            # user can specify that dune rebuilding is off with `None`
            if (
                self._dune_design_elevation is not None
                or self._dune_minimum_elevation is not None
            ):
                self._dune_design_elevation = (
                    self._road_ele + self._relocation_dune_design_height_above_road
                )
                self._dune_minimum_elevation = (
                    self._road_ele + self._relocation_dune_minimum_height_above_road
                )

        else:
            self._road_ele = self._road_ele - (
                barrier3d.RSLR[self._time_index - 1] * 10
            )  # m MHW
            if (
                self._dune_design_elevation is not None
                or self._dune_minimum_elevation is not None
            ):
                self._dune_design_elevation = self._dune_design_elevation - (
                    barrier3d.RSLR[self._time_index - 1] * 10
                )
                self._dune_minimum_elevation = self._dune_minimum_elevation - (
                    barrier3d.RSLR[self._time_index - 1] * 10
                )

        # road cannot be below 0 m MHW (sea level); stop managing!
        if self._road_ele < 0:
            self._drown_break = 1
            print(
                "Roadway drowned in place at {time} years due to SLR - road cannot be below 0 m MHW".format(
                    time=self._time_index - 1
                )
            )

            # an adaptation solution may be to knock down the dunes so that they are small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return
        elif self._drown_break == 1:  # if road drowned from road relocation above
            # an adaptation solution may be to knock down the dunes so that they are small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return

        # when the roadway gets really low in elevation, the dune_design_elevation may not be above the berm; when this
        # happens, we use a design height of 1 m above the berm to keep a dune to protect the roadway and rebuild
        # whenever the dune drops to just above elevation of the berm (0.3 m) -- essentially, we just push the sand back
        self._dune_design_elevation = max(
            self._dune_design_elevation, (barrier3d.BermEl * 10) + 1.0
        )
        # note, Barrier3D adds a random seeded height for the proto/new dune line (20 cm); here we allow the user to
        # specify
        self._dune_minimum_elevation = max(
            self._dune_minimum_elevation,
            (barrier3d.BermEl * 10) + self._absolute_minimum_dune_height,
        )

        # save time series
        self._road_setback_TS[self._time_index - 1] = self._road_setback
        self._road_width_TS[self._time_index - 1] = self._road_width
        self._road_ele_TS[self._time_index - 1] = self._road_ele
        self._dune_design_elevation_TS[
            self._time_index - 1
        ] = self._dune_design_elevation
        self._dune_minimum_elevation_TS[
            self._time_index - 1
        ] = self._dune_minimum_elevation
        self._road_relocated_TS[self._time_index - 1] = road_relocated

        ###############################################################################
        # bulldoze roadway after storms and check for road width drowning
        ###############################################################################

        # bulldoze the road and put bulldozed sand back on the dunes; drown road when a water cell touches either side
        (
            new_dune_domain,  # all in dam
            new_xyz_interior_domain,
            road_overwash_removal,
            self._drown_break,
        ) = bulldoze(
            time_index=self._time_index,
            road_ele=self._road_ele,  # m MHW
            road_width=self._road_width,  # m
            road_setback=self._road_setback,  # m
            xyz_interior_grid=barrier3d.InteriorDomain,  # interior domain from this last time step, dam
            yxz_dune_grid=barrier3d.DuneDomain[
                self._time_index - 1, :, :
            ],  # dune domain from this last time step, dam
            dx=10,
            dy=10,
            dz=10,  # specifies dam for dune and interior domains
            drown_threshold=0,  # 0 m MSL
            percent_water_cells_touching_road=self._percent_water_cells_touching_road,  # fraction cells<drown_threshold
        )
        if self._drown_break == 1:
            # an adaptation solution may be to knock down the dunes so that they are small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return  # exit program

        self._road_overwash_volume[self._time_index - 1] = (
            road_overwash_removal * dm3_to_m3
        )  # convert from dam^3 to m^3

        # update Barrier3D class variables
        new_ave_interior_height = np.average(
            new_xyz_interior_domain[
                new_xyz_interior_domain >= barrier3d.SL
            ]  # all in dam MHW
        )
        barrier3d.h_b_TS[
            -1
        ] = new_ave_interior_height  # slightly altered due to roadway
        barrier3d.InteriorDomain = new_xyz_interior_domain
        barrier3d.DomainTS[self._time_index - 1] = new_xyz_interior_domain

        ###############################################################################
        # Rebuild dunes if knocked down
        ###############################################################################

        # dune management: rebuild dunes!
        if self._dune_design_elevation is None or self._dune_minimum_elevation is None:
            pass
        else:
            # in B3D, dune height is the height above the berm crest
            dune_design_height = self._dune_design_elevation - (barrier3d.BermEl * 10)
            min_dune_height = self._dune_minimum_elevation - (barrier3d.BermEl * 10)

            # if any dune cell in the front row of dunes is less than a minimum threshold -- as measured above the
            # berm crest -- then rebuild the dune (all rows up to dune_design_elevation)
            if np.min(new_dune_domain[:, 0]) < (min_dune_height / 10):  # in dam
                # first document what percentage of the dune field is below this minimum
                dune_cells_below_threshold = np.sum(
                    new_dune_domain < (min_dune_height / 10)
                )
                self._percent_below_min[self._time_index - 1] = (
                    dune_cells_below_threshold / np.size(new_dune_domain) * 100
                )

                new_dune_domain, rebuild_dune_volume = rebuild_dunes(
                    new_dune_domain,  # dam
                    max_dune_height=dune_design_height,  # in m
                    min_dune_height=dune_design_height,  # in m
                    dz=10,  # specifies dune domain is in dam
                    rng=True,  # adds stochasticity to dune height (seeded)
                )
                self._dunes_rebuilt_TS[self._time_index - 1] = 1
                self._rebuild_dune_volume_TS[self._time_index - 1] = (
                    rebuild_dune_volume * dm3_to_m3
                )

        # update Barrier3D class variables
        barrier3d.DuneDomain[self._time_index - 1, :, :] = new_dune_domain

        ###############################################################################
        # Set dune growth rate to zero if dunes rebuilt > natural equilibrium height
        ###############################################################################

        # set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the
        # natural eq. dune height (Dmax)
        new_growth_parameters = set_growth_parameters(
            new_dune_domain,  # in dam
            barrier3d.Dmax,  # in dam
            barrier3d.growthparam,
            original_growth_param=self._original_growth_param,  # use original growth rates for resetting values
        )
        self._growth_params[self._time_index - 1] = copy.deepcopy(new_growth_parameters)

        # update class variables
        barrier3d.growthparam = new_growth_parameters

        return

    @property
    def road_relocation_width(self):
        return self._road_relocation_width

    @road_relocation_width.setter
    def road_relocation_width(self, value):
        self._road_relocation_width = value

    @property
    def road_relocation_setback(self):
        return self._road_relocation_setback

    @road_relocation_setback.setter
    def road_relocation_setback(self, value):
        self._road_relocation_setback = value

    @property
    def drown_break(self):
        return self._drown_break

    @property
    def time_index(self):
        return self._time_index

    @drown_break.setter
    def drown_break(self, value):
        self._drown_break = value

    @property
    def relocation_break(self):
        return self._relocation_break

    @relocation_break.setter
    def relocation_break(self, value):
        self._relocation_break = value

    @property
    def percent_water_cells_touching_road(self):
        return self._percent_water_cells_touching_road

    @percent_water_cells_touching_road.setter
    def percent_water_cells_touching_road(self, value):
        self._percent_water_cells_touching_road = value
