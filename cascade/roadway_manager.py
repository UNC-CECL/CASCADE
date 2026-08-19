"""Manage roadways

This module provides functions for modifying a barrier segment from Barrier3D --
consisting of 1+ rows of dune cells, a separate interior grid, and an idealized
shoreface -- for roadway management decisions, including:

1. overwash removal from the roadway after storms and placement on the dune line,
2. road relocation landward when the dunes migrate over the roadway,
3. dune rebuilding when the dunes fall below a minimum height,
4. optional shoreface/beach nourishment without community dune management.

References
----------

.. [1] Velasquez-Montoya, L., Sciaudone, E. J., Smyre, E., & Overton, M. F. (2021).
       Vulnerability Indicators for Coastal Roadways Based on Barrier Island
       Morphology and Shoreline Change Predictions. Natural Hazards Review, 22(2),
       04021003.
.. [2] Sciaudone, E. J., Velasquez-Montoya, L., Smyre, E. A., & Overton, M. F. (2016).
       Pea Island, North Carolina. Shore & Beach, 84(2), 10.
.. [3] Vinent, O. D., & Moore, L. J. (2015). Barrier island bistability induced
       by biophysical interactions. Nature Climate Change, 5(2), 158-162.

Notes
-----
The alongshore length of a barrier segment in Barrier3D is time-invariant, whereas
the barrier interior width -- and number of cross-shore cells -- varies dynamically
due to storm impacts and SLR.

Because SLR is simulated using a Lagrangian reference frame in Barrier3D, the
roadway and dune elevations are reduced by SLR for each time step.
"""

import copy
import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator

dm3_to_m3 = 1000  # convert from cubic decameters to cubic meters


def shoreface_nourishment(
    x_s,
    x_t,
    nourishment_volume,
    average_barrier_height,
    shoreface_depth,
    beach_width,
):
    """Apply sand to the shoreface without managing or rebuilding the dunes.

    This is the same shoreface-volume formulation used by BeachDuneManager. It is
    included locally so RoadwayManager can nourish a roadway-managed domain without
    turning on the community-management module (which would also filter overwash,
    rebuild community dunes, and control dune migration).

    All inputs must use one consistent length system. RoadwayManager calls this
    function with decameters and ``dam^3/dam`` so it remains consistent with
    Barrier3D.

    Parameters
    ----------
    x_s: float
        Current shoreline position.
    x_t: float
        Current shoreface-toe position.
    nourishment_volume: float
        Nourishment volume per unit alongshore length.
    average_barrier_height: float
        Current average barrier height.
    shoreface_depth: float
        Current shoreface depth.
    beach_width: float
        Current beach width.

    Returns
    -------
    new_x_s: float
        Shoreline position after nourishment.
    new_shoreface_slope: float
        Shoreface slope after nourishment.
    new_beach_width: float
        Beach width after nourishment.
    """

    new_x_s = x_s - (2 * nourishment_volume) / (
        2 * average_barrier_height + shoreface_depth
    )
    new_shoreface_slope = shoreface_depth / (new_x_s - x_t)
    new_beach_width = beach_width + (x_s - new_x_s)

    return new_x_s, new_shoreface_slope, new_beach_width


def beach_width_dune_dynamics(
    current_beach_width,
    beach_width_last_year,
    beach_width_threshold,
    barrier3d,
    time_index,
):
    """Apply the BeachDuneManager beach-width/dune-migration rule.

    Dune migration remains off while a positive managed beach separates the
    shoreline from the dune line. Once beach width reaches the threshold (0 m),
    Barrier3D dune migration turns back on. If erosion has already consumed one or
    more complete 10-m cells, the missed dune migration is applied immediately.
    """

    if current_beach_width <= beach_width_threshold:
        barrier3d.dune_migration_on = True

        if beach_width_last_year > beach_width_threshold:
            cellular_shoreline_change = int(
                math.floor(current_beach_width) / -10
            )

            if cellular_shoreline_change >= 1:
                barrier3d.SCRagg[time_index - 1] = (
                    (barrier3d.x_s_TS[-1] % 1) + cellular_shoreline_change
                ) * -1
                barrier3d.migrate_dunes(
                    shoreline_change_aggregate=barrier3d.SCRagg[time_index - 1],
                    cellular_shoreline_change=cellular_shoreline_change,
                    time_index=time_index - 1,
                )

                (
                    _,
                    _,
                    interior_width_average,
                ) = barrier3d.FindWidths(barrier3d.InteriorDomain, barrier3d.SL)
                barrier3d.InteriorWidth_AvgTS[-1] = interior_width_average
                barrier3d.DomainTS[time_index - 1] = barrier3d.InteriorDomain
            else:
                barrier3d.SCRagg[time_index - 1] = (barrier3d.x_s % 1) * -1

        if current_beach_width < 0:
            current_beach_width = 0
    else:
        barrier3d.dune_migration_on = False

    return current_beach_width


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
    allow_causeway=False,
):
    """
    Remove overwash from roadway and put it back on the adjacent dune. Spreads sand
    evenly across adjacent dune cells.

    Check for width drowning of roadway (i.e., when a water cell touches the roadway).

    Parameters
    ----------
    time_index: int,
        Time index for drowning error message
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; for Barrier3d,
        dz=10, decameters MHW]
    yxz_dune_grid:
        Dune topography [z units specified by dz; for Barrier3d, dz=10, decameters
        above the berm elevation]
    road_ele: float
        Road elevation [m; needs to be in same reference frame as xyz; for Barrier3d,
        decameters MHW]
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
        Elevation threshold for roadway drowning [m; needs to be in same reference
        frame as xyz]
    percent_water_cells_touching_road: float
        Fraction of cells below drown_threshold
    allow_causway: bool
        Whether roadways drowns when surrounded by water [default is allow_causeway=FALSE]

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

    xyz_interior_grid[road_start:road_end, :] = (
        new_road_domain  # update interior domain
    )

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

    if ((seaside_water_cells > percent_water_cells_touching_road) or (
        bayside_water_cells > percent_water_cells_touching_road)) and  allow_causeway==False:
        roadway_drown = True
        # print(
        #     f"Roadway width drowned at {time_index - 1} years, "
        #     f"{percent_water_cells_touching_road * 100.0}% of road borders water"
        # )
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
    """Rebuild dunes if they fall below a minimum elevation after a storm. Note that
    while artificial dune geometry of a given height is constrained by the angle of
    repose of beach, here we just assume the dunes are built to a width capable of
    maintaining dunes of the stated height (default is two grid cells in
    Barrier3D = 20 m).

    If the min and max dune heights differ, a linear gradient is applied from the
    first to last dune row, with small random perturbations.

    From Valasquez et al., (2020), the average elevation of the road along NC-12
    is 1.3 m (NAVD88); they find that in order for the road to not be vulnerable to
    overwash, the dune crest must be higher than 4.3 m (NAVD88), so here the
    default max_dune_height is 3 m. Note that dune height in Barrier3D is measured
    as the height above the dune toe (berm elevation).

    Parameters
    ----------
    yxz_dune_grid: ndarray, shape (ny, nx)
        Dune topography [z units specified by dz; for Barrier3D, dz=10, decameters
        above the berm elevation]
    max_dune_height: float, optional
        Maximum dune height for dune rebuilding [m]
    min_dune_height: float, optional
        Minimum dune height for dune rebuilding [m]
    dz: int, optional
        Vertical discretization of z [default is dz=10, dam]
    rng: bool or np.random.Generator, optional
        If `True`, add random perturbations alongshore to dune height. `rng`
        can also be an object that provides a `uniform` method, like `numpy.random`
        or `numpy.random.Generator`.

    Returns
    -------
    new_dune_domain: ndarray, shape (ny, nx)
        New yxz dune domain in new_road_domain: in units of dx, dy, dz
    rebuild_dune_volume: float
        Volume of sand for dune rebuild, in units of dx*dy*dz
    """
    if rng and isinstance(rng, bool):
        rng = np.random.default_rng(seed=1973)

    ny, nx = yxz_dune_grid.shape

    # convert from m to grid z discretization
    dune_height = np.empty((ny, 2), dtype=float)
    dune_height[:, 0] = max_dune_height / dz
    dune_height[:, 1] = min_dune_height / dz

    # add some random perturbations to dune heights
    if rng:
        dune_height += rng.uniform(high=0.01, size=dune_height.size).reshape((ny, 2))

    interpolate = RegularGridInterpolator((np.arange(ny), [0, nx - 1]), dune_height)

    new_dune_domain = interpolate(
        tuple(np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij"))
    )

    rebuild_dune_volume = np.sum(new_dune_domain - yxz_dune_grid)
    z = 20
    return new_dune_domain, rebuild_dune_volume

def build_interior_dunes(
    b3d, dune_construction_distance =0, max_dune_height=3.0, min_dune_height=2.4, dz=10, rng=True
):
    """Build dunes within the barrier island interior if B3D dunes are below specific minimum elevation.
    Constructing dunes within barrier island interior is based on NCDOT management practices on Ocracoke, as
    the NCDOT was unable to construct dunes in areas outside a 75' right of way on Ocracoke. Here we construct
    a dune line within the first row of the barrier island interior grid that falls within the NCDOT's 75' right of way

    If the min and max dune heights differ, a linear gradient is applied from the
    first to last dune row, with small random perturbations.

    From Valasquez et al., (2020), the average elevation of the road along NC-12
    is 1.3 m (NAVD88); they find that in order for the road to not be vulnerable to
    overwash, the dune crest must be higher than 4.3 m (NAVD88), so here the
    default max_dune_height is 3 m. Note that dune height in Barrier3D is measured
    as the height above the dune toe (berm elevation).

    Parameters
    ----------
    b3d ndarray,
        Current instance of B3D model
    max_dune_height: float, optional
        Maximum dune height for dune rebuilding [m]
    min_dune_height: float, optional
        Minimum dune height for dune rebuilding [m]
    dune_construction_distance: int, optional
        Cell where dunes will be rebuilt
    dz: int, optional
        Vertical discretization of z [default is dz=10, dam]
    rng: bool or np.random.Generator, optional
        If `True`, add random perturbations alongshore to dune height. `rng`
        can also be an object that provides a `uniform` method, like `numpy.random`
        or `numpy.random.Generator`.

    Returns
    -------
    new_dune_domain: ndarray, shape (ny, nx)
        New yxz dune domain in new_road_domain: in units of dx, dy, dz
    rebuild_dune_volume: float
        Volume of sand for dune rebuild, in units of dx*dy*dz
    """

    dune_construction_distance = int(dune_construction_distance)

    if rng and isinstance(rng, bool):
        rng = np.random.default_rng(seed=1973)

    # Find the elevation of dune construction area
    dune_construction_area_elevation = b3d.InteriorDomain[dune_construction_distance,:]

    ny = len(dune_construction_area_elevation)
    nx = 1

    # convert from m to grid z discretization
    dune_height = np.empty((ny, 1), dtype=float)
    dune_height[:, 0] = max_dune_height / dz

    # add some random perturbations to dune heights
    if rng:
        dune_height += rng.uniform(high=0.01, size=dune_height.size).reshape((ny, 1))



    rebuild_dune_volume = np.sum(dune_height - dune_construction_area_elevation)

    b3d.InteriorDomain[dune_construction_distance, :] = dune_height[:,0]

    return dune_height, rebuild_dune_volume





def set_growth_parameters(
    yxz_dune_grid,
    Dmax,
    growthparam,
    original_growth_param=None,
    rmin=0.35,
    rmax=0.85,
):
    """Set dune growth rate to zero for next time step if the dune elevation
    (front row) is larger than the natural eq. dune height (Dmax).

    We understand from modeling work (Duran and Moore, 2013) and from empirical
    evidence (Houser et al., 2015) that dunes can reach a maximum height due to
    negative feedacks between the (cross-shore) wind field and the dune. For
    non-normal wind incidence, it has been suggested that dunes may continue to
    grow in height, albeit at a very slow rate (Davidson Arnott et al., 2018).

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
        New growth parameter array that accounts for human modifications to dune
        height above/below the equilibrium
    """

    ny = np.size(yxz_dune_grid, 0)
    new_growth_param = np.copy(growthparam)
    rng = np.random.default_rng(seed=1973)

    for idune in range(ny):
        if yxz_dune_grid[idune, 0] > Dmax:  # if dune height above dmax, don't grow
            new_growth_param[0, idune] = 0

        else:
            # if dune is now below the Dmax (was formerly above), make sure it has
            # a growth rate either the same as before (if original growth rate
            # provided) or a random number between rmin and rmax
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
    """
    For a given setback distance, check what the road elevation would be if
    relocated at grade.

    If zero MSL, the roadway can't be relocated (height drowned).

    Parameters
    ----------
    time_index: int,
        Time index for drowning error message
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; for Barrier3d,
        dz=10, decameters MHW]
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
        roadway_drown: if the roadway elevation would be zero MSL, the roadway
        can't be constructed or relocated
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

    # if the roadway elevation would be zero MSL, the roadway can't be constructed
    # or relocated
    if road_ele <= 0:
        roadway_drown = 1
        # print(
        #     f"Roadway cannot be relocated at {time_index - 1} years b/c the road "
        #     "would be at or below MSL"
        # )

    return road_ele, roadway_drown

def road_relocation_checks(
    time_index,
    dune_migrated,
    road_setback,
    road_relocation_setback,
    road_relocation_width,
    average_barrier_width,
    relocate_now=False,
):
    """Check whether the roadway should be relocated and whether relocation is viable.

    Relocation can be triggered in either of two ways:

    1. automatically, when dune migration makes ``road_setback < 0``; or
    2. explicitly, when ``relocate_now=True`` for a known historical relocation.

    ``relocate_now=True`` represents an observed hindcast relocation and therefore
    imposes the supplied setback even when the modeled barrier is too narrow. Automatic
    forecast-style relocation still retains the normal width feasibility check.

    Parameters
    ----------
    time_index: int
        Time index for relocation diagnostics.
    dune_migrated: float
        Number of meters the dune migrated in Barrier3D; positive for progradation
        and negative for erosion/landward migration [m].
    road_setback: float
        Road setback from the edge of the interior domain at the previous time step
        [m].
    road_relocation_setback: float
        Target setback for the relocated roadway [m].
    road_relocation_width: float
        Width of the relocated roadway [m].
    average_barrier_width: float
        Average barrier width at the current time step [m].
    relocate_now: bool, optional
        Request relocation during this time step, regardless of whether the modeled
        dune has reached the roadway. Intended for prescribed historical relocation
        years in a hindcast. Default is ``False``.
    Returns
    -------
    road_relocated: bool
        ``True`` when a prescribed hindcast relocation is imposed or when an automatic
        relocation request passes the normal width check.
    road_setback: float
        Updated setback after dune migration and, when successful, relocation.
    relocation_break: bool
        ``True`` when relocation was requested but the island is too narrow.
    """

    relocation_break = False
    road_relocated = False

    # Keep the road fixed in geographic space while the modeled dune line migrates.
    # This update must occur every year, including a prescribed relocation year.
    road_setback = road_setback + dune_migrated

    # Automatic relocation: dunes have migrated landward past the road.
    # Hindcast relocation: the user prescribes that relocation occurred this year.
    relocation_requested = bool(relocate_now) or (road_setback < 0)

    if relocation_requested:
        minimum_required_width = (
            road_relocation_setback + 2 * road_relocation_width
        )

        # An explicit relocate_now request is an observed hindcast boundary
        # condition, so impose the supplied road location. Automatic relocation
        # caused by dune migration retains the normal width feasibility check.
        if relocate_now:
            road_relocated = True
            road_setback = road_relocation_setback
        elif minimum_required_width > average_barrier_width:
            relocation_break = True
            # time_index - 1 because Barrier3D advances the time step at the end of
            # dune_update.
            # print(
            #     "Island is too narrow for roadway relocation at "
            #     f"{time_index - 1} years"
            # )
        else:
            road_relocated = True
            road_setback = road_relocation_setback

    return road_relocated, road_setback, relocation_break

def check_sandbag_need(
        dune_road_distance,
        design_elevation,
        barrier3d,
        sandbag_status,
        threshold_elevation = 0.101,
):
    time_index = barrier3d.time_index -1
    if dune_road_distance == 0:
        min_elev = np.min(barrier3d.DuneDomain[time_index,:,0])
        exceeds_min_dune_threshold = np.min(barrier3d.DuneDomain[time_index,:,0]) < threshold_elevation

        if exceeds_min_dune_threshold == True:
            for width in range(0, barrier3d.DuneWidth):
                for cell in range(0, len(barrier3d.DuneDomain[time_index, :, width])):
                    if barrier3d.DuneDomain[time_index, cell, width] < threshold_elevation:
                        barrier3d._DuneRestart[width][cell] = design_elevation / 10
                        c = 10

        if exceeds_min_dune_threshold == True:
            sandbag_need = True
        elif sandbag_status == True:
            sandbag_need = True
        else:
            sandbag_need = False
        c = 'end'
    elif dune_road_distance != 0:
        # If road is too far away reset to initial threshold rebuild value
        for width in range(0,barrier3d.DuneWidth):
            for cell in range(0,len(barrier3d.DuneDomain[time_index,:,width])):
                barrier3d._DuneRestart[width][cell] = 0.075
        sandbag_need = False

    return sandbag_need

class RoadwayManager:
    """Manage the road!

    Examples
    --------
    # >>> from cascade.roadway_manager import RoadwayManager
    # >>> roadways = RoadwayManager()
    # >>> roadways.request_relocation(100.0)
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
        road_relocation_setback=30,
        allow_causeway=False,
    ):
        """The RoadwayManager module

        Parameters
        ----------
        initial_road_elevation: float, optional
            Initial elevation of the roadway [m MHW]; road relocation elevations
            are built to grade based on setback.
        road_width: int, optional
            Width of roadway [m]. Also used for relocation width.
        road_setback: int, optional
            Setback of roadway from the inital dune line [m]. Also used for
            relocation setback.
        initial_dune_design_elevation: float, optional
            Elevation which dune is originally rebuilt to when road established
            [m MHW]. Used for rebuild dune height.
        initial_dune_minimum_elevation: float, optional
            Elevation threshold which originally triggers rebuilding of dune [m MHW].
            Used for min dune height.
        time_step_count: int, optional
            Number of time steps.
        original_growth_param: optional
            Dune growth parameters from first time step of Barrier3d, before
            human modifications [unitless]
        road_relocation_setback: optional
            Distance from the duneline where the relocated road will be placed [m]
        allow_causway: optional
            Whether roadways drowns when surrounded by water [default is allow_causeway=FALSE]

        """

        self._road_width = road_width
        self._road_setback = road_setback
        self._road_ele = initial_road_elevation
        # can be `None` if user doesn't want to rebuild
        self._dune_design_elevation = initial_dune_design_elevation
        self._dune_minimum_elevation = initial_dune_minimum_elevation
        self._original_growth_param = original_growth_param
        self._nt = time_step_count
        self._drown_break = 0
        self._relocation_break = 0
        self._time_index = 1
        self._absolute_minimum_dune_height = 0.3
        self._percent_water_cells_touching_road = 0.2
        self._allow_causeway = allow_causeway

        # set relocation parameters to original values
        self._road_relocation_width = (
            road_width  # can be updated outside `update` within cascade
        )
        self._road_relocation_setback = (
            road_relocation_setback  # can be updated outside `update` within cascade
        )

        # One-step prescribed relocation request used by hindcast simulations.
        # These fields are deliberately separate from _road_relocation_setback
        # because Cascade.update() may overwrite that normal relocation target
        # immediately before calling this manager.
        self._relocate_now = False
        self._prescribed_relocation_setback = None
        self._prescribed_relocation_elevation = None

        # user can specify that dune rebuilding is off with `None`: mostly for
        # debugging and sensitivity testing
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
        self._interior_dunes_built_TS = np.zeros(self._nt)  # when interior dunes are built (boolean)
        self._interior_dunes_volume_TS = np.zeros(self._nt) # volume (m^3) of sediment used to construct interior dunes
        # total overwash removed from roadway [m^3]
        self._road_overwash_volume = np.zeros(self._nt)
        # Shoreface-only nourishment performed while roadway management remains on.
        # Beach width is supplied explicitly by the simulation and follows the same
        # 0-m dune-migration threshold used by BeachDuneManager.
        self._nourishment_TS = np.zeros(self._nt)
        self._nourishment_volume_TS = np.zeros(self._nt)  # m^3/m
        self._beach_width_TS = [np.nan] * self._nt  # m
        # keep track of what percent of the dune elevations fall below minimum threshold
        self._percent_below_min = [None] * self._nt
        self._growth_params = [
            np.nan
        ] * self._nt  # when dune growth parameters set to zero b/c of rebuild height
        self._growth_params[0] = original_growth_param

        # also keep track of post-storm dune and interior impacts before human
        # modifications
        self._post_storm_dunes = [None] * self._nt
        self._post_storm_interior = [None] * self._nt
        self._post_storm_ave_interior_height = [None] * self._nt

    def update(
        self,
        barrier3d,
        trigger_dune_knockdown,
        relocate_now=None,
    ):
        """Apply roadway management for the current Barrier3D time step.

        Parameters
        ----------
        barrier3d
            Current Barrier3D model instance.
        trigger_dune_knockdown: bool
            Whether to reset dunes when roadway management terminates.
        relocate_now: bool, optional
            Impose a prescribed hindcast road relocation during the current time
            step. When omitted, a one-step request queued with
            :meth:`request_relocation` is used. The prescribed hindcast location
            bypasses event-year width and drowning rejection.
        """
        self._time_index = barrier3d.time_index

        # Consume the queued hindcast request before any possible early return.
        # Clearing it here guarantees that a historical event applies for only
        # one model update.
        queued_relocate_now = bool(self._relocate_now)
        queued_relocation_setback = self._prescribed_relocation_setback
        queued_relocation_elevation = self._prescribed_relocation_elevation
        self._relocate_now = False
        self._prescribed_relocation_setback = None
        self._prescribed_relocation_elevation = None

        if relocate_now is None:
            relocate_now = queued_relocate_now
        else:
            # Keep compatibility with direct calls while allowing a queued request
            # to survive Cascade.update(), which calls this method with two args.
            relocate_now = bool(relocate_now) or queued_relocate_now

        # A prescribed hindcast relocation reactivates roadway management before
        # relocation is evaluated. The request is treated as an observed boundary
        # condition for this one update.
        if relocate_now:
            self._drown_break = 0
            self._relocation_break = 0

        effective_relocation_setback = self._road_relocation_setback
        if relocate_now and queued_relocation_setback is not None:
            effective_relocation_setback = float(queued_relocation_setback)

        if self._original_growth_param is None:
            self._original_growth_param = barrier3d.growthparam

        # save post-storm dune and interior domain before human modifications
        # (essentially a 0.5 yr time step)
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
            self._relocation_break
        ] = road_relocation_checks(
            self._time_index,
            dune_migration,
            self._road_setback,  # current road setback, m
            effective_relocation_setback,  # normal or prescribed target setback, m
            self._road_relocation_width,  # width specified for relocation, m
            average_barrier_width,  # current width, m
            relocate_now=relocate_now,  # enforced historical relocation trigger
        )

        # if road can't be relocated, no longer manage and exit; dune growth
        # parameters reset to original in CASCADE
        if self._relocation_break == 1:
            # an adaptation solution may be to knock down the dunes so that they
            # are small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return

        # if road is relocated, get the new road elevation (built at grade) and
        # update dune elevations which are dependent on the road elevation;
        # otherwise, decrease all elevations (m MHW) this year by the SLR increment
        if road_relocated:
            self._road_width = self._road_relocation_width
            previous_road_elevation = float(self._road_ele)
            grade_road_elevation, grade_drown = get_road_relocation_elevation(
                self._time_index,
                # interior domain from this last time step, dam
                xyz_interior_grid=barrier3d.InteriorDomain,
                road_setback=self._road_setback,  # m
                road_width=self._road_width,  # m
                dx=10,
                dy=10,
                dz=10,  # specifies interior is in dam
            )

            if relocate_now:
                # A prescribed hindcast relocation is imposed even when the modeled
                # grade is at/below MHW. An explicitly prescribed elevation takes
                # priority; otherwise retain the positive pre-relocation roadway
                # elevation as construction fill.
                if queued_relocation_elevation is not None:
                    self._road_ele = float(queued_relocation_elevation)
                elif np.isfinite(grade_road_elevation) and grade_road_elevation > 0:
                    self._road_ele = float(grade_road_elevation)
                else:
                    self._road_ele = max(previous_road_elevation, 0.01)
                self._drown_break = 0
            else:
                self._road_ele = grade_road_elevation
                self._drown_break = grade_drown

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
            # print(
            #     f"Roadway drowned in place at {self._time_index - 1} years due to "
            #     "SLR - road cannot be below 0 m MHW"
            # )

            # an adaptation solution may be to knock down the dunes so that they
            # are small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return
        elif self._drown_break == 1:  # if road drowned from road relocation above
            # an adaptation solution may be to knock down the dunes so that they are
            # small and can easily be overwashed
            if trigger_dune_knockdown:
                barrier3d.DuneDomain[self._time_index - 1, :, :] = barrier3d.DuneDomain[
                    0, :, :
                ]

            return

        # when the roadway gets really low in elevation, the dune_design_elevation
        # may not be above the berm; when this happens, we use a design height of
        # 1 m above the berm to keep a dune to protect the roadway and rebuild
        # whenever the dune drops to just above elevation of the berm (0.3 m) --
        # essentially, we just push the sand back
        self._dune_design_elevation = max(
            self._dune_design_elevation, (barrier3d.BermEl * 10) + 1.0
        )
        # note, Barrier3D adds a random seeded height for the proto/new dune line
        # (20 cm); here we allow the user to specify
        self._dune_minimum_elevation = max(
            self._dune_minimum_elevation,
            (barrier3d.BermEl * 10) + self._absolute_minimum_dune_height,
        )

        # save time series
        self._road_setback_TS[self._time_index - 1] = self._road_setback
        self._road_width_TS[self._time_index - 1] = self._road_width
        self._road_ele_TS[self._time_index - 1] = self._road_ele
        self._dune_design_elevation_TS[self._time_index - 1] = (
            self._dune_design_elevation
        )
        self._dune_minimum_elevation_TS[self._time_index - 1] = (
            self._dune_minimum_elevation
        )
        self._road_relocated_TS[self._time_index - 1] = road_relocated

        ###############################################################################
        # bulldoze roadway after storms and check for road width drowning
        ###############################################################################

        # bulldoze the road and put bulldozed sand back on the dunes; drown road
        # when a water cell touches either side
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
            # interior domain from this last time step, dam
            xyz_interior_grid=barrier3d.InteriorDomain,
            # dune domain from this last time step, dam
            yxz_dune_grid=barrier3d.DuneDomain[self._time_index - 1, :, :],
            dx=10,
            dy=10,
            dz=10,  # specifies dam for dune and interior domains
            drown_threshold=0,  # 0 m MSL
            # fraction cells<drown_threshold
            percent_water_cells_touching_road=self._percent_water_cells_touching_road,
            # During the prescribed hindcast relocation year, permit construction
            # fill/causeway so adjacent modeled water does not cancel the observed
            # event.
            allow_causeway=(self._allow_causeway or bool(relocate_now))
        )
        if self._drown_break == 1:
            # an adaptation solution may be to knock down the dunes so that they
            # are small and can easily be overwashed
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
            # all in dam MHW
            new_xyz_interior_domain[new_xyz_interior_domain >= barrier3d.SL]
        )
        # slightly altered due to roadway
        barrier3d.h_b_TS[-1] = new_ave_interior_height
        barrier3d.InteriorDomain = new_xyz_interior_domain
        barrier3d.DomainTS[self._time_index - 1] = new_xyz_interior_domain

        ###############################################################################
        # Rebuild dunes if knocked down
        ###############################################################################

        # dune management: rebuild dunes!
        # Dune min elevation should likely be altered
        if self._dune_design_elevation is None or self._dune_minimum_elevation is None:
            pass
        else:
            # in B3D, dune height is the height above the berm crest
            dune_design_height = self._dune_design_elevation - (barrier3d.BermEl * 10)
            min_dune_height = self._dune_minimum_elevation - (barrier3d.BermEl * 10)

            # if any dune cell in the front row of dunes is less than a minimum
            # threshold -- as measured above the berm crest -- then rebuild the
            # dune (all rows up to dune_design_elevation)
            if np.min(new_dune_domain[:, 0]) < (min_dune_height / 10) and self._road_setback <= 10:  # in m
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
            elif np.min(new_dune_domain[:, 0]) < (min_dune_height / 10) and self._road_setback <= 20:
                Interior_Dune_Front = (self._road_setback/10)-2
                Interior_Dune_Back = (self._road_setback/10)-1
                new_dune_domain, rebuild_dune_volume = build_interior_dunes(
                    b3d=barrier3d,
                    dune_construction_distance=Interior_Dune_Front,
                    max_dune_height=dune_design_height,
                    min_dune_height=dune_design_height,
                    dz=10,
                    rng=True
                )

                self._interior_dunes_built_TS[self._time_index - 1] = 1
                self._interior_dunes_volume_TS[self._time_index - 1] = (
                    rebuild_dune_volume * dm3_to_m3
                )

        # update Barrier3D class variables
        barrier3d.DuneDomain[self._time_index - 1, :, :] = new_dune_domain

        ###############################################################################
        # Set dune growth rate to zero if dunes rebuilt > natural equilibrium height
        ###############################################################################

        # set dune growth rate to zero for next time step if the dune elevation
        # (front row) is larger than the natural eq. dune height (Dmax)
        new_growth_parameters = set_growth_parameters(
            new_dune_domain,  # in dam
            barrier3d.Dmax,  # in dam
            barrier3d.growthparam,
            # use original growth rates for resetting values
            original_growth_param=self._original_growth_param,
        )
        self._growth_params[self._time_index - 1] = copy.deepcopy(new_growth_parameters)

        # update class variables
        barrier3d.growthparam = new_growth_parameters

        return

    def update_beach_width(
        self,
        barrier3d,
        beach_width_last_year,
    ):
        """Update roadway-domain beach width and dune migration for one model year.

        This duplicates the existing BeachDuneManager rule without applying any
        community overwash filtering or community dune rebuilding. It must be called
        once per simulated year for a roadway domain that is tracking beach width.

        Parameters
        ----------
        barrier3d
            Current Barrier3D instance after its annual physical update.
        beach_width_last_year: float
            Managed beach width from the preceding model year [m].

        Returns
        -------
        float
            Current beach width after annual shoreline change [m].
        """

        beach_width_last_year = float(beach_width_last_year)
        if not np.isfinite(beach_width_last_year) or beach_width_last_year < 0:
            raise ValueError(
                "Previous roadway beach width must be finite and non-negative; "
                f"received {beach_width_last_year}."
            )

        time_index = barrier3d.time_index
        change_in_shoreline = (
            barrier3d.x_s_TS[-1] - barrier3d.x_s_TS[-2]
        ) * 10  # m
        current_beach_width = beach_width_last_year - change_in_shoreline
        current_beach_width = beach_width_dune_dynamics(
            current_beach_width=current_beach_width,
            beach_width_last_year=beach_width_last_year,
            beach_width_threshold=0,
            barrier3d=barrier3d,
            time_index=time_index,
        )

        output_index = time_index - 1
        if 0 <= output_index < self._nt:
            self._beach_width_TS[output_index] = current_beach_width

        return current_beach_width

    def nourish_beach(
        self,
        barrier3d,
        nourishment_volume,
        beach_width,
    ):
        """Nourish a roadway-managed domain without community dune management.

        This method applies only the shoreface/beach nourishment calculation used by
        BeachDuneManager. It intentionally does not:

        * rebuild or otherwise modify ``DuneDomain``;
        * filter or bulldoze overwash;
        * modify roadway dune design/minimum elevations;
        * modify dune growth parameters; or
        * use community dune-rebuilding or overwash-filtering rules.

        Dune migration follows the existing BeachDuneManager beach-width rule:
        migration is off while beach width is above 0 m and turns on when beach
        width reaches 0 m. Nourishment widens the beach, so migration is turned off
        again after a positive nourishment event.

        The simulation script must explicitly choose the domain, year, nourishment
        volume, and current beach width. The updated beach width is returned so the
        script can retain it for the next nourishment event.

        Parameters
        ----------
        barrier3d
            Current Barrier3D model instance after its annual physical update.
        nourishment_volume: float
            Nourishment volume [m^3/m].
        beach_width: float
            Beach width immediately before nourishment [m].

        Returns
        -------
        float
            Updated beach width [m].
        """

        nourishment_volume = float(nourishment_volume)
        beach_width = float(beach_width)

        if not np.isfinite(nourishment_volume) or nourishment_volume < 0:
            raise ValueError(
                "Roadway nourishment volume must be finite and non-negative; "
                f"received {nourishment_volume}."
            )
        if not np.isfinite(beach_width) or beach_width < 0:
            raise ValueError(
                "Roadway beach width must be finite and non-negative; "
                f"received {beach_width}."
            )

        (
            barrier3d.x_s,
            barrier3d.s_sf_TS[-1],
            new_beach_width,
        ) = shoreface_nourishment(
            x_s=barrier3d.x_s,  # dam
            x_t=barrier3d.x_t,  # dam
            nourishment_volume=nourishment_volume / 100,  # m^3/m to dam^3/dam
            average_barrier_height=barrier3d.h_b_TS[-1],  # dam
            shoreface_depth=barrier3d.DShoreface,  # dam
            beach_width=beach_width / 10,  # m to dam
        )

        new_beach_width *= 10  # dam to m
        barrier3d.x_s_TS[-1] = barrier3d.x_s

        # Match BeachDuneManager: nourishment restores a positive beach buffer, so
        # the dune line is pinned until annual erosion reduces beach width to 0 m.
        barrier3d.dune_migration_on = new_beach_width <= 0

        # Keep the full barrier-position bookkeeping internally consistent when this
        # method is called directly from a simulation script.
        barrier3d.x_b_TS[-1] = (
            barrier3d.x_s
            + barrier3d.InteriorWidth_AvgTS[-1]
            + np.size(barrier3d.DuneDomain, 2)
            + (new_beach_width / 10)
        )

        time_index = barrier3d.time_index - 1
        if 0 <= time_index < self._nt:
            self._nourishment_TS[time_index] = 1
            self._nourishment_volume_TS[time_index] = nourishment_volume
            self._beach_width_TS[time_index] = new_beach_width

        return new_beach_width

    def request_relocation(
        self,
        road_setback,
        road_elevation=None,
    ):
        """Queue one prescribed road relocation for the next update.

        Parameters
        ----------
        road_setback: float
            Post-relocation roadway setback from the current modeled dune/interior
            boundary [m].
        road_elevation: float, optional
            Prescribed post-relocation roadway elevation [m MHW]. When omitted, the
            manager uses modeled grade when positive; otherwise it retains the
            positive pre-relocation road elevation as construction fill.

        Notes
        -----
        The request is consumed and cleared the next time ``update`` is called.
        Every queued request is treated as an observed hindcast boundary condition:
        it resets prior roadway break flags and bypasses event-year width and
        drowning rejection without modifying ``cascade.py``.
        """
        road_setback = float(road_setback)

        if not np.isfinite(road_setback) or road_setback < 0:
            raise ValueError(
                "Prescribed road relocation setback must be a finite, "
                f"non-negative value; received {road_setback}."
            )

        if road_elevation is not None:
            road_elevation = float(road_elevation)
            if not np.isfinite(road_elevation) or road_elevation <= 0:
                raise ValueError(
                    "Prescribed relocation road elevation must be finite and "
                    f"greater than zero; received {road_elevation}."
                )

        self._prescribed_relocation_setback = road_setback
        self._prescribed_relocation_elevation = road_elevation
        self._relocate_now = True

        # Reactivate roadway management for the observed historical event so the
        # next normal Cascade.update() reaches and enforces this manager request.
        self._drown_break = 0
        self._relocation_break = 0

    @property
    def road_setback(self):
        """Current modeled roadway setback [m]."""
        return self._road_setback

    @property
    def relocate_now(self):
        """Whether a prescribed relocation is queued for the next update."""
        return self._relocate_now

    @property
    def prescribed_relocation_setback(self):
        """Queued post-relocation setback [m], or ``None``."""
        return self._prescribed_relocation_setback

    @property
    def nourishment_TS(self):
        """Years in which roadway-domain beach nourishment was applied."""
        return self._nourishment_TS

    @property
    def nourishment_volume_TS(self):
        """Roadway-domain nourishment volume time series [m^3/m]."""
        return self._nourishment_volume_TS

    @property
    def beach_width_TS(self):
        """Post-nourishment roadway-domain beach widths [m]."""
        return self._beach_width_TS

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
