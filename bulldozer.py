"""Bulldozer

This module provides functions modifying a 3D grid (X,Y,Z) for human coastal management decision,
including 1) bulldozing sand from roadways and adding to the dune line ....

References
----------

.. [1] NCDOT (North Carolina Department of Transportation), 2008. NC 12 Transportation Improvements. http://www.ncdot.org/projects/
NC12/ (accessed October 24, 2007).
.. [2]
.. [3] P.D. Komar, 1998, Beach processes and sedimentation: Upper Saddle River, New Jersey, Prentice Hall , 544 p.
.. [4] Jaap H. Nienhuis, Jorge Lorenzo Trueba; Simulating barrier island response to sea level rise with the barrier
    island and inlet environment (BRIE) model v1.0 ; Geosci. Model Dev., 12, 4013–4030, 2019;
    https://doi.org/10.5194/gmd-12-4013-2019


Notes
---------

"""
import numpy as np
from scipy.interpolate import interp2d


def bulldoze(
    xyz_interior_grid,
    yxz_dune_grid,
    road_ele=2.0,
    road_width=20,
    road_setback=30,
    dx=10,
    dy=10,
    dz=10,
):
    r"""Remove overwash from roadway and put it back on the dune. Spreads sand evenly across dune cells.

    Parameters
    ----------
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; if dz=10]
    yxz_dune_grid: array
        Dune topography [z units specified by dz; if dz=10, decameters]
    road_ele: float
        Road elevation [m]
    road_width: int
        Width of roadway [m]
    road_setback: int
        Setback distance of roadway from edge of interior domain [m]
    dx: int
        Cross-shore discretization of x [m]
    dy: int
        Alongshore discretization of y [m]
    dz: int
        Vertical discretization of z [m]

    Returns
    -------
    array of float
        new_road_domain: in units of dx, dy, dz
        new_dune_domain: in units of dx, dy, dz

    """

    # road parameters
    road_start = int(
        road_setback / dy
    )  # grid index for start of roadway in interior domain (for B3D, convert to dm)
    road_width = int(road_width / dx)
    road_end = road_start + road_width
    road_ele = (
        road_ele / dz
    )  # convert to units of grid (NOTE: in B3D default simulation, berm is 1.44 m)

    # remove sand from roadway (only account for positive values)
    old_road_domain = xyz_interior_grid[road_start:road_end, :]
    new_road_domain = np.zeros((road_width, np.size(old_road_domain, 1))) + road_ele
    road_overwash_removal = sum(old_road_domain - new_road_domain)
    road_overwash_removal[road_overwash_removal < 0] = 0

    # spread overwash removed from roadway equally over all dune cells
    number_dune_cells = np.size(yxz_dune_grid, 1)
    overwash_to_dune = np.transpose(
        [road_overwash_removal / number_dune_cells] * number_dune_cells
    )
    new_dune_domain = yxz_dune_grid + overwash_to_dune

    xyz_interior_grid[
        road_start:road_end, :
    ] = new_road_domain  # update interior domain

    return new_dune_domain, xyz_interior_grid, road_overwash_removal


def rebuild_dunes(
    yxz_dune_grid, max_dune_height=3.0, min_dune_height=2.4, dz=10, rng=True
):
    r"""Construct artificial dunes after each overwash or inundation event. Along NC-12, artificial dune heights range
    from 2.4 to 4.6 m (NCDOT, 2008; Overton, Judge, and Fisher, 2000; Magliocca et al., 2011). From a more recent pub
    (Valasquez et al., 2020), the average elevation of the road along NC-12 is 1.3 m (NAVD88); they find that in order
    for the road to not be vulnerable to overwash, the dune crest must be higher than 4.3 m (NAVD88), so here the
    default max_dune_height is 3 m (here, dune height is measured as the height above the dune toe).

    Note that while artificial dune geometry of a given height is constrained by the angle of repose of beach, here we
    just assume the dunes are built to a width capable of surviving impacts from a "moderate" storm (25-m wide along
    NC-12, corresponding to a 50-year event; NCDOT, 2008, Overton, Judge, and Fisher, 2000).

    For the given dune grid, we apply a linear gradient from the first to last dune row given the max and min dune
    height, with small random perturbations.

    Parameters
    ----------
    yxz_dune_grid: float or array of float
        Dune topography [z units specified by dz; if dz=10, decameters]
    max_dune_height: float
        Maximum dune height [m]
    min_dune_height: float
        Minimum dune height [m]
    dz: int
        Vertical discretization of z [m]
    rng: boolean
        If True, add random perturbations alongshore to dune height

    Returns
    -------
    new_dune_domain: float or array of float
        New yxz dune domain in z units specified by dz

    """

    # convert from m to grid z discretization
    max_height = max_dune_height / dz
    min_height = min_dune_height / dz
    ny = np.size(yxz_dune_grid, 0)
    nx = np.size(yxz_dune_grid, 1)

    if rng:
        # add some random perturbations to dune heights
        dune_start_max = np.ones([ny]) * (
            max_height + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
        )
        dune_start_min = np.ones([ny]) * (
            min_height + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
        )
    else:
        dune_start_max = np.ones([ny]) * max_height
        dune_start_min = np.ones([ny]) * min_height

    # linearly interpolate (2D) from front row (max height) to back row (min height)
    x = [0, nx - 1]
    y = [np.arange(0, ny, 1)]
    z = np.transpose([dune_start_max, dune_start_min])
    f = interp2d(x, y, z)
    new_dune_domain = f(np.arange(0, nx, 1), np.arange(0, ny, 1))

    return new_dune_domain
