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

    return new_dune_domain, xyz_interior_grid, np.sum(road_overwash_removal)


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
        RNG = np.random.default_rng(seed=1973)

        dune_start_max = np.ones([ny]) * (
            # max_height + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
            max_height
            + (-0.01 + (0.01 - (-0.01)) * RNG.random(ny))
        )
        dune_start_min = np.ones([ny]) * (
            # min_height + (-0.01 + (0.01 - (-0.01)) * np.random.rand(ny))
            min_height
            + (-0.01 + (0.01 - (-0.01)) * RNG.random(ny))
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


def set_growth_parameters(
    yxz_dune_grid,
    Dmax,
    growthparam,
    original_growth_param=None,
    rmin=0.35,
    rmax=0.85,
):
    r"""Set dune growth rate to zero for next time step if the dune elevation (front row) is larger than the natural eq.
    dune height (DMax) -- we know from modeling work (Duran and Moore 2013) that the dune will not grow above
    the natural equilibrium height because of interactions between the wind field and the dune: too high and no sand
    flux

    Parameters
    ----------
    yxz_dune_grid: float or array of float
        Dune topography [decameters]
    Dmax: float
        Maximum natural equilibrium dune height [decameters]
    growthparam: float
        growth parameters from last time step, used in Houser formulation [unitless]
    original_growth_param: float, optional
        growth parameters from first time step, before humans interfered [unitless]
    rmin: float, optional
        Minimum growth rate - used if now original_growth_parm provided [unitless]
    rmax: float, optional
        Maximum growth rate - used if now original_growth_parm provided [unitless]

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
