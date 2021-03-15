"""Bulldozer

This module provides functions modifying a 3D grid (X,Y,Z) for human coastal management decision,
including 1) bulldozing sand from roadways and adding to the dune line ....

References
----------


Notes
---------

"""
import numpy as np


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
    r"""Remove overwash from roadway and put it back on the dune

    Parameters
    ----------
    xyz_interior_grid: array
        Interior barrier island topography [z units specified by dz; if dz=10, decameters]
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
