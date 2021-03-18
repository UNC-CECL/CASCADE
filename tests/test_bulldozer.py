import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal

from bulldozer import bulldoze, rebuild_dunes


def test_bulldoze_volume():
    xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    new_dune_domain, new_xyz_interior_grid, road_overwash_removal = bulldoze(
        xyz_interior_grid,
        yxz_dune_grid,
        road_ele=2.0,
        road_width=20,
        road_setback=30,
        dy=1,
        dx=1,
        dz=1,
    )

    assert (
        road_overwash_removal == np.zeros([20]) + 20
    )  # should equal 1 vertical m of overwash removed per road cell
    assert (
        new_dune_domain == np.zeros([20, 10]) + 2
    )  # should equal 2 vertical m of accretion


def test_rebuild_dunes_interpolation():
    yxz_dune_grid = np.zeros([20, 11])  # dune domain is 11 cells (m) wide
    new_dune_domain = rebuild_dunes(
        yxz_dune_grid, max_dune_height=2.4, min_dune_height=1.4, dz=1, rng=False
    )

    assert_array_almost_equal(
        new_dune_domain, np.array([np.arange(2.4, 1.3, -0.1)] * 20)
    )
