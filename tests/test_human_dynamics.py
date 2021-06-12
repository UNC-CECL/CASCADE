import pytest
import numpy as np
from numpy.testing import assert_array_almost_equal

from roadway_manager import bulldoze, rebuild_dunes, set_growth_parameters
from beach_nourisher import shoreface_nourishment


def test_bulldoze_volume():

    xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    # bulldoze 1 vertical m of sand from a 20 m wide roadway
    new_dune_domain, new_xyz_interior_grid, road_overwash_removal = bulldoze(
        xyz_interior_grid,
        yxz_dune_grid,
        road_ele=2.0,  # in meters
        road_width=20,  # in meters
        road_setback=30,  # in meters
        dy=1,
        dx=1,
        dz=1,
    )
    assert road_overwash_removal == np.sum(
        np.zeros([20]) + (1 * 20)
    )  # should equal 1 vertical m of overwash removed per road cell (20), summed across the domain [m^3]
    assert (
        new_dune_domain == np.zeros([20, 10]) + 2
    )  # should equal 2 vertical m of accretion

    # now decameters
    xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0 / 10
    )  # units of decameters in X (cross-shore) Y (alongshore) Z (vertical)
    yxz_dune_grid = np.zeros([20, 2])  # dune domain is 2 cells (m) wide
    # bulldoze 0.1 vertical decameter of sand from a 2 decameter wide roadway
    new_dune_domain, new_xyz_interior_grid, road_overwash_removal = bulldoze(
        xyz_interior_grid,
        yxz_dune_grid,
        road_ele=2.0,  # in meters
        road_width=20,  # in meters
        road_setback=30,  # in meters
        dy=10,  # specify that grid in decameters
        dx=10,  # specify that grid in decameters
        dz=10,  # specify that grid in decameters
    )
    assert_array_almost_equal(
        road_overwash_removal, np.sum(np.zeros([20]) + (0.1 + 0.1))
    )  # should equal 0.1 vertical decameters of overwash removed per road cell (2), summed across the domain [dm^3]
    assert (
        new_dune_domain == np.zeros([20, 2]) + 1 / 10
    )  # should equal 0.1 vertical decameter of accretion along each dune cell


def test_rebuild_dunes_interpolation():
    yxz_dune_grid = np.zeros([20, 11])  # dune domain is 11 cells (m) wide
    new_dune_domain = rebuild_dunes(
        yxz_dune_grid, max_dune_height=2.4, min_dune_height=1.4, dz=1, rng=False
    )
    assert_array_almost_equal(
        new_dune_domain, np.array([np.arange(2.4, 1.3, -0.1)] * 20)
    )


# @pytest.mark.parametrize(
#     "func", (set_growth_parameters)
# )
def test_growth_params():
    yxz_dune_grid = np.zeros([20, 11]) + 2
    Dmax = 3  # 2  -- NOTE, I don't know how to make this work for multiple values
    growthparam = np.zeros([1, 20]) + 0.7
    rmin = 0.6
    rmax = 0.8

    # growth parameter should only change when yxz_dune_grid is greater than Dmax
    new_growth_parameters = set_growth_parameters(
        yxz_dune_grid,
        Dmax,
        growthparam,
        original_growth_param=None,
        rmin=rmin,
        rmax=rmax,
    )
    assert_array_almost_equal(growthparam, new_growth_parameters)

    # growth parameters should all go to zero
    Dmax = 1
    new_growth_parameters = set_growth_parameters(
        yxz_dune_grid,
        Dmax,
        growthparam,
        original_growth_param=None,
        rmin=rmin,
        rmax=rmax,
    )
    assert np.all(new_growth_parameters) == 0

    # only last grid cell has dune growth parameter go to zero (dune height above dune max)
    yxz_dune_grid = np.zeros([4, 2]) + 2
    yxz_dune_grid[3, :] = 4
    Dmax = 3
    growthparam = np.zeros([1, 4]) + 0.7
    new_growth_parameters = set_growth_parameters(
        yxz_dune_grid,
        Dmax,
        growthparam,
        original_growth_param=None,
        rmin=rmin,
        rmax=rmax,
    )
    assert new_growth_parameters == [0.7, 0.7, 0.7, 0.0]

    # check the scenario where at the previous time step growth rate was zero but now we want to change it back to 0.7
    # because the dune height dropped below Dmax
    yxz_dune_grid = np.zeros([4, 2]) + 2
    Dmax = 3
    original_growth_param = np.zeros([1, 4]) + 0.6
    growthparam = np.zeros([1, 4]) + 0.7
    growthparam[0, 1:3] = 0
    new_growth_parameters = set_growth_parameters(
        yxz_dune_grid,
        Dmax,
        growthparam,
        original_growth_param=original_growth_param,
        rmin=None,
        rmax=None,
    )
    assert new_growth_parameters == [0.7, 0.7, 0.6, 0.6]


# NOTE FROM KA: I don't know how to use pytest, or test it, so leave for now
# @pytest.mark.parametrize("nourishment_volume", (0, 100, 300))  # m^3/m
def test_shoreface_nourishment(nourishment_volume):
    x_t = 0  # m
    x_s = 500  # m
    average_barrier_height = 2  # m
    shoreface_depth = 10  # m
    beach_width = 0  # m

    [new_x_s, s_sf, new_beach_width] = shoreface_nourishment(
        x_s,
        x_t,
        nourishment_volume,
        average_barrier_height,
        shoreface_depth,
        beach_width,
    )

    assert np.floor(new_beach_width) == [0, 14, 42]
