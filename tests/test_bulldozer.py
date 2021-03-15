import pytest
import numpy as np
from scipy.sparse import csr_matrix

from bulldozer import bulldoze


@pytest.mark.parametrize("func", bulldoze)
def test_bulldoze_volume(func):
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
    # new_dune_domain == np.zeros([20, 10]) + 2  # should equal 2 vertical m of accretion
