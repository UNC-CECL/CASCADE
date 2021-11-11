import numpy as np
from numpy.testing import assert_array_almost_equal
from pathlib import Path
from tqdm import tqdm

from cascade.roadway_manager import bulldoze, rebuild_dunes, set_growth_parameters
from cascade.beach_dune_manager import shoreface_nourishment, filter_overwash
from cascade import Cascade

datadir = "../Cascade/B3D_Inputs/"
Path(datadir)


def test_bulldoze_volume():

    # bulldoze 1 vertical m of sand from a 20 m wide roadway
    xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    (
        new_dune_domain,
        new_xyz_interior_grid,
        road_overwash_removal,
        roadway_drown,
    ) = bulldoze(
        time_index=1,
        road_ele=2.0,  # in meters
        road_width=20,  # in meters
        road_setback=30,  # in meters
        xyz_interior_grid=xyz_interior_grid,
        yxz_dune_grid=yxz_dune_grid,
        dx=1,
        dy=1,
        dz=1,
        drown_threshold=0,
    )
    assert road_overwash_removal == np.sum(
        np.zeros([20]) + (1 * 20)
    )  # should equal 1 vertical m of overwash removed per road cell (20), summed across the domain [m^3]
    assert_array_almost_equal(
        new_dune_domain, np.zeros([20, 10]) + 2
    )  # should equal 2 vertical m of accretion

    # now test that decameters works properly; bulldoze 0.1 vertical decameter of sand from a 2 decameter wide roadway
    xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0 / 10
    )  # units of decameters in X (cross-shore) Y (alongshore) Z (vertical)
    yxz_dune_grid = np.zeros([20, 2])  # dune domain is 2 cells (m) wide
    (
        new_dune_domain,
        new_xyz_interior_grid,
        road_overwash_removal,
        roadway_drown,
    ) = bulldoze(
        time_index=1,
        road_ele=2.0,  # in meters
        road_width=20,  # in meters
        road_setback=30,  # in meters
        xyz_interior_grid=xyz_interior_grid,
        yxz_dune_grid=yxz_dune_grid,
        dx=10,
        dy=10,
        dz=10,
        drown_threshold=0,
    )
    assert_array_almost_equal(
        road_overwash_removal, np.floor(np.sum(np.zeros([20]) + (0.1 + 0.1)))
    )  # should equal 0.1 vertical decameters of overwash removed per road cell (2), summed across the domain [dm^3]
    assert_array_almost_equal(
        new_dune_domain, (np.zeros([20, 2]) + 1 / 10)
    )  # should equal 0.1 vertical decameter of accretion along each dune cell

    # NOTE: need to make a test that shows unequal overwash volume removed per grid cell and placement on adjacent dunes


def test_rebuild_dunes_interpolation():
    yxz_dune_grid = np.zeros([20, 11])  # dune domain is 11 cells (m) wide
    new_dune_domain, rebuild_dune_volume = rebuild_dunes(
        yxz_dune_grid, max_dune_height=2.4, min_dune_height=1.4, dz=1, rng=False
    )
    assert_array_almost_equal(
        new_dune_domain, np.array([np.arange(2.4, 1.3, -0.1)] * 20)
    )


def test_growth_params():

    # growth parameter should only change when yxz_dune_grid is greater than Dmax
    yxz_dune_grid = np.zeros([20, 11]) + 2
    Dmax = 3
    growthparam = np.zeros([1, 20]) + 0.7
    rmin = 0.6
    rmax = 0.8
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
    assert all([a == b] for a, b in zip(new_growth_parameters, [0.7, 0.7, 0.7, 0.0]))

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
    assert all([a == b] for a, b in zip(new_growth_parameters, [0.7, 0.7, 0.6, 0.6]))


def test_shoreface_nourishment():

    nourishment_volume = [0, 100, 300]
    new_beach_width = np.zeros(3)
    x_t = 0  # m
    x_s = 500  # m
    average_barrier_height = 2  # m
    shoreface_depth = 10  # m
    beach_width = 0  # m

    for i in range(3):

        [new_x_s, s_sf, new_beach_width[i]] = shoreface_nourishment(
            x_s,
            x_t,
            nourishment_volume[i],
            average_barrier_height,
            shoreface_depth,
            beach_width,
        )

    assert all([a == b] for a, b in zip(np.floor(new_beach_width), [0, 14, 42]))


def test_overwash_filter():
    """
    check that overwash removal volumes are correct and that only overwash is removed (i.e., eroded areas do not
    become deeper)
    """

    barrier_overwash_removed = np.zeros(3)

    post_xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    pre_xyz_interior_grid = post_xyz_interior_grid - 1.0
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    # remove 40% of 1 m^3 of overwash in each cell,
    (
        new_yxz_dune_domain,  # [m], should be 800 m^3 / (10 * 20) added to each cell (4 m -- woah)
        new_xyz_interior_domain,  # [m], should be 2.6 m depth in each cell
        barrier_overwash_removed[0],  # [m^3], should be 0.4 m^3 * 100 * 20 (800 m^3)
    ) = filter_overwash(
        overwash_filter=40,  # remove 40% of overwash deposit
        post_storm_xyz_interior_grid=post_xyz_interior_grid,
        pre_storm_xyz_interior_grid=pre_xyz_interior_grid,
        post_storm_yxz_dune_grid=yxz_dune_grid,
        artificial_maximum_dune_height=0.5,  # maximum height of artificial dunes in dam
    )

    # check that it only removes overwash, doesn't deepen areas where there wasn't overwash or eroded
    post_xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    pre_xyz_interior_grid = post_xyz_interior_grid - 1.0
    post_xyz_interior_grid[
        0, :
    ] = 2.0  # now make it so overwash wasn't deposited in some areas
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    # remove 40% of 1 m^3 of overwash in each cell,
    (
        new_yxz_dune_domain,  # [m], should be 792 m^3 / (10 * 20) added to each cell (3.96 m -- woah)
        new_xyz_interior_domain,  # [m], should be 2.6 m depth in each cell except for row 0, which should be 2.0 m
        barrier_overwash_removed[1],  # [m^3], should be 0.4 m^3 * 99 * 20 (792 m^3)
    ) = filter_overwash(
        overwash_filter=40,  # remove 40% of overwash deposit
        post_storm_xyz_interior_grid=post_xyz_interior_grid,
        pre_storm_xyz_interior_grid=pre_xyz_interior_grid,
        post_storm_yxz_dune_grid=yxz_dune_grid,
        artificial_maximum_dune_height=0.5,  # maximum height of artificial dunes in dam
    )

    post_xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    pre_xyz_interior_grid = post_xyz_interior_grid - 1.0
    post_xyz_interior_grid[
        0, :
    ] = 1.0  # now make it so some places are eroded --> result should be same as above
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    # remove 40% of 1 m^3 of overwash in each cell,
    (
        new_yxz_dune_domain,  # [m], should be 792 m^3 / (10 * 20) added to each cell (3.96 m -- woah)
        new_xyz_interior_domain,  # [m], should be 2.6 m depth in each cell except for row 0, which should be 1.0 m
        barrier_overwash_removed[2],  # [m^3], should be 0.4 m^3 * 99 * 20 (792 m^3)
    ) = filter_overwash(
        overwash_filter=40,  # remove 40% of overwash deposit
        post_storm_xyz_interior_grid=post_xyz_interior_grid,
        pre_storm_xyz_interior_grid=pre_xyz_interior_grid,
        post_storm_yxz_dune_grid=yxz_dune_grid,
        artificial_maximum_dune_height=0.5,  # maximum height of artificial dunes in dam
    )

    assert all(
        [a == b] for a, b in zip(np.ceil(barrier_overwash_removed), [800, 792, 792])
    )
    # assert_array_almost_equal(barrier_overwash_removed, [800, 792, 792])


def test_shoreline_migration():
    """
    As a check on the dynamics in Barrier3D, here we want to see if the dunes migrate when the beach width goes to zero
    and the shoreline surpasses a full cell width (10 m). Indeed, dunes only migrate when humans allow them to
    (beach width = 0), and when the shoreline moves a full cell width -- in year 57 and 65. Note that if one was to
    check the `post_storm_x_s`, they would find that the dunes actually migrated at 56.5 since dune migration
    occurrs prior to any human modifications.
    """

    iB3D = 0
    total_time = 100

    cascade = Cascade(
        datadir,
        name="test",
        sea_level_rise_rate=0.007,
        alongshore_section_count=1,
        time_step_count=total_time,
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=True,
        community_dynamics_module=False,
        dune_design_elevation=3.7,
        nourishment_interval=10,  # yrs
        nourishment_volume=100,
        background_erosion=-1.0,  # m/yr
        overwash_filter=40,
    )

    # Loop for 50 years at a 10 year interval, 100 m^3/m and then 50 years at a 20 year interval with 300 m^3/m
    nt = 50
    for _ in tqdm(range(nt - 1)):
        cascade.update()
        if cascade.road_break or cascade.b3d_break:
            break

    # during the CASCADE initialization, the nourishment interval and volume is specified individually for each
    # barrier3d alongshore cell; so to update these values, we need to specify which barrier3d cell we want to modify
    # (here, we only have one cell)
    cascade.nourishment_interval[iB3D] = 20  # increase to 20 years
    cascade.nourishment_volume[iB3D] = 300  # increase to 300 m^3/m

    for _ in tqdm(range(nt)):
        cascade.update()
        if cascade.road_break or cascade.b3d_break:
            break

    frac_grid_cell = np.array(cascade.barrier3d[iB3D]._x_s_TS) % 1
    diff = np.hstack([0, np.diff(frac_grid_cell)])
    diff[~np.array(cascade.nourishments[iB3D]._dune_migration_on)] = 0
    shoreline_transgressed = diff < 0
    dunes_migrated = cascade.barrier3d[iB3D]._ShorelineChangeTS < 0

    assert all(shoreline_transgressed == dunes_migrated)
