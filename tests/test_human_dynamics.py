from pathlib import Path

import numpy as np
from numpy.testing import assert_array_almost_equal

from cascade import Cascade
from cascade.beach_dune_manager import filter_overwash, shoreface_nourishment
from cascade.roadway_manager import bulldoze, rebuild_dunes, set_growth_parameters

BMI_DATA_DIR = Path(__file__).parent / "cascade_test_human_inputs"
NT = 180
# datadir = "../Cascade/tests/cascade_test_human_inputs/"


def run_cascade_roadway_dynamics():
    # Roadway width drowned at 178 years, 20.0% of road borders water
    cascade = Cascade(
        str(BMI_DATA_DIR) + "/",
        # datadir,
        name="test_roadway_relocation",
        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
        dune_file="pathways-dunes.npy",
        parameter_file="roadway-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=1,
        time_step_count=NT,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=True,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
        # m MHW, average of NC-12 is 1.3 m NAVD88, berm ele is 1.4 m MHW (1.9 m NAVD88)
        road_ele=1.2,
        road_width=20,  # m
        road_setback=20,  # m
        dune_design_elevation=3.2,  # m MHW, rebuild to 2 m dune above the roadway
        # m MHW, allow dune to erode down to 0.5 m above the roadway, v1 = 2.7 m
        dune_minimum_elevation=1.7,
    )

    for _ in range(NT - 1):
        cascade.update()
        if cascade.b3d_break:
            break

    return cascade


def run_cascade_nourishment_dynamics():
    iB3D = 0
    total_time = 100

    cascade = Cascade(
        str(BMI_DATA_DIR) + "/",
        name="test_shoreline_erosion",
        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="b3d_pt45_8750yrs_low-elevations.csv",
        dune_file="pathways-dunes.npy",
        parameter_file="nourishment-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.007,
        sea_level_rise_constant=True,
        background_erosion=-1.0,  # m/yr
        alongshore_section_count=1,
        time_step_count=total_time,
        min_dune_growth_rate=0.25,  # average is 0.45, a low dune growth rate
        max_dune_growth_rate=0.65,
        dune_design_elevation=3.7,  # m MHW, rebuild to 2 m dune above the berm
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=True,
        community_economics_module=False,  # no community dynamics
        nourishment_interval=10,  # yrs
        nourishment_volume=100,
        overwash_filter=40,
        overwash_to_dune=10,
    )

    # Loop for 50 years at a 10 year interval, 100 m^3/m and then 50 years at a 20
    # year interval with 300 m^3/m
    nt = 50
    for _ in range(nt - 1):
        cascade.update()
        if cascade.b3d_break:
            break

    # during the CASCADE initialization, the nourishment interval and volume is
    # specified individually for each barrier3d alongshore cell; so to update these
    # values, we need to specify which barrier3d cell we want to modify
    # (here, we only have one cell)
    cascade.nourishment_interval[iB3D] = 20  # increase to 20 years
    cascade.nourishment_volume[iB3D] = 300  # increase to 300 m^3/m

    for _ in range(nt):
        cascade.update()
        if cascade.b3d_break:
            break

    return cascade


CASCADE_ROADWAY_OUTPUT = run_cascade_roadway_dynamics()
CASCADE_NOURISHMENT_OUTPUT = run_cascade_nourishment_dynamics()


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
    # should equal 1 vertical m of overwash removed per road cell (20), summed across
    # the domain [m^3]
    assert road_overwash_removal == np.sum(np.zeros([20]) + (1 * 20))
    # should equal 2 vertical m of accretion
    assert_array_almost_equal(new_dune_domain, np.zeros([20, 10]) + 2)

    # now test that decameters works properly; bulldoze 0.1 vertical decameter of
    # sand from a 2 decameter wide roadway
    # units of decameters in X (cross-shore) Y (alongshore) Z (vertical)
    xyz_interior_grid = np.zeros([100, 20]) + 3.0 / 10
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
    # should equal 0.1 vertical decameters of overwash removed per road cell (2),
    # summed across the domain [dm^3]
    assert_array_almost_equal(
        road_overwash_removal, np.floor(np.sum(np.zeros([20]) + (0.1 + 0.1)))
    )
    # should equal 0.1 vertical decameter of accretion along each dune cell
    assert_array_almost_equal(new_dune_domain, (np.zeros([20, 2]) + 1 / 10))

    # NOTE: need to make a test that shows unequal overwash volume removed per
    # grid cell and placement on adjacent dunes


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

    # only last grid cell has dune growth parameter go to zero (dune height above
    # dune max)
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

    # check the scenario where at the previous time step growth rate was zero but now
    # we want to change it back to 0.7 because the dune height dropped below Dmax
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
    check overwash volumes removed (includes negative values, which can happen)
    """

    # CASE 1: uniform positive overwash
    barrier_overwash_removed = np.zeros(3)
    new_dunes = []
    expected_dunes = []

    post_xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    pre_xyz_interior_grid = post_xyz_interior_grid - 1.0
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    # remove 40% of 1 m^3 of overwash in each cell,
    (
        new_yxz_dune_domain,  # [m], should be zero
        new_xyz_interior_domain,  # [m], should be 2.6 m depth in each cell
        barrier_overwash_removed[0],  # [m^3], should be 0.4 m^3 * 100 * 20 (800 m^3)
        new_ave_interior_height,  # [m] for the rest
        beach_width_new,
        x_s_new,
        s_sf_new,
    ) = filter_overwash(
        overwash_filter=40,  # remove 40% of overwash deposit
        # leave the rest of the overwash deposit, don't bulldoze to dunes
        overwash_to_dune=0,
        post_storm_xyz_interior_grid=post_xyz_interior_grid,
        pre_storm_xyz_interior_grid=pre_xyz_interior_grid,
        post_storm_yxz_dune_grid=yxz_dune_grid,
        # maximum height of artificial dunes in meters
        artificial_maximum_dune_height=5,
        sea_level=0,  # the rest of the variables are in meters
        barrier_length=20,
        x_s=200,
        x_t=5,
        beach_width=20,
        shoreface_depth=8,
        dune_spread_equal=False,
    )
    new_dunes.append(new_yxz_dune_domain)
    expected_dunes.append(np.zeros([20, 10]))

    # CASE 2: overwash only deposited in some areas
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
        new_yxz_dune_domain,  # [m], should be zero
        # [m], should be 2.6 m depth in each cell except for row 0, which should
        # be 2.0 m
        new_xyz_interior_domain,
        barrier_overwash_removed[1],  # [m^3], should be 0.4 m^3 * 99 * 20 (792 m^3)
        new_ave_interior_height,  # [m] for the rest
        beach_width_new,
        x_s_new,
        s_sf_new,
    ) = filter_overwash(
        overwash_filter=40,  # remove 40% of overwash deposit
        overwash_to_dune=0,
        post_storm_xyz_interior_grid=post_xyz_interior_grid,
        pre_storm_xyz_interior_grid=pre_xyz_interior_grid,
        post_storm_yxz_dune_grid=yxz_dune_grid,
        artificial_maximum_dune_height=5,  # max height of artificial dunes in meters
        sea_level=0,  # the rest of the variables are in meters
        barrier_length=20,
        x_s=200,
        x_t=5,
        beach_width=20,
        shoreface_depth=8,
        dune_spread_equal=False,
    )
    new_dunes.append(new_yxz_dune_domain)
    expected_dunes.append(np.zeros([20, 10]))

    # CASE 3: first row is eroded (channelized) due to overwash, the rest of the
    # grid sees deposition
    post_xyz_interior_grid = (
        np.zeros([100, 20]) + 3.0
    )  # units of meters in X (cross-shore) Y (alongshore) Z (vertical)
    pre_xyz_interior_grid = post_xyz_interior_grid - 1.0
    post_xyz_interior_grid[0, :] = 1.0
    yxz_dune_grid = np.zeros([20, 10])  # dune domain is 10 cells (m) wide
    # remove 50% total of 1 m^3 of overwash in each cell with deposition, and add
    # 50% of 1 m^3 for the row eroded
    # total overwash removed = (0.5 m^3 * 20 cells/row * 99 rows)
    #   - (0.5 m^3 * 20 cells/row * 1 row) = 980 m^3
    # to dune, only 1/5 of that volume == 196 m^3 / (10 * 20) added to each
    # cell == 0.98 m
    (
        new_yxz_dune_domain,
        new_xyz_interior_domain,
        barrier_overwash_removed[2],
        new_ave_interior_height,  # [m] for the rest
        beach_width_new,
        x_s_new,
        s_sf_new,
    ) = filter_overwash(
        overwash_filter=40,  # remove 40% of overwash deposit
        overwash_to_dune=10,
        post_storm_xyz_interior_grid=post_xyz_interior_grid,
        pre_storm_xyz_interior_grid=pre_xyz_interior_grid,
        post_storm_yxz_dune_grid=yxz_dune_grid,
        artificial_maximum_dune_height=5,  # max height of artificial dunes in meters
        sea_level=0,  # the rest of the variables are in meters
        barrier_length=20,
        x_s=200,
        x_t=5,
        beach_width=20,
        shoreface_depth=8,
        dune_spread_equal=False,
    )
    new_dunes.append(new_yxz_dune_domain)
    expected_dunes.append(np.zeros([20, 10]) + 0.98)  # m

    assert (np.round(barrier_overwash_removed) == [800, 792, 980]).all()
    assert all([a == b] for a, b in zip(new_dunes, expected_dunes))
    # assert_array_almost_equal(barrier_overwash_removed, [800, 792, 990])


def test_shoreline_migration():
    """
    As a check on the dynamics in Barrier3D, here we want to see if the dunes
    migrate when the beach width goes to zero and the shoreline surpasses a full
    cell width (10 m). Indeed, dunes only migrate when humans allow them to
    (beach width = 0), and when the shoreline moves a full cell width -- in year
    58 and 65. Note that if one was to check the `post_storm_x_s`, they would find
    that the dunes actually migrated at 57.5 since dune migration occurs prior
    to any human modifications.
    """

    iB3D = 0

    frac_grid_cell = np.array(CASCADE_NOURISHMENT_OUTPUT.barrier3d[iB3D]._x_s_TS) % 1
    diff = np.hstack([0, np.diff(frac_grid_cell)])
    diff[
        ~np.array(CASCADE_NOURISHMENT_OUTPUT.nourishments[iB3D]._dune_migration_on)
    ] = 0
    shoreline_transgressed = diff < 0
    dunes_migrated = CASCADE_NOURISHMENT_OUTPUT.barrier3d[iB3D]._ShorelineChangeTS < 0

    assert np.all(shoreline_transgressed == dunes_migrated)


def test_shoreline_road_relocation():
    """
    Does the roadway relocate when the shoreline eats up the dune?
    """
    iB3D = 0

    dunes_migrated = CASCADE_ROADWAY_OUTPUT.barrier3d[iB3D]._ShorelineChangeTS < 0
    road_relocated = CASCADE_ROADWAY_OUTPUT.roadways[iB3D]._road_relocated_TS > 0
    road_setback_TS = CASCADE_ROADWAY_OUTPUT.roadways[iB3D]._road_setback_TS

    diff_road_setback = np.hstack([0, np.diff(road_setback_TS)])
    road_relocated_based_on_setback = diff_road_setback > 0

    assert np.all(road_relocated_based_on_setback == road_relocated)
    assert np.all(dunes_migrated[road_relocated])
