from pathlib import Path

import numpy as np
from barrier3d import Barrier3dBmi
from numpy.testing import assert_array_almost_equal

from cascade.cascade import Cascade

BMI_DATA_DIR = Path(__file__).parent / "cascade_test_versions_inputs"
# datadir = "../Cascade/tests/cascade_test_versions_inputs/"
NT = 30


def run_cascade_no_human_dynamics():
    cascade = Cascade(
        str(BMI_DATA_DIR) + "/",
        # datadir,
        name="test_coupled_dune_migration",
        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
        dune_file="pathways-dunes.npy",
        parameter_file="barrier3d-default-parameters.yaml",
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
        roadway_management_module=False,
        alongshore_transport_module=False,
        beach_nourishment_module=False,
        community_economics_module=False,  # no community dynamics
    )

    for time_step in range(NT - 1):
        cascade.update()
        if cascade.b3d_break:
            break

    return cascade


def initialize_cascade_no_human_dynamics_ast():
    cascade = Cascade(
        name="test_initialize",
        datadir=str(BMI_DATA_DIR) + "/",
        storm_file="StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
        elevation_file="b3d_pt75_3284yrs_low-elevations.csv",
        dune_file="pathways-dunes.npy",
        parameter_file="ast-barrier3d-parameters.yaml",
        wave_height=1,
        wave_period=7,
        wave_asymmetry=0.8,
        wave_angle_high_fraction=0.2,
        sea_level_rise_rate=0.004,
        sea_level_rise_constant=True,
        background_erosion=0.0,
        alongshore_section_count=6,  # 6 barrier 3d models connected with alongshore sediment transport
        time_step_count=3,
        min_dune_growth_rate=0.55,
        max_dune_growth_rate=0.95,  # rave = 0.75
        num_cores=1,
        roadway_management_module=False,
        alongshore_transport_module=True,
        beach_nourishment_module=False,
        community_economics_module=False,
    )

    return cascade


CASCADE_OUTPUT = run_cascade_no_human_dynamics()
CASCADE_AST_MODEL = initialize_cascade_no_human_dynamics_ast()


def test_barrier3d_versions():
    """
    Check that the barrier3d output in cascade matches the bmi version(s)
    """

    # Barrier3D BMI Version (2.0 and beyond): create an instance of the new BMI class, which is the model
    barrier3d = Barrier3dBmi()
    input_file = "barrier3d-default-parameters.yaml"
    barrier3d.initialize(str(BMI_DATA_DIR / input_file))

    # increase time step
    for time_step in range(NT - 1):
        barrier3d.update()

    assert np.all(
        barrier3d._model._DuneDomain == CASCADE_OUTPUT._barrier3d[0]._DuneDomain
    )  # dune height over time
    assert np.all(
        barrier3d._model._x_s_TS == CASCADE_OUTPUT._barrier3d[0]._x_s_TS
    )  # shoreline change time series
    assert np.all(
        barrier3d._model._ShorelineChangeTS
        == CASCADE_OUTPUT._barrier3d[0]._ShorelineChangeTS
    )  # dune migration time series
    assert np.all(
        barrier3d._model._s_sf_TS == CASCADE_OUTPUT._barrier3d[0]._s_sf_TS
    )  # shoreface slope time series
    assert np.all(barrier3d._model._QowTS == CASCADE_OUTPUT._barrier3d[0]._QowTS)


def test_shoreline_dune_migration():
    """
    As a check on the dynamics in Barrier3D, here we want to see if the dunes migrate correctly for natural simulations,
    when the shoreline surpasses a full cell width (10 m). Barrier3D was originally written for a shoreface toe to start
    at x=0 dam. When coupled with BRIE, the shoreface starts at a non-integer value, so I had to modify the dune
    migration algorithmn within Barrier3D to account for this.
    """
    iB3D = 0
    frac_grid_cell = np.array(CASCADE_OUTPUT.barrier3d[iB3D]._x_s_TS) % 1
    diff = np.hstack([0, np.diff(frac_grid_cell)])
    shoreline_transgressed = diff < 0
    dunes_migrated = CASCADE_OUTPUT.barrier3d[iB3D]._ShorelineChangeTS < 0

    assert np.all(shoreline_transgressed == dunes_migrated)


def test_initialize():
    """
    check that the initial shoreface toe and shoreline are correct between the barrier3d and brie models for a cascade
    model with 6 barrier3d domains
    """

    x_t_TS, x_s_TS, x_b_TS, h_b_TS = ([] for _ in range(4))
    for iB3D in range(CASCADE_AST_MODEL.brie.ny):
        x_t_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].x_t_TS[0] * 10)
        x_s_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].x_s_TS[0] * 10)
        x_b_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].x_b_TS[0] * 10)
        h_b_TS.append(CASCADE_AST_MODEL.barrier3d[iB3D].h_b_TS[0] * 10)

    dt = CASCADE_AST_MODEL.brie.x_t - x_t_TS
    ds = CASCADE_AST_MODEL.brie.x_s - x_s_TS  # this isn't always zero due to rounding
    db = CASCADE_AST_MODEL.brie.x_b - x_b_TS
    dh = CASCADE_AST_MODEL.brie.h_b - h_b_TS

    assert_array_almost_equal([dt, ds, db, dh], np.zeros([4, 6]))


def test_shoreline_variable_exchange_ast():
    """
    check that the brie and barrier3d shoreline, shoreface toe, and back-barrier shorelines are equivalent with each
    update
    """
    CASCADE_AST_MODEL.update()
    CASCADE_AST_MODEL.update()

    x_t_TS, x_s_TS, x_b_TS, h_b_TS = ([] for _ in range(4))
    for iB3D in range(CASCADE_AST_MODEL.brie.ny):
        x_t_TS.append(np.array(CASCADE_AST_MODEL.barrier3d[iB3D].x_t_TS) * 10)
        x_s_TS.append(np.array(CASCADE_AST_MODEL.barrier3d[iB3D].x_s_TS) * 10)
        x_b_TS.append(np.array(CASCADE_AST_MODEL.barrier3d[iB3D].x_b_TS) * 10)
        h_b_TS.append(np.array(CASCADE_AST_MODEL.barrier3d[iB3D].h_b_TS) * 10)

    dt = CASCADE_AST_MODEL.brie._x_t_save - np.array(x_t_TS).astype(int)
    db = CASCADE_AST_MODEL.brie._x_b_save - np.array(x_b_TS).astype(int)
    ds = CASCADE_AST_MODEL.brie._x_s_save - np.array(x_s_TS).astype(int)
    dh = CASCADE_AST_MODEL.brie._h_b_save - np.array(
        h_b_TS
    )  # this isn't always zero; rounding error

    assert_array_almost_equal([dt, ds, db, dh], np.zeros([4, 6, 3]))
