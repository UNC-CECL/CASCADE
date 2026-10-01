"""Roadway functionality regressions, kept separate from the original tests."""

from types import SimpleNamespace

import numpy as np
import pytest

from cascade import Cascade
from cascade.roadway_manager import RoadwayManager
from test_full_barrier3d_roadway_nourishment_relocation import (
    CASCADE_PARAMETERS,
    prepare_inputs,
)

def relocation_test_barrier(average_width_m=100.0, dune_migration_m=0.0):
    """Create the minimum Barrier3D state needed for a manager update."""

    states = 4
    alongshore = 5
    interior = np.full((12, alongshore), 0.2)
    dunes = np.full((states, alongshore, 2), 0.4)
    domain_ts = np.empty(states, dtype=object)
    for index in range(states):
        domain_ts[index] = interior.copy()
    return SimpleNamespace(
        time_index=2,
        growthparam=np.full((1, alongshore), 0.5),
        InteriorDomain=interior,
        DuneDomain=dunes,
        h_b_TS=[0.2, 0.2],
        InteriorWidth_AvgTS=[average_width_m / 10.0],
        ShorelineChangeTS=np.array([0.0, dune_migration_m / 10.0, 0.0, 0.0]),
        RSLR=np.zeros(states),
        BermEl=0.1,
        SL=0.0,
        Dmax=0.5,
        DomainTS=domain_ts,
        x_s=100.0,
        x_t=0.0,
        x_s_TS=[100.0, 100.0],
        x_b_TS=[114.0, 114.0],
        s_sf_TS=[0.01, 0.01],
        DShoreface=1.0,
        dune_migration_on=True,
        SCRagg=np.zeros(states),
    )

def test_default_20_m_setback_trigger_permits_rebuilding_below_trigger():
    roadway = RoadwayManager(
        road_width=10,
        road_setback=10,
        initial_dune_design_elevation=3.0,
        initial_dune_minimum_elevation=2.2,
        time_step_count=4,
    )
    barrier = relocation_test_barrier()
    barrier.DuneDomain[1, :, :] = 0.0

    roadway.update(barrier, trigger_dune_knockdown=False)

    event_index = barrier.time_index - 1
    assert roadway._dunes_rebuilt_TS[event_index]
    assert not roadway.road_dune_rebuild_disabled_TS[event_index]
    np.testing.assert_array_equal(barrier.growthparam, np.full((1, 5), 0.5))

def test_setback_rebuild_rule_uses_advisor_specified_greater_than_condition():
    for road_setback, expected_disabled in (
        (19.999, False),
        (20.0, False),
        (20.001, True),
    ):
        roadway = RoadwayManager(
            road_width=10,
            road_setback=road_setback,
            initial_dune_design_elevation=3.0,
            initial_dune_minimum_elevation=2.2,
            time_step_count=4,
            road_setback_trigger=20.0,
        )
        barrier = relocation_test_barrier()
        barrier.DuneDomain[1, :, :] = 0.0

        roadway.update(barrier, trigger_dune_knockdown=False)

        event_index = barrier.time_index - 1
        assert (
            bool(roadway.road_dune_rebuild_disabled_TS[event_index])
            == expected_disabled
        )
        assert bool(roadway._dunes_rebuilt_TS[event_index]) == (not expected_disabled)
        np.testing.assert_array_equal(barrier.growthparam, np.full((1, 5), 0.5))

def test_setback_rebuild_rule_stops_after_relocation_above_20_m():
    roadway = RoadwayManager(
        road_width=10,
        road_setback=10,
        initial_dune_design_elevation=3.0,
        initial_dune_minimum_elevation=2.2,
        time_step_count=4,
        road_setback_trigger=20.0,
    )
    barrier = relocation_test_barrier()
    barrier.DuneDomain[1:, :, :] = 0.0

    roadway.update(barrier, trigger_dune_knockdown=False)
    assert not roadway.road_dune_rebuild_disabled_TS[1]
    assert roadway._dunes_rebuilt_TS[1]
    np.testing.assert_array_equal(barrier.growthparam, np.full((1, 5), 0.5))

    roadway.request_historical_relocation(30.0)
    barrier.time_index = 3
    roadway.update(barrier, trigger_dune_knockdown=False)

    np.testing.assert_array_equal(barrier.growthparam, np.full((1, 5), 0.5))
    assert roadway._road_setback_TS[2] == pytest.approx(30.0)
    assert roadway.forced_relocation_TS[2]
    assert roadway.road_dune_rebuild_disabled_TS[2]
    assert not roadway._dunes_rebuilt_TS[2]
    assert roadway._dune_design_elevation_TS[2] == pytest.approx(
        roadway._road_ele_TS[2] + 1.3
    )


def test_roadway_nourishment_keeps_original_cascade_back_barrier_formula(tmp_path):
    prepare_inputs(tmp_path)
    parameters = {**CASCADE_PARAMETERS, "time_step_count": 4}
    model = Cascade(str(tmp_path), name="original_back_barrier_geometry", **parameters)
    initial_beach_width = model._initial_beach_width[0]
    model.nourish_now = [1]
    saw_changed_beach_width = False
    for _ in range(3):
        model.update()
        assert not model.b3d_break
        barrier = model.barrier3d[0]
        expected = (
            barrier.x_s
            + barrier.InteriorWidth_AvgTS[-1]
            + np.size(barrier.DuneDomain, 2)
            + initial_beach_width / 10
        )
        assert barrier.x_b_TS[-1] == pytest.approx(expected)
        current_width = model.roadways[0].beach_width[barrier.time_index - 1]
        saw_changed_beach_width |= not np.isclose(current_width, initial_beach_width)
    assert saw_changed_beach_width, "Exercise the formula after nourishment changes beach width"
