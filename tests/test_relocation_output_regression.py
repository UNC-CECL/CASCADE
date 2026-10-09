"""Regression comparison for original and current no-request saved runs.

Local Pea Island outputs are detected automatically. ``CASCADE_ORIGINAL_NPZ``
and ``CASCADE_RELOCATION_NPZ`` can override them for other completed runs.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ORIGINAL_NPZ = (
    WORKSPACE_ROOT
    / "CASCADE/output/baseline_original_roadway_manager"
    / "PEA_1992_2007_OriginalRoadwayManager_Hs2p0_Berm1p7"
    / "PEA_1992_2007_OriginalRoadwayManager_Hs2p0_Berm1p7.npz"
)
DEFAULT_RELOCATION_NPZ = (
    WORKSPACE_ROOT
    / "CASCADE/output/roadway_nourishment/no_requests"
    / "PEA_1992_2007_RoadwayIndependentNourishment_NoRequests_Hs2p0_Berm1p7"
    / "PEA_1992_2007_RoadwayIndependentNourishment_NoRequests_Hs2p0_Berm1p7.npz"
)

BARRIER_OUTPUTS = (
    pytest.param("x_s_TS", id="shoreline-position-x_s_TS"),
    pytest.param("DuneDomain", id="dune-elevation-DuneDomain"),
    pytest.param("DomainTS", id="interior-history-DomainTS"),
    pytest.param("InteriorDomain", id="final-interior-InteriorDomain"),
    pytest.param("InteriorWidth_AvgTS", id="barrier-width-InteriorWidth_AvgTS"),
    pytest.param("ShorelineChangeTS", id="dune-migration-ShorelineChangeTS"),
    pytest.param("QowTS", id="overwash-flux-QowTS"),
    pytest.param("QsfTS", id="shoreface-flux-QsfTS"),
    pytest.param("h_b_TS", id="barrier-height-h_b_TS"),
    pytest.param("s_sf_TS", id="shoreface-slope-s_sf_TS"),
    pytest.param("x_b_TS", id="back-barrier-position-x_b_TS"),
    pytest.param("x_t_TS", id="shoreface-toe-position-x_t_TS"),
)

ORIGINAL_ROADWAY_OUTPUTS = (
    pytest.param("road_setback_TS", id="road-setback-road_setback_TS"),
    pytest.param("road_width_TS", id="road-width-road_width_TS"),
    pytest.param("road_ele_TS", id="road-elevation-road_ele_TS"),
    pytest.param(
        "dune_design_elevation_TS",
        id="dune-design-elevation-dune_design_elevation_TS",
    ),
    pytest.param(
        "dune_minimum_elevation_TS",
        id="dune-minimum-elevation-dune_minimum_elevation_TS",
    ),
    pytest.param("dunes_rebuilt_TS", id="dune-rebuild-events-dunes_rebuilt_TS"),
    pytest.param("road_relocated_TS", id="road-relocation-events-road_relocated_TS"),
    pytest.param(
        "rebuild_dune_volume_TS",
        id="dune-rebuild-volume-rebuild_dune_volume_TS",
    ),
    pytest.param(
        "road_overwash_volume",
        id="road-overwash-volume-road_overwash_volume",
    ),
    pytest.param("percent_below_min", id="dunes-below-minimum-percent_below_min"),
    pytest.param("growth_params", id="dune-growth-parameters-growth_params"),
    pytest.param("post_storm_dunes", id="post-storm-dunes-post_storm_dunes"),
    pytest.param(
        "post_storm_interior",
        id="post-storm-interior-post_storm_interior",
    ),
    pytest.param(
        "post_storm_ave_interior_height",
        id="post-storm-height-post_storm_ave_interior_height",
    ),
)


def attr(obj, name):
    """Read a public property or its saved private attribute."""

    if hasattr(obj, name):
        return getattr(obj, name)
    return getattr(obj, f"_{name}")


def result_path(environment_variable, default_path):
    """Use an environment override or the existing local Pea Island result."""

    value = os.environ.get(environment_variable)
    path = Path(value).expanduser().resolve() if value else default_path
    if not path.is_file():
        if value:
            pytest.fail(f"{environment_variable} does not identify a file: {path}")
        pytest.skip(
            f"Default saved output is unavailable; set {environment_variable}: {path}"
        )
    return path


def load_cascade(path):
    """Load one saved CASCADE object."""

    with np.load(path, allow_pickle=True) as archive:
        return archive["cascade"].item()


def assert_saved_equal(original, relocation, location):
    """Compare saved scalars, arrays, and lists exactly, including NaNs."""

    if original is None or relocation is None:
        assert original is relocation, location
        return

    if isinstance(original, np.ndarray):
        np.testing.assert_array_equal(original, relocation, err_msg=location)
        return

    if isinstance(original, (list, tuple)):
        assert len(original) == len(relocation), location
        for index, (original_item, relocation_item) in enumerate(
            zip(original, relocation)
        ):
            assert_saved_equal(
                original_item,
                relocation_item,
                f"{location}[{index}]",
            )
        return

    if isinstance(original, float) and np.isnan(original):
        assert isinstance(relocation, float) and np.isnan(relocation), location
        return

    assert original == relocation, location


@pytest.fixture(scope="module")
def saved_runs():
    original = load_cascade(result_path("CASCADE_ORIGINAL_NPZ", DEFAULT_ORIGINAL_NPZ))
    relocation = load_cascade(
        result_path("CASCADE_RELOCATION_NPZ", DEFAULT_RELOCATION_NPZ)
    )
    return original, relocation


def test_domain_counts_match(saved_runs):
    original, relocation = saved_runs
    assert len(original.barrier3d) == len(relocation.barrier3d)
    assert len(original.roadways) == len(relocation.roadways)


@pytest.mark.parametrize("output_name", BARRIER_OUTPUTS)
def test_barrier_output_matches_original(saved_runs, output_name):
    """Compare one physical Barrier3D output across every domain."""

    original, relocation = saved_runs
    for domain_index, (original_barrier, relocation_barrier) in enumerate(
        zip(original.barrier3d, relocation.barrier3d)
    ):
        assert_saved_equal(
            attr(original_barrier, output_name),
            attr(relocation_barrier, output_name),
            f"barrier3d[{domain_index}].{output_name}",
        )


@pytest.mark.parametrize("output_name", ORIGINAL_ROADWAY_OUTPUTS)
def test_original_roadway_output_matches(saved_runs, output_name):
    """Compare one original RoadwayManager output across every domain."""

    original, relocation = saved_runs
    for domain_index, (original_roadway, relocation_roadway) in enumerate(
        zip(original.roadways, relocation.roadways)
    ):
        assert_saved_equal(
            attr(original_roadway, output_name),
            attr(relocation_roadway, output_name),
            f"roadways[{domain_index}].{output_name}",
        )


def test_triggered_relocation_matches_original_relocation(saved_runs):
    """Natural relocation flags must reproduce the original relocation record."""

    original, relocation = saved_runs
    for domain_index, (original_roadway, relocation_roadway) in enumerate(
        zip(original.roadways, relocation.roadways)
    ):
        np.testing.assert_array_equal(
            relocation_roadway.triggered_relocation_TS,
            np.asarray(attr(original_roadway, "road_relocated_TS")) != 0,
            err_msg=f"roadways[{domain_index}].triggered_relocation_TS",
        )


def test_no_incomplete_relocation_without_requests(saved_runs):
    _, relocation = saved_runs
    for relocation_roadway in relocation.roadways:
        assert not relocation_roadway.relocation_incomplete_TS.any()


def test_no_historical_relocation_request_without_requests(saved_runs):
    _, relocation = saved_runs
    for relocation_roadway in relocation.roadways:
        assert not relocation_roadway.historical_relocation_requested_TS.any()


def test_no_forced_relocation_without_requests(saved_runs):
    _, relocation = saved_runs
    for relocation_roadway in relocation.roadways:
        assert not relocation_roadway.forced_relocation_TS.any()
