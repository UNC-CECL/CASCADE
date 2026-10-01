"""Full Barrier3D nourishment lifecycle for both management modules.

Unlike the small manager-only test, this file runs CASCADE's normal update loop
and the real ``barrier3d.barrier3d.Barrier3d`` model.  It uses the established
CASCADE test elevation, dune, growth, and parameter inputs.  StormStart is moved
past the end of this controlled run so no storm can overtop the dunes.  One
100 m3/m ``nourish_now`` request is supplied to each scenario; all subsequent
shoreline change and dune-grid migration are calculated by Barrier3D.
"""

from __future__ import annotations

import csv
import json
import shutil
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pytest
import yaml
from matplotlib.animation import FuncAnimation, PillowWriter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

from barrier3d import Barrier3d  # noqa: E402
from cascade import Cascade  # noqa: E402


INPUT_ROOT = SOURCE_ROOT / "tests" / "test_human_dynamics"
INPUT_FILES = (
    "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "b3d_pt45_8750yrs_low-elevations.csv",
    "pathways-dunes.npy",
    "growthparam_1000dam.npy",
    "nourishment-parameters.yaml",
)

TIME_STEP_COUNT = 45
NO_STORM_START = TIME_STEP_COUNT + 1
NOURISHMENT_VOLUME_M3_PER_M = 100.0
NOURISHMENT_TIME_INDEX = 1

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment"
    / "full_barrier3d_nourishment_test"
)
FIGURE_PATH = OUTPUT_DIR / "actual_barrier3d_nourishment_lifecycle.png"
GIF_PATH = OUTPUT_DIR / "actual_barrier3d_nourishment_lifecycle_slow.gif"
CSV_PATH = OUTPUT_DIR / "actual_barrier3d_nourishment_lifecycle.csv"
MANIFEST_PATH = OUTPUT_DIR / "actual_barrier3d_nourishment_test_manifest.json"
NPZ_PATHS = {
    "BeachDuneManager": OUTPUT_DIR / "actual_barrier3d_beach_dune_manager.npz",
    "RoadwayManager": OUTPUT_DIR / "actual_barrier3d_roadway_manager.npz",
}

COMMON_CASCADE_PARAMETERS = {
    "storm_file": "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "elevation_file": "b3d_pt45_8750yrs_low-elevations.csv",
    "dune_file": "pathways-dunes.npy",
    "parameter_file": "nourishment-parameters.yaml",
    "wave_height": 1,
    "wave_period": 7,
    "wave_asymmetry": 0.8,
    "wave_angle_high_fraction": 0.2,
    "sea_level_rise_rate": 0.007,
    "sea_level_rise_constant": True,
    "background_erosion": -1.0,
    "alongshore_section_count": 1,
    "time_step_count": TIME_STEP_COUNT,
    "min_dune_growth_rate": 0.25,
    "max_dune_growth_rate": 0.65,
    "num_cores": 1,
    "alongshore_transport_module": False,
    "community_economics_module": False,
    "nourishment_interval": None,
    "nourishment_volume": NOURISHMENT_VOLUME_M3_PER_M,
    "overwash_filter": 40,
    "overwash_to_dune": 10,
    "outwash_module": False,
    "road_ele": 1.2,
    "road_width": 20,
    "road_setback": 20,
    "dune_design_elevation": 3.2,
    "dune_minimum_elevation": 1.7,
}


def prepare_inputs(run_directory: Path) -> None:
    """Copy repository inputs and defer storms beyond the controlled run."""

    run_directory.mkdir(parents=True)
    for filename in INPUT_FILES:
        shutil.copy2(INPUT_ROOT / filename, run_directory / filename)

    parameter_path = run_directory / "nourishment-parameters.yaml"
    with parameter_path.open() as stream:
        parameters = yaml.safe_load(stream)
    parameters["StormStart"] = NO_STORM_START
    with parameter_path.open("w") as stream:
        yaml.safe_dump(parameters, stream, sort_keys=False)


def run_actual_scenario(run_directory: Path, manager_name: str) -> Cascade:
    """Run one real Barrier3D scenario through all requested timesteps."""

    prepare_inputs(run_directory)
    parameters = dict(COMMON_CASCADE_PARAMETERS)
    parameters["roadway_management_module"] = manager_name == "RoadwayManager"
    parameters["beach_nourishment_module"] = manager_name == "BeachDuneManager"

    cascade = Cascade(
        str(run_directory),
        name=f"actual_barrier3d_{manager_name.lower()}",
        **parameters,
    )
    cascade.nourish_now = [1]
    for _ in range(TIME_STEP_COUNT - 1):
        cascade.update()
        assert not cascade.b3d_break, "Barrier3D drowned before the test completed"

    return cascade


def selected_manager(cascade: Cascade, manager_name: str):
    if manager_name == "BeachDuneManager":
        return cascade.nourishments[0]
    return cascade.roadways[0]


def scenario_series(cascade: Cascade, manager_name: str) -> dict[str, np.ndarray]:
    barrier = cascade.barrier3d[0]
    manager = selected_manager(cascade, manager_name)
    state_count = len(barrier.x_s_TS)
    shoreline_change = np.asarray(barrier.ShorelineChangeTS[:state_count], dtype=float)
    return {
        "time_index": np.arange(state_count),
        "shoreline_m": np.asarray(barrier.x_s_TS, dtype=float) * 10,
        "beach_width_m": np.asarray(manager.beach_width[:state_count], dtype=float),
        "dune_migration_on": np.asarray(
            manager._dune_migration_on[:state_count], dtype=float
        ),
        "shoreline_change_cells": shoreline_change,
        "cumulative_landward_dune_grid_m": np.cumsum(-shoreline_change) * 10,
        "overwash_flux_m3_per_m": np.asarray(barrier.QowTS, dtype=float),
        "storm_count": np.asarray(barrier._StormCount, dtype=int),
        "nourishment": np.asarray(
            manager._nourishment_TS[:state_count], dtype=bool
        ),
        "nourishment_volume_m3_per_m": np.asarray(
            manager._nourishment_volume_TS[:state_count], dtype=float
        ),
    }


@pytest.fixture(scope="module")
def actual_runs(tmp_path_factory):
    root = tmp_path_factory.mktemp("actual_barrier3d_nourishment")
    return {
        manager_name: run_actual_scenario(root / directory, manager_name)
        for manager_name, directory in (
            ("BeachDuneManager", "beach_dune_manager"),
            ("RoadwayManager", "roadway_manager"),
        )
    }


@pytest.fixture(scope="module")
def actual_series(actual_runs):
    return {
        manager_name: scenario_series(cascade, manager_name)
        for manager_name, cascade in actual_runs.items()
    }


@pytest.mark.parametrize("manager_name", ["BeachDuneManager", "RoadwayManager"])
def test_real_barrier3d_runs_with_zero_storms(actual_runs, actual_series, manager_name):
    """Confirm the actual model ran and the controlled storm count stayed zero."""

    barrier = actual_runs[manager_name].barrier3d[0]
    series = actual_series[manager_name]
    assert isinstance(barrier, Barrier3d)
    assert len(series["time_index"]) == TIME_STEP_COUNT
    np.testing.assert_array_equal(series["storm_count"], 0)
    np.testing.assert_allclose(series["overwash_flux_m3_per_m"], 0.0)


@pytest.mark.parametrize("manager_name", ["BeachDuneManager", "RoadwayManager"])
def test_one_real_nourishment_request_is_applied(actual_series, manager_name):
    """The same single nourish_now request and volume must be recorded."""

    series = actual_series[manager_name]
    event_indices = np.flatnonzero(series["nourishment"])
    np.testing.assert_array_equal(event_indices, [NOURISHMENT_TIME_INDEX])
    assert series["nourishment_volume_m3_per_m"][NOURISHMENT_TIME_INDEX] == pytest.approx(
        NOURISHMENT_VOLUME_M3_PER_M
    )
    assert series["beach_width_m"][NOURISHMENT_TIME_INDEX] > series["beach_width_m"][0]
    assert series["dune_migration_on"][NOURISHMENT_TIME_INDEX] == 0


@pytest.mark.parametrize("manager_name", ["BeachDuneManager", "RoadwayManager"])
def test_real_beach_erodes_to_zero_before_migration(actual_series, manager_name):
    """Real model evolution must remove the beach and then permit migration."""

    series = actual_series[manager_name]
    zero_width_indices = np.flatnonzero(np.isclose(series["beach_width_m"], 0.0))
    migration_indices = np.flatnonzero(series["shoreline_change_cells"] < 0)
    assert zero_width_indices.size > 0
    assert migration_indices.size > 0
    first_zero = int(zero_width_indices[0])
    first_migration = int(migration_indices[0])
    assert series["dune_migration_on"][first_zero] == 1
    assert first_migration > first_zero
    assert np.all(series["beach_width_m"][first_zero:] == 0)


@pytest.mark.parametrize("manager_name", ["BeachDuneManager", "RoadwayManager"])
def test_real_barrier3d_applies_discrete_dune_grid_migration(
    actual_runs,
    actual_series,
    manager_name,
):
    """A Barrier3D migration flag must correspond to deleted interior grid rows."""

    barrier = actual_runs[manager_name].barrier3d[0]
    series = actual_series[manager_name]
    migration_indices = np.flatnonzero(series["shoreline_change_cells"] < 0)
    assert migration_indices.size >= 2
    for index in migration_indices:
        moved_cells = int(abs(series["shoreline_change_cells"][index]))
        previous_rows = np.asarray(barrier.DomainTS[index - 1]).shape[0]
        current_rows = np.asarray(barrier.DomainTS[index]).shape[0]
        assert current_rows == previous_rows - moved_cells


def test_road_relocation_does_not_enter_the_controlled_comparison(actual_runs):
    """The established setback must remain uncrossed during this run."""

    cascade = actual_runs["RoadwayManager"]
    roadway = cascade.roadways[0]
    assert not cascade.road_break[0]
    assert not roadway.triggered_relocation_TS.any()
    assert not roadway.relocation_incomplete_TS.any()
    assert not roadway.historical_relocation_requested_TS.any()
    assert not roadway.forced_relocation_TS.any()


def transition_indices(series: dict[str, np.ndarray]) -> tuple[int, int]:
    first_zero = int(np.flatnonzero(np.isclose(series["beach_width_m"], 0.0))[0])
    first_migration = int(np.flatnonzero(series["shoreline_change_cells"] < 0)[0])
    return first_zero, first_migration


def write_csv(series_by_manager: dict[str, dict[str, np.ndarray]]) -> None:
    with CSV_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "manager",
                "time_index",
                "shoreline_m",
                "beach_width_m",
                "dune_migration_on",
                "shoreline_change_cells",
                "cumulative_landward_dune_grid_m",
                "storm_count",
                "overwash_flux_m3_per_m",
                "nourishment",
                "nourishment_volume_m3_per_m",
            ],
        )
        writer.writeheader()
        for manager_name, series in series_by_manager.items():
            for index in series["time_index"]:
                writer.writerow(
                    {
                        key: manager_name if key == "manager" else series[key][index]
                        for key in writer.fieldnames
                    }
                )


def plot_static(series_by_manager: dict[str, dict[str, np.ndarray]]) -> None:
    colors = {"BeachDuneManager": "#238b45", "RoadwayManager": "#d95f0e"}
    figure, axes = plt.subplots(4, 1, figsize=(16, 13), sharex=True)
    fields = (
        ("shoreline_m", "Barrier3D shoreline position", "Position (m)"),
        ("beach_width_m", "Manager-tracked beach width", "Width (m)"),
        (
            "cumulative_landward_dune_grid_m",
            "Actual cumulative Barrier3D dune-grid migration",
            "Landward movement (m)",
        ),
    )
    for manager_name, series in series_by_manager.items():
        x = series["time_index"]
        for axis, (field, title, ylabel) in zip(axes[:3], fields):
            axis.plot(x, series[field], linewidth=2, label=manager_name, color=colors[manager_name])
            axis.set_title(title)
            axis.set_ylabel(ylabel)
            axis.grid(alpha=0.25)
        axes[3].step(
            x,
            series["dune_migration_on"],
            where="post",
            linewidth=2,
            label=manager_name,
            color=colors[manager_name],
        )
        migration_indices = np.flatnonzero(series["shoreline_change_cells"] < 0)
        axes[2].scatter(
            migration_indices,
            series["cumulative_landward_dune_grid_m"][migration_indices],
            marker="D",
            s=38,
            color=colors[manager_name],
            edgecolor="black",
            linewidth=0.5,
        )
        first_zero, first_migration = transition_indices(series)
        axes[1].axvline(first_zero, color=colors[manager_name], linestyle=":", alpha=0.8)
        axes[2].axvline(first_migration, color=colors[manager_name], linestyle=":", alpha=0.8)

    axes[1].axhline(0, color="black", linewidth=1)
    axes[3].set_title("Manager-recorded dune-migration permission")
    axes[3].set_ylabel("State")
    axes[3].set_yticks([0, 1], ["OFF", "ON"])
    axes[3].set_ylim(-0.2, 1.2)
    axes[3].set_xlabel("Actual Barrier3D time index")
    axes[3].grid(alpha=0.25)
    for axis in axes:
        axis.axvline(NOURISHMENT_TIME_INDEX, color="#3182bd", linestyle="--", alpha=0.8)
    axes[0].legend(loc="best")
    figure.suptitle(
        "Full Barrier3D nourishment lifecycle — one real domain, no storms\n"
        "Dashed blue: nourish_now; dotted: beach reaches zero / first real grid migration",
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    figure.savefig(FIGURE_PATH, dpi=190, bbox_inches="tight")
    plt.close(figure)


def plot_gif(series_by_manager: dict[str, dict[str, np.ndarray]]) -> None:
    colors = {"BeachDuneManager": "#238b45", "RoadwayManager": "#d95f0e"}
    key_frames = {NOURISHMENT_TIME_INDEX}
    for series in series_by_manager.values():
        key_frames.update(transition_indices(series))
        key_frames.update(np.flatnonzero(series["shoreline_change_cells"] < 0).tolist())
    frame_indices = []
    for index in range(TIME_STEP_COUNT):
        frame_indices.extend([index] * (3 if index in key_frames else 1))

    figure, axes = plt.subplots(2, 2, figsize=(15, 9))

    def update(frame_number):
        current = frame_indices[frame_number]
        for axis in axes.flat:
            axis.clear()
        visible = slice(0, current + 1)
        for manager_name, series in series_by_manager.items():
            x = series["time_index"][visible]
            color = colors[manager_name]
            axes[0, 0].plot(x, series["shoreline_m"][visible], color=color, linewidth=2, label=manager_name)
            axes[0, 1].plot(x, series["beach_width_m"][visible], color=color, linewidth=2, label=manager_name)
            axes[1, 0].step(x, series["dune_migration_on"][visible], where="post", color=color, linewidth=2, label=manager_name)
            axes[1, 1].plot(x, series["cumulative_landward_dune_grid_m"][visible], color=color, linewidth=2, label=manager_name)
            migrated = np.flatnonzero(series["shoreline_change_cells"][: current + 1] < 0)
            axes[1, 1].scatter(
                migrated,
                series["cumulative_landward_dune_grid_m"][migrated],
                marker="D",
                s=35,
                color=color,
                edgecolor="black",
                linewidth=0.4,
            )

        axes[0, 0].set_title("Barrier3D shoreline position")
        axes[0, 0].set_ylabel("Position (m)")
        axes[0, 1].set_title("Beach width")
        axes[0, 1].set_ylabel("Width (m)")
        axes[0, 1].axhline(0, color="black", linewidth=1)
        axes[1, 0].set_title("Dune-migration permission")
        axes[1, 0].set_ylabel("State")
        axes[1, 0].set_yticks([0, 1], ["OFF", "ON"])
        axes[1, 0].set_ylim(-0.2, 1.2)
        axes[1, 1].set_title("Actual Barrier3D dune-grid movement")
        axes[1, 1].set_ylabel("Cumulative landward movement (m)")
        for axis in axes.flat:
            axis.set_xlim(0, TIME_STEP_COUNT - 1)
            axis.axvline(NOURISHMENT_TIME_INDEX, color="#3182bd", linestyle="--", alpha=0.7)
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        axes[1, 0].set_xlabel("Actual Barrier3D time index")
        axes[1, 1].set_xlabel("Actual Barrier3D time index")

        event_messages = []
        if current == NOURISHMENT_TIME_INDEX:
            event_messages.append("100 m3/m nourishment applied")
        for manager_name, series in series_by_manager.items():
            first_zero, first_migration = transition_indices(series)
            if current == first_zero:
                event_messages.append(f"{manager_name}: beach width = 0, migration ON")
            if series["shoreline_change_cells"][current] < 0:
                cells = int(abs(series["shoreline_change_cells"][current]))
                event_messages.append(f"{manager_name}: Barrier3D moved dune grid {cells} cell(s)")
        subtitle = " | ".join(event_messages) if event_messages else "Barrier3D calculating evolution"
        figure.suptitle(
            f"REAL BARRIER3D — no storms — time index {current}\n{subtitle}",
            fontsize=14,
            fontweight="bold",
        )
        figure.tight_layout(rect=(0, 0, 1, 0.92))
        return []

    animation = FuncAnimation(
        figure,
        update,
        frames=len(frame_indices),
        interval=1000,
        repeat=True,
        blit=False,
    )
    animation.save(GIF_PATH, writer=PillowWriter(fps=1), dpi=100)
    plt.close(figure)


def test_create_full_barrier3d_outputs(actual_runs, actual_series):
    """Write auditable values and plots from the arrays tested above."""

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for manager_name, cascade in actual_runs.items():
        np.savez(NPZ_PATHS[manager_name], cascade=cascade)
    write_csv(actual_series)
    plot_static(actual_series)
    plot_gif(actual_series)

    manifest = {
        "model": "real barrier3d.barrier3d.Barrier3d through Cascade.update",
        "source_tree": str(SOURCE_ROOT),
        "input_root": str(INPUT_ROOT),
        "input_files": list(INPUT_FILES),
        "saved_full_run_npz": {
            manager_name: str(path) for manager_name, path in NPZ_PATHS.items()
        },
        "controlled_parameter_change": {
            "StormStart": NO_STORM_START,
            "reason": "No storm occurs within the 45-state comparison",
        },
        "nourishment_request": {
            "time_index": NOURISHMENT_TIME_INDEX,
            "volume_m3_per_m": NOURISHMENT_VOLUME_M3_PER_M,
            "count": 1,
        },
        "cascade_parameters": {
            **COMMON_CASCADE_PARAMETERS,
        },
        "road_relocation_check": {
            "road_setback_m": COMMON_CASCADE_PARAMETERS["road_setback"],
            "source": "established CASCADE test_human_dynamics roadway input",
            "all_relocation_diagnostics_required_false": True,
        },
        "results": {
            manager_name: {
                "first_zero_beach_width_index": transition_indices(series)[0],
                "first_actual_dune_grid_migration_index": transition_indices(series)[1],
                "actual_migration_indices": np.flatnonzero(
                    series["shoreline_change_cells"] < 0
                ).tolist(),
            }
            for manager_name, series in actual_series.items()
        },
    }
    with MANIFEST_PATH.open("w") as stream:
        json.dump(manifest, stream, indent=2)

    for path in (
        FIGURE_PATH,
        GIF_PATH,
        CSV_PATH,
        MANIFEST_PATH,
        *NPZ_PATHS.values(),
    ):
        assert path.is_file() and path.stat().st_size > 0
