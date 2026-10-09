"""Full-model integration test: roadway nourishment then natural relocation.

The real Barrier3D model runs through Cascade.update on one established
CASCADE numerical test segment. StormStart is moved beyond the experiment,
one 100 m3/m nourish_now request is issued at time index 1, and no shoreline,
dune, road, or relocation state is prescribed after initialization.
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
from matplotlib.patches import Rectangle

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
TIME_STEP_COUNT = 80
NO_STORM_START = TIME_STEP_COUNT + 1
NOURISHMENT_TIME_INDEX = 1
NOURISHMENT_VOLUME_M3_PER_M = 100.0
INITIAL_BEACH_WIDTH_M = 30.0

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment_relocation"
    / "actual_barrier3d_no_storm_integration"
)
FIGURE_PATH = OUTPUT_DIR / "roadway_nourishment_relocation_lifecycle.png"
DOMAIN_FIGURE_PATH = OUTPUT_DIR / "roadway_nourishment_relocation_domain.png"
GIF_PATH = OUTPUT_DIR / "roadway_nourishment_relocation_slow.gif"
CSV_PATH = OUTPUT_DIR / "roadway_nourishment_relocation_lifecycle.csv"
MANIFEST_PATH = OUTPUT_DIR / "roadway_nourishment_relocation_manifest.json"
NPZ_PATH = OUTPUT_DIR / "roadway_nourishment_relocation.npz"

CASCADE_PARAMETERS = {
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
    "roadway_management_module": True,
    "alongshore_transport_module": False,
    "beach_nourishment_module": False,
    "community_economics_module": False,
    "nourishment_interval": None,
    "nourishment_volume": NOURISHMENT_VOLUME_M3_PER_M,
    "outwash_module": False,
    "road_ele": 1.2,
    "road_width": 20,
    "road_setback": 20,
    "dune_design_elevation": 3.2,
    "dune_minimum_elevation": 1.7,
}


def prepare_inputs(run_directory: Path) -> None:
    run_directory.mkdir(parents=True, exist_ok=True)
    for filename in INPUT_FILES:
        shutil.copy2(INPUT_ROOT / filename, run_directory / filename)

    parameter_path = run_directory / "nourishment-parameters.yaml"
    with parameter_path.open() as stream:
        parameters = yaml.safe_load(stream)
    parameters["StormStart"] = NO_STORM_START
    with parameter_path.open("w") as stream:
        yaml.safe_dump(parameters, stream, sort_keys=False)


def run_integration_scenario(run_directory: Path) -> Cascade:
    prepare_inputs(run_directory)
    cascade = Cascade(
        str(run_directory),
        name="roadway_nourishment_then_natural_relocation",
        **CASCADE_PARAMETERS,
    )
    cascade.nourish_now = [1]
    for _ in range(TIME_STEP_COUNT - 1):
        cascade.update()
        assert not cascade.b3d_break, "Barrier3D drowned during integration test"
    return cascade


@pytest.fixture(scope="module")
def integration_run(tmp_path_factory):
    return run_integration_scenario(
        tmp_path_factory.mktemp("roadway_nourishment_relocation")
    )


def test_actual_barrier3d_rebuilding_is_permitted_below_20_m_setback(tmp_path):
    prepare_inputs(tmp_path)
    parameters = {
        **CASCADE_PARAMETERS,
        "road_setback": 10,
        "road_setback_trigger": 20.0,
    }
    cascade = Cascade(
        str(tmp_path),
        name="road_dune_rebuild_management",
        **parameters,
    )

    cascade.update()

    roadway = cascade.roadways[0]
    barrier = cascade.barrier3d[0]
    output_index = barrier.time_index - 1
    assert roadway._road_setback_TS[output_index] < 20.0
    assert not roadway.road_dune_rebuild_disabled_TS[output_index]


def extract_series(cascade: Cascade) -> dict[str, np.ndarray]:
    barrier = cascade.barrier3d[0]
    roadway = cascade.roadways[0]
    count = len(barrier.x_s_TS)
    shoreline_change = np.asarray(barrier.ShorelineChangeTS[:count], dtype=float)
    initial_dune_grid_m = (
        np.floor(barrier.x_s_TS[0] + INITIAL_BEACH_WIDTH_M / 10.0) * 10.0
    )
    dune_grid_m = initial_dune_grid_m + np.cumsum(-shoreline_change) * 10.0
    road_setback_m = np.asarray(roadway._road_setback_TS[:count], dtype=float)
    road_width_m = np.asarray(roadway._road_width_TS[:count], dtype=float)
    road_m = dune_grid_m + road_setback_m
    road_m[road_width_m <= 0] = np.nan

    return {
        "time_index": np.arange(count),
        "shoreline_m": np.asarray(barrier.x_s_TS[:count], dtype=float) * 10.0,
        "beach_width_m": np.asarray(roadway.beach_width[:count], dtype=float),
        "dune_migration_on": np.asarray(
            roadway.dune_migration_on[:count], dtype=bool
        ),
        "shoreline_change_cells": shoreline_change,
        "dune_grid_m": dune_grid_m,
        "road_setback_m": road_setback_m,
        "road_width_m": road_width_m,
        "road_seaward_edge_m": road_m,
        "nourishment": np.asarray(roadway.nourishment_TS[:count], dtype=bool),
        "nourishment_volume_m3_per_m": np.asarray(
            roadway.nourishment_volume_TS[:count], dtype=float
        ),
        "triggered_relocation": np.asarray(
            roadway.triggered_relocation_TS[:count], dtype=bool
        ),
        "incomplete_relocation": np.asarray(
            roadway.relocation_incomplete_TS[:count], dtype=bool
        ),
        "historical_relocation_requested": np.asarray(
            roadway.historical_relocation_requested_TS[:count], dtype=bool
        ),
        "forced_relocation": np.asarray(
            roadway.forced_relocation_TS[:count], dtype=bool
        ),
        "storm_count": np.asarray(barrier._StormCount[:count], dtype=int),
        "overwash_flux_m3_per_m": np.asarray(barrier.QowTS[:count], dtype=float),
    }


@pytest.fixture(scope="module")
def series(integration_run):
    return extract_series(integration_run)


def lifecycle_indices(series: dict[str, np.ndarray]) -> dict[str, int]:
    return {
        "nourishment": int(np.flatnonzero(series["nourishment"])[0]),
        "zero_beach": int(
            np.flatnonzero(np.isclose(series["beach_width_m"], 0))[0]
        ),
        "first_migration": int(
            np.flatnonzero(series["shoreline_change_cells"] < 0)[0]
        ),
        "first_relocation": int(
            np.flatnonzero(series["triggered_relocation"])[0]
        ),
    }


def test_real_barrier3d_runs_with_zero_storms_and_overwash(integration_run, series):
    assert isinstance(integration_run.barrier3d[0], Barrier3d)
    np.testing.assert_array_equal(series["storm_count"], 0)
    np.testing.assert_allclose(series["overwash_flux_m3_per_m"], 0.0)


def test_one_nourishment_event_is_recorded_with_requested_volume(series):
    np.testing.assert_array_equal(
        np.flatnonzero(series["nourishment"]), [NOURISHMENT_TIME_INDEX]
    )
    assert series["nourishment_volume_m3_per_m"][
        NOURISHMENT_TIME_INDEX
    ] == pytest.approx(NOURISHMENT_VOLUME_M3_PER_M)


def test_nourishment_progrades_shoreline_and_widens_beach(series):
    event = NOURISHMENT_TIME_INDEX
    assert series["shoreline_m"][event] < series["shoreline_m"][event - 1]
    assert series["beach_width_m"][event] > series["beach_width_m"][event - 1]


def test_dune_grid_stays_put_until_beach_erodes_away(series):
    indices = lifecycle_indices(series)
    event = indices["nourishment"]
    zero_beach = indices["zero_beach"]
    np.testing.assert_allclose(
        series["dune_grid_m"][event : zero_beach + 1],
        series["dune_grid_m"][event],
    )
    assert not series["dune_migration_on"][event]
    assert series["dune_migration_on"][zero_beach]


def test_model_reactivates_and_calculates_dune_migration(series):
    indices = lifecycle_indices(series)
    assert indices["first_migration"] > indices["zero_beach"]
    assert series["shoreline_change_cells"][indices["first_migration"]] < 0


def test_natural_relocation_occurs_after_nourishment_lifecycle(series):
    indices = lifecycle_indices(series)
    relocations = np.flatnonzero(series["triggered_relocation"])
    assert relocations[0] > indices["first_migration"]
    for relocation in relocations:
        assert series["shoreline_change_cells"][relocation] < 0
        prior_setback = series["road_setback_m"][relocation - 1]
        migration_m = series["shoreline_change_cells"][relocation] * 10.0
        assert prior_setback + migration_m < 0


def test_relocation_moves_road_landward_and_resets_setback(series):
    for relocation in np.flatnonzero(series["triggered_relocation"]):
        assert series["road_setback_m"][relocation] == pytest.approx(
            CASCADE_PARAMETERS["road_setback"]
        )
        assert series["road_seaward_edge_m"][relocation] > series[
            "road_seaward_edge_m"
        ][relocation - 1]


def test_no_forced_historical_or_incomplete_relocation(series):
    assert not series["historical_relocation_requested"].any()
    assert not series["forced_relocation"].any()
    assert not series["incomplete_relocation"].any()


def write_csv(series: dict[str, np.ndarray]) -> None:
    with CSV_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(series))
        writer.writeheader()
        for index in range(len(series["time_index"])):
            writer.writerow(
                {field: values[index] for field, values in series.items()}
            )


def plot_lifecycle(series: dict[str, np.ndarray]) -> None:
    indices = lifecycle_indices(series)
    time = series["time_index"]
    figure, axes = plt.subplots(4, 1, figsize=(16, 14), sharex=True)

    axes[0].plot(time, series["shoreline_m"], color="#08519c", linewidth=2.3)
    axes[0].set_title("Barrier3D shoreline position (smaller = seaward)")
    axes[0].set_ylabel("Position (m)")
    axes[1].plot(time, series["beach_width_m"], color="#d95f0e", linewidth=2.3)
    axes[1].fill_between(
        time, 0, series["beach_width_m"], color="#fdd49e", alpha=0.6
    )
    axes[1].set_title("RoadwayManager beach width")
    axes[1].set_ylabel("Width (m)")
    axes[2].plot(
        time,
        series["dune_grid_m"],
        color="#8c510a",
        linewidth=2.3,
        label="Dune-grid edge",
    )
    axes[2].plot(
        time,
        series["road_seaward_edge_m"],
        color="#3f3f3f",
        linewidth=2.5,
        label="Road seaward edge",
    )
    axes[2].set_title("Calculated dune-grid and roadway positions")
    axes[2].set_ylabel("Absolute position (m)")
    axes[2].legend(loc="best")
    axes[3].step(
        time,
        series["road_setback_m"],
        where="post",
        color="#252525",
        linewidth=2.3,
        label="Road setback",
    )
    axes[3].step(
        time,
        series["dune_migration_on"].astype(int) * 5,
        where="post",
        color="#756bb1",
        linewidth=2,
        label="Dune migration ON (shown at 5 m)",
    )
    relocations = np.flatnonzero(series["triggered_relocation"])
    axes[3].scatter(
        relocations,
        series["road_setback_m"][relocations],
        marker="*",
        s=200,
        color="#e31a1c",
        edgecolor="black",
        zorder=5,
        label="Triggered relocation",
    )
    axes[3].set_title("Road setback, dune-migration state, and relocation")
    axes[3].set_ylabel("State/distance")
    axes[3].set_xlabel("Actual Barrier3D time index")
    axes[3].legend(loc="best", ncol=3)

    event_styles = (
        (indices["nourishment"], "#2171b5", "nourish_now"),
        (indices["zero_beach"], "#d95f0e", "beach = 0"),
        (indices["first_migration"], "#756bb1", "first migration"),
        (indices["first_relocation"], "#e31a1c", "road relocation"),
    )
    for axis in axes:
        axis.grid(alpha=0.25)
        for event, color, _ in event_styles:
            axis.axvline(event, color=color, linestyle="--", alpha=0.75)
    for event, color, label in event_styles:
        axes[0].text(
            event + 0.4,
            axes[0].get_ylim()[1],
            label,
            color=color,
            fontsize=9,
            va="top",
            rotation=90,
        )

    figure.suptitle(
        "RoadwayManager: nourishment → beach erosion → dune migration → road relocation\n"
        "Full Barrier3D model; one controlled segment; zero storms/overtopping",
        fontsize=16,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(FIGURE_PATH, dpi=190, bbox_inches="tight")
    plt.close(figure)


def combined_domain(cascade, series, index):
    barrier = cascade.barrier3d[0]
    dunes = (
        np.asarray(barrier.DuneDomain[index], dtype=float).T + barrier.BermEl
    ) * 10.0
    interior = np.asarray(barrier.DomainTS[index], dtype=float) * 10.0
    domain = np.vstack((dunes, interior))
    y_edges = (
        series["dune_grid_m"][index]
        + np.arange(domain.shape[0] + 1, dtype=float) * 10.0
    )
    x_edges = np.arange(domain.shape[1] + 1, dtype=float) * 10.0
    return x_edges, y_edges, domain


def draw_domain(axis, cascade, series, index):
    x_edges, y_edges, domain = combined_domain(cascade, series, index)
    mesh = axis.pcolormesh(
        x_edges,
        y_edges,
        domain,
        cmap="terrain",
        vmin=-1,
        vmax=3.5,
        shading="flat",
    )
    dune = series["dune_grid_m"][index]
    road = series["road_seaward_edge_m"][index]
    road_width = series["road_width_m"][index]
    axis.axhline(
        dune, color="#8c510a", linewidth=2.2, label="Dune-grid edge"
    )
    if np.isfinite(road) and road_width > 0:
        axis.add_patch(
            Rectangle(
                (0, road),
                x_edges[-1],
                road_width,
                facecolor="#3f3f3f",
                edgecolor="black",
                alpha=0.75,
                label="Road",
            )
        )
    if series["triggered_relocation"][index] and index > 0:
        old_road = series["road_seaward_edge_m"][index - 1]
        old_width = series["road_width_m"][index - 1]
        axis.add_patch(
            Rectangle(
                (0, old_road),
                x_edges[-1],
                old_width,
                facecolor="none",
                edgecolor="#e31a1c",
                linestyle="--",
                linewidth=2.2,
                label="Previous road",
            )
        )
    axis.set_xlim(0, x_edges[-1])
    axis.set_ylim(y_edges[0] - 5, y_edges[-1])
    axis.set_xlabel("Alongshore distance (m)")
    axis.set_ylabel("Absolute cross-shore position (m; landward ↑)")
    return mesh


def plot_domain_snapshots(cascade, series) -> None:
    indices = lifecycle_indices(series)
    snapshots = [
        0,
        indices["nourishment"],
        indices["zero_beach"],
        indices["first_migration"],
        indices["first_relocation"],
    ]
    labels = [
        "Initial",
        "Nourished",
        "Beach reaches zero",
        "First dune migration",
        "Road relocated",
    ]
    figure, axes = plt.subplots(1, 5, figsize=(22, 7))
    mesh = None
    for axis, index, label in zip(axes, snapshots, labels):
        mesh = draw_domain(axis, cascade, series, index)
        axis.set_title(f"{label}\ntime {index}")
    axes[-1].legend(loc="upper right", fontsize=8)
    colorbar = figure.colorbar(mesh, ax=axes, shrink=0.75, pad=0.02)
    colorbar.set_label("Saved elevation (m MHW)")
    figure.suptitle(
        "Saved Barrier3D domain through nourishment and natural road relocation",
        fontsize=16,
        fontweight="bold",
    )
    figure.subplots_adjust(
        top=0.86, bottom=0.13, left=0.05, right=0.91, wspace=0.36
    )
    figure.savefig(DOMAIN_FIGURE_PATH, dpi=180, bbox_inches="tight")
    plt.close(figure)


def plot_gif(cascade, series) -> None:
    indices = lifecycle_indices(series)
    key_events = set(indices.values())
    count = len(series["time_index"])
    frames = []
    for index in range(count):
        frames.extend([index] * (6 if index in key_events else 1))

    figure, (domain_axis, lifecycle_axis) = plt.subplots(
        1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": (1, 1.15)}
    )

    def update(frame_number):
        current = frames[frame_number]
        domain_axis.clear()
        lifecycle_axis.clear()
        draw_domain(domain_axis, cascade, series, current)
        domain_axis.set_title(f"Saved model elevation domain — time {current}")

        visible = slice(0, current + 1)
        time = series["time_index"][visible]
        lifecycle_axis.plot(
            time,
            series["beach_width_m"][visible],
            color="#d95f0e",
            linewidth=2.5,
            label="Beach width",
        )
        lifecycle_axis.step(
            time,
            series["road_setback_m"][visible],
            where="post",
            color="#252525",
            linewidth=2.5,
            label="Road setback",
        )
        lifecycle_axis.step(
            time,
            series["dune_migration_on"][visible].astype(int) * 5,
            where="post",
            color="#756bb1",
            linewidth=2,
            label="Dune migration ON (at 5)",
        )
        passed = np.flatnonzero(
            series["triggered_relocation"][: current + 1]
        )
        if passed.size:
            lifecycle_axis.scatter(
                passed,
                series["road_setback_m"][passed],
                marker="*",
                s=180,
                color="#e31a1c",
                edgecolor="black",
                zorder=5,
                label="Road relocation",
            )
        lifecycle_axis.set_xlim(0, count - 1)
        lifecycle_axis.set_ylim(
            -2, max(55, np.max(series["beach_width_m"]) + 5)
        )
        lifecycle_axis.set_xlabel("Actual Barrier3D time index")
        lifecycle_axis.set_ylabel("Width, setback, or plotted state (m)")
        lifecycle_axis.set_title("Combined functionality lifecycle")
        lifecycle_axis.grid(alpha=0.25)
        lifecycle_axis.legend(loc="upper right")

        if current == indices["nourishment"]:
            message = "NOURISH NOW: shoreline progrades; dune migration OFF"
        elif current == indices["zero_beach"]:
            message = "BEACH WIDTH = 0: dune migration ON"
        elif current == indices["first_migration"]:
            message = "BARRIER3D MOVES DUNE GRID LANDWARD"
        elif current == indices["first_relocation"]:
            message = "ROAD RELOCATED: setback reset to 20 m"
        else:
            message = "Barrier3D calculating no-storm evolution"
        figure.suptitle(
            f"ROADWAY NOURISHMENT + NATURAL RELOCATION — time {current}\n"
            f"{message}",
            fontsize=14,
            fontweight="bold",
        )
        figure.tight_layout(rect=(0, 0, 1, 0.91))
        return []

    animation = FuncAnimation(
        figure,
        update,
        frames=len(frames),
        interval=1000,
        repeat=True,
        blit=False,
    )
    animation.save(GIF_PATH, writer=PillowWriter(fps=1), dpi=100)
    plt.close(figure)


def test_create_auditable_combined_outputs(integration_run, series):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(NPZ_PATH, cascade=integration_run)
    write_csv(series)
    plot_lifecycle(series)
    plot_domain_snapshots(integration_run, series)
    plot_gif(integration_run, series)

    indices = lifecycle_indices(series)
    manifest = {
        "purpose": (
            "Modified RoadwayManager nourishment-plus-relocation integration test"
        ),
        "model": "barrier3d.barrier3d.Barrier3d through Cascade.update",
        "source_tree": str(SOURCE_ROOT),
        "input_root": str(INPUT_ROOT),
        "input_files": list(INPUT_FILES),
        "controlled_parameter_change": {
            "StormStart": NO_STORM_START,
            "reason": "Prevent storms and overtopping during integration test",
        },
        "state_manipulation_after_initialization": "none",
        "nourishment_request": {
            "time_index": NOURISHMENT_TIME_INDEX,
            "volume_m3_per_m": NOURISHMENT_VOLUME_M3_PER_M,
            "count": 1,
        },
        "historical_or_forced_relocation_requests": "none",
        "cascade_parameters": CASCADE_PARAMETERS,
        "results": {
            **indices,
            "triggered_relocation_indices": np.flatnonzero(
                series["triggered_relocation"]
            ).tolist(),
            "incomplete_relocation_indices": np.flatnonzero(
                series["incomplete_relocation"]
            ).tolist(),
        },
    }
    with MANIFEST_PATH.open("w") as stream:
        json.dump(manifest, stream, indent=2)

    for path in (
        FIGURE_PATH,
        DOMAIN_FIGURE_PATH,
        GIF_PATH,
        CSV_PATH,
        MANIFEST_PATH,
        NPZ_PATH,
    ):
        assert path.is_file() and path.stat().st_size > 0
