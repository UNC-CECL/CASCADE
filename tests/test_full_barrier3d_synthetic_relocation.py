"""Actual-Barrier3D demonstration of naturally triggered road relocation.

This is a one-segment synthetic/controlled CASCADE experiment built from the
repository's established ``test_human_dynamics`` inputs.  Shoreline change,
dune-grid migration, storms, overwash, and road relocation are all calculated
through the normal ``Cascade.update`` loop.  The test does not prescribe any
shoreline or dune movement and makes no nourishment or historical relocation
requests.
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
from matplotlib.animation import FuncAnimation, PillowWriter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

from barrier3d import Barrier3d  # noqa: E402
from cascade import Cascade  # noqa: E402


INPUT_ROOT = SOURCE_ROOT / "tests" / "test_human_dynamics"
INPUT_FILES = (
    "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "b3d_pt75_3284yrs_low-elevations.csv",
    "pathways-dunes.npy",
    "growthparam_1000dam.npy",
    "roadway-parameters.yaml",
)
TIME_STEP_COUNT = 180
INITIAL_BEACH_WIDTH_M = 30.0

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_relocation"
    / "actual_barrier3d_synthetic_test"
)
FIGURE_PATH = OUTPUT_DIR / "actual_barrier3d_road_relocation.png"
GIF_PATH = OUTPUT_DIR / "actual_barrier3d_road_relocation_slow.gif"
CSV_PATH = OUTPUT_DIR / "actual_barrier3d_road_relocation.csv"
MANIFEST_PATH = OUTPUT_DIR / "actual_barrier3d_road_relocation_manifest.json"
NPZ_PATH = OUTPUT_DIR / "actual_barrier3d_road_relocation.npz"

CASCADE_PARAMETERS = {
    "storm_file": "StormSeries_1kyrs_VCR_Berm1pt9m_Slope0pt04_01.npy",
    "elevation_file": "b3d_pt75_3284yrs_low-elevations.csv",
    "dune_file": "pathways-dunes.npy",
    "parameter_file": "roadway-parameters.yaml",
    "wave_height": 1,
    "wave_period": 7,
    "wave_asymmetry": 0.8,
    "wave_angle_high_fraction": 0.2,
    "sea_level_rise_rate": 0.004,
    "sea_level_rise_constant": True,
    "background_erosion": 0.0,
    "alongshore_section_count": 1,
    "time_step_count": TIME_STEP_COUNT,
    "min_dune_growth_rate": 0.55,
    "max_dune_growth_rate": 0.95,
    "num_cores": 1,
    "roadway_management_module": True,
    "alongshore_transport_module": False,
    "beach_nourishment_module": False,
    "community_economics_module": False,
    "road_ele": 1.2,
    "road_width": 20,
    "road_setback": 20,
    "dune_design_elevation": 3.2,
    "dune_minimum_elevation": 1.7,
    "outwash_module": False,
}


def prepare_inputs(run_directory: Path) -> None:
    """Copy established inputs without editing them."""

    run_directory.mkdir(parents=True, exist_ok=True)
    for filename in INPUT_FILES:
        shutil.copy2(INPUT_ROOT / filename, run_directory / filename)


def run_scenario(run_directory: Path) -> Cascade:
    """Run the actual one-segment model with no prescribed state changes."""

    prepare_inputs(run_directory)
    cascade = Cascade(
        str(run_directory),
        name="actual_barrier3d_synthetic_road_relocation",
        **CASCADE_PARAMETERS,
    )
    for _ in range(TIME_STEP_COUNT - 1):
        cascade.update()
        if cascade.b3d_break:
            break
    return cascade


@pytest.fixture(scope="module")
def actual_run(tmp_path_factory):
    return run_scenario(tmp_path_factory.mktemp("actual_barrier3d_relocation"))


def extract_series(cascade: Cascade) -> dict[str, np.ndarray]:
    barrier = cascade.barrier3d[0]
    roadway = cascade.roadways[0]
    count = len(barrier.x_s_TS)

    shoreline_change_cells = np.asarray(
        barrier.ShorelineChangeTS[:count], dtype=float
    )
    initial_dune_grid_m = (
        np.floor(barrier.x_s_TS[0] + INITIAL_BEACH_WIDTH_M / 10.0) * 10.0
    )
    dune_grid_m = initial_dune_grid_m + np.cumsum(-shoreline_change_cells) * 10.0

    road_setback_m = np.asarray(roadway._road_setback_TS[:count], dtype=float)
    road_width_m = np.asarray(roadway._road_width_TS[:count], dtype=float)
    managed = road_width_m > 0
    road_seaward_edge_m = np.full(count, np.nan)
    road_seaward_edge_m[managed] = dune_grid_m[managed] + road_setback_m[managed]

    return {
        "time_index": np.arange(count),
        "shoreline_m": np.asarray(barrier.x_s_TS[:count], dtype=float) * 10.0,
        "shoreline_change_cells": shoreline_change_cells,
        "dune_grid_m": dune_grid_m,
        "road_setback_m": road_setback_m,
        "road_width_m": road_width_m,
        "road_seaward_edge_m": road_seaward_edge_m,
        "road_elevation_m_mhw": np.asarray(
            roadway._road_ele_TS[:count], dtype=float
        ),
        "average_barrier_width_m": np.asarray(
            barrier.InteriorWidth_AvgTS[:count], dtype=float
        )
        * 10.0,
        "storm_count": np.asarray(barrier._StormCount[:count], dtype=int),
        "overwash_flux_m3_per_m": np.asarray(barrier.QowTS[:count], dtype=float),
        "dunes_rebuilt": np.asarray(
            roadway._dunes_rebuilt_TS[:count], dtype=bool
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
        "nourishment": np.asarray(roadway.nourishment_TS[:count], dtype=bool),
    }


@pytest.fixture(scope="module")
def series(actual_run):
    return extract_series(actual_run)


def test_uses_actual_barrier3d_and_normal_update_loop(actual_run):
    assert isinstance(actual_run.barrier3d[0], Barrier3d)


def test_relocation_is_natural_and_not_requested(series):
    events = np.flatnonzero(series["triggered_relocation"])
    assert events.size > 0
    assert not series["historical_relocation_requested"].any()
    assert not series["forced_relocation"].any()
    assert not series["incomplete_relocation"].any()
    assert not series["nourishment"].any()


def test_each_relocation_follows_model_calculated_dune_migration(series):
    for event in np.flatnonzero(series["triggered_relocation"]):
        assert event > 0
        assert series["shoreline_change_cells"][event] < 0
        prior_setback = series["road_setback_m"][event - 1]
        migration_m = series["shoreline_change_cells"][event] * 10.0
        assert prior_setback + migration_m < 0


def test_each_successful_relocation_resets_setback_and_moves_road_landward(series):
    target_setback = float(CASCADE_PARAMETERS["road_setback"])
    for event in np.flatnonzero(series["triggered_relocation"]):
        assert series["road_setback_m"][event] == pytest.approx(target_setback)
        assert series["road_seaward_edge_m"][event] > series["road_seaward_edge_m"][
            event - 1
        ]


def draw_cross_shore(axis, series: dict[str, np.ndarray], index: int, event: bool) -> None:
    shoreline = float(series["shoreline_m"][index])
    dune = float(series["dune_grid_m"][index])
    road = float(series["road_seaward_edge_m"][index])
    road_width = float(series["road_width_m"][index])
    barrier_width = float(series["average_barrier_width_m"][index])
    bay = dune + barrier_width

    x_min = min(series["shoreline_m"]) - 20
    x_max = max(
        np.nanmax(series["road_seaward_edge_m"] + series["road_width_m"]) + 35,
        bay + 15,
    )
    axis.axvspan(x_min, shoreline, color="#9ecae1", alpha=0.95, label="Ocean")
    axis.axvspan(shoreline, dune, color="#fdd49e", alpha=0.95, label="Beach")
    axis.axvspan(dune, bay, color="#c7e9c0", alpha=0.95, label="Barrier interior")
    axis.axvspan(bay, x_max, color="#a6bddb", alpha=0.75, label="Back-barrier water")
    axis.axvline(shoreline, color="#08519c", linewidth=2.5, label="Shoreline")
    axis.axvspan(dune - 1.2, dune + 1.2, color="#8c510a", label="Dune grid")
    if np.isfinite(road) and road_width > 0:
        axis.add_patch(
            Rectangle(
                (road, 0.23),
                road_width,
                0.34,
                facecolor="#4d4d4d",
                edgecolor="black",
                linewidth=1.2,
                label="Road",
            )
        )
    if event and index > 0:
        old_road = float(series["road_seaward_edge_m"][index - 1])
        old_width = float(series["road_width_m"][index - 1])
        axis.add_patch(
            Rectangle(
                (old_road, 0.23),
                old_width,
                0.34,
                facecolor="none",
                edgecolor="#cb181d",
                linestyle="--",
                linewidth=2,
                label="Previous road",
            )
        )
        axis.add_patch(
            FancyArrowPatch(
                (old_road + old_width / 2, 0.67),
                (road + road_width / 2, 0.67),
                arrowstyle="-|>",
                mutation_scale=18,
                linewidth=2,
                color="#cb181d",
            )
        )

    axis.set_xlim(x_min, x_max)
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.set_xlabel("Absolute cross-shore position (m; landward →)")
    axis.grid(axis="x", alpha=0.2)


def write_csv(series: dict[str, np.ndarray]) -> None:
    with CSV_PATH.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(series))
        writer.writeheader()
        for index in range(len(series["time_index"])):
            writer.writerow({field: values[index] for field, values in series.items()})


def plot_static(series: dict[str, np.ndarray]) -> None:
    events = np.flatnonzero(series["triggered_relocation"])
    event = int(events[0])
    before = event - 1

    figure = plt.figure(figsize=(16, 12))
    grid = figure.add_gridspec(3, 2, height_ratios=(1, 1, 1.35))
    before_axis = figure.add_subplot(grid[0, 0])
    after_axis = figure.add_subplot(grid[0, 1], sharex=before_axis)
    position_axis = figure.add_subplot(grid[1, :])
    setback_axis = figure.add_subplot(grid[2, :], sharex=position_axis)

    draw_cross_shore(before_axis, series, before, event=False)
    draw_cross_shore(after_axis, series, event, event=True)
    before_axis.set_title(
        f"Immediately before first relocation — time {before}\n"
        f"road setback = {series['road_setback_m'][before]:.0f} m"
    )
    after_axis.set_title(
        f"Successful natural relocation — time {event}\n"
        f"new road setback = {series['road_setback_m'][event]:.0f} m"
    )
    handles, labels = after_axis.get_legend_handles_labels()
    after_axis.legend(handles, labels, loc="upper right", fontsize=8, ncol=2)

    time = series["time_index"]
    position_axis.plot(time, series["dune_grid_m"], color="#8c510a", linewidth=2, label="Actual dune-grid line")
    position_axis.plot(time, series["road_seaward_edge_m"], color="#4d4d4d", linewidth=2.5, label="Road seaward edge")
    position_axis.set_ylabel("Absolute position (m)")
    position_axis.set_title("Model-calculated dune-grid position and resulting road position")
    position_axis.legend(loc="best")
    position_axis.grid(alpha=0.25)

    setback_axis.step(time, series["road_setback_m"], where="post", color="#252525", linewidth=2, label="Road setback")
    migration = series["shoreline_change_cells"] * 10.0
    setback_axis.bar(time, migration, color=np.where(migration < 0, "#d95f0e", "#31a354"), alpha=0.55, label="Barrier3D dune-grid change")
    setback_axis.scatter(events, series["road_setback_m"][events], marker="*", s=190, color="#e31a1c", edgecolor="black", zorder=5, label="triggered_relocation_TS")
    for relocation_index in events:
        for axis in (position_axis, setback_axis):
            axis.axvline(relocation_index, color="#e31a1c", linestyle="--", alpha=0.65)
    setback_axis.axhline(0, color="black", linewidth=0.8)
    setback_axis.set_ylabel("Distance/change (m)")
    setback_axis.set_xlabel("Actual Barrier3D time index")
    setback_axis.set_title("Road setback, actual dune-grid changes, and relocation events")
    setback_axis.legend(loc="best", ncol=3)
    setback_axis.grid(alpha=0.25)

    figure.suptitle(
        "Actual Barrier3D one-segment roadway-relocation test\n"
        "No prescribed shoreline/dune movement; no nourishment; no forced request",
        fontsize=16,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(FIGURE_PATH, dpi=190, bbox_inches="tight")
    plt.close(figure)


def plot_gif(series: dict[str, np.ndarray]) -> None:
    events = set(np.flatnonzero(series["triggered_relocation"]).tolist())
    count = len(series["time_index"])
    selected = set(range(0, count, 3))
    selected.add(count - 1)
    for event in events:
        selected.update(range(max(0, event - 3), min(count, event + 4)))
    frame_indices = []
    for index in sorted(selected):
        frame_indices.extend([index] * (6 if index in events else 1))

    figure, (plan_axis, timeline_axis, setback_axis) = plt.subplots(
        3, 1, figsize=(15, 11), gridspec_kw={"height_ratios": (1.1, 1.2, 1.0)}
    )

    def update(frame_number):
        current = frame_indices[frame_number]
        for axis in (plan_axis, timeline_axis, setback_axis):
            axis.clear()

        draw_cross_shore(plan_axis, series, current, current in events)
        visible = slice(0, current + 1)
        time = series["time_index"][visible]
        timeline_axis.plot(time, series["dune_grid_m"][visible], color="#8c510a", linewidth=2.5, label="Actual dune-grid line")
        timeline_axis.plot(time, series["road_seaward_edge_m"][visible], color="#4d4d4d", linewidth=2.5, label="Road seaward edge")
        timeline_axis.set_xlim(0, count - 1)
        timeline_axis.set_ylim(
            np.nanmin(series["dune_grid_m"]) - 10,
            np.nanmax(series["road_seaward_edge_m"] + series["road_width_m"]) + 10,
        )
        timeline_axis.set_ylabel("Absolute position (m)")
        timeline_axis.set_title("Calculated dune and road positions")
        timeline_axis.legend(loc="upper left")
        timeline_axis.grid(alpha=0.25)

        setback_axis.step(time, series["road_setback_m"][visible], where="post", color="#252525", linewidth=2.5, label="Road setback")
        passed_events = np.array(sorted(event for event in events if event <= current), dtype=int)
        if passed_events.size:
            setback_axis.scatter(passed_events, series["road_setback_m"][passed_events], marker="*", s=180, color="#e31a1c", edgecolor="black", zorder=5, label="Triggered relocation")
        setback_axis.set_xlim(0, count - 1)
        setback_axis.set_ylim(-2, np.nanmax(series["road_setback_m"]) + 5)
        setback_axis.set_ylabel("Setback (m)")
        setback_axis.set_xlabel("Actual Barrier3D time index")
        setback_axis.set_title("Road setback and relocation record")
        setback_axis.legend(loc="upper right")
        setback_axis.grid(alpha=0.25)

        messages = [
            f"storms this step: {series['storm_count'][current]}",
            f"overwash flux: {series['overwash_flux_m3_per_m'][current]:.2f} m³/m",
            f"dune-grid change: {series['shoreline_change_cells'][current] * 10:.0f} m",
            f"road setback: {series['road_setback_m'][current]:.0f} m",
        ]
        if current in events:
            messages.insert(0, "ROAD RELOCATED: dune crossed road; new road built landward")
        figure.suptitle(
            f"ACTUAL BARRIER3D ROAD RELOCATION — time index {current}\n"
            + " | ".join(messages),
            fontsize=13,
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


def test_create_auditable_relocation_outputs(actual_run, series):
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(NPZ_PATH, cascade=actual_run)
    write_csv(series)
    plot_static(series)
    plot_gif(series)

    events = np.flatnonzero(series["triggered_relocation"])
    manifest = {
        "purpose": "Actual-model one-segment natural road-relocation demonstration",
        "model": "barrier3d.barrier3d.Barrier3d through Cascade.update",
        "source_tree": str(SOURCE_ROOT),
        "input_root": str(INPUT_ROOT),
        "input_files": list(INPUT_FILES),
        "state_manipulation_after_initialization": "none",
        "nourishment_requests": "none",
        "historical_or_forced_relocation_requests": "none",
        "cascade_parameters": CASCADE_PARAMETERS,
        "results": {
            "state_count": int(len(series["time_index"])),
            "triggered_relocation_indices": events.tolist(),
            "incomplete_relocation_indices": np.flatnonzero(
                series["incomplete_relocation"]
            ).tolist(),
            "road_management_stopped": bool(actual_run.road_break[0]),
            "barrier3d_drowned": bool(actual_run.b3d_break),
        },
    }
    with MANIFEST_PATH.open("w") as stream:
        json.dump(manifest, stream, indent=2)

    for path in (FIGURE_PATH, GIF_PATH, CSV_PATH, MANIFEST_PATH, NPZ_PATH):
        assert path.is_file() and path.stat().st_size > 0
