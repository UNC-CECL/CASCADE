"""Controlled nourishment lifecycle with no dune-overtopping storms.

The test uses one synthetic barrier segment, deep-copied for the original
BeachDuneManager and modified RoadwayManager. Both public manager ``update``
interfaces are exercised. The synthetic forcing has zero overwash and no
storm-driven dune-line movement. After nourishment, identical prescribed
shoreline retreat erodes the beach until it disappears. Two additional
shoreline-retreat steps then represent Barrier3D dune migration after both
managers have released the dune line.

This is a manager-level integration test, not a real Pea Island simulation.
"""

from __future__ import annotations

import csv
import sys
from copy import deepcopy
from pathlib import Path

import matplotlib
import numpy as np
import pytest
from matplotlib.animation import FuncAnimation, PillowWriter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# Always test the cumulative local source tree, even when pytest is launched
# from the parent workspace where a different ``cascade`` may be installed.
SOURCE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SOURCE_ROOT))

from cascade.beach_dune_manager import BeachDuneManager  # noqa: E402
from cascade.roadway_manager import RoadwayManager  # noqa: E402


WORKSPACE_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment"
    / "controlled_nourishment_test"
)
FIGURE_PATH = OUTPUT_DIR / "zero_overtopping_nourishment_lifecycle.png"
GIF_PATH = OUTPUT_DIR / "zero_overtopping_nourishment_lifecycle_slow.gif"
CSV_PATH = OUTPUT_DIR / "zero_overtopping_nourishment_lifecycle.csv"

NOURISHMENT_VOLUME_M3_PER_M = 100.0
EROSION_INCREMENT_M = 10.0
POST_RELEASE_MIGRATION_STEPS_M = (10.0, 10.0)
ROAD_SETBACK_M = 30.0


class NoOvertoppingSyntheticBarrier:
    """Minimum stable Barrier3D-like state required by both manager updates."""

    def __init__(self, time_step_count=9):
        alongshore_cells = 5
        cross_shore_cells = 12
        self.time_index = 2
        self.x_s = 100.0  # dam
        self.x_t = 0.0  # dam
        self.x_s_TS = [100.0, 100.0]
        self.x_b_TS = [114.0, 114.0]
        self.s_sf_TS = [0.01, 0.01]
        self.h_b_TS = [0.2, 0.2]
        self.InteriorWidth_AvgTS = [12.0, 12.0]
        self.QowTS = [0.0, 0.0]

        self.DShoreface = 1.0
        self.BermEl = 0.0
        self.Dmax = 0.5
        self.SL = 0.0
        self.BayDepth = 0.3
        self.BarrierLength = float(alongshore_cells)
        self.RSLR = np.zeros(time_step_count)
        self.ShorelineChangeTS = np.zeros(time_step_count)
        self.SCRagg = np.zeros(time_step_count)

        self.InteriorDomain = np.full(
            (cross_shore_cells, alongshore_cells), 0.2
        )
        self.PreStorm_InteriorDomain = self.InteriorDomain.copy()
        self.DomainTS = np.empty(time_step_count, dtype=object)
        for index in range(time_step_count):
            self.DomainTS[index] = self.InteriorDomain.copy()

        # Heights can evolve in a real run; this test controls only dune-line
        # position. They remain unchanged here because no Barrier3D storm/growth
        # step is run and neither manager is allowed to rebuild them.
        self.DuneDomain = np.full(
            (time_step_count, alongshore_cells, 2), 0.4
        )
        self.growthparam = np.full((1, alongshore_cells), 0.5)
        self.dune_migration_on = True
        self.migrate_dunes_call_count = 0

    def advance_shoreline_without_overtopping(self, retreat_m):
        """Apply retreat; move the dune line only if migration was already on."""

        migration_was_on = self.dune_migration_on
        self.time_index += 1
        self.x_s += retreat_m / 10.0
        self.x_s_TS.append(self.x_s)
        self.s_sf_TS.append(self.DShoreface / (self.x_s - self.x_t))
        self.h_b_TS.append(self.h_b_TS[-1])
        self.InteriorWidth_AvgTS.append(self.InteriorWidth_AvgTS[-1])
        self.x_b_TS.append(self.x_b_TS[-1])
        self.QowTS.append(0.0)
        self.PreStorm_InteriorDomain = self.InteriorDomain.copy()
        output_index = self.time_index - 1
        self.DomainTS[output_index] = self.InteriorDomain.copy()
        self.DuneDomain[output_index] = self.DuneDomain[output_index - 1].copy()
        # Before release, the manager absorbs shoreline retreat by reducing beach
        # width. After release, this prescribed value represents the Barrier3D
        # dune line migrating landward with the shoreline (negative is erosion).
        self.ShorelineChangeTS[output_index] = (
            -retreat_m / 10.0 if migration_was_on else 0.0
        )

    def FindWidths(self, interior_domain, sea_level):
        del sea_level
        widths = np.full(interior_domain.shape[1], interior_domain.shape[0])
        return widths, widths.copy(), float(np.mean(widths))

    def migrate_dunes(self, **kwargs):
        del kwargs
        self.migrate_dunes_call_count += 1
        raise AssertionError(
            "The dune line moved before the controlled beach-width test ended"
        )


def new_managers(time_step_count):
    beach_manager = BeachDuneManager(
        nourishment_interval=None,
        nourishment_volume=NOURISHMENT_VOLUME_M3_PER_M,
        initial_beach_width=30.0,
        time_step_count=time_step_count,
    )
    beach_manager.overwash_removal = False

    roadway_manager = RoadwayManager(
        initial_road_elevation=2.0,
        road_width=10.0,
        road_setback=ROAD_SETBACK_M,
        initial_dune_design_elevation=3.0,
        initial_dune_minimum_elevation=1.0,
        nourishment_interval=None,
        nourishment_volume=NOURISHMENT_VOLUME_M3_PER_M,
        initial_beach_width=30.0,
        time_step_count=time_step_count,
    )
    return beach_manager, roadway_manager


def record_state(stage, barrier, beach_width):
    shoreline_m = barrier.x_s * 10.0
    return {
        "stage": stage,
        "shoreline_position_m": shoreline_m,
        "beach_width_m": float(beach_width),
        # With coordinates increasing landward, shoreline + beach width is the
        # cross-shore dune-line position while migration is held off.
        # Round the displayed derived coordinate so the intentional 1e-6 m
        # threshold overshoot does not create a misleading plot scale.
        "dune_line_position_m": round(shoreline_m + float(beach_width), 4),
        "dune_migration_on": bool(barrier.dune_migration_on),
        "overwash_flux_m3_per_m": float(barrier.QowTS[-1]),
    }


def create_lifecycle_figure(beach_records, road_records):
    stages = [record["stage"] for record in beach_records]
    x = np.arange(len(stages))
    release_index = next(
        index
        for index, record in enumerate(beach_records[2:], start=2)
        if record["beach_width_m"] == 0
    )
    figure, axes = plt.subplots(4, 1, figsize=(15, 12), sharex=True)

    variables = [
        ("shoreline_position_m", "Shoreline position", "Position (m)"),
        ("beach_width_m", "Beach width", "Width (m)"),
        (
            "dune_line_position_m",
            "Dune-line position (height is not tested)",
            "Position (m)",
        ),
    ]
    for axis, (key, title, ylabel) in zip(axes[:3], variables):
        axis.plot(
            x,
            [record[key] for record in beach_records],
            "o-",
            color="#238b45",
            linewidth=2,
            label="Original BeachDuneManager",
        )
        axis.plot(
            x,
            [record[key] for record in road_records],
            "x--",
            color="#d95f0e",
            linewidth=1.5,
            markersize=8,
            label="Modified RoadwayManager",
        )
        axis.set_title(title)
        axis.set_ylabel(ylabel)
        axis.grid(alpha=0.25)
    axes[0].annotate(
        "Same nourishment applied\nshoreline progrades seaward",
        xy=(1, beach_records[1]["shoreline_position_m"]),
        xytext=(1.6, beach_records[0]["shoreline_position_m"] - 3),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
    )
    axes[1].axhline(0, color="black", linewidth=1)
    axes[1].annotate(
        "Beach is gone",
        xy=(release_index, 0),
        xytext=(release_index - 1.7, 12),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
    )
    axes[2].axhline(
        beach_records[0]["dune_line_position_m"],
        color="#666666",
        linewidth=1,
        alpha=0.6,
    )
    axes[2].annotate(
        "After release, dune line\nmigrates with shoreline",
        xy=(len(stages) - 1, beach_records[-1]["dune_line_position_m"]),
        xytext=(release_index - 1.1, beach_records[-1]["dune_line_position_m"] + 5),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
    )

    axes[3].step(
        x,
        [int(record["dune_migration_on"]) for record in beach_records],
        where="post",
        marker="o",
        color="#238b45",
        linewidth=2,
        label="Original BeachDuneManager",
    )
    axes[3].step(
        x,
        [int(record["dune_migration_on"]) for record in road_records],
        where="post",
        marker="x",
        color="#d95f0e",
        linewidth=1.5,
        linestyle="--",
        label="Modified RoadwayManager",
    )
    axes[3].set_title("Dune-migration permission")
    axes[3].set_ylabel("State")
    axes[3].set_yticks([0, 1], ["OFF", "ON"])
    axes[3].set_ylim(-0.2, 1.2)
    axes[3].grid(alpha=0.25)
    axes[3].annotate(
        "Migration allowed again",
        xy=(release_index, 1),
        xytext=(release_index - 2.2, 0.55),
        arrowprops={"arrowstyle": "->", "color": "#333333"},
    )

    for axis in axes:
        axis.axvline(1, color="#3182bd", linestyle=":", linewidth=1.5)
        axis.axvline(
            release_index,
            color="#756bb1",
            linestyle=":",
            linewidth=1.5,
        )
    axes[0].legend(loc="best")
    axes[3].set_xticks(x, stages, rotation=30, ha="right")
    axes[3].set_xlabel("Controlled test stage")
    figure.suptitle(
        "No-overtopping synthetic nourishment lifecycle\n"
        "Qow = 0 at every stage; manager curves overlap",
        fontsize=16,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    return figure


def create_lifecycle_gif(beach_records, road_records):
    """Animate the same controlled states shown in the static test figure."""

    stages = [record["stage"] for record in beach_records]
    stage_x = np.arange(len(stages))
    release_index = next(
        index
        for index, record in enumerate(beach_records[2:], start=2)
        if record["beach_width_m"] == 0
    )
    # At 0.5 frames/s, hold nourishment and release for six seconds each, then
    # hold the last post-release migration state for four seconds.
    frame_indices = [
        0,
        1,
        1,
        1,
        *range(2, release_index),
        release_index,
        release_index,
        release_index,
        *range(release_index + 1, len(stages) - 1),
        len(stages) - 1,
        len(stages) - 1,
    ]
    figure, axes = plt.subplots(2, 2, figsize=(15, 9))
    beach_axis, road_axis, position_axis, width_axis = axes.flat
    dune_line_m = beach_records[0]["dune_line_position_m"]
    initial_shoreline_m = beach_records[0]["shoreline_position_m"]
    nourished_shoreline_m = beach_records[1]["shoreline_position_m"]
    nourishment_progradation_m = initial_shoreline_m - nourished_shoreline_m
    road_center_m = dune_line_m + ROAD_SETBACK_M
    x_min = 980.0
    x_max = 1075.0

    def draw_plan_view(axis, record, title, show_road, frame_index):
        shoreline_m = record["shoreline_position_m"]
        current_dune_line_m = record["dune_line_position_m"]
        axis.axvspan(x_min, shoreline_m, color="#9ecae1", alpha=0.9, label="Ocean")
        if shoreline_m < current_dune_line_m:
            axis.axvspan(
                shoreline_m,
                current_dune_line_m,
                color="#fdd49e",
                alpha=0.95,
                label="Beach",
            )
        axis.axvspan(
            current_dune_line_m,
            x_max,
            color="#c7e9c0",
            alpha=0.9,
            label="Barrier interior",
        )
        axis.axvline(
            shoreline_m,
            color="#08519c",
            linewidth=3,
            label="Shoreline",
        )
        axis.axvspan(
            current_dune_line_m - 0.8,
            current_dune_line_m + 0.8,
            color="#8c510a",
            alpha=0.95,
            label="Dune line",
        )
        if frame_index == 1:
            # This is a visual guide to the shoreline area gained at the instant
            # of nourishment; it is not intended to track individual sand grains.
            axis.axvspan(
                nourished_shoreline_m,
                initial_shoreline_m,
                facecolor="#ffd92f",
                alpha=0.9,
                hatch="///",
                edgecolor="#e6550d",
                linewidth=2,
                label="New beach from nourishment",
                zorder=3,
            )
            axis.axvline(
                initial_shoreline_m,
                color="#08519c",
                linestyle=":",
                linewidth=2,
                zorder=4,
            )
            axis.annotate(
                f"NOURISHMENT\n100 m³/m\n+{nourishment_progradation_m:.2f} m beach",
                xy=(nourished_shoreline_m + nourishment_progradation_m / 2, 25),
                xytext=(nourished_shoreline_m + nourishment_progradation_m / 2, 36),
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
                color="#7f2704",
                arrowprops={"arrowstyle": "->", "color": "#e6550d", "linewidth": 2},
                bbox={
                    "boxstyle": "round,pad=0.35",
                    "facecolor": "#fff7bc",
                    "edgecolor": "#e6550d",
                    "alpha": 0.96,
                },
                zorder=6,
            )
        if show_road:
            axis.axvspan(
                road_center_m - 1.5,
                road_center_m + 1.5,
                color="#636363",
                alpha=0.95,
                label="Road",
            )
        axis.set_xlim(x_min, x_max)
        axis.set_ylim(0, 50)
        axis.set_yticks([])
        axis.set_xlabel("Cross-shore position (m; landward →)")
        axis.set_title(title, fontweight="bold")
        migration_text = "ON" if record["dune_migration_on"] else "OFF"
        migration_color = "#238b45" if record["dune_migration_on"] else "#cb181d"
        axis.text(
            0.02,
            0.96,
            f"Beach width: {record['beach_width_m']:.2f} m\n"
            f"Dune line: {current_dune_line_m:.2f} m\n"
            f"Migration: {migration_text}\nQow: 0 m³/m",
            transform=axis.transAxes,
            va="top",
            fontsize=10,
            bbox={
                "boxstyle": "round,pad=0.4",
                "facecolor": "white",
                "edgecolor": migration_color,
                "linewidth": 2,
                "alpha": 0.94,
            },
        )
        axis.text(
            shoreline_m,
            2,
            "shoreline",
            rotation=90,
            ha="right",
            va="bottom",
            color="#08519c",
            fontsize=9,
        )
        axis.text(
            current_dune_line_m,
            48,
            "fixed dune line" if frame_index <= release_index else "migrating dune line",
            rotation=90,
            ha="right",
            va="top",
            color="#8c510a",
            fontsize=9,
        )

    def update(frame_number):
        frame_index = frame_indices[frame_number]
        beach_record = beach_records[frame_index]
        road_record = road_records[frame_index]
        for axis in axes.flat:
            axis.clear()

        draw_plan_view(
            beach_axis,
            beach_record,
            "Original BeachDuneManager",
            show_road=False,
            frame_index=frame_index,
        )
        draw_plan_view(
            road_axis,
            road_record,
            "Modified RoadwayManager",
            show_road=True,
            frame_index=frame_index,
        )
        road_axis.legend(loc="lower right", fontsize=8, ncol=2)

        visible = slice(0, frame_index + 1)
        position_axis.plot(
            stage_x[visible],
            [record["shoreline_position_m"] for record in beach_records][visible],
            "o-",
            color="#238b45",
            linewidth=2,
            label="BeachDune shoreline",
        )
        position_axis.plot(
            stage_x[visible],
            [record["shoreline_position_m"] for record in road_records][visible],
            "x--",
            color="#d95f0e",
            linewidth=1.5,
            markersize=8,
            label="Roadway shoreline",
        )
        position_axis.plot(
            stage_x[visible],
            [record["dune_line_position_m"] for record in beach_records][visible],
            "s-",
            color="#8c510a",
            linewidth=2,
            label="Dune-line position (both)",
        )
        if frame_index == 1:
            position_axis.annotate(
                f"Nourish now: 100 m³/m\nshoreline progrades {nourishment_progradation_m:.2f} m",
                xy=(1, nourished_shoreline_m),
                xytext=(2.0, 990),
                fontsize=9,
                fontweight="bold",
                arrowprops={"arrowstyle": "->", "color": "#e6550d", "linewidth": 2},
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "#fff7bc"},
            )
        position_axis.set_xlim(-0.2, len(stages) - 0.8)
        position_axis.set_ylim(980, 1055)
        position_axis.set_xticks(stage_x, stages, rotation=28, ha="right")
        position_axis.set_ylabel("Cross-shore position (m)")
        position_axis.set_title("Position history through current stage")
        position_axis.grid(alpha=0.25)
        position_axis.legend(fontsize=8)

        width_axis.plot(
            stage_x[visible],
            [record["beach_width_m"] for record in beach_records][visible],
            "o-",
            color="#238b45",
            linewidth=2,
            label="BeachDune beach width",
        )
        width_axis.plot(
            stage_x[visible],
            [record["beach_width_m"] for record in road_records][visible],
            "x--",
            color="#d95f0e",
            linewidth=1.5,
            markersize=8,
            label="Roadway beach width",
        )
        width_axis.axhline(0, color="black", linewidth=1)
        width_axis.set_xlim(-0.2, len(stages) - 0.8)
        width_axis.set_ylim(-2, 50)
        width_axis.set_xticks(stage_x, stages, rotation=28, ha="right")
        width_axis.set_ylabel("Beach width (m)")
        width_axis.set_title("Beach-width history through current stage")
        width_axis.grid(alpha=0.25)
        width_axis.legend(fontsize=8)
        migration_text = "ON" if beach_record["dune_migration_on"] else "OFF"
        width_axis.text(
            0.98,
            0.92,
            f"Dune migration: {migration_text}",
            transform=width_axis.transAxes,
            ha="right",
            va="top",
            fontsize=11,
            fontweight="bold",
            color="#238b45" if migration_text == "ON" else "#cb181d",
        )
        if frame_index == release_index:
            width_axis.text(
                0.5,
                0.18,
                "Beach width reached zero:\ndune migration is allowed on the next physical step",
                transform=width_axis.transAxes,
                ha="center",
                fontsize=10,
                bbox={"boxstyle": "round,pad=0.4", "facecolor": "#d9f0d3"},
            )
        elif frame_index > release_index:
            width_axis.text(
                0.5,
                0.18,
                "Migration ON: shoreline and dune line move landward together\n"
                "while beach width remains zero",
                transform=width_axis.transAxes,
                ha="center",
                fontsize=10,
                bbox={"boxstyle": "round,pad=0.4", "facecolor": "#d9f0d3"},
            )

        if frame_index == 1:
            title = (
                "NOURISHMENT APPLIED IN BOTH MODULES — 100 m³/m\n"
                f"Shoreline progradation: {nourishment_progradation_m:.2f} m  |  "
                "dune line remains fixed"
            )
            title_color = "#d94801"
        elif frame_index > release_index:
            title = (
                "DUNE MIGRATION ACTIVE IN BOTH MODULES\n"
                f"Stage: {stages[frame_index]}  |  Qow = 0  |  "
                "shoreline and dune line retreat together"
            )
            title_color = "#238b45"
        else:
            title = (
                "No-overtopping synthetic barrier evolution — both modules\n"
                f"Stage: {stages[frame_index]}  |  Qow = 0  |  "
                "dune height is not evaluated"
            )
            title_color = "black"
        figure.suptitle(title, fontsize=15, fontweight="bold", color=title_color)
        figure.tight_layout(rect=(0, 0, 1, 0.93))
        return []

    animation = FuncAnimation(
        figure,
        update,
        frames=len(frame_indices),
        interval=2000,
        repeat=True,
        blit=False,
    )
    animation.save(GIF_PATH, writer=PillowWriter(fps=0.5), dpi=110)
    plt.close(figure)


def test_no_overtopping_nourishment_then_beach_erosion_releases_dunes():
    """Both managers hold dunes through erosion, then allow equal migration."""

    time_step_count = 9
    original = NoOvertoppingSyntheticBarrier(time_step_count)
    beach_barrier = deepcopy(original)
    road_barrier = deepcopy(original)
    beach_manager, roadway_manager = new_managers(time_step_count)
    original_dunes = original.DuneDomain.copy()

    beach_records = [record_state("Initial", beach_barrier, 30.0)]
    road_records = [record_state("Initial", road_barrier, 30.0)]

    beach_manager.update(
        barrier3d=beach_barrier,
        nourish_now=True,
        rebuild_dune_now=False,
        nourishment_interval=None,
    )
    roadway_manager.update(
        barrier3d=road_barrier,
        trigger_dune_knockdown=False,
        nourish_now=True,
    )
    event_index = beach_barrier.time_index - 1
    beach_records.append(
        record_state("Nourish", beach_barrier, beach_manager.beach_width[event_index])
    )
    road_records.append(
        record_state("Nourish", road_barrier, roadway_manager.beach_width[event_index])
    )

    assert beach_barrier.x_s < original.x_s
    assert road_barrier.x_s < original.x_s
    assert beach_barrier.x_s == pytest.approx(road_barrier.x_s)
    assert not beach_barrier.dune_migration_on
    assert not road_barrier.dune_migration_on
    np.testing.assert_array_equal(beach_barrier.DuneDomain, original_dunes)
    np.testing.assert_array_equal(road_barrier.DuneDomain, original_dunes)

    nourished_width = beach_manager.beach_width[event_index]
    assert roadway_manager.beach_width[event_index] == pytest.approx(nourished_width)
    erosion_steps = [EROSION_INCREMENT_M] * int(nourished_width // EROSION_INCREMENT_M)
    remainder = nourished_width - sum(erosion_steps)
    if remainder > 0:
        # A microscopic overshoot avoids a floating-point residual above zero.
        erosion_steps.append(remainder + 1e-6)
    assert (
        len(erosion_steps) + len(POST_RELEASE_MIGRATION_STEPS_M)
        == time_step_count - 2
    )

    for step, retreat_m in enumerate(erosion_steps, start=1):
        beach_barrier.advance_shoreline_without_overtopping(retreat_m)
        road_barrier.advance_shoreline_without_overtopping(retreat_m)
        beach_manager.update(
            barrier3d=beach_barrier,
            nourish_now=False,
            rebuild_dune_now=False,
            nourishment_interval=None,
        )
        roadway_manager.update(
            barrier3d=road_barrier,
            trigger_dune_knockdown=False,
            nourish_now=False,
        )
        index = beach_barrier.time_index - 1
        stage = f"Erode {step}"
        beach_records.append(
            record_state(stage, beach_barrier, beach_manager.beach_width[index])
        )
        road_records.append(
            record_state(stage, road_barrier, roadway_manager.beach_width[index])
        )

        assert beach_manager.beach_width[index] == pytest.approx(
            roadway_manager.beach_width[index]
        )
        assert beach_barrier.x_s == pytest.approx(road_barrier.x_s)
        assert beach_barrier.dune_migration_on == road_barrier.dune_migration_on
        assert beach_barrier.QowTS[-1] == 0.0
        assert road_barrier.QowTS[-1] == 0.0
        np.testing.assert_array_equal(
            beach_barrier.DuneDomain[index], original_dunes[index]
        )
        np.testing.assert_array_equal(
            road_barrier.DuneDomain[index], original_dunes[index]
        )

        if step < len(erosion_steps):
            assert beach_manager.beach_width[index] > 0
            assert not beach_barrier.dune_migration_on
        else:
            assert beach_manager.beach_width[index] == 0
            assert beach_barrier.dune_migration_on

    release_record_index = len(beach_records) - 1
    original_dune_line_m = beach_records[0]["dune_line_position_m"]
    np.testing.assert_allclose(
        [
            record["dune_line_position_m"]
            for record in beach_records[: release_record_index + 1]
        ],
        original_dune_line_m,
        atol=2e-6,
    )
    np.testing.assert_allclose(
        [
            record["dune_line_position_m"]
            for record in road_records[: release_record_index + 1]
        ],
        original_dune_line_m,
        atol=2e-6,
    )

    # These extra physical steps begin only after the managers have set migration
    # ON. They prescribe equal Barrier3D shoreline/dune retreat to both copies.
    previous_dune_line_m = original_dune_line_m
    for step, retreat_m in enumerate(POST_RELEASE_MIGRATION_STEPS_M, start=1):
        beach_barrier.advance_shoreline_without_overtopping(retreat_m)
        road_barrier.advance_shoreline_without_overtopping(retreat_m)
        beach_manager.update(
            barrier3d=beach_barrier,
            nourish_now=False,
            rebuild_dune_now=False,
            nourishment_interval=None,
        )
        roadway_manager.update(
            barrier3d=road_barrier,
            trigger_dune_knockdown=False,
            nourish_now=False,
        )
        index = beach_barrier.time_index - 1
        stage = f"Migrate {step}"
        beach_records.append(
            record_state(stage, beach_barrier, beach_manager.beach_width[index])
        )
        road_records.append(
            record_state(stage, road_barrier, roadway_manager.beach_width[index])
        )

        assert beach_manager.beach_width[index] == 0
        assert roadway_manager.beach_width[index] == 0
        assert beach_barrier.dune_migration_on
        assert road_barrier.dune_migration_on
        assert beach_barrier.x_s == pytest.approx(road_barrier.x_s)
        assert beach_barrier.ShorelineChangeTS[index] == pytest.approx(-retreat_m / 10)
        assert road_barrier.ShorelineChangeTS[index] == pytest.approx(-retreat_m / 10)
        assert beach_barrier.QowTS[-1] == 0.0
        assert road_barrier.QowTS[-1] == 0.0
        np.testing.assert_array_equal(
            beach_barrier.DuneDomain[index], original_dunes[index]
        )
        np.testing.assert_array_equal(
            road_barrier.DuneDomain[index], original_dunes[index]
        )

        expected_dune_line_m = previous_dune_line_m + retreat_m
        assert beach_records[-1]["dune_line_position_m"] == pytest.approx(
            expected_dune_line_m
        )
        assert road_records[-1]["dune_line_position_m"] == pytest.approx(
            expected_dune_line_m
        )
        previous_dune_line_m = expected_dune_line_m

    np.testing.assert_allclose(
        [record["dune_line_position_m"] for record in beach_records],
        [record["dune_line_position_m"] for record in road_records],
    )
    assert beach_barrier.migrate_dunes_call_count == 0
    assert road_barrier.migrate_dunes_call_count == 0
    assert not roadway_manager.triggered_relocation_TS.any()
    assert not roadway_manager.relocation_incomplete_TS.any()
    assert not roadway_manager.historical_relocation_requested_TS.any()
    assert not roadway_manager.forced_relocation_TS.any()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="") as stream:
        fieldnames = ["manager", *beach_records[0].keys()]
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for manager_name, records in [
            ("Original BeachDuneManager", beach_records),
            ("Modified RoadwayManager", road_records),
        ]:
            for record in records:
                writer.writerow({"manager": manager_name, **record})

    figure = create_lifecycle_figure(beach_records, road_records)
    figure.savefig(FIGURE_PATH, dpi=190, bbox_inches="tight")
    plt.close(figure)
    create_lifecycle_gif(beach_records, road_records)
    assert FIGURE_PATH.is_file() and FIGURE_PATH.stat().st_size > 0
    assert GIF_PATH.is_file() and GIF_PATH.stat().st_size > 0
    assert CSV_PATH.is_file() and CSV_PATH.stat().st_size > 0
