#!/usr/bin/env python3
"""Plot saved full-Barrier3D runs in the controlled-nourishment visual style.

This script is read-only with respect to the simulations. It loads the two
saved CASCADE objects and visualizes their actual Barrier3D shoreline and
``ShorelineChangeTS`` dune-grid movements. It does not calculate, prescribe, or
modify any physical state.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment"
    / "full_barrier3d_nourishment_test"
)
NPZ_PATHS = {
    "Original BeachDuneManager": OUTPUT_DIR
    / "actual_barrier3d_beach_dune_manager.npz",
    "Modified RoadwayManager": OUTPUT_DIR
    / "actual_barrier3d_roadway_manager.npz",
}
PNG_PATH = OUTPUT_DIR / "actual_barrier3d_controlled_nourishment_style.png"
GIF_PATH = OUTPUT_DIR / "actual_barrier3d_controlled_nourishment_style_slow.gif"

COLORS = {
    "Original BeachDuneManager": "#238b45",
    "Modified RoadwayManager": "#d95f0e",
}


def load_cascade(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Full Barrier3D output is missing: {path}")
    with np.load(path, allow_pickle=True) as archive:
        return archive["cascade"].item()


def extract(cascade, label: str) -> dict[str, np.ndarray | float]:
    barrier = cascade.barrier3d[0]
    manager = (
        cascade.nourishments[0]
        if label == "Original BeachDuneManager"
        else cascade.roadways[0]
    )
    shoreline_m = np.asarray(barrier.x_s_TS, dtype=float) * 10
    count = shoreline_m.size
    beach_width_m = np.asarray(manager.beach_width[:count], dtype=float)
    shoreline_change_cells = np.asarray(
        barrier.ShorelineChangeTS[:count], dtype=float
    )
    initial_dune_grid_m = shoreline_m[0] + beach_width_m[0]
    dune_grid_m = initial_dune_grid_m + np.cumsum(-shoreline_change_cells) * 10
    nourishment = np.asarray(manager._nourishment_TS[:count], dtype=bool)
    migration_on = np.asarray(manager._dune_migration_on[:count], dtype=float)
    return {
        "shoreline_m": shoreline_m,
        "beach_width_m": beach_width_m,
        "shoreline_change_cells": shoreline_change_cells,
        "dune_grid_m": dune_grid_m,
        "nourishment": nourishment,
        "migration_on": migration_on,
        "storm_count": np.asarray(barrier._StormCount[:count], dtype=int),
        "qow": np.asarray(barrier.QowTS[:count], dtype=float),
        "road_position_m": (
            np.nan
            if label == "Original BeachDuneManager"
            else initial_dune_grid_m + float(manager._road_setback_TS[0])
        ),
    }


def key_event_text(label: str, values: dict, index: int) -> list[str]:
    messages = []
    if values["nourishment"][index]:
        messages.append("nourish_now = 100 m³/m")
    if index > 0 and values["beach_width_m"][index] == 0:
        if values["beach_width_m"][index - 1] > 0:
            messages.append("beach width reached 0; migration ON")
    if values["shoreline_change_cells"][index] != 0:
        cells = int(abs(values["shoreline_change_cells"][index]))
        direction = "landward" if values["shoreline_change_cells"][index] < 0 else "seaward"
        messages.append(f"Barrier3D moved dune grid {cells} cell {direction}")
    return [f"{label}: {message}" for message in messages]


def draw_plan_view(axis, label: str, values: dict, index: int, x_limits) -> None:
    shoreline = float(values["shoreline_m"][index])
    dune_grid = float(values["dune_grid_m"][index])
    x_min, x_max = x_limits

    axis.axvspan(dune_grid, x_max, color="#c7e9c0", alpha=0.9, label="Barrier interior")
    axis.axvspan(x_min, shoreline, color="#9ecae1", alpha=0.9, label="Ocean")
    if shoreline < dune_grid:
        axis.axvspan(shoreline, dune_grid, color="#fdd49e", alpha=0.95, label="Beach")
    axis.axvline(shoreline, color="#08519c", linewidth=3, label="Shoreline")
    axis.axvspan(
        dune_grid - 0.8,
        dune_grid + 0.8,
        color="#8c510a",
        alpha=0.95,
        label="Actual dune-grid line",
    )
    road_position = float(values["road_position_m"])
    if np.isfinite(road_position):
        axis.axvspan(
            road_position - 1.5,
            road_position + 1.5,
            color="#636363",
            alpha=0.95,
            label="Road",
        )
    if values["nourishment"][index] and index > 0:
        prior_shoreline = float(values["shoreline_m"][index - 1])
        left, right = sorted((shoreline, prior_shoreline))
        axis.axvspan(
            left,
            right,
            facecolor="#ffd92f",
            edgecolor="#e6550d",
            hatch="///",
            linewidth=2,
            alpha=0.9,
            label="Net shoreline progradation in nourishment step",
            zorder=4,
        )

    migration = bool(values["migration_on"][index])
    change = float(values["shoreline_change_cells"][index])
    axis.text(
        0.02,
        0.96,
        f"Beach width: {values['beach_width_m'][index]:.2f} m\n"
        f"Migration permission: {'ON' if migration else 'OFF'}\n"
        f"Barrier3D grid change: {change:g} cell\n"
        f"Storms: {values['storm_count'][index]}  |  Qow: {values['qow'][index]:.1f}",
        transform=axis.transAxes,
        va="top",
        fontsize=9,
        bbox={
            "boxstyle": "round,pad=0.4",
            "facecolor": "white",
            "edgecolor": "#238b45" if migration else "#cb181d",
            "linewidth": 2,
            "alpha": 0.95,
        },
    )
    axis.set_xlim(x_min, x_max)
    axis.set_ylim(0, 50)
    axis.set_yticks([])
    axis.set_xlabel("Cross-shore position (m; landward →)")
    axis.set_title(label, fontweight="bold")


def make_figure(data: dict[str, dict], index: int):
    figure, axes = plt.subplots(2, 2, figsize=(15, 9))
    shoreline_values = np.concatenate([v["shoreline_m"] for v in data.values()])
    dune_values = np.concatenate([v["dune_grid_m"] for v in data.values()])
    road_values = np.array([v["road_position_m"] for v in data.values()], dtype=float)
    finite_roads = road_values[np.isfinite(road_values)]
    x_min = float(min(shoreline_values.min(), dune_values.min()) - 8)
    x_max_candidates = [float(dune_values.max() + 15)]
    if finite_roads.size:
        x_max_candidates.append(float(finite_roads.max() + 8))
    x_limits = (x_min, max(x_max_candidates))

    labels = list(data)
    draw_plan_view(axes[0, 0], labels[0], data[labels[0]], index, x_limits)
    draw_plan_view(axes[0, 1], labels[1], data[labels[1]], index, x_limits)
    axes[0, 1].legend(loc="lower right", fontsize=7, ncol=2)

    time = np.arange(index + 1)
    for label, values in data.items():
        color = COLORS[label]
        axes[1, 0].plot(
            time,
            values["shoreline_m"][: index + 1],
            color=color,
            linewidth=2,
            label=f"{label} shoreline",
        )
        axes[1, 0].plot(
            time,
            values["dune_grid_m"][: index + 1],
            color=color,
            linewidth=1.8,
            linestyle="--",
            label=f"{label} dune grid",
        )
        axes[1, 1].plot(
            time,
            values["beach_width_m"][: index + 1],
            color=color,
            linewidth=2,
            label=label,
        )
        migrated = np.flatnonzero(values["shoreline_change_cells"][: index + 1] != 0)
        axes[1, 0].scatter(
            migrated,
            values["dune_grid_m"][migrated],
            marker="D",
            color=color,
            edgecolor="black",
            linewidth=0.4,
            s=35,
            zorder=5,
        )

    state_count = len(next(iter(data.values()))["shoreline_m"])
    axes[1, 0].set_xlim(0, state_count - 1)
    axes[1, 0].set_ylabel("Cross-shore position (m)")
    axes[1, 0].set_xlabel("Actual Barrier3D time index")
    axes[1, 0].set_title("Actual shoreline and discrete dune-grid positions")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(fontsize=7)
    axes[1, 1].set_xlim(0, state_count - 1)
    axes[1, 1].set_ylim(-2, max(v["beach_width_m"].max() for v in data.values()) + 4)
    axes[1, 1].axhline(0, color="black", linewidth=1)
    axes[1, 1].set_ylabel("Beach width (m)")
    axes[1, 1].set_xlabel("Actual Barrier3D time index")
    axes[1, 1].set_title("Manager-tracked beach width")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(fontsize=8)

    messages = []
    for label, values in data.items():
        messages.extend(key_event_text(label, values, index))
    subtitle = " | ".join(messages) if messages else "All movement calculated by Barrier3D"
    figure.suptitle(
        f"Full Barrier3D controlled nourishment — time index {index}\n{subtitle}",
        fontsize=14,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    return figure


def main() -> None:
    data = {
        label: extract(load_cascade(path), label)
        for label, path in NPZ_PATHS.items()
    }
    state_counts = {len(values["shoreline_m"]) for values in data.values()}
    if len(state_counts) != 1:
        raise ValueError(f"Saved runs have different state counts: {state_counts}")
    state_count = state_counts.pop()

    key_frames = set()
    for values in data.values():
        key_frames.update(np.flatnonzero(values["nourishment"]).tolist())
        key_frames.update(np.flatnonzero(values["shoreline_change_cells"] != 0).tolist())
        zero = np.flatnonzero(np.isclose(values["beach_width_m"], 0.0))
        if zero.size:
            key_frames.add(int(zero[0]))
    frame_indices = []
    for index in range(state_count):
        frame_indices.extend([index] * (3 if index in key_frames else 1))

    final_figure = make_figure(data, state_count - 1)
    final_figure.savefig(PNG_PATH, dpi=190, bbox_inches="tight")
    plt.close(final_figure)

    animation_figure = plt.figure(figsize=(15, 9))

    def update(frame_number):
        animation_figure.clear()
        source = make_figure(data, frame_indices[frame_number])
        source.canvas.draw()
        image = np.asarray(source.canvas.buffer_rgba())
        plt.close(source)
        axis = animation_figure.add_axes((0, 0, 1, 1))
        axis.imshow(image)
        axis.axis("off")
        return []

    animation = FuncAnimation(
        animation_figure,
        update,
        frames=len(frame_indices),
        interval=1000,
        repeat=True,
        blit=False,
    )
    animation.save(GIF_PATH, writer=PillowWriter(fps=1), dpi=100)
    plt.close(animation_figure)
    print(f"Wrote {PNG_PATH}")
    print(f"Wrote {GIF_PATH}")


if __name__ == "__main__":
    main()
