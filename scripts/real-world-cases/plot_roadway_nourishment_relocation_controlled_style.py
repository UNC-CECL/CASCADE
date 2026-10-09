"""Plot the saved nourishment-plus-relocation run in controlled-test style."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyArrowPatch, Rectangle

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_nourishment_relocation"
    / "actual_barrier3d_no_storm_integration"
)
NPZ_PATH = OUTPUT_DIR / "roadway_nourishment_relocation.npz"
FIGURE_PATH = (
    OUTPUT_DIR / "roadway_nourishment_relocation_controlled_style.png"
)
GIF_PATH = (
    OUTPUT_DIR / "roadway_nourishment_relocation_controlled_style_slow.gif"
)
INITIAL_BEACH_WIDTH_M = 30.0


def load_data():
    cascade = np.load(NPZ_PATH, allow_pickle=True)["cascade"][()]
    barrier = cascade.barrier3d[0]
    roadway = cascade.roadways[0]
    count = len(barrier.x_s_TS)
    shoreline_change = np.asarray(barrier.ShorelineChangeTS[:count], dtype=float)
    initial_dune_grid_m = (
        np.floor(barrier.x_s_TS[0] + INITIAL_BEACH_WIDTH_M / 10.0) * 10.0
    )
    dune_grid_m = initial_dune_grid_m + np.cumsum(-shoreline_change) * 10.0
    setback = np.asarray(roadway._road_setback_TS[:count], dtype=float)
    width = np.asarray(roadway._road_width_TS[:count], dtype=float)
    road = dune_grid_m + setback
    road[width <= 0] = np.nan
    values = {
        "time": np.arange(count),
        "shoreline_m": np.asarray(barrier.x_s_TS[:count], dtype=float) * 10.0,
        "beach_width_m": np.asarray(roadway.beach_width[:count], dtype=float),
        "dune_grid_m": dune_grid_m,
        "shoreline_change_cells": shoreline_change,
        "migration_on": np.asarray(
            roadway.dune_migration_on[:count], dtype=bool
        ),
        "road_setback_m": setback,
        "road_width_m": width,
        "road_m": road,
        "nourishment": np.asarray(roadway.nourishment_TS[:count], dtype=bool),
        "nourishment_volume": np.asarray(
            roadway.nourishment_volume_TS[:count], dtype=float
        ),
        "relocation": np.asarray(
            roadway.triggered_relocation_TS[:count], dtype=bool
        ),
        "storm_count": np.asarray(barrier._StormCount[:count], dtype=int),
        "qow": np.asarray(barrier.QowTS[:count], dtype=float),
    }
    return cascade, barrier, values


def event_indices(values):
    return {
        "nourishment": int(np.flatnonzero(values["nourishment"])[0]),
        "zero_beach": int(
            np.flatnonzero(np.isclose(values["beach_width_m"], 0))[0]
        ),
        "first_migration": int(
            np.flatnonzero(values["shoreline_change_cells"] < 0)[0]
        ),
        "first_relocation": int(np.flatnonzero(values["relocation"])[0]),
    }


def plot_limits(barrier, values):
    road_landward = values["road_m"] + values["road_width_m"]
    plan_min = float(np.min(values["shoreline_m"]) - 8)
    plan_max = float(max(np.nanmax(road_landward), np.max(values["dune_grid_m"])) + 20)
    domain_max = max(
        values["dune_grid_m"][index]
        + (2 + np.asarray(barrier.DomainTS[index]).shape[0]) * 10
        for index in range(len(values["time"]))
    )
    return (plan_min, plan_max), (
        float(np.min(values["dune_grid_m"]) - 5),
        float(domain_max),
    )


def draw_plan(axis, values, index, x_limits):
    shoreline = float(values["shoreline_m"][index])
    dune = float(values["dune_grid_m"][index])
    road = float(values["road_m"][index])
    road_width = float(values["road_width_m"][index])
    x_min, x_max = x_limits

    axis.axvspan(x_min, shoreline, color="#9ecae1", alpha=0.95, label="Ocean")
    if shoreline < dune:
        axis.axvspan(
            shoreline, dune, color="#fdd49e", alpha=0.95, label="Beach"
        )
    axis.axvspan(dune, x_max, color="#c7e9c0", alpha=0.92, label="Interior")
    axis.axvline(shoreline, color="#08519c", linewidth=3, label="Shoreline")
    axis.axvspan(
        dune - 0.8,
        dune + 0.8,
        color="#8c510a",
        alpha=0.95,
        label="Dune-grid edge",
    )
    if np.isfinite(road) and road_width > 0:
        axis.axvspan(
            road,
            road + road_width,
            color="#525252",
            alpha=0.92,
            label="Current road",
        )
    if values["nourishment"][index] and index > 0:
        axis.axvline(
            values["shoreline_m"][index - 1],
            color="#2171b5",
            linestyle="--",
            linewidth=2,
            label="Pre-nourishment shoreline",
        )
    if values["relocation"][index] and index > 0:
        old_road = float(values["road_m"][index - 1])
        old_width = float(values["road_width_m"][index - 1])
        axis.axvspan(
            old_road,
            old_road + old_width,
            facecolor="none",
            edgecolor="#e31a1c",
            linestyle="--",
            linewidth=2.2,
            label="Previous road",
        )
        axis.add_patch(
            FancyArrowPatch(
                (old_road + old_width / 2, 0.65),
                (road + road_width / 2, 0.65),
                arrowstyle="-|>",
                mutation_scale=18,
                linewidth=2,
                color="#e31a1c",
            )
        )

    migration = "ON" if values["migration_on"][index] else "OFF"
    status = (
        f"Beach width: {values['beach_width_m'][index]:.2f} m\n"
        f"Migration permission: {migration}\n"
        f"Barrier3D grid change: "
        f"{values['shoreline_change_cells'][index]:.0f} cell\n"
        f"Road setback: {values['road_setback_m'][index]:.0f} m\n"
        f"Storms: {values['storm_count'][index]}  |  "
        f"Qow: {values['qow'][index]:.1f}"
    )
    edge_color = "#e31a1c" if values["relocation"][index] else "#238b45"
    axis.text(
        0.02,
        0.96,
        status,
        transform=axis.transAxes,
        va="top",
        fontsize=10,
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": edge_color,
            "linewidth": 2,
            "alpha": 0.93,
        },
    )
    axis.set_xlim(x_limits)
    axis.set_ylim(0, 1)
    axis.set_yticks([])
    axis.set_xlabel("Cross-shore position (m; landward →)")
    axis.set_title("Controlled cross-shore view", fontweight="bold")


def domain_at(barrier, values, index):
    dunes = (
        np.asarray(barrier.DuneDomain[index], dtype=float).T + barrier.BermEl
    ) * 10.0
    interior = np.asarray(barrier.DomainTS[index], dtype=float) * 10.0
    domain = np.vstack((dunes, interior))
    x_edges = np.arange(domain.shape[1] + 1, dtype=float) * 10.0
    y_edges = (
        values["dune_grid_m"][index]
        + np.arange(domain.shape[0] + 1, dtype=float) * 10.0
    )
    return x_edges, y_edges, domain


def draw_domain(axis, barrier, values, index, y_limits):
    x_edges, y_edges, domain = domain_at(barrier, values, index)
    axis.pcolormesh(
        x_edges,
        y_edges,
        domain,
        cmap="terrain",
        vmin=-1,
        vmax=3.5,
        shading="flat",
    )
    dune = float(values["dune_grid_m"][index])
    road = float(values["road_m"][index])
    road_width = float(values["road_width_m"][index])
    axis.axhline(dune, color="#8c510a", linewidth=2.2, label="Dune-grid edge")
    if np.isfinite(road) and road_width > 0:
        axis.add_patch(
            Rectangle(
                (0, road),
                x_edges[-1],
                road_width,
                facecolor="#525252",
                edgecolor="black",
                alpha=0.78,
                label="Current road",
            )
        )
    if values["relocation"][index] and index > 0:
        old_road = float(values["road_m"][index - 1])
        old_width = float(values["road_width_m"][index - 1])
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
    axis.set_ylim(y_limits)
    axis.set_xlabel("Alongshore distance (m)")
    axis.set_ylabel("Cross-shore position (m; landward ↑)")
    axis.set_title("Saved Barrier3D elevation domain", fontweight="bold")
    axis.legend(loc="upper right", fontsize=8)


def draw_histories(axes, values, index):
    position_axis, management_axis = axes
    visible = slice(0, index + 1)
    time = values["time"][visible]
    all_time = values["time"]
    migrations = np.flatnonzero(
        values["shoreline_change_cells"][: index + 1] < 0
    )
    relocations = np.flatnonzero(values["relocation"][: index + 1])

    position_axis.plot(
        time,
        values["shoreline_m"][visible],
        color="#08519c",
        linewidth=2.4,
        label="Shoreline",
    )
    position_axis.plot(
        time,
        values["dune_grid_m"][visible],
        color="#8c510a",
        linestyle="--",
        linewidth=2.4,
        label="Dune-grid edge",
    )
    position_axis.plot(
        time,
        values["road_m"][visible],
        color="#525252",
        linewidth=2.6,
        label="Road seaward edge",
    )
    if migrations.size:
        position_axis.scatter(
            migrations,
            values["dune_grid_m"][migrations],
            marker="D",
            s=38,
            color="#d95f0e",
            edgecolor="black",
            linewidth=0.5,
            zorder=5,
            label="Barrier3D grid migration",
        )
    if relocations.size:
        position_axis.scatter(
            relocations,
            values["road_m"][relocations],
            marker="*",
            s=170,
            color="#e31a1c",
            edgecolor="black",
            zorder=6,
            label="Road relocation",
        )
    position_axis.set_xlim(0, all_time[-1])
    position_axis.set_ylim(
        min(values["shoreline_m"]) - 5,
        np.nanmax(values["road_m"] + values["road_width_m"]) + 8,
    )
    position_axis.set_xlabel("Actual Barrier3D time index")
    position_axis.set_ylabel("Cross-shore position (m)")
    position_axis.set_title(
        "Actual shoreline, dune-grid, and road positions", fontweight="bold"
    )
    position_axis.grid(alpha=0.25)
    position_axis.legend(loc="best", fontsize=8, ncol=2)

    management_axis.plot(
        time,
        values["beach_width_m"][visible],
        color="#d95f0e",
        linewidth=2.6,
        label="Beach width",
    )
    management_axis.step(
        time,
        values["road_setback_m"][visible],
        where="post",
        color="#252525",
        linewidth=2.4,
        label="Road setback",
    )
    off = ~values["migration_on"]
    management_axis.fill_between(
        all_time,
        0,
        1,
        where=off,
        transform=management_axis.get_xaxis_transform(),
        color="#cbc9e2",
        alpha=0.35,
        label="Dune migration OFF",
    )
    if relocations.size:
        management_axis.scatter(
            relocations,
            values["road_setback_m"][relocations],
            marker="*",
            s=170,
            color="#e31a1c",
            edgecolor="black",
            zorder=6,
            label="Road relocation",
        )
    management_axis.set_xlim(0, all_time[-1])
    management_axis.set_ylim(-2, max(np.max(values["beach_width_m"]) + 5, 55))
    management_axis.set_xlabel("Actual Barrier3D time index")
    management_axis.set_ylabel("Width or setback (m)")
    management_axis.set_title(
        "Manager-tracked beach width and road setback", fontweight="bold"
    )
    management_axis.grid(alpha=0.25)
    management_axis.legend(loc="best", fontsize=8)


def event_message(values, indices, index):
    if index == indices["nourishment"]:
        return "NOURISH NOW: 100 m³/m; shoreline progrades; dune grid stays put"
    if index == indices["zero_beach"]:
        return "BEACH WIDTH = 0: dune migration permission turns ON"
    if index == indices["first_migration"]:
        return "BARRIER3D: first 10 m landward dune-grid migration"
    if values["relocation"][index]:
        return "ROAD RELOCATED: dune crossed road; setback resets to 20 m"
    return "Zero storms and zero overwash; Barrier3D calculates evolution"


def draw_figure(figure, axes, barrier, values, index, limits):
    for axis in axes.flat:
        axis.clear()
    plan_limits, domain_limits = limits
    draw_plan(axes[0, 0], values, index, plan_limits)
    draw_domain(axes[0, 1], barrier, values, index, domain_limits)
    draw_histories((axes[1, 0], axes[1, 1]), values, index)
    indices = event_indices(values)
    figure.suptitle(
        f"RoadwayManager controlled nourishment + relocation — time index {index}\n"
        f"{event_message(values, indices, index)}",
        fontsize=15,
        fontweight="bold",
    )
    figure.tight_layout(rect=(0, 0, 1, 0.91))


def main():
    cascade, barrier, values = load_data()
    limits = plot_limits(barrier, values)
    indices = event_indices(values)
    static_index = indices["first_relocation"]

    figure, axes = plt.subplots(2, 2, figsize=(16, 10))
    draw_figure(figure, axes, barrier, values, static_index, limits)
    figure.savefig(FIGURE_PATH, dpi=190, bbox_inches="tight")
    plt.close(figure)

    key_frames = set(indices.values())
    key_frames.update(np.flatnonzero(values["relocation"]).tolist())
    frame_indices = []
    for index in values["time"]:
        frame_indices.extend([int(index)] * (6 if index in key_frames else 1))

    figure, axes = plt.subplots(2, 2, figsize=(16, 10))

    def update(frame_number):
        draw_figure(
            figure,
            axes,
            barrier,
            values,
            frame_indices[frame_number],
            limits,
        )
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
    print(FIGURE_PATH)
    print(GIF_PATH)


if __name__ == "__main__":
    main()
