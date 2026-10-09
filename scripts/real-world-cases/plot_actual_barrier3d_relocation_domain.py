"""Plot the saved numerical domain used by the road-relocation test.

This is visualization only: it loads the completed NPZ and never calls a model
update or changes model state.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SOURCE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = SOURCE_ROOT.parent
sys.path.insert(0, str(SOURCE_ROOT))

OUTPUT_DIR = (
    WORKSPACE_ROOT
    / "CASCADE"
    / "output"
    / "roadway_relocation"
    / "actual_barrier3d_synthetic_test"
)
NPZ_PATH = OUTPUT_DIR / "actual_barrier3d_road_relocation.npz"
FIGURE_PATH = OUTPUT_DIR / "actual_barrier3d_relocation_domain.png"
GIF_PATH = OUTPUT_DIR / "actual_barrier3d_relocation_domain_slow.gif"
INITIAL_BEACH_WIDTH_M = 30.0

# Water gets a separate blue class. Positive elevations use terrain colors.
TERRAIN = plt.get_cmap("terrain")
COLORS = ["#8ecae6"] + [TERRAIN(value) for value in np.linspace(0.24, 0.95, 18)]
CMAP = ListedColormap(COLORS)
BOUNDS = np.concatenate(([-3.1, 0.0], np.linspace(0.2, 3.6, 18)))
NORM = BoundaryNorm(BOUNDS, CMAP.N)


def load_run():
    return np.load(NPZ_PATH, allow_pickle=True)["cascade"][()]


def geometry(cascade) -> dict[str, np.ndarray]:
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
    return {
        "dune_grid_m": dune_grid_m,
        "road_m": road,
        "road_width_m": width,
        "setback_m": setback,
        "relocation": np.asarray(
            roadway.triggered_relocation_TS[:count], dtype=bool
        ),
        "shoreline_change_m": shoreline_change * 10.0,
        "storm_count": np.asarray(barrier._StormCount[:count], dtype=int),
        "overwash": np.asarray(barrier.QowTS[:count], dtype=float),
    }


def elevation_domain(cascade, values: dict[str, np.ndarray], index: int):
    """Return the exact saved dune and interior elevation cells in metres MHW."""

    barrier = cascade.barrier3d[0]
    dunes = (
        np.asarray(barrier.DuneDomain[index], dtype=float).T + barrier.BermEl
    ) * 10.0
    interior = np.asarray(barrier.DomainTS[index], dtype=float) * 10.0
    domain = np.vstack((dunes, interior))
    cross_shore_edges = (
        values["dune_grid_m"][index]
        + np.arange(domain.shape[0] + 1, dtype=float) * 10.0
    )
    alongshore_edges = np.arange(domain.shape[1] + 1, dtype=float) * 10.0
    return alongshore_edges, cross_shore_edges, domain


def draw_domain(axis, cascade, values, index, show_previous_road=False):
    x_edges, y_edges, domain = elevation_domain(cascade, values, index)
    mesh = axis.pcolormesh(
        x_edges,
        y_edges,
        domain,
        cmap=CMAP,
        norm=NORM,
        shading="flat",
    )
    axis.axhline(
        values["dune_grid_m"][index],
        color="#7f3b08",
        linewidth=2.2,
        label="Seaward dune-grid edge",
    )

    road = values["road_m"][index]
    width = values["road_width_m"][index]
    if np.isfinite(road) and width > 0:
        axis.add_patch(
            Rectangle(
                (x_edges[0], road),
                x_edges[-1] - x_edges[0],
                width,
                facecolor="#404040",
                edgecolor="black",
                linewidth=1.2,
                alpha=0.75,
                label="Current road footprint",
            )
        )
    if show_previous_road and index > 0:
        prior_road = values["road_m"][index - 1]
        prior_width = values["road_width_m"][index - 1]
        axis.add_patch(
            Rectangle(
                (x_edges[0], prior_road),
                x_edges[-1] - x_edges[0],
                prior_width,
                facecolor="none",
                edgecolor="#e31a1c",
                linestyle="--",
                linewidth=2.3,
                label="Previous road footprint",
            )
        )

    axis.set_xlim(x_edges[0], x_edges[-1])
    axis.set_ylim(y_edges[0] - 5, y_edges[-1])
    axis.set_xlabel("Alongshore distance (m)")
    axis.set_ylabel("Absolute cross-shore position (m; landward ↑)")
    event_label = " — ROAD RELOCATED" if values["relocation"][index] else ""
    axis.set_title(
        f"Time {index}{event_label}\n"
        f"dune-grid change {values['shoreline_change_m'][index]:.0f} m; "
        f"road setback {values['setback_m'][index]:.0f} m"
    )
    return mesh


def plot_snapshots(cascade, values):
    events = np.flatnonzero(values["relocation"])
    indices = [0]
    for event in events[:2]:
        indices.extend((int(event - 1), int(event)))

    figure, axes = plt.subplots(2, 3, figsize=(18, 13))
    axes = axes.flat
    mesh = None
    for axis, index in zip(axes, indices):
        mesh = draw_domain(
            axis,
            cascade,
            values,
            index,
            show_previous_road=bool(values["relocation"][index]),
        )

    info_axis = axes[-1]
    info_axis.axis("off")
    info_axis.text(
        0.02,
        0.96,
        "What is displayed",
        fontsize=15,
        fontweight="bold",
        va="top",
    )
    info_axis.text(
        0.02,
        0.87,
        "• Exact saved Barrier3D elevation cells\n"
        "• Two dune rows followed by the interior grid\n"
        "• Gray overlay: roadway footprint\n"
        "• Dashed red: road location before relocation\n"
        "• Brown line: seaward edge of the dune grid\n\n"
        "Input elevation grid:\n"
        "b3d_pt75_3284yrs_low-elevations.csv\n"
        "Initial grid: 26 cross-shore × 50 alongshore cells\n"
        "Cell size: 10 m × 10 m\n"
        "Alongshore length: 500 m",
        fontsize=12,
        va="top",
        linespacing=1.45,
    )
    info_axis.legend(
        handles=[
            Patch(facecolor="#8ecae6", label="Cell elevation ≤ 0 m MHW"),
            Patch(facecolor=TERRAIN(0.55), label="Cell elevation > 0 m MHW"),
            Patch(facecolor="#404040", label="Current road footprint"),
            Patch(facecolor="none", edgecolor="#e31a1c", linestyle="--", label="Previous road footprint"),
        ],
        loc="lower left",
        bbox_to_anchor=(0.02, 0.17),
    )
    colorbar_axis = info_axis.inset_axes([0.04, 0.04, 0.88, 0.045])
    colorbar = figure.colorbar(mesh, cax=colorbar_axis, orientation="horizontal")
    colorbar.set_label("Saved elevation (m MHW)")
    figure.suptitle(
        "CASCADE controlled one-segment numerical test domain\n"
        "Saved Barrier3D model cells before and after natural road relocation",
        fontsize=17,
        fontweight="bold",
    )
    figure.subplots_adjust(top=0.90, bottom=0.07, left=0.07, right=0.90, hspace=0.34, wspace=0.25)
    figure.savefig(FIGURE_PATH, dpi=190, bbox_inches="tight")
    plt.close(figure)


def plot_domain_gif(cascade, values):
    count = len(values["dune_grid_m"])
    events = set(np.flatnonzero(values["relocation"]).tolist())
    selected = set(range(0, count, 3))
    selected.add(count - 1)
    for event in events:
        selected.update(range(max(0, event - 3), min(count, event + 4)))
    frames = []
    for index in sorted(selected):
        frames.extend([index] * (6 if index in events else 1))

    figure, (domain_axis, timeline_axis) = plt.subplots(
        1, 2, figsize=(16, 8), gridspec_kw={"width_ratios": (1, 1.1)}
    )

    def update(frame_number):
        index = frames[frame_number]
        domain_axis.clear()
        timeline_axis.clear()
        draw_domain(
            domain_axis,
            cascade,
            values,
            index,
            show_previous_road=index in events,
        )

        visible = slice(0, index + 1)
        time = np.arange(index + 1)
        timeline_axis.plot(
            time,
            values["dune_grid_m"][visible],
            color="#7f3b08",
            linewidth=2.5,
            label="Dune-grid edge",
        )
        timeline_axis.plot(
            time,
            values["road_m"][visible],
            color="#333333",
            linewidth=2.5,
            label="Road seaward edge",
        )
        passed = np.array(sorted(event for event in events if event <= index), dtype=int)
        if passed.size:
            timeline_axis.scatter(
                passed,
                values["road_m"][passed],
                marker="*",
                s=180,
                color="#e31a1c",
                edgecolor="black",
                zorder=5,
                label="Triggered relocation",
            )
        timeline_axis.set_xlim(0, count - 1)
        timeline_axis.set_ylim(
            np.min(values["dune_grid_m"]) - 10,
            np.nanmax(values["road_m"] + values["road_width_m"]) + 10,
        )
        timeline_axis.set_xlabel("Actual Barrier3D time index")
        timeline_axis.set_ylabel("Absolute cross-shore position (m)")
        timeline_axis.set_title("Model-calculated dune and road positions")
        timeline_axis.grid(alpha=0.25)
        timeline_axis.legend(loc="upper left")

        event_text = "ROAD RELOCATED" if index in events else "Model evolution"
        figure.suptitle(
            f"SAVED BARRIER3D DOMAIN — time {index} — {event_text}\n"
            f"storms: {values['storm_count'][index]} | "
            f"overwash: {values['overwash'][index]:.2f} m³/m | "
            f"dune-grid change: {values['shoreline_change_m'][index]:.0f} m | "
            f"road setback: {values['setback_m'][index]:.0f} m",
            fontsize=13,
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


def main():
    cascade = load_run()
    values = geometry(cascade)
    plot_snapshots(cascade, values)
    plot_domain_gif(cascade, values)
    print(FIGURE_PATH)
    print(GIF_PATH)


if __name__ == "__main__":
    main()
