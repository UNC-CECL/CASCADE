"""Plan-view rendering of a CASCADE initialization surface (t=0).

Composites per-domain Barrier3D topography onto one alongshore canvas, each
domain shifted cross-shore by its BRIE island offset, so the initial island
can be read as a map before a run starts. Nothing here is site-specific --
domain geometry arrives as a DomainGeometry and the file paths and offsets are
supplied by the caller.

Two conventions this module depends on, both Barrier3D's:

- Elevation arrays are stored in decameters on a fixed-size grid cell, and are
  indexed (cross_shore, alongshore) with row 0 seaward.
- The extractor writes topography trimmed to the island interior, so the
  landward water rows are absent from the .npy files. `PlanViewConfig.topo_rows`
  is the untrimmed frame height; the missing rows are refilled with
  `sentinel_water_m` (see RUN_MANIFEST.txt in the extractor's version folder).
"""

import dataclasses

import numpy as np
from matplotlib.colors import FuncNorm
import matplotlib.pyplot as plt


@dataclasses.dataclass(frozen=True)
class PlanViewConfig:
    """Rendering and unit settings for a plan-view initialization figure.

    Attributes:
        topo_rows: Cross-shore rows in the extractor's untrimmed frame.
        sentinel_water_m: Elevation used to refill trimmed water rows.
        cell_size_m: Grid cell size in meters.
        dam_to_m: Decameters-to-meters factor for the stored arrays.
        elev_min_m: Colorbar floor.
        elev_max_m: Colorbar ceiling.
        sea_level_m: Elevation the colormap pivots at.
        sea_level_pos: Colormap position sea level maps to. Below 0.5 pushes
            green onto low above-water terrain instead of burying it in the
            sub-sea-level range.
        ocean_color: Fill for canvas cells no domain covers.
        text_color: Foreground color for titles, labels and ticks.
    """

    topo_rows: int = 200
    sentinel_water_m: float = -3.0
    cell_size_m: float = 10.0
    dam_to_m: float = 10.0
    elev_min_m: float = -1.0
    elev_max_m: float = 4.0
    sea_level_m: float = 0.0
    sea_level_pos: float = 0.35
    ocean_color: str = "#b0cfe8"
    text_color: str = "#1a1a2e"


DEFAULT_PLAN_VIEW = PlanViewConfig()


def pad_cross_shore(elevation_m, config=DEFAULT_PLAN_VIEW):
    """Extends a trimmed interior to the full cross-shore frame.

    Args:
        elevation_m: 2-D (cross_shore, alongshore) elevations in meters.
        config: PlanViewConfig supplying the frame height and water value.

    Returns:
        The array padded (or truncated) to config.topo_rows rows. Rows run
        seaward to landward, so padding is appended on the landward end.
    """
    n_rows = elevation_m.shape[0]
    if n_rows >= config.topo_rows:
        return elevation_m[:config.topo_rows, :]
    water = np.full((config.topo_rows - n_rows, elevation_m.shape[1]),
                    config.sentinel_water_m)
    return np.vstack([elevation_m, water])


def load_domain_grids(elevation_paths, config=DEFAULT_PLAN_VIEW):
    """Loads every domain's elevation array, in meters, on the full frame.

    Args:
        elevation_paths: Padded-order elevation file paths, one per domain.
        config: PlanViewConfig supplying unit and frame settings.

    Returns:
        A list of 2-D arrays in meters, each config.topo_rows tall.
    """
    grids = []
    for path in elevation_paths:
        elevation = np.load(path).astype(float)
        if elevation.ndim == 1:
            elevation = elevation.reshape(-1, 1)
        grids.append(pad_cross_shore(elevation * config.dam_to_m, config))
    return grids


def _warn_if_alongshore_reversed(grids, config, warn_ratio=5.0):
    """Warns when the alongshore axis is reversed within each domain.

    A per-domain alongshore reversal is nearly invisible in a 45 km plan view --
    it reads as roughness -- but it lays every 500 m block down backwards. In
    numbers it is unmistakable: the mean jump across domain seams over the mean
    jump within a domain sits near 1 for a continuous island, and reached 21 on
    the corrected arrays while this function's caller still applied np.fliplr.

    Args:
        grids: The per-domain arrays about to be composited.
        config: PlanViewConfig supplying sentinel_water_m.
        warn_ratio: Ratio above which a warning is printed.

    Returns:
        The seam/inner ratio, or nan when it cannot be computed.
    """
    series = []
    for grid in grids:
        land = (grid > config.sentinel_water_m + 1e-9).sum(axis=0).astype(float)
        if land.size:
            series.append(land)
    if len(series) < 3:
        return float("nan")

    width = min(s.size for s in series)
    stack = np.array([s[:width] for s in series])
    seam = np.abs(stack[:-1, -1] - stack[1:, 0])
    inner = np.abs(np.diff(stack, axis=1))
    if not inner.size or inner.mean() <= 0:
        return float("nan")

    ratio = float(seam.mean() / inner.mean())
    if ratio > warn_ratio:
        print(f"[ALONGSHORE WARNING] plan-view canvas: seam/inner discontinuity "
              f"ratio {ratio:.1f} (expected ~1-2). The alongshore axis is "
              f"reversed within each domain relative to the domain order, so "
              f"every 500 m block is backwards. Either the topography was "
              f"extracted with ALONGSHORE_FLIP = False, or a per-domain "
              f"np.fliplr has been reintroduced in build_canvas.")
    return ratio


def build_canvas(domain_grids, offset_cells, geometry, include_buffers=False,
                 config=DEFAULT_PLAN_VIEW):
    """Composites domains onto one canvas, offset by the island dune line.

    Args:
        domain_grids: Padded-order elevation arrays from load_domain_grids.
        offset_cells: Per-domain cross-shore offset, in whole grid cells.
        geometry: DomainGeometry describing the padded array.
        include_buffers: Whether to draw the buffer domains as well as the
            real ones.
        config: PlanViewConfig supplying the frame height.

    Returns:
        A (canvas, domain_col_starts, cells_per_domain, first_real_idx) tuple.
        Cells no domain covers are NaN. first_real_idx is the position of the
        first real domain within the plotted set.
    """
    if include_buffers:
        plotted_grids, plotted_offsets = domain_grids, offset_cells
        first_real_idx = geometry.start_real_index
    else:
        real = slice(geometry.start_real_index, geometry.end_real_index)
        plotted_grids, plotted_offsets = domain_grids[real], offset_cells[real]
        first_real_idx = 0

    _warn_if_alongshore_reversed(plotted_grids, config)

    canvas_rows = int(max(plotted_offsets)) + config.topo_rows + 5
    total_cols = sum(grid.shape[1] for grid in plotted_grids)
    canvas = np.full((canvas_rows, total_cols), np.nan)

    col_cursor = 0
    domain_col_starts = []
    cells_per_domain = []

    for index, grid in enumerate(plotted_grids):
        n_rows, n_cols = grid.shape
        domain_col_starts.append(col_cursor)
        cells_per_domain.append(n_cols)

        origin = int(plotted_offsets[index])
        row_end = min(origin + n_rows, canvas_rows)
        # No np.fliplr here. It used to reverse the alongshore cells WITHIN each
        # domain, inside this very loop, which put every 500 m block backwards
        # against the ascending domain order. That compensated for the extractor
        # writing the within-domain alongshore order reversed; the extractor now
        # fixes it at load (ALONGSHORE_FLIP), so flipping again double-flips.
        # _warn_if_alongshore_reversed below is the guard.
        canvas[origin:row_end, col_cursor:col_cursor + n_cols] = (
            grid[:row_end - origin, :])
        col_cursor += n_cols

    return canvas, np.array(domain_col_starts), cells_per_domain, first_real_idx


def _elevation_norm(config):
    """Builds the sea-level-shifted normalization for the terrain colormap."""

    def forward(value):
        scaled = np.where(
            value < config.sea_level_m,
            config.sea_level_pos * (value - config.elev_min_m)
            / (config.sea_level_m - config.elev_min_m),
            config.sea_level_pos
            + (1.0 - config.sea_level_pos) * value / config.elev_max_m,
        )
        return np.where(np.isnan(value), np.nan, scaled)

    def inverse(value):
        return np.where(
            value < config.sea_level_pos,
            config.elev_min_m + (value / config.sea_level_pos)
            * (config.sea_level_m - config.elev_min_m),
            (value - config.sea_level_pos) / (1.0 - config.sea_level_pos)
            * config.elev_max_m,
        )

    return FuncNorm((forward, inverse), vmin=config.elev_min_m,
                    vmax=config.elev_max_m)


def plot_canvas(canvas, domain_col_starts, cells_per_domain, first_real_idx,
                geometry, title, ax=None, include_buffers=False,
                xlabel="Domain (S -> N)", colorbar=True,
                config=DEFAULT_PLAN_VIEW):
    """Draws a composited canvas as a plan-view map.

    Args:
        canvas: Canvas from build_canvas.
        domain_col_starts: Column start of each plotted domain.
        cells_per_domain: Alongshore cell count of each plotted domain.
        first_real_idx: Index of the first real domain within the plotted set.
        geometry: DomainGeometry used to label the real domains.
        title: Axes title.
        ax: Axes to draw on. A new figure is created when omitted.
        include_buffers: Whether buffer domains are present, and so whether to
            bracket and label them.
        xlabel: X axis label.
        colorbar: Whether to draw the elevation colorbar. Off for small
            multiples that share one bar with a parent panel.
        config: PlanViewConfig supplying colors and limits.

    Returns:
        The matplotlib Figure the canvas was drawn on.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 5))
    fig = ax.figure

    ax.set_facecolor(config.ocean_color)
    cmap = plt.cm.terrain.copy()
    cmap.set_bad(color=config.ocean_color)

    mesh = ax.pcolormesh(np.ma.masked_invalid(canvas), cmap=cmap,
                         norm=_elevation_norm(config), shading="auto",
                         rasterized=True)

    if colorbar:
        bar = fig.colorbar(mesh, ax=ax, pad=0.01, fraction=0.02)
        bar.set_label("Elevation (m MHW)", color=config.text_color)
        bar.set_ticks([-1, 0, 1, 2, 3, 4])

    # Adaptive so a small multiple (a 7-domain zoom) still gets labels;
    # a fixed step of 5 gives such a panel only two ticks.
    tick_step = max(1, round(geometry.num_real_domains / 18))
    tick_indices = range(0, geometry.num_real_domains, tick_step)
    ax.set_xticks([domain_col_starts[first_real_idx + i]
                   + cells_per_domain[first_real_idx + i] // 2
                   for i in tick_indices])
    ax.set_xticklabels([str(geometry.first_gis_id + i) for i in tick_indices],
                       fontsize=9)
    ax.set_xlim(0, canvas.shape[1])
    ax.set_ylim(0, canvas.shape[0])
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cross-shore cell")
    ax.set_title(title, fontweight="bold", color=config.text_color)

    for index in range(0, geometry.num_real_domains, 10):
        ax.axvline(domain_col_starts[first_real_idx + index] - 0.5,
                   color="#aaaaaa", lw=0.4, alpha=0.5, zorder=2)

    if include_buffers:
        last_real_idx = first_real_idx + geometry.num_real_domains - 1
        real_start = domain_col_starts[first_real_idx]
        real_end = domain_col_starts[last_real_idx] + cells_per_domain[last_real_idx]
        for boundary in (real_start, real_end):
            ax.axvline(boundary, color="#555555", lw=1.2, ls="--", alpha=0.8,
                       zorder=6)
        label = f"buffer\n({geometry.num_buffer_domains} domains, interpolated)"
        for x_mid, y, va in ((real_start / 2, 0.03, "bottom"),
                             ((real_end + canvas.shape[1]) / 2, 0.97, "top")):
            ax.text(x_mid, y, label, transform=ax.get_xaxis_transform(),
                    ha="center", va=va, fontsize=8, color="#333333", zorder=7,
                    bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                              edgecolor="#cccccc", alpha=0.9))

    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    return fig


def plot_initialization(elevation_paths, offset_cells, geometry, title,
                        include_buffers=False, ax=None,
                        config=DEFAULT_PLAN_VIEW, **kwargs):
    """Loads, composites and draws an initialization surface in one call.

    Args:
        elevation_paths: Padded-order elevation file paths, one per domain.
        offset_cells: Per-domain cross-shore offset, in whole grid cells.
        geometry: DomainGeometry describing the padded array.
        title: Axes title.
        include_buffers: Whether to draw the buffer domains.
        ax: Axes to draw on. A new figure is created when omitted.
        config: PlanViewConfig supplying unit and rendering settings.
        **kwargs: Passed through to plot_canvas.

    Returns:
        The matplotlib Figure.
    """
    grids = load_domain_grids(elevation_paths, config)
    canvas, col_starts, cells, first_real = build_canvas(
        grids, offset_cells, geometry, include_buffers, config)
    return plot_canvas(canvas, col_starts, cells, first_real, geometry, title,
                       ax=ax, include_buffers=include_buffers, config=config,
                       **kwargs)
