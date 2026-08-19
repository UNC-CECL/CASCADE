"""Plan-view of the roadway on a CASCADE initialization surface.

Draws NC-12 where CASCADE will actually put it: on the same plan-view canvas
`init_planview` builds, with each domain shifted cross-shore by its island
offset, so the road's setback is read against the island the model runs on
rather than against a map.

The road is drawn per domain as a horizontal bar at `offset + setback`, because
that is literally what CASCADE does -- `road_start = int(setback / cell)` is one
row index applied to every alongshore profile in the domain. A relocated
position, when supplied, is drawn as a second bar so the displacement is
visible against the barrier that has to absorb it.

Nothing here is site-specific: geometry arrives as a DomainGeometry, and the
setbacks, offsets and events are supplied by the caller.
"""

import dataclasses

import numpy as np
import matplotlib.pyplot as plt

from cascade_pipeline.plotting import init_planview


@dataclasses.dataclass(frozen=True)
class IslandSection(object):
    """A named alongshore stretch, drawn as a labelled band.

    Attributes:
        name: Label drawn above the band.
        first_gis: First GIS domain in the stretch.
        last_gis: Last GIS domain in the stretch.
        color: Band fill, or None to label the stretch without shading it.
    """

    name: str
    first_gis: int
    last_gis: int
    color: str = None


@dataclasses.dataclass(frozen=True)
class RoadPlanViewStyle:
    """Colors and sizing for the roadway overlay.

    The three overlay colors carry identity and were checked for colour-vision
    separation (worst adjacent OKLab dE 25.5 against a target of 8; normal
    vision 26.8 against a floor of 15). `relocated` warns on contrast against a
    light surface, so it is always accompanied by a legend entry and never used
    as the only cue.

    Attributes:
        road: Color for the initial road footprint.
        relocated: Color for the post-relocation position.
        drowning: Color for the marker on domains that drown at t=0.
        dune_line: Color for the dune line the setback is measured from.
        bar_linewidth: Line width of a road bar, in points.
        text_color: Foreground for titles and labels.
    """

    road: str = "#B71C1C"
    relocated: str = "#FF8C00"
    drowning: str = "#1565C0"
    dune_line: str = "#1a1a2e"
    bar_linewidth: float = 2.4
    text_color: str = "#1a1a2e"
    section_alpha: float = 0.16
    band_label_size: float = 9.0


DEFAULT_ROAD_STYLE = RoadPlanViewStyle()


def crop_to_data(canvas, pad_rows=6):
    """Finds the rows a canvas actually carries data in.

    The island offsets span kilometres, so a canvas built from them is mostly
    empty. Cropping to the occupied rows is what turns the plan view from a
    thin diagonal ribbon into a readable map; a fixed guess slices the south
    end off, because that is where the offsets are largest.

    Args:
        canvas: The plan-view canvas; uncovered cells are NaN.
        pad_rows: Rows of breathing room to keep either side.

    Returns:
        A (first_row, last_row) tuple, suitable for set_ylim.
    """
    occupied = np.where(~np.isnan(canvas).all(axis=1))[0]
    if not occupied.size:
        return 0, canvas.shape[0]
    return (max(int(occupied[0]) - pad_rows, 0),
            min(int(occupied[-1]) + pad_rows + 1, canvas.shape[0]))


def draw_sections(ax, sections, geometry, col_starts, cells_per_domain,
                  first_real_idx, include_buffers=False,
                  style=DEFAULT_ROAD_STYLE):
    """Shades and labels the named alongshore stretches.

    Args:
        ax: Axes carrying the plan-view canvas.
        sections: IslandSection instances.
        geometry: DomainGeometry describing the padded array.
        col_starts: Canvas column where each plotted domain starts.
        cells_per_domain: Canvas column count for each plotted domain.
        first_real_idx: Position of the first real domain in the plotted set.
        include_buffers: Whether buffers are in the plotted set.
        style: RoadPlanViewStyle.
    """
    pad_offset = 0 if include_buffers else geometry.start_real_index

    def span(gis):
        pad = geometry.gis_to_pad(gis)
        index = pad - pad_offset
        if not 0 <= index < len(col_starts):
            return None
        return col_starts[index], col_starts[index] + cells_per_domain[index]

    for section in sections:
        left, right = span(section.first_gis), span(section.last_gis)
        if left is None or right is None:
            continue
        if section.color:
            ax.axvspan(left[0], right[1], color=section.color,
                       alpha=style.section_alpha, zorder=1, lw=0)
        ax.text((left[0] + right[1]) / 2, 0.965, section.name,
                transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=style.band_label_size, color=style.text_color,
                zorder=10,
                bbox=dict(fc="white", ec="none", alpha=0.82, pad=2))


def road_rows(setbacks_m, offset_cells, geometry, config=None):
    """Converts setbacks into canvas rows, the way bulldoze indexes them.

    Args:
        setbacks_m: Padded per-domain setbacks, in meters.
        offset_cells: Padded per-domain cross-shore offset, in whole cells.
        geometry: DomainGeometry describing the padded array.
        config: PlanViewConfig supplying the cell size.

    Returns:
        A float array of canvas row positions, one per padded domain. Domains
        with a zero setback (no road) are NaN.
    """
    config = config or init_planview.DEFAULT_PLAN_VIEW
    setbacks = np.asarray(setbacks_m, dtype=float)
    offsets = np.asarray(offset_cells, dtype=float)
    rows = offsets + np.floor(setbacks / config.cell_size_m)
    return np.where(setbacks > 0, rows, np.nan)


def overlay_roadway(ax, setbacks_m, offset_cells, geometry, col_starts,
                    cells_per_domain, first_real_idx, relocated_m=None,
                    drowning_gis=(), include_buffers=False,
                    config=None, style=DEFAULT_ROAD_STYLE):
    """Draws the roadway onto an existing plan-view axes.

    Args:
        ax: Axes already carrying an init_planview canvas.
        setbacks_m: Padded per-domain setbacks, in meters.
        offset_cells: Padded per-domain cross-shore offset, in whole cells.
        geometry: DomainGeometry describing the padded array.
        col_starts: Canvas column where each plotted domain starts.
        cells_per_domain: Canvas column count for each plotted domain.
        first_real_idx: Position of the first real domain in the plotted set.
        relocated_m: Optional padded setbacks after relocation, in meters.
        drowning_gis: GIS domains whose road drowns at t=0.
        include_buffers: Whether buffers are in the plotted set.
        config: PlanViewConfig supplying the cell size.
        style: RoadPlanViewStyle.

    Returns:
        The list of legend handles added.
    """
    config = config or init_planview.DEFAULT_PLAN_VIEW
    base_rows = road_rows(setbacks_m, offset_cells, geometry, config)
    moved_rows = (road_rows(relocated_m, offset_cells, geometry, config)
                  if relocated_m is not None else None)

    pad_offset = 0 if include_buffers else geometry.start_real_index
    drowning = set(int(g) for g in drowning_gis)
    handles = []

    for plot_index, col_start in enumerate(col_starts):
        pad = plot_index + pad_offset
        if not 0 <= pad < base_rows.size or np.isnan(base_rows[pad]):
            continue
        col_end = col_start + cells_per_domain[plot_index]
        gis = geometry.first_gis_id + (pad - geometry.start_real_index)

        if moved_rows is not None and not np.isnan(moved_rows[pad]):
            ax.hlines(moved_rows[pad], col_start, col_end,
                      color=style.relocated, lw=style.bar_linewidth, zorder=6)
        ax.hlines(base_rows[pad], col_start, col_end, color=style.road,
                  lw=style.bar_linewidth, zorder=7)
        if gis in drowning:
            ax.plot((col_start + col_end) / 2, base_rows[pad], marker="v",
                    ms=9, color=style.drowning, mec="white", mew=1.1,
                    zorder=9)

    handles.append(plt.Line2D([], [], color=style.road,
                              lw=style.bar_linewidth, label="NC-12 footprint"))
    if moved_rows is not None:
        handles.append(plt.Line2D([], [], color=style.relocated,
                                  lw=style.bar_linewidth,
                                  label="after prescribed relocation"))
    if drowning:
        handles.append(plt.Line2D([], [], lw=0, marker="v", ms=9,
                                  color=style.drowning, mec="white",
                                  label="road drowns at t=0"))
    return handles


def plot_roadway_planview(elevation_paths, offset_cells, setbacks_m, geometry,
                          title, relocated_m=None, drowning_gis=(),
                          sections=(), crop=False, legend=True,
                          include_buffers=False, ax=None, config=None,
                          style=DEFAULT_ROAD_STYLE, **kwargs):
    """Renders the initialization surface with the roadway drawn on it.

    Args:
        elevation_paths: Padded-order elevation file paths, one per domain.
        offset_cells: Padded per-domain cross-shore offset, in whole cells.
        setbacks_m: Padded per-domain setbacks, in meters.
        geometry: DomainGeometry describing the padded array.
        title: Axes title.
        relocated_m: Optional padded setbacks after relocation, in meters.
        drowning_gis: GIS domains whose road drowns at t=0.
        sections: IslandSection instances for the alongshore bands.
        crop: Whether to trim the view to the rows carrying data. The offsets
            span kilometres, so an uncropped canvas is mostly empty.
        legend: Whether to draw the roadway legend.
        include_buffers: Whether to draw the buffer domains.
        ax: Axes to draw on. A new figure is created when omitted.
        config: PlanViewConfig supplying unit and rendering settings.
        style: RoadPlanViewStyle.
        **kwargs: Passed through to init_planview.plot_canvas.

    Returns:
        The matplotlib Figure.
    """
    config = config or init_planview.DEFAULT_PLAN_VIEW
    grids = init_planview.load_domain_grids(elevation_paths, config)
    canvas, col_starts, cells, first_real = init_planview.build_canvas(
        grids, offset_cells, geometry, include_buffers, config)
    figure = init_planview.plot_canvas(
        canvas, col_starts, cells, first_real, geometry, title,
        ax=ax, include_buffers=include_buffers, config=config, **kwargs)

    axes = ax if ax is not None else figure.axes[0]
    handles = overlay_roadway(
        axes, setbacks_m, offset_cells, geometry, col_starts, cells,
        first_real, relocated_m=relocated_m, drowning_gis=drowning_gis,
        include_buffers=include_buffers, config=config, style=style)
    if sections:
        draw_sections(axes, sections, geometry, col_starts, cells, first_real,
                      include_buffers=include_buffers, style=style)
    if crop:
        axes.set_ylim(*crop_to_data(canvas)[::-1] if axes.yaxis_inverted()
                      else crop_to_data(canvas))
    if handles and legend:
        axes.legend(handles=handles, loc="lower left", fontsize=8,
                    framealpha=0.92)
    return figure


def plot_roadway_island(elevation_paths, offset_cells, setbacks_m, geometry,
                        title, relocated_m=None, drowning_gis=(),
                        sections=(), zoom_windows=(), config=None,
                        style=DEFAULT_ROAD_STYLE, figsize=(20, 11),
                        xlabel=None, **kwargs):
    """Renders the island-wide roadway figure, with optional zoom panels.

    The full-island panel is cropped to the occupied rows and drawn with
    `aspect="auto"`, which exaggerates the cross-shore axis. Hatteras is ~45 km
    long and ~1 km wide, so at true aspect the island is a hairline; the
    exaggeration is what makes the road's position readable, and the cross-shore
    axis should be read as indicative rather than to scale.

    Args:
        elevation_paths: Padded-order elevation file paths, one per domain.
        offset_cells: Padded per-domain cross-shore offset, in whole cells.
        setbacks_m: Padded per-domain setbacks, in meters.
        geometry: DomainGeometry describing the padded array.
        title: Figure suptitle.
        relocated_m: Optional padded setbacks after relocation, in meters.
        drowning_gis: GIS domains whose road drowns at t=0.
        sections: IslandSection instances for the alongshore bands.
        zoom_windows: (first_gis, last_gis, label) tuples, each drawn as its
            own panel below the island.
        config: PlanViewConfig supplying unit and rendering settings.
        style: RoadPlanViewStyle.
        figsize: Figure size in inches.
        xlabel: Alongshore axis label.
        **kwargs: Passed through to init_planview.plot_canvas.

    Returns:
        The matplotlib Figure.
    """
    config = config or init_planview.DEFAULT_PLAN_VIEW
    n_zoom = len(zoom_windows)
    figure = plt.figure(figsize=figsize, facecolor="white")
    grid = figure.add_gridspec(2 if n_zoom else 1, max(n_zoom, 1),
                               height_ratios=[2.1, 1.0] if n_zoom else [1.0],
                               hspace=0.30, wspace=0.16)

    island_ax = figure.add_subplot(grid[0, :])
    plot_roadway_planview(
        elevation_paths, offset_cells, setbacks_m, geometry,
        title="", relocated_m=relocated_m, drowning_gis=drowning_gis,
        sections=sections, crop=True, ax=island_ax, config=config,
        style=style, xlabel=xlabel, **kwargs)

    for panel, (first_gis, last_gis, label) in enumerate(zoom_windows):
        zoom_ax = figure.add_subplot(grid[1, panel])
        pads = [geometry.gis_to_pad(g) for g in range(first_gis, last_gis + 1)]
        pads = [p for p in pads if 0 <= p < len(elevation_paths)]
        sub_geometry = dataclasses.replace(
            geometry, num_real_domains=len(pads), num_buffer_domains=0,
            first_gis_id=first_gis)
        plot_roadway_planview(
            [elevation_paths[p] for p in pads],
            np.asarray(offset_cells)[pads],
            np.asarray(setbacks_m)[pads],
            sub_geometry,
            title=f"GIS {first_gis}-{last_gis}  |  {label}",
            relocated_m=(None if relocated_m is None
                         else np.asarray(relocated_m)[pads]),
            drowning_gis=drowning_gis, crop=True, legend=False,
            colorbar=False, ax=zoom_ax, config=config, style=style,
            xlabel="GIS domain", **kwargs)

    figure.suptitle(title, fontsize=14, fontweight="bold",
                    color=style.text_color, y=0.98)
    return figure
