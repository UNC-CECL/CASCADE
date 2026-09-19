#!/usr/bin/env python3
"""Pre-run QC figures for the Hatteras hindcast setup.

WHY THIS MODULE EXISTS
    These three figures answer "does the initial condition look right" before
    a run is started: which way the island is oriented, what the
    initialization surface looks like in plan view, and how much sea level
    rises over the period. They are worth having and cost ~125 lines of
    matplotlib to define, so the definitions live here and the notebook keeps
    one call each.

    Nothing downstream reads their output. Skipping them changes no result --
    which is exactly why they belong out of the file that describes the run.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt

from cascade_pipeline.hindcast import DAM_TO_M
from cascade_pipeline.plotting import init_planview

# `scripts/` is on sys.path already -- cascade_pipeline lives inside it.
from site_layer.hat_figure_style import (
    C, C_1984, C_1997, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style,
    figsize, open_frame,
)

apply_style()

__all__ = ["plot_island_orientation", "plot_initialization_planview",
           "plot_sea_level_rise"]

# The one alongshore axis label the whole project uses. The endpoints this
# constant used to carry ("Cape Point to Pea Island") were the only CORRECT
# pair anywhere in the repo, so they are not lost: ENDPOINT_NOTE states them
# beside each figure's title instead.
GIS_AXIS_LABEL = DOMAIN_AXIS_LABEL
ENDPOINT_NOTE = "domain 1 Cape Point, domain 90 Pea Island"


def plot_island_orientation(offsets_by_year, active_year, geometry):
    """Plots island offset by GIS domain for every period, active one bold.

    Args:
        offsets_by_year: Mapping of start year to a padded offset array (dam).
        active_year: The start year currently selected.
        geometry: DomainGeometry used to slice out the real domains.

    Returns:
        The matplotlib Figure.
    """
    real = slice(geometry.start_real_index, geometry.end_real_index)
    gis_ids = np.arange(geometry.first_gis_id, geometry.last_gis_id + 1)

    fig, (ax_offset, ax_diff) = plt.subplots(
        2, 1, figsize=figsize("double", aspect=0.62), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]}, constrained_layout=True)

    # Two start years are two vintages: the earlier red, the later blue.
    years_sorted = sorted(offsets_by_year)
    vintage = {y: (C_1984 if i == 0 else C_1997)
               for i, y in enumerate(years_sorted)}
    for year in years_sorted:
        is_active = year == active_year
        ax_offset.plot(gis_ids, offsets_by_year[year][real] * DAM_TO_M,
                       color=vintage.get(year, C["BASE"]),
                       lw=1.8 if is_active else 1.0,
                       alpha=1.0 if is_active else 0.45,
                       label=f"{year}" + (" (this run)" if is_active else ""))

    ax_offset.set_ylabel("island offset (m)")
    _title(ax_offset, 0, "cross-shore starting position by domain")
    ax_offset.set_title(ENDPOINT_NOTE, loc="right", fontsize=7.5,
                        color=INK_MUTED)
    ax_offset.legend(frameon=False)

    years = sorted(offsets_by_year)
    if len(years) == 2:
        earlier, later = years
        difference_m = ((offsets_by_year[later][real]
                         - offsets_by_year[earlier][real]) * DAM_TO_M)
        ax_diff.plot(gis_ids, difference_m, color=C["ACCENT"], lw=1.3)
        ax_diff.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 3)))
        ax_diff.set_ylabel(f"{later} minus {earlier} (m)")
        # NOT shoreline change: each year is zeroed on its own most-seaward
        # domain, so this carries a constant offset. Pattern only; section
        # 9.4 rebuilds the real change from the raw transect files. The mean
        # stays on the canvas: these figures are handed back to a notebook and
        # never written to disk, so there is no CAPTIONS.md to hold it.
        _title(ax_diff, 1, "offset-file difference, pattern only")
        ax_diff.set_title(f"mean {difference_m.mean():+.1f} m", loc="right",
                          fontsize=7.5, color=INK_MUTED)

    ax_diff.set_xlabel(GIS_AXIS_LABEL)
    for ax in (ax_offset, ax_diff):
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        open_frame(ax)
    return fig


def plot_initialization_planview(elevation_file_paths, island_offset_dam,
                                 geometry, start_year, config=None,
                                 verbose=True):
    """Draws the initialization surface in plan view, with and without buffers.

    Args:
        elevation_file_paths: Padded list of domain elevation .npy paths.
        island_offset_dam: Padded offsets in decameters.
        geometry: DomainGeometry.
        start_year: Start year, for the titles.
        config: PlanViewConfig, or None for the extractor defaults
            (200 rows, -3.0 m water).
        verbose: Whether to print each canvas shape.

    Returns:
        The matplotlib Figure.
    """
    config = config or init_planview.PlanViewConfig()
    offset_cells = np.round(
        island_offset_dam * DAM_TO_M / config.cell_size_m).astype(int)
    domain_grids = init_planview.load_domain_grids(elevation_file_paths,
                                                   config)

    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.4),
                             constrained_layout=True)
    for ax, with_buffers in zip(axes, (False, True)):
        canvas, col_starts, cells, first_real = init_planview.build_canvas(
            domain_grids, offset_cells, geometry,
            include_buffers=with_buffers, config=config)
        init_planview.plot_canvas(
            canvas, col_starts, cells, first_real, geometry,
            title=("with buffer domains" if with_buffers
                   else "real domains only"),
            ax=ax, include_buffers=with_buffers,
            xlabel=GIS_AXIS_LABEL, config=config)
        if not with_buffers:      # said once, on the upper panel
            ax.set_title(f"{start_year} initialization surface  ·  "
                         f"{ENDPOINT_NOTE}", loc="right", fontsize=7.5,
                         color=INK_MUTED)
        if verbose:
            print(f"{'with buffers' if with_buffers else 'real only  '}: "
                  f"canvas {canvas.shape}")

    return fig


def plot_sea_level_rise(periods, active_year):
    """Plots cumulative RSLR for each period from its own start year.

    Args:
        periods: Mapping of start year to a period config dict.
        active_year: The start year currently selected.

    Returns:
        The matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)

    years_sorted = sorted(periods)
    vintage = {y: (C_1984 if i == 0 else C_1997)
               for i, y in enumerate(years_sorted)}
    for year in years_sorted:
        period = periods[year]
        run_years = period["end_year"] - year
        calendar_years = np.arange(year, period["end_year"] + 1)
        cumulative_m = (calendar_years - year) * period["sea_level_rise_rate"]
        is_active = year == active_year
        ax.plot(calendar_years, cumulative_m,
                color=vintage.get(year, C["BASE"]),
                lw=1.8 if is_active else 1.0,
                alpha=1.0 if is_active else 0.45,
                label=f"{year}–{period['end_year']}, "
                      f"{period['sea_level_rise_rate']} m/yr "
                      f"({cumulative_m[-1]:.2f} m over {run_years} yr)"
                      + (" (this run)" if is_active else ""))

    ax.set_xlabel("calendar year")
    ax.set_ylabel("cumulative RSLR (m)")
    # No panel letter: one panel, nothing to cite it against.
    ax.set_title("relative sea level rise by period", loc="left")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    ax.legend(frameon=False, loc="upper left")
    return fig
