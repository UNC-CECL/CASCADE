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

__all__ = ["plot_island_orientation", "plot_initialization_planview",
           "plot_sea_level_rise"]

GIS_AXIS_LABEL = "GIS domain (S -> N, Cape Point to Pea Island)"


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
        2, 1, figsize=(13, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]})

    for year in sorted(offsets_by_year):
        is_active = year == active_year
        ax_offset.plot(gis_ids, offsets_by_year[year][real] * DAM_TO_M,
                       lw=2.2 if is_active else 1.2,
                       alpha=1.0 if is_active else 0.5,
                       label=f"{year}" + ("  (active)" if is_active else ""))

    ax_offset.set_ylabel("Island offset (m)")
    ax_offset.set_title("Island orientation: cross-shore starting position "
                        "by domain")
    ax_offset.legend()

    years = sorted(offsets_by_year)
    if len(years) == 2:
        earlier, later = years
        difference_m = ((offsets_by_year[later][real]
                         - offsets_by_year[earlier][real]) * DAM_TO_M)
        ax_diff.plot(gis_ids, difference_m, color="#c0392b", lw=1.4)
        ax_diff.axhline(0, color="k", lw=0.8)
        ax_diff.set_ylabel(f"{later} - {earlier} (m)")
        # NOT shoreline change: each year is zeroed on its own most-seaward
        # domain, so this carries a constant offset. Pattern only; section
        # 9.4 rebuilds the real change from the raw transect files.
        ax_diff.set_title(f"Offset-file difference, pattern only -- not "
                          f"shoreline change "
                          f"(mean {difference_m.mean():+.1f} m)")

    ax_diff.set_xlabel(GIS_AXIS_LABEL)
    for ax in (ax_offset, ax_diff):
        ax.grid(alpha=0.3)
    fig.tight_layout()
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

    fig, axes = plt.subplots(2, 1, figsize=(16, 10))
    for ax, with_buffers in zip(axes, (False, True)):
        canvas, col_starts, cells, first_real = init_planview.build_canvas(
            domain_grids, offset_cells, geometry,
            include_buffers=with_buffers, config=config)
        init_planview.plot_canvas(
            canvas, col_starts, cells, first_real, geometry,
            title=f"CASCADE initialization, {start_year} orientation"
                  + ("  |  with buffers" if with_buffers
                     else "  |  real domains only"),
            ax=ax, include_buffers=with_buffers,
            xlabel=GIS_AXIS_LABEL, config=config)
        if verbose:
            print(f"{'with buffers' if with_buffers else 'real only  '}: "
                  f"canvas {canvas.shape}")

    fig.tight_layout()
    return fig


def plot_sea_level_rise(periods, active_year):
    """Plots cumulative RSLR for each period from its own start year.

    Args:
        periods: Mapping of start year to a period config dict.
        active_year: The start year currently selected.

    Returns:
        The matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(9, 4.5))

    for year in sorted(periods):
        period = periods[year]
        run_years = period["end_year"] - year
        calendar_years = np.arange(year, period["end_year"] + 1)
        cumulative_m = (calendar_years - year) * period["sea_level_rise_rate"]
        is_active = year == active_year
        ax.plot(calendar_years, cumulative_m,
                lw=2.2 if is_active else 1.2,
                alpha=1.0 if is_active else 0.5,
                label=f"{year}-{period['end_year']}  "
                      f"{period['sea_level_rise_rate']} m/yr  "
                      f"({cumulative_m[-1]:.2f} m over {run_years} yr)"
                      + ("  (active)" if is_active else ""))

    ax.set_xlabel("Calendar year")
    ax.set_ylabel("Cumulative RSLR (m)")
    ax.set_title("Relative sea level rise by period")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig
