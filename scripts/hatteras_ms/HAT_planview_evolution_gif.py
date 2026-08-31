#!/usr/bin/env python3
"""Animated plan view of the island through a run: elevation, not shoreline.

WHAT WAS MISSING, AND WHY THIS EXISTS
    Every matrix run already writes four GIFs and all four animate the
    SHORELINE -- a one-dimensional cross-shore position per domain, plotted
    against domain number. None of them shows elevation. The plan-view canvas
    that `init_planview` builds, the one the initialization figures use, was
    only ever drawn at t = 0.

    So there was no way to watch the barrier itself evolve: where the interior
    lowers, where overwash reaches, where the island narrows. This draws that
    canvas once per model year.

THE HISTORY IS ALREADY ON DISK -- no re-run is needed
    Barrier3D keeps `DomainTS`, one interior grid per year, and `x_s_TS`, the
    shoreline position per year, for all 120 padded domains. Both survive in
    the run's `.npz`, which is why those files are ~300 MB rather than the
    ~20 KB the shoreline matrix costs. This reads them and plots; it does not
    re-run anything.

TWO THINGS THAT HAVE TO BE RIGHT
    RAGGED GRIDS. `DomainTS[t]` is NOT a fixed shape -- Barrier3D trims water
    rows, so one domain-year is (174, 50) and another (200, 50), and the count
    changes as the island evolves. Each grid is put back on the full frame with
    `pad_cross_shore`, exactly as `load_domain_grids` does for the static
    figure, so the two are directly comparable.

    A FIXED FRAME. `build_canvas` sizes the canvas from the offsets it is
    given, so a per-year canvas is a per-year height and the animation would
    breathe. Offsets are taken against ONE reference for the whole run, and the
    axes get one y limit for every frame, so what moves in the GIF is the
    island and not the camera.

UNITS
    `DomainTS` is in decameters and `x_s_TS` likewise; the plan-view config
    carries `dam_to_m = 10.0` and a 10 m cell, so a decameter of shoreline
    movement is exactly one canvas row. Elevations are converted to meters to
    match the shared colorbar.

Usage:
    python HAT_planview_evolution_gif.py <run_directory> [--fps 3] [--out PATH]

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts",):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import (  # noqa: E402
    HATTERAS_DOMAINS as GEOMETRY, HATTERAS_PERIODS,
    HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)

# Defined the same way HAT_hindcast_1984_2024.py:255 defines it. It is not
# exported by the site config, and the period table's paths are relative to it.
HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
from cascade_pipeline.hindcast import load_island_offset_dam  # noqa: E402
from cascade_pipeline.plotting.init_planview import (  # noqa: E402
    DEFAULT_PLAN_VIEW, build_canvas, pad_cross_shore, plot_canvas)
from cascade_pipeline.plotting.road_planview import overlay_roadway  # noqa: E402

HOLD_FRAMES = 5      # frames held on the final year, so it can be read
ARROW = chr(0x2192)  # a real arrow, not '->'


def _period_start(run_dir, cascade):
    """The hindcast period this run belongs to, for picking its offset file.

    Read from the run NAME first: HATTERAS_PERIODS is keyed on the period start
    year, the run name carries it, and a cascade attribute may not.

    Args:
        run_dir: The run directory, whose name carries the period.
        cascade: The loaded Cascade, used only as a fallback.

    Returns:
        A key present in HATTERAS_PERIODS.

    Raises:
        SystemExit: If no period can be identified.
    """
    for token in Path(run_dir).name.split("_"):
        if token.isdigit() and int(token) in HATTERAS_PERIODS:
            return int(token)
    attribute = getattr(cascade, "_start_year", None)
    if attribute is not None and int(attribute) in HATTERAS_PERIODS:
        return int(attribute)
    raise SystemExit(
        f"cannot tell which hindcast period {Path(run_dir).name} belongs to; "
        f"known starts are {sorted(HATTERAS_PERIODS)}")


def load_history(run_dir):
    """Per-year interior grids and shoreline positions from a run's .npz.

    Args:
        run_dir: A matrix run directory holding exactly one .npz.

    Returns:
        (grids_by_year, offsets_by_year, start_year) where grids_by_year[t] is
        a padded-order list of metre-valued arrays on the full cross-shore
        frame, and offsets_by_year[t] is the matching per-domain row origin.

    Raises:
        SystemExit: If no .npz is present, or it holds no domain history.
    """
    run_dir = Path(run_dir)
    matches = sorted(run_dir.glob("*.npz"))
    if not matches:
        raise SystemExit(
            f"no .npz model state in {run_dir}.\n"
            "The run was made with --no-model-state, and the per-year grids "
            "this needs were never written. Re-run that cell without it.")

    cascade = np.load(matches[0], allow_pickle=True)["cascade"].item()
    barrier3d = cascade._barrier3d
    inner = [getattr(model, "_model", model) for model in barrier3d]

    n_years = len(inner[0].DomainTS)
    config = DEFAULT_PLAN_VIEW

    # THE CANVAS FRAME IS THE ISLAND OFFSET, NOT x_s. This was wrong until
    # 2026-08-31 and the error was visible: x_s is Barrier3D's own shoreline
    # coordinate, and because every domain is a separate Barrier3D it varies by
    # only ~0.55 km across the island. The REAL alongshore geometry lives in
    # the BRIE island-offset file and spans ~6.3 km. Using x_s as the frame
    # compressed the island's diagonal about elevenfold, so the road sat far
    # from the dune line it is measured against and the figure disagreed with
    # HAT_road_island_planview_1984.png, which is built from the offsets.
    #
    # The offset is the STATIC frame; x_s supplies only its CHANGE, so the
    # island is placed where the initialization figures place it and still
    # migrates through the run.
    offset_file = HATTERAS_DATA_BASE / HATTERAS_PERIODS[
        _period_start(run_dir, cascade)]["island_offset_file"]
    base_offset = load_island_offset_dam(offset_file, GEOMETRY)
    x_s_initial = np.array([float(model.x_s_TS[0]) for model in inner])

    # THE ROAD MOVES TOO, and it is the reason to watch this rather than the
    # shoreline GIFs: a setback is measured from the dune line, so a road that
    # never relocates still closes on the ocean as the barrier retreats. The
    # roadway manager keeps _road_setback_TS per domain per year; where a
    # domain carries no road the series is absent and its entry stays 0, which
    # road_rows() renders as NaN rather than as a road at the dune line.
    roadways = getattr(cascade, "_roadways", None) or []

    # WHERE THE ROAD ACTUALLY IS, from this run rather than from a site-wide
    # constant. cascade._roadway_management_module is a per-domain mask of the
    # domains the roadway manager actually manages -- 55 of 90 in the 1984
    # calibBE run, matching that run's own road_management_summary.csv exactly.
    #
    # HATTERAS_FIRST/LAST_ROAD_DOMAIN (9 and 90) is the REACH, not the road:
    # NC-12 is present 9-20, 32-67 and 84-90, with real gaps at 21-31 and
    # 68-83. Drawing the whole reach put road through both gaps. Using the
    # run's own mask also means a scenario with roadway management off draws no
    # road at all, which is correct and which a constant cannot express.
    managed = np.asarray(
        getattr(cascade, "_roadway_management_module", []), dtype=bool)
    if managed.size != len(inner):
        managed = np.zeros(len(inner), dtype=bool)
        for gis in range(HATTERAS_FIRST_ROAD_DOMAIN,
                         HATTERAS_LAST_ROAD_DOMAIN + 1):
            pad = GEOMETRY.gis_to_pad(gis)
            if 0 <= pad < managed.size:
                managed[pad] = True

    setbacks_by_year = []

    grids_by_year, offsets_by_year = [], []
    for year in range(n_years):
        grids_by_year.append([
            pad_cross_shore(np.asarray(model.DomainTS[year], dtype=float)
                            * config.dam_to_m, config)
            for model in inner])
        moved = np.array([float(model.x_s_TS[year]) for model in inner])
        offsets_by_year.append(base_offset + (moved - x_s_initial))
        setbacks_by_year.append(_drawable_setbacks(
            [_setback_at(road, year) for road in roadways]
            if roadways else [0.0] * len(inner), managed))

    start_year = int(getattr(cascade, "_start_year", 0)) or None
    return grids_by_year, offsets_by_year, setbacks_by_year, start_year


# A setback of ZERO means the road sits ON the dune line, not that there is no
# road. road_planview.road_rows() cannot tell those apart -- it returns NaN for
# any setback <= 0 -- so a road whose setback decayed to zero vanished from the
# animation and reappeared if it later relocated seaward.
#
# That is not rare and it is mostly NOT relocation: in the 1984-2004 calibBE
# groin run, 31 of 90 domains cross zero during the run and only five ever
# relocate (GIS 10, 11, 84, 85, 86). The rest are the dune line simply catching
# up with a road that never moved -- GIS 9 decays 40 -> 30 -> 20 -> 10 -> 0 m
# and stays there. Blanking it draws the road as absent exactly when it is most
# exposed, which is backwards.
#
# The road's real extent is a site fact, not something to infer from a setback:
# hatteras_site_config says GIS 9-90 carry NC-12 and "Domains 1-8 (Cape Point)
# have no road in the modelled span". So presence comes from that, and a zero
# setback inside the road reach is nudged just above zero -- floor(eps / 10 m)
# is 0, so the bar lands exactly on the dune line, which is where the road is.
# This keeps road_rows() as the single implementation of the geometry rather
# than reimplementing it here where it could drift.
_ZERO_SETBACK_EPS_M = 1e-6


def _drawable_setbacks(setbacks_m, managed):
    """Setbacks with roadless domains blanked and on-dune roads kept visible.

    Args:
        setbacks_m: Padded per-domain setbacks in metres, one per domain.
        managed: Per-domain boolean mask of the domains carrying a road.

    Returns:
        A float array: 0.0 where there is no road, so road_rows() blanks it,
        and at least _ZERO_SETBACK_EPS_M where there is one so a road whose
        setback has decayed to zero still draws, on the dune line.
    """
    values = np.asarray(setbacks_m, dtype=float)
    drawable = np.zeros_like(values)
    on = np.asarray(managed, dtype=bool)
    drawable[on] = np.maximum(values[on], _ZERO_SETBACK_EPS_M)
    return drawable


def _setback_at(roadway, year):
    """One domain's road setback in metres at a given year, or 0 if roadless."""
    series = getattr(roadway, "_road_setback_TS", None)
    if series is None or len(series) == 0:
        return 0.0
    value = series[min(year, len(series) - 1)]
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", help="matrix run directory")
    parser.add_argument("--fps", type=float, default=3.0)
    parser.add_argument("--out", default=None)
    parser.add_argument("--start-year", type=int, default=None,
                        help="calendar year of frame 0, for the title")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    (grids_by_year, offsets_by_year,
     setbacks_by_year, detected) = load_history(run_dir)
    n_years = len(grids_by_year)

    # The run name carries the period, so the calendar year is recoverable
    # without trusting an attribute that may not exist.
    start_year = args.start_year or detected
    if start_year is None:
        for token in run_dir.name.split("_"):
            if token.isdigit() and len(token) == 4:
                start_year = int(token)
                break

    canvases = [build_canvas(grids, offsets, GEOMETRY)
                for grids, offsets in zip(grids_by_year, offsets_by_year)]
    # ONE y limit for every frame: what moves should be the island.
    frame_rows = max(canvas[0].shape[0] for canvas in canvases)

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10, "axes.linewidth": 0.8, "axes.edgecolor": "#3A4149",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.color": "#3A4149", "ytick.color": "#3A4149",
        "legend.frameon": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    figure = plt.figure(figsize=(16, 6.2))
    axis = figure.add_axes([0.062, 0.135, 0.885, 0.66])
    out = Path(args.out) if args.out else (
        run_dir / f"{run_dir.name}_planview_evolution.gif")

    # Static furniture, drawn once. The run name is provenance, not a title,
    # so it goes in the subtitle at reading weight rather than in bold above
    # the panel where it competed with the map.
    figure.text(0.062, 0.945, "Barrier evolution in plan view",
                fontsize=15, fontweight="bold", color="#15202C",
                ha="left", va="center")
    figure.text(0.062, 0.905,
                f"{run_dir.name}   ·   {GEOMETRY.num_real_domains} domains "
                f"·   10 m cells   ·   elevation relative to MHW",
                fontsize=9.5, color="#5C6874", ha="left", va="center")

    def frame(index):
        year_index = min(index, n_years - 1)
        axis.clear()
        canvas, starts, per_domain, first_real = canvases[year_index]
        plot_canvas(canvas, starts, per_domain, first_real, GEOMETRY,
                    title="", ax=axis, colorbar=False,
                    xlabel=f"Alongshore domain   (south {ARROW} north,  "
                           f"Cape Hatteras {ARROW} Rodanthe)")

        # THE ROAD. Placed by the same road_rows() the static figure uses --
        # offset + floor(setback / cell) -- so the two agree by construction
        # rather than by a second implementation that could drift.
        overlay_roadway(axis, setbacks_by_year[year_index],
                        offsets_by_year[year_index], GEOMETRY,
                        starts, per_domain, first_real)

        axis.set_ylim(0, frame_rows)
        # Cells are an implementation detail; kilometres are what a reader
        # measures the island in.
        ticks = np.arange(0, frame_rows + 1, 50)
        axis.set_yticks(ticks)
        axis.set_yticklabels([f"{t * DEFAULT_PLAN_VIEW.cell_size_m / 1000:g}"
                              for t in ticks])
        axis.set_ylabel("Cross-shore distance (km)")

        label = (f"{start_year + year_index}" if start_year
                 else f"year {year_index}")
        axis.annotate(label, xy=(0.992, 0.94), xycoords="axes fraction",
                      ha="right", va="top", fontsize=26, color="#FFFFFF",
                      fontweight="bold", alpha=0.85, zorder=12)
        axis.annotate(f"year {year_index} of {n_years - 1}",
                      xy=(0.992, 0.80), xycoords="axes fraction",
                      ha="right", va="top", fontsize=9.5, color="#FFFFFF",
                      alpha=0.8, zorder=12)

    # The colorbar is drawn once, outside the frame loop: plot_canvas would
    # otherwise add a new one on every frame and shrink the axes each time.
    frame(0)
    mesh = axis.collections[0]
    cax = figure.add_axes([0.955, 0.135, 0.011, 0.66])
    bar = figure.colorbar(mesh, cax=cax)
    bar.set_label("Elevation (m MHW)", fontsize=9.5)
    bar.set_ticks([-1, 0, 1, 2, 3, 4])
    bar.ax.tick_params(labelsize=8.5)

    # One legend for the road, at figure level so a frame redraw cannot drop
    # it and it never lands on the barrier.
    road_handles = overlay_roadway(axis, setbacks_by_year[0], offsets_by_year[0],
                                   GEOMETRY, *canvases[0][1:])
    if road_handles:
        figure.legend(handles=road_handles, loc="upper right",
                      bbox_to_anchor=(0.947, 0.962), fontsize=9,
                      handlelength=1.8, ncol=len(road_handles))

    animation = FuncAnimation(figure, frame,
                              frames=n_years + HOLD_FRAMES, blit=False)
    animation.save(out, writer=PillowWriter(fps=args.fps))
    plt.close(figure)
    print(f"wrote {out}")
    print(f"  {n_years} model years + {HOLD_FRAMES} held, "
          f"{GEOMETRY.num_real_domains} real domains, frame {frame_rows} rows")


if __name__ == "__main__":
    main()
