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

ROAD PLACEMENT IS VERIFIED AGAINST THE MODEL, NOT ASSERTED
    `bulldoze` indexes the roadway as `road_start = int(road_setback / 10)`
    rows into the interior grid and flattens those rows to `road_ele`. So the
    road it draws is checkable: the bulldozed rows are exactly constant
    alongshore in `DomainTS`, and in the 1984-2004 calibBE groin run
    `floor(setback / 10 m)` lands on that constant row in every domain-year
    tested. The overlay uses the same `road_rows()` the static figure uses, so
    the map, the static figure and the model all index the road identically.

    What the map does NOT show is the dune line the setback is measured from:
    `DomainTS` is the INTERIOR only, and Barrier3D keeps `DuneDomain`
    separately, seaward of interior row 0. A setback of 0 m therefore draws on
    the interior's seaward edge, which is one dune width (20 m) landward of the
    dune crest.

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
from matplotlib.lines import Line2D

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

HOLD_FRAMES = 8      # frames held on the final year, so the totals can be read
ARROW = chr(0x2192)  # a real arrow, not '->'
DOT = chr(0x00B7)
TIMES = chr(0x00D7)
# The relocated colour from road_planview's own style, so a relocation marker
# here and a relocated road bar in the static figure read as the same thing.
RELOCATION_COLOR = "#FF8C00"
# A road the model has stopped managing is still a road on the ground, so it is
# drawn -- greyed, at its last managed position, and never counted again.
ABANDONED_COLOR = "#8A8F94"
HUD_COLOR = "#15202C"
SCALE_BAR_KM = 5.0


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


class RunHistory(object):
    """Everything one run contributes to the animation, already per-year.

    Attributes:
        grids: grids[t] is a padded-order list of metre-valued interior arrays
            on the full cross-shore frame.
        offsets: offsets[t] is the matching per-domain canvas row origin.
        setbacks: setbacks[t] is the per-domain drawable setback in metres,
            zeroed where no road should be drawn that year.
        managed: managed[t] is the per-domain mask of roads the model is still
            managing that year.
        relocations: relocations[t] is the per-domain relocation flag.
        rebuilds: rebuilds[t] is the per-domain dune-rebuild flag.
        prescribed_setback_m: per-domain relocation setback the run was given.
            Where this is 0 a "relocation" cannot move the road landward at
            all, so the event is a rebuild in place, not a retreat.
        road_domains: per-domain mask of domains carrying a road at all.
        start_year: calendar year of frame 0, or None.
    """

    def __init__(self, **fields):
        self.__dict__.update(fields)


def load_history(run_dir):
    """Per-year grids, shoreline offsets and roadway state from a run's .npz.

    Args:
        run_dir: A matrix run directory holding exactly one .npz.

    Returns:
        A RunHistory.

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
    road_domains = np.asarray(
        getattr(cascade, "_roadway_management_module", []), dtype=bool)
    if road_domains.size != len(inner):
        road_domains = np.zeros(len(inner), dtype=bool)
        for gis in range(HATTERAS_FIRST_ROAD_DOMAIN,
                         HATTERAS_LAST_ROAD_DOMAIN + 1):
            pad = GEOMETRY.gis_to_pad(gis)
            if 0 <= pad < road_domains.size:
                road_domains[pad] = True

    prescribed = np.zeros(len(inner), dtype=float)
    given = np.asarray(getattr(cascade, "_road_setback", []), dtype=float)
    if given.size == prescribed.size:
        prescribed = given

    # WHEN THE MODEL STOPS MANAGING A ROAD it returns from RoadwayManager.update
    # BEFORE writing that year's time series, so _road_setback_TS stays 0 for
    # every remaining year. Read literally that draws an abandoned road pinned
    # to the dune line for the rest of the run -- exactly where a drowned road
    # is not. _road_ele_TS is the reliable per-year flag: the manager stops the
    # moment the road elevation would go below 0 m MHW, and leaves 0 behind.
    managed_by_year = []
    for year in range(n_years):
        managed_by_year.append(np.array([
            _managed_at(roadways[pad] if pad < len(roadways) else None, year)
            for pad in range(len(inner))]) & road_domains)

    setbacks_by_year, relocations_by_year, rebuilds_by_year = [], [], []
    grids_by_year, offsets_by_year = [], []
    last_setback = np.zeros(len(inner), dtype=float)

    for year in range(n_years):
        grids_by_year.append([
            pad_cross_shore(np.asarray(model.DomainTS[year], dtype=float)
                            * config.dam_to_m, config)
            for model in inner])
        moved = np.array([float(model.x_s_TS[year]) for model in inner])
        offsets_by_year.append(base_offset + (moved - x_s_initial))

        raw = np.array([_setback_at(road, year) for road in roadways]
                       if roadways else [0.0] * len(inner), dtype=float)
        alive = managed_by_year[year]
        # An abandoned road keeps its last managed setback so it can be drawn
        # greyed where it actually sits, rather than snapping to the dune line.
        last_setback = np.where(alive, raw, last_setback)
        setbacks_by_year.append(_drawable_setbacks(last_setback, road_domains))

        # RELOCATIONS ARE EVENTS, not a state: _road_relocated_TS is a 0/1 flag
        # per domain per year, raised in the year the roadway manager moves the
        # road. Kept per year rather than accumulated so a frame shows what
        # happened THAT year and the running total can be built from it.
        relocations_by_year.append(np.array([
            _relocated_at(road, year) for road in roadways]
            if roadways else [False] * len(inner)) & alive)
        rebuilds_by_year.append(np.array([
            _rebuilt_at(road, year) for road in roadways]
            if roadways else [False] * len(inner)) & alive)

    start_year = int(getattr(cascade, "_start_year", 0)) or None
    return RunHistory(
        grids=grids_by_year, offsets=offsets_by_year,
        setbacks=setbacks_by_year, managed=managed_by_year,
        relocations=relocations_by_year, rebuilds=rebuilds_by_year,
        prescribed_setback_m=prescribed, road_domains=road_domains,
        start_year=start_year)


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


def _drawable_setbacks(setbacks_m, has_road):
    """Setbacks with roadless domains blanked and on-dune roads kept visible.

    Args:
        setbacks_m: Padded per-domain setbacks in metres, one per domain.
        has_road: Per-domain boolean mask of the domains carrying a road.

    Returns:
        A float array: 0.0 where there is no road, so road_rows() blanks it,
        and at least _ZERO_SETBACK_EPS_M where there is one so a road whose
        setback has decayed to zero still draws, on the dune line.
    """
    values = np.asarray(setbacks_m, dtype=float)
    drawable = np.zeros_like(values)
    on = np.asarray(has_road, dtype=bool)
    drawable[on] = np.maximum(values[on], _ZERO_SETBACK_EPS_M)
    return drawable


def _series_at(roadway, name, year, default=0.0):
    """One entry of a roadway time series, or a default when it is absent."""
    series = getattr(roadway, name, None)
    if series is None or len(series) == 0 or year >= len(series):
        return default
    try:
        return float(series[year])
    except (TypeError, ValueError):
        return default


def _managed_at(roadway, year):
    """Whether the model was still managing this road in this model year.

    The roadway manager writes a positive road elevation every year it runs and
    returns without writing once it gives the road up, so a zero elevation is
    the abandonment flag. Year 0 is the initial state, before any update.
    """
    if roadway is None:
        return False
    if year == 0:
        return True
    return _series_at(roadway, "_road_ele_TS", year) > 0.0


def _relocated_at(roadway, year):
    """Whether this domain's road was relocated in this model year."""
    return bool(_series_at(roadway, "_road_relocated_TS", year))


def _rebuilt_at(roadway, year):
    """Whether this domain's dune was rebuilt in this model year."""
    return bool(_series_at(roadway, "_dunes_rebuilt_TS", year))


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


def _road_row(offset, setback_m):
    """The canvas row a setback puts the road on, as road_rows() does it."""
    return offset + np.floor(max(float(setback_m), 0.0)
                             / DEFAULT_PLAN_VIEW.cell_size_m)


def build_relocation_tally(history):
    """Splits relocation events by whether they win the road any clearance.

    A relocation resets the setback to the domain's PRESCRIBED relocation
    setback, and CASCADE has no separate parameter for that: cascade_groin.py
    re-assigns `road_relocation_setback = road_setback`, the domain's STARTING
    setback, every year. So a domain that starts with the road on the dune line
    can only ever be relocated back onto the dune line.

    That is the case at GIS 85 and 86, whose measured 1984 setback floored to
    0 m. Each event still drags the road one 10 m cell landward -- it rides the
    dune toe -- but it ends the year with the same zero clearance it began
    with, so the next cell of retreat re-fires it. The counts are exact: GIS 85
    retreats 72.7 m (7.3 cells) and relocates 7 times, GIS 86 retreats 59.9 m
    (6.0 cells) and relocates 6 times. One relocation per cell.

    A domain with clearance behaves completely differently. GIS 9 holds 40 m,
    absorbs 38 m of retreat over the whole run and never triggers at all.

    So the split is NOT moved against not-moved -- every event moves the road.
    It is relocated-with-clearance against re-pinned-at-the-dune-line, and
    reporting one total for both would say the road retreated eighteen times
    when it bought itself room five times.

    Args:
        history: A RunHistory.

    Returns:
        (moves_by_year, pinned_by_year, cumulative_moves, cumulative_pinned)
        where the first two are per-domain boolean arrays per year and the last
        two are integer running totals per year.
    """
    can_move = np.asarray(history.prescribed_setback_m, dtype=float) > 0.0
    moves = [flags & can_move for flags in history.relocations]
    pinned = [flags & ~can_move for flags in history.relocations]
    return (moves, pinned,
            np.cumsum([int(f.sum()) for f in moves]),
            np.cumsum([int(f.sum()) for f in pinned]))


def draw_timeline(axis, moves_by_year, pinned_by_year, start_year):
    """Draws the static relocation timeline strip below the map.

    A plan view has no time axis, so animating one leaves the reader with no
    sense of WHEN anything happened -- only of what is on screen now. The strip
    carries the whole run at once and the cursor says where in it this frame
    sits, which is what turns a loop into a record.

    Args:
        axis: Axes for the strip.
        moves_by_year: Per-year per-domain masks of relocations that move.
        pinned_by_year: Per-year per-domain masks of 0 m rebuilds in place.
        start_year: Calendar year of frame 0, or None.

    Returns:
        The cursor Line2D, to be moved each frame.
    """
    years = np.arange(len(moves_by_year))
    per_move = np.array([int(f.sum()) for f in moves_by_year])
    per_pinned = np.array([int(f.sum()) for f in pinned_by_year])
    labels = years + start_year if start_year else years
    step = max(1, len(years) // 12)

    axis.bar(years, per_move, width=0.72, color=RELOCATION_COLOR,
             label="relocated, gains clearance", zorder=3)
    axis.bar(years, per_pinned, width=0.72, bottom=per_move, color="white",
             edgecolor=RELOCATION_COLOR, linewidth=0.9, hatch="///",
             label="re-pinned at dune line (0 m setback)", zorder=3)

    ceiling = max(1, int((per_move + per_pinned).max()))
    axis.set_ylim(0, ceiling + 0.6)
    axis.set_xlim(-0.8, len(years) - 0.2)
    axis.set_yticks(range(0, ceiling + 1, max(1, ceiling // 3)))
    axis.set_xticks(years[::step])
    axis.set_xticklabels([str(int(v)) for v in labels[::step]])
    axis.set_ylabel("relocation\nevents", fontsize=8.5, color="#5C6874")
    axis.tick_params(labelsize=8, colors="#5C6874")
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.spines["left"].set_color("#C6CCD2")
    axis.spines["bottom"].set_color("#C6CCD2")
    axis.grid(axis="y", color="#E6EAEE", lw=0.6, zorder=0)

    cursor = axis.axvline(0, color=HUD_COLOR, lw=1.4, zorder=5)
    # No legend here: the map's legend already names both bar colours, and a
    # second copy under the panel only competes with the year axis.
    axis.set_title("Relocation events per model year", fontsize=9,
                   color="#5C6874", loc="left", pad=4)
    return cursor


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", help="matrix run directory")
    parser.add_argument("--fps", type=float, default=3.0)
    parser.add_argument("--out", default=None)
    parser.add_argument("--start-year", type=int, default=None,
                        help="calendar year of frame 0, for the title")
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    history = load_history(run_dir)
    n_years = len(history.grids)

    # The run name carries the period, so the calendar year is recoverable
    # without trusting an attribute that may not exist.
    start_year = args.start_year or history.start_year
    if start_year is None:
        for token in run_dir.name.split("_"):
            if token.isdigit() and len(token) == 4:
                start_year = int(token)
                break

    moves, pinned, cum_moves, cum_pinned = build_relocation_tally(history)
    total_moves = int(cum_moves[-1]) if n_years else 0
    total_pinned = int(cum_pinned[-1]) if n_years else 0
    total_events = total_moves + total_pinned
    n_road_domains = int(history.road_domains.sum())
    total_rebuilds = int(sum(int(f.sum()) for f in history.rebuilds))
    n_moved_domains = int(np.any(np.array(history.relocations), axis=0).sum())
    abandoned = any(bool((history.road_domains & ~alive).any())
                    for alive in history.managed)

    # Domains that have relocated at least once BY a given year, so the map can
    # keep showing where the road has already had to move rather than flashing
    # it for a single frame and losing it.
    ever_moved = np.zeros_like(history.road_domains)
    ever_by_year = []
    for flags in history.relocations:
        ever_moved = ever_moved | flags
        ever_by_year.append(ever_moved.copy())

    canvases = [build_canvas(grids, offsets, GEOMETRY)
                for grids, offsets in zip(history.grids, history.offsets)]
    # ONE y limit for every frame: what moves should be the island.
    frame_rows = max(canvas[0].shape[0] for canvas in canvases)
    total_cols = canvases[0][0].shape[1]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10, "axes.linewidth": 0.8, "axes.edgecolor": "#3A4149",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.color": "#3A4149", "ytick.color": "#3A4149",
        "legend.frameon": False,
        "figure.facecolor": "white", "savefig.facecolor": "white",
    })
    figure = plt.figure(figsize=(16, 7.6))
    map_box = [0.062, 0.345, 0.885, 0.535]
    axis = figure.add_axes(map_box)
    strip = figure.add_axes([0.062, 0.095, 0.885, 0.115])
    out = Path(args.out) if args.out else (
        run_dir / f"{run_dir.name}_planview_evolution.gif")

    # THE CROSS-SHORE AXIS IS STRETCHED. Cells are square, so the distortion is
    # just the pixel aspect of the axes box; stating it means the reader can
    # take the shapes at face value instead of guessing at them.
    exaggeration = ((total_cols / (map_box[2] * figure.get_figwidth()))
                    / (frame_rows / (map_box[3] * figure.get_figheight())))

    # Static furniture, drawn once. The run name is provenance, not a title,
    # so it goes in the subtitle at reading weight rather than in bold above
    # the panel where it competed with the map.
    figure.text(0.062, 0.960, "Barrier evolution in plan view",
                fontsize=15, fontweight="bold", color=HUD_COLOR,
                ha="left", va="center")
    figure.text(0.062, 0.929,
                f"{run_dir.name}   {DOT}   {GEOMETRY.num_real_domains} domains "
                f"{DOT}   10 m cells   {DOT}   elevation relative to MHW   "
                f"{DOT}   cross-shore axis stretched {exaggeration:.1f}{TIMES}",
                fontsize=9.5, color="#5C6874", ha="left", va="center")
    figure.text(0.062, 0.906,
                "interior domain only: the dune line a road setback is "
                "measured from lies one dune width (20 m) seaward of the "
                "drawn island edge, and is not shown",
                fontsize=8.5, color="#8A939C", ha="left", va="center",
                style="italic")

    cursor = draw_timeline(strip, moves, pinned, start_year)

    def frame(index):
        year_index = min(index, n_years - 1)
        axis.clear()
        canvas, starts, per_domain, first_real = canvases[year_index]
        plot_canvas(canvas, starts, per_domain, first_real, GEOMETRY,
                    title="", ax=axis, colorbar=False,
                    xlabel=f"Alongshore domain   (south {ARROW} north,  "
                           f"Cape Hatteras {ARROW} Rodanthe)")

        offsets = history.offsets[year_index]
        setbacks = history.setbacks[year_index]
        alive = history.managed[year_index]

        # THE ROAD. Placed by the same road_rows() the static figure uses --
        # offset + floor(setback / cell) -- so the two agree by construction
        # rather than by a second implementation that could drift. A road the
        # model has given up on is drawn separately, greyed, at the last
        # position it was managed in, so abandonment reads as abandonment
        # rather than as a road that has snapped onto the shoreline.
        overlay_roadway(axis, np.where(alive, setbacks, 0.0), offsets,
                        GEOMETRY, starts, per_domain, first_real)
        for pad in np.flatnonzero(history.road_domains & ~alive):
            plotted = int(pad) - GEOMETRY.start_real_index
            if not 0 <= plotted < len(starts):
                continue
            axis.hlines(_road_row(offsets[pad], setbacks[pad]),
                        starts[plotted], starts[plotted] + per_domain[plotted],
                        color=ABANDONED_COLOR, lw=2.0, linestyle=(0, (3, 2)),
                        zorder=6)

        # RELOCATIONS, marked in the year they happen and kept afterwards as a
        # faint tick. A relocation is the one thing in this animation the model
        # DECIDES rather than suffers, so it is drawn as an event marker above
        # the road rather than as a change of road colour, which would be
        # invisible on a 1-cell bar. An event that cannot move the road -- a
        # prescribed relocation setback of 0 m -- gets a hollow marker, because
        # counting those as retreat is the easiest way to overstate this figure.
        for pad in np.flatnonzero(ever_by_year[year_index]):
            plotted = int(pad) - GEOMETRY.start_real_index
            if not 0 <= plotted < len(starts):
                continue
            centre = starts[plotted] + per_domain[plotted] / 2
            row = _road_row(offsets[pad], setbacks[pad])
            axis.vlines(centre, row + 4, row + 13, color=RELOCATION_COLOR,
                        lw=1.3, alpha=0.75, zorder=10)
            if moves[year_index][pad] or pinned[year_index][pad]:
                filled = bool(moves[year_index][pad])
                axis.plot(centre, row + 15, marker="v", markersize=9,
                          color=RELOCATION_COLOR if filled else "white",
                          markeredgecolor=("white" if filled
                                           else RELOCATION_COLOR),
                          markeredgewidth=0.8 if filled else 1.3,
                          zorder=11, clip_on=False)

        # THE COUNTER, in the empty ocean at lower left where it competes with
        # nothing. Running and final totals both, because "how often did the
        # road have to move" is the question this animation exists to answer,
        # and a per-year flash never answers it.
        this_year = int(moves[year_index].sum() + pinned[year_index].sum())
        running = int(cum_moves[year_index] + cum_pinned[year_index])
        lines = [
            ("NC-12 relocation events", 10.0, "bold", HUD_COLOR, 0.062),
            (f"{running} of {total_events}", 22.0, "bold", RELOCATION_COLOR,
             0.105),
            (f"{int(cum_moves[year_index])} of {total_moves} gain clearance "
             f"from the dune line", 9.5, "normal", HUD_COLOR, 0.058),
            (f"{int(cum_pinned[year_index])} of {total_pinned} re-pinned at "
             f"the dune line, 0 m setback", 9.5, "normal", "#5C6874", 0.058),
            (f"{this_year} this year   {DOT}   {n_moved_domains} of "
             f"{n_road_domains} road domains ever affected",
             8.5, "normal", "#8A939C", 0.0),
        ]
        y = 0.40
        for text, size, weight, colour, drop in lines:
            axis.annotate(text, xy=(0.012, y), xycoords="axes fraction",
                          ha="left", va="top", fontsize=size,
                          fontweight=weight, color=colour, zorder=12)
            y -= drop

        # A SCALE BAR, because the two axes are at different scales and a
        # reader measuring the island off the tick labels alone will get the
        # alongshore distance wrong.
        bar_cells = SCALE_BAR_KM * 1000.0 / DEFAULT_PLAN_VIEW.cell_size_m
        bar_x0 = total_cols * 0.012
        bar_y = frame_rows * 0.045
        axis.hlines(bar_y, bar_x0, bar_x0 + bar_cells, color=HUD_COLOR, lw=2.4,
                    zorder=12)
        axis.annotate(f"{SCALE_BAR_KM:g} km alongshore",
                      xy=(bar_x0 + bar_cells / 2, bar_y + frame_rows * 0.012),
                      ha="center", va="bottom", fontsize=8.5, color=HUD_COLOR,
                      zorder=12)

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
        cursor.set_xdata([year_index, year_index])

    # The colorbar is drawn once, outside the frame loop: plot_canvas would
    # otherwise add a new one on every frame and shrink the axes each time.
    frame(0)
    mesh = axis.collections[0]
    cax = figure.add_axes([0.955, map_box[1], 0.011, map_box[3]])
    bar = figure.colorbar(mesh, cax=cax)
    bar.set_label("Elevation (m MHW)", fontsize=9.5)
    bar.set_ticks([-1, 0, 1, 2, 3, 4])
    bar.ax.tick_params(labelsize=8.5)

    # One legend for the road, at figure level so a frame redraw cannot drop
    # it and it never lands on the barrier.
    road_handles = list(overlay_roadway(
        axis, history.setbacks[0], history.offsets[0], GEOMETRY,
        *canvases[0][1:]))
    if total_moves:
        road_handles.append(Line2D(
            [], [], linestyle="none", marker="v", markersize=9,
            color=RELOCATION_COLOR, markeredgecolor="white",
            label="relocated, gains clearance"))
    if total_pinned:
        road_handles.append(Line2D(
            [], [], linestyle="none", marker="v", markersize=9,
            color="white", markeredgecolor=RELOCATION_COLOR,
            label="re-pinned at dune line (0 m setback)"))
    if total_events:
        road_handles.append(Line2D(
            [], [], color=RELOCATION_COLOR, lw=1.3, alpha=0.75,
            label="relocated earlier in the run"))
    if abandoned:
        road_handles.append(Line2D(
            [], [], color=ABANDONED_COLOR, lw=2.0, linestyle=(0, (3, 2)),
            label="road no longer managed"))
    if road_handles:
        figure.legend(handles=road_handles, loc="upper right",
                      bbox_to_anchor=(0.947, 0.980), fontsize=9,
                      handlelength=1.8, ncol=len(road_handles))

    animation = FuncAnimation(figure, frame,
                              frames=n_years + HOLD_FRAMES, blit=False)
    animation.save(out, writer=PillowWriter(fps=args.fps))
    plt.close(figure)
    print(f"wrote {out}")
    print(f"  {n_years} model years + {HOLD_FRAMES} held, "
          f"{GEOMETRY.num_real_domains} real domains, frame {frame_rows} rows")
    print(f"  {total_events} relocation events on {n_road_domains} road "
          f"domains ({n_moved_domains} ever affected): {total_moves} move the "
          f"road landward, {total_pinned} re-pin it at the dune line")
    print(f"  {total_rebuilds} dune rebuilds; cross-shore axis stretched "
          f"{exaggeration:.2f}x")


if __name__ == "__main__":
    main()
