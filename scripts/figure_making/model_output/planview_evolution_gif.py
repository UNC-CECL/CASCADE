#!/usr/bin/env python3
"""
Animated plan view of the island through a run: elevation, the road, and every relocation.

    python scripts/figure_making/model_output/planview_evolution_gif.py <run_dir> [--out PATH]

Reads one run's .npz (grids, shoreline offsets, roadway state per year) and
its island offset; writes an animated GIF with a relocation timeline strip.
Details: scripts/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# House style, typeface only (site_layer/hat_figure_style.py)
import sys as _sys
from pathlib import Path as _HP
_sys.path.insert(0, str(next(_q for _q in _HP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import apply_style  # noqa: E402
apply_style()
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D

_HERE = Path(__file__).resolve()
# Project root found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / 'pyproject.toml').exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts",):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_DOMAINS as GEOMETRY, HATTERAS_PERIODS,
    HATTERAS_FIRST_ROAD_DOMAIN, HATTERAS_LAST_ROAD_DOMAIN)

# Defined as the hindcast runner defines it; the period table's paths are relative to it
HATTERAS_DATA_BASE = PROJECT_BASE_DIR / "data" / "hatteras_init"
from cascade_pipeline.hindcast import load_island_offset_dam  # noqa: E402
from cascade_pipeline.plotting.init_planview import (  # noqa: E402
    DEFAULT_PLAN_VIEW, build_canvas, pad_cross_shore, plot_canvas)
from cascade_pipeline.plotting.road_planview import overlay_roadway  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HOLD_FRAMES = 8      # frames held on the final year, so the totals can be read
ARROW = chr(0x2192)  # a real arrow, not '->'
DOT = chr(0x00B7)
TIMES = chr(0x00D7)
# The relocation colour from road_planview, so both figures read alike
RELOCATION_COLOR = "#FF8C00"
# An abandoned road stays drawn, greyed at its last managed position
ABANDONED_COLOR = "#8A8F94"
HUD_COLOR = "#15202C"
SCALE_BAR_KM = 5.0
# -----------------------------------------------------------------------------


# The hindcast period a run belongs to, from its name first
def _period_start(run_dir, cascade):
    for token in Path(run_dir).name.split("_"):
        if token.isdigit() and int(token) in HATTERAS_PERIODS:
            return int(token)
    attribute = getattr(cascade, "_start_year", None)
    if attribute is not None and int(attribute) in HATTERAS_PERIODS:
        return int(attribute)
    raise SystemExit(
        f"cannot tell which hindcast period {Path(run_dir).name} belongs to; "
        f"known starts are {sorted(HATTERAS_PERIODS)}")


# Everything one run contributes to the animation, per year
class RunHistory(object):

    def __init__(self, **fields):
        self.__dict__.update(fields)


# Per-year grids, shoreline offsets and roadway state from a run's .npz
def load_history(run_dir):
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

    # The canvas frame is the island offset; x_s supplies only its change (README)
    offset_file = HATTERAS_DATA_BASE / HATTERAS_PERIODS[
        _period_start(run_dir, cascade)]["island_offset_file"]
    base_offset = load_island_offset_dam(offset_file, GEOMETRY)
    x_s_initial = np.array([float(model.x_s_TS[0]) for model in inner])

    # The road moves too: its setback is measured from a retreating dune line
    roadways = getattr(cascade, "_roadways", None) or []

    # Where the road is comes from this run's own roadway mask, not the whole reach (README)
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

    # _road_ele_TS flags management each year; the setback series stops when a road is abandoned
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
        # An abandoned road keeps its last managed setback
        last_setback = np.where(alive, raw, last_setback)
        setbacks_by_year.append(_drawable_setbacks(last_setback, road_domains))

        # Relocations are events: kept per year, totals built from them
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


# A zero setback inside the road reach is the road on the dune line, not no road (README)
_ZERO_SETBACK_EPS_M = 1e-6


# Setbacks with roadless domains blanked and on-dune roads kept visible
def _drawable_setbacks(setbacks_m, has_road):
    values = np.asarray(setbacks_m, dtype=float)
    drawable = np.zeros_like(values)
    on = np.asarray(has_road, dtype=bool)
    drawable[on] = np.maximum(values[on], _ZERO_SETBACK_EPS_M)
    return drawable


# One entry of a roadway time series, or a default when absent
def _series_at(roadway, name, year, default=0.0):
    series = getattr(roadway, name, None)
    if series is None or len(series) == 0 or year >= len(series):
        return default
    try:
        return float(series[year])
    except (TypeError, ValueError):
        return default


# Was the model still managing this road in this year?
def _managed_at(roadway, year):
    if roadway is None:
        return False
    if year == 0:
        return True
    return _series_at(roadway, "_road_ele_TS", year) > 0.0


# Was this domain's road relocated in this year?
def _relocated_at(roadway, year):
    return bool(_series_at(roadway, "_road_relocated_TS", year))


# Was this domain's dune rebuilt in this year?
def _rebuilt_at(roadway, year):
    return bool(_series_at(roadway, "_dunes_rebuilt_TS", year))


# One domain's road setback (m) in a year, 0 if roadless
def _setback_at(roadway, year):
    series = getattr(roadway, "_road_setback_TS", None)
    if series is None or len(series) == 0:
        return 0.0
    value = series[min(year, len(series) - 1)]
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


# The canvas row a setback puts the road on, as road_rows() does
def _road_row(offset, setback_m):
    return offset + np.floor(max(float(setback_m), 0.0)
                             / DEFAULT_PLAN_VIEW.cell_size_m)


# Split relocations by whether they won the road any clearance
def build_relocation_tally(history):
    can_move = np.asarray(history.prescribed_setback_m, dtype=float) > 0.0
    moves = [flags & can_move for flags in history.relocations]
    pinned = [flags & ~can_move for flags in history.relocations]
    return (moves, pinned,
            np.cumsum([int(f.sum()) for f in moves]),
            np.cumsum([int(f.sum()) for f in pinned]))


# The static relocation timeline below the map
def draw_timeline(axis, moves_by_year, pinned_by_year, start_year):
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
    # No legend here
    axis.set_title("Relocation events per model year", fontsize=9,
                   color="#5C6874", loc="left", pad=4)
    return cursor


# Run: load the run, draw every frame, write the GIF
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

    # The calendar year comes from the run name's period
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

    # Domains relocated at least once by a given year, so earlier moves stay visible
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

    # The cross-shore axis is stretched
    exaggeration = ((total_cols / (map_box[2] * figure.get_figwidth()))
                    / (frame_rows / (map_box[3] * figure.get_figheight())))

    # Static furniture, drawn once
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
                           f"Cape Point {ARROW} Rodanthe)")

        offsets = history.offsets[year_index]
        setbacks = history.setbacks[year_index]
        alive = history.managed[year_index]

        # The road
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

        # Relocations: marked in their year, then kept as a faint tick
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

        # The counter, in the empty ocean at lower left
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

        # A scale bar: the two axes are at different scales
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
        # Cells are an implementation detail
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

    # The colorbar is drawn once, outside the frame loop
    frame(0)
    mesh = axis.collections[0]
    cax = figure.add_axes([0.955, map_box[1], 0.011, map_box[3]])
    bar = figure.colorbar(mesh, cax=cax)
    bar.set_label("Elevation (m MHW)", fontsize=9.5)
    bar.set_ticks([-1, 0, 1, 2, 3, 4])
    bar.ax.tick_params(labelsize=8.5)

    # One legend for the road, at figure level
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
