"""Side-by-side animation of NC-12 against a migrating dune line.

WHAT IS DRAWN
    One frame per model year, two panels sharing a year clock and a y-axis:

        left   relocations OFF -- roadway_manager decides on its own
        right  relocations ON  -- the measured historical displacements

    In each panel, per alongshore domain:

        dune line   the modelled shoreline, as displacement from year 0
        road        that line PLUS the domain's current setback
        marker      a domain that relocated in this frame's year

WHY THE ROAD IS DRAWN AS "DUNE LINE PLUS SETBACK"
    `road_setback` is the road's distance LANDWARD of the interior domain's
    seaward edge, and `roadway_manager` decrements it by dune migration every
    year precisely so the road stays geographically put while the dune line
    advances on it. On a landward-positive axis, `dune + setback` reproduces
    that: a road that is not moving traces a FLAT line while the dune line
    climbs toward it, and a relocation shows up as the road stepping up in a
    single frame. A road drawn at an absolute position would hide the very
    mechanism the animation exists to show.

    The y-axis is therefore cross-shore displacement relative to each domain's
    own year-0 dune line, not an absolute cross-shore coordinate. Two domains
    at the same height on the plot are NOT at the same place on the island.

SIGN CONVENTION
    x_s_TS increases landward, and `run.flip_sign_model` turns that into a
    seaward-positive quantity -- which shoreline_gif then draws on an inverted
    axis. This module negates once more instead, so the axis is plainly
    landward-positive and ascending: 0 is the year-0 dune line and larger is
    further landward. Do not pre-flip the matrix before passing it in.

WHERE A LINE STOPS
    A road that drowns stops being drawn. The managed span is dated from
    `_road_ele_TS`, never from the setback, because a setback of exactly
    0.0 m is legitimate -- see `_last_managed`.

Author: Hannah A. Henry, UNC CECL
"""

import io
import os

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from cascade_pipeline.annotations import DEFAULT_ANNOTATIONS
from cascade_pipeline.domains import DEFAULT_DOMAINS
from cascade_pipeline.plotting.shoreline_gif import (
    DEFAULT_GIF_CONFIG,
    _open_file,
    _slug,
)

COLOR_DUNE = "#1f6fb4"
COLOR_ROAD = "#b03030"
COLOR_RELOC = "#f0a202"
COLOR_BAY = "#4a8fbf"
# Deliberately not the star colour: a prescribed move and a module-triggered
# one are different claims and must not share a glyph.
COLOR_PRESCRIBED = "#00B0F0"
COLOR_LAND = "#efe0bd"


def _panel_series(shoreline_m, run, pad_lo, pad_hi):
    """Year-0-referenced dune-line displacement for one arm, LANDWARD-POSITIVE.

    shoreline_gif keeps a seaward-positive axis and then inverts the y-limits
    to put the ocean at the bottom. That works for a single shoreline but reads
    badly here, where the road sits tens of metres landward of the dune line
    and would be drawn at large NEGATIVE numbers. This negates once more, so
    landward is positive and up, the dune line starts at 0 and climbs as it
    retreats, and the road is simply `dune + setback`.

    Args:
        shoreline_m: 2-D [n_years, total_domains] raw x_s_TS matrix, metres.
        run: RunInfo carrying flip_sign_model.
        pad_lo: First padded index to keep.
        pad_hi: One past the last padded index to keep.

    Returns:
        A 2-D array [n_years, pad_hi - pad_lo]; row 0 is identically zero and
        positive means landward of the year-0 dune line.
    """
    raw = np.asarray(shoreline_m, dtype=float)
    flip = -1.0 if run.flip_sign_model else 1.0
    pos = -(flip * raw)[:, pad_lo:pad_hi]
    return pos - pos[0][None, :]


def _referenced(other_m, shoreline_m, run, pad_lo, pad_hi):
    """Puts a second cross-shore series on the dune line's own axis.

    The back-barrier shoreline has to be referenced to the SAME year-0 dune
    line as `_panel_series`, not to its own year-0 position; otherwise both
    curves start at zero and the island is drawn with no width at all.

    Args:
        other_m: 2-D raw x_b_TS matrix in metres, same shape as shoreline_m.
        shoreline_m: The x_s_TS matrix supplying the reference.
        run: RunInfo carrying flip_sign_model.
        pad_lo: First padded index to keep.
        pad_hi: One past the last padded index to keep.

    Returns:
        A 2-D landward-positive array on the dune line's axis.
    """
    flip = -1.0 if run.flip_sign_model else 1.0
    other = -(flip * np.asarray(other_m, dtype=float))[:, pad_lo:pad_hi]
    ref = -(flip * np.asarray(shoreline_m, dtype=float))[0, pad_lo:pad_hi]
    return other - ref[None, :]


def _last_managed(entry):
    """Index of the last year the manager ran, dated from road ELEVATION.

    Not from the setback series: a setback of exactly 0.0 m is legitimate --
    it means the road sits on the dune line, which is where six of the ten
    historical domains start -- so a nonzero test on the setback would draw no
    road at all in precisely the domains this animation exists to show. Road
    elevation has no such ambiguity; the module stops the moment it drops
    below 0 m MHW.

    Args:
        entry: One road_series() value, carrying "elevation".

    Returns:
        The index, or -1 if the manager never ran.
    """
    written = np.flatnonzero(np.asarray(entry["elevation"], dtype=float) != 0)
    return int(written[-1]) if written.size else -1


def _road_matrix(series, road_series, gis_lo, gis_hi, domains, n_years):
    """Builds road position and relocation-flag grids aligned to `series`.

    Args:
        series: Dune-line displacement [n_years, n_domains] for this arm.
        road_series: {gis: {"setback", "relocated", ...}} for this arm.
        gis_lo: First GIS domain in the window.
        gis_hi: Last GIS domain in the window.
        domains: DomainGeometry.
        n_years: Number of animation years.

    Returns:
        A (road, relocated) tuple of [n_years, n_domains] arrays. Domains the
        run did not manage are NaN in `road` and False in `relocated`, so an
        unmanaged domain draws no road rather than a road at the dune line.
    """
    n_dom = gis_hi - gis_lo + 1
    road = np.full((n_years, n_dom), np.nan)
    relocated = np.zeros((n_years, n_dom), dtype=bool)

    for col, gis in enumerate(range(gis_lo, gis_hi + 1)):
        entry = road_series.get(gis)
        if entry is None:
            continue
        setback = np.asarray(entry["setback"], dtype=float)
        reloc = np.asarray(entry["relocated"], dtype=float)
        # The road stops being drawn where the record stops -- a drowned road
        # must leave a gap, not a line frozen at its last position.
        stop = min(n_years, len(setback), _last_managed(entry) + 1)
        road[:stop, col] = series[:stop, col] + setback[:stop]
        m = min(n_years, len(reloc))
        relocated[:m, col] = reloc[:m] > 0
    return road, relocated


def make_road_relocation_gif(
    arm_a, arm_b, road_series_a, road_series_b,
    gis_lo, gis_hi, out_path, back_a=None, back_b=None,
    # A STAR IS ALWAYS THE MODULE, IN BOTH PANELS. It is driven by
    # `_road_relocated_TS`, the RoadwayManager's own counter, and a PRESCRIBED
    # historical move never increments it -- the pipeline applies those as a
    # displacement before the manager updates. So arm B stars in the years its
    # module fired, not in 1989/1999, and the prescribed move shows only as a
    # step in the road. That was silently misleading, hence the separate
    # prescribed marker and these labels.
    label_a="relocations OFF  —  module decides   "
            "[★ = module triggered a move this year]",
    label_b="relocations ON  —  measured moves applied   "
            "[★ = module ALSO fired this year]",
    event_years=None,
    domains=DEFAULT_DOMAINS,
    annotations=DEFAULT_ANNOTATIONS,
    gif_config=DEFAULT_GIF_CONFIG,
    fps=None, stride=None, title=None,
):
    """Writes the two-panel road/dune animation for one alongshore window.

    Args:
        arm_a: (shoreline_m, RunInfo) for the relocations-OFF run.
        arm_b: (shoreline_m, RunInfo) for the relocations-ON run.
        road_series_a: {gis: series dict} for arm A, from the comparison
            script's road_series().
        road_series_b: The same for arm B.
        gis_lo: First GIS domain to draw.
        gis_hi: Last GIS domain to draw.
        out_path: Destination .gif path.
        back_a: Optional raw x_b_TS matrix for arm A, metres. Supplying it
            draws the barrier body between the two shorelines instead of a
            single dune line, which is what makes the plot read as an
            island rather than as a chart.
        back_b: The same for arm B.
        label_a: Panel title for arm A.
        label_b: Panel title for arm B.
        event_years: {gis: year} marking the historical relocation year, drawn
            as a reference tick on the right panel.
        domains: DomainGeometry.
        annotations: Geographic annotation layer (town spans, groins).
        gif_config: Shared GifConfig.
        fps: Frames per second; defaults to gif_config.
        stride: Year stride; defaults to gif_config.
        title: Figure suptitle prefix.

    Returns:
        The written path, or None if the animation was skipped.
    """
    try:
        from PIL import Image
    except Exception as exc:
        print(f"  [GIF] Pillow unavailable ({exc}); skipping animation.")
        return None

    fps = gif_config.fps if fps is None else fps
    stride = gif_config.year_stride if stride is None else stride
    event_years = event_years or {}

    shore_a, run_a = arm_a
    shore_b, run_b = arm_b
    pad_lo = domains.gis_to_pad(gis_lo)
    pad_hi = domains.gis_to_pad(gis_hi) + 1

    series_a = _panel_series(shore_a, run_a, pad_lo, pad_hi)
    series_b = _panel_series(shore_b, run_b, pad_lo, pad_hi)
    bay_a = (None if back_a is None
             else _referenced(back_a, shore_a, run_a, pad_lo, pad_hi))
    bay_b = (None if back_b is None
             else _referenced(back_b, shore_b, run_b, pad_lo, pad_hi))
    n_years = min(series_a.shape[0], series_b.shape[0])
    if n_years < 2:
        print("  [GIF] fewer than 2 model years; skipping animation.")
        return None
    series_a, series_b = series_a[:n_years], series_b[:n_years]
    if bay_a is not None:
        bay_a, bay_b = bay_a[:n_years], bay_b[:n_years]

    road_a, reloc_a = _road_matrix(series_a, road_series_a, gis_lo, gis_hi,
                                   domains, n_years)
    road_b, reloc_b = _road_matrix(series_b, road_series_b, gis_lo, gis_hi,
                                   domains, n_years)

    x = np.arange(gis_lo, gis_hi + 1)

    # One y-axis for both panels, fixed across all frames: a road that jumps
    # in one panel and not the other must be readable as a difference between
    # the arms, not as a difference between two autoscaled axes.
    parts = [series_a.ravel(), series_b.ravel(), road_a.ravel(), road_b.ravel()]
    if bay_a is not None:
        # Include the back-barrier or the island body is drawn clipped, which
        # reads as the sound being part of the barrier.
        parts += [bay_a.ravel(), bay_b.ravel()]
    stack = np.concatenate(parts)
    finite = stack[np.isfinite(stack)]
    ymin, ymax = float(np.min(finite)), float(np.max(finite))
    ypad = (ymax - ymin) * 0.10 or 1.0
    # Landward-positive (see _panel_series), so ocean-at-bottom is the plain
    # ascending axis rather than an inverted one.
    ylim = ((ymin - ypad, ymax + ypad) if gif_config.ocean_at_bottom
            else (ymax + ypad, ymin - ypad))

    n_dom = len(x)
    width = float(np.clip(5.0 + 0.16 * n_dom, 9.0, 18.0))
    year_idx = list(range(0, n_years, max(int(stride), 1)))
    if year_idx[-1] != n_years - 1:
        year_idx.append(n_years - 1)

    frames = []
    for t in year_idx:
        year = run_a.start_year + t
        fig, axes = plt.subplots(1, 2, figsize=(width, 5.2), dpi=110,
                                 sharey=True)
        # Fixed margins (NOT bbox_inches="tight") so every frame is the same
        # size -- mismatched frame dimensions break GIF assembly. Left margin
        # is wide enough for the y-label at every window width.
        fig.subplots_adjust(left=0.105, right=0.985, top=0.855, bottom=0.225,
                            wspace=0.06)
        fig.patch.set_facecolor("white")

        # The last flag marks the arm that CARRIES the prescribed moves.
        panels = ((axes[0], series_a, road_a, reloc_a, label_a, bay_a, False),
                  (axes[1], series_b, road_b, reloc_b, label_b, bay_b, True))
        for ax, series, road, reloc, label, bay, is_prescribed_arm in panels:
            ax.set_facecolor("white")
            ax.set_ylim(*ylim)
            ax.set_xlim(gis_lo - 0.5, gis_hi + 0.5)
            ax.grid(alpha=0.25, lw=0.6)

            for span_label, (d_lo, d_hi) in annotations.town_spans.items():
                if d_hi < gis_lo or d_lo > gis_hi:
                    continue
                ax.axvspan(max(d_lo, gis_lo) - 0.5, min(d_hi, gis_hi) + 0.5,
                           color=annotations.color_town_span, alpha=0.13,
                           zorder=0)

            if bay is not None:
                # The island itself: everything between the ocean dune line
                # and the back-barrier shoreline.
                ax.fill_between(x, series[t], bay[t], color=COLOR_LAND,
                                alpha=0.55, zorder=1, lw=0)
                ax.plot(x, bay[t], color=COLOR_BAY, lw=1.1, zorder=3,
                        label="back-barrier")
            ax.plot(x, series[t], color=COLOR_DUNE, lw=1.9,
                    label="dune line", zorder=3)
            ax.plot(x, road[t], color=COLOR_ROAD, lw=1.7, zorder=4,
                    label="NC-12")
            # The gap between the two lines IS the setback; shading it makes
            # "the dune line is closing on the road" the thing you see.
            ax.fill_between(x, series[t], road[t], where=np.isfinite(road[t]),
                            color=COLOR_ROAD, alpha=0.10, zorder=1)

            fired = np.flatnonzero(reloc[t])
            if fired.size:
                ax.plot(x[fired], road[t][fired], marker="*", ls="none",
                        ms=15, mfc=COLOR_RELOC, mec="k", mew=0.6, zorder=6,
                        label="relocated this year")

            # The PRESCRIBED move, marked only in the arm that carries it and
            # only in its event year. Without this the bottom panel's headline
            # event was invisible: the road simply steps, with no marker at
            # all, while the stars nearby are the module doing something else.
            if is_prescribed_arm:
                for gis, ev_year in event_years.items():
                    if gis_lo <= gis <= gis_hi and year == ev_year:
                        col = gis - gis_lo
                        if np.isfinite(road[t][col]):
                            ax.plot([x[col]], [road[t][col]], marker="D",
                                    ls="none", ms=9, mfc=COLOR_PRESCRIBED,
                                    mec="k", mew=0.8, zorder=7)

            for gis, ev_year in event_years.items():
                if gis_lo <= gis <= gis_hi:
                    ax.axvline(gis, color="0.45", ls=":", lw=0.9, zorder=2,
                               alpha=0.9 if year >= ev_year else 0.35)

            ax.set_title(label, fontsize=10)
            ax.set_xlabel("alongshore domain (GIS)")

        axes[0].set_ylabel("cross-shore displacement since "
                           f"{run_a.start_year} (m)\nlandward ▲")
        handles = [
            Line2D([], [], color=COLOR_DUNE, lw=1.9, label="ocean dune line"),
            Line2D([], [], color=COLOR_BAY, lw=1.1, label="back-barrier"),
            Line2D([], [], color=COLOR_ROAD, lw=1.7, label="NC-12"),
            Line2D([], [], color="none", marker="*", ms=13,
                   mfc=COLOR_RELOC, mec="k", mew=0.6,
                   label="road moved this year (see panel title "
                         "for which kind)"),
            Line2D([], [], color="none", marker="D", ms=8,
                   mfc=COLOR_PRESCRIBED, mec="k", mew=0.8,
                   label="measured 1989/1999 move applied (lower panel only)"),
            Line2D([], [], color="0.45", ls=":", lw=0.9,
                   label="domain that relocated historically "
                         "(bright once its year has passed)"),
        ]
        # Figure-level and below the axes: an in-axes legend sits on top of
        # the road wherever the setback is large, which is most of the island.
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=8,
                   frameon=False, bbox_to_anchor=(0.5, 0.005))

        head = title or f"NC-12 vs the dune line, GIS {gis_lo}-{gis_hi}"
        fig.suptitle(f"{head}   |   {year}", fontsize=12, fontweight="bold")

        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    frames[0].save(out_path, save_all=True, append_images=frames[1:],
                   duration=int(1000 / max(fps, 1)), loop=0)
    print(f"  [GIF] saved ({len(frames)} frames, GIS {gis_lo}-{gis_hi}): "
          f"{out_path}")
    if gif_config.auto_open:
        _open_file(out_path)
    return out_path


def make_all_road_gifs(arm_a, arm_b, road_series_a, road_series_b,
                       windows, out_dir, event_years=None,
                       back_a=None, back_b=None,
                       gif_config=DEFAULT_GIF_CONFIG, **kwargs):
    """Writes one animation per alongshore window.

    Args:
        arm_a: (shoreline_m, RunInfo) for the relocations-OFF run.
        arm_b: (shoreline_m, RunInfo) for the relocations-ON run.
        road_series_a: {gis: series dict} for arm A.
        road_series_b: {gis: series dict} for arm B.
        windows: Sequence of (name, gis_lo, gis_hi) tuples.
        out_dir: Directory for the .gif files.
        event_years: {gis: year} for the historical relocation reference ticks.
        gif_config: Shared GifConfig.
        **kwargs: Forwarded to make_road_relocation_gif.

    Returns:
        A list of written paths, skipped windows omitted.
    """
    written = []
    for name, gis_lo, gis_hi in windows:
        path = os.path.join(out_dir, f"road_relocation_{_slug(name)}.gif")
        result = make_road_relocation_gif(
            arm_a, arm_b, road_series_a, road_series_b, gis_lo, gis_hi, path,
            back_a=back_a, back_b=back_b,
            event_years=event_years, gif_config=gif_config,
            title=f"NC-12 vs the dune line — {name}", **kwargs)
        if result:
            written.append(result)
    return written


# =============================================================================
# Topographic plan view
# =============================================================================
# The line version above abstracts the island to two curves. This one paints
# Barrier3D's actual interior grids, so the animation shows overwash fans,
# the barrier narrowing, and the road sitting on real ground.
#
# GEOMETRY. Every domain carries `DomainTS[t]`, a (cross-shore rows, 50)
# elevation grid in dam MHW whose row 0 is the seaward edge of the interior.
# The domain's absolute cross-shore position is `x_s_TS[t]` (dam), so row r
# lands at `x_s_TS[t] + r`. Alongshore, each domain is 50 cells of 10 m, i.e.
# 500 m, and the domains abut -- so 90 real domains tile 45 km of island.
#
# WHAT IS WATER. The interior grid extends well past the subaerial barrier
# (176 rows where the barrier is only ~32), so most of it is sound. Cells at
# or below 0 m MHW are masked and drawn as water rather than as low land.

WATER_COLOR = "#c8dcea"

# Cross-shore window drawn around the window's reference shoreline, in dam.
# The seaward end is fixed (a few cells of ocean for context); the landward
# end is chosen per window from how far the land actually reaches, because a
# fixed value either clips the wide Tri-Village section or drowns narrow Pea
# Island in an empty frame.
CROSS_SHORE_SEAWARD_DAM = -4
CROSS_SHORE_MIN_ROWS = 40
CROSS_SHORE_HEADROOM_ROWS = 8
CROSS_SHORE_REACH_PERCENTILE = 75


def _window_reference(cascades, gis_lo, gis_hi, domains, years):
    """Cross-shore origin shared by every panel, domain and year.

    Taken as the most SEAWARD shoreline anywhere in the window, over both arms
    and every drawn year, minus a few cells of open ocean. An earlier version
    used the window's mean year-0 shoreline, which put domains seaward of that
    mean at a negative row: `_island_raster` clipped their land to row 0 while
    `_road_track` drew the road at its true negative position, so NC-12
    appeared to float in the sea off Pea Island. A minimum cannot do that.

    Args:
        cascades: Iterable of finished Cascade instances sharing the window.
        gis_lo: First GIS domain in the window.
        gis_hi: Last GIS domain in the window.
        domains: DomainGeometry.
        years: Iterable of year indices that will be drawn.

    Returns:
        The cross-shore dam coordinate of raster row 0.
    """
    pads = [domains.gis_to_pad(g) for g in range(gis_lo, gis_hi + 1)]
    seaward = min(
        float(cas.barrier3d[pad].x_s_TS[min(t, len(cas.barrier3d[pad].x_s_TS) - 1)])
        for cas in cascades for pad in pads for t in years)
    return seaward + CROSS_SHORE_SEAWARD_DAM


def _cross_shore_rows(cascade, gis_lo, gis_hi, domains, years, x_ref):
    """Chooses a landward extent that fits the land in this window.

    Scans every drawn year and takes the furthest landward row that is still
    above 0 m MHW anywhere in the window, then adds headroom. Fixed once for
    the whole animation so the frame does not breathe year to year.

    Args:
        cascade: A finished Cascade instance.
        gis_lo: First GIS domain in the window.
        gis_hi: Last GIS domain in the window.
        domains: DomainGeometry.
        years: Iterable of year indices that will be drawn.
        x_ref: Cross-shore origin from _window_reference.

    Returns:
        The number of cross-shore rows to draw.
    """
    b3d = cascade.barrier3d
    pads = [domains.gis_to_pad(g) for g in range(gis_lo, gis_hi + 1)]
    reach = []
    for pad in pads:
        bb = b3d[pad]
        widest = 0
        for t in years:
            t = min(t, len(bb.DomainTS) - 1)
            interior = np.asarray(bb.DomainTS[t], dtype=float)
            land = np.flatnonzero((interior > 0).any(axis=1))
            if land.size:
                r0 = int(round(float(bb.x_s_TS[t]) - x_ref))
                widest = max(widest, r0 + int(land[-1]))
        reach.append(widest)

    # A PERCENTILE, not the maximum. One wide domain in the window -- a
    # Tri-Village section beside narrow Pea Island, say -- otherwise sets a
    # frame tall enough to make the domains under test unreadable. Clipping
    # the widest domain costs nothing here: the subject is the road corridor,
    # which sits within a few hundred metres of the ocean shoreline.
    tall = float(np.percentile(reach, CROSS_SHORE_REACH_PERCENTILE)) if reach else 0.0
    return int(max(tall, CROSS_SHORE_MIN_ROWS) + CROSS_SHORE_HEADROOM_ROWS)


def _island_raster(cascade, year, gis_lo, gis_hi, domains, n_cross, x_ref):
    """Assembles one year's elevation raster for an alongshore window.

    Args:
        cascade: A finished Cascade instance.
        year: Index into the per-domain time series (0 = first model year).
        gis_lo: First GIS domain in the window.
        gis_hi: Last GIS domain in the window.
        domains: DomainGeometry.
        n_cross: Number of cross-shore rows to draw, from
            _cross_shore_rows.
        x_ref: Cross-shore origin from _window_reference.

    Returns:
        A [n_cross, n_domains * 50] array in metres MHW, NaN where the model
        has no cell.
    """
    b3d = cascade.barrier3d
    pads = [domains.gis_to_pad(g) for g in range(gis_lo, gis_hi + 1)]

    # One reference shoreline for the whole window, so the panel does not
    # shift under the island as the shoreline retreats.
    cell = int(np.shape(b3d[pads[0]].DomainTS[0])[1])
    grid = np.full((n_cross, len(pads) * cell), np.nan)

    for col, pad in enumerate(pads):
        bb = b3d[pad]
        t = min(year, len(bb.DomainTS) - 1)
        interior = np.asarray(bb.DomainTS[t], dtype=float)     # (rows, cell) dam
        r0 = int(round(float(bb.x_s_TS[t]) - x_ref))

        lo = max(r0, 0)
        hi = min(r0 + interior.shape[0], n_cross)
        if hi > lo:
            grid[lo:hi, col * cell:(col + 1) * cell] = (
                interior[lo - r0:hi - r0, :] * 10.0)

        # Dune rows sit seaward of interior row 0, height above the berm.
        dune = np.asarray(bb.DuneDomain[t], dtype=float)        # (cell, width)
        berm = float(bb.BermEl)
        for w in range(dune.shape[1]):
            r = r0 - dune.shape[1] + w
            if 0 <= r < n_cross:
                grid[r, col * cell:(col + 1) * cell] = (dune[:, w] + berm) * 10.0
    return grid


def _road_track(cascade, road_series, year, gis_lo, gis_hi, domains, x_ref):
    """Road cross-shore position for one year, in raster row units.

    Args:
        cascade: A finished Cascade instance.
        road_series: {gis: series dict} for this arm.
        year: Index into the time series.
        gis_lo: First GIS domain in the window.
        gis_hi: Last GIS domain in the window.
        domains: DomainGeometry.
        x_ref: Cross-shore dam coordinate of raster row 0.

    Returns:
        An (x_cells, rows, relocated_mask) tuple. `rows` is NaN wherever the
        domain is unmanaged or its record has stopped, so the road line breaks
        rather than being drawn across a drowned stretch.
    """
    b3d = cascade.barrier3d
    cell = int(np.shape(b3d[domains.gis_to_pad(gis_lo)].DomainTS[0])[1])
    n = (gis_hi - gis_lo + 1)
    xs, rows, flags = [], [], []
    for col, gis in enumerate(range(gis_lo, gis_hi + 1)):
        centre = col * cell + cell / 2.0
        xs.append(centre)
        entry = road_series.get(gis)
        pad = domains.gis_to_pad(gis)
        if entry is None or year > _last_managed(entry):
            rows.append(np.nan)
            flags.append(False)
            continue
        setback = np.asarray(entry["setback"], dtype=float)
        t = min(year, len(setback) - 1)
        # setback is metres landward of the interior's seaward edge; the
        # raster is in 10 m rows, so /10 puts it in the same units.
        rows.append(float(b3d[pad].x_s_TS[t]) - x_ref + setback[t] / 10.0)
        reloc = np.asarray(entry["relocated"], dtype=float)
        flags.append(bool(t < len(reloc) and reloc[t] > 0))
    return np.array(xs), np.array(rows), np.array(flags)


def _land_colormap():
    """matplotlib's `terrain`, restricted to its LAND range.

    Full `terrain` spends its first quarter on blues and cyans for bathymetry.
    Everything drawn here is at or above 0 m MHW -- water is masked out and
    painted separately -- so those blues would render low-lying land in the
    same colour as the sound beside it, which is exactly the distinction the
    animation exists to show. Starting at 0.25 gives the familiar green-to-
    brown-to-white land ramp with nothing wasted below sea level.

    Returns:
        A Colormap whose "bad" colour is the water fill.
    """
    base = plt.get_cmap("terrain")
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "terrain_land", base(np.linspace(0.25, 1.0, 256)))
    cmap.set_bad(WATER_COLOR)
    return cmap


def make_topography_gif(cascade_a, cascade_b, road_series_a, road_series_b,
                        gis_lo, gis_hi, out_path, start_year,
                        label_a="relocations OFF  —  module decides   "
                                "[★ = module triggered a move this year]",
                        label_b="relocations ON  —  measured moves applied   "
                                "[★ = module ALSO fired this year]",
                        event_years=None, vmax_m=3.0,
                        domains=DEFAULT_DOMAINS,
                        gif_config=DEFAULT_GIF_CONFIG,
                        fps=None, stride=None, title=None,
                        planform_note=None):
    """Animated elevation map of the island with NC-12 drawn on it.

    Two panels sharing a colour scale and a year clock. Cells at or below
    0 m MHW are drawn as water, so the barrier's real outline, its overwash
    fans and its narrowing are all visible rather than implied.

    VERTICAL EXAGGERATION. The window is ~640 m across-shore and 500 m per
    domain alongshore, so the aspect is set from the data and stated on the
    figure. Nothing here is to scale in the way a map is; it is a raster of
    model cells.

    Args:
        cascade_a: Finished Cascade for the relocations-OFF run.
        cascade_b: Finished Cascade for the relocations-ON run.
        road_series_a: {gis: series dict} for arm A.
        road_series_b: {gis: series dict} for arm B.
        gis_lo: First GIS domain to draw.
        gis_hi: Last GIS domain to draw.
        out_path: Destination .gif path.
        start_year: Calendar year of model year 0.
        label_a: Panel title for arm A.
        label_b: Panel title for arm B.
        event_years: {gis: year} for the historical relocation reference ticks.
        vmax_m: Top of the elevation colour scale, metres MHW. 3 m rather
            than the barrier's true maximum: almost every cell is under 2 m,
            so a scale topped by the few dune crests washes the island out.
        domains: DomainGeometry.
        gif_config: Shared GifConfig.
        fps: Frames per second; defaults to gif_config.
        stride: Year stride; defaults to gif_config.
        title: Figure suptitle prefix.
        planform_note: Optional caption, e.g. recording that the run's
            shoreline offset does not carry the island's true curvature.

    Returns:
        The written path, or None if the animation was skipped.
    """
    try:
        from PIL import Image
    except Exception as exc:
        print(f"  [GIF] Pillow unavailable ({exc}); skipping animation.")
        return None

    fps = gif_config.fps if fps is None else fps
    stride = gif_config.year_stride if stride is None else stride
    event_years = event_years or {}

    n_years = min(len(cascade_a.barrier3d[0].DomainTS),
                  len(cascade_b.barrier3d[0].DomainTS))
    year_idx = list(range(0, n_years, max(int(stride), 1)))
    if year_idx[-1] != n_years - 1:
        year_idx.append(n_years - 1)

    terrain = _land_colormap()
    # Fixed for the whole animation and shared by both panels, so a colour or
    # a height means the same thing in every frame and in both arms.
    x_ref = _window_reference((cascade_a, cascade_b), gis_lo, gis_hi, domains,
                              year_idx)
    n_cross = max(
        _cross_shore_rows(cascade_a, gis_lo, gis_hi, domains, year_idx, x_ref),
        _cross_shore_rows(cascade_b, gis_lo, gis_hi, domains, year_idx, x_ref))

    n_dom = gis_hi - gis_lo + 1
    width = float(np.clip(5.0 + 0.42 * n_dom, 9.0, 18.0))

    frames = []
    for t in year_idx:
        year = start_year + t
        fig, axes = plt.subplots(2, 1, figsize=(width, 7.4), dpi=110,
                                 sharex=True)
        # bottom raised from 0.105 to clear the new legend strip and the
        # planform note beneath it
        fig.subplots_adjust(left=0.085, right=0.90, top=0.90, bottom=0.155,
                            hspace=0.22)
        fig.patch.set_facecolor("white")

        # Last flag marks the arm carrying the prescribed moves; see
        # COLOR_PRESCRIBED and the label note on make_road_relocation_gif.
        for ax, cas, series, label, is_prescribed_arm in (
                (axes[0], cascade_a, road_series_a, label_a, False),
                (axes[1], cascade_b, road_series_b, label_b, True)):
            grid = _island_raster(cas, t, gis_lo, gis_hi, domains,
                                  n_cross, x_ref)
            masked = np.ma.masked_invalid(np.ma.masked_less_equal(grid, 0.0))
            im = ax.imshow(masked, origin="lower", aspect="auto",
                           cmap=terrain, vmin=0.0, vmax=vmax_m,
                           interpolation="nearest")
            ax.set_facecolor(WATER_COLOR)

            x, rows, fired = _road_track(cas, series, t, gis_lo, gis_hi,
                                         domains, x_ref)
            ax.plot(x, rows, color=COLOR_ROAD, lw=2.0, zorder=5, label="NC-12")
            if fired.any():
                ax.plot(x[fired], rows[fired], marker="*", ls="none", ms=15,
                        mfc=COLOR_RELOC, mec="k", mew=0.6, zorder=6)
            if is_prescribed_arm:
                for gis, ev_year in event_years.items():
                    if gis_lo <= gis <= gis_hi and year == ev_year:
                        col = gis - gis_lo
                        if col < len(rows) and np.isfinite(rows[col]):
                            ax.plot([x[col]], [rows[col]], marker="D",
                                    ls="none", ms=9, mfc=COLOR_PRESCRIBED,
                                    mec="k", mew=0.8, zorder=7)

            cell = grid.shape[1] / n_dom
            ticks = np.arange(n_dom) * cell + cell / 2.0
            ax.set_xticks(ticks[::max(1, n_dom // 12)])
            ax.set_xticklabels(range(gis_lo, gis_hi + 1)[::max(1, n_dom // 12)])
            for gis in event_years:
                if gis_lo <= gis <= gis_hi:
                    ax.axvline((gis - gis_lo) * cell + cell / 2.0, color="0.25",
                               ls=":", lw=0.9, zorder=4, alpha=0.75)
            # Cells are 10 m, so the axis is labelled in metres: "row 23"
            # is not a distance anyone can check against a map.
            step = max(10, int(round(n_cross / 5.0 / 10.0)) * 10)
            rows_at = np.arange(0, n_cross, step)
            ax.set_yticks(rows_at)
            ax.set_yticklabels((rows_at * 10).astype(int))
            ax.set_ylabel("cross-shore (m)\nlandward →", fontsize=9)
            ax.text(0.006, 0.045, "ocean", transform=ax.transAxes,
                    fontsize=7.5, color="0.30", style="italic")
            ax.text(0.006, 0.93, "sound", transform=ax.transAxes,
                    fontsize=7.5, color="0.30", style="italic")
            ax.set_title(label, fontsize=10, pad=4)

        axes[1].set_xlabel("alongshore domain (GIS)")
        cax = fig.add_axes([0.915, 0.155, 0.016, 0.745])
        fig.colorbar(im, cax=cax, label="elevation (m MHW)")

        # This panel had NO legend at all: the star and the dotted lines were
        # drawn unexplained, and a reader had no way to tell whether a star
        # meant "the module decided to move the road" or "history did". The
        # panel titles now carry the per-panel meaning; these entries name the
        # glyphs themselves.
        fig.legend(handles=[
            Line2D([], [], color=COLOR_ROAD, lw=2.0, label="NC-12"),
            Line2D([], [], color="none", marker="*", ms=13,
                   mfc=COLOR_RELOC, mec="k", mew=0.6,
                   label="MODULE triggered a move this year (either panel)"),
            Line2D([], [], color="none", marker="D", ms=8,
                   mfc=COLOR_PRESCRIBED, mec="k", mew=0.8,
                   label="MEASURED 1989/1999 move applied (lower panel)"),
            Line2D([], [], color="0.25", ls=":", lw=0.9,
                   label="domain that relocated historically"),
        ], loc="lower center", ncol=3, fontsize=8, frameon=False,
            bbox_to_anchor=(0.5, 0.038))

        head = title or f"Hatteras topography and NC-12, GIS {gis_lo}-{gis_hi}"
        fig.suptitle(f"{head}   |   {year}", fontsize=12, fontweight="bold")
        if planform_note:
            fig.text(0.085, 0.006, planform_note, fontsize=7.5, color="0.35")

        buf = io.BytesIO()
        fig.savefig(buf, format="png", facecolor=fig.get_facecolor())
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE))

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    frames[0].save(out_path, save_all=True, append_images=frames[1:],
                   duration=int(1000 / max(fps, 1)), loop=0)
    print(f"  [GIF] saved topography ({len(frames)} frames, "
          f"GIS {gis_lo}-{gis_hi}): {out_path}")
    return out_path
