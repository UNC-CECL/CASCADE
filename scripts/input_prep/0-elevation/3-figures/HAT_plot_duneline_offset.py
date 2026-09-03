"""
HAT_plot_duneline_offset.py

The 1984-start DEM with both digitized dune lines and the domain boxes on it,
and the cross-shore distance between the two lines, per domain.

WHAT THIS IS AND IS NOT
-----------------------
This measures ONE thing: how far apart the 1984 and 1997 dune lines are, in
metres, in each of the 90 domains. It says nothing about whether the 1996 ALACE
swath reaches either of them - that is a separate and much harder question, and
it is not asked here.

Nothing is modified. No elevation is written, no product is forked; the
`2009-2014-1996` rasters are read in place.

THE FRAME
---------
Every domain box is axis-aligned, 2000 m in easting by 500 m in northing, and
this pipeline's OCEAN_LOC is "right" - the Atlantic is at increasing easting.
So easting is cross-shore, northing is alongshore, and the separation between
two roughly shore-parallel lines is just a difference in easting at a shared
northing. That is checked at load, not assumed: the script raises if the boxes
are not 2000 x 500.

    offset_m = x_1984 - x_1997     at the same northing

    POSITIVE means the 1984 line lies SEAWARD of the 1997 line, which is the
    sign 13 years of erosion predicts.

Sampled every SAMPLE_SPACING_M along the northing axis of each box, so a domain
contributes up to 500 independent measurements and the per-domain number is
their median, with the quartiles beside it.

A SECOND, ORIENTATION-FREE DISTANCE
-----------------------------------
`nearest_m` is the plain nearest-point distance from each 1984 sample to the
1997 line - no axis, no sign, no assumption about which way the ocean is. It is
reported next to the easting difference as a check on the frame. Where the
island runs obliquely to the grid the two must diverge, because a cross-shore
difference measured along easting is the true separation divided by the cosine
of that obliquity. Large `offset_over_nearest` is not an error; it says the box
axis and the shoreline disagree there, and the number to quote is `nearest_m`.

WHY 1997 AND NOT 2004
---------------------
1997 is one year after the 1996 ALACE flight the DEM's beach comes from, so the
pair brackets the model's 1984 start and the DEM's own vintage. See
`data/hatteras_init/1-barrier3d-domains/raw-duneline-geojson/README.md` for
what each line is and the metadata caveat - 1997 carries `feature_type`,
`method` and `editor`; 1984 carries nothing at all, so "the same feature at
both ends" rests on the numbers rather than on the files.

INPUTS
    D:/Hatteras_GIS/domains.geojson
    data/.../0-elevation/2009-2014-1996/2-resampled-10m/resampled_domain_*.tif
    data/.../1-barrier3d-domains/raw-duneline-geojson/duneline_1984.geojson
    data/.../1-barrier3d-domains/raw-duneline-geojson/duneline_1997.geojson

OUTPUTS (data/hatteras_init/0-elevation/2009-2014-1996-duneline/)
    duneline_offset_by_domain.csv
    figures/HAT_duneline_offset_island.png    the DEM, both lines, the boxes
    figures/HAT_duneline_offset_bydomain.png  the offset, domain by domain

Requires: rasterio, geopandas, shapely, numpy, matplotlib

    python HAT_plot_duneline_offset.py
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from shapely.geometry import LineString, Point


def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from hat_elevation_products import ELEVATION_ROOT  # noqa: E402
# The island mosaic loader, the km axis formatter and the elevation panel are
# imported rather than re-written so this figure and the other 1984-start
# figures cannot drift apart in extent, colour or projection. Importing runs no
# IO beyond resolving the product path.
import HAT_plot_1984_mosaic as m  # noqa: E402


SOURCE_TAG = "2009-2014-1996"
OUT_DIR = ELEVATION_ROOT / f"{SOURCE_TAG}-duneline"
FIG_DIR = OUT_DIR / "figures"
CSV_NAME = "duneline_offset_by_domain.csv"

DUNE_DIR = INIT_ROOT / "1-barrier3d-domains" / "raw-duneline-geojson"
DUNE_LINES = {1984: DUNE_DIR / "duneline_1984.geojson",
              1997: DUNE_DIR / "duneline_1997.geojson"}

# One sample per metre of alongshore, matching the 1 m DEM the rest of the
# 1984-start chain is built on. 500 per domain.
SAMPLE_SPACING_M = 1.0
GRID_10M = 10.0          # the resampled product's cell, for the box fallback

# The row-insert scope table, read ONLY to label a --zoom with how many
# Barrier3D rows the measured offset turns into. Optional; absent is fine.
INSERT_SCOPE_CSV = (INIT_ROOT / "1-barrier3d-domains" / "1984-start"
                    / "row-insert-scope" / "row_insert_scope_by_domain.csv")

# The box shape the easting-is-cross-shore frame depends on. Checked, not
# assumed - see THE FRAME above.
EXPECTED_BOX_M = (2000.0, 500.0)
BOX_TOL_M = 1.0

# 1984 is the older line and the one the model starts from, so it gets the
# emphasis colour; 1997 is the reference. Deliberately NOT the road key from
# HAT_plot_1984_mosaic - these are dune lines, and reusing black/white-dashed
# would read as NC-12 on a figure where NC-12 is absent.
LINE_STYLE = {1984: dict(color="#d62728", linestyle="-", linewidth=1.6),
              1997: dict(color="#1f77b4", linestyle=(0, (4.0, 2.5)),
                         linewidth=1.6)}
LINE_CASING = {1984: dict(color="white", linewidth=3.2),
               1997: dict(color="white", linewidth=3.2)}
LINE_ORDER = [1997, 1984]      # 1984 drawn last, on top

BOX_STYLE = dict(edgecolor="0.25", facecolor="none", linewidth=0.4)
LABEL_EVERY = 5                # label every Nth domain box on the island map

# The island is 46 km long and ~2 km wide. Split into thirds, each panel gets
# its own extent and roughly three times the scale - see fig_island.
N_PANELS = 3
PANEL_PAD_M = 400.0

CELL_M = 10.0                  # the Barrier3D cell, drawn on the bar chart

# The ribbon's baseline: a boxcar over the two lines' mean, alongshore. Long
# enough to keep several domains of shared curve, short enough that the
# island's 6.5 km sweep does not survive it. See fig_ribbon.
BASELINE_WINDOW_M = 2000.0

# (first domain, last domain, short label, what the reach is). Chosen from the
# measured table, not by eye: 62-68 is the largest sustained NEGATIVE run
# (-29.7 to -58.9 m over seven neighbours), 78-85 the largest POSITIVE one (up
# to +70.2 m), and 17-21 is the quietest five-domain run on the island
# (|offset| <= 7.8 m). The control is not optional - without it every figure of
# this kind reads as a discrepancy, and there is no way to see what agreement
# looks like at the same scale.
ZOOM_REACHES = [
    (17, 21, "17-21", "the quietest reach on the island"),
    (62, 68, "62-68", "1984 line LANDWARD of 1997"),
    (78, 85, "78-85", "1984 line SEAWARD of 1997"),
]
# Cross-shore half-width of a zoom panel. The domain box is 2000 m across and
# almost all of it is water and back-barrier; cropping to the dune makes the
# offset a visible fraction of the frame at equal aspect.
ZOOM_HALF_WIDTH_M = 300.0


# =============================================================================
# THE MEASUREMENT
# =============================================================================

def x_at_northings(line_geom, box, ys):
    """
    Easting of a line at each of `ys`, inside one domain box.

    The line is clipped to the box first, then intersected with a horizontal
    segment at each northing. Where a sample crosses the line more than once -
    a hook, or a stretch running momentarily east-west - the MEAN easting is
    taken, the same choice made for the rasterized version of this measurement:
    there is no reason to prefer either crossing, and the mean is the position
    a reader means by "where the line is at this northing".

    Returns an array with NaN at northings the line does not cross. Those are
    NOT interpolated. A domain where the line genuinely runs outside the box
    should report fewer samples, not a fabricated position, and n_1984 / n_1997
    in the CSV are how that shows up.
    """
    minx, _, maxx, _ = box.bounds
    clipped = line_geom.intersection(box)
    out = np.full(ys.size, np.nan)
    if clipped.is_empty:
        return out
    for i, y in enumerate(ys):
        hit = clipped.intersection(LineString([(minx - 1.0, y),
                                               (maxx + 1.0, y)]))
        if hit.is_empty:
            continue
        xs = [p.x for p in getattr(hit, "geoms", [hit])
              if p.geom_type == "Point"]
        if xs:
            out[i] = float(np.mean(xs))
    return out


def measure(gdf, lines):
    """
    Per-domain offset between the two dune lines.

    Returns (rows, samples): one row per domain, and the raw 1 m alongshore
    samples the rows are medians of - northing, and each line's easting - kept
    so the ribbon figure can draw what the medians were computed from rather
    than an interpolation of them.
    """
    rows = []
    s_dom, s_y, s_x84, s_x97 = [], [], [], []
    for _, r in gdf.iterrows():
        box = r.geometry
        minx, miny, maxx, maxy = box.bounds
        if (abs((maxx - minx) - EXPECTED_BOX_M[0]) > BOX_TOL_M
                or abs((maxy - miny) - EXPECTED_BOX_M[1]) > BOX_TOL_M):
            raise SystemExit(
                f"domain {int(r['domain_id'])} box is "
                f"{maxx - minx:.0f} x {maxy - miny:.0f} m, expected "
                f"{EXPECTED_BOX_M[0]:.0f} x {EXPECTED_BOX_M[1]:.0f}.\n"
                f"  The easting-is-cross-shore frame this script measures in "
                f"does not hold for that box. See THE FRAME in the docstring.")

        ys = np.arange(miny + SAMPLE_SPACING_M / 2, maxy, SAMPLE_SPACING_M)
        x84 = x_at_northings(lines[1984], box, ys)
        x97 = x_at_northings(lines[1997], box, ys)

        off = x84 - x97                       # +ve = 1984 seaward of 1997
        both = np.isfinite(off)

        # Orientation-free check: nearest-point distance, 1984 sample to the
        # 1997 line. Unsigned by construction.
        # The 1997 line is used UNCLIPPED here: the nearest point to a 1984
        # sample near a box edge can legitimately lie in the neighbouring
        # domain, and clipping would inflate the distance there.
        near = np.full(ys.size, np.nan)
        for i in np.where(np.isfinite(x84))[0]:
            near[i] = lines[1997].distance(Point(x84[i], ys[i]))

        def q(a, p):
            a = a[np.isfinite(a)]
            return round(float(np.percentile(a, p)), 2) if a.size else ""

        med = q(off, 50)
        nmed = q(near, 50)
        rows.append({
            "domain": int(r["domain_id"]),
            "n_samples": int(ys.size),
            "n_1984": int(np.isfinite(x84).sum()),
            "n_1997": int(np.isfinite(x97).sum()),
            "n_both": int(both.sum()),
            "offset_med_m": med,
            "offset_p25_m": q(off, 25),
            "offset_p75_m": q(off, 75),
            "offset_min_m": q(off, 0),
            "offset_max_m": q(off, 100),
            "offset_mean_m": (round(float(np.nanmean(off)), 2)
                              if both.any() else ""),
            "offset_sd_m": (round(float(np.nanstd(off)), 2)
                            if both.any() else ""),
            "offset_med_cells": (round(med / CELL_M, 2)
                                 if med != "" else ""),
            "nearest_med_m": nmed,
            "offset_over_nearest": (round(abs(med) / nmed, 2)
                                    if med != "" and nmed not in ("", 0.0)
                                    else ""),
            "x1984_med": q(x84, 50),
            "x1997_med": q(x97, 50),
        })
        s_dom.append(np.full(ys.size, int(r["domain_id"])))
        s_y.append(ys)
        s_x84.append(x84)
        s_x97.append(x97)

    samples = {"domain": np.concatenate(s_dom), "y": np.concatenate(s_y),
               "x1984": np.concatenate(s_x84), "x1997": np.concatenate(s_x97)}
    order = np.argsort(samples["y"])
    return rows, {k: v[order] for k, v in samples.items()}


# =============================================================================
# FIGURES
# =============================================================================

def load_lines(dst_crs):
    """Both dune lines, reprojected. UNCLIPPED - this is what gets measured."""
    out = {}
    for yr, p in DUNE_LINES.items():
        g = gpd.read_file(p)
        src_crs = g.crs
        if g.crs is not None and dst_crs is not None:
            g = g.to_crs(dst_crs)
        # EPSG codes only. dst_crs here is the DEM's COMPOUND CRS and its
        # full WKT is ~1200 characters, which buries every other log line.
        print(f"  {yr}: {p.name}  {_epsg(src_crs)} -> {_epsg(dst_crs)}   "
              f"{g.geometry.length.sum() / 1000:.1f} km")
        out[yr] = g
    return out


def clip_for_drawing(lines, footprint):
    """
    The same clip HAT_plot_1984_mosaic.load_roads applies to NC-12, and for the
    same reason: both geojsons run past the 90 domains at the north end, and
    drawn unclipped they show dune line where there is no model domain, which
    reads as coverage that does not exist.

    DRAWING ONLY. The measurement keeps the unclipped lines, because a 1984
    sample near a box edge can legitimately have its nearest 1997 point just
    outside the footprint, and clipping would inflate `nearest_m` there.
    """
    out = {}
    for yr, g in lines.items():
        before = float(g.geometry.length.sum())
        c = gpd.clip(g, footprint)
        after = float(c.geometry.length.sum()) if len(c) else 0.0
        print(f"  {yr}: clipped to domains for drawing, "
              f"{after / 1000:.1f} km of {before / 1000:.1f} km kept")
        out[yr] = c
    return out


def draw_lines(ax, lines, scale=1.0):
    """Casing then line, in LINE_ORDER so 1984 lands on top of 1997."""
    for yr in LINE_ORDER:
        cas = dict(LINE_CASING[yr])
        cas["linewidth"] *= scale
        st = dict(LINE_STYLE[yr])
        st["linewidth"] *= scale
        lines[yr].plot(ax=ax, linestyle="-", alpha=0.9,
                       zorder=6 + LINE_ORDER.index(yr) * 2, **cas)
        lines[yr].plot(ax=ax, zorder=7 + LINE_ORDER.index(yr) * 2, **st)


def line_legend():
    return [Line2D([0], [0], label=f"{yr} dune line", **LINE_STYLE[yr])
            for yr in (1984, 1997)]


def fig_island(elev, extent, gdf, lines, rows):
    """
    The DEM, both dune lines, and the domain boxes.

    THREE PANELS, each a third of the island, rather than one frame. At equal
    aspect the island is 46 km long and about 2 km wide, so a single panel is a
    hair-thin strip in which two lines 10 m apart are one line. Cutting it into
    thirds and giving each panel its own extent triples the scale for free -
    that is the whole reason for the split, and why the panels deliberately do
    NOT share axes.
    """
    vmin, vmax = m.elev_limits(elev)
    ids = gdf["domain_id"].astype(int).to_numpy()
    edges = np.array_split(np.sort(ids), N_PANELS)

    fig, axes = plt.subplots(1, N_PANELS, figsize=(4.6 * N_PANELS, 17.0),
                             constrained_layout=True)
    im = None
    for ax, group in zip(np.atleast_1d(axes), edges):
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        bx = sub.total_bounds
        im = m.panel_elev(ax, elev, extent, vmin, vmax,
                          f"domains {group.min()}-{group.max()}")
        sub.boundary.plot(ax=ax, color=BOX_STYLE["edgecolor"],
                          linewidth=BOX_STYLE["linewidth"], zorder=4)
        for _, r in sub.iterrows():
            d = int(r["domain_id"])
            if d % LABEL_EVERY == 0 or d in (1, ids.max()):
                b = r.geometry.bounds
                ax.text(b[0] - 90, (b[1] + b[3]) / 2, str(d), fontsize=6,
                        ha="right", va="center", color="0.1", zorder=6)
        draw_lines(ax, lines)
        ax.set_xlim(bx[0] - PANEL_PAD_M, bx[2] + PANEL_PAD_M)
        ax.set_ylim(bx[1] - PANEL_PAD_M, bx[3] + PANEL_PAD_M)
        ax.set_aspect("equal")
        m.km_axes(ax, nx=3, ny=8)
        ax.set_xlabel("Easting (km)", fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("Northing (km)", fontsize=9)

    cb = fig.colorbar(im, ax=list(np.atleast_1d(axes)),
                      orientation="horizontal", fraction=0.02, pad=0.01,
                      aspect=55)
    cb.set_label("Elevation (m NAVD88)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    np.atleast_1d(axes)[0].legend(
        handles=line_legend() + [Line2D([0], [0], color=BOX_STYLE["edgecolor"],
                                        lw=0.8,
                                        label="Barrier3D domain, 2000 x 500 m")],
        loc="lower left", fontsize=8, framealpha=0.92)

    med = np.array([r["offset_med_m"] for r in rows
                    if r["offset_med_m"] != ""], float)
    fig.suptitle(f"1984 and 1997 dune lines on the {SOURCE_TAG} DEM",
                 fontsize=13)
    fig.text(0.5, -0.008,
             f"Cross-shore offset between the two lines, median over "
             f"{len(med)} domains: {np.median(med):+.1f} m, positive where the "
             f"1984 line lies seaward. Range {med.min():+.1f} to "
             f"{med.max():+.1f} m. Both lines are clipped to the domain "
             f"footprint for drawing only. Per-domain values in {CSV_NAME}.",
             ha="center", va="top", fontsize=8, color="#444444", wrap=True)

    p = FIG_DIR / "HAT_duneline_offset_island.png"
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


def _smooth(a, win_px):
    """Boxcar over an alongshore series, edges held rather than tapered."""
    k = np.ones(win_px) / win_px
    pad = win_px // 2
    return np.convolve(np.pad(a, pad, mode="edge"), k, mode="same")[
        pad:pad + a.size]


def fig_ribbon(samples, rows, gdf):
    """
    The two lines against a common baseline, with the band between them filled.

    WHY A BASELINE IS NEEDED. Plotted as raw easting the two lines sweep 6.5 km
    across the island's curve, which dwarfs a separation of tens of metres -
    the same reason they are indistinguishable on the map. Subtracting a
    SMOOTHED MIDLINE of the two removes the curve the two share and leaves what
    differs between them, at full 1 m alongshore resolution.

    The baseline is the mean of the two lines, boxcar-smoothed over
    BASELINE_WINDOW_M alongshore. It is a drawing device and carries no claim:
    it is symmetric in the two lines, so it cannot move one relative to the
    other, and the filled band's width is the offset exactly. What the choice
    DOES control is how much of each line's own sinuosity is left in the
    curves - a shorter window flattens both toward the axis, a longer one lets
    shared meanders back in. 2 km keeps four domains of context.
    """
    y = samples["y"]
    x84, x97 = samples["x1984"], samples["x1997"]
    km = (y - y.min()) / 1000.0
    win = max(int(BASELINE_WINDOW_M / SAMPLE_SPACING_M), 3)
    base = _smooth(np.nanmean(np.vstack([x84, x97]), axis=0), win)
    d84, d97 = x84 - base, x97 - base
    off = x84 - x97

    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(15.0, 8.6), sharex=True,
        gridspec_kw=dict(height_ratios=[2.4, 1.0]), constrained_layout=True)

    ax.fill_between(km, d97, d84, where=(off >= 0), interpolate=True,
                    color=LINE_STYLE[1984]["color"], alpha=0.30, linewidth=0,
                    label="1984 seaward of 1997")
    ax.fill_between(km, d97, d84, where=(off < 0), interpolate=True,
                    color=LINE_STYLE[1997]["color"], alpha=0.30, linewidth=0,
                    label="1984 landward of 1997")
    ax.plot(km, d97, **{**LINE_STYLE[1997], "label": "1997 dune line"})
    ax.plot(km, d84, **{**LINE_STYLE[1984], "label": "1984 dune line"})
    ax.axhline(0, color="0.55", linewidth=0.7, zorder=1)
    ax.set_ylabel("cross-shore position, m\n"
                  f"(from a {BASELINE_WINDOW_M / 1000:.0f} km smoothed "
                  f"midline; up = seaward)")
    ax.legend(fontsize=8, loc="upper left", ncol=4, framealpha=0.92)

    bx.fill_between(km, 0, off, where=(off >= 0), interpolate=True,
                    color=LINE_STYLE[1984]["color"], alpha=0.55, linewidth=0)
    bx.fill_between(km, 0, off, where=(off < 0), interpolate=True,
                    color=LINE_STYLE[1997]["color"], alpha=0.55, linewidth=0)
    bx.axhline(0, color="0.2", linewidth=0.9)
    for s in (-CELL_M, CELL_M):
        bx.axhline(s, color="0.2", linewidth=0.7, linestyle=":")
    bx.set_ylabel("1984 minus 1997, m")
    bx.set_xlabel("alongshore distance from the south end of domain 1 (km)")
    bx.text(0.2, CELL_M + 2, "one 10 m Barrier3D cell", fontsize=7,
            color="0.35", va="bottom")

    # The reaches the true-scale zooms cover, marked on both panels so the two
    # figures can be read against each other.
    ymid = {int(r["domain_id"]): r.geometry.bounds[1] + 250.0
            for _, r in gdf.iterrows()}
    for lo, hi, slug, _ in ZOOM_REACHES:
        a = (ymid[lo] - 250 - y.min()) / 1000.0
        b = (ymid[hi] + 250 - y.min()) / 1000.0
        for axis in (ax, bx):
            axis.axvspan(a, b, color="0.85", alpha=0.55, zorder=0)
        # Labelled along the BOTTOM of the upper panel. The top is where the
        # legend sits and where the 1984 excursions run, so a label there is
        # either hidden or hiding something.
        ax.text((a + b) / 2, ax.get_ylim()[0], f"zoom {slug} ", fontsize=7,
                ha="center", va="bottom", color="0.35")

    # Domain numbers on top, since the reader thinks in domains and the CSV is
    # keyed by them, but the axis itself stays metric.
    tx = ax.secondary_xaxis("top")
    ticks = [d for d in sorted(ymid) if d % 10 == 0 or d == 1]
    tx.set_xticks([(ymid[d] - y.min()) / 1000.0 for d in ticks])
    tx.set_xticklabels([str(d) for d in ticks], fontsize=8)
    tx.set_xlabel("domain", fontsize=9)

    med = np.array([r["offset_med_m"] for r in rows
                    if r["offset_med_m"] != ""], float)
    fig.suptitle(
        f"1984 and 1997 dune lines, cross-shore separation along the island   "
        f"|   median {np.median(med):+.1f} m, "
        f"{int((np.abs(med) >= CELL_M).sum())} of {len(med)} domains "
        f"at least one 10 m cell apart", fontsize=12)
    ax.margins(x=0.005)

    p = FIG_DIR / "HAT_duneline_offset_ribbon.png"
    fig.savefig(p, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_zooms(elev, extent, gdf, lines, rows, reaches=None, out=None,
              half_width=None, rows_by_domain=None, suptitle=None):
    """
    True-scale map panels on the reaches where the offset is largest, plus a
    quiet control.

    Equal aspect throughout - nothing is exaggerated. What makes the separation
    visible here is the CROSS-SHORE CROP: each panel is cut to `half_width`
    either side of the local line position instead of the full 2000 m domain
    box, so 50 m of offset is roughly a twelfth of the frame rather than a
    fortieth. The control reach is included so a reader can see what agreement
    looks like at the same scale, and not read every figure of this kind as
    showing a discrepancy.

    `reaches` overrides ZOOM_REACHES so an arbitrary span can be rendered to
    `out` - see --zoom. `rows_by_domain` is an optional {domain: N} mapping; if
    given, each label also carries the number of Barrier3D rows the insert
    would add there, which is what ties this view to row-insert-scope/.
    """
    reaches = reaches or ZOOM_REACHES
    half_width = ZOOM_HALF_WIDTH_M if half_width is None else half_width
    vmin, vmax = m.elev_limits(elev)
    by_dom = {r["domain"]: r for r in rows}

    fig, axes = plt.subplots(1, len(reaches),
                             figsize=(4.4 * len(reaches), 13.0),
                             constrained_layout=True)
    im = None
    for ax, (lo, hi, slug, note) in zip(np.atleast_1d(axes), reaches):
        ids = list(range(lo, hi + 1))
        sub = gdf[gdf["domain_id"].astype(int).isin(ids)]
        bx = sub.total_bounds
        cen = float(np.median([by_dom[d]["x1984_med"] for d in ids
                               if by_dom[d]["x1984_med"] != ""]))
        meds = [by_dom[d]["offset_med_m"] for d in ids]

        im = m.panel_elev(ax, elev, extent, vmin, vmax,
                          f"domains {lo}-{hi}\n{note}\n"
                          f"offset {min(meds):+.0f} to {max(meds):+.0f} m")
        sub.boundary.plot(ax=ax, color=BOX_STYLE["edgecolor"],
                          linewidth=0.6, zorder=4)
        for _, r in sub.iterrows():
            b = r.geometry.bounds
            d = int(r["domain_id"])
            lab = f"{d}   {by_dom[d]['offset_med_m']:+.0f} m"
            if rows_by_domain is not None:
                n = rows_by_domain.get(d, 0)
                lab += f"   →  +{n} row{'' if n == 1 else 's'}"
            ax.text(cen - half_width + 55, (b[1] + b[3]) / 2, lab,
                    fontsize=7.5, ha="left", va="center", color="0.05",
                    zorder=6,
                    bbox=dict(facecolor="white", alpha=0.7, edgecolor="none",
                              boxstyle="square,pad=0.15"))
        draw_lines(ax, lines, scale=1.3)
        ax.set_xlim(cen - half_width, cen + half_width)
        ax.set_ylim(bx[1] - 100.0, bx[3] + 100.0)
        ax.set_aspect("equal")
        m.km_axes(ax, nx=2, ny=8)
        ax.set_xlabel("Easting (km)", fontsize=9)
    np.atleast_1d(axes)[0].set_ylabel("Northing (km)", fontsize=9)
    np.atleast_1d(axes)[0].legend(handles=line_legend(), loc="lower left",
                                  fontsize=8, framealpha=0.92)

    cb = fig.colorbar(im, ax=list(np.atleast_1d(axes)),
                      orientation="horizontal", fraction=0.025, pad=0.01,
                      aspect=50)
    cb.set_label("Elevation (m NAVD88)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    fig.suptitle(suptitle
                 or "The same two lines at true scale, cross-shore crop only",
                 fontsize=13)
    tail = ("" if rows_by_domain is None else
            " The row count beside each is what the 1984 seaward-row insert "
            "would add there: round(offset / 10 m), floored at 0.")
    fig.text(0.5, -0.01,
             f"Equal aspect, no exaggeration. Each panel is cropped to "
             f"{half_width:.0f} m either side of the local line "
             f"position rather than the full 2000 m domain box. Per-domain "
             f"median offset is printed beside each box.{tail}",
             ha="center", va="top", fontsize=8, color="#444444", wrap=True)

    p = Path(out) if out else (FIG_DIR / "HAT_duneline_offset_zooms.png")
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(p, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_by_domain(rows):
    """The requested number, domain by domain, with its spread."""
    dom = np.array([r["domain"] for r in rows])
    med = np.array([r["offset_med_m"] if r["offset_med_m"] != "" else np.nan
                    for r in rows], float)
    p25 = np.array([r["offset_p25_m"] if r["offset_p25_m"] != "" else np.nan
                    for r in rows], float)
    p75 = np.array([r["offset_p75_m"] if r["offset_p75_m"] != "" else np.nan
                    for r in rows], float)

    fig, ax = plt.subplots(figsize=(11.0, 4.2), constrained_layout=True)
    ax.vlines(dom, p25, p75, color="0.75", linewidth=2.4,
              label="interquartile range within the domain", zorder=2)
    ax.scatter(dom, med, s=14, color="#d62728", zorder=3,
               label="median offset")
    ax.axhline(0, color="0.2", linewidth=0.9, zorder=1)
    for s in (-CELL_M, CELL_M):
        ax.axhline(s, color="0.2", linewidth=0.7, linestyle=":", zorder=1)
    ax.text(dom.min(), CELL_M + 1.5, "one 10 m Barrier3D cell", fontsize=7,
            va="bottom", ha="left", color="0.35")

    fin = np.isfinite(med)
    ax.set_xlabel("domain (1 = south)")
    ax.set_ylabel("1984 line minus 1997 line, m\n(+ve = 1984 seaward)")
    ax.set_title(f"Cross-shore offset between the 1984 and 1997 dune lines, "
                 f"per domain   |   island median {np.median(med[fin]):+.1f} m, "
                 f"{int((med[fin] > 0).sum())} of {int(fin.sum())} domains "
                 f"positive", fontsize=10)
    ax.legend(fontsize=8, loc="best")
    ax.margins(x=0.01)

    p = FIG_DIR / "HAT_duneline_offset_bydomain.png"
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return p


# =============================================================================
# ONE ARBITRARY ZOOM
# =============================================================================

def _epsg(crs):
    """A CRS as 'EPSG:nnnn', falling back to its name. Compound CRS WKT runs
    to ~1200 characters and makes the run log unreadable."""
    try:
        code = crs.to_epsg()
        if code:
            return f"EPSG:{code}"
        sub = getattr(crs, "sub_crs_list", None)
        if sub:
            code = sub[0].to_epsg()
            if code:
                return f"EPSG:{code} (compound)"
        return crs.name
    except Exception:
        return str(crs)[:40]


def load_domains():
    """
    The 90 domain boxes, from `domains.geojson` if it is reachable and from the
    resampled rasters if it is not.

    THE GEOJSON LIVES ON D:. That drive is an external disk - not
    version-controlled, not present on another machine, and it can disappear
    mid-session, which is exactly what happened on 2026-09-03. The fallback is
    not an approximation: every `resampled_domain_<N>_filled.tif` carries the
    snapped window this pipeline actually clipped, all 90 come back 2000 x 500 m
    at 200 x 50 cells, and domain 1 and 90's northings reproduce the geojson to
    the metre. The mosaic these figures draw is built from these same rasters,
    so the boxes agree with the pixels by construction.

    Preference order matters: the geojson stays authoritative when it is there,
    so this can never change a result on a machine that has the drive.
    """
    if m.DOMAIN_FILE.exists():
        g = gpd.read_file(m.DOMAIN_FILE).sort_values("domain_id")
        print(f"  {len(g)} domain boxes from {m.DOMAIN_FILE.name}, "
              f"{_epsg(g.crs)}")
        return g

    import re
    import rasterio
    from shapely.geometry import box as _box
    paths = sorted(m.IN_DIR.glob("resampled_domain_*_filled.tif"))
    if not paths:
        raise SystemExit(
            f"\n{m.DOMAIN_FILE} is not reachable and there are no resampled "
            f"rasters in\n    {m.IN_DIR}\nto fall back on. Reconnect the "
            f"drive, or rebuild the 10 m product.")
    ids, geoms, crs = [], [], None
    for p in paths:
        with rasterio.open(p) as s:
            t, crs = s.transform, s.crs
            geoms.append(_box(t.c, t.f - s.height * GRID_10M,
                              t.c + s.width * GRID_10M, t.f))
        ids.append(int(re.search(r"domain_(\d+)_", p.name).group(1)))
    g = gpd.GeoDataFrame({"domain_id": ids}, geometry=geoms, crs=crs)
    print(f"  WARNING: {m.DOMAIN_FILE} not reachable - domain boxes rebuilt "
          f"from {len(g)} resampled rasters instead ({_epsg(g.crs)})")
    return g.sort_values("domain_id").reset_index(drop=True)


def read_rows(path):
    """The per-domain offsets back off disk, so a zoom does not re-measure."""
    if not Path(path).exists():
        raise SystemExit(
            f"\n{path} not found.\n"
            f"  Run this script with no arguments first - the zoom reads the "
            f"table rather than re-measuring, which takes a minute.")
    out = []
    for r in csv.DictReader(open(path)):
        rec = {"domain": int(r["domain"])}
        for k in ("offset_med_m", "x1984_med"):
            rec[k] = float(r[k]) if r[k] not in ("", "nan") else ""
        out.append(rec)
    return out


def read_insert_rows(path):
    """{domain: N} from the row-insert scope table, if it has been built."""
    if not Path(path).exists():
        print(f"  NOTE: {Path(path).name} absent - row counts not labelled")
        return None
    return {int(r["domain"]): int(r["n_rows"]) for r in
            csv.DictReader(open(path))}


def zoom_only(span, out, half_width, note, with_rows):
    """
    Render ONE true-scale zoom for an arbitrary domain span.

    Reads the per-domain table rather than re-running measure(), so this is
    seconds rather than a minute. Everything else - the crop rule, the line
    styling, the equal aspect - is the same code path the three standard
    reaches use, so a custom zoom cannot quietly differ from them.
    """
    lo, hi = (int(x) for x in span.split("-"))
    rows = read_rows(OUT_DIR / CSV_NAME)
    n_by = read_insert_rows(INSERT_SCOPE_CSV) if with_rows else None

    print(f"\nloading the {SOURCE_TAG} mosaic at 10 m ...")
    elev, _surv, extent, _n = m.load_mosaic()
    gdf = load_domains()
    lines = load_lines(gdf.crs)
    drawn = clip_for_drawing(lines, gdf.union_all())

    by = {r["domain"]: r for r in rows}
    meds = [by[d]["offset_med_m"] for d in range(lo, hi + 1) if d in by]
    p = fig_zooms(elev, extent, gdf, drawn, rows,
                  reaches=[(lo, hi, f"{lo}-{hi}",
                            note or f"offset {min(meds):+.0f} to "
                                    f"{max(meds):+.0f} m")],
                  out=out, half_width=half_width,
                  rows_by_domain=n_by)
    # The suptitle is left at its default on purpose: the panel title already
    # names the span, and repeating it reads as two headings for one picture.
    print(f"\n  figure : {p}")
    return p


# =============================================================================
# MAIN
# =============================================================================

def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for p in DUNE_LINES.values():
        if not p.exists():
            raise SystemExit(f"dune line not found: {p}")

    print(f"\nloading the {SOURCE_TAG} mosaic at 10 m ...")
    elev, _surv, extent, n = m.load_mosaic()
    print(f"  {n} domains on one grid, {elev.shape[0]} x {elev.shape[1]}")

    gdf = load_domains()
    print("dune lines:")
    lines = load_lines(gdf.crs)

    print(f"\nmeasuring, {SAMPLE_SPACING_M:.0f} m alongshore sampling ...")
    geoms = {yr: g.union_all() for yr, g in lines.items()}
    rows, samples = measure(gdf, geoms)

    print("\nfor drawing:")
    drawn = clip_for_drawing(lines, gdf.union_all())

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cp = OUT_DIR / CSV_NAME
    with open(cp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for r in rows:
        print(f"  domain {r['domain']:>3}: offset {r['offset_med_m']:>8} m  "
              f"({r['offset_med_cells']:>6} cells)   IQR "
              f"{r['offset_p25_m']:>8} to {r['offset_p75_m']:>8}   "
              f"nearest {r['nearest_med_m']:>7} m   "
              f"n {r['n_both']:>3}/{r['n_samples']}")

    med = np.array([r["offset_med_m"] for r in rows
                    if r["offset_med_m"] != ""], float)
    nea = np.array([r["nearest_med_m"] for r in rows
                    if r["nearest_med_m"] != ""], float)
    print(f"\n{'=' * 78}\n{len(rows)} domains, {len(med)} with both lines")
    print(f"  offset, 1984 minus 1997, m (+ve = 1984 SEAWARD of 1997):")
    print(f"      min {med.min():+.1f}   p25 {np.percentile(med, 25):+.1f}   "
          f"median {np.median(med):+.1f}   p75 {np.percentile(med, 75):+.1f}   "
          f"max {med.max():+.1f}")
    print(f"      in 10 m cells: median {np.median(med) / CELL_M:+.2f}")
    print(f"      {int((med > 0).sum())} domains positive (1984 seaward), "
          f"{int((med < 0).sum())} negative")
    print(f"      |offset| under one 10 m cell in "
          f"{int((np.abs(med) < CELL_M).sum())} of {len(med)} domains")
    print(f"  nearest-point distance, orientation-free: "
          f"median {np.median(nea):.1f} m")

    figs = [fig_ribbon(samples, rows, gdf),
            fig_zooms(elev, extent, gdf, drawn, rows),
            fig_island(elev, extent, gdf, drawn, rows),
            fig_by_domain(rows)]
    print(f"\n  table  : {cp}")
    for f in figs:
        print(f"  figure : {f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Dune-line offset: measurement, table and figures. With "
                    "--zoom, renders one true-scale panel for an arbitrary "
                    "domain span and exits.")
    ap.add_argument("--zoom", default=None, metavar="LO-HI",
                    help="render ONE true-scale zoom for this domain span "
                         "(e.g. 83-87) and exit. Reads the existing table "
                         "rather than re-measuring.")
    ap.add_argument("--zoom-out", default=None,
                    help="output path for --zoom")
    ap.add_argument("--zoom-halfwidth", type=float, default=ZOOM_HALF_WIDTH_M,
                    help=f"cross-shore half-width of the crop, m "
                         f"(default {ZOOM_HALF_WIDTH_M:.0f})")
    ap.add_argument("--zoom-note", default=None,
                    help="subtitle line under the panel title")
    ap.add_argument("--no-row-labels", action="store_true",
                    help="do not annotate each domain with its insert row "
                         "count, even if the scope table is on disk")
    a = ap.parse_args()

    if a.zoom:
        zoom_only(a.zoom, a.zoom_out, a.zoom_halfwidth, a.zoom_note,
                  not a.no_row_labels)
    else:
        main()
