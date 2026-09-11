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
    figures/detail/   HAT_duneline_offset_simple.png   four two-domain pairs,
                                              grey relief, both lines solid,
                                              tight crop, no values
                      HAT_duneline_offset_zooms.png    three reaches of 5-8
                                              domains, same style, wider crop
                      HAT_duneline_offset_zoom_83_87.png   --zoom 83-87
    figures/island/   HAT_duneline_offset_simple_island.png   whole island,
                                              the measured offset as a bar
                                              beside the map (and _mean)
                      HAT_duneline_offset_lines_island.png  maps only, ~5 km
                                              panels cropped to the lines
                                              (and _3panel, 30 per panel)
    figures/offset/   HAT_duneline_offset_ribbon.png   both lines against a
                                              smoothed midline, 1 m sampling
                      HAT_duneline_offset_bydomain.png   the offset per domain
    figures/CAPTIONS.md                       a caption per figure, numbers
                                              filled from the table
    (the terrain-coloured locator HAT_duneline_offset_island.png was retired
    2026-09-08; fig_island_lines at 30 domains per panel replaces it)

STYLE
-----
Every figure in the folder is drawn to the one house style, which lives in
scripts/hat_figure_style.py and is re-exported here: a plain sans face,
8-10 pt type, thin dark-grey axes, a ColorBrewer red/blue pair for the two
lines that survives greyscale and colour-deficient print, panel letters, a
north arrow and a scale bar on the maps, and NO in-figure title sentences or
footnote paragraphs - what a figure needs said goes in figures/CAPTIONS.md,
which write_captions() fills from the same table the figures draw from.

Since 2026-09-10 every figure is also drawn at the width it will be PRINTED,
figsize("double") = 190 mm, so its 8-9 pt type is 8-9 pt on the page rather
than 4 pt after a journal reduces an 12-inch canvas. A panel is then one to
two inches wide, and what fitted on the old canvas does not: the panel letter
moves inside the corner on the island maps, the panel titles carry the domain
span alone with the place and the reading in the caption, the villages are
named vertically beside their bracket, and a figure whose panels all share one
scale gets ONE scale bar rather than one per panel. `save()` writes the PNG at
300 dpi, plus a PDF beside it for the two figures that are lines and bars
rather than shaded relief (offset/).

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
import matplotlib.patheffects as pe
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

# The place names are NOT redefined here. hatteras_site_config owns the
# community spans, the village centres and the end labels for the whole
# project - the same object the shoreline-rate figures annotate from - so
# a town that moves there moves here too, and this figure cannot quietly
# disagree with the rest of the repo about where Avon is.
from hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402


# =============================================================================
# STYLE
# =============================================================================
# One look for every figure. The STYLE block that lived here from 2026-09-04
# (Arial, thin dark-grey axes, the ColorBrewer RdBu poles for the two vintages,
# panel letters, nothing on the canvas that belongs in a caption) is now the
# project-wide standard in scripts/hat_figure_style.py, merged there on
# 2026-09-10 so every figure script can apply it. The names are re-exported
# here because a dozen scripts take both the style and the map loaders from
# this module as `off`.
from hat_figure_style import (  # noqa: E402,F401
    FONT_STACK, INK, INK_MUTED, GRID_C, C_1984, C_1997, C_1984_FILL, C_1997_FILL,
    STYLE_RC, apply_style, _letter, _title, _letter_inside, _north_arrow, _halo,
    figsize, FIG_W_DOUBLE, FIG_H_MAX, DOMAIN_AXIS_LABEL, town_bands, open_frame,
    save,
)
import hat_figure_style as _style  # noqa: E402


SOURCE_TAG = "2009-2014-1996"
OUT_DIR = ELEVATION_ROOT / f"{SOURCE_TAG}-duneline"
FIG_DIR = OUT_DIR / "figures"
# figures/ is sorted by what a figure IS (2026-09-08, Hannah): island/ for the
# whole-island maps, detail/ for the true-scale crops, offset/ for the two
# readings that are not maps. fig_path() is the only way a figure name becomes
# a path, so a figure cannot land at the folder root. Every map in the folder
# is drawn in the simple style - grey relief, both lines solid, red 1984 and
# blue 1997 - since the same date; the terrain-coloured locator and zooms are
# gone (the locator retired, the zooms redrawn through fig_zooms_simple).
FIG_SUBFOLDER = {
    "HAT_duneline_offset_simple_island.png": "island",
    "HAT_duneline_offset_simple_island_mean.png": "island",
    "HAT_duneline_offset_lines_island.png": "island",
    "HAT_duneline_offset_lines_island_3panel.png": "island",
    "HAT_duneline_offset_simple.png": "detail",
    "HAT_duneline_offset_zooms.png": "detail",
    "HAT_duneline_offset_zoom_83_87.png": "detail",
    "HAT_duneline_offset_ribbon.png": "offset",
    "HAT_duneline_offset_bydomain.png": "offset",
}


def fig_path(name):
    """figures/<kind>/<name>, the folder made; raises on a name not in the map."""
    p = FIG_DIR / FIG_SUBFOLDER[name] / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
CSV_NAME = "duneline_offset_by_domain.csv"

DUNE_DIR = INIT_ROOT / "1-barrier3d-domains" / "raw-duneline-geojson"
DUNE_LINES = {1984: DUNE_DIR / "duneline_1984.geojson",
              1997: DUNE_DIR / "duneline_1997.geojson"}

# One sample per metre of alongshore, matching the 1 m DEM the rest of the
# 1984-start chain is built on. 500 per domain.
SAMPLE_SPACING_M = 1.0
GRID_10M = 10.0          # the resampled product's cell, for the box fallback

# The 1984 footprint table, read ONLY to label a --zoom with how many
# Barrier3D rows the measured offset turns into. Optional; absent is fine.
# Since 2026-09-07 this is the SYMMETRIC footprint (HAT_footprint_1984.py):
# n_cells is signed, + rows added, - existing rows removed, trunc(shift/10).
INSERT_SCOPE_CSV = (INIT_ROOT / "1-barrier3d-domains" / "1984-start"
                    / "2-domain-reconstruction-1984" / "2-extent" / "footprint_1984_by_domain.csv")   # step folder since 2026-09-09

# The box shape the easting-is-cross-shore frame depends on. Checked, not
# assumed - see THE FRAME above.
EXPECTED_BOX_M = (2000.0, 500.0)
BOX_TOL_M = 1.0

# 1984 is the older line and the one the model starts from, so it gets the
# emphasis colour; 1997 is the reference. Deliberately NOT the road key from
# HAT_plot_1984_mosaic - these are dune lines, and reusing black/white-dashed
# would read as NC-12 on a figure where NC-12 is absent.
#
# ONE key for every figure in this folder. Both lines SOLID, at the weight the
# simple figure uses, and colour is the only thing that separates them. The
# 1997 line used to be dashed on the DEM figures and solid on the simple ones,
# which meant the same two lines carried two keys across one folder; the simple
# figure is the reference and this is now it. What still varies per figure is
# the WIDTH, through draw_lines(scale=...) - a 46 km locator and a 300 m crop
# cannot carry the same line weight - and the two whole-island figures use the
# same scale as each other, as do the two detail figures.
LINE_STYLE = {1984: dict(color=C_1984, linestyle="-", linewidth=2.0),
              1997: dict(color=C_1997, linestyle="-", linewidth=2.0)}
LINE_CASING = {1984: dict(color="white", linewidth=3.0),
               1997: dict(color="white", linewidth=3.0)}
LINE_ORDER = [1997, 1984]      # 1984 drawn last, on top

# The per-figure line WIDTH, as a multiplier on LINE_STYLE. Paired on purpose:
# the two whole-island figures share one, the two detail figures share the
# other, so a change to either cannot land on one of a pair and not the other.
# 46 km of island next to a 70 m offset puts the two lines inside one line
# width, and a heavier line there merges them further; a 300 m crop has room.
LINE_SCALE_ISLAND = 0.55
LINE_SCALE_DETAIL = 1.0
# The ribbon is a trace, not a map: 46 km of 1 m sampling on one axis, and
# the two lines cross constantly, so they are drawn finer than anywhere
# else. Same key, same colours, both solid - width only.
LINE_SCALE_RIBBON = 0.7

BOX_STYLE = dict(edgecolor="0.30", facecolor="none", linewidth=0.4)
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
    (62, 68, "62-68", "1984 line landward of 1997"),
    (78, 85, "78-85", "1984 line seaward of 1997"),
]
# Cross-shore half-width of a zoom panel. The domain box is 2000 m across and
# almost all of it is water and back-barrier; cropping to the dune makes the
# offset a visible fraction of the frame at equal aspect.
ZOOM_HALF_WIDTH_M = 300.0

# -----------------------------------------------------------------------------
# THE SIMPLE ZOOM  (--simple)
# -----------------------------------------------------------------------------
# A stripped version of the same panels: no elevation values, no colour ramp,
# no colourbar, no coordinate ticks. Both lines SOLID. What is left is the two
# lines, the ridge they sit on, and a scale bar.
#
# WHY THE PAIRS ARE SHORTER THAN THE ZOOM REACHES. The island runs about 7 deg
# oblique to the UTM grid, so a dune line drifts ~130 m in easting per km of
# alongshore. A crop tight enough to show a 50 m offset therefore cannot hold a
# five-domain reach - over 1.5 km the two lines sweep 190-330 m in easting and
# walk straight out of the frame. Two neighbouring domains sweep 150-223 m,
# which fits inside +/-150 m with margin. So each panel is the two-domain pair
# carrying that reach's extreme, not the whole reach.
#
# ONE half-width for all three panels, not one per panel. The control only
# works if it is drawn at exactly the scale of the other two.
SIMPLE_HALF_WIDTH_M = 150.0
# South to north, matching the locator. Each pair's lines were checked to sit
# inside +/-150 m of the pair's median 1984 easting before it was chosen: 3-4
# needs 115 m, 19-20 needs 88, 63-64 needs 112, 79-80 needs 75. 4-5 carries the
# south's single largest offset (+62 m at domain 5) and is NOT used, because
# its lines need 196 m and would have forced a wider crop on all four panels.
# The third field is a ROLE marker, not a description. Where a pair sits and
# which way its offset goes are both DERIVED - the place from
# HATTERAS_ANNOTATIONS, the direction and magnitude from the measured table -
# so a panel title cannot state a direction its own numbers contradict. Only
# the reason a pair is in the figure at all is written by hand.
SIMPLE_REACHES = [
    (3, 4, ""),
    (19, 20, "control"),
    (63, 64, ""),
    (79, 80, ""),
]

# THE ISLAND LOCATOR. Two things side by side per third of the island: a true
# map, and a bar of the measured offset aligned to it row for row.
#
# The map ALONE cannot carry this and no styling fixes that - at equal aspect
# 46 km of island next to a 50 m offset is half a line width, so on the map the
# two lines are one line nearly everywhere. The map is therefore a LOCATOR: it
# says where each detail pair sits and how the island is shaped. The bar beside
# it is what carries the magnitude, and it is the same median that the detail
# panels print, off the same CSV.
ISLAND_PANELS = 3
ISLAND_PAD_M = 500.0
STRIP_WIDTH_RATIO = 0.42       # bar axes width, as a fraction of the map's
OFFSET_POS = C_1984            # 1984 seaward - the sign erosion predicts
OFFSET_NEG = C_1997            # 1984 landward
HIGHLIGHT = INK                # the label and the band marking a detail pair

# Piers and the groin are drawn SEAWARD from the 1984 line, because that is
# where they are. The length is a drawing constant, not a measurement - none of
# these structures has a surveyed length in this repo, and a 46 km panel could
# not show the difference between 200 m and 400 m of pier anyway. They are
# drawn as marks and named in the legend rather than labelled on the map: at
# this scale a label beside the Rodanthe pier lands on the Rodanthe village
# tick and on the 79-80 detail label, and no amount of nudging fixes three
# labels inside one 500 m domain.
STRUCTURE_LEN_M = 380.0

# THE LINES-ONLY LOCATOR (--lines-island). The same whole island with the bar
# strips dropped, and the maps made to carry the offset themselves by
# CROPPING, not by styling. Two things change against the three-panel figure:
#
#   more panels   each covers LINES_ISLAND_DOMAINS domains (~5 km) instead of
#                 30 (~15 km), so at the same panel height a 50 m offset is
#                 about three line widths rather than under one.
#   tight crop    each panel's easting window is the envelope of the two lines
#                 inside that panel's northing window, plus LINES_ISLAND_PAD_M
#                 either side, rather than the whole island. The island runs
#                 ~7 deg oblique to the grid, so over 5 km the lines sweep
#                 ~650 m in easting, and the window is roughly 1 km wide.
#
# Equal aspect is kept, so nothing is stretched: what the reader sees is the
# true separation, just at a scale where it is visible.
LINES_ISLAND_DOMAINS = 10
LINES_ISLAND_PAD_M = 150.0
LINES_ISLAND_LINE_SCALE = 0.9

# The alongshore ruler on the bar strip, in km from the south end of domain 1 -
# the same origin and direction fig_ribbon's x-axis uses, so the two figures
# can be read against each other.
KM_TICK_M = 5000.0
# Kept as a name because the simple figures pass it explicitly, but it is no
# longer a SEPARATE key - LINE_STYLE is now what this used to be, so every
# figure in the folder draws the same two solid lines. Copied rather than
# aliased so a caller mutating one cannot reach the other.
SIMPLE_LINE_STYLE = {yr: dict(LINE_STYLE[yr]) for yr in LINE_STYLE}

# The 1 m gapfilled tiles, read in place. fig_zooms draws the 10 m resampled
# mosaic, which across a 300 m crop is 30 cells wide and renders the dune as a
# staircase. These panels are tight enough to be worth the 1 m source. NOTE the
# tiles carry NO CRS tag; their bounds are checked against the domain boxes
# instead, and a mismatch is fatal rather than silent.
TILE_1M_DIR = ELEVATION_ROOT / SOURCE_TAG / "1-gapfill-1m"
TILE_1M_NAME = "clip_domain_{d}_filled.tif"
TILE_BOUNDS_TOL_M = 1.5

# Relief only, no values. Both of these are DRAWING parameters and neither can
# move either line - they change how the backdrop is shaded, nothing else.
#
#   vert_exag  the dune is ~5 m of relief over ~50 m cross-shore, so at 1:1 the
#              shading is nearly flat grey.
#   SHADE_SMOOTH_M  the shading is computed off a boxcar-smoothed COPY of the
#              DEM. 1 m lidar over a vegetated backdune is speckly at the cell
#              scale, and exaggerating the slope exaggerates the speckle with
#              it until the dune ridge is lost in it. The smoothing is applied
#              to the shading only; the elevation array itself is untouched,
#              and nothing here is measured off the backdrop anyway.
HILLSHADE = dict(azdeg=315.0, altdeg=40.0, vert_exag=2.2)
SHADE_SMOOTH_M = 3.0
NODATA_GREY = "#e8e8e8"

SCALEBAR_M = 50.0              # 5 Barrier3D cells


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


def draw_lines(ax, lines, scale=1.0, style=None):
    """Casing then line, in LINE_ORDER so 1984 lands on top of 1997.

    `style` swaps the per-year style dict - the simple figure draws both lines
    solid. The draw order and the casing are shared, so the two figures cannot
    disagree about which line is on top.
    """
    style = style or LINE_STYLE
    for yr in LINE_ORDER:
        cas = dict(LINE_CASING[yr])
        cas["linewidth"] *= scale
        st = dict(style[yr])
        st["linewidth"] *= scale
        lines[yr].plot(ax=ax, linestyle="-", alpha=0.9,
                       zorder=6 + LINE_ORDER.index(yr) * 2, **cas)
        lines[yr].plot(ax=ax, zorder=7 + LINE_ORDER.index(yr) * 2, **st)


def line_legend():
    return [Line2D([0], [0], label=f"{yr} dune line", **LINE_STYLE[yr])
            for yr in (1984, 1997)]


def fig_island(elev, extent, gdf, lines, rows):
    """
    RETIRED 2026-09-08 (not called): the terrain-coloured locator, superseded
    by fig_island_lines with 30 domains per panel, which shows the same boxes
    and lines on grey relief. Kept so the drawing is on record.

    The DEM, both dune lines, and the domain boxes.

    THREE PANELS, each a third of the island, rather than one frame. At equal
    aspect the island is 46 km long and about 2 km wide, so a single panel is a
    hair-thin strip in which two lines 10 m apart are one line. Cutting it into
    thirds and giving each panel its own extent triples the scale for free -
    that is the whole reason for the split, and why the panels deliberately do
    NOT share axes.

    The three panels DO share one northing span - the longest third, padded -
    so they are the same height, the same scale, and line up top and bottom.
    Each panel's width is then its own easting span at that scale, which is
    why the width ratios are computed rather than equal.
    """
    apply_style()
    vmin, vmax = m.elev_limits(elev)
    ids = gdf["domain_id"].astype(int).to_numpy()
    edges = np.array_split(np.sort(ids), N_PANELS)
    bounds = [gdf[gdf["domain_id"].astype(int).isin(g)].total_bounds
              for g in edges]
    span = max(b[3] - b[1] for b in bounds) + 2 * PANEL_PAD_M
    widths = [(b[2] - b[0] + 2 * PANEL_PAD_M) / span for b in bounds]

    fig, axes = plt.subplots(
        1, N_PANELS, figsize=figsize("double", height=FIG_H_MAX),
        constrained_layout=True, gridspec_kw=dict(width_ratios=widths))
    axes = np.atleast_1d(axes)
    im = None
    for i, (ax, group, bx) in enumerate(zip(axes, edges, bounds)):
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        im = m.panel_elev(ax, elev, extent, vmin, vmax, "")
        _title(ax, i, f"Domains {group.min()}\u2013{group.max()}")
        sub.boundary.plot(ax=ax, color=BOX_STYLE["edgecolor"],
                          linewidth=BOX_STYLE["linewidth"], zorder=4)
        for _, r in sub.iterrows():
            d = int(r["domain_id"])
            if d % LABEL_EVERY == 0 or d in (1, ids.max()):
                b = r.geometry.bounds
                ax.text(b[0] - 80, (b[1] + b[3]) / 2, str(d), fontsize=6.5,
                        ha="right", va="center", color=INK_MUTED, zorder=6)
        draw_lines(ax, lines, scale=LINE_SCALE_ISLAND)
        yc = (bx[1] + bx[3]) / 2
        ax.set_xlim(bx[0] - PANEL_PAD_M, bx[2] + PANEL_PAD_M)
        ax.set_ylim(yc - span / 2, yc + span / 2)
        ax.set_aspect("equal")
        m.km_axes(ax, nx=3, ny=8)
        ax.set_xlabel("Easting (km)")
    axes[0].set_ylabel("Northing (km)")

    cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                      fraction=0.02, pad=0.015, aspect=60)
    cb.set_label("Elevation (m NAVD88)")
    cb.outline.set_linewidth(0.6)
    cb.ax.tick_params(labelsize=8, width=0.6, length=3)

    axes[0].legend(
        handles=line_legend() + [Line2D([0], [0], color=BOX_STYLE["edgecolor"],
                                        lw=0.8,
                                        label="Barrier3D domain, 2000 \u00d7 500 m")],
        loc="lower left")

    p = FIG_DIR / "HAT_duneline_offset_island.png"
    save(fig, p, vector=False, bbox_inches="tight")
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

    THE ALONGSHORE AXIS IS IN DOMAINS, continuously: each 1 m sample sits at
    its own domain's id plus its fraction of that 500 m box, so domain d spans
    d - 0.5 .. d + 0.5 and the village bands land exactly where they do on
    every other alongshore chart in the project. The boxes are not perfectly
    contiguous (502-507 m apart), which is why the position is built per
    sample from its own box rather than from a fixed pitch. Distance in km
    from the south end of domain 1 rides on a second axis along the top -
    the origin the island figure's ruler uses - so the two can still be read
    against each other.
    """
    apply_style()
    y = samples["y"]
    x84, x97 = samples["x1984"], samples["x1997"]
    km = (y - y.min()) / 1000.0
    ymid = {int(r["domain_id"]): r.geometry.bounds[1] + 250.0
            for _, r in gdf.iterrows()}
    dom = samples["domain"].astype(int)
    xd = dom + (y - np.array([ymid[d] for d in dom])) / EXPECTED_BOX_M[1]
    win = max(int(BASELINE_WINDOW_M / SAMPLE_SPACING_M), 3)
    base = _smooth(np.nanmean(np.vstack([x84, x97]), axis=0), win)
    d84, d97 = x84 - base, x97 - base
    off = x84 - x97

    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=figsize("double", aspect=0.58), sharex=True,
        gridspec_kw=dict(height_ratios=[2.3, 1.0]), constrained_layout=True)

    # ---- (a) the two lines about the midline ------------------------------
    ax.fill_between(xd, d97, d84, where=(off >= 0), interpolate=True,
                    color=C_1984_FILL, linewidth=0, zorder=2)
    ax.fill_between(xd, d97, d84, where=(off < 0), interpolate=True,
                    color=C_1997_FILL, linewidth=0, zorder=2)
    for yr, d in ((1997, d97), (1984, d84)):
        st = dict(LINE_STYLE[yr])
        st["linewidth"] = 0.9
        ax.plot(xd, d, zorder=3, **st)
    ax.axhline(0, color=INK_MUTED, linewidth=0.6, zorder=1)
    ax.set_ylabel("Cross-shore position (m)\nabout the smoothed midline")

    # ---- (b) the difference ----------------------------------------------
    bx.fill_between(xd, 0, off, where=(off >= 0), interpolate=True,
                    color=C_1984, alpha=0.75, linewidth=0, zorder=2)
    bx.fill_between(xd, 0, off, where=(off < 0), interpolate=True,
                    color=C_1997, alpha=0.75, linewidth=0, zorder=2)
    bx.axhline(0, color=INK, linewidth=0.7, zorder=3)
    for s_ in (-CELL_M, CELL_M):
        bx.axhline(s_, color=INK_MUTED, linewidth=0.6, linestyle=(0, (3, 2)),
                   zorder=1)
    bx.set_ylabel("1984 minus 1997 (m)")
    bx.set_xlabel(DOMAIN_AXIS_LABEL)

    # Headroom: the village names sit along the top of (a) and the detail
    # brackets along the top of (b), and neither may land on the data.
    ax.set_xlim(dom.min() - 0.5, dom.max() + 0.5)
    for a_, frac in ((ax, 0.14), (bx, 0.30)):
        lo_, hi_ = a_.get_ylim()
        a_.set_ylim(lo_, hi_ + frac * (hi_ - lo_))
    town_bands(ax)
    town_bands(bx, label=False)

    # The reaches the true-scale zooms cover, as brackets along the top of
    # (b): a band would be the same grey as the villages behind it.
    for lo, hi, slug, _ in ZOOM_REACHES:
        bx.plot([lo - 0.5, hi + 0.5], [0.95, 0.95], color=INK_MUTED, lw=0.9,
                solid_capstyle="butt", zorder=4,
                transform=bx.get_xaxis_transform())
        bx.text((lo + hi) / 2, 0.92, f"detail {slug}", fontsize=7,
                ha="center", va="top", color=INK_MUTED,
                transform=bx.get_xaxis_transform())

    # Distance along the top, in km from the south end of domain 1; the ticks
    # are placed by interpolating the samples' own (km, domain) pairs.
    tx = ax.secondary_xaxis("top")
    kt = np.arange(0.0, km.max() + 1e-9, KM_TICK_M / 1000.0)
    tx.set_xticks(np.interp(kt, km, xd))
    tx.set_xticklabels([f"{k:.0f}" for k in kt])
    tx.set_xlabel("Alongshore distance from the south end of domain 1 (km)")
    tx.tick_params(labelsize=8, color=INK, width=0.6, length=3)
    ax.set_xticks([1] + list(range(10, int(dom.max()) + 1, 10)))
    for a_ in (ax, bx):
        open_frame(a_)
        a_.grid(axis="y", zorder=0)
        a_.set_axisbelow(True)
    for i, a_ in enumerate((ax, bx)):
        a_.set_title(_letter(i), loc="left", fontweight="bold")

    from matplotlib.patches import Patch
    fig.legend(loc="outside lower center", ncol=5, frameon=False, handles=[
        Line2D([0], [0], color=C_1984, lw=1.6, label="1984 dune line"),
        Line2D([0], [0], color=C_1997, lw=1.6, label="1997 dune line"),
        Patch(color=C_1984_FILL, label="1984 seaward of 1997"),
        Patch(color=C_1997_FILL, label="1984 landward of 1997"),
        Line2D([0], [0], color=INK_MUTED, lw=0.8, ls=(0, (3, 2)),
               label="\u00b11 Barrier3D cell (10 m)")])

    p = fig_path("HAT_duneline_offset_ribbon.png")
    save(fig, p, bbox_inches="tight")
    plt.close(fig)
    return p


def fig_zooms(elev, extent, gdf, lines, rows, reaches=None, out=None,
              half_width=None, rows_by_domain=None):
    """
    True-scale panels on the reaches where the offset is largest, plus a quiet
    control - since 2026-09-08 drawn by fig_zooms_simple, so they carry the
    same grey relief and the same two solid lines as every other map in the
    folder (Hannah: one style for the whole folder). What this keeps from the
    original is the REACH definition (ZOOM_REACHES, five to eight domains) and
    the wider crop (ZOOM_HALF_WIDTH_M), which is why it is still a separate
    figure from the two-domain pairs.

    Equal aspect throughout - nothing is exaggerated. What makes the separation
    visible is the cross-shore crop: each panel is cut to `half_width` either
    side of the local line position instead of the full 2000 m domain box. The
    control reach is included so a reader can see what agreement looks like at
    the same scale.

    `reaches` overrides ZOOM_REACHES so an arbitrary span can be rendered to
    `out` - see --zoom. `rows_by_domain` is an optional {domain: N} mapping; if
    given, each label also carries the number of Barrier3D rows the 1984
    footprint adds or removes there, which ties this view to 2-domain-reconstruction-1984/.
    """
    # The standard reaches' notes ("the quietest reach on the island") are
    # caption text and are written there; a --zoom-note given by hand is the
    # one thing drawn under a panel title.
    custom = reaches is not None
    reaches = reaches or ZOOM_REACHES
    half_width = ZOOM_HALF_WIDTH_M if half_width is None else half_width
    return fig_zooms_simple(gdf, lines, rows, elev=elev, extent=extent,
                            reaches=[(lo, hi, note if custom else "")
                                     for lo, hi, _slug, note in reaches],
                            out=out or fig_path("HAT_duneline_offset_zooms.png"),
                            half_width=half_width, rows_by_domain=rows_by_domain,
                            tall=True)


def load_1m(gdf, ids):
    """
    The 1 m gapfilled tiles for `ids`, mosaicked onto one array.

    The tiles carry no CRS tag, so this does NOT trust them: each tile's bounds
    are checked against its domain box from `gdf`, and a disagreement over
    TILE_BOUNDS_TOL_M raises. That check is the only thing establishing that
    the raster and the dune lines are in the same frame, so it is not optional.

    Returns (array, extent) with nodata as NaN, or (None, None) if the tiles
    are not on disk - the caller falls back to the 10 m mosaic.
    """
    import rasterio
    paths = {d: TILE_1M_DIR / TILE_1M_NAME.format(d=d) for d in ids}
    missing = [d for d, q in paths.items() if not q.exists()]
    if missing:
        print(f"  NOTE: 1 m tiles absent for domains {missing} - "
              f"falling back to the 10 m mosaic")
        return None, None

    tiles = {}
    for d, q in paths.items():
        box = gdf[gdf["domain_id"].astype(int) == d].total_bounds
        with rasterio.open(q) as src:
            b = src.bounds
            off = max(abs(b.left - box[0]), abs(b.bottom - box[1]),
                      abs(b.right - box[2]), abs(b.top - box[3]))
            if off > TILE_BOUNDS_TOL_M:
                raise SystemExit(
                    f"\n{q.name} does not sit on domain {d}'s box - the worst "
                    f"corner is {off:.1f} m out (tolerance "
                    f"{TILE_BOUNDS_TOL_M} m).\n"
                    f"  tile {[round(v, 1) for v in b]}\n"
                    f"  box  {[round(v, 1) for v in box]}\n"
                    f"  The tiles carry no CRS, so this check is the only "
                    f"thing tying them to the dune lines.")
            if src.res != (1.0, 1.0):
                raise SystemExit(f"{q.name} is {src.res} m, not 1 m")
            a = src.read(1).astype(float)
            a[a == src.nodata] = np.nan
            tiles[d] = (a, b)

    x0 = min(b.left for _, b in tiles.values())
    x1 = max(b.right for _, b in tiles.values())
    y0 = min(b.bottom for _, b in tiles.values())
    y1 = max(b.top for _, b in tiles.values())
    out = np.full((int(round(y1 - y0)), int(round(x1 - x0))), np.nan)
    for a, b in tiles.values():
        r0 = int(round(y1 - b.top))
        c0 = int(round(b.left - x0))
        win = out[r0:r0 + a.shape[0], c0:c0 + a.shape[1]]
        np.copyto(win, a, where=np.isnan(win))
    print(f"  backdrop: {len(tiles)} tiles at 1 m, {out.shape[0]} x "
          f"{out.shape[1]}, {int(np.isfinite(out).sum()):,} valid cells")
    return out, (x0, x1, y0, y1)


def _hillshade(ax, arr, extent, res=1.0):
    """
    Grey relief. No value is readable off this and none is labelled.

    `res` is the array's cell size in metres. It has to be passed, not assumed:
    the detail panels shade the 1 m tiles and the locator shades the 10 m
    mosaic, and a hillshade computed at the wrong cell size gets the slope - and
    so the whole look of the relief - wrong by that factor. The smoothing is
    expressed in metres and converted here for the same reason; below one cell
    it is skipped rather than rounded to a no-op filter.
    """
    from matplotlib.colors import LightSource
    ax.set_facecolor(NODATA_GREY)
    good = np.isfinite(arr)
    filled = np.where(good, arr, np.nanmin(arr))
    win = int(round(SHADE_SMOOTH_M / res))
    if win > 1:
        from scipy.ndimage import uniform_filter
        filled = uniform_filter(filled, size=win, mode="nearest")
    hs = LightSource(azdeg=HILLSHADE["azdeg"],
                     altdeg=HILLSHADE["altdeg"]).hillshade(
        filled, vert_exag=HILLSHADE["vert_exag"], dx=res, dy=res)
    rgba = plt.get_cmap("gray")(0.30 + 0.62 * hs)
    rgba[~good] = matplotlib.colors.to_rgba(NODATA_GREY)
    ax.imshow(rgba, extent=(extent[0], extent[1], extent[2], extent[3]),
              origin="upper", interpolation="bilinear", zorder=1)


def _place_of(lo, hi, ann=HATTERAS_ANNOTATIONS):
    """
    Where GIS domains lo..hi are, in the site's own vocabulary.

    Every name comes from HATTERAS_ANNOTATIONS, so none of them is a place name
    invented for this figure. Most specific first: a community containing the
    pair, then a named shoal zone, then the gap between the two nearest
    communities. A pair inside a community also names the village centre it
    sits on, where the config gives one.
    """
    for name, (a, b) in ann.town_spans.items():
        if a <= lo and hi <= b:
            near = [v for v, g in ann.village_lines.items() if lo <= g <= hi]
            return f"{name}, at {near[0]}" if near else name
    for name, (a, b) in ann.shoal_zones.items():
        if a <= lo and hi <= b:
            return name
    below = [(b, n) for n, (a, b) in ann.town_spans.items() if b < lo]
    above = [(a, n) for n, (a, b) in ann.town_spans.items() if a > hi]
    if below and above:
        return f"{max(below)[1]}–{min(above)[1]}"
    if above:
        return f"south of {min(above)[1]}"
    if below:
        return f"north of {max(below)[1]}"
    return ann.region_name


def _pair_reading(meds):
    """
    What was measured on a pair or reach, as caption text: "1984 seaward by
    30-62 m (3-6 cells)".

    Generated, not written down. The direction word comes from the SIGN of the
    measured medians, so the text cannot read SEAWARD over numbers that are
    negative, and the cell count is the same round(offset / 10 m) the
    row-insert scope uses. Until 2026-09-10 this was the second line of every
    detail panel's title; at the printed width (four panels across 190 mm)
    that line ran into its neighbours, so it now goes under the figure in
    CAPTIONS.md and the panel title is the domain span alone.
    """
    meds = np.asarray(meds, float)
    a = np.abs(meds)
    if (meds > 0).all():
        way = "1984 seaward"
    elif (meds < 0).all():
        way = "1984 landward"
    else:
        way = "1984 both ways"
    if a.max() < CELL_M:
        mag = f"by {a.min():.0f}–{a.max():.0f} m (under one cell)"
    else:
        c0, c1 = round(a.min() / CELL_M), round(a.max() / CELL_M)
        cell = f"{c0:.0f} cell" if c0 == c1 else f"{c0:.0f}–{c1:.0f} cells"
        mag = f"by {a.min():.0f}–{a.max():.0f} m ({cell})"
    return way + " " + mag


def _span_y(gdf, lo, hi):
    """Northing bounds of GIS domains lo..hi, or None if none are on the grid.

    Fractional ids are accepted so a groin on a box boundary resolves, but the
    only thing this figure uses is whole spans.
    """
    sub = gdf[gdf["domain_id"].astype(int).isin(
        range(int(np.floor(lo)), int(np.ceil(hi)) + 1))]
    if sub.empty:
        return None
    b = sub.total_bounds
    return float(b[1]), float(b[3])


def _places(ax, gdf, y0, y1, ann=HATTERAS_ANNOTATIONS):
    """
    The communities, as a bracket in the ocean margin - not a wash.

    A translucent band across the panel would sit on the island and on both
    dune lines, and its colour is a blue close enough to the 1997 line's to be
    read as belonging to it. A bracket at the seaward edge says the same thing
    where nothing else is drawn, and leaves the two lines the only coloured
    marks on the island itself.

    Spans, names and the end labels all come from HATTERAS_ANNOTATIONS. Only
    the parts falling inside this panel's northing window are drawn, and a
    bracket is clipped to that window rather than dropped, so a community
    straddling a panel break appears on both.
    """
    xa, xb = ax.get_xlim()
    xbr = xb - 0.055 * (xb - xa)          # the bracket
    xtx = xbr - 0.02 * (xb - xa)          # community names, landward of it
    # Village names get a column of their own. Both sets run vertically, and
    # a village near the middle of its community (Waves in Tri-Village) put
    # the two names on the same line of text when they shared an x.
    xtv = xbr - 0.11 * (xb - xa)

    for name, (lo, hi) in ann.town_spans.items():
        yy = _span_y(gdf, lo, hi)
        if yy is None or yy[1] < y0 or yy[0] > y1:
            continue
        a, b = max(yy[0], y0), min(yy[1], y1)
        ax.plot([xbr, xbr], [a, b], color=ann.color_town_span, lw=5.0,
                solid_capstyle="butt", zorder=8)
        # Avon is GIS 21-31 and the panels break at 30, so it is drawn in two
        # pieces. The BRACKET is drawn in both - the community really does
        # continue past the break - but the NAME goes only on the panel holding
        # most of it, or a one-domain sliver at the foot of the next panel
        # reads as a second Avon. Centred on the visible piece, never on the
        # clip edge, so a clipped label cannot ride up over the panel title.
        if (b - a) < 0.5 * (yy[1] - yy[0]):
            continue
        # Rotated to run along the bracket: at the printed width a panel is
        # one to two inches across, and a horizontal name would lie over the
        # two lines the panel exists to show.
        ax.text(xtx, (a + b) / 2.0, name, fontsize=8, fontweight="bold",
                ha="right", va="center", rotation=90, color=INK, zorder=11,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                          boxstyle="square,pad=0.15"))

    for name, gid in ann.village_lines.items():
        yy = _span_y(gdf, gid, gid)
        if yy is None:
            continue
        yc = (yy[0] + yy[1]) / 2.0
        if not (y0 <= yc <= y1):
            continue
        ax.plot([xbr - 0.012 * (xb - xa), xbr + 0.012 * (xb - xa)], [yc, yc],
                color=ann.color_village_line, lw=1.0, zorder=9)
        ax.text(xtv, yc, name, fontsize=7, ha="right", va="center",
                rotation=90, color=ann.color_village_line, zorder=11,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                          boxstyle="square,pad=0.1"))


def _x1984(by_dom, gid):
    """Median easting of the 1984 line at GIS id `gid`, fractional allowed.

    A structure sitting on a box boundary (the Buxton groin is at 5.5) has no
    domain of its own, so the two it lies between are averaged. Returns None
    where the table has no line in those domains, and the caller draws nothing
    rather than guessing a position.
    """
    import math
    ids = {math.floor(gid), math.ceil(gid)}
    xs = [by_dom[d]["x1984_med"] for d in ids
          if d in by_dom and by_dom[d]["x1984_med"] != ""]
    return float(np.mean(xs)) if xs else None


def _structures(ax, gdf, by_dom, y0, y1, ann=HATTERAS_ANNOTATIONS):
    """
    The piers and the groin, as seaward marks off the 1984 line.

    Positions come from HATTERAS_ANNOTATIONS. The piers' second field is a
    label height for the domain-axis figures and is ignored here; only the
    domain id is used. Drawn perpendicular to a shore-parallel pair of lines,
    so nothing here can be mistaken for a dune line despite sharing a colour
    family with one.
    """
    for name, spec in ann.piers.items():
        gid = spec[0] if isinstance(spec, (tuple, list)) else spec
        yy, x = _span_y(gdf, gid, gid), _x1984(by_dom, gid)
        if yy is None or x is None:
            continue
        yc = (yy[0] + yy[1]) / 2.0
        if not (y0 <= yc <= y1):
            continue
        ax.plot([x, x + STRUCTURE_LEN_M], [yc, yc], color=ann.color_pier,
                lw=1.8, solid_capstyle="butt", zorder=10)

    for name, gid in ann.groins.items():
        yy, x = _span_y(gdf, gid, gid), _x1984(by_dom, gid)
        if yy is None or x is None:
            continue
        yc = (yy[0] + yy[1]) / 2.0
        if not (y0 <= yc <= y1):
            continue
        ax.plot([x, x + STRUCTURE_LEN_M], [yc, yc], color=ann.color_groin,
                lw=1.8, solid_capstyle="butt", zorder=10)


def _km_axis(bar, y_origin, y0, y1):
    """
    Alongshore distance in km, on the bar strip's left edge.

    It goes on the BAR and not on the map on purpose: the map is pinned to
    equal aspect and its column width is set to its own data aspect, so hanging
    tick labels off it would letterbox it and break the row-for-row alignment
    that is the whole reason the bar sits beside it. The bar is not
    aspect-locked, and it shares the map's northing axis, so a ruler on the bar
    reads correctly against the map.

    Origin is the south end of domain 1 and distance increases north, matching
    fig_ribbon's x-axis.
    """
    lo = np.ceil((y0 - y_origin) / KM_TICK_M) * KM_TICK_M
    if lo <= 0.0:
        # The first panel starts ISLAND_PAD_M south of domain 1, so ceil() here
        # returns -0.0 and the tick formats as "-0". max(-0.0, 0.0) does not
        # fix it - the two compare equal and max keeps the first - so assign.
        lo = 0.0
    hi = (y1 - y_origin)
    vals = np.arange(lo, hi + 1.0, KM_TICK_M)
    bar.set_yticks([y_origin + v for v in vals])
    bar.set_yticklabels([f"{v / 1000:.0f}" for v in vals])
    bar.tick_params(axis="y", labelsize=7, length=2.5, pad=1.5,
                    color="0.55", labelcolor="0.35")


def _end_label(ax, text, at_top, y0, y1):
    """The site's own name for what lies off the end of the reach."""
    xa, xb = ax.get_xlim()
    # right of centre: the panel letter sits in the top-left corner of the
    # narrow island panels and a centred label reached back under it
    ax.text(xa + 0.58 * (xb - xa), y1 if at_top else y0, text,
            fontsize=8, fontstyle="italic", ha="center",
            va="top" if at_top else "bottom", color=INK_MUTED, zorder=12,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                      boxstyle="square,pad=0.22"))


def _scalebar(ax, length_m=SCALEBAR_M, show_cells=None):
    """The house scale bar (hat_figure_style._scalebar), with this module's
    cell size: under 1 km the label also says how many Barrier3D cells.

    `show_cells=False` drops that clause. The panels here are one to two
    inches wide on the page and the bar is a small fraction of one, so
    "500 m (50 cells)" is wider than the panel it sits in; the cell count is
    worth its width on the detail crops and not on the island maps.
    """
    _style._scalebar(ax, length_m=length_m, cell_m=CELL_M,
                     show_cells=show_cells)


# The per-domain statistic the bar strip can draw. Median is the figure's
# default and the number every other figure quotes; the mean is offered
# because a reader may ask for it, and the two disagree exactly where a
# domain's 1 m samples are skewed - a short reach of large offset inside an
# otherwise quiet domain pulls the mean and leaves the median alone.
ISLAND_STATS = {"median": "offset_med_m", "mean": "offset_mean_m"}


def fig_island_simple(elev, extent, gdf, lines, rows, reaches=None,
                      out=None, stat="median"):
    """
    The whole island, and the measured offset beside it.

    `stat` picks the per-domain statistic on the bar - see ISLAND_STATS. The
    default output name carries a suffix for anything but the median, so the
    two versions cannot overwrite each other.

    Each of the three columns is one third of the island as TWO axes sharing a
    northing axis:

      left   a true map at equal aspect - grey relief, both lines solid, the
             detail pairs boxed. It locates. It does NOT show the offset, and
             it cannot: 46 km of island against a 50 m offset is under half a
             line width, so the two lines lie on top of each other nearly
             everywhere on it. That is a property of the scale, not of the
             lines, and the bar beside it exists because of it.

      right  the per-domain median offset as a bar, aligned to the map row for
             row, red where 1984 lies seaward and blue where it lies landward.
             Same numbers the detail panels print, off the same CSV.

    Read together they answer the two halves of the question: the bar says
    where along the island the two lines disagree and by how much, the map says
    what that part of the island looks like and where the detail panels are cut
    from.
    """
    reaches = reaches or SIMPLE_REACHES
    ids = np.sort(gdf["domain_id"].astype(int).to_numpy())
    groups = np.array_split(ids, ISLAND_PANELS)
    by_dom = {r["domain"]: r for r in rows}
    if stat not in ISLAND_STATS:
        raise ValueError(f"stat must be one of {list(ISLAND_STATS)}, "
                         f"not {stat!r}")
    key = ISLAND_STATS[stat]
    med = np.array([r[key] for r in rows if r[key] != ""], float)
    lim = float(np.ceil(np.abs(med).max() / 10.0) * 10.0)

    # THE COLUMN WIDTHS ARE COMPUTED, NOT CHOSEN. Every axes in a one-row
    # figure gets the same height, and each map is pinned to equal aspect, so a
    # map whose column is wider than its own data aspect is letterboxed - it
    # shrinks vertically inside its box and stops lining up with the bar beside
    # it. Alignment row for row is the whole reason the bar sits next to the
    # map, so the width of each map column is set to exactly its own
    # x-span/y-span (the island is 5.7 km wide at the south end and 4.1 at the
    # north, so the three are genuinely different) and the maps then fill their
    # boxes and share a northing axis by construction.
    spans = []
    for group in groups:
        b = gdf[gdf["domain_id"].astype(int).isin(group)].total_bounds
        spans.append(((b[2] - b[0] + 2 * ISLAND_PAD_M),
                      (b[3] - b[1] + 2 * ISLAND_PAD_M)))
    ratios = [dx / dy for dx, dy in spans]
    strip = STRIP_WIDTH_RATIO * float(np.mean(ratios))
    widths = [v for r in ratios for v in (r, strip)]

    # The south end of domain 1, which is where the alongshore ruler starts.
    y_origin = float(gdf[gdf["domain_id"].astype(int) == int(ids.min())]
                     .total_bounds[1])

    apply_style()
    # Printed at the double-column width. The maps are aspect-locked, so the
    # panel height is whatever lets the six columns tile that width with
    # nothing letterboxed; the title row and the legend sit above and below.
    panel_h = (FIG_W_DOUBLE - 0.12 * 2 * ISLAND_PANELS - 0.45) / sum(widths)
    fig, axes = plt.subplots(
        1, 2 * ISLAND_PANELS,
        figsize=figsize("double", height=panel_h + 1.45),
        constrained_layout=True,
        gridspec_kw=dict(width_ratios=widths))

    for i, group in enumerate(groups):
        ax, bar = axes[2 * i], axes[2 * i + 1]
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        bx = sub.total_bounds
        y0, y1 = bx[1] - ISLAND_PAD_M, bx[3] + ISLAND_PAD_M

        # ---- the map -------------------------------------------------------
        _hillshade(ax, elev, extent, res=GRID_10M)
        sub.boundary.plot(ax=ax, color="0.45", linewidth=0.35, zorder=4)
        draw_lines(ax, lines, scale=LINE_SCALE_ISLAND,
                   style=SIMPLE_LINE_STYLE)
        for lo, hi, note in reaches:
            if not (group.min() <= lo <= group.max()):
                continue
            r = gdf[gdf["domain_id"].astype(int).isin(range(lo, hi + 1))]
            rb = r.total_bounds
            # No outline round the pair. At this scale the box is 2 km of a
            # 2 km-wide island, so it enclosed the whole width and read as a
            # feature of the island rather than a crop mark. The label and the
            # band on the bar beside it carry the same information without
            # drawing a rectangle over the only two lines the map has.
            ax.annotate(f"{lo}\u2013{hi}", xy=(rb[0], (rb[1] + rb[3]) / 2),
                        xytext=(-8, 0), textcoords="offset points",
                        fontsize=8, fontweight="bold", ha="right",
                        va="center", color=HIGHLIGHT, zorder=10,
                        bbox=dict(facecolor="white", alpha=0.85,
                                  edgecolor="none",
                                  boxstyle="square,pad=0.18"))
        ax.set_xlim(bx[0] - ISLAND_PAD_M, bx[2] + ISLAND_PAD_M)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        # The span alone: at the printed width the narrowest map column is
        # about 33 mm, and "Domains 61-90" centred over it runs back under
        # the panel letter. The caption says these are GIS domains.
        _title(ax, i, f"{group.min()}\u2013{group.max()}")
        _places(ax, gdf, y0, y1)
        _structures(ax, gdf, by_dom, y0, y1)
        if i == 0:
            _end_label(ax, HATTERAS_ANNOTATIONS.low_end_label, False, y0, y1)
            _north_arrow(ax, x=0.88, y=0.09)
        if i == ISLAND_PANELS - 1:
            _end_label(ax, HATTERAS_ANNOTATIONS.high_end_label, True, y0, y1)
            # One bar for the figure: the three maps share a northing span
            # and a height, so they are at one scale, and the south panel's
            # bottom corner is where the end label goes.
            _scalebar(ax, length_m=2000.0)

        # ---- the offset beside it ------------------------------------------
        for d in group:
            b = gdf[gdf["domain_id"].astype(int) == d].total_bounds
            v = by_dom[d][key]
            if v == "":
                continue
            bar.barh((b[1] + b[3]) / 2, v, height=(b[3] - b[1]) * 0.86,
                     color=OFFSET_POS if v > 0 else OFFSET_NEG,
                     linewidth=0, zorder=3)
        for lo, hi, note in reaches:
            if group.min() <= lo <= group.max():
                r = gdf[gdf["domain_id"].astype(int).isin(range(lo, hi + 1))]
                rb = r.total_bounds
                bar.axhspan(rb[1], rb[3], color=HIGHLIGHT, alpha=0.10,
                            zorder=1)
        bar.axvline(0, color=INK, lw=0.7, zorder=4)
        for c in (-CELL_M, CELL_M):
            bar.axvline(c, color=INK_MUTED, lw=0.6, ls=(0, (3, 2)), zorder=4)
        bar.set_xlim(-lim, lim)
        bar.set_ylim(y0, y1)
        _km_axis(bar, y_origin, y0, y1)
        bar.set_xticks([-lim, 0.0, lim])
        bar.tick_params(axis="x", labelsize=7)
        bar.set_title("Offset (m)", fontsize=9)
        if i == 0:
            bar.text(-0.30, 1.008, "km", transform=bar.transAxes, fontsize=7,
                     ha="center", va="bottom", color=INK_MUTED)
        bar.grid(axis="x", zorder=0)
        bar.set_axisbelow(True)
        for sp in ("top", "right", "left"):
            bar.spines[sp].set_visible(False)
        # Domain numbers live on the bar, not the map - on the map they would
        # sit on top of the two lines, which are the only thing it carries.
        # ... in the margin beyond the bar, not on it: at the printed width
        # a 70 m bar reaches the frame and the number landed on it.
        for d in group:
            if d % 10 == 0 or d in (ids.min(), ids.max()):
                b = gdf[gdf["domain_id"].astype(int) == d].total_bounds
                bar.text(lim * 1.08, (b[1] + b[3]) / 2, str(d), fontsize=7,
                         ha="left", va="center", color=INK_MUTED, zorder=6,
                         clip_on=False)

    # One legend for the whole figure, below the panels. On the map it would
    # have to sit on the island, which is 2 km wide at this scale.
    from matplotlib.patches import Patch
    fig.legend(
        handles=[Line2D([0], [0], label=f"{yr} dune line",
                        **dict(SIMPLE_LINE_STYLE[yr], linewidth=2.0))
                 for yr in (1984, 1997)]
        + [Patch(color=OFFSET_POS, label=f"{stat} offset, 1984 seaward"),
           Patch(color=OFFSET_NEG, label=f"{stat} offset, 1984 landward"),
           Line2D([0], [0], color=INK_MUTED, lw=0.8, ls=(0, (3, 2)),
                  label="\u00b11 Barrier3D cell (10 m)"),
           Line2D([0], [0], color=HATTERAS_ANNOTATIONS.color_town_span,
                  lw=5.0, label="community"),
           Line2D([0], [0], color=HATTERAS_ANNOTATIONS.color_pier, lw=1.8,
                  label="pier"),
           Line2D([0], [0], color=HATTERAS_ANNOTATIONS.color_groin, lw=1.8,
                  label="groin"),
           Patch(facecolor=HIGHLIGHT, alpha=0.10,
                 label="pair in the detail figure")],
        loc="outside lower center", ncol=3, frameon=False)

    suffix = "" if stat == "median" else f"_{stat}"
    q = Path(out) if out else (
        fig_path(f"HAT_duneline_offset_simple_island{suffix}.png"))
    q.parent.mkdir(parents=True, exist_ok=True)
    save(fig, q, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return q


def fig_island_lines(elev, extent, gdf, lines, rows, reaches=None,
                     out=None, per_panel=LINES_ISLAND_DOMAINS,
                     pad_m=LINES_ISLAND_PAD_M):
    """
    The whole island as maps only: no bar strips, the two lines the subject.

    fig_island_simple's map is a LOCATOR and its bar carries the offset,
    because at 15 km per panel the two lines sit inside one line width. This
    figure makes the map itself carry the offset by cutting the island into
    ~5 km panels and cropping each to the strip the lines occupy (see the
    LINES_ISLAND_* constants). Everything stays at equal aspect; it is a zoom,
    not an exaggeration. Communities, village ticks, structures and the detail
    pair labels are kept as location cues; the domain numbers, which lived on
    the bar, move to the landward edge of the map every fifth domain.
    """
    from shapely.geometry import box as _box
    from matplotlib.patches import Patch

    reaches = reaches or SIMPLE_REACHES
    ids = np.sort(gdf["domain_id"].astype(int).to_numpy())
    n_panels = int(np.ceil(len(ids) / per_panel))
    groups = np.array_split(ids, n_panels)
    by_dom = {r["domain"]: r for r in rows}
    med = np.array([r["offset_med_m"] for r in rows
                    if r["offset_med_m"] != ""], float)

    # Each panel's window: northing from the domain boxes, easting from the
    # lines themselves inside that northing window. The lines are what the
    # panel is for, so they set the crop; the boxes only say which domains.
    wins = []
    for group in groups:
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        bx = sub.total_bounds
        y0, y1 = bx[1] - 0.2 * pad_m, bx[3] + 0.2 * pad_m
        band = _box(bx[0] - 5000.0, y0, bx[2] + 5000.0, y1)
        xs = []
        for yr in lines:
            c = gpd.clip(lines[yr], band)
            if len(c):
                cb = c.total_bounds
                xs += [cb[0], cb[2]]
        if not xs:                      # no line here - fall back to the box
            xs = [bx[0], bx[2]]
        x0, x1 = min(xs) - pad_m, max(xs) + pad_m
        wins.append((x0, x1, y0, y1))
    ratios = [(w[1] - w[0]) / (w[3] - w[2]) for w in wins]

    apply_style()
    # Printed at the double-column width. Nine 1 km x 5 km panels in one row
    # of 190 mm are 15 mm wide and 75 mm tall, and a 50 m offset is under a
    # line width again - the thing this figure exists to avoid. Two rows
    # (five over four) give each panel ~2.3x the scale; up to four panels
    # stay in one row. Rows are subfigures because their width ratios differ.
    n_rows = 1 if n_panels <= 4 else 2
    per_row = int(np.ceil(n_panels / n_rows))
    row_ids = [list(range(k, min(k + per_row, n_panels)))
               for k in range(0, n_panels, per_row)]
    row_sum = max(sum(ratios[j] for j in rg) for rg in row_ids)
    # The panel height is whichever binds: the width of the page, or the page
    # itself. At ten domains a panel is 5 km by ~1 km, so it is the page
    # height that binds here and the panels sit in a wide row with margins.
    panel_h = min((FIG_W_DOUBLE - 0.12 * per_row - 0.4) / row_sum,
                  (FIG_H_MAX - 0.45) / n_rows - 0.75)
    fig = plt.figure(
        figsize=figsize("double", height=n_rows * (panel_h + 0.75) + 0.45),
        constrained_layout=True)
    axes = []
    for sf, rg in zip(np.atleast_1d(fig.subfigures(n_rows, 1)), row_ids):
        row_axes = sf.subplots(
            1, len(rg), gridspec_kw=dict(width_ratios=[ratios[j] for j in rg]))
        axes += list(np.atleast_1d(row_axes))
    # A panel is ~20 mm wide on the page, so nothing fits beside a panel
    # letter over it: the letter goes inside the top corner and the title is
    # the domain span alone. The caption says what the numbers are.
    head = "Domains {}\u2013{}" if panel_h * min(ratios) > 1.6 else "{}\u2013{}"
    inside = panel_h * min(ratios) <= 1.6

    for i, (group, win) in enumerate(zip(groups, wins)):
        ax = axes[i]
        x0, x1, y0, y1 = win
        sub = gdf[gdf["domain_id"].astype(int).isin(group)]
        _hillshade(ax, elev, extent, res=GRID_10M)
        sub.boundary.plot(ax=ax, color="0.45", linewidth=0.4, zorder=4)
        draw_lines(ax, lines, scale=LINES_ISLAND_LINE_SCALE,
                   style=SIMPLE_LINE_STYLE)
        for lo, hi, note in reaches:
            if not (group.min() <= lo <= group.max()):
                continue
            r = gdf[gdf["domain_id"].astype(int).isin(range(lo, hi + 1))]
            rb = r.total_bounds
            ax.axhspan(rb[1], rb[3], color=HIGHLIGHT, alpha=0.08, zorder=2)
            ax.annotate(f"{lo}\u2013{hi}", xy=(x0, (rb[1] + rb[3]) / 2),
                        xytext=(4, 0), textcoords="offset points",
                        fontsize=8, fontweight="bold", ha="left",
                        va="center", color=HIGHLIGHT, zorder=10,
                        bbox=dict(facecolor="white", alpha=0.85,
                                  edgecolor="none",
                                  boxstyle="square,pad=0.18"))
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        if inside:
            ax.set_title(head.format(group.min(), group.max()))
            _letter_inside(ax, i)
        else:
            _title(ax, i, head.format(group.min(), group.max()))
        _places(ax, gdf, y0, y1)
        _structures(ax, gdf, by_dom, y0, y1)
        if i == 0:
            _end_label(ax, HATTERAS_ANNOTATIONS.low_end_label, False, y0, y1)
            _north_arrow(ax, x=0.86, y=0.09)
        if i == n_panels - 1:
            _end_label(ax, HATTERAS_ANNOTATIONS.high_end_label, True, y0, y1)
            # One bar: every panel spans the same northing distance at the
            # same height, so they are all at one scale.
            _scalebar(ax, length_m=500.0, show_cells=False)
        # Domain numbers on the landward edge, every fifth. The bar that used
        # to carry them is gone, and on a 1 km-wide crop the landward margin
        # is backdune with nothing else drawn on it.
        for d in group:
            if d % 5 == 0 or d in (ids.min(), ids.max()):
                b = gdf[gdf["domain_id"].astype(int) == d].total_bounds
                # not in the top strip, where the panel letter sits
                if inside and (b[1] + b[3]) / 2 > y0 + 0.88 * (y1 - y0):
                    continue
                ax.text(x0 + 0.03 * (x1 - x0), (b[1] + b[3]) / 2, str(d),
                        fontsize=7, ha="left", va="center", color=INK_MUTED,
                        zorder=6,
                        bbox=dict(facecolor="white", alpha=0.7,
                                  edgecolor="none",
                                  boxstyle="square,pad=0.1"))

    fig.legend(
        handles=[Line2D([0], [0], label=f"{yr} dune line",
                        **dict(SIMPLE_LINE_STYLE[yr], linewidth=2.0))
                 for yr in (1984, 1997)]
        + [Line2D([0], [0], color=HATTERAS_ANNOTATIONS.color_town_span,
                  lw=5.0, label="community"),
           Line2D([0], [0], color=HATTERAS_ANNOTATIONS.color_pier, lw=1.8,
                  label="pier"),
           Line2D([0], [0], color=HATTERAS_ANNOTATIONS.color_groin, lw=1.8,
                  label="groin"),
           Patch(facecolor=HIGHLIGHT, alpha=0.08,
                 label="pair in the detail figure")],
        loc="outside lower center", ncol=6, frameon=False)

    q = Path(out) if out else fig_path("HAT_duneline_offset_lines_island.png")
    q.parent.mkdir(parents=True, exist_ok=True)
    save(fig, q, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return q


def fig_zooms_simple(gdf, lines, rows, elev=None, extent=None, reaches=None,
                     out=None, half_width=None, rows_by_domain=None, tall=False):
    """
    The two lines, and as little else as the picture can carry.

    Same measurement, same frame, same equal aspect and same two colours as
    fig_zooms. REMOVED: the terrain colour ramp and its colourbar, the
    elevation values, the coordinate ticks, and the dashed styling - both lines
    are SOLID here and differ only in colour. ADDED: a 1 m greyscale hillshade
    backdrop, and a scale bar in place of the ticks.

    What is NOT changed is the geometry. Equal aspect, no horizontal
    exaggeration, and the crop is the only thing making the offset visible -
    see SIMPLE_HALF_WIDTH_M for why the pairs are two domains rather than whole
    reaches. The hillshade's vertical exaggeration is a property of the
    BACKDROP only and moves nothing in the map plane.
    """
    apply_style()
    reaches = reaches or SIMPLE_REACHES
    half_width = SIMPLE_HALF_WIDTH_M if half_width is None else half_width
    by_dom = {r["domain"]: r for r in rows}

    # The panel width follows the crop: equal aspect, so a reach of eight
    # domains at +/-300 m is a much taller, narrower panel than a pair at
    # +/-150 m, and a fixed figsize would letterbox one or the other.
    spans = [gdf[gdf["domain_id"].astype(int).isin(range(lo, hi + 1))].total_bounds
             for lo, hi, _ in reaches]
    span_max = max(b[3] - b[1] for b in spans)
    # Printed at the double-column width; the height is whatever equal aspect
    # needs for the panels to tile that width (capped at a page, in which
    # case the reach panels are letterboxed and the gaps between them grow).
    n = len(reaches)
    h_axes = min((FIG_W_DOUBLE - 0.25 * (n + 1)) / (n * 2.0 * half_width / span_max),
                 FIG_H_MAX - 1.1)
    w_axes = h_axes * 2.0 * half_width / span_max
    # A panel narrower than about 40 mm cannot carry "Domains 62-68" centred
    # over it AND a letter beside it; the reach panels are 32 mm.
    head = "Domains {}\u2013{}" if w_axes > 1.5 else "{}\u2013{}"
    fig, axes = plt.subplots(1, n, figsize=figsize("double", height=h_axes + 1.1),
                             constrained_layout=True)
    span = 0.0
    for i, (ax, (lo, hi, note)) in enumerate(zip(np.atleast_1d(axes),
                                                 reaches)):
        ids = list(range(lo, hi + 1))
        sub = gdf[gdf["domain_id"].astype(int).isin(ids)]
        bx = sub.total_bounds
        span = bx[3] - bx[1]
        cen = float(np.median([by_dom[d]["x1984_med"] for d in ids
                               if by_dom[d]["x1984_med"] != ""]))

        # the backdrop covers the whole panel: with the reach panels padded to
        # a common height, the neighbouring domains' tiles are drawn too
        if tall:
            ymid = (bx[1] + bx[3]) / 2
            gb = gdf.bounds
            ids_draw = sorted(int(v) for v in gdf.loc[(gb["maxy"] > ymid - span_max / 2)
                                                        & (gb["miny"] < ymid + span_max / 2),
                                                        "domain_id"])
        else:
            ids_draw = ids
        arr, ext = load_1m(gdf, ids_draw)
        if arr is None:
            arr, ext = elev, extent
        _hillshade(ax, arr, ext)

        sub.boundary.plot(ax=ax, color="0.35", linewidth=0.8, zorder=4)
        draw_lines(ax, lines, scale=LINE_SCALE_DETAIL,
                   style=SIMPLE_LINE_STYLE)
        ax.set_xlim(cen - half_width, cen + half_width)
        if tall:
            # reaches of unequal length: every panel spans the longest one, centred,
            # so the panels come out the same height and their tops line up
            ymid = (bx[1] + bx[3]) / 2
            ax.set_ylim(ymid - span_max / 2, ymid + span_max / 2)
        else:
            ax.set_ylim(bx[1], bx[3])
        ax.set_aspect("equal")

        for _, r in sub.iterrows():
            b = r.geometry.bounds
            d = int(r["domain_id"])
            lab = f"domain {d}\n{by_dom[d]['offset_med_m']:+.0f} m"
            if rows_by_domain is not None:
                n = rows_by_domain.get(d, 0)
                lab += (f"\n{n:+d} row{'' if abs(n) == 1 else 's'}" if n else "\nno rows")
            ax.text(cen - half_width + 12, (b[1] + b[3]) / 2, lab,
                    fontsize=8, ha="left", va="center", color=INK,
                    zorder=8, linespacing=1.35,
                    bbox=dict(facecolor="white", alpha=0.8,
                              edgecolor="none", boxstyle="square,pad=0.25"))

        # The span alone, on ONE line: a second line raised that panel's
        # title above its neighbours' and the row of letters stopped lining
        # up. A note is drawn only where there is one panel to carry it (a
        # --zoom-note); "control" is caption text, like the reading itself
        # (_pair_reading in write_captions).
        _title(ax, i, head.format(lo, hi)
               + (f" \u00b7 {note}" if note and n == 1 else ""))
        if i == 0:
            # One bar and one arrow for the figure: every panel is the same
            # crop width over the same alongshore span, so they share a
            # scale, and the first panel's bottom corner is the only one
            # with no domain label in it.
            _scalebar(ax, show_cells=not tall)
            _north_arrow(ax, x=0.88, y=0.06)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.legend(
        handles=[Line2D([0], [0], label=f"{yr} dune line",
                        **SIMPLE_LINE_STYLE[yr]) for yr in (1984, 1997)],
        loc="outside lower center", ncol=2, frameon=False)

    p = Path(out) if out else fig_path("HAT_duneline_offset_simple.png")
    p.parent.mkdir(parents=True, exist_ok=True)
    save(fig, p, vector=False, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_by_domain(rows):
    """
    The requested number, domain by domain, with its spread.

    One series, so the marker carries the SIGN rather than an identity: red
    where the 1984 line lies seaward, blue where it lies landward, the same
    pair every other figure in the folder uses for the same fact. The
    communities are banded along the axis so a reader can place a domain
    without the map.
    """
    apply_style()
    dom = np.array([r["domain"] for r in rows])
    med = np.array([r["offset_med_m"] if r["offset_med_m"] != "" else np.nan
                    for r in rows], float)
    p25 = np.array([r["offset_p25_m"] if r["offset_p25_m"] != "" else np.nan
                    for r in rows], float)
    p75 = np.array([r["offset_p75_m"] if r["offset_p75_m"] != "" else np.nan
                    for r in rows], float)
    fin = np.isfinite(med)

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.42),
                           constrained_layout=True)
    ax.vlines(dom, p25, p75, color="0.72", linewidth=2.0, zorder=2,
              label="interquartile range within the domain")
    pos, neg = fin & (med >= 0), fin & (med < 0)
    ax.scatter(dom[pos], med[pos], s=16, color=C_1984, zorder=3,
               edgecolor="white", linewidth=0.5,
               label="median offset, 1984 seaward")
    ax.scatter(dom[neg], med[neg], s=16, color=C_1997, zorder=3,
               edgecolor="white", linewidth=0.5,
               label="median offset, 1984 landward")
    ax.axhline(0, color=INK, linewidth=0.7, zorder=1)
    for s_ in (-CELL_M, CELL_M):
        ax.axhline(s_, color=INK_MUTED, linewidth=0.6, linestyle=(0, (3, 2)),
                   zorder=1, label="\u00b11 Barrier3D cell (10 m)")

    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("1984 minus 1997 dune line (m)\nseaward positive")
    ax.set_xlim(dom.min() - 0.5, dom.max() + 0.5)
    # headroom, so the village names along the top clear the tallest bars
    lo_, hi_ = ax.get_ylim()
    ax.set_ylim(lo_, hi_ + 0.14 * (hi_ - lo_))
    town_bands(ax)
    h, l = ax.get_legend_handles_labels()
    fig.legend(h[:4], l[:4], loc="outside lower center", ncol=2, frameon=False)
    ax.set_xticks([1] + list(range(10, int(dom.max()) + 1, 10)))
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)
    open_frame(ax)

    p = fig_path("HAT_duneline_offset_bydomain.png")
    save(fig, p, bbox_inches="tight")
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
        for k in ("offset_med_m", "offset_mean_m", "offset_sd_m",
                  "x1984_med"):
            rec[k] = float(r[k]) if r[k] not in ("", "nan") else ""
        out.append(rec)
    return out


def read_insert_rows(path):
    """{domain: N} from the row-insert scope table, if it has been built."""
    if not Path(path).exists():
        print(f"  NOTE: {Path(path).name} absent - row counts not labelled")
        return None
    return {int(r["domain"]): int(r["n_cells"]) for r in
            csv.DictReader(open(path))}


def simple_only(out, half_width, span=None, island_out=None,
                stat="median"):
    """
    Render the two simple figures alone, off the existing table.

    Reads duneline_offset_by_domain.csv rather than re-measuring, so this is
    seconds. `island_out=False` skips the locator.
    """
    rows = read_rows(OUT_DIR / CSV_NAME)
    gdf = load_domains()
    print("dune lines:")
    lines = load_lines(gdf.crs)
    print("\nfor drawing:")
    drawn = clip_for_drawing(lines, gdf.union_all())

    reaches = None
    if span:
        lo, hi = (int(x) for x in span.split("-"))
        reaches = [(lo, hi, "")]

    # The 10 m mosaic is what the locator draws, so it is loaded unconditionally
    # here - unlike in the detail panels, where it is only a fallback for a
    # missing 1 m tile.
    print(f"\nloading the {SOURCE_TAG} mosaic at 10 m ...")
    elev, _s, extent, _n = m.load_mosaic()

    figs = [fig_zooms_simple(gdf, drawn, rows, elev=elev, extent=extent,
                             reaches=reaches, out=out,
                             half_width=half_width)]
    if island_out is not False:
        figs.append(fig_island_simple(elev, extent, gdf, drawn, rows,
                                      reaches=reaches, out=island_out,
                                      stat=stat))
    figs.append(write_captions(rows, simple_half_width=half_width))
    print()
    for f in figs:
        print(f"  figure : {f}")
    return figs


def lines_island_only(out, per_panel, pad_m):
    """Render the lines-only island figure alone, off the existing table."""
    rows = read_rows(OUT_DIR / CSV_NAME)
    gdf = load_domains()
    print("dune lines:")
    lines = load_lines(gdf.crs)
    print("\nfor drawing:")
    drawn = clip_for_drawing(lines, gdf.union_all())
    print(f"\nloading the {SOURCE_TAG} mosaic at 10 m ...")
    elev, _s, extent, _n = m.load_mosaic()
    f = fig_island_lines(elev, extent, gdf, drawn, rows, out=out,
                         per_panel=per_panel, pad_m=pad_m)
    print(f"\n  figure : {f}")
    return f


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
    # Only a note given by hand goes under the panel title. The offset range
    # used to be generated into it; that is a statistics line, and the
    # per-domain labels inside the panel already carry the numbers.
    p = fig_zooms(elev, extent, gdf, drawn, rows,
                  reaches=[(lo, hi, f"{lo}-{hi}", note or "")],
                  out=out, half_width=half_width,
                  rows_by_domain=n_by)
    print(f"\n  figure : {p}  (offset {min(meds):+.0f} to {max(meds):+.0f} m)")
    return p


# =============================================================================
# MAIN
# =============================================================================

def write_captions(rows, half_width=None, simple_half_width=None):
    """
    A caption per figure, with the numbers filled from the table.

    The figures carry no title sentences or footnote paragraphs - that text
    belongs under the figure in whatever document uses it, and keeping it here
    rather than on the image means it is editable, searchable and cannot go
    stale against the picture without the file's own numbers going stale too.
    """
    half_width = ZOOM_HALF_WIDTH_M if half_width is None else half_width
    shw = SIMPLE_HALF_WIDTH_M if simple_half_width is None else simple_half_width
    med = np.array([r["offset_med_m"] for r in rows
                    if r["offset_med_m"] != ""], float)
    n_cell = int((np.abs(med) >= CELL_M).sum())
    mean = np.array([r["offset_mean_m"] for r in rows
                     if r.get("offset_mean_m", "") != ""], float)
    mean_all = float(np.mean(mean)) if mean.size else float("nan")
    mean_lo = float(mean.min()) if mean.size else float("nan")
    mean_hi = float(mean.max()) if mean.size else float("nan")
    n_pos = int((med > 0).sum())
    stats = (f"Per-domain median offset (1984 minus 1997, seaward positive), "
             f"from `{CSV_NAME}`: island median {np.median(med):+.1f} m, "
             f"range {med.min():+.1f} to {med.max():+.1f} m; {n_cell} of "
             f"{len(med)} domains differ by at least one 10 m Barrier3D "
             f"cell, {n_pos} of {len(med)} are positive.")
    by = {r["domain"]: r for r in rows}

    def reading(lo, hi):
        ms = [by[d]["offset_med_m"] for d in range(lo, hi + 1)
              if d in by and by[d]["offset_med_m"] != ""]
        return _pair_reading(ms) if ms else "no measurement"

    pairs = "; ".join(
        f"{lo}\u2013{hi} ({_place_of(lo, hi)}{'; the control' if note else ''}): "
        f"{reading(lo, hi)}" for lo, hi, note in SIMPLE_REACHES)
    reaches = "; ".join(f"{slug} ({note}): {reading(lo, hi)}"
                        for lo, hi, slug, note in ZOOM_REACHES)
    src = ("Communities, village centres, piers and the groin are the "
           "project's own positions (`hatteras_site_config."
           "HATTERAS_ANNOTATIONS`); structures are drawn seaward off the "
           f"1984 line at a fixed {STRUCTURE_LEN_M:.0f} m, a mark rather than "
           "a surveyed extent.")
    caps = [
        ("HAT_duneline_offset_simple.png",
         f"The 1984 (red) and 1997 (blue) dune lines at true scale on four "
         f"two-domain pairs, one of them a control where the two agree; "
         f"reading south to north, with the median offset on each pair: "
         f"{pairs}. Equal aspect, nothing exaggerated in the map plane: each "
         f"panel is {2 * shw:.0f} m cross-shore, centred on the 1984 line, "
         f"by two 500 m domains alongshore. Grey relief is the 1 m gap-filled "
         f"DEM shaded at {HILLSHADE['vert_exag']:.1f}\u00d7 vertical "
         f"exaggeration; the shading carries no readable elevation. Numbers "
         f"beside each domain are its median offset, positive where 1984 "
         f"lies seaward. Scale bar 50 m = 5 Barrier3D cells."),
        ("HAT_duneline_offset_simple_island.png",
         f"Where the two dune lines are and where they disagree, over the "
         f"whole island in three equal-aspect panels (south at left). Beside "
         f"each map, the per-domain median offset as a bar, aligned row for "
         f"row: red where the 1984 line lies seaward of 1997, blue where it "
         f"lies landward; dashed guides at \u00b110 m, one Barrier3D cell. "
         f"At this scale the two lines coincide within a line width nearly "
         f"everywhere, which is what the bars are for. The ruler beside each "
         f"bar is km north of the south end of domain 1, the same origin the "
         f"ribbon figure uses. Shaded bands mark the pairs shown in the "
         f"detail figure. {stats} {src}"),
        ("HAT_duneline_offset_simple_island_mean.png",
         f"As the previous figure, with the per-domain MEAN offset on the "
         f"bars in place of the median (island mean {mean_all:+.1f} m, "
         f"range {mean_lo:+.1f} to {mean_hi:+.1f} m). The two differ where "
         f"a domain's 1 m samples are skewed: a short stretch of large "
         f"offset inside an otherwise quiet domain moves the mean and not "
         f"the median."),
        ("HAT_duneline_offset_lines_island.png",
         f"The two dune lines over the whole island as maps only, in "
         f"{LINES_ISLAND_DOMAINS}-domain (~5 km) panels reading south to "
         f"north, left to right and then the second row; each panel is "
         f"titled with the domains it spans. Each panel is cropped at equal aspect to the "
         f"envelope of the two lines plus {LINES_ISLAND_PAD_M:.0f} m either "
         f"side, so the separation is visible on the map itself without "
         f"exaggeration. Domain numbers every fifth domain on the landward "
         f"edge; shaded bands mark the pairs in the detail figure. {src}"),
        ("HAT_duneline_offset_lines_island_3panel.png",
         f"As the previous figure, in three 30-domain panels matching the "
         f"map-and-bar figure's layout. At 15 km per panel a 50 m offset is "
         f"about one line width, so the separation reads only where it is "
         f"large. {src}"),
        ("HAT_duneline_offset_ribbon.png",
         f"(a) The 1984 and 1997 dune lines along the island at 1 m "
         f"alongshore sampling, each drawn relative to their common "
         f"{BASELINE_WINDOW_M / 1000:.0f} km boxcar-smoothed midline so the "
         f"island's curvature drops out, seaward positive; the band between "
         f"them is filled red where 1984 lies seaward and blue where it lies "
         f"landward. (b) The difference, 1984 minus 1997, filled by sign; "
         f"dashed guides at \u00b110 m, one Barrier3D cell. The alongshore "
         f"axis is the GIS domain, 1 at Cape Point and 90 at Pea Island, each "
         f"1 m sample placed within its own 500 m box; the axis along the top "
         f"is distance in km from the south end of domain 1, the origin the "
         f"island figure's ruler uses. Grey bands are the communities; the "
         f"brackets over (b) mark the reaches shown in the true-scale reach "
         f"figure ({reaches}). {stats}"),
        ("HAT_duneline_offset_zooms.png",
         f"The two dune lines at true scale on three reaches of five to eight "
         f"domains, south to north, with the median offset on each: "
         f"{reaches}. Equal aspect; each panel is cropped to "
         f"{half_width:.0f} m either side of the local 1984 line rather than "
         f"the full 2000 m domain box. Grey relief is the 1 m gap-filled DEM "
         f"shaded at {HILLSHADE['vert_exag']:.1f}\u00d7 vertical exaggeration "
         f"and carries no readable elevation; the domain boxes are outlined; "
         f"the number beside each domain is its median offset, positive where "
         f"1984 lies seaward. Scale bar 50 m = 5 Barrier3D cells."),
        ("HAT_duneline_offset_zoom_83_87.png",
         f"As the previous figure, for domains 83\u201387: GIS 85 and its "
         f"neighbours, the largest sustained positive run on the island. "
         f"Where the 1984 footprint table is on disk the label also gives "
         f"the number of Barrier3D rows the 1984 start would add (+) or "
         f"remove (\u2212) there, trunc(paired shift / 10 m): a row only once a "
         f"full cell of change is measured, in either direction."),
        ("HAT_duneline_offset_bydomain.png",
         f"Median cross-shore offset between the 1984 and 1997 dune lines in "
         f"each of the 90 domains, 1984 minus 1997 with seaward positive, "
         f"with the interquartile range of the 1 m samples within the "
         f"domain. Red markers: 1984 seaward; blue: 1984 landward. Dashed "
         f"guides at \u00b110 m, one Barrier3D cell. Domain 1 is at Cape "
         f"Point, 90 at Pea Island; grey bands are the communities. {stats}"),
    ]
    q = FIG_DIR / "CAPTIONS.md"
    with open(q, "w", encoding="utf-8") as fh:
        fh.write("# Figure captions\n\n")
        fh.write("Written by `HAT_plot_duneline_offset.py` "
                 "(`write_captions`) from the same table the figures draw "
                 "from. Each heading names the figure's folder under "
                 "`figures/` (island/, detail/, offset/). The figures carry "
                 "no in-image titles or footnotes on purpose; use these "
                 "under them.\n\n")
        for name, text in caps:
            fh.write(f"## `{FIG_SUBFOLDER[name]}/{name}`\n\n{text}\n\n")
    return q


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
            fig_zooms_simple(gdf, drawn, rows, elev=elev, extent=extent),
            fig_island_simple(elev, extent, gdf, drawn, rows),
            fig_island_simple(elev, extent, gdf, drawn, rows, stat="mean"),
            fig_island_lines(elev, extent, gdf, drawn, rows),
            fig_island_lines(elev, extent, gdf, drawn, rows, per_panel=30,
                             out=fig_path("HAT_duneline_offset_lines_island_3panel.png")),
            fig_by_domain(rows)]
    figs.append(write_captions(rows))
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
    ap.add_argument("--simple", action="store_true",
                    help="render ONLY the simplified figure (grey relief, no "
                         "elevation values, both lines solid, tight crop) and "
                         "exit. Reads the existing table rather than "
                         "re-measuring.")
    ap.add_argument("--simple-out", default=None,
                    help="output path for the --simple detail figure")
    ap.add_argument("--simple-island-out", default=None,
                    help="output path for the --simple island locator")
    ap.add_argument("--no-island", action="store_true",
                    help="with --simple, render only the detail panels")
    ap.add_argument("--simple-span", default=None, metavar="LO-HI",
                    help="render the simple figure for ONE domain span "
                         "instead of the three standard pairs")
    ap.add_argument("--simple-stat", choices=sorted(ISLAND_STATS),
                    default="median",
                    help="per-domain statistic on the --simple island "
                         "locator's bar strip (default median)")
    ap.add_argument("--simple-halfwidth", type=float,
                    default=SIMPLE_HALF_WIDTH_M,
                    help=f"cross-shore half-width of a simple panel, m "
                         f"(default {SIMPLE_HALF_WIDTH_M:.0f})")
    ap.add_argument("--lines-island", action="store_true",
                    help="render ONLY the lines-only whole-island figure (no "
                         "bar strips, ~5 km panels cropped to the lines) and "
                         "exit. Reads the existing table.")
    ap.add_argument("--lines-island-out", default=None,
                    help="output path for --lines-island")
    ap.add_argument("--lines-per-panel", type=int,
                    default=LINES_ISLAND_DOMAINS,
                    help=f"domains per panel for --lines-island (default "
                         f"{LINES_ISLAND_DOMAINS})")
    ap.add_argument("--lines-pad", type=float, default=LINES_ISLAND_PAD_M,
                    help=f"easting pad either side of the lines' envelope, m "
                         f"(default {LINES_ISLAND_PAD_M:.0f})")
    ap.add_argument("--captions", action="store_true",
                    help="write figures/CAPTIONS.md from the existing table "
                         "and exit")
    ap.add_argument("--no-row-labels", action="store_true",
                    help="do not annotate each domain with its footprint row "
                         "count (+ added / - removed), even if the table is on disk")
    a = ap.parse_args()

    if a.captions:
        print(f"  captions : {write_captions(read_rows(OUT_DIR / CSV_NAME))}")
    elif a.lines_island:
        lines_island_only(a.lines_island_out, a.lines_per_panel,
                          a.lines_pad)
    elif a.simple:
        simple_only(a.simple_out, a.simple_halfwidth, a.simple_span,
                    island_out=False if a.no_island else a.simple_island_out,
                    stat=a.simple_stat)
    elif a.zoom:
        zoom_only(a.zoom, a.zoom_out, a.zoom_halfwidth, a.zoom_note,
                  not a.no_row_labels)
    else:
        main()
