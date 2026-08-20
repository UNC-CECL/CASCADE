r"""
HAT_road_geojson_map.py
===============================================================================
The 1984 and 2004 NC-12 geojsons drawn on the 2009 DEM, in MAP coordinates, for
the whole modelled island. A reference sheet: find a GIS domain on the real
island, and see where each road line actually runs across it.

  data/hatteras_init/4-mgmt-forcing/road_offset/raster/
      HAT_road_geojson_on_2009_dem.png

WHY THIS IS A MAP AND NOT A DOMAIN-FRAME FIGURE
-----------------------------------------------
Every other road figure in this tree lives in the Barrier3D domain frame, which
is reached through orient -> alongshore flip -> shear -> water trim. That chain
is exactly what makes those figures comparable to the model, and exactly what
makes them useless for checking the model's inputs against the world.

This figure applies NO transform. The DEM is read in its native UTM 18N grid and
the geojson is reprojected onto it with `to_crs`, the same single step
HAT_check_geojson_vs_mask.py uses for the same reason. If the road line and the
road visible in the LiDAR agree here, the source data is right; whether the
DOMAIN arrays are right is a separate question that the domain-frame figures
answer.

WHY THE STRIPS RUN LEFT TO RIGHT
---------------------------------
Hatteras runs very nearly north-south: GIS 1-90 spans 8.0 km east-west and 45.2
km north-south, 5.7:1. Cut into six 15-domain segments at TRUE NORTH and TRUE
SCALE, every segment is portrait (h/w 1.7-2.7) -- the island is simply taller
than it is wide at any segmentation. Stacking portrait panels vertically would
give a figure four feet tall, so the six run SOUTH (left) to NORTH (right)
instead. North is up in every panel and the scale is true in both axes; only the
reading order is unusual, and the panel titles carry it.

WHAT THE COLOURS MEAN, AND THE ONE DISTINCTION ONLY THIS FIGURE CAN DRAW
------------------------------------------------------------------------
The source DEM stores NoData as exactly -10.0 m NAVD88, and everywhere else in
this project that value has already been folded into the water sentinel --
Barrier3D has no representation for "unknown", so by the time the topography is
saved, a LiDAR hole and a genuinely wet cell are the same number. This figure
reads the RAW tif, so it is the one place the two can be told apart, and it
draws them differently: never surveyed is a distinct grey from water.

That matters for reading GIS 78-80, where the roadway relocation fires. Those
domains drown on coverage gaps, not on measured water, and here you can see it.

REQUIREMENTS
------------
  numpy, matplotlib, rasterio, geopandas
===============================================================================
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.patheffects import withStroke

import geopandas as gpd
import rasterio

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[5]
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"

DEM_DIR = (INIT_ROOT / "1-barrier3d-domains" / "2009-raw"
           / "2009-domain-clipresample")
DEM_NAME = "domain_{d}/clip_domain_{d}.tif"          # 1 m, the full source grid

ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"
GEOJSON_FMT = ROADS_ROOT / "raw_offset" / "{year}" / "nc12_{year}.geojson"
OUT_PNG = ROADS_ROOT / "raster" / "HAT_road_geojson_on_2009_dem.png"

EXTRACTOR = (PROJECT_ROOT / "scripts" / "input_prep" / "1-barrier3d-domains"
             / "topography_dunes" / "HAT_dune_topo_extractor.py")
PLACEMENT = (PROJECT_ROOT / "scripts" / "input_prep" / "4-mgmt-forcings"
             / "road_offset" / "1-produce" / "HAT_road_placement_on_domains.py")

DOMAINS = list(range(1, 91))
YEARS = (1984, 2004)
STRIP_SIZE = 15                 # 6 strips; see WHY THE STRIPS RUN LEFT TO RIGHT
PAD_M = 150.0                   # breathing room around each strip's bbox
READ_STRIDE = 3                 # decimated read, ~3 m -- the display resolution
LABEL_EVERY = 5                 # GIS number on every Nth domain box

DPI = 300                       # 3 m/px on the page, so the road is resolvable


def constant_from_extractor(name: str, fallback: float) -> float:
    """
    Read one numeric constant out of the extractor without importing it.

    Importing a 2000-line module for two floats risks its import-time side
    effects; copying the numbers risks them drifting apart. Parsing the
    assignment is neither.
    """
    try:
        src = EXTRACTOR.read_text(encoding="utf-8", errors="replace")
        m = re.search(rf"^{name}\s*=\s*(-?\d+(?:\.\d+)?)", src, re.MULTILINE)
        if m:
            return float(m.group(1))
    except OSError:
        pass
    print(f"[warn] could not read {name} from the extractor; using {fallback}")
    return fallback


MHW_M = constant_from_extractor("MHW_M", 0.36)
RAW_NODATA_MAX_NAVD = constant_from_extractor("RAW_NODATA_MAX_NAVD", -9.0)


def load_placement():
    spec = importlib.util.spec_from_file_location("hat_placement", PLACEMENT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


P = load_placement()                       # palette, sections, rcParams

# Never-surveyed is the PALEST thing on the sheet, deliberately. Most of each
# 2000 m clip box is off-island and therefore NoData, so at any real saturation
# it becomes the largest block of colour in the figure and the island reads as
# mostly unsurveyed. It is absence of information and should recede; the holes
# that matter sit INSIDE the island, where pale against green still reads.
C_NODATA = "#e4e0d9"
C_WATER = "#c3ccd4"         # at or below 0 m MHW
C_BOX = "#5a5a55"

# The two lines coincide except in the relocation blocks, so drawing them at the
# same weight hides 1984 completely. 1984 goes underneath and wide, 2004 dashed
# on top: coincident stretches read as a dashed orange line on a blue casing,
# and the places they genuinely diverge are the places you see two lines.
LW_1984, LW_2004 = 2.8, 1.3


# =============================================================================
# DATA
# =============================================================================

def domain_path(d: int) -> Path:
    return DEM_DIR / DEM_NAME.format(d=d)


def read_bounds() -> dict:
    out = {}
    for d in DOMAINS:
        p = domain_path(d)
        if not p.is_file():
            continue
        with rasterio.open(p) as src:
            out[d] = src.bounds
    if not out:
        raise SystemExit(f"no clip tifs found under {DEM_DIR}")
    return out


def read_domain(d: int):
    """
    One domain at full 1 m, block-reduced to the display resolution.

    Reduced here rather than by rasterio's decimated read, because the two
    quantities need OPPOSITE reductions and rasterio can only apply one per
    read (and refuses `Resampling.min` on reads at all):

        elevation   mean of the SURVEYED cells in the block
        no-data     any() -- a block containing any hole is a hole

    Averaging a hole together with real ground would invent an elevation and
    quietly shrink the coverage gaps, which are the thing this figure exists to
    show honestly. NoData is identified on the RAW NAVD88 values before the
    datum shift, the same order the extractor uses.
    """
    with rasterio.open(domain_path(d)) as src:
        raw = src.read(1).astype(float)
        b, (rx, ry), nodata = src.bounds, src.res, src.nodata

    hole_full = raw <= RAW_NODATA_MAX_NAVD
    if nodata is not None:
        hole_full |= np.isclose(raw, nodata)

    s = max(1, READ_STRIDE)
    h, w = (raw.shape[0] // s) * s, (raw.shape[1] // s) * s
    rb = raw[:h, :w].reshape(h // s, s, w // s, s)
    hb = hole_full[:h, :w].reshape(h // s, s, w // s, s)

    hole = hb.any(axis=(1, 3))
    good = ~hb
    n_good = good.sum(axis=(1, 3))
    total = np.where(good, rb, 0.0).sum(axis=(1, 3))
    elev = np.divide(total, n_good, out=np.full(hole.shape, np.nan),
                     where=n_good > 0) - MHW_M

    # Extent from the TRIMMED size, not the file's, so the image lands on the
    # ground it was read from -- the trim discards up to s-1 rows/cols.
    ext = (b.left, b.left + w * rx, b.top - h * ry, b.top)
    return elev, hole, ext


def load_roads(target_crs) -> dict:
    """Both geojson lines, reprojected onto the DEM's own grid."""
    out = {}
    for year in YEARS:
        p = Path(str(GEOJSON_FMT).format(year=year))
        if not p.is_file():
            print(f"  [warn] no geojson for {year}: {p}")
            continue
        gdf = gpd.read_file(p)
        try:
            out[year] = gdf.to_crs(target_crs)
        except Exception:
            # The tifs carry a COMPOUND CRS (UTM 18N + NAVD88 height). If pyproj
            # declines to transform onto it, fall back to its horizontal part --
            # the vertical component is irrelevant to a plan-view reprojection.
            out[year] = gdf.to_crs("EPSG:26918")
            print(f"  [note] {year}: reprojected to EPSG:26918 "
                  f"(horizontal part of the compound CRS)")
    return out


# =============================================================================
# FIGURE
# =============================================================================

def strips() -> list:
    return [DOMAINS[i:i + STRIP_SIZE]
            for i in range(0, len(DOMAINS), STRIP_SIZE)]


def draw_strip(ax, group, bounds, roads):
    xs = [v for d in group if d in bounds
          for v in (bounds[d].left, bounds[d].right)]
    ys = [v for d in group if d in bounds
          for v in (bounds[d].bottom, bounds[d].top)]
    x0, x1 = min(xs) - PAD_M, max(xs) + PAD_M
    y0, y1 = min(ys) - PAD_M, max(ys) + PAD_M

    ax.set_facecolor(C_WATER)
    for d in group:
        if d not in bounds:
            continue
        elev, hole, ext = read_domain(d)

        # Three states, drawn as three layers rather than one ramp: water and
        # "never surveyed" are states, not small elevations, and the whole point
        # of reading the raw tif is that they can be separated here.
        ax.imshow(np.ma.masked_where(~hole, np.ones_like(elev)),
                  extent=ext, origin="upper", aspect="equal",
                  cmap=matplotlib.colors.ListedColormap([C_NODATA]),
                  vmin=0, vmax=1, interpolation="nearest", zorder=2)
        wet = (~hole) & (elev <= 0.0)
        ax.imshow(np.ma.masked_where(~wet, np.ones_like(elev)),
                  extent=ext, origin="upper", aspect="equal",
                  cmap=matplotlib.colors.ListedColormap([C_WATER]),
                  vmin=0, vmax=1, interpolation="nearest", zorder=2)
        land = np.ma.masked_where(hole | (elev <= 0.0), elev)
        im = ax.imshow(land, extent=ext, origin="upper", aspect="equal",
                       cmap=P.LAND_CMAP, norm=Normalize(0.0, 4.0),
                       interpolation="nearest", zorder=3)

        b = bounds[d]
        ax.add_patch(Rectangle((b.left, b.bottom), b.right - b.left,
                               b.top - b.bottom, fill=False, ec=C_BOX,
                               lw=0.45, alpha=0.55, zorder=6))
        if d % LABEL_EVERY == 0 or d == group[0] or d == group[-1]:
            ax.text(b.left + 60, (b.bottom + b.top) / 2, str(d), ha="left",
                    va="center", fontsize=6.5, color="#22221f", zorder=8,
                    path_effects=[withStroke(linewidth=2.0,
                                             foreground=P.SURFACE)])

    for year, gdf in sorted(roads.items()):
        clipped = gdf.clip_by_rect(x0, y0, x1, y1)
        clipped = clipped[~clipped.is_empty]
        if clipped.empty:
            continue
        first = year == YEARS[0]
        gpd.GeoSeries(clipped).plot(
            ax=ax, color=P.C_YEAR[year],
            linewidth=LW_1984 if first else LW_2004,
            linestyle="-" if first else (0, (3.2, 1.8)),
            zorder=7 if first else 8)

    # No section names on this figure. The strips are only ~3 km wide, so a
    # boxed label sits on the island rather than beside it and hides the DEM
    # and the road lines underneath. The GIS numbers already locate a domain,
    # and the domain-frame figures carry the section bands.
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_edgecolor(P.INK_MUTED)
    ax.set_title(f"GIS {group[0]}–{group[-1]}", loc="left", fontsize=9.5,
                 color=P.INK_SECOND, pad=3)
    return im


def scale_bar(ax, length_m=1000.0):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    xa = x0 + 0.06 * (x1 - x0)
    ya = y0 + 0.035 * (y1 - y0)
    ax.plot([xa, xa + length_m], [ya, ya], color="#0b0b0b", lw=2.4, zorder=10,
            solid_capstyle="butt",
            path_effects=[withStroke(linewidth=4.4, foreground=P.SURFACE)])
    ax.text(xa + length_m / 2, ya + 0.012 * (y1 - y0),
            f"{length_m / 1000:.0f} km", ha="center", va="bottom", fontsize=7,
            color="#0b0b0b", zorder=10,
            path_effects=[withStroke(linewidth=2.2, foreground=P.SURFACE)])


def north_arrow(ax):
    ax.annotate("N", xy=(0.5, 0.965), xytext=(0.5, 0.90),
                xycoords="axes fraction", textcoords="axes fraction",
                ha="center", va="bottom", fontsize=8.5, color="#0b0b0b",
                arrowprops=dict(arrowstyle="-|>", color="#0b0b0b", lw=1.3),
                zorder=11)


def main() -> int:
    print("=" * 84)
    print("NC-12 geojson on the 2009 DEM -- map frame, no transforms")
    print("=" * 84)

    bounds = read_bounds()
    print(f"  domains with a 1 m clip tif : {len(bounds)} of {len(DOMAINS)}")

    with rasterio.open(domain_path(sorted(bounds)[0])) as src:
        crs = src.crs
    roads = load_roads(crs)
    for year, gdf in sorted(roads.items()):
        print(f"  {year} geojson             : {len(gdf)} feature(s), "
              f"reprojected to the DEM grid")

    groups = strips()
    heights, widths = [], []
    for g in groups:
        xs = [v for d in g if d in bounds
              for v in (bounds[d].left, bounds[d].right)]
        ys = [v for d in g if d in bounds
              for v in (bounds[d].bottom, bounds[d].top)]
        widths.append(max(xs) - min(xs) + 2 * PAD_M)
        heights.append(max(ys) - min(ys) + 2 * PAD_M)

    # True scale in both axes means the panel widths are set by the data, not
    # chosen. Width ratios come from each strip's own bbox.
    fig_w = 22.0
    fig_h = fig_w * max(heights) / sum(widths) * 1.16
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(1, len(groups), width_ratios=widths, wspace=0.035,
                          left=0.022, right=0.93, top=0.855, bottom=0.035)

    im = None
    for k, g in enumerate(groups):
        ax = fig.add_subplot(gs[k])
        im = draw_strip(ax, g, bounds, roads) or im
        if k == 0:
            scale_bar(ax)
            north_arrow(ax)
        print(f"  strip GIS {g[0]:>2}-{g[-1]:<2} drawn")

    cax = fig.add_axes([0.939, 0.035, 0.008, 0.82])
    cb = fig.colorbar(im, cax=cax, extend="max")
    cb.set_label("2009 surface elevation (m MHW)", fontsize=8.5)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_edgecolor(P.INK_MUTED)

    fig.text(0.022, 0.985,
             "NC-12 as digitised, on the 2009 DEM — map frame, GIS 1–90",
             fontsize=15, va="top", weight="semibold")
    fig.text(0.022, 0.952,
             "1 m bare-earth LiDAR in its native UTM 18N grid; the geojsons "
             "reprojected onto it with to_crs and NOTHING else applied — no "
             "orient, no alongshore flip, no shear, no water trim.\n"
             "So this sheet checks the SOURCE data against the world; whether "
             "the Barrier3D domain arrays are right is the separate question "
             "the domain-frame figures answer.\n"
             "Strips run SOUTH (left) → NORTH (right). North is up and the "
             "scale is true in both axes in every panel — Hatteras is 45 km "
             "north–south against 8 km east–west, so every true-north segment "
             "is portrait.\n"
             f"Elevations are m MHW (NAVD88 − {MHW_M:.2f} m). Never-surveyed "
             "cells are drawn apart from water: this is the only figure in the "
             "tree reading the raw tif, so it is the only one that can tell "
             "them apart.\n"
             "The two road lines COINCIDE almost everywhere — NC-12 did not "
             "move between these dates outside the relocation blocks, so the "
             "1984 casing under the 2004 dashes is the expected reading. Look "
             "for two separate lines only around GIS 9–16 (1999 inter-village) "
             "and GIS 84–88 (1989 Pea Island).",
             fontsize=9, color=P.INK_SECOND, va="top", linespacing=1.55)

    fig.legend(handles=[
        Line2D([], [], color=P.C_1984, lw=LW_1984,
               label="1984 NC-12 (geojson)"),
        Line2D([], [], color=P.C_2004, lw=LW_2004 + 0.6,
               ls=(0, (3.2, 1.8)), label="2004 NC-12 (geojson, dashed)"),
        Line2D([], [], color=C_WATER, lw=8, label="at or below 0 m MHW"),
        Line2D([], [], color=C_NODATA, lw=8,
               label="never surveyed (raw NoData)"),
        Line2D([], [], color=C_BOX, lw=0.9,
               label="Barrier3D domain clip box (north-up, 2000 × 500 m)"),
    ], loc="upper left", bbox_to_anchor=(0.022, 0.878), ncol=5, fontsize=8.5,
        framealpha=0.0, borderpad=0.4, columnspacing=1.8, handlelength=2.6)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=DPI, facecolor=P.SURFACE)
    plt.close(fig)
    print(f"\n[out] {OUT_PNG}")
    print(f"      {fig_w:.0f} x {fig_h:.1f} in at {DPI} dpi")
    return 0


if __name__ == "__main__":
    sys.exit(main())
