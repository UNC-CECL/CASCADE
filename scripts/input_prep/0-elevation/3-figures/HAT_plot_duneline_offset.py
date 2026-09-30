"""
The 1984-start DEM with both digitized dune lines on it, and the cross-shore distance between them per domain.

    python scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py
    python scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py --simple
    python scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py --zoom 62-68

Measures the 1984-vs-1997 offset into duneline_offset_by_domain.csv and draws
the island, ribbon, zoom and simple figures under
data/hatteras_init/0-elevation/2009-2014-1996-duneline/figures/. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
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


# Walk up until a directory holds data/hatteras_init
def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from site_layer.hat_elevation_products import ELEVATION_ROOT  # noqa: E402
# Mosaic loader, km axes and elevation panel imported from HAT_plot_1984_mosaic, so they cannot drift
import HAT_plot_1984_mosaic as m  # noqa: E402

# The place names are NOT redefined here
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS  # noqa: E402


# House style from site_layer/hat_figure_style.py, re-exported for the scripts that import this one
from site_layer.hat_figure_style import (  # noqa: E402,F401
    FONT_STACK, INK, INK_MUTED, GRID_C, C_1984, C_1997, C_1984_FILL, C_1997_FILL,
    STYLE_RC, apply_style, _letter, _title, _letter_inside, _north_arrow, _halo,
    figsize, FIG_W_DOUBLE, FIG_H_MAX, DOMAIN_AXIS_LABEL, town_bands, open_frame,
    save,
)
from site_layer import hat_figure_style as _style  # noqa: E402


SOURCE_TAG = "2009-2014-1996"
from site_layer.hat_elevation_products import duneline_check_dir, product as _elprod  # noqa: E402
OUT_DIR = duneline_check_dir(SOURCE_TAG)
FIG_DIR = OUT_DIR / "figures"
# Figures/ is sorted by what a figure IS (2026-09-08, Hannah)
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


# figures/<kind>/<name>, the folder made
def fig_path(name):
    p = FIG_DIR / FIG_SUBFOLDER[name] / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
# --- CONFIG ------------------------------------------------------------------
CSV_NAME = "duneline_offset_by_domain.csv"

from site_layer.hat_topo_version import DUNELINE_DIR as DUNE_DIR  # noqa: E402
DUNE_LINES = {1984: DUNE_DIR / "duneline_1984.geojson",
              1997: DUNE_DIR / "duneline_1997.geojson"}

# One sample per alongshore metre, matching the 1 m DEM (500 per domain)
SAMPLE_SPACING_M = 1.0
GRID_10M = 10.0          # the resampled product's cell, for the box fallback

# The 1984 footprint table, read only to label a --zoom with its Barrier3D rows
from site_layer import hat_topo_version as _b3d  # noqa: E402
INSERT_SCOPE_CSV = (_b3d.insert_scope_step("1984-start", "2-extent")
                    / "footprint_1984_by_domain.csv")   # step folder since 2026-09-09

# The box shape the easting-is-cross-shore frame depends on
EXPECTED_BOX_M = (2000.0, 500.0)
BOX_TOL_M = 1.0

# One key for every figure: both lines solid, 1984 the emphasis colour, 1997 the reference
LINE_STYLE = {1984: dict(color=C_1984, linestyle="-", linewidth=2.0),
              1997: dict(color=C_1997, linestyle="-", linewidth=2.0)}
LINE_CASING = {1984: dict(color="white", linewidth=3.0),
               1997: dict(color="white", linewidth=3.0)}
LINE_ORDER = [1997, 1984]      # 1984 drawn last, on top

# The per-figure line WIDTH, as a multiplier on LINE_STYLE
LINE_SCALE_ISLAND = 0.55
LINE_SCALE_DETAIL = 1.0
# The ribbon is a trace, not a map
LINE_SCALE_RIBBON = 0.7

BOX_STYLE = dict(edgecolor="0.30", facecolor="none", linewidth=0.4)
LABEL_EVERY = 5                # label every Nth domain box on the island map

# The island is 46 km long and ~2 km wide
N_PANELS = 3
PANEL_PAD_M = 400.0

CELL_M = 10.0                  # the Barrier3D cell, drawn on the bar chart

# The ribbon's baseline: a 2 km boxcar over the two lines' mean
BASELINE_WINDOW_M = 2000.0

# Zoom reaches from the measured table: largest negative, largest positive, and a quiet control
ZOOM_REACHES = [
    (17, 21, "17-21", "the quietest reach on the island"),
    (62, 68, "62-68", "1984 line landward of 1997"),
    (78, 85, "78-85", "1984 line seaward of 1997"),
]
# Cross-shore half-width of a zoom panel
ZOOM_HALF_WIDTH_M = 300.0

# A stripped version of the same panels
SIMPLE_HALF_WIDTH_M = 150.0
# South to north, matching the locator
SIMPLE_REACHES = [
    (3, 4, ""),
    (19, 20, "control"),
    (63, 64, ""),
    (79, 80, ""),
]

# The island locator: a map per third of the island, with the offset bar beside it
ISLAND_PANELS = 3
ISLAND_PAD_M = 500.0
STRIP_WIDTH_RATIO = 0.42       # bar axes width, as a fraction of the map's
OFFSET_POS = C_1984            # 1984 seaward - the sign erosion predicts
OFFSET_NEG = C_1997            # 1984 landward
HIGHLIGHT = INK                # the label and the band marking a detail pair

# Piers and the groin are drawn SEAWARD from the 1984 line, because that is where they are
STRUCTURE_LEN_M = 380.0

# THE LINES-ONLY LOCATOR (--lines-island)
LINES_ISLAND_DOMAINS = 10
LINES_ISLAND_PAD_M = 150.0
LINES_ISLAND_LINE_SCALE = 0.9

# The alongshore ruler, in km from the south end of domain 1 (as in fig_ribbon)
KM_TICK_M = 5000.0
# The simple figures' line key, a copy of LINE_STYLE
SIMPLE_LINE_STYLE = {yr: dict(LINE_STYLE[yr]) for yr in LINE_STYLE}

# The 1 m gapfilled tiles, read in place
TILE_1M_DIR = _elprod(SOURCE_TAG, check=False).gapfill_1m
TILE_1M_NAME = "clip_domain_{d}_filled.tif"
TILE_BOUNDS_TOL_M = 1.5

# Hillshade drawing parameters: relief only, they cannot move either line
HILLSHADE = dict(azdeg=315.0, altdeg=40.0, vert_exag=2.2)
SHADE_SMOOTH_M = 3.0
NODATA_GREY = "#e8e8e8"

SCALEBAR_M = 50.0              # 5 Barrier3D cells
# -----------------------------------------------------------------------------


# The measurement

# Easting of a line at each of `ys`, inside one domain box
def x_at_northings(line_geom, box, ys):
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


# Per-domain offset between the two dune lines
def measure(gdf, lines):
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

        # Orientation-free check: unsigned nearest-point distance to the unclipped 1997 line
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


# Figures

# Both dune lines, reprojected
def load_lines(dst_crs):
    out = {}
    for yr, p in DUNE_LINES.items():
        g = gpd.read_file(p)
        src_crs = g.crs
        if g.crs is not None and dst_crs is not None:
            g = g.to_crs(dst_crs)
        # EPSG codes only; the compound CRS's full WKT buries the log
        print(f"  {yr}: {p.name}  {_epsg(src_crs)} -> {_epsg(dst_crs)}   "
              f"{g.geometry.length.sum() / 1000:.1f} km")
        out[yr] = g
    return out


# The same clip HAT_plot_1984_mosaic.load_roads applies to NC-12, and for the same reason
def clip_for_drawing(lines, footprint):
    out = {}
    for yr, g in lines.items():
        before = float(g.geometry.length.sum())
        c = gpd.clip(g, footprint)
        after = float(c.geometry.length.sum()) if len(c) else 0.0
        print(f"  {yr}: clipped to domains for drawing, "
              f"{after / 1000:.1f} km of {before / 1000:.1f} km kept")
        out[yr] = c
    return out


# Casing then line, in LINE_ORDER so 1984 lands on top of 1997
def draw_lines(ax, lines, scale=1.0, style=None):
    style = style or LINE_STYLE
    for yr in LINE_ORDER:
        cas = dict(LINE_CASING[yr])
        cas["linewidth"] *= scale
        st = dict(style[yr])
        st["linewidth"] *= scale
        lines[yr].plot(ax=ax, linestyle="-", alpha=0.9,
                       zorder=6 + LINE_ORDER.index(yr) * 2, **cas)
        lines[yr].plot(ax=ax, zorder=7 + LINE_ORDER.index(yr) * 2, **st)


# Legend handles for the two dune lines
def line_legend():
    return [Line2D([0], [0], label=f"{yr} dune line", **LINE_STYLE[yr])
            for yr in (1984, 1997)]


# RETIRED 2026-09-08 (not called)
def fig_island(elev, extent, gdf, lines, rows):
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


# Boxcar over an alongshore series, edges held rather than tapered
def _smooth(a, win_px):
    k = np.ones(win_px) / win_px
    pad = win_px // 2
    return np.convolve(np.pad(a, pad, mode="edge"), k, mode="same")[
        pad:pad + a.size]


# The two lines against a common baseline, with the band between them filled
def fig_ribbon(samples, rows, gdf):
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

    # (a) the two lines about the midline
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

    # (b) the difference
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

    # Headroom for the village names on (a) and the detail brackets on (b)
    ax.set_xlim(dom.min() - 0.5, dom.max() + 0.5)
    for a_, frac in ((ax, 0.14), (bx, 0.30)):
        lo_, hi_ = a_.get_ylim()
        a_.set_ylim(lo_, hi_ + frac * (hi_ - lo_))
    town_bands(ax)
    town_bands(bx, label=False)

    # The reaches the true-scale zooms cover, as brackets along the top of (b)
    for lo, hi, slug, _ in ZOOM_REACHES:
        bx.plot([lo - 0.5, hi + 0.5], [0.95, 0.95], color=INK_MUTED, lw=0.9,
                solid_capstyle="butt", zorder=4,
                transform=bx.get_xaxis_transform())
        bx.text((lo + hi) / 2, 0.92, f"detail {slug}", fontsize=7,
                ha="center", va="top", color=INK_MUTED,
                transform=bx.get_xaxis_transform())

    # Distance along the top, in km from the south end of domain 1
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


# True-scale panels on the reaches where the offset is largest, plus a quiet control - since ...
def fig_zooms(elev, extent, gdf, lines, rows, reaches=None, out=None,
              half_width=None, rows_by_domain=None):
    # A hand-given --zoom-note is the only note drawn under a panel title
    custom = reaches is not None
    reaches = reaches or ZOOM_REACHES
    half_width = ZOOM_HALF_WIDTH_M if half_width is None else half_width
    return fig_zooms_simple(gdf, lines, rows, elev=elev, extent=extent,
                            reaches=[(lo, hi, note if custom else "")
                                     for lo, hi, _slug, note in reaches],
                            out=out or fig_path("HAT_duneline_offset_zooms.png"),
                            half_width=half_width, rows_by_domain=rows_by_domain,
                            tall=True)


# The 1 m gapfilled tiles for `ids`, mosaicked onto one array
def load_1m(gdf, ids):
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


# Grey relief
def _hillshade(ax, arr, extent, res=1.0):
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


# Where GIS domains lo..hi are, in the site's own vocabulary
def _place_of(lo, hi, ann=HATTERAS_ANNOTATIONS):
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


# What was measured on a pair or reach, as caption text
def _pair_reading(meds):
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


# Northing bounds of GIS domains lo..hi, or None if none are on the grid
def _span_y(gdf, lo, hi):
    sub = gdf[gdf["domain_id"].astype(int).isin(
        range(int(np.floor(lo)), int(np.ceil(hi)) + 1))]
    if sub.empty:
        return None
    b = sub.total_bounds
    return float(b[1]), float(b[3])


# The communities, as a bracket in the ocean margin - not a wash
def _places(ax, gdf, y0, y1, ann=HATTERAS_ANNOTATIONS):
    xa, xb = ax.get_xlim()
    xbr = xb - 0.055 * (xb - xa)          # the bracket
    xtx = xbr - 0.02 * (xb - xa)          # community names, landward of it
    # Village names get a column of their own
    xtv = xbr - 0.11 * (xb - xa)

    for name, (lo, hi) in ann.town_spans.items():
        yy = _span_y(gdf, lo, hi)
        if yy is None or yy[1] < y0 or yy[0] > y1:
            continue
        a, b = max(yy[0], y0), min(yy[1], y1)
        ax.plot([xbr, xbr], [a, b], color=ann.color_town_span, lw=5.0,
                solid_capstyle="butt", zorder=8)
        # Avon is GIS 21-31 and the panels break at 30, so it is drawn in two pieces
        if (b - a) < 0.5 * (yy[1] - yy[0]):
            continue
        # Rotated to run along the bracket
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


# Median easting of the 1984 line at GIS id `gid`, fractional allowed
def _x1984(by_dom, gid):
    import math
    ids = {math.floor(gid), math.ceil(gid)}
    xs = [by_dom[d]["x1984_med"] for d in ids
          if d in by_dom and by_dom[d]["x1984_med"] != ""]
    return float(np.mean(xs)) if xs else None


# The piers and the groin, as seaward marks off the 1984 line
def _structures(ax, gdf, by_dom, y0, y1, ann=HATTERAS_ANNOTATIONS):
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


# Alongshore distance in km, on the bar strip's left edge
def _km_axis(bar, y_origin, y0, y1):
    lo = np.ceil((y0 - y_origin) / KM_TICK_M) * KM_TICK_M
    if lo <= 0.0:
        # Assign 0.0: ceil() returns -0.0 here and the tick would read -0
        lo = 0.0
    hi = (y1 - y_origin)
    vals = np.arange(lo, hi + 1.0, KM_TICK_M)
    bar.set_yticks([y_origin + v for v in vals])
    bar.set_yticklabels([f"{v / 1000:.0f}" for v in vals])
    bar.tick_params(axis="y", labelsize=7, length=2.5, pad=1.5,
                    color="0.55", labelcolor="0.35")


# The site's own name for what lies off the end of the reach
def _end_label(ax, text, at_top, y0, y1):
    xa, xb = ax.get_xlim()
    # Right of centre, clear of the panel letter
    ax.text(xa + 0.58 * (xb - xa), y1 if at_top else y0, text,
            fontsize=8, fontstyle="italic", ha="center",
            va="top" if at_top else "bottom", color=INK_MUTED, zorder=12,
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                      boxstyle="square,pad=0.22"))


# The house scale bar (hat_figure_style._scalebar), with this module's cell size
def _scalebar(ax, length_m=SCALEBAR_M, show_cells=None):
    _style._scalebar(ax, length_m=length_m, cell_m=CELL_M,
                     show_cells=show_cells)


# The per-domain statistic the bar strip can draw
ISLAND_STATS = {"median": "offset_med_m", "mean": "offset_mean_m"}


# The whole island, and the measured offset beside it
def fig_island_simple(elev, extent, gdf, lines, rows, reaches=None,
                      out=None, stat="median"):
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

    # THE COLUMN WIDTHS ARE COMPUTED, NOT CHOSEN
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
    # Printed at the double-column width
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

        # the map
        _hillshade(ax, elev, extent, res=GRID_10M)
        sub.boundary.plot(ax=ax, color="0.45", linewidth=0.35, zorder=4)
        draw_lines(ax, lines, scale=LINE_SCALE_ISLAND,
                   style=SIMPLE_LINE_STYLE)
        for lo, hi, note in reaches:
            if not (group.min() <= lo <= group.max()):
                continue
            r = gdf[gdf["domain_id"].astype(int).isin(range(lo, hi + 1))]
            rb = r.total_bounds
            # No outline round the pair; the label and the bar's band carry it
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
        # The span alone: the narrowest map column is about 33 mm
        _title(ax, i, f"{group.min()}\u2013{group.max()}")
        _places(ax, gdf, y0, y1)
        _structures(ax, gdf, by_dom, y0, y1)
        if i == 0:
            _end_label(ax, HATTERAS_ANNOTATIONS.low_end_label, False, y0, y1)
            _north_arrow(ax, x=0.88, y=0.09)
        if i == ISLAND_PANELS - 1:
            _end_label(ax, HATTERAS_ANNOTATIONS.high_end_label, True, y0, y1)
            # One scale bar: the three maps share a scale
            _scalebar(ax, length_m=2000.0)

        # the offset beside it
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
        # Domain numbers in the margin beyond the bar, not on the map
        for d in group:
            if d % 10 == 0 or d in (ids.min(), ids.max()):
                b = gdf[gdf["domain_id"].astype(int) == d].total_bounds
                bar.text(lim * 1.08, (b[1] + b[3]) / 2, str(d), fontsize=7,
                         ha="left", va="center", color=INK_MUTED, zorder=6,
                         clip_on=False)

    # One legend for the whole figure, below the panels
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


# The whole island as maps only
def fig_island_lines(elev, extent, gdf, lines, rows, reaches=None,
                     out=None, per_panel=LINES_ISLAND_DOMAINS,
                     pad_m=LINES_ISLAND_PAD_M):
    from shapely.geometry import box as _box
    from matplotlib.patches import Patch

    reaches = reaches or SIMPLE_REACHES
    ids = np.sort(gdf["domain_id"].astype(int).to_numpy())
    n_panels = int(np.ceil(len(ids) / per_panel))
    groups = np.array_split(ids, n_panels)
    by_dom = {r["domain"]: r for r in rows}
    med = np.array([r["offset_med_m"] for r in rows
                    if r["offset_med_m"] != ""], float)

    # Each panel's window: northing from the domain boxes, easting from the lines
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
    # Printed at the double-column width
    n_rows = 1 if n_panels <= 4 else 2
    per_row = int(np.ceil(n_panels / n_rows))
    row_ids = [list(range(k, min(k + per_row, n_panels)))
               for k in range(0, n_panels, per_row)]
    row_sum = max(sum(ratios[j] for j in rg) for rg in row_ids)
    # The panel height is whichever binds
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
    # A panel is ~20 mm wide on the page, so nothing fits beside a panel letter over it
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
            # One scale bar: all panels share a scale
            _scalebar(ax, length_m=500.0, show_cells=False)
        # Domain numbers on the landward edge, every fifth
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


# The two lines, and as little else as the picture can carry
def fig_zooms_simple(gdf, lines, rows, elev=None, extent=None, reaches=None,
                     out=None, half_width=None, rows_by_domain=None, tall=False):
    apply_style()
    reaches = reaches or SIMPLE_REACHES
    half_width = SIMPLE_HALF_WIDTH_M if half_width is None else half_width
    by_dom = {r["domain"]: r for r in rows}

    # The panel width follows the crop
    spans = [gdf[gdf["domain_id"].astype(int).isin(range(lo, hi + 1))].total_bounds
             for lo, hi, _ in reaches]
    span_max = max(b[3] - b[1] for b in spans)
    # Printed at the double-column width
    n = len(reaches)
    h_axes = min((FIG_W_DOUBLE - 0.25 * (n + 1)) / (n * 2.0 * half_width / span_max),
                 FIG_H_MAX - 1.1)
    w_axes = h_axes * 2.0 * half_width / span_max
    # A narrow panel drops 'Domains' from its title
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

        # The backdrop covers the whole panel
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
            # Reaches of unequal length: every panel spans the longest, so the tops line up
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

        # The span alone, on one line, so the panel titles stay level
        _title(ax, i, head.format(lo, hi)
               + (f" \u00b7 {note}" if note and n == 1 else ""))
        if i == 0:
            # One bar and one arrow for the figure
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


# The requested number, domain by domain, with its spread
def fig_by_domain(rows):
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


# One arbitrary zoom

# A CRS as 'EPSG:nnnn', falling back to its name
def _epsg(crs):
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


# The 90 domain boxes, from `domains.geojson` if it is reachable and from the resampled rasters if it ...
def load_domains():
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


# The per-domain offsets back off disk, so a zoom does not re-measure
def read_rows(path):
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


# {domain
def read_insert_rows(path):
    if not Path(path).exists():
        print(f"  NOTE: {Path(path).name} absent - row counts not labelled")
        return None
    return {int(r["domain"]): int(r["n_cells"]) for r in
            csv.DictReader(open(path))}


# Render the two simple figures alone, off the existing table
def simple_only(out, half_width, span=None, island_out=None,
                stat="median"):
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

    # The 10 m mosaic the locator draws, loaded unconditionally
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


# Render the lines-only island figure alone, off the existing table
def lines_island_only(out, per_panel, pad_m):
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


# Render ONE true-scale zoom for an arbitrary domain span
def zoom_only(span, out, half_width, note, with_rows):
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
    # Only a note given by hand goes under the panel title
    p = fig_zooms(elev, extent, gdf, drawn, rows,
                  reaches=[(lo, hi, f"{lo}-{hi}", note or "")],
                  out=out, half_width=half_width,
                  rows_by_domain=n_by)
    print(f"\n  figure : {p}  (offset {min(meds):+.0f} to {max(meds):+.0f} m)")
    return p


# A caption per figure, with the numbers filled from the table
def write_captions(rows, half_width=None, simple_half_width=None):
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


# Run: measure the offsets, then draw the figures the options ask for
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
