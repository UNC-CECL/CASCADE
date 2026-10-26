"""
Where on Hatteras Island overwash recurs, from the digitised washover-fan outlines.

    python scripts/input_prep/8-overwash-analysis/2-record/overwash_gis_heatmap.py

Reads the overwash_events feature class from the ArcGIS project's Default.gdb,
counts the years each 10 m cell was overwashed, and draws that count on the
island beside the years overwashed per CASCADE domain, in the layout of
overwash_record_heatmap.py. The outlines cover 1984-2004 only, outside the
model windows. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from rasterio.features import rasterize
from rasterio.transform import from_origin
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "8-overwash-analysis"
                       / "1-observations"))

from site_layer.hat_figure_style import (FIG_H_MAX, INK, _north_arrow,  # noqa: E402
                                         apply_style, figsize, save)
from overwash_data import upsert_caption  # noqa: E402
from site_layer.hat_gis_root import GIS_PROJECT  # noqa: E402
from overwash_map_periods import (CLR_BOX, CLR_LAND, CLR_ROAD, PAD_E,  # noqa: E402
                                  PAD_N, PAD_S, draw_island, load_geometry,
                                  scalebar_and_north)
from overwash_record_heatmap import (ZOOM_GAP_IN, ZOOMS, count_cmap,  # noqa: E402
                                     draw_profile, draw_zoom, zoom_window)

from site_layer import hat_overwash as ow  # noqa: E402
FIG_DIR = ow.WHERE
OUT_PNG = FIG_DIR / "overwash_where_polygons_1984_2004.png"

# --- CONFIG ------------------------------------------------------------------
GDB = GIS_PROJECT / "hatteras" / "Default.gdb"
LAYER = "overwash_events"
CELL_M = 10.0          # heat-map grid, metres
# -----------------------------------------------------------------------------


# Inputs

# The overwash polygons, in the domain boxes' CRS
def load_events(crs):
    if not GDB.exists():
        raise SystemExit(f"\n{GDB} is missing.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")      # measured (M) geometry dropped to Z
        ev = gpd.read_file(GDB, layer=LAYER)
    ev = ev[ev["year"].notna()].to_crs(crs)
    ev["year"] = ev["year"].astype(int)
    ev["geometry"] = ev.geometry.force_2d().buffer(0)
    print(f"  {len(ev)} polygons, {ev['year'].nunique()} years "
          f"({ev['year'].min()}-{ev['year'].max()})")
    return ev


# Years overwashed per grid cell; a year's overlapping polygons count once
def year_count_grid(ev, bounds):
    x0, y0, x1, y1 = bounds
    nx = int(np.ceil((x1 - x0) / CELL_M))
    ny = int(np.ceil((y1 - y0) / CELL_M))
    tf = from_origin(x0, y1, CELL_M, CELL_M)
    grid = np.zeros((ny, nx), dtype=np.int16)
    for _, sub in ev.groupby("year"):
        shape = unary_union(sub.geometry.values)
        grid += rasterize([(shape, 1)], out_shape=(ny, nx), transform=tf,
                          fill=0, dtype="uint8")
    return grid, (x0, x0 + nx * CELL_M, y1 - ny * CELL_M, y1)


# Per domain box: (domain, years with overwash, hectares overwashed), and the years
def per_domain(ev, dom):
    rows, years = [], {}
    for _, r in dom.iterrows():
        d = int(r["domain_id"])
        hit = ev[ev.intersects(r.geometry)]
        clipped = hit.geometry.intersection(r.geometry)
        area = unary_union(clipped.values).area if len(hit) else 0.0
        rows.append((d, hit["year"].nunique(), area / 1e4))
        years[d] = sorted(hit["year"].unique().tolist())
    return rows, years


# Caption

def caption_text(ev, rows):
    years = sorted(ev["year"].unique())
    n_dom = sum(1 for r in rows if r[1] > 0)
    top = sorted(rows, key=lambda r: (-r[1], -r[2]))[:5]
    top_s = "; ".join(f"domain {d} ({n} yr, {a:.1f} ha)" for d, n, a in top)
    by_year = ev.groupby("year").size()
    return (
        f"Where overwash recurred on Hatteras Island, {years[0]}–{years[-1]}, "
        f"from the washover-fan outlines (overwash_events polygons) in the "
        f"Hatteras ArcGIS project (Default.gdb, {len(ev)} polygons over "
        f"{len(years)} years: "
        + ", ".join(f"{y} ({by_year[y]})" for y in years) +
        f"). These years lie outside the model windows; the windowed record is "
        f"overwash_where_record_1996_2025.png, in the same layout. (a) The island "
        f"(NC 1:80k coastline, UTM 18N, true scale) with each {CELL_M:.0f} m cell "
        f"shaded by the number of years a polygon covers it; polygons of the same "
        f"year that overlap count once. CASCADE domain boxes in grey, numbered "
        f"every ten; NC-12 as the dark line. (b) The number of years with any "
        f"overwash polygon inside each domain box, level with the island, in the "
        f"same shades. (c–e) Three hotspot reaches at true scale ("
        + "; ".join(f"{t}, domains {lo}–{hi}" for t, lo, hi in ZOOMS) +
        f"), outlined in (a), with the domain boxes numbered. {n_dom} of 90 "
        f"domains were overwashed at least once. Most often: {top_s}. The layer "
        f"carries only the year: its storm, date, source and confidence fields "
        f"are empty, so a year's polygons may come from more than one image or "
        f"storm, and years without polygons may be years without imagery rather "
        f"than years without overwash.")


# Run

def main():
    apply_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    dom, land, road, bounds = load_geometry()
    ev = load_events(dom.crs)

    pad_w = 1500.0
    win = (bounds[0] - pad_w, bounds[1] - PAD_S, bounds[2] + PAD_E, bounds[3] + PAD_N)
    grid, extent = year_count_grid(ev, win)
    rows, years = per_domain(ev, dom)
    vmax = int(max(grid.max(), max(r[1] for r in rows)))
    cols = count_cmap(vmax)
    cmap = ListedColormap(cols)
    norm = BoundaryNorm(np.arange(0.5, vmax + 1.5), vmax)

    # Island | by-domain profile and its reach labels | three zooms stacked
    fig_w = figsize("double")[0]
    fig_h = FIG_H_MAX - 0.4
    TM, BM, LM, RM = 0.30, 0.80, 0.05, 0.05
    PROF_W, PROF_LBL, GAP = 1.0, 1.15, 0.25
    map_h = fig_h - TM - BM
    map_w = map_h * (win[2] - win[0]) / (win[3] - win[1])
    fig = plt.figure(figsize=(fig_w, fig_h))

    def axes_in(x, y, w, h):
        return fig.add_axes([x / fig_w, y / fig_h, w / fig_w, h / fig_h])

    ax_m = axes_in(LM, BM, map_w, map_h)
    x = LM + map_w + GAP
    ax_p = axes_in(x, BM, PROF_W, map_h)
    x += PROF_W + PROF_LBL
    zoom_w = fig_w - RM - x
    zoom_h = (map_h - (len(ZOOMS) - 1) * ZOOM_GAP_IN) / len(ZOOMS)

    draw_island(ax_m, dom, land, road, bounds, {}, "a", "Years overwashed",
                reach_labels=False, pad_w=pad_w)
    ax_m.imshow(np.ma.masked_less(grid, 1), extent=extent, origin="upper",
                cmap=cmap, norm=norm, interpolation="nearest", zorder=2.5)
    # The scale bar to the east of the island, clear of the Cape Point outline
    x0, x1 = ax_m.get_xlim()
    ax_m.set_xlim(x0 + 0.44 * (x1 - x0), x1 + 0.44 * (x1 - x0))
    scalebar_and_north(ax_m)
    ax_m.set_xlim(x0, x1)
    draw_profile(ax_p, dom, rows, cols, bounds, vmax,
                 xlabel="Years with overwash in the domain")

    for k, (title, lo, hi) in enumerate(ZOOMS):
        y = BM + map_h - (k + 1) * zoom_h - k * ZOOM_GAP_IN
        ax_z = axes_in(x, y, zoom_w, zoom_h)
        zw = zoom_window(dom, lo, hi, zoom_h / zoom_w)
        letter = chr(ord("c") + k)
        draw_zoom(ax_z, dom, land, road, {}, zw, letter,
                  f"{title}, domains {lo}–{hi}", grid=(grid, extent, cmap, norm))
        if k == 0:
            _north_arrow(ax_z, x=0.08, y=0.80, length=0.09)
        ax_m.add_patch(mpatches.Rectangle((zw[0], zw[1]), zw[2] - zw[0],
                                          zw[3] - zw[1], fill=False, ec=INK,
                                          lw=0.7, zorder=6))
        ax_m.text(zw[0] + 200, zw[3] - 200, letter, fontsize=7.5,
                  fontweight="bold", ha="left", va="top", color=INK, zorder=6)

    handles = [mpatches.Patch(fc=cols[i], ec="none",
                              label=f"{i + 1} year{'s' if i else ''}")
               for i in range(vmax)]
    handles += [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4,
                               label="CASCADE domain"),
                plt.Line2D([0], [0], color=CLR_ROAD, lw=0.9, label="NC-12")]
    fig.legend(handles=handles, loc="lower left", ncol=7, frameon=False,
               bbox_to_anchor=(LM / fig_w, 0.01), handlelength=1.4,
               handleheight=1.0, columnspacing=1.2,
               title="Years with overwash, 1984–2004", title_fontsize=8,
               alignment="left")

    save(fig, OUT_PNG, close=True)
    print(f"  wrote {OUT_PNG.relative_to(REPO)} (+ .pdf)")
    p = upsert_caption(OUT_PNG.name, "where", caption_text(ev, rows))
    print(f"  wrote {p.relative_to(REPO)}")
    for d, n, a in sorted(rows, key=lambda r: (-r[1], -r[2]))[:10]:
        print(f"    domain {d:2d}: {n} yr, {a:6.1f} ha  {years[d]}")


if __name__ == "__main__":
    main()
