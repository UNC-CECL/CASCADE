"""
Where on Hatteras Island overwash recurs, from the workbook's images over the full model window.

    python scripts/input_prep/8-overwash-analysis/2-record/overwash_record_heatmap.py

The workbook gives overwash present or absent per image and domain, so each
domain box is shaded by the number of images showing it overwashed: on the
island, per domain beside it, and in three true-scale hotspot zooms. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

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

from site_layer.hat_figure_style import (FIG_H_MAX, INK, INK_MUTED,  # noqa: E402
                                         _north_arrow, _scalebar, _title,
                                         apply_style, figsize, open_frame,
                                         save, spines_for_image)
from overwash_data import (FULL, PERIODS, SECTIONS, load_observations,  # noqa: E402
                           period_images, period_label, sources_text,
                           upsert_caption)
from overwash_map_periods import (CLR_BOX, CLR_COAST, CLR_LAND,  # noqa: E402
                                  CLR_ROAD, CLR_WATER, PAD_E, PAD_N, PAD_S,
                                  draw_island, load_geometry,
                                  scalebar_and_north)

from site_layer import hat_overwash as ow  # noqa: E402
FIG_DIR = ow.WHERE

# --- CONFIG ------------------------------------------------------------------
# True-scale zooms on the hotspots: (title, first domain, last domain)
ZOOMS = [("Rodanthe–Pea Island", 82, 90),
         ("Wimble Shoals", 45, 54),
         ("Cape Point", 1, 7)]
ZOOM_GAP_IN = 0.45
# -----------------------------------------------------------------------------


# Inputs

# Images showing overwash, and images reaching the domain, per domain
def per_domain(obs, domains, matrix):
    hits = np.nansum(matrix, axis=0).astype(int)
    seen = (~np.isnan(matrix)).sum(axis=0)
    return [(int(d), int(h), int(s)) for d, h, s in zip(domains, hits, seen)]


# Drawing

# Discrete purples for 1..vmax, darker at the light end than the period maps so one image shows on land
def count_cmap(vmax):
    return [cm.Purples(0.42 + 0.53 * (i - 1) / max(vmax - 1, 1)) for i in range(1, vmax + 1)]


# The window around domains lo..hi, shaped to an axes of the given aspect (h / w)
def zoom_window(dom, lo, hi, aspect):
    x0, y0, x1, y1 = dom[dom["domain_id"].between(lo, hi)].total_bounds
    xc, yc = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    h = (y1 - y0) + 400.0
    w = max(h / aspect, (x1 - x0) + 400.0)
    h = w * aspect
    return xc - w / 2, yc - h / 2, xc + w / 2, yc + h / 2


# One hotspot at true scale: land, domain boxes shaded on land and numbered, NC-12.
# grid = (counts, extent, cmap, norm) draws a cell count instead (overwash_gis_heatmap.py).
def draw_zoom(ax, dom, land, road, fill_of, win, letter, title, grid=None):
    ax.set_facecolor(CLR_WATER)
    for geom in getattr(land, "geoms", [land]):
        ax.add_patch(mpatches.Polygon(np.asarray(geom.exterior.coords), closed=True,
                                      facecolor=CLR_LAND, edgecolor=CLR_COAST,
                                      lw=0.5, zorder=1))
    if grid is not None:
        counts, extent, cmap, norm = grid
        ax.imshow(np.ma.masked_less(counts, 1), extent=extent, origin="upper",
                  cmap=cmap, norm=norm, interpolation="nearest", zorder=2)
    for _, r in dom.iterrows():
        g = r.geometry
        fc = fill_of.get(int(r["domain_id"]))
        if fc is not None:
            on_land = g.intersection(land)
            for part in getattr(on_land, "geoms", [on_land]):
                if part.geom_type == "Polygon" and not part.is_empty:
                    ax.add_patch(mpatches.Polygon(np.asarray(part.exterior.coords),
                                                  closed=True, facecolor=fc,
                                                  edgecolor="none", zorder=2))
        ax.add_patch(mpatches.Polygon(np.asarray(g.exterior.coords), closed=True,
                                      facecolor="none", edgecolor=CLR_BOX,
                                      lw=0.4, zorder=3))
        x0, y0, x1, y1 = g.bounds
        yc = 0.5 * (y0 + y1)
        if win[1] < yc < win[3] and x1 + 120 < win[2]:
            ax.text(x1 + 120, yc, str(int(r["domain_id"])), fontsize=6.5,
                    va="center", ha="left", color=INK_MUTED, zorder=5, clip_on=True)
    road.plot(ax=ax, color=CLR_ROAD, lw=0.8, zorder=4)
    ax.set_xlim(win[0], win[2])
    ax.set_ylim(win[1], win[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    _title(ax, ord(letter) - ord("a"), title)
    _scalebar(ax, 1000.0, show_cells=False)


# Images with overwash per domain, level with the island, reaches named on the right
def draw_profile(ax, dom, rows, cols, bounds, vmax,
                 xlabel="Images with overwash in the domain"):
    ybox = {int(r["domain_id"]): (r.geometry.bounds[1], r.geometry.bounds[3])
            for _, r in dom.iterrows()}
    for d, n, _ in rows:
        if n == 0:
            continue
        y0, y1 = ybox[d]
        ax.barh(0.5 * (y0 + y1), n, height=(y1 - y0) * 0.86, color=cols[n - 1],
                edgecolor="none", zorder=3)
    for _, lo, _, _ in SECTIONS[1:]:
        ax.axhline(ybox[lo][0], color="#d0d0d0", lw=0.6, zorder=1)
    for name, lo, hi, kind in SECTIONS:
        yc = 0.5 * (ybox[lo][0] + ybox[hi][1])
        ax.text(vmax + 0.3, yc, name, fontsize=7,
                va="center", ha="left", linespacing=1.1,
                color=INK if kind == "village" else INK_MUTED,
                fontstyle="normal" if kind == "village" else "italic")
    ax.set_xlim(0, vmax + 0.2)
    ax.set_xticks(range(0, vmax + 1, 2))
    ax.set_ylim(bounds[1] - PAD_S, bounds[3] + PAD_N)
    ax.set_yticks([])
    ax.set_xlabel(xlabel)
    open_frame(ax)
    ax.spines["left"].set_visible(False)
    ax.grid(axis="x", color="#e6e6e6", lw=0.5, zorder=0)
    _title(ax, 1, "")


# Caption

def caption_text(obs, rows, span):
    n_img = len(obs)
    n_poor = int(obs["poor"].sum())
    lo, hi = PERIODS[FULL]
    n_dom = sum(1 for r in rows if r[1] > 0)
    top = sorted(rows, key=lambda r: (-r[1], r[0]))[:6]
    top_s = "; ".join(f"domain {d} ({h} of {s} images)" for d, h, s in top)
    seen = [r[2] for r in rows]
    return (
        f"Where overwash recurred on Hatteras Island over the full model window, "
        f"{span}, from the image record in Hatteras_Overwash_Data.xlsx: the "
        f"{n_img} images taken in its storm years, {lo}–{hi} ("
        f"{sources_text(obs, list(obs.index))}; {n_poor} flagged poor or partial "
        f"and kept), the same images the model comparison uses. Earlier images "
        f"(1984–1996) are left out. Each image records overwash present "
        f"or absent per CASCADE domain. (a) The island (NC 1:80k coastline, "
        f"UTM 18N, true scale) with each domain box, clipped to land, shaded by "
        f"the number of images showing it overwashed; unshaded = none. Domain "
        f"boxes numbered every ten; NC-12 as the dark line. (b) The same count "
        f"per domain, level with the island, in the same shades. (c–e) The three "
        f"most-hit reaches at true scale ("
        + "; ".join(f"{t}, domains {lo}–{hi}" for t, lo, hi in ZOOMS) +
        f"), outlined in (a), with the domain boxes numbered. {n_dom} of 90 "
        f"domains were "
        f"overwashed in at least one image. Most often: {top_s}. Each domain is "
        f"covered by {min(seen)}–{max(seen)} of the {n_img} images (an image may "
        f"not reach every domain), so the counts are comparable along the island. "
        f"The shading is per domain: it says a domain was overwashed in an image, "
        f"not where in the domain or how much of it.")


# Run

def main():
    apply_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    obs, domains, matrix = load_observations()
    keep = period_images(obs, FULL)
    obs, matrix = obs.loc[keep].reset_index(drop=True), matrix[keep]
    span = period_label(FULL)
    out_png = FIG_DIR / f"overwash_where_record_{FULL}.png"
    dom, land, road, bounds = load_geometry()

    rows = per_domain(obs, domains, matrix)
    vmax = max(r[1] for r in rows)
    cols = count_cmap(vmax)
    fill_of = {d: cols[h - 1] for d, h, _ in rows if h > 0}
    print(f"  {len(obs)} images, {span}; most images showing overwash: {vmax}")

    pad_w = 1500.0
    win = (bounds[0] - pad_w, bounds[1] - PAD_S, bounds[2] + PAD_E, bounds[3] + PAD_N)

    # Island | by-domain profile and its reach labels | three zooms stacked
    fig_w = figsize("double")[0]
    fig_h = FIG_H_MAX - 0.4
    TM, BM, LM, RM = 0.30, 1.05, 0.05, 0.05    # two legend rows under the axis label
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

    draw_island(ax_m, dom, land, road, bounds, fill_of, "a", "Images with overwash",
                reach_labels=False, pad_w=pad_w)
    # The scale bar to the east of the island, clear of the Cape Point outline
    x0, x1 = ax_m.get_xlim()
    ax_m.set_xlim(x0 + 0.44 * (x1 - x0), x1 + 0.44 * (x1 - x0))
    scalebar_and_north(ax_m)
    ax_m.set_xlim(x0, x1)
    draw_profile(ax_p, dom, rows, cols, bounds, vmax)

    for k, (title, lo, hi) in enumerate(ZOOMS):
        y = BM + map_h - (k + 1) * zoom_h - k * ZOOM_GAP_IN
        ax_z = axes_in(x, y, zoom_w, zoom_h)
        zw = zoom_window(dom, lo, hi, zoom_h / zoom_w)
        letter = chr(ord("c") + k)
        draw_zoom(ax_z, dom, land, road, fill_of, zw, letter,
                  f"{title}, domains {lo}–{hi}")
        if k == 0:
            _north_arrow(ax_z, x=0.08, y=0.80, length=0.09)
        ax_m.add_patch(mpatches.Rectangle((zw[0], zw[1]), zw[2] - zw[0],
                                          zw[3] - zw[1], fill=False, ec=INK,
                                          lw=0.7, zorder=6))
        ax_m.text(zw[0] + 200, zw[3] - 200, letter, fontsize=7.5,
                  fontweight="bold", ha="left", va="top", color=INK, zorder=6)

    handles = [mpatches.Patch(fc=cols[i], ec="none",
                              label=f"{i + 1} image{'s' if i else ''}")
               for i in range(vmax)]
    handles += [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4,
                               label="CASCADE domain, none"),
                plt.Line2D([0], [0], color=CLR_ROAD, lw=0.9, label="NC-12")]
    fig.legend(handles=handles, loc="lower left", ncol=6, frameon=False,
               bbox_to_anchor=(LM / fig_w, 0.0), handlelength=1.4,
               handleheight=1.0, columnspacing=1.2,
               title=f"Images with overwash, {span} ({len(obs)} images)",
               title_fontsize=8, alignment="left")

    save(fig, out_png, close=True)
    print(f"  wrote {out_png.relative_to(REPO)} (+ .pdf)")
    p = upsert_caption(out_png.name, "where", caption_text(obs, rows, span))
    print(f"  wrote {p.relative_to(REPO)}")
    for d, h, s in sorted(rows, key=lambda r: (-r[1], r[0]))[:10]:
        print(f"    domain {d:2d}: {h} of {s} images")


if __name__ == "__main__":
    main()
