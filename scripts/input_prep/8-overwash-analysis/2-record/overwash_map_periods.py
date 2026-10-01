"""
Where on Hatteras Island overwash was observed, and when, for the two model periods.

    python scripts/input_prep/8-overwash-analysis/2-record/overwash_map_periods.py   # one figure per period
    python scripts/input_prep/8-overwash-analysis/2-record/overwash_map_periods.py --both   # plus the two-period figure

The island with its 90 domains shaded by images showing overwash, beside one
column per image. Details: scripts/input_prep/8-overwash-analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-22
"""

from __future__ import annotations

import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from matplotlib import cm
from shapely.geometry import box as shp_box
from shapely.ops import unary_union

HERE = Path(__file__).resolve().parent
REPO = next(
    _p for _p in HERE.parents
    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
# The stage's shared module and its siblings are found from the repo root
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "8-overwash-analysis"
                       / "1-observations"))

from site_layer.hat_figure_style import (C, FIG_H_MAX, INK, INK_MUTED, _north_arrow,   # noqa: E402
                              _scalebar, _title, apply_style, figsize, save,
                              spines_for_image)
from overwash_data import (OUT_DIR, PERIODS, SECTIONS, assign_capture,   # noqa: E402
                           load_observations, load_storms, upsert_caption)

from site_layer import hat_overwash as ow  # noqa: E402
FIG_DIR = ow.MAP
OUT_PNG = FIG_DIR / "overwash_map_periods.png"

from site_layer import hat_map_layers as _ml  # noqa: E402
DOMAIN_FILE = _ml.DOMAIN_BOXES
COAST_FILE = _ml.NC_COAST
from site_layer.hat_topo_version import road_line_file  # noqa: E402
ROAD_FILE = road_line_file(2008)

# --- CONFIG ------------------------------------------------------------------
CLR_OW = C["ACCENT"]
CLR_WATER = "#eef4f8"
CLR_LAND = "#e9e5dc"
CLR_COAST = "#a8a49a"
CLR_BOX = "#bdbdbd"
CLR_ROAD = C["ROAD"]
CLR_PART_FACE = "#f3f3f3"
CLR_PART_EDGE = "#9a9a9a"
CLR_INK = INK
CLR_MUTED = INK_MUTED

STRIP_IN_MAX = 0.30  # widest an image column gets, inches
COL_MIN_IN = 0.22    # two lines of rotated 7 pt type, the two-line head
COL_MIN_1LINE = 0.12  # one line of it
PAD_S, PAD_N = 1200.0, 1200.0
PAD_W, PAD_E = 7200.0, 3600.0
SCALEBAR_M = 5000.0

# Reach names broken for the sound side of a column-width island panel.
REACH_LABEL = {"Rodanthe–Waves–Salvo": "Rodanthe–\nWaves–Salvo"}
# -----------------------------------------------------------------------------


# Inputs

# The domain boxes, the land around them and NC-12, clipped to one window
def load_geometry():
    if not DOMAIN_FILE.exists():
        raise SystemExit(f"\n{DOMAIN_FILE} is missing.")
    dom = gpd.read_file(DOMAIN_FILE).sort_values("domain_id").reset_index(drop=True)
    b = dom.total_bounds
    win = shp_box(b[0] - PAD_W - 2000, b[1] - PAD_S - 2000,
                  b[2] + PAD_E + 2000, b[3] + PAD_N + 2000)
    coast = gpd.read_file(COAST_FILE).to_crs(dom.crs)
    coast = coast[coast.intersects(win)]
    land = unary_union(coast.geometry.values).intersection(win)
    road = gpd.read_file(ROAD_FILE).to_crs(dom.crs)
    road = gpd.clip(road, win)
    print(f"  {len(dom)} domain boxes, {len(coast)} coast polygons, "
          f"{road.geometry.length.sum() / 1000:.1f} km of NC-12 in the window")
    return dom, land, road, b


# Drawing

# Discrete purples for 1..vmax, the family of C["ACCENT"]; 0 is bare land
def count_cmap(vmax):
    cols = [cm.Purples(0.30 + 0.65 * (i - 1) / max(vmax - 1, 1)) for i in range(1, vmax + 1)]
    return cols


# The island with the domain boxes clipped to land
def draw_island(ax, dom, land, road, bounds, fill_of, letter, title,
                reach_labels=True, pad_w=None, alpha_of=None, pad_e=None,
                pad_s=None, pad_n=None, label_every=10, reach_rotation=0):
    b = bounds
    alpha_of = alpha_of or {}
    if pad_w is None:
        pad_w = PAD_W if reach_labels else 1500.0
    pad_e = PAD_E if pad_e is None else pad_e
    pad_s = PAD_S if pad_s is None else pad_s
    pad_n = PAD_N if pad_n is None else pad_n
    ax.set_facecolor(CLR_WATER)
    for geom in getattr(land, "geoms", [land]):
        ax.add_patch(mpatches.Polygon(np.asarray(geom.exterior.coords),
                                      closed=True, facecolor=CLR_LAND,
                                      edgecolor=CLR_COAST, lw=0.5, zorder=1))
    for _, r in dom.iterrows():
        d = int(r["domain_id"])
        g = r.geometry
        on_land = g.intersection(land)
        if on_land.is_empty:
            on_land = g
        fc = fill_of.get(d, "none")
        for part in getattr(on_land, "geoms", [on_land]):
            if part.is_empty or part.geom_type != "Polygon":
                continue
            ax.add_patch(mpatches.Polygon(np.asarray(part.exterior.coords),
                                          closed=True, facecolor=fc,
                                          alpha=alpha_of.get(d, 1.0),
                                          edgecolor=CLR_BOX, lw=0.3, zorder=3))
        if d == 1 or d % label_every == 0:
            x0, y0, x1, y1 = g.bounds
            ax.plot([x1, x1 + 350], [0.5 * (y0 + y1)] * 2, color=CLR_MUTED,
                    lw=0.5, zorder=4)
            ax.text(x1 + 450, 0.5 * (y0 + y1), str(d), fontsize=7,
                    va="center", ha="left", color=CLR_INK, zorder=5)
    road.plot(ax=ax, color=CLR_ROAD, lw=0.7, zorder=4)

    for name, lo, hi, kind in (SECTIONS if reach_labels else []):
        sub = dom[dom["domain_id"].between(lo, hi)]
        if not len(sub):
            continue
        x0, y0, x1, y1 = sub.total_bounds
        yc = 0.5 * (y0 + y1)
        ax.plot([x0 - 250, x0 - 250], [y0 + 60, y1 - 60], color=CLR_MUTED,
                lw=0.6, zorder=4)
        text = name.replace("\n", " ") if reach_rotation else REACH_LABEL.get(name, name)
        ax.text(x0 - 450, yc, text, fontsize=7,
                va="center", ha="right", zorder=5, rotation=reach_rotation,
                linespacing=1.1,
                color=CLR_INK if kind == "village" else CLR_MUTED,
                fontstyle="normal" if kind == "village" else "italic")

    ax.set_xlim(b[0] - pad_w, b[2] + pad_e)
    ax.set_ylim(b[1] - pad_s, b[3] + pad_n)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    _title(ax, ord(letter) - ord("a"), title)


# The house scale bar and north arrow (hat_figure_style), for a map panel without coordinate ticks
def scalebar_and_north(ax, length_m=SCALEBAR_M):
    _scalebar(ax, length_m, show_cells=False)
    _north_arrow(ax)


# The column-head strings of the strip
def strip_heads(obs, obs_idx, storm_names, one_line=False):
    heads = []
    for i in obs_idx:
        r = obs.loc[i]
        head = f"{r['Imagery_Date']:%d %b %Y}" + ("*" if r["poor"] else "")
        names = storm_names.get(i, [])
        if names:
            shown = names[0] + (f" +{len(names) - 1}" if len(names) > 1 else "")
            head = shown + (", " if one_line else chr(10)) + head
        heads.append((head, bool(names)))
    return heads


# One column per image, level with the island beside it
def draw_strips(ax, dom, obs_idx, obs, matrix, domains, storm_names,
                bounds, letter, title, one_line_heads=False):
    b = bounds
    ybox = {int(r["domain_id"]): (r.geometry.bounds[1], r.geometry.bounds[3])
            for _, r in dom.iterrows()}
    n = len(obs_idx)
    ax.set_facecolor("white")
    heads = strip_heads(obs, obs_idx, storm_names, one_line_heads)
    for k, i in enumerate(obs_idx):
        v = matrix[i]
        for j, d in enumerate(domains):
            y0, y1 = ybox[int(d)]
            if np.isnan(v[j]):
                ax.add_patch(mpatches.Rectangle(
                    (k, y0), 1, y1 - y0, facecolor=CLR_PART_FACE,
                    edgecolor=CLR_PART_EDGE, hatch="....", lw=0, zorder=2))
            elif v[j] >= 0.5:
                ax.add_patch(mpatches.Rectangle(
                    (k, y0), 1, y1 - y0, facecolor=CLR_OW, edgecolor="none",
                    lw=0, zorder=2))
        ax.axvline(k, color="#e0e0e0", lw=0.4, zorder=3)
        head, named = heads[k]
        ax.text(k + 0.5, -0.006, head, rotation=90, ha="center", va="top",
                fontsize=7, linespacing=1.15,
                color=CLR_INK if named else CLR_MUTED,
                transform=mtransforms.blended_transform_factory(ax.transData,
                                                                ax.transAxes),
                zorder=5)
    for _, lo, _, _ in SECTIONS[1:]:
        ax.axhline(ybox[lo][0], color="#bfbfbf", lw=0.6, zorder=3)
    for d in domains:
        if d == 1 or d % 10 == 0:
            y0, y1 = ybox[int(d)]
            ax.text(n + 0.25, 0.5 * (y0 + y1), str(d), fontsize=7,
                    va="center", ha="left", color=CLR_MUTED)
    ax.set_xlim(0, n)
    ax.set_ylim(b[1] - PAD_S, b[3] + PAD_N)
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    for s in ax.spines.values():
        s.set_color("#bbbbbb")
    _title(ax, ord(letter) - ord("a"), title)


# Inches the rotated column heads need below the strip
def head_band_in(obs, per, storm_names, tags, one_line=False):
    longest = 0
    for t in tags:
        for head, _ in strip_heads(obs, per[t]["idx"], storm_names, one_line):
            longest = max(longest, max(len(line) for line in head.split(chr(10))))
    return max(1.45, 0.055 * longest + 0.15)


# Caption

# The CAPTIONS.md entry for the figure of the given periods
def caption_text(tags, per):
    one = len(tags) == 1
    span = {"period1": "1984–2004", "period2": "2004–2024"}
    src = {"period1": "the Hapke and Henderson (2007) delineations",
           "period2": "Google Earth imagery"}
    if one:
        t = tags[0]
        letters = "(a) the island, (b) the images"
        imgs = f"{per[t]['n_img']} images in {span[t]}"
        srcs = f"The images are {src[t]}"
        top = f"Most-hit domains: {per[t]['top']}."
        sides = ("domain numbers every ten on the ocean side, reaches on the "
                 "sound side, villages in black")
    else:
        letters = "(a) 1984–2004 and (c) 2004–2024 the island, (b) and (d) the images"
        imgs = (f"{per['period1']['n_img']} images in 1984–2004, "
                f"{per['period2']['n_img']} in 2004–2024, 2004 counted in both")
        srcs = ("Period 1 images are the Hapke and Henderson (2007) "
                "delineations, Period 2 images are Google Earth")
        top = (f"Most-hit domains: {per['period1']['top']} in 1984–2004; "
               f"{per['period2']['top']} in 2004–2024.")
        sides = "domain numbers every ten on the ocean side"
    return (
        f"Observed overwash on Hatteras Island, "
        f"{span[tags[0]] if one else 'by period'}, on the island outline "
        f"(NC 1:80k coastline, UTM 18N); {letters}. The island: the 90 CASCADE "
        f"domain boxes clipped to land, each shaded by the number of images "
        f"showing it overwashed (unshaded = none; {imgs}); the shade scale is "
        f"the same in both periods' figures. NC-12 as the dark line; {sides}. "
        f"The images: one column per image in date order, "
        f"level with the island beside it, purple where the domain was overwashed "
        f"in that image, dotted where the image does not reach the domain; the "
        f"column head under the strip gives the date (* poor or partial image) "
        f"and the named storm that image is the first to show (+n more), by the "
        f"date rule of the heatmap figure. {srcs}. 2004 belongs to both "
        f"periods. {top}")


# Write one figure's entry into CAPTIONS.md
def write_caption(name, text):
    p = upsert_caption(name, "map", text)
    print(f"  wrote {p.relative_to(REPO)}  ({name})")


# Render

# One figure at double-column width: per period, the island panel and its image strip side by side
def render(tags, out_name, geo, per, cols, vmax, obs, matrix, domains,
           storm_names):
    dom, land, road, bounds = geo
    span = {"period1": "1984–2004", "period2": "2004–2024"}
    n_per = len(tags)
    pair = n_per > 1

    yspan = (bounds[3] - bounds[1]) + PAD_S + PAD_N
    pad_w = PAD_W if not pair else 1500.0
    xspan = (bounds[2] - bounds[0]) + pad_w + PAD_E
    ratio = yspan / xspan
    fig_w = figsize("double")[0]
    LM, G, GG, RM, TM, BM = 0.20, 0.10, 0.70, 0.32, 0.40, 0.08
    # Two periods: the column heads run as one line, and the band below grows
    HEAD = head_band_in(obs, per, storm_names, tags, one_line=pair)
    LG = 0.40 if pair else 0.0                        # a legend row, when two periods
    n_img_tot = sum(per[t]["n_img"] for t in tags)
    chrome = LM + RM + n_per * G + GG * (n_per - 1)

    # The island is as tall as the page allows, giving width back if the columns get too thin
    map_h = FIG_H_MAX - BM - LG - TM - HEAD
    map_w = map_h / ratio
    w_col = min(STRIP_IN_MAX, (fig_w - chrome - n_per * map_w) / n_img_tot)
    col_min = COL_MIN_1LINE if pair else COL_MIN_IN
    if w_col < col_min:
        w_col = col_min
        map_w = (fig_w - chrome - n_img_tot * w_col) / n_per
        map_h = map_w * ratio
    strip_w = {t: w_col * per[t]["n_img"] for t in tags}
    fig_h = BM + LG + HEAD + map_h + TM
    fig = plt.figure(figsize=figsize("double", height=fig_h))

    y = BM + LG + HEAD

    def ax_at(x, w):
        return fig.add_axes([x / fig_w, y / fig_h, w / fig_w, map_h / fig_h])

    letters = iter("abcdefgh")
    x = LM
    for k, t in enumerate(tags):
        ax_m = ax_at(x, map_w)
        x += map_w + G
        ax_s = ax_at(x, strip_w[t])
        x += strip_w[t] + GG
        fill_of = {d: cols[n - 1] for d, n in per[t]["counts"].items() if n > 0}
        draw_island(ax_m, dom, land, road, bounds, fill_of, next(letters), span[t],
                    reach_labels=not pair, pad_w=pad_w)
        if k == 0:
            scalebar_and_north(ax_m)
        draw_strips(ax_s, dom, per[t]["idx"], obs, matrix, domains,
                    storm_names, bounds, next(letters),
                    "images" if pair else f"images {span[t]}",
                    one_line_heads=pair)

    strip_ref = "(b)" if len(tags) == 1 else "(b, d)"
    handles = [mpatches.Patch(fc=CLR_LAND, ec=CLR_BOX, lw=0.4,
                              label="domain, no overwash in the period")]
    handles += [mpatches.Patch(fc=cols[i], ec="none",
                               label=f"{i + 1} image{'s' if i else ''}")
                for i in range(vmax)]
    handles += [plt.Line2D([0], [0], color=CLR_ROAD, lw=0.9, label="NC-12"),
                mpatches.Patch(fc=CLR_OW, ec="none",
                               label=f"overwash in that image {strip_ref}"),
                mpatches.Patch(fc=CLR_PART_FACE, ec=CLR_PART_EDGE, hatch="....",
                               lw=0, label="image does not reach the domain")]
    if pair:
        fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
                   handlelength=1.4, handleheight=1.0, columnspacing=1.6,
                   borderaxespad=0.3)
    else:
        # under the island, in the band the column heads leave free there
        fig.legend(handles=handles, loc="upper left",
                   bbox_to_anchor=(LM / fig_w, (BM + HEAD - 0.05) / fig_h),
                   ncol=1, frameon=False, handlelength=1.4, handleheight=1.0,
                   columnspacing=1.4, borderaxespad=0, labelspacing=0.45)

    out = FIG_DIR / out_name
    save(fig, out, close=True)
    print(f"  wrote {out.relative_to(REPO)} (+ .pdf)")
    write_caption(out_name, caption_text(tags, per))


# Main

# Run: one figure per period; --both adds the two-period figure
def main(argv):
    apply_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    obs, domains, matrix = load_observations()
    storms = assign_capture(load_storms(), obs)
    geo = load_geometry()

    storm_names = {}
    for s in storms:
        if s["capture"] is not None:
            nm = s["name"].title() if s["cat"].startswith("H") else s["name"]
            storm_names.setdefault(s["capture"], []).append(nm)

    per = {}
    for tag in ("period1", "period2"):
        lo, hi = PERIODS[tag]
        idx = [int(i) for i in obs.index[obs["Year"].between(lo, hi)]]
        n_ow = np.nansum(matrix[idx], axis=0)
        counts = {int(d): int(n) for d, n in zip(domains, n_ow)}
        order = np.argsort(-n_ow)[:3]
        top = ", ".join(f"domain {int(domains[k])} ({int(n_ow[k])})" for k in order)
        per[tag] = dict(idx=idx, counts=counts, n_img=len(idx),
                        vmax=int(n_ow.max()), top=top)

    vmax = max(per["period1"]["vmax"], per["period2"]["vmax"], 1)
    cols = count_cmap(vmax)
    common = (geo, per, cols, vmax, obs, matrix, domains, storm_names)

    render(["period1"], "overwash_map_period1.png", *common)
    render(["period2"], "overwash_map_period2.png", *common)
    if "--both" in argv:
        render(["period1", "period2"], "overwash_map_periods.png", *common)


if __name__ == "__main__":
    main(sys.argv[1:])
