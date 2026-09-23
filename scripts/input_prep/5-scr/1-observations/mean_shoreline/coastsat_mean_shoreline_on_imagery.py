"""
coastsat_mean_shoreline_on_imagery.py -- the window mean, on that window's photographs
==============================================================================
The CoastSat mean shoreline for a window (coastsat_mean_shoreline.py) drawn
on the USGS aerial photographs flown inside that window, one panel per
flight year, at a handful of sites along the island. Asked for by Hannah on
2026-09-23 ("the shoreline position imposed on the aerial imagery from that
year").

WHAT IS ON EACH PANEL
    the photograph      the USGS Henderson release (doi 10.5066/P1CXBCDW),
                        read frame by frame from D:\\Hatteras_GIS\\Aerial
                        through the 1984 imagery review's Imagery class (the
                        seamline rule and the film-fringe handling live there,
                        not here); stated accuracy 1.2 m
    the mean shoreline  the window-mean line, ink with a white halo
    +/-1 sd             the within-window scatter of each transect, placed
                        along the transect's own direction and joined
                        alongshore, as a translucent white band, dashed edges
    the positions       (second version only) every satellite position behind
                        the mean, geolocated the same way (origin + chainage *
                        direction), coloured by date on one scale for the
                        whole window; the scale marks each flight date
                        (added 2026-09-23, Hannah: "with and without the dots",
                        "a gradient to show throughout time")

WHAT IT CAN AND CANNOT SHOW
    A photograph is ONE October day; the line is a three-year mean of ~28
    satellite passes. The wet/dry line in a photo is not expected to sit on
    the mean -- it is expected to sit, most days, inside the band. A photo
    edge far outside the band at one site is worth looking at; nothing is
    measured from the photographs here.

SITES (supporting/sites.csv)
    Each window is three domains (1.5 km) alongshore, centred on the site's
    domain, and every panel of a site shares one extent. Cross-shore it runs
    LAND_M landward and SEA_M seaward of the line. The sites are spread along
    the island and include the two piers, which are fixed in all three photos.

OUTPUT   <mean_shoreline_dir(window)>/on_imagery/
    line_and_band/mean_shoreline_<window>_on_imagery_GIS<NN>_<site>.png
    line_and_band/mean_shoreline_<window>_on_imagery_island_1996.png
        the whole island in three north-up segments (GIS 1-30, 31-60, 61-90)
        at one scale on the 1996 photographs, the site windows outlined
    line_and_band/mean_shoreline_<window>_on_imagery_ribbon_1996.png
        panel (a) of mean_shoreline_<window>.png alone on the 1996 photographs,
        the island outline beneath where no frame covers
    with_positions/mean_shoreline_<window>_on_imagery_with_positions_GIS<NN>_<site>.png
        each subfolder with supporting/CAPTIONS.md (no PDFs: raster panels)
    supporting/sites.csv   the windows, both files and the position counts per site
    Also published to output/figures/shoreline/mean_shoreline/<the same subfolders>.

USAGE
    python coastsat_mean_shoreline_on_imagery.py
    python coastsat_mean_shoreline_on_imagery.py --window 1995 1997 --sites 26 79
    python coastsat_mean_shoreline_on_imagery.py --only island   (or sites, ribbon)
    Needs the D: drive; the .venv Python (rasterio).
==============================================================================
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401

import geopandas as gpd  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    DOMAIN_BOXES, mean_shoreline_csv, mean_shoreline_dir,
)

CRS = "EPSG:26918"                        # the line's CRS; the photos are warped to it
DEFAULT_WINDOW = (1995, 1997)
# (centre GIS domain, name). Spread south to north; 26 and 79 are the piers.
DEFAULT_SITES = [(4, "buxton"), (26, "avon_pier"), (40, "central"),
                 (55, "central_north"), (79, "rodanthe_pier"), (88, "mirlo_beach")]
HALF = 1                                  # domains either side of the centre
LAND_M, SEA_M = 300.0, 250.0              # window either side of the line
RES_M = 0.6                               # display resolution of the photographs
C_LINE = fs.INK
# The +/-1 sd band is neutral -- translucent white with a dashed ink edge -- so
# the only colour on a panel is the positions' date scale. It was yellow until
# 2026-09-23, when the dates moved to viridis, whose light end is yellow.
C_BAND = "white"
BAND_EDGE = dict(color=C_LINE, lw=0.5, ls=(0, (3, 2)), alpha=0.9)
# The positions' date scale (Hannah, 2026-09-23: "academic and professional"):
# viridis is perceptually uniform, reads in greyscale and to colour-blind
# readers, and is the scale reviewers expect for an ordered variable. The thin
# ink edge keeps the light end visible on the brightest beach.
CMAP = plt.get_cmap("viridis")
PUBLISH = fs.figure_dir("shoreline", "mean_shoreline")
# One subfolder per version, each with its own supporting/CAPTIONS.md
# (Hannah, 2026-09-23); sites.csv covers both and stays in on_imagery/supporting.
VERSION_DIRS = {False: "line_and_band", True: "with_positions"}

REVIEW_SCRIPT = (_REPO / "scripts" / "input_prep" / "1-barrier3d-domains"
                 / "2-domain-reconstruction-1984" / "3-placement" / "imagery-review"
                 / "HAT_imagery_review_1984.py")


def _review_module():
    """The 1984 imagery review, for its Imagery reader (one copy of it)."""
    spec = importlib.util.spec_from_file_location("imagery_review", REVIEW_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def transects(window):
    """Included transects with their direction, in alongshore order."""
    import coastsat_mean_shoreline as cms
    df = pd.read_csv(mean_shoreline_csv(*window))
    df = df[df["included"] & df["x"].notna()].copy()
    geom = cms.transect_geometry(set(df["transect_id"]))
    df["x0"] = df["transect_id"].map(lambda t: geom[t][0])
    df["y0"] = df["transect_id"].map(lambda t: geom[t][1])
    df["ux"] = df["transect_id"].map(lambda t: geom[t][2])
    df["uy"] = df["transect_id"].map(lambda t: geom[t][3])
    return df.sort_values(["site", "transect_number"]).reset_index(drop=True)


def positions(df, window):
    """Every satellite position behind the included means, on the ground.

    The same window filter as coastsat_mean_shoreline.window_means, so these
    are exactly the positions each mean was taken over.
    """
    import coastsat_mean_shoreline as cms
    from coastsat_lrr import filter_dates, load_timeseries
    lo, hi = "{0}-01-01".format(window[0]), "{0}-12-31".format(window[1])
    out = []
    for t in df.itertuples():
        obs = filter_dates(load_timeseries(str(cms.timeseries_file(t.transect_id))), lo, hi)
        obs = obs.dropna(subset=["chainage_m"])
        ch = obs["chainage_m"].to_numpy(dtype=float)
        out.append(pd.DataFrame({"transect_id": t.transect_id,
                                 "date": obs["date"].dt.date.values,
                                 "year": obs["date"].dt.year.values,
                                 "chainage_m": ch,
                                 "x": t.x0 + ch * t.ux, "y": t.y0 + ch * t.uy}))
    return pd.concat(out, ignore_index=True)


def draw_positions(ax, pos, b, norm):
    """Every position in the window, coloured by its date on one shared scale."""
    x0, y0, x1, y1 = b
    p = pos[pos["x"].between(x0, x1) & pos["y"].between(y0, y1)].sort_values("t")
    ax.scatter(p["x"], p["y"], s=6, c=p["t"], cmap=CMAP, norm=norm,
               edgecolors=C_LINE, linewidths=0.3, zorder=5)
    return p


def extent(df, boxes, centre):
    """North-up window: the three domains alongshore, the line +/- cross-shore."""
    seg = boxes[boxes["gis"].between(centre - HALF, centre + HALF)]
    y0, y1 = seg.total_bounds[1], seg.total_bounds[3]
    pts = df[(df["y"] >= y0) & (df["y"] <= y1)]
    # seaward is +x on this coast (the transects point offshore, ux > 0)
    return (pts["x"].min() - LAND_M, y0, pts["x"].max() + SEA_M, y1)


def draw_line(ax, df, b, lw=1.3, edges=True):
    """The mean line and its band. At island scale (edges=False) the band's
    dashed edges would merge with the line, so only its fill is drawn."""
    x0, y0, x1, y1 = b
    pad = 100.0
    s = df[(df["y"] >= y0 - pad) & (df["y"] <= y1 + pad)]
    sea = np.c_[s["x"] + s["sd_chainage_m"] * s["ux"], s["y"] + s["sd_chainage_m"] * s["uy"]]
    land = np.c_[s["x"] - s["sd_chainage_m"] * s["ux"], s["y"] - s["sd_chainage_m"] * s["uy"]]
    ax.fill(np.r_[sea[:, 0], land[::-1, 0]], np.r_[sea[:, 1], land[::-1, 1]],
            color=C_BAND, alpha=0.35, lw=0, zorder=2)
    if edges:
        for e in (sea, land):
            ax.plot(e[:, 0], e[:, 1], zorder=2, **BAND_EDGE)
    ax.plot(s["x"], s["y"], color=C_LINE, lw=lw, zorder=3,
            path_effects=[pe.withStroke(linewidth=lw + 1.7, foreground="white")])


def read_photos(imagery, b):
    """Each year's photograph for one window, read once for both versions."""
    x0, y0, x1, y1 = b
    out = []
    for im in imagery:
        img = im.read(x0, x1, y0, y1, RES_M, CRS).copy()
        blank = img.max(axis=2) == 0
        img[blank] = 255                                   # no photograph -> white
        out.append((img, float(blank.mean())))
    return out


def panel(ax, im, img, b, df, i):
    x0, y0, x1, y1 = b
    ax.imshow(img, extent=(x0, x1, y0, y1), origin="upper", zorder=0,
              interpolation="bilinear")
    draw_line(ax, df, b)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fs.spines_for_image(ax)
    d = pd.Timestamp(im.date)
    ax.set_title(f"({chr(97 + i)})  {d.day} {d:%B %Y}", loc="left", fontsize=9, pad=4)


def _decimal_year(ts):
    ts = pd.Timestamp(ts)
    start = pd.Timestamp(year=ts.year, month=1, day=1)
    return ts.year + (ts - start).days / (366 if ts.is_leap_year else 365)


def time_bar(fig, axes, norm, imagery, window):
    """The date scale, with each photograph's flight marked by its panel letter."""
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    cb = fig.colorbar(sm, ax=axes, location="bottom", shrink=0.45, aspect=40, pad=0.02)
    cb.set_ticks(list(range(window[0], window[1] + 2)))
    cb.ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v)}"))
    cb.set_label("Date of satellite position")
    cb.outline.set_linewidth(0.5)
    tr = cb.ax.get_xaxis_transform()
    for i, im in enumerate(imagery):
        t = _decimal_year(im.date)
        cb.ax.plot([t], [1.0], marker="v", ms=5, color=C_LINE, clip_on=False,
                   transform=tr, zorder=10)
        cb.ax.text(t, 1.4, f"({chr(97 + i)})", ha="center", va="bottom", fontsize=7.5,
                   transform=tr)
    return cb


def site_figure(df, pos, boxes, imagery, centre, key, window, out_dir, norm):
    """Two versions of one site: the line and band, then the same with the positions."""
    b = extent(df, boxes, centre)
    photos = read_photos(imagery, b)
    w, h = b[2] - b[0], b[3] - b[1]
    n = len(imagery)
    start, end = window
    dates = ", ".join(pd.Timestamp(im.date).strftime("%Y-%m-%d") for im in imagery)
    row = dict(centre_gis=centre, site=key, first_gis=centre - HALF, last_gis=centre + HALF,
               x0=b[0], y0=b[1], x1=b[2], y1=b[3],
               **{f"no_photo_frac_{im.year}": round(f, 3) for im, (_, f) in zip(imagery, photos)})
    for dots in (False, True):
        panel_w = fs.FIG_W_DOUBLE / n * 0.92
        fig, axes = plt.subplots(1, n, layout="constrained",
                                 figsize=(fs.FIG_W_DOUBLE, panel_w * h / w + (1.75 if dots else 0.9)))
        axes = np.atleast_1d(axes)
        for i, (ax, im, (img, _)) in enumerate(zip(axes, imagery, photos)):
            panel(ax, im, img, b, df, i)
            if dots:
                shown = draw_positions(ax, pos, b, norm)
        fs._scalebar(axes[0], 200.0, show_cells=False)
        fs._north_arrow(axes[0], x=0.86, y=0.06)
        handles = [Line2D([], [], color=C_LINE, lw=1.3,
                          path_effects=[pe.withStroke(linewidth=3.0, foreground="white")]),
                   Patch(facecolor="0.9", edgecolor=C_LINE, lw=0.5, ls=(0, (3, 2)))]
        labels = [f"Mean shoreline, {start}–{end} (CoastSat)", "±1 standard deviation"]
        if dots:
            time_bar(fig, list(axes), norm, imagery, window)
            handles.append(Line2D([], [], ls="none", marker="o", ms=3.2,
                                  mfc=CMAP(0.5), mec=C_LINE, mew=0.4))
            labels.append("Individual satellite positions")
        fig.legend(handles, labels, loc="outside lower center", ncol=len(handles), frameon=False)
        fig.suptitle(f"GIS {centre - HALF}–{centre + HALF}", fontsize=9.5)
        what = "_with_positions" if dots else ""
        stem = f"mean_shoreline_{start}_{end}_on_imagery{what}_GIS{centre:02d}_{key}"
        dot_text = (" Dots: every satellite position the mean was taken over, geolocated along "
                    "its transect and coloured by date on the scale below the panels; the "
                    "triangles on the scale mark each photograph's flight date."
                    if dots else "")
        fs.caption(fig, (
            f"The CoastSat mean shoreline for calendar {start}–{end} (ink, white halo) and "
            f"the ±1 standard deviation of the satellite positions behind each transect mean "
            f"(translucent white band with dashed edges, placed along each transect's direction), on the USGS aerial "
            f"photographs flown inside the window ({dates}; doi 10.5066/P1CXBCDW, stated "
            f"accuracy 1.2 m), GIS domains {centre - HALF}–{centre + HALF}, north up, all "
            f"panels at one extent and scale.{dot_text} Each photograph is one autumn day; "
            f"the line is a mean over the window, so the waterline in a photograph is expected "
            f"to fall within the band, not on the line. White is outside the photographs."))
        sub = VERSION_DIRS[dots]
        path = out_dir / sub / f"{stem}.png"
        fs.save(fig, path, vector=False, close=True)
        (PUBLISH / sub).mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, PUBLISH / sub / path.name)
        row["file_with_positions" if dots else "file"] = f"{sub}/{path.name}"
        if dots:
            row.update({f"positions_{y}": int((shown["year"] == y).sum())
                        for y in range(start, end + 1)}, positions_total=len(shown))
    return row


# =============================================================================
# the island overview (Hannah, 2026-09-23: "across the island in 3 vertical
# panels with the 1996 imagery")
# =============================================================================

ISLAND_YEAR = 1996
SEGMENTS = [("GIS 1–30", 1, 30), ("GIS 31–60", 31, 60), ("GIS 61–90", 61, 90)]
ISLAND_RES_M = 5.0                        # ~ the printed pixel at this scale
ISLAND_LAND_M, ISLAND_SEA_M = 900.0, 600.0
ISLAND_PANEL_H_IN = 8.2                   # the segments are 15 km tall; this sets the scale


def island_figure(df, boxes, im, sites, window, out_dir):
    """Three north-up segments side by side at one scale, on one year's photos,
    with the site windows of the zoom figures outlined."""
    ext = []
    for _, lo, hi in SEGMENTS:
        seg = boxes[boxes["gis"].between(lo, hi)]
        y0, y1 = seg.total_bounds[1], seg.total_bounds[3]
        pts = df[df["y"].between(y0, y1)]
        ext.append((pts["x"].min() - ISLAND_LAND_M, y0, pts["x"].max() + ISLAND_SEA_M, y1))
    widths = [b[2] - b[0] for b in ext]
    height = max(b[3] - b[1] for b in ext)
    scale = min(ISLAND_PANEL_H_IN / height, (fs.FIG_W_DOUBLE - 1.3) / sum(widths))   # in per m
    fig = plt.figure(figsize=(fs.FIG_W_DOUBLE, height * scale + 1.0), layout="constrained")
    axes = fig.subplots(1, 3, gridspec_kw={"width_ratios": widths})
    for i, (ax, b, (title, lo, hi)) in enumerate(zip(axes, ext, SEGMENTS)):
        x0, y0, x1, y1 = b
        print(f"    {title} ...")
        img = im.read(x0, x1, y0, y1, ISLAND_RES_M, CRS).copy()
        img[img.max(axis=2) == 0] = 255
        ax.imshow(img, extent=(x0, x1, y0, y1), origin="upper", zorder=0,
                  interpolation="bilinear")
        draw_line(ax, df, b, lw=0.8, edges=False)
        for _, s in sites.iterrows():
            if s["y1"] > y0 and s["y0"] < y1:
                ax.add_patch(matplotlib.patches.Rectangle(
                    (s["x0"], s["y0"]), s["x1"] - s["x0"], s["y1"] - s["y0"], fill=False,
                    ec=C_LINE, lw=0.7, zorder=4,
                    path_effects=[pe.withStroke(linewidth=1.8, foreground="white")]))
                ax.text(s["x0"] + 60, s["y1"] - 60,
                        f"{int(s['first_gis'])}–{int(s['last_gis'])}", fontsize=6.5,
                        va="top", ha="left", color=C_LINE, zorder=5, clip_on=True,
                        bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                                  boxstyle="square,pad=0.15"))
        seg = boxes[boxes["gis"].between(lo, hi)]
        ticks = seg[(seg["gis"] % 5 == 0) | (seg["gis"] == lo)]
        ax.set_yticks([(g.bounds[1] + g.bounds[3]) / 2 for g in ticks.geometry])
        ax.set_yticklabels([str(g) for g in ticks["gis"]], fontsize=7)
        ax.tick_params(axis="y", length=2.5, pad=1.5)
        ax.set_xticks([])
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_aspect("equal")
        fs.spines_for_image(ax)
        ax.set_title(f"({chr(97 + i)})  {title}", loc="left", fontsize=9, pad=4)
    axes[0].set_ylabel(fs.DOMAIN_AXIS_LABEL)
    fs._scalebar(axes[0], 2000.0, show_cells=False)
    fs._north_arrow(axes[0], x=0.80, y=0.07, length=0.03)
    start, end = window
    fig.legend([Line2D([], [], color=C_LINE, lw=0.8,
                       path_effects=[pe.withStroke(linewidth=2.5, foreground="white")]),
                Patch(facecolor="0.9", edgecolor="none"),
                Patch(facecolor="none", edgecolor=C_LINE, lw=0.7)],
               [f"Mean shoreline, {start}–{end} (CoastSat)", "±1 standard deviation",
                "Site windows (zoom figures)"],
               loc="outside lower center", ncol=3, frameon=False)
    d = pd.Timestamp(im.date)
    fig.suptitle(f"USGS aerial photographs, {d.day} {d:%B %Y}", fontsize=9.5)
    fs.caption(fig, (
        f"The CoastSat mean shoreline for calendar {start}–{end} (ink, white halo) and the "
        f"±1 standard deviation of the satellite positions behind each transect mean (white "
        f"band; ~25 m wide, so near the limit of the print at this scale) on the USGS aerial "
        f"photographs of {d:%Y-%m-%d} (doi 10.5066/P1CXBCDW, stated accuracy 1.2 m), the island "
        f"in three north-up segments at one scale: (a) GIS 1–30, Cape Point to Avon; (b) GIS "
        f"31–60; (c) GIS 61–90, the Tri-Village to the north end. Ticks mark every fifth "
        f"Barrier3D domain. Boxes outline the six site windows drawn at full resolution in "
        f"the zoom figures, labelled with their domains. White is outside the photographs."))
    path = out_dir / VERSION_DIRS[False] / f"mean_shoreline_{start}_{end}_on_imagery_island_{im.year}.png"
    fs.save(fig, path, vector=False, close=True)
    (PUBLISH / VERSION_DIRS[False]).mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, PUBLISH / VERSION_DIRS[False] / path.name)
    print(f"  island figure -> {path.name}")


RIBBON_RES_M = 8.0                        # ~ the printed pixel of a 45 km ribbon


def ribbon_figure(df, im, window, out_dir):
    """Panel (a) of the diagnostic, alone, over one year's photographs: the
    same extent and axes as coastsat_mean_shoreline.outline_figure."""
    import coastsat_mean_shoreline as cms
    start, end = window
    ext = cms.ribbon_extent(df)
    n0, n1, e0, e1 = ext
    print("    ribbon ...")
    img = im.read(e0, e1, n0, n1, RIBBON_RES_M, CRS)
    # transparent where no frame covers, so the island outline shows through
    img = np.dstack([img, np.where(img.max(axis=2) == 0, 0, 255).astype(np.uint8)])
    # rows N (top = north) x cols E  ->  rows E (top = east) x cols N
    img = np.transpose(img[::-1], (1, 0, 2))[::-1]
    w = fs.FIG_W_DOUBLE
    fig, ax = plt.subplots(figsize=(w, (e1 - e0) / (n1 - n0) * (w - 0.8) + 1.1),
                           layout="constrained")
    cms.draw_island_outline(ax, ext, edge_on_top="white")
    ax.imshow(img, extent=(n0, n1, e0, e1), origin="upper", zorder=2, interpolation="nearest")
    ax.plot(df["y"], df["x"], "-", color=C_LINE, lw=0.8, zorder=4,
            path_effects=[pe.withStroke(linewidth=2.3, foreground="white")])
    cms.ribbon_axes(ax, ext)
    d = pd.Timestamp(im.date)
    ax.set_title(f"Mean shoreline, {start}–{end}, on the {d.day} {d:%B %Y} photographs")
    fig.legend([Line2D([], [], color=C_LINE, lw=0.8,
                       path_effects=[pe.withStroke(linewidth=2.3, foreground="white")]),
                Line2D([], [], color="white", lw=0.8,
                       path_effects=[pe.withStroke(linewidth=1.8, foreground=fs.INK_MUTED)])],
               [f"Mean shoreline, {start}–{end} (CoastSat)", "Island outline"],
               loc="outside lower center", ncol=2, frameon=False)
    fs.caption(fig, (
        f"The CoastSat mean shoreline for calendar {start}–{end} (ink, white halo) on the USGS "
        f"aerial photographs of {d:%Y-%m-%d} (doi 10.5066/P1CXBCDW), drawn as panel (a) of "
        f"mean_shoreline_{start}_{end}.png but flipped: alongshore across the page, easting "
        f"increasing downward, equal aspect, so the ocean is at the bottom. The ±1 standard deviation band (~25 m) is below "
        f"the resolution of the print and is not drawn. The island outline "
        f"(map_elements/hatteras_outline; its survey date is not recorded) is drawn in white "
        f"over the photographs, and as grey land and pale-blue water where no frame covers."))
    path = out_dir / VERSION_DIRS[False] / f"mean_shoreline_{start}_{end}_on_imagery_ribbon_{im.year}.png"
    fs.save(fig, path, vector=False, close=True)
    shutil.copy2(path, PUBLISH / VERSION_DIRS[False] / path.name)
    print(f"  ribbon figure -> {path.name}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--window", nargs=2, type=int, default=DEFAULT_WINDOW)
    ap.add_argument("--sites", nargs="*", type=int,
                    help="centre GIS domains (default: the six built-in sites)")
    ap.add_argument("--only", choices=("sites", "island", "ribbon"),
                    help="draw only the site figures, the three-segment island "
                         "overview, or the single-panel ribbon")
    a = ap.parse_args(argv)
    window = tuple(a.window)
    sites = ([(c, k) for c, k in DEFAULT_SITES if c in a.sites] + [(c, "site") for c in a.sites
              if c not in dict(DEFAULT_SITES)]) if a.sites else DEFAULT_SITES

    fs.apply_style()
    R = _review_module()
    imagery = []
    for y in range(window[0], window[1] + 1):
        if R._year_files(y):
            imagery.append(R.Imagery(y, CRS))
    if not imagery:
        raise SystemExit(f"no photographs for {window} under {R.AERIAL_ROOT} (is D: on?)")

    df = transects(window)
    pos = positions(df, window)
    pos["t"] = [_decimal_year(d) for d in pos["date"]]
    norm = matplotlib.colors.Normalize(window[0], window[1] + 1)
    print(f"  {len(pos)} satellite positions over {pos['transect_id'].nunique()} transects")
    boxes = gpd.read_file(DOMAIN_BOXES).to_crs(CRS)
    boxes["gis"] = np.arange(1, len(boxes) + 1)   # file order is south -> north, GIS 1-90

    out_dir = mean_shoreline_dir(*window) / "on_imagery"
    out_dir.mkdir(parents=True, exist_ok=True)
    table = fs.support_dir(out_dir) / "sites.csv"
    if a.only in (None, "sites"):
        rows = []
        for c, k in sites:
            print(f"  GIS {c} ({k}) ...")
            rows.append(site_figure(df, pos, boxes, imagery, c, k, window, out_dir, norm))
        new = pd.DataFrame(rows)
        if a.sites and table.is_file():            # a partial run keeps the other sites' rows
            old = pd.read_csv(table)
            new = pd.concat([old[~old["centre_gis"].isin(new["centre_gis"])], new])
        new.sort_values("centre_gis").to_csv(table, index=False)
        print(f"  {len(rows)} site figure(s) -> {out_dir}")
    if a.only in (None, "island", "ribbon"):
        im = next((im for im in imagery if im.year == ISLAND_YEAR), None)
        if im is None:
            print(f"  no {ISLAND_YEAR} photographs in this window; island figures skipped")
        else:
            if a.only in (None, "island"):
                island_figure(df, boxes, im, pd.read_csv(table), window, out_dir)
            if a.only in (None, "ribbon"):
                ribbon_figure(df, im, window, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
