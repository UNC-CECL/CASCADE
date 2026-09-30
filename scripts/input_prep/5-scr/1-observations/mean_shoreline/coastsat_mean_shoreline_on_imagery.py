"""
A window's CoastSat mean shoreline on that window's aerial photographs, one panel per flight year.

    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_on_imagery.py
    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_on_imagery.py --centred-on alace_1996
    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_on_imagery.py --window 2009 2011 --photo-years 2008

Site zooms, an island view and the ribbon; needs the D: drive for the photos. Details: scripts/input_prep/5-scr/1-observations/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
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

# --- CONFIG ------------------------------------------------------------------
CRS = "EPSG:26918"                        # the line's CRS; the photos are warped to it
DEFAULT_WINDOW = (1995, 1997)
# (centre GIS domain, name). Spread south to north; 26 and 79 are the piers.
DEFAULT_SITES = [(4, "buxton"), (26, "avon_pier"), (40, "central"),
                 (55, "central_north"), (79, "rodanthe_pier"), (88, "mirlo_beach")]
HALF = 1                                  # domains either side of the centre
LAND_M, SEA_M = 300.0, 250.0              # window either side of the line
RES_M = 0.6                               # display resolution of the photographs
C_LINE = fs.INK
# The +/-1 sd band is neutral -- translucent white with a dashed ink edge
C_BAND = "white"
BAND_EDGE = dict(color=C_LINE, lw=0.5, ls=(0, (3, 2)), alpha=0.9)
# The positions' date scale (Hannah, 2026-09-23
CMAP = plt.get_cmap("viridis")
PUBLISH = fs.figure_dir("observations", "mean_shoreline")
# One subfolder per version, each with its own supporting/CAPTIONS.md (Hannah, 2026-09-23)
VERSION_DIRS = {"line": "line_and_band", "positions": "with_positions",
                "domains": "with_domains"}
# The domain boxes follow the island outline's convention in the ribbon figure (white, grey edge)
DOMAIN_STYLE = dict(color="white", zorder=4,
                    path_effects=[pe.withStroke(linewidth=1.9, foreground=fs.INK_MUTED)])

# Photographs that are not the USGS Henderson release
USGS_REF = "doi 10.5066/P1CXBCDW, stated accuracy 1.2 m"
PHOTO_SOURCES = {
    # D:\Hatteras_GIS\Aerial\2008, 2008_IOCM_NaturalColorImagery_J1129187_metadata.xml
    2008: dict(date="2008-03-26", label="26–27 March 2008", short="NOAA NGS orthomosaic",
               ref="the NOAA NGS colour orthomosaic flown 2008-03-26/27 (InPort 48695)"),
}
# -----------------------------------------------------------------------------


# The flight date as a panel title reads it
def photo_label(im):
    if im.year in PHOTO_SOURCES:
        return PHOTO_SOURCES[im.year]["label"]
    d = pd.Timestamp(im.date)
    return f"{d.day} {d:%B %Y}"


# Was this photograph flown inside the window? By date since 2026-09-29, when windows stopped being ...
def _inside(im, window):
    d = PHOTO_SOURCES.get(im.year, {}).get("date", im.date)
    try:
        t = pd.Timestamp(d)
    except (TypeError, ValueError):
        return window.lo.year <= im.year <= window.hi.year
    if len(str(d)) == 4:                    # a bare year: inside if the year is
        return window.lo.year <= im.year <= window.hi.year
    return window.lo <= t <= window.hi


# 'the USGS aerial photographs 
def photo_ref(imagery, window):
    usgs = [im for im in imagery if im.year not in PHOTO_SOURCES]
    parts = []
    if usgs:
        dates = ", ".join(pd.Timestamp(im.date).strftime("%Y-%m-%d") for im in usgs)
        inside = all(_inside(im, window) for im in usgs)
        parts.append(f"the USGS aerial photographs flown inside the window ({dates}; {USGS_REF})"
                     if inside and len(usgs) > 1 else
                     f"the USGS aerial photographs of {dates} ({USGS_REF})")
    parts += [PHOTO_SOURCES[im.year]["ref"] for im in imagery if im.year in PHOTO_SOURCES]
    return " and ".join(parts)


# The sentence on what a photograph can show against a window mean
def photo_timing(imagery, window):
    outside = [im for im in imagery if not _inside(im, window)]
    if not outside:
        return ("Each photograph is one autumn day; the line is a mean over the window, so "
                "the waterline in a photograph is expected to fall within the band, not on "
                "the line.")
    names = ", ".join(photo_label(im) for im in outside)
    return (f"No georeferenced photographs were flown inside the window; {names} is the "
            f"nearest, outside it, so the waterline can sit off the band through real change "
            f"over that gap as well as through the day's water level and waves.")


REVIEW_SCRIPT = (_REPO / "scripts" / "input_prep" / "1-barrier3d-domains"
                 / "2-domain-reconstruction-1984" / "3-placement" / "imagery-review"
                 / "HAT_imagery_review_1984.py")


# The 1984 imagery review, for its Imagery reader (one copy of it)
def _review_module():
    spec = importlib.util.spec_from_file_location("imagery_review", REVIEW_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# Included transects with their direction, in alongshore order
def transects(window):
    import coastsat_mean_shoreline as cms
    df = pd.read_csv(mean_shoreline_csv(*window.key))
    df = df[df["included"] & df["x"].notna()].copy()
    geom = cms.transect_geometry(set(df["transect_id"]))
    df["x0"] = df["transect_id"].map(lambda t: geom[t][0])
    df["y0"] = df["transect_id"].map(lambda t: geom[t][1])
    df["ux"] = df["transect_id"].map(lambda t: geom[t][2])
    df["uy"] = df["transect_id"].map(lambda t: geom[t][3])
    return df.sort_values(["site", "transect_number"]).reset_index(drop=True)


# Every satellite position behind the included means, on the ground
def positions(df, window):
    import coastsat_mean_shoreline as cms
    from coastsat_lrr import load_timeseries
    out = []
    for t in df.itertuples():
        obs = window.clip(load_timeseries(str(cms.timeseries_file(t.transect_id))))
        obs = obs.dropna(subset=["chainage_m"])
        ch = obs["chainage_m"].to_numpy(dtype=float)
        out.append(pd.DataFrame({"transect_id": t.transect_id,
                                 "date": obs["date"].dt.date.values,
                                 "year": obs["date"].dt.year.values,
                                 "chainage_m": ch,
                                 "x": t.x0 + ch * t.ux, "y": t.y0 + ch * t.uy}))
    return pd.concat(out, ignore_index=True)


# Every position in the window, coloured by its date on one shared scale
def draw_positions(ax, pos, b, norm):
    x0, y0, x1, y1 = b
    p = pos[pos["x"].between(x0, x1) & pos["y"].between(y0, y1)].sort_values("t")
    ax.scatter(p["x"], p["y"], s=6, c=p["t"], cmap=CMAP, norm=norm,
               edgecolors=C_LINE, linewidths=0.3, zorder=5)
    return p


# North-up window
def extent(df, boxes, centre):
    seg = boxes[boxes["gis"].between(centre - HALF, centre + HALF)]
    y0, y1 = seg.total_bounds[1], seg.total_bounds[3]
    pts = df[(df["y"] >= y0) & (df["y"] <= y1)]
    # seaward is +x on this coast (the transects point offshore, ux > 0)
    return (pts["x"].min() - LAND_M, y0, pts["x"].max() + SEA_M, y1)


# The mean line and its band
def draw_line(ax, df, b, lw=1.3, edges=True):
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


# The Barrier3D domain boxes crossing the window, each labelled at its landward (west) side, clear of ...
def draw_domains(ax, boxes, b, lw=0.8, labels=True):
    x0, y0, x1, y1 = b
    bd = boxes.bounds
    seg = boxes[(bd["maxy"] > y0) & (bd["miny"] < y1)]
    for gis, g in zip(seg["gis"], seg.geometry):
        xs, ys = g.exterior.xy
        ax.plot(xs, ys, lw=lw, **DOMAIN_STYLE)
        if labels:
            gy0, gy1 = max(g.bounds[1], y0), min(g.bounds[3], y1)
            ax.text(x0 + 0.04 * (x1 - x0), (gy0 + gy1) / 2, f"GIS {gis}", fontsize=7.5,
                    va="center", ha="left", color=C_LINE, zorder=5, clip_on=True,
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                              boxstyle="square,pad=0.2"))


# Legend handle for the domain boxes
def _domain_handle(lw=0.8):
    return Line2D([], [], lw=lw, **{k: v for k, v in DOMAIN_STYLE.items() if k != "zorder"})


# Each year's photograph for one window, read once for both versions
def read_photos(imagery, b):
    x0, y0, x1, y1 = b
    out = []
    for im in imagery:
        img = im.read(x0, x1, y0, y1, RES_M, CRS).copy()
        blank = img.max(axis=2) == 0
        img[blank] = 255                                   # no photograph -> white
        out.append((img, float(blank.mean())))
    return out


# One photo panel with the line on it
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
    ax.set_title(f"({chr(97 + i)})  {photo_label(im)}", loc="left", fontsize=9, pad=4)


# A timestamp as a decimal year
def _decimal_year(ts):
    ts = pd.Timestamp(ts)
    start = pd.Timestamp(year=ts.year, month=1, day=1)
    return ts.year + (ts - start).days / (366 if ts.is_leap_year else 365)


# The date scale, with each photograph's flight marked by its panel letter
def time_bar(fig, axes, norm, imagery, window):
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=norm)
    cb = fig.colorbar(sm, ax=axes, location="bottom", shrink=0.45, aspect=40, pad=0.02)
    cb.set_ticks([y for y in range(int(norm.vmin), int(norm.vmax) + 1)
                  if norm.vmin <= y <= norm.vmax])
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


# Three versions of one site
def site_figure(df, pos, boxes, imagery, centre, key, window, out_dir, norm):
    b = extent(df, boxes, centre)
    photos = read_photos(imagery, b)
    w, h = b[2] - b[0], b[3] - b[1]
    n = len(imagery)
    row = dict(centre_gis=centre, site=key, first_gis=centre - HALF, last_gis=centre + HALF,
               x0=b[0], y0=b[1], x1=b[2], y1=b[3],
               **{f"no_photo_frac_{im.year}": round(f, 3) for im, (_, f) in zip(imagery, photos)})
    for version, sub in VERSION_DIRS.items():
        dots, doms = version == "positions", version == "domains"
        panel_w = fs.FIG_W_DOUBLE / n * 0.92
        fig, axes = plt.subplots(1, n, layout="constrained",
                                 figsize=(fs.FIG_W_DOUBLE, panel_w * h / w + (1.75 if dots else 0.9)))
        axes = np.atleast_1d(axes)
        for i, (ax, im, (img, _)) in enumerate(zip(axes, imagery, photos)):
            panel(ax, im, img, b, df, i)
            if dots:
                shown = draw_positions(ax, pos, b, norm)
            if doms:
                draw_domains(ax, boxes, b)
        fs._scalebar(axes[0], 200.0, show_cells=False)
        fs._north_arrow(axes[0], x=0.86, y=0.06)
        handles = [Line2D([], [], color=C_LINE, lw=1.3,
                          path_effects=[pe.withStroke(linewidth=3.0, foreground="white")]),
                   Patch(facecolor="0.9", edgecolor=C_LINE, lw=0.5, ls=(0, (3, 2)))]
        labels = [f"Mean shoreline, {window.span} (CoastSat)", "±1 standard deviation"]
        if dots:
            time_bar(fig, list(axes), norm, imagery, window)
            handles.append(Line2D([], [], ls="none", marker="o", ms=3.2,
                                  mfc=CMAP(0.5), mec=C_LINE, mew=0.4))
            labels.append("Individual satellite positions")
        if doms:
            handles.append(_domain_handle())
            labels.append("Barrier3D domains")
        fig.legend(handles, labels, loc="outside lower center", ncol=len(handles), frameon=False)
        fig.suptitle(f"GIS {centre - HALF}–{centre + HALF}", fontsize=9.5)
        what = "" if version == "line" else f"_{sub}"
        stem = f"mean_shoreline_{window.label}_on_imagery{what}_GIS{centre:02d}_{key}"
        dot_text = (" Dots: every satellite position the mean was taken over, geolocated along "
                    "its transect and coloured by date on the scale below the panels; the "
                    f"{'triangles' if n > 1 else 'triangle'} on the scale "
                    f"{'mark each' if n > 1 else 'marks the'} photograph's flight date."
                    if dots else
                    " White boxes with grey edges: the Barrier3D model domains (500 m "
                    "alongshore), labelled with their GIS numbers; the panel's top and "
                    "bottom edges are the outer domain boundaries."
                    if doms else "")
        fs.caption(fig, (
            f"The CoastSat mean shoreline for {window.described} (ink, white halo) and "
            f"the ±1 standard deviation of the satellite positions behind each transect mean "
            f"(translucent white band with dashed edges, placed along each transect's direction), on "
            f"{photo_ref(imagery, window)}, GIS domains {centre - HALF}–{centre + HALF}, north up"
            f"{', all panels at one extent and scale' if n > 1 else ''}.{dot_text} "
            f"{photo_timing(imagery, window)} White is outside the photographs."))
        path = out_dir / sub / f"{stem}.png"
        fs.save(fig, path, vector=False, close=True)
        (PUBLISH / sub).mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, PUBLISH / sub / path.name)
        row["file" if version == "line" else f"file_{sub}"] = f"{sub}/{path.name}"
        if dots:
            row.update({f"positions_{y}": int((shown["year"] == y).sum())
                        for y in range(window.lo.year, window.hi.year + 1)},
                       positions_total=len(shown))
    return row


# The island overview (Hannah, 2026-09-23

ISLAND_YEAR = 1996                        # the default; --island-year overrides
SEGMENTS = [("GIS 1–30", 1, 30), ("GIS 31–60", 31, 60), ("GIS 61–90", 61, 90)]
ISLAND_RES_M = 5.0                        # ~ the printed pixel at this scale
ISLAND_LAND_M, ISLAND_SEA_M = 900.0, 600.0
ISLAND_PANEL_H_IN = 8.2                   # the segments are 15 km tall; this sets the scale


# Three north-up segments side by side at one scale, on one year's photos, with the site windows of ...
def island_figure(df, boxes, im, sites, window, out_dir):
    ext = []
    for _, lo, hi in SEGMENTS:
        seg = boxes[boxes["gis"].between(lo, hi)]
        y0, y1 = seg.total_bounds[1], seg.total_bounds[3]
        pts = df[df["y"].between(y0, y1)]
        ext.append((pts["x"].min() - ISLAND_LAND_M, y0, pts["x"].max() + ISLAND_SEA_M, y1))
    widths = [b[2] - b[0] for b in ext]
    height = max(b[3] - b[1] for b in ext)
    scale = min(ISLAND_PANEL_H_IN / height, (fs.FIG_W_DOUBLE - 1.3) / sum(widths))   # in per m
    imgs = []
    for b, (title, _, _) in zip(ext, SEGMENTS):
        print(f"    {title} ...")
        img = im.read(b[0], b[2], b[1], b[3], ISLAND_RES_M, CRS).copy()
        img[img.max(axis=2) == 0] = 255
        imgs.append(img)
    for version in ("line", "domains"):
        _island_version(df, boxes, im, sites, window, out_dir, ext, imgs, height * scale,
                        version)


# One island figure (domains or plain)
def _island_version(df, boxes, im, sites, window, out_dir, ext, imgs, panel_h, version):
    doms = version == "domains"
    widths = [b[2] - b[0] for b in ext]
    fig = plt.figure(figsize=(fs.FIG_W_DOUBLE, panel_h + 1.0), layout="constrained")
    axes = fig.subplots(1, 3, gridspec_kw={"width_ratios": widths})
    for i, (ax, b, img, (title, lo, hi)) in enumerate(zip(axes, ext, imgs, SEGMENTS)):
        x0, y0, x1, y1 = b
        ax.imshow(img, extent=(x0, x1, y0, y1), origin="upper", zorder=0,
                  interpolation="bilinear")
        draw_line(ax, df, b, lw=0.8, edges=False)
        if doms:
            draw_domains(ax, boxes, b, lw=0.5, labels=False)
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
    handles = [Line2D([], [], color=C_LINE, lw=0.8,
                      path_effects=[pe.withStroke(linewidth=2.5, foreground="white")]),
               Patch(facecolor="0.9", edgecolor="none"),
               Patch(facecolor="none", edgecolor=C_LINE, lw=0.7)]
    labels = [f"Mean shoreline, {window.span} (CoastSat)", "±1 standard deviation",
              "Site windows (zoom figures)"]
    if doms:
        handles.append(_domain_handle(0.5))
        labels.append("Barrier3D domains")
    fig.legend(handles, labels, loc="outside lower center", ncol=len(handles), frameon=False)
    short = PHOTO_SOURCES.get(im.year, {}).get("short", "USGS aerial photographs")
    fig.suptitle(f"{short}, {photo_label(im)}", fontsize=9.5)
    fs.caption(fig, (
        f"The CoastSat mean shoreline for {window.described} (ink, white halo) and the "
        f"±1 standard deviation of the satellite positions behind each transect mean (white "
        f"band; ~25 m wide, so near the limit of the print at this scale) on "
        f"{photo_ref([im], window)}, the island "
        f"in three north-up segments at one scale: (a) GIS 1–30, Cape Point to Avon; (b) GIS "
        f"31–60; (c) GIS 61–90, the Tri-Village to the north end. Ticks mark every fifth "
        f"Barrier3D domain. Ink boxes outline the six site windows drawn at full resolution "
        f"in the zoom figures, labelled with their domains."
        + (" White boxes with grey edges: every Barrier3D model domain (500 m alongshore), "
           "numbered by the ticks." if doms else "")
        + ("" if _inside(im, window) else
           f" The photograph is from outside the window ({photo_label(im)}).")
        + " White is outside the photographs."))
    sub = VERSION_DIRS[version]
    what = "_with_domains" if doms else ""
    path = out_dir / sub / f"mean_shoreline_{window.label}_on_imagery{what}_island_{im.year}.png"
    fs.save(fig, path, vector=False, close=True)
    (PUBLISH / sub).mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, PUBLISH / sub / path.name)
    print(f"  island figure -> {sub}/{path.name}")


RIBBON_RES_M = 8.0                        # ~ the printed pixel of a 45 km ribbon


# Panel (a) of the diagnostic, alone, over one year's photographs
def ribbon_figure(df, im, window, out_dir):
    import coastsat_mean_shoreline as cms
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
    ax.set_title(f"Mean shoreline, {window.span}, on the {photo_label(im)} photographs")
    fig.legend([Line2D([], [], color=C_LINE, lw=0.8,
                       path_effects=[pe.withStroke(linewidth=2.3, foreground="white")]),
                Line2D([], [], color="white", lw=0.8,
                       path_effects=[pe.withStroke(linewidth=1.8, foreground=fs.INK_MUTED)])],
               [f"Mean shoreline, {window.span} (CoastSat)", "Island outline"],
               loc="outside lower center", ncol=2, frameon=False)
    fs.caption(fig, (
        f"The CoastSat mean shoreline for {window.described} (ink, white halo) on "
        f"{photo_ref([im], window)}, drawn as panel (a) of "
        f"mean_shoreline_{window.label}.png: alongshore across the page, easting "
        f"increasing downward, equal aspect, so the ocean is at the bottom. The ±1 standard deviation band (~25 m) is below "
        f"the resolution of the print and is not drawn. The island outline "
        f"(map_elements/hatteras_outline; its survey date is not recorded) is drawn in white "
        f"over the photographs, and as grey land and pale-blue water where no frame covers."))
    path = out_dir / VERSION_DIRS["line"] / f"mean_shoreline_{window.label}_on_imagery_ribbon_{im.year}.png"
    fs.save(fig, path, vector=False, close=True)
    shutil.copy2(path, PUBLISH / VERSION_DIRS["line"] / path.name)
    print(f"  ribbon figure -> {path.name}")


# Run: the window's line, then the site, island and ribbon figures
def main(argv=None) -> int:
    import coastsat_mean_shoreline as cms
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    how = ap.add_mutually_exclusive_group()
    how.add_argument("--window", nargs=2, type=int, metavar=("START", "END"),
                     help="calendar years (the default, 1995 1997)")
    how.add_argument("--window-dates", nargs=2, metavar=("FIRST", "LAST"),
                     help="a window built by coastsat_mean_shoreline.py --window-dates")
    how.add_argument("--centred-on", choices=sorted(cms.SURVEY_ANCHORS),
                     help="a window built by coastsat_mean_shoreline.py --centred-on")
    ap.add_argument("--half-width", type=int, default=cms.DEFAULT_HALF_WIDTH_YEARS)
    ap.add_argument("--sites", nargs="*", type=int,
                    help="centre GIS domains (default: the six built-in sites)")
    ap.add_argument("--only", choices=("sites", "island", "ribbon"),
                    help="draw only the site figures, the three-segment island "
                         "overview, or the single-panel ribbon")
    ap.add_argument("--photo-years", nargs="+", type=int,
                    help="photograph years to draw on (default: the window's own years); "
                         "for a window with none, the nearest year, e.g. 2008 for 2009-2011")
    ap.add_argument("--island-year", type=int,
                    help=f"photographs for the island and ribbon figures (default: "
                         f"{ISLAND_YEAR} if drawn, else the first photo year)")
    a = ap.parse_args(argv)
    if a.centred_on:
        window = cms.Window.centred_on(a.centred_on, a.half_width)
    elif a.window_dates:
        window = cms.Window.from_dates(*a.window_dates)
    else:
        window = cms.Window.from_years(*(a.window or DEFAULT_WINDOW))
    sites = ([(c, k) for c, k in DEFAULT_SITES if c in a.sites] + [(c, "site") for c in a.sites
              if c not in dict(DEFAULT_SITES)]) if a.sites else DEFAULT_SITES

    fs.apply_style()
    R = _review_module()
    imagery = []
    years = a.photo_years or range(window.lo.year, window.hi.year + 1)
    for y in years:
        if R._year_files(y):
            im = R.Imagery(y, CRS)
            if y in PHOTO_SOURCES:
                im.date = PHOTO_SOURCES[y]["date"]
            imagery.append(im)
        elif a.photo_years:
            raise SystemExit(f"no photographs for {y} under {R.AERIAL_ROOT} (is D: on?)")
    if not imagery:
        raise SystemExit(f"no photographs for {window.label} under {R.AERIAL_ROOT} (is D: on?); "
                         f"a window with none takes --photo-years")

    df = transects(window)
    pos = positions(df, window)
    pos["t"] = [_decimal_year(d) for d in pos["date"]]
    # the date scale spans the window, widened to reach any photograph outside it
    flights = [_decimal_year(im.date) for im in imagery]
    lo_t = _decimal_year(window.lo)
    hi_t = _decimal_year(window.hi + pd.Timedelta(days=1))
    if window.calendar:     # whole years, as drawn since 09-23
        norm = matplotlib.colors.Normalize(min(lo_t, int(min(flights))),
                                           max(hi_t, int(max(flights)) + 1))
    else:                   # a date window: the scale is the window itself
        norm = matplotlib.colors.Normalize(min(lo_t, min(flights)), max(hi_t, max(flights)))
    print(f"  {len(pos)} satellite positions over {pos['transect_id'].nunique()} transects")
    boxes = gpd.read_file(DOMAIN_BOXES).to_crs(CRS)
    boxes["gis"] = np.arange(1, len(boxes) + 1)   # file order is south -> north, GIS 1-90

    out_dir = mean_shoreline_dir(*window.key) / "on_imagery"
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
        want = a.island_year or (ISLAND_YEAR if any(i.year == ISLAND_YEAR for i in imagery)
                                 else imagery[0].year)
        im = next((im for im in imagery if im.year == want), None)
        if im is None:
            print(f"  no {want} photographs drawn; island figures skipped")
        else:
            if a.only in (None, "island"):
                island_figure(df, boxes, im, pd.read_csv(table), window, out_dir)
            if a.only in (None, "ribbon"):
                ribbon_figure(df, im, window, out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
