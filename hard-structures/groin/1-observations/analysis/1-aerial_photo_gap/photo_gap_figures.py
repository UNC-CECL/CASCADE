#!/usr/bin/env python3
"""
The aerial-photo gap figures, drawn from the tables photo_gap.py writes.

    python photo_gap_figures.py   ->  figures/photo_gap_<n>_<what>.png (and the GIF)

The photos alone, in the order of the analysis: (1) the wet/dry lines, (2) how one
survey becomes a position per domain, (3) every domain's position in every survey,
(4) the two domains either side of the groins and the gap between them (and the same
animated), (5) the gap's rate in each stage, (6) where along the coast it built out and
retreated, and the gap across both candidate model boundaries (GIS 5/6 and 6/7). Captions go to
figures/supporting/CAPTIONS.md. CoastSat, and how the two records compare, are in
../2-coastsat_gap/ and ../3-photos_vs_coastsat/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.colors import ListedColormap, Normalize  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
from photo_gap import BOUNDARIES, ERAS, PHOTOS, STRIP_DOMAINS, STRIP_M, TABLES  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1997, INK, INK_MUTED, MAP_TEXT_DARK, _letter_inside, _title, apply_style, figsize, north_dart,
    open_frame, record_caption, save, scale_bar_km, spines_for_image)

# --- CONFIG ------------------------------------------------------------------
FIGS = HERE / "figures"
SIDE_COL = {"north": C_1997, "south": C["ADDED"]}
SIDE = {"north": "north side (updrift), GIS 6", "south": "south side (downdrift), GIS 5"}
EVENTS = ((1969, "built"), (1975, "repaired"), (1994, "Gordon"), (1995, "last repair"), (2003, "Isabel"),
          (2017, "Buxton fill"))
LOSSES = (2015, 2021)                # structural losses reported (../../data/structure_history.md)
WETDRY_FILE = "wet_dry_shorelines_groin.geojson"   # in the sanity check's INPUT_DIR
GROINS = HERE.parent.parent / "data" / "gis" / "groins_hatteras.geojson"
YEAR_CMAP = "viridis"
# -----------------------------------------------------------------------------


def tab(name):
    return pd.read_csv(TABLES / f"{name}.csv")


def out(name):
    return FIGS / f"photo_gap_{name}.png"


def chart(ax):
    ax.grid(axis="y")
    open_frame(ax)


# Dotted event lines, named along the top
def events(ax):
    for yr, lab in EVENTS:
        ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":", zorder=1)
        ax.text(yr + (-0.3 if yr == 1994 else 0.3), 0.985, lab, transform=ax.get_xaxis_transform(), fontsize=6.5,
                color=INK_MUTED, va="top", ha="right" if yr == 1994 else "left")


# Geometry for the method figures: the same files and functions HAT_geometric_distance_sanity_check.py uses
def geometry():
    import HAT_geometric_distance_sanity_check as h
    doms = h.load_domain_polygons(h.DOMAIN_POLYGONS_FILE)
    datum = h.load_single_geom(h.DATUM_LINE_FILE)
    lines = {k: v for k, v in h.load_shorelines(os.path.join(h.INPUT_DIR, WETDRY_FILE)).items()}
    groins = [np.asarray(f["geometry"]["coordinates"]) for f in json.load(open(GROINS))["features"]]
    return h, doms, datum, lines, groins


def year_colours(lines):
    keys = sorted(lines, key=lambda k: lines[k][1])
    yrs = np.array([lines[k][1] for k in keys])
    cmap = plt.get_cmap(YEAR_CMAP)
    norm = Normalize(yrs.min(), yrs.max())
    return keys, {k: cmap(0.08 + 0.82 * norm(lines[k][1])) for k in keys}, cmap, norm


def domain_boxes(ax, doms, label=True, ids=None, x_label=None):
    for d, poly in doms.items():
        if ids and d not in ids:
            continue
        x, y = poly.exterior.xy
        ax.fill(x, y, facecolor="0.95" if d % 2 else "0.90", edgecolor="0.72", lw=0.5, zorder=0)
        if label:
            cx = poly.bounds[0] + 60 if x_label is None else x_label
            ax.text(cx, poly.centroid.y, f"GIS {d}", fontsize=6.5, color=INK_MUTED, va="center", zorder=6,
                    clip_on=True)


def groin_lines(ax, groins, lw=2.2):
    for c in groins:
        ax.plot(c[:, 0], c[:, 1], color=C["GROIN"], lw=lw, solid_capstyle="butt", zorder=7)


def map_frame(ax, xlim, ylim):
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)


# 1. The wet/dry lines themselves
def fig_lines():
    h, doms, datum, lines, groins = geometry()
    keys, col, cmap, norm = year_colours(lines)
    b = np.array([p.bounds for p in doms.values()])
    xlim, ylim = (b[:, 0].min() + 1050, b[:, 2].max() - 50), (b[:, 1].min() - 150, b[:, 3].max() + 150)
    gx = np.concatenate([c[:, 0] for c in groins])
    gy = np.concatenate([c[:, 1] for c in groins])
    zx, zy = (gx.mean() - 450, gx.mean() + 450), (gy.mean() - 600, gy.mean() + 600)
    fig = plt.figure(figsize=figsize("double", height=6.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=(1, 1.05, 0.08))
    for i, (lim_x, lim_y, title) in enumerate(((xlim, ylim, "All 24 surveys, GIS 2-12"),
                                               (zx, zy, "Close-up at the groins"))):
        ax = fig.add_subplot(gs[0, i])
        domain_boxes(ax, doms, x_label=lim_x[0] + 40)
        for k in keys:
            g = lines[k][0]
            for part in getattr(g, "geoms", [g]):
                x, y = part.xy
                ax.plot(x, y, color=col[k], lw=0.9 if i == 0 else 1.2, zorder=3)
        groin_lines(ax, groins)
        map_frame(ax, lim_x, lim_y)
        scale_bar_km(ax, length_m=500 if i == 0 else 200, segments=1, unit="m", x=0.55, y=0.04,
                     text_kw=MAP_TEXT_DARK)
        if i == 0:
            ax.add_patch(Rectangle((zx[0], zy[0]), zx[1] - zx[0], zy[1] - zy[0], fill=False, edgecolor=INK,
                                   lw=0.9, ls="--", zorder=8))
            north_dart(ax, (lim_x[1] - 250, lim_y[0] + 700), arrow_m=300, text_kw=MAP_TEXT_DARK)
            ax.text(lim_x[1] - 30, lim_y[0] + 2200, "ocean →", ha="right", fontsize=7, color=INK_MUTED)
        _letter_inside(ax, i)
        ax.set_title(title, loc="center")
    cax = fig.add_subplot(gs[0, 2])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=ListedColormap(cmap(np.linspace(0.08, 0.90, 256))))
    fig.colorbar(sm, cax=cax, label="survey year")
    fig.legend(handles=[Line2D([], [], color=C["GROIN"], lw=2.2, label="Buxton groins"),
                        Patch(facecolor="0.92", edgecolor="0.72", label="model domains (GIS), 500 m alongshore")],
               loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("1_wetdry_lines")
    save(fig, png, close=True)
    yrs = sorted({int(lines[k][1]) for k in keys})
    record_caption(png, (
        "**Tests:** nothing; this is the raw record the rest of the folder measures. **How to read:** UTM 18N, "
        "true scale, north up, ocean to the right. Each line is one wet/dry line digitized from a dated aerial "
        f"photo, {len(keys)} surveys in {len(yrs)} years from {yrs[0]} to {yrs[-1]} (two in 1972), coloured by "
        "date (`wet_dry_shorelines_groin.geojson` in the 1967 rig's inputs, the file the position table is "
        "made from). Grey boxes: the model's GIS domains 2-12, each 500 m alongshore. Red: the four Buxton "
        "groin lines. (a) The whole stretch the table covers; dashed box: the extent of (b). (b) The groins "
        "close up. **Shows:** the coast retreats landward (left) through time almost everywhere, and the lines "
        "fan out most just south of the groins, where the downdrift side lost the most ground."))


# 2. How one survey becomes a position per domain
def fig_measure():
    h, doms, datum, lines, groins = geometry()
    key = "wetdry_1996"
    line = lines[key][0]
    dx_datum = datum.bounds[0]
    show = (4, 5, 6, 7)
    b = np.array([doms[d].bounds for d in show])
    xlim, ylim = (b[:, 0].min() + 1000, b[:, 2].max() + 80), (b[:, 1].min() - 30, b[:, 3].max() + 30)
    fig = plt.figure(figsize=figsize("double", height=5.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=(1.25, 1))

    ax = fig.add_subplot(gs[0, 0])
    domain_boxes(ax, doms, ids=show, x_label=xlim[0] + 30)
    x, y = line.xy
    ax.plot(x, y, color=INK_MUTED, lw=0.8, zorder=2)
    angles = {}
    for d in show:
        clip = line.intersection(doms[d])
        parts = getattr(clip, "geoms", [clip])
        pts = [p.interpolate(i / (h.N_SAMPLES - 1), normalized=True) for p in parts for i in range(h.N_SAMPLES)]
        px, py = np.array([p.x for p in pts]), np.array([p.y for p in pts])
        c = SIDE_COL["north"] if d == 6 else SIDE_COL["south"] if d == 5 else INK
        for p in parts:
            ax.plot(*p.xy, color=c, lw=1.8, zorder=4)
        ax.plot(px, py, "o", ms=2.2, color=c, mec="white", mew=0.3, zorder=5)
        for xx, yy in zip(px[::4], py[::4]):
            ax.annotate("", xy=(xlim[1], yy), xytext=(xx, yy),
                        arrowprops=dict(arrowstyle="->", color=C["LOCATOR"], lw=0.5, alpha=0.8), zorder=3)
        dist = dx_datum - px.mean()
        ax.text(xlim[1] - 20, doms[d].centroid.y + 95, f"mean {dist:,.0f} m", ha="right", fontsize=6.6,
                color=c, zorder=6, bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.8))
        e = np.polyfit(py, px, 1)[0]                       # metres east per metre north along the line
        angles[d] = np.degrees(np.arctan(e))
    groin_lines(ax, groins)
    map_frame(ax, xlim, ylim)
    ax.text(xlim[0] + 30, ylim[1] - 120, "teal arrows: east to the datum line,\n~7.4 km offshore, off the panel",
            fontsize=6.6, color=C["LOCATOR"], va="top",
            bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85), zorder=8)
    scale_bar_km(ax, length_m=200, segments=1, unit="m", x=0.06, y=0.07, text_kw=MAP_TEXT_DARK)
    _letter_inside(ax, 0)
    ax.set_title("One survey (Oct 1996), GIS 4-7", loc="center")

    # (b) the steps, in words, with this survey's numbers for GIS 5 and 6
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off")
    _letter_inside(ax, 1)
    ax.set_title("From a line to a number", loc="center")
    d67 = h.per_domain_mean_distance(lines["wetdry_1967"][0], doms, datum)
    d96 = h.per_domain_mean_distance(line, doms, datum)
    ch = {d: d67[d] - d96[d] for d in (5, 6)}
    steps = [
        ("1  Clip", "cut the survey's line at each domain's box"),
        ("2  Sample", f"{h.N_SAMPLES} points evenly along each clipped piece"),
        ("3  Measure", "each point's distance to the offshore datum line,\na straight north-south line ~7.4 km "
                       "east,\nso the distance runs due east"),
        ("4  Average", f"the mean of the {h.N_SAMPLES} distances is the domain's\nposition in that survey"),
        ("5  Difference", "subtract it from the same domain's 1967 value:\nchange since 1967 "
                          "(here + seaward)"),
    ]
    yy = 0.90
    for head, body in steps:
        ax.text(0.0, yy, head, fontsize=7.3, fontweight="bold", va="top", transform=ax.transAxes)
        ax.text(0.30, yy, body, fontsize=6.9, va="top", transform=ax.transAxes, linespacing=1.25)
        yy -= 0.05 + 0.055 * (body.count("\n") + 1)
    ax.text(0.0, yy - 0.02, "This survey", fontsize=7.3, fontweight="bold", va="top", transform=ax.transAxes)
    for k, d in enumerate((6, 5)):
        ax.text(0.30, yy - 0.02 - 0.06 * k,
                f"GIS {d}:  {d96[d]:,.0f} m  (1967: {d67[d]:,.0f} m)  →  {ch[d]:+.0f} m",
                fontsize=6.9, va="top", transform=ax.transAxes,
                color=SIDE_COL["north"] if d == 6 else SIDE_COL["south"])
    png = out("2_measuring_position")
    save(fig, png, close=True)
    ang = (angles[5] + angles[6]) / 2
    record_caption(png, (
        "**Tests:** nothing; this shows how `HAT_geometric_distance_sanity_check.py` turns a digitized line "
        "into the per-domain numbers in `Change_from_wetdry_1967_D2_D12.csv`. **How to read:** (a) The October "
        "1996 wet/dry line over GIS 4-7 at true scale, north up; each domain's clipped piece is drawn heavy "
        f"(blue GIS 6, gold GIS 5) with its {h.N_SAMPLES} sample points; teal arrows (every fourth point) run "
        "east toward the offshore datum line, which lies about 7.4 km offshore, off the panel; the label is "
        "the mean distance for that domain. (b) The five steps, and this survey's numbers for the two domains "
        "either side of the groins. **Shows:** a domain's position is one number per survey, the mean "
        "east-west distance from its piece of line to a fixed line offshore, and the table holds that number's "
        "change since 1967 (stored landward positive; the figures here flip it to seaward positive). The datum "
        "runs due north-south, so the distance is measured due east rather than straight offshore: the "
        f"shoreline at GIS 5-6 runs about {abs(ang):.0f}° from north, so an east-west change is about "
        f"{1 / np.cos(np.radians(ang)):.2f} times the shore-normal change. That scales both sides about alike and "
        "does not move the gap's timing."))


# 3. Every domain's position, every survey
def fig_domains():
    t = pd.read_csv(PHOTOS).set_index("Domain_ID")
    cols = [c for c in t.columns if c.startswith("change_from_wetdry_1967_wetdry_")]
    lab = [c.replace("change_from_wetdry_1967_wetdry_", "")[:-2].replace("_", " ") for c in cols]
    m = -t[cols].to_numpy()                                # seaward positive
    order = np.argsort([float(s.split()[0]) + (0.5 if s.endswith("august") else 0) for s in lab])
    m, lab = m[:, order], [lab[i] for i in order]
    doms = t.index.to_numpy()
    fig, ax = plt.subplots(figsize=figsize("double", height=4.4), constrained_layout=True)
    lim = np.nanmax(np.abs(m))
    im = ax.imshow(m[::-1], cmap="RdBu", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m[::-1][i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{round(v) + 0:.0f}", ha="center", va="center", fontsize=4.8,
                        color="white" if abs(v) > 0.6 * lim else INK)
    ax.set_xticks(range(len(lab)), [s.replace(" july", "a").replace(" august", "b") for s in lab], rotation=90,
                  fontsize=6.5)
    ax.set_yticks(range(len(doms)), [f"GIS {d}" for d in doms[::-1]], fontsize=6.5)
    i5 = list(doms[::-1]).index(5)
    ax.axhline(i5 - 0.5, color=C["GROIN"], lw=1.6)
    ax.text(len(lab) - 0.4, i5 - 0.5, " groins", color=C["GROIN"], fontsize=6.8, va="center", ha="left",
            clip_on=False)
    for s in ax.spines.values():
        s.set_visible(False)
    fig.colorbar(im, ax=ax, label="change since 1967 (m)\n+ seaward, − landward", shrink=0.9)
    ax.set_xlabel("Survey (1972a July, 1972b August)")
    ax.set_title("Where did each domain's shoreline go, survey by survey?", loc="center")
    png = out("3_positions_by_domain")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** nothing; this is the position table every hindcast fit reads, drawn whole. **How to read:** "
        "each cell is one domain (rows, north at the top) in one survey (columns, in date order), the change in "
        "its mean distance to the datum since the 1967 survey (figure 2), seaward positive (blue) or landward "
        "(red), in metres; blank: that survey's line does not reach the domain. The red line is the GIS 5/6 "
        "boundary, where the southernmost groin sits (the field lies just inside GIS 6), so GIS 5 is the "
        "downdrift side and GIS 6 the updrift side. Source: `Change_from_wetdry_1967_D2_D12.csv`, sign flipped. **Shows:** the whole stretch "
        "retreats, but not evenly: GIS 4-5, just south of the groins, lose far more ground than GIS 6-8 just "
        "north of them, and the step in colour across the red line is the groin's effect. The next figure "
        "takes those two rows, GIS 5 and GIS 6, and their difference is the gap."))


# 4. Each side and the gap between them
def fig_sides():
    s = tab("photo_gap_annual").set_index("year")
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True, constrained_layout=True,
                             gridspec_kw=dict(height_ratios=(1.4, 1)))
    ax = axes[0]
    ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
    ax.fill_between(s.index, s["south_m"], s["north_m"], color=C["ACCENT_FILL"], alpha=0.5, lw=0)
    for side in ("north", "south"):
        ax.plot(s.index, s[f"{side}_m"], color=SIDE_COL[side], lw=1.3, marker="s", ms=3.2, label=SIDE[side])
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1965, 2026)
    ax.set_ylabel("Shoreline change since\n1967 (m), + seaward")
    events(ax)
    chart(ax)
    _title(ax, 0, "Where was each side of the groins?")

    ax = axes[1]
    ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
    ax.fill_between(s.index, 0, s["gap_m"], color=C["ACCENT_FILL"], alpha=0.5, lw=0)
    ax.plot(s.index, s["gap_m"], color=C["ACCENT"], lw=1.3, marker="s", ms=3.2)
    peak = int(s["gap_m"].idxmax())
    last = int(s.index.max())
    for yr, (dx, dy, ha) in ((peak, (-5.0, 8, "right")), (1996, (2.5, -28, "left")), (last, (1.0, 18, "left"))):
        g = s.loc[yr, "gap_m"]
        ax.annotate(f"{g:.0f} m\n{yr}", xy=(yr, g), xytext=(yr + dx, g + dy), ha=ha, va="center", fontsize=6.8,
                    color=C["ACCENT"], arrowprops=dict(arrowstyle="-", color=C["ACCENT"], lw=0.5))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_ylim(-5, 190)
    ax.set_ylabel("Gap (m)\nnorth minus south")
    ax.set_xlabel("Year")
    events(ax)
    chart(ax)
    _title(ax, 1, "How wide was the gap between them?")

    handles = [Line2D([], [], color=SIDE_COL[k], lw=1.3, marker="s", ms=3.2, label=SIDE[k]) for k in SIDE] + [
        Patch(facecolor=C["ACCENT_FILL"], alpha=0.5, label="the gap: the groin's effect"),
        Patch(facecolor=C["ADDED_FILL"], alpha=0.5, label="structural losses reported, 2015-2021")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("4_two_sides")
    save(fig, png, close=True)
    g96, gp, gl = s.loc[1996, "gap_m"], s.loc[peak, "gap_m"], s.loc[last, "gap_m"]
    record_caption(png, (
        "**Tests:** nothing; this is the photo record of the groin's effect, read on its own. **How to read:** "
        f"{len(s)} survey years of dated aerial-photo wet/dry lines, 1967-{last} (two surveys in 1972, "
        "averaged). (a) Each side's position as change since 1967, seaward positive: the mean distance from the "
        "wet/dry line to a fixed offshore datum line over the GIS domain, from "
        "`Change_from_wetdry_1967_D2_D12.csv` with its sign flipped. (b) The gap, north minus south: rising means "
        "the north side is holding while the south side retreats. Dotted lines: the structure's recorded events "
        "(`../../data/structure_history.md`); orange: the structural losses reported for 2015-2021. "
        f"**Shows:** both sides retreated; the gap is the north side retreating less. It grew fast after the "
        f"groins were built, peaked at {gp:.0f} m in {peak} (partly the south side's loss to Gordon), stood at "
        f"{g96:.0f} m in 1996 and fell to {gl:.0f} m by {last}, most of that fall around the 2017 fill and the "
        "losses. A photo is one day, so single surveys (1995, 2004) can sit well off the trend."))


# What is happening in each stretch of the record, for the animation
STORY = ((1969, "1967: the first survey. Both sides start at zero; there is no gap yet."),
         (1978, "Groins built 1969. Sand drifting south is caught on the north side, so the south side\n"
                "is starved and retreats faster. The gap grows quickly."),
         (1995.5, "The gap keeps growing, more slowly. Hurricane Gordon (1994) cuts the south side;\n"
                  "the last repair is in 1995, and the gap is at its widest."),
         (2015, "No more repairs. The gap stops growing and slowly shrinks: the north side keeps\n"
                "retreating while the south side holds."),
         (2024, "The groins fall apart (2015-2021) and the 2017 fill is placed across them. Sand gets past:\n"
                "the south side gains, the north side keeps retreating, and the gap shrinks."))


# 4, animated: the two shorelines and the gap through the photo record
def gif_sides():
    from matplotlib.animation import FuncAnimation, PillowWriter
    s = tab("photo_gap_annual").set_index("year")
    years = np.arange(1967, s.index.max() + 0.001, 0.25)
    north = np.interp(years, s.index, s["north_m"])
    south = np.interp(years, s.index, s["south_m"])
    frames = list(range(len(years))) + [len(years) - 1] * 25          # hold the last frame
    fig = plt.figure(figsize=figsize("double", height=4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=(0.55, 1), left=0.07, right=0.98, top=0.80, bottom=0.12, wspace=0.28)
    ax, at = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    story = fig.text(0.5, 0.975, "", ha="center", va="top", fontsize=8.2)
    stamp = fig.text(0.07, 0.86, "", fontsize=11, fontweight="bold")

    at.set_xlim(1965, 2025)
    at.set_ylim(-215, 20)
    at.axhline(0, color=INK_MUTED, lw=0.6)
    at.set_ylabel("Shoreline change since 1967 (m)\n+ seaward")
    at.set_xlabel("Year")
    chart(at)
    for yr, lab in ((1969, "built"), (1995, "last repair"), (2017, "fill")):
        at.axvline(yr, color=INK_MUTED, lw=0.6, ls=":")
        at.text(yr + 0.4, 15, lab, fontsize=6.5, color=INK_MUTED, va="top")
    at.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0)
    at.plot(s.index, s["north_m"], "s", color=SIDE_COL["north"], ms=3, alpha=0.35)
    at.plot(s.index, s["south_m"], "s", color=SIDE_COL["south"], ms=3, alpha=0.35)
    ln_n, = at.plot([], [], color=SIDE_COL["north"], lw=1.6, label="north side (updrift)")
    ln_s, = at.plot([], [], color=SIDE_COL["south"], lw=1.6, label="south side (downdrift)")
    at.legend(frameon=False, fontsize=7, loc="lower left")
    at.set_title("Each side's position; the shaded band is the gap", fontsize=8.5)
    dyn = []

    def draw(k):
        for a in dyn:
            a.remove()
        dyn.clear()
        yr, n, so = years[k], north[k], south[k]
        # plan view, true scale
        ax.clear()
        y = np.linspace(0, 1000, 401)
        x = np.where(y >= 500, n, so)
        ax.fill_betweenx(y, -320, x, color=C["ADDED_FILL"], lw=0)
        ax.fill_betweenx(y, x, 140, color=C["WATER"], alpha=0.35, lw=0)
        ax.plot(x, y, color=INK, lw=1.6)
        ax.axvline(0, color=INK_MUTED, lw=1.0, ls=(0, (4, 2)))
        ax.text(3, 990, "1967", fontsize=6.5, color=INK_MUTED, va="top")
        ax.axhline(500, color="white", lw=1.0)
        if yr >= 1969:
            ax.plot([-90, 70], [500, 500], color=C["GROIN"], lw=3.0, solid_capstyle="butt")
            ax.text(72, 515, "groins", fontsize=6.5, color=C["GROIN"])
        ax.annotate("", xy=(105, 650), xytext=(105, 900), arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.1))
        ax.text(97, 775, "drift", rotation=90, ha="right", va="center", fontsize=6.5, color=INK_MUTED)
        for yy, xx, side in ((700, n, "north"), (200, so, "south")):
            if xx < -2:
                ax.annotate("", xy=(xx, yy), xytext=(0, yy), arrowprops=dict(arrowstyle="->", color=SIDE_COL[side],
                                                                              lw=1.2))
            ax.text(min(xx, 0) - 6, yy + 18, f"{xx:.0f} m", ha="right", fontsize=7.5, color=SIDE_COL[side],
                    fontweight="bold")
        gap = n - so
        if abs(gap) > 2:
            ax.plot([n, n], [500, 330], color=SIDE_COL["north"], lw=1.0, ls=(0, (2, 1.5)))
            ax.annotate("", xy=(so, 360), xytext=(n, 360), arrowprops=dict(arrowstyle="<->", color=C["ACCENT"], lw=1.4))
        ax.text(-310, 390, f"gap\n{gap:.0f} m", fontsize=9, color=C["ACCENT"], fontweight="bold", va="center")
        ax.text(-310, 760, "GIS 6\nnorth", fontsize=7, va="center")
        ax.text(-310, 110, "GIS 5\nsouth", fontsize=7, va="center")
        ax.set_xlim(-320, 140)
        ax.set_ylim(0, 1000)
        ax.set_aspect("equal")
        ax.set_yticks([0, 500, 1000], ["0", "500", "1000 m"])
        ax.set_xlabel("Cross-shore (m), ocean →")
        open_frame(ax)
        ax.set_title("Map view (true scale)", fontsize=8.5)
        # time series up to this year
        m = years <= yr
        ln_n.set_data(years[m], north[m])
        ln_s.set_data(years[m], south[m])
        dyn.append(at.fill_between(years[m], south[m], north[m], color=C["ACCENT_FILL"], alpha=0.5, lw=0))
        dyn.append(at.annotate("", xy=(yr, so), xytext=(yr, n),
                               arrowprops=dict(arrowstyle="<->", color=C["ACCENT"], lw=1.2)))
        stamp.set_text(f"{int(yr)}")
        story.set_text(next(t for end, t in STORY if yr < end or end == STORY[-1][0]))
        return []

    anim = FuncAnimation(fig, draw, frames=frames, blit=False)
    gif = FIGS / "photo_gap_4_two_sides.gif"
    anim.save(gif, writer=PillowWriter(fps=12), dpi=110)
    plt.close(fig)
    print("wrote", gif)


# 5. The gap's rate in each stage
def fig_eras():
    e = tab("photo_era_rates")
    s = tab("photo_gap_annual").set_index("year")["gap_m"]
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.2), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=(1.5, 1)))
    ax = axes[0]
    ax.plot(s.index, s.values, ls="none", marker="s", ms=3.6, color=C["ACCENT"], zorder=4)
    for (lo, hi), r in zip(ERAS, e.itertuples()):
        w = s.loc[lo:hi]
        t = np.array([w.index.min(), w.index.max()], float)
        ax.plot(t, w.mean() + r.rate_m_yr * (t - w.index.to_numpy(float).mean()), color=INK, lw=1.2)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1965, 2026)
    ax.set_ylabel("Gap (m)\nnorth minus south")
    ax.set_xlabel("Year")
    events(ax)
    chart(ax)
    _title(ax, 0, "The gap and a line through each stage")
    ax = axes[1]
    x = np.arange(len(e))
    ax.bar(x, e["rate_m_yr"], yerr=e["se_m_yr"], color=C["ACCENT"], width=0.6, capsize=2.5, error_kw=dict(lw=0.8))
    for xi, r in zip(x, e.itertuples()):
        ax.text(xi, r.rate_m_yr + np.sign(r.rate_m_yr) * (r.se_m_yr + 0.4), f"{r.rate_m_yr:+.1f}", ha="center",
                va="bottom" if r.rate_m_yr >= 0 else "top", fontsize=7)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    names = ("building", "slow growth\nto the peak", "after the\nlast repair")
    ax.set_xticks(x, [f"{n}\n{era}\n({k} surveys)" for n, era, k in zip(names, e["era"], e["n_years"])], fontsize=6.6)
    ax.set_ylim(-4, 15)
    ax.set_ylabel("Rate of gap change (m/yr)\n± 1 SE; + = widening")
    chart(ax)
    _letter_inside(ax, 1)
    ax.set_title("Widening or narrowing?", loc="center")
    png = out("5_stage_rates")
    save(fig, png, close=True)
    b, g, a = e.itertuples()
    record_caption(png, (
        "**Tests:** how fast the gap widened or narrowed in each stage of the groin's life, from the photos alone. "
        "**How to read:** (a) The photo gap (north minus south, seaward positive) with a least-squares line through "
        f"the surveys in each stage: building ({b.era}), slow growth to the {g.era[-4:]} peak ({g.era}), and after "
        f"the last repair to the year before the fill ({a.era}). (b) Each line's slope ± 1 standard error, with the "
        "number of survey years behind it. Dotted lines: the structure's recorded events. "
        f"**Shows:** the gap built at {b.rate_m_yr:+.1f} m/yr in the first decade, slowed to {g.rate_m_yr:+.1f} "
        f"± {g.se_m_yr:.1f} m/yr up to the peak, and narrowed at {a.rate_m_yr:+.1f} ± {a.se_m_yr:.1f} m/yr after "
        "the last repair. The photos are sparse after 1997 (one survey between 1997 and 2008, in 2004), so the "
        "timing of the turn is better read from CoastSat (`../2-coastsat_gap/`)."))


# Where the candidate boundaries and the groins sit, in metres of northing from the southern groin
def placement_frame():
    h, doms, datum, lines, groins = geometry()
    south = min(c[:, 1].min() for c in groins)
    bounds = {f"{s}/{n}": doms[n].bounds[1] - south for s, n in BOUNDARIES}
    groin_y = sorted({round(c[:, 1].min() - south) for c in groins})
    mids = {d: doms[d].centroid.y - south for d in range(STRIP_DOMAINS[0], STRIP_DOMAINS[1] + 1)}
    return south, bounds, groin_y, mids


def placement_axes(ax, bounds, groin_y, mids, label=True):
    for name, x in bounds.items():
        ax.axvline(x, color=INK, lw=0.9, ls="--", zorder=1)
        if label:
            ax.text(x, 0.985, f" GIS {name}", transform=ax.get_xaxis_transform(), fontsize=6.8, va="top")
    for x in groin_y:
        ax.plot(x, 0, transform=ax.get_xaxis_transform(), marker="^", ms=6, color=C["GROIN"], clip_on=False,
                zorder=6)
    for d, x in mids.items():
        ax.text(x, 0.03, f"GIS {d}", transform=ax.get_xaxis_transform(), fontsize=6.5, color=INK_MUTED,
                ha="center")


# 6. Where along the coast the shoreline built out and where it retreated
def fig_placement():
    r = tab("photo_alongshore_change")
    g = tab("photo_boundary_gaps").set_index("year")
    south, bounds, groin_y, mids = placement_frame()
    show = (1978, 1989, 1995, 2004, 2016, 2023)
    cmap = plt.get_cmap(YEAR_CMAP)
    norm = Normalize(1967, 2023)
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.2), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=(1.4, 1)))
    ax = axes[0]
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    for yr in show:
        s = r[r["year"] == yr].sort_values("y_mid")
        ax.plot(s["y_mid"] - south, s["change_m"], color=cmap(0.08 + 0.82 * norm(yr)), lw=1.4, marker="o", ms=2.6,
                label=str(yr), zorder=3)
    placement_axes(ax, bounds, groin_y, mids)
    ax.set_xlabel("Northing from the southern groin (m), south → north")
    ax.set_ylabel("Change since 1967 (m)\n+ seaward (built out)")
    ax.legend(title="survey", frameon=False, fontsize=6.6, title_fontsize=6.8, ncol=2, loc="upper left")
    chart(ax)
    _title(ax, 0, "Where did the shoreline build out, and where did it retreat? (100 m strips)")

    ax = axes[1]
    for (s, n), col in zip(BOUNDARIES, (C["ACCENT"], INK)):
        k = f"gap_{s}_{n}_m"
        ax.plot(g.index, g[k], color=col, lw=1.2, marker="s", ms=3.4, label=f"across GIS {s}/{n}: GIS {n} minus GIS {s}")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
    ax.set_xlim(1965, 2026)
    ax.set_ylabel("Gap across the\nboundary (m)")
    ax.set_xlabel("Year")
    events(ax)
    chart(ax)
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    _title(ax, 1, "How big is the gap across each candidate boundary?")
    png = out("6_placement")
    save(fig, png, close=True)
    p89 = r[r["year"] == 1989].set_index("y_mid")["change_m"]
    peak_x = p89.idxmax() - south
    zero = p89[(p89.index - south > -100) & (p89.index - south < 450)]
    cross = zero.index[np.argmin(np.abs(zero.values))] - south
    pre = g.loc[1978:1997]
    record_caption(png, (
        "**Tests:** whether the observed pattern of building out and retreat puts the groin's effect at the "
        "GIS 5/6 boundary, where the model groin sits, or at 6/7. **How to read:** (a) Each line is one survey's "
        "change since 1967 along the coast, measured as the position table is (mean distance to the offshore "
        f"datum line) but in {STRIP_M:.0f} m strips of northing instead of 500 m domains, seaward positive; x is "
        "northing from the southern groin. Dashed lines: the two candidate boundaries; red triangles: the groins "
        "(the two northern ones nearly coincide). (b) The gap across each boundary from the domain table: the "
        "north domain's change minus the south domain's. Dotted lines: the structure's recorded events; orange: "
        "the structural losses reported for 2015-2021. **Shows:** the shoreline does not step at one place. "
        "Through the groins' working years it rises steadily across the field, from the deepest retreat in "
        f"GIS 4-5 through zero near {cross:.0f} m (inside GIS 6, between the groins) to the most building-out at "
        f"{peak_x:.0f} m (1989), at the 6/7 boundary just north of the northern groins, then falls off through "
        "GIS 7 to near zero in GIS 8. So the starved side begins at the southern groin and the fillet sits "
        "against the northern groins: one domain-boundary groin can catch only one of the two. Across 5/6 the "
        f"gap is {pre['gap_5_6_m'].min():.0f}-{pre['gap_5_6_m'].max():.0f} m over 1978-1997; across 6/7 it is "
        f"{pre['gap_6_7_m'].min():.0f}-{pre['gap_6_7_m'].max():.0f} m, smaller because GIS 6's mean averages the "
        "ramp inside the field."))


def main():
    apply_style()
    for f in (fig_lines, fig_measure, fig_domains, fig_sides, gif_sides, fig_eras, fig_placement):
        f()
        print(f"drew {f.__name__}")


if __name__ == "__main__":
    main()
