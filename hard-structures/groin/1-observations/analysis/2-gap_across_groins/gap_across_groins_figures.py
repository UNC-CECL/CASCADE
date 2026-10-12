#!/usr/bin/env python3
"""
The gap-across-groins figures, drawn from the tables coastsat_gap_across_groins.py writes.

    python gap_across_groins_figures.py   ->  figures/gap_across_groins_<n>_<what>.png

One figure per question, numbered in reading order: where the gap is measured, which
side moved, when the gap stopped widening, how certain that year is, the rate in each
era, whether the break survives other fits, and where the photos and CoastSat disagree.
Each panel title is the question it answers; captions (what it tests, how to read it,
what it shows) go to figures/supporting/CAPTIONS.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import geopandas as gpd  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
from coastsat_gap_across_groins import (  # noqa: E402
    BREAK_RANGE, ERAS, FIT, GROINS, N_BOOT, N_BOOT_CHECK, NEAR_M, PHOTOS, STEP_YEAR, TABLES, DOWNDRIFT, UPDRIFT, era_rate)
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, INK, INK_MUTED, MAP_TEXT_DARK, _letter_inside, _title, apply_style, figsize,
    map_label, north_dart, open_frame, place_label, record_caption, save, scale_bar_km,
    spines_for_image, water_label)
from site_layer.hat_map_layers import ISLAND_OUTLINE  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
FIGS = HERE / "figures"
DOMAINS = REPO / "data" / "hatteras_init" / "5-scr" / "2-transect-frame" / "transect_domains" / "HAT_domains.json"
MEAN_LINE = (REPO / "data" / "hatteras_init" / "5-scr" / "1-observations" / "mean_shoreline"
             / "1995-10-12_1997-10-12" / "transect_means_1995-10-12_1997-10-12.csv")
UTM = "EPSG:32618"
COL = {"domain": C_1997, "near_field": C["ACCENT"], "photos": C_1984}
SIDE_COL = {"updrift": C_1997, "downdrift": C["ADDED"]}
NAME = {"domain": "whole domains (GIS 6 vs GIS 5)", "near_field": f"within {NEAR_M:.0f} m of the groins",
        "photos": "aerial photos (wet/dry line)"}
SIDE = {"updrift": "north side (updrift)", "downdrift": "south side (downdrift)"}
CHECK = {"as fitted": "baseline, 1984-2016", "from 1988": "start in 1988\n(sparse early images out)",
         "1995 left out": "1995 peak left out", "to 2013": "end in 2013"}
EVENTS = {1994: "Gordon", 1995: "last repair", 2003: "Isabel", 2017: "Buxton fill"}
GAP_LABEL = "Gap across the groins (m)\n+ = north side seaward"
MAP_PAD_M = 250.0
# -----------------------------------------------------------------------------

WHAT = ("The gap is the north (updrift) shoreline position minus the south (downdrift) one, seaward "
        "positive, so a rising gap means the north side is holding while the south side retreats, "
        "which is what a working groin does. CoastSat values are calendar-year means per transect "
        f"(at least 3 images), each taken about its own {FIT[0]}-{FIT[1]} mean and averaged per side "
        "in years where at least half the side's transects have data.")


def tab(name):
    return pd.read_csv(TABLES / f"{name}.csv")


# The question centred above; the letter beside it, or inside the corner on a narrow panel
def question(ax, i, text, inside=False):
    if inside:
        _letter_inside(ax, i)
        ax.set_title(text, loc="center")
    else:
        _title(ax, i, text)


# Dotted event lines, named once along the top
def events(ax, label=True, years=None):
    lo, hi = ax.get_xlim()
    for yr, lab in EVENTS.items():
        if (years and yr not in years) or not lo <= yr <= hi:
            continue
        ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":", zorder=1)
        if label:
            ax.text(yr + (-0.25 if yr == 1994 else 0.25), 0.985, lab, transform=ax.get_xaxis_transform(),
                    fontsize=6.5, color=INK_MUTED, va="top", ha="right" if yr == 1994 else "left")


# The years after the 2017 fill, left out of every fit
def fill_shade(ax):
    ax.axvspan(FIT[1] + 0.5, ax.get_xlim()[1], color=C["BASE_FILL"], alpha=0.45, lw=0, zorder=0)


FILL_HANDLE = Patch(facecolor=C["BASE_FILL"], alpha=0.45, label="after the 2017 fill (not fitted)")


def chart(ax):
    ax.grid(axis="y")
    open_frame(ax)


def out(name):
    return FIGS / f"gap_across_groins_{name}.png"


# Each side's wet/dry position per survey year, change since 1967, seaward positive
def photo_sides():
    t = pd.read_csv(PHOTOS).set_index("Domain_ID")
    out = {}
    for c in t.columns:
        m = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", c)
        if m:
            out.setdefault(int(m.group(1)), []).append(-t.loc[[UPDRIFT, DOWNDRIFT], c].to_numpy())
    return pd.DataFrame({y: np.mean(v, axis=0) for y, v in out.items()}, index=["north", "south"]).T.sort_index()


# 0. What the gap is, and how each record measures it
def fig_what_gap():
    s = photo_sides()
    n96, s96 = s.loc[1996, "north"], s.loc[1996, "south"]
    fig = plt.figure(figsize=figsize("double", height=5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=(0.62, 1), height_ratios=(1, 0.95))

    # (a) plan view at true scale: GIS 5 and 6, the 1967 and 1996 shorelines, both measurements
    ax = fig.add_subplot(gs[:, 0])
    y = np.linspace(0, 1000, 401)
    x96 = np.where(y >= 500, n96, s96)
    ax.fill_betweenx(y, -320, x96, color=C["ADDED_FILL"], lw=0, zorder=0)
    ax.fill_betweenx(y, x96, 230, color=C["WATER"], alpha=0.35, lw=0, zorder=0)
    ax.plot(x96, y, color=INK, lw=1.4, zorder=3)
    ax.axvline(0, color=INK_MUTED, lw=1.0, ls=(0, (4, 2)), zorder=3)
    ax.axhline(500, color="white", lw=1.0, zorder=1)
    ax.plot([-90, 70], [500, 500], color=C["GROIN"], lw=3.0, solid_capstyle="butt", zorder=4)
    for k in range(13):                                     # CoastSat transects: 6 in GIS 6, 13 in GIS 5
        yy = 20 + k * 36
        ax.plot([-300, 60], [yy, yy], color=SIDE_COL["downdrift"], lw=0.5, alpha=0.8, zorder=2)
    for k in range(6):
        yy = 540 + k * 80
        ax.plot([-300, 60], [yy, yy], color=SIDE_COL["updrift"], lw=0.5, alpha=0.8, zorder=2)
    ax.axvline(160, color=C["LOCATOR"], lw=1.0, ls=":", zorder=3)
    ax.text(166, 980, "offshore\ndatum line", fontsize=6.3, color=C["LOCATOR"], va="top")
    for yy, xx in ((600, n96), (250, s96)):                 # the photo distance to the datum
        ax.annotate("", xy=(xx, yy), xytext=(160, yy),
                    arrowprops=dict(arrowstyle="<->", color=C["LOCATOR"], lw=0.8), zorder=5)
    for yy, xx, side in ((690, n96, "updrift"), (180, s96, "downdrift")):   # retreat since 1967
        ax.annotate("", xy=(xx, yy), xytext=(0, yy), arrowprops=dict(arrowstyle="->", color=INK, lw=1.0), zorder=5)
        ax.text(xx - 6 if xx > -100 else xx / 2, yy + 14, f"{xx:.0f} m", ha="right" if xx > -100 else "center",
                fontsize=7, zorder=6)
    ax.plot([n96, n96], [500, 360], color=SIDE_COL["updrift"], lw=1.0, ls=(0, (2, 1.5)), zorder=4)
    ax.annotate("", xy=(s96, 380), xytext=(n96, 380),
                arrowprops=dict(arrowstyle="<->", color=C["ACCENT"], lw=1.3), zorder=6)
    ax.text((s96 + n96) / 2, 392, f"gap {n96 - s96:.0f} m", ha="center", fontsize=7.5, color=C["ACCENT"],
            fontweight="bold", zorder=6)
    ax.annotate("", xy=(110, 760), xytext=(110, 960), arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.1))
    ax.text(100, 860, "net drift", rotation=90, ha="right", va="center", fontsize=6.5, color=INK_MUTED)
    ax.text(-310, 750, "GIS 6\nnorth\n(updrift)", fontsize=7, va="center")
    ax.text(-310, 110, "GIS 5\nsouth\n(downdrift)", fontsize=7, va="center")
    ax.text(4, 990, "1967", fontsize=6.5, color=INK_MUTED, va="top")
    ax.text(n96 - 4, 990, "1996", fontsize=6.5, va="top", ha="right")
    ax.text(75, 508, "groins", fontsize=6.5, color=C["GROIN"])
    ax.set_xlim(-320, 230)
    ax.set_ylim(0, 1000)
    ax.set_aspect("equal")
    ax.set_xlabel("Cross-shore (m), ocean to the right")
    ax.set_yticks([0, 500, 1000], ["0", "500", "1000 m"])
    open_frame(ax)
    question(ax, 0, "What is the gap?", inside=True)

    # (b) each side's position over time from the photos; the gap is the space between them
    ax = fig.add_subplot(gs[0, 1])
    ax.fill_between(s.index, s["south"], s["north"], color=C["ACCENT_FILL"], alpha=0.5, lw=0)
    for side, key in (("updrift", "north"), ("downdrift", "south")):
        ax.plot(s.index, s[key], color=SIDE_COL[side], lw=1.3, marker="s", ms=3.2, label=SIDE[side])
    ax.annotate("", xy=(1996, s96), xytext=(1996, n96),
                arrowprops=dict(arrowstyle="<->", color=C["ACCENT"], lw=1.2))
    ax.text(1996.8, (n96 + s96) / 2, f"gap in 1996\n= {n96:.0f} − ({s96:.0f})\n= {n96 - s96:.0f} m",
            fontsize=6.8, color=C["ACCENT"], va="center")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1965, 2025)
    ax.set_ylabel("Shoreline change since\n1967 (m), + seaward")
    ax.legend(frameon=False, fontsize=6.8, loc="lower left")
    chart(ax)
    question(ax, 1, "Both sides retreat; the gap is the space between")

    # (c) how each record measures it, and how to read it
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    rows = [
        ("Aerial photos", "24 dated surveys, 1967-2023; one day each. Each side's position is the mean\n"
                          "distance from the wet/dry line to a fixed offshore datum line (25 points per\n"
                          "domain), as change since 1967. Gap = south's retreat − north's retreat."),
        ("CoastSat", "satellite shorelines from 1984, many images a year. Each transect's yearly mean\n"
                     "(≥ 3 images), about its own 1984-2016 mean, averaged over the side (6 transects\n"
                     "in GIS 6, 13 in GIS 5). Gap = north − south. With no 1967 survey, only its\n"
                     "changes compare, so it is moved onto the photos' 1967 baseline."),
        ("Reading it", "rising: the north holds while the south retreats (the groin trapping sand);\n"
                       "flat: the groin holds the gap but adds nothing; falling: the two sides converge\n"
                       "(sand getting past the groin, or a storm or fill on one side)."),
    ]
    yy = 1.0
    for head, body in rows:
        ax.text(0.0, yy, head, fontsize=7.3, fontweight="bold", va="top", transform=ax.transAxes)
        ax.text(0.0, yy - 0.075, body, fontsize=6.8, va="top", transform=ax.transAxes, linespacing=1.25)
        yy -= 0.12 + 0.07 * (body.count("\n") + 1)
    question(ax, 2, "How each record measures it")

    handles = [Line2D([], [], color=INK_MUTED, lw=1.0, ls=(0, (4, 2)), label="1967 shoreline"),
               Line2D([], [], color=INK, lw=1.4, label="1996 shoreline (photos)"),
               Line2D([], [], color=C["GROIN"], lw=3.0, label="groins (simplified)"),
               Line2D([], [], color=SIDE_COL["updrift"], lw=0.6, label="CoastSat transects, north"),
               Line2D([], [], color=SIDE_COL["downdrift"], lw=0.6, label="CoastSat transects, south"),
               Line2D([], [], color=C["LOCATOR"], lw=0.8, label="photo distance to the datum")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=6.8)
    png = out("0_what_the_gap_is")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** nothing; this figure explains the quantity the rest of the figures, and every groin fit, "
        "use. **How to read:** (a) Plan view at true scale, north up, ocean to the right, simplified: GIS 5 and "
        "GIS 6 drawn as straight 500 m domains with the groin field at the face between them (in reality the "
        "field lies just inside GIS 6), the 1967 shoreline dashed, and the 1996 shoreline at each domain's "
        f"photo-measured position: GIS 6 {n96:.0f} m and GIS 5 {s96:.0f} m from 1967. The gap (purple) is the "
        "difference between the two sides' positions. Teal: what the photos measure, the distance from the "
        "wet/dry line to a fixed offshore datum line; thin lines: the CoastSat transects whose positions are "
        "averaged per side (spacing schematic). (b) Each side's position from the dated aerial-photo "
        "wet/dry surveys, change since 1967, seaward positive (`../1-wetdry_photo_positions/"
        "Change_from_wetdry_1967_D2_D12.csv`; two 1972 surveys averaged); the shaded band between them is "
        "the gap. (c) How the photo and CoastSat records each measure it, and how to read a change in it. "
        f"**Shows:** the gap is not new beach: both sides retreated, the north side by {-n96:.0f} m and the "
        f"south side by {-s96:.0f} m by 1996, and the {n96 - s96:.0f} m between them is the groin's effect."))


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


# 0b. The same, animated: the two shorelines and the gap through the photo record
def gif_what_gap():
    from matplotlib.animation import FuncAnimation, PillowWriter
    s = photo_sides()
    years = np.arange(1967, s.index.max() + 0.001, 0.25)
    north = np.interp(years, s.index, s["north"])
    south = np.interp(years, s.index, s["south"])
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
    at.axvspan(2015, 2021, color=C["ADDED_FILL"], alpha=0.5, lw=0)
    at.plot(s.index, s["north"], "s", color=SIDE_COL["updrift"], ms=3, alpha=0.35)
    at.plot(s.index, s["south"], "s", color=SIDE_COL["downdrift"], ms=3, alpha=0.35)
    ln_n, = at.plot([], [], color=SIDE_COL["updrift"], lw=1.6, label=SIDE["updrift"])
    ln_s, = at.plot([], [], color=SIDE_COL["downdrift"], lw=1.6, label=SIDE["downdrift"])
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
        for yy, xx, side in ((700, n, "updrift"), (200, so, "downdrift")):
            if xx < -2:
                ax.annotate("", xy=(xx, yy), xytext=(0, yy), arrowprops=dict(arrowstyle="->", color=SIDE_COL[side],
                                                                              lw=1.2))
            ax.text(min(xx, 0) - 6, yy + 18, f"{xx:.0f} m", ha="right", fontsize=7.5, color=SIDE_COL[side],
                    fontweight="bold")
        gap = n - so
        if abs(gap) > 2:
            ax.plot([n, n], [500, 330], color=SIDE_COL["updrift"], lw=1.0, ls=(0, (2, 1.5)))
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
    gif = FIGS / "gap_across_groins_0_what_the_gap_is.gif"
    anim.save(gif, writer=PillowWriter(fps=12), dpi=110)
    plt.close(fig)
    print("wrote", gif)


# 1. Where each gap is measured
def fig_map():
    used = tab("transects_used")
    g = json.load(open(GROINS))
    dom = json.load(open(DOMAINS))
    mean = pd.read_csv(MEAN_LINE)
    xs = pd.concat([used["x0"], used["x1"]])
    ys = pd.concat([used["y0"], used["y1"]])
    xlim = (xs.min() - MAP_PAD_M, xs.max() + MAP_PAD_M)
    ylim = (ys.min() - MAP_PAD_M, ys.max() + MAP_PAD_M)
    fig = plt.figure(figsize=figsize("double", height=4.3), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=(1, 1, 0.5))
    axes = [fig.add_subplot(gs[0, k]) for k in range(2)]
    titles = ("Which transects make up GIS 6 and GIS 5?",
              f"Which transects lie within {NEAR_M:.0f} m of the groins?")
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        for f in dom["features"]:
            d = f["properties"]["domain_id"]
            if d not in (4, 5, 6, 7):
                continue
            ring = np.asarray(f["geometry"]["coordinates"][0])
            ax.fill(ring[:, 0], ring[:, 1], facecolor="0.95" if d % 2 else "0.90", edgecolor="0.7",
                    lw=0.5, zorder=0)
            cy = np.clip(ring[:, 1].mean(), *ylim)
            map_label(ax, xlim[0] + 0.12 * (xlim[1] - xlim[0]), cy, f"GIS {d}", text_kw=MAP_TEXT_DARK,
                      fontsize=7.5)
        m = mean[mean["included"].astype(bool)].sort_values("y")
        ax.plot(m["x"], m["y"], color=INK, lw=1.0, zorder=3)
        sel = used[used["gap"] == name]
        for side, s in sel.groupby("side"):
            for _, r in s.iterrows():
                ax.plot([r.x0, r.x1], [r.y0, r.y1], color=SIDE_COL[side], lw=1.6, zorder=4,
                        solid_capstyle="butt")
        for f in g["features"]:
            c = np.asarray(f["geometry"]["coordinates"])
            ax.plot(c[:, 0], c[:, 1], color=C["GROIN"], lw=2.4, zorder=6, solid_capstyle="butt")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        spines_for_image(ax)
        w, h = xlim[1] - xlim[0], ylim[1] - ylim[0]
        # net drift runs south here (option A waves); an arrow along the coast, offshore
        ax.annotate("", xy=(xlim[0] + 0.90 * w, ylim[0] + 0.30 * h), xytext=(xlim[0] + 0.90 * w, ylim[0] + 0.66 * h),
                    arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.2), zorder=7)
        ax.text(xlim[0] + 0.935 * w, ylim[0] + 0.48 * h, "net longshore drift", rotation=90, ha="left",
                va="center", fontsize=6.5, color=INK_MUTED)
        scale_bar_km(ax, length_m=500, segments=2, unit="m", x=0.08, y=0.05, text_kw=MAP_TEXT_DARK)
        north_dart(ax, (xlim[0] + 0.88 * w, ylim[0] + 0.13 * h), arrow_m=130, text_kw=MAP_TEXT_DARK)
        question(ax, i, titles[i], inside=True)

    # (c) where on the island
    ax = fig.add_subplot(gs[0, 2])
    o = gpd.read_file(ISLAND_OUTLINE).to_crs(UTM)
    o.plot(ax=ax, facecolor="0.86", edgecolor="0.55", lw=0.4)
    cx, cy = (xlim[0] + xlim[1]) / 2, (ylim[0] + ylim[1]) / 2
    ax.add_patch(Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0], fill=False,
                           edgecolor=C["LOCATOR"], lw=1.4, zorder=5))
    ax.plot(cx, cy, marker="s", ms=7, mfc="none", mec=C["LOCATOR"], mew=1.4, zorder=6)
    bx = o.total_bounds
    ax.set_xlim(bx[0] - 2000, bx[2] + 16000)
    ax.set_ylim(bx[1] - 2000, bx[3] + 2000)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    spines_for_image(ax)
    place_label(ax, cx - 2500, cy + 4500, "Buxton", text_kw=MAP_TEXT_DARK, fontsize=7.5, ha="right")
    place_label(ax, bx[0] + 15000, bx[1] + 4000, "Cape Point", text_kw=MAP_TEXT_DARK, fontsize=6.5)
    water_label(ax, bx[2] + 8000, cy + 30000, "Atlantic Ocean", text_kw=MAP_TEXT_DARK, rotation=90, fontsize=6.5)
    question(ax, 2, "Where on Hatteras?", inside=True)

    handles = [Line2D([], [], color=SIDE_COL["updrift"], lw=1.6, label="transects, north side (updrift)"),
               Line2D([], [], color=SIDE_COL["downdrift"], lw=1.6, label="transects, south side (downdrift)"),
               Line2D([], [], color=C["GROIN"], lw=2.4, label="Buxton groins"),
               Line2D([], [], color=INK, lw=1.0, label="CoastSat mean shoreline, Oct 1995-Oct 1997"),
               Patch(facecolor="0.92", edgecolor="0.7", label="model domains (GIS)"),
               Line2D([], [], color=C["LOCATOR"], lw=1.4, label="extent of (a) and (b)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=7.5)
    n = used.groupby(["gap", "side"]).size()
    png = out("1_transect_map")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** which CoastSat transects define each of the two gaps used in this analysis. "
        "**How to read:** UTM 18N, north up, ocean to the right; each coloured line is one CoastSat "
        "transect, blue on the north (updrift) side of the groins and gold on the south (downdrift) side; "
        "red bars are the four groins of the Buxton field; grey boxes are the model domains. "
        f"(a) The domain gap, the quantity the model is fitted on: every transect in GIS 6 "
        f"({n['domain', 'updrift']}) and GIS 5 ({n['domain', 'downdrift']}), the two domains the model's "
        f"groin sits between. (b) The near-field gap, a check: transects whose landward end lies within "
        f"{NEAR_M:.0f} m north of the northernmost groin ({n['near_field', 'updrift']}) or south of the "
        f"southernmost ({n['near_field', 'downdrift']}). (c) Locator: Hatteras Island, with the extent of "
        "(a) and (b) in teal. The black line is the CoastSat mean shoreline over the 1996 DEM window; the "
        "arrow is the direction of the wave-climate net longshore drift at Buxton (southward), which makes "
        "the north side updrift. **Shows:** the southernmost GIS 6 transect lies just south of the "
        "southern groin, so the domain gap is not purely across the field; the near-field gap has no such "
        "transect and gives the same break year (figure 3). The fourth groin line is the landward "
        "anti-flanking extension of the northern groin."))


# 2. Which side moved
def fig_sides():
    a = tab("gap_annual")
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8), sharex=True,
                             constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        s = a[a["gap"] == name]
        for side in ("updrift", "downdrift"):
            ax.plot(s["year"], s[f"{side}_anom_m"], color=SIDE_COL[side], lw=1.2, marker="o", ms=2.8)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1983, 2027)
        fill_shade(ax)
        events(ax)
        chart(ax)
        ax.set_ylabel(f"Shoreline position (m)\n+ seaward, about {FIT[0]}-{FIT[1]} mean")
        question(ax, i, f"Which side of the groins moved? {NAME[name][0].upper()}{NAME[name][1:]}")
    axes[-1].set_xlabel("Year")
    handles = [Line2D([], [], color=SIDE_COL[s], lw=1.2, marker="o", ms=2.8, label=SIDE[s])
               for s in ("updrift", "downdrift")] + [FILL_HANDLE]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = out("2_side_positions")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether a change in the gap comes from the north side holding (trapping) or from the "
        "south side moving. **How to read:** CoastSat annual shoreline position of each side, seaward "
        f"positive, each transect taken about its own {FIT[0]}-{FIT[1]} mean and averaged per side; the "
        "gap in figure 3 is blue minus gold. Dotted lines: Hurricane Gordon (1994), the last repair "
        "(1995), Hurricane Isabel (2003), the Buxton fill (2017); grey: years after the fill, left out of "
        f"every fit. (a) Whole domains, GIS 6 and GIS 5. (b) Transects within {NEAR_M:.0f} m of the "
        "groin field. **Shows:** the 1991-1995 widening of the gap is mostly the south side retreating "
        "(about 98 m in GIS 5 from 1984 to 1995, against about 24 m on the north side), with Gordon in "
        "that stretch, so the 1995 peak is partly storm loss downdrift rather than trapping. After 2017 "
        "the south side gains as the fill sand moves south past the groins."))


# 3. When the gap stopped widening
def fig_gap():
    a, h, ph, br = tab("gap_annual"), tab("hinge_fit"), tab("photo_gap"), tab("gap_breakpoint")
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.0), sharex=True,
                             constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        s, f = a[a["gap"] == name], h[h["gap"] == name]
        b = br.set_index("gap").loc[name]
        ax.axvspan(b.break_ci90_lo, b.break_ci90_hi, color=COL[name], alpha=0.10, lw=0, zorder=0)
        ax.plot(f["year"], f["fit_m"], color=COL[name], lw=3.0, alpha=0.35, solid_capstyle="round")
        ax.plot(s["year"], s["gap_m"], color=COL[name], lw=1.2, marker="o", ms=3)
        if name == "domain":
            p = ph[ph["year"] >= 1984]
            ax.plot(p["year"], p["shifted_m"], ls="none", marker="s", ms=4, color=COL["photos"],
                    mec="white", mew=0.5, zorder=5)
        ax.axvline(b.break_year, color=COL[name], lw=0.9, ls="--")
        ax.text(b.break_year + 0.3, 0.05, f"break {b.break_year}\n90% range {b.break_ci90_lo}-{b.break_ci90_hi}",
                transform=ax.get_xaxis_transform(), fontsize=7, color=COL[name])
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1983, 2027)
        fill_shade(ax)
        events(ax, years=(1995, 2003, 2017))
        chart(ax)
        ax.set_ylabel(GAP_LABEL)
        question(ax, i, f"When did the gap stop widening? {NAME[name][0].upper()}{NAME[name][1:]}")
    axes[-1].set_xlabel("Year")
    handles = [Line2D([], [], color=COL["domain"], lw=1.2, marker="o", ms=3, label="CoastSat gap, annual mean"),
               Line2D([], [], color=INK_MUTED, lw=3.0, alpha=0.5, label=f"best one-break fit, {FIT[0]}-{FIT[1]}"),
               Patch(facecolor=COL["domain"], alpha=0.15, label="90% range of the break year"),
               Line2D([], [], ls="none", marker="s", ms=4, color=COL["photos"], mec="white",
                      label="aerial-photo gap, aligned to CoastSat"),
               FILL_HANDLE]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = out("3_gap_and_break")
    save(fig, png, close=True)
    d, n = br.set_index("gap").loc["domain"], br.set_index("gap").loc["near_field"]
    record_caption(png, (
        "**Tests:** in which year the groin stopped holding the north side ahead of the south, i.e. when "
        f"the gap stopped widening. **How to read:** {WHAT} Pale thick line: a continuous line with one "
        f"change of slope (a hinge), fitted on {FIT[0]}-{FIT[1]}; dashed line: the year of that change, "
        f"with its 90% range from {N_BOOT} residual bootstraps shaded. Grey: after the 2017 fill, not "
        f"fitted. (a) Whole domains: {d.rate_before_m_yr:+.1f} m/yr to {d.break_year}, then "
        f"{d.rate_after_m_yr:+.1f} m/yr. Red squares: the gap from the dated aerial-photo wet/dry lines "
        "(GIS 5 minus GIS 6 change since 1967, landward positive, the same sense), moved up or down to "
        "the CoastSat mean over the shared pre-fill years because the two datums differ, so only the "
        f"shape compares. (b) Within {NEAR_M:.0f} m of the field: {n.rate_before_m_yr:+.1f} then "
        f"{n.rate_after_m_yr:+.1f} m/yr, break {n.break_year}. **Shows:** both gaps stop widening in "
        "1995, the year of the last repair, not at Hurricane Isabel (2003); what follows is a slow "
        "decline, not a collapse."))


# 4. How certain the break year is
def fig_break():
    prof, boots = tab("break_profile"), tab("break_bootstrap")
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.2), constrained_layout=True)
    ax = axes[0]
    for name in ("domain", "near_field"):
        p = prof[prof["gap"] == name]
        ax.plot(p["break_year"], p["sse_over_best"], color=COL[name], lw=1.4)
    ax.set_ylabel("Fit error, relative to\nthe best year (1 = best)")
    ax.set_xlabel("Candidate break year")
    events(ax, years=(1995, 2003))
    chart(ax)
    question(ax, 0, "Which break year fits best?", inside=True)
    ax = axes[1]
    years = np.arange(BREAK_RANGE[0], BREAK_RANGE[1] + 1)
    for k, name in enumerate(("domain", "near_field")):
        share = boots[boots["gap"] == name]["break_year"].value_counts(normalize=True)
        ax.bar(years + (k - 0.5) * 0.4, 100 * share.reindex(years, fill_value=0), width=0.4, color=COL[name])
    ax.axvspan(2001.5, 2005.5, color=C["BASE_FILL"], alpha=0.6, lw=0, zorder=0)
    ax.text(2003.5, 0.97, "model's\nfailure step", transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=6.5, color=INK_MUTED)
    ax.set_ylabel(f"Share of {N_BOOT:,} bootstrap fits (%)")
    ax.set_xlabel("Break year")
    chart(ax)
    question(ax, 1, "How certain is that year?", inside=True)
    handles = [Patch(color=COL[n], label=f"CoastSat, {NAME[n]}") for n in ("domain", "near_field")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    png = out("4_break_year")
    save(fig, png, close=True)
    br = tab("gap_breakpoint").set_index("gap")
    record_caption(png, (
        "**Tests:** how sharply the data pick out the break year, and whether the timing the model "
        "assumes (a failure around 2004) is compatible with them. **How to read:** (a) For every "
        f"candidate break year, the misfit (sum of squared residuals) of the one-break fit to the "
        f"{FIT[0]}-{FIT[1]} annual gap, divided by the misfit of the best year, so 1 is the best fit and "
        "higher is worse. (b) The best break year in each of "
        f"{N_BOOT} residual bootstraps, as a share of all bootstraps. Grey band: 2002-2005, around the "
        "model's 2004 failure step. **Shows:** both gaps fit best with a 1995 break; 2002-2005 holds "
        f"{100 * br.loc['domain', 'boot_share_2002_2005']:.0f}% (whole domains) and "
        f"{100 * br.loc['near_field', 'boot_share_2002_2005']:.0f}% (near field) of the bootstraps, "
        f"while 1997 or earlier holds {100 * br.loc['domain', 'boot_share_le_1997']:.0f}% and "
        f"{100 * br.loc['near_field', 'boot_share_le_1997']:.0f}%."))


# 5. The gap's rate in each era
def fig_eras():
    e = tab("era_rates")
    fig, ax = plt.subplots(figsize=figsize("single", height=3.4), constrained_layout=True)
    eras = [f"{lo}-{hi}" for lo, hi in ERAS]
    for k, src in enumerate(("domain", "near_field", "photos")):
        s = e[e["source"] == src].set_index("era").reindex(eras)
        x = np.arange(len(eras)) + (k - 1) * 0.22
        ok = s["rate_m_yr"].notna().to_numpy()
        lab = NAME[src] if src == "photos" else f"CoastSat, {NAME[src]}"
        ax.errorbar(x[ok], s["rate_m_yr"][ok], yerr=s["se_m_yr"][ok], fmt="s" if src == "photos" else "o",
                    color=COL[src], ms=4.5, capsize=2.5, lw=1.0, label=lab)
        for xi, (_, r) in zip(x, s.iterrows()):
            if not np.isfinite(r.rate_m_yr):
                ax.text(xi, 0.4, "too few\nsurveys", ha="center", va="bottom", fontsize=6.3, color=COL[src])
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(range(len(eras)), [f"before the last\nrepair\n{eras[0]}", f"repair to\nIsabel\n{eras[1]}",
                                     f"Isabel to\nthe fill\n{eras[2]}"])
    ax.set_ylabel("Rate of gap change (m/yr)\n+ = widening")
    chart(ax)
    ax.set_title("Was the gap widening or narrowing?")
    ax.legend(frameon=False, fontsize=6.8, loc="upper right")
    png = out("5_era_rates")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether the gap widened (groin trapping), held, or narrowed in each era of the "
        "structure's history, and whether CoastSat and the aerial photos agree. **How to read:** the "
        "least-squares rate of the gap in each era, ± 1 standard error; above zero the north side is "
        "gaining on the south. Eras: up to the 1995 last repair, from 1996 to Hurricane Isabel (2003), and "
        "from 2004 to the 2017 fill. CoastSat on annual means; photos on the dated wet/dry surveys inside "
        "each era (only 1996 and 1997 fall in 1996-2003, too few for a rate). **Shows:** both sources "
        "show the gap widening before the repair. From 1996 to 2003 CoastSat shows a slight narrowing "
        "that is not distinguishable from zero. After 2004 they disagree: CoastSat is flat while the "
        "photos fall, a fall that starts from the single high 2004 survey (figure 7)."))


# 6. Whether the break survives other fits
def fig_robust():
    r = tab("gap_robustness")
    checks = list(dict.fromkeys(r["check"]))
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.0), sharey=True,
                             constrained_layout=True)
    for k, name in enumerate(("domain", "near_field")):
        s = r[r["gap"] == name].set_index("check").loc[checks]
        y = np.arange(len(checks)) + (k - 0.5) * 0.25
        axes[0].errorbar(s["break_year"], y, xerr=[s["break_year"] - s["break_ci90_lo"],
                                                   s["break_ci90_hi"] - s["break_year"]],
                         fmt="o", color=COL[name], ms=4, capsize=2.5, lw=1.0)
        axes[1].errorbar(s["step_m"], y, xerr=s["step_se_m"], fmt="o", color=COL[name], ms=4,
                         capsize=2.5, lw=1.0)
    axes[0].axvspan(2001.5, 2005.5, color=C["BASE_FILL"], alpha=0.6, lw=0, zorder=0)
    axes[0].text(2003.5, 0.02, "model's\nfailure step", transform=axes[0].get_xaxis_transform(),
                 ha="center", va="bottom", fontsize=6.5, color=INK_MUTED)
    axes[0].set_xlim(1988, 2008)
    axes[0].set_xticks(range(1990, 2009, 5))
    axes[0].set_xlabel("Break year, with 90% bootstrap range")
    axes[1].axvline(0, color=INK_MUTED, lw=0.6)
    axes[1].set_xlabel(f"Size of an extra drop at {STEP_YEAR} (m, ± 1 SE)")
    axes[0].set_yticks(range(len(checks)), [CHECK.get(c, c) for c in checks])
    axes[0].invert_yaxis()
    for i, (ax, t) in enumerate(zip(axes, ("Does the break year hold under other fits?",
                                           f"Is there a separate drop at {STEP_YEAR}?"))):
        ax.grid(axis="x")
        open_frame(ax)
        question(ax, i, t, inside=True)
    handles = [Line2D([], [], color=COL[n], marker="o", ms=4, ls="none", label=f"CoastSat, {NAME[n]}")
               for n in ("domain", "near_field")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    png = out("6_robustness")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether the 1995 break depends on choices in the fit, and whether the data also "
        f"need a sudden drop at {STEP_YEAR}, the year the model's groin fails. **How to read:** each row "
        f"is one version of the fit: the baseline ({FIT[0]}-{FIT[1]}); starting in 1988, which drops the "
        "sparse early Landsat years (5-10 images a year); leaving out the 1995 peak; ending in 2013. "
        f"(a) The best break year with its 90% range from {N_BOOT_CHECK} residual bootstraps; grey: "
        f"2002-2005. (b) The size of a level step at {STEP_YEAR} fitted alongside the free break, ± 1 "
        "standard error; a step that matters would sit clearly away from zero. **Shows:** the break stays "
        "at 1994-1996 in every version. The step is negative but within two standard errors of zero in "
        "all eight fits, and adding it lowers the BIC below the break-only fit in only one (whole domains "
        "with 1995 left out), by 0.1 (tables/gap_robustness.csv)."))


# 7. Where the photos and CoastSat disagree
def fig_photos():
    a, ph = tab("gap_annual"), tab("photo_gap")
    d = a[a["gap"] == "domain"]
    p = ph[ph["year"] >= 1984]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.2), sharex=True,
                             constrained_layout=True, gridspec_kw=dict(height_ratios=(1.6, 1)))
    ax = axes[0]
    ax.plot(d["year"], d["gap_m"], color=COL["domain"], lw=1.2, marker="o", ms=3)
    ax.plot(p["year"], p["shifted_m"], ls="none", marker="s", ms=4.5, color=COL["photos"], mec="white",
            mew=0.5, zorder=5)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1983, 2027)
    fill_shade(ax)
    events(ax, years=(1995, 2003, 2017))
    chart(ax)
    ax.set_ylabel(GAP_LABEL)
    question(ax, 0, "Do the aerial photos and CoastSat agree? Whole domains")
    ax = axes[1]
    q = p.dropna(subset=["photo_minus_coastsat_m"])
    cols = [COL["photos"] if y == 2004 else C["BASE"] for y in q["year"]]
    ax.bar(q["year"], q["photo_minus_coastsat_m"], width=0.7, color=cols)
    v = q.set_index("year").loc[2004, "photo_minus_coastsat_m"]
    ax.text(2004.6, v, "2004 survey", fontsize=7, color=COL["photos"], va="top")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    fill_shade(ax)
    chart(ax)
    ax.set_ylabel("Photo minus\nCoastSat (m)")
    ax.set_xlabel("Year")
    question(ax, 1, "By how much does each photo differ?")
    handles = [Line2D([], [], color=COL["domain"], lw=1.2, marker="o", ms=3, label="CoastSat gap, annual mean"),
               Line2D([], [], ls="none", marker="s", ms=4.5, color=COL["photos"], mec="white",
                      label="aerial-photo gap, aligned to CoastSat"), FILL_HANDLE]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = out("7_photos_vs_coastsat")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether the aerial-photo record, which the model's current failure timing rests on, "
        "agrees with CoastSat, and if not, where. **How to read:** (a) The CoastSat annual gap across the "
        "groins (whole domains) and the gap from each dated wet/dry survey, the photos moved up or down "
        f"to the CoastSat mean over the shared {FIT[0]}-{FIT[1]} years because their datums differ, so "
        "only the shape compares. (b) Photo minus CoastSat in each photo year; bars near zero mean the "
        f"two agree. **Shows:** before the fill they agree to within 13 m in every year but one (within "
        f"16 m after it). The 2004 survey "
        f"sits {v:+.0f} m above CoastSat's 2004 mean, the largest difference before the fill. That single "
        "survey is the basis for reading the gap as held until 2004, and the date that pinned the "
        "blocking strength in the 2026-10-05 calibration fit. A photo is one day; a CoastSat value is a "
        "year's mean."))


# 8. The groin's state in 1996, the start of the hindcast
def fig_state_1996():
    a, ph, e, br = tab("gap_annual"), tab("photo_gap"), tab("era_rates"), tab("gap_breakpoint")
    d = a[a["gap"] == "domain"].set_index("year")
    ph = ph.set_index("year")
    shift = (ph["shifted_m"] - ph["gap_change_since_1967_m"]).mean()    # CoastSat onto the photos' 1967 datum
    cs = d["gap_m"] - shift
    photos = ph["gap_change_since_1967_m"]
    peak_yr = int(photos.idxmax())
    b = br.set_index("gap").loc["domain"]
    fig = plt.figure(figsize=figsize("double", height=6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=(1.35, 1))

    # (a) the gap over the groin's whole life, on the photos' 1967 datum
    ax = fig.add_subplot(gs[0, :])
    for (lo, hi), name in (((1967, 1978), "building"), ((1978, 1995), "slow growth"),
                           ((1995, 2016), "after the last repair")):
        ax.axvspan(lo, hi, color=C["BASE_FILL"] if name == "slow growth" else "white", alpha=0.5, lw=0, zorder=0)
        ax.text((lo + hi) / 2, 0.03, name, transform=ax.get_xaxis_transform(), ha="center", fontsize=6.8,
                color=INK_MUTED)
    ax.axvspan(2015, 2021, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
    ax.plot(cs.index, cs.values, color=COL["domain"], lw=1.1, marker="o", ms=2.8)
    ax.plot(photos.index, photos.values, ls="none", marker="s", ms=4.5, color=COL["photos"], mec="white",
            mew=0.5, zorder=5)
    ax.axvline(1996, color=INK, lw=1.4)
    g96 = photos.loc[1996]
    ax.annotate(f"1996: gap {g96:.0f} m (photos), {cs.loc[1996]:.0f} m (CoastSat)\n"
                f"peak {photos.max():.0f} m in {peak_yr}", xy=(1996, g96), xytext=(1998.5, 30),
                fontsize=7, arrowprops=dict(arrowstyle="-", color=INK, lw=0.6), va="bottom")
    for yr, lab in ((1969, "built"), (1975, "repaired"), (1995, "last repair"), (2003, "Isabel"),
                    (2017, "Buxton fill")):
        ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":", zorder=1)
        ax.text(yr + 0.3, 0.985, lab, transform=ax.get_xaxis_transform(), fontsize=6.5, color=INK_MUTED,
                va="top")
    ax.text(2019.2, 0.80, "decay\nreported", transform=ax.get_xaxis_transform(), ha="center", fontsize=6.5,
            color=INK_MUTED)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1965, 2027)
    ax.set_ylim(-10, 215)
    ax.set_ylabel("Gap since 1967 (m)\n+ = north side seaward")
    ax.set_xlabel("Year")
    chart(ax)
    question(ax, 0, "Where did the gap stand in 1996, in the groin's whole life?")

    # (b) the rate of gap change in each stage
    ax = fig.add_subplot(gs[1, 0])
    r_build, se_build, _ = era_rate(photos, 1967, 1978)
    ce = e[e["source"] == "domain"].set_index("era")
    rows = [("1967-1978\nbuilding\n(photos)", r_build, se_build, COL["photos"]),
            ("1984-1995\nup to the\nlast repair", ce.loc["1984-1995", "rate_m_yr"], ce.loc["1984-1995", "se_m_yr"],
             COL["domain"]),
            ("1996-2003\nfrom the\nhindcast start", ce.loc["1996-2003", "rate_m_yr"],
             ce.loc["1996-2003", "se_m_yr"], INK),
            ("2004-2016\nIsabel to\nthe fill", ce.loc["2004-2016", "rate_m_yr"], ce.loc["2004-2016", "se_m_yr"],
             COL["domain"])]
    x = np.arange(len(rows))
    ax.bar(x, [r[1] for r in rows], yerr=[r[2] for r in rows], color=[r[3] for r in rows], width=0.6,
           capsize=2.5, error_kw=dict(lw=0.8))
    for xi, (_, r, se, _) in zip(x, rows):
        ax.text(xi, r + np.sign(r) * (se + 0.4) if r else 0.4, f"{r:+.1f}", ha="center",
                va="bottom" if r >= 0 else "top", fontsize=7)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(x, [r[0] for r in rows], fontsize=6.6)
    ax.set_ylim(-6, 15.5)
    ax.set_ylabel("Rate of gap change (m/yr)\n± 1 SE; + = widening")
    chart(ax)
    question(ax, 1, "Was the gap still widening in 1996?", inside=True)

    # (c) the state in 1996, in words
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    _letter_inside(ax, 2)
    ax.set_title("The groin in 1996", loc="center")
    lines = [
        ("Structure", "south groin repaired 1995 (last repair);\nnorth groin gone since 1972;\n"
                      "maintenance stops after 1999"),
        ("Gap", f"near its full size: {g96:.0f} m since 1967,\n{100 * g96 / photos.max():.0f}% of the "
                f"{photos.max():.0f} m peak ({peak_yr})"),
        ("Trend", f"no longer widening: the break is {b.break_year}\n(90% range {b.break_ci90_lo}-"
                  f"{b.break_ci90_hi}); after it {b.rate_after_m_yr:+.1f} m/yr"),
        ("Not failed", "the gap fades slowly; no separate\ndrop at Isabel (2003-04)"),
        ("Caveat", "the 1995 peak is partly Gordon's\n(1994) loss on the south side"),
    ]
    y = 0.86
    for head, body in lines:
        ax.text(0.02, y, head, fontsize=7.4, fontweight="bold", va="top", transform=ax.transAxes)
        ax.text(0.30, y, body, fontsize=7.2, va="top", transform=ax.transAxes)
        y -= 0.035 + 0.058 * (body.count("\n") + 1)
    ax.text(0.02, 0.0, "In short: a full fillet, held in place,\nno longer growing and starting to fade.",
            fontsize=7.4, va="bottom", style="italic", transform=ax.transAxes)

    handles = [Line2D([], [], color=COL["domain"], lw=1.1, marker="o", ms=2.8,
                      label="CoastSat gap, annual mean (moved onto the photos' 1967 datum)"),
               Line2D([], [], ls="none", marker="s", ms=4.5, color=COL["photos"], mec="white",
                      label="aerial-photo gap (wet/dry line), change since 1967"),
               Patch(facecolor=C["ADDED_FILL"], alpha=0.5, label="structural losses reported, 2015-2021")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("8_state_in_1996")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** what state the Buxton groin was in at 1996, the start of the calibration (hindcast) "
        "period, from the observations alone. **How to read:** (a) The gap across the groins (GIS 6 "
        "seaward of GIS 5) over the groin's whole life: the dated aerial-photo wet/dry surveys as change "
        "since 1967 (squares), and the CoastSat annual gap (circles) moved onto the same datum by its "
        f"mean offset from the photos over the shared {FIT[0]}-{FIT[1]} years ({-shift:.0f} m), so the "
        "two read on one scale. Black line: 1996. Dotted lines: the structure's recorded events "
        "(`../../data/structure_history.md`); orange shading: the structural losses reported for 2015-2021. "
        "(b) The least-squares rate of gap change in each stage, ± 1 standard error: the photos for "
        "1967-1978 (CoastSat starts in 1984), CoastSat for the rest (tables/era_rates.csv). (c) The "
        f"state in 1996 in words, from (a), (b), figures 2-6 and the structure record. **Shows:** in "
        f"1996 the gap stood at {g96:.0f} m (photos), {100 * g96 / photos.max():.0f}% of its {photos.max():.0f} m "
        f"peak in {peak_yr}, a year after the last repair. It had built at {r_build:+.1f} m/yr in the "
        f"first decade and was still widening at {ce.loc['1984-1995', 'rate_m_yr']:+.1f} m/yr to 1995, "
        f"but from 1996 to 2003 it changed by {ce.loc['1996-2003', 'rate_m_yr']:+.1f} ± "
        f"{ce.loc['1996-2003', 'se_m_yr']:.1f} m/yr, not distinguishable from zero. So in 1996 the groin "
        "held a full-sized gap but no longer widened it, and nothing in the record or the data marks a "
        "collapse before the 2015-2021 losses."))


def main():
    apply_style()
    for f in (fig_what_gap, gif_what_gap, fig_map, fig_sides, fig_gap, fig_break, fig_eras, fig_robust, fig_photos, fig_state_1996):
        f()
        print(f"drew {f.__name__}")


if __name__ == "__main__":
    main()
