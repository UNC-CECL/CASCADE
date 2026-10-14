#!/usr/bin/env python3
"""
The photos-vs-CoastSat figures, drawn from the tables photos_vs_coastsat.py writes
and the two single-source analyses' tables.

    python photos_vs_coastsat_figures.py   ->  figures/photos_vs_coastsat_<n>_<what>.png

One figure per step, in the order of the analysis: (1) what the gap is and how each
record measures it, (2) putting the two records on one scale, (3) the gap from both and
how far each photo is from CoastSat, (4) whether their era rates agree, (5) what
the two together say about the groin in 1996, and (6) the gap across both candidate model
boundaries (GIS 5/6 and 6/7) on one scale. The plotting helpers come from
../2-coastsat_gap/coastsat_gap_figures.py. Captions go to figures/supporting/CAPTIONS.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(REPO / "scripts"), str(HERE), str(HERE.parent / "2-coastsat_gap"),
                str(HERE.parent / "1-aerial_photo_gap")]
from coastsat_gap_across_groins import BOUNDARIES, ERAS, FIT  # noqa: E402
from coastsat_gap_figures import (  # noqa: E402
    COL, FILL_HANDLE, GAP_LABEL, NAME, SIDE, SIDE_COL, chart, events, fill_shade, question)
from photo_gap import era_rate  # noqa: E402
from photos_vs_coastsat import KEYS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
FIGS = HERE / "figures"
TABLES = HERE / "tables"
PHOTO_TABLES = HERE.parent / "1-aerial_photo_gap" / "tables"
COASTSAT_TABLES = HERE.parent / "2-coastsat_gap" / "tables"
C_CAL, C_TEST = C_1984, C_1997      # calibration red, test blue
# -----------------------------------------------------------------------------


def tab(name):
    for d in (TABLES, COASTSAT_TABLES, PHOTO_TABLES):
        if (d / f"{name}.csv").exists():
            return pd.read_csv(d / f"{name}.csv")
    raise FileNotFoundError(name)


def out(name):
    return FIGS / f"photos_vs_coastsat_{name}.png"


# Each side's wet/dry position per survey year, change since 1967, seaward positive
def photo_sides():
    s = tab("photo_gap_annual").set_index("year")
    return s.rename(columns={"north_m": "north", "south_m": "south"})[["north", "south"]]


# 1. What the gap is, and how each record measures it
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
    png = out("1_what_the_gap_is")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** nothing; this figure explains the quantity every groin figure and fit uses, and how "
        "the photos and CoastSat each measure it. **How to read:** (a) Plan view at true scale, north up, ocean to the right, simplified: GIS 5 and "
        "GIS 6 drawn as straight 500 m domains with the groin field at the face between them (in reality the "
        "field lies just inside GIS 6), the 1967 shoreline dashed, and the 1996 shoreline at each domain's "
        f"photo-measured position: GIS 6 {n96:.0f} m and GIS 5 {s96:.0f} m from 1967. The gap (purple) is the "
        "difference between the two sides' positions. Teal: what the photos measure, the distance from the "
        "wet/dry line to a fixed offshore datum line; thin lines: the CoastSat transects whose positions are "
        "averaged per side (spacing schematic). (b) Each side's position from the dated aerial-photo "
        "wet/dry surveys, change since 1967, seaward positive (`../1-aerial_photo_gap/tables/"
        "photo_gap_annual.csv`; two 1972 surveys averaged); the shaded band between them is "
        "the gap. (c) How the photo and CoastSat records each measure it, and how to read a change in it. "
        f"**Shows:** the gap is not new beach: both sides retreated, the north side by {-n96:.0f} m and the "
        f"south side by {-s96:.0f} m by 1996, and the {n96 - s96:.0f} m between them is the groin's effect."))


# The structure's recorded events, for the whole-life figures
LIFE_EVENTS = ((1969, "built"), (1975, "repaired"), (1994, "Gordon"), (1995, "last repair"), (2003, "Isabel"),
               (2017, "Buxton fill"))
LOSSES = (2015, 2021)


def life_events(ax):
    for yr, lab in LIFE_EVENTS:
        ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":", zorder=1)
        ax.text(yr + (-0.3 if yr == 1994 else 0.3), 0.985, lab, transform=ax.get_xaxis_transform(), fontsize=6.5,
                color=INK_MUTED, va="top", ha="right" if yr == 1994 else "left")


# 2. Putting the two records on one scale
def fig_align():
    ph = tab("photo_gap_annual").set_index("year")
    cs = tab("coastsat_sides_on_1967_datum").set_index("year")
    a = tab("gap_annual")
    an = a[a["gap"] == "domain"].set_index("year")
    shift = {"north": (cs["north_m"] - an["updrift_anom_m"]).mean(), "south": (cs["south_m"] - an["downdrift_anom_m"]).mean()}
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.4), sharex=True, constrained_layout=True)
    for i, ax in enumerate(axes):
        ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
        for side, key, col in (("updrift", "north", "updrift_anom_m"), ("downdrift", "south", "downdrift_anom_m")):
            y = an[col] if i == 0 else cs[f"{key}_m"]
            ax.plot(y.index, y.values, color=SIDE_COL[side], lw=1.2, marker="o", ms=2.6, zorder=3)
            ax.plot(ph.index, ph[f"{key}_m"], ls="none", marker="s", ms=4.2, color=SIDE_COL[side], mec="white",
                    mew=0.6, zorder=4)
        if i == 1:
            ax.fill_between(cs.index, cs["south_m"], cs["north_m"], color=C["ACCENT_FILL"], alpha=0.5, lw=0)
            for key, side in (("north", "updrift"), ("south", "downdrift")):
                y1 = cs[f"{key}_m"].iloc[0]
                ax.text(cs.index[0] - 0.8, y1, f"CoastSat {key}\nmoved {shift[key]:+.0f} m", ha="right", va="center",
                        fontsize=6.6, color=SIDE_COL[side])
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1965, 2027)
        ax.set_ylabel("Shoreline position (m)\n+ seaward")
        life_events(ax)
        chart(ax)
    question(axes[0], 0, "Before: each record has its own zero")
    question(axes[1], 1, "After: CoastSat moved onto the photos' 1967 datum")
    axes[1].set_xlabel("Year")
    handles = [Line2D([], [], color=SIDE_COL["updrift"], lw=1.2, label=SIDE["updrift"] + ", GIS 6"),
               Line2D([], [], color=SIDE_COL["downdrift"], lw=1.2, label=SIDE["downdrift"] + ", GIS 5"),
               Line2D([], [], color=INK_MUTED, lw=1.2, marker="o", ms=2.6, label="CoastSat, annual mean"),
               Line2D([], [], ls="none", marker="s", ms=4.2, color=INK_MUTED, mec="white",
                      label="aerial photos, change since 1967"),
               Patch(facecolor=C["ACCENT_FILL"], alpha=0.5, label="the gap between the sides"),
               Patch(facecolor=C["ADDED_FILL"], alpha=0.5, label="structural losses reported, 2015-2021")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("2_one_scale")
    save(fig, png, close=True)
    shared = ph.index.intersection(an.index)
    shared = shared[(shared >= FIT[0]) & (shared <= FIT[1])]
    record_caption(png, (
        "**Tests:** nothing; this shows how the two records are put on one scale before they are compared. "
        "**How to read:** squares: each side's position in the dated aerial-photo surveys, change since the 1967 "
        "survey (`../1-aerial_photo_gap/`); lines: each side's CoastSat annual mean (`../2-coastsat_gap/` "
        f"figure 3). (a) As each record comes: CoastSat is centred on its own {FIT[0]}-{FIT[1]} mean, so its zero "
        "has nothing to do with 1967 and the lines sit far from the squares. (b) Each CoastSat side moved by one "
        f"constant, the mean difference from the same side's photos over the {len(shared)} photo years in "
        f"{FIT[0]}-{FIT[1]} (labelled at each line's start: north {shift['north']:+.0f} m, south {shift['south']:+.0f} m). The shift "
        "changes only the level, never the shape, and the gap between the two shifted sides is the CoastSat "
        "gap on the photos' scale. **Shows:** CoastSat starts in 1984, 15 years into the groin's life, so its "
        "level has to come from the photos; once shifted, the annual lines run through the photo squares on "
        "both sides. How closely is figure 3."))


# 3. The gap from both records, and how far each photo is from CoastSat
def fig_two_sides():
    ph = tab("photo_gap_annual").set_index("year")
    cs = tab("coastsat_sides_on_1967_datum").set_index("year")
    res = tab("photo_minus_coastsat_by_side").set_index("year")
    spread = tab("coastsat_one_day_spread").set_index("year")
    sd_cols = [f"{k}_one_day_sd_m" for k in ("north", "south", "gap")]
    band = 2 * float(np.median(spread[sd_cols].to_numpy()))       # 2 x the typical one-day SD, all three pooled
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True, constrained_layout=True)
    ax = axes[0]
    ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
    ax.fill_between(cs.index, 0, cs["gap_m"], color=C["ACCENT_FILL"], alpha=0.5, lw=0)
    ax.plot(cs.index, cs["gap_m"], color=C["ACCENT"], lw=1.2, marker="o", ms=2.6, zorder=3)
    ax.plot(ph.index, ph["gap_m"], ls="none", marker="s", ms=4.2, color=C["ACCENT"], mec="white", mew=0.6, zorder=4)
    for yr, (dx, dy, ha) in ((1996, (-5.0, 22, "right")), (int(cs.index.max()), (-1.0, -30, "right"))):
        g = cs.loc[yr, "gap_m"]
        ax.annotate(f"{g:.0f} m\n{yr}", xy=(yr, g), xytext=(yr + dx, g + dy), ha=ha, va="center", fontsize=6.8,
                    color=C["ACCENT"], arrowprops=dict(arrowstyle="-", color=C["ACCENT"], lw=0.5, shrinkA=0,
                                                       shrinkB=2, relpos=(1, 0.5) if yr == 1996 else (0.5, 1)))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_ylim(-5, 190)
    ax.set_ylabel("Gap (m)\nnorth minus south")
    life_events(ax)
    chart(ax)
    ax.set_xlim(1965, 2027)
    question(ax, 0, "How wide was the gap, in both records?")

    # (c) photo minus CoastSat, against the scatter one day shows about its year's mean
    ax = axes[1]
    ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
    ax.axhspan(-band, band, color=C["BASE_FILL"], alpha=0.7, lw=0, zorder=0)
    ax.text(1985.5, band - 2, f"± {band:.0f} m: what one day's\nscatter about a year's mean explains",
            fontsize=6.5, color=INK_MUTED, va="top", ha="right")
    keys = (("north", SIDE_COL["updrift"], "o", -0.25), ("gap", C["ACCENT"], "s", 0.0),
            ("south", SIDE_COL["downdrift"], "o", 0.25))
    for key, col, mk, dx in keys:
        v = res[f"{key}_photo_minus_coastsat_m"]
        outside = v.abs() > band
        ax.plot(res.index + dx, v, ls="none", marker=mk, ms=4.2, color=col, mec="white", mew=0.5, zorder=4)
        ax.plot(res.index[outside] + dx, v[outside], ls="none", marker=mk, ms=4.2, color=col, mec=INK, mew=0.9,
                zorder=5)
    v = res.loc[2004, "gap_photo_minus_coastsat_m"]
    ax.annotate(f"2004 gap {v:+.0f} m", xy=(2004, v), xytext=(2006.5, v), fontsize=6.8, color=C["ACCENT"],
                va="center", arrowprops=dict(arrowstyle="-", color=C["ACCENT"], lw=0.5))
    ax.set_ylim(-52, 62)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_ylabel("Photo minus\nCoastSat (m)")
    ax.set_xlabel("Year")
    life_events(ax)
    chart(ax)
    question(ax, 1, "How far is each photo from CoastSat?")

    handles = [Line2D([], [], color=SIDE_COL["updrift"], ls="none", marker="o", ms=4.2, mec="white",
                      label="(b) " + SIDE["updrift"] + ", GIS 6"),
               Line2D([], [], color=SIDE_COL["downdrift"], ls="none", marker="o", ms=4.2, mec="white",
                      label="(b) " + SIDE["downdrift"] + ", GIS 5"),
               Line2D([], [], color=INK_MUTED, lw=1.2, marker="o", ms=2.6,
                      label="CoastSat, annual mean, moved onto the photos (figure 2)"),
               Line2D([], [], ls="none", marker="s", ms=4.2, color=INK_MUTED, mec="white",
                      label="aerial photos, dated surveys"),
               Patch(facecolor=C["ACCENT_FILL"], alpha=0.5, label="the gap: the groin's effect"),
               Patch(facecolor=C["ADDED_FILL"], alpha=0.5, label="structural losses reported, 2015-2021"),
               Line2D([], [], ls="none", marker="s", ms=4.2, color=C["ACCENT"], mec="white", label="(b) the gap"),
               Line2D([], [], ls="none", marker="o", ms=4.2, color=INK_MUTED, mec=INK, mew=0.9,
                      label="(b) outlined: outside the grey band")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("3_gap_and_agreement")
    save(fig, png, close=True)
    shared = ph.index.intersection(cs.index)
    shared = shared[(shared >= FIT[0]) & (shared <= FIT[1])]
    shift = {k: (ph.loc[shared, f"{k}_m"] - cs.loc[shared, f"{k}_m"]).mean() for k in ("north", "south")}
    last = int(cs.index.max())
    pre, post = res.index <= FIT[1], res.index > FIT[1]
    keys = ("north", "south", "gap")
    inside = {k: (res[f"{k}_photo_minus_coastsat_m"].abs() <= band) for k in keys}
    rms = {k: float(np.sqrt((res.loc[pre, f"{k}_photo_minus_coastsat_m"] ** 2).mean())) for k in keys}
    after = {k: res.loc[post, f"{k}_photo_minus_coastsat_m"].mean() for k in keys}
    z04 = res.loc[2004, "gap_photo_minus_coastsat_m"] / spread.loc[2004, "gap_one_day_sd_m"]
    r04 = res.loc[2004]
    record_caption(png, (
        "**Tests:** whether the two records agree, side by side and for the gap, once on one scale. **How to "
        "read:** (a) The gap, north minus south, from the aerial photos (squares) and from CoastSat annual "
        "means (line), both on the photos' 1967 datum (figure 2), so its size since 1967 reads directly. "
        "Dotted lines: the structure's recorded events; orange: the ""structural losses reported for 2015-2021. (b) Each photo minus the CoastSat annual mean on that datum, for "
        "each side (circles) and the gap (squares). Grey band: ± 2 times the typical scatter of one CoastSat "
        f"image day about its year's mean ({band / 2:.0f} m, the median over every year of the side means' "
        "per-date standard deviation), the difference a one-day photo would show against an annual mean if the "
        "two records agreed exactly. Inside the band is agreement; points outside are outlined. "
        f"**Shows:** before the fill the two records agree: {int(inside['north'][pre].sum())} of {pre.sum()} "
        f"north, {int(inside['south'][pre].sum())} of {pre.sum()} south and {int(inside['gap'][pre].sum())} of "
        f"{pre.sum()} gap differences sit inside the band, with RMS differences of {rms['north']:.0f}, "
        f"{rms['south']:.0f} and {rms['gap']:.0f} m (the shift makes the mean zero by construction, so only the "
        f"scatter tests agreement). The one exception is the 2004 gap, {r04['gap_photo_minus_coastsat_m']:+.0f} m "
        f"(north {r04['north_photo_minus_coastsat_m']:+.0f}, south {r04['south_photo_minus_coastsat_m']:+.0f}), "
        f"though 2004, after Isabel, is CoastSat's most scattered year, and against that year's own scatter it "
        f"is {z04:.1f} SD. After the fill both sides of the photos sit landward of CoastSat (north "
        f"{after['north']:+.0f} m, south {after['south']:+.0f} m on average; {int((~inside['north'][post]).sum())} "
        f"north and {int((~inside['south'][post]).sum())} south outside the band) while every gap difference stays "
        f"inside ({after['gap']:+.0f} m on average): a shift common to both sides, not a disagreement about the "
        "groin; a wet/dry line and a satellite waterline could part that way on a wider, re-graded fill beach, "
        "which this figure cannot confirm. Per-year spreads are in tables/coastsat_one_day_spread.csv. The annual record fills in what the photos skip: "
        f"the gap was {cs.loc[1996, 'gap_m']:.0f} m in 1996 by CoastSat ({ph.loc[1996, 'gap_m']:.0f} m in the "
        "photos), slipped slowly with no step at Isabel, where the single 2004 photo sits well above CoastSat, "
        f"and fell to {cs.loc[last, 'gap_m']:.0f} m by {last} after the fill. The shared shift moves the "
        f"CoastSat north side by {shift['north']:+.0f} m and the south side by {shift['south']:+.0f} m."))


# 4. The gap's rate in each era
def fig_eras():
    e = tab("era_rates_both")
    fig, ax = plt.subplots(figsize=figsize("single", height=3.4), constrained_layout=True)
    eras = [f"{lo}-{hi}" for lo, hi in ERAS]
    for k, src in enumerate(("domain", "photos")):
        s = e[e["source"] == src].set_index("era").reindex(eras)
        x = np.arange(len(eras)) + (k - 0.5) * 0.22
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
    ax.set_title("Do the two records give the same rates?")
    ax.legend(frameon=False, fontsize=6.8, loc="upper right")
    png = out("4_era_rates")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether the gap widened (groin trapping), held, or narrowed in each era of the "
        "structure's history, and whether CoastSat and the aerial photos agree. **How to read:** the "
        "least-squares rate of the gap in each era, ± 1 standard error; above zero the north side is "
        "gaining on the south. Eras: up to the 1995 last repair, from 1996 to Hurricane Isabel (2003), and "
        "from 2004 to the 2017 fill. CoastSat on annual means; photos on the dated wet/dry surveys inside "
        "each era (only 1996 and 1997 fall in 1996-2003, too few for a rate); CoastSat on the whole "
        "domains, the quantity the model is fitted on. **Shows:** both sources "
        "show the gap widening before the repair. From 1996 to 2003 CoastSat shows a slight narrowing "
        "that is not distinguishable from zero. After 2004 they disagree: CoastSat is flat while the "
        "photos fall, a fall that starts from the single high 2004 survey (figure 3)."))


# 5. The groin's state in 1996, the start of the hindcast
def fig_state_1996():
    a, ph, e, br = tab("gap_annual"), tab("photo_vs_coastsat"), tab("era_rates_both"), tab("gap_breakpoint")
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
    png = out("5_state_in_1996")
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
        f"state in 1996 in words, from (a), (b), `../2-coastsat_gap/` figures 4-8 and the structure record. **Shows:** in "
        f"1996 the gap stood at {g96:.0f} m (photos), {100 * g96 / photos.max():.0f}% of its {photos.max():.0f} m "
        f"peak in {peak_yr}, a year after the last repair. It had built at {r_build:+.1f} m/yr in the "
        f"first decade and was still widening at {ce.loc['1984-1995', 'rate_m_yr']:+.1f} m/yr to 1995, "
        f"but from 1996 to 2003 it changed by {ce.loc['1996-2003', 'rate_m_yr']:+.1f} ± "
        f"{ce.loc['1996-2003', 'se_m_yr']:.1f} m/yr, not distinguishable from zero. So in 1996 the groin "
        "held a full-sized gap but no longer widened it, and nothing in the record or the data marks a "
        "collapse before the 2015-2021 losses."))


# 6. The gap across both candidate model boundaries, from both records, and what a run would be scored on
def fig_placement():
    pg = tab("photo_boundary_gaps").set_index("year")
    on = tab("boundary_gaps_on_1967_datum").set_index("year")
    tg = pd.read_csv(TABLES / "placement_targets.csv", index_col=0)
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8), sharex=True, constrained_layout=True)
    for i, ((s, n), ax) in enumerate(zip(BOUNDARIES, axes)):
        k = f"gap_{s}_{n}_m"
        col = C["ACCENT"] if (s, n) == (5, 6) else INK
        ax.axvspan(KEYS[0], KEYS[1], color=C_CAL, alpha=0.08, lw=0, zorder=0)
        ax.axvspan(KEYS[1], KEYS[2], color=C_TEST, alpha=0.08, lw=0, zorder=0)
        ax.axvspan(*LOSSES, color=C["ADDED_FILL"], alpha=0.5, lw=0, zorder=0)
        ax.plot(on.index, on[f"coastsat_{k}"], color=col, lw=1.2, marker="o", ms=2.6, zorder=3)
        ax.plot(pg.index, pg[k], ls="none", marker="s", ms=4.2, color=col, mec="white", mew=0.6, zorder=4)
        b = tg[f"{s}/{n}"]
        for y in KEYS:
            v = float(b[f"coastsat_{y}_m"])
            ax.annotate(f"{v:.0f}", xy=(y, v), xytext=(0, 7), textcoords="offset points", ha="center", fontsize=6.8,
                        color=col, fontweight="bold")
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1965, 2027)
        ax.set_ylim(-15, 175)
        ax.set_ylabel(f"Gap across GIS {s}/{n} (m)\nGIS {n} minus GIS {s}")
        life_events(ax)
        chart(ax)
        question(ax, i, f"The gap across GIS {s}/{n}" + (" (where the model groin sits now)" if (s, n) == (5, 6) else ""))
    axes[-1].set_xlabel("Year")
    handles = [Line2D([], [], color=INK_MUTED, lw=1.2, marker="o", ms=2.6, label="CoastSat, annual, moved onto the photos"),
               Line2D([], [], ls="none", marker="s", ms=4.2, color=INK_MUTED, mec="white", label="aerial photos"),
               Patch(facecolor=C_CAL, alpha=0.25, label=f"calibration {KEYS[0]}-{KEYS[1]}"),
               Patch(facecolor=C_TEST, alpha=0.25, label=f"test {KEYS[1]}-{KEYS[2]}"),
               Patch(facecolor=C["ADDED_FILL"], alpha=0.5, label="structural losses reported, 2015-2021")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    png = out("6_placement")
    save(fig, png, close=True)
    a, b = tg["5/6"].astype(float), tg["6/7"].astype(float)
    record_caption(png, (
        "**Tests:** what the observed gap is across each candidate boundary for the model groin, GIS 5/6 (where it "
        "sits now) or GIS 6/7, so a run with the groin at either can be scored on its own boundary. **How to "
        "read:** each panel is the north domain's position minus the south domain's, since 1967, seaward "
        "positive: squares from the aerial photos (`../1-aerial_photo_gap/` figure 6), the line from CoastSat "
        "(`../2-coastsat_gap/` figure 9) moved onto the photos by the mean difference over their shared "
        f"{FIT[0]}-{FIT[1]} years ({a['shift_m']:+.0f} m at 5/6, {b['shift_m']:+.0f} m at 6/7; RMS difference after "
        f"the shift {a['photo_minus_coastsat_rms_m']:.0f} and {b['photo_minus_coastsat_rms_m']:.0f} m). Bold "
        f"numbers: the CoastSat gap at {KEYS[0]}, {KEYS[1]} and {KEYS[2]}, the calibration and test keys (red and "
        "blue shading); `tables/placement_targets.csv` holds them with the era rates. **Shows:** the groin's "
        f"effect is across 5/6. The 5/6 gap built to {a['photo_peak_m']:.0f} m ({int(a['photo_peak_year'])}) and "
        f"was {a[f'coastsat_{KEYS[0]}_m']:.0f} m in {KEYS[0]}, {a[f'coastsat_{KEYS[1]}_m']:.0f} m in {KEYS[1]} and "
        f"{a[f'coastsat_{KEYS[2]}_m']:.0f} m in {KEYS[2]}. The 6/7 gap peaked at {b['photo_peak_m']:.0f} m "
        f"({int(b['photo_peak_year'])}) and was {b[f'coastsat_{KEYS[0]}_m']:.0f}, {b[f'coastsat_{KEYS[1]}_m']:.0f} "
        f"and {b[f'coastsat_{KEYS[2]}_m']:.0f} m: it did not widen while the groins worked "
        f"({b['coastsat_rate_1984_1995_m_yr']:+.1f} ± {b['coastsat_rate_1984_1995_se']:.1f} m/yr over 1984-1995, "
        f"against {a['coastsat_rate_1984_1995_m_yr']:+.1f} ± {a['coastsat_rate_1984_1995_se']:.1f} at 5/6), and over "
        f"the calibration it fell {abs(b[f'coastsat_change_{KEYS[0]}_{KEYS[1]}_m']):.0f} m against "
        f"{abs(a[f'coastsat_change_{KEYS[0]}_{KEYS[1]}_m']):.0f} m at 5/6; over the test it fell "
        f"{abs(b[f'coastsat_change_{KEYS[1]}_{KEYS[2]}_m']):.0f} m against {abs(a[f'coastsat_change_{KEYS[1]}_{KEYS[2]}_m']):.0f} m. "
        "Scored on its own boundary, a groin at 6/7 has a small target that never widened while the groins "
        "worked; scored on 5/6, it has to reproduce the larger gap across a boundary it no longer sits on."))


def main():
    apply_style()
    for f in (fig_what_gap, fig_align, fig_two_sides, fig_eras, fig_state_1996, fig_placement):
        f()
        print(f"drew {f.__name__}")


if __name__ == "__main__":
    main()
