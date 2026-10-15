#!/usr/bin/env python3
"""
The CoastSat gap figures, drawn from the tables coastsat_gap_across_groins.py writes.

    python coastsat_gap_figures.py   ->  figures/coastsat_gap_<n>_<what>.png

CoastSat only. One figure per step, in the order of the analysis: (1) the transects,
(2) one transect's images to annual means, (3) transects to one position per side,
(4) the two sides together, (5) the gap and when it stopped widening, (6) which side
moved when it did, (7) how certain that year is, (8) how the rate in each era is fitted,
(9) whether the break survives other fits, (10) each transect's rate along GIS 4-8 and the gap across both candidate
model boundaries (GIS 5/6 and 6/7). Each panel title is
the question it answers; captions (what it tests, how to read it, what it shows) go
to figures/supporting/CAPTIONS.md. The photos are drawn in ../1-aerial_photo_gap/,
the two compared in ../3-photos_vs_coastsat/, which imports the helpers here.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import json
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
    BOUNDARIES, BREAK_RANGE, ERAS, FIT, GROINS, MIN_OBS_YEAR, MIN_SHARE, N_BOOT, N_BOOT_CHECK, NEAR_M, PROFILE_DOMAINS,
    PROFILE_ERAS, STEP_YEAR, TABLES, era_rate, images)
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
    return FIGS / f"coastsat_gap_{name}.png"


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
        "transect and gives the same break year (figure 5). The fourth groin line is the landward "
        "anti-flanking extension of the northern groin."))


# 2. How one transect becomes a value per year
def fig_transect():
    per = tab("transect_annual")
    d = per[(per["gap"] == "domain") & (per["side"] == "downdrift")]
    tids = sorted(d["transect_id"].unique())
    tid = tids[len(tids) // 2]                                # the middle GIS 5 transect
    a = d[d["transect_id"] == tid].set_index("year")
    img = images(tid)
    ref = float(a["fit_mean_m"].iloc[0])
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.0), sharex=True, constrained_layout=True,
                             gridspec_kw=dict(height_ratios=(1.6, 1)))
    ax = axes[0]
    yr = img["date"].dt.year + (img["date"].dt.dayofyear - 0.5) / 365.25
    ax.plot(yr, img["chainage_m"], ls="none", marker="o", ms=1.8, color=INK_MUTED, alpha=0.45, zorder=2)
    used, short = a[a["used"]], a[~a["used"]]
    ax.plot(used.index + 0.5, used["annual_mean_m"], color=SIDE_COL["downdrift"], lw=1.3, marker="o", ms=3.2,
            zorder=4)
    ax.plot(short.index + 0.5, short["annual_mean_m"], ls="none", marker="o", ms=3.6, mfc="white",
            mec=SIDE_COL["downdrift"], zorder=4)
    ax.axhline(ref, color=INK, lw=0.8, ls="--")
    ax.text(1984.3, ref, f"its {FIT[0]}-{FIT[1]} mean, {ref:.0f} m: the zero for this transect", fontsize=6.6,
            va="bottom", color=INK)
    ax.set_xlim(1983, 2027)
    ax.set_ylabel("Position along the\ntransect (m, + seaward)")
    fill_shade(ax)
    chart(ax)
    question(ax, 0, f"One transect in GIS 5 ({tid.replace('usa_NC_0032_', '')}): images to annual means")
    ax = axes[1]
    ax.bar(a.index + 0.5, a["n_images"], width=0.8, color=[C["BASE"] if u else "0.85" for u in a["used"]])
    ax.axhline(MIN_OBS_YEAR - 0.5, color=INK_MUTED, lw=0.6, ls=":")
    ax.set_ylabel("Images\nthat year")
    ax.set_xlabel("Year")
    fill_shade(ax)
    chart(ax)
    question(ax, 1, "How many images is each annual mean?")
    handles = [Line2D([], [], ls="none", marker="o", ms=2.5, color=INK_MUTED, alpha=0.6, label="one satellite image"),
               Line2D([], [], color=SIDE_COL["downdrift"], lw=1.3, marker="o", ms=3.2,
                      label=f"calendar-year mean (at least {MIN_OBS_YEAR} images)"),
               Line2D([], [], ls="none", marker="o", ms=3.6, mfc="white", mec=SIDE_COL["downdrift"],
                      label=f"fewer than {MIN_OBS_YEAR} images: not used"), FILL_HANDLE]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("2_one_transect")
    save(fig, png, close=True)
    n = a.loc[FIT[0]:FIT[1], "n_images"]
    record_caption(png, (
        "**Tests:** nothing; this shows the first step of the analysis on one transect. **How to read:** (a) "
        f"Grey dots: the shoreline position CoastSat found on transect {tid} (the middle of the 13 GIS 5 "
        "transects) in each satellite image, as distance along the transect, seaward positive. Gold: the "
        f"calendar-year mean, kept where at least {MIN_OBS_YEAR} images fall in the year (hollow: fewer, not "
        f"used). Dashed: the transect's own mean over {FIT[0]}-{FIT[1]}; every later figure uses each annual mean "
        "minus this value (an anomaly), so transects with different starting positions can be averaged. (b) "
        "The number of images behind each annual mean. Grey shading: after the 2017 fill, left out of every "
        f"fit. **Shows:** a single image scatters tens of metres about its year's mean; the yearly mean "
        f"averages that out, from {int(n.min())} to {int(n.max())} images a year over {FIT[0]}-{FIT[1]}, fewest "
        "before 1999."))


# 3. How the transects become one position per side
def fig_to_sides():
    per = tab("transect_annual")
    a = tab("gap_annual")
    a = a[a["gap"] == "domain"].set_index("year")
    d = per[(per["gap"] == "domain") & per["used"]]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True, constrained_layout=True)
    for i, (ax, side, gis) in enumerate(zip(axes, ("updrift", "downdrift"), (6, 5))):
        s = d[d["side"] == side]
        for _, g in s.groupby("transect_id"):
            ax.plot(g["year"], g["anom_m"], color=SIDE_COL[side], lw=0.5, alpha=0.35, zorder=2)
        ax.plot(a.index, a[f"{side}_anom_m"], color=SIDE_COL[side], lw=1.8, marker="o", ms=2.8, zorder=4)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1983, 2027)
        fill_shade(ax)
        events(ax)
        chart(ax)
        ax.set_ylabel(f"Position (m), + seaward\nabout {FIT[0]}-{FIT[1]} mean")
        question(ax, i, f"GIS {gis}, the {SIDE[side]}: {s['transect_id'].nunique()} transects and their mean")
    axes[-1].set_xlabel("Year")
    handles = [Line2D([], [], color=INK_MUTED, lw=0.5, alpha=0.6, label="one transect's annual anomaly"),
               Line2D([], [], color=INK_MUTED, lw=1.8, marker="o", ms=2.8,
                      label=f"side mean (years with at least {MIN_SHARE:.0%} of the transects)"), FILL_HANDLE]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    png = out("3_transects_to_sides")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** nothing; this shows the second step, from transects to one position per side. **How to "
        "read:** thin lines: each transect's annual anomaly (figure 2) in GIS 6 (a, north, updrift) and GIS 5 "
        "(b, south, downdrift); thick line: the side's position, the mean of its transects in each year where "
        f"at least {MIN_SHARE:.0%} of them have an annual mean. Dotted lines: Hurricane Gordon (1994), the last "
        "repair (1995), Hurricane Isabel (2003), the Buxton fill (2017); grey: after the fill, not fitted. "
        "**Shows:** before the fill the transects on a side move together, so the mean is a fair summary of the "
        "side. After it one GIS 6 line parts from the rest and rises with GIS 5: transect 0077, the "
        "southernmost GIS 6 transect, which lies just south of the southern groin (figure 1) and so takes the "
        "fill sand moving south past the field. The next figure puts the two side means together."))



# 4. Which side moved
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
    png = out("4_side_positions")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether a change in the gap comes from the north side holding (trapping) or from the "
        "south side moving. **How to read:** CoastSat annual shoreline position of each side, seaward "
        f"positive, each transect taken about its own {FIT[0]}-{FIT[1]} mean and averaged per side; the "
        "gap in figure 5 is blue minus gold. Dotted lines: Hurricane Gordon (1994), the last repair "
        "(1995), Hurricane Isabel (2003), the Buxton fill (2017); grey: years after the fill, left out of "
        f"every fit. (a) Whole domains, GIS 6 and GIS 5. (b) Transects within {NEAR_M:.0f} m of the "
        "groin field. **Shows:** the 1991-1995 widening of the gap is mostly the south side retreating "
        "(about 98 m in GIS 5 from 1984 to 1995, against about 24 m on the north side), with Gordon in "
        "that stretch, so the 1995 peak is partly storm loss downdrift rather than trapping. After 2017 "
        "the south side gains as the fill sand moves south past the groins."))


# 5. When the gap stopped widening
def fig_gap():
    a, h, br = tab("gap_annual"), tab("hinge_fit"), tab("gap_breakpoint")
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.0), sharex=True,
                             constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        s, f = a[a["gap"] == name], h[h["gap"] == name]
        b = br.set_index("gap").loc[name]
        ax.axvspan(b.break_ci90_lo, b.break_ci90_hi, color=COL[name], alpha=0.10, lw=0, zorder=0)
        ax.plot(f["year"], f["fit_m"], color=COL[name], lw=3.0, alpha=0.35, solid_capstyle="round")
        ax.plot(s["year"], s["gap_m"], color=COL[name], lw=1.2, marker="o", ms=3)
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
               FILL_HANDLE]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = out("5_gap_and_break")
    save(fig, png, close=True)
    d, n = br.set_index("gap").loc["domain"], br.set_index("gap").loc["near_field"]
    record_caption(png, (
        "**Tests:** in which year the groin stopped holding the north side ahead of the south, i.e. when "
        f"the gap stopped widening. **How to read:** {WHAT} Pale thick line: a continuous line with one "
        f"change of slope (a hinge), fitted on {FIT[0]}-{FIT[1]}; dashed line: the year of that change, "
        f"with its 90% range from {N_BOOT} residual bootstraps shaded. Grey: after the 2017 fill, not "
        f"fitted. (a) Whole domains: {d.rate_before_m_yr:+.1f} m/yr to {d.break_year}, then "
        f"{d.rate_after_m_yr:+.1f} m/yr. (b) Within {NEAR_M:.0f} m of the field: {n.rate_before_m_yr:+.1f} then "
        f"{n.rate_after_m_yr:+.1f} m/yr, break {n.break_year}. **Shows:** both gaps stop widening in "
        "1995, the year of the last repair, not at Hurricane Isabel (2003); what follows is a slow "
        "decline, not a collapse. The photo gap on the same axes is in "
        "`../3-photos_vs_coastsat/`."))


# 6. What the 1995 break is made of: which side moved, against the coast either side
def fig_break_sides():
    pos = tab("domain_positions").set_index("year")
    rates = tab("domain_era_rates")
    dcol = {4: SIDE_COL["downdrift"], 5: SIDE_COL["downdrift"], 6: SIDE_COL["updrift"], 7: "0.35", 8: "0.62"}
    dlab = {4: "GIS 4", 5: "GIS 5, south side (downdrift)", 6: "GIS 6, north side (updrift)", 7: "GIS 7", 8: "GIS 8"}
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.2), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=(1.3, 1)))
    ax = axes[0]
    for lo, hi, name in ((1991, 1995.5, "Gordon loss"), (1995.5, 2003.5, "south side recovers")):
        ax.axvspan(lo, hi, color=C["BASE_FILL"] if name == "Gordon loss" else C["ADDED_FILL"], alpha=0.35, lw=0,
                   zorder=0)
        ax.text((lo + hi) / 2, 0.03, name, transform=ax.get_xaxis_transform(), ha="center", fontsize=6.6,
                color=INK_MUTED)
    for d in (8, 7, 4, 5, 6):
        main = d in (5, 6)
        ax.plot(pos.index, pos[f"gis{d}_m"], color=dcol[d], lw=1.6 if main else 0.8, marker="o" if main else None,
                ms=2.8, ls="-" if main or d == 7 else (0, (3, 2)) if d == 4 else (0, (1, 1.2)), alpha=1 if main else 0.85, label=dlab[d],
                zorder=4 if main else 3)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1983, 2027)
    fill_shade(ax)
    events(ax)
    chart(ax)
    ax.set_ylabel(f"Position (m), + seaward\nabout {FIT[0]}-{FIT[1]} mean")
    question(ax, 0, "Which side moved when the gap stopped widening?")
    ax.legend(frameon=False, fontsize=6.6, loc="lower left", ncol=2)

    ax = axes[1]
    eras = list(dict.fromkeys(rates["era"]))
    doms = (4, 5, 6, 7, 8)
    w = 0.15
    for k, d in enumerate(doms):
        s = rates[rates["domain"] == d].set_index("era").reindex(eras)
        x = np.arange(len(eras)) + (k - 2) * w
        ax.bar(x, s["rate_m_yr"], width=w * 0.9, yerr=s["se_m_yr"], color=dcol[d],
               alpha=1 if d in (5, 6) else 0.5, capsize=1.5, error_kw=dict(lw=0.6), label=f"GIS {d}")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    names = ("before", "Gordon,\nto the repair", "after the repair,\nto Isabel", "after Isabel,\nto the fill")
    ax.set_xticks(range(len(eras)), [f"{n}\n{e}" for n, e in zip(names, eras)], fontsize=6.6)
    ax.set_ylabel("Rate (m/yr) ± 1 SE\n+ = building out")
    chart(ax)
    question(ax, 1, "How fast did each domain move in each era?")
    ax.legend(frameon=False, fontsize=6.4, ncol=5, loc="upper left")
    png = out("6_break_sides")
    save(fig, png, close=True)
    r = rates.set_index(["domain", "era"])["rate_m_yr"]
    se = rates.set_index(["domain", "era"])["se_m_yr"]
    g = tab("gap_annual")
    gd = g[g["gap"] == "domain"].set_index("year")["gap_m"]
    record_caption(png, (
        "**Tests:** whether the gap stopped widening in 1995 because the groin stopped holding the north side, or "
        "because the south side moved. **How to read:** (a) Each domain's CoastSat position, the mean of its "
        f"transects about their own {FIT[0]}-{FIT[1]} means (figure 3): GIS 6, the north (updrift) side, in blue; "
        "GIS 5, the south (downdrift) side, in gold; GIS 4 (dashed gold) further south and GIS 7-8 (grey) further "
        "north (GIS 7 solid, GIS 8 dotted), for comparison. Shading: 1991-1995, the years of Hurricane Gordon's (1994) loss on the south side; "
        "1996-2003, from the last repair (1995) to Hurricane Isabel (2003). (b) Each domain's least-squares rate in "
        "each era, ± 1 standard error. **Shows:** the 1995 break is the south side, not the north. Through "
        f"1991-1995 GIS 5 retreated {r[5, '1991-1995']:+.1f} ± {se[5, '1991-1995']:.1f} m/yr (GIS 4 "
        f"{r[4, '1991-1995']:+.1f}) while GIS 6 held ({r[6, '1991-1995']:+.1f} ± {se[6, '1991-1995']:.1f}), as GIS 7-8 "
        f"did ({r[7, '1991-1995']:+.1f}, {r[8, '1991-1995']:+.1f}): the gap "
        f"widened to {gd[1995]:.0f} m above its 1984-2016 mean in 1995 by storm loss downdrift. From 1996 to 2003 "
        f"the south side recovered ({r[5, '1996-2003']:+.1f} ± {se[5, '1996-2003']:.1f}, GIS 4 "
        f"{r[4, '1996-2003']:+.1f}) while GIS 6 still held ({r[6, '1996-2003']:+.1f} ± {se[6, '1996-2003']:.1f}): "
        "the gap narrowed because the south side came back, not because the groin let the north side go. From "
        f"2004 both sides retreated (GIS 6 {r[6, '2004-2016']:+.1f}, GIS 5 {r[5, '2004-2016']:+.1f} m/yr), as the "
        f"coast north of the field did (GIS 7 {r[7, '2004-2016']:+.1f}, GIS 8 {r[8, '2004-2016']:+.1f}), so the gap "
        "held; GIS 6 still lost ground a little more slowly than its northern neighbours. Throughout, GIS 6 moves "
        "with GIS 7-8: in CoastSat the north side behaves like the coast north of it, and the groin's mark on the "
        "gap is in how the south side moves (starved, cut by Gordon, recovering). The 1995 annual means mix months "
        "before and after the repair, whose month is not recorded."))


# 7. How certain the break year is
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
    png = out("7_break_year")
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


# 8. How the gap's rate in each era is fitted
def fig_eras():
    e = tab("era_rates")
    g = tab("gap_annual")
    eras = [f"{lo}-{hi}" for lo, hi in ERAS]
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.3), constrained_layout=True,
                             gridspec_kw=dict(width_ratios=(1.5, 1)))
    ax = axes[0]
    s = g[g["gap"] == "domain"].set_index("year")["gap_m"]
    ax.plot(s.index, s.values, ls="none", marker="o", ms=3, color=COL["domain"], zorder=4)
    r = e[e["source"] == "domain"].set_index("era")
    for (lo, hi), era in zip(ERAS, eras):
        w = s.loc[lo:hi]
        t = np.array([lo, hi], float)
        ax.plot(t, w.mean() + r.loc[era, "rate_m_yr"] * (t - w.index.to_numpy(float).mean()), color=INK, lw=1.3)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1983, 2027)
    fill_shade(ax)
    events(ax, years=(1995, 2003, 2017))
    chart(ax)
    ax.set_ylabel(GAP_LABEL)
    ax.set_xlabel("Year")
    question(ax, 0, "A straight line through each era (whole domains)")
    ax = axes[1]
    for k, src in enumerate(("domain", "near_field")):
        q = e[e["source"] == src].set_index("era").reindex(eras)
        x = np.arange(len(eras)) + (k - 0.5) * 0.22
        ax.errorbar(x, q["rate_m_yr"], yerr=q["se_m_yr"], fmt="o", color=COL[src], ms=4.5, capsize=2.5, lw=1.0,
                    label=f"CoastSat, {NAME[src]}")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(range(len(eras)), [f"to the last\nrepair\n{eras[0]}", f"repair to\nIsabel\n{eras[1]}",
                                     f"Isabel to\nthe fill\n{eras[2]}"], fontsize=6.6)
    ax.set_ylabel("Rate of gap change (m/yr)\n± 1 SE; + = widening")
    chart(ax)
    question(ax, 1, "Each line's slope")
    ax.legend(frameon=False, fontsize=6.3, loc="lower left")
    png = out("8_era_rates")
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** whether the gap widened (groin trapping), held, or narrowed in each era of the "
        "structure's history. **How to read:** (a) The whole-domain annual gap (figure 5) with a least-squares "
        "line through the years of each era: up to the 1995 last repair, from 1996 to Hurricane Isabel (2003), "
        "and from 2004 to the 2017 fill. (b) Each line's slope ± 1 standard error, for both gap definitions. "
        "**Shows:** both gaps widened before the repair. From 1996 to 2003 they narrowed slightly, not "
        "distinguishably from zero, and from 2004 to the fill they were flat. The photo rates on the same eras "
        "are compared in `../3-photos_vs_coastsat/`."))


# 9. Whether the break survives other fits
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
    png = out("9_robustness")
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


# 10. Where along the coast the shoreline retreated fastest in each era, and the gap across both boundaries
def fig_placement():
    r = tab("alongshore_rates")
    b = tab("boundary_gaps").set_index("year")
    g = json.load(open(GROINS))
    south = min(c[1] for f in g["features"] for c in f["geometry"]["coordinates"])
    groin_y = sorted({round(min(c[1] for c in f["geometry"]["coordinates"]) - south) for f in g["features"]})
    dom = {f["properties"]["domain_id"]: np.asarray(f["geometry"]["coordinates"][0]) for f in json.load(open(DOMAINS))["features"]}
    bounds = {f"{s}/{n}": dom[n][:, 1].min() - south for s, n in BOUNDARIES}
    mids = {d: dom[d][:, 1].mean() - south for d in PROFILE_DOMAINS}
    era_col = dict(zip([f"{lo}-{hi}" for lo, hi in PROFILE_ERAS], (COL["domain"], C["ACCENT"], INK, INK_MUTED)))
    era_name = dict(zip(era_col, ("to the last repair", "repair to Isabel", "Isabel to the fill", "after the fill")))
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=6.2), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=(1.4, 1)))
    ax = axes[0]
    for era, col in era_col.items():
        s = r[r["era"] == era].sort_values("from_south_groin_m")
        ax.fill_between(s["from_south_groin_m"], s["rate_m_yr"] - s["se_m_yr"], s["rate_m_yr"] + s["se_m_yr"],
                        color=col, alpha=0.12, lw=0)
        ax.plot(s["from_south_groin_m"], s["rate_m_yr"], color=col, lw=1.3, marker="o", ms=2.6,
                ls="--" if era.startswith("2017") else "-", label=f"{era}, {era_name[era]}")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    for name, x in bounds.items():
        ax.axvline(x, color=INK, lw=0.9, ls="--", zorder=1)
        ax.text(x, 0.985, f" GIS {name}", transform=ax.get_xaxis_transform(), fontsize=6.8, va="top")
    for x in groin_y:
        ax.plot(x, 0, transform=ax.get_xaxis_transform(), marker="^", ms=6, color=C["GROIN"], clip_on=False, zorder=6)
    for d, x in mids.items():
        ax.text(x, 0.03, f"GIS {d}", transform=ax.get_xaxis_transform(), fontsize=6.5, color=INK_MUTED, ha="center")
    ax.set_xlabel("Northing of the transect from the southern groin (m), south → north")
    ax.set_ylabel("Shoreline rate (m/yr)\n+ = building out")
    ax.legend(frameon=False, fontsize=6.6, loc="lower right", bbox_to_anchor=(1.0, 0.07))
    chart(ax)
    question(ax, 0, "How fast did each transect move, era by era?")

    ax = axes[1]
    for (s, n), col in zip(BOUNDARIES, (C["ACCENT"], INK)):
        ax.plot(b.index, b[f"gap_{s}_{n}_m"], color=col, lw=1.2, marker="o", ms=2.8,
                label=f"across GIS {s}/{n}: GIS {n} minus GIS {s}")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1983, 2027)
    fill_shade(ax)
    events(ax)
    chart(ax)
    ax.set_ylabel(f"Gap (m), each side\nabout its {FIT[0]}-{FIT[1]} mean")
    ax.set_xlabel("Year")
    question(ax, 1, "How did the gap across each candidate boundary change?")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = out("10_placement")
    save(fig, png, close=True)
    w = r[r["era"] == f"{ERAS[0][0]}-{ERAS[0][1]}"].set_index("transect_id")
    south_of = w[w["from_south_groin_m"] < 0]
    north_of = w[w["from_south_groin_m"] > 0]
    by_dom = north_of.groupby("domain")["rate_m_yr"].mean()
    e = {f"{s}/{n}": era_rate(b[f"gap_{s}_{n}_m"], *ERAS[0]) for s, n in BOUNDARIES}
    record_caption(png, (
        "**Tests:** whether the CoastSat record puts the groin's effect at the GIS 5/6 boundary, where the model "
        "groin sits, or at 6/7. **How to read:** (a) Each point is one CoastSat transect in GIS 4-8, placed by the "
        "northing of its landward end from the southern groin; y is its least-squares rate over the annual means "
        "of each era (shading ± 1 standard error), seaward positive. Dashed verticals: the two candidate "
        "boundaries; red triangles: the groins. (b) The gap across each boundary: the north domain's mean "
        "position minus the south domain's, each side about its own 1984-2016 mean, so only its changes read "
        "(the size since 1967 is in `../3-photos_vs_coastsat/`). Grey: after the 2017 fill. **Shows:** over "
        f"{ERAS[0][0]}-{ERAS[0][1]}, while the groins worked, the rate steps at the southern groin: "
        f"{south_of['rate_m_yr'].min():.1f} to {south_of['rate_m_yr'].max():.1f} m/yr south of it (GIS 4-5 and "
        f"transect 0077) against {north_of['rate_m_yr'].min():.1f} to {north_of['rate_m_yr'].max():.1f} m/yr north "
        f"of it, with no second step at 6/7 (GIS 6 {by_dom.get(6, np.nan):.1f}, GIS 7 {by_dom.get(7, np.nan):.1f}, "
        f"GIS 8 {by_dom.get(8, np.nan):.1f} m/yr on average). The gap across 5/6 widened at "
        f"{e['5/6'][0]:+.1f} ± {e['5/6'][1]:.1f} m/yr over those years and across 6/7 at {e['6/7'][0]:+.1f} ± "
        f"{e['6/7'][1]:.1f} m/yr. In each later era the step at the southern groin reverses, GIS 4-5 building out "
        "or retreating more slowly than the field (the 5/6 gap narrowing), and 6/7 still shows no step. CoastSat starts in 1984, after the fillet had built (photos, "
        "`../1-aerial_photo_gap/` figure 6), so these are the rates of a working groin, not the build-up."))


def main():
    apply_style()
    for f in (fig_map, fig_transect, fig_to_sides, fig_sides, fig_gap, fig_break_sides, fig_break, fig_eras, fig_robust,
              fig_placement):
        f()
        print(f"drew {f.__name__}")


if __name__ == "__main__":
    main()
