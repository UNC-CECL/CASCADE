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
    BREAK_RANGE, ERAS, FIT, GROINS, N_BOOT, N_BOOT_CHECK, NEAR_M, STEP_YEAR, TABLES)
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


def main():
    apply_style()
    for f in (fig_map, fig_sides, fig_gap, fig_break, fig_eras, fig_robust, fig_photos):
        f()
        print(f"drew {f.__name__}")


if __name__ == "__main__":
    main()
