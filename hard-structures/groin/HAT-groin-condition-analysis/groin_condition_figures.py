#!/usr/bin/env python3
"""
The groin-condition figures, drawn from the tables coastsat_groin_gap.py writes.

    python groin_condition_figures.py   ->  figures/groin_condition_<n>_<what>.png

One figure per question, numbered in reading order: where the gap is measured,
how each side moved, the gap and its break, how sure the break year is, the
rate in each era, whether the break survives the checks, where the photos and
CoastSat disagree, and the failure schedules that follow. Captions go to
figures/supporting/CAPTIONS.md.

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
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(HERE))
from coastsat_groin_gap import (  # noqa: E402
    BREAK_RANGE, ERAS, FIT, GROINS, N_BOOT, N_BOOT_CHECK, NEAR_M, STEP_YEAR, TABLES)
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, INK, INK_MUTED, _title, apply_style, figsize, map_label, north_dart,
    open_frame, record_caption, save, scale_bar_km, spines_for_image)

# --- CONFIG ------------------------------------------------------------------
FIGS = HERE / "figures"
DOMAINS = REPO / "data" / "hatteras_init" / "5-scr" / "2-transect-frame" / "transect_domains" / "HAT_domains.json"
MEAN_LINE = (REPO / "data" / "hatteras_init" / "5-scr" / "1-observations" / "mean_shoreline"
             / "1995-10-12_1997-10-12" / "transect_means_1995-10-12_1997-10-12.csv")
COL = {"domain": C_1997, "near_field": C["ACCENT"], "photos": C_1984}
SIDE_COL = {"updrift": C_1997, "downdrift": C["ADDED"]}
NAME = {"domain": "GIS 6 minus GIS 5", "near_field": f"{NEAR_M:.0f} m either side of the field",
        "photos": "wet/dry photos"}
EVENTS = {1994: "Gordon", 1995: "last repair", 2003: "Isabel", 2017: "Buxton fill"}
F_EXAMPLE = 0.6                       # the pinned post-failure fraction, for the schedule figure
MAP_PAD_M = 250.0
# -----------------------------------------------------------------------------

WHAT = ("Gap = updrift (north) minus downdrift (south) shoreline position, seaward positive, so a "
        "rising gap means the north side holding while the south retreats. CoastSat: calendar-year "
        "means per transect (>= 3 images), each about its own "
        f"{FIT[0]}-{FIT[1]} mean, averaged per side in years with >= half the side's transects.")


def tab(name):
    return pd.read_csv(TABLES / f"{name}.csv")


# Dotted event lines, named once along the top
def events(ax, label=True, years=None):
    lo, hi = ax.get_xlim()
    for yr, lab in EVENTS.items():
        if (years and yr not in years) or not lo <= yr <= hi:
            continue
        ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":", zorder=1)
        if label:
            ax.text(yr + (-0.25 if yr == 1994 else 0.25), 1.0, lab, transform=ax.get_xaxis_transform(),
                    fontsize=6.5, color=INK_MUTED, va="top", ha="right" if yr == 1994 else "left")


# The years after the 2017 fill, left out of every fit
def fill_shade(ax):
    ax.axvspan(FIT[1] + 0.5, ax.get_xlim()[1], color=C["BASE_FILL"], alpha=0.45, lw=0, zorder=0)


def chart(ax):
    ax.grid(axis="y")
    open_frame(ax)


def out(name):
    return FIGS / f"groin_condition_{name}.png"


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
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=5.3), constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        for f in dom["features"]:
            d = f["properties"]["domain_id"]
            if d not in (4, 5, 6, 7):
                continue
            ring = np.asarray(f["geometry"]["coordinates"][0])
            ax.fill(ring[:, 0], ring[:, 1], facecolor="0.95" if d % 2 else "0.90", edgecolor="0.7",
                    lw=0.5, zorder=0)
            cx, cy = ring[:, 0].mean(), np.clip(ring[:, 1].mean(), *ylim)
            map_label(ax, xlim[0] + 0.12 * (xlim[1] - xlim[0]), cy, f"GIS {d}", color=INK,
                      fontsize=7.5, path_effects=[])
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
        scale_bar_km(ax, length_m=500, segments=2, unit="m", x=0.08, y=0.04)
        north_dart(ax, (xlim[1] - 0.14 * (xlim[1] - xlim[0]), ylim[1] - 0.07 * (ylim[1] - ylim[0])),
                   arrow_m=160)
        _title(ax, i, NAME[name])
    handles = [Line2D([], [], color=SIDE_COL["updrift"], lw=1.6, label="updrift transects (north)"),
               Line2D([], [], color=SIDE_COL["downdrift"], lw=1.6, label="downdrift transects (south)"),
               Line2D([], [], color=C["GROIN"], lw=2.4, label="Buxton groins"),
               Line2D([], [], color=INK, lw=1.0, label="CoastSat mean shoreline, Oct 1995-Oct 1997")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False, fontsize=7.5)
    n = used.groupby(["gap", "side"]).size()
    png = out("1_transect_map")
    save(fig, png, close=True)
    record_caption(png, (
        "The CoastSat transects behind each gap, UTM 18N, north up, ocean to the right. "
        f"(a) The fit quantity: every transect in GIS 6 (updrift, {n['domain', 'updrift']}) and GIS 5 "
        f"(downdrift, {n['domain', 'downdrift']}), the domains the model's groin sits between. (b) The "
        "condition check: transects whose landward end lies within "
        f"{NEAR_M:.0f} m north of the northernmost groin ({n['near_field', 'updrift']}) or south of the "
        f"southernmost ({n['near_field', 'downdrift']}). "
        "Grey bands: the GIS domains. Black line: the CoastSat mean shoreline over the 1996 DEM "
        "window. The fourth groin line is the landward anti-flanking extension of the northern groin. "
        "The southernmost GIS 6 transect lies just south of the southern groin, so the domain gap is "
        "not purely across the field; the near-field gap in (b) has no such transect and gives the "
        "same break year."))


# 2. How each side moved
def fig_sides():
    a = tab("gap_annual")
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True,
                             constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        s = a[a["gap"] == name]
        for side in ("updrift", "downdrift"):
            ax.plot(s["year"], s[f"{side}_anom_m"], color=SIDE_COL[side], lw=1.2, marker="o", ms=2.8)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1983, 2027)
        fill_shade(ax)
        events(ax, label=(i == 0))
        chart(ax)
        ax.set_ylabel("Position about the\n1984-2016 mean (m), seaward ▲")
        _title(ax, i, NAME[name])
    handles = [Line2D([], [], color=SIDE_COL[s], lw=1.2, marker="o", ms=2.8, label=f"{s} side")
               for s in ("updrift", "downdrift")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    png = out("2_side_positions")
    save(fig, png, close=True)
    record_caption(png, (
        "Each side of the groin on its own, CoastSat annual means. Each transect is taken about its "
        f"{FIT[0]}-{FIT[1]} mean, then averaged per side. The 1991-1995 widening of the gap is mostly "
        "the downdrift (south) side retreating, by about 98 m in GIS 5 from 1984 to 1995, against "
        "about 24 m updrift; Hurricane Gordon (1994) falls in that stretch. After 2017 the downdrift "
        "side gains: the Buxton fill sand moving south past the groin. Grey: years after the 2017 "
        "fill, left out of every fit. (a) GIS 6 and GIS 5. (b) Transects within "
        f"{NEAR_M:.0f} m of the groin field."))


# 3. The gap, its hinge, and the photos
def fig_gap():
    a, h, ph, br = tab("gap_annual"), tab("hinge_fit"), tab("photo_gap"), tab("gap_breakpoint")
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8), sharex=True,
                             constrained_layout=True)
    for i, (ax, name) in enumerate(zip(axes, ("domain", "near_field"))):
        s, f = a[a["gap"] == name], h[h["gap"] == name]
        b = br.set_index("gap").loc[name]
        ax.plot(f["year"], f["fit_m"], color=COL[name], lw=3.0, alpha=0.35, solid_capstyle="round")
        ax.plot(s["year"], s["gap_m"], color=COL[name], lw=1.2, marker="o", ms=3)
        if name == "domain":
            p = ph[ph["year"] >= 1984]
            ax.plot(p["year"], p["shifted_m"], ls="none", marker="s", ms=4, color=COL["photos"],
                    mec="white", mew=0.5, zorder=5)
        ax.axvline(b.break_year, color=COL[name], lw=0.9, ls="--")
        ax.text(b.break_year + 0.3, 0.06, f"break {b.break_year}\n(90%: {b.break_ci90_lo}-{b.break_ci90_hi})",
                transform=ax.get_xaxis_transform(), fontsize=7, color=COL[name])
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1983, 2027)
        fill_shade(ax)
        events(ax, label=(i == 0), years=(1995, 2003, 2017))
        chart(ax)
        ax.set_ylabel("Gap across the groin (m)\nupdrift seaward ▲")
        _title(ax, i, NAME[name])
    handles = [Line2D([], [], color=COL["domain"], lw=1.2, marker="o", ms=3, label="CoastSat gap, annual"),
               Line2D([], [], color=INK_MUTED, lw=3.0, alpha=0.5, label="one-break hinge, 1984-2016"),
               Line2D([], [], ls="none", marker="s", ms=4, color=COL["photos"], mec="white",
                      label="wet/dry photo gap (shifted)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = out("3_gap_and_break")
    save(fig, png, close=True)
    d, n = br.set_index("gap").loc["domain"], br.set_index("gap").loc["near_field"]
    record_caption(png, (
        f"The gap across the Buxton groin, 1984-2026. {WHAT} Pale thick lines: a continuous one-break "
        f"hinge fitted on {FIT[0]}-{FIT[1]}; dashed: its break year, with the 90% range from "
        f"{N_BOOT} residual bootstraps. (a) GIS 6 minus GIS 5: rises {d.rate_before_m_yr:+.1f} m/yr to "
        f"{d.break_year}, then {d.rate_after_m_yr:+.1f} m/yr. Red squares: the wet/dry photo gap "
        "(GIS 5 minus GIS 6 change since 1967, landward positive, the same sense), shifted to the "
        "CoastSat series' mean over the shared pre-fill years, so only its shape compares. (b) "
        f"Transects within {NEAR_M:.0f} m of the field: {n.rate_before_m_yr:+.1f} then "
        f"{n.rate_after_m_yr:+.1f} m/yr, break {n.break_year}. Grey: after the 2017 fill, not fitted."))


# 4. How sure the break year is
def fig_break():
    prof, boots = tab("break_profile"), tab("break_bootstrap")
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=3.0), constrained_layout=True)
    ax = axes[0]
    for name in ("domain", "near_field"):
        p = prof[prof["gap"] == name]
        ax.plot(p["break_year"], p["sse_over_best"], color=COL[name], lw=1.4)
    ax.set_ylabel("Misfit / best misfit")
    ax.set_xlabel("Break year")
    events(ax, years=(1995, 2003))
    chart(ax)
    _title(ax, 0, "Misfit by break year")
    ax = axes[1]
    years = np.arange(BREAK_RANGE[0], BREAK_RANGE[1] + 1)
    for k, name in enumerate(("domain", "near_field")):
        share = boots[boots["gap"] == name]["break_year"].value_counts(normalize=True)
        ax.bar(years + (k - 0.5) * 0.4, share.reindex(years, fill_value=0), width=0.4, color=COL[name])
    ax.axvspan(2001.5, 2005.5, color=C["BASE_FILL"], alpha=0.6, lw=0, zorder=0)
    ax.text(2003.5, 0.97, "2002-05", transform=ax.get_xaxis_transform(), ha="center", va="top",
            fontsize=7, color=INK_MUTED)
    ax.set_ylabel("Share of bootstraps")
    ax.set_xlabel("Break year")
    chart(ax)
    _title(ax, 1, "Break year over bootstraps")
    handles = [Patch(color=COL[n], label=f"CoastSat, {NAME[n]}") for n in ("domain", "near_field")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    png = out("4_break_year")
    save(fig, png, close=True)
    br = tab("gap_breakpoint").set_index("gap")
    record_caption(png, (
        f"How well the break year is pinned. (a) Misfit (sum of squares) of a one-break hinge on "
        f"the {FIT[0]}-{FIT[1]} annual gap, for each candidate break year, over the best; the minimum "
        f"is 1995 in both. (b) The best break year in {N_BOOT} residual bootstraps. Grey band: "
        "2002-2005, the timing of the model's failure step, which holds "
        f"{100 * br.loc['domain', 'boot_share_2002_2005']:.0f}% (GIS 6 minus GIS 5) and "
        f"{100 * br.loc['near_field', 'boot_share_2002_2005']:.0f}% (near field) of the bootstraps; "
        f"1997 or earlier holds {100 * br.loc['domain', 'boot_share_le_1997']:.0f}% and "
        f"{100 * br.loc['near_field', 'boot_share_le_1997']:.0f}%."))


# 5. The gap's rate in each era
def fig_eras():
    e = tab("era_rates")
    fig, ax = plt.subplots(figsize=figsize("single", height=3.1), constrained_layout=True)
    eras = [f"{lo}-{hi}" for lo, hi in ERAS]
    for k, src in enumerate(("domain", "near_field", "photos")):
        s = e[e["source"] == src].set_index("era").reindex(eras)
        x = np.arange(len(eras)) + (k - 1) * 0.22
        ok = s["rate_m_yr"].notna().to_numpy()
        ax.errorbar(x[ok], s["rate_m_yr"][ok], yerr=s["se_m_yr"][ok], fmt="s" if src == "photos" else "o",
                    color=COL[src], ms=4.5, capsize=2.5, lw=1.0, label=NAME[src])
        for xi, (_, r) in zip(x, s.iterrows()):
            if not np.isfinite(r.rate_m_yr):
                ax.text(xi, 0.4, f"{int(r.n_years)}\nsurveys", ha="center", va="bottom", fontsize=6.5,
                        color=COL[src])
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(range(len(eras)), ["to the last\nrepair", "after it,\nto Isabel", "after Isabel,\nto the fill"])
    for i, lab in enumerate(eras):
        ax.text(i, -0.2, lab, transform=ax.get_xaxis_transform(), ha="center", fontsize=7, color=INK_MUTED)
    ax.set_ylabel("Gap rate (m/yr), widening ▲")
    chart(ax)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    png = out("5_era_rates")
    save(fig, png, close=True)
    record_caption(png, (
        "The gap's least-squares rate in three eras, ± 1 standard error: to the 1995 last repair, "
        "from 1996 to Isabel (2003), and from 2004 to the 2017 fill. CoastSat on annual means; photos "
        "on the dated wet/dry surveys inside each era (only 1996 and 1997 fall in 1996-2003, too few "
        "for a rate). Both sources show the widening before the repair. After 2004 they disagree: "
        "CoastSat is flat, the photos fall, a fall that starts from the single high 2004 survey."))


# 6. Does the break survive the checks
def fig_robust():
    r = tab("gap_robustness")
    checks = list(dict.fromkeys(r["check"]))
    fig, axes = plt.subplots(1, 2, figsize=figsize("double", height=2.8), sharey=True,
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
    axes[0].set_xlim(1988, 2008)
    axes[0].set_xticks(range(1990, 2009, 5))
    axes[0].set_xlabel("Break year (90% bootstrap range)")
    axes[1].axvline(0, color=INK_MUTED, lw=0.6)
    axes[1].set_xlabel(f"Level step at {STEP_YEAR} (m, ± 1 SE)")
    axes[0].set_yticks(range(len(checks)), checks)
    axes[0].invert_yaxis()
    for i, (ax, t) in enumerate(zip(axes, ("Break year", f"A separate drop at {STEP_YEAR}?"))):
        ax.grid(axis="x")
        open_frame(ax)
        _title(ax, i, t)
    handles = [Line2D([], [], color=COL[n], marker="o", ms=4, ls="none", label=f"CoastSat, {NAME[n]}")
               for n in ("domain", "near_field")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    png = out("6_robustness")
    save(fig, png, close=True)
    record_caption(png, (
        "Whether the 1995 break survives changes to the fit. Rows: the fit as made (1984-2016); "
        "starting in 1988, which drops the sparse early Landsat years (5-10 images a year); leaving "
        "out the 1995 peak; ending in 2013. (a) Best break year with its 90% range from "
        f"{N_BOOT_CHECK} residual bootstraps; grey: 2002-2005. (b) The size of a level step at "
        f"{STEP_YEAR}, the model's failure step, fitted beside the free break; none differs from zero "
        "by two standard errors, and adding it lowers the BIC below the break-only fit in one of the "
        "eight fits only (GIS 6 minus GIS 5 with 1995 left out), by 0.1 (tables/gap_robustness.csv)."))


# 7. Where the photos and CoastSat disagree
def fig_photos():
    a, ph = tab("gap_annual"), tab("photo_gap")
    d = a[a["gap"] == "domain"]
    p = ph[ph["year"] >= 1984]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.0), sharex=True,
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
    ax.set_ylabel("Gap (m)\nupdrift seaward ▲")
    _title(ax, 0, "CoastSat and the photos, GIS 6 minus GIS 5")
    ax = axes[1]
    q = p.dropna(subset=["photo_minus_coastsat_m"])
    cols = [COL["photos"] if y == 2004 else C["BASE"] for y in q["year"]]
    ax.bar(q["year"], q["photo_minus_coastsat_m"], width=0.7, color=cols)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    fill_shade(ax)
    chart(ax)
    ax.set_ylabel("Photo minus\nCoastSat (m)")
    _title(ax, 1, "Difference in each photo year")
    handles = [Line2D([], [], color=COL["domain"], lw=1.2, marker="o", ms=3, label="CoastSat, annual mean"),
               Line2D([], [], ls="none", marker="s", ms=4.5, color=COL["photos"], mec="white",
                      label="wet/dry photo (shifted)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    png = out("7_photos_vs_coastsat")
    save(fig, png, close=True)
    v = ph.set_index("year").loc[2004, "photo_minus_coastsat_m"]
    record_caption(png, (
        "The photo record against CoastSat, GIS 6 minus GIS 5. (a) The CoastSat annual gap and the "
        "wet/dry photo gap, the photos shifted to the CoastSat mean over the shared 1984-2016 years "
        "(their datums differ, so only shape compares). (b) Photo minus CoastSat in each photo year. "
        f"The 2004 survey (red) sits {v:+.0f} m from CoastSat's 2004 mean, the largest pre-fill "
        "difference; it is the survey behind the reading that the gap held until 2004, and the date "
        "that pinned b in the 2026-10-05 calibration fit. A photo is one day; a CoastSat value is a "
        "year's mean."))


# 8. The failure schedules that follow
def fig_schedules():
    yrs = np.arange(1984, 2026)
    sched = {
        "pinned now: full to 2004, then f": np.where(yrs < 2004, 1.0, F_EXAMPLE),
        "instant from 1996": np.where(yrs < 1996, 1.0, F_EXAMPLE),
        "ramp 1996 to 2003": np.interp(yrs, [1995, 2003], [1.0, F_EXAMPLE]),
    }
    style = {0: (C["BASE"], "-"), 1: (C_1997, "--"), 2: (C["ACCENT"], ":")}
    fig, ax = plt.subplots(figsize=figsize("double", height=2.8), constrained_layout=True)
    ax.axvspan(1996, 2009, color=C_1997, alpha=0.08, lw=0)
    ax.axvspan(2009, 2025, color=C["ADDED"], alpha=0.08, lw=0)
    ax.text(2002.5, 1.13, "calibration 1996-2009", ha="center", fontsize=7, color=INK_MUTED)
    ax.text(2017, 1.13, "test 2009-2025", ha="center", fontsize=7, color=INK_MUTED)
    for k, (lab, v) in enumerate(sched.items()):
        c, ls = style[k]
        ax.step(yrs, v, where="post", color=c, ls=ls, lw=1.6, label=lab)
    ax.set_ylim(0, 1.2)
    ax.set_xlim(1984, 2025)
    ax.set_ylabel("Blocking strength / b")
    events(ax, label=False, years=(1995, 2003))
    chart(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    png = out("8_candidate_schedules")
    save(fig, png, close=True)
    record_caption(png, (
        "The three failure schedules for the refit, as the fraction of the full blocking strength b "
        f"the groin applies each model year, drawn with the pinned post-failure fraction f = {F_EXAMPLE} "
        "(the refit fits b and f for each). Grey: the schedule in the model now, full strength to the "
        "2004 step. Blue dashed: failure from 1996, the first year after the 1995 last repair, where "
        "CoastSat puts the break. Purple dotted: a linear ramp from 1996 to Isabel (2003). Shading: the "
        "calibration and test windows. Dotted lines: the 1995 repair and Isabel."))


def main():
    apply_style()
    for f in (fig_map, fig_sides, fig_gap, fig_break, fig_eras, fig_robust, fig_photos, fig_schedules):
        f()
        print(f"drew {f.__name__}")


if __name__ == "__main__":
    main()
