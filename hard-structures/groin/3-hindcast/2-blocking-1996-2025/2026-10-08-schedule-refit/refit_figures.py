#!/usr/bin/env python3
"""
Figures for the schedule refit: the three schedules, the b-f misfit maps, and the gap trajectories.

    python refit_figures.py   ->  figures/schedule_refit_<n>_<what>.png

Reads grid_scores.csv and coastsat_gap_change.csv from schedule_refit.py score and the
runs' shoreline matrices. Best cell per schedule = lowest annual CoastSat RMSE.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
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
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import schedule_refit as sr  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, INK, INK_MUTED, _title, apply_style, figsize, open_frame, record_caption,
    save)

# --- CONFIG ------------------------------------------------------------------
FIGS = HERE / "figures"
LABEL = {"instant2004": "failure at 2004 (as pinned)",
         "instant1996": "failure from 1996",
         "ramp1996": "ramp 1996 to 2003"}
STYLE = {"instant2004": (C["BASE"], "-"), "instant1996": (C_1997, "--"), "ramp1996": (C["ACCENT"], "-."),
         "pinned": (C["ADDED"], "-")}
PRODUCT_STYLE = ((C_1997, "-"), (C["ACCENT"], "--"), (C["ADDED"], ":"))
PINNED = ("instant2004", 0.6, 0.6)
# -----------------------------------------------------------------------------


def scores():
    return pd.read_csv(HERE / "grid_scores.csv")


def best(t, period=sr.CALIBRATION, metric="rmse_annual"):
    g = t[(t.period == period) & (t.schedule != "none")]
    return {s: sub.nsmallest(1, metric).iloc[0] for s, sub in g.groupby("schedule")}


# 1. Misfit over b and f, per schedule: annual CoastSat (top) and photo dates (bottom)
def fig_maps(t):
    cal = t[(t.period == sr.CALIBRATION) & (t.schedule != "none")]
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.6), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, (metric, lab) in enumerate((("rmse_annual", "annual CoastSat"), ("rmse_photos", "photo dates"))):
        vmax = np.nanpercentile(cal[metric], 95)
        for col, s in enumerate(sr.SCHEDULES):
            ax = axes[row, col]
            p = cal[cal.schedule == s].pivot(index="b", columns="f", values=metric)
            im = ax.imshow(p.values, origin="lower", cmap="Blues_r", vmin=0, vmax=vmax, aspect="auto")
            for i in range(p.shape[0]):
                for j in range(p.shape[1]):
                    v = p.values[i, j]
                    if np.isfinite(v):
                        ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=6,
                                color="white" if v < 0.45 * vmax else INK)
            i, j = np.unravel_index(np.nanargmin(p.values), p.shape)
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, ec=C_1984, lw=1.6))
            ax.set_xticks(range(p.shape[1]), [f"{v:g}" for v in p.columns])
            ax.set_yticks(range(p.shape[0]), [f"{v:g}" for v in p.index])
            if row == 1:
                ax.set_xlabel("f, post-failure fraction")
            if col == 0:
                ax.set_ylabel(f"b, blocking strength\n({lab} RMSE)")
            _title(ax, 3 * row + col, LABEL[s] if row == 0 else "")
        fig.colorbar(im, ax=axes[row, :], shrink=0.85, label=f"{lab} RMSE (m)")
    png = FIGS / "schedule_refit_1_misfit_maps.png"
    save(fig, png, close=True)
    base = t[(t.period == sr.CALIBRATION) & (t.schedule == "none")].iloc[0]
    record_caption(png, (
        "Calibration 1996-2009 misfit of the blocking groin over strength b and post-failure fraction "
        "f, for the three failure schedules (columns). Top: RMSE of the modelled GIS 5|6 gap change "
        "against the annual CoastSat gap change, 1996-2008, model mid-year against the calendar-year "
        "mean, both relative to the start. Bottom: RMSE at the three wet/dry photo dates (1997, 2004, "
        "2008), the 2026-10-05 score. Red box: the best cell. Each row shares one colour scale "
        f"(capped at the 95th percentile). No groin: {base.rmse_annual:.0f} m annual, "
        f"{base.rmse_photos:.0f} m photos. Under instant failure from 1996 the groin is at b x f for "
        "the whole window, so only the product is constrained (equal products, equal runs)."))


# Model gap change through a run, seaward-positive, from its shoreline matrix
def model_gap(path, start):
    x = np.load(path)
    g = x[:, sr.DOWN] - x[:, sr.UP]
    return start + np.arange(len(g)), g - g[0]


def photos_change(start):
    from score_instant_grid import observed_series
    p = observed_series()
    return p - np.interp(start, p.index, p.values)


# 2/3. Gap trajectories for the best cell per schedule, the pinned cell and no groin
def fig_trajectories(t, period, picks, name, letter_title):
    obs = pd.read_csv(HERE / "coastsat_gap_change.csv")
    obs = obs[obs.period == period]
    from site_layer.hatteras_site_config import run_years
    end = period + run_years(period)
    o = obs[(obs.year >= period - 1) & (obs.year <= end)]
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6), constrained_layout=True)
    ax.plot(o.year + 0.5, o.coastsat_change_m, color=C_1997, lw=0, marker="o", ms=4, alpha=0.9,
            zorder=5)
    ph = photos_change(period)
    ph = ph[(ph.index >= period - 1) & (ph.index <= end)]
    ax.plot(ph.index + 0.5, ph.values, ls="none", marker="s", ms=5, color=C_1984, mec="white",
            mew=0.5, zorder=6)
    yrs, g = model_gap(sr.baseline_file(period), period)
    ax.plot(yrs, g, color=INK_MUTED, lw=1.0, ls=":")
    handles = [Line2D([], [], color=C_1997, marker="o", ms=4, lw=0, label="CoastSat, annual mean"),
               Line2D([], [], color=C_1984, marker="s", ms=5, lw=0, mec="white", label="wet/dry photos"),
               Line2D([], [], color=INK_MUTED, ls=":", label="no groin")]
    for k, (s, b, f, lab, style) in enumerate(picks):
        p = sr.shoreline_file(s, period, b, f)
        if p is None:
            continue
        c, ls = style
        yrs, g = model_gap(p, period)
        ax.plot(yrs, g, color=c, ls=ls, lw=1.5)
        handles.append(Line2D([], [], color=c, ls=ls, lw=1.5, label=lab))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    for yr, lab in ((1995, "last repair"), (2003, "Isabel"), (2017, "Buxton fill")):
        if period - 1 <= yr <= end:
            ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":")
            ax.text(yr + 0.2, 1.0, lab, transform=ax.get_xaxis_transform(), fontsize=6.5,
                    color=INK_MUTED, va="top")
    ax.set_xlim(period - 1, end + 0.5)
    ax.set_ylabel("Gap change since the start (m)\nupdrift seaward ▲")
    ax.grid(axis="y")
    open_frame(ax)
    ax.set_title(letter_title, loc="center")
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False, fontsize=7)
    png = FIGS / name
    save(fig, png, close=True)
    return png


# 0. The three failure schedules, as the model applies them (moved from the condition analysis 2026-10-08)
def fig_schedules(f_example=0.6):
    from cascade.groin import _scheduled_strength
    yrs = np.arange(1984, 2026)
    fig, ax = plt.subplots(figsize=figsize("double", height=2.9), constrained_layout=True)
    ax.axvspan(sr.CALIBRATION, sr.TEST, color=C_1997, alpha=0.08, lw=0)
    ax.axvspan(sr.TEST, 2025, color=C["ADDED"], alpha=0.08, lw=0)
    ax.text((sr.CALIBRATION + sr.TEST) / 2, 1.12, f"calibration {sr.CALIBRATION}-{sr.TEST}", ha="center",
            fontsize=7, color=INK_MUTED)
    ax.text((sr.TEST + 2025) / 2, 1.12, f"test {sr.TEST}-2025", ha="center", fontsize=7, color=INK_MUTED)
    for s, (delay, mode, ramp) in sr.SCHEDULES.items():
        v = [_scheduled_strength(1.0, y, sr.INSTALL + delay, mode, f_example, ramp) for y in yrs]
        c, ls = STYLE[s]
        ax.step(yrs, v, where="post", color=c, ls=ls, lw=1.6, label=LABEL[s])
    for yr, lab in ((1995, "last repair"), (2003, "Isabel")):
        ax.axvline(yr, color=INK_MUTED, lw=0.6, ls=":")
        ax.text(yr - 0.25, 0.04, lab, ha="right", fontsize=6.5, color=INK_MUTED)
    ax.set_ylim(0, 1.2)
    ax.set_xlim(1984, 2025)
    ax.set_ylabel("Blocking strength,\nfraction of full b")
    ax.set_xlabel("Model year")
    ax.grid(axis="y")
    open_frame(ax)
    ax.set_title("When does the groin weaken, under each schedule tested?")
    ax.legend(frameon=False, fontsize=7.5, loc="lower left", bbox_to_anchor=(0.0, 0.08))
    png = FIGS / "schedule_refit_0_failure_schedules.png"
    save(fig, png, close=True)
    record_caption(png, (
        "**Tests:** nothing; this is the set-up the refit compares. **How to read:** the fraction of the "
        "full blocking strength b that the groin applies in each model year under each failure "
        "schedule, computed with the model's own schedule function (`cascade.groin`), drawn with the "
        f"post-failure fraction f = {f_example} for illustration (the refit fits b and f for each "
        "schedule). Grey: the schedule pinned in the model, full strength until the 2004 step. Blue "
        "dashed: failure from 1996, the first year after the 1995 last repair, where CoastSat puts the "
        "end of trapping (`hard-structures/groin/1-observations/gap_across_groins/`). Purple: a "
        "linear ramp from full strength in 1995 to the floor in 2003 (Isabel). Shading: the calibration "
        "and test periods. **Shows:** the schedules differ only inside the calibration period; from 2004 "
        "on all three apply b × f, so test and forward runs depend on the product alone."))


def main():
    apply_style()
    t = scores()
    fig_schedules()
    fig_maps(t)
    bc = best(t)
    picks = [(s, r.b, r.f, f"{LABEL[s]}: best, b {r.b:g} f {r.f:g}", STYLE[s]) for s, r in bc.items()]
    if PINNED not in [(s, r.b, r.f) for s, r in bc.items()]:
        picks.append((*PINNED, "failure at 2004: pinned now, b 0.6 f 0.6", STYLE["pinned"]))
    png = fig_trajectories(t, sr.CALIBRATION, picks, "schedule_refit_2_calibration_gap.png",
                           "Calibration 1996-2009: gap across the groin")
    record_caption(png, (
        "The GIS 5|6 gap change through the calibration run, seaward positive (updrift holding is "
        "up), each relative to its start. Lines: model at 1 Jan of each year; no groin dotted; the "
        "best cell of each schedule on the annual CoastSat score; the pinned cell (failure at 2004, "
        "b 0.6 f 0.6, amber). Blue dots: CoastSat calendar-year means, relative to the mean over the "
        "DEM-centred start window (Oct 1995-Oct 1997), plotted mid-year. Red squares: the wet/dry "
        "photo gap relative to its value interpolated to 1996, plotted mid-year."))
    # After 2009 every schedule has failed, so a test run depends on b x f only: one line per product
    cells = [(s, r.b, r.f, f"{LABEL[s]} best") for s, r in bc.items()] + [(*PINNED, "the pin")]
    by_product = {}
    for s, b, f, who in cells:
        by_product.setdefault(round(b * f, 3), []).append((s, b, f, who))
    tests = [(*members[0][:3], f"b × f = {bf:g}: " + ", ".join(dict.fromkeys(m[3] for m in members)),
              PRODUCT_STYLE[k % len(PRODUCT_STYLE)])
             for k, (bf, members) in enumerate(sorted(by_product.items(), reverse=True))]
    if any(sr.shoreline_file(s, sr.TEST, b, f) for s, b, f, *_ in tests):
        png = fig_trajectories(t, sr.TEST, tests, "schedule_refit_3_test_gap.png",
                               "Test 2009-2025: gap across the groin (not fitted)")
        record_caption(png, (
            "The same on the 2009-2025 test period, which no cell was fitted to. Every schedule has "
            "failed by 2009, so a test run depends only on b x f: one line per product, naming the "
            "calibration cells it stands for (their runs are identical). CoastSat relative to the mean over "
            "Aug 2008-Aug 2010; the Buxton 2017 fill lands on GIS 6-16 and its sand moved south past "
            "the groin, which the model cannot do (the +32 m model jump at 2018)."))


if __name__ == "__main__":
    main()
