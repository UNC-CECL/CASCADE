#!/usr/bin/env python3
"""Summary figures documenting the 2026-08-30 calibration, and how it was tested.

WHY THESE EXIST
    The reasoning behind the source/sink and groin calibrations lives in prose,
    in comments in scripts/hatteras_site_config.py. That is the right home for
    the conclusions, but two things were not recorded anywhere at all:

      * the BE convergence sequence. `--overwrite` REPLACES a run's row in
        run_index.csv, so only the last pass survives there, and
        convergence_history.json still carries 2026-08-24 baselines from the
        pre-restructure topography. The pass-by-pass numbers existed only as
        text typed into a comment.
      * the three-target comparison. That the groin's fitted M is set by the
        CHOICE OF TARGET rather than by the data is the methodological result
        of the exercise, and nothing on disk showed it.

    These figures are the durable record of both.

WHAT EACH ONE SHOWS
    fig_three_targets.png    the same 61 sweep cells scored three ways. The
                             fillet says M = 95, the D1-D12 profile says M = 0,
                             D4-D8 demeaned says M = 60. Identical model runs.
    fig_Mf_identifiability.png  D4-D8 demeaned RMSE over the (M, f) grid, with
                             iso-M*f contours. Built to TEST the claim that
                             "only the product M*f is identified" -- and it
                             refuted it: corr(RMSE, M*f) = -0.07 against
                             +0.61 for M and -0.49 for f, and equal-product
                             cells score 10.4 to 12.5 m. But the REPLACEMENT
                             claim ("M and f each weakly constrained") was
                             also wrong: per GROIN_PLAN.md the invariant is
                             period-1 cumulative trapping, M(15.5 + 4.5f).
    fig_top_profiles.png     the top cells and the no-groin baseline against
                             the observed change profile, fit window marked.
    fig_be_convergence.png   interior RMSE per calibration pass, both periods,
                             with the GIS 90 re-solve marked.
    fig_period2_and_bug.png  why period 2 is not fitted, and what the
                             topography-product bug was worth.

A CAVEAT THAT NO LONGER APPLIES, WITHDRAWN 2026-08-31
    This said the D1-D12 panel came from the 40-year window on the PRE-FIX
    topography while the other two were period 1 on the corrected one, so
    the three were not a controlled comparison. True when written; not true
    now. The fullperiod sweep was re-run 2026-08-30 18:20, five hours AFTER
    the worker topography fix (562c75c, 13:01), and this figure was rebuilt
    at 23:52 from those cells.

    VERIFIED, not assumed: re-running cell M60_f0.50 through the worker
    reproduces its stored result to 8.3e-05 on rates of ~2.9 m/yr -- the
    SAME noise floor a period-2 cell shows against itself (1.0e-04 between
    two re-runs), which the 1984-2024 window inherits because it contains
    period 2. A wrong-island run would differ in the first or second
    decimal, not the fifth.

    All three panels are on the same corrected topography.

STYLE, 2026-09-11
    Under the project house style (`scripts/hat_figure_style.py`), which
    replaced this file's own INK/MUTED/ACCENT/FOIL palette and its local
    rcParams block. Three consequences worth knowing before reading an older
    copy of these images side by side with a new one:

      * The canvases were 8.6-15 in wide and are now a 190 mm printed column,
        so the type is the size it claims to be on a page.
      * The two calibration periods are drawn in the house VINTAGE pair -- the
        earlier period red, the later blue -- in every panel that shows both.
        They were pink/teal here and pink meant "the answer" elsewhere in the
        same figure set.
      * The suptitles, the italic per-panel verdicts and the footnote
        paragraphs are off the canvas and in CAPTIONS.md beside the images.
        The five captions carry every sentence they used to, so nothing in the
        argument was dropped to make room.

Usage:
    python HAT_calibration_summary_figures.py

Writes output/groin_sweep/figures/ (untracked; the PDFs beside the PNGs are).
Reasoning and results: CALIBRATION_FIGURES.md, beside this file.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
for _p in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from HAT_fullperiod_target import observed_change_profile  # noqa: E402
from hat_figure_style import (apply_style, C, C_1984, C_1997,  # noqa: E402
                              INK, INK_MUTED, caption, error_cmap, figsize,
                              open_frame, save, _title)

SWEEP = (PROJECT_BASE_DIR / "output" / "groin_sweep" / "1984_2004_edgeBE"
         / "sweep_results.jsonl")
# The LIVE full-period sweep, not the 2026-08-28 archive. This pointed into
# `superseded_20260828/`, which is a "do not use for analysis" tree, and
# that one line was the only thing keeping 12 GB of superseded output
# undeletable. The live file carries the same 16 columns and the same 43 rows.
FULLPERIOD = (PROJECT_BASE_DIR / "output" / "groin_sweep"
              / "fullperiod_1984_2024" / "results.csv")
OUT = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

PINNED_BE1 = -42.6
FIT_DOMAINS = list(range(4, 9))
SHOW_DOMAINS = list(range(1, 13))
PERIOD_YEARS = 20.0

# Semantic colours, from the house palette. ACCENT is the target that works
# and the answer it gives; BASE is a target that fails or a baseline; REF is a
# reference construction laid over the data (the iso-product curves, a
# tolerance line); the VINTAGE pair is the two calibration periods.
ANSWER, FOIL, REF = C["ACCENT"], C["BASE"], C["REF"]
BAND = "0.94"

# Five top cells as one family: a ramp of the accent, so they read as
# variations of the same thing rather than five unrelated series. viridis was
# used here until 2026-09-11 and put the best cell in the same green as the
# reference curves in the neighbouring figure.
TOP_CMAP = LinearSegmentedColormap.from_list(
    "hat_accent_ramp", [C["ACCENT_FILL"], C["ACCENT"]])


def load_cells():
    """Sweep cells at the calibrated be1, with all three scores attached."""
    rows = [json.loads(l) for l in open(SWEEP, encoding="utf-8") if l.strip()]
    d = pd.DataFrame(rows)
    d = d[(d.be1 == PINNED_BE1) & d.differential_err.notna()].copy()

    obs = np.array([observed_change_profile(1984, 2004, FIT_DOMAINS)[k]
                    for k in FIT_DOMAINS])
    obs_dm = obs - obs.mean()

    def d48(r):
        # Rates are m/yr; the observed target is CHANGE over the period, so the
        # model side is scaled by the period length. Demeaned because a uniform
        # level offset in the groin's neighbourhood belongs to the source/sink
        # term, not to the groin -- what the groin must get right is the shape.
        m = np.array([r[f"rate_D{k}"] for k in FIT_DOMAINS]) * PERIOD_YEARS
        m = m - m.mean()
        return float(np.sqrt(((m - obs_dm) ** 2).mean()))

    d["score_d48"] = d.apply(d48, axis=1)
    return d, obs_dm


def fig_three_targets(d):
    """The same cells, scored three ways, with each minimum marked."""
    fp = pd.read_csv(FULLPERIOD)
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", aspect=0.36),
                             constrained_layout=True)
    panels = [
        (axes[0], d.groupby("M").differential_err.min(),
         "fillet, D5−D6", "|error| (m/yr)", FOIL),
        (axes[1], fp.groupby("M").rmse_m.min(),
         "profile, D1−D12", "RMSE (m)", FOIL),
        (axes[2], d.groupby("M").score_d48.min(),
         "D4−D8, demeaned", "RMSE (m)", ANSWER),
    ]
    bests = []
    for i, (ax, series, title, ylab, colour) in enumerate(panels):
        ax.plot(series.index, series.values, marker="o", ms=3.0, lw=1.4,
                color=colour)
        best = series.idxmin()
        bests.append(best)
        ax.axvline(best, color=colour, lw=0.8, ls=(0, (3, 2)))
        ax.annotate(f"min at M = {best:g}", xy=(best, series.min()),
                    xytext=(4, 10), textcoords="offset points",
                    fontsize=7.5, color=colour)
        _title(ax, i, title)
        ax.set_xlabel("groin trapping M (m/yr)")
        ax.set_ylabel(ylab)
        ax.grid(axis="y")
        ax.set_axisbelow(True)
        open_frame(ax)

    caption(fig,
            "The same model runs, scored three ways. (a) The fillet, the "
            "D5−D6 scalar: no admissible M matches it, and fitting it anyway "
            "rails f at the grid bound. (b) The D1−D12 change profile: ranks "
            "M = {b:g} best and monotonically, because D2−D4 is Cape Point "
            "accretion the parameterisation does not represent and D6−D7 is an "
            "erosion trough peaking one domain north of the structure, which a "
            "groin actively worsens; a groin signal of about 17 m is swamped. "
            "(c) D4−D8 demeaned, the target that works: a clean interior "
            "minimum at M = {c:g}, rising on both sides through M = 160. "
            "Demeaning removes the level offset the source/sink term owns, and "
            "D1 is excluded because the cape's 81−104 m change is about five "
            "times the groin's signal. The fitted groin is therefore set by "
            "the target, not by the data. (a) and (c) are period 1 at "
            "be1 = {be:g} on the corrected topography, {n} cells; (b) is the "
            "40-year continuous window on the SAME corrected topography, "
            "re-run 2026-08-30 after the worker's topography fix and verified "
            "2026-08-31 to reproduce to 8e-05."
            .format(b=bests[1], c=bests[2], be=PINNED_BE1, n=len(d)))

    p = OUT / "fig_three_targets.png"
    print("  {}".format(save(fig, p, close=True)[0].name))


def fig_identifiability(d):
    """Is the valley along constant M*f? Tests the config's claim."""
    piv = d.pivot_table(index="fraction", columns="M", values="score_d48")
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.50),
                           constrained_layout=True)
    mesh = ax.pcolormesh(piv.columns, piv.index, piv.values,
                         cmap=error_cmap(), shading="auto")
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label("D4−D8 demeaned RMSE (m)")
    cb.outline.set_linewidth(0.6)

    # Iso-M*f curves. If the valley follows these, only the product is fitted.
    Mg = np.linspace(max(piv.columns.min(), 1), piv.columns.max(), 200)
    for prod in (30, 40, 50, 60):
        f = prod / Mg
        keep = (f >= piv.index.min()) & (f <= piv.index.max())
        ax.plot(Mg[keep], f[keep], color=REF, lw=1.0, ls=(0, (4, 2)))
        if keep.any():
            ax.annotate(f"M·f={prod}", xy=(Mg[keep][-1], f[keep][-1]),
                        fontsize=7, color=REF, ha="right",
                        bbox=dict(facecolor="white", alpha=0.75,
                                  edgecolor="none",
                                  boxstyle="square,pad=0.12"))

    best = d.loc[d.score_d48.idxmin()]
    ax.plot(best.M, best.fraction, marker="*", ms=13, color=REF,
            markeredgecolor="white", markeredgewidth=0.7, zorder=5,
            label=f"best cell, M {best.M:g}, f {best.fraction:g}")
    ax.plot(60, 0.6, marker="o", ms=8, color="none", markeredgecolor=REF,
            markeredgewidth=1.6, zorder=5,
            label="the pair taken forward, M 60, f 0.6")
    ax.legend(loc="upper right")
    ax.set_xlabel("groin trapping M (m/yr)")
    ax.set_ylabel("deterioration floor f")
    ax.set_title("D4−D8 demeaned error over the (M, f) grid", loc="left")

    caption(fig,
            "Built to test the claim that only the product M·f is identified, "
            "and it refutes it: the dashed iso-product curves cut across the "
            "valley rather than following it. Equal-product cells are not "
            "interchangeable — at M·f of about 40 the error spans 10.4 to "
            "12.5 m, wider than the 3.8 m the groin buys — and the "
            "correlations are −0.07 for M·f against +0.61 for M and −0.49 for "
            "f, so the target responds to M and f separately. The replacement "
            "claim, that M and f are each only weakly constrained, is wrong "
            "too: the invariant is period-1 cumulative trapping, "
            "M(15.5 + 4.5f). Period 1 mostly precedes the 1996−2003 "
            "deterioration ramp, so f moves it by only 29% across its whole "
            "range, while period 2 is 20·M·f and f = 0 gives zero there. M is "
            "therefore set from period 1 and f from the 1967 rig and period 2, "
            "and because the best M is 50 at f = 1.0 but 70 at f = 0.6, a poor "
            "score at M = 50, f = 0.6 means M was too low for that f, not that "
            "f = 0.6 is wrong.")

    p = OUT / "fig_Mf_identifiability.png"
    print("  {}".format(save(fig, p, close=True)[0].name))


def fig_top_profiles(d, obs_dm):
    """Top cells and the no-groin baseline against the observed profile."""
    # LANDWARD-POSITIVE, so erosion is UP and this panel reads as a plan view,
    # matching the gifs and the other profile figures. Both sources are
    # SEAWARD-positive, so both are negated -- at the PLOTTING layer only.
    # score_d48 is computed upstream and is unaffected.
    # fig_three_targets needs no flip: it plots error against M, not a profile.
    obs_full = -np.array([observed_change_profile(1984, 2004, SHOW_DOMAINS)[k]
                          for k in SHOW_DOMAINS])
    fit_idx = [SHOW_DOMAINS.index(k) for k in FIT_DOMAINS]

    def centred(v):
        return v - np.asarray(v)[fit_idx].mean()

    top = d.nsmallest(5, "score_d48")
    nog = d[d.M == 0]
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.46),
                           constrained_layout=True)
    ax.axvspan(min(FIT_DOMAINS) - 0.5, max(FIT_DOMAINS) + 0.5,
               color=BAND, zorder=0)
    ax.annotate("fit window, D4−D8", xy=(6, 0.02),
                xycoords=("data", "axes fraction"),
                ha="center", fontsize=7.5, color=INK_MUTED)
    ax.axvline(5.5, color=INK_MUTED, lw=0.8, ls=(0, (4, 2)), zorder=2)
    ax.annotate("Buxton groin", xy=(5.5, 0.95),
                xycoords=("data", "axes fraction"),
                rotation=90, ha="right", va="top", fontsize=7.5,
                color=INK_MUTED)

    ax.plot(SHOW_DOMAINS, centred(obs_full), marker="s", ms=4.0, lw=1.8,
            ls="--", color=INK, label="observed, 1984 to 2004", zorder=6)
    if len(nog):
        v = [-nog.iloc[0][f"rate_D{k}"] * PERIOD_YEARS for k in SHOW_DOMAINS]
        ax.plot(SHOW_DOMAINS, centred(v), lw=1.4, ls=":", color=FOIL,
                label=f"no groin, {nog.score_d48.iloc[0]:.1f} m", zorder=4)
    # Darkest first: `top` is sorted best to worst, and the caption says the
    # best cell is the darkest.
    for colour, (_, r) in zip(TOP_CMAP(np.linspace(1.0, 0.15, len(top))),
                              top.iterrows()):
        v = [-r[f"rate_D{k}"] * PERIOD_YEARS for k in SHOW_DOMAINS]
        ax.plot(SHOW_DOMAINS, centred(v), lw=1.3, color=colour,
                label=f"M {r.M:g}, f {r.fraction:g}, {r.score_d48:.1f} m",
                zorder=5)

    ax.set_xlabel("GIS domain")
    ax.set_ylabel("shoreline change 1984 to 2004 (m)\ndemeaned over the fit"
                  " window; positive is landward")
    ax.set_xticks(SHOW_DOMAINS)
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    # Outside the axes: seven entries inside the panel covered the observed
    # line across D1-D3, which is the part of the profile the caption is about.
    fig.legend(loc="outside lower center", ncol=4, frameon=False, fontsize=7.5)
    ax.set_title("The five best cells and the no-groin baseline", loc="left")

    caption(fig,
            "The groin improves the fit but does not reproduce the shape. "
            "Inside the fit window the observed profile peaks at D6, dips at "
            "D7 and peaks again at D8; every cell draws a smooth monotonic "
            "rise, and the error gain comes from matching the overall D4−D8 "
            "slope. Read the residual as the split between what the groin "
            "explains and what the source/sink calibration absorbs, not as a "
            "successful shape fit. Everything is demeaned over the fit window "
            "and plotted landward-positive, so erosion is up and the panel "
            "reads as a plan view; both sources are seaward-positive and are "
            "negated at the plotting layer only. The five cells are drawn in "
            "one colour family because they are variations of one thing, with "
            "the best cell darkest. Absolute values here are not comparable "
            "with the production ranking figure: this script approximates "
            "change as rate times 20 years from the sweep results, while the "
            "production script reads the shoreline matrix directly. The "
            "ranking and the shape of the curves are unaffected, and the "
            "production figure is the one to quote.")

    p = OUT / "fig_top_profiles.png"
    print("  {}".format(save(fig, p, close=True)[0].name))


# Transcribed from the pass-by-pass table in hatteras_site_config.py. NOT
# derivable from run_index.csv: --overwrite replaces a run's row, so only the
# final pass survives there. This figure is the durable record.
BE_PASSES = ["zeroBE", "edgeBE", "calib\npass 0", "pass 1", "pass 2",
             "+GIS 90\nre-solve", "final"]
BE_RMSE = {1984: [1.422, 1.219, 0.721, 0.547, 0.527, 0.556, 0.517],
           2004: [2.124, 1.794, 0.763, 0.615, 0.583, 0.603, 0.563]}


def fig_be_convergence():
    """Interior RMSE per calibration pass, both periods."""
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.42),
                           constrained_layout=True)
    x = np.arange(len(BE_PASSES))
    for period, colour, mark in ((1984, C_1984, "o"), (2004, C_1997, "s")):
        ax.plot(x, BE_RMSE[period], marker=mark, ms=4.0, lw=1.6, color=colour,
                label=f"{period} to {period + 20}")
        ax.annotate(f"{BE_RMSE[period][-1]:.3f}",
                    xy=(x[-1], BE_RMSE[period][-1]),
                    xytext=(6, -2), textcoords="offset points",
                    fontsize=7.5, color=colour)
    ax.axvspan(4.5, 5.5, color=BAND, zorder=0)
    ax.annotate("edges gained,\ninterior gave back", xy=(5, 1.55),
                ha="center", fontsize=7.5, color=INK_MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels(BE_PASSES, fontsize=7.5)
    ax.set_ylabel("interior RMSE against the CoastSat rate (m/yr)")
    ax.set_xlabel("calibration stage")
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    open_frame(ax)
    ax.legend(loc="upper right")
    ax.set_title("Source/sink calibration, pass by pass", loc="left")

    caption(fig,
            "The source/sink calibration converged in five passes in both "
            "periods, from a base run with the road and the beach/dune "
            "management on and the groin off. The stopping rule, recorded in "
            "convergence_history.json, is that a pass must buy more than 5% of "
            "the standing error; it was met from pass 2 to pass 3. The shaded "
            "stage is the GIS 90 re-solve, which cost the interior 0.03 while "
            "buying 0.61 at the edge, about 20 to 1, so the edge/interior loop "
            "contracts rather than oscillating. These numbers are TRANSCRIBED "
            "from the pass-by-pass table in the site config, not re-derived: "
            "--overwrite replaces a run's row in run_index.csv, so only the "
            "final pass survives there, and convergence_history.json still "
            "carries 2026-08-24 baselines from the pre-restructure topography. "
            "This figure is their durable record, and if one of them was "
            "mistyped it is mistyped here too.")

    p = OUT / "fig_be_convergence.png"
    print("  {}".format(save(fig, p, close=True)[0].name))


def fig_period2_and_bug(d):
    """Why period 2 is unfittable, and what the topography bug was worth."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=figsize("double", aspect=0.42),
                                 constrained_layout=True)

    reach_lo = float(d.differential_m_yr.min())
    reach_hi = float(d.differential_m_yr.max())
    a1.axhspan(reach_lo, reach_hi, color=BAND, zorder=0)
    a1.annotate("reachable by the module,\ntrapping at or above zero",
                xy=(0.5, (reach_lo + reach_hi) / 2), ha="center",
                fontsize=7.5, color=INK_MUTED)
    a1.errorbar([0], [3.46], yerr=[0.70], marker="o", ms=5, lw=1.4,
                color=C_1984, capsize=3, label="observed, 1984 to 2004")
    a1.errorbar([1], [-3.85], yerr=[0.76], marker="s", ms=5, lw=1.4,
                color=C_1997, capsize=3, label="observed, 2004 to 2024")
    a1.axhline(0, color=INK, lw=0.8)
    a1.set_xlim(-0.6, 1.6)
    a1.set_xticks([0, 1])
    a1.set_xticklabels(["1984 to 2004", "2004 to 2024"])
    a1.set_ylabel("fillet trend, D5−D6 (m/yr)")
    _title(a1, 0, "period 2 against the module's range")
    a1.grid(axis="y")
    a1.set_axisbelow(True)
    open_frame(a1)
    a1.legend(loc="lower left")

    dom = np.arange(1, 13)
    drift = np.zeros(12)
    drift[3] = 0.1258
    a2.bar(dom, drift, color=FOIL, width=0.65)
    a2.axhline(0.005, color=REF, lw=1.0, ls=(0, (4, 2)))
    a2.annotate("guard tolerance 0.005", xy=(12, 0.007), ha="right",
                fontsize=7.5, color=REF)
    a2.annotate("after the fix, exactly zero\nat every domain", xy=(8, 0.075),
                ha="center", fontsize=7.5, color=INK_MUTED)
    a2.set_xlabel("GIS domain")
    a2.set_ylabel("|sweep − published| (m/yr)")
    a2.set_xticks(dom)
    _title(a2, 1, "the topography-product bug")
    a2.grid(axis="y")
    a2.set_axisbelow(True)
    open_frame(a2)

    caption(fig,
            "(a) Why period 2 is not fitted. A groin whose trapping is bounded "
            "at or above zero can only WIDEN the updrift-to-downdrift gap. "
            "Period 1's observed gap widens, period 2's narrows, so no M at or "
            "above zero reaches it and there is no joint fit; the 2026-08-30 "
            "joint fit railing at M = 160 was not a failure to repair, because "
            "fitting period 2 is the wrong thing to attempt. The shaded band "
            "is the range the sweep cells actually span. The groin still RUNS "
            "in period 2: GroinCallback carries an absolute calendar timeline, "
            "so no period-specific configuration exists, and running it there "
            "is right for consistency of the structure's timeline rather than "
            "because it explains that period's shoreline. (b) What the "
            "topography-product bug was worth, at 25 times the guard "
            "tolerance: the sweep worker resolved topography without naming a "
            "product, so the 2004-start default answered and 1984 cells were "
            "built on the 2004 island. The single non-zero domain is the "
            "measured drift before the fix; after it, every domain matches the "
            "published run exactly.")

    p = OUT / "fig_period2_and_bug.png"
    print("  {}".format(save(fig, p, close=True)[0].name))


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    d, obs_dm = load_cells()
    print(f"{len(d)} cells at be1 = {PINNED_BE1}")
    fig_three_targets(d)
    fig_identifiability(d)
    fig_top_profiles(d, obs_dm)
    fig_be_convergence()
    fig_period2_and_bug(d)


if __name__ == "__main__":
    main()
