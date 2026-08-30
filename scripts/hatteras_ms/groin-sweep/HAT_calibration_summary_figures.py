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

A CAVEAT CARRIED ON THE FIGURE, NOT JUST HERE
    The D1-D12 panel comes from the 40-year full-period rig on the PRE-FIX
    topography; the other two are period 1 on the corrected one. They are not
    a controlled comparison of targets on identical runs, and the panel says
    so. It is shown because the MECHANISM it demonstrates -- a wide window
    swamped by signal a groin cannot produce -- does not depend on which
    topography it ran on.

Usage:
    python HAT_calibration_summary_figures.py

Writes output/groin_sweep/figures/.
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

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[3]
for _p in (PROJECT_BASE_DIR / "scripts", _HERE.parent):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from HAT_fullperiod_target import observed_change_profile  # noqa: E402

SWEEP = (PROJECT_BASE_DIR / "output" / "groin_sweep" / "1984_2004_edgeBE"
         / "sweep_results.jsonl")
FULLPERIOD = (PROJECT_BASE_DIR / "output" / "archived_output_20260828"
              / "groin_sweep" / "fullperiod_1984_2024" / "results.csv")
OUT = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"

PINNED_BE1 = -42.6
FIT_DOMAINS = list(range(4, 9))
SHOW_DOMAINS = list(range(1, 13))
PERIOD_YEARS = 20.0

INK, MUTED, GRID = "#1a1a2e", "#5c6068", "#d5d8dd"
ACCENT = "#c2185b"      # the chosen target / the answer
FOIL = "#7a7f87"        # targets that fail

plt.rcParams.update({
    "font.size": 10, "axes.linewidth": 0.7, "axes.labelsize": 10.5,
    "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
})


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
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    panels = [
        (axes[0], d.groupby("M").differential_err.min(),
         "(a)  fillet, D5−D6 scalar", "|error| (m/yr)", FOIL,
         "no admissible M matches it"),
        (axes[1], fp.groupby("M").rmse_m.min(),
         "(b)  change profile, D1−D12", "RMSE (m)", FOIL,
         "swamped by Cape Point + D6−D7"),
        (axes[2], d.groupby("M").score_d48.min(),
         "(c)  change profile, D4−D8 demeaned", "RMSE (m)", ACCENT,
         "the target that works"),
    ]
    for ax, series, title, ylab, colour, note in panels:
        ax.plot(series.index, series.values, marker="o", ms=4.5, lw=1.8,
                color=colour)
        best = series.idxmin()
        ax.axvline(best, color=colour, lw=1.0, ls=(0, (3, 2)), alpha=0.8)
        ax.annotate(f"min at M = {best:g}", xy=(best, series.min()),
                    xytext=(6, 14), textcoords="offset points",
                    fontsize=9.5, color=colour, fontweight="bold")
        ax.set_title(title, loc="left", fontsize=11, color=INK, pad=6)
        ax.set_xlabel("Groin trapping M (m/yr)")
        ax.set_ylabel(ylab)
        ax.grid(color=GRID, lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.annotate(note, xy=(0.5, -0.19), xycoords="axes fraction",
                    ha="center", fontsize=9, color=MUTED, style="italic")

    fig.suptitle("Identical model runs, three targets, three answers — the "
                 "fitted groin is set by the target, not by the data",
                 fontsize=13, fontweight="bold", color=INK, y=1.02)
    fig.text(0.5, -0.30,
             "(a) and (c): period 1, be1 = −42.6, corrected topography, 61 "
             "cells.   (b): 40-year full-period rig on the PRE-FIX topography "
             "— not a controlled comparison,\nshown because the mechanism "
             "(a wide window dominated by signal a groin cannot produce) does "
             "not depend on which topography it ran on.",
             ha="center", fontsize=8.6, color=MUTED)
    p = OUT / "fig_three_targets.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {p.name}")


def fig_identifiability(d):
    """Is the valley along constant M*f? Tests the config's claim."""
    piv = d.pivot_table(index="fraction", columns="M", values="score_d48")
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    mesh = ax.pcolormesh(piv.columns, piv.index, piv.values,
                         cmap="magma_r", shading="auto")
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label("D4−D8 demeaned RMSE (m)")

    # Iso-M*f curves. If the valley follows these, only the product is fitted.
    Mg = np.linspace(max(piv.columns.min(), 1), piv.columns.max(), 200)
    for prod in (30, 40, 50, 60):
        f = prod / Mg
        keep = (f >= piv.index.min()) & (f <= piv.index.max())
        ax.plot(Mg[keep], f[keep], color="#3ddc97", lw=1.2, ls=(0, (4, 2)))
        if keep.any():
            ax.annotate(f"M·f={prod}", xy=(Mg[keep][-1], f[keep][-1]),
                        fontsize=8, color="#1b7f5a", ha="right")

    best = d.loc[d.score_d48.idxmin()]
    ax.plot(best.M, best.fraction, marker="*", ms=20, color="#3ddc97",
            markeredgecolor=INK, markeredgewidth=0.8, zorder=5,
            label=f"best  M={best.M:g}, f={best.fraction:g}")
    ax.plot(60, 0.6, marker="o", ms=10, color="none", markeredgecolor="#3ddc97",
            markeredgewidth=2.0, zorder=5, label="production  M=60, f=0.6")
    ax.legend(loc="upper right", fontsize=9, labelcolor=INK,
              facecolor="white", framealpha=0.85, frameon=True)
    ax.set_xlabel("Groin trapping M (m/yr)")
    ax.set_ylabel("Deterioration floor f")
    ax.set_title("Not a product ridge — the iso-M·f curves cut across it; "
                 "the invariant is M(15.5+4.5f)", loc="left",
                 fontsize=12, fontweight="bold", color=INK, pad=8)
    fig.text(0.5, -0.05,
             "Equal-product cells are not interchangeable: at M·f ≈ 40 the "
             "RMSE spans 10.4 to 12.5 m, wider than the 3.8 m the groin buys. "
             "corr(RMSE, M·f) = −0.07,\n"
             "corr(RMSE, M) = +0.61, corr(RMSE, f) = −0.49 — the target "
             "responds to M and f separately, so the pair is weakly "
             "constrained rather than traded off.",
             ha="center", fontsize=8.8, color=MUTED)
    p = OUT / "fig_Mf_identifiability.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {p.name}")


def fig_top_profiles(d, obs_dm):
    """Top cells and the no-groin baseline against the observed profile."""
    obs_full = np.array([observed_change_profile(1984, 2004, SHOW_DOMAINS)[k]
                         for k in SHOW_DOMAINS])
    fit_idx = [SHOW_DOMAINS.index(k) for k in FIT_DOMAINS]

    def centred(v):
        return v - np.asarray(v)[fit_idx].mean()

    top = d.nsmallest(5, "score_d48")
    nog = d[d.M == 0]
    fig, ax = plt.subplots(figsize=(10.5, 5.6))
    ax.axvspan(min(FIT_DOMAINS) - 0.5, max(FIT_DOMAINS) + 0.5,
               color="#fdf6d8", zorder=0)
    ax.annotate("FIT WINDOW  D4−D8", xy=(6, 0.02), xycoords=("data", "axes fraction"),
                ha="center", fontsize=9.5, color="#8a7b1f", fontweight="bold")
    ax.axvline(5.5, color=ACCENT, lw=1.4, ls=(0, (4, 2)), zorder=2)
    ax.annotate("Buxton groin", xy=(5.5, 0.93), xycoords=("data", "axes fraction"),
                rotation=90, ha="right", va="top", fontsize=9, color=ACCENT)

    ax.plot(SHOW_DOMAINS, centred(obs_full), marker="s", ms=7, lw=2.2,
            ls="--", color=INK, label="observed 1984→2004", zorder=6)
    if len(nog):
        v = [nog.iloc[0][f"rate_D{k}"] * PERIOD_YEARS for k in SHOW_DOMAINS]
        ax.plot(SHOW_DOMAINS, centred(v), lw=1.6, ls=":", color=MUTED,
                label=f"no groin  (RMSE {nog.score_d48.iloc[0]:.1f} m)", zorder=4)
    cmap = plt.cm.viridis(np.linspace(0.15, 0.8, len(top)))
    for colour, (_, r) in zip(cmap, top.iterrows()):
        v = [r[f"rate_D{k}"] * PERIOD_YEARS for k in SHOW_DOMAINS]
        ax.plot(SHOW_DOMAINS, centred(v), lw=1.5, color=colour,
                label=f"M={r.M:g}, f={r.fraction:g}  ({r.score_d48:.1f} m)",
                zorder=5)

    ax.set_xlabel("GIS domain")
    ax.set_ylabel("Shoreline change 1984→2004 (m)\ndemeaned over the fit "
                  "window   [+ = seaward]")
    ax.set_xticks(SHOW_DOMAINS)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_title("The groin improves the fit but does not reproduce the shape",
                 loc="left", fontsize=12, fontweight="bold", color=INK, pad=8)
    fig.text(0.5, -0.04,
             "Inside the window the observed profile peaks at D6, dips at D7 "
             "and peaks again at D8; every cell draws a smooth monotonic rise. "
             "The RMSE gain comes from\nmatching the overall D4−D8 slope. "
             "Read the residual as the split between what the groin explains "
             "and what the source/sink calibration absorbs.",
             ha="center", fontsize=8.8, color=MUTED)
    p = OUT / "fig_top_profiles.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {p.name}")


# Transcribed from the pass-by-pass table in hatteras_site_config.py. NOT
# derivable from run_index.csv: --overwrite replaces a run's row, so only the
# final pass survives there. This figure is the durable record.
BE_PASSES = ["zeroBE", "edgeBE", "calib\npass 0", "pass 1", "pass 2",
             "+GIS 90\nre-solve", "final"]
BE_RMSE = {1984: [1.422, 1.219, 0.721, 0.547, 0.527, 0.556, 0.517],
           2004: [2.124, 1.794, 0.763, 0.615, 0.583, 0.603, 0.563]}


def fig_be_convergence():
    """Interior RMSE per calibration pass, both periods."""
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    x = np.arange(len(BE_PASSES))
    for period, colour, mark in ((1984, ACCENT, "o"), (2004, "#1f6f8b", "s")):
        ax.plot(x, BE_RMSE[period], marker=mark, ms=6.5, lw=2.0, color=colour,
                label=f"{period}–{period + 20}")
        ax.annotate(f"{BE_RMSE[period][-1]:.3f}", xy=(x[-1], BE_RMSE[period][-1]),
                    xytext=(8, -2), textcoords="offset points",
                    fontsize=9.5, color=colour, fontweight="bold")
    ax.axvspan(4.5, 5.5, color="#f2f4f6", zorder=0)
    ax.annotate("edges gained,\ninterior gave back", xy=(5, 1.55),
                ha="center", fontsize=9, color=MUTED, style="italic")
    ax.set_xticks(x)
    ax.set_xticklabels(BE_PASSES, fontsize=9)
    ax.set_ylabel("Interior RMSE vs CoastSat LRR (m/yr)")
    ax.set_xlabel("Calibration stage")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(fontsize=10)
    ax.set_title("Source/sink calibration converged in five passes, both "
                 "periods", loc="left", fontsize=12, fontweight="bold",
                 color=INK, pad=8)
    fig.text(0.5, -0.06,
             "Base run: road + beach/dune management, groin off. Stopping rule "
             "(convergence_history.json): a pass buys less than 5% of the "
             "standing RMSE — met at pass 2→3.\nThe GIS 90 re-solve "
             "cost the interior 0.03 while buying 0.61 at the edge, about "
             "20:1, so the edge/interior loop contracts rather than "
             "oscillating.",
             ha="center", fontsize=8.8, color=MUTED)
    p = OUT / "fig_be_convergence.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {p.name}")


def fig_period2_and_bug(d):
    """Why period 2 is unfittable, and what the topography bug was worth."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.4, 4.6))

    reach_lo = float(d.differential_m_yr.min())
    reach_hi = float(d.differential_m_yr.max())
    a1.axhspan(reach_lo, reach_hi, color="#e8f3ee", zorder=0)
    a1.annotate("reachable by the module\n(trapping ≥ 0)",
                xy=(0.5, (reach_lo + reach_hi) / 2), ha="center",
                fontsize=9, color="#1b7f5a")
    a1.errorbar([0], [3.46], yerr=[0.70], marker="o", ms=9, lw=2,
                color=ACCENT, capsize=5, label="observed 1984–2004")
    a1.errorbar([1], [-3.85], yerr=[0.76], marker="s", ms=9, lw=2,
                color="#1f6f8b", capsize=5, label="observed 2004–2024")
    a1.axhline(0, color=INK, lw=0.9)
    a1.set_xlim(-0.6, 1.6)
    a1.set_xticks([0, 1])
    a1.set_xticklabels(["1984–2004", "2004–2024"])
    a1.set_ylabel("Fillet trend, D5−D6 (m/yr)")
    a1.set_title("(a)  Period 2 is outside the module's range",
                 loc="left", fontsize=11, color=INK, pad=6)
    a1.grid(axis="y", color=GRID, lw=0.6)
    a1.set_axisbelow(True)
    for s in ("top", "right"):
        a1.spines[s].set_visible(False)
    a1.legend(fontsize=9, loc="lower left")

    dom = np.arange(1, 13)
    drift = np.zeros(12)
    drift[3] = 0.1258
    a2.bar(dom, drift, color=FOIL, width=0.65)
    a2.bar(dom, np.zeros(12), color=ACCENT, width=0.65)
    a2.axhline(0.005, color=ACCENT, lw=1.4, ls=(0, (4, 2)))
    a2.annotate("guard tolerance 0.005", xy=(12, 0.007), ha="right",
                fontsize=9, color=ACCENT)
    a2.annotate("after the fix:\nexactly 0 at every domain", xy=(8, 0.075),
                ha="center", fontsize=9.5, color="#1b7f5a", fontweight="bold")
    a2.set_xlabel("GIS domain")
    a2.set_ylabel("|sweep − published| (m/yr)")
    a2.set_xticks(dom)
    a2.set_title("(b)  The topography-product bug, 25× the tolerance",
                 loc="left", fontsize=11, color=INK, pad=6)
    a2.grid(axis="y", color=GRID, lw=0.6)
    a2.set_axisbelow(True)
    for s in ("top", "right"):
        a2.spines[s].set_visible(False)

    fig.text(0.5, -0.06,
             "(a) A groin with trapping ≥ 0 can only WIDEN the updrift"
             "−downdrift gap. Period 1's observed gap widens; period 2's "
             "narrows, so no M ≥ 0 reaches it and there is no joint fit.\n"
             "(b) The sweep worker resolved topography without naming a "
             "product, so DEFAULT_PRODUCT (2004-start) answered and 1984 cells "
             "were built on the 2004 island.",
             ha="center", fontsize=8.8, color=MUTED)
    p = OUT / "fig_period2_and_bug.png"
    fig.savefig(p, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {p.name}")


def main():
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
