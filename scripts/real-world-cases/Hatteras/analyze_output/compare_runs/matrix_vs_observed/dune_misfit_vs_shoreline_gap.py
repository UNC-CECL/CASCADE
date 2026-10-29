#!/usr/bin/env python3
"""
Where the full-management runs' dunes miss the lidar, against where their shorelines miss CoastSat.

    python scripts/real-world-cases/Hatteras/analyze_output/compare_runs/matrix_vs_observed/dune_misfit_vs_shoreline_gap.py

One panel per window (1996-2009, 2009-2025), one point per interior GIS domain: the
mean model minus survey median dune crest over the surveys after the start (from the
run's dune_crest_vs_lidar.csv) against the model minus CoastSat shoreline position
change over the window (both smoothed as the v2 dune crest overview draws them).
Writes output/raw_runs/matrix/figures/dune_misfit_vs_shoreline_gap_full_management.png.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-28
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

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
for _path in (_HAT_DIR, _HAT_DIR / "cascade_pipeline" / "plotting"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
from hatteras_paths import WORKSPACE as _REPO  # noqa: E402

import dune_crest_comparison as dc  # noqa: E402
from cascade_pipeline.run_layout import TABLES  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    INK, INK_MUTED, _title, apply_style, open_frame, record_caption, save)
from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, SCORE_INTERIOR_GIS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = _REPO / "output" / "raw_runs" / "matrix"
OUT = RAW / "figures"
RUNS = {(1996, 2009): "1996_2009/edgeBE/HAT_1996_2009_edgeBE_offsetmetres_road_bdm_groin",
        (2009, 2025): "2009_2025/edgeBE/HAT_2009_2025_edgeBE_offsetmetres_road_bdm_nourish_groin"}
# A domain's dunes count as too tall or too low beyond the overview's white band
MATCH_M = 0.5
TALL_FILL = "#F4A582"    # RdBu light red: dunes too tall and shoreline seaward, the model short of erosion
LOW_FILL = "#92C5DE"     # RdBu light blue: dunes too low and shoreline landward, the model over-eroding
# -----------------------------------------------------------------------------


# Per interior domain: mean dune misfit after the start survey, and the model minus CoastSat shoreline change
def domain_pairs(run_dir, window):
    table = pd.read_csv(run_dir / TABLES / dc.TABLE)
    misfit = table[table.survey > window[0]].groupby("gis").model_minus_survey_m.mean()
    model, observed = dc.position_change(run_dir, window)
    lo, hi = SCORE_INTERIOR_GIS
    return pd.DataFrame({"misfit": misfit, "gap": model - observed}).loc[lo:hi].dropna()


def in_village(gis):
    return any(lo <= gis <= hi for lo, hi in HATTERAS_ANNOTATIONS.town_spans.values())


def main():
    apply_style()
    pairs = {w: domain_pairs(RAW / r, w) for w, r in RUNS.items()}
    xmax = max(np.abs(p.misfit).max() for p in pairs.values()) * 1.08
    ymax = max(np.abs(p.gap).max() for p in pairs.values()) * 1.08
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), sharex=True, sharey=True, layout="constrained")
    stats = []
    for i, (ax, (window, df)) in enumerate(zip(axes, pairs.items())):
        # The two quadrants the misfits would share if the dunes and the shoreline erred together
        ax.fill_between([MATCH_M, xmax], 0, ymax, color=TALL_FILL, alpha=0.35, lw=0, zorder=0)
        ax.fill_between([-xmax, -MATCH_M], -ymax, 0, color=LOW_FILL, alpha=0.35, lw=0, zorder=0)
        ax.text(xmax * 0.97, ymax * 0.95, "dunes too tall,\nshoreline too seaward", ha="right", va="top",
                fontsize=8.5, color=INK_MUTED)
        ax.text(-xmax * 0.97, -ymax * 0.95, "dunes too low,\nshoreline too landward", ha="left",
                va="bottom", fontsize=8.5, color=INK_MUTED)
        ax.axhline(0, color=INK_MUTED, lw=0.8, zorder=1)
        ax.axvspan(-MATCH_M, MATCH_M, color="0.93", lw=0, zorder=0)

        village = df.index.map(in_village)
        ax.scatter(df.misfit[~village], df.gap[~village], s=26, fc="white", ec=INK, lw=0.9, zorder=3)
        ax.scatter(df.misfit[village], df.gap[village], s=26, fc=INK, ec=INK, lw=0.9, zorder=3)
        slope, icept = np.polyfit(df.misfit, df.gap, 1)
        xs = np.array([df.misfit.min(), df.misfit.max()])
        ax.plot(xs, icept + slope * xs, color=INK, lw=1.0, ls="--", zorder=2)
        # The domains outside the white band, by GIS number
        for gis, row in df[np.abs(df.misfit) > MATCH_M].iterrows():
            ax.annotate(str(gis), (row.misfit, row.gap), xytext=(3, 3), textcoords="offset points",
                        fontsize=6.5, color=INK_MUTED, zorder=4)

        tall, low = df[df.misfit > MATCH_M], df[df.misfit < -MATCH_M]
        r = df.misfit.corr(df.gap)
        rho = df.misfit.corr(df.gap, method="spearman")
        ax.text(0.02, 0.98,
                f"Pearson r = {r:+.2f}, Spearman ρ = {rho:+.2f}, {len(df)} domains\n"
                f"too tall: {len(tall)} domains, {(tall.gap > 0).sum()} seaward, mean gap {tall.gap.mean():+.1f} m\n"
                f"too low: {len(low)} domains, {(low.gap < 0).sum()} landward, mean gap {low.gap.mean():+.1f} m",
                transform=ax.transAxes, ha="left", va="top", fontsize=8,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8"))
        _title(ax, i, f"{window[0]}–{window[1]}, full management")
        ax.set_xlabel("Mean model − survey median dune crest (m)\n← model lower     model higher →")
        open_frame(ax)
        stats.append((window, r, rho, len(tall), (tall.gap > 0).sum(), tall.gap.mean(),
                      len(low), (low.gap < 0).sum(), low.gap.mean()))
    axes[0].set_ylabel("Model − CoastSat shoreline position change (m)\n← model landward     model seaward →")
    axes[0].set_xlim(-xmax, xmax)
    axes[0].set_ylim(-ymax, ymax)
    fig.legend(handles=[Line2D([], [], ls="", marker="o", mfc=INK, mec=INK, label="village domain"),
                        Line2D([], [], ls="", marker="o", mfc="white", mec=INK, label="outside the villages"),
                        Line2D([], [], color=INK, lw=1.0, ls="--", label="least-squares fit")],
               loc="outside lower center", ncol=3, frameon=False)

    png = OUT / "dune_misfit_vs_shoreline_gap_full_management.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Dune crest misfit against shoreline position misfit, per interior GIS domain (GIS "
        f"{SCORE_INTERIOR_GIS[0]}–{SCORE_INTERIOR_GIS[1]}), for the edgeBE full-management matrix runs: (a) "
        f"{RUNS[(1996, 2009)].split('/')[-1]}, (b) {RUNS[(2009, 2025)].split('/')[-1]}. x: the mean of model "
        "minus lidar median dune crest over the surveys after the start (the start survey is the model's "
        "starting topography). y: the model's shoreline position change over the window minus CoastSat's, "
        "both 7-domain LOWESS, seaward positive, as drawn over the v2 dune crest overview. Grey band: "
        f"within ±{MATCH_M} m, the overview's white band; GIS numbers mark the domains outside it. Red "
        "quadrant: dunes too tall and shoreline too seaward (the model short of erosion); blue: dunes too "
        "low and shoreline too landward (the model over-eroding). Filled points are village domains. "
        "Neighbouring domains are not independent (the shoreline is smoothed over seven), so r overstates "
        "the evidence."))
    for window, r, rho, nt, ts, tg, nl, ll, lg in stats:
        print(f"  {window[0]}-{window[1]}  r {r:+.2f}  rho {rho:+.2f}  tall {nt} ({ts} seaward, {tg:+.1f} m)"
              f"  low {nl} ({ll} landward, {lg:+.1f} m)")
    print(f"  -> {png}")


if __name__ == "__main__":
    main()
