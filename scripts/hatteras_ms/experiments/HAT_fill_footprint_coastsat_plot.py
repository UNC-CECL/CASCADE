"""
Net change 2009-2025 with the fills on their reported and on their CoastSat-observed footprints.

    python scripts/hatteras_ms/experiments/HAT_fill_footprint_coastsat_plot.py [--positions]

Reads tables/per_domain_net_change.csv and scores.csv written by HAT_fill_footprint_coastsat.py;
draws the alongshore net change for both runs against the test target, and the
difference between the runs with each fill's two footprints marked. --positions draws the
shoreline positions instead: the 2009 start, observed 2025 and both runs' 2025, island-wide
and zoomed on the three fill reaches.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-06
"""
from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path[:0] = [str(PROJECT_ROOT / "scripts"), str(_HERE.parent),
                str(PROJECT_ROOT / "scripts" / "input_prep" / "5-scr" / "3-rates" / "coastsat" / "net_change")]
import HAT_fill_footprint_coastsat as X  # noqa: E402
from coastsat_net_change import smooth_like_model  # noqa: E402
from site_layer.hat_figure_style import (C, C_1984, C_1997, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,  # noqa: E402
                                         figsize, mark_offaxis, offaxis_clause, open_frame, record_caption,
                                         save, structures, town_bands)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = X.EXP_DIR / "figures" / "fill_footprint_reported_vs_coastsat_net_change_2009_2025.png"
OUT_POS = X.EXP_DIR / "figures" / "fill_footprint_reported_vs_coastsat_positions_2009_2025.png"
ZERO_GIS = 76                   # the offset build's zero domain, position 0 here
REACHES = (("Buxton", 1, 20), ("Avon", 17, 32), ("Rodanthe", 76, 90))
Y_HALF = 100.0                  # the fixed net-change axis every run figure uses
XLIM = (0.5, 90.5)
# -----------------------------------------------------------------------------


# Cross-shore position of the model start, seaward +, relative to ZERO_GIS
def start_position():
    D = HATTERAS_DOMAINS
    m0 = np.load(next(X.MATRIX_RUN.glob("*_shoreline_matrix.npy")))[0][D.start_real_index:D.end_real_index]
    p = pd.Series(-m0, index=pd.RangeIndex(1, 91))
    return p - p.loc[ZERO_GIS]


def positions():
    apply_style()
    per = pd.read_csv(X.EXP_DIR / "tables" / "per_domain_net_change.csv", index_col=0)
    p0 = start_position()
    lines = [("2009 start (model input)", p0, C_1984, "-", 1.4, None),
             ("Observed 2025 (CoastSat)", p0 + per.observed_m, C_1997, "-", 1.6, "o"),
             ("Model 2025, reported footprints", p0 + per.model_reported_m, C["BASE"], (0, (4, 2)), 1.3, None),
             ("Model 2025, CoastSat footprints", p0 + per.model_coastsat_m, C["ACCENT"], "-", 1.4, None)]
    observed = X.observed_ranges()

    fig, axes = plt.subplots(3, 1, figsize=figsize("double", height=7.6), constrained_layout=True)
    for i, (ax, (name, lo, hi)) in enumerate(zip(axes, REACHES)):
        g = np.arange(lo, hi + 1)
        base = np.polyval(np.polyfit(g, p0.loc[lo:hi].values, 1), g)
        for label, y, col, ls, lw, mk in lines:
            ax.plot(g, y.loc[lo:hi].values - base, color=col, ls=ls, lw=lw, marker=mk, ms=2.5,
                    zorder=4 if col == C["ACCENT"] else 3)
        ax.axhline(0, color=INK_MUTED, lw=0.5, zorder=1)
        y0, y1 = ax.get_ylim()
        span = y1 - y0
        ax.set_ylim(y0 - 0.2 * span, y1)
        for p in HATTERAS_NOURISHMENT_PROJECTS:
            for gis, frac, col in ((p.gis_domains, 0.11, C["BASE"]), (observed[(p.name, p.year)], 0.04, C["ACCENT"])):
                if lo <= min(gis) <= hi:
                    ax.plot([min(gis) - 0.45, max(gis) + 0.45], [y0 - 0.2 * span + frac * 1.2 * span] * 2,
                            color=col, lw=3, solid_capstyle="butt", zorder=5)
        ax.set_xlim(lo - 0.5, hi + 0.5)
        ax.set_xticks(g)
        ax.set_ylabel("Position from reach baseline (m)")
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, i, f"{name}, GIS {lo}–{hi}")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)

    handles = [Line2D([], [], color=col, ls=ls, lw=lw, marker=mk, ms=3, label=label)
               for label, _, col, ls, lw, mk in lines]
    handles += [Patch(color=C["BASE"], label="Reported footprint"), Patch(color=C["ACCENT"], label="CoastSat footprint")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS as A
    for ax, (_, lo, hi) in zip(axes, REACHES):
        inside = lambda pos: lo - 0.5 < pos < hi + 0.5  # noqa: E731
        structures(ax, spans=dataclasses.replace(
            A, groins={k: v for k, v in A.groins.items() if inside(v)},
            piers={k: v for k, v in A.piers.items() if inside(v[0])}))
        ax.set_xlim(lo - 0.5, hi + 0.5)
    save(fig, OUT_POS, close=True)
    record_caption(OUT_POS, (
        "Cross-shore shoreline position per 500 m domain in the three fill reaches, positive seaward, measured "
        "from a straight baseline fitted to the 2009 start over each reach (the island's planform curvature, "
        "kilometres across the island, would otherwise hide differences of tens of metres). The 2009 start is the "
        "model's initial shoreline, built from the CoastSat mean over 2008-08-17 to 2010-08-17. Observed 2025 is "
        "that start plus the observed net change to the 2025-08-17 ±6 month CoastSat mean (unsmoothed domain "
        "means). Model 2025 is the 1 Jan 2025 state of the 2009–2025 test run (per-domain source/sink set 1, "
        "blocking groin, full management, relocations off) with the fills on their reported footprints (dashed "
        "grey) and on the CoastSat-observed footprints (purple); volumes are the same in both. Where the two runs "
        "agree the purple line hides the dashed one. Bars at the bottom mark each fill's reported (grey) and "
        "CoastSat (purple) footprint; Buxton's two fills share GIS 6–16 reported, with CoastSat ranges 8–13 "
        "(2017) and 7–9 (2022). Each panel has its own y range."))
    print(OUT_POS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positions", action="store_true")
    if ap.parse_args().positions:
        positions()
        return
    apply_style()
    per = pd.read_csv(X.EXP_DIR / "tables" / "per_domain_net_change.csv", index_col=0)
    sc = pd.read_csv(X.EXP_DIR / "tables" / "scores.csv")
    interior = sc[sc.area.str.startswith("interior")].set_index("variant")
    g = per.index.values
    mod = {v: smooth_like_model(per[f"model_{v}_m"]) for v in X.VARIANTS}
    observed = X.observed_ranges()

    fig, (a, b) = plt.subplots(2, 1, figsize=figsize("double", height=5.6), sharex=True,
                               constrained_layout=True, gridspec_kw=dict(height_ratios=(1.5, 1)))

    # (a) net change, both runs against the target
    a.scatter(g, per.observed_m, s=5, color=C["BASE_FILL"], lw=0, zorder=1)
    a.plot(g, per.observed_lowess7_m, color=INK, lw=1.6, zorder=4)
    a.plot(g, mod["reported"], color=C["BASE"], lw=1.3, ls=(0, (4, 2)), zorder=3)
    a.plot(g, mod["coastsat"], color=C["ACCENT"], lw=1.4, zorder=5)
    a.axhline(0, color=INK_MUTED, lw=0.5, zorder=2)
    a.set_ylim(-Y_HALF, Y_HALF)
    off = [("the observed change", mark_offaxis(a, g, per.observed_lowess7_m, Y_HALF)),
           ("the model, reported footprints", mark_offaxis(a, g, mod["reported"], Y_HALF, C["BASE"])),
           ("the model, CoastSat footprints", mark_offaxis(a, g, mod["coastsat"], Y_HALF, C["ACCENT"]))]
    a.set_ylabel("Net shoreline change (m)")
    a.grid(axis="y")
    open_frame(a)
    _title(a, 0, "Net shoreline change, 2009–2025")

    # (b) what moving the footprints changed, raw per domain, with both footprints marked
    diff = per.coastsat_minus_reported_m
    b.bar(g, diff, width=0.8, color=C["ADDED"], lw=0, zorder=3)
    b.axhline(0, color=INK, lw=0.6, zorder=4)
    lim = max(10.0, float(diff.abs().max()) * 1.15)
    b.set_ylim(-lim, lim * 1.35)
    top = lim * 1.35
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        for gis, y, col in ((p.gis_domains, top * 0.90, C["BASE"]),
                            (observed[(p.name, p.year)], top * 0.78, C["ACCENT"])):
            b.plot([min(gis) - 0.45, max(gis) + 0.45], [y, y], color=col, lw=3, solid_capstyle="butt", zorder=5)
    b.set_ylabel("CoastSat − reported (m)")
    b.set_xlabel(DOMAIN_AXIS_LABEL)
    b.set_xlim(*XLIM)
    b.set_xticks(range(1, 91, 5))
    b.grid(axis="y")
    open_frame(b)
    _title(b, 1, "Change made by moving the fills, unsmoothed")
    town_bands(a, strip=0.06)
    town_bands(b, strip=0.06, label=False)
    structures(a)

    handles = [Line2D([], [], color=INK, lw=1.6, label="Observed (CoastSat, smoothed)"),
               Line2D([], [], color=C["BASE_FILL"], marker="o", ms=3, ls="none", label="Observed, per domain"),
               Line2D([], [], color=C["BASE"], lw=1.3, ls=(0, (4, 2)),
                      label=f"Model, reported footprints (RMSE {interior.loc['reported', 'rmse_m']:.1f} m)"),
               Line2D([], [], color=C["ACCENT"], lw=1.4,
                      label=f"Model, CoastSat footprints (RMSE {interior.loc['coastsat', 'rmse_m']:.1f} m)"),
               Patch(color=C["BASE"], label="Reported footprint"),
               Patch(color=C["ACCENT"], label="CoastSat footprint")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    save(fig, OUT, close=True)

    fills = "; ".join(f"{p.name} {p.year} GIS {min(p.gis_domains)}–{max(p.gis_domains)} reported, "
                      f"{min(observed[(p.name, p.year)])}–{max(observed[(p.name, p.year)])} CoastSat"
                      for p in HATTERAS_NOURISHMENT_PROJECTS)
    s = interior
    record_caption(OUT, (
        "2009–2025 test run (per-domain source/sink set 1, blocking groin, full management, relocations off), "
        "run twice: with each fill on its reported footprint, and on the range where CoastSat saw the "
        "shoreline move after it (nourishment/4-extent-checks/coastsat). Volumes are the same in both, so a narrower "
        f"footprint carries more sand per metre. Footprints: {fills}. "
        "(a) Net shoreline change, positive seaward. The model is the 1 Jan 2025 state minus the start, "
        "the target the 2025-08-17 ±6 month CoastSat mean minus the 2008-08-17 to 2010-08-17 mean. Lines pass "
        "through a 7-domain LOWESS with GIS 1–10 raw; grey points are the unsmoothed observed domain values. "
        "Away from the fills the two runs are identical and the purple line hides the dashed one. "
        "Legend RMSE is over the interior, GIS 2–89. Bias is "
        f"{s.loc['reported', 'bias_m']:+.1f} m (reported) and {s.loc['coastsat', 'bias_m']:+.1f} m (CoastSat); "
        f"r is {s.loc['reported', 'r']:.2f} and {s.loc['coastsat', 'r']:.2f}. "
        "(b) The CoastSat-footprint run minus the reported-footprint run, per domain, unsmoothed (positive: more seaward with the CoastSat footprints). "
        "Bars at the top mark each fill's reported (grey) and CoastSat (purple) footprint."
        + offaxis_clause(off, Y_HALF)))
    print(OUT)


if __name__ == "__main__":
    main()
