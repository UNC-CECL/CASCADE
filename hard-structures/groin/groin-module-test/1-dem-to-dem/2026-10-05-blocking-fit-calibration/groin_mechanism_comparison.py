#!/usr/bin/env python3
"""
The old dipole groin against the pinned blocking groin, on the same DEM-to-DEM setup.

    python groin_mechanism_comparison.py   ->  figures/groin_dipole_vs_blocking.png + tables/

What each groin applies per year at its two flanks, the GIS 5-6 gap against the
photos, and the net change along GIS 1-20 against CoastSat, for 1996-2009 and
2009-2025. Every run: full management, the solved edgeBE ends, relocations off,
failure instant from the 2004 step.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
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

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
sys.path[:0] = [str(HERE.parents[1] / "0-solver-audit" / "2026-09-29-option-a-real-planform"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]
from score_instant_grid import observed_series  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style, figsize, open_frame,
    record_caption, save)
from site_layer.hat_observed_rates import net_change_domain_csv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW = REPO / "output" / "raw_runs"
STUDY = RAW / "experiments" / "groin" / "2026-10-05-blocking-fit-dem-to-dem"
PERIODS = ((1996, 2009), (2009, 2025))
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1          # GIS 6 / GIS 5 in the padded array
PAD = 15                                   # padded index of GIS 1 is PAD
ZOOM = (1, 20)
PROFILE_M, DOMAIN_M = 18.77, 500.0         # active profile height (runner) and domain length
# label, member folder (None = matrix no-groin), colour, line style
SERIES = (
    ("no groin", None, C["BASE"], "-"),
    ("dipole M 60, f 0.6 (pinned 2026-08-30)", "dipole_M60_f0.6", C["ADDED"], "--"),
    ("dipole M 12, f 0.3 (best dipole, 2026-09-29)", "dipole_M12_f0.3", C["ADDED"], ":"),
    ("blocking b 0.6, f 0.6 (pinned 2026-10-05)", "b0.60_f0.6", C["ACCENT"], "-"),
)
LW = 1.4
# -----------------------------------------------------------------------------


def window(p):
    return f"{p[0]}_{p[1]}"


def run_dir(period, member):
    if member is None:
        fill = "_nourish" if period[0] == 2009 else ""
        return RAW / "matrix" / window(period) / "edgeBE" / (
            f"HAT_{window(period)}_edgeBE_offsetmetres_road_bdm{fill}_nogroin")
    hits = sorted((STUDY / member / window(period) / "edgeBE").glob("*/"))
    return hits[-1] if hits else None


def matrix(rd):
    return np.load(next(rd.glob("*_shoreline_matrix.npy")))


# Gap change since the period start, + = GIS 6 gains on GIS 5 (x grows landward)
def gap(x):
    g = x[:, DOWN] - x[:, UP]
    return g - g[0]


# Seaward shoreline change the groin applied each year at GIS 6 (+) and GIS 5 (-)
def applied(rd):
    f = rd / "tables" / "groin_diagnostics.csv"
    if not f.is_file():
        return None
    d = pd.read_csv(f)
    return d["model_year"].to_numpy(), -d["applied_dx_updrift_m"].to_numpy(), \
        -d["applied_dx_downdrift_m"].to_numpy()


def obs_gap_changes(obs, period):
    g0 = np.interp(period[0], obs.index, obs.values)
    yrs = [y for y in obs.index if period[0] < y <= period[1]]
    return np.array(yrs), np.array([obs[y] for y in yrs]) - g0


def main():
    apply_style()
    obs = observed_series()
    fig, axes = plt.subplots(4, 1, figsize=figsize("double", height=9.2), constrained_layout=True,
                             gridspec_kw=dict(height_ratios=(1, 1.1, 1, 1)))
    rows = []

    # (a) what each groin does per year
    ax = axes[0]
    for label, member, col, ls in SERIES[1:]:
        for period in PERIODS:
            rd = run_dir(period, member)
            a = applied(rd) if rd else None
            if a is None:
                continue
            yr, up, down = a
            ax.plot(yr, up, color=col, ls=ls, lw=LW)
            ax.plot(yr, down, color=col, ls=ls, lw=LW)
            vol = np.mean(np.abs(up)) * PROFILE_M * DOMAIN_M
            rows.append(dict(groin=label, period=window(period),
                             mean_gain_gis6_m_yr=float(np.mean(up)),
                             mean_change_gis5_m_yr=float(np.mean(down)),
                             trapped_m3_yr_at_gis6=float(vol)))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.axvline(2004, color=INK_MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(2004.2, ax.get_ylim()[1] * 0.92, "2004: failure", fontsize=6.5, color=INK_MUTED,
            va="top")
    ax.axvline(2009, color=INK, lw=0.6)
    ax.set_xlim(1996, 2025)
    ax.set_ylabel("Applied change\n(m/yr, seaward +)")
    ax.text(1996.3, ax.get_ylim()[1] * 0.75, "GIS 6 (north of groin) held seaward", fontsize=6.5,
            color=INK_MUTED)
    ax.text(1996.3, ax.get_ylim()[0] * 0.8, "GIS 5 (south) held landward", fontsize=6.5,
            color=INK_MUTED)
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 0, "What the groin adds to each flank every model year")

    # (b) the gap at the groin against the photos
    ax = axes[1]
    for label, member, col, ls in SERIES:
        for period in PERIODS:
            rd = run_dir(period, member)
            if rd is None:
                continue
            g = gap(matrix(rd))
            ax.plot(period[0] + np.arange(g.size), g, color=col, ls=ls, lw=LW)
    for period in PERIODS:
        yrs, o = obs_gap_changes(obs, period)
        ax.plot(yrs, o, "o", ms=4.5, mfc=INK, mec="white", mew=0.6, zorder=6)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.axvline(2009, color=INK, lw=0.6)
    ax.axvline(2017.5, color=INK_MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(2017.6, ax.get_ylim()[1] * 0.9, "Buxton fill", fontsize=6.5, color=INK_MUTED, va="top")
    ax.set_xlim(1996, 2025)
    ax.set_ylabel("GIS 5-6 gap change\nfrom period start (m)")
    ax.set_xlabel("Year (calibration run 1996-2008 | test run 2009-2024, each from its own start)")
    ax.grid(axis="y")
    open_frame(ax)
    _title(ax, 1, "The step at the groin against the photos (dots)")

    # (c, d) net change along GIS 1-20
    for i, period in enumerate(PERIODS):
        ax = axes[2 + i]
        gis = np.arange(ZOOM[0], ZOOM[1] + 1)
        o = pd.read_csv(net_change_domain_csv(*period), index_col=0)["net_change_m"]
        for label, member, col, ls in SERIES:
            rd = run_dir(period, member)
            if rd is None:
                continue
            x = matrix(rd)
            ch = -(x[-1] - x[0])[PAD + gis - 1]
            ax.plot(gis, ch, color=col, ls=ls, lw=LW)
            rows.append(dict(groin=label, period=window(period), net_gis5_m=float(ch[4]),
                             net_gis6_m=float(ch[5]), obs_gis5_m=float(o[5]), obs_gis6_m=float(o[6])))
        ax.plot(gis, o.loc[gis].values, "o", ms=4, mfc=INK, mec="white", mew=0.6, zorder=6)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.axvline(5.5, color=C["GROIN"], lw=0.8)
        ax.text(5.6, ax.get_ylim()[0] * 0.9, "groin", fontsize=6.5, color=C["GROIN"], va="bottom")
        ax.set_xlim(*ZOOM)
        ax.set_xticks(gis)
        ax.set_ylabel("Net change (m,\nseaward +)")
        ax.grid(axis="y")
        open_frame(ax)
        role = "Calibration" if i == 0 else "Test"
        _title(ax, 2 + i, f"{role} {period[0]}-{period[1]}: end minus start, CoastSat (dots)")
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)

    handles = [Line2D([], [], color=c, ls=ls, lw=LW, label=lab) for lab, _, c, ls in SERIES]
    handles.append(Line2D([], [], ls="", marker="o", ms=4.5, mfc=INK, mec="white", label="observed"))
    fig.legend(handles=handles, loc="outside upper center", ncol=3, fontsize=6.8, frameon=False)

    (HERE / "tables").mkdir(exist_ok=True)
    pd.DataFrame(rows).round(3).to_csv(HERE / "tables" / "groin_dipole_vs_blocking.csv", index=False)
    png = HERE / "figures" / "groin_dipole_vs_blocking.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "The old dipole groin against the pinned blocking groin on the same DEM-to-DEM setup: full "
        "management, the solved edgeBE ends (+1.4981 / +10.5659 m/yr), relocations off, failure "
        "instant from the 2004 step. Calibration runs 1996-2008, test runs 2009-2024. "
        "(a) The shoreline change each groin adds every model year, seaward positive: above zero at "
        "GIS 6, below at GIS 5. The dipole adds a fixed +M / -M (times f after 2004) whatever the "
        "shoreline does. The blocking groin cancels a fraction b of the change BRIE's alongshore step "
        "would move across the GIS 5|6 face that year, so its size follows the shoreline and BRIE's "
        "local diffusivity. (b) The gap between GIS 5 and GIS 6, as change from each period's start; "
        "dots are the wet/dry photo gaps from the same start. (c, d) End minus start shoreline along "
        "GIS 1-20 against the CoastSat net-change target (domain means). The red line is the groin."))
    print(f"wrote {png.relative_to(REPO)}")
    print(pd.DataFrame(rows).round(2).to_string())


if __name__ == "__main__":
    main()
