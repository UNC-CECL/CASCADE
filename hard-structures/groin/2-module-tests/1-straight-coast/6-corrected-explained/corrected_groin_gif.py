#!/usr/bin/env python3
"""
How the corrected groin works, year by year: one groin on a straight coast at -20 degrees.

    python corrected_groin_gif.py [--years 60] [--theta -20]   ->  figures/corrected_groin_how_it_works.gif

Reuses straight_coast_test.py's emulator and module 6, unchanged. Top: the shoreline either side of
the groin, the drift arriving, the sand the groin holds back that year (moved from the downdrift
cell to the updrift cell), and the sand BRIE's diffusion carries back across the groin face, which
BRIE treats as open. Bottom: the two flows by year, and the gap across the groin, which slows as
the return flow catches up, against the 25 m real-sized field with 110 m groins
(../../6-long-term/tables/, if present).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, figsize, open_frame  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = HERE / "figures" / "corrected_groin_how_it_works.gif"
REFERENCE = HERE.parents[1] / "6-long-term" / "tables" / "long_term_by_year.csv"
FPS = 6
HOLD_FRAMES = 12                          # repeats of the last frame before the loop restarts
WINDOW = 10                               # domains either side of the groin drawn
# -----------------------------------------------------------------------------


# The sand the groin holds back in the coming year, m of shoreline, from the shoreline x (seaward +)
def held_back(t, theta0, x):
    eta = -x
    step = np.tan(np.deg2rad(theta0)) * sc.DY_M + (np.roll(eta, -1) - eta)
    approach = np.degrees(np.arctan2(step[sc.HI], sc.DY_M))       # the face just updrift of the groin
    q = abs(float(sc.drift(t, approach)[0]))
    return sc.BLOCKING_FRACTION * q * sc.DT_YR / (sc.DY_M * t["depth"]), q


# The sand BRIE's diffusion carries back across the groin face in the coming year, m of shoreline (approximate:
# the face's r times the step across it, for the two Crank-Nicolson halves)
def carried_back(t, theta0, x):
    eta = -x
    step = np.tan(np.deg2rad(theta0)) * sc.DY_M + (np.roll(eta, -1) - eta)
    theta = np.degrees(np.arctan2(step, sc.DY_M))
    r = max(0.0, float(t["coast_diff"][sc.angle_index(np.array([theta[sc.LO]]), t["climl"])][0]) * sc.DT_YR / 2 / sc.DY_M ** 2)
    return 2 * r * (x[sc.HI] - x[sc.LO])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=60)
    ap.add_argument("--theta", type=int, default=-20)
    a = ap.parse_args()
    apply_style()
    t = sc.brie_tables()
    th = a.theta
    assert float(sc.drift(t, th)[0]) < 0, "drawn for drift toward lower index (updrift on the right)"
    sc.PROFILE_YEARS = tuple(range(1, a.years + 1))
    _, prof = sc.run(t, th, "corrected", a.years)
    shore = {0: np.zeros(sc.NY), **prof}                          # seaward +, updrift at higher index
    held, drift_in, gap, back = [], [], [], []
    for y in range(a.years + 1):
        h, q = held_back(t, th, shore[y])
        held.append(h)
        back.append(carried_back(t, th, shore[y]))
        drift_in.append(q)
        gap.append(shore[y][sc.HI] - shore[y][sc.LO])
    ref = None
    if REFERENCE.exists():
        r = pd.read_csv(REFERENCE)
        ref = r[(r.source == "reference") & (r.theta0_deg == th) & (r.year <= a.years)]

    dom = np.arange(sc.NY) - sc.LO - 0.5                          # domain centres from the groin face
    near = np.abs(dom) <= WINDOW + 1
    lim = max(np.abs(shore[a.years][near]).max(), 10) * 1.25
    fig = plt.figure(figsize=figsize("double", height=5.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.35, 1))

    def draw(year):
        fig.clear()
        gs = fig.add_gridspec(2, 2, height_ratios=(1.35, 1))
        ax = fig.add_subplot(gs[0, :])
        x = shore[year]
        ax.fill_between(dom, 0, x, where=x >= 0, step="mid", color=C["ADDED_FILL"], lw=0)
        ax.fill_between(dom, 0, x, where=x < 0, step="mid", color=C["WATER"], alpha=0.6, lw=0)
        ax.step(dom, x, where="mid", color=INK, lw=1.4)
        ax.axhline(0, color=INK_MUTED, lw=0.6, ls=(0, (3, 2)))
        ax.axvline(0, color=C["GROIN"], lw=2.2)
        ax.set_xlim(-WINDOW - 0.5, WINDOW + 0.5)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("Domains from the groin (500 m each; updrift right)")
        ax.set_ylabel("Shoreline change (m, seaward +)")
        ax.text(WINDOW * 0.65, lim * 0.35, "built out (updrift)", color=INK_MUTED, fontsize=7, ha="center")
        ax.text(-WINDOW * 0.65, -lim * 0.45, "cut back (downdrift)", color=INK_MUTED, fontsize=7, ha="center")
        # Step 1: the drift arriving from the right
        ax.annotate("", xy=(0.55, 0.92), xytext=(0.97, 0.92), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.4))
        ax.text(0.76, 0.84, f"1. drift arriving: {drift_in[year] / 1e6:.2f} M m$^3$/yr", transform=ax.transAxes,
                ha="center", va="top", fontsize=7.5, color=INK)
        # Step 2: the share held back, moved from the downdrift cell to the updrift cell
        ax.annotate("", xy=(0.6, lim * 0.62), xytext=(-0.6, lim * 0.62),
                    arrowprops=dict(arrowstyle="-|>", color=C["GROIN"], lw=1.3, connectionstyle="arc3,rad=-0.5"))
        ax.text(-WINDOW * 0.5, lim * 0.92, f"2. held back this year: {held[year]:.0f} m\n"
                               f"({sc.BLOCKING_FRACTION:.0%} of the drift, downdrift cell $\\to$ updrift cell)",
                ha="center", va="top", fontsize=7.5, color=C["GROIN"])
        ax.annotate("", xy=(-0.6, -lim * 0.62), xytext=(0.6, -lim * 0.62),
                    arrowprops=dict(arrowstyle="-|>", color=C["LATE"], lw=1.3, connectionstyle="arc3,rad=-0.5"))
        ax.text(WINDOW * 0.45, -lim * 0.4, f"BRIE carries about {back[year]:.0f} m back\nacross the groin, which it treats as open",
                ha="center", va="top", fontsize=7.5, color=C["LATE"])
        ax.text(0.99, 0.03, "3. BRIE's leak near the groin is put back each year,\nso no sand is made or lost",
                transform=ax.transAxes, fontsize=6.8, color=INK_MUTED, va="bottom", ha="right")
        open_frame(ax)
        ax.set_title(f"Year {year}", loc="left", fontsize=9)
        # Bottom left: sand held back each year
        ax = fig.add_subplot(gs[1, 0])
        ax.plot(range(year + 1), held[:year + 1], color=C["GROIN"], lw=1.6, label="held back by the groin")
        ax.plot(range(year + 1), back[:year + 1], color=C["LATE"], lw=1.6, label="carried back by BRIE")
        ax.plot([year, year], [held[year], back[year]], "o", color=INK, ms=3)
        ax.set_xlim(0, a.years)
        ax.set_ylim(0, max(held) * 1.25)
        ax.set_xlabel("Year")
        ax.set_ylabel("Sand per year (m of shoreline)")
        ax.set_title("The groin holds back about the same each year", loc="left", fontsize=7.5)
        ax.legend(loc="lower right", fontsize=6.5, frameon=False)
        ax.grid(axis="y")
        open_frame(ax)
        # Bottom right: the gap across the groin
        ax = fig.add_subplot(gs[1, 1])
        if ref is not None:
            ax.plot(ref.year, ref.gap_m, color=C["BASE"], lw=4, alpha=0.45, label="real-sized field, 110 m groins")
        ax.plot(range(year + 1), gap[:year + 1], color=INK, lw=1.6, label="corrected groin")
        ax.plot(year, gap[year], "o", color=INK, ms=4)
        ax.set_xlim(0, a.years)
        top = max(max(gap), ref.gap_m.max() if ref is not None else 0) * 1.15
        ax.set_ylim(0, top)
        ax.set_xlabel("Year")
        ax.set_ylabel("Gap across the groin (m)")
        ax.set_title("The gap slows as BRIE's return flow catches up", loc="left", fontsize=7.5)
        ax.legend(loc="lower right", fontsize=6.5, frameon=False)
        ax.grid(axis="y")
        open_frame(ax)
        fig.suptitle(f"How the corrected groin works: one groin, coast at {th:+d}$^\\circ$", fontsize=10)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(OUT), dpi=120):
        for year in range(a.years + 1):
            draw(year)
            writer.grab_frame()
        for _ in range(HOLD_FRAMES):
            writer.grab_frame()
    plt.close(fig)
    print(f"wrote {OUT.relative_to(sc.REPO)}; held back {held[0]:.0f} m in year 1, {held[-1]:.1f} m by year {a.years}; "
          f"carried back {back[1]:.0f} m in year 1, {back[-1]:.0f} m by year {a.years}; gap {gap[-1]:.0f} m")


if __name__ == "__main__":
    main()
