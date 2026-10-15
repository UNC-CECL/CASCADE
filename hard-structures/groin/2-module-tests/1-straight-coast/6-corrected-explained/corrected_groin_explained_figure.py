#!/usr/bin/env python3
"""
How the corrected groin works, in one figure: a sketch of the three flows and the balance they reach.

    python corrected_groin_explained_figure.py [--years 60] [--theta -20]
        ->  figures/corrected_groin_explained.png

(a) a plan-view sketch, not to scale; (b) the sand the groin holds back and the sand BRIE carries
back across the groin each year; (c) the gap across the groin. (b) and (c) come from module 6 in
straight_coast_test.py, computed as corrected_groin_gif.py computes them; both are imported unchanged.

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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
import corrected_groin_gif as cg  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = HERE / "figures" / "corrected_groin_explained.png"
SAND, OCEAN = C["ADDED_FILL"], C["WATER"]
HELD_C, BACK_C = C["GROIN"], C["LATE"]
# -----------------------------------------------------------------------------


# (a) Plan view: land above, ocean below, updrift on the right; the groin between two cells
def sketch(ax):
    x = np.array([0, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 10])
    y = np.array([0, 0, 0.35, 0.35, 0.8, 0.8, -1.1, -1.1, -0.55, -0.55, 0, 0])   # seaward is down
    ax.fill_between(x, y, 3.3, color=SAND, lw=0)
    ax.fill_between(x, -4.3, y, color=OCEAN, lw=0)
    ax.plot([0, 10], [0, 0], color=INK_MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.plot(x, y, color=INK, lw=1.3)
    ax.plot([5, 5], [1.8, -2.6], color=C["GROIN"], lw=3, solid_capstyle="butt")
    ax.text(5, 2.0, "groin", ha="center", va="bottom", fontsize=7.5, color=C["GROIN"])
    ax.text(7.6, 1.6, "land", fontsize=7, color=INK_MUTED)
    ax.text(0.3, -3.9, "ocean", fontsize=7, color=INK_MUTED)
    ax.text(1.4, 0.95, "downdrift:\ncut back", ha="center", fontsize=7, color=INK_MUTED)
    ax.text(8.6, -1.05, "updrift:\nbuilt out", ha="center", fontsize=7, color=INK_MUTED)
    ax.plot([4.0, 5.75], [0.8, 0.8], color=INK, lw=0.6, ls=(0, (1, 1.5)))           # the downdrift shoreline, carried across
    ax.annotate("", xy=(5.6, 0.8), xytext=(5.6, -1.1), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
    ax.text(5.75, -0.15, "gap", fontsize=7.5, color=INK, va="center")
    # 1. the drift arriving
    ax.annotate("", xy=(5.3, -3.4), xytext=(9.7, -3.4), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
    ax.text(7.5, -3.2, "1. waves bring sand\ntoward the groin", ha="center", va="bottom", fontsize=7.5, color=INK)
    # 2. held back: downdrift cell to updrift cell
    ax.annotate("", xy=(5.6, -1.6), xytext=(4.4, -1.6),
                arrowprops=dict(arrowstyle="-|>", color=HELD_C, lw=1.6, connectionstyle="arc3,rad=0.6"))
    ax.text(3.2, -2.55, "2. the groin holds back\n60% of it each year", ha="center", va="center", fontsize=7.5,
            color=HELD_C)
    # 3. carried back by BRIE
    ax.annotate("", xy=(4.4, 1.25), xytext=(5.6, 1.25),
                arrowprops=dict(arrowstyle="-|>", color=BACK_C, lw=1.6, connectionstyle="arc3,rad=0.6"))
    ax.text(7.3, 2.45, "3. BRIE treats the groin as open\nand carries sand back across it", ha="center",
            va="center", fontsize=7.5, color=BACK_C)
    ax.text(0.9, 3.05, "4. a correction puts back BRIE's leak near the groin, so no sand is made or lost",
            fontsize=6.8, color=INK_MUTED, va="top")
    ax.set_xlim(0, 10)
    ax.set_ylim(-4.3, 3.3)
    ax.set_xticks([])
    ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_color(INK_MUTED)
    _letter_inside(ax, 0)
    ax.set_title("What happens each year (sketch, not to scale)", loc="center", fontsize=8.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=60)
    ap.add_argument("--theta", type=int, default=-20)
    a = ap.parse_args()
    apply_style()
    t = sc.brie_tables()
    th = a.theta
    sc.PROFILE_YEARS = tuple(range(1, a.years + 1))
    _, prof = sc.run(t, th, "corrected", a.years)
    shore = {0: np.zeros(sc.NY), **prof}
    yrs = np.arange(a.years + 1)
    held = np.array([cg.held_back(t, th, shore[y])[0] for y in yrs])
    back = np.array([cg.carried_back(t, th, shore[y]) for y in yrs])
    gap = np.array([shore[y][sc.HI] - shore[y][sc.LO] for y in yrs])

    fig = plt.figure(figsize=figsize("double", height=5.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.45, 1))
    sketch(fig.add_subplot(gs[0, :]))
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(yrs[1:], held[1:], color=HELD_C, lw=1.8)
    ax.plot(yrs[1:], back[1:], color=BACK_C, lw=1.8)
    ax.text(a.years, held[-1] + 3, "held back by the groin", color=HELD_C, fontsize=7, ha="right", va="bottom")
    ax.text(a.years, back[-1] - 3, "carried back by BRIE", color=BACK_C, fontsize=7, ha="right", va="top")
    ax.set_xlim(0, a.years)
    ax.set_ylim(0, held.max() * 1.25)
    ax.set_xlabel("Year")
    ax.set_ylabel("Sand per year\n(m of shoreline)")
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 1)
    ax.set_title("The return flow catches up", loc="center", fontsize=8.5)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(yrs, gap, color=INK, lw=1.8)
    ax.set_xlim(0, a.years)
    ax.set_ylim(0, gap.max() * 1.2)
    ax.set_xlabel("Year")
    ax.set_ylabel("Gap across the groin (m)")
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 2)
    ax.set_title("So the gap levels off", loc="center", fontsize=8.5)
    fig.suptitle("How the corrected groin works", fontsize=10)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    save(fig, OUT, dpi=250, close=True)
    record_caption(OUT, (
        "**Tests:** nothing new; this figure explains how the corrected trapping groin (straight-coast module 6) "
        "works. **How to read:** (a) a plan-view sketch, not to scale: land above, ocean below, the updrift side "
        "on the right, the dashed line the coast before the groin. Each year: (1) the waves bring the net longshore "
        f"drift toward the groin; (2) the groin holds back {sc.BLOCKING_FRACTION:.0%} of it, moving that sand from "
        "the cell just downdrift to the cell just updrift; (3) BRIE's alongshore diffusion does not know the groin "
        "is there and carries sand back across the groin face, more as the step across it grows; (4) a correction "
        "returns what BRIE's solve creates or destroys near the groin, so total sand is conserved. (b) and (c): "
        f"the model, one groin on a straight coast at {th:+d} degrees, option A waves. (b) The sand held back and "
        "the sand carried back each year, as metres of shoreline in one 500 m cell (the carried-back flow is "
        "approximate: the groin face's diffusion number times the step). (c) The gap between the updrift and "
        f"downdrift cells. **Shows:** the groin holds back nearly the same every year ({held[1]:.0f} m in year 1, "
        f"{held[-1]:.1f} m in year {a.years}), while what BRIE carries back grows with the gap ({back[1]:.0f} m in "
        f"year 1, {back[-1]:.0f} m in year {a.years}); the rest spreads along the coast. The gap therefore grows "
        f"fast at first ({gap[1]:.0f} m after year 1) and slowly later ({gap[a.years // 2]:.0f} m at year "
        f"{a.years // 2}, {gap[-1]:.0f} m at year {a.years})."))
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
