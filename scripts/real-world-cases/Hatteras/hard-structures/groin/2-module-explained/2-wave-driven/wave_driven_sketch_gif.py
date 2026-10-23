#!/usr/bin/env python3
"""
How the wave-driven groin works, animated in the style of wave_driven_explained_figure.py.

    python wave_driven_sketch_gif.py [--years 60] [--theta -20]
        ->  figures/wave_driven_2_explained.gif

Top: the plan-view sketch, land above and ocean below, with the shoreline each year taken from
module 6 in straight_coast_test.py (stretched across the shore so the change is visible; not to
scale). The red arrow (sand the groin holds back) and the blue arrow (sand BRIE carries back across
the groin) are drawn with widths in proportion to the sand they move. Bottom: the two flows and the
gap, drawn in year by year. The flows are computed as wave_driven_60_years_gif.py computes them; both
scripts are imported unchanged.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import sys
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402

RUN = Path(__file__).resolve().parents[2] / "3-module-tests" / "2-straight-coast" / "1-run"
sys.path.insert(0, str(RUN))   # the straight-coast test
import straight_coast_test as sc  # noqa: E402
import wave_driven_60_years_gif as cg  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, figsize, open_frame  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = workspace_twin(HERE) / "figures" / "wave_driven_2_explained.gif"
SAND, OCEAN = C["ADDED_FILL"], C["WATER"]
HELD_C, BACK_C = C["GROIN"], C["LATE"]
WINDOW = 8                                # domains either side of the groin drawn
FPS = 6
HOLD_FRAMES = 12
ARROW_LW_PER_M = 0.09                     # arrow width per metre of shoreline moved per year
# -----------------------------------------------------------------------------


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
    shore = {0: np.zeros(sc.NY), **prof}                       # seaward +, updrift at higher index
    yrs = np.arange(a.years + 1)
    held = np.array([cg.held_back(t, th, shore[y])[0] for y in yrs])
    back = np.array([cg.carried_back(t, th, shore[y]) for y in yrs])
    gap = np.array([shore[y][sc.HI] - shore[y][sc.LO] for y in yrs])

    cells = np.arange(sc.LO - WINDOW, sc.HI + WINDOW + 1)
    edges = np.r_[cells - sc.LO - 1, cells[-1] - sc.LO]       # cell edges, the groin face at 0
    lim = max(np.abs(shore[a.years][cells]).max(), 10) * 2.2   # metres of shoreline change shown either side
    fig = plt.figure(figsize=figsize("double", height=5.4), constrained_layout=True)

    def draw(year):
        fig.clear()
        gs = fig.add_gridspec(2, 2, height_ratios=(1.45, 1))
        ax = fig.add_subplot(gs[0, :])
        y_shore = -shore[year][cells]                          # seaward is down
        xs = np.repeat(edges, 2)[1:-1]
        ys = np.repeat(y_shore, 2)
        ax.fill_between(xs, ys, lim, color=SAND, lw=0)
        ax.fill_between(xs, -lim, ys, color=OCEAN, lw=0)
        ax.plot([edges[0], edges[-1]], [0, 0], color=INK_MUTED, lw=0.7, ls=(0, (3, 2)))
        ax.plot(xs, ys, color=INK, lw=1.3)
        ax.plot([0, 0], [0.6 * lim, -0.7 * lim], color=C["GROIN"], lw=3, solid_capstyle="butt")
        ax.text(0, 0.63 * lim, "groin", ha="center", va="bottom", fontsize=7.5, color=C["GROIN"])
        # the gap: downdrift shoreline carried across, and the arrow to the updrift shoreline
        y_down, y_up = -shore[year][sc.LO], -shore[year][sc.HI]
        if gap[year] > 1:
            ax.plot([-0.6, 0.75], [y_down, y_down], color=INK, lw=0.6, ls=(0, (1, 1.5)))
            ax.annotate("", xy=(0.6, y_down), xytext=(0.6, y_up), arrowprops=dict(arrowstyle="<->", color=INK, lw=0.9))
            ax.text(0.75, (y_down + y_up) / 2, f"gap {gap[year]:.0f} m", fontsize=7.5, color=INK, va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=1.5))
        # 1. the drift arriving
        ax.annotate("", xy=(1.0, -0.86 * lim), xytext=(edges[-1] - 0.3, -0.86 * lim),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6))
        ax.text((1.0 + edges[-1]) / 2, -0.82 * lim, "1. waves bring sand toward the groin", ha="center",
                va="bottom", fontsize=7.5, color=INK)
        # 2. held back, and 3. carried back, widths in proportion to the sand moved
        ax.annotate("", xy=(0.8, -0.45 * lim), xytext=(-0.8, -0.45 * lim),
                    arrowprops=dict(arrowstyle="-|>", color=HELD_C, lw=max(0.8, held[year] * ARROW_LW_PER_M),
                                    mutation_scale=14, connectionstyle="arc3,rad=0.6"))
        ax.text(-1.4, -0.62 * lim, f"2. the groin holds back {held[year]:.0f} m\nof shoreline this year",
                ha="right", va="center", fontsize=7.5, color=HELD_C)
        if back[year] > 0.5:
            ax.annotate("", xy=(-0.8, 0.42 * lim), xytext=(0.8, 0.42 * lim),
                        arrowprops=dict(arrowstyle="-|>", color=BACK_C, lw=max(0.8, back[year] * ARROW_LW_PER_M),
                                        mutation_scale=14, connectionstyle="arc3,rad=0.6"))
        ax.text(1.4, 0.62 * lim, f"3. BRIE carries {back[year]:.0f} m back across it\n(it treats the groin as open)",
                ha="left", va="center", fontsize=7.5, color=BACK_C)
        ax.text(edges[0] + 0.2, 0.93 * lim, "4. a correction puts back BRIE's leak near the groin, "
                "so no sand is made or lost", fontsize=6.8, color=INK_MUTED, va="top")
        ax.text(edges[0] + 0.2, -0.95 * lim, "ocean", fontsize=7, color=INK_MUTED, va="bottom")
        ax.text(edges[-1] - 0.2, 0.93 * lim, "land", fontsize=7, color=INK_MUTED, va="top", ha="right")
        ax.set_xlim(edges[0], edges[-1])
        ax.set_ylim(-lim, lim)
        ax.set_xticks([])
        ax.set_yticks([])
        for s_ in ax.spines.values():
            s_.set_color(INK_MUTED)
        ax.set_title(f"Year {year}: downdrift (left) cut back, updrift (right) built out  "
                     "(sketch, not to scale)", loc="left", fontsize=8.5)
        # Bottom left: the two flows
        ax = fig.add_subplot(gs[1, 0])
        ax.plot(yrs[1:year + 1], held[1:year + 1], color=HELD_C, lw=1.8)
        ax.plot(yrs[1:year + 1], back[1:year + 1], color=BACK_C, lw=1.8)
        if year >= 1:
            ax.plot([year, year], [held[year], back[year]], "o", color=INK, ms=3)
        ax.text(a.years, held.max() * 1.07, "held back by the groin", color=HELD_C, fontsize=7, ha="right", va="bottom")
        ax.text(a.years, back.max() * 0.9, "carried back by BRIE", color=BACK_C, fontsize=7, ha="right", va="top")
        ax.set_xlim(0, a.years)
        ax.set_ylim(0, held.max() * 1.3)
        ax.set_xlabel("Year")
        ax.set_ylabel("Sand per year\n(m of shoreline)")
        ax.set_title("The return flow catches up", loc="center", fontsize=8.5)
        ax.grid(axis="y")
        open_frame(ax)
        # Bottom right: the gap
        ax = fig.add_subplot(gs[1, 1])
        ax.plot(yrs[:year + 1], gap[:year + 1], color=INK, lw=1.8)
        ax.plot(year, gap[year], "o", color=INK, ms=3)
        ax.set_xlim(0, a.years)
        ax.set_ylim(0, gap.max() * 1.2)
        ax.set_xlabel("Year")
        ax.set_ylabel("Gap across the groin (m)")
        ax.set_title("So the gap levels off", loc="center", fontsize=8.5)
        ax.grid(axis="y")
        open_frame(ax)
        fig.suptitle(f"How the wave-driven groin works (coast at {th:+d}$^\\circ$, option A waves)", fontsize=10)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, str(OUT), dpi=120):
        for year in range(a.years + 1):
            draw(year)
            writer.grab_frame()
        for _ in range(HOLD_FRAMES):
            writer.grab_frame()
    plt.close(fig)
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
