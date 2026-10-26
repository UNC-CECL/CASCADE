#!/usr/bin/env python3
"""
How the fixed-shift groin worked, in two figures.

    python fixed_shift_explained_figure.py   ->  figures/fixed_shift_1_in_the_model_year.png
                                            figures/fixed_shift_2_explained.png

(1) where the fixed-shift groin sat in a model year: Barrier3D (with its background erosion) gives each domain its
shoreline change, the fixed-shift groin edits two domains of it, and BRIE moves the result alongshore; with an
illustrative example of the edit. (2) a plan-view sketch of what it did each year, and, from the
straight-coast runs, the fixed amount it moved against the sand BRIE's solve made or lost at the bend it
created. The companion of wave_driven_explained_figure.py. Reads the fixed-shift groin's by-year table from
straight_coast_test.py; nothing is rerun.

The fixed-shift groin (cascade/groin.py GroinCallback, removed from the runner on 2026-10-13 in favour of the wave-driven
groin; in git before 563a542e) was a callback on the sandbox Cascade (cascade/cascade_groin.py,
`_groin_callback`). Each model year, after every Barrier3D domain had run and before BRIE's alongshore step,
it added M of seaward change to the updrift domain (GIS 6) and M of landward change to the downdrift one
(GIS 5) in x_s_dt, the shoreline change each domain hands to BRIE. It never read the waves.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
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
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

RUN = Path(__file__).resolve().parents[2] / "3-module-tests" / "2-straight-coast" / "1-run"
sys.path.insert(0, str(RUN))   # the straight-coast test
import straight_coast_test as sc  # noqa: E402  (puts scripts/ on the path)
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
FIGS = workspace_twin(HERE) / "figures"
BY_YEAR = workspace_twin(RUN) / "tables" / "straight_coast_fixed_shift_by_year.csv"
SAND, OCEAN = C["ADDED_FILL"], C["WATER"]
MOVED_C, LEAK_C = C["ADDED"], C["GROIN"]
ANGLES = ((-20, "coast turned against the drift (−20°)", C["LATE"]),
          (10, "coast turned with the drift (+10°)", C["EARLY"]))
EXAMPLE_BE = 1.5                     # m/yr of retreat in every domain, the illustrative Barrier3D change in (1b)
# -----------------------------------------------------------------------------


def box(ax, x, y, w, h, text, fc, ec=INK, fs=7.5, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK,
            fontweight="bold" if bold else "normal", linespacing=1.3)


def arrow(ax, x0, y0, x1, y1, col=INK, lw=1.4):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=col, lw=lw))


# 1. Where the fixed-shift groin sat in a model year
def fig_model_year():
    fig = plt.figure(figsize=figsize("double", height=4.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=(1.25, 1))
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    box(ax, 0.2, 7.3, 8.6, 2.2, "1. Barrier3D runs in every domain\nstorms, overwash, dunes and its "
        "background erosion (BE)\n→ each domain's shoreline change this year", "#eef3f8")
    box(ax, 0.2, 4.0, 8.6, 2.5, "2. The fixed-shift groin edits two domains of that change\nGIS 6 (north, updrift): "
        f"{sc.M:g} m seaward\nGIS 5 (south, downdrift): {sc.M:g} m landward\nthe same {sc.M:g} m every year; "
        "it never reads the waves", "#fbeed3", ec=MOVED_C, bold=False)
    box(ax, 0.2, 0.5, 8.6, 2.7, "3. BRIE moves the shoreline alongshore\nit smooths the edited shoreline "
        "between neighbouring domains;\nit does not know a groin is there, and at the bend the\nfixed shift makes "
        "its step makes or loses sand", "#eef3f8")
    arrow(ax, 4.5, 7.3, 4.5, 6.55)
    arrow(ax, 4.5, 4.0, 4.5, 3.25)
    ax.annotate("", xy=(8.9, 8.4), xytext=(8.9, 1.8),
                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.0, connectionstyle="arc3,rad=0.45"))
    ax.text(9.95, 5.1, "next\nyear", fontsize=6.8, color=INK_MUTED, va="center", ha="right")
    _letter_inside(ax, 0)
    ax.set_title("One model year, in order", loc="center", fontsize=8.5)

    ax = fig.add_subplot(gs[0, 1])
    gis = np.arange(2, 10)
    before = -np.full(gis.size, EXAMPLE_BE)                  # seaward +
    after = before.copy()
    after[gis == 6] += sc.M
    after[gis == 5] -= sc.M
    w = 0.38
    ax.bar(gis - w / 2, before, width=w, color="0.75", label="from Barrier3D (example)")
    ax.bar(gis + w / 2, after, width=w, color=[MOVED_C if g in (5, 6) else "0.45" for g in gis],
           label="to BRIE, after the fixed-shift groin")
    ax.axvline(5.5, color=C["GROIN"], lw=2.2)
    ax.text(5.6, -sc.M * 1.3, "groin (GIS 5/6)", color=C["GROIN"], fontsize=6.8, va="bottom")
    ax.annotate(f"+{sc.M:g} m", xy=(6 + w / 2, after[gis == 6][0]), xytext=(7.0, sc.M * 0.85), fontsize=7,
                color=MOVED_C, arrowprops=dict(arrowstyle="-", color=MOVED_C, lw=0.6))
    ax.annotate(f"−{sc.M:g} m", xy=(5 + w / 2, after[gis == 5][0]), xytext=(3.0, -sc.M * 1.05), fontsize=7,
                color=MOVED_C, arrowprops=dict(arrowstyle="-", color=MOVED_C, lw=0.6))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(gis, [f"GIS {g}" for g in gis], fontsize=6.5)
    ax.set_xlabel("south ← domain → north")
    ax.set_ylabel("Shoreline change this year (m)\n+ seaward")
    ax.set_ylim(-sc.M * 1.35, sc.M * 1.35)
    ax.legend(frameon=False, fontsize=6.6, loc="upper left", bbox_to_anchor=(0.54, 1.0))
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 1)
    ax.set_title("What step 2 does to the numbers", loc="center", fontsize=8.5)
    fig.suptitle("Where the fixed-shift groin sat in the model", fontsize=10)
    png = workspace_twin(FIGS) / "fixed_shift_1_in_the_model_year.png"
    save(fig, png, dpi=250, close=True)
    record_caption(png, (
        "**Tests:** nothing; this figure explains where the fixed-shift groin (the source/sink dipole, "
        "the runner's groin until 2026-10-13) acted in a CASCADE model year. **How to read:** (a) the order "
        "within one year: every Barrier3D domain runs first, including its background erosion (BE), and hands "
        "BRIE its shoreline change for the year; the fixed-shift groin then edits that change in the two domains either side "
        "of the groin, adding M seaward to GIS 6 (north, updrift) and M landward to GIS 5 (south, downdrift); BRIE "
        "then moves the shoreline alongshore. The fixed-shift groin was a callback on the sandbox Cascade "
        "(`cascade/cascade_groin.py`, `_groin_callback`), not part of the BE rates, though it entered the same "
        "number BE ends up in. (b) an illustrative example of the edit: grey, a uniform "
        f"{EXAMPLE_BE:g} m of retreat from Barrier3D in every domain (made up for the example); dark, what BRIE "
        f"receives, the same except GIS 5 and GIS 6, moved {sc.M:g} m apart (gold). **Shows:** the fixed-shift groin is a "
        "fixed transfer between two neighbouring domains each year, set by one number (M) and blind to the waves; "
        "everything else that happens to the groin's shoreline is BRIE's doing (figure fixed_shift_2)."))


# (a) of figure 2: plan view, land above, ocean below, updrift on the right
def sketch(ax):
    x = np.array([0, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 10])
    y = np.array([0, 0, 0.3, 0.3, 0.7, 0.7, -1.0, -1.0, -0.45, -0.45, 0, 0])
    ax.fill_between(x, y, 3.3, color=SAND, lw=0)
    ax.fill_between(x, -4.3, y, color=OCEAN, lw=0)
    ax.plot([0, 10], [0, 0], color=INK_MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.plot(x, y, color=INK, lw=1.3)
    ax.plot([5, 5], [1.8, -2.6], color=C["GROIN"], lw=3, solid_capstyle="butt")
    ax.text(5, 2.0, "groin", ha="center", va="bottom", fontsize=7.5, color=C["GROIN"])
    ax.text(7.8, 1.6, "land", fontsize=7, color=INK_MUTED)
    ax.text(0.3, -3.9, "ocean", fontsize=7, color=INK_MUTED)
    ax.text(1.3, 0.9, "GIS 5 (south):\n−M each year", ha="center", fontsize=7, color=INK_MUTED)
    ax.text(8.6, -1.05, "GIS 6 (north):\n+M each year", ha="center", fontsize=7, color=INK_MUTED)
    # the drift, drawn but not read
    ax.annotate("", xy=(5.3, -3.4), xytext=(9.7, -3.4),
                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.4, ls="--"))
    ax.text(7.5, -3.2, "the waves' drift: not read", ha="center", va="bottom", fontsize=7.5, color=INK_MUTED)
    # 1. the fixed transfer
    ax.annotate("", xy=(5.6, -1.5), xytext=(4.4, -1.5),
                arrowprops=dict(arrowstyle="-|>", color=MOVED_C, lw=1.8, connectionstyle="arc3,rad=0.6"))
    ax.text(3.0, -2.45, f"1. moves a fixed M = {sc.M:g} m\nfrom GIS 5 to GIS 6", ha="center", va="center",
            fontsize=7.5, color=MOVED_C)
    # 2. BRIE spreads it
    ax.annotate("", xy=(4.4, 1.2), xytext=(5.6, 1.2),
                arrowprops=dict(arrowstyle="-|>", color=C["LATE"], lw=1.6, connectionstyle="arc3,rad=0.6"))
    ax.text(7.2, 2.4, "2. BRIE smooths it along the coast,\nunaware of the groin", ha="center", va="center",
            fontsize=7.5, color=C["LATE"])
    # 3. the leak at the bend
    ax.text(5.32, 0.15, "*", ha="center", va="center", fontsize=16, fontweight="bold", color=LEAK_C)
    ax.text(0.9, 3.05, "3. at the bend, BRIE's step makes or loses sand (*), and nothing puts it back",
            fontsize=6.9, color=LEAK_C, va="top")
    ax.set_xlim(0, 10)
    ax.set_ylim(-4.3, 3.3)
    ax.set_xticks([])
    ax.set_yticks([])
    for s_ in ax.spines.values():
        s_.set_color(INK_MUTED)
    _letter_inside(ax, 0)
    ax.set_title("What happens each year (sketch, not to scale)", loc="center", fontsize=8.5)


# 2. How it worked, and why the sand did not balance
def fig_explained():
    b = pd.read_csv(BY_YEAR)
    fig = plt.figure(figsize=figsize("double", height=5.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.35, 1))
    sketch(fig.add_subplot(gs[0, :]))
    ax = fig.add_subplot(gs[1, 0])
    yrs = np.arange(1, b["year"].max() + 1)
    ax.plot(yrs, np.full(yrs.size, sc.M), color=MOVED_C, lw=2.0)
    ax.text(yrs[-1], sc.M - 1.5, f"moved by the fixed-shift groin ({sc.M:g} m)", color=MOVED_C, fontsize=7, ha="right",
            va="top")
    vals = {}
    for th, lab, col in ANGLES:
        s = b[b["theta0_deg"] == th].set_index("year")["solve_error_m"]
        per = s.diff().fillna(s)
        vals[th] = (per, s)
        ax.plot(per.index, per.values, color=col, lw=1.5, label=f"made (+) or lost (−) by BRIE, {lab}")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(0, yrs[-1])
    ax.set_xlabel("Year")
    ax.set_ylabel("Sand per year\n(m of shoreline)")
    ax.legend(frameon=False, fontsize=6.2, loc="center left", bbox_to_anchor=(0.15, 0.55))
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 1)
    ax.set_title("BRIE's leak outgrows the fixed-shift groin", loc="center", fontsize=8.5)
    ax = fig.add_subplot(gs[1, 1])
    for th, lab, col in ANGLES:
        s = vals[th][1]
        ax.plot(s.index, s.values, color=col, lw=1.8)
        ax.text(s.index[-1], s.values[-1], f" {s.values[-1]:+,.0f} m", color=col, fontsize=7, va="center")
    ax.axhspan(-30, 30, color=C["BASE_FILL"], alpha=0.6, lw=0)
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.text(24, -60, "balanced (±30 m)", fontsize=6.5, color=INK_MUTED, va="top")
    ax.set_xlim(0, yrs[-1] + 9)
    ax.set_xlabel("Year")
    ax.set_ylabel("Sand made or lost so far\n(m of summed shoreline)")
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 2)
    ax.set_title("So the reach gains or loses sand", loc="center", fontsize=8.5)
    fig.suptitle("How the fixed-shift groin worked", fontsize=10)
    png = workspace_twin(FIGS) / "fixed_shift_2_explained.png"
    save(fig, png, dpi=250, close=True)
    p20, s20 = vals[-20]
    p10, s10 = vals[10]
    f20 = b[(b["theta0_deg"] == -20)].set_index("year")["face_angle_deg"]
    f10 = b[(b["theta0_deg"] == 10)].set_index("year")["face_angle_deg"]
    record_caption(png, (
        "**Tests:** nothing new; this figure explains how the fixed-shift groin (the source/sink dipole) "
        "worked and why its sand did not balance. The companion of `wave_driven_2_explained.png` "
        "(`../2-wave-driven/figures/`). **How to read:** (a) a plan-view sketch, not to scale: land above, "
        "ocean below, the updrift side (GIS 6, north) on the right, the dashed line the coast before the groin. "
        f"Each year: (1) the fixed-shift groin moves a fixed M = {sc.M:g} m of shoreline from the domain south of the groin "
        "to the one north of it, whatever the waves are doing; (2) BRIE's alongshore step smooths the result "
        "along the coast, not knowing a groin is there; (3) at the bend this makes, BRIE's step makes or loses "
        "sand, because each domain takes its alongshore rate from its own boundary (`../../3-module-tests/"
        "2-straight-coast/3-sand-balance-at-a-step/figures/brie_sand_balance_at_a_step.png`), and nothing puts "
        "it back. (b) and (c), from "
        "the straight-coast runs (BRIE's alongshore step alone, option A waves, 50 years): (b) the sand the fixed-shift groin "
        "moves each year (it balances: GIS 6 gains what GIS 5 loses) against the sand BRIE's step makes or loses "
        "each year, on a coast turned against the drift and one turned with it; (c) the running total over the "
        "61-domain reach; grey, the ±30 m band the wave-driven groin stays inside. **Shows:** the fixed-shift groin's own "
        "transfer is balanced, but the bend it builds keeps steepening (the face angle goes from −20 to "
        f"{f20.iloc[-1]:.0f} degrees on the first coast, from +10 to {f10.iloc[-1]:.0f} on the second), and BRIE's "
        f"leak grows with it: by year 50 it makes {p20.iloc[-1]:+.0f} m a year on the first and "
        f"{p10.iloc[-1]:+.0f} m a year on the second, more than the fixed-shift groin's {sc.M:g} m, and the reach has gained "
        f"{s20.iloc[-1]:+,.0f} m and lost {abs(s10.iloc[-1]):,.0f} m. What the fixed-shift groin built at Buxton therefore "
        "depended mostly on the coast's angle at the groin, not on the drift it is meant to trap; that is why it "
        "was replaced by the wave-driven groin."))


def main():
    apply_style()
    fig_model_year()
    fig_explained()
    print("wrote", FIGS.relative_to(sc.REPO))


if __name__ == "__main__":
    main()
