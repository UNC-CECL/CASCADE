#!/usr/bin/env python3
"""
How the six groin modules differ: what each one reads, and how much sand it moves.

    python straight_coast_module_rules_figure.py   ->  comparison/figures/straight_coast_module_rules.png

Top row: a plan-view sketch of one module per panel, the two cells either side of the groin and
the rule that moves their shorelines each year. Bottom row: the sand each module moves in its
first year on a straight coast, against the coast's orientation, read from the module tables
that straight_coast_test.py writes (run that first).

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

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "1-straight-coast"))   # the test, kept in place
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
OUT = Path(__file__).resolve().parent / "comparison" / "figures" / "straight_coast_module_rules.png"
SAND, OCEAN = C["ADDED_FILL"], C["WATER"]
PER_ROW = 4                          # modules per row of sketches
MAX_THETA = 30                       # beyond the zero-drift angle the drift reverses and the sides swap
# Per module: what it reads, the updrift/downdrift labels, the sketch's updrift gain and downdrift loss
SKETCH = {
    "dipole": dict(reads="reads nothing", up="+M", down="$-$M", gain=1.0, loss=1.0, drift=False),
    "pinned": dict(reads="reads the step across the groin\n$\\times$ each cell's own diffusivity",
                   up="+f$_{block}$ r$_{up}$ s", down="$-$f$_{block}$ r$_{down}$ s", gain=1.35, loss=0.7, drift=True),
    "conserving": dict(reads="reads the step across the groin\n$\\times$ the groin face's diffusivity",
                       up="+f$_{block}$ r$_{face}$ s", down="$-$f$_{block}$ r$_{face}$ s", gain=1.0, loss=1.0, drift=True),
    "drift": dict(reads="reads the wave climate's\nnet drift at the groin",
                  up="+f$_{block}$ Q$_{net}$", down="$-$f$_{block}$ Q$_{net}$", gain=1.0, loss=1.0, drift=True),
    "hybrid": dict(reads="reads the net drift, and\ncorrects BRIE's r mismatch",
                   up="+f$_{block}$ Q$_{net}$\n+ f$_{block}$ $\\Delta$r s", down="$-$f$_{block}$ Q$_{net}$", gain=1.25, loss=1.0, drift=True),
    "corrected": dict(reads="reads the net drift, and\nreturns BRIE's leak\nnear the groin",
                      up="+f$_{block}$ Q$_{net}$", down="$-$f$_{block}$ Q$_{net}$", gain=1.0, loss=1.0, drift=True),
    "length": dict(reads="as corrected, scaled by\ng = (1$-$a/L)$^{2/3}$: zero once\nthe fillet reaches the tip",
                   up="+f$_{block}$ g Q$_{net}$", down="$-$f$_{block}$ g Q$_{net}$", gain=0.8, loss=0.8, drift=True),
}
# -----------------------------------------------------------------------------


# Plan view: land above, ocean below, the updrift side on the right, the groin between two cells
def sketch(ax, name, spec, i):
    s = SKETCH[name]
    x = np.array([0, 4, 4, 5, 5, 6, 6, 10])
    y = np.array([0, 0, s["loss"], s["loss"], -s["gain"], -s["gain"], 0, 0])
    ax.fill_between(x, y, 3.2, color=SAND, lw=0)
    ax.fill_between(x, -4, y, color=OCEAN, lw=0)
    ax.fill_between([5, 6], [-s["gain"]] * 2, 0, color=spec["col"], alpha=0.45, lw=0)
    ax.fill_between([4, 5], 0, [s["loss"]] * 2, color="white", alpha=0.75, lw=0)
    ax.plot([0, 10], [0, 0], color=INK_MUTED, lw=0.7, ls=(0, (3, 2)))
    ax.plot(x, y, color=INK, lw=1.2)
    for xc in (3, 4, 6, 7):
        ax.plot([xc, xc], [0.25, -0.25], color=INK_MUTED, lw=0.5)
    ax.plot([5, 5], [1.2, -2.2], color=C["GROIN"], lw=2.2, solid_capstyle="butt")
    ax.text(6.15, -s["gain"] / 2, s["up"], ha="left", va="center", fontsize=6.0, color=INK)
    ax.text(3.85, s["loss"] / 2, s["down"], ha="right", va="center", fontsize=6.0, color=INK)
    if s["drift"]:
        ax.annotate("", xy=(5.15, -2.9), xytext=(9.6, -2.9),
                    arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3))
        ax.annotate("", xy=(1.0, -2.9), xytext=(4.85, -2.9),
                    arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.7, ls=(0, (2, 1.5))))
        ax.text(7.4, -3.15, "drift", ha="center", va="top", fontsize=6.5, color=INK)
        ax.text(2.9, -3.15, "1$-$f$_{block}$ passes", ha="center", va="top", fontsize=6.5, color=INK_MUTED)
    else:
        ax.annotate("", xy=(1.0, -2.9), xytext=(9.6, -2.9),
                    arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.7, ls=(0, (2, 1.5))))
        ax.text(7.4, -3.15, "drift (not read)", ha="center", va="top", fontsize=6.5, color=INK_MUTED)
    ax.text(5, 2.25, s["reads"], ha="center", va="top", fontsize=6.6, color=INK)
    ax.text(9.8, -1.8, "updrift", ha="right", va="top", fontsize=6.3, color=INK_MUTED, style="italic")
    ax.text(0.2, -1.8, "downdrift", ha="left", va="top", fontsize=6.3, color=INK_MUTED, style="italic")
    ax.set_xlim(0, 10)
    ax.set_ylim(-4, 3.2)
    ax.set_xticks([])
    ax.set_yticks([])
    for side in ax.spines.values():
        side.set_color(INK_MUTED)
        side.set_linewidth(0.6)
    _letter_inside(ax, i)
    ax.set_title(spec["label"], loc="center", color=spec["col"], fontweight="bold")


def first_year():
    out = {}
    for name, spec in sc.MODULES.items():
        f = sc.module_dir(spec) / "tables" / f"straight_coast_{spec['folder']}_by_year.csv"
        d = pd.read_csv(f)
        out[name] = d[(d.year == 1) & (d.theta0_deg <= MAX_THETA)].sort_values("theta0_deg")
    return out


def main():
    apply_style()
    t = sc.brie_tables()
    th = np.arange(-40, 41)
    q = sc.drift(t, th)
    k = np.flatnonzero(np.diff(np.sign(q)) != 0)
    zero_q = float(np.interp(0, q[k[0]:k[0] + 2][::-1], th[k[0]:k[0] + 2][::-1]))
    yr1 = first_year()
    lo = min(v.applied_up_m.min() for v in yr1.values())
    hi = max(v.applied_up_m.max() for v in yr1.values())
    pad = 0.08 * (hi - lo)

    n = len(sc.MODULES)
    blocks = int(np.ceil(n / PER_ROW))
    fig, axes = plt.subplots(2 * blocks, PER_ROW, figsize=figsize("double", height=4.4 * blocks),
                             constrained_layout=True)
    for j, (name, spec) in enumerate(sc.MODULES.items()):
        row, col = 2 * (j // PER_ROW), j % PER_ROW
        sketch(axes[row, col], name, spec, j)
        ax = axes[row + 1, col]
        ax.axhspan(lo - pad, 0, color=C["BASE_FILL"], alpha=0.6, lw=0)
        for other, ospec in sc.MODULES.items():
            if other != name:
                ax.plot(yr1[other].theta0_deg, yr1[other].applied_up_m, color="0.78", lw=0.9)
        d = yr1[name]
        ax.plot(d.theta0_deg, d.applied_up_m, color=spec["col"], lw=1.8, marker="o", ms=2.5)
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.axvline(0, color=INK_MUTED, lw=0.5, ls=(0, (2, 2)))
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel("Orientation ($^\\circ$)")
        if col == 0:
            ax.set_ylabel("Sand moved to the updrift cell\nin year 1 (m of shoreline)")
            ax.text(-38, lo - pad * 0.4, "trapped on the downdrift side", fontsize=6.3, color=INK_MUTED,
                    va="bottom")
        else:
            ax.tick_params(labelleft=False)
        ax.grid(axis="y")
        open_frame(ax)
        _letter_inside(ax, j + n)
    for j in range(n, blocks * PER_ROW):                     # the empty slots in the last row
        axes[2 * (j // PER_ROW), j % PER_ROW].axis("off")
        axes[2 * (j // PER_ROW) + 1, j % PER_ROW].axis("off")
    fig.suptitle("What does each groin version read, and how much sand does it move?")
    save(fig, OUT, dpi=250, close=True)

    p = {n: yr1[n].set_index("theta0_deg").applied_up_m for n in yr1}
    pin = pd.read_csv(sc.module_dir(sc.MODULES["pinned"]) / "tables"
                      / f"straight_coast_{sc.MODULES['pinned']['folder']}_by_year.csv")
    bent = pin[(pin.theta0_deg == 20) & (pin.year == 10)].iloc[0]
    record_caption(OUT, (
        "**Tests:** nothing new; this figure explains the seven groin modules compared in the straight-coast "
        "test, by what each one reads and how much sand it moves. **How to read:** (a-g) Plan-view sketches, "
        "not to scale: land (sand) above, ocean (blue) below, the updrift side on the right and the net drift "
        "arriving from the right. The red line is the groin, between two 500 m cells; the dashed line is the "
        "straight coast before the module acts. Each year the module moves the updrift cell's shoreline "
        "seaward (shaded) and the downdrift cell's landward (pale) by the amounts marked. (a) The source/sink "
        f"dipole moves a fixed M = {sc.M:g} m per year whatever the coast and waves do. (b) The pinned trapping "
        "version reads BRIE's own alongshore exchange across the groin, the shoreline step s between the two "
        "cells times the diffusivity r, and cancels a fraction f$_{block}$ of it (the blocking fraction); each cell uses its own r, so once the "
        "coast bends at the groin the two sides move by different amounts. (c) The conserving version is the "
        "same with the groin face's r on both sides, so the two amounts always match. (d) The drift version "
        f"cancels a fraction f$_{{block}}$ of the wave climate's net longshore drift Q$_{{net}}$ arriving at the groin, "
        f"converted to shoreline metres over one 500 m cell and the {t['depth']:.1f} m active depth; f$_{{block}}$ = {sc.BLOCKING_FRACTION} "
        "in (b)-(g). (e) The hybrid is the drift version plus the pinned version's row correction: the updrift "
        "cell also gets f$_{block}$ times the difference between the two cells' r, times the step s, which is what the "
        "pinned version moves beyond the conserving one. (f) The corrected version is the drift version plus a leak correction within five cells of the groin: each year it gives back what BRIE's row-scaled solve creates or destroys there, so sand is conserved without changing BRIE. Its first-year transfer is the drift version's. (g) The corrected version with a groin length L (110 m): the same, with f$_{block}$ scaled by (1 $-$ a/L)$^{2/3}$, a the updrift cell's build-out, so it passes all the drift once the fillet reaches the groin's tip; on a straight start a = 0, so its first-year transfer is the corrected version's. (h-n) What each module moves in its first year on a straight coast, against the coast's "
        "orientation (its own colour; the others in grey). Below zero (shaded) the module builds the "
        f"downdrift side instead. Orientations stop at +{MAX_THETA} degrees: the net drift reaches zero at "
        f"{zero_q:+.0f} degrees and reverses beyond it, which swaps the two sides. **Shows:** "
        "the versions differ in what drives them. The dipole moves the same 12 m at every orientation. The pinned "
        "and conserving versions read the step across the groin, which on a straight coast is the coast's own "
        "tilt, so they move nothing at 0 degrees (no step, and BRIE's diffusivity is zero there), build the "
        f"updrift side on coasts turned against the drift ({p['pinned'].loc[-30]:+.0f} m at -30 degrees) and the "
        f"downdrift side on coasts turned with it ({p['pinned'].loc[20]:+.0f} m at +20 degrees). The two are "
        "identical in year 1 because a straight coast has one r everywhere; they part once a fillet bends the "
        f"coast (pinned at +20 degrees, year 10: r$_{{up}}$ {bent.r_lo:.2f} vs r$_{{down}}$ {bent.r_hi:.2f}). "
        f"Only the drift version follows the drift itself: {p['drift'].loc[0]:+.0f} m at 0 degrees, falling to "
        f"zero near {zero_q:+.0f} degrees where the drift stops."))
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
