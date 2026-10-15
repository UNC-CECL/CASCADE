#!/usr/bin/env python3
"""
The dipole (version 1) and the corrected groin (version 6) side by side, year by year, on two straight coasts.

    python straight_coast_v1_vs_v6_gif.py   ->  figures/straight_coast_v1_vs_v6.gif

Reruns both modules in straight_coast_test.py's emulator of BRIE's alongshore step (about 20 s), keeping
every year's shoreline, and animates them as plan views: land below the shoreline, ocean above, the
groin in red, and each panel's running total of sand gained or lost. The GIF is git-ignored; rerun to make it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
"""
from __future__ import annotations

import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, figsize  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = Path(__file__).resolve().parent / "figures" / "straight_coast_v1_vs_v6.gif"
PAIR = (("dipole", "Version 1: dipole",
         f"moves a fixed {sc.M:g} m a year north,\nwhatever the waves do"),
        ("corrected", "Version 6: corrected",
         f"holds back {sc.BLOCKING_FRACTION:.0%} of the drift arriving,\nand puts back BRIE's lost or made sand"))
ANGLES = ((-20, "coast turned against the drift (−20°)"), (10, "coast turned with the drift (+10°)"))
YEARS = 50
CELLS = 12                          # cells either side of the groin drawn
HOLD = 20                           # frames the last year is held
# -----------------------------------------------------------------------------


# Every year's shoreline (seaward +, updrift at positive offsets) and the reach's total sand change
def runs():
    sc.PROFILE_YEARS = tuple(range(1, YEARS + 1))            # keep every year, not just the five sampled
    t = sc.brie_tables()
    out = {}
    for name, *_ in PAIR:
        for a, _ in ANGLES:
            rows, prof = sc.run(t, a, name, YEARS)
            out[name, a] = dict(
                prof=np.array([prof[y] for y in range(1, YEARS + 1)]),
                total=rows.reach_change_m.to_numpy())
    return out


def main():
    apply_style()
    data = runs()
    offs = np.arange(sc.NY) - sc.LO                       # as the tables' domain_from_groin; the groin face is at +0.5
    keep = (offs >= -CELLS) & (offs <= CELLS + 1)
    lim = max(np.abs(d["prof"][:, keep]).max() for d in data.values()) * 1.1
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.82, bottom=0.1, hspace=0.32, wspace=0.08)
    stamp = fig.text(0.5, 0.975, "", ha="center", va="top", fontsize=12, fontweight="bold")
    for j, (name, title, rule) in enumerate(PAIR):
        fig.text(0.1 + 0.445 * j + 0.22, 0.925, title, ha="center", fontsize=9.5, fontweight="bold",
                 color=sc.MODULES[name]["col"])
        fig.text(0.1 + 0.445 * j + 0.22, 0.905, rule, ha="center", va="top", fontsize=7, color=INK_MUTED)

    def draw(k):
        yr = min(k, YEARS - 1)
        stamp.set_text(f"Year {yr + 1}")
        for i, (a, label) in enumerate(ANGLES):
            for j, (name, *_) in enumerate(PAIR):
                ax = axes[i, j]
                ax.clear()
                d = data[name, a]
                y = d["prof"][yr][keep]
                x = offs[keep]
                xs = np.repeat(np.r_[x - 0.5, x[-1] + 0.5], 2)[1:-1]
                ys = np.repeat(y, 2)
                ax.fill_between(xs, -lim, ys, color=C["ADDED_FILL"], lw=0)
                ax.fill_between(xs, ys, lim, color=C["WATER"], alpha=0.45, lw=0)
                ax.plot(xs, ys, color=INK, lw=1.3)
                ax.axhline(0, color=INK_MUTED, lw=0.8, ls=(0, (4, 2)))
                ax.plot([0.5, 0.5], [-0.25 * lim, 0.55 * lim], color=C["GROIN"], lw=2.4, solid_capstyle="butt")
                ax.annotate("", xy=(1.5, 0.8 * lim), xytext=(CELLS - 0.5, 0.8 * lim),
                            arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.0))
                ax.text(CELLS / 2 + 0.5, 0.86 * lim, "drift", ha="center", fontsize=6.5, color=INK_MUTED)
                ax.text(-CELLS + 0.3, 0.86 * lim, "south (downdrift)", fontsize=6.5, color=INK_MUTED)
                ax.text(CELLS - 0.3, -0.92 * lim, "north (updrift)", fontsize=6.5, color=INK_MUTED, ha="right")
                tot = d["total"][yr]
                ok = abs(tot) <= 30
                ax.text(-CELLS + 0.3, -0.92 * lim,
                        f"sand made or lost: {tot:+,.0f} m" + (" (balanced)" if ok else ""),
                        fontsize=7.2, color=C["REF"] if ok else C["GROIN"], fontweight="bold")
                ax.set_xlim(-CELLS - 0.5, CELLS + 1.5)
                ax.set_ylim(-lim, lim)
                ax.set_title(label, fontsize=8, loc="left")
                if j == 0:
                    ax.set_ylabel("Shoreline (m)\n+ seaward")
                if i == 1:
                    ax.set_xlabel("Cells from the groin (500 m each)")
        return []

    anim = FuncAnimation(fig, draw, frames=YEARS + HOLD, blit=False)
    anim.save(OUT, writer=PillowWriter(fps=5), dpi=110)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
