#!/usr/bin/env python3
"""
Animate the straight-coast test: the four groin modules year by year at six orientations.

    python straight_coast_gif.py [--years 50]   ->  comparison/figures/straight_coast_over_time.gif

Reuses the emulator and modules of straight_coast_test.py, recording every year instead of
five. Each panel keeps one y range for the whole animation, and the text in each panel is the
total sand change across the reach so far, per module (summed shoreline metres, seaward +).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK_MUTED, _title, apply_style, figsize, open_frame)

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = HERE / "comparison" / "figures" / "straight_coast_over_time.gif"
FPS = 5
HOLD_FRAMES = 6                          # repeats of the last frame before the loop restarts
# -----------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--years", type=int, default=50)
    a = ap.parse_args()
    apply_style()
    sc.PROFILE_YEARS = tuple(range(1, a.years + 1))
    t = sc.brie_tables()

    # profiles[(module, theta0)][year] and the reach total by year
    profiles, totals = {}, {}
    for name in sc.MODULES:
        for th in sc.SHOWN:
            frame, prof = sc.run(t, th, name, a.years)
            prof[0] = np.zeros(sc.NY)
            profiles[(name, th)] = prof
            totals[(name, th)] = np.r_[0.0, frame.reach_change_m.to_numpy()]
    lims = {th: max(5.0, np.ceil(max(np.abs(profiles[(n, th)][y]).max() for n in sc.MODULES
                                     for y in range(a.years + 1)) / 10) * 10 * 1.08)
            for th in sc.SHOWN}

    dom = np.arange(sc.NY) - sc.LO - 0.5
    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.4), sharex=True)
    fig.subplots_adjust(left=0.09, right=0.98, top=0.83, bottom=0.1, hspace=0.38, wspace=0.3)
    handles = [Line2D([], [], color=s["col"], lw=1.6, label=s["label"]) for s in sc.MODULES.values()]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.94), ncol=4,
               fontsize=7, frameon=False)

    def draw(year):
        for i, (ax, th) in enumerate(zip(axes.flat, sc.SHOWN)):
            ax.clear()
            ax.axvline(0, color=C["GROIN"], lw=0.8)
            ax.axhline(0, color=INK_MUTED, lw=0.5)
            lines = []
            for name, spec in sc.MODULES.items():
                ax.plot(dom, profiles[(name, th)][year], color=spec["col"], lw=1.5)
                lines.append((spec["col"], totals[(name, th)][year]))
            ax.text(0.03, 0.97, "total sand (m)", transform=ax.transAxes, fontsize=6.0,
                    color=INK_MUTED, va="top")
            for k, (col, tot) in enumerate(lines):
                ax.text(0.03, 0.97 - 0.075 * (k + 1), f"{tot:+.0f}", transform=ax.transAxes,
                        fontsize=6.3, color=col, va="top")
            ax.set_xlim(-12, 12)
            ax.set_ylim(-lims[th], lims[th])
            ax.tick_params(labelsize=6.5)
            ax.grid(axis="y")
            open_frame(ax)
            q = float(sc.drift(t, th)[0])
            _title(ax, i, f"coast at {th:+d}$^\\circ$")
            ax.text(0.97, 0.04, f"drift {abs(q) / 1e3:.0f}k m$^3$/yr", transform=ax.transAxes,
                    fontsize=6.3, color=INK_MUTED, ha="right")
        for ax in axes[:, 0]:
            ax.set_ylabel("Shoreline change\n(m, seaward +)")
        for ax in axes[1]:
            ax.set_xlabel("Domains from the groin (updrift right)", fontsize=7)
        fig.suptitle(f"Year {year}", fontsize=10, y=0.99)

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
