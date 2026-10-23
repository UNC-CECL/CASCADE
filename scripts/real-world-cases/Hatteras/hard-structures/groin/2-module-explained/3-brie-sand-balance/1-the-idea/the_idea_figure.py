#!/usr/bin/env python3
"""
The bookkeeping problem behind BRIE's sand error, with no model in it.

    python the_idea_figure.py   ->  figures/brie_sand_the_idea.png

Two sections and the boundary between them. Every handover of sand across the boundary is written down twice,
once by each section; BRIE lets each work out the amount with its own rate, so where the two rates differ the
books disagree and the difference is sand made or lost. Illustrative numbers.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
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

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from common import *  # noqa: E402,F403  (BRIE runs, drawing helpers, constants; puts scripts/ on the path)
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
FIGS = workspace_twin(HERE) / "figures"

# -----------------------------------------------------------------------------


# 1. The idea, with no model in it
def fig_idea():
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
    titles = ("1  What should happen: sand that leaves A\narrives in B, and both write down the same amount",
              "2  BRIE: each side works out the amount\nitself: its OWN rate × the slope between them",
              "3  If the two sides' rates differ, they write\ndown different amounts for the same handover")
    for k, ax in enumerate(axes):
        blank(ax)
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.set_title(titles[k], fontsize=12.5, pad=10)
        for j, (x, name) in enumerate(((0.4, "section A"), (5.6, "section B"))):
            ax.add_patch(FancyBboxPatch((x, 2.0), 4.0, 4.2, boxstyle="round,pad=0.05,rounding_size=0.3", fc=SAND,
                                        ec=SAND_EDGE, lw=1.5))
            ax.text(x + 2.0, 1.3, name, ha="center", fontsize=11, color=INK_MUTED)
        ax.plot([5.0, 5.0], [1.6, 6.6], color=INK, lw=2, ls=(0, (4, 3)))
        ax.text(5.0, 1.0, "boundary", ha="center", fontsize=9.5, color=INK_MUTED)
        arrow(ax, 3.4, 3.4, 6.6, 3.4, col=INK, lw=3, ms=24)
        ax.text(5.0, 7.6, "sand moves A → B", ha="center", fontsize=10, color=INK_MUTED)
    a, b_, c = axes
    bubble(a, 2.4, 8.4, "A writes:\n10 left me", fs=11)
    bubble(a, 7.6, 8.4, "B writes:\n10 arrived", fs=11)
    a.text(5, 0.2, "−10 + 10 = 0:  nothing made or lost", ha="center", fontsize=12, color=GOOD, fontweight="bold")
    bubble(b_, 2.4, 8.4, "A: rate 2 × slope 5\n= 10 left me", fs=10.5)
    bubble(b_, 7.6, 8.4, "B: rate 2 × slope 5\n= 10 arrived", fs=10.5)
    b_.text(5, 0.2, "same rate on both sides: the books agree", ha="center", fontsize=12, color=GOOD,
            fontweight="bold")
    bubble(c, 2.4, 8.4, "A: rate 2 × slope 5\n= 10 left me", fs=10.5, ec=BAD)
    bubble(c, 7.6, 8.4, "B: rate 0.8 × slope 5\n= only 4 arrived", fs=10.5, ec=BAD)
    c.text(5, 0.2, "−10 + 4 = −6:  6 vanished on the way", ha="center", fontsize=12, color=BAD, fontweight="bold")
    fig.suptitle("The bookkeeping problem, without any model (illustrative numbers)", fontsize=15,
                 fontweight="bold", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    png = workspace_twin(FIGS) / "brie_sand_the_idea.png"
    save(fig, png, dpi=200, close=True)
    record_caption(png, (
        "**Tests:** nothing; the idea behind BRIE's sand error, without the model. **How to read:** two "
        "neighbouring sections of coast and the boundary between them. Every handover of sand across the boundary "
        "is written down twice, once by each section. Left: both write the same amount, so the island's total is "
        "unchanged. Middle: BRIE lets each section work out the amount with its own rate; when the two rates "
        "are the same the books still agree. Right: when the rates differ, one side books more than the other, "
        "and the difference is sand made or lost. The numbers are illustrative. **Shows:** sand can only be made "
        "or lost if the two sides of a boundary use different rates, which is what BRIE does at a bend "
        "(../3-one-year/)."))

def main():
    apply_style()
    FIGS.mkdir(parents=True, exist_ok=True)
    fig_idea()


if __name__ == "__main__":
    main()
