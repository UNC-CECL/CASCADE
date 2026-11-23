"""
Schematic of the three collision designs on two idealized dune rows (2026-11-19), drawn as the model holds them:
two 10 m dune cells on the berm, back row landward, front row seaward, ocean on the right. Pre-storm rows are the
grey dashed outline, post-storm rows filled in the design's colour; panels stacked top to bottom.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_collision_design_schematic.py

- Lower both rows (red): both rows lowered in proportion to their height (cascade/dune_collision.py "proportional")
- Lower front row (yellow): only the front row lowered ("front")
- Face retreat (blue): the face retreats landward through the front row ("lateral", erode_laterally); the row
  keeps its height until the sand it holds is used up, then drops to the toe; the back row is next
Writes collision_designs_schematic.png in 3-collision-design/0-all-designs/ (collision_designs.py).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
DESIGNS = (("Lower both rows", "#D7301F", "Both rows lowered in proportion to their height"),
           ("Lower front row", "#F2C511", "Only the front row lowered"),
           ("Face retreat", "#2166AC", "Face retreats landward through the front row;\n"
                                       "the row keeps its height until its sand is gone"))
BERM = 1.0                     # m, berm above sea level
HEIGHT = 3.0                   # m, pre-storm dune height above the berm, both rows
LOWER = 0.55                   # fraction of height kept after lowering
ROW = 10.0                     # m, one dune cell cross-shore
BACK_X = 25.0                  # m, landward edge of the back row; front row right behind it
RETREAT = 6.0                  # m, face retreat into the front row
X0, XMAX = 0.0, 65.0             # m, drawn cross-shore extent
OUT = CD.ALL_DIR / "collision_designs_schematic.png"
# -----------------------------------------------------------------------------


def main():
    front_x = BACK_X + ROW
    toe = front_x + ROW
    fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True)
    for ax, (name, col, note) in zip(axes, DESIGNS):
        # ground: berm, then beach sloping to the sea; ocean below sea level
        ax.fill([X0, toe, XMAX, XMAX, X0], [BERM, BERM, 0, 0, 0], color="0.88", lw=0)
        ax.plot([X0, toe, XMAX], [BERM, BERM, 0], color="k", lw=1.0)
        ax.axhspan(-0.6, 0, color="#BFD7EA", lw=0)
        ax.text(XMAX - 1, -0.3, "ocean", ha="right", va="center", fontsize=8, color="#2166AC")

        back_h = HEIGHT * LOWER if name == "Lower both rows" else HEIGHT
        front_h = HEIGHT * LOWER if name != "Face retreat" else HEIGHT
        front_w = ROW - RETREAT if name == "Face retreat" else ROW
        ax.add_patch(Rectangle((BACK_X, BERM), ROW, back_h, fc=col, ec="k", lw=1.0))
        ax.add_patch(Rectangle((front_x, BERM), front_w, front_h, fc=col, ec="k", lw=1.0))
        if name == "Face retreat":   # what the face has eroded so far, still counted in the row's height
            ax.add_patch(Rectangle((front_x + front_w, BERM), RETREAT, HEIGHT, fc="none", ec=col,
                                   hatch="///", lw=0))
            y = BERM + HEIGHT * 0.5
            ax.annotate("", xy=(front_x + front_w + 0.3, y), xytext=(toe + 7, y),
                        arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
        else:
            xs = (BACK_X, front_x) if name == "Lower both rows" else (front_x,)
            for x0 in xs:
                ax.annotate("", xy=(x0 + ROW / 2, BERM + HEIGHT * LOWER + 0.05),
                            xytext=(x0 + ROW / 2, BERM + HEIGHT - 0.05),
                            arrowprops=dict(arrowstyle="->", color="k", lw=1.4))
        for x0 in (BACK_X, front_x):  # pre-storm rows
            ax.add_patch(Rectangle((x0, BERM), ROW, HEIGHT, fc="none", ec="0.45", lw=1.3, ls=(0, (4, 2))))
        ax.text(BACK_X + ROW / 2, BERM + HEIGHT + 0.2, "back\nrow", ha="center", va="bottom", fontsize=9)
        ax.text(front_x + ROW / 2, BERM + HEIGHT + 0.2, "front\nrow", ha="center", va="bottom", fontsize=9)

        ax.set_title(name, fontsize=12, fontweight="bold", loc="left")
        ax.text(0.0, 0.93, note, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_xlim(X0, XMAX)
        ax.set_ylim(-0.6, BERM + HEIGHT + 3.0)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    axes[0].plot([], [], color="0.45", lw=1.3, ls=(0, (4, 2)), label="Before storm")
    axes[0].legend(loc="upper right", frameon=False, fontsize=8)
    axes[-1].set_xlabel("landward  ←            →  seaward", fontsize=9)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(OUT)


if __name__ == "__main__":
    main()
