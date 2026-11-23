"""
Schematic of the two starting dune-row designs (2026-11-19), drawn as the model holds them: two 10 m dune cells on
the berm, ocean on the right, over one idealized lidar dune profile and its search window; panels stacked.

    python scripts/real-world-cases/Hatteras/dune-collision/4-hindcast/3-on-vs-off/HAT_dune_start_schematic.py

- Equal rows (green): both rows the window's maximum (dune-topo v3)
- Low back row (purple): front row the window's maximum, back row the window's median
  (input_prep/1-barrier3d-domains/1-extraction/HAT_dune_topo_maxmed.py; HAT_dune_collision.py MAXMED_START)
Writes dune_start_schematic.png in 3-collision-design/0-all-designs/ (collision_designs.py).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import collision_designs as CD  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
STARTS = (("Equal rows", "#1A9850", "Both rows set to the window's maximum", "max"),
          ("Low back row", "#6A3D9A", "Front row the window's maximum, back row its median", "median"))
BERM = 1.0                     # m, berm above sea level
ROW = 10.0                     # m, one dune cell cross-shore
BACK_X = 25.0                  # m, landward edge of the back row; front row right behind it
WIN = (10.0, 50.0)             # m, dune search window on the lidar profile
X0, XMAX = 0.0, 65.0           # m, drawn cross-shore extent
OUT = CD.ALL_DIR / "dune_start_schematic.png"
# -----------------------------------------------------------------------------


def lidar():
    """Idealized lidar dune profile above the berm, m: a tall foredune and a lower, wider back dune."""
    x = np.linspace(X0, XMAX, 1301)
    z = (3.2 * np.exp(-((x - 41) / 4.5) ** 2) + 1.4 * np.exp(-((x - 26) / 6) ** 2)
         + 0.25 * np.sin(x / 2.3) * np.exp(-((x - 33) / 12) ** 2))
    return x, np.clip(z, 0, None)


def main():
    x, z = lidar()
    inwin = (x >= WIN[0]) & (x <= WIN[1])
    zmax, zmed = float(z[inwin].max()), float(np.median(z[inwin]))
    front_x = BACK_X + ROW
    toe = front_x + ROW

    fig, axes = plt.subplots(2, 1, figsize=(6, 6.4), sharex=True)
    for ax, (name, col, note, back) in zip(axes, STARTS):
        ax.fill([X0, toe, XMAX, XMAX, X0], [BERM, BERM, 0, 0, 0], color="0.88", lw=0)
        ax.plot([X0, toe, XMAX], [BERM, BERM, 0], color="k", lw=1.0)
        ax.axhspan(-0.6, 0, color="#BFD7EA", lw=0)
        ax.text(XMAX - 1, -0.3, "ocean", ha="right", va="center", fontsize=8, color="#2166AC")

        back_h = zmax if back == "max" else zmed
        ax.add_patch(Rectangle((BACK_X, BERM), ROW, back_h, fc=col, ec="k", lw=1.0, alpha=0.85))
        ax.add_patch(Rectangle((front_x, BERM), ROW, zmax, fc=col, ec="k", lw=1.0, alpha=0.85))
        ground = np.where(x < toe, BERM, BERM - (x - toe) * BERM / (XMAX - toe))
        ax.plot(x, ground + z, color="k", lw=1.3, label="Lidar dune profile")
        ax.fill_between(WIN, BERM, BERM + zmax + 0.8, color="0.5", alpha=0.08, lw=0)
        ax.text(WIN[0] + 0.5, BERM + zmax + 0.6, "search window", va="top", fontsize=8, color="0.35")
        for v, lab in ((zmax, "max"), (zmed, "median")):
            ax.plot(WIN, [BERM + v] * 2, color="0.35", lw=1.0, ls=(0, (4, 2)))
            ax.text(WIN[0] - 0.5, BERM + v, lab, ha="right", va="center", fontsize=8, color="0.35")
        ax.text(BACK_X + ROW / 2, BERM + back_h + 0.15, "back\nrow", ha="center", va="bottom", fontsize=9)
        ax.text(front_x + ROW / 2, BERM + zmax + 0.15, "front\nrow", ha="center", va="bottom", fontsize=9)

        ax.set_title(name, fontsize=12, fontweight="bold", loc="left")
        ax.text(0.0, 0.93, note, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_xlim(X0, XMAX)
        ax.set_ylim(-0.6, BERM + zmax + 2.6)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    axes[0].legend(loc="upper right", frameon=False, fontsize=8)
    axes[-1].set_xlabel("landward  ←            →  seaward", fontsize=9)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print(OUT)


if __name__ == "__main__":
    main()
