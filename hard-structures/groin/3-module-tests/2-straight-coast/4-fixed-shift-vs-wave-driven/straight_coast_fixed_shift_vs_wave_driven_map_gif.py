#!/usr/bin/env python3
"""
The fixed-shift groin and the wave-driven groin as maps, year by year, on coasts of many angles.

    python straight_coast_fixed_shift_vs_wave_driven_map_gif.py   ->  figures/straight_coast_fixed_shift_vs_wave_driven_map.gif

The animated companion of straight_coast_fixed_shift_vs_wave_driven_map_figure.py: true-scale, north-up maps of each coast at
its own angle, from -30 to +30 degrees, the fixed-shift groin on the top row and the wave-driven groin below, each year
for 50 years. Each map shows the straight start coast, the shoreline that year, the sand built out (gold) and
cut back (red) since the start, the groin, the net drift (arrow length scaled to its size on that coast) and
the running total of sand made or lost over the 61-cell reach. Reruns both modules in straight_coast_test.py's
emulator of BRIE's alongshore step to keep every year (about a minute). The GIF is git-ignored; rerun to make it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import C, C_1984_FILL, INK, INK_MUTED, apply_style  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = Path(__file__).resolve().parent / "figures" / "straight_coast_fixed_shift_vs_wave_driven_map.gif"
VERSIONS = (("dipole", "Fixed-shift groin", "moves a fixed 12 m a year, whatever the waves do"),
            ("corrected", "Wave-driven groin", "holds back 60% of the drift arriving, sand balanced"))
ANGLES = (-30, -20, -10, 0, 10, 20, 30)
YEARS = 50
CELLS = 3                            # cells drawn either side of the groin
HALF_WIDTH_M = 1150.0                # map half-width east-west about the groin, the same on every map
GROIN_DRAWN_M = 250.0                # length drawn; in the model the groin is a face between two cells
BALANCED_M = 30.0                    # the sand-balance band used in the comparison figure
HOLD = 20                            # frames the last year is held
# -----------------------------------------------------------------------------


def runs():
    sc.PROFILE_YEARS = tuple(range(1, YEARS + 1))            # keep every year, not just the five sampled
    t = sc.brie_tables()
    out = {}
    for name, *_ in VERSIONS:
        for a in ANGLES:
            rows, prof = sc.run(t, a, name, YEARS)
            out[name, a] = dict(prof=np.array([prof[y] for y in range(1, YEARS + 1)]),
                                total=rows.reach_change_m.to_numpy(),
                                drift=abs(float(sc.drift(t, a)[0])))
    return out


# A shoreline as a north-up polyline: cell k (offset from the groin, + north) spans northing (k - 1) to k
# cells, so the groin face (between offsets 0 and 1) sits at 0; each cell's shoreline moves as one
def shoreline(theta, offs, seaward):
    tan = np.tan(np.deg2rad(theta))
    xs, ys = [], []
    for k, s in zip(offs, seaward):
        y0, y1 = (k - 1) * sc.DY_M, k * sc.DY_M
        xs += [-tan * y0 + s, -tan * y1 + s]
        ys += [y0, y1]
    return np.array(xs), np.array(ys)


def main():
    apply_style()
    data = runs()
    offs = np.arange(sc.NY) - sc.LO
    keep = (offs >= -CELLS + 1) & (offs <= CELLS)
    qmax = max(d["drift"] for d in data.values())
    fig, axes = plt.subplots(2, len(ANGLES), figsize=(11.0, 6.4))
    fig.subplots_adjust(left=0.06, right=0.995, top=0.84, bottom=0.08, wspace=0.04, hspace=0.14)
    stamp = fig.text(0.5, 0.985, "", ha="center", va="top", fontsize=13, fontweight="bold")
    fig.text(0.5, 0.94, "Coast angle: negative turned against the drift, positive turned with it; the drift "
             "fades toward +32° and reverses beyond. True scale, north up, 1.5 km either side of the groin; each 500 m "
             "cell's shoreline moves as one. BRIE's alongshore step alone, option A waves.", ha="center", wrap=True, va="top", fontsize=8.5, color=INK_MUTED)
    for i, (name, title, rule) in enumerate(VERSIONS):
        y = 0.83 - i * 0.375 - 0.18
        fig.text(0.012, y, title, rotation=90, ha="center", va="center", fontsize=9.5, fontweight="bold",
                 color=C["ADDED"] if name == "dipole" else INK)
        fig.text(0.032, y, rule, rotation=90, ha="center", va="center", fontsize=7, color=INK_MUTED)
    fig.legend(handles=[Line2D([], [], color=INK_MUTED, lw=1.0, ls=(0, (4, 2)), label="start: a straight coast"),
                        Line2D([], [], color=INK, lw=1.5, label="shoreline this year"),
                        Patch(facecolor=C["ADDED"], alpha=0.75, label="built out since the start"),
                        Patch(facecolor=C_1984_FILL, alpha=0.95, label="cut back since the start"),
                        Line2D([], [], color=C["GROIN"], lw=2.6, label="groin"),
                        Line2D([], [], color=INK_MUTED, lw=1.2, marker=">", ms=4, label="net drift (length = size)")],
               loc="lower center", ncol=6, frameon=False, fontsize=7.5)

    def draw(k):
        yr = min(k, YEARS - 1)
        stamp.set_text(f"Year {yr + 1}")
        for i, (name, *_) in enumerate(VERSIONS):
            for j, a in enumerate(ANGLES):
                ax = axes[i, j]
                ax.clear()
                d = data[name, a]
                x0, y0 = shoreline(a, offs[keep], np.zeros(keep.sum()))
                x1, y1 = shoreline(a, offs[keep], d["prof"][yr][keep])
                xmin, xmax = -HALF_WIDTH_M, HALF_WIDTH_M
                ax.fill_betweenx(y1, xmin, x1, color=C["ADDED_FILL"], alpha=0.55, lw=0, zorder=0)
                ax.fill_betweenx(y1, x1, xmax, color=C["WATER"], alpha=0.45, lw=0, zorder=0)
                ax.fill_betweenx(y0, x0, x1, where=x1 >= x0, color=C["ADDED"], alpha=0.75, lw=0, zorder=1)
                ax.fill_betweenx(y0, x1, x0, where=x1 < x0, color=C_1984_FILL, alpha=0.95, lw=0, zorder=1)
                ax.plot(x0, y0, color=INK_MUTED, lw=0.9, ls=(0, (4, 2)), zorder=3)
                ax.plot(x1, y1, color=INK, lw=1.4, zorder=4)
                nx, ny = np.cos(np.deg2rad(a)), np.sin(np.deg2rad(a))
                ax.plot([-50 * nx, GROIN_DRAWN_M * nx], [-50 * ny, GROIN_DRAWN_M * ny], color=C["GROIN"], lw=2.4,
                        solid_capstyle="butt", zorder=5)
                # net drift offshore, running south along the coast, length scaled to its size here
                tan = np.tan(np.deg2rad(a))
                half = 0.9 * sc.DY_M * d["drift"] / qmax
                yc, off = 1.6 * sc.DY_M, 420.0
                if half > 20:
                    ax.annotate("", xy=(-tan * (yc - half) + off, yc - half), xytext=(-tan * (yc + half) + off, yc + half),
                                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.1, clip_on=True),
                                annotation_clip=True, zorder=6)
                tot = d["total"][yr]
                ok = abs(tot) <= BALANCED_M
                ax.text(0.5, 0.015, f"sand made/lost\n{tot:+,.0f} m", transform=ax.transAxes, ha="center",
                        va="bottom", fontsize=7, fontweight="bold", color=C["REF"] if ok else C["GROIN"],
                        bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=8)
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(-CELLS * sc.DY_M, CELLS * sc.DY_M)
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(f"{a:+d}°" if a else "0°", fontsize=9)
                if i == 1 and j == 0:
                    ax.plot([-HALF_WIDTH_M + 100, -HALF_WIDTH_M + 600], [CELLS * sc.DY_M - 200] * 2, color=INK, lw=2, zorder=8)
                    ax.text(-HALF_WIDTH_M + 350, CELLS * sc.DY_M - 160, "500 m", ha="center", fontsize=6.5, zorder=8)
        return []

    anim = FuncAnimation(fig, draw, frames=YEARS + HOLD, blit=False)
    anim.save(OUT, writer=PillowWriter(fps=6), dpi=100)
    draw(YEARS - 1)
    fig.savefig(OUT.with_suffix(".png").with_name(OUT.stem + "_year50_preview.png"), dpi=100)
    plt.close(fig)
    print(f"wrote {OUT.relative_to(sc.REPO)}")


if __name__ == "__main__":
    main()
