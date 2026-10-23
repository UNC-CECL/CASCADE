#!/usr/bin/env python3
"""
The fixed-shift groin and the wave-driven groin as maps, year by year, on coasts of many angles.

    python straight_coast_fixed_shift_vs_wave_driven_map_gif.py   ->  figures/straight_coast_fixed_shift_vs_wave_driven_map.gif
    python straight_coast_fixed_shift_vs_wave_driven_map_gif.py --module fixed-shift
                                                 ->  ../5-fixed-shift/figures/straight_coast_fixed_shift_map.gif
    python straight_coast_fixed_shift_vs_wave_driven_map_gif.py --module wave-driven
                                                 ->  ../6-wave-driven/figures/straight_coast_wave_driven_map.gif

The animated companion of straight_coast_fixed_shift_vs_wave_driven_map_figure.py: true-scale, north-up maps of each coast at
its own angle, from -30 to +30 degrees, the fixed-shift groin on the top row and the wave-driven groin below, each year
for 50 years. Each map shows the straight start coast, the shoreline that year, the sand built out (gold) and
cut back (red) since the start, the groin and the net drift (arrow length scaled to its size on that coast),
under one wave climate, printed in the header. Reruns both modules in straight_coast_test.py's
emulator of BRIE's alongshore step to keep every year (about a minute). --module animates one version alone,
one row. The GIF and its year-50 preview are git-ignored; rerun to make them.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-16
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import argparse
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
from map_geometry import draw_coast, shoreline  # noqa: E402
from site_layer.hat_figure_style import C, C_1984_FILL, INK, INK_MUTED, apply_style  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = workspace_twin(HERE) / "figures" / "straight_coast_fixed_shift_vs_wave_driven_map.gif"
VERSIONS = (("dipole", "Fixed-shift groin", "moves a fixed 12 m a year,\nwhatever the waves do"),
            ("corrected", "Wave-driven groin", "holds back 60% of the drift\narriving, sand balanced"))
ALONE = {                            # --module: the version animated alone, and its own folder
    "fixed-shift": ("dipole", workspace_twin(HERE.parent) / "5-fixed-shift" / "figures" / "straight_coast_fixed_shift_map.gif"),
    "wave-driven": ("corrected", workspace_twin(HERE.parent) / "6-wave-driven" / "figures" / "straight_coast_wave_driven_map.gif"),
}
ANGLES = (-30, -20, -10, 0, 10, 20, 30)
YEARS = 50
CELLS = 3                            # cells in view either side of the groin
DRAWN = CELLS + 3                    # cells drawn, so a turned coast still fills the view
FILLS = dict(land=dict(facecolor=C["ADDED_FILL"], alpha=0.55), sea=dict(facecolor=C["WATER"], alpha=0.45),
             built=dict(facecolor=C["ADDED"], alpha=0.75), cut=dict(facecolor=C_1984_FILL, alpha=0.95))
HALF_WIDTH_M = 1150.0                # map half-width east-west about the groin, the same on every map
GROIN_DRAWN_M = 250.0                # length drawn; in the model the groin is a face between two cells
HOLD = 20                            # frames the last year is held
HEADER = (                           # the lines under the year stamp: (text, colour)
    ("Waves on every coast:  Hs {Hs:g} m  ·  Tp {Tp:g} s  ·  asymmetry {asym:g}  ·  "
     "high-angle fraction {ahf:g}".format(**sc.CLIMATE), INK),
    ("Coast angle: negative turned against the drift, positive turned with it; the drift fades toward +32° "
     "and reverses beyond", INK_MUTED),
    ("True scale, north up, 1.5 km either side of the groin  ·  each 500 m cell's shoreline moves as one  ·  "
     "BRIE's alongshore step alone", INK_MUTED),
)
# Layout in inches, so the header sits just above the maps whatever the number of rows
FIG_W = 11.0
H_HEAD, H_TITLE, H_PANEL, H_GAP, H_FOOT = 1.0, 0.22, 1.85, 0.12, 0.42
# -----------------------------------------------------------------------------


def runs(versions):
    sc.PROFILE_YEARS = tuple(range(1, YEARS + 1))            # keep every year, not just the five sampled
    t = sc.brie_tables()
    out = {}
    for name, *_ in versions:
        for a in ANGLES:
            rows, prof = sc.run(t, a, name, YEARS)
            out[name, a] = dict(prof=np.array([prof[y] for y in range(1, YEARS + 1)]),
                                drift=abs(float(sc.drift(t, a)[0])))
    return out


def main(module="both"):
    apply_style()
    if module == "both":
        versions, out = VERSIONS, OUT
    else:
        versions, out = tuple(v for v in VERSIONS if v[0] == ALONE[module][0]), ALONE[module][1]
    n = len(versions)
    data = runs(versions)
    offs = np.arange(sc.NY) - sc.LO
    keep = (offs >= -DRAWN + 1) & (offs <= DRAWN)
    qmax = max(d["drift"] for d in data.values())
    fig_h = H_FOOT + n * H_PANEL + (n - 1) * (H_GAP + H_TITLE) + H_TITLE + H_HEAD
    fig, axes = plt.subplots(n, len(ANGLES), figsize=(FIG_W, fig_h), squeeze=False)
    fig.subplots_adjust(left=0.06, right=0.995, top=1 - (H_HEAD + H_TITLE) / fig_h, bottom=H_FOOT / fig_h,
                        wspace=0.04, hspace=(H_GAP + H_TITLE) / H_PANEL)
    stamp = fig.text(0.5, 1 - 0.1 / fig_h, "", ha="center", va="top", fontsize=13, fontweight="bold")
    for k, (line, col) in enumerate(HEADER):
        fig.text(0.5, 1 - (0.46 + 0.17 * k) / fig_h, line, ha="center", va="top", fontsize=8.5 if k == 0 else 8,
                 color=col, fontweight="bold" if k == 0 else "normal")
    for i, (name, title, rule) in enumerate(versions):
        box = axes[i, 0].get_position()
        y = (box.y0 + box.y1) / 2
        fig.text(0.011, y, title, rotation=90, ha="center", va="center", fontsize=9.5, fontweight="bold",
                 color=C["ADDED"] if name == "dipole" else INK)
        fig.text(0.034, y, rule, rotation=90, ha="center", va="center", fontsize=7, linespacing=1.3, color=INK_MUTED)
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
        for i, (name, *_) in enumerate(versions):
            for j, a in enumerate(ANGLES):
                ax = axes[i, j]
                ax.clear()
                d = data[name, a]
                x0, y0 = shoreline(a, offs[keep], np.zeros(keep.sum()), sc.DY_M)
                x1, y1 = shoreline(a, offs[keep], d["prof"][yr][keep], sc.DY_M)
                xmin, xmax = -HALF_WIDTH_M, HALF_WIDTH_M
                draw_coast(ax, a, offs[keep], d["prof"][yr][keep], sc.DY_M, **FILLS)
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
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(-CELLS * sc.DY_M, CELLS * sc.DY_M)
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])
                if i == 0:
                    ax.set_title(f"{a:+d}°" if a else "0°", fontsize=9)
                if i == n - 1 and j == 0:
                    ax.plot([-HALF_WIDTH_M + 100, -HALF_WIDTH_M + 600], [CELLS * sc.DY_M - 200] * 2, color=INK, lw=2, zorder=8)
                    ax.text(-HALF_WIDTH_M + 350, CELLS * sc.DY_M - 160, "500 m", ha="center", fontsize=6.5, zorder=8)
        return []

    anim = FuncAnimation(fig, draw, frames=YEARS + HOLD, blit=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out, writer=PillowWriter(fps=6), dpi=100)
    draw(YEARS - 1)
    fig.savefig(out.with_suffix(".png").with_name(out.stem + "_year50_preview.png"), dpi=100)
    plt.close(fig)
    print(f"wrote {out.relative_to(sc.REPO)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--module", choices=("both", "fixed-shift", "wave-driven"), default="both")
    main(ap.parse_args().module)
