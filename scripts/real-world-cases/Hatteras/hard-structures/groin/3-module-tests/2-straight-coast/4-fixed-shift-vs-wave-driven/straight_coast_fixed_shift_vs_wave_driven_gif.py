#!/usr/bin/env python3
"""
The fixed-shift groin and the wave-driven groin side by side, year by year, on two straight coasts.

    python straight_coast_fixed_shift_vs_wave_driven_gif.py   ->  figures/straight_coast_fixed_shift_vs_wave_driven.gif
    python straight_coast_fixed_shift_vs_wave_driven_gif.py --module fixed-shift
                                                 ->  ../5-fixed-shift/figures/straight_coast_fixed_shift.gif
    python straight_coast_fixed_shift_vs_wave_driven_gif.py --module wave-driven
                                                 ->  ../6-wave-driven/figures/straight_coast_wave_driven.gif

Reruns both modules in straight_coast_test.py's emulator of BRIE's alongshore step (about 20 s), keeping
every year's shoreline, and animates them as plan views: land below the shoreline, ocean above, the
groin in red, and each panel's running total of sand gained or lost. --module animates one version alone.
The GIF is git-ignored; rerun to make it.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, figsize  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = workspace_twin(HERE) / "figures" / "straight_coast_fixed_shift_vs_wave_driven.gif"
PAIR = (("dipole", "Fixed-shift groin",
         f"moves a fixed {sc.M:g} m a year north,\nwhatever the waves do"),
        ("corrected", "Wave-driven groin",
         f"holds back {sc.BLOCKING_FRACTION:.0%} of the drift arriving,\nand puts back BRIE's lost or made sand"))
ANGLES = ((-20, "coast turned against the drift (−20°)"), (10, "coast turned with the drift (+10°)"))
ALONE = {                           # --module: the version animated alone, and its own folder
    "fixed-shift": ("dipole", workspace_twin(HERE.parent) / "5-fixed-shift" / "figures" / "straight_coast_fixed_shift.gif"),
    "wave-driven": ("corrected", workspace_twin(HERE.parent) / "6-wave-driven" / "figures" / "straight_coast_wave_driven.gif"),
}
YEARS = 50
CELLS = 12                          # cells either side of the groin drawn
HOLD = 20                           # frames the last year is held
# -----------------------------------------------------------------------------


# Every year's shoreline (seaward +, updrift at positive offsets) and the reach's total sand change
def runs(pair):
    sc.PROFILE_YEARS = tuple(range(1, YEARS + 1))            # keep every year, not just the five sampled
    t = sc.brie_tables()
    out = {}
    for name, *_ in pair:
        for a, _ in ANGLES:
            rows, prof = sc.run(t, a, name, YEARS)
            out[name, a] = dict(
                prof=np.array([prof[y] for y in range(1, YEARS + 1)]),
                total=rows.reach_change_m.to_numpy())
    return out


def main(module="both"):
    apply_style()
    if module == "both":
        pair, out = PAIR, OUT
    else:
        pair, out = tuple(p for p in PAIR if p[0] == ALONE[module][0]), ALONE[module][1]
    data = runs(pair)
    offs = np.arange(sc.NY) - sc.LO                       # as the tables' domain_from_groin; the groin face is at +0.5
    keep = (offs >= -CELLS) & (offs <= CELLS + 1)
    lim = max(np.abs(d["prof"][:, keep]).max() for d in data.values()) * 1.1
    n = len(pair)
    fig, axes = plt.subplots(2, n, figsize=figsize("double" if n == 2 else "single", height=5.4), sharex=True,
                             sharey=True, squeeze=False)
    left = 0.1 if n == 2 else 0.18
    fig.subplots_adjust(left=left, right=0.98, top=0.82, bottom=0.1, hspace=0.32, wspace=0.08)
    stamp = fig.text(0.5, 0.975, "", ha="center", va="top", fontsize=12, fontweight="bold")
    for j, (name, title, rule) in enumerate(pair):
        xc = left + (0.98 - left) * (j + 0.5) / n
        fig.text(xc, 0.925, title, ha="center", fontsize=9.5, fontweight="bold", color=sc.MODULES[name]["col"])
        fig.text(xc, 0.905, rule, ha="center", va="top", fontsize=7, color=INK_MUTED)

    def draw(k):
        yr = min(k, YEARS - 1)
        stamp.set_text(f"Year {yr + 1}")
        for i, (a, label) in enumerate(ANGLES):
            for j, (name, *_) in enumerate(pair):
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
    out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out, writer=PillowWriter(fps=5), dpi=110)
    plt.close(fig)
    print(f"wrote {out.relative_to(sc.REPO)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--module", choices=("both", "fixed-shift", "wave-driven"), default="both")
    main(ap.parse_args().module)
