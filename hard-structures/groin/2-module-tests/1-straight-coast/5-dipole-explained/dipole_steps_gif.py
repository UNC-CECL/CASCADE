#!/usr/bin/env python3
"""
The steps figure animated: the model year by year, step by step, from a straight coast, for either groin.

    python dipole_steps_gif.py                       ->  figures/dipole_3_steps.gif
    python dipole_steps_gif.py --module corrected    ->  ../6-corrected-explained/figures/corrected_groin_steps.gif
    options: --years 50 --theta -20 --slow 5

The animated companion of dipole_steps_figure.py. Top: the steps of one call of Cascade.update() (Barrier3D,
groin hook, BRIE, management), the active one lit. Below: the shoreline around the groin, the straight start
dashed, and the running total of sand made or lost over the 61-cell reach. The first `slow` years show three
frames each (after Barrier3D, after the groin hook, after BRIE); later years one frame each, after BRIE.

The yearly loop repeats straight_coast_test.run() for the two modules so the state after the groin's edit can
be drawn; it is checked against run()'s own profiles at the last year. The GIFs are git-ignored; rerun to make
them (under a minute each).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.linalg import spsolve  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "1-run"))   # the test, in ../1-run/
import straight_coast_test as sc  # noqa: E402  (puts scripts/ on the path)
from site_layer.hat_figure_style import C, INK, INK_MUTED, apply_style, open_frame  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
CELLS = 5                            # cells drawn either side of the groin
OUT = {"dipole": HERE / "figures" / "dipole_3_steps.gif",
       "corrected": HERE.parent / "6-corrected-explained" / "figures" / "corrected_groin_steps.gif"}
NAME = {"dipole": "the dipole groin (version 1)", "corrected": "the corrected groin (version 6)"}
HOOK = {"dipole": f"moves a fixed {sc.M:g} m\nGIS 5 → GIS 6",
        "corrected": f"holds back {sc.BLOCKING_FRACTION:.0%} of the drift\narriving; cancels BRIE's leak"}
BOX_C = ("#eef3f8", "#fbeed3", "#e6f0e9", "#f2f2f2")
HOLD = 18                            # frames the last year is held
BALANCED_M = 30.0
# -----------------------------------------------------------------------------


# Each year's shoreline (seaward +) at the start, after the groin's edit, and after BRIE; as straight_coast_test.run()
def history(t, theta0, module, years):
    slope = np.tan(np.deg2rad(theta0))
    eta = np.zeros(sc.NY)
    up, down = (sc.LO, sc.HI) if float(sc.drift(t, theta0)[0]) > 0 else (sc.HI, sc.LO)
    out = []
    for _ in range(years):
        step = slope * sc.DY_M + (np.roll(eta, -1) - eta)
        theta = 180.0 * np.arctan2(step, sc.DY_M) / np.pi
        r = np.maximum(0.0, t["coast_diff"][sc.angle_index(theta, t["climl"])] * sc.DT_YR / 2 / sc.DY_M ** 2)
        dx = np.zeros(sc.NY)
        if module == "dipole":
            dx[up] -= sc.M
            dx[down] += sc.M
        else:
            approach = theta[sc.HI] if up == sc.HI else theta[sc.LO - 1]
            q = sc.BLOCKING_FRACTION * float(sc.drift(t, approach)[0]) * sc.DT_YR / (sc.DY_M * t["depth"])
            dx[sc.LO] -= q
            dx[sc.HI] += q
            dx += sc.leak_correction(r, step)
        rhs = eta + r * (np.roll(eta, -1) - 2 * eta + np.roll(eta, 1)) + dx
        dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
        new = spsolve(csr_matrix((dv, (t["di"], t["dj"]))), rhs)
        out.append((-eta.copy(), -(eta + dx), -new.copy()))
        eta = new
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--module", choices=("dipole", "corrected"), default="dipole")
    ap.add_argument("--years", type=int, default=50)
    ap.add_argument("--theta", type=int, default=-20)
    ap.add_argument("--slow", type=int, default=5, help="years shown step by step")
    a = ap.parse_args()
    apply_style()
    t = sc.brie_tables()
    hist = history(t, a.theta, a.module, a.years)
    sc.PROFILE_YEARS = (a.years,)
    _, prof = sc.run(t, a.theta, a.module, a.years)
    assert np.allclose(hist[-1][2], prof[a.years], atol=1e-6), "the loop does not match straight_coast_test.run()"
    offs = np.arange(sc.NY) - sc.LO
    keep = (offs >= -CELLS + 1) & (offs <= CELLS)
    x = offs[keep]
    xs = np.repeat(np.r_[x - 0.5, x[-1] + 0.5], 2)[1:-1]
    lim = 1.25 * max(np.abs(s[keep]).max() for h in hist for s in h)
    # frames: (year index, stage) with stage 0 after Barrier3D, 1 after the groin, 2 after BRIE
    frames = [(y, s) for y in range(min(a.slow, a.years)) for s in (0, 1, 2)]
    frames += [(y, 2) for y in range(a.slow, a.years)] + [(a.years - 1, 2)] * HOLD

    fig = plt.figure(figsize=(8.6, 6.0), dpi=95)
    gs = fig.add_gridspec(2, 1, height_ratios=(0.26, 1), left=0.09, right=0.98, top=0.9, bottom=0.08, hspace=0.06)
    top, ax = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    fig.suptitle(f"One model year at a time with {NAME[a.module]}, from a straight coast (at {a.theta:+d}°)",
                 fontsize=10.5)
    stamp = fig.text(0.09, 0.925, "", fontsize=10, fontweight="bold", va="bottom")
    titles = ("1  Barrier3D", "2  Groin hook", "3  BRIE", "4  Management")
    bodies = ("storms, dunes, BE\n(none in this test)", HOOK[a.module], "smooths the coast\nalongshore",
              "roads, sandbags,\nnourishment")

    def draw(k):
        y, s = frames[k]
        start, edited, after = hist[y]
        shown = (start, edited, after)[s]
        before = (start, start, edited)[s]
        stamp.set_text(f"Year {y + 1}" + ("" if y >= a.slow else f": after step {s + 1}"))
        top.clear()
        top.set_xlim(0, 4)
        top.set_ylim(0, 1)
        top.axis("off")
        for i in range(4):
            lit = i == s
            top.add_patch(FancyBboxPatch((i + 0.06, 0.06), 0.82, 0.86, boxstyle="round,pad=0.01,rounding_size=0.05",
                                         fc=BOX_C[i], ec=C["ADDED"] if lit else "0.7", lw=2.4 if lit else 0.8))
            top.text(i + 0.47, 0.8, titles[i], ha="center", va="top", fontsize=8.5, fontweight="bold",
                     color=INK if lit else INK_MUTED)
            top.text(i + 0.47, 0.5, bodies[i], ha="center", va="top", fontsize=7, color=INK if lit else INK_MUTED)
            if i < 3:
                top.annotate("", xy=(i + 1.05, 0.5), xytext=(i + 0.9, 0.5),
                             arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.0))
        ax.clear()
        ya = np.repeat(shown[keep], 2)
        ax.fill_between(xs, -lim * 2, ya, color=C["ADDED_FILL"], lw=0)
        ax.fill_between(xs, ya, lim * 2, color=C["WATER"], alpha=0.45, lw=0)
        ax.axhline(0, color=INK_MUTED, lw=1.0, ls=(0, (4, 2)))
        if y < a.slow and s > 0:
            ax.plot(xs, np.repeat(before[keep], 2), color=C["ADDED"], lw=1.0, ls=(0, (2, 1.5)))
        ax.plot(xs, ya, color=INK, lw=1.8)
        ax.plot([0.5, 0.5], [-lim * 0.8, lim * 0.9], color=C["GROIN"], lw=3, solid_capstyle="butt")
        ax.text(0.62, lim * 0.86, "groin", color=C["GROIN"], fontsize=8, va="top")
        ax.text(-CELLS + 0.7, -lim * 0.92, "← GIS 5 side (south, downdrift)", fontsize=7.5, color=INK_MUTED)
        ax.text(CELLS + 0.3, -lim * 0.92, "GIS 6 side (north, updrift) →", fontsize=7.5, color=INK_MUTED, ha="right")
        tot = float(after.sum()) if s == 2 else float(start.sum())
        ok = abs(tot) <= BALANCED_M
        ax.text(-CELLS + 0.7, lim * 0.88, f"sand made or lost so far: {tot:+,.0f} m" + (" (balanced)" if ok else ""),
                fontsize=8.5, fontweight="bold", color=C["REF"] if ok else C["GROIN"], va="top",
                bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.85))
        ax.text(CELLS + 0.3, lim * 0.72, "dashed grey: the straight start" +
                ("\ndashed gold: before this step" if y < a.slow and s > 0 else ""), fontsize=7, color=INK_MUTED,
                ha="right", va="top")
        ax.set_xlim(-CELLS + 0.5, CELLS + 0.5)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([])
        ax.set_xlabel("each step is one 500 m domain; each domain's shoreline moves as one", fontsize=7.5,
                      color=INK_MUTED)
        ax.set_ylabel("Shoreline (m), + seaward")
        open_frame(ax)
        return []

    from PIL import Image
    imgs, durations = [], []
    n_main = len(frames) - HOLD
    for k in range(n_main):
        draw(k)
        fig.canvas.draw()
        imgs.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB"))
        y, s = frames[k]
        durations.append(1100 if y < a.slow else 140)
    durations[-1] = 3500                               # hold the last year
    plt.close(fig)
    out = OUT[a.module]
    out.parent.mkdir(parents=True, exist_ok=True)
    imgs[0].save(out, save_all=True, append_images=imgs[1:], duration=durations, loop=0, optimize=True)
    print("wrote", out, f"| {len(imgs)} frames")


if __name__ == "__main__":
    main()
