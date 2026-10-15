#!/usr/bin/env python3
"""
Where the wave-driven groin sits in a CASCADE model year, and what it does to the numbers.

    python wave_driven_in_the_model_year_figure.py [--theta -20] [--year 10]
        ->  figures/wave_driven_1_in_the_model_year.png

The companion of ../1-fixed-shift/figures/fixed_shift_1_in_the_model_year.png, panel for panel. (a) the order within
one year: Barrier3D in every domain, the groin's edit, BRIE's alongshore step. (b) the edit's two parts in
one year of the straight-coast run (straight_coast_test.py in ../../3-module-tests/2-straight-coast/1-run/:
BRIE's alongshore step alone, option A waves, strength 0.6), on top of a made-up uniform Barrier3D retreat:
the drift held back across the groin, and the leak correction that cancels the sand BRIE's step would make
or lose at the groin's bends. The yearly loop repeats straight_coast_test.run() for the wave-driven groin so
the two parts can be drawn apart; it is checked against run()'s own profile. Reruns one coast (seconds).

The wave-driven groin is cascade/groin.py CorrectedGroinCallback, the runner's groin since 2026-10-13. Like
the fixed-shift groin it sits on the sandbox Cascade's `_groin_callback` (cascade/cascade_groin.py) and edits x_s_dt,
the shoreline change each domain hands to BRIE. On a straight coast the groin-free shoreline has no bends,
so the runner's footprint-based correction equals the rig's leak_correction().

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

RUN = Path(__file__).resolve().parents[2] / "3-module-tests" / "2-straight-coast" / "1-run"
sys.path.insert(0, str(RUN))   # the straight-coast test
import straight_coast_test as sc  # noqa: E402  (puts scripts/ on the path)
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
OUT = HERE / "figures" / "wave_driven_1_in_the_model_year.png"
GIS = np.arange(2, 10)               # domains drawn; the groin is the face between GIS 5 (sc.LO) and GIS 6 (sc.HI)
EXAMPLE_BE = 1.5                     # m/yr of retreat in every domain, the illustrative Barrier3D change in (b)
HELD_C, CORR_C = C["ADDED"], C["LATE"]
# -----------------------------------------------------------------------------


def box(ax, x, y, w, h, text, fc, ec=INK, fs=7.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=INK, linespacing=1.3)


def arrow(ax, x0, y0, x1, y1, col=INK, lw=1.4):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=col, lw=lw))


# The groin's edit in `year`, split into the drift held back and the leak correction (seaward +), as run() does it
def edit_parts(t, theta0, year):
    slope = np.tan(np.deg2rad(theta0))
    eta = np.zeros(sc.NY)
    up = sc.LO if float(sc.drift(t, theta0)[0]) > 0 else sc.HI
    for y in range(1, year + 1):
        step = slope * sc.DY_M + (np.roll(eta, -1) - eta)
        theta = 180.0 * np.arctan2(step, sc.DY_M) / np.pi
        r = np.maximum(0.0, t["coast_diff"][sc.angle_index(theta, t["climl"])] * sc.DT_YR / 2 / sc.DY_M ** 2)
        approach = theta[sc.HI] if up == sc.HI else theta[sc.LO - 1]
        q = sc.BLOCKING_FRACTION * float(sc.drift(t, approach)[0]) * sc.DT_YR / (sc.DY_M * t["depth"])
        held = np.zeros(sc.NY)
        held[sc.LO] -= q
        held[sc.HI] += q
        corr = sc.leak_correction(r, step)
        if y == year:
            parts = (-held, -corr)
        rhs = eta + r * (np.roll(eta, -1) - 2 * eta + np.roll(eta, 1)) + held + corr
        dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
        eta = spsolve(csr_matrix((dv, (t["di"], t["dj"]))), rhs)
    return parts, -eta


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--theta", type=int, default=-20)
    ap.add_argument("--year", type=int, default=10)
    a = ap.parse_args()
    apply_style()
    t = sc.brie_tables()
    (held, corr), end = edit_parts(t, a.theta, a.year)
    sc.PROFILE_YEARS = (a.year,)
    _, prof = sc.run(t, a.theta, "corrected", a.year)
    assert np.allclose(end, prof[a.year], atol=1e-6), "the loop does not match straight_coast_test.run()"
    idx = sc.LO + (GIS - 5)
    held, corr = held[idx], corr[idx]
    leak = -float(corr.sum())                         # what BRIE's step would have made (+) or lost (−)
    q = float(held[GIS == 6][0])

    fig = plt.figure(figsize=figsize("double", height=4.6), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=(1.25, 1))
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    box(ax, 0.2, 7.3, 8.6, 2.2, "1. Barrier3D runs in every domain\nstorms, overwash, dunes and its "
        "background erosion (BE)\n→ each domain's shoreline change this year", "#eef3f8")
    box(ax, 0.2, 3.6, 8.6, 3.0, "2. The wave-driven groin edits that change, in two parts\n"
        "held back: a share of the drift arriving at the groin,\nread from the waves each year, "
        "moved GIS 5 → GIS 6\nleak correction: within 10 domains of the groin, cancels\nthe sand "
        "BRIE's step would make or lose at its bends", "#fbeed3", ec=HELD_C)
    box(ax, 0.2, 0.5, 8.6, 2.4, "3. BRIE moves the shoreline alongshore, unchanged\nit smooths the edited "
        "shoreline between neighbouring domains;\nwith step 2's correction the island's sand "
        "balances", "#eef3f8")
    arrow(ax, 4.5, 7.3, 4.5, 6.65)
    arrow(ax, 4.5, 3.6, 4.5, 2.95)
    ax.annotate("", xy=(8.9, 8.4), xytext=(8.9, 1.7),
                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=1.0, connectionstyle="arc3,rad=0.45"))
    ax.text(9.95, 5.1, "next\nyear", fontsize=6.8, color=INK_MUTED, va="center", ha="right")
    _letter_inside(ax, 0)
    ax.set_title("One model year, in order", loc="center", fontsize=8.5)

    ax = fig.add_subplot(gs[0, 1])
    before = -np.full(GIS.size, EXAMPLE_BE)
    after = before + held + corr
    w = 0.21
    ax.bar(GIS - 1.5 * w, before, width=w, color="0.75", label="from Barrier3D (example)")
    ax.bar(GIS - 0.5 * w, held, width=w, color=HELD_C, label="held back")
    ax.bar(GIS + 0.5 * w, corr, width=w, color=CORR_C, label="leak correction")
    ax.bar(GIS + 1.5 * w, after, width=w, color="0.35", label="to BRIE, after the groin")
    lim = 1.3 * max(np.abs(held).max(), np.abs(after).max())
    ax.axvline(5.5, color=C["GROIN"], lw=2.2)
    ax.text(5.6, -lim * 0.97, "groin (GIS 5/6)", color=C["GROIN"], fontsize=6.8, va="bottom")
    ax.annotate(f"+{q:.0f} m", xy=(6 - 0.5 * w, q), xytext=(6.9, q * 0.55), fontsize=7, color=HELD_C,
                arrowprops=dict(arrowstyle="-", color=HELD_C, lw=0.6))
    ax.annotate(f"−{q:.0f} m", xy=(5 - 0.5 * w, -q), xytext=(2.6, -q * 1.02), fontsize=7, color=HELD_C,
                arrowprops=dict(arrowstyle="-", color=HELD_C, lw=0.6))
    ax.text(9.6, -lim * 0.62, (f"correction sums to {-leak:+.0f} m:\nBRIE's step would make\n{leak:+.0f} m "
            "at the bends").replace("-", "−"), fontsize=6.6, color=CORR_C, ha="right", va="top")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xticks(GIS, [f"GIS {g}" for g in GIS], fontsize=6.5)
    ax.set_xlabel("south ← domain → north")
    ax.set_ylabel("Shoreline change this year (m)\n+ seaward")
    ax.set_ylim(-lim, 1.35 * lim)                     # headroom for the legend
    ax.legend(frameon=False, fontsize=6.4, loc="upper right")
    ax.grid(axis="y")
    open_frame(ax)
    _letter_inside(ax, 1)
    ax.set_title(f"What step 2 does (year {a.year}, coast at {a.theta:+d}°)".replace("-", "−"),
                 loc="center", fontsize=8.5)
    fig.suptitle("Where the wave-driven groin sits in the model", fontsize=10)
    save(fig, OUT, dpi=250, close=True)
    record_caption(OUT, (
        "**Tests:** nothing new; this figure explains where the wave-driven groin ("
        "`CorrectedGroinCallback`, the runner's groin since 2026-10-13) acts in a CASCADE model year. The "
        "companion of `../1-fixed-shift/figures/fixed_shift_1_in_the_model_year.png`. **How to read:** (a) the order "
        "within one year: every Barrier3D domain runs first, including its background erosion (BE), and hands "
        "BRIE its shoreline change for the year; the groin then edits that change, in the same place the fixed-shift groin "
        "did (`_groin_callback` on the sandbox Cascade); BRIE then moves the shoreline alongshore, unchanged. "
        f"(b) the edit in year {a.year} of the straight-coast run with the coast at {a.theta:+d}° (BRIE's "
        f"alongshore step alone, option A waves, strength {sc.BLOCKING_FRACTION:g}), drawn on GIS labels as if "
        f"it were Buxton: grey, a uniform {EXAMPLE_BE:g} m of retreat from Barrier3D in every domain (made up "
        f"for the example); gold, the drift held back, {q:.1f} m moved from GIS 5 to GIS 6, read from BRIE's "
        "own drift at the angle of the face the drift arrives on; blue, the leak correction near the groin, "
        f"which sums to {-leak:+.1f} m and cancels the {leak:+.1f} m BRIE's step would otherwise make at the "
        "groin's bends (zero in year 1, when there are no bends yet); dark, what BRIE receives. At Buxton the "
        "strength is 0.16 (1996-2009) or 0.10 (2009-2017), so the amounts are smaller. **Shows:** unlike the "
        "fixed-shift groin's 12 m, the amount held back follows the waves, and the correction keeps the island's "
        "sand constant, so the groin moves sand but never makes or loses it."))
    print("wrote", OUT, f"| held {q:.2f} m | correction {-leak:+.2f} m")


if __name__ == "__main__":
    main()
