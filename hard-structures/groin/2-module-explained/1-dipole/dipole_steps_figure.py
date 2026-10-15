#!/usr/bin/env python3
"""
One model year, step by step, from a straight coast, for the dipole groin (version 1) or the corrected groin (version 6).

    python dipole_steps_figure.py                       ->  figures/dipole_3_steps.png
    python dipole_steps_figure.py --module corrected    ->  ../2-corrected/figures/corrected_groin_steps.png
    options: --years 50 --theta -20

Top: the order of operations in one call of Cascade.update(): Barrier3D in every domain (storms, overwash,
dunes, background erosion), then the groin hook, then BRIE's alongshore transport, then the management
modules. The dipole sat on the sandbox Cascade's `_groin_callback` (cascade/cascade_groin.py, before
563a542e); the corrected groin uses the same hook. Bottom: the first year of the straight-coast run
(straight_coast_test.py in ../../3-module-tests/2-straight-coast/1-run/: BRIE's alongshore step alone), from
a straight coast, at the shoreline after each step, then the same year repeated. The groin's edit is
computed as straight_coast_test.run() computes it; the profiles are run()'s own. Reruns one coast (a few
seconds). One script for both modules so the two figures match panel for panel.

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

HERE = Path(__file__).resolve().parent
RUN = HERE.parents[1] / "3-module-tests" / "2-straight-coast" / "1-run"
sys.path.insert(0, str(RUN))   # the straight-coast test
import straight_coast_test as sc  # noqa: E402  (puts scripts/ on the path)
from site_layer.hat_figure_style import (  # noqa: E402
    C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame, record_caption, save)

# --- CONFIG ------------------------------------------------------------------
CELLS = 4                            # cells drawn either side of the groin
STEP_C = {"b3d": "#eef3f8", "groin": "#fbeed3", "brie": "#e6f0e9", "mgmt": "#f2f2f2"}
OUT = {"dipole": HERE / "figures" / "dipole_3_steps.png",
       "corrected": HERE.parent / "2-corrected" / "figures" / "corrected_groin_steps.png"}
NAME = {"dipole": "the dipole groin (version 1)", "corrected": "the corrected groin (version 6)"}
# -----------------------------------------------------------------------------


# The groin's edit to year 1's shoreline change, seaward +, as straight_coast_test.run() makes it;
# for the corrected groin also its two parts (the sand held back, the leak correction)
def year_one_edit(t, theta0, module):
    eta = np.zeros(sc.NY)
    step = np.tan(np.deg2rad(theta0)) * sc.DY_M + (np.roll(eta, -1) - eta)
    theta = np.degrees(np.arctan2(step, sc.DY_M))
    r = np.maximum(0.0, t["coast_diff"][sc.angle_index(theta, t["climl"])] * sc.DT_YR / 2 / sc.DY_M ** 2)
    up, down = (sc.LO, sc.HI) if float(sc.drift(t, theta0)[0]) > 0 else (sc.HI, sc.LO)
    dx, leak = np.zeros(sc.NY), np.zeros(sc.NY)
    if module == "dipole":
        dx[up] -= sc.M
        dx[down] += sc.M
    else:
        approach = theta[sc.HI] if up == sc.HI else theta[sc.LO - 1]
        q = sc.BLOCKING_FRACTION * float(sc.drift(t, approach)[0]) * sc.DT_YR / (sc.DY_M * t["depth"])
        dx[sc.LO] -= q
        dx[sc.HI] += q
        leak = sc.leak_correction(r, step)
        dx += leak
    return -dx, -leak


def box(ax, x, y, w, h, title, body, fc, ec=INK, lw=1.0):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12", fc=fc, ec=ec, lw=lw))
    ax.text(x + w / 2, y + h - 0.18, title, ha="center", va="top", fontsize=7.6, fontweight="bold", color=INK)
    ax.text(x + w / 2, y + h - 0.62, body, ha="center", va="top", fontsize=6.6, color=INK, linespacing=1.3)


def flow(ax, module, held):
    ax.set_xlim(0, 20)
    ax.set_ylim(-0.9, 2.6)
    ax.axis("off")
    if module == "dipole":
        hook = (f"the dipole edits that change:\nGIS 6 {sc.M:g} m seaward\nGIS 5 {sc.M:g} m landward\n"
                "same every year,\nwaves not read")
        brie = "alongshore transport:\nsmooths the shoreline\nbetween domains, unaware\nof the groin; leaks sand\nat the step"
    else:
        hook = (f"the corrected groin holds back\n{sc.BLOCKING_FRACTION:.0%} of the drift arriving\n"
                f"(year 1: {held:.0f} m GIS 5 → GIS 6),\nand cancels BRIE's leak\nnear the groin in advance")
        brie = ("alongshore transport:\nsmooths the shoreline\nbetween domains, unaware\nof the groin; its leak is\n"
                "already cancelled")
    steps = [
        ("1  Barrier3D", "every domain at once:\nstorms, overwash, dunes,\nshoreline and shoreface,\n"
                         "background erosion (BE)\n→ each domain's change", STEP_C["b3d"], INK, 1.0),
        ("2  Groin hook", hook, STEP_C["groin"], C["ADDED"], 2.0),
        ("3  BRIE", brie, STEP_C["brie"], INK, 1.0),
        ("4  Management", "road manager\nsandbag manager\nbeach & dune manager\n(nourishment,\ndune rebuilds)",
         STEP_C["mgmt"], INK_MUTED, 1.0),
    ]
    w, gap, x = 4.1, 0.85, 0.35
    for k, (title, body, fc, ec, lw) in enumerate(steps):
        box(ax, x, 0.0, w, 2.45, title, body, fc, ec, lw)
        if k < len(steps) - 1:
            ax.annotate("", xy=(x + w + gap - 0.08, 1.2), xytext=(x + w + 0.08, 1.2),
                        arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3))
        x += w + gap
    ax.annotate("", xy=(2.4, -0.05), xytext=(x - gap - 2.0, -0.05),
                arrowprops=dict(arrowstyle="-|>", color=INK_MUTED, lw=0.9, connectionstyle="arc3,rad=-0.18"))
    ax.text(10, -0.3, "then the next year starts again at step 1", ha="center", va="center", fontsize=6.6,
            color=INK_MUTED)
    _letter_inside(ax, 0)
    ax.set_title("One call of Cascade.update(): the order within a model year", loc="center", fontsize=8.5)


def strip(ax, offs, before, after, i, title, note, lim, changed=None, show_before=True):
    keep = (offs >= -CELLS + 1) & (offs <= CELLS)
    x = offs[keep]
    xs = np.repeat(np.r_[x - 0.5, x[-1] + 0.5], 2)[1:-1]
    yb, ya = np.repeat(before[keep], 2), np.repeat(after[keep], 2)
    ax.fill_between(xs, -lim * 2, ya, color=C["ADDED_FILL"], lw=0)
    ax.fill_between(xs, ya, lim * 2, color=C["WATER"], alpha=0.45, lw=0)
    if show_before:
        ax.plot(xs, yb, color=INK_MUTED, lw=1.0, ls=(0, (3, 2)))
    ax.plot(xs, ya, color=INK, lw=1.5)
    if changed is not None:
        for o in changed:
            ax.annotate("", xy=(o, after[offs == o][0]), xytext=(o, before[offs == o][0]),
                        arrowprops=dict(arrowstyle="-|>", color=C["ADDED"], lw=1.6, shrinkA=0, shrinkB=0))
    ax.plot([0.5, 0.5], [-lim * 0.75, lim * 0.95], color=C["GROIN"], lw=2.4, solid_capstyle="butt")
    ax.text(-CELLS + 0.7, -lim * 0.93, "← GIS 5", fontsize=6.3, color=INK_MUTED)
    ax.text(CELLS + 0.3, -lim * 0.93, "GIS 6 →", fontsize=6.3, color=INK_MUTED, ha="right")
    ax.set_xlabel(note, fontsize=6.4, linespacing=1.25, color=INK)
    ax.set_xlim(-CELLS + 0.5, CELLS + 0.5)
    ax.set_ylim(-lim, lim)
    ax.set_xticks([])
    if i == 1:
        ax.set_ylabel("Shoreline (m)\n+ seaward")
    elif i < 5:
        ax.set_yticklabels([])
    else:
        ax.tick_params(axis="y", labelsize=6.5)
    open_frame(ax)
    _letter_inside(ax, i)
    ax.set_title(title, loc="center", fontsize=8)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--module", choices=("dipole", "corrected"), default="dipole")
    ap.add_argument("--years", type=int, default=50, help="how long the last panel repeats the year")
    ap.add_argument("--theta", type=int, default=-20)
    a = ap.parse_args()
    apply_style()
    sc.PROFILE_YEARS = (1, a.years)
    t = sc.brie_tables()
    rows, prof = sc.run(t, a.theta, a.module, a.years)
    offs = np.arange(sc.NY) - sc.LO                     # updrift (north) at positive offsets; groin face at +0.5
    start = np.zeros(sc.NY)                             # year 1 starts from the straight coast
    edit, _ = year_one_edit(t, a.theta, a.module)
    edited = start + edit
    end = prof[1]
    late = prof[a.years]
    made = float(end.sum() - edit.sum())                # sand BRIE's step made in year 1 (the edit itself sums to 0)
    total = float(rows.reach_change_m.iloc[-1])
    held = float(edit[offs == 1][0])
    keep = (offs >= -CELLS + 1) & (offs <= CELLS)
    lim = 1.7 * max(np.abs(edited[keep]).max(), np.abs(end[keep]).max())
    lim_late = 1.35 * np.abs(late[keep]).max()
    if a.module == "dipole":
        hook_note = f"the dipole: GIS 6 +{sc.M:g} m,\nGIS 5 −{sc.M:g} m (dashed:\nthe straight coast)"
        late_note = f"the same year repeated:\nthe groin shape (dashed:\nthe straight start);\nsand made {total:+,.0f} m"
    else:
        hook_note = f"the groin holds back {held:.0f} m:\nGIS 6 +{held:.0f} m, GIS 5 −{held:.0f} m\n(dashed: the straight coast)"
        late_note = f"the same year repeated:\nthe groin shape (dashed:\nthe straight start);\nsand made {total:+,.0f} m"
    fig = plt.figure(figsize=figsize("double", height=6.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 5, height_ratios=(0.7, 1))
    flow(fig.add_subplot(gs[0, :]), a.module, held)
    strip(fig.add_subplot(gs[1, 0]), offs, start, start, 1, "Start: a straight coast",
          "year 1 begins\nstraight", lim, show_before=False)
    strip(fig.add_subplot(gs[1, 1]), offs, start, start, 2, "After step 1",
          "Barrier3D: no change\nin this test (at Buxton,\nstorms and BE act here)", lim)
    strip(fig.add_subplot(gs[1, 2]), offs, start, edited, 3, "After step 2", hook_note, lim, changed=(0, 1))
    strip(fig.add_subplot(gs[1, 3]), offs, edited, end, 4, "After step 3",
          "BRIE smooths it\nalongshore (dashed:\nafter step 2)", lim)
    strip(fig.add_subplot(gs[1, 4]), offs, start, late, 5, f"After {a.years} years", late_note, lim_late)
    fig.suptitle(f"One model year with {NAME[a.module]}, from a straight coast (at {a.theta:+d}°), and the shape "
                 f"it builds by repeating", fontsize=9.5)
    out = OUT[a.module]
    out.parent.mkdir(parents=True, exist_ok=True)
    save(fig, out, dpi=250, close=True)
    common = (
        "**How to read:** (a) the order within one call of `Cascade.update()`: Barrier3D in every domain (storms, "
        "overwash, dunes, background erosion), then the groin hook (`_groin_callback`), then BRIE's alongshore "
        "transport, then the management modules; then the next year. (b)-(e) the first year of the "
        f"straight-coast run, from a straight coast at {a.theta:+d} degrees (BRIE's alongshore step alone, option "
        "A waves): the shoreline at the start, after Barrier3D (unchanged: the test has none; at Buxton storms "
        "and background erosion act here), after the groin's edit (arrows), and after BRIE smooths it alongshore "
        f"(dashed: the shoreline before each step). (f) the same year repeated for {a.years} years, on its own "
        "scale (dashed: the straight start). Red: the groin, the face between GIS 5 and GIS 6. ")
    if a.module == "dipole":
        text = ("**Tests:** nothing new; this figure shows where the source/sink dipole groin (version 1) acted in a "
                "CASCADE model year, and the shape it builds, starting from a straight coast. The dipole sat on the "
                "sandbox Cascade's hook (`cascade/cascade_groin.py`, before 563a542e). " + common +
                f"**Shows:** each year the dipole cuts a {sc.M:g} m notch and adds a {sc.M:g} m bump either side "
                f"of the groin, and BRIE spreads them; in year 1 BRIE makes no sand ({made:+.1f} m), but as the step "
                f"steepens it does, and after {a.years} years the reach has gained {total:+,.0f} m of sand that was "
                "never added (`dipole_2_explained.png`). The companion of "
                "`../2-corrected/figures/corrected_groin_steps.png`.")
    else:
        text = ("**Tests:** nothing new; this figure shows where the corrected groin (version 6, the runner's groin "
                "since 2026-10-13) acts in a CASCADE model year, and the shape it builds, starting from a straight "
                "coast; the same layout as `../1-dipole/figures/dipole_3_steps.png`. " + common +
                f"**Shows:** the groin holds back {sc.BLOCKING_FRACTION:.0%} of the drift arriving: in year 1, "
                f"{held:.0f} m of shoreline moved from GIS 5 to GIS 6, more than the dipole's fixed {sc.M:g} m, "
                "because it is read from the waves. Its second part, cancelling BRIE's leak near the groin, is zero "
                "in year 1 (the coast is still straight, so there is no bend to leak at) and grows as the step "
                f"does. BRIE then spreads the notch and bump, and after {a.years} years the shape is a fillet north "
                f"of the groin and an equal cut south of it, with {total:+,.0f} m of sand made: balanced, where "
                "the dipole made 539 m on the same coast.")
    record_caption(out, text)
    print("wrote", out, "| held/moved yr1:", round(held, 2), "| made yr1:", round(made, 2), "| total:", round(total, 1))


if __name__ == "__main__":
    main()
