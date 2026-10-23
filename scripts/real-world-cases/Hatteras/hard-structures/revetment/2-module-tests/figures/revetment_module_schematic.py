"""
How the simple revetment works, as two schematics: four cross-sections and the yearly loop.

    python revetment_module_schematic.py

Draws only; reads no model output. Writes revetment_how_it_works.png and
revetment_yearly_loop.png beside this script, PDFs and captions in supporting/.

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
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent
REPO = next(p for p in HERE.parents if (p / "cascade").is_dir())
sys.path.insert(0, str(REPO / "scripts"))
from site_layer.hat_figure_style import (C, INK, INK_MUTED, apply_style,  # noqa: E402
                                         caption, figsize, save, _title)

# --- CONFIG ------------------------------------------------------------------
SAND = "#e9dcc0"
BERM = 1.4          # m above MHW, schematic
DUNE = 2.0          # dune height above the berm
FLOOR = 1.2         # the wall's crest floor above the berm
GAP0 = 30           # m, the starting gap
DUNE_X = 40         # dune toe at build
# -----------------------------------------------------------------------------


# Cross-shore profile: shoreface, beach, dune, interior (ocean to the left)
def profile(shore_x, dune_x, dune_h=DUNE):
    x = np.array([-60, shore_x, dune_x, dune_x + 8, dune_x + 16, dune_x + 24, 160])
    z = np.array([-3.5, 0, BERM, BERM + dune_h, BERM + dune_h, BERM + 0.3, BERM + 0.3])
    return x, z


# Sand, water and an optional ghost of an earlier profile
def draw_profile(ax, shore_x, dune_x, dune_h=DUNE, ghost=None):
    ax.axhspan(-4, 0, color=C["WATER"], alpha=0.5, lw=0, zorder=0)
    if ghost is not None:
        gx, gz = profile(*ghost)
        ax.plot(gx, gz, color=INK_MUTED, lw=0.8, ls=":", zorder=3)
    x, z = profile(shore_x, dune_x, dune_h)
    ax.fill_between(x, -4, z, color=SAND, lw=0, zorder=1)
    ax.plot(x, z, color=INK, lw=0.9, zorder=2)
    ax.set(xlim=(-30, 120), ylim=(-2.5, 5.2))
    ax.axis("off")


# The wall on the seaward dune face, solid while standing, dashed once failed
def draw_wall(ax, dune_x, failed=False):
    poly = np.array([[dune_x - 2, BERM - 0.6], [dune_x + 6, BERM + DUNE],
                     [dune_x + 9, BERM + DUNE], [dune_x + 2, BERM - 0.6]])
    if failed:
        ax.fill(*poly.T, fill=False, ec=INK_MUTED, ls="--", lw=0.9, zorder=4)
    else:
        ax.fill(*poly.T, color=C["ROAD"], zorder=4)
    ax.axhline(BERM + FLOOR, xmin=0.42, xmax=0.62, color=C["ACCENT"], lw=0.9, ls="--")


# A labelled bracket between two cross-shore positions
def bracket(ax, x0, x1, y, text):
    ax.annotate("", xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="<->", color=C["ACCENT"], lw=0.9))
    ax.text((x0 + x1) / 2, y + 0.25, text, ha="center", va="bottom",
            color=C["ACCENT"], fontsize=7)


# Figure 1: four cross-sections
def how_it_works():
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", aspect=0.55))
    a, b, c, d = axes.flat

    # (a) no wall: shoreline and dune retreat together
    draw_profile(a, 20, DUNE_X + 20, ghost=(0, DUNE_X))
    a.annotate("", xy=(DUNE_X + 30, BERM + DUNE + 0.9), xytext=(DUNE_X + 8, BERM + DUNE + 0.9),
               arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
    a.text(DUNE_X + 19, BERM + DUNE + 1.1, "shoreline and dune move together",
           ha="center", fontsize=7)
    a.text(-25, -2.3, "ocean", color=INK_MUTED, fontsize=7)
    _title(a, 0, "no wall: the island rolls back")

    # (b) wall built: dune line fixed, crest floored, gap = the beach
    draw_profile(b, 0, DUNE_X)
    draw_wall(b, DUNE_X)
    bracket(b, 0, DUNE_X - 2, -0.9, f"gap {GAP0} m")
    b.text(DUNE_X + 27, BERM + FLOOR + 0.05, "crest floor", color=C["ACCENT"],
           fontsize=7, va="top")
    b.text(DUNE_X + 4, BERM + DUNE + 0.5, "wall: dune line can't move",
           ha="center", fontsize=7)
    _title(b, 1, "wall built")

    # (c) years later: shoreline retreats, wall stays, gap shrinks, overwash stops
    draw_profile(c, 22, DUNE_X, ghost=(0, DUNE_X))
    draw_wall(c, DUNE_X)
    bracket(c, 22, DUNE_X - 2, -0.9, "gap 8 m")
    c.annotate("", xy=(DUNE_X - 1, BERM + 0.9), xytext=(5, BERM + 0.9),
               arrowprops=dict(arrowstyle="->", color=C["REF"], lw=1.1))
    c.text(4, BERM + 1.15, "storm", color=C["REF"], fontsize=7)
    c.text(DUNE_X + 30, BERM + DUNE + 0.5, "no overwash while\nthe crest tops the runup",
           ha="center", fontsize=7)
    _title(c, 2, "shoreline retreats, wall holds")

    # (d) gap reaches 0: the repaired wall holds the shoreline, with added sand
    draw_profile(d, DUNE_X - 2, DUNE_X, ghost=(0, DUNE_X))
    xs = np.linspace(-30, DUNE_X - 2, 200)
    held = np.interp(xs, *profile(DUNE_X - 2, DUNE_X))
    would = np.interp(xs, np.array([-60, DUNE_X + 30]), np.array([-3.5, 0]))
    d.fill_between(xs, would, held, where=held > would, color=C["ADDED_FILL"],
                   ec=C["ADDED"], lw=0.6, zorder=2)
    draw_wall(d, DUNE_X)
    d.text(5, 1.0, "sand added each year\nto hold the shoreline", color=C["ADDED"],
           fontsize=7, ha="center")
    d.text(DUNE_X + 30, BERM + DUNE + 0.5, "wall repaired:\nnever fails", ha="center", fontsize=7)
    _title(d, 3, "gap reaches 0: shoreline held at the wall")

    fig.tight_layout()
    caption(fig, (
        "How the simple revetment works, in cross-section (schematic, not to scale; "
        "ocean to the left). (a) With no wall, Barrier3D moves the dune line with the "
        "shoreline, so the island rolls back. (b) The wall fixes the dune line and "
        f"floors the dune crest (dashed purple, {FLOOR} m above the berm); the beach "
        "in front is bookkeeping, the gap from the shoreline to the wall, starting at "
        f"Barrier3D's implicit beach width ({GAP0} m in the test). (c) The shoreline "
        "keeps retreating and the gap shrinks by the same amount each year; storms "
        "that don't top the crest no longer overwash. (d) The wall is repaired, so it "
        "never fails: once the gap reaches 0, the shoreline is put back at the wall "
        "each year through CASCADE's nourishment path, and the sand that takes "
        "(orange) is tallied as a hidden nourishment. A real wall adds no sand; the "
        "beach in front lowers instead, which the model has no beach to show. "
        "Dotted: the profile at build."))
    out = save(fig, workspace_twin(HERE) / "revetment_how_it_works.png", close=True)
    print("wrote", out[0])


# One rounded box with text
def box(ax, x, y, text, w=0.27, h=0.2, color=None):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc=color or "white", ec=INK, lw=0.8))
    ax.text(x, y, text, ha="center", va="center", fontsize=7.2)


# An arrow between two points, with an optional label
def arrow(ax, p0, p1, text=None, offset=(0, 0.03)):
    ax.annotate("", xy=p1, xytext=p0, arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
    if text:
        ax.text((p0[0] + p1[0]) / 2 + offset[0], (p0[1] + p1[1]) / 2 + offset[1],
                text, ha="center", fontsize=6.8, color=INK_MUTED)


# Figure 2: what the module does each model year
def yearly_loop():
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.3))
    ax.set(xlim=(0, 1), ylim=(0.2, 0.92))
    ax.axis("off")
    box(ax, 0.17, 0.68, "before the year\n\nhold the dune line\n(dune migration off)\nfloor the dune crest",
        h=0.42, color=C["ACCENT_FILL"])
    box(ax, 0.5, 0.68, "CASCADE year\n\nstorms, dune growth,\nshoreline moves\n(BRIE + Barrier3D, unchanged)",
        h=0.42, color=C["BASE_FILL"])
    box(ax, 0.83, 0.68, "after the year\n\ngap = gap − shoreline retreat\nif gap < 0: put the shoreline\nback at the wall (nourishment\npath), tally the sand",
        h=0.42, color=C["ACCENT_FILL"])
    arrow(ax, (0.305, 0.68), (0.365, 0.68))
    arrow(ax, (0.635, 0.68), (0.695, 0.68))
    ax.plot([0.83, 0.83, 0.17], [0.47, 0.3, 0.3], color=INK, lw=0.8)
    arrow(ax, (0.17, 0.3), (0.17, 0.465))
    ax.text(0.5, 0.32, "next year (the wall is repaired, never fails)", fontsize=6.8,
            color=INK_MUTED, ha="center")
    caption(fig, (
        "What the simple revetment does each model year. It runs around "
        "cascade.update() like the nourishment schedule and edits nothing in CASCADE, "
        "BRIE or Barrier3D: before the year it switches dune migration off and floors "
        "the dune crest; after it, it moves the gap by the shoreline change and, once "
        "the shoreline has reached the wall, puts it back with the sand that takes."))
    out = save(fig, workspace_twin(HERE) / "revetment_yearly_loop.png", close=True)
    print("wrote", out[0])


# Draw both
def main() -> None:
    apply_style()
    how_it_works()
    yearly_loop()


if __name__ == "__main__":
    main()
