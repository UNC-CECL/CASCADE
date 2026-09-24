r"""
HAT_plot_overwash_fix_explained.py -- a picture of the Barrier3D route_overwash bug
==============================================================================
Hannah asked for a figure that shows the bug plainly. Five panels:

  (a) the question the model means to ask: from the cell carrying overwash sand,
      look at the nine cells LANDWARD of it, down its own column
  (b) what the code actually looks at: row and column swapped, so a strip
      ACROSS the island somewhere else
  (c) on a narrow island (fewer rows than columns) that strip lies off the
      grid: the read is outside the island's memory
  (d) every compared run's score with the bug and fixed
  (e) what fixing it changes along the island in the run it moved most,
      natural 2010-2024 (experiments/2026-09-24-overwash-fix)

The grids in (a)-(c) are schematic (a small domain, not to scale); the
indexing is the real one from barrier3d.py line 1092.

Output: output/raw_runs/experiments/2026-09-24-overwash-fix/figures/
        route_overwash_bug_explained.png
==============================================================================
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Rectangle  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_overwash_fix_check as X  # noqa: E402
import HAT_offset_scale_wave_tuning as C  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C as COL, INK, INK_MUTED, DOMAIN_AXIS_LABEL, _title, apply_style, figsize,
    open_frame, record_caption, save)

def title(ax, i, text):
    """Letter and title as one left-aligned line: the house _title centres the
    title, which collides with the letter on these narrow panels."""
    ax.set_title(f"({'abcde'[i]})  {text}", loc="left", fontsize=9.5)


OUT = X.EXP_DIR / "figures" / "route_overwash_bug_explained.png"
LAND, WATER = "#e8dcc0", "#a8c8e0"
CELL = "#1a1a1a"
GOOD, BAD = COL["REF"], COL["EARLY"]       # green: what is meant; red: what is read


def draw_grid(ax, n_rows, n_cols, water_from_row):
    for r in range(n_rows):
        for c in range(n_cols):
            ax.add_patch(Rectangle((c, r), 1, 1, facecolor=LAND if r < water_from_row else WATER,
                                   edgecolor="white", lw=0.4))
    ax.set_xlim(-0.5, n_cols + 0.5)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def mark_cell(ax, r, c, color, lw=0, fill=True, alpha=0.85, ls="-"):
    ax.add_patch(Rectangle((c, r), 1, 1, facecolor=color if fill else "none",
                           edgecolor=color, lw=lw, alpha=alpha, ls=ls))


def axis_labels(ax, n_rows, n_cols):
    ax.annotate("", xy=(n_cols * 0.62, -1.3), xytext=(n_cols * 0.38, -1.3),
                arrowprops=dict(arrowstyle="->", color=INK_MUTED, lw=0.8), annotation_clip=False)
    ax.text(n_cols / 2, -2.0, "columns: along the coast", ha="center", va="bottom",
            fontsize=7, color=INK_MUTED)
    ax.annotate("", xy=(-1.3, n_rows * 0.7), xytext=(-1.3, n_rows * 0.3),
                arrowprops=dict(arrowstyle="->", color=INK_MUTED, lw=0.8), annotation_clip=False)
    ax.text(-2.0, n_rows / 2, "rows: across the island\n(ocean → sound)", ha="right",
            va="center", fontsize=7, color=INK_MUTED, rotation=90)


D, I = 3, 12                # the cell carrying sand: row d, column i
ROWS, COLS, WATER_ROW = 20, 26, 16


def panel_meant(ax):
    draw_grid(ax, ROWS, COLS, WATER_ROW)
    for r in range(D + 1, D + 10):
        mark_cell(ax, r, I, GOOD)
    mark_cell(ax, D, I, CELL)
    ax.add_patch(FancyArrowPatch((I + 0.5, D - 2.2), (I + 0.5, D + 0.1), arrowstyle="-|>",
                                 mutation_scale=9, color=CELL, lw=1.0))
    ax.text(I + 0.5, D - 2.4, "overwash sand", ha="center", va="bottom", fontsize=7)
    ax.text(I + 1.6, D + 5.5, "the 9 cells\nlandward\n(same column)", color=GOOD, fontsize=7,
            va="center", fontweight="bold")
    axis_labels(ax, ROWS, COLS)
    ax.set_ylim(ROWS + 0.5, -3.5)
    title(ax, 0, "What the model means to check")


def panel_read(ax):
    draw_grid(ax, ROWS, COLS, WATER_ROW)
    for r in range(D + 1, D + 10):
        mark_cell(ax, r, I, GOOD, lw=1.0, fill=False, alpha=1, ls="--")
    for c in range(D + 1, D + 10):
        mark_cell(ax, I, c, BAD)
    mark_cell(ax, D, I, CELL)
    ax.text(D + 5.5, I + 1.6, "what it actually reads:\nrow 12, columns 4–12",
            color=BAD, fontsize=7, ha="center", va="top", fontweight="bold")
    ax.text(I + 1.6, D + 3.0, "should be here", color=GOOD, fontsize=7, va="center")
    ax.set_ylim(ROWS + 0.5, -3.5)
    title(ax, 1, "What the code actually checks")


def panel_offgrid(ax):
    """A narrow island: fewer rows (8) than the column number of the cell (12),
    so the swapped read, row 12, lies below the end of the grid."""
    rows, water = 8, 6
    draw_grid(ax, rows, COLS, water)
    for r in range(D + 1, rows):
        mark_cell(ax, r, I, GOOD, lw=1.0, fill=False, alpha=1, ls="--")
    ax.add_patch(Rectangle((D + 1, I), 9, 1, facecolor=BAD, alpha=0.18, lw=0))
    for c in range(D + 1, D + 10):
        ax.add_patch(Rectangle((c, I), 1, 1, facecolor="none", edgecolor=BAD, lw=1.0, ls="--"))
    mark_cell(ax, D, I, CELL)
    ax.plot([-0.5, COLS + 0.5], [rows, rows], color=INK, lw=1.2)
    ax.text(COLS + 0.4, rows - 0.1, "end of this\nisland's grid\n(8 rows)", ha="left",
            va="bottom", fontsize=7, color=INK)
    ax.text(D + 5.5, I + 1.5, "reads row 12, which does not exist:\nwhatever is in memory "
            "here → junk values,\nand sometimes Windows kills the run",
            color=BAD, fontsize=7, ha="center", va="top", fontweight="bold")
    ax.set_ylim(ROWS + 0.5, -3.5)
    title(ax, 2, "On a narrow island: off the map")


RUN_LABELS = {
    "natural_baseline_1996": "natural 1996–2010",
    "managed_baseline_1996": "managed 1996–2010",
    "natural_ridge_1996": "natural 1996–2010, ridge",
    "natural_baseline_2010": "natural 2010–2024",
    "managed_baseline_2010": "managed 2010–2024",
    "managed_highangle0.4_2010": "managed 2010–2024, high-angle 0.4",
    "div10_managed_1996": "÷10 managed 1996–2010",
    "div10_managed_2010": "÷10 managed 2010–2024",
}


def panel_scores(ax):
    """Each run's interior RMSE with the bug and fixed."""
    import json
    members = list(RUN_LABELS)
    for y, m in enumerate(members):
        period = X.MEMBERS[m][0]
        pdir, udir = X.run_dir(m, "patched"), X.twin(m)
        rm = lambda d: float(json.loads(next(d.glob("*_run_metadata.json")).read_text(
            encoding="utf-8"))["index row"]["rmse_interior_m_yr"])
        if udir is not None:
            ax.plot(rm(udir), y, "o", ms=8, color=BAD, zorder=3)
        else:
            ax.text(rm(pdir) + 0.15, y, "with the bug: crashed (year 13)", color=BAD,
                    fontsize=7, va="center")
        ax.plot(rm(pdir), y, "o", ms=4.5, mfc="white", mec=GOOD, mew=1.4, zorder=4)
    ax.set_yticks(range(len(members)))
    ax.set_yticklabels([RUN_LABELS[m] for m in members], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0.8, 6.3)
    ax.set_xlabel("Interior RMSE against CoastSat (m/yr)")
    ax.grid(axis="x")
    open_frame(ax)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], ls="none", marker="o", ms=8, color=BAD, label="With the bug"),
                       Line2D([], [], ls="none", marker="o", ms=4.5, mfc="white", mec=GOOD,
                              mew=1.4, label="Fixed")],
              frameon=False, fontsize=7, loc="lower right")
    title(ax, 3, "Scores with the bug and fixed")


def panel_effect(ax):
    member = "natural_baseline_2010"
    p = C.interior(C.run_rates(X.run_dir(member, "patched")))
    u = C.interior(C.run_rates(X.twin(member)))
    ax.plot(u.index, u.values, color=BAD, lw=2.2, label="With the bug (runs made so far)")
    ax.plot(p.index, p.values, color=GOOD, lw=1.1, ls="--", label="Fixed")
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    big = (p - u).abs().idxmax()
    ax.annotate(f"largest change anywhere: GIS {big}, {float((p - u)[big]):+.1f} m/yr",
                (big, float(p[big])), xytext=(big + 6, 1.1), textcoords="data",
                fontsize=7, arrowprops=dict(arrowstyle="->", color=INK, lw=0.6))
    ax.set_ylim(top=2.2)
    ax.set_xlim(1, 90)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Modelled shoreline\nchange rate (m/yr)")
    ax.grid(axis="y")
    open_frame(ax)
    ax.legend(frameon=False, fontsize=7, loc="upper right", ncol=2)
    title(ax, 4, "Along the island: the run the fix changed most (natural 2010–2024)")


def main():
    apply_style()
    fig = plt.figure(figsize=figsize("double", height=9.0), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 0.8])
    panel_meant(fig.add_subplot(gs[0, 0]))
    panel_read(fig.add_subplot(gs[0, 1]))
    panel_offgrid(fig.add_subplot(gs[1, 0]))
    panel_scores(fig.add_subplot(gs[1, 1]))
    panel_effect(fig.add_subplot(gs[2, :]))
    save(fig, OUT, dpi=250, close=True)
    record_caption(OUT, (
        "The Barrier3D route_overwash indexing bug. The island grid of one domain: rows "
        "run across the island from the ocean (top) to the sound (bottom, blue), columns "
        "along the coast. (a) When overwash sand reaches a cell (black), the model means to "
        "ask whether any of the nine cells landward of it in the same column (green) is "
        "above sea level: if so the sand keeps moving inland, if not it is spread into the "
        "sound. (b) Line 1092 of barrier3d.py swaps row and column, so it reads a strip "
        "across the island elsewhere (red): row = the cell's column number, columns = its "
        "row number plus 1 to 9. (c) When the island has fewer rows than the column "
        "number, that strip lies beyond the grid: the compiled code reads whatever is in "
        "memory there and occasionally crashes the run. Grids are schematic, not to scale. "
        "(d) Interior RMSE against CoastSat for the eight runs compared, with the bug (red) "
        "and fixed (green, open): they coincide; the one run that crashed with the bug "
        "completes when fixed. (e) The natural 2010-2024 baseline along the island with the "
        "bug (red) and fixed (green dashed), the run the fix changed most: the lines "
        "coincide except at a few domains. experiments/2026-09-24-overwash-fix/NOTE.md has "
        "the numbers."))
    print(OUT)


if __name__ == "__main__":
    main()
