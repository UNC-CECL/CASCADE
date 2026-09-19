#!/usr/bin/env python3
"""Why GIS 11 drowns at a 30 m standard relocation setback, and not at 77 m.

THE QUESTION THIS ANSWERED, AND WHAT WAS DECIDED
    Introducing a flat 30 m relocation target (`relocation_setback_m`) drowned
    NC-12 at GIS 11 in all eight 1984-2004 reloc-arm runs. Nothing drowned
    under the previous per-domain measured targets. The figure exists to settle
    whether 97 m is a physically implausible place to put a road, or whether it
    is merely one cell past a threshold.

    It is the latter, and that is the point of the third panel: the criterion
    takes whole-cell values only, and cells 7 and 8 are 0% wet against cell 9
    at 24%. There is no gradual approach to failure to read a margin off.

    DECIDED 2026-09-01: the standard is 20 m, and the whole matrix runs at it.
    A 20 m target lands the 1999 event at 87 m, which is cell 8, and all eight
    drownings go away -- confirmed by re-running the twelve reloc arms, where
    relocation counts moved only 26 -> 28 and no new drowning appeared
    anywhere. But 87 m clears the threshold by ONE CELL, so 20 m is not robust
    to different forcing; it is the dry side of a cliff, not a margin.

    THE UNDERLYING COUPLING WAS NOT FIXED. `_apply_relocation` still adds a
    surveyed displacement to a modelled position, so any future change to the
    emergent rule will move where the historical events land. Anchoring it to
    an absolute setback (initial measured setback + cumulative displacement)
    would decouple the two and let the standard be chosen on its own merits.
    That remains open.

THE CHAIN, WHICH IS NOT WHAT IT LOOKS LIKE
    The standard did NOT push the road into the bay directly. It raised GIS
    11's emergent relocation from 10 m to 30 m in 1993, so the setback stood at
    20 m rather than 0 m when the historical 1999 event fired. That event is
    stored as a DISPLACEMENT, and `_apply_relocation` adds it to whatever the
    model's current setback is:

        measured target:   0 m + 77 m  =  77 m
        30 m standard:    20 m + 77 m  =  97 m

    So a change to the emergent rule moved where a PRESCRIBED historical
    relocation lands. The measured displacements were surveyed against the real
    road; adding them to a modelled road inherits the model's drift.

WHAT DROWNS IT
    `bulldoze` drowns a road when more than 20% of the cells in the row
    BORDERING the road are at or below 0 m MHW. At GIS 11 that criterion has a
    cliff between 80 m and 90 m of setback -- 0% wet against 24% wet, one cell
    apart. 77 m clears it; 97 m does not.

    The road occupies `int(setback / 10)` and the cell behind it, and the
    bordering row is the one after that, so the axis is drawn in the same whole
    cells the model indexes in rather than in smooth metres.

Usage:
    python HAT_gis11_relocation_drown_figure.py [--out PATH]

Reads output/comparisons/relocation/standard_setback/GIS11_profiles.npz, the
per-domain extract taken before the superseded runs were deleted.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# HOUSE STYLE: one typeface and one palette across every figure in this
# project. See scripts/site_layer/hat_figure_style.py and figure_making/STYLE.md. The root is
# found by searching upward (ORGANIZATION.md rule 5). This file drew in
# matplotlib's defaults until 2026-09-17 -- it never called apply_style().
import sys as _sys
from pathlib import Path as _HP
_sys.path.insert(0, str(next(_q for _q in _HP(__file__).resolve().parents
                             if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer.hat_figure_style import (apply_style, figsize,  # noqa: E402
                              FIG_W_DOUBLE, record_caption, save)
apply_style()
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

_HERE = Path(__file__).resolve()
# Anchored by SEARCHING UPWARD for the project root rather than by
# counting parent directories (2026-09-13). A counted depth is correct
# only while the file stays where it was written, and these moved into
# subfolders of hatteras_ms. Six files here already did it this way.
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / 'pyproject.toml').exists())
from site_layer import hat_figure_style as _hs  # noqa: E402
PROFILES = _hs.COMPARISONS_ROOT / "relocation" / "standard_setback" / "GIS11_profiles.npz"
# The manuscript copy, with the other figures by subject (2026-09-18): the
# same panels, the headline and its two lines moved to supporting/CAPTIONS.md.
# Written only when --out is not given.
PUBLISHED = _hs.figure_dir("management") / "gis11_relocation_drown.png"
# ON since the layout was redrawn for the house-style column (2026-09-18).
# Before that the labels, title and legend overlapped, and it was held off so
# a broken figure could not land in output/figures/.
PUBLISH = True

CELL_M = 10.0
ROAD_CELLS = 2          # 20 m road, as road_width_TS carries all run
DROWN_THRESHOLD_M = 0.0    # bulldoze drown_threshold, 0 m MSL
WET_LIMIT = 0.20           # percent_water_cells_touching_road

MEASURED_SETBACK_M = 77.0   # where the 1999 event landed under measured targets
STANDARD_SETBACK_M = 97.0   # where it landed under the 30 m standard
EVENT_YEAR_INDEX = 15       # 1999
DROWN_YEAR_INDEX = 17       # 2001, the year the road is given up

COLOR_OK = "#1B7F4B"
COLOR_DROWN = "#B71C1C"
COLOR_LAND = "#C9B27C"
COLOR_WATER = "#3E6E9E"
INK = "#15202C"
MUTED = "#5C6874"


def wet_fraction(grid, setback_m):
    """Fraction of the bordering row at or below the drown threshold.

    Mirrors `bulldoze`: the road starts at `int(setback / 10)`, runs
    `ROAD_CELLS` cells, and the row checked is the one after it.

    Args:
        grid: (cross_shore, alongshore) interior elevations in m MHW.
        setback_m: Road setback in metres.

    Returns:
        (fraction_wet, bordering_row_index, bordering_row) or (nan, idx, None)
        when the bordering row is past the end of the interior.
    """
    start = int(setback_m / CELL_M)
    border_index = start + ROAD_CELLS + 1
    if border_index >= grid.shape[0]:
        return float("nan"), border_index, None
    row = grid[border_index, :]
    return float((row <= DROWN_THRESHOLD_M).mean()), border_index, row


def draw_profile(axis, grid, setbacks, year_label):
    """One cross-shore profile with the candidate road positions on it.

    Args:
        axis: Axes to draw on.
        grid: (cross_shore, alongshore) interior elevations in m MHW.
        setbacks: (setback_m, label, colour) tuples to mark.
        year_label: Calendar year for the panel title.
    """
    rows = np.arange(grid.shape[0]) * CELL_M
    mean = grid.mean(axis=1)
    low = np.percentile(grid, 10, axis=1)
    high = np.percentile(grid, 90, axis=1)

    keep = rows <= 200
    axis.fill_between(rows[keep], low[keep], high[keep], color=COLOR_LAND,
                      alpha=0.45, lw=0, zorder=2,
                      label="alongshore 10th-90th percentile")
    axis.plot(rows[keep], mean[keep], color=INK, lw=1.6, zorder=4,
              label="alongshore mean elevation")
    axis.axhline(0.0, color=COLOR_WATER, lw=1.0, ls="--", zorder=3)

    span = axis.get_ylim()
    for setback_m, label, colour in setbacks:
        start = int(setback_m / CELL_M) * CELL_M
        fraction, border_index, row = wet_fraction(grid, setback_m)
        # The road itself, as the two cells bulldoze flattens.
        axis.add_patch(Rectangle(
            (start, -1.4), ROAD_CELLS * CELL_M, 3.9, facecolor=colour,
            alpha=0.16, edgecolor=colour, lw=1.4, zorder=5))
        # The two candidate positions are 20 m apart on a 200 m axis, so a
        # label above each band overlaps its neighbour. Rotated inside the band
        # each label sits on the thing it names and cannot collide; at the
        # 190 mm column only the distance fits, and the legend says which
        # target each colour is (2026-09-18).
        axis.annotate(label, xy=(start + ROAD_CELLS * CELL_M / 2, 2.42),
                      ha="center", va="top", fontsize=7.5, color=colour,
                      fontweight="bold", rotation=90, zorder=7)
        # The row the drowning test actually reads.
        axis.plot([border_index * CELL_M], [row.mean() if row is not None
                                            else 0.0],
                  marker="v", ms=6, color=colour, mec="white", mew=0.8,
                  zorder=8, clip_on=False)
        axis.annotate(f"{fraction * 100:.0f}% wet",
                      xy=(border_index * CELL_M + 5,
                          row.mean() if row is not None else 0.0),
                      ha="left", va="center", fontsize=7, color=colour,
                      fontweight="bold", zorder=8,
                      bbox=dict(facecolor="white", alpha=0.8, edgecolor="none",
                                pad=0.5))

    axis.set_xlim(0, 200)
    axis.set_ylim(-1.4, 2.5)
    axis.set_xlabel("Distance landward of the dune line (m)")
    axis.set_title(year_label, loc="left", fontweight="bold")
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.grid(axis="y", color="#EDF0F3", lw=0.6, zorder=0)


def draw_cliff(axis, grid):
    """Wet fraction of the bordering row against setback, with the 20% limit.

    This is the panel that answers the question: the criterion is not a slope,
    it is a step between two adjacent cells, so 77 m and 97 m sit either side
    of it and 87 m -- where a 20 m standard would land -- sits on the safe side
    by one cell.

    Args:
        axis: Axes to draw on.
        grid: (cross_shore, alongshore) interior elevations in m MHW.
    """
    setbacks = np.arange(30, 151, 10, dtype=float)
    fractions = np.array([wet_fraction(grid, s)[0] for s in setbacks])
    colours = [COLOR_DROWN if f > WET_LIMIT else COLOR_OK for f in fractions]

    # A 0% bar draws nothing, which reads as "not evaluated" rather than
    # "evaluated and dry". Every setback gets a visible stub.
    axis.bar(setbacks, np.maximum(fractions * 100, 1.4), width=7.0,
             color=colours, zorder=3)
    axis.axhline(WET_LIMIT * 100, color=INK, lw=1.0, ls="--", zorder=4,
                 label="20% wet: above it the road drowns")

    # Three callouts within 20 m of each other on a 130 m axis. Stacked
    # horizontally they ran off the axis into panel (b) at the 190 mm column
    # (2026-09-18); rotated, each sits in its own bar's column above the limit
    # line, where the three dry cells leave the panel empty.
    height = 27
    for setback_m, label, colour in (
            (MEASURED_SETBACK_M, "77 m, cell 7: measured", COLOR_OK),
            (87.0, "87 m, cell 8: 20 m standard", COLOR_OK),
            (STANDARD_SETBACK_M, "97 m, cell 9: 30 m standard", COLOR_DROWN)):
        # SNAP TO THE CELL, not to the metre. bulldoze indexes the road at
        # int(setback / 10), so 77 m and 70 m are the SAME road position and
        # the same bar. Pointing the callout at its raw metre value drops it
        # between two bars and invites the reader to interpolate a criterion
        # that only ever takes whole-cell values.
        snapped = int(setback_m / CELL_M) * CELL_M
        axis.annotate(label, xy=(snapped, height + 3), ha="center",
                      va="bottom", fontsize=7, color=colour, rotation=90,
                      fontweight="bold", zorder=7)
        axis.plot([snapped], [height], marker="o", ms=4, color=colour,
                  mec="white", mew=0.8, zorder=7)
        axis.plot([snapped, snapped], [0, height], color=colour, lw=0.8,
                  ls=":", alpha=0.75, zorder=6)

    axis.set_xlim(25, 155)
    axis.set_ylim(0, 100)
    axis.set_xlabel("Road setback (m behind the dune line)")
    axis.set_ylabel("Bordering row wet (%)")
    axis.set_title("(c) The drowning test", loc="left", fontweight="bold")
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.spines["left"].set_color("#C6CCD2")
    axis.spines["bottom"].set_color("#C6CCD2")
    axis.grid(axis="y", color="#EDF0F3", lw=0.7, zorder=0)
    axis.tick_params(labelsize=8.5, colors=MUTED)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    if not PROFILES.exists():
        raise SystemExit(
            f"{PROFILES} not found. It is the extract taken from the reloc-arm "
            f"runs before the superseded archive was deleted; regenerate it by "
            f"re-running the arms at relocation_setback_m: measured.")
    data = np.load(PROFILES)

    # HOUSE STYLE (2026-09-18): apply_style() only, no local rcParams; the
    # headline and its two lines are the caption, not the canvas.
    figure = plt.figure(figsize=figsize("double", height=3.1))
    grid_spec = figure.add_gridspec(1, 3, wspace=0.34,
                                    width_ratios=(1, 1, 1.1),
                                    left=0.07, right=0.99,
                                    top=0.92, bottom=0.32)

    marks = ((MEASURED_SETBACK_M, "77 m", COLOR_OK),
             (STANDARD_SETBACK_M, "97 m", COLOR_DROWN))
    grid_event = data[f"standard_domain_t{EVENT_YEAR_INDEX}"]
    grid_drown = data[f"standard_domain_t{DROWN_YEAR_INDEX}"]

    ax_event = figure.add_subplot(grid_spec[0, 0])
    draw_profile(ax_event, grid_event, marks,
                 "(a) 1999, the relocation fires")
    ax_event.set_ylabel("Elevation (m MHW)")

    ax_drown = figure.add_subplot(grid_spec[0, 1], sharey=ax_event)
    draw_profile(ax_drown, grid_drown, marks,
                 "(b) 2001, the road is given up")

    ax_cliff = figure.add_subplot(grid_spec[0, 2])
    draw_cliff(ax_cliff, grid_drown)

    headline = ("GIS 11: a standard relocation setback moves where a "
                "PRESCRIBED relocation lands")
    lede = ("The historical 1999 event is stored as a displacement, so it "
            "is added to the model's current setback, not to the road's "
            "surveyed position.")
    chain = ("Raising the emergent target 10 m → 30 m left the setback "
             "at 20 m rather than 0 m in 1999, so 0 + 77 = 77 m became "
             "20 + 77 = 97 m — two cells further back, across the "
             "drowning threshold.")
    handles = [
        Line2D([], [], color=INK, lw=1.6, label="alongshore mean elevation"),
        Line2D([], [], color=COLOR_LAND, lw=6, alpha=0.6,
               label="alongshore 10th-90th percentile"),
        Line2D([], [], color=COLOR_WATER, lw=1.0, ls="--",
               label="0 m MHW, the drowning threshold"),
        Line2D([], [], marker="v", ms=6, lw=0, color=INK, mec="white",
               label="row the drowning test reads"),
        Rectangle((0, 0), 1, 1, facecolor=COLOR_OK, alpha=0.3,
                  edgecolor=COLOR_OK, label="road at 77 m, measured target"),
        Rectangle((0, 0), 1, 1, facecolor=COLOR_DROWN, alpha=0.3,
                  edgecolor=COLOR_DROWN, label="road at 97 m, 30 m standard"),
        Line2D([], [], color=INK, lw=1.0, ls="--",
               label="(c) 20% wet: above it the road drowns"),
    ]
    figure.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.5,
                  frameon=False, bbox_to_anchor=(0.5, 0.0), handlelength=1.8)

    out = Path(args.out) if args.out else (
        PROFILES.parent / "HAT_GIS11_relocation_drown.png")
    # 300 dpi, the house savefig resolution. It was 170, which is a screen
    # export: the geometry was already right at 190 mm, the pixels were not.
    figure.savefig(out, dpi=300)
    if PUBLISH and args.out is None:
        save(figure, PUBLISHED)
        record_caption(PUBLISHED, (
            f"{headline.replace('PRESCRIBED', 'prescribed')}. {lede} {chain} "
            f"(a) and (b) are the GIS 11 cross-shore profile in 1999, "
            f"the year the historical relocation fires, and in 2001, the "
            f"year the road is given up, with the 77 m and 97 m road "
            f"positions marked. (c) is the drowning test by setback, "
            f"in the whole 10 m cells the model indexes in: the bordering row "
            f"goes from 0% wet to 24% wet in one cell, against a 20% limit, "
            f"so the criterion is a cliff rather than a slope. The project "
            f"settled on a 20 m standard, which lands the 1999 event at 87 m, "
            f"one cell on the dry side."))
        print(f"wrote {PUBLISHED}")
    plt.close(figure)
    print(f"wrote {out}")

    for label, grid in (("1999", grid_event), ("2001", grid_drown)):
        for setback_m in (MEASURED_SETBACK_M, 87.0, STANDARD_SETBACK_M):
            fraction, index, _ = wet_fraction(grid, setback_m)
            verdict = "DROWN" if fraction > WET_LIMIT else "ok"
            print(f"  {label}  setback {setback_m:5.0f} m  "
                  f"bordering row {index:3d}  {fraction * 100:5.1f}% wet  "
                  f"{verdict}")


if __name__ == "__main__":
    main()
