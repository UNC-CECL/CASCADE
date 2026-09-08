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

Reads output/comparisons/relocation_standard_setback/GIS11_profiles.npz, the
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
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
PROFILES = (PROJECT_BASE_DIR / "output" / "comparisons"
            / "relocation_standard_setback" / "GIS11_profiles.npz")

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
    axis.axhline(0.0, color=COLOR_WATER, lw=1.2, ls="--", zorder=3)
    axis.annotate("0 m MHW  (the drowning threshold)", xy=(2, 0.06),
                  fontsize=8, color=COLOR_WATER, va="bottom", zorder=6)

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
        # each label sits on the thing it names and cannot collide.
        axis.annotate(label, xy=(start + ROAD_CELLS * CELL_M / 2, 2.42),
                      ha="center", va="top", fontsize=8.5, color=colour,
                      fontweight="bold", rotation=90, zorder=7)
        # The row the drowning test actually reads.
        axis.plot([border_index * CELL_M], [row.mean() if row is not None
                                            else 0.0],
                  marker="v", ms=9, color=colour, mec="white", mew=1.0,
                  zorder=8, clip_on=False)
        axis.annotate(f"{fraction * 100:.0f}% wet",
                      xy=(border_index * CELL_M, (row.mean() if row is not None
                                                  else 0.0) - 0.16),
                      ha="center", va="top", fontsize=8.5, color=colour,
                      fontweight="bold", zorder=8)

    axis.set_xlim(0, 200)
    axis.set_ylim(-1.4, 2.5)
    axis.set_xlabel("Distance landward of the dune line (m)")
    axis.set_title(year_label, fontsize=10.5, color=INK, fontweight="bold",
                   loc="left")
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.spines["left"].set_color("#C6CCD2")
    axis.spines["bottom"].set_color("#C6CCD2")
    axis.grid(axis="y", color="#EDF0F3", lw=0.7, zorder=0)
    axis.tick_params(labelsize=8.5, colors=MUTED)


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
    axis.axhline(WET_LIMIT * 100, color=INK, lw=1.3, ls="--", zorder=4)
    axis.annotate("above this line the road drowns",
                  xy=(28, WET_LIMIT * 100 - 2.5), fontsize=8.5, color=INK,
                  ha="left", va="top", zorder=8)

    # Three callouts within 20 m of each other on a 130 m axis: stack them at
    # different heights and push the text off the marker, or they merge.
    for setback_m, label, colour, height, align, dx in (
            (MEASURED_SETBACK_M, "77 m = cell 7   measured", COLOR_OK,
             46, "right", -5),
            (87.0, "87 m = cell 8   20 m standard", COLOR_OK,
             62, "right", -5),
            (STANDARD_SETBACK_M, "97 m = cell 9   30 m standard",
             COLOR_DROWN, 78, "right", -5)):
        # SNAP TO THE CELL, not to the metre. bulldoze indexes the road at
        # int(setback / 10), so 77 m and 70 m are the SAME road position and
        # the same bar. Pointing the callout at its raw metre value drops it
        # between two bars and invites the reader to interpolate a criterion
        # that only ever takes whole-cell values.
        snapped = int(setback_m / CELL_M) * CELL_M
        axis.annotate(label, xy=(snapped + dx, height), ha=align,
                      va="center", fontsize=8.5, color=colour,
                      fontweight="bold", zorder=7)
        axis.plot([snapped], [height], marker="o", ms=6, color=colour,
                  mec="white", mew=1.0, zorder=7)
        axis.plot([snapped, snapped], [0, height], color=colour, lw=0.9,
                  ls=":", alpha=0.75, zorder=6)

    axis.set_xlim(25, 155)
    axis.set_ylim(0, 100)
    axis.set_xlabel("Road setback (m behind the dune line)")
    axis.set_ylabel("Bordering row at or below 0 m MHW (%)")
    axis.set_title("The drowning criterion is a cliff, not a slope",
                   fontsize=10.5, color=INK, fontweight="bold", loc="left")
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

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10, "figure.facecolor": "white",
        "savefig.facecolor": "white", "legend.frameon": False,
    })
    figure = plt.figure(figsize=(14.5, 5.6))
    grid_spec = figure.add_gridspec(1, 3, wspace=0.24,
                                    left=0.055, right=0.985,
                                    top=0.735, bottom=0.115)

    marks = ((MEASURED_SETBACK_M, "77 m\nmeasured", COLOR_OK),
             (STANDARD_SETBACK_M, "97 m\n30 m standard", COLOR_DROWN))
    grid_event = data[f"standard_domain_t{EVENT_YEAR_INDEX}"]
    grid_drown = data[f"standard_domain_t{DROWN_YEAR_INDEX}"]

    ax_event = figure.add_subplot(grid_spec[0, 0])
    draw_profile(ax_event, grid_event, marks,
                 "1999  -  the year the historical relocation fires")
    ax_event.set_ylabel("Elevation (m MHW)")

    ax_drown = figure.add_subplot(grid_spec[0, 1], sharey=ax_event)
    draw_profile(ax_drown, grid_drown, marks,
                 "2001  -  the year the road is given up")

    ax_cliff = figure.add_subplot(grid_spec[0, 2])
    draw_cliff(ax_cliff, grid_drown)

    figure.text(0.055, 0.945,
                "GIS 11: a standard relocation setback moves where a "
                "PRESCRIBED relocation lands",
                fontsize=14, fontweight="bold", color=INK, ha="left")
    figure.text(0.055, 0.895,
                "The historical 1999 event is stored as a displacement, so it "
                "is added to the model's current setback, not to the road's "
                "surveyed position.",
                fontsize=9.5, color=MUTED, ha="left")
    figure.text(0.055, 0.862,
                "Raising the emergent target 10 m → 30 m left the setback "
                "at 20 m rather than 0 m in 1999, so 0 + 77 = 77 m became "
                "20 + 77 = 97 m — two cells further back, across the "
                "drowning threshold.",
                fontsize=9.5, color=MUTED, ha="left")

    handles = [
        Line2D([], [], color=INK, lw=1.6, label="alongshore mean elevation"),
        Line2D([], [], color=COLOR_LAND, lw=7, alpha=0.6,
               label="alongshore 10th-90th percentile"),
        Line2D([], [], color=COLOR_WATER, lw=1.2, ls="--", label="0 m MHW"),
        Line2D([], [], marker="v", ms=8, lw=0, color=INK, mec="white",
               label="row the drowning test reads"),
    ]
    figure.legend(handles=handles, loc="upper right", ncol=4, fontsize=8.5,
                  bbox_to_anchor=(0.985, 0.845), handlelength=1.7)

    out = Path(args.out) if args.out else (
        PROFILES.parent / "HAT_GIS11_relocation_drown.png")
    figure.savefig(out, dpi=170)
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
