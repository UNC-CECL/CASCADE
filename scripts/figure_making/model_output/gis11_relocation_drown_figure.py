#!/usr/bin/env python3
"""
Why GIS 11 drowns at a 30 m standard relocation setback and not at 77 m: the wet-cell criterion is a cliff.

    python scripts/figure_making/model_output/gis11_relocation_drown_figure.py [--out PATH]

Reads the saved GIS 11 profiles (comparisons/relocation/standard_setback/);
writes the figure there and the manuscript copy to output/figures/5-results/.
Details: scripts/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# House style (site_layer/hat_figure_style.py), applied at import
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
# Project root found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / 'pyproject.toml').exists())
from site_layer import hat_figure_style as _hs  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
PROFILES = _hs.COMPARISONS_ROOT / "relocation" / "standard_setback" / "GIS11_profiles.npz"
# The manuscript copy, written only when --out is not given
PUBLISHED = _hs.figure_dir("results") / "gis11_relocation_drown.png"
# On since the layout was redrawn for the house-style column (2026-09-18)
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
# -----------------------------------------------------------------------------


# Fraction of the road's bordering row at or below the drown threshold, as bulldoze sees it
def wet_fraction(grid, setback_m):
    start = int(setback_m / CELL_M)
    border_index = start + ROAD_CELLS + 1
    if border_index >= grid.shape[0]:
        return float("nan"), border_index, None
    row = grid[border_index, :]
    return float((row <= DROWN_THRESHOLD_M).mean()), border_index, row


# One cross-shore profile with the candidate road positions on it
def draw_profile(axis, grid, setbacks, year_label):
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
        # Labels rotated inside each band: the two positions are too close for labels above
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


# Wet fraction against setback, with the 20% limit: the cliff
def draw_cliff(axis, grid):
    setbacks = np.arange(30, 151, 10, dtype=float)
    fractions = np.array([wet_fraction(grid, s)[0] for s in setbacks])
    colours = [COLOR_DROWN if f > WET_LIMIT else COLOR_OK for f in fractions]

    # Every setback gets a visible stub, so 0% reads as evaluated and dry
    axis.bar(setbacks, np.maximum(fractions * 100, 1.4), width=7.0,
             color=colours, zorder=3)
    axis.axhline(WET_LIMIT * 100, color=INK, lw=1.0, ls="--", zorder=4,
                 label="20% wet: above it the road drowns")

    # Callouts rotated, each in its own bar's column above the limit line
    height = 27
    for setback_m, label, colour in (
            (MEASURED_SETBACK_M, "77 m, cell 7: measured", COLOR_OK),
            (87.0, "87 m, cell 8: 20 m standard", COLOR_OK),
            (STANDARD_SETBACK_M, "97 m, cell 9: 30 m standard", COLOR_DROWN)):
        # Snap to the cell: bulldoze places the road at int(setback / 10)
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


# Run: load the profiles, draw the three panels, save both copies
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    if not PROFILES.exists():
        raise SystemExit(
            f"{PROFILES} not found. It is the extract taken from the reloc-arm "
            f"runs before the superseded archive was deleted; regenerate it by "
            f"re-running the arms at relocation_setback_m: measured.")
    data = np.load(PROFILES)

    # House style: apply_style() only; the headline and notes are the caption
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
    # 300 dpi, the house savefig resolution
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
