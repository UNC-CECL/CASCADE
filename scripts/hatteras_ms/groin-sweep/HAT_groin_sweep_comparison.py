#!/usr/bin/env python3
"""Cross-reference figures spanning every groin sweep at once.

`HAT_groin_sweep_figures.py` draws one sweep at a time, into that sweep's own
directory. That is the right place for a diagnostic and the wrong place for a
comparison: answering "does edgeBE put the optimum where zeroBE does?" or
"do the two periods agree about f?" currently means opening four folders and
holding four colour scales in your head.

This file draws the four sweeps together.

WHAT THE COMPARISON IS FOR
    The sweep grid is the same in every cell -- same M values, same f values,
    same observed target per period -- so the four surfaces are directly
    comparable and the interesting content is where they DIFFER:

    across presets (zeroBE vs edgeBE)
        Same period, same groin, different background erosion. If the optimum
        moves, the fitted groin is absorbing background erosion rather than
        describing the structure.

    across periods (1984-2004 vs 2004-2024)
        Same structure, different window. Period 1 straddles the 1996 repair
        and the 2003 storm and is the only window that can separate M from f;
        period 2 sits entirely past the ramp and sees only the product M*f.
        Disagreement here is the scientific result, not a defect.

WHY THE PANELS DO NOT SHARE A COLOUR SCALE
    Period 1's fillet error spans roughly 0-40 m and period 2's roughly 43-95
    m, because period 2's observed fillet is NEGATIVE (-43.2 m: the fillet
    relaxed) and no M >= 0 can build a negative fillet, so every cell carries
    at least that much error. Forcing one scale would flatten period 1 --
    where the actual optimum lives -- into a single colour. Each panel is
    scaled to its own range and the numbers are given on the panel, so the
    comparison is read from the annotations rather than from the hue.

TIES ARE DRAWN, NOT RESOLVED
    In both period-2 sweeps the whole f = 0 row scores identically: a fully
    deteriorated groin traps nothing, so M has no effect and seven cells tie
    to within 5e-5 m. Marking one of them as "best" would report a fitted M
    that is really just whichever cell sorted first. Tied sets are drawn as
    open circles and labelled as unconstrained.

Usage:
    python HAT_groin_sweep_comparison.py
    python HAT_groin_sweep_comparison.py --top-n 3

Writes to output/groin_sweep/figures/:
    comparison_surfaces.png   the four M-f error surfaces side by side
    comparison_optima.png     every sweep's optimum in one (M, f) plane
    comparison_profiles.png   each sweep's best LRR curve against CoastSat

Sweeps that have not run yet are drawn as labelled placeholders rather than
skipped, so a missing panel reads as "not swept" instead of silently
shrinking the figure.

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
PROJECT_BASE_DIR = _HERE.parents[3]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from HAT_groin_sweep_config import (  # noqa: E402
    END_YEAR,
    F_VALUES,
    M_VALUES,
    OBSERVED_FILLET_M,
    PERIODS,
    PRESETS,
)

# Imported, never copied. The single-sweep figures and these comparisons must
# reduce a sweep the SAME way -- same be1 profiling, same tie rule, same
# ranking metric -- or the two figure sets would disagree about which cell won.
from HAT_groin_sweep_figures import (  # noqa: E402
    GROIN_COLOR,
    OBSERVED_COLOR,
    RANK_METRIC,
    REACH_METRIC,
    _cell_label,
    _footnote,
    _profile_axis,
    _tie_note,
    load_scored,
    observed_curve,
    profile_be1,
    rate_curve,
    tied_best,
)

OUTPUT_DIR = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# One colour per sweep, stable across all three figures so a reader who learns
# "orange is 1984 edgeBE" on one figure keeps it on the next.
SWEEP_COLORS = {
    (1984, "zeroBE"): "#1565C0",
    (1984, "edgeBE"): "#FF8C00",
    (2004, "zeroBE"): "#2E7D32",
    (2004, "edgeBE"): "#B71C1C",
}
SWEEP_MARKERS = {"zeroBE": "o", "edgeBE": "s"}


def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def collect():
    """Loads every sweep that has results.

    Returns:
        {(period, preset): surface}, where surface is one row per (M, f) with
        be1 profiled out. Sweeps with no results are absent from the dict.
    """
    found = {}
    for period in PERIODS:
        for preset in PRESETS:
            frame = load_scored(period, preset)
            if frame is None or frame.empty:
                continue
            surface = profile_be1(frame)
            if (surface["M"] > 0).any():
                found[(period, preset)] = surface
    return found


def _placeholder(axis, period, preset):
    """Marks a panel whose sweep has not produced results yet."""
    axis.text(0.5, 0.5, f"{period}-{END_YEAR[period]}  {preset}\nnot swept yet",
              ha="center", va="center", fontsize=11, color="#999999",
              transform=axis.transAxes)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_edgecolor("#CCCCCC")


# =============================================================================
# FIGURE 1 -- the four surfaces
# =============================================================================

def fig_surfaces(surfaces):
    """The four M-f fillet-error surfaces in one 2x2 block."""
    plt = _matplotlib()

    figure, axes = plt.subplots(len(PERIODS), len(PRESETS), figsize=(14, 9))
    for row, period in enumerate(PERIODS):
        for col, preset in enumerate(PRESETS):
            axis = axes[row][col]
            surface = surfaces.get((period, preset))
            if surface is None:
                _placeholder(axis, period, preset)
                continue

            groin = surface[surface["M"] > 0]
            grid = groin.pivot(index="fraction", columns="M",
                               values=RANK_METRIC)
            mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                                   shading="nearest", cmap="viridis_r")
            figure.colorbar(mesh, ax=axis, label="|model - obs| fillet (m)")

            best, tied = tied_best(groin)
            if len(tied) > 1:
                axis.plot(tied["M"], tied["fraction"], marker="o",
                          markersize=8, color="white",
                          markeredgecolor=GROIN_COLOR, markeredgewidth=1.4,
                          linestyle="none", zorder=6,
                          label=f"{len(tied)} tied (M free)")
            else:
                axis.plot(best["M"], best["fraction"], marker="*",
                          markersize=20, color=GROIN_COLOR,
                          markeredgecolor="white", markeredgewidth=1.2,
                          linestyle="none", zorder=6,
                          label=_cell_label(best))

            axis.set_title(
                f"{period}-{END_YEAR[period]}  {preset}\n"
                f"observed fillet {OBSERVED_FILLET_M[period]:+.1f} m   |   "
                f"best err {best[RANK_METRIC]:.2f} m   |   "
                f"reach RMSE {best[REACH_METRIC]:.2f} m/yr",
                fontsize=9.5)
            axis.set_xlabel("M (m/yr)")
            axis.set_ylabel("deterioration floor f")
            axis.legend(loc="upper right", fontsize=8)

    figure.suptitle("Groin sweep error surfaces -- all periods and presets",
                    fontsize=13)
    figure.tight_layout(rect=(0, 0.05, 1, 0.965))
    _footnote(
        figure,
        "Colour scales are PER PANEL, not shared: period 2's observed fillet "
        "is negative (-43.2 m, a relaxing fillet), which no M >= 0 can build, "
        "so every period-2 cell carries at least that much error and a shared "
        "scale would flatten period 1 into one colour. Compare the annotated "
        "numbers, not the hues.", width=170)

    path = OUTPUT_DIR / "comparison_surfaces.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


# =============================================================================
# FIGURE 2 -- every optimum in one plane, plus the numbers
# =============================================================================

def fig_optima(surfaces):
    """Each sweep's optimum and valley floor in a single (M, f) plane."""
    plt = _matplotlib()

    figure, (axis, table_axis) = plt.subplots(
        2, 1, figsize=(12, 10), gridspec_kw=dict(height_ratios=[3, 2]))

    rows = []
    for (period, preset), surface in sorted(surfaces.items()):
        groin = surface[surface["M"] > 0]
        colour = SWEEP_COLORS[(period, preset)]
        marker = SWEEP_MARKERS[preset]
        label = f"{period}-{END_YEAR[period]} {preset}"

        best, tied = tied_best(groin)
        tie = _tie_note(tied)

        # The valley floor shows the SHAPE of each sweep's constraint, which is
        # what makes two sweeps comparable even when their optima coincide: a
        # flat floor means the axis is unconstrained, a steep one means it bites.
        grid = groin.pivot(index="fraction", columns="M", values=RANK_METRIC)
        floor_M = [grid.columns[int(np.nanargmin(grid.loc[f].values))]
                   if np.isfinite(grid.loc[f].values).any() else np.nan
                   for f in grid.index]
        axis.plot(floor_M, grid.index, color=colour, alpha=0.35, linewidth=1.4,
                  zorder=2)

        if tie:
            axis.plot(tied["M"], tied["fraction"], marker=marker, markersize=9,
                      color="white", markeredgecolor=colour, markeredgewidth=2,
                      linestyle="none", zorder=5,
                      label=f"{label}  ({len(tied)} tied, M free)")
        else:
            axis.plot(best["M"], best["fraction"], marker=marker,
                      markersize=15, color=colour, markeredgecolor="white",
                      markeredgewidth=1.4, linestyle="none", zorder=6,
                      label=f"{label}  {_cell_label(best)}")

        rows.append([
            f"{period}-{END_YEAR[period]}", preset,
            "tied" if tie else f"{best['M']:g}",
            f"{best['fraction']:.2f}",
            f"{best['fillet_m']:+.1f}",
            f"{OBSERVED_FILLET_M[period]:+.1f}",
            f"{best[RANK_METRIC]:.2f}",
            f"{best[REACH_METRIC]:.2f}",
            "-" if pd.isna(best.get("be1")) else f"{best['be1']:g}",
            f"{len(tied)}" if tie else "1",
        ])

    axis.set_xlim(min(m for m in M_VALUES if m > 0) - 8, max(M_VALUES) + 8)
    axis.set_ylim(min(F_VALUES) - 0.12, max(F_VALUES) + 0.12)
    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.set_title("Where each sweep puts the groin\n"
                   "faint lines are valley floors: flat = that axis is "
                   "unconstrained", fontsize=11)
    axis.grid(alpha=0.25)
    axis.legend(loc="best", fontsize=8.5)

    table_axis.axis("off")
    table = table_axis.table(
        cellText=rows,
        colLabels=["period", "preset", "M", "f", "fillet", "observed",
                   "err (m)", "reach RMSE", "be1", "tied"],
        loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    for col in range(10):
        table[0, col].set_facecolor("#EEEEEE")
        table[0, col].set_text_props(weight="bold")

    figure.tight_layout(rect=(0, 0.05, 1, 1))
    _footnote(
        figure,
        "M reads 'tied' where a whole row of the grid scores alike -- in both "
        "period-2 sweeps the fillet relaxed, so the score is minimised by "
        "trapping nothing and every M at f = 0 is equally good. That is a "
        "statement about identifiability, not a fitted value. be1 is the "
        "background-erosion rate profiled out of the surface; it is swept "
        "only in the 1984 edgeBE sweep.", width=175)

    path = OUTPUT_DIR / "comparison_optima.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


# =============================================================================
# FIGURE 3 -- best profile per sweep
# =============================================================================

def fig_profiles(surfaces):
    """Each sweep's best LRR curve against its own CoastSat target."""
    plt = _matplotlib()

    figure, axes = plt.subplots(len(PERIODS), len(PRESETS), figsize=(14, 9),
                               sharex=True)
    for row, period in enumerate(PERIODS):
        for col, preset in enumerate(PRESETS):
            axis = axes[row][col]
            surface = surfaces.get((period, preset))
            if surface is None:
                _placeholder(axis, period, preset)
                continue

            groin = surface[surface["M"] > 0]
            best, tied = tied_best(groin)
            gis, model = rate_curve(best)
            _, observed = observed_curve(period)

            axis.plot(gis, observed, marker="o", color=OBSERVED_COLOR,
                      linewidth=2.0, label="CoastSat observed", zorder=4)
            axis.plot(gis, model, marker="s",
                      color=SWEEP_COLORS[(period, preset)], linewidth=2.0,
                      label=f"model {_cell_label(best)}", zorder=3)
            _profile_axis(axis, period)

            bias = float(best.get("bias_window", float("nan")))
            axis.set_title(
                f"{period}-{END_YEAR[period]}  {preset}   "
                f"RMSE {best[REACH_METRIC]:.2f}  bias {bias:+.2f} m/yr",
                fontsize=9.5)
            axis.legend(loc="best", fontsize=8)

    figure.suptitle("Best-fit shoreline change against CoastSat -- all sweeps",
                    fontsize=13)
    figure.tight_layout(rect=(0, 0.05, 1, 0.965))
    _footnote(
        figure,
        "Each panel is scored against its OWN period's CoastSat target, so "
        "compare a curve to the black line beside it rather than across "
        "panels. A large offset with the right shape is a background-erosion "
        "problem (see bias); a right level with the wrong shape at D5/D6 is "
        "the groin.", width=170)

    path = OUTPUT_DIR / "comparison_profiles.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.parse_args()

    surfaces = collect()
    print("=" * 72)
    print("GROIN SWEEP COMPARISON")
    print("=" * 72)
    for period in PERIODS:
        for preset in PRESETS:
            key = (period, preset)
            if key not in surfaces:
                print(f"  {period}-{END_YEAR[period]} {preset:<8} not swept yet")
                continue
            groin = surfaces[key][surfaces[key]["M"] > 0]
            best, tied = tied_best(groin)
            tie = _tie_note(tied)
            print(f"  {period}-{END_YEAR[period]} {preset:<8} "
                  f"{len(groin):>3} cells   "
                  f"{'TIED (M free)' if tie else _cell_label(best):<28} "
                  f"err {best[RANK_METRIC]:>6.2f} m")

    if not surfaces:
        print("\n  no sweep has results yet; nothing to compare")
        return 1

    for path in (fig_surfaces(surfaces), fig_optima(surfaces),
                 fig_profiles(surfaces)):
        print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
