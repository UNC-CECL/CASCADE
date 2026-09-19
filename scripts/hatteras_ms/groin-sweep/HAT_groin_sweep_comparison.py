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
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents
                        if (_p / "pyproject.toml").exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in "
        f"scripts/hatteras_ms/groin-sweep/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hat_figure_style import (apply_style, C, C_1984, C_1997,  # noqa: E402
                              INK, INK_MUTED, error_cmap, figsize,
                              open_frame, save, _title)
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
# Colour is the PERIOD, marker is the preset (2026-09-11). Four invented hues
# were in use here -- a blue, an orange, a green and a dark red -- which spent
# two colours on a distinction the marker already carries, and put the vintage
# red on one arbitrary sweep. The house pair means the same thing on every
# figure in the project: red is the earlier period, blue the later.
SWEEP_COLORS = {
    (1984, "zeroBE"): C_1984,
    (1984, "edgeBE"): C_1984,
    (2004, "zeroBE"): C_1997,
    (2004, "edgeBE"): C_1997,
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
    axis.text(0.5, 0.5, f"{period} to {END_YEAR[period]}, {preset}\n"
              "not swept yet",
              ha="center", va="center", fontsize=8, color=INK_MUTED,
              transform=axis.transAxes)
    axis.set_xticks([])
    axis.set_yticks([])
    for spine in axis.spines.values():
        spine.set_edgecolor("0.85")


# =============================================================================
# FIGURE 1 -- the four surfaces
# =============================================================================

def fig_surfaces(surfaces):
    """The four M-f fillet-error surfaces in one 2x2 block."""
    plt = _matplotlib()

    apply_style()
    figure, axes = plt.subplots(len(PERIODS), len(PRESETS),
                                figsize=figsize("double", aspect=0.66),
                                constrained_layout=True)
    notes = []
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
                                   shading="nearest", cmap=error_cmap())
            cb = figure.colorbar(mesh, ax=axis)
            cb.set_label("|modelled − observed| fillet (m)")
            cb.outline.set_linewidth(0.6)

            best, tied = tied_best(groin)
            if len(tied) > 1:
                axis.plot(tied["M"], tied["fraction"], marker="o",
                          markersize=5.5, color="none",
                          markeredgecolor=C["ACCENT"], markeredgewidth=1.2,
                          linestyle="none", zorder=6,
                          label=f"{len(tied)} tied, M free")
            else:
                axis.plot(best["M"], best["fraction"], marker="*",
                          markersize=13, color=C["ACCENT"],
                          markeredgecolor="white", markeredgewidth=0.8,
                          linestyle="none", zorder=6,
                          label=_cell_label(best))

            _title(axis, row * len(PRESETS) + col,
                   f"{period} to {END_YEAR[period]}, {preset}")
            axis.set_xlabel("M (m/yr)")
            axis.set_ylabel("deterioration floor f")
            axis.legend(loc="upper right", fontsize=7)
            notes.append(
                "({}) observed fillet {:+.1f} m, best error {:.2f} m, reach "
                "RMSE {:.2f} m/yr.".format(
                    chr(ord("a") + row * len(PRESETS) + col),
                    OBSERVED_FILLET_M[period], best[RANK_METRIC],
                    best[REACH_METRIC]))

    _footnote(
        figure,
        "The fillet-error surface over the (M, f) grid for every swept period "
        "and source/sink preset, dark worse, with the best cell marked. "
        + " ".join(notes) +
        " The grey scales are PER PANEL, not shared: period 2's observed "
        "fillet is negative, a relaxing fillet, which no trapping at or above "
        "zero can build, so every period-2 cell carries at least that much "
        "error and one shared scale would flatten period 1 into a single tone. "
        "Compare the numbers quoted here, not the tones.", width=170)

    return save(figure, OUTPUT_DIR / "comparison_surfaces.png", close=True)[0]


# =============================================================================
# FIGURE 2 -- every optimum in one plane, plus the numbers
# =============================================================================

def fig_optima(surfaces):
    """Each sweep's optimum and valley floor in a single (M, f) plane."""
    plt = _matplotlib()

    apply_style()
    figure, (axis, table_axis) = plt.subplots(
        2, 1, figsize=figsize("double", aspect=0.85),
        gridspec_kw=dict(height_ratios=[3, 2]), constrained_layout=True)

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
        axis.plot(floor_M, grid.index, color=colour, alpha=0.35, linewidth=1.0,
                  zorder=2)

        if tie:
            axis.plot(tied["M"], tied["fraction"], marker=marker,
                      markersize=5.5, color="none", markeredgecolor=colour,
                      markeredgewidth=1.2, linestyle="none", zorder=5,
                      label=f"{label}, {len(tied)} tied, M free")
        else:
            axis.plot(best["M"], best["fraction"], marker=marker,
                      markersize=8, color=colour, markeredgecolor="white",
                      markeredgewidth=0.8, linestyle="none", zorder=6,
                      label=f"{label}, {_cell_label(best)}")

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
    _title(axis, 0, "where each sweep puts the groin")
    axis.grid()
    axis.set_axisbelow(True)
    open_frame(axis)
    axis.legend(loc="best", fontsize=7)

    table_axis.axis("off")
    table = table_axis.table(
        cellText=rows,
        colLabels=["period", "preset", "M", "f", "fillet", "observed",
                   "err (m)", "reach RMSE", "be1", "tied"],
        loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(7)
    table.scale(1, 1.35)
    for col in range(10):
        table[0, col].set_facecolor("0.93")
        table[0, col].set_text_props(weight="bold")
    _title(table_axis, 1, "every optimum, as numbers")

    _footnote(
        figure,
        "(a) Each sweep's optimum in one (M, f) plane, coloured by period — "
        "red for 1984 to 2004, blue for 2004 to 2024, as everywhere in this "
        "project — and marked by preset. The faint line behind each is that "
        "sweep's valley floor, the best M at every f, which is what makes two "
        "sweeps comparable even when their optima coincide: a flat floor means "
        "that axis is unconstrained, a steep one means it bites. (b) The same "
        "four sweeps as numbers. "
        "M reads 'tied' where a whole row of the grid scores alike -- in both "
        "period-2 sweeps the fillet relaxed, so the score is minimised by "
        "trapping nothing and every M at f = 0 is equally good. That is a "
        "statement about identifiability, not a fitted value. be1 is the "
        "background-erosion rate profiled out of the surface; it is swept "
        "only in the 1984 edgeBE sweep.", width=175)

    return save(figure, OUTPUT_DIR / "comparison_optima.png", close=True)[0]


# =============================================================================
# FIGURE 3 -- best profile per sweep
# =============================================================================

def fig_profiles(surfaces):
    """Each sweep's best LRR curve against its own CoastSat target."""
    plt = _matplotlib()

    apply_style()
    figure, axes = plt.subplots(len(PERIODS), len(PRESETS),
                                figsize=figsize("double", aspect=0.66),
                                sharex=True, constrained_layout=True)
    notes = []
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

            axis.plot(gis, observed, marker="o", markersize=3.0,
                      color=OBSERVED_COLOR, linewidth=1.4,
                      label="observed, CoastSat", zorder=4)
            axis.plot(gis, model, marker="s", markersize=3.0,
                      color=SWEEP_COLORS[(period, preset)], linewidth=1.4,
                      label=f"modelled, {_cell_label(best)}", zorder=3)
            _profile_axis(axis, period)

            bias = float(best.get("bias_window", float("nan")))
            _title(axis, row * len(PRESETS) + col,
                   f"{period} to {END_YEAR[period]}, {preset}")
            axis.legend(loc="best", fontsize=7)
            notes.append(
                "({}) reach RMSE {:.2f} m/yr, bias {:+.2f}.".format(
                    chr(ord("a") + row * len(PRESETS) + col),
                    best[REACH_METRIC], bias))

    _footnote(
        figure,
        "Each sweep's best-scoring cell against the observed shoreline change "
        "rate, one panel per period and preset, with the modelled curve in its "
        "period's colour. " + " ".join(notes) + " "
        "Each panel is scored against its OWN period's CoastSat target, so "
        "compare a curve to the black line beside it rather than across "
        "panels. A large offset with the right shape is a background-erosion "
        "problem (see bias); a right level with the wrong shape at D5/D6 is "
        "the groin.", width=170)

    return save(figure, OUTPUT_DIR / "comparison_profiles.png", close=True)[0]


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
