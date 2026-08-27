#!/usr/bin/env python3
"""Every scenario, every preset, both periods, on one page.

THE FIGURE
    A 2 x 3 grid: periods as rows (1984-2004 above 2004-2024), source/sink
    presets as columns (zeroBE, edgeBE, calibBE -- left to right in order of
    increasing correction). Each panel draws one line per management scenario
    against the section 8 CoastSat target.

    So the two contrasts read on different axes: scanning ACROSS a row shows
    what the source/sink term does, and the spread WITHIN a panel shows what
    management does. The target is the same heavy black line in every panel,
    which is what makes the across-row read a skill comparison rather than
    just a shape comparison.

WHAT IS ON THE Y AXIS
    Shoreline change rate, m/yr, (+) seaward -- read from each run's own
    `*_shoreline_change_rate.csv`. That file is written by the pipeline from
    the same array section 12 scores, so this figure and the reported skill
    numbers cannot disagree about what a run did.

SHARED SCALES
    x is shared down each column: GIS domain 1-90, the whole island.

    y is shared across ALL SIX PANELS by default, so a change rate has the
    same height everywhere on the page and the two periods can be compared
    directly by eye. That is the whole point of the figure: if the axes
    differed, a period-2 line that looked steeper than a period-1 line might
    only be a different scale, and every amplitude read would need a glance
    at the tick labels first.

    The cost is small here, and was measured rather than assumed. Period 1
    spans 6.40 m/yr across every run and the target, period 2 spans 8.23, and
    the two together span 8.33 -- so on a common axis period 1 still occupies
    77% of the height. There is no meaningful squashing to trade away.
    `--y-per-row` restores an independent range per period for the case where
    one period's detail has to be read closely.

COLOUR
    A single-hue sequential ramp ordered by management intensity: natural
    (lightest) through to full_management (darkest). The ordering is in the
    colour, so "more management" reads as "darker" without consulting the
    legend, and a single hue stays legible under the common colour-vision
    deficiencies. The observed target is black and heavier than any model
    line, so it never competes with a scenario for attention.

    Relocation arms are the SAME colour as their non-reloc twin, dashed. They
    sit almost exactly on top of it -- the two differ in the fifth decimal of
    mean bias -- and drawing them as a distinct colour would imply a
    separation that is not there. Dashed-over-solid shows the overlap
    honestly, and any real divergence would immediately stand out.

Usage:
    python scripts/hatteras_ms/HAT_scenario_grid.py
    python scripts/hatteras_ms/HAT_scenario_grid.py --no-reloc --out FIG.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live in scripts/hatteras_ms/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
for _path in (SCRIPTS_DIR, _HERE.parent):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline.annotations import (                    # noqa: E402
    add_geographic_annotations,
)
from cascade_pipeline.coastsat_loess import (                 # noqa: E402
    CoastSatDataset,
    LoessConfig,
    build_coastsat_series,
)
from cascade_pipeline.hindcast import build_target_table      # noqa: E402
from hatteras_site_config import (                            # noqa: E402
    HATTERAS_ANNOTATIONS,
    HATTERAS_DOMAINS,
)

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
COASTSAT_BASE = SCRIPTS_DIR / "input_prep" / "5-scr" / "CoastSat"
DEFAULT_OUT = (PROJECT_BASE_DIR / "output" / "comparisons"
               / "scenario_grid_by_preset.png")

PERIODS = ((1984, 2004), (2004, 2024))
PRESETS = ("zeroBE", "edgeBE", "calibBE")

# Section 8's settings, matching the runner, so the target drawn here is the
# curve the runs were scored against rather than a second opinion.
LOESS_CONFIG = LoessConfig(window_domains=(7, 10), skip_southern_domains=10)

# WHICH MODEL COLUMN THESE PANELS DRAW. The observed curve on every panel
# is a CoastSat LRR -- a per-transect OLS slope through the period -- so
# the model side is read from lrr_m_yr, the run's matching OLS slope
# through its annual states, rather than change_rate_m_yr, which is a net
# displacement over a span. Set to "change_rate_m_yr" to redraw a
# pre-2026-08-22 version of this figure.
RATE_COLUMN = "lrr_m_yr"
RATE_LABEL = ("LRR" if RATE_COLUMN == "lrr_m_yr"
              else "endpoint difference")

TARGET_WINDOW = 10

# Ordered by management intensity -- the order the colour ramp encodes.
SCENARIO_ORDER = ("natural", "beachdune_only", "roadway_only",
                  "full_no_fill", "full_management")
SCENARIO_LABEL = {
    "natural": "natural (no management)",
    "beachdune_only": "beach/dune only",
    "roadway_only": "roadway only",
    "full_no_fill": "full, no nourishment",
    "full_management": "full management",
}

TARGET_COLOUR = "black"
RAMP = plt.get_cmap("YlGnBu")
# Starts at 0.35, not 0.0: the pale end of any sequential map disappears
# against white, and the lightest scenario still has to be readable.
SCENARIO_COLOUR = {
    name: RAMP(0.35 + 0.62 * i / (len(SCENARIO_ORDER) - 1))
    for i, name in enumerate(SCENARIO_ORDER)
}


def classify(run_name):
    """Maps a run directory name to (scenario, relocations).

    Reads the switch tokens rather than matching whole names, because the
    token set is not the same in both periods -- period 2 carries a
    nourishment token that period 1 has no reason to.

    Returns:
        (scenario_key, reloc_bool), or (None, None) for a run this figure
        does not draw (anything with the groin attached).
    """
    tokens = run_name.split("_")
    if "groin" in tokens:            # "nogroin" is its own token, so this is
        return None, None            # the groin-attached arm only
    road = "road" in tokens
    bdm = "bdm" in tokens
    reloc = "reloc" in tokens

    if not road and not bdm:
        scenario = "natural"
    elif road and not bdm:
        scenario = "roadway_only"
    elif not road and bdm:
        scenario = "beachdune_only"
    elif "nonourish" in tokens:
        scenario = "full_no_fill"
    else:
        scenario = "full_management"
    return scenario, reloc


def load_runs(period_start, period_end, preset):
    """Every nogroin run for one period/preset, keyed by (scenario, reloc).

    Returns:
        {(scenario, reloc): Series indexed by GIS domain}, empty if the
        preset directory does not exist.
    """
    preset_dir = RAW_RUNS / f"{period_start}_{period_end}" / preset
    if not preset_dir.is_dir():
        return {}

    series = {}
    for run_dir in sorted(preset_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        scenario, reloc = classify(run_dir.name)
        if scenario is None:
            continue
        hits = sorted(run_dir.glob("*_shoreline_change_rate.csv"))
        if not hits:
            print(f"  ! no rate CSV in {run_dir.name}")
            continue
        frame = pd.read_csv(hits[0]).set_index("gis_domain")
        # lrr_m_yr where the run has it, change_rate_m_yr otherwise.
        # The target on these axes is a CoastSat LRR, so the model side
        # has to be one too; a run written before the column existed
        # still plots, and says so.
        if RATE_COLUMN in frame.columns:
            series[(scenario, reloc)] = frame[RATE_COLUMN]
        else:
            print(f"  ! {run_dir.name}: no {RATE_COLUMN}; falling back "
                  f"to change_rate_m_yr (endpoint difference)")
            series[(scenario, reloc)] = frame["change_rate_m_yr"]
    return series


def load_target(period_start):
    """The section 8 CoastSat target for one period, as a Series by domain."""
    csv_path = COASTSAT_BASE / f"{period_start}_{period_start + 20}" \
        / "transect_lrr_full.csv"
    built = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {period_start}",
                         period_start=period_start, csv_path=str(csv_path))],
        period_start, LOESS_CONFIG)
    if not built:
        raise FileNotFoundError(f"CoastSat transects failed to load: "
                                f"{csv_path}")
    table = build_target_table(built[0], LOESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    return pd.Series(np.asarray(table["target_lrr_m_yr"], dtype=float),
                     index=np.asarray(table["gis_domain"], dtype=int))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--no-reloc", action="store_true",
                        help="omit the dashed relocation arms")
    parser.add_argument("--no-annotations", action="store_true",
                        help="omit villages, piers, groin and shoal zones")
    parser.add_argument("--y-per-row", action="store_true",
                        help="give each period its own y range instead of "
                             "one shared across all panels (default shared, "
                             "so periods compare directly)")
    args = parser.parse_args()

    print("=" * 74)
    print("SCENARIO GRID: management x source/sink x period")
    print("=" * 74)

    targets = {start: load_target(start) for start, _ in PERIODS}

    figure, axes = plt.subplots(
        len(PERIODS), len(PRESETS), figsize=(19, 9.5),
        sharex=True, squeeze=False)

    drawn_scenarios = set()
    any_reloc = False
    row_ranges = []      # per-row value arrays, pooled below if y is shared

    for row, (start, end) in enumerate(PERIODS):
        row_axes = axes[row]
        row_values = []

        for column, preset in enumerate(PRESETS):
            ax = row_axes[column]
            runs = load_runs(start, end, preset)
            print(f"  {start}-{end} {preset:<8} {len(runs)} run(s)")

            target = targets[start]
            ax.plot(target.index, target.values, color=TARGET_COLOUR,
                    linewidth=2.6, zorder=5,
                    label="CoastSat target (LOESS 10-domain)")
            row_values.append(target.values)

            for scenario in SCENARIO_ORDER:
                for reloc in (False, True):
                    if reloc and args.no_reloc:
                        continue
                    rates = runs.get((scenario, reloc))
                    if rates is None:
                        continue
                    ax.plot(
                        rates.index, rates.values,
                        color=SCENARIO_COLOUR[scenario],
                        linewidth=1.7 if not reloc else 1.5,
                        linestyle="--" if reloc else "-",
                        alpha=0.95, zorder=3 + reloc,
                        label=SCENARIO_LABEL[scenario] if not reloc else None)
                    row_values.append(rates.values)
                    drawn_scenarios.add(scenario)
                    any_reloc = any_reloc or reloc

            ax.axhline(0.0, color="0.55", linewidth=0.8, zorder=1)
            if row == 0:
                ax.set_title(preset, fontsize=13, fontweight="bold", pad=8)
            if column == 0:
                ax.set_ylabel(f"{start}–{end}\n"
                              f"shoreline change rate, {RATE_LABEL} "
                              f"(m/yr)", fontsize=11)
            if row == len(PERIODS) - 1:
                ax.set_xlabel("GIS domain (1 = south / Cape Point  →  "
                              "90 = north / Oregon Inlet)", fontsize=10)
            ax.set_xlim(HATTERAS_DOMAINS.first_gis_id,
                        HATTERAS_DOMAINS.last_gis_id)
            ax.grid(True, alpha=0.25, linewidth=0.6)

        stacked = np.concatenate([np.asarray(v, dtype=float)
                                  for v in row_values])
        row_ranges.append(stacked[np.isfinite(stacked)])

        if not args.no_annotations:
            for ax in row_axes:
                add_geographic_annotations(ax, HATTERAS_ANNOTATIONS)

    # Applied after every panel is drawn, because the shared case needs the
    # pooled range of the whole figure and cannot be set row by row.
    def apply_limits(target_axes, values):
        finite = values[np.isfinite(values)]
        if not finite.size:
            return
        pad = 0.08 * (finite.max() - finite.min())
        for ax in target_axes:
            ax.set_ylim(finite.min() - pad, finite.max() + pad)

    if args.y_per_row:
        for row_axes, values in zip(axes, row_ranges):
            apply_limits(row_axes, values)
    else:
        apply_limits(axes.ravel(), np.concatenate(row_ranges))

    handles = [Line2D([], [], color=TARGET_COLOUR, linewidth=2.6,
                      label="CoastSat target (observed)")]
    handles += [Line2D([], [], color=SCENARIO_COLOUR[s], linewidth=1.9,
                       label=SCENARIO_LABEL[s])
                for s in SCENARIO_ORDER if s in drawn_scenarios]
    if any_reloc:
        handles.append(Line2D([], [], color="0.35", linewidth=1.5,
                              linestyle="--",
                              label="+ historical relocations (1989, 1999)"))
    figure.legend(handles=handles, loc="lower center", ncol=len(handles),
                  frameon=False, fontsize=10.5,
                  bbox_to_anchor=(0.5, -0.005))

    figure.suptitle(
        "Hatteras Island hindcast: management scenarios by source/sink preset",
        fontsize=15, fontweight="bold", y=0.985)
    scale_note = ("y axis shared across all panels — rates compare directly "
                  "between periods and presets"
                  if not args.y_per_row else
                  "y axis independent per period — row amplitudes are NOT "
                  "directly comparable")
    figure.text(
        0.5, 0.945,
        f"columns share the observed target; {scale_note}",
        ha="center", fontsize=9.5, style="italic", color="0.35")

    figure.tight_layout(rect=[0, 0.045, 1, 0.935])
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"\n  saved -> {out_path}")
    print("=" * 74)


if __name__ == "__main__":
    main()
