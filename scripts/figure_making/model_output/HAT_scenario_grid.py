#!/usr/bin/env python3
"""Every scenario, every preset, both periods, on one page.

THE FIGURE
    Periods as rows and source/sink presets as columns, left to right in
    order of increasing correction. Each panel draws one line per management
    scenario against the section 8 CoastSat target.

    NEITHER DIMENSION IS FIXED. The rows come from PERIOD_STARTS through
    HATTERAS_PERIODS -- the canonical chain is 1996-2010 and 2010-2024 since
    2026-09-17 -- and a preset only gets a column if be_rates() has it solved
    for EVERY period drawn. calibBE is solved for 1984 and 2004 only, so on
    the current chain it is dropped rather than drawn as a column that can
    never be filled.

    A cell with no run on disk draws nothing, so the script PRINTS the
    missing (period, preset, scenario) combinations: a sparse grid should
    read as runs not yet done, not as a result.

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
    python scripts/figure_making/model_output/HAT_scenario_grid.py
    python scripts/figure_making/model_output/HAT_scenario_grid.py --no-reloc --out FIG.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

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
                              FIG_W_DOUBLE, DOMAIN_AXIS_LABEL,
                              record_caption, save)
apply_style()          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402

_HERE = Path(__file__).resolve()
# Anchored by SEARCHING UPWARD for the project root rather than by
# counting parent directories (2026-09-13). A counted depth is correct
# only while the file stays where it was written, and these moved into
# subfolders of hatteras_ms. Six files here already did it this way.
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / 'pyproject.toml').exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live under scripts/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
# HAT_run_all lives in scripts/hatteras_ms/, named outright since this file
# moved to scripts/figure_making/model_output/ (2026-09-18); it used to be
# found as this file's grandparent.
for _path in (SCRIPTS_DIR, SCRIPTS_DIR / "hatteras_ms",
              SCRIPTS_DIR / "hatteras_ms" / "groin-sweep"):
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
from cascade_pipeline.run_layout import resolve             # noqa: E402
from cascade_pipeline.run_registry import preset_dir_for      # noqa: E402
# THE RUN DRIVER'S OWN GUARD. Which (period, scenario) pairs are distinct
# runs is decided by HAT_run_all.scenario_applies -- full_no_fill only exists
# where a fill is actually scheduled -- and this figure asks it rather than
# keeping a second copy of the rule that could disagree (2026-09-17).
# Importing is safe: HAT_run_all does its work under a __main__ guard.
from HAT_run_all import scenario_applies                      # noqa: E402
from site_layer.hatteras_site_config import (                            # noqa: E402
    HATTERAS_ANNOTATIONS,
    HATTERAS_DOMAINS,
    HATTERAS_PERIODS,
    resolve_be_preset,
)

RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
# Resolved through hat_observed_rates.py (2026-09-18), not typed.
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE  # noqa: E402
DEFAULT_OUT = (PROJECT_BASE_DIR / "output" / "comparisons"
               / "scenario_grid" / "scenario_grid_by_preset.png")
# The manuscript copy, with the other figures by subject (2026-09-18). Written
# only from a default run: an --out or any flagged variant is a working figure
# and must not overwrite it.
PUBLISHED = (PROJECT_BASE_DIR / "output" / "figures" / "shoreline"
             / "scenario_grid.png")

# THE CANONICAL CHAIN, 1996 -> 2010 -> 2024 (Hannah, 2026-09-17). Ends come
# from HATTERAS_PERIODS, so changing PERIOD_STARTS moves the whole figure.
PERIOD_STARTS = (1996, 2010)
PERIODS = tuple((st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS)

# ONLY PRESETS SOLVED FOR EVERY PERIOD DRAWN. calibBE is solved for 1984 and
# 2004 only, so on the new chain it would be a column that can never be
# filled -- not a gap in the runs but a preset that does not exist for those
# windows. Asking be_rates() is what decides, so this cannot go stale.
_WANTED = ("zeroBE", "edgeBE", "calibBE")


def _solved_everywhere(name):
    try:
        _canonical, by_period = resolve_be_preset(name)
    except Exception:
        return False
    return all(start in by_period for start, _end in PERIODS)


PRESETS = tuple(p for p in _WANTED if _solved_everywhere(p))
_DROPPED = tuple(p for p in _WANTED if p not in PRESETS)

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
    # Resolved rather than joined: runs forced off the calibration wave
    # climate sit under an arm component this join had no slot for. The
    # default arm is the calibration one, which is what this grid draws.
    preset_dir = preset_dir_for(RAW_RUNS, (period_start, period_end), preset)
    if not preset_dir.is_dir():
        return {}

    series = {}
    for run_dir in sorted(preset_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        scenario, reloc = classify(run_dir.name)
        if scenario is None:
            continue
        rate_csv = resolve(run_dir, "rate_csv", run_dir.name)
        if not rate_csv.is_file():
            print(f"  ! no rate CSV in {run_dir.name}")
            continue
        frame = pd.read_csv(rate_csv).set_index("gis_domain")
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
    # The end year comes from HATTERAS_PERIODS, not from start + 20: the older
    # pair happened to be 20-year windows, the canonical chain is 1996-2010
    # and 2010-2024, both 14 (2026-09-17).
    end = HATTERAS_PERIODS[period_start]["end_year"]
    csv_path = COASTSAT_BASE / f"{period_start}_{end}" / "transect_lrr_full.csv"
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
        len(PERIODS), len(PRESETS), figsize=figsize("double", height=3.74),
        sharex=True, squeeze=False)

    drawn_scenarios = set()
    drawn_cells = set()
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
                    drawn_cells.add((start, preset, scenario))
                    any_reloc = any_reloc or reloc

            ax.axhline(0.0, color="0.55", linewidth=0.8, zorder=1)
            if row == 0:
                ax.set_title(preset, fontweight="bold", pad=4)
            if column == 0:
                # short: the full wording is one shared y label below
                ax.set_ylabel(f"{start}–{end}", fontweight="bold")
            # one shared x label for the whole grid, set after the loop
            ax.set_xlim(HATTERAS_DOMAINS.first_gis_id,
                        HATTERAS_DOMAINS.last_gis_id)
            ax.grid(True, alpha=0.25, linewidth=0.6)

        stacked = np.concatenate([np.asarray(v, dtype=float)
                                  for v in row_values])
        row_ranges.append(stacked[np.isfinite(stacked)])

        if not args.no_annotations:
            for ax in row_axes:
                # bands and lines only: six panels cannot each carry the
                # eight annotation names at this width
                add_geographic_annotations(ax, HATTERAS_ANNOTATIONS,
                                           label=False)

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
    # The strip the legend needs depends on how many scenarios actually had
    # runs, so it is measured rather than fixed: with only full_management
    # on disk the key is one row, not three.
    _leg_rows = math.ceil(len(handles) / 3)
    figure.legend(handles=handles, loc="lower center", ncol=3,
                  frameon=False, bbox_to_anchor=(0.5, 0.004))

    # The title and the y-axis note used to be drawn here. Both are caption
    # material under figure_making/STYLE.md, and on a 190 mm six-panel grid they
    # were also the two widest things on the page (2026-09-17).
    # WHAT IS MISSING, SAID OUT LOUD. A cell with no run draws nothing, and a
    # near-empty grid looks like a result rather than an absence of runs
    # (2026-09-17: the move to 1996/2010 left most scenarios unrun).
    if _DROPPED:
        print(f"  presets not solved for {[p[0] for p in PERIODS]}, column "
              f"dropped: {', '.join(_DROPPED)}")
    # A cell only counts if the driver says it is a DISTINCT run: a scenario
    # that collapses onto another in this period is not a gap, and reporting
    # it as one would leave the figure permanently claiming missing work.
    _expected, _degenerate = [], []
    for (start, _end) in PERIODS:
        for scen in SCENARIO_LABEL:
            applies, why = scenario_applies(start, scen)
            if not applies:
                _degenerate.append(f"{start}/{scen} -- {why}")
                continue
            for preset in PRESETS:
                _expected.append((start, preset, scen))
    _missing = [f"{start}/{preset}/{scen}"
                for (start, preset, scen) in _expected
                if (start, preset, scen) not in drawn_cells]
    for note in _degenerate:
        print(f"  not a cell: {note}")
    if _missing:
        print(f"  {len(_missing)} of {len(_expected)} cells that SHOULD exist "
              f"have no run on disk:")
        for cell in _missing:
            print(f"    {cell}")
    else:
        print(f"  all {len(_expected)} cells present")

    figure.supxlabel(DOMAIN_AXIS_LABEL)
    figure.supylabel(f"shoreline change rate, {RATE_LABEL} (m/yr)")

    figure.tight_layout(rect=[0.015, 0.075 + 0.042 * _leg_rows, 1, 0.99])
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    # NO TIGHT BBOX: it trims to content, and with the labels hanging outside
    # the axes this saved at 14.20 in wide however figsize() was set -- nearly
    # double the column. tight_layout above already reserves the margins.
    figure.savefig(out_path, dpi=300)
    print(f"\n  saved -> {out_path}")
    _variant = args.no_reloc or args.no_annotations or args.y_per_row
    if out_path.resolve() == DEFAULT_OUT.resolve() and not _variant:
        save(figure, PUBLISHED)
        periods = " and ".join(f"{a}–{b}" for a, b in PERIODS)
        record_caption(PUBLISHED, (
            f"Modelled shoreline change rate by domain for every management "
            f"scenario, {periods}, under each source/sink preset. Periods are "
            f"rows and presets ({', '.join(PRESETS)}) columns, in order of "
            f"increasing correction. Each panel draws one line per scenario "
            f"against the CoastSat LOESS target in black; scenarios are a "
            f"single-hue ramp ordered by management intensity, and each "
            f"relocation arm is dashed in its non-relocation twin's colour "
            f"because the two overlap at this scale. The y axis is shared "
            f"across every panel. Reading across a row shows what the "
            f"source/sink term does; the spread within a panel shows what "
            f"management does. The model rate is lrr_m_yr, the OLS slope the "
            f"runs are scored with; rates are m/yr, seaward positive. An "
            f"empty panel is a run not yet made, not a result."))
        print(f"  saved -> {PUBLISHED}")
    print("=" * 74)


if __name__ == "__main__":
    main()
