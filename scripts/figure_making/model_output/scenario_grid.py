#!/usr/bin/env python3
"""
Every management scenario, every solved preset, both periods, on one page against the CoastSat target.

    python scripts/figure_making/model_output/scenario_grid.py [--out PATH]

Periods as rows, presets as columns (only presets solved for every period);
reads the nogroin matrix runs. Writes output/comparisons/scenario_grid/ and the
manuscript copy to output/figures/5-results/. Details: scripts/figure_making/model_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# House style (site_layer/hat_figure_style.py), applied at import
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
# Project root found by searching upward (ORGANIZATION.md rule 5)
PROJECT_BASE_DIR = next(_p for _p in _HERE.parents if (_p / 'pyproject.toml').exists())
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no "
        f"pyproject.toml. This file expects to live under scripts/.")
SCRIPTS_DIR = PROJECT_BASE_DIR / "scripts"
# HAT_run_all lives in scripts/hatteras_ms/
for _path in (SCRIPTS_DIR, SCRIPTS_DIR / "hatteras_ms",
              SCRIPTS_DIR / "hatteras_ms" / "groin-sweep"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cascade_pipeline.annotations import (                    # noqa: E402
    add_geographic_annotations,
)
from cascade_pipeline.coastsat_lowess import (                 # noqa: E402
    CoastSatDataset,
    LowessConfig,
    build_coastsat_series,
)
from cascade_pipeline.hindcast import build_target_table      # noqa: E402
from cascade_pipeline.run_layout import resolve             # noqa: E402
from cascade_pipeline.run_registry import preset_dir_for      # noqa: E402
# Which scenarios are distinct runs is the run driver's rule, asked, not copied
from HAT_run_all import scenario_applies                      # noqa: E402
from site_layer.hatteras_site_config import (                            # noqa: E402
    HATTERAS_ANNOTATIONS,
    HATTERAS_DOMAINS,
    HATTERAS_PERIODS,
    resolve_be_preset,
)

# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
# Resolved through hat_observed_rates.py (2026-09-18), not typed.
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE  # noqa: E402
from site_layer import hat_figure_style as _hs  # noqa: E402
DEFAULT_OUT = _hs.COMPARISONS_ROOT / "scenario_grid" / "scenario_grid_by_preset.png"
# The manuscript copy, written only from a default run
PUBLISHED = _hs.figure_dir("results") / "scenario_grid.png"

# The canonical chain, 1996 -> 2010 -> 2024; ends come from HATTERAS_PERIODS
PERIOD_STARTS = (1996, 2010)
PERIODS = tuple((st, HATTERAS_PERIODS[st]["end_year"]) for st in PERIOD_STARTS)

# Only presets solved for every period drawn; be_rates() decides
_WANTED = ("zeroBE", "edgeBE", "calibBE")


# Is this preset's end rate solved for every period drawn?
def _solved_everywhere(name):
    try:
        _canonical, by_period = resolve_be_preset(name)
    except Exception:
        return False
    return all(start in by_period for start, _end in PERIODS)


PRESETS = tuple(p for p in _WANTED if _solved_everywhere(p))
_DROPPED = tuple(p for p in _WANTED if p not in PRESETS)

# Section 8's LOWESS settings, as the runner scores
LOWESS_CONFIG = LowessConfig(window_domains=(7,), skip_southern_domains=10)

# Model column drawn: lrr_m_yr, the OLS slope matching the CoastSat LRR
RATE_COLUMN = "lrr_m_yr"
RATE_LABEL = ("LRR" if RATE_COLUMN == "lrr_m_yr"
              else "endpoint difference")

TARGET_WINDOW = 7   # 10 until 2026-09-28, with the runner

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
# Ramp starts at 0.35 so the lightest scenario stays readable on white
SCENARIO_COLOUR = {
    name: RAMP(0.35 + 0.62 * i / (len(SCENARIO_ORDER) - 1))
    for i, name in enumerate(SCENARIO_ORDER)
}
# -----------------------------------------------------------------------------


# (scenario, relocations) from a run name's switch tokens
def classify(run_name):
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


# Every nogroin run of one period and preset, keyed by (scenario, reloc)
def load_runs(period_start, period_end, preset):
    # Resolved through run_registry, never joined by hand
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
        # lrr_m_yr where the run has it, change_rate_m_yr otherwise (and says so)
        if RATE_COLUMN in frame.columns:
            series[(scenario, reloc)] = frame[RATE_COLUMN]
        else:
            print(f"  ! {run_dir.name}: no {RATE_COLUMN}; falling back "
                  f"to change_rate_m_yr (endpoint difference)")
            series[(scenario, reloc)] = frame["change_rate_m_yr"]
    return series


# The section 8 CoastSat target for one period, by domain
def load_target(period_start):
    # End year from HATTERAS_PERIODS, not start + 20
    end = HATTERAS_PERIODS[period_start]["end_year"]
    csv_path = COASTSAT_BASE / f"{period_start}_{end}" / "transect_lrr_full.csv"
    built = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {period_start}",
                         period_start=period_start, csv_path=str(csv_path))],
        period_start, LOWESS_CONFIG)
    if not built:
        raise FileNotFoundError(f"CoastSat transects failed to load: "
                                f"{csv_path}")
    table = build_target_table(built[0], LOWESS_CONFIG, HATTERAS_DOMAINS,
                               TARGET_WINDOW)
    return pd.Series(np.asarray(table["target_lrr_m_yr"], dtype=float),
                     index=np.asarray(table["gis_domain"], dtype=int))


# Run: load runs and targets, draw the grid, name what is missing, save both copies
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
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
                    label=f"CoastSat target (LOWESS {TARGET_WINDOW}-domain)")
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
                # Bands and lines only: no room for the annotation names in six panels
                add_geographic_annotations(ax, HATTERAS_ANNOTATIONS,
                                           label=False)

    # Applied after every panel: the shared range pools the whole figure
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
    # Legend strip measured from how many scenarios had runs
    _leg_rows = math.ceil(len(handles) / 3)
    figure.legend(handles=handles, loc="lower center", ncol=3,
                  frameon=False, bbox_to_anchor=(0.5, 0.004))

    # Title and axis note are the caption; missing cells are named, not left blank (README)
    if _DROPPED:
        print(f"  presets not solved for {[p[0] for p in PERIODS]}, column "
              f"dropped: {', '.join(_DROPPED)}")
    # Only distinct runs count as missing; collapsed scenarios are not gaps
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
    # No tight bbox: it trimmed to content and doubled the width
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
            f"against the CoastSat LOWESS target in black; scenarios are a "
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
