#!/usr/bin/env python3
"""
Figures for the hindcast sensitivity sweep: skill per axis, alongshore rates per cell, road outcomes.

    python scripts/sensitivity_analysis/plot_sensitivity.py --start-year 1996 --preset edgeBE
    python scripts/sensitivity_analysis/plot_sensitivity.py --start-year 1996 --quantity position --reference total
    python scripts/sensitivity_analysis/plot_sensitivity.py --start-year 1996 --quantity position --reference net_change
    python scripts/sensitivity_analysis/plot_sensitivity.py --start-year 1996 --quantity position --reference net_change --no-groin
    python scripts/sensitivity_analysis/plot_sensitivity.py --start-year 1996 --circularity

Reads the sweep manifest, run_index.csv and each cell's rates; scores against
the production target (interior GIS 2-89). Since 2026-10-26 the folder
output/raw_runs/sensitivity/figures/<start>_<end>_<preset>/ holds the position-change
set (--quantity position --reference net_change: the run figure's CoastSat net change
and LOWESS-smoothed model, skill in metres); the rate set goes to its rate/ subfolder.
Details: scripts/sensitivity_analysis/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

# House style (site_layer/hat_figure_style.py), applied at import
import sys as _sys
from pathlib import Path as _HP
_sys.path.insert(0, str(next(_q for _q in _HP(__file__).resolve().parents
                             if (_q / "hatteras_paths.py").exists())))
from site_layer.hat_figure_style import (apply_style, figsize,  # noqa: E402
                              FIG_W_DOUBLE, town_bands)
apply_style()
import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_BASE_DIR  # noqa: E402
if not (PROJECT_BASE_DIR / "data" / "hatteras_init").is_dir():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (_HAT_DIR,
              _HAT_DIR / "hatteras_ms" / "hindcast",
              _HAT_DIR / "sensitivity_analysis"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS, HATTERAS_PERIODS)
from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, LowessConfig, build_coastsat_series, scale_coastsat_series)
from cascade_pipeline.shoreline import compute_change_rate  # noqa: E402
from cascade_pipeline.hindcast import build_target_table  # noqa: E402
from cascade_pipeline.run_layout import resolve as resolve_run_file  # noqa: E402
from cascade_pipeline.run_registry import (  # noqa: E402
    MATRIX_KIND, find_run_dir, scan_runs, sweep_family)
from cascade_pipeline.plotting.rate_comparison import (  # noqa: E402
    DEFAULT_RATE_COMPARISON)
from cascade_pipeline.plotting.net_change_comparison import smooth_model  # noqa: E402
from site_layer.hat_observed_rates import net_change_dir  # noqa: E402
from hindcast_sensitivity import SWEEPS, normalise  # noqa: E402
from HAT_hindcast_config import field_default  # noqa: E402

# The value this sweep's base runs used for one setting
def sweep_base_value(setting):
    return field_default(setting)

# Reading order, most informative first; also sets the filename prefixes
AXIS_ORDER = ("wave_height", "wave_period", "wave_asymmetry",
              "wave_angle_high_fraction", "relocation_setback")

# --- CONFIG ------------------------------------------------------------------
OUT_ROOT = PROJECT_BASE_DIR / "output" / "calibration" / "sensitivity"
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
# Figures beside the cells they describe (since 2026-09-28)
FIGURES_ROOT = RAW_RUNS / "sensitivity" / "figures"
# -----------------------------------------------------------------------------
# Observed rates from where the runner reads them
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT as COASTSAT_BASE_DIR  # noqa: E402

# The runner's LOWESS target; check_target_matches() asserts it per run
LOWESS_CONFIG = LowessConfig()
TARGET_WINDOW = max(LOWESS_CONFIG.window_domains)


# The "<start>_<end>" path component of a period, from the period table
def period_component(start_year):
    start_year = int(start_year)
    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    return f"{start_year}_{end_year}"

# Model curves: a warm sequential ramp, since the observed layer owns the blues
MODEL_RAMP = plt.get_cmap("YlOrRd")
RAMP_LO, RAMP_HI = 0.30, 0.92          # skip the near-white and near-black ends
BASELINE_COLOR = "#111111"
INK, INK_MUTED = "#222222", "#666666"

# The values the alongshore figure draws per axis, the current setting included (advisor, 2026-10-27);
# an axis not named here draws every cell
ALONGSHORE_VALUES = {"wave_height": (2.0, 3.0, 4.0, 5.0)}


# Reading the sweep

# Every completed cell for one period, preset and groin state, newest wins
def load_cells(start_year, preset, groin=True):
    manifest = OUT_ROOT / f"sensitivity_{start_year}.jsonl"
    if not manifest.exists():
        raise FileNotFoundError(f"no sweep manifest at {manifest}")

    rows = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if row.get("preset") != preset or not row.get("ok"):
            continue
        # Groin-on and groin-off cells share values; mixed, the newest would win (2026-10-28)
        if row.get("groin", True) != groin:
            continue
        # The manifest's absolute path goes stale; re-resolve from the run name
        run_dir = Path(row["detail"])
        family = sweep_family(run_dir.name)
        if not family:
            # A 2026-09-01..16 cell filed under the matrix run's own name: skip it
            continue
        # Resolved through the registry; the run name is the stable handle
        try:
            run_dir = find_run_dir(RAW_RUNS, run_dir.name,
                                   period_component(row["start_year"]),
                                   row["preset"], kind="sensitivity", tag=family)
        except FileNotFoundError:
            if not run_dir.is_dir():
                raise
        rows.append(dict(sweep=row["sweep"], setting=row["setting"],
                         value=row["value"], run_dir=run_dir,
                         run_name=run_dir.name, kind="sensitivity", tag=family,
                         base_name=baseline_name(run_dir.name), groin=groin))
    if not rows:
        return pd.DataFrame(columns=["sweep", "setting", "value", "sort_key",
                                     "run_dir", "run_name", "kind", "tag",
                                     "base_name", "groin", "key", "base_key"])

    frame = pd.DataFrame(rows)
    # Keyed on (run_name, kind, tag): a cell and its baseline are two rows
    frame["key"] = list(zip(frame.run_name, frame.kind, frame.tag))
    frame["base_key"] = list(zip(frame.base_name, [MATRIX_KIND] * len(frame),
                                 [""] * len(frame)))
    frame["norm"] = frame["value"].map(normalise)
    # `measured` sorts last as its own category, not a fake number
    frame["sort_key"] = frame["norm"].map(
        lambda v: np.inf if v is None else v)
    frame = frame.drop_duplicates(subset=["sweep", "norm"], keep="last")
    return frame.sort_values(["sweep", "sort_key"]).reset_index(drop=True)


# The values of one axis whose cells failed, for one period, preset and groin state
def failed_values(start_year, preset, sweep, groin=True):
    manifest = OUT_ROOT / f"sensitivity_{start_year}.jsonl"
    failed = set()
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        if (row.get("preset") == preset and row.get("sweep") == sweep
                and row.get("groin", True) == groin and not row.get("ok")):
            failed.add(normalise(row["value"]))
    return failed


# The matrix run a cell is a departure from
def baseline_name(run_name):
    stem, _, token = run_name.rpartition("_")
    if token.startswith("wave") or token.startswith("rset"):
        return stem
    return run_name


# Every run's row, indexed by (run_name, kind, tag)
def load_index():
    frame = scan_runs(RAW_RUNS)
    text_columns = {"run_name", "kind", "tag", "status", "timestamp",
                    "source_sink_preset", "scenario", "arm", "git_commit",
                    "topo_product", "topo_dune_version",
                    "island_offset_version", "rate_estimator",
                    "be_values_digest"}
    for column in frame.columns:
        if column not in text_columns:
            try:
                frame[column] = pd.to_numeric(frame[column])
            except (ValueError, TypeError):
                pass
    keys = list(zip(frame.run_name, frame.kind, frame.tag))
    repeated = sorted({k for k in keys if keys.count(k) > 1})
    if repeated:
        raise ValueError(
            f"run_index.csv has more than one row for {repeated}. The "
            f"(run_name, kind, tag) key should make that impossible; the index "
            f"needs repairing, not deduplicating.")
    return frame.set_index(["run_name", "kind", "tag"])


# Per-domain modelled LRR for one run
def model_rates(run_dir, run_name):
    # Resolved, not joined: the rate CSV's path depends on the run layout
    path = resolve_run_file(run_dir, "rate_csv", run_name)
    frame = pd.read_csv(path)
    return frame[["gis_domain", "lrr_m_yr"]]


# Per-domain modelled shoreline position change, end minus start (m); `smooth` as the run figure draws it
def model_position_change(run_dir, run_name, smooth=False):
    matrix = np.load(resolve_run_file(run_dir, "matrix", run_name))
    change = compute_change_rate(matrix, span_years=1, flip_sign=True)
    if smooth:
        change = smooth_model(change)
    geometry = HATTERAS_DOMAINS.for_total(len(change))   # the run's own buffers, not the config's
    real = change[geometry.start_real_index:geometry.end_real_index]
    gis = np.arange(HATTERAS_DOMAINS.first_gis_id,
                    HATTERAS_DOMAINS.last_gis_id + 1)
    return pd.DataFrame({"gis_domain": gis, "lrr_m_yr": real})


# The full-record rate behind the projected change (2026-09-19 advisor target)
LONG_TERM_WINDOW = (1996, 2024)


# The observed layer in metres
def position_layers(start_year, reference):
    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    span = end_year - start_year
    if reference == "net_change":
        # The DEM-to-DEM target the run's shoreline_position_change figure draws (net_change_comparison)
        tag = f"{start_year}_{end_year}"
        series = build_coastsat_series(
            [CoastSatDataset(
                label=f"CoastSat net change ({start_year}-{end_year})", period_start=start_year,
                csv_path=str(net_change_dir(start_year, end_year)
                             / f"transect_net_change_{tag}.csv"),
                rate_col="net_change_m")],
            active_period_start=start_year, lowess_config=LOWESS_CONFIG,
            domains=HATTERAS_DOMAINS)
        active = next(cs for cs in series if cs["active"])
        target = build_target_table(active, LOWESS_CONFIG, HATTERAS_DOMAINS,
                                    TARGET_WINDOW)
        label = (f"CoastSat net change {start_year}–{end_year}, {TARGET_WINDOW}-domain "
                 f"LOWESS (domain means D1–{LOWESS_CONFIG.skip_southern_domains})")
        return series, target, label
    if reference == "total":
        series, _ = coastsat_layers(start_year)
        fit = f"{start_year}–{end_year}"
        name = "total change"
    else:
        lo, hi = LONG_TERM_WINDOW
        series = build_coastsat_series(
            [CoastSatDataset(
                label=f"CoastSat LRR ({lo}-{hi})", period_start=lo,
                csv_path=str(COASTSAT_BASE_DIR / f"{lo}_{hi}"
                             / "transect_lrr_full.csv"))],
            active_period_start=lo, lowess_config=LOWESS_CONFIG,
            domains=HATTERAS_DOMAINS)
        fit = f"{lo}–{hi}"
        name = "projected change"
    active = next(cs for cs in series if cs["active"])
    target = build_target_table(active, LOWESS_CONFIG, HATTERAS_DOMAINS,
                                TARGET_WINDOW)
    target = target.assign(target_lrr_m_yr=target.target_lrr_m_yr * span)
    scaled = scale_coastsat_series(series, span, active=True if reference ==
                                   "projected" else None)
    label = (f"CoastSat {name}, LRR {fit} × {span} yr, {TARGET_WINDOW}-domain "
             f"LOWESS (domain means D1–{LOWESS_CONFIG.skip_southern_domains})")
    return scaled, target, label


# Assert every run was scored against the target this figure draws
def check_target_matches(index, keys):
    expected = f"CoastSat LOWESS {TARGET_WINDOW}-domain"
    for name, kind, tag in keys:
        row = index.loc[(name, kind, tag)]
        # Resolved in the arm the row claims; find_run_dir raises rather than skipping
        directory = find_run_dir(
            RAW_RUNS, name, (int(row.start_year), int(row.end_year)),
            row.source_sink_preset, kind, tag)
        meta = resolve_run_file(directory, "metadata_json", name)
        # No metadata means no recorded target: the only allowed skip
        if not meta.exists():
            continue
        target = json.loads(meta.read_text(encoding="utf-8")).get(
            "skill", {}).get("target")
        if target and target != expected:
            raise ValueError(
                f"{name} was scored against {target!r} but these figures draw "
                f"{expected!r}; the plotted curve would not be the one its "
                f"RMSE refers to")


# The observed layer

# The CoastSat series and the spliced target, built the hindcast's way
def coastsat_layers(start_year):
    datasets = [
        CoastSatDataset(
            label="CoastSat LRR (1984-2004)", period_start=1984,
            csv_path=str(COASTSAT_BASE_DIR / "1984_2004"
                         / "transect_lrr_full.csv")),
        CoastSatDataset(
            label="CoastSat LRR (2004-2024)", period_start=2004,
            csv_path=str(COASTSAT_BASE_DIR / "2004_2024"
                         / "transect_lrr_full.csv")),
        # The calibration and test windows, built from the period table as
        # the runner builds them (section 8.1). Literal 1996_2010 / 2010_2024
        # here drew the 1996 sweep against a window it was not scored on and
        # had no 2009 dataset at all (2026-10-15).
        *[CoastSatDataset(
            label=f"CoastSat LRR ({period_component(start).replace('_', '-')})",
            period_start=start,
            csv_path=str(COASTSAT_BASE_DIR / period_component(start)
                         / "transect_lrr_full.csv"))
          for start in (1996, 2009)],
    ]
    series = build_coastsat_series(
        datasets, active_period_start=start_year, lowess_config=LOWESS_CONFIG,
        domains=HATTERAS_DOMAINS)
    active = next((cs for cs in series if cs["active"]), None)
    if active is None:
        raise RuntimeError(f"no CoastSat dataset starts at {start_year}")
    target = build_target_table(
        active, LOWESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)
    return series, target


# House style

# One place, applied once at import, so every figure in the folder matches

HOUSE_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8,
    "figure.titlesize": 11,
    # Two spines, not four
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#222222",
    "text.color": "#222222",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.6,
    "ytick.minor.width": 0.6,
    "axes.grid": True,
    "grid.color": "#C8CDD2",
    "grid.linewidth": 0.5,
    "grid.alpha": 0.55,
    "axes.axisbelow": True,          # data over the grid, never under it
    "legend.frameon": False,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.04,
}
plt.rcParams.update(HOUSE_STYLE)

PANEL_LETTERS = "abcdefghijklmnopqrstuvwxyz"


# A lower-case panel letter outside the axes, top left
def panel_label(ax, letter, dx=-0.02):
    ax.text(dx, 1.04, f"({letter})", transform=ax.transAxes,
            fontsize=9.5, fontweight="bold", va="bottom", ha="right",
            color="#222222")


# The per-axes cleanup every panel gets
def tidy(ax, minor=True):
    if minor:
        ax.minorticks_on()
        ax.tick_params(which="minor", length=2)
        ax.grid(which="minor", visible=False)
    ax.tick_params(which="major", length=3.5)


# The community zones (grey, named at the top) and shoal zones (hatched gold, named at the
# foot), from the site config, under the data; call after the x limits are set
def site_zones(ax):
    town_bands(ax, where="top", fontsize=7.5)
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    for name, (lo, hi) in HATTERAS_ANNOTATIONS.shoal_zones.items():
        ax.axvspan(lo - 0.5, hi + 0.5, facecolor=HATTERAS_ANNOTATIONS.color_shoal, alpha=0.10,
                   zorder=0, hatch="///", edgecolor=HATTERAS_ANNOTATIONS.color_shoal,
                   linewidth=0)
        ax.text((lo + hi) / 2.0, 0.03, name, transform=trans, ha="center", va="bottom",
                fontsize=7.5, color="#7A5800", style="italic", zorder=9,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8))


# Figures

# The groin state as a title suffix: the calibrated groin-on set is unmarked
def groin_text(groin):
    return "" if groin else ", groin off"


# A swept value as it appears in a legend or colourbar tick
def value_label(value):
    norm = normalise(value)
    return "measured" if norm is None else f"{norm:g}"


# `count` colours from the sequential model ramp, light to dark
def ramp_colors(count):
    if count == 1:
        return [MODEL_RAMP(RAMP_HI)]
    return [MODEL_RAMP(RAMP_LO + (RAMP_HI - RAMP_LO) * i / (count - 1))
            for i in range(count)]


# The x-axis label for one swept axis, with its unit
def axis_xlabel(sweep):
    units = SWEEPS[sweep]["units"]
    label = SWEEPS[sweep]["label"]
    return f"{label} ({units})" if units else label


# Interior RMSE and bias against the swept value, one column per axis
# `scores` maps a run key to (rmse, bias) in `unit`; None reads the index's rate skill
def plot_skill_overview(cells, index, start_year, preset, out_dir, scores=None,
                        unit="m/yr", note=""):
    sweeps = [s for s in AXIS_ORDER if s in set(cells["sweep"])]
    if not sweeps:
        return None
    if scores is None:
        def scores(key):
            return index.loc[key].rmse_interior_m_yr, index.loc[key].mean_bias_interior_m_yr
    unit_tex = {"m/yr": "m yr$^{-1}$"}.get(unit, unit)
    fig, axes = plt.subplots(
        2, len(sweeps),
        figsize=figsize("double",
                        height=5.6 * FIG_W_DOUBLE / (2.55 * len(sweeps) + 0.9)),
        squeeze=False,
        sharex="col", constrained_layout=True)

    footnotes = []
    for col, sweep in enumerate(sweeps):
        block = cells[cells.sweep == sweep]
        base_rmse, base_bias = scores(block.iloc[0].base_key)

        numeric = block[np.isfinite(block.sort_key)]
        # The calibration value is a point on the curve, not just a reference line
        default = normalise(sweep_base_value(SWEEPS[sweep]["setting"]))
        x = list(numeric.sort_key.to_numpy(dtype=float))
        rmse = [scores(k)[0] for k in numeric.key]
        bias = [scores(k)[1] for k in numeric.key]
        if default is not None:
            x.append(default)
            rmse.append(base_rmse)
            bias.append(base_bias)
        order = np.argsort(x)
        x = np.asarray(x)[order]
        rmse = np.asarray(rmse)[order]
        bias = np.asarray(bias)[order]

        for row, (values, base_value, label) in enumerate((
                (rmse, base_rmse, f"Interior RMSE ({unit})"),
                (bias, base_bias, f"Interior mean bias ({unit})"))):
            ax = axes[row][col]
            ax.plot(x, values, color=MODEL_RAMP(RAMP_HI), lw=1.5,
                    marker="o", ms=3.8, mec="white", mew=0.6, zorder=3,
                    label="swept value")
            ax.axhline(base_value, color=BASELINE_COLOR, lw=0.9, ls=(0, (4, 2)),
                       zorder=2, label="calibration run")
            if default is not None:
                ax.plot([default], [base_value], marker="D", ms=5.5,
                        color=BASELINE_COLOR, zorder=4,
                        label="calibration value")
            if row == 1:
                ax.axhline(0.0, color=INK_MUTED, lw=0.7, ls=":", zorder=1)

            # Floor the y span on the run's RMSE so a flat response looks flat
            floor_span = 0.12 * abs(base_rmse)
            lo, hi = ax.get_ylim()
            if hi - lo < floor_span:
                mid = 0.5 * (lo + hi)
                ax.set_ylim(mid - floor_span / 2, mid + floor_span / 2)

            # A margin so an endpoint calibration diamond is not clipped
            ax.margins(x=0.07)
            tidy(ax)
            panel_label(ax, PANEL_LETTERS[row * len(sweeps) + col])
            if col == 0:
                ax.set_ylabel(label)
            if row == 0:
                ax.set_title(SWEEPS[sweep]["label"], pad=14)
                ax.tick_params(labelbottom=False)
            else:
                units = SWEEPS[sweep]["units"]
                ax.set_xlabel(units if units else "fraction")

        # `measured` goes in the caption strip, not on the axes
        for key in block[~np.isfinite(block.sort_key)].key:
            footnotes.append(
                f"{SWEEPS[sweep]['label']}, measured: "
                f"RMSE {scores(key)[0]:.3f} {unit_tex}")

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.085), frameon=False)

    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    caption = (f"Parameter sensitivity, {start_year}–{end_year}, {preset}"
               f"{groin_text(bool(cells.groin.iloc[0]))}. "
               f"Interior domains (GIS 2–89).{note}")
    if footnotes:
        caption += "  " + "; ".join(footnotes) + "."
    fig.suptitle(caption, y=1.035, fontsize=10)

    path = out_dir / "01_skill_overview.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# One axis of the sweep on `ax`: the observed layer, then one ramp-coloured line per value
# Returns the (handle, text) value-legend entries and how many of them are drawn lines
def draw_alongshore(ax, cells, sweep, start_year, preset, cs_series, target,
                    position=None, shown=None):
    block = cells[cells.sweep == sweep]
    base = block.iloc[0].base_name
    groin = bool(block.iloc[0].groin)
    # `shown`: the values to draw; one with no completed cell is reported, not drawn
    absent = []
    if shown is not None:
        wanted = {normalise(v) for v in shown}
        absent = sorted(wanted - set(block.norm)
                        - {normalise(sweep_base_value(SWEEPS[sweep]["setting"]))})
        if absent:
            print(f"  {sweep}: no completed cell at {', '.join(f'{v:g}' for v in absent)}; "
                  f"not drawn")
        block = block[block.norm.isin(wanted)]

    # The observed layer is one curve: the scoring target
    skip = LOWESS_CONFIG.skip_southern_domains
    if position is None:
        load = model_rates
    else:
        smooth = position[0] == "net_change"

        def load(run_dir, run_name):
            return model_position_change(run_dir, run_name, smooth=smooth)
    active = next(cs for cs in cs_series if cs["active"])
    south = np.asarray(active["transect_domains"]) <= skip
    ax.scatter(
        np.asarray(active["transect_along_coast"])[south]
        / HATTERAS_DOMAINS.domain_spacing_m + HATTERAS_DOMAINS.first_gis_id,
        np.asarray(active["transect_rates"])[south],
        color=DEFAULT_RATE_COMPARISON.raw_color, s=7, alpha=0.6,
        linewidths=0, zorder=2,
        label=f"CoastSat transects, D1–{skip}")
    # The target line runs through the raw-mean domains too
    ax.plot(target.gis_domain, target.target_lrr_m_yr, color="#08306B", lw=1.8,
            zorder=5,
            label=(position[1] if position else
                   f"CoastSat LRR, {TARGET_WINDOW}-domain LOWESS "
                   f"(domain means D1–{skip})"))

    # The baseline from the registry, not the cell's sibling folder
    base_dir = find_run_dir(RAW_RUNS, base, period_component(start_year),
                            preset, kind=MATRIX_KIND)
    base_rates = load(base_dir, base)
    # The current setting, named with its value
    current = normalise(sweep_base_value(SWEEPS[sweep]["setting"]))
    units = SWEEPS[sweep]["units"]
    unit = f" {units.split()[0]}" if units else ""
    # (sort key, curve, legend text); the current setting is one value among the rest
    curves = [(np.inf if current is None else current, base_rates,
               f"{value_label(current)}{unit} (current)")]
    for _, cell in block.iterrows():
        curves.append((cell.sort_key, load(cell.run_dir, cell.run_name),
                       f"{value_label(cell.value)}{unit}"))
    curves.sort(key=lambda c: c[0])

    # One ramp colour per value, in value order, every line the same weight (2026-10-27)
    entries = []
    for color, (_, rates, text) in zip(ramp_colors(len(curves)), curves):
        line, = ax.plot(rates.gis_domain, rates.lrr_m_yr, color=color, lw=1.4,
                        alpha=0.95, zorder=6)
        entries.append((line, text))
    # A value whose run failed is named in the legend, greyed, with no line
    # Every failed value on an unfiltered axis; only the selected ones on a filtered one
    failed = failed_values(start_year, preset, sweep, groin) - set(block.norm) - {current}
    if shown is not None:
        failed &= set(absent)
    for value in sorted(failed):
        entries.append((Line2D([], [], color="none"), f"{value:g}{unit}: run failed"))

    ax.axhline(0.0, color=INK_MUTED, lw=0.7, ls="--", zorder=1)
    ax.set_xlim(HATTERAS_DOMAINS.first_gis_id - 0.5,
                HATTERAS_DOMAINS.last_gis_id + 0.5)
    site_zones(ax)
    ax.set_ylabel("Shoreline change rate, LRR (m/yr)" if position is None
                  else "Shoreline position change (m)")
    tidy(ax)
    return entries, len(curves)


# The value legend to the right of `ax`: one coloured line per value, failed runs greyed
def value_legend(ax, sweep, entries, drawn):
    values = ax.legend([e[0] for e in entries], [e[1] for e in entries],
                       title=SWEEPS[sweep]["label"], loc="upper left",
                       bbox_to_anchor=(1.01, 1.0), frameon=False, handlelength=2.4,
                       alignment="left")
    for text in values.get_texts()[drawn:]:
        text.set_color(INK_MUTED)
        text.set_style("italic")


# The title tail naming the observed quantity
def observed_text(position):
    if position is None:
        return ""
    if position[0] == "net_change":
        return ": position change vs CoastSat net change"
    return f": position change vs {position[0]} change"


# Modelled LRR per domain for every cell of one axis, over the observed layer
def plot_alongshore(cells, sweep, start_year, preset, cs_series, target,
                    out_dir, position=None, shown=None):
    block = cells[cells.sweep == sweep]
    if block.empty:
        return None
    groin = bool(block.iloc[0].groin)
    fig, ax = plt.subplots(figsize=(9.6, 4.5), constrained_layout=True)
    entries, drawn = draw_alongshore(ax, cells, sweep, start_year, preset, cs_series,
                                     target, position=position, shown=shown)
    ax.set_xlabel("Alongshore position (GIS domain, south → north)")

    # Legend below the axes: the observed curve, then the transect dots
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: "LOWESS" not in labels[i])
    reference = ax.legend([handles[i] for i in order], [labels[i] for i in order],
                          loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
                          frameon=False, handlelength=2.4, columnspacing=2.0)
    ax.add_artist(reference)

    # Value legend on the right, in value order
    value_legend(ax, sweep, entries, drawn)

    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    ax.set_title(f"{SWEEPS[sweep]['label']} sensitivity, "
                 f"{start_year}–{end_year}, {preset}{groin_text(groin)}"
                 + observed_text(position), pad=8)

    path = out_dir / f"{AXIS_ORDER.index(sweep) + 2:02d}_alongshore_{sweep}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# Every swept axis stacked on one alongshore axis, a panel per wave parameter (2026-10-28)
def plot_alongshore_panels(cells, start_year, preset, cs_series, target, out_dir,
                           position=None):
    sweeps = [s for s in AXIS_ORDER if s in set(cells.sweep)]
    if len(sweeps) < 2:
        return None
    groin = bool(cells.iloc[0].groin)
    fig, axes = plt.subplots(len(sweeps), 1, figsize=(9.6, 2.9 * len(sweeps)),
                             sharex=True, constrained_layout=True)
    towns = set(HATTERAS_ANNOTATIONS.town_spans)
    for i, (ax, sweep) in enumerate(zip(axes, sweeps)):
        entries, drawn = draw_alongshore(ax, cells, sweep, start_year, preset, cs_series,
                                         target, position=position,
                                         shown=ALONGSHORE_VALUES.get(sweep))
        value_legend(ax, sweep, entries, drawn)
        panel_label(ax, PANEL_LETTERS[i])
        # Community names on the top panel only; the grey bands stay on every panel
        if i:
            for text in ax.texts:
                if text.get_text() in towns:
                    text.set_visible(False)
    axes[-1].set_xlabel("Alongshore position (GIS domain, south → north)")

    # The observed layer once, below the last panel
    handles, labels = axes[-1].get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: "LOWESS" not in labels[i])
    fig.legend([handles[i] for i in order], [labels[i] for i in order],
               loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=2,
               frameon=False, handlelength=2.4, columnspacing=2.0)

    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    fig.suptitle(f"Wave climate sensitivity, {start_year}–{end_year}, {preset}"
                 f"{groin_text(groin)}{observed_text(position)}")

    path = out_dir / "06_alongshore_all_parameters.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


# Road outcomes against the relocation target
def plot_relocation_outcomes(cells, index, start_year, preset, out_dir):
    block = cells[cells.sweep == "relocation_setback"]
    if block.empty:
        return None
    base_row = index.loc[block.iloc[0].base_key]
    default = normalise(field_default("relocation_setback_m"))

    labels = [value_label(v) for v in block.value] + [f"{default:g}"]
    drowned = [index.loc[k].roads_drowned for k in block.key] + \
              [base_row.roads_drowned]
    blocked = [index.loc[k].roads_reloc_blocked for k in block.key] + \
              [base_row.roads_reloc_blocked]
    keys = [np.inf if normalise(v) is None else normalise(v)
            for v in block.value] + [default]
    order = np.argsort(keys)
    labels = [labels[i] for i in order]
    drowned = [drowned[i] for i in order]
    blocked = [blocked[i] for i in order]
    is_base = [abs(keys[i] - default) < 1e-9 for i in order]

    # The relocation panel kept short: it is usually all zeros
    fig, axes = plt.subplots(2, 1, figsize=(6.2, 4.2), sharex=True,
                             height_ratios=[2, 1], constrained_layout=True)
    for row, (ax, values, label) in enumerate((
            (axes[0], drowned, "NC-12 domains drowned"),
            (axes[1], blocked, "Relocations blocked"))):
        colors = [BASELINE_COLOR if b else MODEL_RAMP(RAMP_HI) for b in is_base]
        ax.bar(range(len(values)), values, color=colors, width=0.6, zorder=3)
        top = max(max(values), 1)
        for i, v in enumerate(values):
            ax.annotate(f"{int(v)}", (i, v), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=8.5)
        # Headroom so the value labels never touch the axes frame or the title.
        ax.set_ylim(0, top * 1.35 + 0.35)
        # Counts are integers; 0.5 of a drowned road is not a thing.
        ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=4))
        ax.set_ylabel(label)
        ax.grid(axis="x", visible=False)
        tidy(ax, minor=False)
        # Panel letter further out for the long y-labels
        panel_label(ax, PANEL_LETTERS[row], dx=-0.115)
        if not any(values):
            ax.annotate("none at any value", (0.5, 0.5),
                        xycoords="axes fraction", ha="center", va="center",
                        fontsize=8.5, color=INK_MUTED, style="italic")

    axes[1].set_xticks(range(len(labels)))
    axes[1].set_xticklabels(labels)
    for i, tick in enumerate(axes[1].get_xticklabels()):
        if is_base[i]:
            tick.set_fontweight("bold")
            tick.set_color(BASELINE_COLOR)
    if True in is_base:
        axes[1].annotate("calibration", (is_base.index(True), -0.36),
                         xycoords=("data", "axes fraction"), ha="center",
                         va="top", fontsize=7.5, color=BASELINE_COLOR)
    axes[1].set_xlabel("Relocation target (m behind the dune line)",
                       labelpad=18)
    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    fig.suptitle(f"Road outcomes vs relocation target, "
                 f"{start_year}–{end_year}, {preset}", fontsize=10)

    path = out_dir / "08_relocation_outcomes.png"
    fig.savefig(path)
    plt.close(fig)
    return path


# The Hs skill curve under both presets, which shows the circularity
def plot_circularity(start_year, index, out_dir):
    curves = {}
    for preset in ("calibBE", "edgeBE"):
        cells = load_cells(start_year, preset)
        block = cells[cells.sweep == "wave_height"] if not cells.empty \
            else cells
        if len(block) < 2:
            continue
        curves[preset] = (
            block.sort_key.to_numpy(dtype=float),
            np.array([index.loc[k].rmse_interior_m_yr for k in block.key]),
            index.loc[block.iloc[0].base_key],
        )
    if len(curves) < 2:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0),
                             constrained_layout=True)
    # Fixed colours per preset
    preset_colors = {"calibBE": "#B24502", "edgeBE": "#08519C"}
    default_hs = float(sweep_base_value("hs"))

    # Panel (b): each curve as a rise above its own minimum, to compare shapes
    for preset, (x, rmse, base_row) in curves.items():
        # The calibration Hs is a point on this curve
        xs = np.append(x, default_hs)
        ys = np.append(rmse, base_row.rmse_interior_m_yr)
        order = np.argsort(xs)
        xs, ys = xs[order], ys[order]
        for ax, values in ((axes[0], ys), (axes[1], ys - ys.min())):
            ax.plot(xs, values, color=preset_colors[preset], lw=1.6,
                    marker="o", ms=3.8, mec="white", mew=0.6,
                    label=preset, zorder=3)
            k = int(np.argmin(np.abs(xs - default_hs)))
            ax.plot([xs[k]], [values[k]], marker="D", ms=6,
                    color=preset_colors[preset],
                    markeredgecolor=BASELINE_COLOR, markeredgewidth=1.0,
                    zorder=4, label="_nolegend_")

    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    for ax, ylab, ttl in (
            (axes[0], "Interior RMSE (m/yr)", "Absolute skill"),
            (axes[1], "RMSE above each curve's own minimum (m/yr)",
             "Basin shape, magnitude removed")):
        ax.axvline(default_hs, color=INK_MUTED, lw=0.7, ls=":", zorder=1)
        ax.set_xlabel("Significant wave height, H$_s$ (m)")
        ax.set_ylabel(ylab)
        ax.set_title(ttl, pad=8)
        ax.margins(x=0.05)
        tidy(ax)
    axes[1].set_ylim(-0.02, 0.75)

    # Annotated at the top of panel (a), clear of both curves
    axes[0].annotate(f"calibration H$_s$ = {default_hs:g} m",
                     xy=(default_hs, 1.0), xycoords=("data", "axes fraction"),
                     xytext=(5, -4), textcoords="offset points",
                     fontsize=8, color=INK_MUTED, ha="left", va="top")

    # Legend below the axes; the diamond is explained in the caption
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, -0.07), frameon=False)
    fig.suptitle("Apparent optimum H$_s$ depends on the background-erosion "
                 f"treatment, {start_year}–{end_year}\n"
                 "Diamonds mark the calibration H$_s$",
                 fontsize=10, y=1.09)
    panel_label(axes[0], "a", dx=-0.10)
    panel_label(axes[1], "b", dx=-0.10)

    path = out_dir.parent / "_comparisons" / f"hs_circularity_{start_year}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return path


# One tidy CSV behind the figures
def write_summary(cells, index, start_year, preset, out_dir):
    rows = []
    for _, cell in cells.iterrows():
        row = index.loc[cell.key]
        base = index.loc[cell.base_key]
        rows.append(dict(
            start_year=start_year, preset=preset, axis=cell.sweep,
            setting=cell.setting, value=value_label(cell.value),
            run_name=cell.run_name, kind=cell.kind, tag=cell.tag,
            rmse_interior_m_yr=row.rmse_interior_m_yr,
            mean_bias_interior_m_yr=row.mean_bias_interior_m_yr,
            baseline_rmse_interior_m_yr=base.rmse_interior_m_yr,
            delta_rmse=row.rmse_interior_m_yr - base.rmse_interior_m_yr,
            roads_drowned=row.roads_drowned,
            roads_reloc_blocked=row.roads_reloc_blocked))
    frame = pd.DataFrame(rows)
    path = out_dir / "summary.csv"
    frame.to_csv(path, index=False)
    return path, frame


# Interior RMSE and mean bias (m) of the smoothed position change against the net-change target
def position_scores(cells, start_year, preset, target):
    obs = target.set_index("gis_domain")["target_lrr_m_yr"].loc[2:89]
    keys = {cell.key: (cell.run_dir, cell.run_name) for _, cell in cells.iterrows()}
    for base in set(cells.base_name):
        keys[(base, MATRIX_KIND, "")] = (find_run_dir(
            RAW_RUNS, base, period_component(start_year), preset, kind=MATRIX_KIND), base)
    out = {}
    for key, (run_dir, run_name) in keys.items():
        model = model_position_change(run_dir, run_name, smooth=True).set_index("gis_domain")
        err = model["lrr_m_yr"].loc[2:89] - obs
        out[key] = (float(np.sqrt((err ** 2).mean())), float(err.mean()))
    return out


# The position-change summary: one row per cell, skill in metres
def write_position_summary(cells, scores, start_year, preset, out_dir):
    rows = []
    for _, cell in cells.iterrows():
        rmse, bias = scores[cell.key]
        base_rmse, _ = scores[cell.base_key]
        rows.append(dict(
            start_year=start_year, preset=preset, axis=cell.sweep,
            setting=cell.setting, value=value_label(cell.value),
            run_name=cell.run_name, rmse_interior_m=rmse, mean_bias_interior_m=bias,
            baseline_rmse_interior_m=base_rmse, delta_rmse_m=rmse - base_rmse))
    path = out_dir / "summary.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


# Run: the figures for one period and preset, and the summary CSV
def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--start-year", type=int, default=1984,
                        choices=sorted(HATTERAS_PERIODS))
    parser.add_argument("--preset", default="edgeBE")
    parser.add_argument("--circularity", action="store_true",
                        help="also draw the calibBE-vs-edgeBE Hs panel")
    parser.add_argument("--out-dir", default=None)
    parser.add_argument("--quantity", choices=("rate", "position"),
                        default="rate",
                        help="position: model end-minus-start change (m) "
                             "against CoastSat LRR x span; alongshore figures "
                             "only, under figures/position_change/<reference>/")
    parser.add_argument("--reference", choices=("total", "projected", "net_change"),
                        default="total",
                        help="with --quantity position: the window's own LRR "
                             "(total), the 1996-2024 LRR (projected), or the "
                             "CoastSat net change the run figure draws (net_change: "
                             "the main folder, with skill and summary in metres)")
    parser.add_argument("--no-groin", dest="groin", action="store_false",
                        help="the groin-off cells, into <start>_<end>_<preset>_nogroin/")
    args = parser.parse_args()

    # One directory per (period, preset)
    end_year = HATTERAS_PERIODS[args.start_year].get(
        "end_year", args.start_year + 20)
    root = Path(args.out_dir) if args.out_dir else FIGURES_ROOT
    net_change = args.quantity == "position" and args.reference == "net_change"
    if args.quantity == "position" and not net_change:
        root = root / "position_change" / args.reference
    out_dir = root / (f"{args.start_year}_{end_year}_{args.preset}"
                      + ("" if args.groin else "_nogroin"))
    # The main folder holds the net-change position set; the rate set sits under rate/ (2026-10-26)
    if args.quantity == "rate":
        out_dir = out_dir / "rate"
    out_dir.mkdir(parents=True, exist_ok=True)

    index = load_index()
    cells = load_cells(args.start_year, args.preset, groin=args.groin)
    if cells.empty:
        print(f"no completed cells for {args.start_year} / {args.preset}"
              f"{groin_text(args.groin)}")
        return 1
    known = set(index.index)
    missing = sorted(set(cells.key) - known)
    if missing:
        raise ValueError(
            f"{len(missing)} cell(s) are in the manifest but not in "
            f"run_index.csv: {missing[:5]}. The index may have lost rows to a "
            f"concurrent write; re-run those cells before plotting.")
    absent_base = sorted(set(cells.base_key) - known)
    if absent_base:
        raise ValueError(
            f"baseline run(s) {absent_base} are not in run_index.csv; every "
            f"cell is drawn against its calibration-arm matrix run.")
    check_target_matches(index, list(cells.key))

    if args.quantity == "position":
        cs_series, target, obs_label = position_layers(args.start_year,
                                                       args.reference)
        written = [plot_alongshore(cells, sweep, args.start_year, args.preset,
                                   cs_series, target, out_dir,
                                   position=(args.reference, obs_label),
                                   shown=ALONGSHORE_VALUES.get(sweep))
                   for sweep in AXIS_ORDER if sweep in set(cells.sweep)]
        written.append(plot_alongshore_panels(cells, args.start_year, args.preset,
                                              cs_series, target, out_dir,
                                              position=(args.reference, obs_label)))
        if net_change:
            scores = position_scores(cells, args.start_year, args.preset, target)
            written.insert(0, plot_skill_overview(
                cells, index, args.start_year, args.preset, out_dir,
                scores=scores.__getitem__, unit="m",
                note=" Shoreline position change against the CoastSat net change, "
                     "both 7-domain LOWESS."))
            written.append(write_position_summary(cells, scores, args.start_year,
                                                  args.preset, out_dir))
        print(f"{len(cells)} cells  |  {args.start_year}  |  {args.preset}  |  "
              f"position change vs {args.reference}")
        for path in [p for p in written if p is not None]:
            print(f"  wrote {Path(path).name}")
        return 0

    cs_series, target = coastsat_layers(args.start_year)

    written = []
    written.append(plot_skill_overview(cells, index, args.start_year,
                                       args.preset, out_dir))
    for sweep in [s for s in AXIS_ORDER if s in set(cells.sweep)]:
        written.append(plot_alongshore(cells, sweep, args.start_year,
                                       args.preset, cs_series, target, out_dir,
                                       shown=ALONGSHORE_VALUES.get(sweep)))
    written.append(plot_alongshore_panels(cells, args.start_year, args.preset,
                                          cs_series, target, out_dir))
    written.append(plot_relocation_outcomes(cells, index, args.start_year,
                                            args.preset, out_dir))
    if args.circularity:
        written.append(plot_circularity(args.start_year, index, out_dir))
    summary_path, summary = write_summary(cells, index, args.start_year,
                                          args.preset, out_dir)
    written.append(summary_path)

    print(f"{len(cells)} cells  |  {args.start_year}  |  {args.preset}")
    for path in [p for p in written if p is not None]:
        print(f"  wrote {Path(path).name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
