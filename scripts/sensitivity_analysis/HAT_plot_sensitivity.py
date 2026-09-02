#!/usr/bin/env python3
"""Figures for the hindcast parameter sensitivity sweep.

WHY THIS REPLACED plot_sensitivity_vs_coastsat.py
    That script read a "timestamped session folder" produced by the two
    hand-built sensitivity scripts that have since been deleted, and it needed
    SESSION_DIR and two CoastSat CSV paths retyped at the top of the file before
    every use. Neither the folder layout nor those scripts exist any more.

    This reads the run registry instead. `HAT_hindcast_sensitivity.py` writes a
    manifest line per cell recording the parameter, its value and the directory
    the run landed in; `run_index.csv` carries the skill metrics; each run
    directory carries its own per-domain rates. Nothing is retyped and nothing
    is passed between the two scripts except the manifest.

THE OBSERVED LAYER IS NOT REDRAWN HERE
    The CoastSat scatter and the LOESS curves come from
    `cascade_pipeline.plotting.rate_comparison.plot_coastsat_overlay`, the same
    call every other figure in the repo makes, and the target table comes from
    `build_target_table`. A second implementation here would be free to drift on
    styling, on the southern splice, and on which window is the reference.

WHAT THE MODEL IS SCORED AGAINST, AND WHY THE FIGURE SHOWS THREE THINGS
    The target is a HYBRID, not one curve: GIS 1..skip_southern_domains are raw
    per-domain means (LOESS is suppressed near Oregon Inlet, where boundary
    effects dominate and smoothing would hide the gradient), and the rest is the
    LOESS 10-domain reference. So the alongshore panels carry the raw transect
    scatter, both LOESS windows, and the spliced target -- reporting RMSE
    against a hybrid while plotting only one of its halves would misstate what
    the number means.

INTERIOR METRICS, NOT ISLAND-WIDE
    Every skill number here is the interior one (GIS 2-89). The two end domains
    carry imposed background-erosion rates under edgeBE and calibBE alike, so
    island-wide skill partly scores the boundary condition rather than the model.

Usage:
    python HAT_plot_sensitivity.py --start-year 1984 --preset edgeBE
    python HAT_plot_sensitivity.py --start-year 2004 --preset edgeBE
    python HAT_plot_sensitivity.py --start-year 1984 --circularity

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.patheffects as mpatheffects
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_BASE_DIR = _HERE.parents[2]
if not (PROJECT_BASE_DIR / "pyproject.toml").exists():
    raise RuntimeError(
        f"CASCADE repo root not found: {PROJECT_BASE_DIR} has no pyproject.toml.")
for _path in (PROJECT_BASE_DIR / "scripts",
              PROJECT_BASE_DIR / "scripts" / "hatteras_ms",
              PROJECT_BASE_DIR / "scripts" / "sensitivity_analysis"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_PERIODS  # noqa: E402
from cascade_pipeline.coastsat_loess import (  # noqa: E402
    CoastSatDataset, LoessConfig, build_coastsat_series)
from cascade_pipeline.hindcast import build_target_table  # noqa: E402
from cascade_pipeline.run_registry import (  # noqa: E402
    CALIBRATION_ARM, find_run_dir)
from cascade_pipeline.plotting.rate_comparison import (  # noqa: E402
    DEFAULT_RATE_COMPARISON, plot_coastsat_overlay)
from HAT_hindcast_sensitivity import SWEEPS, normalise  # noqa: E402
from HAT_hindcast_config import field_default  # noqa: E402

# Reading order, most informative first. Alphabetical put the relocation axis
# -- the one that barely moves shoreline skill -- in the leftmost panel and the
# first file, which is the opposite of how these should be read. Used for both
# the panel columns and the numeric filename prefixes, so the figure and a
# directory listing agree.
AXIS_ORDER = ("wave_height", "wave_period", "wave_asymmetry",
              "wave_angle_high_fraction", "relocation_setback")

OUT_ROOT = PROJECT_BASE_DIR / "output" / "sensitivity_analysis"
RAW_RUNS = PROJECT_BASE_DIR / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
COASTSAT_BASE_DIR = (PROJECT_BASE_DIR / "scripts" / "input_prep"
                     / "5-scr" / "CoastSat")

# Must match section 8.1 of HAT_hindcast_1984_2024.py. These are the LoessConfig
# defaults, so the two agree by construction rather than by copying -- but the
# run metadata records the target it actually used, and check_target_matches()
# below asserts against that rather than trusting this line.
LOESS_CONFIG = LoessConfig()
TARGET_WINDOW = max(LOESS_CONFIG.window_domains)

# The model curves are a SEQUENTIAL ramp: the swept values are ordered
# magnitudes, not categories, so one hue light-to-dark is the correct encoding
# and a categorical cycle would imply the values are unrelated. Warm on purpose
# -- the observed layer owns the blues (#5BA3C9 scatter, #6BAED6 / #08519C
# LOESS), so a blue model ramp would collide with the thing it is measured
# against.
MODEL_RAMP = plt.get_cmap("YlOrRd")
RAMP_LO, RAMP_HI = 0.30, 0.92          # skip the near-white and near-black ends
BASELINE_COLOR = "#111111"
# The current setting, on the alongshore panels only. It needs a hue that is in
# NEITHER family already on those axes -- the warm model ramp or the blue
# observed layer -- because it is the one curve a reader goes looking for.
# Black was too close to the navy target and, worse, the legend never said
# which value it was: the calibration cell is skipped by the sweep (it IS the
# baseline), so it never appears on the colourbar either.
CURRENT_COLOR = "#D01C8B"
INK, INK_MUTED = "#222222", "#666666"


# =============================================================================
# READING THE SWEEP
# =============================================================================

def load_cells(start_year, preset):
    """Every completed cell for one period and preset, newest wins.

    The manifest is append-only, so re-running a cell adds a second line rather
    than replacing the first. Keeping the LAST line per (sweep, value) is what
    makes a re-run authoritative -- and it has to be, because a cell re-run
    after the resize_interior_domain fix supersedes the one that failed.

    Args:
        start_year: Period start year.
        preset: Source/sink preset the cells ran under.

    Returns:
        DataFrame with sweep, setting, value, sort_key, run_dir, run_name.
    """
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
        run_dir = Path(row["detail"])
        rows.append(dict(sweep=row["sweep"], setting=row["setting"],
                         value=row["value"], run_dir=run_dir,
                         run_name=run_dir.name))
    if not rows:
        return pd.DataFrame(columns=["sweep", "setting", "value", "sort_key",
                                     "run_dir", "run_name"])

    frame = pd.DataFrame(rows)
    frame["norm"] = frame["value"].map(normalise)
    # `measured` (None) is not a point on the numeric axis; it sorts last and is
    # drawn as its own category rather than being given a fake number.
    frame["sort_key"] = frame["norm"].map(
        lambda v: np.inf if v is None else v)
    frame = frame.drop_duplicates(subset=["sweep", "norm"], keep="last")
    return frame.sort_values(["sweep", "sort_key"]).reset_index(drop=True)


def baseline_name(run_name):
    """The matrix run a cell is a departure from.

    Every cell name is its baseline plus one trailing token, so the baseline is
    the name with that token removed. Derived rather than rebuilt from switches
    because the token is the only difference by construction -- see the run-name
    tokens in cascade_pipeline.hindcast.

    Args:
        run_name: A sensitivity cell's directory name.

    Returns:
        The baseline run name.

    Raises:
        ValueError: If the name carries no sensitivity token, which would mean
            the caller handed over a matrix run and the "baseline" returned
            would silently be a different scenario.
    """
    stem, _, token = run_name.rpartition("_")
    if not (token.startswith("wave") or token.startswith("rset")):
        raise ValueError(
            f"{run_name} has no sensitivity token; its last token is {token!r}")
    return stem


def load_index():
    """run_index.csv for the calibration arm, indexed by run name.

    FILTERED TO ONE ARM, NOT DEDUPLICATED. A run name describes the SCENARIO
    and an arm describes the FORCING, so one name legitimately appears once per
    arm -- and the index is keyed on (run_name, Hs_m, arm) for that reason.
    This function previously collapsed those with
    `drop_duplicates(keep="last")`, which does not pick the calibration run: it
    picks whichever row sorts last. On 2026-09-01 that made
    `HAT_1984_2004_edgeBE_road_bdm_groin` resolve to the `waveHs3_probe` arm --
    a Newton probe at Hs = 3.0 -- and that row is the BASELINE every edgeBE
    wave cell is drawn against, so each panel measured its sweep against a
    different experiment. Silently: the row exists and reads cleanly.

    Every sweep cell is a calibration-arm run, so the baselines are too. Any
    other arm is a different experiment and is dropped rather than ranked.

    Returns:
        DataFrame indexed by run_name, one row per name.

    Raises:
        ValueError: If a name is still duplicated within the calibration arm,
            which the index key makes impossible and so means a fault in the
            index rather than something to pick a winner from.
    """
    frame = pd.read_csv(RUN_INDEX)
    if "arm" in frame.columns:
        frame = frame[frame["arm"].fillna(CALIBRATION_ARM) == CALIBRATION_ARM]
    repeated = sorted(frame["run_name"][frame["run_name"].duplicated()].unique())
    if repeated:
        raise ValueError(
            f"run_index.csv has more than one {CALIBRATION_ARM}-arm row for "
            f"{repeated}. The (run_name, Hs_m, arm) key should make that "
            f"impossible; the index needs repairing, not deduplicating.")
    return frame.set_index("run_name")


def model_rates(run_dir, run_name):
    """Per-domain modelled LRR for one run.

    Args:
        run_dir: The run's directory.
        run_name: The run's name, which prefixes its files.

    Returns:
        DataFrame with gis_domain and lrr_m_yr.
    """
    path = run_dir / f"{run_name}_shoreline_change_rate.csv"
    frame = pd.read_csv(path)
    return frame[["gis_domain", "lrr_m_yr"]]


def check_target_matches(index, run_names):
    """Assert every run was scored against the target this figure draws.

    The skill numbers plotted here come from run_index.csv, and the curve comes
    from build_target_table with TARGET_WINDOW. If a run was scored against a
    different window, the panel would show a curve that is not the one its RMSE
    refers to. Each run records its target in run_metadata.json, so this is
    checkable rather than assumed.
    """
    expected = f"CoastSat LOESS {TARGET_WINDOW}-domain"
    for name in run_names:
        row = index.loc[name]
        # RESOLVED, NOT JOINED BY HAND, AND IN THE ARM THE ROW CLAIMS. This
        # join had no slot for the arm component, so a run filed under one
        # resolved to a path that does not exist -- and the `continue` below
        # then skipped the check in silence, which is the opposite of what a
        # check is for. find_run_dir raises instead, naming the arms the run
        # IS under. Passing row.arm rather than letting it default is what
        # makes this figure read the run it is scoring: the index carries the
        # arm for exactly this, and as of the 2026-09-01 backfill every row's
        # arm is the directory it is actually in.
        directory = find_run_dir(
            RAW_RUNS, name, (int(row.start_year), int(row.end_year)),
            row.source_sink_preset, row.arm)
        meta = directory / f"{name}_run_metadata.json"
        # A run predating the metadata field records no target and cannot be
        # checked. That is the only thing this skip is allowed to mean now.
        if not meta.exists():
            continue
        target = json.loads(meta.read_text(encoding="utf-8")).get(
            "skill", {}).get("target")
        if target and target != expected:
            raise ValueError(
                f"{name} was scored against {target!r} but these figures draw "
                f"{expected!r}; the plotted curve would not be the one its "
                f"RMSE refers to")


# =============================================================================
# THE OBSERVED LAYER
# =============================================================================

def coastsat_layers(start_year):
    """The CoastSat series and the spliced target, built the hindcast's way.

    Args:
        start_year: Period start year, which selects the active dataset.

    Returns:
        (cs_series, target_frame).
    """
    datasets = [
        CoastSatDataset(
            label="CoastSat LRR (1984-2004)", period_start=1984,
            csv_path=str(COASTSAT_BASE_DIR / "1984_2004"
                         / "transect_lrr_full.csv")),
        CoastSatDataset(
            label="CoastSat LRR (2004-2024)", period_start=2004,
            csv_path=str(COASTSAT_BASE_DIR / "2004_2024"
                         / "transect_lrr_full.csv")),
    ]
    series = build_coastsat_series(
        datasets, active_period_start=start_year, loess_config=LOESS_CONFIG,
        domains=HATTERAS_DOMAINS)
    active = next((cs for cs in series if cs["active"]), None)
    if active is None:
        raise RuntimeError(f"no CoastSat dataset starts at {start_year}")
    target = build_target_table(
        active, LOESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)
    return series, target


# =============================================================================
# HOUSE STYLE
# =============================================================================
# One place, applied once at import. Set here rather than per-figure so every
# panel in the directory carries the same type sizes and rules -- a set of
# figures read side by side is the unit that has to look consistent, not any
# one of them.

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
    # Two spines, not four. A box around the data adds a rule the reader has to
    # look past on every panel and encodes nothing.
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


def panel_label(ax, letter, dx=-0.02):
    """Puts a lower-case panel letter outside the axes, top left.

    Outside, so it can never land on a data point. Journals want the letter on
    every panel of a multi-panel figure and want it in a fixed place -- reading
    a ten-panel grid means finding (f) without hunting for it.
    """
    ax.text(dx, 1.04, f"({letter})", transform=ax.transAxes,
            fontsize=9.5, fontweight="bold", va="bottom", ha="right",
            color="#222222")


def tidy(ax, minor=True):
    """The per-axes cleanup every panel gets."""
    if minor:
        ax.minorticks_on()
        ax.tick_params(which="minor", length=2)
        ax.grid(which="minor", visible=False)
    ax.tick_params(which="major", length=3.5)


# =============================================================================
# FIGURES
# =============================================================================

def value_label(value):
    """A swept value as it appears in a legend or colourbar tick."""
    norm = normalise(value)
    return "measured" if norm is None else f"{norm:g}"


def ramp_colors(count):
    """`count` colours from the sequential model ramp, light to dark."""
    if count == 1:
        return [MODEL_RAMP(RAMP_HI)]
    return [MODEL_RAMP(RAMP_LO + (RAMP_HI - RAMP_LO) * i / (count - 1))
            for i in range(count)]


def axis_xlabel(sweep):
    """The x-axis label for one swept axis, with its unit."""
    units = SWEEPS[sweep]["units"]
    label = SWEEPS[sweep]["label"]
    return f"{label} ({units})" if units else label


def plot_skill_overview(cells, index, start_year, preset, out_dir):
    """Interior RMSE and bias against the swept value, one column per axis.

    Two ROWS rather than two y-scales on one axes: RMSE and bias are different
    measures on different scales, and overlaying them on twin axes is the one
    chart form that reliably misleads, because the crossing point of the two
    curves is an artifact of the scaling choice.
    """
    sweeps = [s for s in AXIS_ORDER if s in set(cells["sweep"])]
    if not sweeps:
        return None
    fig, axes = plt.subplots(
        2, len(sweeps), figsize=(2.55 * len(sweeps) + 0.9, 5.6), squeeze=False,
        sharex="col", constrained_layout=True)

    footnotes = []
    for col, sweep in enumerate(sweeps):
        block = cells[cells.sweep == sweep]
        base = baseline_name(block.iloc[0].run_name)
        base_row = index.loc[base]

        numeric = block[np.isfinite(block.sort_key)]
        # The calibration value is a POINT ON THE CURVE, not just a reference
        # line. Without it the series jumps straight across its own reference
        # -- the Hs cells run 2.0 then 3.0 and the segment between them hides
        # the 2.5 the whole sweep is measured against.
        default = normalise(field_default(SWEEPS[sweep]["setting"]))
        x = list(numeric.sort_key.to_numpy(dtype=float))
        rmse = [index.loc[n].rmse_interior_m_yr for n in numeric.run_name]
        bias = [index.loc[n].mean_bias_interior_m_yr for n in numeric.run_name]
        if default is not None:
            x.append(default)
            rmse.append(base_row.rmse_interior_m_yr)
            bias.append(base_row.mean_bias_interior_m_yr)
        order = np.argsort(x)
        x = np.asarray(x)[order]
        rmse = np.asarray(rmse)[order]
        bias = np.asarray(bias)[order]

        for row, (values, base_value, label) in enumerate((
                (rmse, base_row.rmse_interior_m_yr, "Interior RMSE (m yr$^{-1}$)"),
                (bias, base_row.mean_bias_interior_m_yr,
                 "Interior mean bias (m yr$^{-1}$)"))):
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

            # An axis the model barely responds to autoscales to its own noise
            # and reads as a cliff. Floor the span on the run's own RMSE so a
            # flat response LOOKS flat, at whatever scale that run lives at.
            floor_span = 0.12 * abs(base_row.rmse_interior_m_yr)
            lo, hi = ax.get_ylim()
            if hi - lo < floor_span:
                mid = 0.5 * (lo + hi)
                ax.set_ylim(mid - floor_span / 2, mid + floor_span / 2)

            # Without this the calibration diamond is half-clipped by the
            # spine on any axis whose calibration value is an endpoint of the
            # swept range -- which is three of the five.
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

        # `measured` has no numeric x. It is reported in the caption strip
        # rather than dropped, and rather than floated inside the axes where it
        # can land on the curve.
        for name in block[~np.isfinite(block.sort_key)].run_name:
            footnotes.append(
                f"{SWEEPS[sweep]['label']}, measured: "
                f"RMSE {index.loc[name].rmse_interior_m_yr:.3f} "
                f"m yr$^{{-1}}$")

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3,
               bbox_to_anchor=(0.5, -0.035), frameon=False)

    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    caption = (f"Parameter sensitivity, {start_year}–{end_year}, {preset}. "
               f"Interior domains (GIS 2–89).")
    if footnotes:
        caption += "  " + "; ".join(footnotes) + "."
    fig.suptitle(caption, y=1.035, fontsize=10)

    path = out_dir / "01_skill_overview.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_alongshore(cells, sweep, start_year, preset, cs_series, target,
                    out_dir):
    """Modelled LRR per domain for every cell of one axis, over the observed layer.

    THE LEGEND IS NOT IN THE DATA AREA. A 16-entry legend box placed inside
    these axes sat directly on the target curve and the southern transect
    scatter, which is the densest and most important corner of the panel. The
    swept values become a discrete colourbar on the right -- correct anyway,
    because they are an ordered ramp rather than unrelated categories -- and
    the four reference curves get a horizontal legend beneath the axes.
    """
    block = cells[cells.sweep == sweep]
    if block.empty:
        return None
    base = baseline_name(block.iloc[0].run_name)
    colors = ramp_colors(len(block))

    fig = plt.figure(figsize=(9.6, 4.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 2, width_ratios=[60, 1], wspace=0.02)
    ax = fig.add_subplot(grid[0, 0])
    cax = fig.add_subplot(grid[0, 1])

    plot_coastsat_overlay(
        ax, cs_series, LOESS_CONFIG, DEFAULT_RATE_COMPARISON,
        x_transform=lambda along_m: (along_m / HATTERAS_DOMAINS.domain_spacing_m
                                     + HATTERAS_DOMAINS.first_gis_id),
    )
    # Dashed on purpose: over D11-90 this IS the LOESS-10 curve already drawn
    # underneath, so a second solid line would just thicken it. Dashed, the
    # reader sees them coincide there and separate over D1-10, which is where
    # the target stops being a LOESS curve at all.
    ax.plot(target.gis_domain, target.target_lrr_m_yr, color="#08306B", lw=1.8,
            ls=(0, (5, 2)), zorder=5,
            label=f"Target (raw D1–{LOESS_CONFIG.skip_southern_domains}, "
                  f"LOESS {TARGET_WINDOW} elsewhere)")

    base_dir = block.iloc[0].run_dir.parent / base
    base_rates = model_rates(base_dir, base)
    # Named with its VALUE, not just "calibration run". On the Hs panel this is
    # the line the reader is looking for -- where the current setting sits among
    # the alternatives -- and "calibration run" does not answer that.
    current = normalise(field_default(SWEEPS[sweep]["setting"]))
    units = SWEEPS[sweep]["units"]
    current_text = "measured" if current is None else f"{current:g}"
    if units:
        current_text += f" {units.split()[0]}"
    # Identity is carried by colour AND dash pattern AND weight, so the line
    # survives a greyscale print and red-green colour blindness alike.
    ax.plot(base_rates.gis_domain, base_rates.lrr_m_yr, color=CURRENT_COLOR,
            lw=2.6, ls=(0, (7, 1.6)), zorder=8,
            path_effects=[mpatheffects.withStroke(linewidth=4.6,
                                                  foreground="white")],
            label=f"Model at the CURRENT setting ({current_text})")

    for color, (_, cell) in zip(colors, block.iterrows()):
        rates = model_rates(cell.run_dir, cell.run_name)
        ax.plot(rates.gis_domain, rates.lrr_m_yr, color=color, lw=1.0,
                alpha=0.95, zorder=6)

    ax.axhline(0.0, color=INK_MUTED, lw=0.7, ls="--", zorder=1)
    ax.set_xlim(HATTERAS_DOMAINS.first_gis_id - 0.5,
                HATTERAS_DOMAINS.last_gis_id + 0.5)
    ax.set_xlabel("Alongshore position (GIS domain, south → north)")
    ax.set_ylabel("Shoreline change rate, LRR (m yr$^{-1}$)")
    tidy(ax)

    # Reference curves only -- the swept values are on the colourbar. Below the
    # axes, so it cannot cover data at any y-limit.
    handles, labels = ax.get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: labels[i].startswith("Co"))
    ax.legend([handles[i] for i in order], [labels[i] for i in order],
              loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              frameon=False, handlelength=2.4, columnspacing=2.0)

    # Discrete colourbar: one swatch per cell, because the swept values are not
    # evenly spaced and the colours are assigned by rank. A continuous bar would
    # imply a linear mapping from value to colour that is not what was drawn.
    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(len(colors) + 1) - 0.5
    bar = mpl.colorbar.ColorbarBase(
        cax, cmap=cmap, norm=mcolors.BoundaryNorm(bounds, cmap.N),
        ticks=np.arange(len(colors)), spacing="uniform")
    bar.set_ticklabels([value_label(v) for v in block.value])
    bar.ax.tick_params(length=0, labelsize=7.5)
    bar.outline.set_linewidth(0.6)
    bar.outline.set_edgecolor("#333333")
    units = SWEEPS[sweep]["units"]
    bar.set_label(f"{SWEEPS[sweep]['label']}" + (f" ({units})" if units else ""),
                  fontsize=8.5, labelpad=8)

    end_year = HATTERAS_PERIODS[start_year].get("end_year", start_year + 20)
    ax.set_title(f"{SWEEPS[sweep]['label']} sensitivity, "
                 f"{start_year}–{end_year}, {preset}", pad=8)

    path = out_dir / f"{AXIS_ORDER.index(sweep) + 2:02d}_alongshore_{sweep}.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def plot_relocation_outcomes(cells, index, start_year, preset, out_dir):
    """Road outcomes against the relocation target.

    A bar chart, not a line: the question is a count at each candidate value,
    and the interesting feature is a THRESHOLD between two adjacent values, not
    a slope. `measured` sits beside the numeric bars as its own category because
    it is a different policy, not a larger number.
    """
    block = cells[cells.sweep == "relocation_setback"]
    if block.empty:
        return None
    base = baseline_name(block.iloc[0].run_name)
    base_row = index.loc[base]
    default = normalise(field_default("relocation_setback_m"))

    labels = [value_label(v) for v in block.value] + [f"{default:g}"]
    drowned = [index.loc[n].roads_drowned for n in block.run_name] + \
              [base_row.roads_drowned]
    blocked = [index.loc[n].roads_reloc_blocked for n in block.run_name] + \
              [base_row.roads_reloc_blocked]
    keys = [np.inf if normalise(v) is None else normalise(v)
            for v in block.value] + [default]
    order = np.argsort(keys)
    labels = [labels[i] for i in order]
    drowned = [drowned[i] for i in order]
    blocked = [blocked[i] for i in order]
    is_base = [abs(keys[i] - default) < 1e-9 for i in order]

    # The blocked-relocation panel is usually all zeros -- real information,
    # but it does not need half the figure to say it.
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
        # Further out than the default: this figure's y-labels are long and
        # the short lower panel brings them closer to the letter.
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


def plot_circularity(start_year, index, out_dir):
    """The Hs skill curve under both presets, which is what makes the circularity visible.

    calibBE fits a per-domain background-erosion rate against the CoastSat LRR
    at the calibration Hs, so its RMSE minimum says as much about the fit as
    about the wave climate. edgeBE imposes a rate on the two end domains only,
    so its interior is generated by the physics. Two series, so both are named
    in a legend and neither is identified by colour alone.
    """
    curves = {}
    for preset in ("calibBE", "edgeBE"):
        cells = load_cells(start_year, preset)
        block = cells[cells.sweep == "wave_height"] if not cells.empty \
            else cells
        if len(block) < 2:
            continue
        curves[preset] = (
            block.sort_key.to_numpy(dtype=float),
            np.array([index.loc[n].rmse_interior_m_yr for n in block.run_name]),
            index.loc[baseline_name(block.iloc[0].run_name)],
        )
    if len(curves) < 2:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0),
                             constrained_layout=True)
    # Two entities, fixed colours: the preset owns the colour, so adding a third
    # preset later cannot repaint these two.
    preset_colors = {"calibBE": "#B24502", "edgeBE": "#08519C"}
    default_hs = float(field_default("hs"))

    # PANEL (b) EXISTS BECAUSE OF THE MAGNITUDE GAP. calibBE fits 48 domains to
    # the target and sits around 0.5-1.2; edgeBE fits two and sits around
    # 1.2-3.0. On one shared axis the calibBE basin is squashed into the bottom
    # third and the reader cannot compare the two SHAPES -- which is the entire
    # question. (b) re-plots each curve as a rise above its own minimum, so the
    # location and width of each basin can be read directly.
    for preset, (x, rmse, base_row) in curves.items():
        # The calibration Hs is a POINT ON THIS CURVE. Without it calibBE draws
        # a straight 2.0-to-3.0 segment straight over the top of its own
        # minimum at 2.5 -- hiding a feature the panel exists to show.
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
            (axes[0], "Interior RMSE (m yr$^{-1}$)", "Absolute skill"),
            (axes[1], "RMSE above each curve's own minimum (m yr$^{-1}$)",
             "Basin shape, magnitude removed")):
        ax.axvline(default_hs, color=INK_MUTED, lw=0.7, ls=":", zorder=1)
        ax.set_xlabel("Significant wave height, H$_s$ (m)")
        ax.set_ylabel(ylab)
        ax.set_title(ttl, pad=8)
        ax.margins(x=0.05)
        tidy(ax)
    axes[1].set_ylim(-0.02, 0.75)

    # Annotated at the TOP of panel (a), which is empty at the calibration Hs --
    # both curves are near their minima there, well below. At the foot it
    # grazed the calibBE diamond; moving the legend out from under the title
    # freed this corner up.
    axes[0].annotate(f"calibration H$_s$ = {default_hs:g} m",
                     xy=(default_hs, 1.0), xycoords=("data", "axes fraction"),
                     xytext=(5, -4), textcoords="offset points",
                     fontsize=8, color=INK_MUTED, ha="left", va="top")

    # Below the axes, where it cannot cover a curve at any y-limit. The diamond
    # is explained once in the caption instead of doubling every series entry.
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


def write_summary(cells, index, start_year, preset, out_dir):
    """One tidy CSV behind the figures, so a number can be read not measured off a line."""
    rows = []
    for _, cell in cells.iterrows():
        row = index.loc[cell.run_name]
        base = index.loc[baseline_name(cell.run_name)]
        rows.append(dict(
            start_year=start_year, preset=preset, axis=cell.sweep,
            setting=cell.setting, value=value_label(cell.value),
            run_name=cell.run_name,
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--start-year", type=int, default=1984,
                        choices=sorted(HATTERAS_PERIODS))
    parser.add_argument("--preset", default="edgeBE")
    parser.add_argument("--circularity", action="store_true",
                        help="also draw the calibBE-vs-edgeBE Hs panel")
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()

    # One directory per (period, preset). Flat output put 22 files with long
    # names in one listing, where the only way to find the pair you wanted to
    # compare was to read every name.
    end_year = HATTERAS_PERIODS[args.start_year].get(
        "end_year", args.start_year + 20)
    root = Path(args.out_dir) if args.out_dir else OUT_ROOT / "figures"
    out_dir = root / f"{args.start_year}_{end_year}_{args.preset}"
    out_dir.mkdir(parents=True, exist_ok=True)

    index = load_index()
    cells = load_cells(args.start_year, args.preset)
    if cells.empty:
        print(f"no completed cells for {args.start_year} / {args.preset}")
        return 1
    known = set(index.index)
    missing = sorted(set(cells.run_name) - known)
    if missing:
        raise ValueError(
            f"{len(missing)} cell(s) are in the manifest but not in "
            f"run_index.csv: {missing[:5]}. The index may have lost rows to a "
            f"concurrent write; re-run those cells before plotting.")
    check_target_matches(index, list(cells.run_name))

    cs_series, target = coastsat_layers(args.start_year)

    written = []
    written.append(plot_skill_overview(cells, index, args.start_year,
                                       args.preset, out_dir))
    for sweep in [s for s in AXIS_ORDER if s in set(cells.sweep)]:
        written.append(plot_alongshore(cells, sweep, args.start_year,
                                       args.preset, cs_series, target, out_dir))
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
