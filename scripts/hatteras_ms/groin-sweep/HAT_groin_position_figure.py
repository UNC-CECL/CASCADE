#!/usr/bin/env python3
"""
Observed against modelled shoreline position across the groin, per period.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_position_figure.py
    python scripts/hatteras_ms/groin-sweep/HAT_groin_position_figure.py --preset zeroBE

Where the shoreline started, where it ended, and the fit through CoastSat. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
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

from site_layer.hat_figure_style import (apply_style, C, INK, INK_MUTED,  # noqa: E402
                              caption, figsize, open_frame, save, _title)
from HAT_groin_sweep_config import (  # noqa: E402
    GROIN_SWEEP_ROOT,
    END_YEAR,
    FIT_DOMAINS_GIS,
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
    PERIODS,
    PRESETS,
    combo_dir_name,
    sweep_output_dir,
)
from HAT_groin_sweep_figures import (  # noqa: E402
    GROIN_COLOR,
    OBSERVED_COLOR,
    MODEL_COLOR,
    RANK_METRIC,
    _cell_label,
    _footnote,
    load_scored,
    profile_be1,
    tied_best,
)

# --- CONFIG ------------------------------------------------------------------
CHAINAGE_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
                / "1-observations" / "shoreline_rates_by_era" / "output"
                / "groin_analysis_chainage_all.csv")

# The groin field's real footprint, from HAT_groin_shoreline_analysis_v2.py's metadata
GROIN_FIELD_NORTHING = (3901373.14, 3901788.79)

ZOOM_DOMAINS = (4, 8)
NOGROIN_COLOR = C["BASE"]   # the run without the modification under test
MIN_OBS_PER_DOMAIN = 20
# -----------------------------------------------------------------------------


# matplotlib on the Agg backend, imported when needed
def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


# Observations

# The shoreline-position observations, CoastSat only
def load_chainage():
    if not CHAINAGE_CSV.exists():
        raise FileNotFoundError(
            f"shoreline chainage not found at {CHAINAGE_CSV}. It is produced "
            f"by HAT_groin_shoreline_analysis_v2.py in "
            f"hard-structures/groin/1-observations/shoreline_rates_by_era/.")
    frame = pd.read_csv(
        CHAINAGE_CSV,
        usecols=["domain", "decimal_year", "chainage_m", "alongshore_m",
                 "source", "transect_id"])
    # CoastSat only: other sources measure a different feature
    return frame[frame["source"] == "coastsat"].copy()


# Start and end shoreline position per domain, by OLS over the period
def fitted_positions(chainage, period, domains):
    lo, hi = float(period), float(END_YEAR[period])
    window = chainage[(chainage["decimal_year"] >= lo)
                      & (chainage["decimal_year"] <= hi)]
    start, end, counts = {}, {}, {}
    for domain in domains:
        rows = window[window["domain"] == domain]
        if len(rows) < MIN_OBS_PER_DOMAIN:
            continue
        slope, intercept = np.polyfit(rows["decimal_year"],
                                      rows["chainage_m"], 1)
        start[domain] = float(intercept + slope * lo)
        end[domain] = float(intercept + slope * hi)
        counts[domain] = len(rows)
    return start, end, counts


# Per-transect start and end positions, for the zoom panel
def transect_positions(chainage, period, domains):
    lo, hi = float(period), float(END_YEAR[period])
    window = chainage[(chainage["decimal_year"] >= lo)
                      & (chainage["decimal_year"] <= hi)
                      & chainage["domain"].between(*domains)]
    rows = []
    for transect, grp in window.groupby("transect_id"):
        if len(grp) < MIN_OBS_PER_DOMAIN:
            continue
        slope, intercept = np.polyfit(grp["decimal_year"], grp["chainage_m"], 1)
        rows.append(dict(alongshore_m=float(grp["alongshore_m"].iloc[0]),
                         domain=float(grp["domain"].iloc[0]),
                         start_m=float(intercept + slope * lo),
                         end_m=float(intercept + slope * hi)))
    return pd.DataFrame(rows).sort_values("alongshore_m")


# Alongshore extent of each domain, measured from the observations
def domain_alongshore_bounds(chainage, domains):
    window = chainage[chainage["domain"].between(*domains)]
    grouped = window.groupby("domain")["alongshore_m"].agg(["min", "max"])
    return {int(d): (float(r["min"]), float(r["max"]))
            for d, r in grouped.iterrows()}


# Model

# Modelled start->end shoreline change per domain, seaward-positive
def model_change(period, preset, combo, domains):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as geometry
    path = sweep_output_dir(period, preset) / combo / "shoreline_matrix.npy"
    if not path.exists():
        return None
    matrix = np.load(path)
    # x_s increases LANDWARD; negate so + is seaward, matching chainage.
    change = -(matrix[-1] - matrix[0])
    return {d: float(change[geometry.gis_to_pad(d)]) for d in domains}


# The fitted cell for one sweep, or None if it has not been scored
def best_cell(period, preset):
    frame = load_scored(period, preset)
    if frame is None or frame.empty:
        return None
    surface = profile_be1(frame)
    groin = surface[surface["M"] > 0]
    if groin.empty:
        return None
    best, tied = tied_best(groin)
    return best


# The M = 0 combination paired with `best`, at the same be1
def baseline_combo(period, preset, best):
    be1 = None if pd.isna(best.get("be1")) else float(best["be1"])
    return combo_dir_name(0.0, be1, 0.0)


# Figure

# Draws the two-period position figure for one preset
def draw(preset, chainage):
    plt = _matplotlib()

    domains = list(FIT_DOMAINS_GIS)
    apply_style()
    figure, axes = plt.subplots(
        len(PERIODS), 2, figsize=figsize("double", aspect=0.66),
        gridspec_kw=dict(width_ratios=[1.5, 1]), constrained_layout=True)
    notes = []

    drew_any = False
    for row, period in enumerate(PERIODS):
        reach_axis, zoom_axis = axes[row]
        best = best_cell(period, preset)
        if best is None:
            for axis in (reach_axis, zoom_axis):
                axis.text(0.5, 0.5,
                          f"{period} to {END_YEAR[period]}, {preset}\n"
                          "not swept yet",
                          ha="center", va="center", fontsize=8,
                          color=INK_MUTED, transform=axis.transAxes)
                axis.set_xticks([]); axis.set_yticks([])
            continue
        drew_any = True

        start, end, counts = fitted_positions(chainage, period, domains)
        have = sorted(start)
        x = np.array(have, dtype=float)
        obs_start = np.array([start[d] for d in have])
        obs_end = np.array([end[d] for d in have])

        fitted = model_change(period, preset, best["combo"], have)
        nogroin = model_change(period, preset,
                               baseline_combo(period, preset, best), have)

        # Reach panel
        reach_axis.plot(x, obs_start, marker="o", markersize=2.8,
                        linestyle="--", color=INK_MUTED, linewidth=1.2,
                        label=f"observed at {period}, the start", zorder=3)
        reach_axis.plot(x, obs_end, marker="o", markersize=2.8,
                        color=OBSERVED_COLOR, linewidth=1.8,
                        label=f"observed at {END_YEAR[period]}, the end",
                        zorder=5)
        if nogroin is not None:
            reach_axis.plot(x, obs_start + np.array([nogroin[d] for d in have]),
                            marker="^", markersize=2.8, color=NOGROIN_COLOR,
                            linewidth=1.4, linestyle=":",
                            label="modelled end, groin off", zorder=4)
        if fitted is not None:
            reach_axis.plot(x, obs_start + np.array([fitted[d] for d in have]),
                            marker="s", markersize=2.8, color=MODEL_COLOR,
                            linewidth=1.6,
                            label=f"modelled end, {_cell_label(best)}",
                            zorder=6)

        # DOMAIN-COORDINATE shading belongs to the reach panel ONLY
        reach_axis.axvspan(GROIN_UPDRIFT_GIS - 0.5, GROIN_UPDRIFT_GIS + 0.5,
                           color="0.90", zorder=0)
        reach_axis.axvspan(GROIN_DOWNDRIFT_GIS - 0.5,
                           GROIN_DOWNDRIFT_GIS + 0.5,
                           color="0.94", zorder=0)
        for axis in (reach_axis, zoom_axis):
            axis.grid(axis="y")
            axis.set_axisbelow(True)
            open_frame(axis)

        misfit = (np.nan if fitted is None else
                  float(np.sqrt(np.mean(
                      (obs_start + np.array([fitted[d] for d in have])
                       - obs_end) ** 2))))
        _title(reach_axis, row * 2, f"{period} to {END_YEAR[period]}, the reach")
        reach_axis.set_xlabel("GIS domain (south → north)")
        reach_axis.set_ylabel("shoreline position (m from datum)\n"
                              "positive is seaward")
        reach_axis.set_xticks(domains)
        reach_axis.legend(loc="best", fontsize=7)
        notes.append("({}) {} to {}: end-position RMSE {:.1f} m.".format(
            chr(ord("a") + row * 2), period, END_YEAR[period], misfit))

        # Zoom panel: transect resolution

        # X IS REAL ALONGSHORE DISTANCE, NOT THE DOMAIN ID
        tr = transect_positions(chainage, period, ZOOM_DOMAINS)
        bounds = domain_alongshore_bounds(chainage, ZOOM_DOMAINS)
        if not tr.empty:
            zoom_axis.plot(tr["alongshore_m"], tr["start_m"], linestyle="--",
                           color=INK_MUTED, linewidth=1.0,
                           label=f"observed at {period}", zorder=3)
            zoom_axis.plot(tr["alongshore_m"], tr["end_m"],
                           color=OBSERVED_COLOR, linewidth=1.4,
                           label=f"observed at {END_YEAR[period]}", zorder=5)
        zoom_have = [d for d in have
                     if ZOOM_DOMAINS[0] <= d <= ZOOM_DOMAINS[1] and d in bounds]
        if fitted is not None and zoom_have:
            for index, domain in enumerate(zoom_have):
                lo_m, hi_m = bounds[domain]
                zoom_axis.hlines(start[domain] + fitted[domain], lo_m, hi_m,
                                 color=MODEL_COLOR, linewidth=1.8, zorder=6,
                                 label="modelled end, one 500 m cell each"
                                 if index == 0 else None)
            for domain in zoom_have:
                zoom_axis.axvline(bounds[domain][0], color="0.90",
                                  linewidth=0.6, zorder=0)
        for domain, shade in ((GROIN_UPDRIFT_GIS, "0.90"),
                              (GROIN_DOWNDRIFT_GIS, "0.94")):
            if domain in bounds:
                zoom_axis.axvspan(*bounds[domain], color=shade, zorder=0)
        _title(zoom_axis, row * 2 + 1,
               f"the structure, D{ZOOM_DOMAINS[0]} to D{ZOOM_DOMAINS[1]}")
        if bounds:
            lo_m = min(v[0] for v in bounds.values())
            hi_m = max(v[1] for v in bounds.values())
            zoom_axis.set_xlim(lo_m - 60, hi_m + 60)
        zoom_axis.set_xlabel("alongshore distance (m, south → north)")
        zoom_axis.legend(loc="best", fontsize=7)

    _footnote(
        figure,
        "Observed shoreline POSITION at the start and end of each hindcast "
        "period against the model at the fitted groin, for the {} preset. The "
        "left panels are the whole D1\u2013D12 fit window; the right panels "
        "zoom to the structure at transect resolution, where the observations "
        "come every 62 m or so and the model has one value per 500 m cell. "
        "{} "
        "Model and observed START at the same position BY CONSTRUCTION: the "
        .format(preset, " ".join(notes)) +
        "".join([]) +
        "model's cross-shore origin is Barrier3D's own, so it is plotted as "
        "observed start + model change (as build_shoreline_target does). Only "
        "the separation at the END year is informative. Positions are OLS fits "
        "to CoastSat chainage over each period, evaluated at the endpoints -- "
        "the same estimator the sweep is scored on. The groin field (4 "
        "structures) occupies D6; its ~190 m fillet is narrower than one 500 m "
        "model cell, which is what the right-hand panels show.", width=185)

    if not drew_any:
        plt.close(figure)
        return None
    figure_dir = GROIN_SWEEP_ROOT / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    return save(figure, figure_dir / f"position_{preset}.png", close=True)[0]


# Run: one figure per period and preset
def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--preset", choices=PRESETS, action="append",
                        help="repeatable; default every preset")
    args = parser.parse_args()

    chainage = load_chainage()
    print("=" * 72)
    print("SHORELINE POSITION FIGURE")
    print("=" * 72)
    print(f"  {len(chainage):,} CoastSat chainage observations loaded")

    wrote = []
    for preset in (args.preset or list(PRESETS)):
        path = draw(preset, chainage)
        if path is None:
            print(f"  {preset:<8} no scored sweep yet -- skipped")
            continue
        print(f"  {preset:<8} -> {path}")
        wrote.append(path)
    return 0 if wrote else 1


if __name__ == "__main__":
    sys.exit(main())
