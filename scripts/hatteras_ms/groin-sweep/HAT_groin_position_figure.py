#!/usr/bin/env python3
"""Observed vs modelled SHORELINE POSITION across the groin, per period.

Every other figure in this directory plots a RATE. This one plots position:
where the shoreline started, where it ended, and where the model put it at the
fitted (M, f). A rate figure can hide a run that gets the trend right from the
wrong place; a position figure cannot.

WHERE THE OBSERVATIONS COME FROM
    `HAT-groin-gis-analysis/.../groin_analysis_chainage_all.csv` -- 904k
    shoreline observations carrying `chainage_m`, the cross-shore distance from
    the project's offshore datum line, already mapped to CASCADE domains.
    Chainage is SEAWARD-POSITIVE: verified on 2026-08-23 by differencing the
    1984 and 2004 domain means and correlating against the published CoastSat
    LRR (+0.74). The model's x_s is landward-positive, so it is negated.

    POSITIONS ARE FITTED, NOT BINNED. Two year-bins differenced give a noisy
    endpoint rate -- at D1 that route gave -1.6 m/yr against a published LRR of
    -4.2. Instead an OLS line is fitted to chainage against decimal year across
    the whole period and evaluated at both ends, so the plotted start and end
    are consistent with the LRR that the sweep is scored on. Same estimator,
    same answer.

WHY THE START LINES COINCIDE -- read this before reading the figure
    The model's cross-shore origin is Barrier3D's own, not a real datum, so the
    model cannot be placed on a surveyed axis independently. Following
    `build_shoreline_target` in cascade_pipeline.hindcast, the model is plotted
    as OBSERVED START + MODEL CHANGE. The consequence is unavoidable and worth
    stating plainly: model and observed start at the same place BY
    CONSTRUCTION. Only the separation at the END year carries information. The
    figure annotates this rather than letting a reader mistake a shared start
    for a validated initial condition.

THE ZOOM PANEL, AND WHY IT IS NOT A FAILURE
    The groin field -- four structures, northing 3901373-3901789 -- sits
    entirely inside D6, which spans 3901298-3901798. Its fillet is ~190 m
    wide while a model domain is 500 m, so the model CANNOT resolve it: its
    dipole is 500 m wide by construction. The right panel plots the
    observations at transect resolution against the model's 500 m steps so the
    mismatch is visible as a resolution limit rather than being averaged away.
    A reader who sees the model miss a 190 m notch should know the model was
    never able to draw one.

Usage:
    python HAT_groin_position_figure.py
    python HAT_groin_position_figure.py --preset zeroBE

Writes to output/groin_sweep/figures/:
    position_<preset>.png

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
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

CHAINAGE_CSV = (PROJECT_BASE_DIR / "hard-structures" / "groin"
                / "HAT-groin-gis-analysis" / "shoreline_output_coastsat"
                / "groin_analysis_chainage_all.csv")

# The groin field's real footprint, from HAT_groin_shoreline_analysis_v2.py's
# metadata. Drawn on the zoom panel so the structure's extent is visible
# rather than implied by a single line.
GROIN_FIELD_NORTHING = (3901373.14, 3901788.79)

ZOOM_DOMAINS = (4, 8)
NOGROIN_COLOR = "#777777"
MIN_OBS_PER_DOMAIN = 20


def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


# =============================================================================
# OBSERVATIONS
# =============================================================================

def load_chainage():
    """The shoreline-position observations, CoastSat only.

    Returns:
        A DataFrame with domain, decimal_year, chainage_m, alongshore_m.

    Raises:
        FileNotFoundError: If the GIS analysis output is absent.
    """
    if not CHAINAGE_CSV.exists():
        raise FileNotFoundError(
            f"shoreline chainage not found at {CHAINAGE_CSV}. It is produced "
            f"by HAT_groin_shoreline_analysis_v2.py in "
            f"hard-structures/groin/HAT-groin-gis-analysis/.")
    frame = pd.read_csv(
        CHAINAGE_CSV,
        usecols=["domain", "decimal_year", "chainage_m", "alongshore_m",
                 "source", "transect_id"])
    # CoastSat only: nc_state and wet_dry are different features measured to
    # different definitions, and mixing them into one fitted line would put
    # the source's offset into the trend.
    return frame[frame["source"] == "coastsat"].copy()


def fitted_positions(chainage, period, domains):
    """Start and end shoreline position per domain, by OLS over the period.

    The same estimator the sweep is scored on (see `compute_lrr`), so the
    endpoints of this line and the LRR are the same statement about the same
    data. Differencing two year-bins instead gives an endpoint rate that
    disagrees with the published LRR by a factor of two at the noisiest
    domains.

    Args:
        chainage: Observation frame from `load_chainage`.
        period: 1984 or 2004.
        domains: Iterable of GIS domain ids.

    Returns:
        (start, end, n_obs) dicts keyed by domain, in metres seaward-positive.
        Domains with too few observations are absent from all three.
    """
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


def transect_positions(chainage, period, domains):
    """Per-transect start and end positions, for the zoom panel.

    Same OLS treatment as `fitted_positions` but without the domain averaging,
    so a fillet narrower than a domain survives.

    Returns:
        A DataFrame with alongshore_m, start_m, end_m, one row per transect.
    """
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


def domain_alongshore_bounds(chainage, domains):
    """Alongshore extent of each domain, measured from the observations.

    Derived from the transects rather than from the domain polygons because
    the zoom panel plots observations on an alongshore axis and the model's
    cells must be drawn on the SAME axis. Taking the extent from the data that
    is being plotted keeps the two aligned even where a domain's transect
    coverage is partial.

    Args:
        chainage: Observation frame from `load_chainage`.
        domains: (first, last) inclusive domain range.

    Returns:
        {domain: (min_alongshore_m, max_alongshore_m)}.
    """
    window = chainage[chainage["domain"].between(*domains)]
    grouped = window.groupby("domain")["alongshore_m"].agg(["min", "max"])
    return {int(d): (float(r["min"]), float(r["max"]))
            for d, r in grouped.iterrows()}


# =============================================================================
# MODEL
# =============================================================================

def model_change(period, preset, combo, domains):
    """Modelled start->end shoreline change per domain, seaward-positive.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".
        combo: Combination directory name.
        domains: GIS domain ids wanted.

    Returns:
        A dict of domain -> change in metres, or None if the run is absent.
    """
    from hatteras_site_config import HATTERAS_DOMAINS as geometry
    path = sweep_output_dir(period, preset) / combo / "shoreline_matrix.npy"
    if not path.exists():
        return None
    matrix = np.load(path)
    # x_s increases LANDWARD; negate so + is seaward, matching chainage.
    change = -(matrix[-1] - matrix[0])
    return {d: float(change[geometry.gis_to_pad(d)]) for d in domains}


def best_cell(period, preset):
    """The fitted cell for one sweep, or None if it has not been scored."""
    frame = load_scored(period, preset)
    if frame is None or frame.empty:
        return None
    surface = profile_be1(frame)
    groin = surface[surface["M"] > 0]
    if groin.empty:
        return None
    best, tied = tied_best(groin)
    return best


def baseline_combo(period, preset, best):
    """The M = 0 combination paired with `best`, at the same be1."""
    be1 = None if pd.isna(best.get("be1")) else float(best["be1"])
    return combo_dir_name(0.0, be1, 0.0)


# =============================================================================
# FIGURE
# =============================================================================

def draw(preset, chainage):
    """Draws the two-period position figure for one preset."""
    plt = _matplotlib()

    domains = list(FIT_DOMAINS_GIS)
    figure, axes = plt.subplots(
        len(PERIODS), 2, figsize=(15, 9.5),
        gridspec_kw=dict(width_ratios=[1.5, 1]))

    drew_any = False
    for row, period in enumerate(PERIODS):
        reach_axis, zoom_axis = axes[row]
        best = best_cell(period, preset)
        if best is None:
            for axis in (reach_axis, zoom_axis):
                axis.text(0.5, 0.5,
                          f"{period}-{END_YEAR[period]} {preset}\nnot swept yet",
                          ha="center", va="center", color="#999999",
                          transform=axis.transAxes)
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

        # --- reach panel -------------------------------------------------
        reach_axis.plot(x, obs_start, marker="o", linestyle="--",
                        color="#999999", linewidth=1.6,
                        label=f"observed {period} (start)", zorder=3)
        reach_axis.plot(x, obs_end, marker="o", color=OBSERVED_COLOR,
                        linewidth=2.4,
                        label=f"observed {END_YEAR[period]} (end)", zorder=5)
        if nogroin is not None:
            reach_axis.plot(x, obs_start + np.array([nogroin[d] for d in have]),
                            marker="^", color=NOGROIN_COLOR, linewidth=1.8,
                            linestyle=":", label="model end, groin OFF",
                            zorder=4)
        if fitted is not None:
            reach_axis.plot(x, obs_start + np.array([fitted[d] for d in have]),
                            marker="s", color=MODEL_COLOR, linewidth=2.2,
                            label=f"model end, {_cell_label(best)}", zorder=6)

        # DOMAIN-COORDINATE shading belongs to the reach panel ONLY. The zoom
        # panel's x axis is alongshore METRES, so a span drawn at 5.5-6.5
        # lands six metres from the origin and drags the axis back to zero --
        # which is exactly what it did before this was split.
        reach_axis.axvspan(GROIN_UPDRIFT_GIS - 0.5, GROIN_UPDRIFT_GIS + 0.5,
                           color=GROIN_COLOR, alpha=0.10, zorder=0)
        reach_axis.axvspan(GROIN_DOWNDRIFT_GIS - 0.5,
                           GROIN_DOWNDRIFT_GIS + 0.5,
                           color="#1565C0", alpha=0.10, zorder=0)
        reach_axis.grid(alpha=0.25)
        zoom_axis.grid(alpha=0.25)

        misfit = (np.nan if fitted is None else
                  float(np.sqrt(np.mean(
                      (obs_start + np.array([fitted[d] for d in have])
                       - obs_end) ** 2))))
        reach_axis.set_title(
            f"{period}-{END_YEAR[period]}  {preset}   "
            f"end-position RMSE {misfit:.1f} m", fontsize=10.5)
        reach_axis.set_xlabel("GIS domain (south to north)")
        reach_axis.set_ylabel("shoreline position (m from datum, + seaward)")
        reach_axis.set_xticks(domains)
        reach_axis.legend(loc="best", fontsize=8)

        # --- zoom panel: transect resolution ------------------------------
        # X IS REAL ALONGSHORE DISTANCE, NOT THE DOMAIN ID. Plotting transects
        # against their integer domain puts every transect in a cell on the
        # same x, so a domain's whole spread draws as one vertical spike and
        # the fillet -- the entire point of this panel -- is unreadable.
        tr = transect_positions(chainage, period, ZOOM_DOMAINS)
        bounds = domain_alongshore_bounds(chainage, ZOOM_DOMAINS)
        if not tr.empty:
            zoom_axis.plot(tr["alongshore_m"], tr["start_m"], linestyle="--",
                           color="#999999", linewidth=1.4,
                           label=f"observed {period}", zorder=3)
            zoom_axis.plot(tr["alongshore_m"], tr["end_m"],
                           color=OBSERVED_COLOR, linewidth=2.0,
                           label=f"observed {END_YEAR[period]}", zorder=5)
        zoom_have = [d for d in have
                     if ZOOM_DOMAINS[0] <= d <= ZOOM_DOMAINS[1] and d in bounds]
        if fitted is not None and zoom_have:
            for index, domain in enumerate(zoom_have):
                lo_m, hi_m = bounds[domain]
                zoom_axis.hlines(start[domain] + fitted[domain], lo_m, hi_m,
                                 color=MODEL_COLOR, linewidth=2.6, zorder=6,
                                 label="model end (500 m cells)"
                                 if index == 0 else None)
            for domain in zoom_have:
                zoom_axis.axvline(bounds[domain][0], color="#DDDDDD",
                                  linewidth=0.8, zorder=0)
        for domain, colour in ((GROIN_UPDRIFT_GIS, GROIN_COLOR),
                               (GROIN_DOWNDRIFT_GIS, "#1565C0")):
            if domain in bounds:
                zoom_axis.axvspan(*bounds[domain], color=colour, alpha=0.10,
                                  zorder=0)
        zoom_axis.set_title(
            f"groin zoom D{ZOOM_DOMAINS[0]}-D{ZOOM_DOMAINS[1]}: "
            f"observations at ~62 m, model at 500 m", fontsize=10)
        if bounds:
            lo_m = min(v[0] for v in bounds.values())
            hi_m = max(v[1] for v in bounds.values())
            zoom_axis.set_xlim(lo_m - 60, hi_m + 60)
        zoom_axis.set_xlabel("alongshore distance (m, south to north)")
        zoom_axis.legend(loc="best", fontsize=8)

    figure.suptitle(
        f"Shoreline position: observed start and end vs model at the fitted "
        f"groin -- {preset}", fontsize=13)
    figure.tight_layout(rect=(0, 0.065, 1, 0.965))
    _footnote(
        figure,
        "Model and observed START at the same position BY CONSTRUCTION: the "
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
    figure_dir = PROJECT_BASE_DIR / "output" / "groin_sweep" / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    path = figure_dir / f"position_{preset}.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
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
