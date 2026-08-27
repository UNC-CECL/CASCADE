#!/usr/bin/env python3
"""Per-period diagnostic figures for one groin sweep.

The sweep orchestrator writes numbers; the joint fit draws the two-period
surface. Neither draws the per-period diagnostics that say WHETHER a winning
cell is any good: what the M-f error surface looks like on its own, and what
the winning shoreline-change curve looks like against CoastSat. This file
draws those, reading only what the sweep already wrote.

WHAT IS PLOTTED, AND WHY THAT AND NOT SOMETHING ELSE
    The heatmap carries TWO panels, not one, because the sweep's ranking
    metric and the reach-scale fit are different questions and can disagree:

        fillet_err   |modelled fillet - observed fillet| at the groin pair.
                     This is what the sweep RANKS on. The fillet saturates
                     (~18 m by 1990 at M = 40, where BRIE's diffusion balances
                     the trapping), so its LEVEL carries the information about
                     M while its slope -- which is what `differential` measures
                     -- carries almost none.
        rmse_window  RMSE of the modelled LRR against CoastSat across D1-D12.
                     A whole-reach number, dominated by domains the groin
                     never touches.

    If the two panels pick different cells that is a RESULT worth seeing, not
    a defect to average away: it means the groin parameters that reproduce the
    local fillet are not the ones that reproduce the reach.

    The profile figures plot the per-domain LRR, because that is the quantity
    the sweep is scored against -- the figure and the ranking then agree. The
    alternative (end-of-run shoreline position) shows the fillet's shape more
    directly but has no observed counterpart to overlay.

THE NOTCH CAVEAT, DRAWN RATHER THAN FOOTNOTED
    `observed_fillet_m` in the sweep config derives the observed fillet by
    fitting a regional trend across the fit window EXCLUDING the groin pair
    and taking the pair's departure from it. That anomaly is ONE-SIGNED --
    both D5 and D6 sit ABOVE the trend in 1984-2004 (+0.79 and +1.48 m/yr).
    A volume-neutral source/sink dipole must put a NOTCH at the downdrift
    domain, and the model does at every cell on the grid.

    So the observed shape is not a groin dipole, and the scalar fillet is the
    best available summary of a feature whose SHAPE the model cannot
    reproduce. Every profile figure that hits this case says so on its face,
    because a reader looking at a mismatched D5 should be told it is a
    property of the target and not a bad fit.

M = 0 IS NOT A COLUMN
    The M = 0 cells are the paired baselines the fillet is differenced
    against, so they have no fillet by definition and would be a blank column
    on the ranking panel. They are reported as a reference line in the panel
    subtitle instead, which is also the honest way to read them: the number to
    beat, not a candidate.

Usage:
    python HAT_groin_sweep_figures.py                     # every swept cell
    python HAT_groin_sweep_figures.py --period 2004 --preset zeroBE
    python HAT_groin_sweep_figures.py --top-n 8

Writes to output/groin_sweep/<start>_<end>_<preset>/figures/:
    heatmap.png              fillet error and reach RMSE over the M-f grid
    best_fit_profile.png     winning cell's LRR against CoastSat
    top_n_profiles.png       the best N cells on the same axes
    period2_surface.png      2004-2024 only: the M*f ridge, with contours

Author: Hannah A. Henry, UNC CECL
"""

from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
# The guard below is what makes a future move fail here, loudly, instead of
# resolving to scripts/scripts and surfacing as a missing data file several
# imports deeper.
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
    COASTSAT_DIR,
    END_YEAR,
    FIT_DOMAINS_GIS,
    F_VALUES,
    GROIN_DOWNDRIFT_GIS,
    GROIN_UPDRIFT_GIS,
    M_VALUES,
    OBSERVED_FILLET_M,
    OBSERVED_LRR,
    PERIODS,
    PRESETS,
    sweep_output_dir,
)

# Shared with HAT_groin_joint_fit.py so a groin is the same colour in every
# figure the sweep produces.
MODEL_COLOR = "#FF8C00"
GROIN_COLOR = "#B71C1C"
OBSERVED_COLOR = "#1A1A1A"

RANK_METRIC = "fillet_err"
REACH_METRIC = "rmse_window"


def _matplotlib():
    """Imports pyplot with a headless backend.

    Deferred rather than imported at module scope so `--list` and the argument
    parsing stay usable on a machine where matplotlib is missing.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


# =============================================================================
# LOADING
# =============================================================================

def load_scored(period, preset):
    """Loads one sweep's scored results.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".

    Returns:
        A DataFrame of scored rows, or None if that sweep has no CSV yet.

    Raises:
        ValueError: If the CSV predates the fillet-size rescore. Ranking such
            a frame would silently fall back to the differential -- the
            fillet's SLOPE -- which is the metric that railed at both grid
            edges. Re-running the orchestrator over a finished sweep is free
            (every cell resumes from disk) and adds the column, so this is a
            fixable error rather than a reason to plot the weaker metric.
    """
    path = sweep_output_dir(period, preset) / "sweep_results.csv"
    if not path.exists():
        return None
    frame = pd.read_csv(path)
    frame = frame[frame["differential_err"].notna()].copy()
    if frame.empty:
        return None
    if RANK_METRIC not in frame.columns:
        raise ValueError(
            f"{path} has no '{RANK_METRIC}' column, so it predates the "
            f"fillet-size rescore. Re-run:\n"
            f"    python HAT_groin_sweep.py --period {period} "
            f"--preset {preset}\n"
            f"Every cell is already on disk, so it only re-collates.")
    return frame


def profile_be1(frame):
    """Reduces rows to one per (M, f) by keeping the best-scoring be1.

    be1 is a nuisance axis here: it exists in the 1984 edgeBE sweep only, and
    the figure's subject is (M, f). Profiling it out -- keeping, for each
    (M, f), the be1 that scored best -- is what the joint fit does, so the
    surface drawn here and the surface fitted there are the same reduction.
    The winning be1 is carried along so it can be annotated rather than lost.

    Args:
        frame: Scored rows for one period/preset.

    Returns:
        A DataFrame with one row per (M, f), sorted by the ranking metric.
    """
    ranked = frame.sort_values(RANK_METRIC, na_position="last")
    return (ranked.drop_duplicates(subset=["M", "fraction"], keep="first")
            .reset_index(drop=True))


def rate_curve(row):
    """The modelled per-domain LRR for one row, over the fit window.

    Args:
        row: One row of a scored frame.

    Returns:
        (gis, rate) arrays over FIT_DOMAINS_GIS.
    """
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    rate = np.array([float(row[f"rate_D{int(g)}"]) for g in gis])
    return gis, rate


def observed_curve(period):
    """The CoastSat per-domain LRR for one period, over the fit window."""
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    rate = np.array([OBSERVED_LRR[period][int(g)] for g in gis])
    return gis, rate


def observed_anomaly_is_one_signed(period, trend_order=2):
    """Whether both groin domains depart the regional trend the SAME way.

    Mirrors the trend removal in `observed_fillet_m`. When true, the observed
    pair is not a dipole and no volume-neutral source/sink can reproduce its
    shape -- which the profile figures state on their face.

    Args:
        period: 1984 or 2004.
        trend_order: Polynomial order for the regional trend. Must match the
            default in `observed_fillet_m` or the two disagree about the same
            data.

    Returns:
        (is_one_signed, updrift_anomaly, downdrift_anomaly) in m/yr.
    """
    gis, rate = observed_curve(period)
    keep = ~np.isin(gis, [GROIN_DOWNDRIFT_GIS, GROIN_UPDRIFT_GIS])
    trend = np.polyval(np.polyfit(gis[keep], rate[keep], trend_order), gis)
    anomaly = rate - trend
    up = float(anomaly[list(gis).index(float(GROIN_UPDRIFT_GIS))])
    down = float(anomaly[list(gis).index(float(GROIN_DOWNDRIFT_GIS))])
    return (up * down > 0.0), up, down


# The whole island, not just the fit window. Loaded lazily and cached: this
# reads the transect file, and the module is imported by the comparison and
# position scripts too.
_OBSERVED_FULL = {}


def observed_curve_full(period):
    """CoastSat per-domain LRR across ALL 90 domains.

    `observed_curve` covers D1-D12 because that is the window the sweep is
    SCORED on. Nothing about the data stops at D12 -- both periods have
    transects in all 90 domains -- so the full reach is available for showing
    where a groin fitted on 12 domains leaves the other 78.

    Args:
        period: 1984 or 2004.

    Returns:
        (gis, rate) arrays over the domains that have transects.
    """
    if period not in _OBSERVED_FULL:
        import pandas as _pd
        from cascade_pipeline.coastsat_loess import compute_domain_means
        path = (COASTSAT_DIR / f"{period}_{END_YEAR[period]}"
                / "transect_lrr_full.csv")
        frame = _pd.read_csv(path)
        gis, means = compute_domain_means(
            frame["domain_number"].values, frame["lrr_m_yr"].values, 1, 90)
        _OBSERVED_FULL[period] = (np.asarray(gis, dtype=float),
                                  np.asarray(means, dtype=float))
    return _OBSERVED_FULL[period]


def model_curve_full(period, preset, combo):
    """One cell's modelled LRR across all 90 domains.

    Read from the cell's own `shoreline_change_rate.csv` rather than from the
    rate_D1..rate_D12 columns of sweep_results.csv, which carry the fit window
    only. `lrr_m_yr` is the column, not `change_rate_m_yr`: the sweep is scored
    on the OLS slope, and rate_D* is that same quantity.

    Returns:
        (gis, rate) arrays, or (None, None) if the cell has no rate file.
    """
    path = sweep_output_dir(period, preset) / combo / "shoreline_change_rate.csv"
    if not path.exists():
        return None, None
    frame = pd.read_csv(path)
    return (frame["gis_domain"].to_numpy(dtype=float),
            frame["lrr_m_yr"].to_numpy(dtype=float))


# =============================================================================
# SHARED AXIS FURNITURE
# =============================================================================

def _mark_groin(axis):
    """Draws the groin between the downdrift and updrift domains."""
    axis.axvline((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0,
                 color=GROIN_COLOR, linewidth=2.0, zorder=1)


def _shade_pair(axis):
    """Shades the two domains the fillet is measured across."""
    axis.axvspan(GROIN_UPDRIFT_GIS - 0.5, GROIN_UPDRIFT_GIS + 0.5,
                 color=GROIN_COLOR, alpha=0.10, zorder=0)
    axis.axvspan(GROIN_DOWNDRIFT_GIS - 0.5, GROIN_DOWNDRIFT_GIS + 0.5,
                 color="#1565C0", alpha=0.10, zorder=0)


def _profile_axis(axis, period):
    """Applies the shared labelling of a per-domain LRR panel."""
    _shade_pair(axis)
    _mark_groin(axis)
    axis.axhline(0.0, color="#888888", linewidth=0.8, zorder=1)
    axis.set_xlabel("GIS domain (south to north)")
    axis.set_ylabel("shoreline change rate (m/yr, + = seaward)")
    axis.set_xticks(list(FIT_DOMAINS_GIS))
    axis.grid(alpha=0.25)
    axis.text(GROIN_UPDRIFT_GIS + 0.15, axis.get_ylim()[1], " updrift",
              color=GROIN_COLOR, fontsize=8, va="top")
    axis.text(GROIN_DOWNDRIFT_GIS - 0.15, axis.get_ylim()[1], "downdrift ",
              color="#1565C0", fontsize=8, va="top", ha="right")


def _footnote(figure, text, width=150):
    """Writes a wrapped footnote under a figure.

    Wrapped with `textwrap` at a fixed column rather than matplotlib's
    `wrap=True`, which measures against the figure edge and clips the last
    word of a long line.
    """
    figure.text(0.01, 0.005, "\n".join(textwrap.wrap(text, width)),
                fontsize=7.5, color="#444444", va="bottom")


def _notch_note(figure, period):
    """Writes the one-signed-anomaly caveat under a profile figure, if it applies."""
    one_signed, up, down = observed_anomaly_is_one_signed(period)
    if not one_signed:
        return
    _footnote(
        figure,
        f"Observed D{GROIN_UPDRIFT_GIS} and D{GROIN_DOWNDRIFT_GIS} both "
        f"depart the regional trend the same way ({up:+.2f} and {down:+.2f} "
        f"m/yr), so the observed pair is not a dipole. A volume-neutral "
        f"source/sink must notch D{GROIN_DOWNDRIFT_GIS}; the mismatch there "
        f"is a property of the target, not of this cell.")


def reach_panel_is_bias_driven(groin, threshold=-0.7):
    """Whether the reach RMSE panel is tracking bias rather than groin skill.

    THE FAILURE THIS CATCHES, measured on the 1984 zeroBE sweep. With
    background erosion switched off the whole reach is far too accretional
    (bias +2.50 m/yr against CoastSat), and `rmse_window` is dominated by that
    offset rather than by shape. Trapping sediment adds net erosion, so
    raising M shaves the bias and the reach RMSE falls MONOTONICALLY with M --
    corr(M, bias) = -0.96 -- until it rails at the largest M on the grid.

    A reader seeing that panel rail at M = 110 would reasonably conclude the
    reach fit wants an enormous groin. It does not: it wants background
    erosion, and M is the only knob on the grid that can supply any. Flagged
    on the figure so the railing cannot be quoted as a groin result.

    Args:
        groin: Scored rows with M > 0 for one period/preset.
        threshold: Correlation at or below which bias is called dominant.

    Returns:
        (is_bias_driven, correlation, mean_absolute_bias).
    """
    if "bias_window" not in groin.columns or len(groin) < 3:
        return False, float("nan"), float("nan")
    correlation = float(groin["M"].corr(groin["bias_window"]))
    mean_bias = float(groin["bias_window"].abs().mean())
    mean_rmse = float(groin[REACH_METRIC].abs().mean())
    # Both conditions matter: a strong correlation with a NEGLIGIBLE bias is
    # just a well-fit reach responding mildly to M.
    dominant = mean_rmse > 0 and (mean_bias / mean_rmse) > 0.7
    return (correlation <= threshold and dominant), correlation, mean_bias


def tied_best(groin, column=RANK_METRIC, rel_tol=0.01):
    """The best cell and every cell statistically tied with it.

    WHY THIS IS NOT `idxmin`. In 2004-2024 the observed fillet is NEGATIVE
    (-43.2 m: the fillet RELAXED across the window), and no M >= 0 can build a
    negative fillet. The whole f = 0 row therefore scores identically -- a
    fully deteriorated groin traps nothing, so M has no effect at all -- and on
    the 2026-08-23 sweep seven cells from M = 40 to M = 110 tied to within
    5e-5 m. `idxmin` picks one of them arbitrarily and the figure then
    announces "best: M=95", which is not a fitted value: it is whichever tied
    cell numpy happened to reach first.

    Reporting the tie is the honest version, because the tie IS the result --
    it says this period cannot constrain M.

    Args:
        groin: Scored rows with M > 0.
        column: Metric to rank on.
        rel_tol: Fraction of the best score within which a cell counts as tied.

    Returns:
        (best_row, tied_frame). `tied_frame` always contains at least the best
        row, so `len(tied) > 1` is the test for an unidentified parameter.
    """
    ranked = groin.sort_values(column)
    best = ranked.iloc[0]
    threshold = abs(float(best[column])) * rel_tol
    tied = ranked[(ranked[column] - float(best[column])).abs() <= threshold]
    return best, tied


def _tie_note(tied, column=RANK_METRIC):
    """One line describing a tie, or None when the best cell is unique."""
    if len(tied) <= 1:
        return None
    m_span = (tied["M"].min(), tied["M"].max())
    f_span = (tied["fraction"].min(), tied["fraction"].max())
    free = []
    if m_span[0] != m_span[1]:
        free.append(f"M spans {m_span[0]:g}-{m_span[1]:g}")
    if f_span[0] != f_span[1]:
        free.append(f"f spans {f_span[0]:g}-{f_span[1]:g}")
    return (f"{len(tied)} cells tied within 1% of the best score"
            + (f" ({'; '.join(free)})" if free else "")
            + " -- not a fitted optimum")


def _cell_label(row):
    """Short human label for one cell, with be1 only where it was swept."""
    label = f"M={row['M']:g}, f={row['fraction']:.2f}"
    if not pd.isna(row.get("be1")):
        label += f", be1={row['be1']:g}"
    return label


# =============================================================================
# FIGURES
# =============================================================================

def fig_heatmap(period, preset, surface, out_dir):
    """Two-panel M-f error surface: the ranking metric and the reach fit.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".
        surface: One row per (M, f), from `profile_be1`.
        out_dir: Directory to write into.

    Returns:
        The written path.
    """
    plt = _matplotlib()

    groin = surface[surface["M"] > 0]
    baseline = surface[surface["M"] == 0]
    baseline_rmse = (float(baseline[REACH_METRIC].min())
                     if not baseline.empty else np.nan)

    bias_driven, bias_corr, mean_bias = reach_panel_is_bias_driven(groin)
    reach_subtitle = (f"no-groin baseline {baseline_rmse:.3f} m/yr"
                      if np.isfinite(baseline_rmse)
                      else "no M = 0 baseline on disk")
    if bias_driven:
        reach_subtitle += (f"   |   BIAS-DRIVEN: mean |bias| {mean_bias:.2f} "
                           f"m/yr, corr(M, bias) = {bias_corr:+.2f}")

    panels = [
        (RANK_METRIC, "|modelled - observed| fillet size (m)",
         f"RANKED ON THIS. Observed fillet "
         f"{OBSERVED_FILLET_M[period]:.1f} m"),
        (REACH_METRIC, "LRR RMSE over D1-D12 (m/yr)", reach_subtitle),
    ]

    figure, axes = plt.subplots(1, 2, figsize=(14, 5.4))
    for axis, (column, label, subtitle) in zip(axes, panels):
        grid = groin.pivot(index="fraction", columns="M", values=column)
        mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                               shading="nearest", cmap="viridis_r")
        figure.colorbar(mesh, ax=axis, label=label)

        best, tied = tied_best(groin, column)
        panel_tie = _tie_note(tied, column)
        axis.plot(tied["M"], tied["fraction"], marker="*", markersize=20,
                  color=GROIN_COLOR, markeredgecolor="white",
                  markeredgewidth=1.2, linestyle="none", zorder=5,
                  label=(f"{len(tied)} tied cells" if panel_tie
                         else f"best: {_cell_label(best)}"))

        rails = []
        if best["M"] in (min(m for m in M_VALUES if m > 0), max(M_VALUES)):
            rails.append("M")
        if best["fraction"] in (min(F_VALUES), max(F_VALUES)):
            rails.append("f")
        title = column if not rails else f"{column}   (railed on {', '.join(rails)})"

        axis.set_title(f"{title}\n{subtitle}", fontsize=10)
        axis.set_xlabel("groin trapping rate M (m/yr)")
        axis.set_ylabel("deterioration floor f")
        axis.legend(loc="upper right", fontsize=8)

    # Whether the two panels agree is the point of drawing both.
    pick_rank = groin.loc[groin[RANK_METRIC].idxmin()]
    pick_reach = groin.loc[groin[REACH_METRIC].idxmin()]
    agree = (pick_rank["M"] == pick_reach["M"]
             and pick_rank["fraction"] == pick_reach["fraction"])
    verdict = ("Both panels pick the same cell."
               if agree else
               f"The panels DISAGREE: fillet picks {_cell_label(pick_rank)}, "
               f"reach RMSE picks {_cell_label(pick_reach)}. The groin "
               f"parameters that reproduce the local fillet are not the ones "
               f"that best fit the reach.")
    if bias_driven:
        verdict += (
            f" Read the right panel with care: mean |bias| is "
            f"{mean_bias:.2f} m/yr and corr(M, bias) = {bias_corr:+.2f}, so "
            f"reach RMSE falls with M because trapping shaves a whole-reach "
            f"offset, not because a larger groin fits better. Its rail is "
            f"asking for background erosion, which M is the only knob on this "
            f"grid able to supply.")

    figure.suptitle(
        f"Groin sweep error surface -- {period}-{END_YEAR[period]}  {preset}",
        fontsize=12)
    figure.tight_layout(rect=(0, 0.075, 1, 0.96))
    _footnote(figure, verdict, width=185)

    path = out_dir / "heatmap.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def fig_best_fit_profile(period, preset, surface, out_dir):
    """The winning cell's per-domain LRR against CoastSat."""
    plt = _matplotlib()

    groin = surface[surface["M"] > 0]
    best = groin.loc[groin[RANK_METRIC].idxmin()]
    gis, model = rate_curve(best)
    _, observed = observed_curve(period)

    figure, (axis, full_axis) = plt.subplots(
        2, 1, figsize=(11, 9), gridspec_kw=dict(height_ratios=[1, 1]))
    axis.plot(gis, observed, marker="o", color=OBSERVED_COLOR, linewidth=2.0,
              label="CoastSat observed", zorder=4)
    axis.plot(gis, model, marker="s", color=MODEL_COLOR, linewidth=2.0,
              label=f"model: {_cell_label(best)}", zorder=3)
    _profile_axis(axis, period)

    # --- full reach ------------------------------------------------------
    # The fit window is 12 of 90 domains. A cell that matches the fillet says
    # nothing on its own about the other 78, and the groin's own influence is
    # measured (not fitted) out to a few km -- so this panel is where an
    # emergent extent can be checked against the observations that were never
    # part of the objective.
    gis_full, obs_full = observed_curve_full(period)
    mod_gis, mod_full = model_curve_full(period, preset, best["combo"])
    full_axis.axvspan(min(FIT_DOMAINS_GIS) - 0.5, max(FIT_DOMAINS_GIS) + 0.5,
                      color="#FFD54F", alpha=0.22, zorder=0,
                      label="fit window (scored)")
    full_axis.plot(gis_full, obs_full, color=OBSERVED_COLOR, linewidth=1.6,
                   label="CoastSat observed", zorder=4)
    if mod_full is not None:
        full_axis.plot(mod_gis, mod_full, color=MODEL_COLOR, linewidth=1.6,
                       label=f"model: {_cell_label(best)}", zorder=3)
        common = np.intersect1d(gis_full, mod_gis)
        obs_i = np.array([obs_full[list(gis_full).index(g)] for g in common])
        mod_i = np.array([mod_full[list(mod_gis).index(g)] for g in common])
        outside = ~np.isin(common, FIT_DOMAINS_GIS)
        rmse_in = float(np.sqrt(np.nanmean((mod_i[~outside] - obs_i[~outside]) ** 2)))
        rmse_out = float(np.sqrt(np.nanmean((mod_i[outside] - obs_i[outside]) ** 2)))
        full_axis.set_title(
            f"full reach D1-D90   |   RMSE inside fit window {rmse_in:.2f}, "
            f"outside {rmse_out:.2f} m/yr", fontsize=10)
    full_axis.axhline(0.0, color="#888888", linewidth=0.8, zorder=1)
    full_axis.axvline((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0,
                      color=GROIN_COLOR, linewidth=1.6, zorder=2)
    full_axis.set_xlabel("GIS domain (south to north)")
    full_axis.set_ylabel("shoreline change rate (m/yr, + = seaward)")
    full_axis.grid(alpha=0.25)
    full_axis.legend(loc="best", fontsize=8)

    axis.set_title(
        f"Best-fit groin cell -- {period}-{END_YEAR[period]}  {preset}\n"
        f"fillet {best['fillet_m']:.1f} m vs observed "
        f"{OBSERVED_FILLET_M[period]:.1f} m "
        f"(err {best[RANK_METRIC]:.2f} m)   |   "
        f"D1-D12 RMSE {best[REACH_METRIC]:.3f} m/yr",
        fontsize=11)
    axis.legend(loc="best", fontsize=9)
    figure.tight_layout(rect=(0, 0.075, 1, 1))
    _notch_note(figure, period)

    path = out_dir / "best_fit_profile.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def fig_top_n_profiles(period, preset, surface, out_dir, n=5):
    """The best N cells on one set of axes, against CoastSat.

    Shows how tightly the ranking discriminates: curves that lie on top of
    each other mean the metric cannot tell those cells apart, which is a
    statement about identifiability rather than about the fit.
    """
    plt = _matplotlib()

    groin = surface[surface["M"] > 0].sort_values(RANK_METRIC)
    top = groin.head(n)
    gis, observed = observed_curve(period)

    figure, axis = plt.subplots(figsize=(10, 5.6))
    axis.plot(gis, observed, marker="o", color=OBSERVED_COLOR, linewidth=2.5,
              label="CoastSat observed", zorder=5)

    shades = plt.get_cmap("autumn")(np.linspace(0.0, 0.7, len(top)))
    for colour, (_, row) in zip(shades, top.iterrows()):
        _, model = rate_curve(row)
        axis.plot(gis, model, marker=".", color=colour, linewidth=1.5,
                  label=f"{_cell_label(row)}  (err {row[RANK_METRIC]:.2f} m)",
                  zorder=3)
    _profile_axis(axis, period)

    spread = float(top[RANK_METRIC].max() - top[RANK_METRIC].min())
    axis.set_title(
        f"Top {len(top)} cells by fillet error -- "
        f"{period}-{END_YEAR[period]}  {preset}\n"
        f"error spread across these cells: {spread:.2f} m",
        fontsize=11)
    axis.legend(loc="best", fontsize=8)
    figure.tight_layout(rect=(0, 0.075, 1, 1))
    _notch_note(figure, period)

    path = out_dir / "top_n_profiles.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


def fig_period2_surface(period, preset, surface, out_dir):
    """The second period's error surface with constant-M*f contours drawn on.

    WHY THIS FIGURE EXISTS. 2004-2024 sits entirely past the 2003 end of the
    deterioration ramp, so its cumulative trapping is 20*M*f: only the PRODUCT
    is identifiable and the surface is a valley running along a hyperbola, not
    a bowl with a minimum. Drawing the constant-M*f contours makes that
    visible -- if the valley floor follows a contour, the non-identifiability
    is shown rather than asserted in a caption.

    Args:
        period: Must be the second period; the figure is meaningless for the
            first, which straddles the ramp.
        preset: "edgeBE" or "zeroBE".
        surface: One row per (M, f), from `profile_be1`.
        out_dir: Directory to write into.

    Returns:
        The written path.
    """
    plt = _matplotlib()

    groin = surface[surface["M"] > 0]
    grid = groin.pivot(index="fraction", columns="M", values=RANK_METRIC)

    figure, axis = plt.subplots(figsize=(9.5, 5.8))
    mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                           shading="nearest", cmap="viridis_r")
    figure.colorbar(mesh, ax=axis,
                    label="|modelled - observed| fillet size (m)")

    # Constant-M*f hyperbolae, anchored on the products the GRID spans rather
    # than on the best cell. Anchoring on the best cell silently drew nothing
    # here: period 2's optimum sits at f = 0, so every product was 0 and each
    # contour was skipped as non-positive -- the one figure whose entire point
    # is to show the M*f ridge came out with no ridge on it.
    m_axis = np.linspace(min(grid.columns), max(grid.columns), 200)
    best, tied = tied_best(groin)
    all_products = (groin["M"] * groin["fraction"])
    all_products = all_products[all_products > 0]
    products = ([float(all_products.quantile(q)) for q in (0.25, 0.5, 0.75)]
                if not all_products.empty else [])
    for product in products:
        if product <= 0:
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            f_axis = product / m_axis
        visible = (f_axis >= min(F_VALUES)) & (f_axis <= max(F_VALUES))
        if not visible.any():
            continue
        axis.plot(m_axis[visible], f_axis[visible], linestyle="--",
                  color="white", linewidth=1.2, alpha=0.9, zorder=4)
        axis.annotate(f"M·f = {product:.0f}",
                      xy=(m_axis[visible][-1], f_axis[visible][-1]),
                      color="white", fontsize=8, va="bottom", ha="right",
                      zorder=4)

    # The valley floor: for each f, the M this period likes best.
    floor_M = [grid.columns[int(np.nanargmin(grid.loc[f].values))]
               if np.isfinite(grid.loc[f].values).any() else np.nan
               for f in grid.index]
    axis.plot(floor_M, grid.index, marker="o", color=GROIN_COLOR,
              linewidth=1.5, label="valley floor (best M at each f)", zorder=5)

    # Draw the whole tied set, not one arbitrary member of it. A single star
    # on a tie reads as a fitted value; a row of them reads as what it is.
    tie_note = _tie_note(tied)
    if tie_note:
        axis.plot(tied["M"], tied["fraction"], marker="o", markersize=9,
                  color="white", markeredgecolor=GROIN_COLOR,
                  markeredgewidth=1.4, linestyle="none", zorder=6,
                  label=f"{len(tied)} tied cells (M unconstrained)")
    else:
        axis.plot(best["M"], best["fraction"], marker="*", markersize=20,
                  color="white", markeredgecolor=GROIN_COLOR,
                  markeredgewidth=1.4, linestyle="none", zorder=6,
                  label=f"best: {_cell_label(best)}")

    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.set_title(
        f"{period}-{END_YEAR[period]} {preset}: only the PRODUCT M·f is "
        f"identifiable\n"
        f"the run sits entirely past the 2003 ramp, so cumulative trapping is "
        f"20·M·f",
        fontsize=11)
    axis.legend(loc="upper right", fontsize=8)

    figure.tight_layout(rect=(0, 0.075, 1, 1))
    note = ("A valley floor that tracks a dashed contour means this period "
            "cannot separate M from f: read it as a bound on the product, not "
            "as a fitted M. Period 1 straddles the 1996-2003 ramp and is where "
            "the separation comes from.")
    if tie_note:
        note += (f" {tie_note.capitalize()}: the observed fillet here is "
                 f"{OBSERVED_FILLET_M[period]:+.1f} m, which no M >= 0 can "
                 f"build, so the score is minimised by trapping NOTHING and "
                 f"every M scores alike once f = 0.")
    _footnote(figure, note, width=165)

    path = out_dir / "period2_surface.png"
    figure.savefig(path, dpi=150, facecolor="white")
    plt.close(figure)
    return path


# =============================================================================
# DRIVER
# =============================================================================

def figures_for(period, preset, top_n):
    """Draws every figure for one swept period/preset.

    Args:
        period: 1984 or 2004.
        preset: "edgeBE" or "zeroBE".
        top_n: How many cells the top-N profile figure carries.

    Returns:
        A list of written paths, empty if that sweep has no results yet.
    """
    frame = load_scored(period, preset)
    if frame is None:
        print(f"  {period}-{END_YEAR[period]} {preset:<8} no sweep_results.csv "
              f"-- skipped")
        return []

    surface = profile_be1(frame)
    n_groin = int((surface["M"] > 0).sum())
    if n_groin == 0:
        print(f"  {period}-{END_YEAR[period]} {preset:<8} only M = 0 baselines "
              f"scored -- nothing to plot")
        return []

    out_dir = sweep_output_dir(period, preset) / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    written = [
        fig_heatmap(period, preset, surface, out_dir),
        fig_best_fit_profile(period, preset, surface, out_dir),
        fig_top_n_profiles(period, preset, surface, out_dir, top_n),
    ]
    if period == PERIODS[1]:
        written.append(fig_period2_surface(period, preset, surface, out_dir))

    best, tied = tied_best(surface[surface["M"] > 0])
    tie = _tie_note(tied)
    print(f"  {period}-{END_YEAR[period]} {preset:<8} {n_groin:>3} cells   "
          f"best {_cell_label(best)}   fillet err {best[RANK_METRIC]:.2f} m   "
          f"-> {out_dir}")
    if tie:
        print(f"     ^ {tie}")
    return written


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--period", type=int, choices=PERIODS, action="append",
                        help="repeatable; default every period")
    parser.add_argument("--preset", choices=PRESETS, action="append",
                        help="repeatable; default every preset")
    parser.add_argument("--top-n", type=int, default=5,
                        help="cells in the top-N profile figure (default 5)")
    args = parser.parse_args()

    periods = args.period or list(PERIODS)
    presets = args.preset or list(PRESETS)

    print("=" * 72)
    print("GROIN SWEEP FIGURES")
    print("=" * 72)

    written = []
    for period in periods:
        for preset in presets:
            written.extend(figures_for(period, preset, args.top_n))

    print("=" * 72)
    print(f"  {len(written)} figures written")
    return 0 if written else 1


if __name__ == "__main__":
    sys.exit(main())
