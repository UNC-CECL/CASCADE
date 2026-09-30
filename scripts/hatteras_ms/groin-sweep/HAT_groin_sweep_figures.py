#!/usr/bin/env python3
"""
Per-period diagnostic figures for one groin sweep: whether a winning cell is worth believing.

    python scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_figures.py
    python scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_figures.py --period 2004 --preset zeroBE

Written into the sweep's own folder. Details: scripts/hatteras_ms/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
# parents[3]: this file is in hatteras_ms/groin-sweep/; the guard makes a move fail loudly here
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
                              caption, error_cmap, figsize, open_frame,
                              save, _title)
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

# --- CONFIG ------------------------------------------------------------------
# Shared with HAT_groin_joint_fit.py so a groin is the same colour in every figure the sweep produces
MODEL_COLOR = C["ACCENT"]
GROIN_COLOR = INK_MUTED
OBSERVED_COLOR = INK

RANK_METRIC = "fillet_err"
REACH_METRIC = "rmse_window"
# -----------------------------------------------------------------------------


# Imports pyplot with a headless backend
def _matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


# Loading

# Loads one sweep's scored results
def load_scored(period, preset):
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


# Reduces rows to one per (M, f) by keeping the best-scoring be1
def profile_be1(frame):
    ranked = frame.sort_values(RANK_METRIC, na_position="last")
    return (ranked.drop_duplicates(subset=["M", "fraction"], keep="first")
            .reset_index(drop=True))


# The modelled per-domain LRR for one row, over the fit window
def rate_curve(row):
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    rate = np.array([float(row[f"rate_D{int(g)}"]) for g in gis])
    return gis, rate


# The CoastSat per-domain LRR for one period, over the fit window
def observed_curve(period):
    gis = np.array(FIT_DOMAINS_GIS, dtype=float)
    rate = np.array([OBSERVED_LRR[period][int(g)] for g in gis])
    return gis, rate


# Whether both groin domains depart the regional trend the SAME way
def observed_anomaly_is_one_signed(period, trend_order=2):
    gis, rate = observed_curve(period)
    keep = ~np.isin(gis, [GROIN_DOWNDRIFT_GIS, GROIN_UPDRIFT_GIS])
    trend = np.polyval(np.polyfit(gis[keep], rate[keep], trend_order), gis)
    anomaly = rate - trend
    up = float(anomaly[list(gis).index(float(GROIN_UPDRIFT_GIS))])
    down = float(anomaly[list(gis).index(float(GROIN_DOWNDRIFT_GIS))])
    return (up * down > 0.0), up, down


# The whole island, not just the fit window
_OBSERVED_FULL = {}


# CoastSat per-domain LRR across ALL 90 domains
def observed_curve_full(period):
    if period not in _OBSERVED_FULL:
        import pandas as _pd
        from cascade_pipeline.coastsat_lowess import compute_domain_means
        path = (COASTSAT_DIR / f"{period}_{END_YEAR[period]}"
                / "transect_lrr_full.csv")
        frame = _pd.read_csv(path)
        gis, means = compute_domain_means(
            frame["domain_number"].values, frame["lrr_m_yr"].values, 1, 90)
        _OBSERVED_FULL[period] = (np.asarray(gis, dtype=float),
                                  np.asarray(means, dtype=float))
    return _OBSERVED_FULL[period]


# One cell's modelled LRR across all 90 domains
def model_curve_full(period, preset, combo):
    path = sweep_output_dir(period, preset) / combo / "shoreline_change_rate.csv"
    if not path.exists():
        return None, None
    frame = pd.read_csv(path)
    return (frame["gis_domain"].to_numpy(dtype=float),
            frame["lrr_m_yr"].to_numpy(dtype=float))


# Shared axis furniture

# Draws the groin between the downdrift and updrift domains
def _mark_groin(axis):
    axis.axvline((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0,
                 color=GROIN_COLOR, linewidth=2.0, zorder=1)


# Shades the two domains the fillet is measured across
def _shade_pair(axis):
    # Two greys, not two hues: the bands say where; colour is for what is plotted
    axis.axvspan(GROIN_UPDRIFT_GIS - 0.5, GROIN_UPDRIFT_GIS + 0.5,
                 color="0.90", zorder=0)
    axis.axvspan(GROIN_DOWNDRIFT_GIS - 0.5, GROIN_DOWNDRIFT_GIS + 0.5,
                 color="0.94", zorder=0)


# Applies the shared labelling of a per-domain LRR panel
def _profile_axis(axis, period):
    _shade_pair(axis)
    _mark_groin(axis)
    axis.axhline(0.0, color=INK_MUTED, linewidth=0.8, linestyle=(0, (4, 3)),
                 zorder=1)
    axis.set_xlabel("GIS domain (south → north)")
    axis.set_ylabel("shoreline change rate (m/yr)\npositive is seaward")
    axis.set_xticks(list(FIT_DOMAINS_GIS))
    axis.grid(axis="y")
    axis.set_axisbelow(True)
    open_frame(axis)
    axis.text(GROIN_UPDRIFT_GIS + 0.15, axis.get_ylim()[1], " updrift",
              color=INK_MUTED, fontsize=7, va="top")
    axis.text(GROIN_DOWNDRIFT_GIS - 0.15, axis.get_ylim()[1], "downdrift ",
              color=INK_MUTED, fontsize=7, va="top", ha="right")


# Registers `text` as the figure's caption
def _footnote(figure, text, width=150):
    caption(figure, text)


# Adds a sentence to whatever caption the figure already carries
def _append_caption(figure, text):
    existing = getattr(figure, "_hat_caption", "")
    caption(figure, f"{existing} {text}".strip())


# Writes the one-signed-anomaly caveat under a profile figure, if it applies
def _notch_note(figure, period):
    one_signed, up, down = observed_anomaly_is_one_signed(period)
    if not one_signed:
        return
    _append_caption(
        figure,
        f"Observed D{GROIN_UPDRIFT_GIS} and D{GROIN_DOWNDRIFT_GIS} both "
        f"depart the regional trend the same way ({up:+.2f} and {down:+.2f} "
        f"m/yr), so the observed pair is not a dipole. A volume-neutral "
        f"source/sink must notch D{GROIN_DOWNDRIFT_GIS}; the mismatch there "
        f"is a property of the target, not of this cell.")


# Whether the reach RMSE panel is tracking bias rather than groin skill
def reach_panel_is_bias_driven(groin, threshold=-0.7):
    if "bias_window" not in groin.columns or len(groin) < 3:
        return False, float("nan"), float("nan")
    correlation = float(groin["M"].corr(groin["bias_window"]))
    mean_bias = float(groin["bias_window"].abs().mean())
    mean_rmse = float(groin[REACH_METRIC].abs().mean())
    # Both conditions: a strong correlation with little bias is just a mild response
    dominant = mean_rmse > 0 and (mean_bias / mean_rmse) > 0.7
    return (correlation <= threshold and dominant), correlation, mean_bias


# The best cell and every cell statistically tied with it
def tied_best(groin, column=RANK_METRIC, rel_tol=0.01):
    ranked = groin.sort_values(column)
    best = ranked.iloc[0]
    threshold = abs(float(best[column])) * rel_tol
    tied = ranked[(ranked[column] - float(best[column])).abs() <= threshold]
    return best, tied


# One line describing a tie, or None when the best cell is unique
def _tie_note(tied, column=RANK_METRIC):
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


# Short human label for one cell, with be1 only where it was swept
def _cell_label(row):
    label = f"M={row['M']:g}, f={row['fraction']:.2f}"
    if not pd.isna(row.get("be1")):
        label += f", be1={row['be1']:g}"
    return label


# Figures

# Two-panel M-f error surface
def fig_heatmap(period, preset, surface, out_dir):
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

    panel_notes = []
    panels = [
        (RANK_METRIC, "|modelled − observed| fillet size (m)",
         f"RANKED ON THIS. Observed fillet "
         f"{OBSERVED_FILLET_M[period]:.1f} m"),
        (REACH_METRIC, "LRR RMSE over D1-D12 (m/yr)", reach_subtitle),
    ]

    apply_style()
    figure, axes = plt.subplots(1, 2, figsize=figsize("double", aspect=0.42),
                                constrained_layout=True)
    for i, (axis, (column, label, subtitle)) in enumerate(zip(axes, panels)):
        grid = groin.pivot(index="fraction", columns="M", values=column)
        mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                               shading="nearest", cmap=error_cmap())
        cb = figure.colorbar(mesh, ax=axis)
        cb.set_label(label)
        cb.outline.set_linewidth(0.6)

        best, tied = tied_best(groin, column)
        panel_tie = _tie_note(tied, column)
        axis.plot(tied["M"], tied["fraction"], marker="*", markersize=13,
                  color=C["ACCENT"], markeredgecolor="white",
                  markeredgewidth=0.8, linestyle="none", zorder=5,
                  label=(f"{len(tied)} tied cells" if panel_tie
                         else f"best: {_cell_label(best)}"))

        rails = []
        if best["M"] in (min(m for m in M_VALUES if m > 0), max(M_VALUES)):
            rails.append("M")
        if best["fraction"] in (min(F_VALUES), max(F_VALUES)):
            rails.append("f")
        short = ("fillet size error" if column == RANK_METRIC
                 else "reach rate error")
        _title(axis, i, short if not rails
               else f"{short}, railed on {', '.join(rails)}")
        axis.set_xlabel("groin trapping rate M (m/yr)")
        axis.set_ylabel("deterioration floor f")
        axis.legend(loc="upper right", fontsize=7)
        panel_notes.append(f"({chr(ord('a') + i)}) {subtitle}.")

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

    _footnote(figure,
              "The groin sweep's error surface over the (M, f) grid for "
              "{period} to {end}, {preset}, scored two ways. Dark is worse and "
              "the marked cells are the best on each score. (a) is the fillet "
              "size error, which the sweep RANKS on; (b) is the error in the "
              "shoreline change rate over the whole D1\u2013D12 reach. {notes} "
              "{verdict}"
              .format(period=period, end=END_YEAR[period], preset=preset,
                      notes=" ".join(panel_notes), verdict=verdict))

    return save(figure, out_dir / "heatmap.png", close=True)[0]


# The winning cell's per-domain LRR against CoastSat
def fig_best_fit_profile(period, preset, surface, out_dir):
    plt = _matplotlib()

    groin = surface[surface["M"] > 0]
    best = groin.loc[groin[RANK_METRIC].idxmin()]
    gis, model = rate_curve(best)
    _, observed = observed_curve(period)

    apply_style()
    figure, (axis, full_axis) = plt.subplots(
        2, 1, figsize=figsize("double", aspect=0.80),
        gridspec_kw=dict(height_ratios=[1, 1]), constrained_layout=True)
    axis.plot(gis, observed, marker="o", markersize=3.4,
              color=OBSERVED_COLOR, linewidth=1.6,
              label="observed, CoastSat", zorder=4)
    axis.plot(gis, model, marker="s", markersize=3.4, color=MODEL_COLOR,
              linewidth=1.6, label=f"modelled, {_cell_label(best)}", zorder=3)
    _profile_axis(axis, period)

    # Full reach

    # The fit window is 12 of 90 domains
    reach_note = ""
    gis_full, obs_full = observed_curve_full(period)
    mod_gis, mod_full = model_curve_full(period, preset, best["combo"])
    full_axis.axvspan(min(FIT_DOMAINS_GIS) - 0.5, max(FIT_DOMAINS_GIS) + 0.5,
                      color="0.92", zorder=0, label="the scored fit window")
    full_axis.plot(gis_full, obs_full, color=OBSERVED_COLOR, linewidth=1.2,
                   label="observed, CoastSat", zorder=4)
    if mod_full is not None:
        full_axis.plot(mod_gis, mod_full, color=MODEL_COLOR, linewidth=1.2,
                       label=f"modelled, {_cell_label(best)}", zorder=3)
        common = np.intersect1d(gis_full, mod_gis)
        obs_i = np.array([obs_full[list(gis_full).index(g)] for g in common])
        mod_i = np.array([mod_full[list(mod_gis).index(g)] for g in common])
        outside = ~np.isin(common, FIT_DOMAINS_GIS)
        rmse_in = float(np.sqrt(np.nanmean((mod_i[~outside] - obs_i[~outside]) ** 2)))
        rmse_out = float(np.sqrt(np.nanmean((mod_i[outside] - obs_i[outside]) ** 2)))
        reach_note = ("RMSE inside the fit window is {:.2f} m/yr and "
                      "outside it {:.2f}.".format(rmse_in, rmse_out))
    _title(full_axis, 1, "the whole reach, D1 to D90")
    full_axis.axhline(0.0, color=INK_MUTED, linewidth=0.8,
                      linestyle=(0, (4, 3)), zorder=1)
    full_axis.axvline((GROIN_DOWNDRIFT_GIS + GROIN_UPDRIFT_GIS) / 2.0,
                      color=INK_MUTED, linewidth=0.8, zorder=2)
    full_axis.set_xlabel("GIS domain (south → north)")
    full_axis.set_ylabel("shoreline change rate (m/yr)\npositive is seaward")
    full_axis.grid(axis="y")
    full_axis.set_axisbelow(True)
    open_frame(full_axis)
    full_axis.legend(loc="best", fontsize=7)

    _title(axis, 0, "the fit window, domain by domain")
    axis.legend(loc="best", fontsize=7)

    _footnote(figure,
              "The best-scoring groin cell for {period} to {end}, {preset}: "
              "{cell}. Its modelled fillet is {fil:.1f} m against an observed "
              "{obs:.1f} m, an error of {err:.2f} m, and its error over the "
              "D1\u2013D12 reach is {reach:.3f} m/yr. (a) The scored fit "
              "window domain by domain, with the two shaded bands marking the "
              "domains either side of the structure. (b) The same cell over "
              "the whole 90-domain reach. The fit window is 12 of those 90, so "
              "a cell that matches the fillet says nothing on its own about "
              "the other 78; the groin's own influence is MEASURED, not "
              "fitted, out to a few kilometres, which makes this panel the "
              "place to check an emergent extent against observations that "
              "were never part of the objective. {reach_note}"
              .format(period=period, end=END_YEAR[period], preset=preset,
                      cell=_cell_label(best), fil=best["fillet_m"],
                      obs=OBSERVED_FILLET_M[period], err=best[RANK_METRIC],
                      reach=best[REACH_METRIC], reach_note=reach_note))
    _notch_note(figure, period)

    return save(figure, out_dir / "best_fit_profile.png", close=True)[0]


# The best N cells on one set of axes, against CoastSat
def fig_top_n_profiles(period, preset, surface, out_dir, n=5):
    plt = _matplotlib()

    groin = surface[surface["M"] > 0].sort_values(RANK_METRIC)
    top = groin.head(n)
    gis, observed = observed_curve(period)

    apply_style()
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.46),
                                constrained_layout=True)
    axis.plot(gis, observed, marker="o", markersize=3.4, color=OBSERVED_COLOR,
              linewidth=1.8, label="observed, CoastSat", zorder=5)

    # One colour family, darkest is best
    from matplotlib.colors import LinearSegmentedColormap
    ramp = LinearSegmentedColormap.from_list(
        "hat_accent_ramp", [C["ACCENT_FILL"], C["ACCENT"]])
    shades = ramp(np.linspace(1.0, 0.2, len(top)))
    for colour, (_, row) in zip(shades, top.iterrows()):
        _, model = rate_curve(row)
        axis.plot(gis, model, marker=".", markersize=3.0, color=colour,
                  linewidth=1.2,
                  label=f"{_cell_label(row)}, {row[RANK_METRIC]:.2f} m",
                  zorder=3)
    _profile_axis(axis, period)

    spread = float(top[RANK_METRIC].max() - top[RANK_METRIC].min())
    _title(axis, 0, f"the top {len(top)} cells by fillet error")
    axis.legend(loc="best", fontsize=7)

    _footnote(figure,
              "The {n} best-scoring cells for {period} to {end}, {preset}, on "
              "one set of axes against the observations, with the best cell "
              "darkest. Their scores span {spread:.2f} m. This is a statement "
              "about identifiability rather than about fit: curves lying on "
              "top of one another mean the metric cannot tell those cells "
              "apart. The two shaded bands are the domains either side of the "
              "structure."
              .format(n=len(top), period=period, end=END_YEAR[period],
                      preset=preset, spread=spread))
    _notch_note(figure, period)

    return save(figure, out_dir / "top_n_profiles.png", close=True)[0]


# The second period's error surface with constant-M*f contours drawn on
def fig_period2_surface(period, preset, surface, out_dir):
    plt = _matplotlib()

    groin = surface[surface["M"] > 0]
    grid = groin.pivot(index="fraction", columns="M", values=RANK_METRIC)

    apply_style()
    figure, axis = plt.subplots(figsize=figsize("double", aspect=0.48),
                                constrained_layout=True)
    mesh = axis.pcolormesh(grid.columns, grid.index, grid.values,
                           shading="nearest", cmap=error_cmap())
    cb = figure.colorbar(mesh, ax=axis)
    cb.set_label("|modelled − observed| fillet size (m)")
    cb.outline.set_linewidth(0.6)

    # Constant-M*f hyperbolae, anchored on the products the GRID spans rather than on the best cell
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
                  color=C["REF"], linewidth=1.0, zorder=4)
        axis.annotate(f"M·f = {product:.0f}",
                      xy=(m_axis[visible][-1], f_axis[visible][-1]),
                      color=C["REF"], fontsize=7, va="bottom", ha="right",
                      zorder=4,
                      bbox=dict(facecolor="white", alpha=0.75,
                                edgecolor="none",
                                boxstyle="square,pad=0.12"))

    # The valley floor: for each f, the M this period likes best.
    floor_M = [grid.columns[int(np.nanargmin(grid.loc[f].values))]
               if np.isfinite(grid.loc[f].values).any() else np.nan
               for f in grid.index]
    axis.plot(floor_M, grid.index, marker="o", markersize=3.0,
              color=C["ACCENT"], linewidth=1.2,
              label="valley floor: the best M at each f", zorder=5)

    # Draw the whole tied set, not one arbitrary member of it
    tie_note = _tie_note(tied)
    if tie_note:
        axis.plot(tied["M"], tied["fraction"], marker="o", markersize=5.5,
                  color="none", markeredgecolor=C["ACCENT"],
                  markeredgewidth=1.2, linestyle="none", zorder=6,
                  label=f"{len(tied)} tied cells, M unconstrained")
    else:
        axis.plot(best["M"], best["fraction"], marker="*", markersize=13,
                  color=C["ACCENT"], markeredgecolor="white",
                  markeredgewidth=0.8, linestyle="none", zorder=6,
                  label=f"best: {_cell_label(best)}")

    axis.set_xlabel("groin trapping rate M (m/yr)")
    axis.set_ylabel("deterioration floor f")
    axis.set_title("Only the product of the pair is identifiable here",
                   loc="left")
    axis.legend(loc="upper right", fontsize=7)

    note = ("The error surface for {} to {}, {}, with contours of constant "
            "M·f drawn on. This window sits entirely past the 2003 end of the "
            "deterioration ramp, so its cumulative trapping is 20·M·f and only "
            "the PRODUCT is identifiable: the surface is a valley running "
            "along a hyperbola rather than a bowl with a minimum. "
            .format(period, END_YEAR[period], preset)
            + "A valley floor that tracks a dashed contour means this period "
            "cannot separate M from f: read it as a bound on the product, not "
            "as a fitted M. Period 1 straddles the 1996-2003 ramp and is where "
            "the separation comes from.")
    if tie_note:
        note += (f" {tie_note.capitalize()}: the observed fillet here is "
                 f"{OBSERVED_FILLET_M[period]:+.1f} m, which no M >= 0 can "
                 f"build, so the score is minimised by trapping NOTHING and "
                 f"every M scores alike once f = 0.")
    _footnote(figure, note, width=165)

    return save(figure, out_dir / "period2_surface.png", close=True)[0]


# Driver

# Draws every figure for one swept period/preset
def figures_for(period, preset, top_n):
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


# Run: every swept cell's figures
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
