r"""
HAT_metres_1_offset_units_plot.py -- the figures for the offset-scale study
==============================================================================
Reads tables/all_runs.csv (written by `HAT_metres_1_offset_units.py score`),
each scored run's tables/shoreline_change_rate.csv, and the CoastSat target the
runner scores against. Writes, under figures/ in the study folder:

  wave_height/rmse_bias_vs_wave_height_by_offset_scale_1996_2010.png
      the offset each scale hands BRIE; interior RMSE and bias against Hs
  wave_angle/rmse_bias_vs_high_angle_fraction_by_asymmetry_1996_2010.png
      RMSE and bias against the high-angle fraction, one line per
      asymmetry, one column per offset scale
  combined/best_rmse_by_offset_scale_1996_2010.png
      the RMSE range each scale reaches in each sweep, best run labelled
  combined/bias_vs_rmse_all_runs_1996_2010.png
      every scored run of both sweeps, bias against RMSE
  combined/alongshore_rate_best_runs_vs_coastsat_1996_2010.png
      the best run of each scale against the CoastSat target, GIS 1-90
  combined/variance_explained_by_offset_scale_1996_2010.png
      the share of the observed alongshore variation each scale explains
  combined/spread_vs_correlation_all_runs_1996_2010.png
      each run's alongshore spread against its correlation with CoastSat
  alongshore_sensitivity/<scale>/rate_and_position_change_by_<parameter>_<scale>_1996_2010.png
      one parameter moved (Hs, high-angle fraction, asymmetry), the others at
      their defaults: (a) rate vs the CoastSat LRR target, (b) position change
      vs the observed CoastSat change 1996 -> 2010

Scores are the runner's (interior GIS 2-89, LRR, CoastSat LOESS 10-domain
target), plus the alongshore-variation scores `score` adds from the same
target (study.coastsat_target, checked there against the runner's RMSE).
==============================================================================
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_metres_1_offset_units as study  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1997, INK, INK_MUTED, DOMAIN_AXIS_LABEL, SMOOTH_RAMP, _title,
    apply_style, figsize, open_frame, record_caption, save, structures,
    support_dir, town_bands)
from site_layer.hat_topo_version import INIT_ROOT, offset_file  # noqa: E402
from site_layer.hatteras_site_config import (  # noqa: E402
    HATTERAS_DOMAINS, SCORE_INTERIOR_GIS)
from cascade_pipeline.hindcast import build_island_offset  # noqa: E402

FIG_DIR = study.STUDY_DIR / "figures"
WINDOW = "1996_2010"

SCALES = tuple(study.SCALES)          # div10, metres, metres-detrended
SCALE_STYLE = {
    "div10":            dict(color=C["BASE"],   label="Offset ÷ 10 (calibrated runs)"),
    "metres":           dict(color=C["ADDED"],  label="Offset in metres"),
    "metres-detrended": dict(color=C["ACCENT"], label="Offset in metres, trend removed"),
}
# short enough to clear the panel letter in a three-column row
PANEL_TITLE = {"div10": "Offset ÷ 10", "metres": "Offset in metres",
               "metres-detrended": "Metres, trend removed"}
SOURCE_STYLE = {
    "duneline":  dict(ls="-",  filled=True,  label="Dune-line offset"),
    "shoreline": dict(ls="--", filled=False, label="Shoreline offset"),
}
SWEEP_LABEL = {"wave_height": "Wave height tuned",
               "wave_angle": "Wave angles tuned"}
SWEEP_MARKER = {"wave_height": "o", "wave_angle": "^"}
HIGH_ANGLE_LABEL = "Fraction of high-angle waves (> 45°)"


def load_runs():
    df = pd.read_csv(study.TABLES_DIR / "all_runs.csv")
    df["scored"] = df["status"] == "scored"
    return df


def run_rates(row):
    return study.run_rates(study.STUDY_DIR / row.run_dir)


# One definition of the target and the interior, shared with `score`.
coastsat_target = study.coastsat_target
interior = study.interior


def flat_line_rmse(target):
    """RMSE of predicting the observed interior mean at every domain: the
    score a model with no alongshore pattern at all would get."""
    t = interior(target)
    return float(np.sqrt(((t - t.mean()) ** 2).mean()))


NULL_STYLE = dict(color=INK, lw=0.9, ls="--")
NULL_LABEL = "Flat line at the observed mean"


def log_rmse_axis(ax):
    ax.set_yscale("log")
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator([1, 2, 5, 10]))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def scale_handles():
    return [Line2D([], [], color=SCALE_STYLE[s]["color"], lw=1.6,
                   label=SCALE_STYLE[s]["label"]) for s in SCALES]


def caption_common():
    return ("1996-2010 hindcast, no background erosion, full management, no "
            "groin, relocations off. Scores are the modelled LRR shoreline-change "
            "rate against the CoastSat LRR target (LOESS, 10 domains) over the "
            "interior domains GIS 2-89. The tuned parameters are fitted on the "
            "window they are scored on, so a best value is a band, not a "
            "calibrated value.")


# =============================================================================
# 1. WAVE HEIGHT
# =============================================================================

def fig_wave_height(df):
    d = df[df.sweep == "wave_height"]
    fig = plt.figure(figsize=figsize("double", height=6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15])
    ax_a, ax_b, ax_c = (fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]),
                        fig.add_subplot(gs[1, 1]))

    path = INIT_ROOT / offset_file(1996, "padded", HATTERAS_DOMAINS.total_domains,
                                   source="duneline")
    lo, hi = HATTERAS_DOMAINS.start_real_index, HATTERAS_DOMAINS.end_real_index
    gis = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
    for scale in SCALES:
        arr = build_island_offset(path, HATTERAS_DOMAINS,
                                  mode=study.SCALES[scale])[lo:hi]
        ax_a.plot(gis, (arr - arr.mean()) / 1000.0, lw=1.4,
                  color=SCALE_STYLE[scale]["color"])
    ax_a.axhline(0, color=INK_MUTED, lw=0.6)
    ax_a.set_xlim(gis[0], gis[-1])
    ax_a.set_xlabel(DOMAIN_AXIS_LABEL)
    ax_a.set_ylabel("Offset given to BRIE,\nrelative to its mean (km)")
    ax_a.grid(axis="y")
    open_frame(ax_a)
    _title(ax_a, 0, "Island offset as given to the model (dune-line source)")

    for ax, col, ylab, i, title in (
            (ax_b, "rmse_interior_m_yr", "Interior RMSE (m/yr)", 1, "Error against CoastSat"),
            (ax_c, "mean_bias_interior_m_yr", "Interior mean bias (m/yr)", 2, "Bias against CoastSat")):
        for scale in SCALES:
            for source, s in SOURCE_STYLE.items():
                g = d[(d.offset_scale == scale) & (d.offset_source == source) & d.scored]
                g = g.sort_values("Hs_m")
                colr = SCALE_STYLE[scale]["color"]
                ax.plot(g.Hs_m, g[col], ls=s["ls"], marker="o", ms=4, lw=1.3,
                        color=colr, mec=colr, mfc=colr if s["filled"] else "white")
        # drowned cells, at the foot of the panel, so a missing point says why
        # drowned cells, at the foot of the panel in the scale's colour, one
        # row per scale, so a missing point says which run and why
        drowned = d[~d.scored]
        for _, g in drowned.iterrows():
            ax.plot(g.Hs_m, 0.03 + 0.045 * SCALES.index(g.offset_scale), marker="x", ms=5,
                    mew=1.4, color=SCALE_STYLE[g.offset_scale]["color"],
                    transform=ax.get_xaxis_transform(), clip_on=False)
        if col.startswith("rmse"):
            log_rmse_axis(ax)
        else:
            ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlabel("Significant wave height, Hs (m)")
        ax.set_ylabel(ylab)
        # label every half metre plus 0.75; the other Hs run values (0.6,
        # 0.65) are unlabelled minor ticks, too close to label
        run_hs = sorted(d.Hs_m.unique())
        major = [h for h in run_hs if np.isclose(h * 2, round(h * 2)) or np.isclose(h, 0.75)]
        ax.set_xticks(major)
        ax.set_xticks([h for h in run_hs if h not in major], minor=True)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
        ax.grid(axis="y", which="both")
        open_frame(ax)
        _title(ax, i, title)

    handles = scale_handles() + [
        Line2D([], [], color=INK_MUTED, ls=s["ls"], marker="o", ms=4,
               mfc=INK_MUTED if s["filled"] else "white", mec=INK_MUTED,
               label=s["label"]) for s in SOURCE_STYLE.values()] + [
        Line2D([], [], color=INK, ls="none", marker="x", ms=6,
               label="Barrier drowned (no score)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)

    png = FIG_DIR / "wave_height" / f"rmse_bias_vs_wave_height_by_offset_scale_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Island-offset scale against significant wave height. (a) The alongshore "
        "offset BRIE receives under each scale for the dune-line source, each "
        "with its own mean removed: divided by ten (the units error every "
        "calibrated run carries), in metres as measured (the island's lean "
        "included), and in metres with the linear trend removed. (b) Interior "
        "RMSE, log scale, and (c) interior mean bias against Hs. Solid lines and "
        "filled markers: dune-line offset; dashed and open: shoreline-derived "
        "offset. Crosses at the foot, in the colour of the offset scale, mark Hs "
        "values where a run drowned and has no score (every scale at Hs 0.5 and 0.6 m; the trend-removed "
        "shoreline-offset run at 3.0 m). Wave asymmetry 0.7, high-angle fraction 0.1, Tp 8 s. "
        + caption_common()))
    return png


# =============================================================================
# 2. WAVE ANGLE
# =============================================================================

def fig_wave_angle(df, null):
    d = df[(df.sweep == "wave_angle") & df.scored]
    control = d[(d.offset_scale == "div10")
                & (d.wave_asymmetry == study.DEFAULT_ASYMMETRY)
                & (d.wave_angle_high_fraction == study.DEFAULT_HIGH_FRACTION)]
    asyms = sorted(d.wave_asymmetry.unique())
    ramp = dict(zip(asyms, SMOOTH_RAMP))

    fig, axes = plt.subplots(2, 3, figsize=figsize("double", height=5.2),
                             sharex=True, sharey="row", constrained_layout=True)
    for j, scale in enumerate(SCALES):
        g = d[d.offset_scale == scale]
        for i, (col, ylab) in enumerate((("rmse_interior_m_yr", "Interior RMSE (m/yr)"),
                                         ("mean_bias_interior_m_yr", "Interior mean bias (m/yr)"))):
            ax = axes[i, j]
            for a in asyms:
                h = g[g.wave_asymmetry == a].sort_values("wave_angle_high_fraction")
                ax.plot(h.wave_angle_high_fraction, h[col], marker="o", ms=3.5,
                        lw=1.3, color=ramp[a])
            if not control.empty:
                ax.axhline(control[col].iloc[0], color=C["BASE"], lw=0.9, ls=":")
            if i == 0:
                ax.axhline(null, **NULL_STYLE)
            if i == 1:
                ax.axhline(0, color=INK_MUTED, lw=0.6)
                ax.set_xlabel(HIGH_ANGLE_LABEL)
            if j == 0:
                ax.set_ylabel(ylab)
            ax.set_xticks(sorted(d.wave_angle_high_fraction.unique()))
            ax.grid(axis="y", which="both")
            open_frame(ax)
            if i == 0:
                _title(ax, j, PANEL_TITLE[scale])
    # linear: the whole sweep spans 1.0-2.5 m/yr, too narrow for log ticks
    axes[0, 0].set_ylim(bottom=0.9)

    handles = [Line2D([], [], color=ramp[a], marker="o", ms=3.5, lw=1.3,
                      label=f"{a:g}") for a in asyms]
    handles.append(Line2D([], [], color=C["BASE"], ls=":", lw=0.9,
                          label="Calibrated settings, offset ÷ 10"))
    handles.append(Line2D([], [], label=NULL_LABEL, **NULL_STYLE))
    fig.legend(handles=handles, title="Wave asymmetry", loc="outside lower center",
               ncol=4, frameon=False)

    png = FIG_DIR / "wave_angle" / f"rmse_bias_vs_high_angle_fraction_by_asymmetry_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Wave-angle tuning for each island-offset scale at Hs 1.0 m, Tp 8 s, "
        "dune-line offset. Top row: interior RMSE (log scale, shared); bottom "
        "row: interior mean bias. Columns: offset divided by ten (the "
        "calibrated runs' units error), in metres, and in metres with the "
        "linear trend removed. Line shade is wave asymmetry, light (0.5, "
        "waves balanced between the two directions) to dark (0.8); the dotted "
        "grey line is the offset-÷-10 run at the calibrated wave angles "
        "(asymmetry 0.7, high-angle fraction 0.1); the dashed black line is the "
        "RMSE of a flat line at the observed mean rate, i.e. of a model with no "
        "alongshore pattern at all. " + caption_common()))
    return png


# =============================================================================
# 3. BEST ACHIEVABLE
# =============================================================================

def best_label(row):
    if row.sweep == "wave_height":
        return f"Hs {row.Hs_m:g} m"
    return f"asymmetry {row.wave_asymmetry:g}, high-angle {row.wave_angle_high_fraction:g}"


def fig_best(df, null):
    d = df[df.scored & (df.offset_source == "duneline")]
    control_best = d[d.offset_scale == "div10"].rmse_interior_m_yr.min()
    rows = [(scale, sweep) for scale in SCALES for sweep in ("wave_height", "wave_angle")]

    fig, ax = plt.subplots(figsize=figsize("double", height=3.4), constrained_layout=True)
    table = []
    for k, (scale, sweep) in enumerate(rows):
        y = len(rows) - 1 - k - (len(SCALES) - 1 - SCALES.index(scale)) * 0.0
        g = d[(d.offset_scale == scale) & (d.sweep == sweep)]
        if g.empty:
            continue
        colr = SCALE_STYLE[scale]["color"]
        lo, hi = g.rmse_interior_m_yr.min(), g.rmse_interior_m_yr.max()
        best = g.loc[g.rmse_interior_m_yr.idxmin()]
        ax.plot([lo, hi], [y, y], color=colr, lw=2.2, alpha=0.45, solid_capstyle="butt")
        ax.plot(lo, y, marker=SWEEP_MARKER[sweep], ms=7, color=colr)
        ax.annotate(f"{lo:.2f}  ({best_label(best)})", (hi, y), xytext=(6, 0),
                    textcoords="offset points", va="center", fontsize=7.5, color=INK)
        table.append(dict(offset_scale=scale, sweep=sweep, n_runs=len(g),
                          rmse_best=lo, rmse_worst=hi,
                          bias_at_best=best.mean_bias_interior_m_yr,
                          best_settings=best_label(best), best_run=best.run_dir))
    ax.axvline(control_best, color=C["BASE"], lw=0.9, ls=":")
    ax.axvline(null, **NULL_STYLE)
    ax.set_yticks([len(rows) - 1 - k for k in range(len(rows))])
    ax.set_yticklabels([f"{SCALE_STYLE[s]['label']} — {SWEEP_LABEL[w].lower()}"
                        for s, w in rows])
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([1, 2, 5, 10, 20]))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel("Interior RMSE (m/yr), log scale")
    ax.grid(axis="x", which="both")
    open_frame(ax)
    ax.set_xlim(right=ax.get_xlim()[1] * 4)   # room for the labels

    handles = [Line2D([], [], color=INK_MUTED, marker=SWEEP_MARKER[w], ls="none", ms=6,
                      label=f"Best run, {SWEEP_LABEL[w].lower()}") for w in SWEEP_MARKER]
    handles.append(Line2D([], [], color=C["BASE"], ls=":", lw=0.9,
                          label="Best run, offset ÷ 10"))
    handles.append(Line2D([], [], label=NULL_LABEL, **NULL_STYLE))
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)

    png = FIG_DIR / "combined" / f"best_rmse_by_offset_scale_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(table).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The error each island-offset scale reaches when the wave climate is "
        "tuned, dune-line offset. Each bar spans the interior RMSE from the "
        "best to the worst scored run of one sweep: wave height tuned (Hs "
        "0.75-3.0 m at the calibrated wave angles; circles) or wave angles "
        "tuned (asymmetry 0.5-0.8 by high-angle fraction 0.1-0.4, and to 0.5 for "
        "metres, at Hs 1.0 m; "
        "triangles). The marker and label give the best run and its settings; "
        "the dotted line is the best offset-÷-10 run of either sweep, and the "
        "dashed line the RMSE of a flat line at the observed mean rate (a "
        "model with no alongshore pattern). Drowned "
        "runs are not in the bars. " + caption_common()))
    return png


# =============================================================================
# 4. BIAS AGAINST RMSE
# =============================================================================

def fig_scatter(df):
    d = df[df.scored]
    fig, ax = plt.subplots(figsize=figsize("single", height=3.6), constrained_layout=True)
    for (scale, sweep, source), g in d.groupby(["offset_scale", "sweep", "offset_source"]):
        colr = SCALE_STYLE[scale]["color"]
        ax.scatter(g.mean_bias_interior_m_yr, g.rmse_interior_m_yr, s=16,
                   marker=SWEEP_MARKER[sweep], linewidths=0.8, edgecolors=colr,
                   facecolors=colr if SOURCE_STYLE[source]["filled"] else "white",
                   zorder=3)
    ax.axvline(0, color=INK_MUTED, lw=0.6)
    log_rmse_axis(ax)
    ax.set_xlabel("Interior mean bias (m/yr)")
    ax.set_ylabel("Interior RMSE (m/yr)")
    ax.grid(which="both")
    open_frame(ax)
    handles = scale_handles() + [
        Line2D([], [], color=INK_MUTED, marker=SWEEP_MARKER[w], ls="none", ms=5,
               label=SWEEP_LABEL[w]) for w in SWEEP_MARKER] + [
        Line2D([], [], color=INK_MUTED, marker="o", ls="none", ms=5, mfc="white",
               label="Shoreline offset (open)")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False,
               fontsize=7)

    png = FIG_DIR / "combined" / f"bias_vs_rmse_all_runs_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Every scored run of both sweeps, interior mean bias against interior "
        "RMSE (log scale). Colour: island-offset scale. Circles: wave height "
        "tuned; triangles: wave angles tuned. Filled: dune-line offset; open: "
        "shoreline-derived offset (wave-height sweep only). A good run sits "
        "low and on the zero-bias line. " + caption_common()))
    return png


# =============================================================================
# 5. ALONGSHORE, BEST RUNS
# =============================================================================

def fig_alongshore(df, target):
    d = df[df.scored & (df.offset_source == "duneline")]
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6), constrained_layout=True)
    gis = target.index.values
    ax.plot(gis, target.values, color=C_1997, lw=2.0, label="CoastSat LRR (LOESS, 10 domains)",
            zorder=4)
    rows = []
    for scale in SCALES:
        best = d[d.offset_scale == scale].sort_values("rmse_interior_m_yr").iloc[0]
        rates = run_rates(best)
        # the runner's RMSE, reproduced from the drawn curves: a check that
        # this figure shows what was scored
        rmse = float(np.sqrt(((interior(rates) - interior(target)) ** 2).mean()))
        if not np.isclose(rmse, best.rmse_interior_m_yr, rtol=1e-3):
            raise ValueError(f"{best.run_dir}: RMSE from the drawn curves {rmse:.4f} "
                             f"!= runner's {best.rmse_interior_m_yr:.4f}")
        ax.plot(rates.index, rates.values, color=SCALE_STYLE[scale]["color"], lw=1.3,
                label=f"{SCALE_STYLE[scale]['label']}: {best_label(best)}", zorder=5)
        rows.append(dict(offset_scale=scale, run_dir=best.run_dir,
                         settings=best_label(best), rmse=best.rmse_interior_m_yr))
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(gis[0], gis[-1])
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("Shoreline change rate, LRR (m/yr)")
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax)
    fig.legend(loc="outside lower center", ncol=2, frameon=False, fontsize=7.5)
    structures(ax)

    png = FIG_DIR / "combined" / f"alongshore_rate_best_runs_vs_coastsat_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "Where the error is: the modelled shoreline-change rate (LRR) along the "
        "island for the lowest-RMSE run of each island-offset scale, across "
        "both sweeps, dune-line offset, against the CoastSat target (blue). "
        "The settings of each run are in the legend. Positive is accretion. "
        "Villages shaded; the Buxton groin and the two piers are marked. "
        + caption_common()))
    return png


# =============================================================================
# 6. ALONGSHORE VARIATION EXPLAINED
# =============================================================================

VE_FLOOR = -1.0   # the axis stops here; a worse run is marked at the edge


def fig_variance_explained(df):
    d = df[df.scored & (df.offset_source == "duneline")]
    rows = [(scale, sweep) for scale in SCALES for sweep in ("wave_height", "wave_angle")]
    fig, ax = plt.subplots(figsize=figsize("double", height=3.4), constrained_layout=True)
    table = []
    for k, (scale, sweep) in enumerate(rows):
        y = len(rows) - 1 - k
        g = d[(d.offset_scale == scale) & (d.sweep == sweep)]
        if g.empty:
            continue
        colr = SCALE_STYLE[scale]["color"]
        best = g.loc[g.variance_explained.idxmax()]
        hi, lo = best.variance_explained, g.variance_explained.min()
        ax.plot([max(lo, VE_FLOOR), hi], [y, y], color=colr, lw=2.2, alpha=0.45,
                solid_capstyle="butt")
        if lo < VE_FLOOR:
            ax.plot(VE_FLOOR, y, marker="<", ms=5, color=colr, clip_on=False)
        # a best run below the floor is drawn at the floor, its value in the label
        ax.plot(max(hi, VE_FLOOR), y, marker=SWEEP_MARKER[sweep], ms=7, color=colr,
                clip_on=False)
        ax.annotate(f"{100 * hi:.0f}%  ({best_label(best)})", (max(hi, VE_FLOOR), y),
                    xytext=(9, 0), textcoords="offset points", va="center",
                    fontsize=7.5, color=INK)
        table.append(dict(offset_scale=scale, sweep=sweep, n_runs=len(g),
                          variance_explained_best=hi, variance_explained_worst=lo,
                          pattern_variance_explained_at_best=best.pattern_variance_explained,
                          rmse_at_best=best.rmse_interior_m_yr,
                          bias_at_best=best.mean_bias_interior_m_yr,
                          best_settings=best_label(best), best_run=best.run_dir))
    ax.axvline(0, **NULL_STYLE)
    ax.set_yticks([len(rows) - 1 - k for k in range(len(rows))])
    ax.set_yticklabels([f"{SCALE_STYLE[s]['label']} — {SWEEP_LABEL[w].lower()}"
                        for s, w in rows])
    ax.set_xlim(VE_FLOOR, 1.0)
    ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax.set_xlabel("Share of the observed alongshore variation explained")
    ax.grid(axis="x")
    open_frame(ax)
    handles = [Line2D([], [], color=INK_MUTED, marker=SWEEP_MARKER[w], ls="none", ms=6,
                      label=f"Best run, {SWEEP_LABEL[w].lower()}") for w in SWEEP_MARKER]
    handles.append(Line2D([], [], label=NULL_LABEL + " (0%)", **NULL_STYLE))
    handles.append(Line2D([], [], color=INK_MUTED, marker="<", ls="none", ms=5,
                          label="Worst run beyond −100%"))
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)

    png = FIG_DIR / "combined" / f"variance_explained_by_offset_scale_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(table).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "How much of the observed alongshore variation in shoreline-change rate "
        "each island-offset scale explains when the wave climate is tuned, "
        "dune-line offset. The score is 1 - sum((model - observed)^2) / "
        "sum((observed - observed mean)^2) over GIS 2-89: 100% is a perfect "
        "match, 0% (dashed line) is no better than a flat line at the observed "
        "mean rate, and below 0% is worse than that line. Bias counts against "
        "the score. Each bar spans the best to the worst scored run of one "
        "sweep (circles: wave height tuned; triangles: wave angles tuned); the "
        "marker and label give the best run. Bars are cut at -100%; an arrow "
        "marks a sweep whose worst run is lower. " + caption_common()))
    return png


# =============================================================================
# 7. SPREAD AGAINST PLACEMENT
# =============================================================================

def fig_spread_vs_placement(df):
    d = df[df.scored & (df.offset_source == "duneline")]
    fig, ax = plt.subplots(figsize=figsize("single", height=3.6), constrained_layout=True)
    # Lines of equal pattern skill (bias removed): skill = 2 r s - s^2 with s
    # the sd ratio, so r = (skill + s^2) / (2 s).
    s = np.geomspace(0.1, 20, 400)
    for k, style in ((0.0, dict(NULL_STYLE)),
                     (0.2, dict(color=INK_MUTED, lw=0.6, ls=":")),
                     (0.4, dict(color=INK_MUTED, lw=0.6, ls=":"))):
        r = (k + s ** 2) / (2 * s)
        ok = r <= 1.0
        ax.plot(s[ok], r[ok], **style, zorder=1)
        if k > 0:
            # at the curve's lowest point, s = r = sqrt(k)
            ax.annotate(f"{100 * k:.0f}%", (np.sqrt(k), np.sqrt(k)), xytext=(0, -9),
                        textcoords="offset points", ha="center", fontsize=6.5,
                        color=INK_MUTED)
    for (scale, sweep), g in d.groupby(["offset_scale", "sweep"]):
        colr = SCALE_STYLE[scale]["color"]
        ax.scatter(g.sd_ratio, g.r_alongshore, s=16, marker=SWEEP_MARKER[sweep],
                   color=colr, edgecolors=colr, linewidths=0.6, zorder=3)
    ax.axvline(1.0, color=INK_MUTED, lw=0.6)
    ax.plot(1.0, 1.0, marker="*", ms=10, color=C_1997, zorder=4, clip_on=False)
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([0.2, 0.5, 1, 2, 5, 10]))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(0.2, 15)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Modelled ÷ observed alongshore spread (sd)")
    ax.set_ylabel("Correlation with CoastSat, r")
    ax.grid(which="major")
    open_frame(ax)
    handles = scale_handles() + [
        Line2D([], [], color=INK_MUTED, marker=SWEEP_MARKER[w], ls="none", ms=5,
               label=SWEEP_LABEL[w]) for w in SWEEP_MARKER] + [
        Line2D([], [], label="Pattern no better than a flat line", **NULL_STYLE),
        Line2D([], [], color=C_1997, marker="*", ls="none", ms=9, label="Perfect match")]
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False,
               fontsize=7)

    png = FIG_DIR / "combined" / f"spread_vs_correlation_all_runs_{WINDOW}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "Why a run explains little of the alongshore variation: too little "
        "spread, or spread in the wrong places. Every scored dune-line run of "
        "both sweeps (circles: wave height tuned; triangles: wave angles "
        "tuned). x: the standard deviation of the modelled rate along GIS 2-89 "
        "divided by the observed one (log scale; 1 = the observed amount of "
        "variation). y: the correlation of the two alongshore series. The "
        "curves are lines of equal pattern skill with bias removed, 2rs - s^2: "
        "dashed, 0% (no better than a flat line); dotted, 20% and 40%. The star "
        "is a perfect match. " + caption_common()))
    return png


# =============================================================================
# 8. ALONGSHORE SENSITIVITY: RATE AND POSITION CHANGE, ONE PARAMETER AT A TIME
# =============================================================================

# Observed position change over the window, from 5-scr: mean CoastSat position
# over calendar 2010 minus calendar 1996, seaward positive, LOESS-smoothed at
# 10 domains to match the rate target's window.
OBSERVED_CHANGE = (INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change"
                   / WINDOW / "smoothed" / "tables" / "domain_smoothed.csv")
RUN_YEARS = 14
SCALE_CMAP = {"div10": "Greys", "metres": "Oranges", "metres-detrended": "Purples"}

# One parameter moved, the other two at the calibration defaults
# (Hs 2.5 is not a default here: the wave-angle sweep runs at Hs 1.0).
PARAMETERS = {
    "wave_height": dict(sweep="wave_height", column="Hs_m",
                        fixed={"wave_asymmetry": 0.7, "wave_angle_high_fraction": 0.1},
                        label="Hs (m)", fixed_text="asymmetry 0.7, high-angle fraction 0.1"),
    "high_angle_fraction": dict(sweep="wave_angle", column="wave_angle_high_fraction",
                                fixed={"Hs_m": 1.0, "wave_asymmetry": 0.7},
                                label="High-angle fraction",
                                fixed_text="Hs 1.0 m, asymmetry 0.7"),
    "asymmetry": dict(sweep="wave_angle", column="wave_asymmetry",
                      fixed={"Hs_m": 1.0, "wave_angle_high_fraction": 0.1},
                      label="Asymmetry", fixed_text="Hs 1.0 m, high-angle fraction 0.1"),
}


def observed_change():
    t = pd.read_csv(OBSERVED_CHANGE)
    return t[t.window_domains == 10].set_index("domain_number")["observed_m"]


def run_table(row):
    return pd.read_csv(study.STUDY_DIR / row.run_dir / "tables" / "shoreline_change_rate.csv"
                       ).set_index("gis_domain")


def fig_alongshore_sensitivity(df, target, obs_change, scale, param):
    spec = PARAMETERS[param]
    g = df[(df.sweep == spec["sweep"]) & (df.offset_scale == scale)
           & (df.offset_source == "duneline")]
    for col, val in spec["fixed"].items():
        g = g[np.isclose(g[col], val)]
    g = g.sort_values(spec["column"])
    values = g[spec["column"]].tolist()
    cmap = plt.get_cmap(SCALE_CMAP[scale])
    shade = dict(zip(values, cmap(np.linspace(0.35, 0.95, len(values)))))

    fig, (ax_r, ax_p) = plt.subplots(2, 1, figsize=figsize("double", height=5.6),
                                     sharex=True, constrained_layout=True)
    gis = target.index.values
    ax_r.plot(gis, target.values, color=C_1997, lw=2.2, zorder=6)
    ax_p.plot(obs_change.index, obs_change.values, color=C_1997, lw=2.2, zorder=6)
    handles = [Line2D([], [], color=C_1997, lw=2.2, label="CoastSat (LOESS, 10 domains)")]
    table = []
    for _, row in g.iterrows():
        v = row[spec["column"]]
        if row.scored:
            t = run_table(row)
            ax_r.plot(t.index, t.lrr_m_yr, color=shade[v], lw=1.2, zorder=4)
            ax_p.plot(t.index, t.change_rate_m_yr * RUN_YEARS, color=shade[v], lw=1.2, zorder=4)
            lab = (f"{v:g}  ({100 * row.variance_explained:.0f}% of rate variation, "
                   f"bias {row.mean_bias_interior_m_yr:+.2f} m/yr)")
            table.append(dict(value=v, run_dir=row.run_dir,
                              variance_explained=row.variance_explained,
                              bias_m_yr=row.mean_bias_interior_m_yr))
        else:
            lab = f"{v:g}  (barrier drowned, no output)"
        handles.append(Line2D([], [], color=shade[v], lw=1.6 if row.scored else 0,
                              marker=None if row.scored else "x", label=lab))
    for ax in (ax_r, ax_p):
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(gis[0], gis[-1])
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(ax is ax_r))
    ax_r.set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    ax_p.set_ylabel("Shoreline position change,\n2010 minus 1996 (m)")
    ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax_r, 0, "Rate")
    _title(ax_p, 1, "Position change")
    fig.legend(handles=handles, title=f"{SCALE_STYLE[scale]['label']}: {spec['label']}",
               loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    structures(ax_p, label=True)
    structures(ax_r, label=False)

    png = (FIG_DIR / "alongshore_sensitivity" / scale
           / f"rate_and_position_change_by_{param}_{scale}_{WINDOW}.png")
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(table).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        f"Alongshore response to {spec['label'].lower()} with the island offset "
        f"{SCALE_STYLE[scale]['label'].lower()} (dune-line source), the other "
        f"wave parameters held at {spec['fixed_text']}, Tp 8 s. Light to dark: "
        f"increasing {spec['label'].lower()}; the legend gives each run's share "
        "of the observed alongshore rate variation explained (GIS 2-89) and "
        "its interior mean bias. (a) Modelled LRR shoreline-change rate "
        "against the CoastSat LRR target. (b) Modelled shoreline position "
        "change, end of 2010 minus start of 1996, against the observed CoastSat "
        "change, the mean position over calendar 2010 minus that over calendar "
        "1996 (5-scr/3-rates/coastsat/total_change/1996_2010, smoothed at 10 "
        "domains). Seaward positive. Villages shaded; the Buxton groin and the "
        "two piers marked. " + caption_common()))
    return png


def alongshore_sensitivity_figures(df, target):
    obs = observed_change()
    return [fig_alongshore_sensitivity(df, target, obs, scale, param)
            for scale in SCALES for param in PARAMETERS]


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    df = load_runs()
    target = coastsat_target()
    null = flat_line_rmse(target)
    print(f"flat-line RMSE {null:.3f}")
    for f in (fig_wave_height(df), fig_wave_angle(df, null), fig_best(df, null),
              fig_scatter(df), fig_alongshore(df, target),
              fig_variance_explained(df), fig_spread_vs_placement(df),
              *alongshore_sensitivity_figures(df, target)):
        print(f.relative_to(study.STUDY_DIR))


if __name__ == "__main__":
    main()
