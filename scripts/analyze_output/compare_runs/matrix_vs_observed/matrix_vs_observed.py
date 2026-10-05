#!/usr/bin/env python3
"""
Every matrix run against the observations: rate and position change, and start/end positions.

    python scripts/analyze_output/compare_runs/matrix_vs_observed/matrix_vs_observed.py

Two figures per run, written into each run's figures/vs_observed/ and to
output/comparisons/matrix_vs_observed/, plus per-scenario figures, scores.csv
and y_bounds.txt. Other scripts import its loaders (matrix_runs, rates, score).
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-03
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

_REPO = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "hatteras_ms"))
sys.path.insert(0, str(_REPO / "scripts" / "hatteras_ms" / "experiments"))

import HAT_metres_1_offset_units as common  # noqa: E402
from cascade_pipeline.hindcast import build_shoreline_target  # noqa: E402
from cascade_pipeline.run_registry import find_run_dir, load_run_index  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,
    figsize, mark_offaxis, open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_topo_version import (  # noqa: E402
    INIT_ROOT, RAW_OFFSET_DIR, dune_line_for_year)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_PERIODS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = _REPO / "output" / "raw_runs"
OUT = COMPARISONS_ROOT / "matrix_vs_observed"
RUN_SUBDIR = Path("figures") / "vs_observed"
PRESETS = ("edgeBE", "zeroBE")
RATE_YLIM = (-7.5, 7.5)
POSITION_YLIM = None        # set in main()
ENDS_YLIM = None            # set in main()

OBSERVED = dict(color=INK, lw=2.0)
MODEL_ONE = dict(color="#2166ac", lw=1.4)
# start/end figure: shoreline blue and dune red for the observations, the model black
C_START = INK_MUTED
C_MODEL_END = INK
C_COASTSAT = "#2166ac"
C_DUNE = "#b2182b"

SCEN_ORDER = ("natural", "beachdune_only", "roadway_only", "full_no_fill", "full_management")
SCEN_STYLE = {
    "natural":         dict(color="#1b7f6b", lw=1.3, ls="-"),
    "beachdune_only":  dict(color="#7b3294", lw=1.3, ls="-"),
    "roadway_only":    dict(color="#e08214", lw=1.3, ls="-"),
    "full_no_fill":    dict(color="#2166ac", lw=1.3, ls=(0, (4, 2))),
    "full_management": dict(color="#2166ac", lw=1.6, ls="-"),
}
SCEN_LABEL = {"natural": "Natural", "beachdune_only": "Beach and dune management only",
              "roadway_only": "Road management only", "full_no_fill": "Full management, no fills",
              "full_management": "Full management"}
PRESET_TEXT = {"edgeBE": "end rates solved on CoastSat (edgeBE)",
               "zeroBE": "no imposed end rates (zeroBE)"}
OPTION_A = ("Option A matrix on the adopted model, the 1996-2015 and 2010-2026 windows (re-run "
            "2026-10-03): Barrier3D hatteras/adopted (the three overwash fixes, per-cell dune "
            "ceilings), storms v3_split12_trim24 (events split at 12 h gaps, 24 h around each "
            "peak), edgeBE ends solved on that setup against each window's own 7-domain target; "
            "island offset in metres, Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle fraction 0.5, no groin.")
RATE_CAPTION = (
    OPTION_A + " (a) Modelled OLS shoreline-change rate against the CoastSat LRR scoring "
    "target (7-domain LOWESS, raw means GIS 1-10). (b) Modelled position change over the "
    "window (endpoint rate x the window's years) against the observed CoastSat change, mean position over "
    "the last calendar year minus the first, smoothed at 7 domains. The model lines are "
    "smoothed the same way (7-domain LOWESS, GIS 1-10 raw). Seaward positive; "
    "scores over the interior GIS 2-89. A value beyond the axis (GIS 1, where both the "
    "observation and the held end domain run off it) is marked with a triangle at the edge.")
SMOOTH = 7                  # smoothing width in domains, target and observed (10 until 2026-09-28)
# -----------------------------------------------------------------------------


# The window's end year, from HATTERAS_PERIODS
def end_year(p):
    return int(HATTERAS_PERIODS[int(p)]["end_year"])


# The window's length in years
def years(p):
    return end_year(p) - int(p)


# '1996_2015' from a start year
def window(p):
    return f"{p}_{end_year(p)}"


# The smoothed total-change table for a window
def _coastsat_change(period):
    f = (INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change" / window(period)
         / "smoothed" / "tables" / "domain_smoothed.csv")
    return pd.read_csv(f)


# CoastSat position change per domain, seaward +; window 0 is the unsmoothed mean
def observed_change(period, window_domains=SMOOTH):
    d = _coastsat_change(period)
    return d[d.window_domains == window_domains].set_index("domain_number")["observed_m"]


# The option A matrix runs from the run index, with their folders
def matrix_runs():
    idx = load_run_index(RAW_RUNS / "run_index.csv")
    m = idx[(idx["kind"] == "matrix") & idx["run_name"].str.contains("offsetmetres")].copy()
    m["period"] = m["start_year"].astype(int)
    # The configured windows only; the index also holds the retired 14-yr ones
    m = m[m["end_year"].astype(int) == m["period"].map(end_year)].copy()
    # The index stores the flag as text, and "False" is truthy: parse it.
    m["relocations_enabled"] = m["relocations_enabled"].astype(str).str.lower().eq("true")
    m["run_dir"] = [find_run_dir(RAW_RUNS, r.run_name, window(r.period), r.source_sink_preset)
                    for r in m.itertuples()]
    return m


# A run's per-domain rate table
def rates(run_dir):
    return pd.read_csv(run_dir / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")


# A model series smoothed as the observation is (LOWESS 7 domains, GIS 1-10 raw)
def smoothed(series):
    return common.smooth_like_target(series)


# Last annual shoreline minus the first, seaward + (x_s grows landward, so the sign flips)
def model_end_change(run_dir):
    m = np.load(next(run_dir.glob("*_shoreline_matrix.npy")))
    D = HATTERAS_DOMAINS
    return pd.Series(-(m[-1] - m[0])[D.start_real_index:D.end_real_index],
                     index=np.arange(D.first_gis_id, D.last_gis_id + 1))


# The runner's end-year dune-line target as a change, seaward +
def dune_end_change(run_dir, period):
    m = np.load(next(run_dir.glob("*_shoreline_matrix.npy")))
    _, change = build_shoreline_target(m[0], period, end_year(period), HATTERAS_DOMAINS,
                                       RAW_OFFSET_DIR)
    if change is None:     # no dune line for the end year
        return None
    D = HATTERAS_DOMAINS
    return pd.Series(-np.asarray(change, float),
                     index=np.arange(D.first_gis_id, D.last_gis_id + 1))


# Interior bias and RMSE of model minus observation
def score(model, obs):
    r = (common.interior(model) - common.interior(obs)).dropna()
    return float(r.mean()), float(np.sqrt((r ** 2).mean()))


# Symmetric y limits from the widest interior value (GIS 2-89), rounded up to 10
def _sym(vals):
    vals = [common.interior(v) for v in vals]
    half = float(np.ceil(np.nanmax(np.abs(np.concatenate([np.asarray(v, float) for v in vals])))
                         / 10.0) * 10.0)
    return (-half, half)


# Triangles at the axis edge for a model line beyond it (the held end domains)
def _mark(ax, series, ylim, color):
    mark_offaxis(ax, series.index.to_numpy(float), series.to_numpy(float), ylim[1], color=color)


# Figures

# Shared axis: zero line, GIS 1-90, grid, village bands
def _axis(ax, label_towns):
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1, 90)
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=label_towns)


# Limits, labels and titles for the rate and position panels
def frame(axes, period):
    ax_r, ax_p = axes
    for ax in axes:
        _axis(ax, ax is ax_r)
    ax_r.set_ylim(*RATE_YLIM)
    ax_p.set_ylim(*POSITION_YLIM)
    ax_r.set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    ax_p.set_ylabel(f"Shoreline position change,\n{end_year(period)} minus {period} (m)")
    ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax_r, 0, "Rate")
    _title(ax_p, 1, "Position change")
    structures(ax_p, label=True)
    structures(ax_r, label=False)


# Readable run description for legends and captions
def run_title(r):
    reloc = ", relocations on" if r.relocations_enabled else ""
    return f"{SCEN_LABEL[r.scenario]}{reloc}, {r.source_sink_preset}, {window(r.period).replace('_', '–')}"


# File-name stem for a run's figures
def stem(r):
    return f"{r.source_sink_preset}_{r.scenario}{'_reloc' if r.relocations_enabled else ''}_{window(r.period)}"


# Rate and position change for one run, into the run and the comparison tree
def fig_rate_position(r, target, obs):
    rt = rates(r.run_dir)
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.2), sharex=True,
                             constrained_layout=True)
    ax_r, ax_p = axes
    ax_r.plot(target.index, target.values, zorder=6, **OBSERVED)
    ax_p.plot(obs.index, obs.values, zorder=6, **OBSERVED)
    rate_m, pos_m = smoothed(rt.lrr_m_yr), smoothed(rt.change_rate_m_yr * years(r.period))
    ax_r.plot(rate_m.index, rate_m.values, zorder=5, **MODEL_ONE)
    ax_p.plot(pos_m.index, pos_m.values, zorder=5, **MODEL_ONE)
    frame(axes, r.period)
    _mark(ax_r, rate_m, RATE_YLIM, MODEL_ONE["color"])
    _mark(ax_p, pos_m, POSITION_YLIM, MODEL_ONE["color"])
    _mark(ax_r, target, RATE_YLIM, OBSERVED["color"])
    _mark(ax_p, obs, POSITION_YLIM, OBSERVED["color"])
    br, er = score(rate_m, target)
    bp, ep = score(pos_m, obs)
    fig.legend(handles=[
        Line2D([], [], label="CoastSat (LOWESS, 7 domains)", **OBSERVED),
        Line2D([], [], label=f"Model, smoothed alike: rate bias {br:+.2f} m/yr, RMSE {er:.2f}; "
                             f"position bias {bp:+.1f} m, RMSE {ep:.1f}", **MODEL_ONE)],
        title=run_title(r), loc="outside lower center", ncol=1, frameon=False, fontsize=7)
    caption = f"{run_title(r)}. Run {r.run_name}. " + RATE_CAPTION
    pngs = [r.run_dir / RUN_SUBDIR / "rate_and_position_change.png",
            OUT / "rate_and_position_change" / window(r.period) / r.source_sink_preset
            / f"rate_and_position_{stem(r)}.png"]
    for png in pngs:
        save(fig, png, dpi=300, close=False)
        record_caption(png, caption)
    plt.close(fig)
    return dict(rate_bias_m_yr=br, rate_rmse_m_yr=er,
                position_bias_m=bp, position_rmse_m=ep), pngs[1]


# Start and end positions relative to the start line, for one run
def fig_start_end(r, cs_raw, cs_smooth, dune):
    model = model_end_change(r.run_dir)
    s, e = r.period, end_year(r.period)
    # No dune line for the end year: the red line and its scores are left out
    has_dune = dune is not None
    dv0, dv1 = (dune_line_for_year(s), dune_line_for_year(e)) if has_dune else (None, None)
    fig, ax = plt.subplots(figsize=figsize("double", height=3.6), constrained_layout=True)
    _axis(ax, True)
    ax.axhline(0, color=C_START, lw=2.2, zorder=3)
    ax.plot(cs_smooth.index, cs_smooth.values, color=C_COASTSAT, lw=1.0, alpha=0.45, zorder=4)
    ax.plot(cs_raw.index, cs_raw.values, color=C_COASTSAT, lw=1.3, zorder=5,
            marker="o", ms=2.2)
    if has_dune:
        ax.plot(dune.index, dune.values, color=C_DUNE, lw=1.3, zorder=5, marker="s", ms=2.0)
    ax.plot(model.index, model.values, color=C_MODEL_END, lw=1.8, zorder=6)
    ax.set_ylim(*ENDS_YLIM)
    _mark(ax, model, ENDS_YLIM, C_MODEL_END)
    ax.set_ylabel(f"Position relative to the\n{s} start line (m, seaward +)")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    structures(ax, label=True)
    bc, ec = score(model, cs_raw)
    bd, ed = score(model, dune) if has_dune else (np.nan, np.nan)
    handles = [
        Line2D([], [], color=C_START, lw=2.2, label=f"Start position, {s} (model year 0)"),
        Line2D([], [], color=C_MODEL_END, lw=1.8, label=f"Modelled end position, {e}"),
        Line2D([], [], color=C_COASTSAT, lw=1.3, marker="o", ms=2.2,
               label=f"CoastSat end, {e} (domain means; faint: 7-domain LOWESS); "
                     f"model bias {bc:+.1f} m, RMSE {ec:.1f}"),
    ] + ([Line2D([], [], color=C_DUNE, lw=1.3, marker="s", ms=2.0,
                 label=f"Dune-line end ({dv0} → {dv1} change); model bias {bd:+.1f} m, "
                       f"RMSE {ed:.1f}")] if has_dune else [])
    fig.legend(handles=handles, title=run_title(r), loc="outside lower center", ncol=1,
               frameon=False, fontsize=7)
    caption = (
        f"{run_title(r)}. Run {r.run_name}. Start and end shoreline positions, drawn relative "
        f"to the model's {s} start line (zero) because the island's own position varies by "
        "~6 km along the reach while the changes are tens of metres. Black: the modelled "
        f"{e} position. Blue: the observed CoastSat change added to the start (mean position "
        "over the last calendar year minus the first, per domain; faint line smoothed at 7 "
        + (f"domains). Red: the digitised dune-line change, the {dv0} line to the {dv1} line, "
           "the runner's own end-year target; its survey interval is not the calendar window "
           "and is not rescaled. " if has_dune else
           f"domains). No dune line was digitized for {e}, so none is drawn. ")
        + "Seaward positive; scores over the interior GIS 2-89. "
        + OPTION_A)
    pngs = [r.run_dir / RUN_SUBDIR / "start_and_end_positions.png",
            OUT / "start_and_end_positions" / window(r.period) / r.source_sink_preset
            / f"start_and_end_positions_{stem(r)}.png"]
    for png in pngs:
        save(fig, png, dpi=300, close=False)
        record_caption(png, caption)
    plt.close(fig)
    return dict(end_bias_vs_coastsat_m=bc, end_rmse_vs_coastsat_m=ec,
                end_bias_vs_duneline_m=bd, end_rmse_vs_duneline_m=ed), pngs[1]


# Every scenario of one window and preset on one figure
def fig_scenarios(runs, period, preset, target, obs):
    sub = runs[(runs.period == period) & (runs.source_sink_preset == preset)
               & ~runs.relocations_enabled]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8), sharex=True,
                             constrained_layout=True)
    ax_r, ax_p = axes
    ax_r.plot(target.index, target.values, zorder=6, **OBSERVED)
    ax_p.plot(obs.index, obs.values, zorder=6, **OBSERVED)
    handles = [Line2D([], [], label="CoastSat (LOWESS, 7 domains)", **OBSERVED)]
    drawn = []
    for scen in SCEN_ORDER:
        hit = sub[sub.scenario == scen]
        if hit.empty:
            continue
        rt = rates(hit.iloc[0].run_dir)
        rate_m = smoothed(rt.lrr_m_yr)
        pos_m = smoothed(rt.change_rate_m_yr * years(period))
        ax_r.plot(rate_m.index, rate_m.values, zorder=4, **SCEN_STYLE[scen])
        ax_p.plot(pos_m.index, pos_m.values, zorder=4, **SCEN_STYLE[scen])
        drawn.append((rate_m, pos_m, SCEN_STYLE[scen]["color"]))
        _, er = score(rate_m, target)
        bp, _ = score(pos_m, obs)
        handles.append(Line2D([], [], **SCEN_STYLE[scen],
                              label=f"{SCEN_LABEL[scen]}: rate RMSE {er:.2f} m/yr, "
                                    f"position bias {bp:+.1f} m"))
    frame(axes, period)
    for rate_m, pos_m, c in drawn + [(target, obs, OBSERVED["color"])]:
        _mark(ax_r, rate_m, RATE_YLIM, c)
        _mark(ax_p, pos_m, POSITION_YLIM, c)
    fig.legend(handles=handles, title=f"{PRESET_TEXT[preset]}, {window(period).replace('_', '–')}",
               loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = (OUT / "rate_and_position_change" / window(period)
           / f"scenarios_rate_and_position_{preset}_{window(period)}.png")
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"Every management scenario, {PRESET_TEXT[preset]}, {window(period).replace('_', '-')}, "
        "relocation arms left out (within 0.001 m/yr of their twins). " + RATE_CAPTION))
    return png


# Run: load every run and observation, fix the y axes, draw every figure, write scores
def main():
    global POSITION_YLIM, ENDS_YLIM
    apply_style()
    runs = matrix_runs()
    periods = sorted(runs.period.unique())
    cs = {p: (observed_change(p, 0), observed_change(p, SMOOTH)) for p in periods}
    dunes = {r.run_name: dune_end_change(r.run_dir, r.period) for r in runs.itertuples()}
    models = {r.run_name: model_end_change(r.run_dir) for r in runs.itertuples()}

    # One y range per panel type, across every run and observation
    POSITION_YLIM = _sym([cs[p][1] for p in periods]
                         + [rates(r.run_dir).change_rate_m_yr * years(r.period) for r in runs.itertuples()])
    ENDS_YLIM = _sym([cs[p][0] for p in periods] + [v for v in dunes.values() if v is not None]
                     + list(models.values()))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "y_bounds.txt").write_text(
        f"rate panel: {RATE_YLIM[0]:g} to +{RATE_YLIM[1]:g} m/yr, as the runner's real-domains rate "
        "figure (rerender_run_figures.py --ylim-real=-7.5,7.5); its with-buffers figure is +/-10\n"
        f"position-change panel: {POSITION_YLIM[0]:g} to +{POSITION_YLIM[1]:g} m = max |change| over "
        "the interior GIS 2-89 of every run and the smoothed CoastSat change, rounded up to "
        "10 m; values beyond (the held end domains) are marked at the edge\n"
        f"start/end positions: {ENDS_YLIM[0]:g} to +{ENDS_YLIM[1]:g} m = max |position relative to "
        "the start| over the interior GIS 2-89 of every run's end, the CoastSat domain means and "
        "the dune-line change where there is one, rounded up to 10 m\n", encoding="utf-8")

    # Per-run figures, then per-scenario figures, per window
    rows, written = [], []
    for period in periods:
        target = common.coastsat_target(period)
        cs_raw, cs_smooth = cs[period]
        for r in runs[runs.period == period].itertuples():
            a, pa = fig_rate_position(r, target, cs_smooth)
            b, pb = fig_start_end(r, cs_raw, cs_smooth, dunes[r.run_name])
            rows.append(dict(run_name=r.run_name, window=window(period),
                             preset=r.source_sink_preset, scenario=r.scenario,
                             relocations=bool(r.relocations_enabled), **a, **b))
            written += [pa, pb]
        for preset in PRESETS:
            written.append(fig_scenarios(runs, period, preset, target, cs_smooth))
    # Scores table and summary
    s = pd.DataFrame(rows).sort_values(["window", "preset", "scenario", "relocations"])
    s.to_csv(OUT / "scores.csv", index=False)
    print(s.drop(columns="run_name").to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
    print((OUT / "y_bounds.txt").read_text(encoding="utf-8"))
    print(f"{2 * len(rows)} figures written into the runs' {RUN_SUBDIR.as_posix()}/, "
          f"{len(written)} under {OUT}")


if __name__ == "__main__":
    main()
