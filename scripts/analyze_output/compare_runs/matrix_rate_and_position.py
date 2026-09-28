#!/usr/bin/env python3
"""
matrix_rate_and_position.py
==============================================================================
The option A matrix, rate AND position change, against CoastSat (Hannah,
2026-09-27: "Where are these position plots?"). The runner draws only the
rate figures into each run's figures/, and model_vs_observed/ is rates only,
so the matrix had no position-change figure. The layout is the one the wave
experiments use (HAT_metres_2_wave_sensitivity_plot.fig_alongshore):

    (a) rate      the model's OLS rate (lrr_m_yr) against the CoastSat LRR
                  target, 10-domain LOESS with the southern 10 raw, as scored
    (b) position  the model's position change over the window, endpoint rate
                  x 14 yr (last annual shoreline minus first), against the
                  observed CoastSat change: the mean position over the last
                  calendar year minus that over the first, smoothed at 10
                  domains (5-scr/3-rates/coastsat/total_change/<w>/smoothed)

Seaward positive. GIS 1-90, the interior GIS 2-89 is what the scores cover.

WRITES
    each matrix run's own figures/rate_and_position_change.png (22 runs; run
        folders are untracked, so these live beside the run)
    output/comparisons/matrix_rate_and_position/
        <window>/scenarios_rate_and_position_<preset>_<window>.png
            every scenario of one preset on one pair of panels
        <window>/by_scenario/rate_and_position_<preset>_<scenario>_<window>.png
            the same figure as the run's own, filed together
        scores.csv   bias and RMSE per run, rate (m/yr) and position (m)

USAGE
    python scripts/analyze_output/compare_runs/matrix_rate_and_position.py
==============================================================================
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
from cascade_pipeline.run_registry import find_run_dir, load_run_index  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, INK_MUTED, _title, apply_style,
    figsize, open_frame, record_caption, save, structures, town_bands)
from site_layer.hat_topo_version import INIT_ROOT  # noqa: E402

RAW_RUNS = _REPO / "output" / "raw_runs"
OUT = COMPARISONS_ROOT / "matrix_rate_and_position"
YEARS = 14
# ONE Y AXIS FOR EVERY FIGURE (Hannah, 2026-09-27: "all of the figures ... should
# share a consistent y axis so we can compare between them"). The rate panel is
# fixed (below); the runner's own rate figures in each run folder are fixed to
# match (rerender_run_figures.py --ylim=-10,10 --ylim-real=-7.5,7.5). The position
# panel is symmetric and set in main() from the widest run or observation,
# rounded up to 10 m, and written to y_bounds.txt.
# +/-7.5 since the same day, matching the runner's real-domains figure
# (--ylim-real): this panel draws the real reach and the LOESS target only.
RATE_YLIM = (-7.5, 7.5)
POSITION_YLIM = None
PRESETS = ("edgeBE", "zeroBE")
OBSERVED = dict(color=INK, lw=2.0)
MODEL_ONE = dict(color="#2166ac", lw=1.4)
# One colour per scenario, fixed, so a scenario reads the same in every figure.
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
COMMON_CAPTION = (
    "Option A matrix (2026-09-27): island offset in metres, Hs 2.0 m, Tp 7.5 s, asymmetry "
    "0.6, high-angle fraction 0.5, no groin. (a) Modelled OLS shoreline-change rate against "
    "the CoastSat LRR scoring target (10-domain LOESS, raw means GIS 1-10). (b) Modelled "
    "position change over the window (endpoint rate x 14 yr) against the observed CoastSat "
    "change, mean position over the last calendar year minus the first, smoothed at 10 "
    "domains. Seaward positive; scores over the interior GIS 2-89.")


def window(p):
    return f"{p}_{p + YEARS}"


def observed_change(period):
    f = (INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change" / window(period)
         / "smoothed" / "tables" / "domain_smoothed.csv")
    d = pd.read_csv(f)
    return d[d.window_domains == 10].set_index("domain_number")["observed_m"]


def matrix_runs():
    idx = load_run_index(RAW_RUNS / "run_index.csv")
    m = idx[(idx["kind"] == "matrix") & idx["run_name"].str.contains("offsetmetres")].copy()
    m["period"] = m["start_year"].astype(int)
    # The index stores the flag as text, and "False" is truthy: parse it.
    m["relocations_enabled"] = m["relocations_enabled"].astype(str).str.lower().eq("true")
    m["run_dir"] = [find_run_dir(RAW_RUNS, r.run_name, window(r.period), r.source_sink_preset)
                    for r in m.itertuples()]
    return m


def rates(run_dir):
    return pd.read_csv(run_dir / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")


def score(model, obs):
    r = (common.interior(model) - common.interior(obs)).dropna()
    return float(r.mean()), float(np.sqrt((r ** 2).mean()))


def frame(axes, period):
    ax_r, ax_p = axes
    for ax in axes:
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1, 90)
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(ax is ax_r))
    ax_r.set_ylim(*RATE_YLIM)
    ax_p.set_ylim(*POSITION_YLIM)
    ax_r.set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    ax_p.set_ylabel(f"Shoreline position change,\n{period + YEARS} minus {period} (m)")
    ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax_r, 0, "Rate")
    _title(ax_p, 1, "Position change")
    structures(ax_p, label=True)
    structures(ax_r, label=False)


def fig_one(r, target, obs, pngs):
    rt = rates(r.run_dir)
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.2), sharex=True,
                             constrained_layout=True)
    ax_r, ax_p = axes
    ax_r.plot(target.index, target.values, zorder=6, **OBSERVED)
    ax_p.plot(obs.index, obs.values, zorder=6, **OBSERVED)
    ax_r.plot(rt.index, rt.lrr_m_yr, zorder=5, **MODEL_ONE)
    ax_p.plot(rt.index, rt.change_rate_m_yr * YEARS, zorder=5, **MODEL_ONE)
    frame(axes, r.period)
    br, er = score(rt.lrr_m_yr, target)
    bp, ep = score(rt.change_rate_m_yr * YEARS, obs)
    reloc = ", relocations on" if r.relocations_enabled else ""
    fig.legend(handles=[
        Line2D([], [], label="CoastSat (LOESS, 10 domains)", **OBSERVED),
        Line2D([], [], label=f"Model: rate bias {br:+.2f} m/yr, RMSE {er:.2f}; "
                             f"position bias {bp:+.1f} m, RMSE {ep:.1f}", **MODEL_ONE)],
        title=f"{SCEN_LABEL[r.scenario]}{reloc}, {r.source_sink_preset}, "
              f"{window(r.period).replace('_', '–')}",
        loc="outside lower center", ncol=1, frameon=False, fontsize=7)
    caption = (f"{SCEN_LABEL[r.scenario]}{reloc}, {PRESET_TEXT[r.source_sink_preset]}, "
               f"{window(r.period).replace('_', '-')}. Run {r.run_name}. " + COMMON_CAPTION)
    for png in pngs:
        save(fig, png, dpi=300, close=False)
        record_caption(png, caption)
    plt.close(fig)
    return dict(run_name=r.run_name, window=window(r.period), preset=r.source_sink_preset,
                scenario=r.scenario, relocations=bool(r.relocations_enabled),
                rate_bias_m_yr=br, rate_rmse_m_yr=er, position_bias_m=bp, position_rmse_m=ep)


def fig_scenarios(runs, period, preset, target, obs):
    sub = runs[(runs.period == period) & (runs.source_sink_preset == preset)
               & ~runs.relocations_enabled.astype(bool)]
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.8), sharex=True,
                             constrained_layout=True)
    ax_r, ax_p = axes
    ax_r.plot(target.index, target.values, zorder=6, **OBSERVED)
    ax_p.plot(obs.index, obs.values, zorder=6, **OBSERVED)
    handles = [Line2D([], [], label="CoastSat (LOESS, 10 domains)", **OBSERVED)]
    for scen in SCEN_ORDER:
        hit = sub[sub.scenario == scen]
        if hit.empty:
            continue
        rt = rates(hit.iloc[0].run_dir)
        ax_r.plot(rt.index, rt.lrr_m_yr, zorder=4, **SCEN_STYLE[scen])
        ax_p.plot(rt.index, rt.change_rate_m_yr * YEARS, zorder=4, **SCEN_STYLE[scen])
        br, er = score(rt.lrr_m_yr, target)
        bp, ep = score(rt.change_rate_m_yr * YEARS, obs)
        handles.append(Line2D([], [], **SCEN_STYLE[scen],
                              label=f"{SCEN_LABEL[scen]}: rate RMSE {er:.2f} m/yr, "
                                    f"position bias {bp:+.1f} m"))
    frame(axes, period)
    fig.legend(handles=handles, title=f"{PRESET_TEXT[preset]}, {window(period).replace('_', '–')}",
               loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    png = OUT / window(period) / f"scenarios_rate_and_position_{preset}_{window(period)}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"Every management scenario, {PRESET_TEXT[preset]}, {window(period).replace('_', '-')}, "
        "relocation arms left out (within 0.001 m/yr of their twins). " + COMMON_CAPTION))
    return png


def position_bound(runs):
    vals = [observed_change(p).abs().max() for p in runs.period.unique()]
    vals += [(rates(r.run_dir).change_rate_m_yr * YEARS).abs().max() for r in runs.itertuples()]
    half = float(np.ceil(max(vals) / 10.0) * 10.0)
    return (-half, half)


def main():
    global POSITION_YLIM
    apply_style()
    runs = matrix_runs()
    POSITION_YLIM = position_bound(runs)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "y_bounds.txt").write_text(
        f"rate panel: {RATE_YLIM[0]:g} to +{RATE_YLIM[1]:g} m/yr, as the runner's real-domains rate "
        "figure (rerender_run_figures.py --ylim-real=-7.5,7.5); its with-buffers figure is +/-10\n"
        f"position panel: {POSITION_YLIM[0]:g} to +{POSITION_YLIM[1]:g} m = max |position change| "
        "over every matrix run (GIS 1-90) and both observed windows, rounded up to 10 m\n",
        encoding="utf-8")
    rows, written = [], []
    for period in sorted(runs.period.unique()):
        target = common.coastsat_target(period)
        obs = observed_change(period)
        for r in runs[runs.period == period].itertuples():
            reloc = "_reloc" if r.relocations_enabled else ""
            by = (OUT / window(period) / "by_scenario"
                  / f"rate_and_position_{r.source_sink_preset}_{r.scenario}{reloc}_{window(period)}.png")
            rows.append(fig_one(r, target, obs, [r.run_dir / "figures" / "rate_and_position_change.png", by]))
            written.append(by)
        for preset in PRESETS:
            written.append(fig_scenarios(runs, period, preset, target, obs))
    s = pd.DataFrame(rows).sort_values(["window", "preset", "scenario", "relocations"])
    s.to_csv(OUT / "scores.csv", index=False)
    print(s.to_string(index=False, float_format=lambda v: f"{v:+.2f}"))
    print(f"\n{len(rows)} run figures written into each run's figures/, {len(written)} under {OUT}")


if __name__ == "__main__":
    main()
