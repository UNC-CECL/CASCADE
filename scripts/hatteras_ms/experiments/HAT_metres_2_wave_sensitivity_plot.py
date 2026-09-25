#!/usr/bin/env python3
r"""
HAT_metres_2_wave_sensitivity_plot.py -- the figures for the natural-scenario wave sensitivity
==============================================================================
Reads tables/all_runs.csv and tables/stage2_selection.csv (written by
HAT_metres_2_wave_sensitivity.py), each scored run's tables/shoreline_change_rate.csv,
the CoastSat LRR target and the observed CoastSat position change. Writes,
under output/raw_runs/experiments/2026-09-24-metres-2-wave-sensitivity/figures/:

  stage1/scores_by_<parameter>_both_periods.png
      share of alongshore variation explained, bias and RMSE against the
      parameter, both periods on one axis (cross-period consistency)
  alongshore/<period>/rate_and_position_change_by_<parameter>_<period>.png
      the modelled rate and position change along the island for each value
      against CoastSat
  management/natural_vs_full_management_baseline_both_periods.png
      the baseline under the natural scenario and under full management
  stage2/grid_<p1>_x_<p2>_both_periods.png
      the grid as lines: the score against the first parameter, one line per
      value of the second, a panel per period
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
import HAT_metres_2_wave_sensitivity as study  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    C, C_1984, C_1997, INK, INK_MUTED, DOMAIN_AXIS_LABEL, _title, apply_style,
    figsize, open_frame, record_caption, save, structures, support_dir,
    town_bands)
from site_layer.hat_topo_version import INIT_ROOT  # noqa: E402

FIG = study.STUDY_DIR / "figures"
common = study.common
# the two windows: the earlier red, the later blue (the house vintage pair)
PERIOD_STYLE = {1996: dict(color=C_1984, label="1996–2010"),
                2010: dict(color=C_1997, label="2010–2024")}
OBSERVED = dict(color=INK, lw=2.0)
NULL_STYLE = dict(color=INK, lw=0.9, ls="--")
PARAM_CMAP = {"wave_height": "Greys", "high_angle": "Purples",
              "asymmetry": "Oranges", "wave_period": "Greens"}
COMMON_CAPTION = ("Natural scenario (no road, beach or dune management, no fills, no "
                  "relocations, no groin), no imposed background erosion, island offset in "
                  "metres (dune line). Baseline Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle "
                  "fraction 0.45. Scores are the modelled LRR shoreline-change rate against "
                  "each window's CoastSat LRR target (LOESS, 10 domains) over the interior "
                  "domains GIS 2-89. Every value is fitted on the window it is scored on, "
                  "so a best value is a band, not a calibrated value.")


def no_score_marker(status):
    """A drowned barrier is a model result; a crash is not (Barrier3D's
    route_overwash access violation), so they get different marks."""
    return ("x", "drowned") if "drowned" in str(status) else ("D", "crashed")


def load():
    t = pd.read_csv(study.TABLES_DIR / "all_runs.csv")
    t["scored"] = t.status == "scored"
    return t


def baseline_mask(t, except_setting=None):
    keep = [k for k in study.BASELINE if k != except_setting]
    return np.logical_and.reduce([np.isclose(t[k], study.BASELINE[k]) for k in keep])


SCENARIO_LS = {study.SCENARIO: "-", study.MANAGED: "--"}
SCENARIO_LABEL = {study.SCENARIO: "Natural", study.MANAGED: "Full management"}


def one_at_a_time(t, group, period, scenario=study.SCENARIO):
    """The baseline and every stage-1 value of one parameter, in one period."""
    setting = study.PARAMS[group][0]
    p = t[(t.period_start == period) & (t.scenario == scenario)]
    return p[baseline_mask(p, setting)].drop_duplicates(setting).sort_values(setting)


def observed_change(period):
    f = (INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change"
         / study.window(period) / "smoothed" / "tables" / "domain_smoothed.csv")
    d = pd.read_csv(f)
    return d[d.window_domains == 10].set_index("domain_number")["observed_m"]


def run_table(row):
    return pd.read_csv(study.STUDY_DIR / row.run_dir / "tables" / "shoreline_change_rate.csv"
                       ).set_index("gis_domain")


def flat_null(period):
    return float(common.interior(common.coastsat_target(period)).std(ddof=0))


# =============================================================================
# STAGE 1: scores against each parameter, both periods
# =============================================================================

def fig_scores(t, group):
    setting, values, label = study.PARAMS[group]
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=2.9),
                             constrained_layout=True)
    panels = (("variance_explained", "Share of alongshore\nvariation explained"),
              ("mean_bias_interior_m_yr", "Interior mean bias (m/yr)"),
              ("rmse_interior_m_yr", "Interior RMSE (m/yr)"))
    rows = []
    scenarios = [sc for sc in (study.SCENARIO, study.MANAGED)
                 if len(one_at_a_time(t, group, study.PERIODS[0], sc)) > 1]
    for period, st in PERIOD_STYLE.items():
        for sc in scenarios:
            line = one_at_a_time(t, group, period, sc)
            ok = line[line.scored]
            for ax, (col, _) in zip(axes, panels):
                ax.plot(ok[setting], ok[col], marker="o" if sc == study.SCENARIO else "s",
                        ms=4, lw=1.4, ls=SCENARIO_LS[sc], color=st["color"],
                        mfc=st["color"] if sc == study.SCENARIO else "white")
                for _, d in line[~line.scored].iterrows():
                    m, _ = no_score_marker(d.status)
                    ax.plot(d[setting], 0.03, marker=m, ms=6 if m == "x" else 4.5, mew=1.4,
                            mfc="white" if m == "D" else None, color=st["color"],
                            transform=ax.get_xaxis_transform(), clip_on=False)
            rows += line.assign(parameter=group)[["parameter", "period", "scenario", setting,
                                                  "status", "variance_explained",
                                                  "mean_bias_interior_m_yr",
                                                  "rmse_interior_m_yr"]].to_dict("records")
        axes[2].axhline(flat_null(period), color=st["color"], lw=0.8, ls=":")
    axes[0].axhline(0, **NULL_STYLE)
    axes[1].axhline(0, color=INK_MUTED, lw=0.6)
    axes[0].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    for i, (ax, (_, ylab)) in enumerate(zip(axes, panels)):
        ax.axvline(study.BASELINE[setting], color=INK_MUTED, lw=0.8, ls=":")
        ax.set_xlabel(label)
        ax.set_ylabel(ylab)
        ax.grid(axis="y")
        open_frame(ax)
        _title(ax, i, "")
    handles = [Line2D([], [], marker="o", ms=4, lw=1.4, **st) for st in PERIOD_STYLE.values()]
    if len(scenarios) > 1:
        handles += [Line2D([], [], color=INK_MUTED, ls=SCENARIO_LS[sc],
                           marker="o" if sc == study.SCENARIO else "s", ms=4,
                           mfc=INK_MUTED if sc == study.SCENARIO else "white",
                           label=SCENARIO_LABEL[sc]) for sc in scenarios]
    handles += [Line2D([], [], label="Flat line at the observed mean", **NULL_STYLE),
                Line2D([], [], color=INK_MUTED, ls=":", label="Baseline value"),
                Line2D([], [], color=INK, ls="none", marker="x", label="Barrier drowned"),
                Line2D([], [], color=INK, ls="none", marker="D", ms=4.5, mfc="white",
                       label="Model crashed")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False, fontsize=7)
    png = FIG / "stage1" / f"scores_by_{group}_both_periods.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        f"Sensitivity to {label.lower()}, the other three wave parameters at the "
        "baseline, both windows (red 1996-2010, blue 2010-2024). (a) Share of the "
        "observed alongshore variation in shoreline-change rate explained: 1 - "
        "sum((model - observed)^2) / sum((observed - mean)^2); 0% (dashed) is a flat "
        "line at the observed mean. (b) Interior mean bias. (c) Interior RMSE, with "
        "each window's flat-line RMSE dotted in its colour. Solid lines and filled "
        "circles: natural scenario; dashed lines and open squares, where present: "
        "full management (road, beach and dune management, the historical fills; no "
        "relocations, no groin). The dotted vertical line is the baseline value; at the foot, crosses mark values where the barrier "
        "drowned and open diamonds values where the model crashed (a silent access "
        "violation in Barrier3D's overwash routing, not a model result); neither has a "
        "score. " + COMMON_CAPTION))
    return png


# =============================================================================
# ALONGSHORE: rate and position change, each value, each period
# =============================================================================

def fig_alongshore(t, group, period, target, obs, scenario=study.SCENARIO):
    setting, values, label = study.PARAMS[group]
    line = one_at_a_time(t, group, period, scenario)
    shade = dict(zip(line[setting], plt.get_cmap(PARAM_CMAP[group])(
        np.linspace(0.35, 0.95, len(line)))))
    fig, (ax_r, ax_p) = plt.subplots(2, 1, figsize=figsize("double", height=5.6),
                                     sharex=True, constrained_layout=True)
    ax_r.plot(target.index, target.values, zorder=6, **OBSERVED)
    ax_p.plot(obs.index, obs.values, zorder=6, **OBSERVED)
    handles = [Line2D([], [], label="CoastSat (LOESS, 10 domains)", **OBSERVED)]
    for _, row in line.iterrows():
        v = row[setting]
        if row.scored:
            rt = run_table(row)
            ax_r.plot(rt.index, rt.lrr_m_yr, color=shade[v], lw=1.2, zorder=4)
            ax_p.plot(rt.index, rt.change_rate_m_yr * 14, color=shade[v], lw=1.2, zorder=4)
            lab = (f"{v:g}{'  (baseline)' if np.isclose(v, study.BASELINE[setting]) else ''}"
                   f"  {100 * row.variance_explained:.0f}%, bias {row.mean_bias_interior_m_yr:+.2f} m/yr")
            handles.append(Line2D([], [], color=shade[v], lw=1.6, label=lab))
        else:
            m, why = no_score_marker(row.status)
            handles.append(Line2D([], [], color=shade[v], lw=0, marker=m, mfc="white" if m == "D" else None,
                                  label=f"{v:g}  ({'barrier drowned' if why == 'drowned' else 'model crashed'}, no output)"))
    for ax in (ax_r, ax_p):
        ax.axhline(0, color=INK_MUTED, lw=0.6)
        ax.set_xlim(1, 90)
        ax.grid(axis="y")
        open_frame(ax)
        town_bands(ax, label=(ax is ax_r))
    ax_r.set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    ax_p.set_ylabel(f"Shoreline position change,\n{period + 14} minus {period} (m)")
    ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    _title(ax_r, 0, "Rate")
    _title(ax_p, 1, "Position change")
    fig.legend(handles=handles, title=f"{label}, {study.window(period).replace('_', '-')}"
               + ("" if scenario == study.SCENARIO else ", full management"),
               loc="outside lower center", ncol=2, frameon=False, fontsize=7)
    structures(ax_p, label=True)
    structures(ax_r, label=False)
    tag = "" if scenario == study.SCENARIO else "_full_management"
    png = (FIG / "alongshore" / study.window(period)
           / f"rate_and_position_change_by_{group}{tag}_{study.window(period)}.png")
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"Alongshore response to {label.lower()}, {study.window(period).replace('_', '-')}, "
        "the other wave parameters at the baseline. Light to dark: increasing value; the "
        "legend gives each run's share of the observed alongshore rate variation "
        "explained and its interior bias. (a) Modelled LRR rate against the CoastSat LRR "
        "target (black). (b) Modelled shoreline position change, end minus start of the "
        "window, against the observed CoastSat change: the mean position over the last "
        "calendar year minus that over the first (5-scr/3-rates/coastsat/total_change, "
        "smoothed at 10 domains). Seaward positive. "
        + ("" if scenario == study.SCENARIO else
           "FULL MANAGEMENT here (road, beach and dune management, the historical "
           "fills; no relocations, no groin), not the natural scenario. ")
        + COMMON_CAPTION))
    return png


# =============================================================================
# MANAGEMENT: the baseline natural against full management
# =============================================================================

def fig_management(t, targets, obs):
    fig, axes = plt.subplots(2, 2, figsize=figsize("double", height=5.4), sharex=True,
                             constrained_layout=True)
    styles = {study.SCENARIO: dict(color=C["ACCENT"], lw=1.4, label="Natural"),
              study.MANAGED: dict(color=C["ADDED"], lw=1.4, label="Full management")}
    rows = []
    for j, period in enumerate(study.PERIODS):
        ax_r, ax_p = axes[0, j], axes[1, j]
        ax_r.plot(targets[period].index, targets[period].values, **OBSERVED)
        ax_p.plot(obs[period].index, obs[period].values, **OBSERVED)
        p = t[(t.period_start == period) & baseline_mask(t)
              & t.group.isin(["baseline", "baseline_full_management"])]
        for _, row in p.iterrows():
            if not row.scored:
                continue
            rt = run_table(row)
            ax_r.plot(rt.index, rt.lrr_m_yr, **styles[row.scenario])
            ax_p.plot(rt.index, rt.change_rate_m_yr * 14, **styles[row.scenario])
            rows.append(dict(period=row.period, scenario=row.scenario,
                             variance_explained=row.variance_explained,
                             bias=row.mean_bias_interior_m_yr, rmse=row.rmse_interior_m_yr))
        for ax in (ax_r, ax_p):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(1, 90)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r))
        _title(ax_r, j, f"Rate, {study.window(period).replace('_', '-')}")
        _title(ax_p, 2 + j, f"Position change, {study.window(period).replace('_', '-')}")
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
    axes[0, 0].set_ylabel("LRR (m/yr)")
    axes[1, 0].set_ylabel("End minus start (m)")
    handles = [Line2D([], [], label="CoastSat", **OBSERVED),
               *[Line2D([], [], **s) for s in styles.values()]]
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False)
    png = FIG / "management" / "natural_vs_full_management_baseline_both_periods.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        "The baseline wave climate (Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle 0.45) "
        "under the natural scenario (purple) and under full management (orange: road "
        "management, beach and dune management, the historical fills; relocations off; "
        "no groin), against CoastSat (black), in each window. Top: LRR rate; bottom: "
        "position change, end minus start. No imposed background erosion, offset in "
        "metres. Scores are in the CSV under supporting/."))
    return png


# =============================================================================
# STAGE 2: the grid as lines
# =============================================================================

# The two grids: the one the stage-2 rule chose (both windows), and the Hs x
# Tp grid the first rule chose, stopped, then finished for 1996-2010 only.
GRIDS = (
    ("stage2_selection.csv", study.PERIODS,
     "the two parameters whose range over the runs that survived in 1996-2010 moved "
     "the share of alongshore variation explained the most (tables/stage2_selection.csv)"),
    ("stage2_selection_first_rule.csv", (1996,),
     "the pair the first stage-2 rule chose (tables/stage2_selection_first_rule.csv); "
     "stopped when that rule was replaced, then finished for 1996-2010 only because "
     "its best cell was the best natural 1996-2010 run of the study"),
)


def fig_grid(t, selection="stage2_selection.csv", periods=study.PERIODS, why=""):
    sel_path = study.TABLES_DIR / selection
    if not sel_path.is_file():
        return None
    sel = pd.read_csv(sel_path)
    chosen = list(sel[sel.chosen_for_grid].parameter)
    g1, g2 = chosen
    s1, s2 = study.PARAMS[g1][0], study.PARAMS[g2][0]
    v1 = [float(x) for x in sel.set_index("parameter").loc[g1, "grid_values"].split(",")]
    v2 = [float(x) for x in sel.set_index("parameter").loc[g2, "grid_values"].split(",")]
    others = [k for k in study.BASELINE if k not in (s1, s2)]
    fig, axes = plt.subplots(2, len(periods), sharex=True, squeeze=False,
                             figsize=figsize("double" if len(periods) > 1 else "single",
                                             height=5.2),
                             constrained_layout=True)
    ramp = dict(zip(v2, plt.get_cmap(PARAM_CMAP[g2])(np.linspace(0.35, 0.95, len(v2)))))
    rows = []
    for j, period in enumerate(periods):
        p = t[(t.period_start == period) & (t.scenario == study.SCENARIO)]
        p = p[np.logical_and.reduce([np.isclose(p[k], study.BASELINE[k]) for k in others])]
        p = p.drop_duplicates([s1, s2])
        for b in v2:
            line = p[np.isclose(p[s2], b) & p[s1].isin(v1)].sort_values(s1)
            ok = line[line.scored]
            for i, col in enumerate(("variance_explained", "mean_bias_interior_m_yr")):
                axes[i, j].plot(ok[s1], ok[col], marker="o", ms=3.5, lw=1.3, color=ramp[b])
                for _, d in line[~line.scored].iterrows():
                    m, _ = no_score_marker(d.status)
                    axes[i, j].plot(d[s1], 0.03, marker=m, ms=5 if m == "x" else 4, color=ramp[b],
                                    mfc="white" if m == "D" else None,
                                    transform=axes[i, j].get_xaxis_transform(), clip_on=False)
            rows += line[[s1, s2, "period", "status", "variance_explained",
                          "mean_bias_interior_m_yr", "rmse_interior_m_yr"]].to_dict("records")
        _title(axes[0, j], j, study.window(period).replace("_", "-"))
        axes[1, j].set_xlabel(study.PARAMS[g1][2])
    axes[0, 0].set_ylabel("Share of alongshore\nvariation explained")
    axes[1, 0].set_ylabel("Interior mean bias (m/yr)")
    for ax in axes.flat:
        ax.grid(axis="y")
        open_frame(ax)
    for ax in axes[0]:
        ax.axhline(0, **NULL_STYLE)
        ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    for ax in axes[1]:
        ax.axhline(0, color=INK_MUTED, lw=0.6)
    handles = [Line2D([], [], color=ramp[b], marker="o", ms=3.5, label=f"{b:g}") for b in v2]
    handles.append(Line2D([], [], label="Flat line at the observed mean", **NULL_STYLE))
    fig.legend(handles=handles, title=study.PARAMS[g2][2], loc="outside lower center",
               ncol=len(handles) if len(periods) > 1 else 3, frameon=False, fontsize=7)
    which = ("both_periods" if len(periods) > 1
             else study.window(periods[0]))
    png = FIG / "stage2" / f"grid_{g1}_x_{g2}_{which}.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        f"Stage 2: {study.PARAMS[g1][2].lower()} against "
        f"{study.PARAMS[g2][2].lower()}, {why or 'the pair the stage-2 rule chose'}; "
        "the other two at the baseline. Top: share of the observed alongshore variation "
        "explained (0% dashed = a flat line at the observed mean); bottom: interior mean "
        "bias. Columns: the window(s). Lines: values of the second parameter, light to "
        "dark. At the foot: crosses drowned, open diamonds crashed. " + COMMON_CAPTION))
    return png


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    apply_style()
    t = load()
    targets = {p: common.coastsat_target(p) for p in study.PERIODS}
    obs = {p: observed_change(p) for p in study.PERIODS}
    out = [fig_scores(t, g) for g in study.PARAMS]
    out += [fig_alongshore(t, g, p, targets[p], obs[p])
            for p in study.PERIODS for g in study.PARAMS]
    if (t.group.str.startswith(study.MANAGED_PREFIX)).any():
        out += [fig_alongshore(t, g, p, targets[p], obs[p], study.MANAGED)
                for p in study.PERIODS for g in study.PARAMS]
    out.append(fig_management(t, targets, obs))
    for selection, periods, why in GRIDS:
        grid = fig_grid(t, selection, periods, why)
        if grid:
            out.append(grid)
    for f in out:
        print(f.relative_to(study.STUDY_DIR))


if __name__ == "__main__":
    main()
