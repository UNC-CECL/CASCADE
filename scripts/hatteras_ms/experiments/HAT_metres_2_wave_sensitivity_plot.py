#!/usr/bin/env python3
"""
The figures for the natural-scenario wave sensitivity.

    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity_plot.py

Stage-1 scores, alongshore profiles, management, the stage-2 grids and the
best settings. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
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

# --- CONFIG ------------------------------------------------------------------
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
                  "each window's CoastSat LRR target (LOWESS, 7 domains) over the interior "
                  "domains GIS 2-89. Every value is fitted on the window it is scored on, "
                  "so a best value is a band, not a calibrated value.")
# -----------------------------------------------------------------------------


# A marker for a run with no score: drowned or crashed
def no_score_marker(status):
    return ("x", "drowned") if "drowned" in str(status) else ("D", "crashed")


# The study's run table, with a scored flag
def load():
    t = pd.read_csv(study.TABLES_DIR / "all_runs.csv")
    t["scored"] = t.status == "scored"
    return t


# Rows at the baseline in every setting but one
def baseline_mask(t, except_setting=None):
    keep = [k for k in study.BASELINE if k != except_setting]
    return np.logical_and.reduce([np.isclose(t[k], study.BASELINE[k]) for k in keep])


SCENARIO_LS = {study.SCENARIO: "-", study.MANAGED: "--"}
SCENARIO_LABEL = {study.SCENARIO: "Natural", study.MANAGED: "Full management"}


# The baseline and every stage-1 value of one parameter, in one period
def one_at_a_time(t, group, period, scenario=study.SCENARIO):
    setting = study.PARAMS[group][0]
    p = t[(t.period_start == period) & (t.scenario == scenario)]
    return p[baseline_mask(p, setting)].drop_duplicates(setting).sort_values(setting)


# Observed CoastSat position change per domain for a window
def observed_change(period):
    f = (INIT_ROOT / "5-scr" / "3-rates" / "coastsat" / "total_change"
         / study.window(period) / "smoothed" / "tables" / "domain_smoothed.csv")
    d = pd.read_csv(f)
    return d[d.window_domains == 7].set_index("domain_number")["observed_m"]


# A table row's rate table
def run_table(row):
    return pd.read_csv(study.STUDY_DIR / row.run_dir / "tables" / "shoreline_change_rate.csv"
                       ).set_index("gis_domain")


# A window's flat-line RMSE: the target's interior spread
def flat_null(period):
    return float(common.interior(common.coastsat_target(period)).std(ddof=0))


# Stage 1: scores against each parameter, both periods

# Share explained, bias and RMSE against one parameter, both periods
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


# Alongshore: rate and position change, each value, each period

# Rate and position change along the island for each value of one parameter
def fig_alongshore(t, group, period, target, obs, scenario=study.SCENARIO):
    setting, values, label = study.PARAMS[group]
    line = one_at_a_time(t, group, period, scenario)
    shade = dict(zip(line[setting], plt.get_cmap(PARAM_CMAP[group])(
        np.linspace(0.35, 0.95, len(line)))))
    fig, (ax_r, ax_p) = plt.subplots(2, 1, figsize=figsize("double", height=5.6),
                                     sharex=True, constrained_layout=True)
    ax_r.plot(target.index, target.values, zorder=6, **OBSERVED)
    ax_p.plot(obs.index, obs.values, zorder=6, **OBSERVED)
    handles = [Line2D([], [], label="CoastSat (LOWESS, 7 domains)", **OBSERVED)]
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
           / ("natural" if scenario == study.SCENARIO else "full_management")
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
        "smoothed at 7 domains). Seaward positive. "
        + ("" if scenario == study.SCENARIO else
           "FULL MANAGEMENT here (road, beach and dune management, the historical "
           "fills; no relocations, no groin), not the natural scenario. ")
        + COMMON_CAPTION))
    return png


# Management: the baseline natural against full management

# The baseline, natural against full management
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


# Stage 2: the grid as lines

# The two grids: the stage-2 rule's, and the first rule's Hs x Tp (1996-2010 only)
GRIDS = (
    ("stage2_selection.csv", study.PERIODS,
     "the two parameters whose range over the runs that survived in 1996-2010 moved "
     "the share of alongshore variation explained the most (tables/stage2_selection.csv)"),
    ("stage2_selection_first_rule.csv", (1996,),
     "the pair the first stage-2 rule chose (tables/stage2_selection_first_rule.csv); "
     "stopped when that rule was replaced, then finished for 1996-2010 only because "
     "its best cell was the best natural 1996-2010 run of the study"),
)


# A stage-2 grid as lines, a panel per period
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
    others = [k for k in study.STAGE2_BASELINE if k not in (s1, s2)]
    fig, axes = plt.subplots(2, len(periods), sharex=True, squeeze=False,
                             figsize=figsize("double" if len(periods) > 1 else "single",
                                             height=5.2),
                             constrained_layout=True)
    ramp = dict(zip(v2, plt.get_cmap(PARAM_CMAP[g2])(np.linspace(0.35, 0.95, len(v2)))))
    rows = []
    for j, period in enumerate(periods):
        p = t[(t.period_start == period) & (t.scenario == study.SCENARIO)]
        p = p[np.logical_and.reduce([np.isclose(p[k], study.STAGE2_BASELINE[k]) for k in others])]
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


# Best settings: each period's own, and one set for both (added 2026-09-24)

KEYS = ["hs", "wave_period_s", "wave_asymmetry", "wave_angle_high_fraction"]
SCEN_LINE = {study.SCENARIO: dict(color=C["ACCENT"], lw=1.4, label="Natural"),
             study.MANAGED: dict(color=C["ADDED"], lw=1.4, label="Full management")}


# A setting, as text
def settings_text(r):
    return (f"Hs {r.hs:g}, Tp {r.wave_period_s:g}, asym {r.wave_asymmetry:g}, "
            f"high-angle {r.wave_angle_high_fraction:g}")


# The setting with the highest share explained, per scenario and period
def best_per_period(t):
    out = {}
    for sc in (study.SCENARIO, study.MANAGED):
        for p in study.PERIODS:
            x = t[(t.scenario == sc) & (t.period_start == p) & t.scored]
            out[(sc, p)] = x.loc[x.variance_explained.idxmax()]
    return out


# The one setting best across both periods, per scenario
def best_both_periods(t):
    flat = pd.read_csv(study.TABLES_DIR / "observed_targets.csv").set_index("period")[
        "flat_line_rmse_m_yr"]
    out, table = {}, []
    for sc in (study.SCENARIO, study.MANAGED):
        x = t[(t.scenario == sc) & t.scored].copy()
        x["rel"] = x.rmse_interior_m_yr / x.period.map(flat)
        per = {p: x[x.period_start == p].drop_duplicates(KEYS).set_index(KEYS)
               for p in study.PERIODS}
        both = per[study.PERIODS[0]][["rel"]].join(per[study.PERIODS[1]][["rel"]],
                                                   lsuffix="_a", rsuffix="_b", how="inner")
        both["mean_rel"] = (both.rel_a + both.rel_b) / 2
        key = both.mean_rel.idxmin()
        for p in study.PERIODS:
            r = per[p].loc[key].copy()
            for k, v in zip(KEYS, key):
                r[k] = v
            out[(sc, p)] = r
        table.append(dict(scenario=sc, **dict(zip(KEYS, key)),
                          mean_rmse_over_flat=float(both.mean_rel.min()),
                          candidates=len(both)))
    return out, pd.DataFrame(table)


# One best-settings figure per scenario and rule, drawn for the screen
BEST_COLOR = {study.SCENARIO: "#1b7f6b", study.MANAGED: "#c2571a"}
BEST_NAME = {study.SCENARIO: "Natural", study.MANAGED: "Full management"}
smooth_like_target = common.smooth_like_target   # one implementation, shared


# A best-settings figure, drawn for the screen
def fig_best(targets, obs, picks, scenario, rule):
    with plt.rc_context(SCREEN_RC):
        return _fig_best(targets, obs, picks, scenario, rule)


# The best-settings figure itself
def _fig_best(targets, obs, picks, scenario, rule):
    col = BEST_COLOR[scenario]
    fig, axes = plt.subplots(2, 2, figsize=(16, 10.5), sharex=True, constrained_layout=True)
    rows = []
    for j, period in enumerate(study.PERIODS):
        r = picks[(scenario, period)]
        rt = pd.read_csv(study.STUDY_DIR / r.run_dir / "tables" / "shoreline_change_rate.csv"
                         ).set_index("gis_domain")
        rate, change = rt.lrr_m_yr, rt.change_rate_m_yr * 14
        ax_r, ax_p = axes[0, j], axes[1, j]
        for ax, o, m in ((ax_r, targets[period], rate), (ax_p, obs[period], change)):
            ax.plot(o.index, o.values, color=INK, lw=2.6, zorder=5)
            ax.plot(m.index, m.values, color=col, lw=1.5, zorder=4)
            sm = smooth_like_target(m)
            ax.plot(sm.index, sm.values, color=col, lw=2.4, ls=(0, (5, 3)), alpha=0.5,
                    zorder=4)
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(1, 90)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r), fontsize=11)
        structures(ax_p, label=True, label_pt=11)
        structures(ax_r, label=False)
        w = study.window(period).replace("_", "–")
        _title(ax_r, j, f"{w}\n{settings_text(r)}\n"
                        f"{100 * r.variance_explained:+.0f}% explained, "
                        f"bias {r.mean_bias_interior_m_yr:+.2f} m/yr")
        _title(ax_p, 2 + j, "")
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
        rows.append(dict(period=r.period, scenario=scenario, **{k: r[k] for k in KEYS},
                         variance_explained=r.variance_explained,
                         bias_m_yr=r.mean_bias_interior_m_yr,
                         rmse_m_yr=r.rmse_interior_m_yr, run_dir=r.run_dir))
    axes[0, 0].set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    axes[1, 0].set_ylabel("Shoreline position change,\nend minus start (m)")
    handles = [Line2D([], [], color=INK, lw=2.6, label="CoastSat (LOWESS, 7 domains)"),
               Line2D([], [], color=col, lw=1.5, label=f"Model, {BEST_NAME[scenario].lower()}"),
               Line2D([], [], color=col, lw=2.4, ls=(0, (5, 3)), alpha=0.5,
                      label="Model smoothed like CoastSat (LOWESS, 7 domains)")]
    what = ("the best wave settings found for each window" if rule == "per_period"
            else "one wave setting for both windows")
    fig.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
               title=f"{BEST_NAME[scenario]}: {what}")
    tag = "full_management" if scenario == study.MANAGED else "natural"
    png = FIG / "best" / rule / f"best_{'settings_by_period' if rule == 'per_period' else 'shared_settings'}_{tag}.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    head = ("The best wave settings found for each window separately: the run with the "
            "highest share of alongshore variation explained among every scored run of the "
            "study (the one-at-a-time sweeps, the Hs x high-angle grid in both windows, the "
            "Hs x Tp grid in 1996-2010 only, and the hand-picked combos). No search covered "
            "all four parameters jointly. In 2010-2024 the best run is the least bad: no "
            "setting beats a flat line there. " if rule == "per_period" else
            "One wave setting for both windows: among settings run in both windows, the one "
            "with the lowest mean of RMSE divided by each window's flat-line RMSE "
            "(tables/best_settings_both_periods.csv). The Hs x Tp grid was run in 1996-2010 "
            "only and so cannot be picked. The choice is dominated by 2010-2024, where every "
            "setting is far worse than a flat line. ")
    record_caption(png, (
        f"{BEST_NAME[scenario]}. " + head + "Top: the modelled LRR rate along the island "
        "against the CoastSat LRR target (black); bottom: the modelled position change, end "
        "minus start of the window, against the observed CoastSat change (mean position over "
        "the last calendar year minus the first, smoothed at 7 domains). Solid colour: the "
        "model per domain; faint dashed: the model smoothed as the target is (LOWESS over 10 "
        "domains, the southern 10 left raw), for comparing like with like. Settings and "
        "scores (interior GIS 2-89) above each column. Offset in metres (dune line), zeroBE, "
        "no groin, no relocations."))
    return png


# The best-settings figures, per rule and scenario
def best_figures(t, targets, obs):
    per = best_per_period(t)
    joint, table = best_both_periods(t)
    table.to_csv(study.TABLES_DIR / "best_settings_both_periods.csv", index=False)
    return [fig_best(targets, obs, picks, sc, rule)
            for rule, picks in (("per_period", per), ("shared", joint))
            for sc in (study.SCENARIO, study.MANAGED)]


# The asymmetry x high-angle 2x2 (added 2026-09-25)

# The corners drawn: the old baseline and the combination
CORNERS = [  # (asymmetry, high-angle, style)
    (0.8, 0.45, dict(color=C["BASE"], lw=2.0, label="Old baseline: asym 0.8, high-angle 0.45")),
    (0.7, 0.4, dict(color=C["ACCENT"], lw=2.4, label="Combination: asym 0.7, high-angle 0.4")),
]
# Drawn for reading on screen, not for the page: larger canvas and type.
SCREEN_RC = {"font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13,
             "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
             "legend.title_fontsize": 13}
BUXTON_ZOOM = (1, 16)


# One asymmetry x high-angle corner's run, or None
def corner_row(t, period, scenario, a, f):
    x = t[t.scored & (t.period_start == period) & (t.scenario == scenario)
          & np.isclose(t.hs, 1.0) & np.isclose(t.wave_period_s, 8.0)
          & np.isclose(t.wave_asymmetry, a) & np.isclose(t.wave_angle_high_fraction, f)]
    return None if x.empty else x.iloc[0]


# The 2x2 combination figure, drawn for the screen
def fig_combo(t, period, target, obs):
    with plt.rc_context(SCREEN_RC):
        return _fig_combo(t, period, target, obs)


# The combination figure itself
def _fig_combo(t, period, target, obs):
    w = study.window(period).replace("_", "-")
    fig = plt.figure(figsize=(16, 10.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[3, 3, 1.35])
    axes = [[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(2)]
    cols = [(study.SCENARIO, "Natural"), (study.MANAGED, "Full management"),
            (study.MANAGED, "Close-up")]
    rows, handles = [], []
    for j, (sc, name) in enumerate(cols):
        ax_r, ax_p = axes[0][j], axes[1][j]
        zoom = j == 2
        ax_r.plot(target.index, target.values, **{**OBSERVED, "lw": 2.6})
        ax_p.plot(obs.index, obs.values, **{**OBSERVED, "lw": 2.6})
        for a, f, st in CORNERS:
            r = corner_row(t, period, sc, a, f)
            if r is None:
                continue
            rt = run_table(r)
            ax_r.plot(rt.index, rt.lrr_m_yr, **{k: v for k, v in st.items() if k != "label"})
            ax_p.plot(rt.index, rt.change_rate_m_yr * 14,
                      **{k: v for k, v in st.items() if k != "label"})
            if not zoom:
                rows.append(dict(period=w, scenario=sc, wave_asymmetry=a,
                                 wave_angle_high_fraction=f,
                                 variance_explained=r.variance_explained,
                                 bias_m_yr=r.mean_bias_interior_m_yr,
                                 rmse_m_yr=r.rmse_interior_m_yr, run_dir=r.run_dir))
        lim = BUXTON_ZOOM if zoom else (1, 90)
        for ax in (ax_r, ax_p):
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(*lim)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r), fontsize=11)
        if zoom:
            for ax, key in ((ax_r, None), (ax_p, None)):
                lines = [l.get_ydata() for l in ax.get_lines()]
                xs = [l.get_xdata() for l in ax.get_lines()]
                vals = np.concatenate([np.asarray(y)[(np.asarray(x) >= lim[0]) & (np.asarray(x) <= lim[1])]
                                       for x, y in zip(xs, lines) if len(np.atleast_1d(x)) > 2])
                pad = 0.08 * (vals.max() - vals.min())
                ax.set_ylim(vals.min() - pad, vals.max() + pad)
        structures(ax_p, label=not zoom, label_pt=11)
        structures(ax_r, label=False)
        _title(ax_r, j, name)
        _title(ax_p, 3 + j, "")
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL if not zoom else "GIS domain")
    axes[0][0].set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    axes[1][0].set_ylabel(f"Shoreline position change,\n{period + 14} minus {period} (m)")
    scores = {(r["scenario"], r["wave_asymmetry"], r["wave_angle_high_fraction"]): r for r in rows}
    handles = [Line2D([], [], **{**OBSERVED, "lw": 2.6}, label="CoastSat (LOWESS, 7 domains)")]
    for a, f, st in CORNERS:
        n = scores.get((study.SCENARIO, a, f)); m = scores.get((study.MANAGED, a, f))
        tail = (f"   natural {100 * n['variance_explained']:+.0f}%, managed "
                f"{100 * m['variance_explained']:+.0f}%") if n and m else ""
        handles.append(Line2D([], [], color=st["color"], lw=st["lw"], label=st["label"] + tail))
    fig.legend(handles=handles, loc="outside lower center", ncol=1, frameon=False,
               title=f"Hs 1.0 m, Tp 8 s, {w}: share of the alongshore variation explained")
    png = FIG / "combos" / f"old_baseline_vs_asym0.7_highangle0.4_{study.window(period)}.png"
    save(fig, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        f"The old baseline (asymmetry 0.8, high-angle 0.45; grey) against Hannah's "
        f"combination (asymmetry 0.7, high-angle 0.4; purple), Hs 1.0 m and Tp 8 s, {w}, "
        "against CoastSat (black). Top: the modelled LRR rate against the "
        "CoastSat LRR target; bottom: the modelled position change, end minus start, against "
        "the observed CoastSat change. Left, natural; middle, full management; right, the "
        "full-management runs close up on Buxton (GIS 1-16), where the observed rate dips to "
        "-2.1 m/yr at GIS 7; the combination dips to more than twice that at GIS 6. "
        "Legend: share of the alongshore rate variation explained, interior "
        "GIS 2-89. Offset in metres (dune line), zeroBE, no groin, no relocations."))
    return png


# Run: every figure
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
    out += best_figures(t, targets, obs)
    out += [fig_combo(t, p, targets[p], obs[p]) for p in study.PERIODS]
    for selection, periods, why in GRIDS:
        grid = fig_grid(t, selection, periods, why)
        if grid:
            out.append(grid)
    for f in out:
        print(f.relative_to(study.STUDY_DIR))


if __name__ == "__main__":
    main()
