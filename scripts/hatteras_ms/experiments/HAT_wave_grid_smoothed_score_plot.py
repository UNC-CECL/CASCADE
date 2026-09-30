"""Figures for wave-climate/2026-09-25-wave-grid-smoothed-score: the best settings on the
smoothed score, for each window and shared, one figure per scenario.

Two sources, drawn by the same code:
    grid    this study's runs (tables/all_runs.csv), once the sweep has run
    step2   the 2026-09-24 step-2 runs rescored on the smoothed output
            (tables/step2_rescored_smoothed.csv): drawn first, 2026-09-25,
            while the grid was still running (Hannah asked for the figures)

Writes, under output/raw_runs/experiments/wave-climate/2026-09-25-wave-grid-smoothed-score/figures/:
    best/<source>/per_period/best_by_period_<scenario>_<source>.png
    best/<source>/shared/best_shared_<scenario>_<source>.png
each with supporting/ (PDF, CAPTIONS.md, the data CSV).

Reads tables only; it never rebuilds the run index, so it is safe to run
while the sweep is going.

    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score_plot.py [--source step2|grid|both]

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_wave_grid_smoothed_score as study  # noqa: E402
import HAT_metres_2_wave_sensitivity_plot as p2  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    INK, INK_MUTED, DOMAIN_AXIS_LABEL, _title, apply_style, open_frame, record_caption,
    save, structures, support_dir, town_bands)

common = study.common
FIG = study.STUDY_DIR / "figures"
KEYS = list(study.KEYS)
COLOR = {"natural": "#1b7f6b", "full_management": "#c2571a"}
DARK = {"natural": "#0b4f43", "full_management": "#8a3410"}
# Which score ranks the runs (2026-09-27, Hannah: pick on the UNSMOOTHED model,
# as the runner's own RMSE and every matrix run are scored). "raw" draws the
# per-domain model as a dark line with no smoothed curve; figures go to
# figures/best/<source>_raw/.
SCORE = "smoothed"
OTHER = {"smoothed": "raw", "raw": "smoothed"}
NAME = {"natural": "Natural", "full_management": "Full management"}
SOURCE_NAME = {"grid": "four-parameter grid (this study)",
               "step2": "2026-09-24 step-2 runs, rescored on the smoothed output"}


def load(source):
    if source == "grid":
        t = pd.read_csv(study.TABLES_DIR / "all_runs.csv")
        t = t[t.status == "scored"].copy()
        t["run_path"] = [study.STUDY_DIR / d for d in t.run_dir]
    else:
        t = pd.read_csv(study.TABLES_DIR / "step2_rescored_smoothed.csv")
        t["run_path"] = [study.step2.STUDY_DIR / d for d in t.step2_run_dir]
    return t


def picks(t):
    """Best per window x scenario (smoothed share explained), and the shared
    setting per scenario (lowest mean smoothed RMSE / flat line, run in both)."""
    per, shared, table = {}, {}, []
    if SCORE == "raw":
        _, flat = study.targets()
        t = t.copy()
        t["raw_rmse_over_flat"] = t.raw_rmse_m_yr / t.period_start.map(flat)
    rel = f"{SCORE}_rmse_over_flat"
    for sc, p in product(study.SCENARIOS, study.PERIODS):
        x = t[(t.scenario == sc) & (t.period_start == p)]
        if not x.empty:
            per[(sc, p)] = x.loc[x[f"{SCORE}_variance_explained"].idxmax()]
    for sc in study.SCENARIOS:
        x = t[t.scenario == sc].drop_duplicates(["period_start", *KEYS])
        a = x[x.period_start == study.PERIODS[0]].set_index(KEYS)
        b = x[x.period_start == study.PERIODS[1]].set_index(KEYS)
        both = a[[rel]].join(b[[rel]],
                                                   lsuffix="_a", rsuffix="_b", how="inner")
        if both.empty:
            continue
        both["mean"] = both.mean(axis=1)
        key = both["mean"].idxmin()
        for p, frame in ((study.PERIODS[0], a), (study.PERIODS[1], b)):
            r = frame.loc[key].copy()
            for k, v in zip(KEYS, key):
                r[k] = v
            shared[(sc, p)] = r
        table.append(dict(scenario=sc, **dict(zip(KEYS, key)),
                          score=SCORE, mean_rmse_over_flat=float(both["mean"].min()),
                          candidates=len(both)))
    return per, shared, pd.DataFrame(table)


def settings_text(r):
    return (f"Hs {r['hs']:g} m, Tp {r['wave_period_s']:g} s, asym {r['wave_asymmetry']:g}, "
            f"high-angle {r['wave_angle_high_fraction']:g}")


def fig(targets, obs, chosen, scenario, rule, source):
    with plt.rc_context(p2.SCREEN_RC):
        return _fig(targets, obs, chosen, scenario, rule, source)


def _fig(targets, obs, chosen, scenario, rule, source):
    col = COLOR[scenario]
    f, axes = plt.subplots(2, 2, figsize=(16, 10.5), sharex=True, constrained_layout=True)
    rows = []
    for j, period in enumerate(study.PERIODS):
        r = chosen.get((scenario, period))
        ax_r, ax_p = axes[0, j], axes[1, j]
        w = study.window(period).replace("_", "–")
        if r is None:
            _title(ax_r, j, f"{w}\nno run yet")
            continue
        rt = pd.read_csv(Path(r["run_path"]) / "tables" / "shoreline_change_rate.csv"
                         ).set_index("gis_domain")
        for ax, o, m in ((ax_r, targets[period], rt.lrr_m_yr),
                         (ax_p, obs[period], rt.change_rate_m_yr * 14)):
            ax.plot(o.index, o.values, color=INK, lw=2.6, zorder=5)
            if SCORE == "raw":
                ax.plot(m.index, m.values, color=DARK[scenario], lw=2.0, zorder=4)
            else:
                ax.plot(m.index, m.values, color=col, lw=1.5, zorder=4)
                sm = common.smooth_like_target(m)
                ax.plot(sm.index, sm.values, color=col, lw=2.4, ls=(0, (5, 3)), alpha=0.5, zorder=4)
            ax.axhline(0, color=INK_MUTED, lw=0.6)
            ax.set_xlim(1, 90)
            ax.grid(axis="y")
            open_frame(ax)
            town_bands(ax, label=(ax is ax_r), fontsize=11)
        structures(ax_p, label=True, label_pt=11)
        structures(ax_r, label=False)
        _title(ax_r, j, f"{w}\n{settings_text(r)}\n"
                        f"{SCORE} {100 * r[SCORE + '_variance_explained']:+.0f}% explained "
                        f"({OTHER[SCORE]} {100 * r[OTHER[SCORE] + '_variance_explained']:+.0f}%), "
                        f"bias {r['bias_m_yr']:+.2f} m/yr")
        _title(ax_p, 2 + j, "")
        ax_p.set_xlabel(DOMAIN_AXIS_LABEL)
        rows.append(dict(period=w, scenario=scenario, **{k: r[k] for k in KEYS},
                         smoothed_variance_explained=r["smoothed_variance_explained"],
                         raw_variance_explained=r["raw_variance_explained"],
                         smoothed_rmse_m_yr=r["smoothed_rmse_m_yr"], bias_m_yr=r["bias_m_yr"],
                         run_path=str(r["run_path"])))
    axes[0, 0].set_ylabel("Shoreline change rate,\nLRR (m/yr)")
    axes[1, 0].set_ylabel("Shoreline position change,\nend minus start (m)")
    handles = [Line2D([], [], color=INK, lw=2.6, label="CoastSat (LOWESS, 10 domains)")]
    if SCORE == "raw":
        handles.append(Line2D([], [], color=DARK[scenario], lw=2.0,
                              label=f"Model, {NAME[scenario].lower()} (per domain, scored)"))
    else:
        handles += [Line2D([], [], color=col, lw=1.5, label=f"Model, {NAME[scenario].lower()}"),
                    Line2D([], [], color=col, lw=2.4, ls=(0, (5, 3)), alpha=0.5,
                           label="Model smoothed like CoastSat (scored)")]
    what = ("best settings for each window" if rule == "per_period"
            else "one setting for both windows")
    f.legend(handles=handles, loc="outside lower center", ncol=3, frameon=False,
             title=f"{NAME[scenario]}: {what}, from the {SOURCE_NAME[source]}")
    stem = f"best_{'by_period' if rule == 'per_period' else 'shared'}_{scenario}_{source}"
    png = FIG / "best" / (source + ("_raw" if SCORE == "raw" else "")) / rule / f"{stem}.png"
    save(f, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    model = ("unsmoothed (per-domain) model" if SCORE == "raw"
             else "model smoothed like the CoastSat target (LOWESS over 10 domains, the "
                  "southern 10 raw)")
    rule_text = (f"The run with the highest share of alongshore variation explained by the "
                 f"{model}, in each window. " if rule == "per_period" else
                 f"Among settings run in both windows, the one with the lowest mean of "
                 f"RMSE ({SCORE} model) divided by each window's flat-line RMSE. ")
    lines = ("Dark line: the model per domain, as scored." if SCORE == "raw" else
             "Solid colour: the model per domain; faint dashed: the smoothed model that is scored.")
    record_caption(png, (
        f"{NAME[scenario]}, {SOURCE_NAME[source]}. {rule_text}Score: 1 - SSE/SST of the "
        f"{model} against the CoastSat LRR target, interior GIS 2-89; the other score in "
        "brackets. Top: LRR rate; bottom: position change, end minus start, against the "
        f"observed CoastSat change. {lines} Offset in metres (dune line), no groin, "
        "no relocations."))
    return png


TOP_N = 5
TOP_COLORS = ["#1b1b8f", "#2f7fc1", "#3aa39a", "#d08c1f", "#b83a5e"]


def fig_top(targets, t, source):
    """The top TOP_N runs by smoothed score in each window x scenario, drawn as
    the smoothed profiles that were scored (added 2026-09-25: Hannah asked why
    the best figures looked the same as the raw-scored ones)."""
    with plt.rc_context(p2.SCREEN_RC):
        f, axes = plt.subplots(2, 2, figsize=(16, 11), sharex=True, constrained_layout=True)
        rows = []
        for i, sc in enumerate(study.SCENARIOS):
            for j, period in enumerate(study.PERIODS):
                ax = axes[i, j]
                x = t[(t.scenario == sc) & (t.period_start == period)].drop_duplicates(
                    ["period_start", *KEYS]).copy()
                x["rank_raw"] = x.raw_variance_explained.rank(ascending=False).astype(int)
                top = x.nlargest(TOP_N, f"{SCORE}_variance_explained")
                o = targets[period]
                ax.plot(o.index, o.values, color=INK, lw=3.0, zorder=6, label="CoastSat (LOWESS, 10 domains)")
                for k, (_, r) in enumerate(top.iterrows()):
                    rt = pd.read_csv(Path(r.run_path) / "tables" / "shoreline_change_rate.csv"
                                     ).set_index("gis_domain")
                    sm = rt.lrr_m_yr if SCORE == "raw" else common.smooth_like_target(rt.lrr_m_yr)
                    ax.plot(sm.index, sm.values, color=TOP_COLORS[k], lw=2.4 if k == 0 else 1.6,
                            zorder=5 - k * 0.1,
                            label=(f"#{k + 1}  {settings_text(r)}   "
                                   f"{100 * r[SCORE + '_variance_explained']:+.0f}% "
                                   + (f"(smoothed {100 * r.smoothed_variance_explained:+.0f}%)"
                                      if SCORE == "raw" else
                                      f"(raw {100 * r.raw_variance_explained:+.0f}%, raw rank #{r.rank_raw})")))
                    rows.append(dict(scenario=sc, period=study.window(period), rank=k + 1,
                                     raw_rank=int(r.rank_raw), **{kk: r[kk] for kk in KEYS},
                                     smoothed_variance_explained=r.smoothed_variance_explained,
                                     raw_variance_explained=r.raw_variance_explained,
                                     run_path=str(r.run_path)))
                ax.axhline(0, color=INK_MUTED, lw=0.6)
                ax.set_xlim(1, 90)
                ax.grid(axis="y")
                open_frame(ax)
                town_bands(ax, label=True, fontsize=11)
                structures(ax, label=False)
                _title(ax, 2 * i + j, f"{NAME[sc]}, {study.window(period).replace('_', '–')}")
                ax.legend(loc="lower left", fontsize=10, frameon=True, framealpha=0.9)
                if j == 0:
                    ax.set_ylabel(("Shoreline change rate" if SCORE == "raw"
                                   else "Smoothed shoreline change rate") + ",\nLRR (m/yr)")
                if i == 1:
                    ax.set_xlabel(DOMAIN_AXIS_LABEL)
        png = (FIG / "best" / (source + ("_raw" if SCORE == "raw" else ""))
               / f"top{TOP_N}_{SCORE}_profiles_{source}.png")
        save(f, png, dpi=300, close=True)
    pd.DataFrame(rows).to_csv(support_dir(png.parent) / f"{png.stem}.csv", index=False)
    record_caption(png, (
        f"The top {TOP_N} runs by smoothed score in each window and scenario, from the "
        f"{SOURCE_NAME[source]}, each drawn as the smoothed profile that was scored (LOWESS "
        "over 10 domains, the southern 10 raw) against the CoastSat LRR target (black). "
        "Legend: settings, smoothed share of alongshore variation explained, the raw "
        "(unsmoothed) score and the run's rank under the raw score. Interior GIS 2-89."))
    return png


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=("step2", "grid", "both"), default="both")
    ap.add_argument("--score", choices=("smoothed", "raw"), default="smoothed")
    a = ap.parse_args()
    global SCORE
    SCORE = a.score
    apply_style()
    targets = {p: common.coastsat_target(p) for p in study.PERIODS}
    obs = {p: p2.observed_change(p) for p in study.PERIODS}
    sources = ("step2", "grid") if a.source == "both" else (a.source,)
    for source in sources:
        if source == "grid" and not (study.TABLES_DIR / "all_runs.csv").is_file():
            continue
        t = load(source)
        per, shared, table = picks(t)
        table.to_csv(study.TABLES_DIR / f"best_shared_{source}{'_raw' if SCORE == 'raw' else ''}.csv",
                     index=False)
        print(fig_top(targets, t, source).relative_to(study.STUDY_DIR))
        for rule, chosen in (("per_period", per), ("shared", shared)):
            for sc in study.SCENARIOS:
                print(fig(targets, obs, chosen, sc, rule, source).relative_to(study.STUDY_DIR))


if __name__ == "__main__":
    main()
