"""
The 2010-2026 test run against five alternative CoastSat targets, beside the calibration one.

    python scripts/analyze_output/compare_runs/hindcast_vs_observed/alternative_targets.py

alt 1: the 1996-2026 LRR against the modelled 2010-2026 rate.
alt 2: model start + 16 x the 2010-2026 LRR, against the modelled end position.
alt 3: the 1996 start + 30 x the 1996-2026 LRR, against the modelled end position.
alt 4: CoastSat mean position over calendar 2025, against the modelled end position.
alt 5: model start + 16 x the 1996-2026 LRR, against the modelled end position.
Writes output/comparisons/alternative_targets/2010_2026/. Details: README.md there.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-05
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

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as _REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR / "analyze_output" / "compare_runs" / "matrix_vs_observed"))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "5-scr" / "lib"))
import matrix_vs_observed as mvo  # noqa: E402  (run loaders, model end change, scoring)
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
from coastsat_vs_duneline import load_chainage  # noqa: E402

from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK_MUTED, _title, apply_style,
    compare_header, figsize, mark_offaxis, offaxis_clause, open_frame,
    record_caption, save, structures, town_bands)
from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT  # noqa: E402
from site_layer.hat_topo_version import INIT_ROOT  # noqa: E402

common = mvo.common

# --- CONFIG ------------------------------------------------------------------
PERIOD = 2010                       # the test run
LONG_WINDOW = (1996, 2026)          # the long-term LRR (alts 1, 3, 5)
LONG_START = 1996                   # alt 3's start shoreline
OBS_YEAR = 2025                     # alt 4: calendar 2025, the last full year (Hannah, 2026-10-05)
# 16 x for alts 2 and 5, the run length (the methodology text says 15; Hannah chose 16, 2026-10-05)
SCALE_YEARS = mvo.years(PERIOD)
LONG_YEARS = LONG_WINDOW[1] - LONG_WINDOW[0]
# The CoastSat means within +/-1 yr of each start DEM (the shoreline offsets, v2)
START_WINDOWS = {2010: "2008-08-17_2010-08-17", 1996: "1995-10-12_1997-10-12"}
MEAN_ROOT = INIT_ROOT / "5-scr" / "1-observations" / "mean_shoreline"
OUT = COMPARISONS_ROOT / "alternative_targets" / f"{PERIOD}_{mvo.end_year(PERIOD)}"
FIG_PRESET = "zeroBE"               # no source/sink anywhere (methodology note 3b)
FIG_SCENARIO = "full_management"
POS_HALF = 100.0                    # m, as target_comparison
RATE_HALF = 7.5                     # m/yr, as matrix_vs_observed
N_DOMAINS = 90

C_OBS = "#92c5de"                   # observed pale and thick, model dark and thin
C_MODEL = "#2166ac"
C_CURRENT = INK_MUTED               # the calibration-style target, for reference
LW_OBS, LW_MODEL = 2.8, 1.2

_S10 = START_WINDOWS[PERIOD].replace("_", " to ")
_S96 = START_WINDOWS[LONG_START].replace("_", " to ")
# target key -> (quantity, short label, long description)
TARGETS = {
    "current_rate": ("rate", "Calibration method: CoastSat LRR 2010–2026",
                     "the calibration procedure's target, the CoastSat LRR fitted on 2010-2026"),
    "alt1_rate_1996_2026": ("rate", "Alt 1: CoastSat LRR 1996–2026",
                            "the CoastSat LRR fitted on 1996-2026"),
    "alt2_start_plus_16x_lrr_2010_2026": (
        "position", f"Alt 2: 2010 start + {SCALE_YEARS} × LRR 2010–2026",
        f"the 2010 start shoreline (CoastSat mean {_S10}) moved by {SCALE_YEARS} x the "
        "2010-2026 CoastSat LRR"),
    "alt3_1996_start_plus_30x_lrr_1996_2026": (
        "position", f"Alt 3: 1996 start + {LONG_YEARS} × LRR 1996–2026",
        f"the 1996 start shoreline (CoastSat mean {_S96}) moved by {LONG_YEARS} x the "
        "1996-2026 CoastSat LRR, drawn relative to the 2010 start"),
    "alt4_mean_2025": (
        "position", "Alt 4: CoastSat mean position, 2025",
        f"the mean CoastSat position over calendar {OBS_YEAR}, relative to the 2010 start"),
    "alt5_start_plus_16x_lrr_1996_2026": (
        "position", f"Alt 5: 2010 start + {SCALE_YEARS} × LRR 1996–2026",
        f"the 2010 start shoreline moved by {SCALE_YEARS} x the 1996-2026 CoastSat LRR"),
}
# Panel order and titles on the per-run figure
PANELS = (
    ("alt1_rate_1996_2026", "Alt 1: LRR 1996–2026 against the modelled 2010–2026 rate"),
    ("alt2_start_plus_16x_lrr_2010_2026", f"Alt 2: 2010 start + {SCALE_YEARS} × LRR 2010–2026"),
    ("alt3_1996_start_plus_30x_lrr_1996_2026", f"Alt 3: 1996 start + {LONG_YEARS} × LRR 1996–2026"),
    ("alt4_mean_2025", f"Alt 4: mean CoastSat position, {OBS_YEAR}"),
    ("alt5_start_plus_16x_lrr_1996_2026", f"Alt 5: 2010 start + {SCALE_YEARS} × LRR 1996–2026"),
)
# -----------------------------------------------------------------------------


# A domain LRR summary as a GIS-indexed series
def domain_lrr(window):
    d = pd.read_csv(COASTSAT_LRR_ROOT / "{}_{}".format(*window) / "domain_lrr_summary.csv")
    return d.set_index(d["domain_number"].astype(int))["mean_lrr"].reindex(range(1, N_DOMAINS + 1))


# Per-transect CoastSat mean around a start DEM; transects left out of the offset are NaN
def start_means(year):
    w = START_WINDOWS[year]
    t = pd.read_csv(MEAN_ROOT / w / f"transect_means_{w}.csv")
    t = t[t["domain_number"].between(1, N_DOMAINS)].set_index("transect_id")
    return t["mean_chainage_m"].where(t["included"].astype(bool)), t["domain_number"].astype(int)


# (mean chainage, n) over one calendar year
def year_mean(df, year):
    sel = df.loc[df["date"].dt.year == year, "chainage"] if df is not None else pd.Series(dtype=float)
    return (float(sel.mean()), int(sel.size)) if sel.size else (np.nan, 0)


# Per transect: both start positions, the 2025 mean, and the alt 3 and alt 4 changes
def transects():
    p10, dom = start_means(PERIOD)
    p96, _ = start_means(LONG_START)
    lrr = pd.read_csv(COASTSAT_LRR_ROOT / "{}_{}".format(*LONG_WINDOW) / "transect_lrr_full.csv"
                      ).set_index("transect_id")["lrr_m_yr"]
    rows = []
    for tid in p10.index:
        p, n = year_mean(load_chainage(tid), OBS_YEAR)
        rows.append(dict(transect_id=tid, position_2025_m=p, n_2025=n))
    t = pd.DataFrame(rows).set_index("transect_id")
    t.insert(0, "domain_number", dom)
    t["position_start_2010_m"] = p10
    t["position_start_1996_m"] = p96
    t["lrr_1996_2026_m_yr"] = lrr.reindex(t.index)
    t["alt3_end_position_m"] = t["position_start_1996_m"] + LONG_YEARS * t["lrr_1996_2026_m_yr"]
    t["alt3_change_vs_2010_start_m"] = t["alt3_end_position_m"] - t["position_start_2010_m"]
    t["alt4_change_vs_2010_start_m"] = t["position_2025_m"] - t["position_start_2010_m"]
    return t.reset_index()


# Every target per domain, raw (domain means) and smoothed as the scoring target is
def targets(t):
    g = t.groupby("domain_number")
    idx = range(1, N_DOMAINS + 1)
    raw = {
        "current_rate": domain_lrr((PERIOD, mvo.end_year(PERIOD))),
        "alt1_rate_1996_2026": domain_lrr(LONG_WINDOW),
        "alt2_start_plus_16x_lrr_2010_2026": domain_lrr((PERIOD, mvo.end_year(PERIOD))) * SCALE_YEARS,
        "alt3_1996_start_plus_30x_lrr_1996_2026": g["alt3_change_vs_2010_start_m"].mean().reindex(idx),
        "alt4_mean_2025": g["alt4_change_vs_2010_start_m"].mean().reindex(idx),
        "alt5_start_plus_16x_lrr_1996_2026": domain_lrr(LONG_WINDOW) * SCALE_YEARS,
    }
    # A fitted rate is smoothed as the scoring target; a position as the model is
    long_t, own_t = common.coastsat_target(*LONG_WINDOW), common.coastsat_target(PERIOD)
    smooth = {
        "current_rate": own_t,
        "alt1_rate_1996_2026": long_t,
        "alt2_start_plus_16x_lrr_2010_2026": own_t * SCALE_YEARS,
        "alt3_1996_start_plus_30x_lrr_1996_2026": mvo.smoothed(raw["alt3_1996_start_plus_30x_lrr_1996_2026"]),
        "alt4_mean_2025": mvo.smoothed(raw["alt4_mean_2025"]),
        "alt5_start_plus_16x_lrr_1996_2026": long_t * SCALE_YEARS,
    }
    for d in (raw, smooth):
        for k in d:
            d[k] = d[k].rename_axis("gis_domain").astype(float)
    return raw, smooth


# A run's modelled LRR and end-minus-start position change, raw and smoothed
def model_series(r):
    rate = mvo.rates(r.run_dir)["lrr_m_yr"]
    pos = mvo.model_end_change(r.run_dir)
    return ({"rate": rate, "position": pos},
            {"rate": mvo.smoothed(rate), "position": mvo.smoothed(pos)})


# Interior bias, RMSE and r of model minus target
def skill(model, obs):
    m, o = common.interior(model), common.interior(obs)
    ok = m.notna() & o.notna()
    d = (m - o)[ok]
    return dict(n=int(ok.sum()), bias=float(d.mean()), rmse=float(np.sqrt((d ** 2).mean())),
                r=float(np.corrcoef(m[ok], o[ok])[0, 1]))


# Figures

# Alongshore axis: zero line, GIS 1-90, grid, village bands
def _axis(ax, half, label_towns):
    ax.axhline(0, color=INK_MUTED, lw=0.6)
    ax.set_xlim(1, N_DOMAINS)
    ax.set_ylim(-half, half)
    ax.grid(axis="y")
    open_frame(ax)
    town_bands(ax, label=label_towns)


# One panel: the alternative target, the model, and on the rate panel the calibration target
def _panel(ax, i, title, obs, model, half, unit, s, ref=None):
    _axis(ax, half, i == 0)
    handles = [Line2D([], [], color=C_OBS, lw=LW_OBS, label=TARGETS[obs.name][1])]
    fmt = ".2f" if unit == "m/yr" else ".1f"
    handles.append(Line2D([], [], color=C_MODEL, lw=LW_MODEL,
                          label=f"CASCADE: bias {s['bias']:+{fmt}} {unit}, "
                                f"RMSE {s['rmse']:{fmt}}, r {s['r']:.2f}"))
    if ref is not None:
        ax.plot(ref.index, ref.values, color=C_CURRENT, lw=0.9, ls=(0, (4, 2)), zorder=4)
        handles.append(Line2D([], [], color=C_CURRENT, lw=0.9, ls=(0, (4, 2)),
                              label=TARGETS[ref.name][1]))
    ax.plot(obs.index, obs.values, color=C_OBS, lw=LW_OBS, zorder=5)
    ax.plot(model.index, model.values, color=C_MODEL, lw=LW_MODEL, zorder=6)
    off = [("observed", mark_offaxis(ax, obs.index, obs.values, half, color=C_OBS)),
           ("model", mark_offaxis(ax, model.index, model.values, half, color=C_MODEL))]
    _title(ax, i, title)
    ax.legend(handles=handles, loc="lower center", fontsize=6.5, frameon=False)
    return offaxis_clause(off, half, unit)


# The five alternatives for one run, stacked
def figure(r, tgt, model, smoothed):
    fig, axes = plt.subplots(len(PANELS), 1, figsize=figsize("double", height=9.4), sharex=True,
                             constrained_layout=True)
    clauses = []
    for i, (key, title) in enumerate(PANELS):
        qty = TARGETS[key][0]
        unit, half = ("m", POS_HALF) if qty == "position" else ("m/yr", RATE_HALF)
        obs = tgt[key].rename(key)
        ref = tgt["current_rate"].rename("current_rate") if qty == "rate" else None
        clauses.append(_panel(axes[i], i, title, obs, model[qty], half, unit,
                              skill(model[qty], obs), ref))
        axes[i].set_ylabel("Position change,\nfrom 2010 start (m)" if qty == "position"
                           else "Change rate,\nLRR (m/yr)")
    structures(axes[-1], label=True)
    for ax in axes[:-1]:
        structures(ax, label=False)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    reading = "7-domain LOWESS, GIS 1-10 raw" if smoothed else "domain means"
    compare_header(fig, f"{mvo.run_title(r)}: CASCADE against five alternative CoastSat "
                        f"targets ({reading}, both sides)")
    tag = "smoothed" if smoothed else "domain_means"
    png = OUT / tag / f"alternative_targets_{mvo.stem(r)}_{tag}.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        f"{mvo.run_title(r)}. Run {r.run_name}, no source/sink in any domain. The 2010-2026 "
        "test run against five alternative observation targets: pale line the alternative, "
        "dark line CASCADE. (a) Alt 1: " + TARGETS["alt1_rate_1996_2026"][2] + ", against the "
        "modelled OLS rate over 2010-2026; dashed grey, the calibration procedure's target. "
        "(b-e) Positions, each drawn as metres seaward of the 2010 start shoreline (the "
        f"CoastSat mean {_S10}, the DEM-centred mean that built the model's offset), against "
        "the modelled end (1 Jan 2026) minus start. "
        "(b) Alt 2: " + TARGETS["alt2_start_plus_16x_lrr_2010_2026"][2] + ". "
        "(c) Alt 3: " + TARGETS["alt3_1996_start_plus_30x_lrr_1996_2026"][2] + ". "
        "(d) Alt 4: " + TARGETS["alt4_mean_2025"][2] + "; calendar 2025 rather than +/-1 yr "
        "of the run end, because CoastSat stops at 2026-01-13, so its centre is about half a "
        "year before the model end. "
        "(e) Alt 5: " + TARGETS["alt5_start_plus_16x_lrr_1996_2026"][2] + ". "
        f"Alts 2 and 5 scale by {SCALE_YEARS} yr, the run length (the methodology text says 15). "
        f"Both sides are {reading}. Seaward positive; scores over the interior GIS 2-89."
        + "".join(clauses) + " " + mvo.OPTION_A))
    return png


# The targets alone: what each alternative asks of the model
def targets_figure(smooth):
    fig, axes = plt.subplots(2, 1, figsize=figsize("double", height=5.4), sharex=True,
                             constrained_layout=True)
    styles = {"current_rate": dict(color=C_CURRENT, lw=1.0, ls=(0, (4, 2))),
              "alt1_rate_1996_2026": dict(color="#1b7f6b", lw=1.4),
              "alt2_start_plus_16x_lrr_2010_2026": dict(color="#e08214", lw=1.4),
              "alt3_1996_start_plus_30x_lrr_1996_2026": dict(color="#7b3294", lw=1.4),
              "alt4_mean_2025": dict(color="#2166ac", lw=1.4),
              "alt5_start_plus_16x_lrr_1996_2026": dict(color="#b2182b", lw=1.4)}
    clauses = []
    for i, (qty, half, unit) in enumerate((("rate", RATE_HALF, "m/yr"), ("position", POS_HALF, "m"))):
        ax = axes[i]
        _axis(ax, half, i == 0)
        handles = []
        for k, (q, label, _) in TARGETS.items():
            if q != qty:
                continue
            ax.plot(smooth[k].index, smooth[k].values, zorder=5, **styles[k])
            clauses.append(offaxis_clause(
                [(label, mark_offaxis(ax, smooth[k].index, smooth[k].values, half,
                                      color=styles[k]["color"]))], half, unit))
            handles.append(Line2D([], [], label=label, **styles[k]))
        ax.legend(handles=handles, loc="lower center", fontsize=6.5, frameon=False)
        _title(ax, i, "Rate" if qty == "rate" else "Position change from the 2010 start")
        ax.set_ylabel("Change rate (m/yr)" if qty == "rate"
                      else "Position change,\nfrom 2010 start (m)")
        structures(ax, label=(i == 1))
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    compare_header(fig, "The 2010-2026 CoastSat targets: calibration method and alternatives "
                        "(7-domain LOWESS, GIS 1-10 raw)")
    png = OUT / "targets_2010_2026.png"
    save(fig, png, dpi=300, close=True)
    record_caption(png, (
        "The observation targets the 2010-2026 test run can be scored against, no model. "
        "(a) Rate: " + "; ".join(f"{v[1]} = {v[2]}" for v in TARGETS.values() if v[0] == "rate")
        + ". (b) Position change from the 2010 start shoreline: "
        + "; ".join(f"{v[1]} = {v[2]}" for v in TARGETS.values() if v[0] == "position")
        + ". All smoothed at 7 domains, GIS 1-10 raw. Seaward positive." + "".join(clauses)))
    return png


# Run: targets, every 2010-2026 matrix run scored, figures for the zeroBE full-management arm
def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tables").mkdir(exist_ok=True)
    t = transects()
    t.round(4).to_csv(OUT / "tables" / "transect_positions.csv", index=False)
    raw, smooth = targets(t)

    runs = mvo.matrix_runs()
    runs = runs[runs.period == PERIOD]
    rows, dom, written = [], {}, []
    for r in runs.itertuples():
        m_raw, m_smooth = model_series(r)
        for smoothed, tgt, mod in ((False, raw, m_raw), (True, smooth, m_smooth)):
            for key, (qty, _, _) in TARGETS.items():
                rows.append(dict(run_name=r.run_name, preset=r.source_sink_preset,
                                 scenario=r.scenario, relocations=bool(r.relocations_enabled),
                                 target=key, quantity=qty, smoothed=smoothed,
                                 **skill(mod[qty], tgt[key])))
        if (r.scenario == FIG_SCENARIO and r.source_sink_preset == FIG_PRESET
                and not r.relocations_enabled):
            for smoothed, tgt, mod in ((False, raw, m_raw), (True, smooth, m_smooth)):
                written.append(figure(r, tgt, mod, smoothed))
                tag = "smoothed" if smoothed else "raw"
                dom[f"model_{FIG_PRESET}_lrr_m_yr_{tag}"] = mod["rate"]
                dom[f"model_{FIG_PRESET}_position_change_m_{tag}"] = mod["position"]
    written.append(targets_figure(smooth))

    for k, (qty, _, _) in TARGETS.items():
        u = "m_yr" if qty == "rate" else "m"
        dom[f"{k}_{u}_raw"] = raw[k]
        dom[f"{k}_{u}_smoothed"] = smooth[k]
    d = pd.DataFrame(dom).reindex(range(1, N_DOMAINS + 1)).rename_axis("gis_domain")
    d.round(3).to_csv(OUT / "tables" / f"domain_values_{PERIOD}_{mvo.end_year(PERIOD)}.csv")
    s = pd.DataFrame(rows).sort_values(["smoothed", "target", "preset", "scenario", "relocations"])
    s.round(3).to_csv(OUT / "scores.csv", index=False)
    runs[["run_name", "source_sink_preset", "scenario", "relocations_enabled"]].to_csv(
        OUT / "runs_used.csv", index=False)

    n25 = t.loc[t["n_2025"] > 0, "n_2025"]
    print(f"2025 passes per transect: median {n25.median():.0f}, min {n25.min()}, "
          f"{int((t['n_2025'] == 0).sum())} of {len(t)} transects with none; "
          f"{int(t['alt3_change_vs_2010_start_m'].isna().sum())} without an alt 3 value")
    fm = s[(s.scenario == FIG_SCENARIO) & ~s.relocations]
    print(fm.drop(columns=["run_name", "scenario", "relocations"]).to_string(
        index=False, float_format=lambda v: f"{v:+.2f}"))
    for p in written:
        print("wrote   ", Path(p).relative_to(_REPO))


if __name__ == "__main__":
    main()
