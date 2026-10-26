"""
Which observation should CASCADE be graded against? Both candidate targets beside the model.

    python scripts/analyze_output/compare_runs/hindcast_vs_observed/target_comparison.py
    python ... --units rate                              # change_rate/ in m/yr
    python ... --coastsat-target total [--units rate]

CoastSat and dune-line targets and three model sets, as net change over each
14-yr window (or as a rate); loaders from rate_windows.py. Writes
output/comparisons/target_comparison/. Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-03
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import rate_windows as rw  # noqa: E402  (loaders, runs, style constants)

_REPO = rw._REPO
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
from coastsat_vs_duneline import beach_width_handles, shade_beach_width  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    COMPARISONS_ROOT, DOMAIN_AXIS_LABEL, INK, _title, apply_style, caption,
    mark_offaxis,
    figsize, save,
)


# --- CONFIG ------------------------------------------------------------------
obs = rw.obs
ROOT_DIR = COMPARISONS_ROOT / "target_comparison"
# CoastSat target: projected = 1996-2024 LRR x 14 in both windows (in use); total_change = each window's own
CS_MODES = {"projected": "projected", "total": "total_change",
            # the pre-2026-09-21 names, kept so old commands still run
            "full": "projected", "subperiod": "total_change"}
CS_MODE_NOUN = {"projected": "Projected shoreline change",
                "total": "Total shoreline change"}
# The same, in rate units
CS_MODE_NOUN_RATE = {"projected": "Long-term shoreline change rate",
                     "total": "Shoreline change rate"}
# Method strings carry the fit window; 'total' is filled in per window
CS_MODE_METHOD = {"projected": "CoastSat LRR 1996–2024 × 14 yr",
                  "total": "CoastSat LRR {}–{} × {} yr"}
CS_MODE_METHOD_RATE = {"projected": "CoastSat LRR 1996–2024",
                       "total": "CoastSat LRR {}–{}"}
# Canonical mode for each alias
CS_CANON = {"full": "projected", "subperiod": "total",
            "projected": "projected", "total": "total"}
FULL_WINDOW = (1996, 2024)
# The end-domain solve against the 1996-2024 LRR on the current setup (history in README)
FULL_SOLVE_DIR = rw.RAW_RUNS / "experiments" / "end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-split12"
CS_MODE = "projected"
OUT_DIR = ROOT_DIR / CS_MODES[CS_MODE]
# 'net' = metres over the window, 'rate' = m/yr under <mode>/change_rate/
UNITS = "net"
RATE_SUBDIR = "change_rate"
Y_HALF_RATE = 8.0          # m/yr, the model_vs_observed rate figures' range
# The zeroBE set: no source/sink term anywhere, ends included
UNSOLVED = "unsolved"
UNSOLVED_PRESET = "zeroBE"
# Window set (rate_windows.WINDOW_SETS) -> {window: (zeroBE run, arm)}, the same matrix cell
UNSOLVED_RUN_SETS = {
    "current": {
        (1996, 2015): ("HAT_1996_2015_zeroBE_offsetmetres_road_bdm_nourish_nogroin", "calibration"),
        (2010, 2026): ("HAT_2010_2026_zeroBE_offsetmetres_road_bdm_nourish_nogroin", "calibration"),
    },
    "14yr": {   # the option A matrix, 2026-09-27
        (1996, 2010): ("HAT_1996_2010_zeroBE_offsetmetres_road_bdm_nogroin", "calibration"),
        (2010, 2024): ("HAT_2010_2024_zeroBE_offsetmetres_road_bdm_nourish_nogroin", "calibration"),
    },
}
ALL_MODEL_SETS = {"coastsat": "ends_solved_on_coastsat",
                  rw.MAIN_DUNE: "ends_solved_on_duneline",
                  UNSOLVED: "ends_unsolved"}
# Without a dune line: one figure per window, both runs against CoastSat
NO_DUNE_FOLDER = "runs_vs_coastsat"
MODEL_LABEL = {"coastsat": "ends solved on CoastSat",
               rw.MAIN_DUNE: "ends solved on the dune line",
               UNSOLVED: "ends not solved"}
# How each set's run is named in a caption
MODEL_CLAUSE = {
    "coastsat": "the edgeBE run with its ends solved on CoastSat",
    rw.MAIN_DUNE: "the edgeBE run with its ends solved on the dune line",
    UNSOLVED: ("the zeroBE run of the same matrix cell, which carries no "
               "source/sink term in any domain including the two ends"),
}
LW = 1.1
# Fixed y range (+/- m) on every figure here; anything beyond is named and marked
Y_HALF_M = 100.0
_COL_NAME = {"coastsat_target_m": "CoastSat target", "dune_target_m": "the total dune line change",
             "model_ends_solved_on_coastsat_m": "CoastSat-solved run",
             "model_ends_solved_on_duneline_m": "dune-solved run",
             "model_ends_unsolved_m": "unsolved run"}
LW_MODEL = 1.4
Y_LABEL = "Net change in position (m)"
Y_LABEL_RATE = "Change rate (m/yr)"
PAIR_KEYS = (("coastsat", "coastsat_target_m", rw.C_CS_TARGET,
              None,
              "CASCADE, ends solved on CoastSat"),
             (rw.MAIN_DUNE, "dune_target_m", rw.C_DUNE_TARGET,
              "Total dune line change (measured, scaled to 14 yr)",
              "CASCADE, ends solved on the dune line"))
LS_MODEL = "-"
# Observed pale and thick, model dark and thin
PALE = {"coastsat": "#92c5de", rw.MAIN_DUNE: "#f4a582"}
LW_TARGET_PALE = 2.8
LW_MODEL_DARK = 1.1
# window -> (domain frame, LOWESS target frame) in projected mode
CS_SOURCE = {}   # window -> (domain frame, LOWESS target frame) when full-period
# -----------------------------------------------------------------------------


# Point this module and rate_windows at one window set
def select_windows(name):
    global WINDOWS, UNSOLVED_RUNS, MODEL_SETS
    rw.select_windows(name)
    UNSOLVED_RUNS = dict(UNSOLVED_RUN_SETS[name])
    WINDOWS = list(UNSOLVED_RUNS)
    MODEL_SETS = {k: v for k, v in ALL_MODEL_SETS.items() if rw.DUNE_LINE or k != rw.MAIN_DUNE}


select_windows(rw.WINDOW_SET)


# True for metres over the window, False for m/yr
def _net():
    return UNITS == "net"


# The unit every number on a figure is in
def _u():
    return "m" if _net() else "m/yr"


# A score in the figure's unit: 0.1 m, or 0.01 m/yr
def _fmt(v, signed=False):
    d = 1 if _net() else 2
    return f"{v:+.{d}f}" if signed else f"{v:.{d}f}"


# ' multiplied by 14 yr' on the metres figures, nothing on the rate ones
def _x14(prefix=" multiplied by"):
    return f"{prefix} 14 yr" if _net() else ""


# The target mode as a filename token, '_rate' added in rate units
def _stem_mode():
    return CS_MODES[CS_MODE] + ("" if _net() else "_rate")


# The y-axis label for the current units
def y_label():
    return Y_LABEL if _net() else Y_LABEL_RATE


# The y tick spacing for a given half-range
def y_tick(half):
    return (20.0 if half > 60 else 10.0) if _net() else 2.0


# What the CoastSat target is called, in the current mode and units
def cs_noun():
    return (CS_MODE_NOUN if _net() else CS_MODE_NOUN_RATE)[CS_MODE]


# The CoastSat method string for one window, fit window filled in
def cs_method(window):
    m = (CS_MODE_METHOD if _net() else CS_MODE_METHOD_RATE)[CS_MODE]
    return m.format(*window, window[1] - window[0]) if CS_MODE == "total" else m


# Legend label for the CoastSat target, named by its fit window
def cs_label():
    return ("CoastSat target — " + cs_noun().lower()
            + (f" ({cs_method(None)})" if CS_MODE == "projected"
               else " (each window's own CoastSat LRR"
               + (" × its years" if _net() else "") + ")"))


# Caption clause describing the CoastSat rate
def cs_clause():
    return ("the full-period 1996–2024 linear regression rate (the same rate in both "
            "windows)" if CS_MODE == "projected" else
            "the window's own linear regression rate (the runner's scoring series)")


# Legend label for the dune-line target
def dune_label():
    return ("Total dune line change (measured, scaled to 14 yr)" if _net()
            else "Dune line change rate (measured over the survey interval)")


# Panel title for the dune-line target in one window
def _dune_title(window):
    return (f"Total dune line change, {window[0]}–{window[1]} (measured, scaled to 14 yr)"
            if _net() else
            f"Dune line change rate, {window[0]}–{window[1]} (measured over the survey interval)")


# Legend label for the target fill
def _target_legend():
    return ("Target, net change over 14 yr (seaward / landward)" if _net()
            else "Target, change rate (seaward / landward)")


# Figure header naming both targets and the dune-line dates
def _targets_line(o):
    m = o.meta
    return ("CoastSat target: " + cs_method(o.window) + "   ·   "
            f"dune line: {m['start_date']} → {m['end_date']}"
            + (" (assumed)" if bool(m.get("end_date_assumed")) else "")
            + f",  {float(m['interval_yr']):.1f} yr,  "
            + (f"measured, scaled to {o.window[1] - o.window[0]} yr" if _net()
               else "measured rate"))


# What the y axis is, in caption words
def _quantity():
    return ("net change in shoreline position over the model window" if _net()
            else "shoreline change rate over the model window")


# What the model line is, in caption words
def _model_quantity():
    return ("net change over the window" if _net() else
            "endpoint rate (last annual shoreline minus first, over the run years)")


# Caption words naming every value beyond the y range
def over_note(frames_cols, half):
    hits = []
    for df, cols, label in frames_cols:
        for col in cols:
            v = df.set_index("domain_number")[col]
            for g, x in v[v.abs() > half].items():
                hits.append(f"{label}{_COL_NAME.get(col, col)} {_fmt(x, True)} "
                            f"{_u()} at GIS {g}")
    return (f" Beyond ±{half:g} {_u()}, off the axis and marked with a triangle at "
            "the edge: " + "; ".join(hits) + "."
            if hits else "")


# window -> (run, tag): the converged step of the 1996-2024 CoastSat edge solve
def full_period_runs():
    log = pd.read_csv(FULL_SOLVE_DIR / "loop_log.csv")
    runs = {}
    for w in rw.WINDOWS:
        hit = log[(log["window"] == w[0]) & log["run_tag"].notna()]
        if w not in WINDOWS or hit.empty:
            runs[w] = None
            continue
        tag = str(hit.iloc[-1]["run_tag"])
        # Tag is <theme>/<study>/<reading>/step<k>: strip the whole study tag
        study_tag = FULL_SOLVE_DIR.relative_to(rw.RAW_RUNS / "experiments").as_posix()
        step = tag[len(study_tag) + 1:] if tag.startswith(study_tag + "/") else tag.split("/", 1)[1]
        d = next((FULL_SOLVE_DIR / step / "{}_{}".format(*w) / "edgeBE").glob("HAT_*"))
        runs[w] = (d.name, tag)
    return runs


# The three model sets' runs (smoothing_scale.py reads the same); set CS_MODE first
def load_model_sets():
    models = {}
    for key in MODEL_SETS:
        if key == UNSOLVED:
            loaded = [rw.load_model(w, UNSOLVED_RUNS.get(w), key, UNSOLVED_PRESET)
                      for w in rw.WINDOWS]
            mdfs, rows = [m for m, _ in loaded], [r for _, r in loaded]
        elif key == "coastsat" and CS_MODE == "projected":
            runs = full_period_runs()
            loaded = [rw.load_model(w, runs[w], key) for w in rw.WINDOWS]
            mdfs, rows = [m for m, _ in loaded], [r for _, r in loaded]
        else:
            mdfs, rows = rw.load_models(key)
        models[key] = ({w: m for w, m in zip(rw.WINDOWS, mdfs)},
                       [r for r in rows if tuple(int(x) for x in r["window"].split("_")) in WINDOWS])
    return models


# Per domain: both targets raw and LOWESS, and each model set's change
def window_values(o, mdfs):
    s, e = o.window
    years = e - s
    df = pd.DataFrame({"domain_number": np.arange(1, rw.N + 1)})
    df["model_years"] = years
    cs, cs_t = CS_SOURCE.get(o.window, (o.coastsat, o.coastsat_target))
    df["coastsat_lrr_m_yr"] = cs["mean_lrr"].to_numpy(float)
    df["coastsat_target_m"] = df["coastsat_lrr_m_yr"] * years
    df["coastsat_target_lowess_m"] = cs_t["target_lrr_m_yr"].to_numpy(float) * years
    if o.endpoint is not None:
        df["dune_interval_yr"] = o.meta["interval_yr"]
        df["dune_rate_m_yr"] = o.endpoint["mean_lrr"].to_numpy(float)
        df["dune_measured_change_m"] = df["dune_rate_m_yr"] * o.meta["interval_yr"]
        df["dune_target_m"] = df["dune_rate_m_yr"] * years
        df["dune_target_lowess_m"] = o.endpoint_target["target_lrr_m_yr"].to_numpy(float) * years
        df["target_difference_m"] = df["coastsat_target_m"] - df["dune_target_m"]
    for key, mdf in mdfs.items():
        col = f"model_{MODEL_SETS[key]}_m"
        df[col] = (mdf["change_rate_m_yr"].to_numpy(float) * years
                   if mdf is not None else np.nan)
        # Smoothed exactly as the target, for the smoothed scores and figures
        df[col[:-2] + "_lowess_m"] = (
            rw.smooth_like_target(df.set_index("domain_number")[col]).to_numpy()
            if mdf is not None else np.nan)
    return df


# Every metres column over the model years (the x 14 undone), names kept
def to_rate(df):
    out = df.copy()
    for c in [c for c in out.columns if c.endswith("_m")]:
        out[c] = out[c] / out["model_years"]
    return out


# Rate-frame columns renamed for the CSV: *_m -> *_m_yr
def _table_names(df):
    return df.rename(columns={c: c[:-2] + "_m_yr" for c in df.columns if c.endswith("_m")})


# Interior bias, RMSE and r of every model set against every target
def skill_rows(window, df):
    lo, hi = rw.INTERIOR
    sel = df[df["domain_number"].between(lo, hi)]
    rows = []
    for key, folder in MODEL_SETS.items():
        for target, col in (("coastsat", "coastsat_target_m"),
                            ("coastsat_lowess", "coastsat_target_lowess_m"),
                            ("duneline", "dune_target_m"),
                            ("duneline_lowess", "dune_target_lowess_m")):
            if col not in sel:    # no dune line for this window
                continue
            # A smoothed target is scored against the model smoothed the same way
            m = sel[f"model_{folder}" + ("_lowess_m" if target.endswith("_lowess") else "_m")]
            r = (m - sel[col]).dropna()
            ok = m.notna() & sel[col].notna()
            rows.append({"window": "{}_{}".format(*window), "model_ends": folder,
                         "target": target, "n": int(len(r)),
                         "bias_m": round(float(r.mean()), 3),
                         "rmse_m": round(float(np.sqrt((r ** 2).mean())), 3),
                         "r": round(float(np.corrcoef(m[ok], sel[col][ok])[0, 1]), 3)})
    return rows


# {(set, window): (GIS 1 term, GIS 90 term, nonzero domains)} from the run index
def end_values(rows_by_key):
    idx = rw.load_run_index(rw.RUN_INDEX)
    out = {}
    for key, rows in rows_by_key.items():
        for r in rows:
            kind, tag = rw.legacy_arm_to_kind_tag(r["arm"])
            h = idx[(idx["run_name"] == r["run_name"]) & (idx["kind"] == kind)
                    & (idx["tag"] == tag)].iloc[0]
            w = tuple(int(x) for x in r["window"].split("_"))
            out[(key, w)] = (float(h["be_rate_gis1_m_yr"]), float(h["be_rate_gis90_m_yr"]),
                             int(h["be_nonzero_domains"]))
    return out


# One panel: both targets, the beach-width band between them, the model line
def draw(ax, o, df, folder, half, label):
    obs.draw_panel(ax, df.assign(mean_lrr=np.nan, std_lrr=0.0), half,
                   label=label, std=False)
    x = df["domain_number"].to_numpy(float)
    cs, du = df["coastsat_target_m"].to_numpy(float), df["dune_target_m"].to_numpy(float)
    shade_beach_width(ax, x, cs, du)
    ax.plot(x, du, color=rw.C_DUNE_TARGET, lw=LW, zorder=11)
    ax.plot(x, cs, color=rw.C_CS_TARGET, lw=LW, zorder=11)
    ax.plot(x, df[f"model_{folder}_m"], color=INK, lw=LW_MODEL, zorder=12)
    mark_offaxis(ax, x, du, half, color=rw.C_DUNE_TARGET)
    mark_offaxis(ax, x, cs, half, color=rw.C_CS_TARGET)
    mark_offaxis(ax, x, df[f"model_{folder}_m"], half, color=INK)
    obs.draw_shoals(ax, label=label)
    fills = obs.fills_in(*o.window)
    if fills:
        obs.draw_fills(ax, fills, half)


# Lift the centred title clear of the fill bars
def _pad_title(ax, window):
    if obs.fills_in(*window):
        ax.set_title(ax.get_title(loc="center"), loc="center", pad=20)


# Both windows for one model set
def figure(observations, frames, key, folder, half, skill_df):
    fig, axes = plt.subplots(len(WINDOWS), 1, sharex=True, sharey=True,
                             constrained_layout=True, figsize=figsize("double", height=5.6))
    for i, (ax, o) in enumerate(zip(axes, observations)):
        draw(ax, o, frames[o.window], folder, half, label=(i == 0))
        ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
        _title(ax, i, "{}, {}–{} ({})".format(
            cs_noun(), *o.window, cs_method(o.window)))
        _pad_title(ax, o.window)
    axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(y_label(), fontsize=9)
    handles = [Line2D([], [], color=rw.C_CS_TARGET, lw=LW, label=cs_label()),
               Line2D([], [], color=rw.C_DUNE_TARGET, lw=LW, label=dune_label()),
               Line2D([], [], color=INK, lw=LW_MODEL,
                      label=f"CASCADE ({MODEL_LABEL[key]})")] + beach_width_handles()
    # Two columns: at three, the long target labels ran off the page
    fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
    sk = skill_df[skill_df["model_ends"] == folder]
    def _s(w, t):
        r = sk[(sk["window"] == w) & (sk["target"] == t)].iloc[0]
        return f"bias {_fmt(r['bias_m'], True)} {_u()}, RMSE {_fmt(r['rmse_m'])} {_u()}"
    stats = " ".join(
        f"{w.replace('_', '–')}: against CoastSat {_s(w, 'coastsat')}; against the "
        f"dune line {_s(w, 'duneline')}." for w in ("1996_2010", "2010_2024"))
    metas = {o.window: o.meta for o in observations}
    dates = "; ".join(f"{m['start_date']} to {m['end_date']} ({m['interval_yr']:.1f} yr) "
                      f"for {w[0]}–{w[1]}" for w, m in metas.items())
    caption(fig, (
        "The two candidate targets and the CASCADE hindcast as "
        + _quantity().replace("the model", "each model") + " by GIS domain (1 at Cape "
        "Point, 90 at Pea Island), seaward positive, domain means. Blue: the "
        f"CoastSat target, {cs_clause()}{_x14(', multiplied by')}. Red: the dune-line target, the measured net "
        f"change between the digitized lines ({dates}; the 2023 date assumed) "
        f"divided by its interval{_x14(' and multiplied by')}. Black: the model's own "
        f"{_model_quantity()}, {MODEL_CLAUSE[key]}, "
        "full management, groin off. The space between the two targets is the "
        "beach-width change they imply: solid grey where the beach widened, hatched "
        "where it narrowed. Interior GIS 2–89, model minus target: " + stats
        + " The LOWESS-smoothed scores, as the runs are graded, are in "
        f"tables/skill.csv. One y axis, ±{half:g} {_u()}, the same on every figure here."
        + over_note([(frames[o.window], ["coastsat_target_m", "dune_target_m",
                                         f"model_{folder}_m"], "{}–{} ".format(*o.window))
                     for o in observations], half)))
    # Stem carries the target mode and the model set, so no two figures share a name
    out = save(fig, OUT_DIR / folder
               / f"target_comparison_{_stem_mode()}_{folder}_1996_2010_2024")
    plt.close(fig)
    return out


# Each target with the run solved on it, one figure per window
def paired_figure(observations, frames, half, skill_df, ends, smoothed=False):
    sfx = "_lowess" if smoothed else ""
    tgt = {"coastsat": "coastsat" + sfx, rw.MAIN_DUNE: "duneline" + sfx}
    pair_keys = {k for k, *_ in PAIR_KEYS}
    for (key, w), (_, _, n) in ends.items():
        if key in pair_keys and n != 2:
            raise SystemExit(f"{key} {w}: {n} nonzero source/sink domains, expected "
                             "the two ends only; the caption would be wrong")
    sk = skill_df.set_index(["window", "model_ends", "target"])
    names = {"coastsat": ("CoastSat target", "coastsat", "ends_solved_on_coastsat"),
             rw.MAIN_DUNE: ("Total dune line change", "duneline", "ends_solved_on_duneline")}
    out = []
    for o in observations:
        w = "{}_{}".format(*o.window)
        df = frames[o.window]
        x = df["domain_number"].to_numpy(float)
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,
                                 figsize=figsize("double", height=5.6))
        for i, (ax, (key, col, _, _, _)) in enumerate(zip(axes, PAIR_KEYS)):
            fill_col = col.replace("_m", "_lowess_m") if smoothed else col
            obs.draw_panel(ax, df.assign(mean_lrr=df[fill_col], std_lrr=0.0), half,
                           label=(i == 0), std=False)
            if smoothed:
                raw = df[col].to_numpy(float)
                ax.scatter(x, raw, s=9, lw=0, alpha=0.8, zorder=11,
                           c=np.where(raw < 0, obs.C_ERODE, obs.C_ACCRETE))
            obs.draw_shoals(ax, label=(i == 0))
            ax.plot(x, df[f"model_{MODEL_SETS[key]}_m"], color=INK, lw=LW_MODEL, zorder=12)
            mark_offaxis(ax, x, df[fill_col], half)
            mark_offaxis(ax, x, df[f"model_{MODEL_SETS[key]}_m"], half, color=INK)
            ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
            _title(ax, i, (f"{cs_noun()}, {o.window[0]}–{o.window[1]} "
                           f"({cs_method(o.window)}), and its run"
                           if key == "coastsat" else
                           _dune_title(o.window) + ", and its run"))
        fills = obs.fills_in(*o.window)
        if fills:
            obs.draw_fills(axes[0], fills, half)
            _pad_title(axes[0], o.window)
        axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
        fig.supylabel(y_label(), fontsize=9)
        handles = [(Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                    Line2D([], [], color=obs.C_ERODE, lw=1.0))]
        labels = [("Target as graded (raw GIS 1–10, LOWESS beyond)"
                   if smoothed else _target_legend())]
        if smoothed:
            handles.append((Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=3, lw=0),
                            Line2D([], [], color=obs.C_ERODE, marker="o", ms=3, lw=0)))
            labels.append("Raw domain means")
        handles.append(Line2D([], [], color=INK, lw=LW_MODEL))
        labels.append("CASCADE, ends solved on that target")
        fig.legend(handles=handles, labels=labels, loc="outside lower center",
                   ncol=len(handles), frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        (c1, c90, _), (d1, d90, _) = ends[("coastsat", o.window)], ends[(rw.MAIN_DUNE, o.window)]
        fig.suptitle(_targets_line(o) + chr(10)
                     + "{}–{}: source/sink correction at the end domains (GIS 1 and 90) "
                       "only; none on GIS 2–89".format(*o.window) + chr(10)
                     + f"end terms GIS 1 / 90 (m/yr): (a) {c1:+.1f} / {c90:+.1f}   "
                       f"(b) {d1:+.1f} / {d90:+.1f}", fontsize=8.5)
        m = o.meta
        caption(fig, (
            f"{o.window[0]}–{o.window[1]}: each candidate target with the CASCADE run "
            "calibrated to it, as " + _quantity() + " "
            "by GIS domain (1 at Cape Point, 90 at Pea Island), seaward positive; "
            + ("the targets SMOOTHED as the runs are graded: the raw domain means over "
               f"GIS 1–10 and a {rw.TARGET_WINDOW}-domain LOWESS of the transect values beyond, drawn as the "
               "fill, with the raw domain means as dots. " if smoothed else "domain means. ")
            + "(a) The CoastSat target, " + cs_clause() + _x14(" x") + ", "
            "as the fill (blue seaward, red landward), and the edgeBE run whose "
            "two end domains were solved against it (black). (b) The dune-line target, "
            f"the measured net change between the digitized lines ({m['start_date']} to "
            f"{m['end_date']}, {m['interval_yr']:.1f} yr"
            + (", the end date assumed" if m['end_date_assumed'] else "")
            + ") divided by its interval" + _x14(" and x") + ", as the fill, and the edgeBE run whose "
            "end domains were solved against the dune line (black). The model lines are "
            "each run's own " + ("net change" if _net() else "endpoint rate") + ", unchanged; the gap between line and fill is the "
            "misfit. THE ONLY SOURCE/SINK CORRECTION IN EITHER RUN IS THE BOUNDARY TERM AT "
            "GIS 1 AND GIS 90, whose values are in the figure title; every domain from GIS "
            "2 to 89 carries none, so the interior is the model's own response. Full "
            "management, groin off. Interior GIS 2–89, model minus its own "
            + ("smoothed " if smoothed else "") + "target: (a) "
            "{} {u} bias, {} {u} RMSE; (b) {} {u} bias, {} {u} RMSE. The y axis "
            "(±{:g} {u}) is the same on every {e}figure in target_comparison.{} Scores against the other "
            "target, raw and smoothed, are in tables/skill.csv.".format(
                _fmt(sk.loc[(w, "ends_solved_on_coastsat", tgt["coastsat"]), "bias_m"], True),
                _fmt(sk.loc[(w, "ends_solved_on_coastsat", tgt["coastsat"]), "rmse_m"]),
                _fmt(sk.loc[(w, "ends_solved_on_duneline", tgt[rw.MAIN_DUNE]), "bias_m"], True),
                _fmt(sk.loc[(w, "ends_solved_on_duneline", tgt[rw.MAIN_DUNE]), "rmse_m"]), half,
                over_note([(df, ["coastsat_target_m", "model_ends_solved_on_coastsat_m",
                                 "dune_target_m", "model_ends_solved_on_duneline_m"], "")],
                          half), u=_u(), e="" if _net() else "m/yr ")))
        out += save(fig, OUT_DIR / ("paired_smoothed" if smoothed else "paired")
                    / (f"target_and_own_run_{_stem_mode()}_{w}"
                       f"{'_smoothed' if smoothed else ''}"))
        plt.close(fig)
    return out


# The unsolved zeroBE run against both targets, one figure per window
def unsolved_figure(observations, frames, half, skill_df, ends, smoothed=False):
    sfx = "_lowess" if smoothed else ""
    col_model = f"model_{MODEL_SETS[UNSOLVED]}_m"
    for w_, (_, _, n) in ((k[1], v) for k, v in ends.items() if k[0] == UNSOLVED):
        if n != 0:
            raise SystemExit(f"unsolved {w_}: {n} nonzero source/sink domains, "
                             "expected none; the caption would be wrong")
    sk = skill_df.set_index(["window", "model_ends", "target"])
    out = []
    for o in observations:
        w = "{}_{}".format(*o.window)
        df = frames[o.window]
        x = df["domain_number"].to_numpy(float)
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,
                                 figsize=figsize("double", height=5.6))
        for i, (ax, (key, col, _, _, _)) in enumerate(zip(axes, PAIR_KEYS)):
            fill_col = col.replace("_m", "_lowess_m") if smoothed else col
            obs.draw_panel(ax, df.assign(mean_lrr=df[fill_col], std_lrr=0.0), half,
                           label=(i == 0), std=False)
            if smoothed:
                raw = df[col].to_numpy(float)
                ax.scatter(x, raw, s=9, lw=0, alpha=0.8, zorder=11,
                           c=np.where(raw < 0, obs.C_ERODE, obs.C_ACCRETE))
            obs.draw_shoals(ax, label=(i == 0))
            ax.plot(x, df[col_model], color=INK, lw=LW_MODEL, zorder=12)
            mark_offaxis(ax, x, df[fill_col], half)
            mark_offaxis(ax, x, df[col_model], half, color=INK)
            ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
            _title(ax, i, (f"{cs_noun()}, {o.window[0]}–{o.window[1]} "
                           f"({cs_method(o.window)})"
                           if key == "coastsat" else _dune_title(o.window)))
        fills = obs.fills_in(*o.window)
        if fills:
            obs.draw_fills(axes[0], fills, half)
            _pad_title(axes[0], o.window)
        axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
        fig.supylabel(y_label(), fontsize=9)
        handles = [(Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                    Line2D([], [], color=obs.C_ERODE, lw=1.0))]
        labels = [("Target as graded (raw GIS 1–10, LOWESS beyond)"
                   if smoothed else _target_legend())]
        if smoothed:
            handles.append((Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=3, lw=0),
                            Line2D([], [], color=obs.C_ERODE, marker="o", ms=3, lw=0)))
            labels.append("Raw domain means")
        handles.append(Line2D([], [], color=INK, lw=LW_MODEL))
        labels.append("CASCADE, no source/sink anywhere (the same run in both panels)")
        fig.legend(handles=handles, labels=labels, loc="outside lower center",
                   ncol=len(handles), frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        fig.suptitle(_targets_line(o) + chr(10)
                     + "{}–{}: NO source/sink correction in any domain, the two ends "
                     "included".format(*o.window) + chr(10)
                     + "the same run in both panels; neither target was fitted",
                     fontsize=8.5)
        m = o.meta
        caption(fig, (
            f"{o.window[0]}–{o.window[1]}: the UNCALIBRATED CASCADE run against both "
            "candidate targets, as " + _quantity() + " "
            "by GIS domain (1 at Cape Point, 90 at Pea Island), seaward positive; "
            + ("the targets SMOOTHED as the runs are graded: the raw domain means over "
               f"GIS 1–10 and a {rw.TARGET_WINDOW}-domain LOWESS of the transect values beyond, drawn as the "
               "fill, with the raw domain means as dots. " if smoothed else "domain means. ")
            + "(a) The CoastSat target, " + cs_clause() + _x14(" x") + ", as the fill (blue "
            "seaward, red landward). (b) The dune-line target, the measured net change "
            f"between the digitized lines ({m['start_date']} to {m['end_date']}, "
            f"{m['interval_yr']:.1f} yr"
            + (", the end date assumed" if m['end_date_assumed'] else "")
            + ") divided by its interval" + _x14(" and x") + ", as the fill. THE BLACK LINE IS THE "
            "SAME RUN IN BOTH PANELS: the zeroBE arm of the matrix cell, which carries NO "
            "source/sink term in ANY domain, the two ends included, so every one of the 90 "
            "domains is the model's own response and neither target was fitted. Full "
            "management, groin off. The gap between line and fill is the misfit. Interior "
            "GIS 2–89, model minus " + ("smoothed " if smoothed else "") + "target: (a) "
            "{} {u} bias, {} {u} RMSE; (b) {} {u} bias, {} {u} RMSE. The y axis "
            "(±{:g} {u}) is the same on every {e}figure in target_comparison.{} Scores for every "
            "model set against every target, raw and smoothed, are in tables/skill.csv.".format(
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "coastsat" + sfx), "bias_m"], True),
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "coastsat" + sfx), "rmse_m"]),
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "duneline" + sfx), "bias_m"], True),
                _fmt(sk.loc[(w, MODEL_SETS[UNSOLVED], "duneline" + sfx), "rmse_m"]), half,
                over_note([(df, ["coastsat_target_m", "dune_target_m", col_model], "")],
                          half), u=_u(), e="" if _net() else "m/yr ")))
        out += save(fig, OUT_DIR / MODEL_SETS[UNSOLVED]
                    / (f"unsolved_run_and_targets_{_stem_mode()}_{w}"
                       f"{'_smoothed' if smoothed else ''}"))
        plt.close(fig)
    return out


# No dune line: the edgeBE and zeroBE runs against the CoastSat target, one figure per window
def runs_vs_coastsat_figure(observations, frames, half, skill_df, ends, smoothed=False):
    sfx = "_lowess" if smoothed else ""
    rows = (("coastsat", 2), (UNSOLVED, 0))   # (model set, nonzero source/sink domains expected)
    for (key, w_), (_, _, n) in ends.items():
        expect = dict(rows).get(key)
        if expect is not None and n != expect:
            raise SystemExit(f"{key} {w_}: {n} nonzero source/sink domains, expected "
                             f"{expect}; the caption would be wrong")
    sk = skill_df.set_index(["window", "model_ends", "target"])
    fill_col = "coastsat_target" + sfx + "_m"
    out = []
    for o in observations:
        w = "{}_{}".format(*o.window)
        years = o.window[1] - o.window[0]
        df = frames[o.window]
        x = df["domain_number"].to_numpy(float)
        fig, axes = plt.subplots(2, 1, sharex=True, sharey=True, constrained_layout=True,
                                 figsize=figsize("double", height=5.6))
        c1, c90, _ = ends[("coastsat", o.window)]
        titles = (f"edgeBE, ends solved on CoastSat (GIS 1 / 90: {c1:+.1f} / {c90:+.1f} m/yr)",
                  "zeroBE, no source/sink in any domain")
        for i, (ax, (key, _)) in enumerate(zip(axes, rows)):
            col_model = f"model_{MODEL_SETS[key]}" + ("_lowess_m" if smoothed else "_m")
            obs.draw_panel(ax, df.assign(mean_lrr=df[fill_col], std_lrr=0.0), half,
                           label=(i == 0), std=False)
            if smoothed:
                raw = df["coastsat_target_m"].to_numpy(float)
                ax.scatter(x, raw, s=9, lw=0, alpha=0.8, zorder=11,
                           c=np.where(raw < 0, obs.C_ERODE, obs.C_ACCRETE))
            obs.draw_shoals(ax, label=(i == 0))
            ax.plot(x, df[col_model], color=INK, lw=LW_MODEL, zorder=12)
            mark_offaxis(ax, x, df[fill_col], half)
            mark_offaxis(ax, x, df[col_model], half, color=INK)
            ax.yaxis.set_major_locator(MultipleLocator(y_tick(half)))
            _title(ax, i, titles[i])
        fills = obs.fills_in(*o.window)
        if fills:
            obs.draw_fills(axes[0], fills, half)
            _pad_title(axes[0], o.window)
        axes[-1].set_xlabel(DOMAIN_AXIS_LABEL)
        fig.supylabel(y_label(), fontsize=9)
        handles = [(Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                    Line2D([], [], color=obs.C_ERODE, lw=1.0))]
        labels = [("CoastSat target as graded (raw GIS 1–10, LOWESS beyond)"
                   if smoothed else cs_label())]
        if smoothed:
            handles.append((Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=3, lw=0),
                            Line2D([], [], color=obs.C_ERODE, marker="o", ms=3, lw=0)))
            labels.append("Raw domain means")
        handles.append(Line2D([], [], color=INK, lw=LW_MODEL))
        labels.append("CASCADE, full management, no groin"
                      + (", smoothed the same way" if smoothed else ""))
        fig.legend(handles=handles, labels=labels, loc="outside lower center",
                   ncol=len(handles), frameon=False,
                   handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})
        fig.suptitle(f"{cs_noun()}, {o.window[0]}–{o.window[1]} ({cs_method(o.window)})",
                     fontsize=8.5)

        def _s(key):
            r = sk.loc[(w, MODEL_SETS[key], "coastsat" + sfx)]
            return (f"{_fmt(r['bias_m'], True)} {_u()} bias, {_fmt(r['rmse_m'])} {_u()} RMSE, "
                    f"r {r['r']:.2f}")
        caption(fig, (
            f"{o.window[0]}–{o.window[1]}: the two CASCADE runs of the full-management matrix "
            "cell against the CoastSat target, as " + _quantity() + " by GIS domain (1 at "
            "Cape Point, 90 at Pea Island), seaward positive; "
            + ("the target SMOOTHED as the runs are graded: the raw domain means over GIS "
               f"1–10 and a {rw.TARGET_WINDOW}-domain LOWESS of the transect values beyond, "
               "drawn as the fill, with the raw domain means as dots; the model lines are "
               "smoothed the same way. " if smoothed else "domain means. ")
            + "The fill is the CoastSat target, " + cs_clause()
            + (f", multiplied by {years} yr" if _net() else "")
            + " (blue seaward, red landward). (a) The edgeBE run, whose two end domains "
            "(GIS 1 and 90) carry a source/sink term solved against this target, values in "
            "the panel title; GIS 2–89 carry none. (b) The zeroBE run, no source/sink term in "
            "any domain, the ends included. The model lines are each run's own "
            + ("net change" if _net() else "endpoint rate")
            + (", smoothed as the target is. " if smoothed else ", unchanged. ") + "Full management, "
            "groin off. No dune line was digitized for this window's end year, so the "
            "dune-line target is not drawn. Interior GIS 2–89, model minus "
            + ("smoothed " if smoothed else "") + f"target: (a) {_s('coastsat')}; "
            f"(b) {_s(UNSOLVED)}. The y axis (±{half:g} {_u()}) is the same on every figure "
            "here." + over_note([(df, ["coastsat_target_m", "model_ends_solved_on_coastsat_m",
                                       "model_ends_unsolved_m"], "")], half)))
        out += save(fig, OUT_DIR / NO_DUNE_FOLDER
                    / f"runs_vs_coastsat_{_stem_mode()}_{w}{'_smoothed' if smoothed else ''}")
        plt.close(fig)
    return out


# Run: parse mode and units, build the tables, draw every figure
def main() -> int:
    global CS_MODE, OUT_DIR, UNITS
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", choices=sorted(UNSOLVED_RUN_SETS), default=rw.WINDOW_SET,
                    help="current (default): the configured windows, CoastSat only. "
                         "14yr: 1996 -> 2010 -> 2024 with the dune line, under 14yr_windows/.")
    ap.add_argument("--coastsat-target", choices=sorted(CS_MODES), default=None,
                    help="total (default on current): each window's own LRR. projected "
                         "(default on 14yr): the 1996-2024 LRR x 14 yr, carried onto "
                         "windows it was not fitted on. 'full' and 'subperiod' are the "
                         "pre-2026-09-21 aliases.")
    ap.add_argument("--units", choices=("net", "rate"), default="net",
                    help="net (default): metres over the window. rate: the same "
                         "in m/yr, written to <mode>/change_rate/.")
    args = ap.parse_args()
    select_windows(args.windows)
    CS_MODE = CS_CANON[args.coastsat_target
                       or ("total" if args.windows == "current" else "projected")]
    if CS_MODE == "projected" and not rw.DUNE_LINE:
        raise SystemExit("projected needs the 1996-2024 LRR end solve, which exists on the "
                         "14yr windows only; use --coastsat-target total")
    UNITS = args.units
    sub = "" if args.windows == "current" else f"{args.windows}_windows"
    OUT_DIR = ROOT_DIR / sub / CS_MODES[CS_MODE]
    if not _net():
        OUT_DIR = OUT_DIR / RATE_SUBDIR
    apply_style()
    observations = [rw.Observation(w) for w in WINDOWS]
    models = load_model_sets()
    if CS_MODE == "projected":
        full = (obs.load_window(*FULL_WINDOW), rw.load_coastsat_target(FULL_WINDOW))
        CS_SOURCE.update({o.window: full for o in observations})
    frames, skill = {}, []
    for o in observations:
        df = window_values(o, {k: m[0][o.window] for k, m in models.items()})
        if not _net():
            df = to_rate(df)
        frames[o.window] = df
        skill += skill_rows(o.window, df)
    skill_df = pd.DataFrame(skill)

    tables = OUT_DIR / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    for w, df in frames.items():
        (df if _net() else _table_names(df)).round(3).to_csv(
            tables / "domain_values_{}_{}.csv".format(*w), index=False)
    (skill_df if _net() else skill_df.rename(
        columns={"bias_m": "bias_m_yr", "rmse_m": "rmse_m_yr"})).to_csv(
        tables / "skill.csv", index=False)
    pd.DataFrame([dict(r, model_ends=MODEL_SETS[k]) for k, (_, rows) in models.items()
                  for r in rows]).to_csv(OUT_DIR / "runs_used.csv", index=False)

    # One fixed y range for every figure here; anything beyond it goes in the caption
    half = Y_HALF_M if _net() else Y_HALF_RATE
    written = []
    ends = end_values({k: rows for k, (_, rows) in models.items()})
    if not rw.DUNE_LINE:
        for smoothed in (False, True):
            written += runs_vs_coastsat_figure(observations, frames, half, skill_df, ends,
                                               smoothed=smoothed)
        print(skill_df[skill_df["target"] == "coastsat"].to_string(index=False))
        print(f"\ny axis +/-{half:g} {_u()}")
        for p in written:
            print("wrote   ", Path(p).relative_to(_REPO))
        return 0
    for key, folder in MODEL_SETS.items():
        written += figure(observations, frames, key, folder, half, skill_df)
    written += paired_figure(observations, frames, half, skill_df, ends)
    written += unsolved_figure(observations, frames, half, skill_df, ends)
    if CS_MODE == "projected":   # the smoothed versions, full-period only (Hannah)
        written += paired_figure(observations, frames, half, skill_df, ends, smoothed=True)
        written += unsolved_figure(observations, frames, half, skill_df, ends, smoothed=True)

    print(skill_df[skill_df["target"].isin(["coastsat", "duneline"])].to_string(index=False))
    print(f"\ny axis +/-{half:g} {_u()}")
    for p in written:
        print("wrote   ", Path(p).relative_to(_REPO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
