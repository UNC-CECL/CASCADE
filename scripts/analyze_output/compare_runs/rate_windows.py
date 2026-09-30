"""
Model against observation, every hindcast window: CoastSat shoreline and dune line, each with its end-solved run.

    python scripts/analyze_output/compare_runs/rate_windows.py
    python ... --no-sensitivity        # the main level only

Draws the observation in each reading with the model over it, edgeBE, full
management, groin off; writes output/comparisons/model_vs_observed/.
target_comparison.py, smoothing_scale.py and smoothed_lowess7_with_cascade.py
import its loaders. Details: scripts/analyze_output/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import MultipleLocator  # noqa: E402

from cascade_pipeline.coastsat_lowess import (  # noqa: E402
    CoastSatDataset, LowessConfig, build_coastsat_series, compute_domain_means,
    lowess_smooth_transect_to_domains)
from cascade_pipeline.hindcast import build_target_table  # noqa: E402
from cascade_pipeline.run_registry import (  # noqa: E402
    find_run_dir, legacy_arm_to_kind_tag, load_run_index)
from site_layer.hat_observed_rates import (  # noqa: E402
    coastsat_endpoint_csv, dune_endpoint_csv, lrr_csv)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    DOMAIN_AXIS_LABEL, INK, INK_MUTED, apply_style, caption, figsize, save,
    _title,
)
# The 5-scr producers (panel drawing, survey dates), found through scr_paths
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
import coastsat_lrr_windows as obs  # noqa: E402
import coastsat_vs_duneline as dune  # noqa: E402
from site_layer.hat_figure_style import COMPARISONS_ROOT  # noqa: E402


# --- CONFIG ------------------------------------------------------------------
RAW_RUNS = _REPO / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
OUT_DIR = COMPARISONS_ROOT / "model_vs_observed"
PRESET = "edgeBE"
# window -> (run, arm): the option A matrix, ends solved on CoastSat; no metres run for 1984 or 2004
MATRIX_RUNS = {
    (1984, 2004): None,
    (1996, 2010): ("HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin", "calibration"),
    (2004, 2024): None,
    (2010, 2024): ("HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin", "calibration"),
}
WINDOWS = list(MATRIX_RUNS)
# The dune-line end-domain solve on the current setup (history in README)
DUNE_SOLVE_DIR = RAW_RUNS / "experiments" / "end-domain-boundaries/2026-09-29-ends-solved-on-duneline-split12"
# Where each model set's ends were solved
MODEL_SETS = ("coastsat", "dune-mean3", "dune-raw")   # where the ends were solved
MAIN_DUNE = "dune-mean3"
# True: the dune solve is current (option A, 2026-09-27), so its sets are drawn
DUNE_SOLVE_CURRENT = True   # re-solved under option A, 2026-09-27
DRAWN_SETS = MODEL_SETS if DUNE_SOLVE_CURRENT else ("coastsat",)
DUNE_FIG_SET = MAIN_DUNE if DUNE_SOLVE_CURRENT else "coastsat"
INTERIOR = (2, 89)          # the domains the index scores, GIS 2-89
N = obs.N_DOMAINS
Y_LABEL = "Change rate (m/yr)"
# Black: the observation already carries two hues and a fill
C_MODEL = INK
C_CS_TARGET = "#2166ac"     # the house shoreline blue (coastsat_vs_duneline C_LRR)
C_DUNE_TARGET = "#b2182b"   # the house dune red (coastsat_vs_duneline C_DUNE)
LS_SECOND = (0, (4, 2))     # the second model line where two share a panel
NO_RUN_NOTE = "model not yet run for this window"
# The scoring target's LOWESS width, as the runner builds it (10 until 2026-09-28)
TARGET_WINDOW = 7
LOWESS_CONFIG = LowessConfig(window_domains=(TARGET_WINDOW,),
                           skip_southern_domains=10)
SKIP = LOWESS_CONFIG.skip_southern_domains
TARGET_OUTLINE_LW = 0.8    # the edge of the target's fill
RAW_DOT_PT2 = 4.0          # the per-domain means as dots: marker area, ~2 pt across
# variant -> (observation, reading, model column, file stem); OUTPUT_FOLDER says where
VARIANTS = {
    "coastsat/means":             ("coastsat", "means",          "lrr_m_yr",         "model_vs_shoreline_means"),
    "coastsat/lowess":             ("coastsat", "lowess",          "lrr_m_yr",         "model_vs_shoreline_smoothed"),
    "duneline/endpoint":          ("duneline", "endpoint",       "change_rate_m_yr", "model_vs_duneline_netchange"),
    "duneline/endpoint-lowess":    ("duneline", "endpoint-lowess", "change_rate_m_yr", "model_vs_duneline_netchange_smoothed"),
    "both":                       ("both",     "both",           "lrr_m_yr",         "model_vs_shoreline_and_duneline_rate"),
    "both-netchange":             ("both",     "both-netchange", "change_rate_m_yr", "model_vs_shoreline_and_duneline_netchange"),
    "sensitivity/mixed-estimator": ("duneline", "endpoint",      "lrr_m_yr",         "model_ols_vs_duneline_netchange"),
}
OUTPUT_FOLDER = {
    "coastsat/means":              "vs_shoreline/domain_means",
    "coastsat/lowess":              "vs_shoreline/smoothed",
    "duneline/endpoint":           "vs_duneline/endpoint_net_change",
    "duneline/endpoint-lowess":     "vs_duneline/net_change_smoothed",
    "both":                        "vs_shoreline_and_duneline/change_rate",
    "both-netchange":              "vs_shoreline_and_duneline/net_change",
    "sensitivity/mixed-estimator": "sensitivity/mixed-estimator",
}
# Variants drawn in metres: each line x the window's calendar span
NET_CHANGE_VARIANTS = ("both-netchange",)
Y_LABEL_NET = "Net change in position (m)"
NET_PAD_M = 5.0            # the metres bound: largest |change| + this, up to NET_STEP_M
NET_STEP_M = 5.0
COASTSAT_VARIANTS = ("coastsat/means", "coastsat/lowess")
DUNELINE_VARIANTS = ("duneline/endpoint", "duneline/endpoint-lowess")
# Both-panel model lines use the endpoint rate, like both observations
BOTH_COLS = {"coastsat": "change_rate_m_yr", "dune-mean3": "change_rate_m_yr",
             "dune-raw": "change_rate_m_yr"}
# What is drawn: (root under OUT_DIR, variant, model sets in drawing order)
MAIN_PLAN = (
    [("", v, ["coastsat"]) for v in COASTSAT_VARIANTS]
    + [("", v, [DUNE_FIG_SET]) for v in DUNELINE_VARIANTS]
    + [("", v, list(dict.fromkeys(["coastsat", DUNE_FIG_SET])))
       for v in ("both", "both-netchange")]
    + [("", "sensitivity/mixed-estimator", [DUNE_FIG_SET])]
)
SENSITIVITY_PLAN = (
    [("sensitivity/ends-swapped", v, [MAIN_DUNE]) for v in COASTSAT_VARIANTS]
    + [("sensitivity/ends-swapped", v, ["coastsat"]) for v in DUNELINE_VARIANTS]
    + [("sensitivity/dune-raw-solve", v, ["dune-raw"]) for v in DUNELINE_VARIANTS]
) if DUNE_SOLVE_CURRENT else []
TARGET_LABEL = (f"{TARGET_WINDOW}-domain LOWESS (raw means D1–{SKIP})")
TARGET_CLAUSE = (f"a {TARGET_WINDOW}-domain LOWESS of the transect rates north of "
                 f"domain {SKIP}, and the raw domain means over domains 1–{SKIP} "
                 "where the Oregon Inlet boundary dominates")
# -----------------------------------------------------------------------------


# Import a script from the input-prep tree by path
def _import_by_path(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# An empty frame of GIS 1-90
def _full():
    return pd.DataFrame({"domain_number": np.arange(1, N + 1)})


# The dune solve's reading ('mean3', 'raw') from a model-set name, or None
def _smooth_of(model_set):
    return model_set.split("-", 1)[1] if model_set.startswith("dune-") else None


# window -> (run, arm) from the dune solve's solved.csv
def dune_solved_runs(smooth):
    table = pd.read_csv(DUNE_SOLVE_DIR / "solved.csv")
    table = table[table["smooth"] == smooth]
    runs = {}
    for w in WINDOWS:
        hit = table[table["window"] == "{}_{}".format(*w)]
        if hit.empty:
            runs[w] = None
            continue
        r = hit.iloc[-1]
        runs[w] = (str(r["run_name"]), str(r["tag"]))
    return runs


# window -> (run, arm) for one model set
def runs_for(model_set):
    if model_set not in MODEL_SETS:
        raise ValueError(f"model set {model_set!r}; have {MODEL_SETS}")
    if model_set == "coastsat":
        return dict(MATRIX_RUNS)
    return dune_solved_runs(_smooth_of(model_set))


# Both per-domain estimators for one run, and its runs_used.csv row
def load_model(window, spec, model_set, preset=None):
    period = "{}_{}".format(*window)
    if spec is None:
        return None, {"window": period, "model_ends": model_set, "run_name": "",
                      "arm": "", "run_dir": "", "note": NO_RUN_NOTE}
    run_name, arm = spec
    run_dir = find_run_dir(RAW_RUNS, run_name, window, preset or PRESET, arm)
    rates = pd.read_csv(run_dir / "tables" / "shoreline_change_rate.csv")
    rates = rates.rename(columns={"gis_domain": "domain_number"})
    df = _full().merge(rates[["domain_number", "change_rate_m_yr", "lrr_m_yr"]],
                       on="domain_number", how="left")
    row = {"window": period, "model_ends": model_set, "run_name": run_name,
           "arm": arm, "run_dir": str(run_dir.relative_to(_REPO)), "note": ""}
    if RUN_INDEX.is_file():
        # The index is keyed on (run_name, kind, tag); translate the legacy arm name
        idx = load_run_index(RUN_INDEX)
        kind, tag = legacy_arm_to_kind_tag(arm)
        hit = idx[(idx["run_name"] == run_name) & (idx["kind"] == kind)
                  & (idx["tag"] == tag)]
        if len(hit) == 1:
            r = hit.iloc[0]
            for col in ("timestamp", "git_commit", "topo_product",
                        "topo_dune_version", "island_offset_version",
                        "scenario", "groin_enabled", "source_sink_preset",
                        "rate_estimator", "mean_bias_interior_m_yr",
                        "rmse_interior_m_yr"):
                row[col] = r.get(col, "")
        else:
            row["note"] = f"{len(hit)} index rows match (run_name, arm)"
    return df, row


# Every window's run for one model set
def load_models(model_set):
    runs = runs_for(model_set)
    loaded = [load_model(w, runs[w], model_set) for w in WINDOWS]
    return [m for m, _ in loaded], [r for _, r in loaded]


# A LOWESS target table as a GIS 1-90 frame
def _target_frame(series):
    table = build_target_table(series, LOWESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)
    table = table.rename(columns={"gis_domain": "domain_number"})
    return _full().merge(table[["domain_number", "target_lrr_m_yr", "source"]],
                         on="domain_number", how="left")


# The CoastSat scoring target: raw means D1-10, TARGET_WINDOW LOWESS beyond
def load_coastsat_target(window):
    start, end = window
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {start}-{end}", period_start=start,
                         csv_path=str(lrr_csv(start, end)))],
        active_period_start=start, lowess_config=LOWESS_CONFIG,
        domains=HATTERAS_DOMAINS)
    return _target_frame(series[0])


# The stored dune-line endpoint (rate per domain) and its survey metadata
def load_dune_endpoint(window):
    start, end = window
    dom = pd.read_csv(dune_endpoint_csv(start, end, "domain"))
    tr = pd.read_csv(dune_endpoint_csv(start, end, "transect"))
    df = _full().merge(dom[["domain_number", "mean_rate_m_yr"]]
                       .rename(columns={"mean_rate_m_yr": "mean_lrr"}),
                       on="domain_number", how="left")
    df["std_lrr"] = 0.0
    first = tr.iloc[0]
    meta = {"window": "{}_{}".format(*window),
            "start_vintage": int(first["start_vintage"]),
            "end_vintage": int(first["end_vintage"]),
            "start_date": first["start_date"], "end_date": first["end_date"],
            "start_date_assumed": bool(first["start_date_assumed"]),
            "end_date_assumed": bool(first["end_date_assumed"]),
            "interval_yr": float(first["interval_yr"]),
            "n_domains": int(df["mean_lrr"].notna().sum())}
    return df, meta


# CoastSat net change at the dune-line dates, raw and as a target
def load_coastsat_endpoint(window):
    start, end = window
    dom = pd.read_csv(coastsat_endpoint_csv(start, end, "domain"))
    ddf = _full().merge(dom[["domain_number", "mean_rate_m_yr"]]
                        .rename(columns={"mean_rate_m_yr": "mean_lrr"}),
                        on="domain_number", how="left")
    ddf["std_lrr"] = 0.0
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat endpoint {start}-{end}", period_start=start,
                         csv_path=str(coastsat_endpoint_csv(start, end, "transect")),
                         rate_col="rate_m_yr")],
        active_period_start=start, lowess_config=LOWESS_CONFIG,
        domains=HATTERAS_DOMAINS)
    return ddf, _target_frame(series[0])


# The dune-line endpoint per transect, given the scoring target's treatment
def load_dune_endpoint_target(window, meta):
    start, end = window
    t = (pd.read_csv(dune_endpoint_csv(start, end, "transect"))
         .rename(columns={"domain_number": "domain_id", "rate_m_yr": "rate"})
         .sort_values(["domain_id", "line_id"]).reset_index(drop=True))
    rank = t.groupby("domain_id").cumcount()
    n = t.groupby("domain_id")["domain_id"].transform("count")
    sp = HATTERAS_DOMAINS.domain_spacing_m
    along = ((t["domain_id"] - HATTERAS_DOMAINS.first_gis_id) * sp
             + (rank + 0.5) * (sp / n)).to_numpy(dtype=float)
    dom = t["domain_id"].to_numpy(dtype=int)
    rate = t["rate"].to_numpy(dtype=float)
    gis_x, smoothed, frac = lowess_smooth_transect_to_domains(
        along, rate, dom, TARGET_WINDOW, domains=HATTERAS_DOMAINS)
    print(f"  LOWESS applied: window={TARGET_WINDOW} domains  frac={frac:.3f}  "
          f"(dune line {start}-{end}, {len(t)} transects)")
    smooth = dict(zip(gis_x, smoothed))
    raw_x, raw_y = compute_domain_means(dom, rate, HATTERAS_DOMAINS.first_gis_id, SKIP)
    raw = dict(zip(raw_x, raw_y))
    rows = [(g, raw.get(g, np.nan), f"raw mean (D1-{SKIP})") if g <= SKIP
            else (g, smooth.get(g, np.nan), f"LOWESS {TARGET_WINDOW}-dom")
            for g in range(1, N + 1)]
    return pd.DataFrame(rows, columns=["domain_number", "target_lrr_m_yr", "source"])


# Everything observed for one window
class Observation:

    def __init__(self, window):
        self.window = window
        self.coastsat = obs.load_window(*window)
        self.coastsat_target = load_coastsat_target(window)
        self.endpoint, self.meta = load_dune_endpoint(window)
        self.endpoint_target = load_dune_endpoint_target(window, self.meta)
        self.cs_endpoint, self.cs_endpoint_target = load_coastsat_endpoint(window)

    def frames(self, reading):
        """(line/dots frame, target frame or None) for one reading."""
        return {
            "means":          (self.coastsat, None),
            "lowess":          (self.coastsat, self.coastsat_target),
            "endpoint":       (self.endpoint, None),
            "endpoint-lowess": (self.endpoint, self.endpoint_target),
        }[reading]

    # the scoring targets, for tables/skill.csv
    TARGETS = (("coastsat_lowess", lambda o: o.coastsat_target["target_lrr_m_yr"]),
               ("endpoint_raw",   lambda o: o.endpoint["mean_lrr"]),
               ("endpoint_lowess", lambda o: o.endpoint_target["target_lrr_m_yr"]),
               ("cs_endpoint_raw",   lambda o: o.cs_endpoint["mean_lrr"]),
               ("cs_endpoint_lowess", lambda o: o.cs_endpoint_target["target_lrr_m_yr"]))


# The window's calendar span in years
def _span(window):
    return window[1] - window[0]


# Rate -> drawn quantity: 1 for a rate, the span for net change
def _scale(variant, window):
    return float(_span(window)) if variant in NET_CHANGE_VARIANTS else 1.0


# The y-axis label for a variant
def _y_label(variant):
    return Y_LABEL_NET if variant in NET_CHANGE_VARIANTS else Y_LABEL


# One rate half-range for every panel (largest |rate| + 1, rounded up)
def shared_bounds(observations, model_sets):
    frames = [f for o in observations for f in (o.coastsat, o.endpoint, o.cs_endpoint)]
    half = obs.shared_bounds(frames)
    for mdfs, _ in model_sets.values():
        for df in mdfs:
            if df is not None:
                for col in ("change_rate_m_yr", "lrr_m_yr"):
                    m = float(np.nanmax(df[col].abs()))
                    half = max(half, float(math.ceil(m + obs.Y_PAD_M)))
    return half


# The metres half-range for the net-change panels
def shared_bounds_net(observations, model_sets):
    extreme = 0.0
    for i, o in enumerate(observations):
        span = _span(o.window)
        for f in (o.endpoint, o.cs_endpoint):
            extreme = max(extreme, float(np.nanmax(f["mean_lrr"].abs())) * span)
        for mdfs, _ in model_sets.values():
            if mdfs[i] is not None:
                extreme = max(extreme,
                              float(np.nanmax(mdfs[i]["change_rate_m_yr"].abs())) * span)
    return float(math.ceil((extreme + NET_PAD_M) / NET_STEP_M) * NET_STEP_M)


# A y tick giving four to eight intervals across +/-half
def _net_tick(half):
    for step in (5, 10, 20, 25, 50, 100):
        if 2 * half / step <= 8:
            return step
    return 200


# Bias and RMSE of model minus observation over GIS 2-89
def skill(obs_series, mdf, col):
    if mdf is None:
        return np.nan, np.nan, 0
    lo, hi = INTERIOR
    o = pd.Series(np.asarray(obs_series, dtype=float), index=np.arange(1, N + 1))
    m = mdf.set_index("domain_number")[col]
    r = (m.loc[lo:hi] - o.loc[lo:hi]).dropna()
    return float(r.mean()), float(np.sqrt((r ** 2).mean())), int(len(r))


# The model line
def draw_model(ax, df, col, ls="-", scale=1.0):
    ax.plot(df["domain_number"], df[col] * scale, color=C_MODEL, lw=1.3, ls=ls, zorder=8)


# Per-domain means as sign-coloured dots over the target fill
def draw_raw_dots(ax, odf):
    x = odf["domain_number"].to_numpy(dtype=float)
    y = odf["mean_lrr"].to_numpy(dtype=float)
    cols = np.where(y < 0, obs.C_ERODE, obs.C_ACCRETE)
    ax.scatter(x, y, s=RAW_DOT_PT2, c=cols, linewidths=0, zorder=6)


# The 'no run yet' note on an empty panel
def note_no_run(ax, pt):
    ax.text(0.5, 0.93, NO_RUN_NOTE, transform=ax.transAxes, ha="center",
            va="top", fontsize=pt, color=INK_MUTED, style="italic", zorder=9)


# The observation in one reading: line and fill, or target fill and dots
def _draw_observed(ax, reading, ddf, tdf, half, **panel_kw):
    if tdf is None:
        obs.draw_panel(ax, ddf, half, std=(reading == "means"), **panel_kw)
    else:
        obs.draw_panel(ax, ddf, half, std=False, line_lw=0.0,
                       fill_y=tdf["target_lrr_m_yr"].to_numpy(dtype=float),
                       fill_outline_lw=TARGET_OUTLINE_LW, **panel_kw)
        draw_raw_dots(ax, ddf)


# The empty frame, then both targets as lines (scaled for net change)
def _draw_both(ax, o: Observation, half, scale=1.0, **panel_kw):
    blank = _full().assign(mean_lrr=np.nan, std_lrr=0.0)
    obs.draw_panel(ax, blank, half, std=False, line_lw=0.0, **panel_kw)
    ax.plot(o.cs_endpoint_target["domain_number"],
            o.cs_endpoint_target["target_lrr_m_yr"] * scale,
            color=C_CS_TARGET, lw=1.2, zorder=6)
    ax.plot(o.endpoint_target["domain_number"],
            o.endpoint_target["target_lrr_m_yr"] * scale,
            color=C_DUNE_TARGET, lw=1.2, zorder=6)
    if scale != 1.0:
        ax.yaxis.set_major_locator(MultipleLocator(_net_tick(half)))


# One window on one axes: observation, then model line(s); True if any model drawn
def _panel(ax, o: Observation, variant, models, model_keys, half, **panel_kw):
    observation, reading, col, _ = VARIANTS[variant]
    scale = _scale(variant, o.window)
    if observation == "both":
        _draw_both(ax, o, half, scale=scale, **panel_kw)
    else:
        ddf, tdf = o.frames(reading)
        _draw_observed(ax, reading, ddf, tdf, half, **panel_kw)
    i = WINDOWS.index(o.window)
    drawn = False
    for k, key in enumerate(model_keys):
        mdf = models[key][0][i]
        if mdf is not None:
            c = BOTH_COLS[key] if observation == "both" else col
            draw_model(ax, mdf, c, ls="-" if k == 0 else LS_SECOND, scale=scale)
            drawn = True
    return drawn


# Legend clause naming where a model line's ends were solved
def ends_text(model_set):
    sm = _smooth_of(model_set)
    return ", ends solved on CoastSat" if sm is None else \
        f", ends solved on the dune line ({sm})"


# Caption clause naming where a model line's ends were solved
def _ends_clause(model_set):
    sm = _smooth_of(model_set)
    if sm is None:
        return " with the two end domains solved against the CoastSat target"
    return (" with the two end domains solved against the dune-line change "
            f"({sm} reading, endpoint estimator; "
            f"{DUNE_SOLVE_DIR.relative_to(RAW_RUNS).as_posix()}) instead of CoastSat")


# The model estimator in legend words
def _estimator_label(variant):
    if variant in NET_CHANGE_VARIANTS:
        return "net change (last annual shoreline minus first)"
    if VARIANTS[variant][0] == "both":
        return "endpoint rate"
    return "OLS rate" if VARIANTS[variant][2] == "lrr_m_yr" else "endpoint rate"


# One legend entry per row
def add_legend(fig, variant, model_keys):
    observation, reading, _, _ = VARIANTS[variant]
    base = f"modelled shoreline, {_estimator_label(variant)}: edgeBE, full management, no groin"
    dot_pair = (Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=2.3, lw=0),
                Line2D([], [], color=obs.C_ERODE, marker="o", ms=2.3, lw=0))
    fill_pair = (Line2D([], [], color=obs.C_ACCRETE_FILL, lw=6),
                 Line2D([], [], color=obs.C_ERODE_FILL, lw=6))
    line_pair = (Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                 Line2D([], [], color=obs.C_ERODE, lw=1.0))
    if observation == "both" and variant in NET_CHANGE_VARIANTS:
        handles = [Line2D([], [], color=C_CS_TARGET, lw=1.2),
                   Line2D([], [], color=C_DUNE_TARGET, lw=1.2)]
        labels = ["CoastSat shoreline, net change at the dune-line dates, measured, "
                  f"scaled to the window's years, {TARGET_LABEL}",
                  "dune line, net change, measured, scaled to the window's years, "
                  "the same treatment"]
    elif observation == "both":
        handles = [Line2D([], [], color=C_CS_TARGET, lw=1.2),
                   Line2D([], [], color=C_DUNE_TARGET, lw=1.2)]
        labels = [f"CoastSat shoreline, net change at the dune-line dates, {TARGET_LABEL}",
                  "dune line, net change, the same treatment"]
    elif reading == "means":
        handles = [line_pair,
                   Line2D([], [], color=INK_MUTED, lw=0.5, ls=(0, (1, 1.6)))]
        labels = ["observed CoastSat domain mean LRR (accreting / eroding)",
                  "observed ±1 std across the domain's transects"]
    elif reading == "lowess":
        handles = [dot_pair, fill_pair]
        labels = ["observed CoastSat domain mean LRR (accreting / eroding)",
                  f"scoring target: {TARGET_LABEL}"]
    elif reading == "endpoint":
        handles = [line_pair]
        labels = ["observed dune-line change, two surveys (seaward / landward)"]
    else:   # endpoint-lowess
        handles = [dot_pair, fill_pair]
        labels = ["observed dune-line change, two surveys, per domain (seaward / landward)",
                  f"smoothed as the scoring target: {TARGET_LABEL}"]
    for k, key in enumerate(model_keys):
        handles.append(Line2D([], [], color=C_MODEL, lw=1.3,
                              ls="-" if k == 0 else LS_SECOND))
        labels.append(base + ends_text(key))
    fig.legend(handles, labels, loc="outside lower center", ncol=1, frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})


# Caption clause with each window's dune-line vintages, dates and interval
def _dates_clause(metas):
    parts = []
    for m in metas:
        note = []
        if m["start_date_assumed"]:
            note.append(f"the {m['start_vintage']} date assumed")
        if m["end_date_assumed"]:
            note.append(f"the {m['end_vintage']} date assumed")
        parts.append(
            f"{m['window'].replace('_', '–')}: the {m['start_vintage']} line "
            f"({m['start_date']}) to the {m['end_vintage']} line ({m['end_date']}), "
            f"{m['interval_yr']:.2f} yr"
            + (f" ({'; '.join(note)}, mid-year)" if note else ""))
    return "; ".join(parts)


# Caption clause naming each window's run
def _runs_clause(rows):
    return "; ".join(f"{r['window'].replace('_', '–')}: {r['run_name']} "
                     f"({r['arm']} arm)" for r in rows if r["run_name"])


# Caption text describing the observation in one reading
def _observed_clause(observation, reading, metas):
    if reading == "both-netchange":
        return (
            " The two coloured lines are the two observations as net change in "
            "position, each measured between the same two dune-line dates, divided by "
            "that survey interval and multiplied by the window's calendar years, so "
            "they sit on the model's span (the same scaling as target_comparison), "
            f"then given the scoring target's treatment ({TARGET_CLAUSE}). Blue is the "
            "CoastSat shoreline: the mean satellite position within six months of "
            "each dune-line date, differenced (5-scr/3-rates/coastsat/endpoint). Red "
            "is the digitised dune line, end line minus start line "
            f"(5-scr/3-rates/duneline/endpoint). Vintages, dates and survey intervals: "
            f"{_dates_clause(metas)}. Where the interval is shorter than the window the "
            "scaling assumes the same rate over the missing years. The gap between "
            "the two lines is beach-width change, which the model, whose shoreline "
            "is a dune line behind a fixed berm, cannot represent. The same figure "
            "in m/yr is under change_rate/.")
    if observation == "both":
        return (
            " The two coloured lines are the two observations, BOTH AS NET CHANGE "
            "between the same two dune-line dates over the survey interval, each "
            f"given the scoring target's treatment ({TARGET_CLAUSE}). Blue is the "
            "CoastSat shoreline: the mean satellite position within six months of "
            "each dune-line date, differenced (5-scr/3-rates/coastsat/endpoint). Red "
            "is the digitised dune line, end line minus start line "
            f"(5-scr/3-rates/duneline/endpoint). Vintages and dates: {_dates_clause(metas)}. "
            "The gap between them is beach-width change, which the model, whose "
            "shoreline is a dune line behind a fixed berm, cannot represent. The "
            "CoastSat LRR, the model's scoring target, is drawn under vs_shoreline/.")
    if reading == "means":
        return (
            " The observed line is the mean linear regression rate of the CoastSat "
            "transects inside each 500 m domain, blue and filled where the shoreline "
            "moved seaward, red where it moved landward; the dotted lines are ±1 "
            "standard deviation across those transects.")
    if reading == "lowess":
        return (
            " The filled shape is the scoring target the model is graded against, "
            "blue where the shoreline moved seaward and red where it moved landward: "
            f"{TARGET_CLAUSE}, as built by cascade_pipeline.hindcast.build_target_table. "
            "The dots in the same colours are the unsmoothed mean linear regression "
            "rate of the CoastSat transects inside each 500 m domain, one per domain.")
    if reading == "endpoint-lowess":
        return (
            " The filled shape is the dune-line change given the treatment the "
            "CoastSat scoring target gets: the digitised dune line at the window's end "
            "vintage minus the line at its start vintage, per 100 m transect (one "
            "station each, matched between the two lines), divided by the interval "
            f"between the two survey dates, then {TARGET_CLAUSE}, blue where the dune "
            "line moved seaward and red where it moved landward. The dots in the same "
            "colours are the unsmoothed per-domain means (five transects each). "
            f"Vintages and dates: {_dates_clause(metas)}.")
    return (   # endpoint
        " The observed line is the digitised dune line at the window's end vintage "
        "minus the line at its start vintage, per 500 m domain (mean over its ~5 "
        "transects, first station per transect, as the hindcast's end-year target "
        "loader reads the raw offsets), divided by the interval between the two "
        "survey dates, blue and filled where the dune line moved seaward and red "
        f"where it moved landward. Vintages and dates: {_dates_clause(metas)}.")


# The full caption for one figure
def caption_text(windows, rows_by_key, metas, half, grid, variant, model_keys):
    observation, reading, col, _ = VARIANTS[variant]
    wins = ", ".join(f"{a}–{b}" for a, b in windows)
    net = variant in NET_CHANGE_VARIANTS
    estimator = (
        "its net change, the run's last annual shoreline minus its first (the "
        "endpoint rate x the run years, 14 in 1996–2010 and 2010–2024)"
        if net else
        "the endpoint rate, the run's last annual shoreline minus its first over "
        "the run years, like both observations here"
        if observation == "both" else
        "the endpoint rate, the run's last annual shoreline minus its first over "
        "the run years, the like-for-like estimator for two surveys"
        if col == "change_rate_m_yr" else
        "the linear regression rate, the OLS slope over the run's annual "
        "shorelines, the estimator the CoastSat comparison and the run index use")
    what = {"coastsat": "Observed CoastSat shoreline change rate",
            "duneline": "Observed dune-line change",
            "both": "The two scoring targets"}[observation]
    quantity = "net change in shoreline position" if net else "shoreline change rate"
    head = (f"{what} and modelled {quantity} by GIS domain (1 at Cape "
            f"Point, 90 at Pea Island) for {wins}"
            + (": the 1984-start period in the left column, the 1996-start period "
               "in the right, the earlier window of each above the later." if grid
               else "."))
    observed = _observed_clause(observation, reading, metas)
    preset = ("from the edgeBE source/sink preset under full management with the "
              "groin off")
    if len(model_keys) == 2:
        a, b = model_keys
        models = (
            f" The two black lines are the modelled shoreline change of the same "
            f"window, {estimator}, {preset}, each against the target its two end "
            f"domains were solved on: solid{_ends_clause(a)} "
            f"({_runs_clause(rows_by_key[a])}); dashed{_ends_clause(b)} "
            f"({_runs_clause(rows_by_key[b])}).")
    else:
        key = model_keys[0]
        models = (
            f" The black line is the modelled shoreline change of the same window, "
            f"{estimator}, {preset}{_ends_clause(key)}: {_runs_clause(rows_by_key[key])}.")
    missing = sorted({r["window"].replace("_", "–") for k in model_keys
                      for r in rows_by_key[k] if not r["run_name"]})
    body = observed + models
    if missing:
        body += (f" No run exists yet for {', '.join(missing)}; that panel shows "
                 "the observation alone.")
    if observation != "coastsat":
        body += (" A dune line and a shoreline are different features, so a gap "
                 "between the two curves is beach-width change as much as model "
                 "misfit.")
    body += (" Village spans are shaded; the solid hairline is the Buxton groin and "
             "the dotted hairlines are the Avon and Rodanthe piers. ")
    if net:
        body += (f"The y axis is held at ±{half:g} m on every panel, the largest "
                 "|net change| over both observations and every modelled curve of every "
                 f"window plus {NET_PAD_M:g} m rounded up to {NET_STEP_M:g} m, so the "
                 "panels are directly comparable; the 1984–2004 and 2004–2024 windows "
                 "span 20 yr and the other two 14 yr, so equal rates draw larger there.")
    else:
        body += (f"The y axis is held at ±{half:g} m/yr on every panel, the largest "
                 "|rate| over every observed reading and modelled curve of every window "
                 "plus 1 m rounded up, so the panels are directly comparable.")
    return head + body


# The sensitivity arm as a filename token, '' for the main level
def _stem_tag(root):
    return root.rsplit("/", 1)[-1] if root else ""


# The file stem for a variant on a level
def _stem(variant, root):
    stem = VARIANTS[variant][3]
    tag = _stem_tag(root)
    return f"{stem}_{tag}" if tag else stem


# Save and close
def _save(fig, folder, stem):
    out = save(fig, folder / stem, vector=True)
    plt.close(fig)
    return out


# One window, one variant
def single_figure(o: Observation, variant, models, model_keys, half, folder, root=""):
    start, end = o.window
    stem = _stem(variant, root)
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    if not _panel(ax, o, variant, models, model_keys, half):
        note_no_run(ax, 7.5)
    ax.set_title(f"{start}–{end}", loc="center")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(_y_label(variant))
    add_legend(fig, variant, model_keys)
    i = WINDOWS.index(o.window)
    rows_by_key = {k: [models[k][1][i]] for k in model_keys}
    caption(fig, caption_text([o.window], rows_by_key, [o.meta], half, grid=False,
                              variant=variant, model_keys=model_keys))
    return _save(fig, folder, f"{stem}_{start}_{end}")


# All four windows in a 2 x 2 grid, one variant
def grid_figure(observations, variant, models, model_keys, half, folder, root=""):
    stem = _stem(variant, root)
    chains = obs._chains(WINDOWS)
    assert len(chains) == 2 and all(len(c) == 2 for c in chains), chains
    by_w = {o.window: o for o in observations}
    cells = [(r, c, chain[r]) for r in range(2) for c, chain in enumerate(chains)]
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True,
                             figsize=figsize("double", height=5.0),
                             constrained_layout=True)
    for i, (r, c, w) in enumerate(cells):
        ax = axes[r, c]
        if not _panel(ax, by_w[w], variant, models, model_keys, half,
                      label=(i == 0), label_pt=obs.STRUCTURE_LABEL_PT_GRID):
            note_no_run(ax, 6.5)
        _title(ax, i, "{}–{}".format(*w))
        if c > 0:
            ax.tick_params(labelleft=False)
    for ax in axes[-1, :]:
        ax.set_xlabel(DOMAIN_AXIS_LABEL)
    fig.supylabel(_y_label(variant), fontsize=9)
    add_legend(fig, variant, model_keys)
    ordered = [w for _, _, w in cells]
    rows_by_key = {k: [models[k][1][WINDOWS.index(w)] for w in ordered] for k in model_keys}
    caption(fig, caption_text(ordered, rows_by_key, [by_w[w].meta for w in ordered],
                              half, grid=True, variant=variant, model_keys=model_keys))
    return _save(fig, folder, f"{stem}_grid")


# Model-set name as a column token
def _tag(model_set):
    return model_set.replace("-", "")     # coastsat, dunemean3, duneraw


# domain_rates_<w>.csv and skill.csv
def write_tables(observations, models, tables_dir):
    tables_dir.mkdir(parents=True, exist_ok=True)
    skill_rows = []
    for i, o in enumerate(observations):
        tab = _full()
        for name, getter in Observation.TARGETS:
            tab[f"{name}_m_yr"] = getter(o).to_numpy()
        for key, (mdfs, _) in models.items():
            mdf = mdfs[i]
            if mdf is None:
                continue
            tag = _tag(key)
            tab[f"model_endpoint_{tag}_m_yr"] = mdf["change_rate_m_yr"].to_numpy()
            tab[f"model_lrr_{tag}_m_yr"] = mdf["lrr_m_yr"].to_numpy()
            for name, _g in Observation.TARGETS:
                est = "endpoint" if "endpoint" in name else "lrr"
                tab[f"resid_{est}_{tag}_vs_{name}_m_yr"] = (
                    tab[f"model_{est}_{tag}_m_yr"] - tab[f"{name}_m_yr"])
            for est, col in (("endpoint", "change_rate_m_yr"), ("lrr", "lrr_m_yr")):
                for name, getter in Observation.TARGETS:
                    b, e, n = skill(getter(o), mdf, col)
                    skill_rows.append({"window": o.meta["window"], "model_ends": key,
                                       "model_estimator": est, "target": name,
                                       "interval_yr": round(o.meta["interval_yr"], 2),
                                       "n_interior": n, "bias_m_yr": b, "rmse_m_yr": e})
        tab.to_csv(tables_dir / "domain_rates_{}_{}.csv".format(*o.window), index=False)
    skill_df = pd.DataFrame(skill_rows)
    skill_df.to_csv(tables_dir / "skill.csv", index=False)
    return skill_df


# Each target against the runs solved on it, in its own estimator
def fair_rows(skill_df):
    s = skill_df
    return s[((s.model_ends == "coastsat") & (s.target == "coastsat_lowess")
              & (s.model_estimator == "lrr"))
             | ((s.model_ends == DUNE_FIG_SET) & (s.target == "endpoint_lowess")
                & (s.model_estimator == "endpoint"))]


# Run: load observations and runs, fix the y ranges, write tables, draw every figure
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--no-sensitivity", action="store_true",
                    help="draw the main level only")
    args = ap.parse_args(argv)
    plan = MAIN_PLAN + ([] if args.no_sensitivity else SENSITIVITY_PLAN)

    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    observations = [Observation(w) for w in WINDOWS]
    models = {key: load_models(key) for key in DRAWN_SETS}
    half = shared_bounds(observations, models)
    half_net = shared_bounds_net(observations, models)

    # runs_used.csv: one row per window per model set, with the dune line's dates
    prov = []
    for key, (_, rows) in models.items():
        for o, r in zip(observations, rows):
            prov.append({**r, **{k: v for k, v in o.meta.items() if k != "window"}})
    pd.DataFrame(prov).to_csv(OUT_DIR / "runs_used.csv", index=False)
    (OUT_DIR / "y_bounds.txt").write_text(
        f"y axis on every panel: -{half:g} to +{half:g} m/yr\n"
        f"= ceil(max |rate| + {obs.Y_PAD_M:g}) over every observed reading (CoastSat "
        "means, dune-line endpoint) and both estimators of every model set "
        f"({', '.join(DRAWN_SETS)}), " + ", ".join("{}-{}".format(*w) for w in WINDOWS)
        + "\n(the CoastSat std lines are not in the bound)\n"
        f"net-change panels (vs_shoreline_and_duneline/net_change): -{half_net:g} to "
        f"+{half_net:g} m\n= ceil to {NET_STEP_M:g} m of (max |rate x window years| + "
        f"{NET_PAD_M:g}) over both endpoint observations and every model set's "
        "endpoint rate\n", encoding="utf-8")

    # Tables, then every figure in the plan
    skill_df = write_tables(observations, models, OUT_DIR / "tables")

    written = []
    for root, variant, keys in plan:
        folder = OUT_DIR / root / OUTPUT_FOLDER[variant]
        h = half_net if variant in NET_CHANGE_VARIANTS else half
        for o in observations:
            written += single_figure(o, variant, models, keys, h, folder, root)
        written += grid_figure(observations, variant, models, keys, h, folder, root)

    print(f"y bounds  +/-{half:g} m/yr, net change +/-{half_net:g} m")
    for o in observations:
        m = o.meta
        print(f"{m['window']}  dune line {m['start_vintage']} ({m['start_date']}) -> "
              f"{m['end_vintage']} ({m['end_date']})  {m['interval_yr']:.2f} yr")
    for key, (_, rows) in models.items():
        for r in rows:
            print(f"{r['window']}  {key:<11}  {r['run_name'] or '(no run)'}  {r['arm']}")
    print()
    print(fair_rows(skill_df)[["window", "model_ends", "model_estimator", "target",
                               "bias_m_yr", "rmse_m_yr"]]
          .to_string(index=False, float_format=lambda v: f"{v:+.3f}"))
    for p in written:
        print("wrote    ", p.relative_to(_REPO))


if __name__ == "__main__":
    main()
