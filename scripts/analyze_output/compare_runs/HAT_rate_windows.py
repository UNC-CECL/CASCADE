"""
HAT_rate_windows.py
==============================================================================
The four hindcast windows, the model against an observation: the CoastSat
waterline and the digitised dune line, each with the run whose end domains
were solved on it drawn over it. edgeBE, full management, groin OFF.

ONE FOLDER, ONE SCRIPT (Hannah, 2026-09-17: "more organized and clear, and
potentially condensed, especially with the naming")
    Until 2026-09-17 this was two scripts and two output trees,
    observed_vs_modeled_windows/ (CoastSat, 09-15) and
    duneline_vs_modeled_windows/ (the dune line, 09-16), with two
    vocabularies for one axis (domain_means / loess_target beside
    endpoint/raw / endpoint/loess), five stems that all said "model", and the
    dune tree carrying four copies of its own layout under sensitivity/.
    Now: one tree under output/comparisons/model_vs_observed/ (named
    rate_windows/ until 2026-09-18), organised by OBSERVATION and then
    READING, one naming rule for every file, and each sensitivity drawn once.

WHY IT IS HERE AND NOT IN 5-scr
    The observed-only figure lives with the observations
    (data/hatteras_init/5-scr/3-rates/coastsat/lrr/, drawn with the panel code of
    scripts/input_prep/5-scr/CoastSat/coastsat_lrr_windows.py). Once a run's
    curve is on the panel the figure spans runs, and every cross-run figure
    is filed under output/comparisons/ (output/README.md). This script
    IMPORTS the observed drawing from the 5-scr producer rather than copying
    it, so the two cannot drift apart in how the observation is drawn.

THE RUNS (Hannah, 2026-09-15)
    1984-2004   HAT_1984_2004_edgeBE_road_bdm_nogroin, the version-pair/v2 arm
                (topography v2, as asked; the calibration arm is on v1)
    1996-2010   HAT_1996_2010_edgeBE_road_bdm_nogroin, calibration arm
                (offsets v1: the re-digitized 1997 line)
    2004-2024   HAT_2004_2024_edgeBE_road_bdm_nourish_nogroin, calibration arm
                (full management in this period includes the nourishment)
    2010-2024   HAT_2010_2024_edgeBE_road_bdm_nourish_nogroin, calibration arm
                (run 2026-09-16 once the 2009 dune line gave a 2010 offset
                and the end domains were solved; nourishment on, as in
                2004-2024)
    Run folders are resolved through cascade_pipeline.run_registry with the
    arm named explicitly; the registry raises, listing the arms a run IS
    under, rather than guessing.

    Those are the MATRIX runs: two end domains solved against CoastSat, so
    against the CoastSat target the model meets the observation at GIS 1 and
    90 by construction. A dune-line target deserves the same (Hannah,
    2026-09-17: "a fair comparison"), so the dune-line figures draw the runs
    whose ends were solved on the DUNE LINE (experiments/2026-09-16-dune-
    edgesolve/solved.csv, the three-domain-mean reading). Three model sets,
    named for where their ends were solved:
        coastsat     the matrix
        dune-mean3   the dune solve, mean of the end domain and its two
                     inward neighbours (the main dune-line set)
        dune-raw     the dune solve, the end domain's own value (sensitivity)

THE OBSERVATIONS
    coastsat    the CoastSat transect LRR, an OLS slope through ~250
                satellite dates per transect. Two readings:
        means        the per-domain means as the 5-scr figure draws them
                     (sign-coloured line, fill, +/-1 std)
        loess        the SCORING TARGET as the fill: a 10-domain LOESS of the
                     transect rates north of D10 and the raw means over
                     D1-10, as cascade_pipeline.hindcast.build_target_table
                     makes it for the runner; the means as dots over it
    duneline    the digitised dune line, the feature CASCADE's shoreline
                actually is (a dune line behind a fixed berm). Vintages from
                hat_topo_version.DUNE_LINE_FOR_YEAR (1996 reads the 1997
                line, 2010 the 2009, 2024 the 2023); stations from
                2-brie-offset/raw_offsets/<vintage>_duneline_offset_raw.csv
                read as the hindcast's end-year target loader reads them;
                seaward positive. Survey dates from
                duneline_vs_coastsat.KNOWN_SURVEY_DATES; 2023 has no known
                flight date and is centred on 2023-07-01, flagged in every
                caption that uses it. READ FROM the stored product
                5-scr/3-rates/duneline/endpoint/<window>/ (2026-09-18), not
                computed here, so this figure and the stored numbers cannot
                disagree. Two readings:
        endpoint        two surveys differenced, per domain, over the interval
        endpoint-loess  the same per transect, then the scoring target's
                        treatment (LOESS frac 0.111 vs CoastSat's 0.110)
    A third reading, lrr (an OLS through every dune line in the window), was
    RETIRED 2026-09-18 with 3-rates/duneline_lrr (Hannah: "these should not be
    lrr, they would just be endpoint, we are tracking net change").
    both        the two observations as lines on one panel, no fill, BOTH
                AS NET CHANGE between the same two dune-line dates (2026-09-18,
                Hannah): CoastSat blue (3-rates/coastsat/endpoint, the mean
                position within +/-6 months of each date, differenced) and the
                dune line red, each given the scoring target's LOESS treatment;
                a model line per solve, all as the endpoint rate. The CoastSat
                LRR, the model's actual scoring target, is in vs_shoreline/.

THE MODEL LINE
    lrr_m_yr           the OLS slope over the run's annual shorelines, the
                       estimator CoastSat and the run index use
    change_rate_m_yr   the endpoint rate, last annual shoreline minus first,
                       the like-for-like estimator for a two-survey reading
    Each reading is paired with its own estimator; the one pairing that
    mixes them (endpoint observation, OLS model) is kept as a sensitivity.

NAMING   model_vs_<feature>_<reading>_<start>_<end>.png, and _grid for the
    2 x 2 by period (1984-start left, 1996-start right, earlier window above).
    The feature is shoreline (CoastSat) or duneline, the reading means,
    smoothed or netchange (2026-09-18, Hannah: the old coastsat_ / duneline_
    endpoint_ / both_ stems never said a model was being compared). The
    internal variant keys (coastsat/means, ...) are unchanged; OUTPUT_FOLDER
    maps them to the folders below.
    Every figure folder keeps its PDFs and CAPTIONS.md under supporting/,
    as hat_figure_style.save() and caption() put them.

OUTPUT   output/comparisons/model_vs_observed/
    vs_shoreline/domain_means/   model_vs_shoreline_means_<w>.png     ends solved
    vs_shoreline/smoothed/       model_vs_shoreline_smoothed_<w>.png  on CoastSat
    vs_duneline/net_change/      model_vs_duneline_netchange_<w>.png  ends solved
    vs_duneline/net_change_smoothed/
                        model_vs_duneline_netchange_smoothed_<w>.png  on the dune
                                                                      line (mean3)
    vs_shoreline_and_duneline/   model_vs_shoreline_and_duneline_<w>.png
                                                    both targets, both solves
    tables/             domain_rates_<w>.csv   every reading, every model set,
                                               the residual against each
                        skill.csv              bias and RMSE, GIS 2-89, per
                                               window x model set x estimator
                                               x target
    runs_used.csv       one row per window per model set: run, arm, folder,
                        timestamp, commit, topography, offsets, and the dune
                        line's vintages, dates and interval
    y_bounds.txt        the shared y range and the rule behind it
    sensitivity/
        ends-swapped/       vs_shoreline/* on the dune-solved runs,
                            vs_duneline/* on the CoastSat-solved runs: each
                            target against the OTHER solve
        dune-raw-solve/     vs_duneline/* on the raw-reading dune solve
        mixed-estimator/    model_ols_vs_duneline_netchange_<w>.png: the
                            net-change observation against the model's OLS rate

USAGE
    python scripts/analyze_output/compare_runs/HAT_rate_windows.py
    python ... --no-sensitivity        # the main level only

Author: Hannah A. Henry, UNC CECL
==============================================================================
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

from cascade_pipeline.coastsat_loess import (  # noqa: E402
    CoastSatDataset, LoessConfig, build_coastsat_series, compute_domain_means,
    loess_smooth_transect_to_domains)
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


def _import_by_path(name, path):
    """A script in the input-prep tree, not a package; its drawing or its
    readers are what is reused."""
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# The 5-scr producer: panel drawing, colours, chains, the CoastSat reader.
obs = _import_by_path(
    "coastsat_lrr_windows",
    _REPO / "scripts" / "input_prep" / "5-scr" / "CoastSat" / "coastsat_lrr_windows.py")
# The dune-line producer: survey dates by vintage, the per-domain position.
dune = _import_by_path(
    "duneline_vs_coastsat",
    _REPO / "scripts" / "input_prep" / "5-scr" / "duneline_vs_coastsat"
    / "duneline_vs_coastsat.py")

RAW_RUNS = _REPO / "output" / "raw_runs"
RUN_INDEX = RAW_RUNS / "run_index.csv"
from site_layer.hat_figure_style import COMPARISONS_ROOT  # noqa: E402
OUT_DIR = COMPARISONS_ROOT / "model_vs_observed"

PRESET = "edgeBE"
# window -> (run_name, arm): the matrix, ends solved on CoastSat
MATRIX_RUNS = {
    (1984, 2004): ("HAT_1984_2004_edgeBE_road_bdm_nogroin", "version-pair/v2"),
    (1996, 2010): ("HAT_1996_2010_edgeBE_road_bdm_nogroin", "calibration"),
    (2004, 2024): ("HAT_2004_2024_edgeBE_road_bdm_nourish_nogroin", "calibration"),
    (2010, 2024): ("HAT_2010_2024_edgeBE_road_bdm_nourish_nogroin", "calibration"),
}
WINDOWS = list(MATRIX_RUNS)
# The dune-line end-domain solve. 2026-09-18: re-solved on the re-digitized
# lines (1984-2004 carried over from the 09-16 solve, whose lines did not
# change); each row of solved.csv names its own run tag, so a carried-over
# row points back into the 09-16 experiment.
DUNE_SOLVE_DIR = RAW_RUNS / "experiments" / "2026-09-18-dune-edgesolve"
MODEL_SETS = ("coastsat", "dune-mean3", "dune-raw")   # where the ends were solved
MAIN_DUNE = "dune-mean3"

INTERIOR = (2, 89)          # the domains the index scores, GIS 2-89
N = obs.N_DOMAINS
Y_LABEL = "Change rate (m/yr)"

# Black, not the site config's model orange (Hannah, 2026-09-15): the observed
# line already carries two hues and a fill, and a third hue on top read as
# noise. Black sits over both fills and survives greyscale.
C_MODEL = INK
C_CS_TARGET = "#2166ac"     # the house shoreline blue (duneline_vs_coastsat C_LRR)
C_DUNE_TARGET = "#b2182b"   # the house dune red (duneline_vs_coastsat C_DUNE)
LS_SECOND = (0, (4, 2))     # the second model line where two share a panel
NO_RUN_NOTE = "model not yet run for this window"

# The scoring target, as the runner builds it (HAT_hindcast_1984_2024.py
# section 8): one 10-domain LOESS window, raw means over D1-D10.
TARGET_WINDOW = 10
LOESS_CONFIG = LoessConfig(window_domains=(TARGET_WINDOW,),
                           skip_southern_domains=10)
SKIP = LOESS_CONFIG.skip_southern_domains
TARGET_OUTLINE_LW = 0.8    # the edge of the target's fill
RAW_DOT_PT2 = 4.0          # the per-domain means as dots: marker area, ~2 pt across

# variant -> (observation, reading, model column, file stem). The variant is
# the internal key; OUTPUT_FOLDER says where it is written. Folder and stem
# name the comparison in words (Hannah, 2026-09-18): the feature, not the
# data source, and "netchange", not "endpoint".
#   coastsat  means           line + std + fill        lrr_m_yr
#             loess           target fill + dots       lrr_m_yr
#   duneline  endpoint        line + fill              change_rate_m_yr
#             endpoint-loess  target fill + dots       change_rate_m_yr
#   both      both            two target lines         per solve (BOTH_COLS)
VARIANTS = {
    "coastsat/means":             ("coastsat", "means",          "lrr_m_yr",         "model_vs_shoreline_means"),
    "coastsat/loess":             ("coastsat", "loess",          "lrr_m_yr",         "model_vs_shoreline_smoothed"),
    "duneline/endpoint":          ("duneline", "endpoint",       "change_rate_m_yr", "model_vs_duneline_netchange"),
    "duneline/endpoint-loess":    ("duneline", "endpoint-loess", "change_rate_m_yr", "model_vs_duneline_netchange_smoothed"),
    "both":                       ("both",     "both",           "lrr_m_yr",         "model_vs_shoreline_and_duneline"),
    "sensitivity/mixed-estimator": ("duneline", "endpoint",      "lrr_m_yr",         "model_ols_vs_duneline_netchange"),
}
OUTPUT_FOLDER = {
    "coastsat/means":              "vs_shoreline/domain_means",
    "coastsat/loess":              "vs_shoreline/smoothed",
    "duneline/endpoint":           "vs_duneline/net_change",
    "duneline/endpoint-loess":     "vs_duneline/net_change_smoothed",
    "both":                        "vs_shoreline_and_duneline",
    "sensitivity/mixed-estimator": "sensitivity/mixed-estimator",
}
COASTSAT_VARIANTS = ("coastsat/means", "coastsat/loess")
DUNELINE_VARIANTS = ("duneline/endpoint", "duneline/endpoint-loess")
# The "both" panel draws each solve in its own target's estimator.
# Since 2026-09-18 both observations in both/ are NET CHANGE between the same
# two dates (Hannah), so every model line there is the endpoint rate too.
BOTH_COLS = {"coastsat": "change_rate_m_yr", "dune-mean3": "change_rate_m_yr",
             "dune-raw": "change_rate_m_yr"}

# What is drawn: (root under OUT_DIR, variant, model sets in drawing order).
# The main level pairs each target with the runs solved on it.
MAIN_PLAN = (
    [("", v, ["coastsat"]) for v in COASTSAT_VARIANTS]
    + [("", v, [MAIN_DUNE]) for v in DUNELINE_VARIANTS]
    + [("", "both", ["coastsat", MAIN_DUNE]),
       ("", "sensitivity/mixed-estimator", [MAIN_DUNE])]
)
SENSITIVITY_PLAN = (
    [("sensitivity/ends-swapped", v, [MAIN_DUNE]) for v in COASTSAT_VARIANTS]
    + [("sensitivity/ends-swapped", v, ["coastsat"]) for v in DUNELINE_VARIANTS]
    + [("sensitivity/dune-raw-solve", v, ["dune-raw"]) for v in DUNELINE_VARIANTS]
)


def _full():
    return pd.DataFrame({"domain_number": np.arange(1, N + 1)})


def _smooth_of(model_set):
    return model_set.split("-", 1)[1] if model_set.startswith("dune-") else None


def ends_text(model_set):
    """The legend clause naming where a model line's ends were solved."""
    sm = _smooth_of(model_set)
    return ", ends solved on CoastSat" if sm is None else \
        f", ends solved on the dune line ({sm})"


def _ends_clause(model_set):
    """The caption clause for the same."""
    sm = _smooth_of(model_set)
    if sm is None:
        return " with the two end domains solved against the CoastSat target"
    return (" with the two end domains solved against the dune-line change "
            f"({sm} reading, endpoint estimator; experiments/2026-09-16-dune-"
            "edgesolve) instead of CoastSat")


# -----------------------------------------------------------------------------
# the runs
# -----------------------------------------------------------------------------
def dune_solved_runs(smooth):
    """window -> (run_name, arm) from the dune solve's solved.csv, where the
    arm is the experiment tag the registry translates."""
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


def runs_for(model_set):
    if model_set not in MODEL_SETS:
        raise ValueError(f"model set {model_set!r}; have {MODEL_SETS}")
    if model_set == "coastsat":
        return dict(MATRIX_RUNS)
    return dune_solved_runs(_smooth_of(model_set))


def load_model(window, spec, model_set):
    """Both per-domain estimators for one run, and the provenance row for
    runs_used.csv. (None, row) where no run exists."""
    period = "{}_{}".format(*window)
    if spec is None:
        return None, {"window": period, "model_ends": model_set, "run_name": "",
                      "arm": "", "run_dir": "", "note": NO_RUN_NOTE}
    run_name, arm = spec
    run_dir = find_run_dir(RAW_RUNS, run_name, window, PRESET, arm)
    rates = pd.read_csv(run_dir / "tables" / "shoreline_change_rate.csv")
    rates = rates.rename(columns={"gis_domain": "domain_number"})
    df = _full().merge(rates[["domain_number", "change_rate_m_yr", "lrr_m_yr"]],
                       on="domain_number", how="left")
    row = {"window": period, "model_ends": model_set, "run_name": run_name,
           "arm": arm, "run_dir": str(run_dir.relative_to(_REPO)), "note": ""}
    if RUN_INDEX.is_file():
        # The index is keyed on (run_name, kind, tag) since 2026-09-16; the
        # arm named above is the legacy spelling and is translated.
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


def load_models(model_set):
    """[frame or None per window], [provenance row per window]."""
    runs = runs_for(model_set)
    loaded = [load_model(w, runs[w], model_set) for w in WINDOWS]
    return [m for m, _ in loaded], [r for _, r in loaded]


# -----------------------------------------------------------------------------
# the observations
# -----------------------------------------------------------------------------
def _target_frame(series):
    table = build_target_table(series, LOESS_CONFIG, HATTERAS_DOMAINS, TARGET_WINDOW)
    table = table.rename(columns={"gis_domain": "domain_number"})
    return _full().merge(table[["domain_number", "target_lrr_m_yr", "source"]],
                         on="domain_number", how="left")


def load_coastsat_target(window):
    """The per-domain scoring target: raw means D1-10, the 10-domain LOESS
    beyond, from build_target_table on the window's transect_lrr_full.csv."""
    start, end = window
    series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat {start}-{end}", period_start=start,
                         csv_path=str(lrr_csv(start, end)))],
        active_period_start=start, loess_config=LOESS_CONFIG,
        domains=HATTERAS_DOMAINS)
    return _target_frame(series[0])


def load_dune_endpoint(window):
    """Two surveys differenced per domain, seaward positive, from the stored
    product 3-rates/duneline/endpoint/<window>/ (duneline_endpoint.py builds
    it; 2026-09-18). Returns a frame (domain_number / mean_lrr / std_lrr, the
    columns obs.draw_panel expects, holding the RATE in m/yr) and the
    vintages, dates and interval it was built from."""
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


def load_coastsat_endpoint(window):
    """CoastSat NET CHANGE at the dune-line dates (3-rates/coastsat/endpoint,
    the rate over the survey interval): (domain frame, target frame). The
    target goes through the SAME builder as the LRR target, only reading the
    endpoint rate column, so the two differ in the estimator alone."""
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
        active_period_start=start, loess_config=LOESS_CONFIG,
        domains=HATTERAS_DOMAINS)
    return ddf, _target_frame(series[0])


def load_dune_endpoint_target(window, meta):
    """The two-survey rate per transect (from the stored product), then the
    scoring target's treatment. Returns domain_number / target_lrr_m_yr /
    source."""
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
    gis_x, smoothed, frac = loess_smooth_transect_to_domains(
        along, rate, dom, TARGET_WINDOW, domains=HATTERAS_DOMAINS)
    print(f"  LOESS applied: window={TARGET_WINDOW} domains  frac={frac:.3f}  "
          f"(dune line {start}-{end}, {len(t)} transects)")
    smooth = dict(zip(gis_x, smoothed))
    raw_x, raw_y = compute_domain_means(dom, rate, HATTERAS_DOMAINS.first_gis_id, SKIP)
    raw = dict(zip(raw_x, raw_y))
    rows = [(g, raw.get(g, np.nan), f"raw mean (D1-{SKIP})") if g <= SKIP
            else (g, smooth.get(g, np.nan), f"LOESS {TARGET_WINDOW}-dom")
            for g in range(1, N + 1)]
    return pd.DataFrame(rows, columns=["domain_number", "target_lrr_m_yr", "source"])


class Observation:
    """Everything observed for one window: the CoastSat means and target,
    the dune line in its three readings."""

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
            "loess":          (self.coastsat, self.coastsat_target),
            "endpoint":       (self.endpoint, None),
            "endpoint-loess": (self.endpoint, self.endpoint_target),
        }[reading]

    # the scoring targets, for tables/skill.csv
    TARGETS = (("coastsat_loess", lambda o: o.coastsat_target["target_lrr_m_yr"]),
               ("endpoint_raw",   lambda o: o.endpoint["mean_lrr"]),
               ("endpoint_loess", lambda o: o.endpoint_target["target_lrr_m_yr"]),
               ("cs_endpoint_raw",   lambda o: o.cs_endpoint["mean_lrr"]),
               ("cs_endpoint_loess", lambda o: o.cs_endpoint_target["target_lrr_m_yr"]))


def shared_bounds(observations, model_sets):
    """One half-range for every panel: the 5-scr rule (largest |rate| plus
    1 m, rounded up) over every observed reading and both estimators of
    every model set drawn."""
    frames = [f for o in observations for f in (o.coastsat, o.endpoint, o.cs_endpoint)]
    half = obs.shared_bounds(frames)
    for mdfs, _ in model_sets.values():
        for df in mdfs:
            if df is not None:
                for col in ("change_rate_m_yr", "lrr_m_yr"):
                    m = float(np.nanmax(df[col].abs()))
                    half = max(half, float(math.ceil(m + obs.Y_PAD_M)))
    return half


def skill(obs_series, mdf, col):
    """bias and RMSE of model - observation over GIS 2-89."""
    if mdf is None:
        return np.nan, np.nan, 0
    lo, hi = INTERIOR
    o = pd.Series(np.asarray(obs_series, dtype=float), index=np.arange(1, N + 1))
    m = mdf.set_index("domain_number")[col]
    r = (m.loc[lo:hi] - o.loc[lo:hi]).dropna()
    return float(r.mean()), float(np.sqrt((r ** 2).mean())), int(len(r))


# -----------------------------------------------------------------------------
# drawing
# -----------------------------------------------------------------------------
def draw_model(ax, df, col, ls="-"):
    ax.plot(df["domain_number"], df[col], color=C_MODEL, lw=1.3, ls=ls, zorder=8)


def draw_raw_dots(ax, odf):
    """The per-domain means as sign-coloured dots over the target's fill."""
    x = odf["domain_number"].to_numpy(dtype=float)
    y = odf["mean_lrr"].to_numpy(dtype=float)
    cols = np.where(y < 0, obs.C_ERODE, obs.C_ACCRETE)
    ax.scatter(x, y, s=RAW_DOT_PT2, c=cols, linewidths=0, zorder=6)


def note_no_run(ax, pt):
    ax.text(0.5, 0.93, NO_RUN_NOTE, transform=ax.transAxes, ha="center",
            va="top", fontsize=pt, color=INK_MUTED, style="italic", zorder=9)


def _draw_observed(ax, reading, ddf, tdf, half, **panel_kw):
    """The observation: the 5-scr panel as it draws itself (means, with the
    std lines), a line with its fill (endpoint), or with a target frame the
    target as the fill and the per-domain values as dots over it."""
    if tdf is None:
        obs.draw_panel(ax, ddf, half, std=(reading == "means"), **panel_kw)
    else:
        obs.draw_panel(ax, ddf, half, std=False, line_lw=0.0,
                       fill_y=tdf["target_lrr_m_yr"].to_numpy(dtype=float),
                       fill_outline_lw=TARGET_OUTLINE_LW, **panel_kw)
        draw_raw_dots(ax, ddf)


def _draw_both(ax, o: Observation, half, **panel_kw):
    """Axes, bands and structures from the observed panel drawn empty, then
    the two targets as lines. The model lines go on afterwards."""
    blank = _full().assign(mean_lrr=np.nan, std_lrr=0.0)
    obs.draw_panel(ax, blank, half, std=False, line_lw=0.0, **panel_kw)
    ax.plot(o.cs_endpoint_target["domain_number"], o.cs_endpoint_target["target_lrr_m_yr"],
            color=C_CS_TARGET, lw=1.2, zorder=6)
    ax.plot(o.endpoint_target["domain_number"], o.endpoint_target["target_lrr_m_yr"],
            color=C_DUNE_TARGET, lw=1.2, zorder=6)


def _panel(ax, o: Observation, variant, models, model_keys, half, **panel_kw):
    """One window on one axes: the observation in the variant's reading and
    the model line(s) over it. Returns whether any model line was drawn."""
    observation, reading, col, _ = VARIANTS[variant]
    if observation == "both":
        _draw_both(ax, o, half, **panel_kw)
    else:
        ddf, tdf = o.frames(reading)
        _draw_observed(ax, reading, ddf, tdf, half, **panel_kw)
    i = WINDOWS.index(o.window)
    drawn = False
    for k, key in enumerate(model_keys):
        mdf = models[key][0][i]
        if mdf is not None:
            c = BOTH_COLS[key] if observation == "both" else col
            draw_model(ax, mdf, c, ls="-" if k == 0 else LS_SECOND)
            drawn = True
    return drawn


# -----------------------------------------------------------------------------
# legends
# -----------------------------------------------------------------------------
def _estimator_label(variant):
    if VARIANTS[variant][0] == "both":
        return "endpoint rate"
    return "OLS rate" if VARIANTS[variant][2] == "lrr_m_yr" else "endpoint rate"


TARGET_LABEL = (f"{TARGET_WINDOW}-domain LOESS (raw means D1–{SKIP})")


def add_legend(fig, variant, model_keys):
    """One entry per row: the model label carries the estimator and the
    solve, and two abreast ran past the page edge at 190 mm (2026-09-16)."""
    observation, reading, _, _ = VARIANTS[variant]
    base = f"modelled shoreline, {_estimator_label(variant)}: edgeBE, full management, no groin"
    dot_pair = (Line2D([], [], color=obs.C_ACCRETE, marker="o", ms=2.3, lw=0),
                Line2D([], [], color=obs.C_ERODE, marker="o", ms=2.3, lw=0))
    fill_pair = (Line2D([], [], color=obs.C_ACCRETE_FILL, lw=6),
                 Line2D([], [], color=obs.C_ERODE_FILL, lw=6))
    line_pair = (Line2D([], [], color=obs.C_ACCRETE, lw=1.0),
                 Line2D([], [], color=obs.C_ERODE, lw=1.0))
    if observation == "both":
        handles = [Line2D([], [], color=C_CS_TARGET, lw=1.2),
                   Line2D([], [], color=C_DUNE_TARGET, lw=1.2)]
        labels = [f"CoastSat shoreline, net change at the dune-line dates, {TARGET_LABEL}",
                  "dune line, net change, the same treatment"]
    elif reading == "means":
        handles = [line_pair,
                   Line2D([], [], color=INK_MUTED, lw=0.5, ls=(0, (1, 1.6)))]
        labels = ["observed CoastSat domain mean LRR (accreting / eroding)",
                  "observed ±1 std across the domain's transects"]
    elif reading == "loess":
        handles = [dot_pair, fill_pair]
        labels = ["observed CoastSat domain mean LRR (accreting / eroding)",
                  f"scoring target: {TARGET_LABEL}"]
    elif reading == "endpoint":
        handles = [line_pair]
        labels = ["observed dune-line change, two surveys (seaward / landward)"]
    else:   # endpoint-loess
        handles = [dot_pair, fill_pair]
        labels = ["observed dune-line change, two surveys, per domain (seaward / landward)",
                  f"smoothed as the scoring target: {TARGET_LABEL}"]
    for k, key in enumerate(model_keys):
        handles.append(Line2D([], [], color=C_MODEL, lw=1.3,
                              ls="-" if k == 0 else LS_SECOND))
        labels.append(base + ends_text(key))
    fig.legend(handles, labels, loc="outside lower center", ncol=1, frameon=False,
               handler_map={tuple: HandlerTuple(ndivide=None, pad=0.3)})


# -----------------------------------------------------------------------------
# captions
# -----------------------------------------------------------------------------
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


def _runs_clause(rows):
    return "; ".join(f"{r['window'].replace('_', '–')}: {r['run_name']} "
                     f"({r['arm']} arm)" for r in rows if r["run_name"])


TARGET_CLAUSE = (f"a {TARGET_WINDOW}-domain LOESS of the transect rates north of "
                 f"domain {SKIP}, and the raw domain means over domains 1–{SKIP} "
                 "where the Oregon Inlet boundary dominates")


def _observed_clause(observation, reading, metas):
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
    if reading == "loess":
        return (
            " The filled shape is the scoring target the model is graded against, "
            "blue where the shoreline moved seaward and red where it moved landward: "
            f"{TARGET_CLAUSE}, as built by cascade_pipeline.hindcast.build_target_table. "
            "The dots in the same colours are the unsmoothed mean linear regression "
            "rate of the CoastSat transects inside each 500 m domain, one per domain.")
    if reading == "endpoint-loess":
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


def caption_text(windows, rows_by_key, metas, half, grid, variant, model_keys):
    observation, reading, col, _ = VARIANTS[variant]
    wins = ", ".join(f"{a}–{b}" for a, b in windows)
    estimator = (
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
    head = (f"{what} and modelled shoreline change rate by GIS domain (1 at Cape "
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
             "the dotted hairlines are the Avon and Rodanthe piers. The y axis is "
             f"held at ±{half:g} m/yr on every panel, the largest |rate| over every "
             "observed reading and modelled curve of every window plus 1 m rounded "
             "up, so the panels are directly comparable.")
    return head + body


# -----------------------------------------------------------------------------
# figures
# -----------------------------------------------------------------------------
def _save(fig, folder, stem):
    out = save(fig, folder / stem, vector=True)
    plt.close(fig)
    return out


def single_figure(o: Observation, variant, models, model_keys, half, folder):
    start, end = o.window
    stem = VARIANTS[variant][3]
    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.40),
                           constrained_layout=True)
    if not _panel(ax, o, variant, models, model_keys, half):
        note_no_run(ax, 7.5)
    ax.set_title(f"{start}–{end}", loc="center")
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel(Y_LABEL)
    add_legend(fig, variant, model_keys)
    i = WINDOWS.index(o.window)
    rows_by_key = {k: [models[k][1][i]] for k in model_keys}
    caption(fig, caption_text([o.window], rows_by_key, [o.meta], half, grid=False,
                              variant=variant, model_keys=model_keys))
    return _save(fig, folder, f"{stem}_{start}_{end}")


def grid_figure(observations, variant, models, model_keys, half, folder):
    stem = VARIANTS[variant][3]
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
    fig.supylabel(Y_LABEL, fontsize=9)
    add_legend(fig, variant, model_keys)
    ordered = [w for _, _, w in cells]
    rows_by_key = {k: [models[k][1][WINDOWS.index(w)] for w in ordered] for k in model_keys}
    caption(fig, caption_text(ordered, rows_by_key, [by_w[w].meta for w in ordered],
                              half, grid=True, variant=variant, model_keys=model_keys))
    return _save(fig, folder, f"{stem}_grid")


# -----------------------------------------------------------------------------
# tables
# -----------------------------------------------------------------------------
def _tag(model_set):
    return model_set.replace("-", "")     # coastsat, dunemean3, duneraw


def write_tables(observations, models, tables_dir):
    """domain_rates_<w>.csv: every reading of the observation, both
    estimators of every model set, the residual against each. skill.csv:
    bias and RMSE over GIS 2-89 per window x model set x estimator x target."""
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


def fair_rows(skill_df):
    """Each target scored against the runs solved on it, with its own
    estimator: the three columns of the README's skill table."""
    s = skill_df
    return s[((s.model_ends == "coastsat") & (s.target == "coastsat_loess")
              & (s.model_estimator == "lrr"))
             | ((s.model_ends == MAIN_DUNE) & (s.target == "endpoint_loess")
                & (s.model_estimator == "endpoint"))]


# -----------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--no-sensitivity", action="store_true",
                    help="draw the main level only")
    args = ap.parse_args(argv)
    plan = MAIN_PLAN + ([] if args.no_sensitivity else SENSITIVITY_PLAN)

    apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    observations = [Observation(w) for w in WINDOWS]
    models = {key: load_models(key) for key in MODEL_SETS}
    half = shared_bounds(observations, models)

    # provenance: one row per window per model set, with the dune line's
    # vintages, dates and interval beside the run
    prov = []
    for key, (_, rows) in models.items():
        for o, r in zip(observations, rows):
            prov.append({**r, **{k: v for k, v in o.meta.items() if k != "window"}})
    pd.DataFrame(prov).to_csv(OUT_DIR / "runs_used.csv", index=False)
    (OUT_DIR / "y_bounds.txt").write_text(
        f"y axis on every panel: -{half:g} to +{half:g} m/yr\n"
        f"= ceil(max |rate| + {obs.Y_PAD_M:g}) over every observed reading (CoastSat "
        "means, dune-line endpoint) and both estimators of every model set "
        f"({', '.join(MODEL_SETS)}), " + ", ".join("{}-{}".format(*w) for w in WINDOWS)
        + "\n(the CoastSat std lines are not in the bound)\n", encoding="utf-8")

    skill_df = write_tables(observations, models, OUT_DIR / "tables")

    written = []
    for root, variant, keys in plan:
        folder = OUT_DIR / root / OUTPUT_FOLDER[variant]
        for o in observations:
            written += single_figure(o, variant, models, keys, half, folder)
        written += grid_figure(observations, variant, models, keys, half, folder)

    print(f"y bounds  +/-{half:g} m/yr")
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
