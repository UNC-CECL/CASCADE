"""Dune line vs shoreline as the island offset, re-run at the option A waves (2026-09-28).

Asked by Hannah on 2026-09-28: does the 09-25 result (shoreline offset matches
the dune line under management, beats it in the natural run) hold at the
waves adopted on 09-27?

    offsets   duneline (1996/duneline/v1) and shoreline (1996/shoreline/v1),
              both metres. There is no 2010 shoreline offset, so 1996-2010 only
    waves     option A: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5 --
              one setting, no sweep (Hannah's choice)
    ends      zeroBE, as on 09-25: option A's edgeBE ends were solved on the
              dune-line offset and would favour it
    scope     natural and full management; 4 runs, relocations and groins off
    score     RAW share of the alongshore variation explained, interior GIS 2-89
              (the score option A was chosen on); smoothed, bias and r beside it

Everything but the settings and the headline score is the 09-25 driver
(HAT_offset_source_comparison), pointed at this study's folder.

WHERE: output/raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py score
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_offset_source_comparison as base  # noqa: E402

common = base.common
TAG = "island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a"
STUDY_DIR = base.grid.RAW_RUNS / "experiments" / TAG
OPTION_A = {"hs": 2.0, "wave_period_s": 7.5, "wave_asymmetry": 0.6}
HEADLINE = 0.5

# The 09-25 functions read these module globals at call time.
base.TAG, base.STUDY_DIR = TAG, STUDY_DIR
base.TABLES_DIR, base.LOGS_DIR, base.FIG = (STUDY_DIR / "tables", STUDY_DIR / "logs",
                                            STUDY_DIR / "figures")
base.BASE, base.HIGH_ANGLE, base.HEADLINE = OPTION_A, (HEADLINE,), HEADLINE

# ---------------------------------------------------------------------------
# EACH OFFSET GRADED ON ITS OWN FEATURE (Hannah, 2026-09-28). A run started
# from the dune line is scored against the dune line's change, a run started
# from the shoreline against the shoreline's; full management, both periods.
#   dune line   net change in metres: the observed dune-line endpoint change
#               (1997->2009, 2009->2023 as surveyed, 11.6 and 14.1 yr) against
#               the model's endpoint change over its 14 calendar years. The
#               interval mismatch is reported, not corrected
#               ([[cascade-period-is-the-calendar-year]]).
#   shoreline   the CoastSat LRR target against the model's LRR, as before.
# The two scores are on different targets and different estimators (Hannah's
# choice), so they rank each offset against its own feature; they are not a
# head-to-head.
PERIODS = (1996, 2010)
OWN_SCENARIO = "full_management"
# The research group's alongshore smoothing range is 7 domains (Hannah,
# 2026-09-28), for the CoastSat target and the dune line alike. The runner and
# the 09-25 helpers still smooth at 10; this study sets 7 here and leaves them.
LOWESS_DOMAINS = 7
SKIP_SOUTHERN = 10
LONG_WINDOW = "1996_2024"


def coastsat_target7(start, window=None):
    """The CoastSat LRR target built as the runner builds it, at 7 domains.
    `window` fits the LRR on another span ("1996_2024" for the long-term rate);
    default the period's own."""
    from site_layer.hat_observed_rates import COASTSAT_LRR_ROOT
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_lowess import (CoastSatDataset, LowessConfig,
                                                 build_coastsat_series)
    window = window or period_label(start)
    ds = CoastSatDataset(label=f"CoastSat LRR ({window.replace('_', '-')})",
                         period_start=start,
                         csv_path=str(COASTSAT_LRR_ROOT / window / "transect_lrr_full.csv"))
    cfg = LowessConfig(window_domains=(LOWESS_DOMAINS,), skip_southern_domains=SKIP_SOUTHERN)
    cs = build_coastsat_series([ds], active_period_start=start, lowess_config=cfg,
                               domains=HATTERAS_DOMAINS)[0]
    return build_target_table(cs, cfg, HATTERAS_DOMAINS, LOWESS_DOMAINS).set_index(
        "gis_domain")["target_lrr_m_yr"]


def smooth7(series):
    """A per-domain series smoothed as the target is: LOWESS over 7 domains,
    the southern 10 left raw (common.smooth_like_target, at 7)."""
    import numpy as np
    import pandas as pd
    from statsmodels.nonparametric.smoothers_lowess import lowess
    x = series.index.to_numpy(dtype=float)
    y = series.to_numpy(dtype=float)
    ok = np.isfinite(y)
    out = pd.Series(np.nan, index=series.index)
    out[ok] = lowess(y[ok], x[ok], frac=LOWESS_DOMAINS / len(x), return_sorted=False)
    raw = series.index <= SKIP_SOUTHERN
    out[raw] = series[raw]
    return out
OWN_SETTINGS = {**OPTION_A, "wave_angle_high_fraction": HEADLINE}
YEARS = 14


def period_label(start):
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    return f"{start}_{HATTERAS_PERIODS[start]['end_year']}"


def own_log(src, start):
    return (STUDY_DIR / "logs" / period_label(start) / base.group(src, OWN_SCENARIO)
            / f"{base.grid.label(OWN_SETTINGS)}.log")


def own_env(src, start):
    e = base.grid.run_env("x", OWN_SCENARIO, start, OWN_SETTINGS)
    e["HAT_ISLAND_OFFSET_SOURCE"] = src
    e["HAT_RUN_TAG"] = f"{TAG}/runs/{base.group(src, OWN_SCENARIO)}"
    return e


def own_launch(cell):
    import subprocess
    import time
    src, start = cell
    log = own_log(src, start)
    if base.grid.finished(log):
        return
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    p = subprocess.run([sys.executable, str(base.grid.HINDCAST)], env=own_env(src, start),
                       cwd=str(base.grid.PROJECT_ROOT), capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=base.grid.RUN_TIMEOUT_S)
    log.write_text((p.stdout or "") + "\n--- STDERR ---\n" + (p.stderr or ""), encoding="utf-8")
    what = "done" if p.returncode == 0 else f"FAILED ({base.step2.stop_reason(log)})"
    print(f"{what} {src} {period_label(start)} in {(time.perf_counter() - t0) / 60:.1f} min",
          flush=True)


def cmd_run_2010(a):
    """The 2010-2024 full-management pair; 1996's pair is already in runs/."""
    from concurrent.futures import ThreadPoolExecutor
    base.grid.check_barrier3d()
    common.keep_awake()
    todo = [(src, 2010) for src in base.SOURCES if not base.grid.finished(own_log(src, 2010))]
    print(f"{len(todo)} to run", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(own_launch, todo))
    return cmd_grade()


def own_run_dir(src, start):
    import json
    root = (STUDY_DIR / "runs" / base.group(src, OWN_SCENARIO) / period_label(start)
            / base.grid.PRESET)
    hits = []
    for d in sorted(p for p in root.glob("*") if p.is_dir()):
        md = json.loads(next(d.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
        wave = md["wave climate"]
        if not (abs(float(wave["wave_angle_high_frac"]) - HEADLINE) < 1e-9):
            continue
        got = str(md["identity"]["island_offset_version"])
        if not got.startswith(src + "/"):
            raise ValueError(f"{d}: ran on offset {got!r}, filed as {src}")
        hits.append(d)
    if len(hits) != 1:
        raise SystemExit(f"{root}: expected one option A run, found {len(hits)}")
    return hits[0]


def duneline_change(start):
    import pandas as pd
    from site_layer.hat_observed_rates import DUNELINE_ENDPOINT_ROOT
    t = pd.read_csv(DUNELINE_ENDPOINT_ROOT / period_label(start) / "domain_endpoint_summary.csv")
    return t.set_index("domain_number")["mean_change_m"]


def own_profiles():
    """{(source, start): (observed, model, run_dir)} on each offset's own feature."""
    import pandas as pd
    out = {}
    for src in base.SOURCES:
        for start in PERIODS:
            d = own_run_dir(src, start)
            rt = pd.read_csv(d / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
            if src == "duneline":
                # smoothed as the CoastSat target is: LOWESS over 7 domains, the
                # southern 10 raw (Hannah, 2026-09-28); the model stays raw
                out[(src, start)] = (smooth7(duneline_change(start)),
                                     rt.change_rate_m_yr * YEARS, d)
            else:
                out[(src, start)] = (coastsat_target7(start) * YEARS, rt.lrr_m_yr * YEARS, d)
                # PROJECTED shoreline change: the 1996-2024 LRR x 14 yr, the same
                # observed profile in both periods, against the same model runs
                out[("shoreline_projected", start)] = (
                    coastsat_target7(1996, LONG_WINDOW) * YEARS, rt.lrr_m_yr * YEARS, d)
    return out


def cmd_grade(_=None):
    import pandas as pd
    rows = []
    for (src, start), (obs, mod, d) in own_profiles().items():
        s = common.alongshore_scores(mod, obs)
        mi, oi = common.interior(mod), common.interior(obs.reindex(mod.index))
        rows.append(dict(source=src.split("_")[0], period=period_label(start),
                         scenario=OWN_SCENARIO,
                         graded_on={
                             "duneline": "dune-line net change (m), LOWESS 7 domains",
                             "shoreline": "total shoreline change (m): CoastSat LRR of the "
                                          "same period x 14 yr, LOWESS 7 domains",
                             "shoreline_projected": "projected shoreline change (m): CoastSat "
                                                    "LRR 1996-2024 x 14 yr, LOWESS 7 domains",
                         }[src],
                         variance_explained=s["variance_explained"], r=s["r_alongshore"],
                         bias=float((mi - oi).mean()),
                         run_dir=str(d.relative_to(STUDY_DIR)).replace("\\", "/")))
    t = pd.DataFrame(rows)
    (STUDY_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(STUDY_DIR / "tables" / "graded_on_own_feature.csv", index=False)
    print(t.drop(columns="run_dir").round(3).to_string(index=False))
    return 0


def cmd_plot_own(_=None):
    """This study's figures through the shared house_figures (09-25 driver):
    (a) full management 1996-2010, (b) full management 2010-2024."""
    import pandas as pd
    panels = []
    for start in PERIODS:
        a, b = period_label(start).split("_")
        panels.append(dict(
            label=f"Full management {a}–{b}", start=int(a), end=int(b),
            rates={src: pd.read_csv(own_run_dir(src, start) / "tables"
                                    / "shoreline_change_rate.csv").set_index("gis_domain")
                   for src in base.SOURCES}))
    note = ("(a) full management 1996-2010, (b) full management 2010-2024; option A "
            "waves (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5); island offset in "
            "metres (duneline/v1; shoreline/v1 of each start); no source/sink correction at "
            "the ends (zeroBE); relocations and groins off.")
    for p in base.house_figures(panels, base.FIG, "No source/sink correction at the ends",
                                note, "full_management"):
        print(p.relative_to(STUDY_DIR))
    return 0


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=4)
    sub.add_parser("score")
    r = sub.add_parser("run-2010")
    r.add_argument("--jobs", type=int, default=2)
    sub.add_parser("grade")
    sub.add_parser("plot-own")
    a = ap.parse_args()
    return {"run": base.cmd_run, "score": base.cmd_score,
            "run-2010": cmd_run_2010, "grade": cmd_grade, "plot-own": cmd_plot_own}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
