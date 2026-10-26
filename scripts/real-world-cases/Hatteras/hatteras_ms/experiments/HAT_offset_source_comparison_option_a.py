"""
Does the 09-25 offset-source result hold at the option A waves?

    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py run --jobs 4
    python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py score

The 09-25 driver at option A, 1996-2010; run-2010, grade and plot-own add the
2010-2024 pair and grade each offset on its own feature. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_offset_source_comparison as base  # noqa: E402

common = base.common
# --- CONFIG ------------------------------------------------------------------
TAG = "island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a"
STUDY_DIR = base.grid.RAW_RUNS / "experiments" / TAG
OPTION_A = {"hs": 2.0, "wave_period_s": 7.5, "wave_asymmetry": 0.6}
HEADLINE = 0.5

# The 09-25 functions read these module globals at call time.
base.TAG, base.STUDY_DIR = TAG, STUDY_DIR
base.TABLES_DIR, base.LOGS_DIR, base.FIG = (STUDY_DIR / "tables", STUDY_DIR / "logs",
                                            STUDY_DIR / "figures")
base.BASE, base.HIGH_ANGLE, base.HEADLINE = OPTION_A, (HEADLINE,), HEADLINE

# Each offset graded on its own feature: dune line on net change, shoreline on LRR
PERIODS = (1996, 2010)
OWN_SCENARIO = "full_management"
# Smoothing at 7 domains, the research group's range; the 09-25 helpers keep 10
LOWESS_DOMAINS = 7
SKIP_SOUTHERN = 10
LONG_WINDOW = "1996_2024"
# -----------------------------------------------------------------------------


# The CoastSat LRR target built as the runner builds it, at 7 domains
def coastsat_target7(start, window=None):
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


# A per-domain series smoothed as the target is
def smooth7(series):
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


# A period's window label, e.g. 1996_2010
def period_label(start):
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    return f"{start}_{HATTERAS_PERIODS[start]['end_year']}"


# The log file for one own-feature run
def own_log(src, start):
    return (STUDY_DIR / "logs" / period_label(start) / base.group(src, OWN_SCENARIO)
            / f"{base.grid.label(OWN_SETTINGS)}.log")


# The runner's environment for one own-feature run
def own_env(src, start):
    e = base.grid.run_env("x", OWN_SCENARIO, start, OWN_SETTINGS)
    e["HAT_ISLAND_OFFSET_SOURCE"] = src
    e["HAT_RUN_TAG"] = f"{TAG}/runs/{base.group(src, OWN_SCENARIO)}"
    return e


# One hindcast run in a subprocess, its output logged; skipped if already finished
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


# The 2010-2024 full-management pair
def cmd_run_2010(a):
    from concurrent.futures import ThreadPoolExecutor
    base.grid.check_barrier3d()
    common.keep_awake()
    todo = [(src, 2010) for src in base.SOURCES if not base.grid.finished(own_log(src, 2010))]
    print(f"{len(todo)} to run", flush=True)
    with ThreadPoolExecutor(max_workers=a.jobs) as pool:
        list(pool.map(own_launch, todo))
    return cmd_grade()


# The one option A run of a source and period, checked against the offset it read
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


# Observed dune-line net change per domain (m) over a period
def duneline_change(start):
    import pandas as pd
    from site_layer.hat_observed_rates import DUNELINE_ENDPOINT_ROOT
    t = pd.read_csv(DUNELINE_ENDPOINT_ROOT / period_label(start) / "domain_endpoint_summary.csv")
    return t.set_index("domain_number")["mean_change_m"]


# (observed, model, run dir) per (source, start), each offset on its own feature
def own_profiles():
    import pandas as pd
    out = {}
    for src in base.SOURCES:
        for start in PERIODS:
            d = own_run_dir(src, start)
            rt = pd.read_csv(d / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
            if src == "duneline":
                # Smoothed as the CoastSat target is; the model stays raw
                out[(src, start)] = (smooth7(duneline_change(start)),
                                     rt.change_rate_m_yr * YEARS, d)
            else:
                out[(src, start)] = (coastsat_target7(start) * YEARS, rt.lrr_m_yr * YEARS, d)
                # Projected shoreline change: the 1996-2024 LRR x 14 yr, against the same runs
                out[("shoreline_projected", start)] = (
                    coastsat_target7(1996, LONG_WINDOW) * YEARS, rt.lrr_m_yr * YEARS, d)
    return out


# Score each offset against its own feature and write the table
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


# This study's figures through the shared house_figures (09-25 driver)
def cmd_plot_own(_=None):
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


# Run: the chosen subcommand
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
