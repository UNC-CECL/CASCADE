"""Four-parameter wave grid, scored on the smoothed model (2026-09-25).

Designed with Hannah on 2026-09-25, after the 2026-09-24 step-2 study
(one parameter at a time plus two 2-D grids) left no search over all four
wave parameters together, and none at all under full management:

    score     share of the alongshore variation explained, 1 - SSE/SST, with
              the MODEL SMOOTHED LIKE THE COASTSAT TARGET (LOWESS over 10
              domains, the southern 10 raw: common.smooth_like_target),
              interior GIS 2-89, against the window's CoastSat LRR target.
              Bias, RMSE, correlation and the raw (unsmoothed) score beside it.
    coarse    Hs 0.75 1 1.5 2  x  Tp 7 8 10  x  asym 0.5 0.7 0.9
              x  high-angle 0.3 0.45 0.55  = 108 per period x scenario
              (Hs 0.65 and Tp 12 left out: they drowned the barrier)
    scope     natural and full management, 1996-2010 first, then 2010-2024
    refine    per period x scenario, a 3x3x3x3 grid at half the coarse step
              around the best coarse-or-refine setting (the midpoints to the
              neighbouring coarse values), launched automatically
    cross     the top 5 per period x scenario that lack a run in the other
              window are run there, so the shared pick has candidates
    shared    one setting for both windows: the lowest mean of smoothed RMSE
              / that window's flat-line RMSE, among settings run in both

Every run here is made on the Barrier3D route_overwash fix (checked at
start). The step-2 runs are not reused as grid cells (most predate the fix);
they are rescored on the smoothed output into tables/step2_rescored_smoothed.csv
for comparison.

WHERE: output/raw_runs/experiments/wave-climate/2026-09-25-wave-grid-smoothed-score/
    README.md, tables/, figures/, logs/<phase>_<scenario>/<period>/<settings>.log
    runs/<phase>_<scenario>/<period>/zeroBE/<run_name>/   (on disk only)
    phase = coarse, refine, cross

USAGE (from the project root):
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py run all
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py run coarse --periods 1996
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py score
    python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py rescore-step2
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
import HAT_metres_1_offset_units as common  # noqa: E402
import HAT_metres_2_wave_sensitivity as step2  # noqa: E402

HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RAW_RUNS = PROJECT_ROOT / "output" / "raw_runs"
TAG = "wave-climate/2026-09-25-wave-grid-smoothed-score"
STUDY_DIR = RAW_RUNS / "experiments" / TAG
TABLES_DIR = STUDY_DIR / "tables"
LOGS_DIR = STUDY_DIR / "logs"
RUN_TIMEOUT_S = 3600

PERIODS = (1996, 2010)
SCENARIOS = ("natural", "full_management")
PRESET = "zeroBE"
KEYS = ("hs", "wave_period_s", "wave_asymmetry", "wave_angle_high_fraction")
COARSE = {"hs": (0.75, 1.0, 1.5, 2.0),
          "wave_period_s": (7.0, 8.0, 10.0),
          "wave_asymmetry": (0.5, 0.7, 0.9),
          "wave_angle_high_fraction": (0.3, 0.45, 0.55)}
CROSS_TOP = 5


def window(period):
    return f"{period}_{period + 14}"


def label(s):
    return step2.settings_label(s)


def group(phase, scenario):
    return f"{phase}_{scenario}"


def log_path(phase, scenario, period, s):
    return LOGS_DIR / group(phase, scenario) / window(period) / f"{label(s)}.log"


def run_env(phase, scenario, period, s):
    import os
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1",
        "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": PRESET,
        "HAT_SCENARIO": scenario,
        "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "false",
        "HAT_OFFSET_MODE": "metres",
        "HAT_ISLAND_OFFSET_SOURCE": "duneline",
        "HAT_RUN_KIND": "experiment",
        "HAT_RUN_TAG": f"{TAG}/runs/{group(phase, scenario)}",
        "HAT_OVERWRITE": "false",
        "HAT_SAVE_MODEL_STATE": "false",
        "HAT_MAKE_GIFS": "false",
        "MPLBACKEND": "Agg",
        "PYTHONIOENCODING": "utf-8",
    })
    env.update({step2.ENV[k]: f"{v}" for k, v in s.items()})
    return env


def finished(log):
    """A clean finish or a drowned barrier is a result; anything else re-runs."""
    if not log.is_file():
        return False
    text = log.read_text(encoding="utf-8", errors="replace")
    index_race = "rebuild_run_index" in text and "PermissionError" in text
    return "Traceback" not in text or "Model stopped at year" in text or index_race


def launch(cell):
    phase, scenario, period, s = cell
    name = f"{group(phase, scenario)} {window(period)} {label(s)}"
    log = log_path(phase, scenario, period, s)
    if finished(log):
        return True
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    proc = subprocess.run([sys.executable, str(HINDCAST)],
                          env=run_env(phase, scenario, period, s),
                          cwd=str(PROJECT_ROOT), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=RUN_TIMEOUT_S)
    log.write_text((proc.stdout or "") + "\n--- STDERR ---\n" + (proc.stderr or ""),
                   encoding="utf-8")
    minutes = (time.perf_counter() - t0) / 60
    if proc.returncode != 0:
        print(f"FAILED {name} after {minutes:.1f} min: {step2.stop_reason(log)}", flush=True)
        return False
    print(f"done {name} in {minutes:.1f} min", flush=True)
    return True


def run_cells(cells, jobs):
    todo = [c for c in cells if not finished(log_path(*c))]
    print(f"{len(cells)} cells, {len(cells) - len(todo)} already run, {len(todo)} to run, "
          f"{jobs} at a time", flush=True)
    if todo:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(launch, todo))


def check_barrier3d():
    from cascade_pipeline.run_registry import barrier3d_provenance
    b = barrier3d_provenance()
    if not b.get("route_overwash_fix"):
        raise SystemExit(f"refusing to run: Barrier3D without the route_overwash fix ({b})")
    print(f"BARRIER3D = {b['branch']} {b['commit'][:7]} (route_overwash fix)", flush=True)


# =============================================================================
# CELLS
# =============================================================================

def coarse_cells(periods, scenarios):
    return [("coarse", sc, p, dict(zip(KEYS, v)))
            for p in periods for sc in scenarios for v in product(*(COARSE[k] for k in KEYS))]


def refine_values(key, best):
    """The best value and the midpoints to its coarse neighbours (half the
    coarse step), clipped at the ends of the coarse range."""
    c = list(COARSE[key])
    if best in c:
        i = c.index(best)
        lo = c[i - 1] if i > 0 else None
        hi = c[i + 1] if i < len(c) - 1 else None
    else:                                   # a refine value: its own bracket
        lo = max((v for v in c if v < best), default=None)
        hi = min((v for v in c if v > best), default=None)
    out = {round(best, 4)}
    if lo is not None:
        out.add(round((lo + best) / 2, 4))
    if hi is not None:
        out.add(round((best + hi) / 2, 4))
    return sorted(out)


def refine_cells(t):
    cells = []
    for p, sc in product(PERIODS, SCENARIOS):
        x = t[(t.period_start == p) & (t.scenario == sc) & (t.status == "scored")]
        if x.empty:
            continue
        b = x.loc[x.smoothed_variance_explained.idxmax()]
        grid = [refine_values(k, float(b[k])) for k in KEYS]
        print(f"refine {sc} {window(p)} around {label({k: b[k] for k in KEYS})} "
              f"({100 * b.smoothed_variance_explained:+.1f}%): "
              + " x ".join(str(g) for g in grid), flush=True)
        cells += [("refine", sc, p, dict(zip(KEYS, v))) for v in product(*grid)]
    return drop_already_run(cells, t)


def cross_cells(t):
    """The top CROSS_TOP settings per period x scenario, run in the other window."""
    cells = []
    for p, sc in product(PERIODS, SCENARIOS):
        other = PERIODS[1 - PERIODS.index(p)]
        x = t[(t.period_start == p) & (t.scenario == sc) & (t.status == "scored")]
        for _, r in x.nlargest(CROSS_TOP, "smoothed_variance_explained").iterrows():
            cells.append(("cross", sc, other, {k: float(r[k]) for k in KEYS}))
    return drop_already_run(cells, t)


def drop_already_run(cells, t):
    """Skip a cell whose settings already have a result (any phase) in that
    period and scenario: a refine grid overlaps the coarse one at its centre."""
    have = {(int(r.period_start), r.scenario, label({k: r[k] for k in KEYS}))
            for _, r in t.iterrows() if r.status != "not run"}
    out, seen = [], set()
    for c in cells:
        key = (c[2], c[1], label(c[3]))
        if key in have or key in seen:
            continue
        seen.add(key)
        out.append(c)
    return out


# =============================================================================
# SCORE
# =============================================================================

def score_run(run_dir, target):
    rates = common.run_rates(run_dir)
    raw = common.alongshore_scores(rates, target)
    sm = common.smooth_like_target(rates)
    s = common.alongshore_scores(sm, target)
    mi, oi = common.interior(rates), common.interior(target)
    return {"smoothed_variance_explained": s["variance_explained"],
            "smoothed_pattern_variance_explained": s["pattern_variance_explained"],
            "smoothed_rmse_m_yr": s["_rmse"],
            "smoothed_r": s["r_alongshore"], "smoothed_sd_ratio": s["sd_ratio"],
            "bias_m_yr": float((mi - oi).mean()),
            "raw_variance_explained": raw["variance_explained"],
            "raw_rmse_m_yr": raw["_rmse"], "raw_r": raw["r_alongshore"]}


def targets():
    t = {p: common.coastsat_target(p) for p in PERIODS}
    flat = {p: float(common.interior(t[p]).std(ddof=0)) for p in PERIODS}
    return t, flat


def cmd_score(_=None, quiet=False):
    import pandas as pd
    from cascade_pipeline.run_registry import load_run_index, rebuild_run_index
    rebuild_run_index(RAW_RUNS)
    index = load_run_index(RAW_RUNS / "run_index.csv")
    index = index[index["tag"].astype(str).str.startswith(TAG + "/")
                  & (index["status"] == "current")]
    tg, flat = targets()
    by_key = {}
    for _, r in index.iterrows():
        run_dir = (RAW_RUNS / "experiments" / r["tag"] / f"{r['start_year']}_{r['end_year']}"
                   / r["source_sink_preset"] / r["run_name"])
        md = json.loads((run_dir / f"{r['run_name']}_run_metadata.json").read_text(encoding="utf-8"))
        w = md["wave climate"]
        s = {"hs": float(w["wave_height_m"]), "wave_period_s": float(w["wave_period_s"]),
             "wave_asymmetry": float(w["wave_asymmetry"]),
             "wave_angle_high_fraction": float(w["wave_angle_high_frac"])}
        by_key[(r["tag"].split("/")[-1], int(r["start_year"]), label(s))] = (run_dir, md)
    rows = []
    for log in sorted(LOGS_DIR.glob("*/*/*.log")):
        g, win = log.parts[-3], log.parts[-2]
        phase, scenario = g.split("_", 1)
        period = int(win[:4])
        vals = [float(v) for v in log.stem.replace("Hs", "").replace("period", "").replace(
            "asymmetry", "").replace("highangle", "").split("_")]
        s = dict(zip(KEYS, vals))
        rec = {"phase": phase, "scenario": scenario, "period": win.replace("_", "-"),
               "period_start": period, **s}
        hit = by_key.get((g, period, label(s)))
        if hit and finished(log):
            run_dir, md = hit
            if md["scenario"]["shoreline offset"] != "metres":
                raise ValueError(f"{run_dir}: offset {md['scenario']['shoreline offset']!r}")
            rec.update(status="scored", **score_run(run_dir, tg[period]),
                       barrier3d_fix=md.get("identity", {}).get("barrier3d_route_overwash_fix"),
                       run_dir=str(run_dir.relative_to(STUDY_DIR)).replace("\\", "/"))
        else:
            rec.update(status=step2.stop_reason(log) if log.is_file() else "not run")
        rec["log"] = str(log.relative_to(STUDY_DIR)).replace("\\", "/")
        rows.append(rec)
    t = pd.DataFrame(rows)
    if t.empty:
        return t
    t["smoothed_rmse_over_flat"] = t.smoothed_rmse_m_yr / t.period_start.map(flat) \
        if "smoothed_rmse_m_yr" in t else np.nan
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    t.to_csv(TABLES_DIR / "all_runs.csv", index=False)
    pd.DataFrame([{"period": window(p).replace("_", "-"),
                   "observed_mean_m_yr": float(common.interior(tg[p]).mean()),
                   "flat_line_rmse_m_yr": flat[p]} for p in PERIODS]
                 ).to_csv(TABLES_DIR / "observed_targets.csv", index=False)
    best_tables(t)
    if not quiet:
        print(t.groupby(["phase", "scenario", "period", "status"]).size().to_string())
    return t


def best_tables(t):
    import pandas as pd
    s = t[t.status == "scored"]
    if s.empty:
        return
    per = [s[(s.period_start == p) & (s.scenario == sc)]
           .nlargest(1, "smoothed_variance_explained").assign(rule="per_period")
           for p, sc in product(PERIODS, SCENARIOS)]
    per = pd.concat([x for x in per if not x.empty])
    shared = []
    for sc in SCENARIOS:
        x = s[s.scenario == sc].drop_duplicates(["period_start", *KEYS])
        a = x[x.period_start == PERIODS[0]].set_index(list(KEYS))
        b = x[x.period_start == PERIODS[1]].set_index(list(KEYS))
        both = a[["smoothed_rmse_over_flat"]].join(b[["smoothed_rmse_over_flat"]],
                                                   lsuffix="_a", rsuffix="_b", how="inner")
        if both.empty:
            continue
        both["mean"] = both.mean(axis=1)
        key = both["mean"].idxmin()
        for p, frame in ((PERIODS[0], a), (PERIODS[1], b)):
            r = frame.loc[[key]].reset_index().assign(
                rule="shared", shared_mean_rmse_over_flat=float(both["mean"].min()),
                shared_candidates=len(both))
            shared.append(r)
    out = pd.concat([per, *shared], ignore_index=True)
    cols = ["rule", "scenario", "period", *KEYS, "smoothed_variance_explained",
            "smoothed_rmse_m_yr", "bias_m_yr", "smoothed_r", "raw_variance_explained",
            "smoothed_rmse_over_flat", "shared_mean_rmse_over_flat", "shared_candidates",
            "phase", "run_dir"]
    out[[c for c in cols if c in out]].to_csv(TABLES_DIR / "best_settings.csv", index=False)


def cmd_rescore_step2(_=None):
    """The 2026-09-24 step-2 runs scored the same way, for comparison."""
    import pandas as pd
    t = pd.read_csv(step2.TABLES_DIR / "all_runs.csv")
    tg, flat = targets()
    rows = []
    for _, r in t[t.status == "scored"].iterrows():
        rows.append({**{k: r[k] for k in ("group", "scenario", "period", "period_start", *KEYS)},
                     **score_run(step2.STUDY_DIR / r.run_dir, tg[int(r.period_start)]),
                     "step2_run_dir": r.run_dir})
    out = pd.DataFrame(rows)
    out["smoothed_rmse_over_flat"] = out.smoothed_rmse_m_yr / out.period_start.map(flat)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES_DIR / "step2_rescored_smoothed.csv", index=False)
    print(f"{len(out)} step-2 runs rescored -> tables/step2_rescored_smoothed.csv")
    for (p, sc), x in out.groupby(["period", "scenario"]):
        b = x.loc[x.smoothed_variance_explained.idxmax()]
        print(f"  {p} {sc:16s} best smoothed {100 * b.smoothed_variance_explained:+6.1f}% "
              f"(raw {100 * b.raw_variance_explained:+6.1f}%) at {label({k: b[k] for k in KEYS})}")
    return 0


# =============================================================================
# MAIN
# =============================================================================

def cmd_run(a):
    check_barrier3d()
    common.keep_awake()
    periods = a.periods
    if a.phase in ("coarse", "all"):
        for p in periods:                        # 1996-2010 first (Hannah)
            print(f"=== coarse {window(p)}", flush=True)
            run_cells(coarse_cells([p], SCENARIOS), a.jobs)
            cmd_score(quiet=True)
    if a.phase in ("refine", "all"):
        print("=== refine", flush=True)
        run_cells(refine_cells(cmd_score(quiet=True)), a.jobs)
    if a.phase in ("cross", "all"):
        print("=== cross", flush=True)
        run_cells(cross_cells(cmd_score(quiet=True)), a.jobs)
    cmd_score()
    return 0


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("phase", choices=("all", "coarse", "refine", "cross"))
    r.add_argument("--periods", nargs="+", type=int, choices=PERIODS, default=list(PERIODS))
    r.add_argument("--jobs", type=int, default=6)
    sub.add_parser("score")
    sub.add_parser("rescore-step2")
    a = p.parse_args()
    if a.cmd == "score":
        cmd_score()
        return 0
    if a.cmd == "rescore-step2":
        return cmd_rescore_step2()
    return cmd_run(a)


if __name__ == "__main__":
    sys.exit(main())
