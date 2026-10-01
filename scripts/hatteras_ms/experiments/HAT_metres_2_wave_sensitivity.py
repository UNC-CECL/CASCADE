#!/usr/bin/env python3
"""
Wave-climate sensitivity with the offset in metres, natural scenario, both windows.

    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py run stage1 --jobs 6
    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py run stage1 --scenario full_management
    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py score
    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py score-window --start 2010 --end 2020
    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py run stage2 --jobs 6
    python scripts/hatteras_ms/experiments/HAT_metres_2_wave_sensitivity.py run combo --hs 1 --tp 8 --asym 0.7 --ahf 0.4

Stage 1 moves one parameter at a time around the baseline; stage 2 is a
5 x 5 grid over the two that matter most. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "hatteras_ms" / "experiments"))
# Shared with the offset-scale study
import HAT_metres_1_offset_units as common  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RAW_RUNS = PROJECT_ROOT / "output" / "raw_runs"
STUDY_TAG = "wave-climate/2026-09-24-metres-2-wave-sensitivity"
STUDY_DIR = RAW_RUNS / "experiments" / STUDY_TAG
TABLES_DIR = STUDY_DIR / "tables"
LOGS_DIR = STUDY_DIR / "logs"
RUN_TIMEOUT_S = 3600

PERIODS = (1996, 2010)
PRESET = "zeroBE"
SCENARIO = "natural"
MANAGED = "full_management"

# setting -> environment variable the runner reads (HAT_ + field name)
ENV = {"hs": "HAT_HS", "wave_period_s": "HAT_WAVE_PERIOD_S",
       "wave_asymmetry": "HAT_WAVE_ASYMMETRY",
       "wave_angle_high_fraction": "HAT_WAVE_ANGLE_HIGH_FRACTION"}
# The baseline stays at asymmetry 0.8 until re-centring on 0.7 is run (README)
BASELINE = {"hs": 1.0, "wave_period_s": 8.0, "wave_asymmetry": 0.8,
            "wave_angle_high_fraction": 0.45}
# Stage 2 (the two grids) was chosen and run around the first baseline, asymmetry 0.8, and stays there
STAGE2_BASELINE = {**BASELINE, "wave_asymmetry": 0.8}
# folder name -> (setting, stage-1 values, label)
PARAMS = {
    "wave_height": ("hs", (0.65, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0),
                    "Significant wave height, Hs (m)"),
    "high_angle": ("wave_angle_high_fraction", (0.1, 0.2, 0.3, 0.4, 0.45, 0.5, 0.55),
                   "Fraction of high-angle waves (> 45°)"),
    "asymmetry": ("wave_asymmetry", (0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
                  "Wave asymmetry (share of waves from the left)"),
    "wave_period": ("wave_period_s", (6.0, 7.0, 8.0, 10.0, 12.0),
                    "Wave period, Tp (s)"),
}
SETTING_TO_PARAM = {v[0]: k for k, v in PARAMS.items()}
GRID_SIZE = 5
# -----------------------------------------------------------------------------


# A start year's window label
def window(period):
    return f"{period}_{period + 14}"


# A setting's label, every parameter spelled
def settings_label(s):
    return (f"Hs{common.num(s['hs'])}_period{common.num(s['wave_period_s'])}"
            f"_asymmetry{common.num(s['wave_asymmetry'])}"
            f"_highangle{common.num(s['wave_angle_high_fraction'])}")


# A cell's log file
def log_path(group, period, s, scenario=SCENARIO):
    return LOGS_DIR / group / window(period) / f"{settings_label(s)}.log"


# Cells

# Why a run has no score
def stop_reason(log):
    import re
    text = log.read_text(encoding="utf-8", errors="replace")
    if "Model stopped at year" in text or "Traceback" in text:
        return common.stop_reason(log)
    years = re.findall(r"(\d+)/14 \[", text)
    return (f"process crashed in year {years[-1] if years else '?'} with no Python error "
            f"(likely Barrier3D route_overwash access violation)")


# Full-management sweep: the stage-1 values, filed as full_management_<parameter>/
MANAGED_PREFIX = "full_management_"


# The scenario a group runs
def scenario_of(group):
    return (MANAGED if group == "baseline_full_management"
            or group.startswith(MANAGED_PREFIX) else SCENARIO)


# Stage 1 under full management
def stage1_managed_cells():
    cells = []
    for period in PERIODS:
        cells.append(("baseline_full_management", period, dict(BASELINE), MANAGED))
        for group, (setting, values, _) in PARAMS.items():
            for v in values:
                if np.isclose(v, BASELINE[setting]):
                    continue                     # that cell IS the baseline
                cells.append((MANAGED_PREFIX + group, period,
                              {**BASELINE, setting: float(v)}, MANAGED))
    return cells


# Stage 1: the baseline and one parameter at a time
def stage1_cells():
    cells = []
    for period in PERIODS:
        cells.append(("baseline", period, dict(BASELINE), SCENARIO))
        cells.append(("baseline_full_management", period, dict(BASELINE), MANAGED))
        for group, (setting, values, _) in PARAMS.items():
            for v in values:
                if np.isclose(v, BASELINE[setting]):
                    continue                     # that cell IS the baseline
                cells.append((group, period, {**BASELINE, setting: float(v)}, SCENARIO))
    return cells


# The environment one run reads
def run_env(group, period, s, scenario):
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
        "HAT_RUN_TAG": f"{STUDY_TAG}/runs/{group}",
        "HAT_OVERWRITE": "false",
        "HAT_SAVE_MODEL_STATE": "false",
        "HAT_MAKE_GIFS": "false",
        "MPLBACKEND": "Agg",
        "PYTHONIOENCODING": "utf-8",
    })
    env.update({ENV[k]: f"{v}" for k, v in s.items()})
    return env


# One run, skipped if already a result
def launch(cell, dry_run=False):
    group, period, s, scenario = cell
    label = f"{group} {window(period)} {settings_label(s)}"
    log = log_path(group, period, s)
    if log.is_file():
        text = log.read_text(encoding="utf-8", errors="replace")
        # A clean finish, a drowned barrier or an index-lock failure is a result: skip it
        index_race = "rebuild_run_index" in text and "PermissionError" in text
        if "Traceback" not in text or "Model stopped at year" in text or index_race:
            print(f"skip {label}: already run", flush=True)
            return True
    if dry_run:
        print(f"would run: {label}")
        return True
    log.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    proc = subprocess.run([sys.executable, str(HINDCAST)],
                          env=run_env(group, period, s, scenario),
                          cwd=str(PROJECT_ROOT), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", timeout=RUN_TIMEOUT_S)
    log.write_text((proc.stdout or "") + "\n--- STDERR ---\n" + (proc.stderr or ""),
                   encoding="utf-8")
    minutes = (time.perf_counter() - t0) / 60
    if proc.returncode != 0:
        print(f"FAILED {label} after {minutes:.1f} min: {stop_reason(log)}", flush=True)
        return False
    print(f"done {label} in {minutes:.1f} min", flush=True)
    return True


# Run cells in parallel, keeping the machine awake
def run_cells(cells, jobs, dry_run):
    common.keep_awake()
    print(f"{len(cells)} runs, {jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        ok = list(pool.map(lambda c: launch(c, dry_run), cells))
    print(f"{sum(ok)} of {len(cells)} succeeded (a drowned barrier counts as not succeeded)")
    return 0


# Stage 2

# Selection rule: the range over runs that survived in 1996-2010 (README)
SELECTION_PERIOD = 1996


# The two parameters with the largest effect, and 5 values for each
def choose_grid():
    import pandas as pd
    t = pd.read_csv(TABLES_DIR / "all_runs.csv")
    t = t[(t.scenario == SCENARIO) & t.group.isin(["baseline", *PARAMS])
          & (t.period_start == SELECTION_PERIOD) & (t.status == "scored")]
    rows, best = [], {}
    for group, (setting, values, _) in PARAMS.items():
        others = [k for k in STAGE2_BASELINE if k != setting]
        mask = np.logical_and.reduce([np.isclose(t[k], STAGE2_BASELINE[k]) for k in others])
        line = t[mask].drop_duplicates(setting).set_index(setting)["variance_explained"]
        best[group] = float(line.idxmax())
        rows.append(dict(parameter=group, setting=setting,
                         n_survived=len(line), n_values=len(values),
                         effect=float(line.max() - line.min()),
                         best_value=best[group], best_variance_explained=float(line.max()),
                         rule=f"range of variance explained over runs that survived in "
                              f"{SELECTION_PERIOD}-{SELECTION_PERIOD + 14}"))
    sel = pd.DataFrame(rows).sort_values("effect", ascending=False).reset_index(drop=True)
    pair = list(sel.parameter[:2])
    grid = {}
    for group in pair:
        values = list(PARAMS[group][1])
        i = values.index(best[group])
        lo = min(max(0, i - GRID_SIZE // 2), len(values) - GRID_SIZE)
        grid[group] = values[lo:lo + GRID_SIZE]
    sel["chosen_for_grid"] = sel.parameter.isin(pair)
    sel["grid_values"] = [", ".join(f"{v:g}" for v in grid[p]) if p in grid else ""
                          for p in sel.parameter]
    sel.to_csv(TABLES_DIR / "stage2_selection.csv", index=False)
    print(sel.to_string(index=False))
    return pair, grid


# A named grid
def grid_cells(pair, values, periods):
    import pandas as pd
    group = f"grid_{pair[0]}_x_{pair[1]}"
    s1, s2 = PARAMS[pair[0]][0], PARAMS[pair[1]][0]
    t = pd.read_csv(TABLES_DIR / "all_runs.csv")
    t = t[t.scenario == SCENARIO]
    cells = []
    for period in periods:
        done = {settings_label({k: r[k] for k in STAGE2_BASELINE}) for _, r in
                t[t.period_start == period].iterrows()}
        for v1 in values[0]:
            for v2 in values[1]:
                st = {**STAGE2_BASELINE, s1: float(v1), s2: float(v2)}
                if settings_label(st) not in done:
                    cells.append((group, period, st, SCENARIO))
    print(f"grid {group}: {len(values[0])} x {len(values[1])} in {list(periods)}, "
          f"{len(cells)} new runs (the rest already run)")
    return cells


# The stage-2 grid, minus cells stage 1 already ran
def stage2_cells():
    pair, grid = choose_grid()
    group = f"grid_{pair[0]}_x_{pair[1]}"
    s1, s2 = PARAMS[pair[0]][0], PARAMS[pair[1]][0]
    done = {(c[1], settings_label(c[2])) for c in stage1_cells() if c[3] == SCENARIO}
    cells = []
    for period in PERIODS:
        for v1 in grid[pair[0]]:
            for v2 in grid[pair[1]]:
                s = {**STAGE2_BASELINE, s1: float(v1), s2: float(v2)}
                if (period, settings_label(s)) in done:
                    continue                     # already a stage-1 run
                cells.append((group, period, s, SCENARIO))
    print(f"grid {group}: {len(grid[pair[0]])} x {len(grid[pair[1]])} per period, "
          f"{len(cells)} new runs (the rest are stage-1 runs)")
    return cells


# Score

# Score every logged run and write the tables
def cmd_score(a):
    import pandas as pd
    from cascade_pipeline.run_registry import load_run_index, rebuild_run_index

    rebuild_run_index(RAW_RUNS)
    index = load_run_index(RAW_RUNS / "run_index.csv")
    index = index[index["tag"].astype(str).str.startswith(STUDY_TAG + "/")
                  & (index["status"] == "current")]
    targets = {p: common.coastsat_target(p) for p in PERIODS}
    run_targets = {p: common.coastsat_target(p, lowess_domains=common.RUN_TARGET_DOMAINS) for p in PERIODS}
    runs = []
    for _, r in index.iterrows():
        run_dir = (RAW_RUNS / "experiments" / r["tag"]
                   / f"{r['start_year']}_{r['end_year']}" / r["source_sink_preset"]
                   / r["run_name"])
        md = json.loads((run_dir / f"{r['run_name']}_run_metadata.json").read_text(encoding="utf-8"))
        w = md["wave climate"]
        runs.append((r, run_dir, md, {"hs": float(w["wave_height_m"]),
                                      "wave_period_s": float(w["wave_period_s"]),
                                      "wave_asymmetry": float(w["wave_asymmetry"]),
                                      "wave_angle_high_fraction": float(w["wave_angle_high_frac"])}))
    records = []
    for log in sorted(LOGS_DIR.glob("*/*/*.log")):
        group, win = log.parts[-3], log.parts[-2]
        # crash_check logs are diagnostics, not sweep cells
        if group.startswith("crash_check") or log.stem.endswith("_boundscheck"):
            continue
        period = int(win[:4])
        vals = log.stem.replace("Hs", "").replace("period", "").replace(
            "asymmetry", "").replace("highangle", "").split("_")
        s = dict(zip(("hs", "wave_period_s", "wave_asymmetry", "wave_angle_high_fraction"),
                     map(float, vals)))
        scenario = scenario_of(group)
        rec = {"group": group, "period": win.replace("_", "-"), "period_start": period,
               "scenario": scenario, **s, "source_sink_preset": PRESET, "groin": False,
               "offset_mode": "metres",
               "log": str(log.relative_to(STUDY_DIR)).replace("\\", "/")}
        match = [x for x in runs if x[0]["tag"] == f"{STUDY_TAG}/runs/{group}"
                 and int(x[0]["start_year"]) == period and x[3] == s]
        if match:
            r, run_dir, md, _ = match[0]
            if md["scenario"]["shoreline offset"] != "metres":
                raise ValueError(f"{run_dir}: offset mode {md['scenario']['shoreline offset']!r}")
            # Rebuilt at the runs' own window, all four numbers must reproduce the runner's
            for k, v in common.rate_skill(run_dir, run_targets[period]).items():
                if not np.isclose(v, float(r[k]), rtol=1e-3, atol=1e-6):
                    raise ValueError(f"{run_dir}: {k} against the rebuilt target does not "
                                     f"match the runner's")
            # Then scored against the current target (common.SMOOTH_DOMAINS)
            skill = common.rate_skill(run_dir, targets[period])
            sc = common.alongshore_scores(common.run_rates(run_dir), targets[period])
            assert np.isclose(sc.pop("_rmse"), skill["rmse_interior_m_yr"])
            # Which Barrier3D: a run without this field ran on the unfixed model
            fix = md["identity"].get("barrier3d_route_overwash_fix")
            rec["barrier3d_route_overwash_fix"] = bool(fix[0] if isinstance(fix, list) else fix)                 if fix is not None else False
            rec.update(status="scored", island_offset_version=md["identity"]["island_offset_version"],
                       mean_bias_interior_m_yr=skill["mean_bias_interior_m_yr"],
                       rmse_interior_m_yr=skill["rmse_interior_m_yr"], **sc,
                       roads_drowned=r["roads_drowned"], run_name=r["run_name"],
                       run_dir=str(run_dir.relative_to(STUDY_DIR)).replace("\\", "/"))
        else:
            rec["status"] = stop_reason(log)
        records.append(rec)
    cols = ["group", "period", "period_start", "scenario", "hs", "wave_period_s",
            "wave_asymmetry", "wave_angle_high_fraction", "source_sink_preset", "groin",
            "offset_mode", "island_offset_version", "barrier3d_route_overwash_fix", "status",
            "mean_bias_interior_m_yr",
            "rmse_interior_m_yr", "variance_explained", "pattern_variance_explained",
            "r_alongshore", "model_sd_m_yr", "sd_ratio", "roads_drowned", "run_name",
            "run_dir", "log"]
    out = pd.DataFrame(records).reindex(columns=cols).sort_values(
        ["period_start", "group", "hs", "wave_period_s", "wave_asymmetry",
         "wave_angle_high_fraction"])
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES_DIR / "all_runs.csv", index=False)
    obs = []
    for p, t in targets.items():
        i = common.interior(t)
        obs.append(dict(period=window(p).replace("_", "-"),
                        target=f"CoastSat LRR, LOWESS {common.SMOOTH_DOMAINS} domains", domains="GIS 2-89",
                        mean_m_yr=i.mean(), sd_m_yr=i.std(ddof=0),
                        flat_line_rmse_m_yr=i.std(ddof=0)))
    pd.DataFrame(obs).to_csv(TABLES_DIR / "observed_targets.csv", index=False)
    print(out.groupby(["period", "group", "status"]).size().to_string())
    return 0


# Score a window inside a period

# Score the 2010-2024 runs on 2010-2020, clear of the 2021 CoastSat step
def cmd_score_window(a):
    import pandas as pd
    from cascade_pipeline.shoreline import compute_lrr
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as G

    start, end = a.start, a.end
    n_states = end - start + 1
    target = common.coastsat_target(start, end)
    full_target = common.coastsat_target(start)
    t = pd.read_csv(TABLES_DIR / "all_runs.csv")
    t = t[(t.period_start == start) & (t.status == "scored")]
    real = slice(G.start_real_index, G.end_real_index)
    gis = np.arange(G.first_gis_id, G.first_gis_id + G.num_real_domains)
    rows = []
    for _, r in t.iterrows():
        run_dir = STUDY_DIR / r.run_dir
        m = np.load(next(run_dir.glob("*_shoreline_matrix.npy")))
        table = common.run_rates(run_dir)
        # the order check: the full record must give the run's own table
        full, _ = compute_lrr(m)
        full = pd.Series(full[real], index=gis)
        if not np.allclose(full.values, table.reindex(gis).values, atol=1e-9):
            raise ValueError(f"{run_dir}: compute_lrr on the full record does not "
                             f"reproduce the run's own rate table; domain order unknown")
        sub_lrr, _ = compute_lrr(m[:n_states], span_years=end - start)
        rates = pd.Series(sub_lrr[real], index=gis)
        sc = common.alongshore_scores(rates, target)
        mi, oi = common.interior(rates), common.interior(target)
        rows.append({**{k: r[k] for k in ("group", "scenario", "hs", "wave_period_s",
                                           "wave_asymmetry", "wave_angle_high_fraction")},
                     "window": f"{start}-{end}",
                     "bias_m_yr": float((mi - oi).mean()),
                     "rmse_m_yr": sc.pop("_rmse"),
                     "variance_explained": sc["variance_explained"],
                     "pattern_variance_explained": sc["pattern_variance_explained"],
                     "r_alongshore": sc["r_alongshore"], "sd_ratio": sc["sd_ratio"],
                     "full_window_bias_m_yr": r.mean_bias_interior_m_yr,
                     "full_window_rmse_m_yr": r.rmse_interior_m_yr,
                     "full_window_variance_explained": r.variance_explained,
                     "model_mean_m_yr": float(mi.mean()), "observed_mean_m_yr": float(oi.mean()),
                     "run_dir": r.run_dir})
    out = pd.DataFrame(rows).sort_values(["scenario", "group", "hs", "wave_period_s",
                                          "wave_asymmetry", "wave_angle_high_fraction"])
    path = TABLES_DIR / f"window_{start}_{end}.csv"
    out.to_csv(path, index=False)
    i = common.interior(target)
    print(f"{len(out)} runs scored on {start}-{end} -> {path.name}")
    print(f"observed {start}-{end}: mean {i.mean():+.3f} m/yr, flat-line RMSE {i.std(ddof=0):.3f}"
          f"  (full period: mean {common.interior(full_target).mean():+.3f}, "
          f"flat-line RMSE {common.interior(full_target).std(ddof=0):.3f})")
    return 0


# Run: the subcommand asked for
def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("stage", choices=("stage1", "stage2", "grid", "combo"))
    r.add_argument("--pair", nargs=2, choices=tuple(PARAMS),
                   help="grid only: the two parameters")
    r.add_argument("--values1", nargs="+", type=float, help="grid only")
    r.add_argument("--values2", nargs="+", type=float, help="grid only")
    r.add_argument("--periods", nargs="+", type=int, choices=PERIODS, default=list(PERIODS),
                   help="grid only")
    r.add_argument("--scenario", choices=(SCENARIO, MANAGED), default=SCENARIO,
                   help="stage1 only: natural (default) or the full_management sweep")
    r.add_argument("--hs", type=float, help="combo only (default: the baseline's)")
    r.add_argument("--tp", type=float, help="combo only")
    r.add_argument("--asym", type=float, help="combo only")
    r.add_argument("--ahf", type=float, help="combo only")
    r.add_argument("--scenarios", nargs="+", choices=(SCENARIO, MANAGED),
                   default=[SCENARIO, MANAGED], help="combo only")
    r.add_argument("--jobs", type=int, default=6)
    r.add_argument("--dry-run", action="store_true")
    s = sub.add_parser("score")
    w = sub.add_parser("score-window", help="score one period's runs on a window inside it")
    w.add_argument("--start", type=int, default=2010)
    w.add_argument("--end", type=int, default=2020)
    a = p.parse_args()
    if a.cmd == "score":
        return cmd_score(a)
    if a.cmd == "score-window":
        return cmd_score_window(a)
    if a.stage == "stage2":
        cells = stage2_cells()
    elif a.stage == "combo":
        # One hand-picked setting of all four, filed under runs/combos/
        s = {**BASELINE, **{k: v for k, v in (("hs", a.hs), ("wave_period_s", a.tp),
                                              ("wave_asymmetry", a.asym),
                                              ("wave_angle_high_fraction", a.ahf))
                            if v is not None}}
        cells = [(MANAGED_PREFIX + "combos" if sc == MANAGED else "combos", period, s, sc)
                 for sc in a.scenarios for period in a.periods]
    elif a.stage == "grid":
        cells = grid_cells(a.pair, (a.values1, a.values2), a.periods)
    else:
        cells = stage1_cells() if a.scenario == SCENARIO else stage1_managed_cells()
    return run_cells(cells, a.jobs, a.dry_run)


if __name__ == "__main__":
    sys.exit(main())
