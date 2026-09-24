#!/usr/bin/env python3
r"""
natural_wave_sensitivity.py -- wave-climate sensitivity, natural scenario, offset in metres
==============================================================================
THE QUESTION (Hannah, 2026-09-24, by interview). With the island offset in
metres, how does the model's alongshore shoreline change respond to each wave
parameter when nothing human acts on the island, in both canonical windows?

THE DESIGN (every choice Hannah's)
    scenario     natural: no road management, no beach/dune management, no
                 fills, no relocations; no groin; zeroBE; offset in metres
                 (dune line, CURRENT v1)
    periods      1996-2010 and 2010-2024, each against its own CoastSat LRR
    baseline     Hs 1.0 m, Tp 8 s, asymmetry 0.8, high-angle fraction 0.45
                 (the best metres point of experiments/2026-09-24-island-
                 offset-scale-wave-tuning, found under full management)
    stage 1      one parameter at a time around the baseline:
                   wave_height  Hs    0.65 0.75 1.0 1.25 1.5 2.0 2.5 3.0
                   high_angle   ahf   0.1 0.2 0.3 0.4 0.45 0.5 0.55
                   asymmetry    asym  0.3 0.4 0.5 0.6 0.7 0.8 0.9
                   wave_period  Tp    6 7 8 10 12
                 plus the baseline under full_management, per period
    stage 2      a 5 x 5 grid per period over the two parameters whose range
                 moves the share of alongshore variation explained the most
                 (mean over both periods; a drowned run counts as the worst
                 score of its period). Each axis: the 5 stage-1 values centred
                 on that parameter's best value, shifted inward at the ends.
                 The other two stay at the baseline. Cells already run are
                 reused, not re-run.

LAYOUT (output/raw_runs/sensitivity/2026-09-24-natural-waves/)
    README.md, tables/, figures/, logs/<group>/<period>/<settings>.log
    <group>/<period>/zeroBE/<run_name>/     group = baseline, wave_height,
                                            high_angle, asymmetry, wave_period,
                                            baseline_full_management,
                                            grid_<p1>_x_<p2>
    The run folders stay on disk only: this scenario's run name is ~85
    characters and sits twice in each path, past Windows' 260, so git cannot
    index them (the same choice as the offset-scale study).

USAGE
    python natural_wave_sensitivity.py run stage1 [--jobs 6] [--dry-run]
    python natural_wave_sensitivity.py run stage1 --scenario full_management
    python natural_wave_sensitivity.py score
    python natural_wave_sensitivity.py run stage2 [--jobs 6] [--dry-run]
==============================================================================
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
# Shared with the offset-scale study: the target, the alongshore scores, the
# drowning reader, the number spelling and the keep-awake.
import HAT_offset_scale_wave_tuning as common  # noqa: E402

HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
RAW_RUNS = PROJECT_ROOT / "output" / "raw_runs"
STUDY_TAG = "2026-09-24-natural-waves"
STUDY_DIR = RAW_RUNS / "sensitivity" / STUDY_TAG
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
BASELINE = {"hs": 1.0, "wave_period_s": 8.0, "wave_asymmetry": 0.8,
            "wave_angle_high_fraction": 0.45}
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


def window(period):
    return f"{period}_{period + 14}"


def settings_label(s):
    return (f"Hs{common.num(s['hs'])}_period{common.num(s['wave_period_s'])}"
            f"_asymmetry{common.num(s['wave_asymmetry'])}"
            f"_highangle{common.num(s['wave_angle_high_fraction'])}")


def log_path(group, period, s, scenario=SCENARIO):
    return LOGS_DIR / group / window(period) / f"{settings_label(s)}.log"


# =============================================================================
# CELLS
# =============================================================================

def stop_reason(log):
    """Why a run has no score: drowned, crashed, or the error it raised.

    A run that dies with no Python traceback and no drowning is the silent
    access violation in Barrier3D's jitted route_overwash (an out-of-bounds
    read once a domain has prograded; memory note of 2026-09-11, where it was
    found with the groin). Seen here with no groin, in 2010-2024 natural runs.
    """
    import re
    text = log.read_text(encoding="utf-8", errors="replace")
    if "Model stopped at year" in text or "Traceback" in text:
        return common.stop_reason(log)
    years = re.findall(r"(\d+)/14 \[", text)
    return (f"process crashed in year {years[-1] if years else '?'} with no Python error "
            f"(likely Barrier3D route_overwash access violation)")


# THE FULL-MANAGEMENT SWEEP (Hannah, 2026-09-24, after stage 1 showed
# management halves the 2010-2024 bias): the same stage-1 values under
# full_management, filed as full_management_<parameter>/ beside the natural
# folders. Its baseline is baseline_full_management/, already run.
MANAGED_PREFIX = "full_management_"


def scenario_of(group):
    return (MANAGED if group == "baseline_full_management"
            or group.startswith(MANAGED_PREFIX) else SCENARIO)


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
        "HAT_RUN_KIND": "sensitivity",
        "HAT_RUN_TAG": f"{STUDY_TAG}/{group}",
        "HAT_OVERWRITE": "false",
        "HAT_SAVE_MODEL_STATE": "false",
        "HAT_MAKE_GIFS": "false",
        "MPLBACKEND": "Agg",
        "PYTHONIOENCODING": "utf-8",
    })
    env.update({ENV[k]: f"{v}" for k, v in s.items()})
    return env


def launch(cell, dry_run=False):
    group, period, s, scenario = cell
    label = f"{group} {window(period)} {settings_label(s)}"
    log = log_path(group, period, s)
    if log.is_file():
        text = log.read_text(encoding="utf-8", errors="replace")
        # A clean finish, or a drowned barrier, is a result: re-running either
        # would only hit the runner's existing-folder guard and overwrite the
        # log with that refusal (found when the sweep was paused 2026-09-24).
        # ...and a run that finished and only failed to replace the shared
        # run_index.csv (a Windows lock between parallel runs, 2026-09-24):
        # its outputs are complete and scored.
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


def run_cells(cells, jobs, dry_run):
    common.keep_awake()
    print(f"{len(cells)} runs, {jobs} at a time", flush=True)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        ok = list(pool.map(lambda c: launch(c, dry_run), cells))
    print(f"{sum(ok)} of {len(cells)} succeeded (a drowned barrier counts as not succeeded)")
    return 0


# =============================================================================
# STAGE 2
# =============================================================================

# THE SELECTION RULE, and why it changed (2026-09-24). The first rule took
# each parameter's range of variation explained over BOTH periods and counted
# a run without a score as the worst score of its period. It chose Hs x Tp
# for the wrong reasons: 2010-2024 natural runs explain -800% to -2800% at
# every setting (a ~-4.5 m/yr bias), so its ranges measure how badly a
# setting fails; the drowned Tp 12 became each period's worst score; and two
# crashes were scored as worst although the rule named drownings only. That
# selection is kept in tables/stage2_selection_first_rule.csv and its 12
# finished runs in grid_wave_height_x_wave_period/. Hannah then chose: the
# range over runs that SURVIVED in 1996-2010 only, best value from the same
# runs; the grid still covers both periods.
SELECTION_PERIOD = 1996


def choose_grid():
    """The two parameters with the largest effect, and 5 values for each."""
    import pandas as pd
    t = pd.read_csv(TABLES_DIR / "all_runs.csv")
    t = t[(t.scenario == SCENARIO) & t.group.isin(["baseline", *PARAMS])
          & (t.period_start == SELECTION_PERIOD) & (t.status == "scored")]
    rows, best = [], {}
    for group, (setting, values, _) in PARAMS.items():
        others = [k for k in BASELINE if k != setting]
        mask = np.logical_and.reduce([np.isclose(t[k], BASELINE[k]) for k in others])
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


def stage2_cells():
    pair, grid = choose_grid()
    group = f"grid_{pair[0]}_x_{pair[1]}"
    s1, s2 = PARAMS[pair[0]][0], PARAMS[pair[1]][0]
    done = {(c[1], settings_label(c[2])) for c in stage1_cells() if c[3] == SCENARIO}
    cells = []
    for period in PERIODS:
        for v1 in grid[pair[0]]:
            for v2 in grid[pair[1]]:
                s = {**BASELINE, s1: float(v1), s2: float(v2)}
                if (period, settings_label(s)) in done:
                    continue                     # already a stage-1 run
                cells.append((group, period, s, SCENARIO))
    print(f"grid {group}: {len(grid[pair[0]])} x {len(grid[pair[1]])} per period, "
          f"{len(cells)} new runs (the rest are stage-1 runs)")
    return cells


# =============================================================================
# SCORE
# =============================================================================

def cmd_score(a):
    import pandas as pd
    from cascade_pipeline.run_registry import load_run_index, rebuild_run_index

    rebuild_run_index(RAW_RUNS)
    index = load_run_index(RAW_RUNS / "run_index.csv")
    index = index[index["tag"].astype(str).str.startswith(STUDY_TAG + "/")
                  & (index["status"] == "current")]
    targets = {p: common.coastsat_target(p) for p in PERIODS}
    runs = []
    for _, r in index.iterrows():
        run_dir = (RAW_RUNS / "sensitivity" / r["tag"]
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
        match = [x for x in runs if x[0]["tag"] == f"{STUDY_TAG}/{group}"
                 and int(x[0]["start_year"]) == period and x[3] == s]
        if match:
            r, run_dir, md, _ = match[0]
            if md["scenario"]["shoreline offset"] != "metres":
                raise ValueError(f"{run_dir}: offset mode {md['scenario']['shoreline offset']!r}")
            sc = common.alongshore_scores(common.run_rates(run_dir), targets[period])
            if not np.isclose(sc.pop("_rmse"), float(r["rmse_interior_m_yr"]), rtol=1e-3):
                raise ValueError(f"{run_dir}: RMSE against the rebuilt target does not "
                                 f"match the runner's")
            # which Barrier3D: recorded by the runner since 2026-09-24; a run
            # without the field predates it and ran on the unfixed model
            fix = md["identity"].get("barrier3d_route_overwash_fix")
            rec["barrier3d_route_overwash_fix"] = bool(fix[0] if isinstance(fix, list) else fix)                 if fix is not None else False
            rec.update(status="scored", island_offset_version=md["identity"]["island_offset_version"],
                       mean_bias_interior_m_yr=float(r["mean_bias_interior_m_yr"]),
                       rmse_interior_m_yr=float(r["rmse_interior_m_yr"]), **sc,
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
                        target="CoastSat LRR, LOESS 10 domains", domains="GIS 2-89",
                        mean_m_yr=i.mean(), sd_m_yr=i.std(ddof=0),
                        flat_line_rmse_m_yr=i.std(ddof=0)))
    pd.DataFrame(obs).to_csv(TABLES_DIR / "observed_targets.csv", index=False)
    print(out.groupby(["period", "group", "status"]).size().to_string())
    return 0


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    p = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("stage", choices=("stage1", "stage2"))
    r.add_argument("--scenario", choices=(SCENARIO, MANAGED), default=SCENARIO,
                   help="stage1 only: natural (default) or the full_management sweep")
    r.add_argument("--jobs", type=int, default=6)
    r.add_argument("--dry-run", action="store_true")
    s = sub.add_parser("score")
    a = p.parse_args()
    if a.cmd == "score":
        return cmd_score(a)
    if a.stage == "stage2":
        cells = stage2_cells()
    else:
        cells = stage1_cells() if a.scenario == SCENARIO else stage1_managed_cells()
    return run_cells(cells, a.jobs, a.dry_run)


if __name__ == "__main__":
    sys.exit(main())
