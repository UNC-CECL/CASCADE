"""
Why does a storm series with longer events drown the barrier?

    python scripts/hatteras_ms/experiments/HAT_storm_max_duration.py build
    python scripts/hatteras_ms/experiments/HAT_storm_max_duration.py run --workers 4
    python scripts/hatteras_ms/experiments/HAT_storm_max_duration.py cause
    python scripts/hatteras_ms/experiments/HAT_storm_max_duration.py diagnose
    python scripts/hatteras_ms/experiments/HAT_storm_max_duration.py compare

Rebuilds the storm series at longer caps with the builder's own functions and
runs them through the unchanged runner. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
"""
from __future__ import annotations

import argparse
import ast
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# --- CONFIG ------------------------------------------------------------------
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast_1984_2024.py"
BUILDER = PROJECT_ROOT / "scripts" / "input_prep" / "3-env-forcings" / "3-storms" / "historical_storm_creation_v3_HAT.py"
TAG = "storms-and-overwash/2026-09-28-storm-max-duration"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "matrix"

WINDOWS = ((1996, 2010), (2010, 2024))
VARIANTS = {"72": 72, "96": 96, "120": 120, "240": 240, "nocap": 10 ** 6, "72trim": 72}
RUN_VARIANTS = ("96", "120", "240", "72trim")      # nocap == 240 for 1996-2024 (longest event 193 h)
SCENARIOS = ("natural", "full_management")
CONTROLS = {
    (1996, "natural"): "HAT_1996_2010_edgeBE_offsetmetres_noroad_nobdm_nogroin",
    (1996, "full_management"): "HAT_1996_2010_edgeBE_offsetmetres_road_bdm_nogroin",
    (2010, "natural"): "HAT_2010_2024_edgeBE_offsetmetres_noroad_nobdm_nogroin",
    (2010, "full_management"): "HAT_2010_2024_edgeBE_offsetmetres_road_bdm_nourish_nogroin",
}
# the builder's settings (historical_storm_creation_v3_HAT.py, "inputs" block)
BEACH_SLOPE, BERM, MHW, GROUPING, MIN_DUR = 0.06, 1.7, 0.36, 24, 8
# -----------------------------------------------------------------------------


# A window's folder name, start_end
def wtag(w):
    return f"{w[0]}_{w[1]}"


# A variant's storm file
def storm_file(w, variant):
    return EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_{variant}.npy"


# The builder's own functions

# The builder's load_data, calculate_r2_percent and create_storms, without its module code
def builder_functions():
    tree = ast.parse(BUILDER.read_text(encoding="utf-8"))
    keep = [n for n in tree.body if isinstance(n, ast.FunctionDef)
            and n.name in ("find_time_gaps", "load_data", "calculate_r2_percent", "create_storms")]
    ns = {"np": np, "pd": pd, "os": os}
    from datetime import datetime, timedelta
    ns.update(datetime=datetime, timedelta=timedelta)
    exec(compile(ast.Module(body=keep, type_ignores=[]), str(BUILDER), "exec"), ns)
    return ns


# The builder's merged water-level and wave record for one window
def merged_record(fn, w):
    from site_layer import hat_env_forcings as env
    df = fn["load_data"](start_time=f"{w[0]}-01-01 00:00:00", end_time=f"{w[1]}-12-31 23:00:00",
                         water_levels_file=str(env.DUCK_GAUGE_FILE), wis_file=str(env.WIS_FILE),
                         t_name_water="t", water_name="v", t_name_wis="time",
                         waveHs_name="waveHs", waveTp_name="waveTp")
    df["R2"] = fn["calculate_r2_percent"](df["Hs"], df["Tp"], BEACH_SLOPE)
    df["TWL"] = df["water_level"] + df["R2"]
    return df


# Cut each event longer than `limit` to its `limit` hours around the peak
def trim_long_events(df, nocap_summary, w, limit=72):
    above = df[df["TWL"] > BERM]
    rows = []
    for _, ev in nocap_summary.iterrows():
        hrs = above[(above.index >= ev.StartTime) & (above.index <= ev.EndTime)]
        if len(hrs) > limit:
            k = int(np.argmax(hrs["TWL"].values))
            lo = min(max(0, k - limit // 2), len(hrs) - limit)
            hrs = hrs.iloc[lo:lo + limit]
        peak = hrs["TWL"].idxmax()
        rows.append(dict(calendar_year=ev.calendar_year, StartTime=hrs.index[0], EndTime=hrs.index[-1],
                         Rhigh=(hrs["TWL"].max() - MHW) / 10, Rlow=(hrs["TWL"].min() - MHW) / 10,
                         period=hrs.loc[peak, "Tp"], duration=len(hrs),
                         trimmed_from=int(ev.duration) if ev.duration > limit else 0))
    s = pd.DataFrame(rows)
    s["time"] = s["calendar_year"] - w[0] + 1
    return s


# Build every variant's series, checking the 72 h one against the committed file
def build():
    fn = builder_functions()
    for w in WINDOWS:
        out = EXP_DIR / "storms" / wtag(w)
        out.mkdir(parents=True, exist_ok=True)
        df = merged_record(fn, w)
        for v, mx in VARIANTS.items():
            if v == "72trim":
                continue
            fn["create_storms"](df_merged=df, berm_elevation=BERM, weather_grouping=GROUPING, MHW=MHW,
                                min_storm_dur=MIN_DUR, max_storm_dur=mx, save_dfs=True,
                                save_dir=str(out), save_name=f"{wtag(w)}_storms_v3_{v}",
                                window_start_year=w[0])
        nocap = pd.read_csv(out / f"{wtag(w)}_storms_v3_nocap_summary.csv", parse_dates=["StartTime", "EndTime"])
        s = trim_long_events(df, nocap, w)
        s.to_csv(out / f"{wtag(w)}_storms_v3_72trim_summary.csv", index=False)
        np.save(storm_file(w, "72trim"), s[["time", "Rhigh", "Rlow", "period", "duration"]].to_numpy())
        # THE CHECK: the 72 h rebuild is the committed series, or nothing here means anything
        from site_layer import hat_env_forcings as env
        committed = np.load(env.storm_series_file(*w))
        rebuilt = np.load(storm_file(w, "72"))
        same = committed.shape == rebuilt.shape and np.array_equal(committed, rebuilt)
        print(f"{wtag(w)}: 72 h rebuild == committed series: {same}")
        if not same:
            raise SystemExit("the rebuilt 72 h series differs from the committed one; stopping")
        for v in VARIANTS:
            a = np.load(storm_file(w, v))
            print(f"  {v:7s} {len(a):4d} storms, max duration {a[:, 4].max():4.0f} h, "
                  f"max Rhigh {a[:, 1].max() * 10:.2f} m MHW, storm-hours {a[:, 4].sum():6.0f}")


# Runs

# Subprocess entry: the unchanged runner with this period's storm file swapped in
def _launch(start, storm_path):
    import runpy
    sys.path.insert(0, str(HINDCAST.parent))
    from site_layer import hatteras_site_config as sc
    sc.HATTERAS_PERIODS[int(start)]["storm_file"] = str(storm_path)
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


# A run group's experiment tag
def member_tag(variant, scenario):
    return f"{TAG}/runs/{variant}_{scenario}"


# Pre-fix Barrier3D worktree (ce36866), for the cause test on the 72 h and 240 h series
PREFIX_WORKTREE = PROJECT_ROOT.parent / "Barrier3D-prefix-ce36866"


# Run one member in its own process, logged
def run_member(job, prefix=False):
    w, variant, scenario = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    who = f"prefix_{variant}" if prefix else variant
    log = logs / f"{who}_{scenario}_{wtag(w)}.log"
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": scenario,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": member_tag(who, scenario),
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "PYTHONIOENCODING": "utf-8"})
    if prefix:
        env["PYTHONPATH"] = str(PREFIX_WORKTREE) + os.pathsep + env.get("PYTHONPATH", "")
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), str(storm_file(w, variant))],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    text = log.read_text(encoding="utf-8", errors="replace")
    stopped = [ln.strip() for ln in text.splitlines() if "b3d_break" in ln or "DROWNED" in ln]
    b3d = next((ln.strip() for ln in text.splitlines() if ln.startswith("BARRIER3D =")), "")
    rec = dict(window=wtag(w), variant=who, scenario=scenario, returncode=p.returncode, barrier3d=b3d,
               minutes=round((time.time() - t0) / 60, 1), stopped=stopped,
               storm_file_in_log=str(storm_file(w, variant).name) in text)
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{wtag(w)} {who:12s} {scenario:16s} exit {p.returncode} {rec['minutes']} min "
          f"{'; '.join(stopped) if stopped else 'ran to the end'}")
    return rec


# Every variant x scenario x window, in parallel
def run(workers):
    for w in WINDOWS:
        for v in RUN_VARIANTS:
            if not storm_file(w, v).exists():
                raise SystemExit(f"missing {storm_file(w, v)}; run `build` first")
    jobs = [(w, v, s) for w in WINDOWS for v in RUN_VARIANTS for s in SCENARIOS]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(run_member, jobs))


# Run the 72 h and 240 h natural series on the pre-fix Barrier3D
def cause(workers):
    if not (PREFIX_WORKTREE / "barrier3d").exists():
        raise SystemExit(f"no pre-fix worktree at {PREFIX_WORKTREE}")
    jobs = [(w, v, "natural") for w in WINDOWS for v in ("72", "240")]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(lambda j: run_member(j, prefix=True), jobs))


# Diagnosis

# A member's run, or None
def run_dir(w, variant, scenario):
    base = EXP_DIR / "runs" / f"{variant}_{scenario}" / wtag(w) / "edgeBE"
    hits = list(base.glob("HAT_*")) if base.exists() else []
    return hits[0] if hits else None


# Per run: whether it drowned, where, when and how, and that year's storms
def diagnose():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rows = []
    for w in WINDOWS:
        summ72 = pd.read_csv(EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_72_summary.csv")
        for v in RUN_VARIANTS:
            summ = pd.read_csv(EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_{v}_summary.csv")
            for s in SCENARIOS:
                d = run_dir(w, v, s)
                if d is None:
                    rows.append(dict(window=wtag(w), variant=v, scenario=s, status="no run"))
                    continue
                npz = d / f"{d.name}.npz"
                if not npz.exists():
                    rows.append(dict(window=wtag(w), variant=v, scenario=s, status="no model state"))
                    continue
                c = np.load(npz, allow_pickle=True)["cascade"][0]
                years_run = len(c.barrier3d[0].x_s_TS) - 1
                drowned = [(p, b) for p, b in enumerate(c.barrier3d) if getattr(b, "drown_break", 0) == 1]
                base = dict(window=wtag(w), variant=v, scenario=s, years_run=years_run,
                            years_planned=w[1] - w[0], n_drowned=len(drowned))
                if not drowned:
                    rows.append(dict(base, status="ran to the end"))
                    continue
                for p, b in drowned:
                    t = len(b.x_s_TS) - 1                          # the year that broke
                    interior = np.asarray(b.InteriorDomain)
                    mech = ("width: the shoreline ate every interior row" if interior.shape[0] <= 0
                            else "height: every front-row cell at or below sea level")
                    yr = w[0] + t - 1
                    yv, y72 = summ[summ.calendar_year == yr], summ72[summ72.calendar_year == yr]
                    xs = np.asarray(b.x_s_TS) * 10
                    gis = DOM.pad_to_gis(p) if DOM.gis_to_pad(1) <= p <= DOM.gis_to_pad(90) else None
                    rows.append(dict(
                        base, status="drowned", pad=p, gis=gis if gis else f"buffer (pad {p})",
                        model_year=t, calendar_year=yr, mechanism=mech,
                        interior_rows_left=int(interior.shape[0]),
                        shoreline_retreat_that_year_m=float(xs[-1] - xs[-2]) if len(xs) > 1 else np.nan,
                        overwash_that_year_m3_per_m=float(np.asarray(b.QowTS)[-1]),
                        storms_that_year=len(yv), longest_h=int(yv.duration.max()) if len(yv) else 0,
                        max_rhigh_m=float(yv.Rhigh.max() * 10) if len(yv) else np.nan,
                        storm_hours=int(yv.duration.sum()),
                        storm_hours_in_72h_series=int(y72.duration.sum()),
                        max_rhigh_in_72h_series_m=float(y72.Rhigh.max() * 10) if len(y72) else np.nan))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "drowning.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 40):
        print(t.to_string(index=False))


# Each variant against its matrix control: change, overwash, skill, hit rate
def compare():
    import importlib.util
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    path = PROJECT_ROOT / "scripts" / "input_prep" / "8-overwash-analysis" / "4-vs-model" / "overwash_vs_model.py"
    spec = importlib.util.spec_from_file_location("overwash_vs_model", path)
    ovm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ovm)
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    rp = slice(DOM.gis_to_pad(1), DOM.gis_to_pad(90) + 1)

    def one(run_path, w, summary_csv):
        name = run_path.name
        c = np.load(run_path / f"{name}.npz", allow_pickle=True)["cascade"][0]
        meta = json.loads((run_path / f"{name}_run_metadata.json").read_text(encoding="utf-8"))
        xs = np.array([np.asarray(b.x_s_TS) for b in c.barrier3d]).T * 10
        q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
        st = pd.read_csv(summary_csv, parse_dates=["EndTime"])
        st["dated"] = st.EndTime - ovm.GRACE
        largest = st.loc[st.groupby("time").Rhigh.idxmax()].set_index("time").dated
        rows = []
        for _, im in ovm.images_in(w, obs).iterrows():
            o = obs[obs.Obs_ID == im.Obs_ID].set_index("domain").overwash
            for gis in range(1, 91):
                pp = DOM.gis_to_pad(gis)
                m = any(q[t, pp] > 0 and im["from"] < largest[t] <= im.date
                        for t in range(1, q.shape[0]) if t in largest)
                rows.append(dict(gis=gis, observed=o.get(gis, np.nan), model=int(m)))
        sc = ovm.scores(pd.DataFrame(rows))
        sk = meta.get("skill", {})
        return dict(mean_net_change_m=float(-(xs[-1] - xs[0])[rp].mean()),
                    domain_years_overwash=int((q[1:, rp] > 0).sum()),
                    total_overwash_m3_per_m=float(q[1:, rp].sum()),
                    obs_hit_rate=sc["hit_rate"], obs_both=sc["both"], obs_observed_only=sc["observed_only"],
                    obs_model_only=sc["model_only"],
                    rmse_interior_m_yr=float(sk.get("rmse_interior_m_yr", "nan")),
                    mean_bias_interior_m_yr=float(sk.get("mean_bias_interior_m_yr", "nan")))

    out = []
    for w in WINDOWS:
        for s in SCENARIOS:
            ctrl = MATRIX / wtag(w) / "edgeBE" / CONTROLS[(w[0], s)]
            summ = EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_72_summary.csv"
            out.append(dict(window=wtag(w), scenario=s, variant="72 (matrix)", **one(ctrl, w, summ)))
            for v in RUN_VARIANTS:
                d = run_dir(w, v, s)
                summ = EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_{v}_summary.csv"
                out.append(dict(window=wtag(w), scenario=s, variant=v, **one(d, w, summ)))
    t = pd.DataFrame(out)
    t.to_csv(EXP_DIR / "tables" / "comparison.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(t.round(3).to_string(index=False))


# Run: the subprocess entry, or the action asked for
def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2], sys.argv[3])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["build", "run", "cause", "diagnose", "compare"])
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    {"build": build, "run": lambda: run(a.workers), "cause": lambda: cause(a.workers),
     "diagnose": diagnose, "compare": compare}[a.action]()


if __name__ == "__main__":
    main()
