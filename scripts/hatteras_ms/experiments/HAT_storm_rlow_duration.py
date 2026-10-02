"""
Do Rlow and duration, defined as Barrier3D expects them, improve the model?

    python scripts/hatteras_ms/experiments/HAT_storm_rlow_duration.py build
    python scripts/hatteras_ms/experiments/HAT_storm_rlow_duration.py run
    python scripts/hatteras_ms/experiments/HAT_storm_rlow_duration.py score

Four split12 storm series (control, rlow, dureq, rlow_dureq), run managed and
natural in both windows, scored against the overwash imagery and CoastSat. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from __future__ import annotations

import argparse
import contextlib
import io
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
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_max_duration as MD  # noqa: E402
import HAT_storm_length_selection as S  # noqa: E402
import HAT_storm_event_splitting as ES  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-and-overwash/2026-10-01-rlow-and-duration"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
VARIANTS = ("split12", "rlow", "dureq", "rlow_dureq")
SCENARIOS = ("full_management", "natural")
TRIM_H = 24              # the adopted cap, kept wherever duration is not dureq
SPLIT_GAP_H = 12         # the adopted split rule
FLOW_EXPONENT = 1.5      # Barrier3D gap discharge goes as Rexcess^1.5 (Vdune * Rexcess)
WATCH = [("Fran 1996", "1996-09-05", "1996-09-07", 1996), ("Isabel 2003", "2003-09-18", "2003-09-19", 1996),
         ("Irene 2011", "2011-08-27", "2011-08-28", 2010), ("Sandy 2012", "2012-10-27", "2012-10-30", 2010),
         ("Mar 2018 nor'easter", "2018-03-03", "2018-03-05", 2010)]
# -----------------------------------------------------------------------------


# A window's folder name, start_end
def wtag(w):
    return f"{w[0]}_{w[1]}"


# A variant's storm file and its summary CSV (split12 is the adopted file)
def storm_paths(w, v):
    if v == "split12":
        from site_layer import hat_env_forcings as env
        f = env.storm_series_file(*w, variant="v3_split12_trim24")
    else:
        f = EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_split12_{v}.npy"
    return f, f.with_name(f.stem + "_summary.csv")


# Stockdon (2006) total swash, the S the builder's R2 is built from
def swash(hs, tp, slope=MD.BEACH_SLOPE):
    l0 = 9.81 * tp ** 2 / (2 * np.pi)
    s_inc = 0.75 * slope * np.sqrt(hs * l0)
    s_ig = 0.06 * np.sqrt(hs * l0)
    return np.sqrt(s_inc ** 2 + s_ig ** 2)


# Hours at the peak level that move the same Rexcess^1.5 flow over the berm as the whole event
def equivalent_duration(twl):
    ex = np.clip(twl - MD.BERM, 0, None)
    return max(1, int(round(float((ex ** FLOW_EXPONENT).sum() / ex.max() ** FLOW_EXPONENT))))


# One event row; the variant decides Rlow and the duration
def event_row(hrs, w, v):
    full = len(hrs)
    kept, trimmed = hrs, 0
    if full > TRIM_H:
        k = int(np.argmax(hrs["TWL"].values))
        lo = min(max(0, k - TRIM_H // 2), full - TRIM_H)
        kept, trimmed = hrs.iloc[lo:lo + TRIM_H], full
    # MSSM (multivariateSeaStorm.m): hourly Rlow = TWL - S/2, the storm's Rlow is its maximum
    rlow = (hrs["TWL"] - hrs["S"] / 2).max() if "rlow" in v else kept["TWL"].min()
    if "dureq" in v:
        dur, trimmed = equivalent_duration(hrs["TWL"].values), 0
    else:
        hrs, dur = kept, len(kept)
    peak = hrs["TWL"].idxmax()
    start = hrs.index[0]
    return dict(calendar_year=start.year, StartTime=start, EndTime=hrs.index[-1],
                Rhigh=(hrs["TWL"].max() - MD.MHW) / 10, Rlow=(rlow - MD.MHW) / 10,
                period=hrs.loc[peak, "Tp"], duration=dur, hours_above_berm=full, trimmed_from=trimmed,
                time=start.year - w[0] + 1)


# One variant's series from the untrimmed 24 h systems, split as split12
def series(df, systems, w, v):
    above = df[df["TWL"] > MD.BERM]
    rows = []
    for _, ev in systems.iterrows():
        hrs = above[(above.index >= ev.StartTime) & (above.index <= ev.EndTime)]
        for g in ES._pieces(hrs.index, SPLIT_GAP_H, MD.MIN_DUR):
            rows.append(event_row(hrs.iloc[g], w, v))
    return pd.DataFrame(rows).sort_values("StartTime").reset_index(drop=True)


# Write one series as .npy, .csv and its summary
def _save(s, w, v):
    npy, summ = storm_paths(w, v)
    npy.parent.mkdir(parents=True, exist_ok=True)
    cols = ["time", "Rhigh", "Rlow", "period", "duration"]
    s.to_csv(summ, index=False)
    s[cols].to_csv(npy.with_suffix(".csv"), index=False)
    np.save(npy, s[cols].to_numpy())


# Build the three test series, after the control rebuilds the adopted file exactly
def build():
    fn = MD.builder_functions()
    rows = []
    for w in MD.WINDOWS:
        out = EXP_DIR / "storms" / wtag(w)
        out.mkdir(parents=True, exist_ok=True)
        with contextlib.redirect_stdout(io.StringIO()):
            df = MD.merged_record(fn, w)
            fn["create_storms"](df_merged=df, berm_elevation=MD.BERM, weather_grouping=MD.GROUPING, MHW=MD.MHW,
                                min_storm_dur=MD.MIN_DUR, max_storm_dur=10 ** 6, save_dfs=True, save_dir=str(out),
                                save_name=f"{wtag(w)}_systems_untrimmed", window_start_year=w[0])
        df["S"] = swash(df["Hs"], df["Tp"])
        systems = pd.read_csv(out / f"{wtag(w)}_systems_untrimmed_summary.csv", parse_dates=["StartTime", "EndTime"])
        # THE CHECK: the control rule must reproduce the adopted split12 series exactly
        cols = ["time", "Rhigh", "Rlow", "period", "duration"]
        mine = series(df, systems, w, "split12")[cols].to_numpy()
        adopted = np.load(storm_paths(w, "split12")[0])
        same = mine.shape == adopted.shape and np.allclose(mine, adopted, rtol=0, atol=1e-12)
        print(f"{wtag(w)}: rebuilt split12 == adopted v3_split12_trim24: {same}")
        if not same:
            raise SystemExit("the rebuild differs from the adopted series; stopping")
        last = w[1] - 1 if w[0] == 1996 else w[1]      # the calendar years the run spends
        for v in VARIANTS:
            s = series(df, systems, w, v)
            if v != "split12":
                _save(s, w, v)
            s = s[s.calendar_year <= last]
            gap = 10 * (s.Rhigh - s.Rlow)
            row = dict(window=wtag(w), variant=v, events=len(s), storm_hours=int(s.duration.sum()),
                       max_duration_h=int(s.duration.max()), median_duration_h=float(s.duration.median()),
                       median_rlow_m_mhw=round(10 * s.Rlow.median(), 2), max_rlow_m_mhw=round(10 * s.Rlow.max(), 2),
                       median_rhigh_minus_rlow_m=round(float(gap.median()), 2))
            for lab, a, b, ws in WATCH:
                if ws == w[0]:
                    hit = s[(s.EndTime >= a) & (s.StartTime <= pd.Timestamp(b) + pd.Timedelta(days=1))]
                    if len(hit):
                        e = hit.loc[hit.Rhigh.idxmax()]
                        row[lab] = f"Rh {10 * e.Rhigh:.2f} / Rl {10 * e.Rlow:.2f} / {int(e.duration)} h"
                    else:
                        row[lab] = "absent"
            rows.append(row)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "series.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        print(t.to_string(index=False))


# A variant's finished run, or None
def run_dir(w, v, scenario):
    base = EXP_DIR / "runs" / f"{v}_{scenario}" / wtag(w) / "edgeBE"
    # complete runs only: the metadata is written last
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# Run one job in its own process, logged, via this file's _launch entry
def launch(job):
    w, v, scenario = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": scenario,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{v}_{scenario}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{v}_{scenario}_{wtag(w)}.log"
    with open(log, "w", encoding="utf-8") as fh:
        # Relative to data/hatteras_init (through ..), or the runner fails writing its metadata
        rel = os.path.relpath(storm_paths(w, v)[0], PROJECT_ROOT / "data" / "hatteras_init")
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), rel],
                           env=env, cwd=MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    b3d = next((ln.strip() for ln in log.read_text(encoding="utf-8", errors="replace").splitlines()
                if ln.startswith("BARRIER3D =")), "")
    rec = dict(window=wtag(w), variant=v, scenario=scenario, returncode=p.returncode, barrier3d=b3d,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{wtag(w)} {v:10s} {scenario:16s} exit {p.returncode} {rec['minutes']} min  {b3d}", flush=True)


# Score every run: overwash against the imagery, shoreline skill, which regime Barrier3D used
def score():
    import HAT_dune_ceiling_per_domain as P
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    rows = []
    for scenario in SCENARIOS:
        for w in MD.WINDOWS:
            for v in VARIANTS:
                rd = run_dir(w, v, scenario)
                if rd is None:
                    print(f"  missing {wtag(w)} {v} {scenario}")
                    continue
                c = S.load_state(rd)
                summ = pd.read_csv(storm_paths(w, v)[1])
                cells = S.overwash_cells(c, w, summ, obs, ovm)
                b3d = [c.barrier3d[p] for p in pads]
                inun = sum(int(b._InundationCount) for b in b3d)
                runup = sum(int(b._RunUpCount) for b in b3d)
                rows.append(dict(scenario=scenario, window=wtag(w), storms=v,
                                 **S.overwash_scores(cells, 0.0), **P.shoreline_lowess7(rd, w),
                                 overwash_total_m3_per_m=float(sum(np.asarray(b.QowTS)[1:].sum() for b in b3d)),
                                 inundation_storms=inun, runup_storms=runup,
                                 inundation_share=inun / max(inun + runup, 1)))
                print(f"  scored {wtag(w)} {v} {scenario}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    cols = ["scenario", "window", "storms", "POD", "POFD", "PSS", "timing_r", "space_r",
            "rmse_interior_m_yr", "bias_interior_m_yr", "overwash_total_m3_per_m", "inundation_share"]
    with pd.option_context("display.width", 240, "display.max_columns", 20):
        print(t[[x for x in cols if x in t.columns]].round(2).to_string(index=False))


# Run: the subprocess entry, or the action asked for
def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        MD._launch(*sys.argv[2:4])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["build", "run", "score"])
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    if a.action == "build":
        build()
    elif a.action == "run":
        jobs = [(w, v, sc) for sc in SCENARIOS for w in MD.WINDOWS for v in VARIANTS if run_dir(w, v, sc) is None]
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
