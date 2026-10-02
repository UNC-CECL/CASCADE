"""
Does storm run-up from each domain's own beach slope improve the overwash skill?

    python scripts/hatteras_ms/experiments/HAT_storm_runup_by_slope.py slopes
    python scripts/hatteras_ms/experiments/HAT_storm_runup_by_slope.py build
    python scripts/hatteras_ms/experiments/HAT_storm_runup_by_slope.py run
    python scripts/hatteras_ms/experiments/HAT_storm_runup_by_slope.py score

Measures the foreshore slope per domain on each period's lidar, rebuilds every
storm's Rhigh/Rlow per domain from it (same events, same hours), sets each
domain's series in-process, and scores overwash and shoreline. Details: scripts/hatteras_ms/experiments/README.md.

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
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "8-overwash-analysis" / "4-vs-model"))
import HAT_storm_max_duration as MD  # noqa: E402
import HAT_storm_event_splitting as ES  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-and-overwash/2026-10-01-runup-by-slope"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
# Both windows read the 2009 lidar: the 1996 ALACE beach graft gives 19 unmeasurable domains and slopes to 0.72
DEM_PRODUCT = {1996: "2009-2014", 2010: "2009-2014"}
SLOPE_RANGE = (0.01, 0.16)  # Stockdon et al. (2006)'s foreshore slopes; outside it the fit is extrapolated
SMOOTH_DOMAINS = 5          # running median alongshore, so one scarped profile does not set a domain's storms
MHW_NAVD = 0.36
MIN_ROWS = 20            # rows with a clean MHW-to-berm profile a domain needs
TRIM_H, SPLIT_GAP_H = 24, 12            # the adopted series
VARIANTS = ("control", "measured", "pattern")
SCENARIOS = ("full_management", "natural")
COLS = ["time", "Rhigh", "Rlow", "period", "duration"]
# -----------------------------------------------------------------------------


def wtag(w):
    return f"{w[0]}_{w[1]}"


# Foreshore slope

# Median MHW-to-berm slope of one 1 m domain tile, corrected for the shoreline's angle to the rows
def domain_slope(tif):
    import rasterio
    with rasterio.open(tif) as r:
        a = r.read(1).astype(float)
        a[a == r.nodata] = np.nan
    xs, ss = [], []
    for row in a:
        v = np.where(np.isfinite(row))[0]
        if len(v) < 10 or row[v[-1]] > MHW_NAVD + 0.3:      # no wet edge on the ocean side
            continue
        j = v[-1]
        while j > v[0] and row[j] < MHW_NAVD:
            j -= 1
        k = j
        while k > v[0] and row[k] < MD.BERM:
            k -= 1
        if j - k > 0:
            xs.append(j)
            ss.append((MD.BERM - MHW_NAVD) / (j - k))
    if len(ss) < MIN_ROWS:
        return np.nan, np.nan, len(ss)
    theta = np.degrees(np.arctan(abs(np.polyfit(np.arange(len(xs)), xs, 1)[0])))
    return float(np.median(ss) / np.cos(np.radians(theta))), float(theta), len(ss)


def slopes():
    root = PROJECT_ROOT / "data" / "hatteras_init" / "0-elevation"
    rows = []
    for start, prod in DEM_PRODUCT.items():
        for gis in range(1, 91):
            s, th, n = domain_slope(root / prod / "1-gapfill-1m" / f"clip_domain_{gis}_filled.tif")
            rows.append(dict(start=start, product=prod, gis=gis, slope=s, shore_angle_deg=th, rows=n))
    t = pd.DataFrame(rows)
    for start, g in t.groupby("start"):
        s = g.set_index("gis").slope.clip(*SLOPE_RANGE)
        s = s.rolling(SMOOTH_DOMAINS, center=True, min_periods=1).median().fillna(s.median())
        t.loc[g.index, "slope_filled"] = s.values
        t.loc[g.index, "slope_pattern"] = s.values * MD.BEACH_SLOPE / s.median()
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "beach_slope_by_domain.csv", index=False)
    print(t.groupby("start")[["slope", "slope_filled", "slope_pattern"]].describe(percentiles=[0.1, 0.5, 0.9]).round(3).T)
    print("domains without a measurement:", t[t.slope.isna()][["start", "gis"]].values.tolist())


def slope_table():
    return pd.read_csv(EXP_DIR / "tables" / "beach_slope_by_domain.csv")


# Storm series per domain

# Stockdon R2 for one slope
def r2(hs, tp, beta):
    l0 = 9.81 * tp ** 2 / (2 * np.pi)
    sq = np.sqrt(hs * l0)
    return 1.1 * (0.35 * beta * sq + np.sqrt((0.75 * beta * sq) ** 2 + (0.06 * sq) ** 2) / 2)


# The adopted events as (kept hours, calendar year) pairs, found on the island series
def events(w):
    fn = MD.builder_functions()
    out = EXP_DIR / "storms" / wtag(w)
    out.mkdir(parents=True, exist_ok=True)
    with contextlib.redirect_stdout(io.StringIO()):
        df = MD.merged_record(fn, w)
        fn["create_storms"](df_merged=df, berm_elevation=MD.BERM, weather_grouping=MD.GROUPING, MHW=MD.MHW,
                            min_storm_dur=MD.MIN_DUR, max_storm_dur=10 ** 6, save_dfs=True, save_dir=str(out),
                            save_name=f"{wtag(w)}_systems_untrimmed", window_start_year=w[0])
    systems = pd.read_csv(out / f"{wtag(w)}_systems_untrimmed_summary.csv", parse_dates=["StartTime", "EndTime"])
    above = df[df["TWL"] > MD.BERM]
    evs = []
    for _, ev in systems.iterrows():
        hrs = above[(above.index >= ev.StartTime) & (above.index <= ev.EndTime)]
        for g in ES._pieces(hrs.index, SPLIT_GAP_H, MD.MIN_DUR):
            h = hrs.iloc[g]
            if len(h) > TRIM_H:
                k = int(np.argmax(h["TWL"].values))
                lo = min(max(0, k - TRIM_H // 2), len(h) - TRIM_H)
                h = h.iloc[lo:lo + TRIM_H]
            evs.append(h)
    evs.sort(key=lambda h: h.index[0])
    return df, evs


# One domain's series: the same events and hours, run-up from its slope
def series_for(df, evs, w, beta):
    rows = []
    for h in evs:
        twl = df.loc[h.index, "water_level"] + r2(df.loc[h.index, "Hs"], df.loc[h.index, "Tp"], beta)
        peak = h["TWL"].idxmax()                  # the island peak dates the wave period, as adopted
        rows.append([h.index[0].year - w[0] + 1, (twl.max() - MD.MHW) / 10, (twl.min() - MD.MHW) / 10,
                     h.loc[peak, "Tp"], len(h)])
    return np.array(rows, dtype=float)


def series_file(w, v):
    return EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_by_domain_{v}.npz"


def build():
    from site_layer import hat_env_forcings as env
    st = slope_table()
    rows = []
    for w in MD.WINDOWS:
        df, evs = events(w)
        adopted = np.load(env.storm_series_file(*w, variant=env.DEFAULT_STORM_VARIANT))
        sl = st[st.start == w[0]].set_index("gis")
        for v in VARIANTS:
            per = {}
            for gis in range(1, 91):
                beta = MD.BEACH_SLOPE if v == "control" else sl.at[gis, "slope_filled" if v == "measured" else "slope_pattern"]
                per[str(gis)] = series_for(df, evs, w, beta)
            if v == "control":
                # THE CHECK: slope 0.06 everywhere must be the adopted series, for every domain
                same = all(a.shape == adopted.shape and np.allclose(a, adopted, rtol=0, atol=1e-12) for a in per.values())
                print(f"{wtag(w)}: control == adopted {env.DEFAULT_STORM_VARIANT} in all 90 domains: {same}")
                if not same:
                    raise SystemExit("the per-domain rebuild differs from the adopted series; stopping")
            np.savez(series_file(w, v), **per)
            rh = np.array([10 * per[str(g)][:, 1].max() for g in range(1, 91)])
            yearly = np.array([pd.Series(per[str(g)][:, 1]).groupby(per[str(g)][:, 0]).max().median() * 10
                               for g in range(1, 91)])
            rows.append(dict(window=wtag(w), variant=v, max_rhigh_min=rh.min(), max_rhigh_max=rh.max(),
                             typical_year_max_p10=np.percentile(yearly, 10), typical_year_max_median=np.median(yearly),
                             typical_year_max_p90=np.percentile(yearly, 90)))
    t = pd.DataFrame(rows).round(2)
    t.to_csv(EXP_DIR / "tables" / "series.csv", index=False)
    print(t.to_string(index=False))


# Runs

# Child process: set each domain's own series after the model is built, then run the unchanged hindcast
def _launch(start, storm_path, variant):
    import cascade_pipeline.hindcast as H
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    w = next(x for x in MD.WINDOWS if x[0] == int(start))
    per = np.load(series_file(w, variant))
    orig = H.build_cascade

    def build_cascade(*a, **k):
        c = orig(*a, **k)
        for gis in range(1, 91):
            c.barrier3d[DOM.gis_to_pad(gis)].StormSeries = per[str(gis)]
        print(f"RUNUP_BY_SLOPE = {variant}: storm series set per domain for GIS 1-90", flush=True)
        return c

    H.build_cascade = build_cascade
    sys.argv = [sys.argv[0]]
    MD._launch(start, storm_path)


def run_dir(w, v, sc):
    base = EXP_DIR / "runs" / f"{v}_{sc}" / wtag(w) / "edgeBE"
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    from site_layer import hat_env_forcings as env
    w, v, sc = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    e = dict(os.environ)
    e.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
              "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": sc,
              "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{v}_{sc}",
              "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
              "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    rel = os.path.relpath(env.storm_series_file(*w, variant=env.DEFAULT_STORM_VARIANT),
                          PROJECT_ROOT / "data" / "hatteras_init")
    t0 = time.time()
    log = logs / f"{v}_{sc}_{wtag(w)}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), rel, v],
                           env=e, cwd=MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    text = log.read_text(encoding="utf-8", errors="replace")
    b3d = next((ln.strip() for ln in text.splitlines() if ln.startswith("BARRIER3D =")), "")
    rec = dict(window=wtag(w), variant=v, scenario=sc, returncode=p.returncode, barrier3d=b3d,
               series_set="RUNUP_BY_SLOPE" in text, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{wtag(w)} {v:9s} {sc:16s} exit {p.returncode} set={rec['series_set']} {rec['minutes']} min  {b3d}",
          flush=True)


# Scores

def score():
    import HAT_dune_ceiling_per_domain as P
    import storms_vs_overwash as SV
    from site_layer import hat_overwash as ow
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    zone = lambda g: "cape_point" if g <= 6 else ("community" if (7 <= g <= 8 or 21 <= g <= 31 or 68 <= g <= 83)  # noqa: E731
                                                  else "road")
    rows, by_dom = [], []
    for sc in SCENARIOS:
        for w in MD.WINDOWS:
            for v in VARIANTS:
                rd = run_dir(w, v, sc)
                if rd is None:
                    print(f"  missing {wtag(w)} {v} {sc}")
                    continue
                c = np.load(rd / f"{rd.name}.npz", allow_pickle=True)["cascade"][0]
                q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
                orig = SV.ovm.load_qow
                SV.ovm.load_qow = lambda window, arm: (q, rd.name)
                try:
                    d = SV.cells(w, obs)
                finally:
                    SV.ovm.load_qow = orig
                sk = SV.skill(d)
                z = d.assign(z=d.gis.map(zone))
                fa = z[z.outcome == "false_alarm"].groupby("z").size()
                ms = z[z.outcome == "miss"].groupby("z").size()
                per_dom = d.groupby("gis").outcome.value_counts().unstack(fill_value=0)
                for gis, r in per_dom.iterrows():
                    by_dom.append(dict(scenario=sc, window=wtag(w), variant=v, gis=gis, **r.to_dict()))
                rows.append(dict(scenario=sc, window=wtag(w), variant=v, hit_rate=sk["hit_rate"],
                                 false_alarm_rate=sk["false_alarm_rate"], skill=sk["skill"], hits=sk["hits"],
                                 misses=sk["misses"], false_alarms=sk["false_alarms"],
                                 fa_cape=int(fa.get("cape_point", 0)), fa_village=int(fa.get("community", 0)),
                                 fa_road=int(fa.get("road", 0)), miss_cape=int(ms.get("cape_point", 0)),
                                 miss_village=int(ms.get("community", 0)), miss_road=int(ms.get("road", 0)),
                                 space_r=_space_r(d), **P.shoreline_lowess7(rd, w)))
                print(f"  scored {wtag(w)} {v} {sc}", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    pd.DataFrame(by_dom).to_csv(EXP_DIR / "tables" / "outcomes_by_domain.csv", index=False)
    with pd.option_context("display.width", 280, "display.max_columns", 30):
        print(t.round(3).to_string(index=False))


# Alongshore r between the share of images each domain overwashed, observed and modelled
def _space_r(d):
    g = d.groupby("gis").agg(o=("observed", "mean"), m=("model", "mean"))
    return float(np.corrcoef(g.o, g.m)[0, 1]) if g.o.std() > 0 and g.m.std() > 0 else np.nan


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:5])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["slopes", "build", "run", "score"])
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    if a.action == "slopes":
        slopes()
    elif a.action == "build":
        build()
    elif a.action == "run":
        jobs = [(w, v, sc) for sc in SCENARIOS for w in MD.WINDOWS for v in VARIANTS if run_dir(w, v, sc) is None]
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
