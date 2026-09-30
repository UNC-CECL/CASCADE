r"""
HAT_storm_height_test.py -- are the storm water levels too low for Hatteras?
==============================================================================
WHY (storms-and-overwash/2026-09-28-dune-ceiling-per-domain): once the dunes
match the 2009 lidar, the low spots overwash where Irene really did. But
Irene also overwashed 47 of the 60 domains with dunes of 4.3 m and up, and the
model manages 10: a modelled Irene Rhigh of 3.80 m MHW cannot clear them.
The storm series takes water level from the Duck gauge (8651370), about 80 km
north of the reach, and adds Stockdon (2006) R2% run-up on WIS ST63228 waves
with one beach slope, 0.06.

PART 1 -- observations (no model runs)
    gauges   peak water level (m above MHW) at Duck against the gauges in or
             near the reach, for the named storms 1996-2024. The ocean-side
             record inside the reach is the Cape Hatteras Fishing Pier
             (8654400, historic); Oregon Inlet Marina (8652587) and USCG
             Station Hatteras (8654467) sit inside the inlets, on the sound
             side. Fetched from the NOAA CO-OPS API into this folder's data/.
    slope    the foreshore slope the run-up should use, measured from the
             domains' 10 m elevation profiles (first land cell to the dune
             toe), against the 0.06 in the storm builder
PART 2 -- sensitivity runs
    storm variants of the trim24 series: run-up slope 0.08 and 0.10
    (Stockdon recomputed, events re-found), and Rhigh/Rlow raised by 0.25 and
    0.5 m (a local surge Duck does not see; the events are unchanged)
    x dune ceilings uniform 5.5 m NAVD88 and per-cell (the two candidates)
    x managed, both windows. Controls: the trim24 runs of those two ceilings.

NOTHING IN THE MAIN CODE CHANGES (storm files and ceilings are swapped in
each run's own process, as in the earlier experiments).

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-height/

USAGE
    python HAT_storm_height_test.py gauges
    python HAT_storm_height_test.py slope
    python HAT_storm_height_test.py build
    python HAT_storm_height_test.py run [--workers 6]
    python HAT_storm_height_test.py score
==============================================================================

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_length_selection as S  # noqa: E402
import HAT_dune_ceiling_rebuild as E  # noqa: E402
import HAT_dune_ceiling_per_domain as P  # noqa: E402

TAG = "storms-and-overwash/2026-09-28-storm-height"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
DATA = EXP_DIR / "data"

STATIONS = {"8651370": "Duck (ocean, the storm file's gauge)",
            "8654400": "Cape Hatteras Fishing Pier (ocean, in the reach)",
            "8652226": "Jennette's Pier (ocean, Nags Head)",
            "8652587": "Oregon Inlet Marina (sound side)",
            "8654467": "USCG Station Hatteras (sound side)",
            "8653215": "Rodanthe, Pamlico Sound (sound side)"}
EVENTS = {"Fran 1996": "1996-09-05", "Bonnie 1998": "1998-08-27", "Dennis 1999": "1999-09-04",
          "Floyd 1999": "1999-09-16", "Isabel 2003": "2003-09-18", "Ophelia 2005": "2005-09-14",
          "Nov 2006 nor'easter": "2006-11-22", "Nor'Ida 2009": "2009-11-13", "Earl 2010": "2010-09-03",
          "Irene 2011": "2011-08-27", "Sandy 2012": "2012-10-29", "Arthur 2014": "2014-07-04",
          "Joaquin 2015": "2015-10-03", "Matthew 2016": "2016-10-09", "Mar 2018 nor'easter": "2018-03-04",
          "Florence 2018": "2018-09-14", "Dorian 2019": "2019-09-06", "Isaias 2020": "2020-08-04",
          "Ian 2022": "2022-09-30", "Lee 2023": "2023-09-16"}


# --- PART 1: gauges ------------------------------------------------------------------

def fetch(station, begin, end):
    """Hourly heights (historical stations: 'hourly_height'; otherwise
    'water_level'), m above MHW. Cached under data/."""
    DATA.mkdir(parents=True, exist_ok=True)
    cache = DATA / f"{station}_{begin}_{end}.csv"
    if cache.exists():
        df = pd.read_csv(cache)
        df["t"] = pd.to_datetime(df["t"], errors="coerce")
        return df.dropna(subset=["t"])
    for product in ("hourly_height", "water_level"):
        url = ("https://api.tidesandcurrents.noaa.gov/api/prod/datagetter?"
               f"product={product}&application=CASCADE_Hatteras&begin_date={begin}&end_date={end}"
               f"&datum=MHW&station={station}&time_zone=GMT&units=metric&format=csv")
        try:
            raw = urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception:
            continue
        if raw.lower().startswith("date time"):
            df = pd.read_csv(io.StringIO(raw))
            df = df.rename(columns={df.columns[0]: "t", df.columns[1]: "v"})[["t", "v"]]
            df["t"] = pd.to_datetime(df["t"], errors="coerce")
            df["v"] = pd.to_numeric(df["v"], errors="coerce")
            df = df.dropna(subset=["t"])
            if len(df):
                df.to_csv(cache, index=False)
                return df
    pd.DataFrame(columns=["t", "v"]).to_csv(cache, index=False)
    return pd.DataFrame(columns=["t", "v"])


def gauges():
    rows = []
    for ev, day in EVENTS.items():
        d = pd.Timestamp(day)
        b, e = (d - pd.Timedelta(days=3)).strftime("%Y%m%d"), (d + pd.Timedelta(days=3)).strftime("%Y%m%d")
        row = dict(event=ev, date=day)
        for st in STATIONS:
            df = fetch(st, b, e)
            v = df["v"].dropna() if len(df) else pd.Series(dtype=float)
            row[st] = float(v.max()) if len(v) else np.nan
        rows.append(row)
        print(f"  {ev}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "gauge_peaks_m_above_mhw.csv", index=False)
    show = t.rename(columns={k: v.split(" (")[0] for k, v in STATIONS.items()})
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(show.round(2).to_string(index=False))
    for st in ("8654400", "8652587", "8654467", "8653215", "8652226"):
        both = t[["8651370", st]].dropna()
        if len(both):
            d = both[st] - both["8651370"]
            print(f"{STATIONS[st]:48s} minus Duck: median {d.median():+.2f} m over {len(both)} events "
                  f"(range {d.min():+.2f} to {d.max():+.2f})")


# --- PART 2: sensitivity runs ------------------------------------------------------

VARIANTS = ("slope0p10", "plus0p25", "plus0p50")
CEILINGS = ("uniform5p5", "cell")


def storm_file(w, v):
    return EXP_DIR / "storms" / S.wtag(w) / f"{S.wtag(w)}_storms_{v}.npy"


def summary_csv(w, v):
    return EXP_DIR / "storms" / S.wtag(w) / f"{S.wtag(w)}_storms_{v}_summary.csv"


def build():
    """slope0p10: the builder's chain with run-up slope 0.10, events re-found
    (every event kept), trimmed to 24 h. plusX: the trim24 series with Rhigh
    and Rlow raised X m (the same events)."""
    fn = S.MD.builder_functions()
    for w in S.WINDOWS:
        out = EXP_DIR / "storms" / S.wtag(w)
        out.mkdir(parents=True, exist_ok=True)
        df = S.MD.merged_record(fn, w)
        df["R2"] = fn["calculate_r2_percent"](df["Hs"], df["Tp"], 0.10)
        df["TWL"] = df["water_level"] + df["R2"]
        fn["create_storms"](df_merged=df, berm_elevation=S.MD.BERM, weather_grouping=S.MD.GROUPING, MHW=S.MD.MHW,
                            min_storm_dur=S.MD.MIN_DUR, max_storm_dur=10 ** 6, save_dfs=True, save_dir=str(out),
                            save_name=f"{S.wtag(w)}_storms_slope0p10full", window_start_year=w[0])
        full = pd.read_csv(out / f"{S.wtag(w)}_storms_slope0p10full_summary.csv", parse_dates=["StartTime", "EndTime"])
        tr = S.MD.trim_long_events(df, full, w, limit=24)
        tr.to_csv(summary_csv(w, "slope0p10"), index=False)
        np.save(storm_file(w, "slope0p10"), tr[["time", "Rhigh", "Rlow", "period", "duration"]].to_numpy())
        base = pd.read_csv(S.summary_csv(w, "trim24"))
        for v, off in (("plus0p25", 0.25), ("plus0p50", 0.50)):
            b = base.copy()
            b["Rhigh"] += off / 10
            b["Rlow"] += off / 10
            b.to_csv(summary_csv(w, v), index=False)
            np.save(storm_file(w, v), b[["time", "Rhigh", "Rlow", "period", "duration"]].to_numpy())
        for v in ("trim24",) + VARIANTS:
            a = np.load(S.storm_file(w, v) if v == "trim24" else storm_file(w, v))
            print(f"{S.wtag(w)} {v:10s} {len(a):4d} storms  Rhigh median {np.median(a[:, 1]) * 10:.2f}  "
                  f"max {a[:, 1].max() * 10:.2f} m MHW  storm-hours {a[:, 4].sum():.0f}")


def group(v, ceiling):
    return f"{v}_{ceiling}_full_management"


def run_dir(w, v, ceiling):
    if v == "trim24":
        return (E.run_dir(w, "trim24", 5.5, "current", "full_management") if ceiling == "uniform5p5"
                else P.run_dir(w, "trim24", "cell", "full_management"))
    base = EXP_DIR / "runs" / group(v, ceiling) / S.wtag(w) / "edgeBE"
    hits = sorted(base.glob("HAT_*")) if base.exists() else []
    return hits[-1] if hits else None


def _launch(start, storm_path, ceiling):
    if ceiling == "uniform5p5":
        E._launch(start, storm_path, 5.5, "current")
    else:
        P._launch(start, storm_path, "cell")


def launch(job):
    w, v, ceiling = job
    g = group(v, ceiling)
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{g}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    with open(logs / f"{g}_{S.wtag(w)}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), str(storm_file(w, v)), ceiling],
                           env=env, cwd=S.MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(window=S.wtag(w), variant=v, ceiling=ceiling, returncode=p.returncode,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{S.wtag(w)} {g:40s} exit {p.returncode} {rec['minutes']} min", flush=True)
    return rec


def run(workers):
    jobs = [(w, v, c) for w in S.WINDOWS for v in VARIANTS for c in CEILINGS]
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(launch, [j for j in jobs if run_dir(*j) is None]))


def score():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    P.install_cell_growth()
    S.MATRIX = P.ARCHIVED_MATRIX
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    lidar = E.crest_series(S.load_state(S.MATRIX / "2010_2024/edgeBE" / S.CONTROLS[(2010, "full_management")]), pads)[0]
    low = lidar < np.percentile(lidar, 33)
    rows = []
    for w in S.WINDOWS:
        for c in CEILINGS:
            for v in ("trim24",) + VARIANTS:
                rd = run_dir(w, v, c)
                if rd is None:
                    print(f"  missing {S.wtag(w)} {v} {c}")
                    continue
                summ = pd.read_csv(S.summary_csv(w, v) if v == "trim24" else summary_csv(w, v))
                cells = S.overwash_cells(S.load_state(rd), w, summ, obs, ovm)
                row = dict(window=S.wtag(w), ceiling=c, storms=v, **S.overwash_scores(cells, 0.0),
                           **P.shoreline_lowess7(rd, w))
                if w[0] == 2010:
                    ir = cells[cells.obs_id == "OBS-018"].set_index("gis").reindex(range(1, 91))
                    m, o = (ir.model_m3_per_m > 0).values, ir.observed.values
                    for name, sel in (("low", low), ("rest", ~low)):
                        k = sel & ~np.isnan(o)
                        row[f"irene_{name}"] = f"{int(m[k].sum())}/{int(np.nansum(o[k]))} obs"
                rows.append(row)
                print(f"  scored {S.wtag(w)} {c} {v}", flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(EXP_DIR / "tables" / "sensitivity_scores.csv", index=False)
    cols = ["window", "ceiling", "storms", "POD", "POFD", "PSS", "timing_r", "space_r",
            "rmse_interior_m_yr", "bias_interior_m_yr", "irene_low", "irene_rest"]
    with pd.option_context("display.width", 220, "display.max_columns", 20):
        print(t[[x for x in cols if x in t.columns]].round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:5])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["gauges", "build", "run", "score"])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    {"gauges": gauges, "build": build, "run": lambda: run(a.workers), "score": score}[a.action]()


if __name__ == "__main__":
    main()
