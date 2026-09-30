r"""
HAT_storm_event_splitting.py -- should back-to-back storms be separate events?
==============================================================================
WHY (Hannah, 2026-09-29: "test splitting the events"). The storm builder joins
above-berm spells less than 24 h apart into one event (weather_grouping = 24).
At Hatteras the berm is overtopped at most high tides during an active spell,
so two storms a week apart chain into one long event: Edouard + Fran 1996
(one event, 119 h above the berm) and Jose + Maria 2017 (186 h). The adopted
series (v3_trim24) then keeps the 24 h around the event's single highest peak,
so Fran and Jose are not in the model at all, and 1996-2024 loses 56 spells of
>= 8 h above the berm inside merged events (41 peak above 2 m MHW).

THE CANDIDATES (every one trims each event to 24 h, as adopted)
    trim24   the adopted series (hindcast_storms/*_v3_trim24): the control.
             Rebuilt here from the builder's functions and checked identical.
    g12      the builder with weather_grouping = 12 h: the one-number change.
             It recovers Fran and Jose, but a merged storm's tidal fragments
             shorter than 8 h then fall under the minimum-duration rule and
             are dropped, so it has FEWER events and storm-hours than trim24.
    split12  the 24 h grouping kept to define a weather system, then the system
             split wherever the water stays below the berm for >= 12 h; a
             piece shorter than 8 h is folded into the piece before it (or
             after, for the first), so no hour the adopted series counts is
             lost. Each piece is then trimmed to 24 h and dated by its start.

RUNS: managed (full_management), both windows, edgeBE, the site config's end
rates (not re-solved, as in the trim-length check), the unchanged runner with
the period's storm file swapped in its own process (HAT_storm_max_duration).
SCORES: as the trim-length check -- overwash against the imagery (POD, POFD,
PSS, timing and space r), interior RMSE/bias against LOWESS-7, total overwash.

NOTHING IN THE MAIN CODE CHANGES.

    python HAT_storm_event_splitting.py build   # series + what each recovers
    python HAT_storm_event_splitting.py run     # 6 runs, 5 at a time
    python HAT_storm_event_splitting.py score

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-29-event-splitting/
==============================================================================
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

TAG = "storms-and-overwash/2026-09-29-event-splitting"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
WORKTREE = PROJECT_ROOT.parent / "Barrier3D"          # hatteras/adopted is checked out there
VARIANTS = ("trim24", "g12", "split12")
TRIM_H = 24
SPLIT_GAP_H = 12
# storms to report on: (label, first day, last day, window start)
WATCH = [("Edouard 1996", "1996-08-30", "1996-09-02", 1996), ("Fran 1996", "1996-09-05", "1996-09-07", 1996),
         ("Dennis 1999", "1999-08-29", "1999-09-05", 1996), ("Isabel 2003", "2003-09-18", "2003-09-19", 1996),
         ("Jose 2017", "2017-09-18", "2017-09-20", 2010), ("Maria 2017", "2017-09-26", "2017-09-27", 2010),
         ("Mar 2018 nor'easter", "2018-03-03", "2018-03-05", 2010)]


def wtag(w):
    return f"{w[0]}_{w[1]}"


def storm_paths(w, v):
    """(npy, summary csv) for a variant."""
    if v == "trim24":
        from site_layer import hat_env_forcings as env
        f = env.storm_series_file(*w, variant="v3_trim24")
    else:
        f = EXP_DIR / "storms" / wtag(w) / f"{wtag(w)}_storms_v3_{v}.npy"
    return f, f.with_name(f.stem + "_summary.csv")


def _pieces(hours, gap_h, min_h):
    """Split one system's above-berm hours (a sorted DatetimeIndex) at gaps
    >= gap_h; fold a piece shorter than min_h into its neighbour."""
    gaps = np.diff(hours.values).astype("timedelta64[h]").astype(int)
    cuts = np.where(gaps >= gap_h)[0] + 1
    pieces = [list(p) for p in np.split(np.arange(len(hours)), cuts)]
    i = 0
    while len(pieces) > 1 and i < len(pieces):
        if len(pieces[i]) < min_h:
            j = i - 1 if i > 0 else i + 1
            pieces[j] = sorted(pieces[j] + pieces[i])
            del pieces[i]
            i = 0
            continue
        i += 1
    return pieces


def _event_row(hrs, w, trimmed_from):
    """One event from its kept above-berm hours, computed as the builder does."""
    peak = hrs["TWL"].idxmax()
    start = hrs.index[0]
    return dict(calendar_year=start.year, StartTime=start, EndTime=hrs.index[-1],
                Rhigh=(hrs["TWL"].max() - MD.MHW) / 10, Rlow=(hrs["TWL"].min() - MD.MHW) / 10,
                period=hrs.loc[peak, "Tp"], duration=len(hrs), trimmed_from=trimmed_from,
                time=start.year - w[0] + 1)


def split_series(df, systems, w, gap_h):
    """systems: the untrimmed 24 h-grouped events (the builder's 'full' run).
    gap_h None = no split (reproduces trim24)."""
    above = df[df["TWL"] > MD.BERM]
    rows = []
    for _, ev in systems.iterrows():
        hrs = above[(above.index >= ev.StartTime) & (above.index <= ev.EndTime)]
        groups = [list(range(len(hrs)))] if gap_h is None else _pieces(hrs.index, gap_h, MD.MIN_DUR)
        for g in groups:
            h = hrs.iloc[g]
            trimmed = 0
            if len(h) > TRIM_H:
                k = int(np.argmax(h["TWL"].values))
                lo = min(max(0, k - TRIM_H // 2), len(h) - TRIM_H)
                trimmed, h = len(h), h.iloc[lo:lo + TRIM_H]
            rows.append(_event_row(h, w, trimmed))
    s = pd.DataFrame(rows).sort_values("StartTime").reset_index(drop=True)
    if gap_h is None:                       # the builder dates by the SYSTEM's start year
        s["calendar_year"] = systems.calendar_year.values
        s["time"] = s.calendar_year - w[0] + 1
    return s


def _save(s, w, v):
    npy, summ = storm_paths(w, v)
    npy.parent.mkdir(parents=True, exist_ok=True)
    s.to_csv(summ, index=False)
    arr = s[["time", "Rhigh", "Rlow", "period", "duration"]].to_numpy()
    s[["time", "Rhigh", "Rlow", "period", "duration"]].to_csv(npy.with_suffix(".csv"), index=False)
    np.save(npy, arr)
    return arr


def build():
    fn = MD.builder_functions()
    rows = []
    for w in MD.WINDOWS:
        out = EXP_DIR / "storms" / wtag(w)
        out.mkdir(parents=True, exist_ok=True)
        with contextlib.redirect_stdout(io.StringIO()):
            df = MD.merged_record(fn, w)
            # the untrimmed 24 h-grouped systems
            fn["create_storms"](df_merged=df, berm_elevation=MD.BERM, weather_grouping=MD.GROUPING, MHW=MD.MHW,
                                min_storm_dur=MD.MIN_DUR, max_storm_dur=10 ** 6, save_dfs=True, save_dir=str(out),
                                save_name=f"{wtag(w)}_systems_untrimmed", window_start_year=w[0])
            # g12: the builder itself, 12 h grouping, trim 24
            fn["create_storms"](df_merged=df, berm_elevation=MD.BERM, weather_grouping=12, MHW=MD.MHW,
                                min_storm_dur=MD.MIN_DUR, max_storm_dur=TRIM_H, save_dfs=True, save_dir=str(out),
                                save_name=f"{wtag(w)}_storms_v3_g12", window_start_year=w[0], long_events="trim")
        systems = pd.read_csv(out / f"{wtag(w)}_systems_untrimmed_summary.csv", parse_dates=["StartTime", "EndTime"])
        # THE CHECK: no split must reproduce the adopted series exactly
        rebuilt = split_series(df, systems, w, None)
        committed = np.load(storm_paths(w, "trim24")[0])
        mine = rebuilt[["time", "Rhigh", "Rlow", "period", "duration"]].to_numpy()
        same = mine.shape == committed.shape and np.allclose(mine, committed, rtol=0, atol=1e-12)
        print(f"{wtag(w)}: rebuilt trim24 == adopted v3_trim24: {same}")
        if not same:
            raise SystemExit("the rebuild differs from the adopted series; stopping")
        _save(split_series(df, systems, w, SPLIT_GAP_H), w, "split12")
        last = w[1] - 1 if w[0] == 1996 else w[1]      # the calendar years the run spends
        for v in VARIANTS:
            s = pd.read_csv(storm_paths(w, v)[1], parse_dates=["StartTime", "EndTime"])
            s = s[s.calendar_year <= last]
            row = dict(window=wtag(w), variant=v, events=len(s), storm_hours=int(s.duration.sum()),
                       trimmed=int((s.get("trimmed_from", 0) > 0).sum()),
                       max_rhigh_m_mhw=round(10 * s.Rhigh.max(), 2))
            for lab, a, b, ws in WATCH:
                if ws == w[0]:
                    hit = s[(s.EndTime >= a) & (s.StartTime <= pd.Timestamp(b) + pd.Timedelta(days=1))]
                    row[lab] = f"{10 * hit.Rhigh.max():.2f}" if len(hit) else "absent"
            rows.append(row)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "series.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print(t.to_string(index=False))


def run_dir(w, v):
    base = EXP_DIR / "runs" / f"{v}_full_management" / wtag(w) / "edgeBE"
    # complete runs only: the metadata is written last
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    w, v = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{v}_full_management",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{v}_{wtag(w)}.log"
    with open(log, "w", encoding="utf-8") as fh:
        # The runner records the storm file relative to data/hatteras_init and
        # fails at the end of the run on a path outside it (2026-09-29: four
        # runs lost their metadata and shoreline matrix that way). A relative
        # path with ".." resolves to the same file and satisfies it.
        rel = os.path.relpath(storm_paths(w, v)[0], PROJECT_ROOT / "data" / "hatteras_init")
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), rel],
                           env=env, cwd=MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    b3d = next((ln.strip() for ln in log.read_text(encoding="utf-8", errors="replace").splitlines()
                if ln.startswith("BARRIER3D =")), "")
    rec = dict(window=wtag(w), variant=v, returncode=p.returncode, barrier3d=b3d,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{wtag(w)} {v:8s} exit {p.returncode} {rec['minutes']} min  {b3d}", flush=True)


def score():
    import HAT_dune_ceiling_per_domain as P
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    rows = []
    for w in MD.WINDOWS:
        for v in VARIANTS:
            rd = run_dir(w, v)
            if rd is None:
                print(f"  missing {wtag(w)} {v}")
                continue
            c = S.load_state(rd)
            summ = pd.read_csv(storm_paths(w, v)[1])
            cells = S.overwash_cells(c, w, summ, obs, ovm)
            npy = np.load(storm_paths(w, v)[0])
            rows.append(dict(window=wtag(w), storms=v, events=len(npy), storm_hours=int(npy[:, 4].sum()),
                             **S.overwash_scores(cells, 0.0), **P.shoreline_lowess7(rd, w),
                             overwash_total_m3_per_m=float(np.array(
                                 [np.asarray(b.QowTS)[1:].sum() for b in c.barrier3d])[pads].sum())))
            print(f"  scored {wtag(w)} {v}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    cols = ["window", "storms", "events", "storm_hours", "POD", "POFD", "PSS", "timing_r", "space_r",
            "rmse_interior_m_yr", "bias_interior_m_yr", "overwash_total_m3_per_m"]
    with pd.option_context("display.width", 240, "display.max_columns", 20):
        print(t[[x for x in cols if x in t.columns]].round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        MD._launch(*sys.argv[2:4])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["build", "run", "score"])
    a = ap.parse_args()
    if a.action == "build":
        build()
    elif a.action == "run":
        jobs = [(w, v) for w in MD.WINDOWS for v in VARIANTS if run_dir(w, v) is None]
        with ThreadPoolExecutor(5) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
