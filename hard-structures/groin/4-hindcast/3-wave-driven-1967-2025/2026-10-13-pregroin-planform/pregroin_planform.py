#!/usr/bin/env python3
"""
The corrected groin on the coast as it was before the groins: the 1967 wet/dry shoreline as the start.

    python pregroin_planform.py build                -> inputs/: the start offset, its table and check
    python pregroin_planform.py run [--streams 4]    -> no groin and strengths 0.10, 0.16, 0.25
    python pregroin_planform.py gif                  -> figures/pregroin_planform_panels.gif
    python pregroin_planform.py storms               -> inputs/: storms for 1968-1995
    python pregroin_planform.py run --full           -> the maintained life, 1968-1995, strengths 0-0.35
    python pregroin_planform.py score                -> full_life_scores.csv: the gap against the photo dates
    python pregroin_planform.py gif-full             -> figures/pregroin_full_life_panels.gif

The 1967 wet/dry line (August 1967, two years before the groins) is measured on the model's 100 m transects
the way every offset is (step1_line_to_raw_offsets.intersect). It covers about GIS 2-12, so the start is the
1996 shoreline offset (shoreline/v4) everywhere else, and over the covered domains the 1996 values plus
(1967 - 1996) less that difference's straight line between the covered ends, so the two join without a step.
Everything else is the calibration run's (1996 topography, storms 1997-2009, RSLR 0.002, edgeBE ends, no
fills): only the planform around Buxton is 1967's. The child process swaps the period's offset file for the
built one, and the corrected groin in, and drives the unchanged runner.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

import argparse
import os
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
CG_DIR = HERE.parent / "2026-10-09-corrected-groin"
sys.path[:0] = [str(CG_DIR)]
import corrected_groin as cgm  # noqa: E402

REPO = cgm.REPO
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "2-brie-offset" / "1-produce"))

# --- CONFIG ------------------------------------------------------------------
PERIOD = 1996                               # the calibration run's setup
STUDY = "groin/2026-10-13-pregroin-planform"
WETDRY = REPO / "hard-structures" / "groin" / "1-observations" / "data" / "gis" / "wet_dry_groin.geojson"
START_LINE_YEAR, TARGET_LINE_YEAR = 1967, 1978   # before the groins; after their build-up
MIN_CROSSINGS = 3                           # transects a domain needs on the 1967 line to be covered
PAD = 15                                    # buffer domains before GIS 1 in the padded file
OFFSET_OUT = HERE / "inputs" / "Island_Shoreline_Offsets_1996_with_1967_GIS2-12_PADDED_120.csv"
STRENGTHS = (0.0, 0.10, 0.16, 0.25)         # 0 is no groin: the corrected groin at 0 reproduces it exactly
LOGS = HERE / "logs"
# The full maintained life, `--full`: storms 1968-1995 (installed 1969, last repaired 1995), its own study
FULL_FIRST, FULL_LAST = 1968, 1995
FULL_RSLR = 0.004                           # the 1984 period's Duck fit (0.00391 m/yr over 1984-2004)
FULL_STUDY = "groin/2026-10-13-pregroin-full-life"
FULL_STRENGTHS = (0.0, 0.16, 0.20, 0.25, 0.30, 0.35)
STORMS_OUT = HERE / "inputs" / "storms_1968_1995_v3_split12_trim24.npy"
STAND_IN_FROM = 1984                        # 1968-1983 replay the recorded 1984-1999; 1984-1995 are their own years
# -----------------------------------------------------------------------------


def wetdry_line(year):
    import geopandas as gpd
    g = gpd.read_file(WETDRY)
    return g[g["year"] == year].geometry.iloc[0]


# Per-domain mean offset (m from the datum, smaller = seaward) of a line, and the transect rows behind it
def domain_offsets(line):
    import step1_line_to_raw_offsets as s1
    rows = s1.intersect(s1.load_transects(), line)
    ok = rows.dropna(subset=["ORIG_LEN"])
    n = ok.groupby("domain_id")["ORIG_LEN"].count()
    means = ok.groupby("domain_id")["ORIG_LEN"].mean()
    return means[n >= MIN_CROSSINGS], rows


def cmd_build(_a):
    import geopandas as gpd
    from site_layer import hat_topo_version as tv
    # Check the measurement on a line that has a raw file: the 1996 shoreline must come back exactly
    raw96 = pd.read_csv(tv.RAW_OFFSET_DIR / "1996-01-01_1996-12-31_shoreline_offset_raw.csv")
    line96 = gpd.read_file(REPO / "data" / "hatteras_init" / "5-scr" / "1-observations" / "mean_shoreline"
                           / "1996-01-01_1996-12-31" / "shoreline_mean_1996-01-01_1996-12-31.geojson").geometry.iloc[0]
    d96_check, _ = domain_offsets(line96)
    d96 = raw96.drop_duplicates("LineID").groupby("domain_id")["ORIG_LEN"].mean()
    worst = float((d96_check - d96.loc[d96_check.index]).abs().max())
    print(f"check: the 1996 shoreline re-measured here vs its raw file, worst domain {worst:.2e} m")
    assert worst < 1e-6

    d67, rows67 = domain_offsets(wetdry_line(START_LINE_YEAR))
    d78, _ = domain_offsets(wetdry_line(TARGET_LINE_YEAR))
    gis = d67.index.to_numpy()
    print(f"{START_LINE_YEAR} wet/dry line covers GIS {gis.min()}-{gis.max()} ({len(gis)} domains)")
    padded = pd.read_csv(tv.offset_file(PERIOD, "padded", source="shoreline"))
    col = padded.columns[0]
    v = padded[col].to_numpy(dtype=float)
    # The 1967 shape over the covered domains, joined to 1996 at both ends: add (1967 - 1996) less its end-to-end line
    diff = (d67 - d96.loc[gis]).to_numpy()
    join = np.interp(gis, [gis[0], gis[-1]], [diff[0], diff[-1]])
    new = v.copy()
    new[PAD + gis - 1] = v[PAD + gis - 1] + diff - join
    OFFSET_OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({col: new}).to_csv(OFFSET_OUT, index=False)
    # The 1978 line in the same frame, with its own end-to-end line taken off the same way: the target shape
    g78 = [g for g in gis if g in d78.index]
    diff78 = (d78.loc[g78] - d67.loc[g78]).to_numpy()
    join78 = np.interp(g78, [g78[0], g78[-1]], [diff78[0], diff78[-1]])
    t = pd.DataFrame({"gis": gis, "offset_1996_m": v[PAD + gis - 1], f"wetdry_{START_LINE_YEAR}_m": d67.to_numpy(),
                      f"wetdry_{TARGET_LINE_YEAR}_m": d78.reindex(gis).to_numpy(),
                      "start_offset_m": new[PAD + gis - 1]})
    t[f"shape_change_{START_LINE_YEAR}_{TARGET_LINE_YEAR}_m"] = np.nan
    t.loc[t.gis.isin(g78), f"shape_change_{START_LINE_YEAR}_{TARGET_LINE_YEAR}_m"] = diff78 - join78
    t.round(2).to_csv(HERE / "inputs" / "start_offset_table.csv", index=False)
    rows67.to_csv(HERE / "inputs" / f"wetdry_{START_LINE_YEAR}_transects.csv", index=False)
    print(t.round(1).to_string(index=False))
    for name, arr in (("1996 start", v[PAD + gis - 1]), ("1967 start", new[PAD + gis - 1])):
        st = -np.diff(arr)
        print(f"seaward step to the next domain, {name}: "
              + " ".join(f"{a}->{a + 1}: {s:.0f}" for a, s in zip(gis[:-1], st)))
    print("wrote", OFFSET_OUT.relative_to(REPO))


# 1968-1995 storms on one 1-based model-year index: recorded 1984-1995 on their own years, and the recorded
# 1984-1999 replayed for 1968-1983, before the storm record starts
def cmd_storms(_a):
    from site_layer.hat_env_forcings import storm_series_file
    src = np.load(storm_series_file(1984, 2004))             # column 0: model year, 1 = 1984
    parts = []
    for k, year in enumerate(range(FULL_FIRST, FULL_LAST + 1), start=1):
        source = year if year >= STAND_IN_FROM else STAND_IN_FROM + (year - FULL_FIRST)
        rows = src[src[:, 0] == source - 1983].copy()
        rows[:, 0] = k
        parts.append(rows)
    out = np.vstack(parts)
    np.save(STORMS_OUT, out)
    print(f"wrote {STORMS_OUT.relative_to(REPO)}: {int(out[:, 0].max())} years, {len(out)} storms "
          f"({len(src[src[:, 0] <= 21])} in the 21-year source)")


def spec(s):
    return f"1969:{s:g}"


def tag(s, full=False):
    return f"{FULL_STUDY if full else STUDY}/periods_{cgm.periods_label(spec(s))}"


def shoreline(s, full=False):
    if not full:
        return cgm.shoreline_file(tag(s), PERIOD, s, 1.0)
    hits = sorted((cgm.RAW / "experiments" / tag(s, True) / cgm.dg.member(s, 1.0))
                  .glob(f"{PERIOD}_*/edgeBE/*/*_shoreline_matrix.npy"))
    return hits[-1] if hits else None


# Child: the period's offset file swapped for the built one (and, for the full life, its years, storms and
# sea level), then corrected_groin's own launch
def _launch():
    from site_layer import hatteras_site_config as hsc
    p = hsc.HATTERAS_PERIODS[PERIOD]
    p["island_offset_file"] = str(OFFSET_OUT)                # absolute: Path / abs is abs
    if os.environ.get("PREGROIN_FULL") == "1":
        # end_year stays 2009, so the period's label and its CoastSat report target still load (unused here);
        # only the years stepped change, so last_model_year's end_year check is relaxed for this period
        # The storm file relative to the runner's data root: its metadata records the path relative to that root
        storms = os.path.relpath(STORMS_OUT, REPO / "data" / "hatteras_init")
        p.update(first_model_year=FULL_FIRST, last_model_year=FULL_LAST,
                 storm_file=storms, sea_level_rise_rate=FULL_RSLR, enable_nourishment=False)
        strict = hsc.last_model_year
        hsc.last_model_year = lambda y: FULL_LAST if int(y) == PERIOD else strict(y)
    cgm._launch()


def run(s, full=False):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "CORRECTED_GROIN_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(PERIOD),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(s), "HAT_GROIN_DETERIORATION_FRACTION": "1.0",
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{tag(s, full)}/{cgm.dg.member(s, 1.0)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
        "CORRECTED_GROIN_SCHEDULE": "periods", "CORRECTED_GROIN_OFF": spec(s),
        "CORRECTED_GROIN_REFERENCE": "tracked", "PREGROIN_FULL": "1" if full else "0",
    })
    LOGS.mkdir(exist_ok=True)
    t0 = time.time()
    with open(LOGS / f"{'full_' if full else ''}s{s:g}.log", "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    f = shoreline(s, full)
    if f is not None:
        (f.parent / "pregroin_planform.txt").write_text(
            f"patched in-process by pregroin_planform.py: island offset {OFFSET_OUT.name} (the 1996 shoreline with "
            f"the {START_LINE_YEAR} wet/dry shape over GIS 2-12), not the metadata's shoreline/v4; the corrected "
            f"groin at a constant strength {s:g} from 1969 (tracked reference)"
            + (f"; years {FULL_FIRST}-{FULL_LAST}, storms {STORMS_OUT.name} (1968-1983 replay 1984-1999), RSLR "
               f"{FULL_RSLR}, no fills; the 1996 topography" if full else "")
            + ". The runner's report does not show these.")
    print(f"{time.strftime('%H:%M:%S')}  s {s:g}: exit {code} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def cmd_run(a):
    full = a.full
    strengths = FULL_STRENGTHS if full else STRENGTHS
    todo = [s for s in strengths if shoreline(s, full) is None]
    print(f"{len(strengths) - len(todo)} of {len(strengths)} done; running {len(todo)} in {a.streams} streams")
    jobs, failed = queue.Queue(), []
    for s in todo:
        jobs.put(s)

    def worker():
        while True:
            try:
                s = jobs.get_nowait()
            except queue.Empty:
                return
            if run(s, full) or shoreline(s, full) is None:
                failed.append(s)

    threads = [threading.Thread(target=worker) for _ in range(a.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}")
    x = np.load(shoreline(0.0, full))
    off = pd.read_csv(OFFSET_OUT).iloc[:, 0].to_numpy()
    i = np.arange(PAD + 1, PAD + 12)
    print(f"start check: alongshore steps over GIS 2-12, run row 0 vs the built offset, worst "
          f"{np.abs(np.diff(x[0, i]) - np.diff(off[i])).max():.1f} m")


def cmd_gif(_a):
    import groin_panels_gif as GP
    from site_layer.hat_figure_style import apply_style
    apply_style()
    t = pd.read_csv(HERE / "inputs" / "start_offset_table.csv").set_index("gis")
    gis = GP.GIS
    idx = gis + PAD - 1
    model = {s: np.load(shoreline(s))[:, idx] for s in STRENGTHS}
    n = len(model[0.0])
    x0_raw = model[0.0][0]
    shape = t[f"shape_change_{START_LINE_YEAR}_{TARGET_LINE_YEAR}_m"].reindex(gis).to_numpy()
    target_raw = x0_raw + shape                     # offsets and the model's x are both landward +
    lines = GP.rotate([model[s] for s in STRENGTHS] + [x0_raw, target_raw], x0_raw)
    *runs_m, x0, target = lines
    runs = [("No groin" if s == 0 else f"Corrected groin, strength {s:g}", m, False,
             None if s == 0 else (lambda k: True)) for s, m in zip(STRENGTHS, runs_m)]
    out = HERE / "figures" / "pregroin_planform_panels.gif"
    out.parent.mkdir(exist_ok=True)
    GP.render(runs, x0, target, f"{START_LINE_YEAR} shoreline, before the groins (start)",
              f"{TARGET_LINE_YEAR} observed shoreline shape (overall retreat removed)",
              [f"year {k}" for k in range(n)],
              f"The corrected groin on the pre-groin coast: {START_LINE_YEAR} shoreline at GIS 2-12, "
              f"groin in from year 0",
              lambda k: f"storms of {PERIOD + 1}-{PERIOD + k}" if k else "start", out)


# A wet/dry line's change in shape since 1967 over GIS 2-12 (landward +), its end-to-end line taken off as in build
def shape_change(year, gis):
    d67, _ = domain_offsets(wetdry_line(START_LINE_YEAR))
    dy, _ = domain_offsets(wetdry_line(year))
    g = [x for x in gis if x in d67.index and x in dy.index]
    diff = (dy.loc[g] - d67.loc[g]).to_numpy()
    out = pd.Series(diff - np.interp(g, [g[0], g[-1]], [diff[0], diff[-1]]), index=g)
    return out.reindex(gis).to_numpy()


# Full life: the gap across the groins against the wet/dry photo dates, 1970-1995 (row k is the end of 1967 + k)
def cmd_score(_a):
    from score_instant_grid import observed_series
    photos = observed_series()
    photos = photos - photos.loc[START_LINE_YEAR]
    yrs = [y for y in photos.index if START_LINE_YEAR < y <= FULL_LAST]
    rows = []
    for s in FULL_STRENGTHS:
        m = np.load(shoreline(s, True))
        g = m[:, cgm.sr.DOWN] - m[:, cgm.sr.UP]
        g = g - g[0]
        mv = np.array([g[y - START_LINE_YEAR] for y in yrs])
        ov = photos.loc[yrs].to_numpy()
        no95 = [i for i, y in enumerate(yrs) if y != 1995]
        rows.append(dict(strength=s, rmse_photos=np.sqrt(np.mean((mv - ov) ** 2)),
                         rmse_photos_no1995=np.sqrt(np.mean((mv[no95] - ov[no95]) ** 2)),
                         bias=np.mean(mv - ov), **{f"gap_{y}": v for y, v in zip(yrs, mv)}))
    t = pd.DataFrame(rows)
    t.round(2).to_csv(HERE / "full_life_scores.csv", index=False)
    print("observed gap change since 1967 (photos):", " ".join(f"{y}:{v:.0f}" for y, v in photos.loc[yrs].items()))
    print(t.round(1).to_string(index=False))


def cmd_gif_full(_a):
    import groin_panels_gif as GP
    from site_layer.hat_figure_style import apply_style
    apply_style()
    gis = GP.GIS
    idx = gis + PAD - 1
    model = {s: np.load(shoreline(s, True))[:, idx] for s in FULL_STRENGTHS}
    x0_raw = model[0.0][0]
    target_raw = x0_raw + shape_change(FULL_LAST, gis)
    shown = [0.0, 0.16, 0.25, 0.35]
    *runs_m, x0, target = GP.rotate([model[s] for s in shown] + [x0_raw, target_raw], x0_raw)
    years = list(range(START_LINE_YEAR, FULL_LAST + 1))
    runs = [("No groin" if s == 0 else f"Corrected groin, strength {s:g} from 1969", m, False,
             None if s == 0 else (lambda k: years[k] >= 1969)) for s, m in zip(shown, runs_m)]
    out = HERE / "figures" / "pregroin_full_life_panels.gif"
    GP.render(runs, x0, target, f"{START_LINE_YEAR} shoreline, before the groins (start)",
              f"{FULL_LAST} observed shoreline shape (overall retreat removed)", [str(y) for y in years],
              f"The maintained groin, {START_LINE_YEAR}-{FULL_LAST}: {START_LINE_YEAR} shoreline at GIS 2-12",
              lambda k: ("groins built 1969" if years[k] >= 1969 else "before the groins")
              + ("; storms replayed from 1984-1999" if years[k] < 1984 else "; recorded storms"), out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("build")
    sub.add_parser("storms")
    r = sub.add_parser("run")
    r.add_argument("--streams", type=int, default=4)
    r.add_argument("--full", action="store_true", help=f"the maintained life, {FULL_FIRST}-{FULL_LAST}")
    sub.add_parser("gif")
    sub.add_parser("score")
    sub.add_parser("gif-full")
    a = ap.parse_args()
    {"build": cmd_build, "storms": cmd_storms, "run": cmd_run, "gif": cmd_gif, "score": cmd_score,
     "gif-full": cmd_gif_full}[a.cmd](a)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch()
    else:
        main()
