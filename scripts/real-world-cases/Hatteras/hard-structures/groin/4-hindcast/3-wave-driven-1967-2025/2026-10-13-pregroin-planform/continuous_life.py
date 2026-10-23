#!/usr/bin/env python3
"""
The groin's whole life in one run, 1967-2025, with no reset to the observed shoreline at 1996 or 2009.

    python continuous_life.py storms                 -> inputs/: storms for 1968-2025
    python continuous_life.py run [--streams 4]      -> no groin and the stepped strength histories
    python continuous_life.py score                  -> continuous_scores.csv; the gap by era
    python continuous_life.py figure                 -> figures/continuous_gap_1967_2025.png
    python continuous_life.py gif                    -> figures/continuous_life_panels.gif
    python continuous_life.py fit                    -> continuous_fit.csv: the strength era by era, each held
    python continuous_life.py followups              -> continuous_followups.csv: repair, storm seeds, ramp-in

Do the strengths fitted window by window (about 0.25 maintained to 1995, 0.16 for 1996-2009, 0.10 for
2010-2016, off from the 2017 fill) join up as one history? Each window's fit restarts from that window's
observed shoreline; this run does not. Start: the 1967 wet/dry shoreline at GIS 2-12 joined to the 1996
offset elsewhere (pregroin_planform.py build). Storms: 1968-1983 replay the recorded 1984-1999, 1984-1995
from the 1984-2004 series, 1996-2025 the 1996-2025 series the calibration and test windows were sliced from.
RSLR 0.004 throughout (the runner takes one rate). Fills on (Rodanthe 2014, Buxton 2017 and 2022, Avon 2022;
the 1971 and 1973 Buxton fills are not in the model's list). The 1996 topography from the start.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-13
"""
from __future__ import annotations

# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
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

import pregroin_planform as pp  # noqa: E402

cgm, REPO, HERE = pp.cgm, pp.REPO, pp.HERE

# --- CONFIG ------------------------------------------------------------------
FIRST, LAST = 1968, 2025                    # storm years; row 0 is the 1967 start
RSLR = 0.004
STUDY = "groin-wave-driven/2026-10-13-continuous-1967-2025"
STORMS = workspace_twin(HERE) / "inputs" / "storms_1968_2025_v3_split12_trim24.npy"
STAND_IN_FROM = 1984                        # 1968-1983 replay the recorded 1984-1999
SPECS = {                                   # label: strength from each year on
    "no groin": "1969:0",
    "stepped, 0.25 maintained": "1969:0.25,1996:0.16,2010:0.1,2017:0",
    "stepped, 0.20 maintained": "1969:0.2,1996:0.16,2010:0.1,2017:0",
}
ERAS = (("photos 1970-1995", 1968, 1995), ("CoastSat 1997-2009", 1997, 2009),
        ("CoastSat 2010-2016", 2010, 2016), ("CoastSat 2017-2025", 2017, 2025))
LOGS = workspace_twin(HERE) / "logs" / "continuous"
# -----------------------------------------------------------------------------


def cmd_storms(_a):
    from site_layer.hat_env_forcings import storm_series_file
    early = np.load(storm_series_file(1984, 2004))           # model year 1 = 1984
    late = np.load(storm_series_file(1996, 2025))            # model year 1 = 1996
    parts = []
    for k, year in enumerate(range(FIRST, LAST + 1), start=1):
        if year >= 1996:
            rows = late[late[:, 0] == year - 1995].copy()
        else:
            source = year if year >= STAND_IN_FROM else STAND_IN_FROM + (year - FIRST)
            rows = early[early[:, 0] == source - 1983].copy()
        rows[:, 0] = k
        parts.append(rows)
    out = np.vstack(parts)
    np.save(STORMS, out)
    print(f"wrote {STORMS.relative_to(REPO)}: {int(out[:, 0].max())} years, {len(out)} storms")


def strength0(spec):
    return cgm.parse_periods(spec)[0][1]


# A run is filed by its strength history, or a short name for a long one (a ramp), under the study, or under
# the study suffixed with a storm variant ("seed1": 1968-1983 resampled instead of replayed)
def tag(spec, name=None, variant=""):
    study = f"{STUDY}-{variant}" if variant else STUDY
    return f"{study}/periods_{name or cgm.periods_label(spec)}"


def shoreline(spec, name=None, variant=""):
    hits = sorted((cgm.RAW / "experiments" / tag(spec, name, variant) / cgm.dg.member(strength0(spec), 1.0))
                  .glob(f"{pp.PERIOD}_*/edgeBE/*/*_shoreline_matrix.npy"))
    return hits[-1] if hits else None


def storms_file(variant=""):
    return STORMS if not variant else STORMS.with_name(f"storms_1968_2025_{variant}_v3_split12_trim24.npy")


# The default storms with 1968-1983 drawn instead as whole recorded years from 1984-2004, at random (seed)
def build_resampled(seed):
    from site_layer.hat_env_forcings import storm_series_file
    early = np.load(storm_series_file(1984, 2004))
    base = np.load(STORMS)
    rng = np.random.default_rng(seed)
    parts = []
    for k, year in enumerate(range(FIRST, STAND_IN_FROM), start=1):
        rows = early[early[:, 0] == int(rng.integers(1, 22))].copy()
        rows[:, 0] = k
        parts.append(rows)
    parts.append(base[base[:, 0] > STAND_IN_FROM - FIRST])
    out = np.vstack(parts)
    f = storms_file(f"seed{seed}")
    np.save(f, out)
    print(f"wrote {f.name}: {int(out[:, 0].max())} years, {len(out)} storms", flush=True)


# Child: the 1967 offset, the 1968-2025 years, storms and sea level, fills on; then corrected_groin's launch
def _launch():
    from site_layer import hatteras_site_config as hsc
    p = hsc.HATTERAS_PERIODS[pp.PERIOD]
    # end_year stays 2009, so the period's label and CoastSat report target load (unused); the storm file is
    # given relative to the runner's data root, where its metadata records it
    p.update(island_offset_file=str(pp.OFFSET_OUT), first_model_year=FIRST, last_model_year=LAST,
             storm_file=os.path.relpath(os.environ.get("CONTINUOUS_STORMS", str(STORMS)),
                                        REPO / "data" / "hatteras_init"), sea_level_rise_rate=RSLR,
             enable_nourishment=True)
    strict = hsc.last_model_year
    hsc.last_model_year = lambda y: LAST if int(y) == pp.PERIOD else strict(y)
    cgm._launch()


def run(label, spec, name=None, variant=""):
    s = strength0(spec)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "CORRECTED_GROIN_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(pp.PERIOD),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(s), "HAT_GROIN_DETERIORATION_FRACTION": "1.0",
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{tag(spec, name, variant)}/{cgm.dg.member(s, 1.0)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
        "CORRECTED_GROIN_SCHEDULE": "periods", "CORRECTED_GROIN_OFF": spec,
        "CORRECTED_GROIN_REFERENCE": "tracked", "CONTINUOUS_STORMS": str(storms_file(variant)),
    })
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(LOGS / f"{(variant + '_') if variant else ''}{name or cgm.periods_label(spec)}.log", "w",
              encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    f = shoreline(spec, name, variant)
    if f is not None:
        (f.parent / "continuous_life.txt").write_text(
            f"patched in-process by continuous_life.py: island offset {pp.OFFSET_OUT.name} (the 1996 shoreline with "
            f"the 1967 wet/dry shape over GIS 2-12), years {FIRST}-{LAST}, storms {storms_file(variant).name} "
            f"({'1968-1983 resampled from 1984-2004, ' + variant if variant else '1968-1983 replay 1984-1999'}), RSLR {RSLR}, fills on, the 1996 topography; the corrected groin with strength {spec} "
            "(tracked reference). The runner's report, labels and metadata show the 1996 period's instead.")
    print(f"{time.strftime('%H:%M:%S')}  {label}: exit {code} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def cmd_run(a):
    todo = [(lab, sp) for lab, sp in SPECS.items() if shoreline(sp) is None]
    print(f"{len(SPECS) - len(todo)} of {len(SPECS)} done; running {len(todo)} in {a.streams} streams", flush=True)
    jobs, failed = queue.Queue(), []
    for c in todo:
        jobs.put(c)

    def worker():
        while True:
            try:
                lab, sp = jobs.get_nowait()
            except queue.Empty:
                return
            if run(lab, sp) or shoreline(sp) is None:
                failed.append(lab)

    threads = [threading.Thread(target=worker) for _ in range(a.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}")


# The model gap across the groins (GIS 6 seaward of GIS 5) by calendar year, change since 1967
def model_gap(spec, name=None, variant=""):
    m = np.load(shoreline(spec, name, variant))
    g = m[:, cgm.sr.DOWN] - m[:, cgm.sr.UP]
    return pd.Series(g - g[0], index=range(FIRST - 1, FIRST - 1 + len(g)))


def observed():
    from score_instant_grid import observed_series
    photos = observed_series()
    photos = photos - photos.loc[pp.START_LINE_YEAR]
    cgm.SETUP = "calendar"
    coastsat = cgm.observed_change(pp.PERIOD)                # change since the calendar-1996 mean
    return photos, coastsat


def cmd_score(_a):
    photos, cs = observed()
    rows = []
    for label, spec in SPECS.items():
        g = model_gap(spec)
        row = dict(run=label, spec=spec)
        for name, y0, y1 in ERAS:
            if name.startswith("photos"):
                yrs = [y for y in photos.index if y0 <= y <= y1]
                d = g.loc[yrs].to_numpy() - photos.loc[yrs].to_numpy()
            else:                                            # CoastSat as change since 1996, the model likewise
                yrs = [y for y in cs.index if y0 <= y <= y1]
                d = (g.loc[yrs] - g.loc[1996]).to_numpy() - cs.loc[yrs].to_numpy()
            row[f"rmse {name}"] = float(np.sqrt(np.mean(d ** 2)))
        yrs = [y for y in photos.index if y <= LAST]
        row["rmse photos, all dates"] = float(np.sqrt(np.mean((g.loc[yrs].to_numpy() - photos.loc[yrs].to_numpy()) ** 2)))
        for y in (1978, 1995, 2009, 2016, 2025):
            row[f"gap {y}"] = float(g.loc[y])
        rows.append(row)
    t = pd.DataFrame(rows)
    t.round(2).to_csv(workspace_twin(HERE) / "continuous_scores.csv", index=False)
    pd.set_option("display.width", 250)
    print(t.drop(columns="spec").round(1).to_string(index=False))
    print("observed photos:", " ".join(f"{y}:{v:.0f}" for y, v in photos.items() if y <= LAST))


# The observed gap level since 1967 by year: the photo dates, and annual CoastSat placed on the photos by its
# mean offset over the years both have (2004 left out: the photos sit 42 m above CoastSat there)
def observed_levels():
    photos, cs = observed()
    shared = [y for y in photos.index if y in cs.index and y != 2004]
    placed = cs + float(np.mean(photos.loc[shared].to_numpy() - cs.loc[shared].to_numpy()))
    return photos, placed


# RMS error of a run's gap level against every photo date and CoastSat year in y0..y1
def level_rmse(spec, y0, y1, name=None, variant=""):
    photos, placed = observed_levels()
    g = model_gap(spec, name, variant)
    obs = pd.concat([photos[(photos.index >= y0) & (photos.index <= y1)],
                     placed[(placed.index >= y0) & (placed.index <= y1)]])
    return float(np.sqrt(np.mean((g.loc[obs.index].to_numpy() - obs.to_numpy()) ** 2)))


def history(s1, s2=None, s3=None):
    s2 = s1 if s2 is None else s2
    s3 = s2 if s3 is None else s3
    return f"1969:{s1:g},1996:{s2:g},2010:{s3:g},2017:0"


def run_specs(specs, streams=4):
    todo = [sp for sp in dict.fromkeys(specs) if shoreline(sp) is None]
    jobs, failed = queue.Queue(), []
    for sp in todo:
        jobs.put(sp)

    def worker():
        while True:
            try:
                sp = jobs.get_nowait()
            except queue.Empty:
                return
            if run(sp, sp) or shoreline(sp) is None:
                failed.append(sp)

    threads = [threading.Thread(target=worker) for _ in range(streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if failed:
        raise SystemExit(f"failed: {failed}")


# One era: run each candidate, score its level over the era, extend the grid once if the best is on an open edge
def fit_era(name, make_spec, grid, y0, y1, ceiling=None, log=None):
    grid = sorted({round(v, 3) for v in grid if v >= 0 and (ceiling is None or v <= ceiling + 1e-9)})
    for attempt in range(2):
        run_specs([make_spec(v) for v in grid])
        scores = {v: level_rmse(make_spec(v), y0, y1) for v in grid}
        best = min(scores, key=scores.get)
        step = grid[1] - grid[0] if len(grid) > 1 else 0.05
        edge = None
        if best == grid[-1] and (ceiling is None or best + step <= ceiling + 1e-9):
            edge = best + step
        elif best == grid[0] and best - step >= 0:
            edge = best - step
        if edge is None or attempt:
            break
        grid = sorted(set(grid) | {round(edge, 3)})
    for v in sorted(scores):
        print(f"  {name}: {v:g} -> {scores[v]:.1f} m{'  <- best' if v == best else ''}", flush=True)
        log.append(dict(era=name, strength=v, spec=make_spec(v), rmse_level_m=round(scores[v], 2), best=v == best))
    return best


def cmd_fit(a):
    log = []
    print("1) maintained, 1969-1995, scored on its settled years 1978-1995", flush=True)
    s1 = fit_era("1969-1995", lambda v: history(v), (0.25, 0.30, 0.35, 0.40), 1978, 1995, log=log)
    print("2) 1996-2009, holding 1969-1995 at", s1, flush=True)
    s2 = fit_era("1996-2009", lambda v: history(s1, v), [s1 - d for d in (0, 0.05, 0.10, 0.15)], 1996, 2009,
                 ceiling=s1, log=log)
    print("3) 2010-2016, holding", s1, "and", s2, flush=True)
    s3 = fit_era("2010-2016", lambda v: history(s1, s2, v), [s2 - d for d in (0, 0.05, 0.10, 0.15)], 2010, 2016,
                 ceiling=s2, log=log)
    best = history(s1, s2, s3)
    print(f"fitted history: {best}", flush=True)
    for name, y0, y1 in (("1969-1995", 1970, 1995), ("1996-2009", 1996, 2009), ("2010-2016", 2010, 2016),
                         ("2017-2025", 2017, 2025), ("1967-2025", 1968, 2025)):
        print(f"  level RMSE {name}: {level_rmse(best, y0, y1):.1f} m", flush=True)
    pd.DataFrame(log).to_csv(workspace_twin(HERE) / "continuous_fit.csv", index=False)
    (workspace_twin(HERE) / "continuous_fit_best.txt").write_text(best + "\n")


# Jobs (label, spec, name, variant), run four at a time
def run_jobs(jobs, streams=4):
    todo = [j for j in jobs if shoreline(*j[1:]) is None]
    print(f"{len(jobs) - len(todo)} of {len(jobs)} done; running {len(todo)} in {streams} streams", flush=True)
    q, failed = queue.Queue(), []
    for j in todo:
        q.put(j)

    def worker():
        while True:
            try:
                j = q.get_nowait()
            except queue.Empty:
                return
            if run(*j) or shoreline(*j[1:]) is None:
                failed.append(j[0])

    threads = [threading.Thread(target=worker) for _ in range(streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    if failed:
        raise SystemExit(f"failed: {failed}")


# A linear ramp from 0 to s over 1969..end, then the fitted later eras
def ramp(s, end, s2, s3):
    n = end - 1969 + 1
    steps = ",".join(f"{y}:{s * (y - 1968) / n:.4f}" for y in range(1969, end + 1))
    return f"{steps},1996:{s2:g},2010:{s3:g},2017:0"


# After the fit: (1) did the 1995 repair strengthen the groin? (2) does the maintained strength hold under other
# pre-1984 storms? (3) does a ramp-in fix the too-fast build without moving the plateau?
def cmd_followups(a):
    s1, s2, s3 = (v for _, v in cgm.parse_periods((workspace_twin(HERE) / "continuous_fit_best.txt").read_text().strip())[:3])
    fitted = history(s1, s2, s3)
    for seed in (1, 2):
        if not storms_file(f"seed{seed}").is_file():
            build_resampled(seed)
    repair = [(f"1996-2009 at {v:g}", history(s1, v, s3), None, "") for v in (0.30, 0.35)]
    storms = [(f"maintained {v:g}, seed{seed}", history(v), None, f"seed{seed}")
              for seed in (1, 2) for v in (0.20, 0.25, 0.30)]
    ramps = [(f"ramp to {s1:g} by {end}", ramp(s1, end, s2, s3), f"ramp{end}", "") for end in (1974, 1978)]
    run_jobs(repair + storms + ramps, a.streams)
    rows = []
    for label, spec, name, variant in [("fitted", fitted, None, "")] + repair:
        rows.append(dict(test="repair", run=label, **{f"rmse {e}": level_rmse(spec, *yy, name, variant) for e, yy in
                                                      (("1996-2009", (1996, 2009)), ("2010-2016", (2010, 2016)))}))
    for label, spec, name, variant in [(f"maintained {v:g}, replay", history(v), None, "") for v in (0.20, 0.25, 0.30)] + storms:
        rows.append(dict(test="storms", run=label, **{"rmse 1978-1995": level_rmse(spec, 1978, 1995, name, variant)}))
    for label, spec, name, variant in [("fitted, no ramp", fitted, None, "")] + ramps:
        rows.append(dict(test="ramp", run=label, **{f"rmse {e}": level_rmse(spec, *yy, name, variant) for e, yy in
                                                    (("1970-1977 build", (1970, 1977)), ("1978-1995", (1978, 1995)),
                                                     ("1970-1995", (1970, 1995)))}))
    t = pd.DataFrame(rows)
    t.round(2).to_csv(workspace_twin(HERE) / "continuous_followups.csv", index=False)
    pd.set_option("display.width", 250)
    for test, part in t.groupby("test", sort=False):
        print(f"\n== {test} ==")
        print(part.dropna(axis=1, how="all").drop(columns="test").round(1).to_string(index=False))


def cmd_figure(_a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import (C, INK, INK_MUTED, _letter_inside, apply_style, figsize, open_frame,
                                             record_caption, save)
    apply_style()
    photos, placed = observed_levels()
    s1, s2, s3 = (v for _, v in cgm.parse_periods((workspace_twin(HERE) / "continuous_fit_best.txt").read_text().strip())[:3])
    fitted = ramp(s1, 1978, s2, s3)
    runs = [  # label, spec, name, colour, line style, width
        ("no groin", SPECS["no groin"], None, C["BASE"], "-", 1.4),
        ("window fits run continuously (0.25, then 0.16, 0.10, off)", SPECS["stepped, 0.25 maintained"], None,
         C["ACCENT"], (0, (3, 2)), 1.4),
        (f"fitted history (ramp to {s1:g} by 1978, {s2:g} to 2009, {s3:g} to 2016, off)", fitted, "ramp1978",
         INK, "-", 2.2),
    ]
    fig, (ax, axs) = plt.subplots(2, 1, figsize=figsize("double", height=5.0), sharex=True, constrained_layout=True,
                                  gridspec_kw=dict(height_ratios=(3, 1.1)))
    for label, spec, name, col, ls, lw in runs:
        g = model_gap(spec, name)
        ax.plot(g.index, g.values, color=col, ls=ls, lw=lw, label=f"model, {label}")
        if spec != SPECS["no groin"]:
            yrs = np.arange(1966, 2027)
            axs.step(yrs, [cgm.strength_for(spec, y) if y >= 1969 else 0 for y in yrs], where="post", color=col,
                     ls=ls, lw=lw)
    ax.plot(placed.index, placed.values, "o-", color=C["LOCATOR"], lw=0.8, ms=3, mfc="white",
            label="observed, annual CoastSat (placed on the photos)")
    ax.plot(photos.index, photos.values, "s", color=C["LOCATOR"], ms=5, label="observed, wet/dry photo dates")
    for a_ in (ax, axs):
        for x in (1969, 1995, 2017, 2022):
            a_.axvline(x, color=INK_MUTED, lw=0.6, ls=(0, (1, 2)))
        a_.axvspan(1966, 1984, color=C["BASE_FILL"], alpha=0.35, lw=0)
        a_.grid(axis="y")
        open_frame(a_)
    for x, nm in ((1969, "built"), (1995, "last repair"), (2017, "2017 fill;\ngroin off"), (2022, "2022 fill")):
        ax.text(x + 0.4, 178, nm, fontsize=6.3, color=INK_MUTED, va="top")
    ax.text(1975, -40, "storms replayed from 1984-1999", ha="center", fontsize=6.3, color=INK_MUTED)
    ax.axhline(0, color=INK_MUTED, lw=0.5)
    ax.set_xlim(1966, 2026)
    ax.set_ylim(-45, 180)
    ax.set_ylabel("Gap across the groins\nsince 1967 (m)")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center", ncol=2, fontsize=6.3, frameon=False)
    axs.set_ylim(0, 0.32)
    axs.set_ylabel("Strength")
    axs.set_xlabel("Year")
    _letter_inside(ax, 0)
    _letter_inside(axs, 1)
    fig.suptitle("The groin's life in one continuous run, 1967-2025")
    png = workspace_twin(HERE) / "figures" / "continuous_gap_1967_2025.png"
    save(fig, png, dpi=250, close=True)
    e = {k: level_rmse(fitted, *yy, "ramp1978") for k, yy in
         (("all", (1968, 2025)), ("build", (1970, 1977)), ("1978-1995", (1978, 1995)), ("1996-2009", (1996, 2009)),
          ("2010-2016", (2010, 2016)), ("2017-2025", (2017, 2025)))}
    w = {k: level_rmse(SPECS["stepped, 0.25 maintained"], *yy) for k, yy in
         (("all", (1968, 2025)), ("1996-2009", (1996, 2009)))}
    record_caption(png, (
        "**Tests:** whether one continuous history of the corrected groin's strength reproduces the groin's whole "
        "life, 1967-2025, with no reset to the observed shoreline. **How to read:** (a) the change since 1967 in the "
        "gap across the Buxton groins (GIS 6 seaward of GIS 5): one model run per line from the 1967 wet/dry "
        "shoreline (GIS 2-12, joined to the 1996 offset elsewhere); teal squares, the wet/dry photo dates; teal "
        "circles, annual CoastSat, which measures a different line, placed on the photos by its mean offset over "
        "the years both have (2004 left out). (b) The strength each groin run uses, as the share of the arriving "
        "drift held back. Shading: storms before 1984 replay 1984-1999. Model: the 1996 topography, RSLR 0.004, "
        "fills on (2014, 2017, 2022). **Shows:** the strengths fitted window by window (purple), run as one history, "
        f"drop the gap by about 25 m when the strength falls in 1996 and stay 30-40 m low ({w['1996-2009']:.0f} m "
        f"error for 1996-2009): those fits describe how a reset fillet relaxes, not how strong the groin was. Fitted "
        f"era by era on the continuous run (black), the groin ramps to {s1:g} by 1978, holds {s2:g} to 2009 and "
        f"{s3:g} to 2016, and is off from the 2017 fill; it follows the record to {e['all']:.0f} m over 1967-2025 "
        f"({e['build']:.0f} m in the 1970-1977 build-up, {e['1978-1995']:.0f} m for 1978-1995, {e['1996-2009']:.0f} m "
        f"for 1996-2009, {e['2010-2016']:.0f} m for 2010-2016, {e['2017-2025']:.0f} m after the fill). The maintained "
        f"{s1:g} is the minimum under three different pre-1984 storm sequences (continuous_followups.csv)."))
    print("wrote", png.relative_to(REPO), {k: round(v, 1) for k, v in e.items()})


def cmd_gif(_a):
    import groin_panels_gif as GP
    from site_layer.hat_figure_style import apply_style
    apply_style()
    gis = GP.GIS
    idx = gis + pp.PAD - 1
    model = {lab: np.load(shoreline(sp))[:, idx] for lab, sp in SPECS.items()}
    x0_raw = model["no groin"][0]
    target_raw = x0_raw + pp.shape_change(2023, gis)
    shown = ["no groin", "stepped, 0.25 maintained"]
    *runs_m, x0, target = GP.rotate([model[lab] for lab in shown] + [x0_raw, target_raw], x0_raw)
    years = list(range(FIRST - 1, LAST + 1))
    runs = [("No groin, fills on", runs_m[0], True, None),
            ("Groin 0.25 to 1995, 0.16 to 2009, 0.10 to 2016, off from 2017; fills on", runs_m[1], True,
             (lambda k: 1969 <= years[k] < 2017))]

    def status(k):
        y = years[k]
        s = cgm.strength_for(SPECS["stepped, 0.25 maintained"], y)
        return (("before the groins" if y < 1969 else f"groin strength {s:g}" if s else "groin off")
                + ("; storms replayed from 1984-1999" if y < 1984 else ""))

    GP.render(runs, x0, target, "1967 shoreline, before the groins (start)",
              "2023 observed shoreline shape (overall retreat removed)", [str(y) for y in years],
              "The groin's whole life in one run, 1967-2025, from the 1967 shoreline at GIS 2-12", status,
              workspace_twin(HERE) / "figures" / "continuous_life_panels.gif", fills_from=lambda k: years[k])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("storms", "score", "figure", "gif", "fit"):
        sub.add_parser(name)
    for name in ("run", "followups"):
        r = sub.add_parser(name)
        r.add_argument("--streams", type=int, default=4)
    a = ap.parse_args()
    {"storms": cmd_storms, "run": cmd_run, "score": cmd_score, "figure": cmd_figure, "gif": cmd_gif,
     "fit": cmd_fit, "followups": cmd_followups}[a.cmd](a)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch()
    else:
        main()
