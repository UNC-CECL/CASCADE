"""
Dune toes from the lidar, and a collision-regime dune erosion rule calibrated on the 1996-2009 surveys.

    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py toes
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py twl
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py calibrate
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py empirical
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py run
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py score
    python scripts/real-world-cases/Hatteras/dune-storm-response/HAT_dune_collision.py plot

Barrier3D lowers a dune only when Rhigh overtops it. Sallenger (2000) puts a storm in the collision regime
when the water reaches the dune toe but not the crest (D_low < R_high < D_high), where the real dunes lose
sand to scarping. The magnitude here is the wave-impact rule of Larson et al. (2004),
dV = 4 Cs (R - z_toe)^2 t / T, summed hour by hour over the hours the total water level stands above the
toe, and taken as a crest lowering over the model's two-cell (20 m) dune: dH = dV / 20 m.

toes       per survey, domain and 10 m cell, on the 1996-start dune windows: the crest as
           HAT_survey_dune_crests.py finds it, the beach start (first cell above 0.5 m MHW) and the toe,
           the point between them farthest below the straight line from the beach start to the crest
           (the maximum-distance-below-line definition). Also the start's own arrays, the toe the model
           would use. Writes dune_toes.csv
twl        the hourly total water level the storm series is built from (Duck gauge + Stockdon R2% on WIS
           ST63228, slope 0.06, as historical_storm_creation_v3_HAT.py), checked against the series' peaks.
           Writes twl_hourly_1996_2009.csv
calibrate  for each pair of consecutive surveys, the cells the water never overtopped between them: the
           observed crest change against the exposure sum((TWL - toe)^2 dt / Tp) from the cell's toe (the
           median over the surveys that pass the toe screen), and the Cs that fits. Writes
           collision_calibration_cells.csv, collision_calibration.csv and collision_calibration.png
empirical  the same cells in the terms of Long et al. (2014): crest lowering over pre-storm crest (dz/D_high)
           against R_high/D_high and against beach width, and the fit dz/D_high = a max(0, 1 - W / W0) on the
           storm intervals. Writes collision_empirical.csv, collision_empirical_fit.csv and collision_empirical.png
run        the managed edgeBE 1996-2009 hindcast with the groin, as the headline matrix run, with collision
           erosion added after each model year, in-process only: every dune cell loses
           Cs * 4 sum((TWL - toe)^2 dt / Tp) / 20 m over that calendar year's hours with the water between its
           toe and its crest, down to no lower than the toe; both dune rows scaled together, as the forced
           crests are. Optionally the eroded sand goes onto the beach, moving the shoreline seaward as
           CASCADE's shoreface nourishment does (dx_s = 2 V / (2 h_b + D_sf)), passed on to BRIE:
               larson_lost     Cs 1.7e-4 (Dennis fit, = Larson's laboratory value), sand lost
               larson_beach    Cs 1.7e-4, sand to the beach
               strong_lost     Cs 4e-4 (what the Dennis hotspots, GIS 11 and 86, need), sand lost
               strong_beach    Cs 4e-4, sand to the beach
           or by the beach-width rule from empirical: once in each calendar year whose highest water reaches
           the cell's toe but not its crest, the crest loses a max(0, 1 - W / W0) of its height (m NAVD88),
           W the cell's beach width in the last survey flown before that year's highest water:
               width_lost          sand lost
               width_beach         sand to the beach
               width_storms_beach  sand to the beach, only in years whose highest water reaches 3.1 m NAVD88,
                                   the lowest R_high of the two storm intervals the rule is fit on
           A threshold rule (no loss below an exposure, steeper above) fit Dennis no better than the linear
           one (r 0.45 against 0.50 on domain medians; best threshold 0), so it is not run. The control is
           the forced-dune-crests experiment's control run (same code, nothing added)
score      island skill, per-GIS LRR against CoastSat, and the model's crests against the lidar's in each
           survey's model state. Writes collision_runs_summary.csv, collision_runs_by_gis.csv and
           collision_runs_crests.csv
plot       every run's end-minus-start shoreline position change against the CoastSat net change, in the
           layout of the run's shoreline_position_change figure, for the island, Buxton and Rodanthe; the
           runs named after it (and the control), or every run.
           Writes collision_runs_position_<view>.png

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-26
"""
from __future__ import annotations

import argparse
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
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-26-dune-collision"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
START = "1996-start"
SURVEY_CRESTS = (PROJECT_ROOT / "data" / "hatteras_init" / "1-barrier3d-domains" / START
                 / "survey-crests" / "survey_dune_crests.csv")
STORM_SUMMARY = (PROJECT_ROOT / "data" / "hatteras_init" / "3-env-forcings" / "3-storms" / "hindcast_storms"
                 / "1997_2009" / "1997_2009_storms_v3_240_summary.csv")
MHW_NAVD88_M = 0.36
TWL_START, TWL_END = "1996-01-01 00:00:00", "2009-12-31 23:00:00"
BEACH_SLOPE = 0.06               # as the storm series
DUNE_WIDTH_M = 20.0              # Barrier3D's two 10 m dune cells
MIN_BEACH_CELLS = 2              # toe needs a beach start at least this many cells seaward of the crest
# The calibration's toe: per cell, the median over the surveys of the toes that pass the screen. Not each
# survey's own: the 1998 grid puts 38% of its toes on the wet beach (median 1.7 m, some below 0 m NAVD88),
# and it is the survey before Dennis
TOE_SURVEYS_EXCLUDED = (1998,)
TOE_MIN_M_NAVD88 = 1.0           # above the beach start (0.86 m NAVD88)
TOE_BELOW_CREST_M = 0.5          # and at least this far below the crest
# Fit on the storm intervals without the post-Isabel repairs (2001-2004) inside them; all intervals as well
FIT_INTERVALS = ("1996-1997", "1997-1998", "1998-1999", "1999-2000", "2000-2001", "2004-2005", "2005-2009")
# Empirical scaling: fit on the intervals with a storm signal in the crests (Dennis; 2004-05), each less its
# median change on beaches wider than this
EMPIRICAL_INTERVALS = ("1998-1999", "2004-2005")
BACKGROUND_BEACH_WIDTH_M = 60.0
# Model runs
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
CONTROL_DIR = (PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes" / "2026-10-26-forced-dune-crests"
               / "runs" / "control")
RUN_START, RUN_END, PRESET, SCENARIO = 1996, 2009, "edgeBE", "full_management"
# The rule, its Cs (Larson) or the lowest year-maximum TWL it acts in (beach width, m NAVD88; None: any year
# the water reaches the toe), and whether the eroded sand goes onto the beach
VARIANTS = {"larson_lost": ("larson", 1.7e-4, False),
            "larson_beach": ("larson", 1.7e-4, True),
            "strong_lost": ("larson", 4e-4, False),
            "strong_beach": ("larson", 4e-4, True),
            "width_lost": ("width", None, False),
            "width_beach": ("width", None, True),
            "width_storms_beach": ("width", 3.1, True)}
WORKERS = 4
REPORT_GIS = range(5, 21)
# -----------------------------------------------------------------------------


# Toes

# Beach start, crest and toe per profile (m MHW, ocean first): the crest as HAT_survey_dune_crests.crests finds it
def crest_and_toe(z, i0, i1, X):
    n = z.shape[0]
    out = {k: np.full(n, np.nan) for k in ("crest", "toe", "beach_start_idx", "crest_idx", "toe_idx")}
    above = z > X.BEACH_START_THR_M
    sb = np.where(above.any(axis=1), above.argmax(axis=1), -1)
    for i in range(n):
        if sb[i] < 0:
            continue
        a = max(i0, int(sb[i])) if X.CLIP_WINDOW_TO_BEACH else i0
        b = min(i1, z.shape[1])
        if b <= a:
            continue
        win = z[i, a:b]
        if (win <= X.NODATA_SENTINEL_M + 1e-9).any():
            continue
        valid = win > X.SENTINEL_WATER_M + 1e-9
        if not valid.any():
            continue
        c = a + int(np.argmax(np.where(valid, win, -np.inf)))
        out["crest"][i], out["crest_idx"][i], out["beach_start_idx"][i] = z[i, c], c, sb[i]
        s = int(sb[i])
        if c - s < MIN_BEACH_CELLS:
            continue
        seg = z[i, s:c + 1]
        if (seg <= X.NODATA_SENTINEL_M + 1e-9).any():
            continue
        line = seg[0] + (seg[-1] - seg[0]) * np.arange(len(seg)) / (len(seg) - 1)
        k = int(np.argmax(line - seg))
        if k in (0, len(seg) - 1):
            k = 0 if seg[0] >= seg[-1] else k      # a straight or convex face: the toe is the beach start
        out["toe"][i], out["toe_idx"][i] = seg[k], s + k
    return out


def toes():
    sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
    import HAT_survey_dune_crests as S
    import json
    import rasterio
    X, tv = S.X, S.tv
    prod, surveys = S.STARTS[START]
    windows = json.loads(tv.picks_dir(START).joinpath(
        f"HAT_dune_search_windows_{tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
    frames = {g: S.frame(START, g, windows) for g in range(1, 91)}
    clips = S.product(prod).gapfill_1m
    grids = {g: S.clip_grid(clips / f"clip_domain_{g}_filled.tif") for g in range(1, 91)}
    offsets = S.SURVEY_VERTICAL_OFFSET_M.get(START, {})
    rows = []

    def add(survey, g, z, i0, i1):
        r = crest_and_toe(z, i0, i1, X)
        for cell in range(50):
            rows.append(dict(survey=survey, gis=g, cell=cell,
                             crest_m_navd88=r["crest"][cell] + MHW_NAVD88_M,
                             toe_m_navd88=r["toe"][cell] + MHW_NAVD88_M,
                             beach_width_m=(r["crest_idx"][cell] - r["beach_start_idx"][cell]) * 10,
                             toe_to_crest_m=(r["crest_idx"][cell] - r["toe_idx"][cell]) * 10))

    # The start's own arrays: the dunes the model starts from
    for g in range(1, 91):
        shear, c0, i0, i1 = frames[g]
        z = S.profiles(np.load(tv.npy_dirs(START)[0] / f"domain_{g}.npy").astype(float), shear, c0)
        add("start", g, z, i0, i1)
    print("start arrays done", flush=True)
    for year, (kind, path, flown, state, note) in surveys.items():
        files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
        srcs = [rasterio.open(f) for f in files]
        try:
            unit = "m" if kind == "mosaic" else S.unit_of(srcs, START, grids)[0]
            scale = S.US_FT if unit == "ft" else 1.0
            for g in range(1, 91):
                shear, c0, i0, i1 = frames[g]
                on_clip = (S.mosaic_on_clip(srcs[0], grids[g]) if kind == "mosaic"
                           else S.raster_on_clip(srcs, grids[g]))
                z = S.profiles(S.downsample(on_clip * scale + offsets.get(year, 0.0)), shear, c0)
                add(year, g, z, i0, i1)
        finally:
            for s in srcs:
                s.close()
        print(f"{year} done", flush=True)
    t = pd.DataFrame(rows)
    # The crest script's spike screen: keep its status beside each survey cell
    sc = pd.read_csv(SURVEY_CRESTS)[["survey", "gis", "cell", "status", "crest_m_navd88"]]
    sc["survey"] = sc.survey.astype(str)
    t["survey"] = t.survey.astype(str)
    t = t.merge(sc.rename(columns={"crest_m_navd88": "crest_script_m_navd88"}), on=["survey", "gis", "cell"],
                how="left")
    t.loc[t.survey == "start", "status"] = "ok"
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    t.round(3).to_csv(EXP_DIR / "dune_toes.csv", index=False)
    same = (t.crest_m_navd88 - t.crest_script_m_navd88).abs()
    print(f"crest check against survey_dune_crests.csv: max |difference| {same.max():.4f} m")
    ok = t[(t.status == "ok") & t.toe_m_navd88.notna()]
    print(ok.groupby("survey")[["toe_m_navd88", "crest_m_navd88", "beach_width_m"]].median().round(2))
    print(f"wrote {EXP_DIR / 'dune_toes.csv'}")


# Hourly total water level

def twl():
    from site_layer import hat_env_forcings as env
    w = pd.read_csv(env.DUCK_GAUGE_FILE, index_col="t", parse_dates=True)["v"]
    wis = pd.read_csv(env.WIS_FILE, index_col="time", parse_dates=True)[["waveHs", "waveTp"]]
    idx = pd.date_range(TWL_START, TWL_END, freq="h")
    df = pd.DataFrame({"water_level": w.reindex(idx), "Hs": wis.waveHs.reindex(idx),
                       "Tp": wis.waveTp.reindex(idx)}).dropna()
    # Stockdon et al. (2006) R2%, as historical_storm_creation_v3_HAT.calculate_r2_percent
    L0 = 9.81 * df.Tp ** 2 / (2 * np.pi)
    setup = 0.35 * BEACH_SLOPE * np.sqrt(df.Hs * L0)
    swash = np.sqrt((0.75 * BEACH_SLOPE * np.sqrt(df.Hs * L0)) ** 2 + (0.06 * np.sqrt(df.Hs * L0)) ** 2)
    df["R2"] = 1.1 * (setup + swash / 2)
    df["TWL"] = df.water_level + df.R2
    # Check: each storm's peak in the series is the peak of these hours
    ev = pd.read_csv(STORM_SUMMARY, parse_dates=["StartTime", "EndTime"])
    peak = np.array([df.TWL.loc[a:b].max() for a, b in zip(ev.StartTime, ev.EndTime)])
    diff = peak - (ev.Rhigh * 10 + MHW_NAVD88_M)
    print(f"{len(df)} hours; storm peaks against the series: max |difference| {np.nanmax(np.abs(diff)):.4f} m "
          f"over {len(ev)} events")
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    df.round(4).to_csv(EXP_DIR / "twl_hourly_1996_2009.csv", index_label="time")
    print(f"wrote {EXP_DIR / 'twl_hourly_1996_2009.csv'}")


# Calibration

# The survey rows that passed the crest script's screen
def survey_table():
    t = pd.read_csv(EXP_DIR / "dune_toes.csv")
    t = t[(t.survey != "start") & (t.status == "ok")].copy()
    t["survey"] = t.survey.astype(int)
    return t


# Per cell, the median toe over the surveys whose toe passes the screen (m NAVD88)
def cell_toes(t):
    good = t[~t.survey.isin(TOE_SURVEYS_EXCLUDED) & (t.toe_m_navd88 > TOE_MIN_M_NAVD88)
             & (t.toe_m_navd88 < t.crest_m_navd88 - TOE_BELOW_CREST_M)]
    return good.groupby(["gis", "cell"]).toe_m_navd88.median().rename("toe_cell_m_navd88")


# First and last day flown, from "YYYY-MM" or "YYYY-MM to YYYY-MM"
def flown_span(flown):
    parts = flown.split(" to ")
    first = pd.Timestamp(parts[0] + "-01")
    last = pd.Timestamp(parts[-1] + "-01") + pd.offsets.MonthEnd(0) + pd.Timedelta(hours=23)
    return first, last


def calibrate():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()

    t = survey_table()
    toe_cell = cell_toes(t)
    print(f"per-cell toe: {len(toe_cell)} cells, median {toe_cell.median():.2f}, "
          f"10-90% {toe_cell.quantile(0.1):.2f}-{toe_cell.quantile(0.9):.2f} m NAVD88")
    flown = pd.read_csv(SURVEY_CRESTS).groupby("survey").flown.first()
    h = pd.read_csv(EXP_DIR / "twl_hourly_1996_2009.csv", index_col="time", parse_dates=True)
    years = sorted(t.survey.unique())
    cells, summary = [], []
    for a, b in zip(years[:-1], years[1:]):
        t0, t1 = flown_span(flown[a])[1], flown_span(flown[b])[0]
        hh = h.loc[t0:t1]
        pa = t[t.survey == a].set_index(["gis", "cell"])
        pb = t[t.survey == b].set_index(["gis", "cell"])
        j = pa.join(pb[["crest_m_navd88"]], rsuffix="_b", how="inner").join(toe_cell, how="inner").dropna(
            subset=["toe_cell_m_navd88", "crest_m_navd88", "crest_m_navd88_b"])
        twl_max = hh.TWL.max()
        # Not overtopped: the water never reached the earlier crest
        j = j[j.crest_m_navd88 > twl_max].copy()
        # Exposure, m^2: sum over hours of (TWL - toe)^2 * 3600 s / Tp
        toe = j.toe_cell_m_navd88.to_numpy()[:, None]
        excess = np.clip(hh.TWL.to_numpy()[None, :] - toe, 0, None)
        j["exposure_m2"] = (excess ** 2 * 3600.0 / hh.Tp.to_numpy()[None, :]).sum(axis=1)
        j["hours_above_toe"] = (excess > 0).sum(axis=1)
        j["crest_loss_m"] = j.crest_m_navd88 - j.crest_m_navd88_b
        j["interval"] = f"{a}-{b}"
        j["years"] = (t1 - t0).days / 365.25
        cells.append(j.reset_index())
        summary.append(dict(interval=f"{a}-{b}", years=round(j.years.iloc[0], 2), twl_max_m_navd88=round(twl_max, 2),
                            cells=len(j), toe_median=j.toe_cell_m_navd88.median(),
                            hours_above_toe_median=j.hours_above_toe.median(),
                            exposure_median_m2=j.exposure_m2.median(),
                            crest_loss_median_m=j.crest_loss_m.median()))
    c = pd.concat(cells, ignore_index=True)
    s = pd.DataFrame(summary).set_index("interval")

    # Fit dH = (4 Cs / W) E - g * years: Cs the collision coefficient, g the net growth between surveys
    A = np.column_stack([4 * c.exposure_m2 / DUNE_WIDTH_M, -c.years])
    cs_all, g_all = np.linalg.lstsq(A, c.crest_loss_m, rcond=None)[0]
    # The same on domain medians per interval, which damps the lidar's cell-to-cell noise
    d = c.groupby(["interval", "gis"]).agg(exposure_m2=("exposure_m2", "median"), years=("years", "first"),
                                           crest_loss_m=("crest_loss_m", "median")).reset_index()
    A = np.column_stack([4 * d.exposure_m2 / DUNE_WIDTH_M, -d.years])
    cs_every, _ = np.linalg.lstsq(A, d.crest_loss_m, rcond=None)[0]
    fit = d.interval.isin(FIT_INTERVALS).to_numpy()
    cs_dom, g_dom = np.linalg.lstsq(A[fit], d.crest_loss_m[fit], rcond=None)[0]
    # Per interval, the growth term held at the island fit
    s["cs_interval"] = [
        float(np.sum(4 * x.exposure_m2 / DUNE_WIDTH_M * (x.crest_loss_m + g_dom * x.years))
              / np.sum((4 * x.exposure_m2 / DUNE_WIDTH_M) ** 2)) if (x.exposure_m2 > 0).any() else np.nan
        for _, x in d.groupby("interval", sort=False)]
    c["predicted_loss_m"] = 4 * cs_dom * c.exposure_m2 / DUNE_WIDTH_M - g_dom * c.years
    d["predicted_loss_m"] = 4 * cs_dom * d.exposure_m2 / DUNE_WIDTH_M - g_dom * d.years
    r_cells = np.corrcoef(c.predicted_loss_m, c.crest_loss_m)[0, 1]
    r_dom = np.corrcoef(d.predicted_loss_m, d.crest_loss_m)[0, 1]

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(s.round(3))
    print(f"\ncells:   Cs {cs_all:.2e}  growth {g_all:+.3f} m/yr   (n={len(c)})")
    print(f"domains: Cs {cs_dom:.2e}  growth {g_dom:+.3f} m/yr   r {r_dom:.2f} (cells r {r_cells:.2f}, n={len(d)}); "
          f"every interval, 2001-2004 too: Cs {cs_every:.2e}")
    print(f"Larson et al. (2004) laboratory Cs: 1.7e-04")
    c.round(4).to_csv(EXP_DIR / "collision_calibration_cells.csv", index=False)
    s.assign(cs_domains=cs_dom, growth_domains_m_yr=g_dom, r_domains=r_dom).round(6).to_csv(
        EXP_DIR / "collision_calibration.csv")

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axs[0]
    for name, x in d.groupby("interval", sort=False):
        ax.plot(4 * x.exposure_m2 / DUNE_WIDTH_M, x.crest_loss_m, "o", ms=3.5, alpha=0.8, label=name)
    xs = np.linspace(0, (4 * d.exposure_m2 / DUNE_WIDTH_M).max(), 50)
    ax.plot(xs, cs_dom * xs, color="k", lw=1.2, label=f"fit, Cs = {cs_dom:.1e}")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("4 Σ (TWL − toe)² Δt / T / 20 m   (m)")
    ax.set_ylabel("Crest lowering between surveys (m)")
    ax.set_title("Domain medians, cells never overtopped", loc="left")
    ax.legend(fontsize=7, ncol=2)
    ax = axs[1]
    toe0 = t[t.survey == years[0]].groupby("gis").toe_m_navd88.median()
    crest0 = t[t.survey == years[0]].groupby("gis").crest_m_navd88.median()
    ax.plot(toe0.index, crest0, color="0.5", lw=1.2, label=f"crest, {years[0]} survey")
    ax.plot(toe0.index, toe0, color="#2166AC", lw=1.6, label=f"toe, {years[0]} survey")
    st = pd.read_csv(EXP_DIR / "dune_toes.csv")
    st = st[st.survey == "start"].groupby("gis").toe_m_navd88.median()
    ax.plot(st.index, st, color="#2166AC", lw=1.0, ls="--", label="toe, model start arrays")
    ax.axhline(1.7, color="0.3", lw=0.8, ls=":", label="Barrier3D berm, 1.7 m")
    ax.set_xlabel("GIS domain")
    ax.set_ylabel("m NAVD88")
    ax.set_title("Dune toe and crest by domain (median)", loc="left")
    ax.legend(fontsize=7)
    fig.tight_layout()
    out = EXP_DIR / "collision_calibration.png"
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Left: per domain and survey interval, the median crest lowering of the 10 m cells the total water "
        "level never overtopped between the two surveys, against the Larson et al. (2004) wave-impact exposure "
        "4 Σ (TWL − z_toe)² Δt / T over the hours above the earlier survey's toe, per 20 m of dune width; line: "
        f"the least-squares Cs with a net growth term ({g_dom:+.3f} m/yr) held for every interval. "
        "Right: domain-median toe and crest, the toe as the point farthest below the line from the beach start "
        f"(0.5 m MHW) to the crest. Experiment {TAG}."))
    print(f"wrote {out}")


# Empirical scaling, in the terms of Long et al. (2014): dz / D_high against R_high / D_high, and the beach width

def empirical():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    from site_layer.hat_figure_style import apply_style, record_caption
    apply_style()

    c = pd.read_csv(EXP_DIR / "collision_calibration_cells.csv")
    r_high = pd.read_csv(EXP_DIR / "collision_calibration.csv").set_index("interval").twl_max_m_navd88
    c = c[c.interval != "2001-2004"].copy()
    c["r_high"] = c.interval.map(r_high)
    c["x"] = c.r_high / c.crest_m_navd88
    c["y"] = c.crest_loss_m / c.crest_m_navd88
    # Less each interval's background (growth, repair, lidar offset): the median on the wide beaches
    c["y_net"] = c.y - c.interval.map(c[c.beach_width_m > BACKGROUND_BEACH_WIDTH_M].groupby("interval").y.median())
    st = c[c.interval.isin(EMPIRICAL_INTERVALS) & (c.y.abs() < 1.5)]

    def rule(w, a, w0):
        return a * np.clip(1 - w / w0, 0, None)

    (a, w0), _ = curve_fit(rule, st.beach_width_m, st.y_net, p0=[0.4, 40.0])
    r_cells = np.corrcoef(rule(st.beach_width_m, a, w0), st.y_net)[0, 1]
    dm = st.assign(p=rule(st.beach_width_m, a, w0)).groupby(["interval", "gis"])[["p", "y_net", "x"]].median()
    r_dom, r_dom_x = np.corrcoef(dm.p, dm.y_net)[0, 1], np.corrcoef(dm.x, dm.y_net)[0, 1]
    print(f"dz/D_high = {a:.3f} max(0, 1 - W / {w0:.1f} m)   r {r_cells:.2f} cells, {r_dom:.2f} domain medians "
          f"(R_high/D_high alone: {r_dom_x:.2f}); fit on {', '.join(EMPIRICAL_INTERVALS)}, n={len(st)}")

    widths = [0, 20, 30, 40, 60, 400]
    tab = c.groupby(["interval", pd.cut(c.beach_width_m, widths), pd.cut(c.x, [0, .6, .7, .8, .9, 1])],
                    observed=True).y.agg(["size", "median"]).reset_index()
    tab.columns = ["interval", "beach_width_m", "r_high_over_d_high", "cells", "dz_over_d_high_median"]
    tab.round(4).to_csv(EXP_DIR / "collision_empirical.csv", index=False)
    pd.DataFrame([dict(a=a, w0_m=w0, r_cells=r_cells, r_domains=r_dom, r_domains_rhigh_only=r_dom_x,
                       intervals=" ".join(EMPIRICAL_INTERVALS), cells=len(st))]).round(4).to_csv(
        EXP_DIR / "collision_empirical_fit.csv", index=False)

    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    shades = ["#7F2704", "#D94801", "#FD8D3C", "#FDBE85", "#BDBDBD"]
    xb = np.arange(0.45, 1.0, 0.1)
    ax = axs[0]
    for (lo, hi), col in zip(zip(widths[:-1], widths[1:]), shades):
        s = st[(st.beach_width_m > lo) & (st.beach_width_m <= hi)]
        m = s.groupby(pd.cut(s.x, np.arange(0.4, 1.01, 0.1)), observed=False).y
        ax.plot(xb, m.median(), "o-", color=col, lw=1.6, ms=5,
                label=f"beach {lo}–{hi} m" if hi < 400 else f"beach > {lo} m")
    ax.axhspan(-0.1, 0.1, color="0.9", zorder=0, label="Long et al. (2014): within 10%")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("R_high / D_high")
    ax.set_ylabel("Crest lowering / D_high")
    ax.set_title("Dennis and 2004–05, cells in collision", loc="left")
    ax.legend(fontsize=7)
    ax = axs[1]
    wb = np.arange(5, 100, 10)
    for iv, x in c.groupby("interval"):
        storm = iv in EMPIRICAL_INTERVALS
        m = x.groupby(pd.cut(x.beach_width_m, np.arange(0, 101, 10)), observed=False).y_net.median()
        ax.plot(wb, m, "o-" if storm else "-", color="#D94801" if storm else "0.7", lw=1.6 if storm else 1.0,
                ms=4, label=f"{iv} (R_high {r_high[iv]:.2f} m)", zorder=3 if storm else 2)
    ws = np.linspace(0, 100, 200)
    ax.plot(ws, rule(ws, a, w0), color="k", lw=1.4, label=f"{a:.2f} max(0, 1 − W / {w0:.0f} m)")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("Beach width, 0.5 m MHW contour to the crest (m)")
    ax.set_ylabel("Crest lowering / D_high, less background")
    ax.set_title("Every survey interval", loc="left")
    ax.legend(fontsize=6.5)
    fig.tight_layout()
    out = EXP_DIR / "collision_empirical.png"
    fig.savefig(out, dpi=200, facecolor="white")
    record_caption(out, (
        "Dune crest change across the 1996–2009 lidar surveys in the normalisation of Long et al. (2014), for the "
        "10 m cells the total water level (Duck gauge + Stockdon R2%, slope 0.06) never overtopped between two "
        "surveys. Left: median crest lowering over pre-storm crest elevation (NAVD88) against R_high / D_high, "
        "per beach-width class, for the two storm intervals (1998–99, Dennis; 2004–05); shaded, the collision-"
        "regime bound Long et al. report. Right: the same, less each interval's median on beaches wider than "
        f"{BACKGROUND_BEACH_WIDTH_M:g} m, against beach width (beach start at 0.5 m MHW to the crest); storm "
        f"intervals orange, the rest grey; line, the least-squares fit on the storm intervals (r {r_cells:.2f} "
        f"cells, {r_dom:.2f} domain medians). Experiment {TAG}."))
    print(f"wrote {out}")


# Model runs

# Crest lowering (m) by Larson et al. (2004) over the hours in the collision regime: above the toe, below the crest
def larson_loss(hours, cs):
    twl, tp = hours.TWL.to_numpy(), hours.Tp.to_numpy()

    def loss(g, crest, z_toe):
        reach = (twl[None, :] > z_toe[:, None]) & (twl[None, :] < crest[:, None])
        excess = np.where(reach, twl[None, :] - z_toe[:, None], 0.0)
        return cs * 4 * (excess ** 2 * 3600.0 / tp[None, :]).sum(axis=1) / DUNE_WIDTH_M
    return loss


# Crest lowering (m) by the beach-width rule, once, where the year's highest water is between the toe and the crest
def width_loss(hours, width, a, w0, min_twl):
    r_high = hours.TWL.max()

    def loss(g, crest, z_toe):
        if min_twl is not None and r_high < min_twl:
            return np.zeros_like(crest)
        w = width[g][: len(crest)]
        hit = (r_high > z_toe) & (r_high < crest)
        return np.where(hit, a * np.clip(1 - w / w0, 0, None) * crest, 0.0)
    return loss


# Collision erosion over one calendar year for every real domain, applied to the state just written, the
# lowering per cell from loss(gis, crest, toe); returns per GIS the mean crest lowering (m), the sand removed
# (m^3/m) and the shoreline move (m, + seaward)
def collide(cascade, dom, toe, loss, to_beach):
    out = {}
    for g in range(1, 91):
        b = cascade.barrier3d[dom.gis_to_pad(g)]
        t = b._time_index - 1
        rows = b._DuneDomain[t]                                       # dam above the berm, (cells, 2)
        top = rows.max(axis=1)
        crest = (top + b._BermEl) * 10 + MHW_NAVD88_M                 # m NAVD88
        z_toe = toe[g][: len(crest)]
        lowering = np.nan_to_num(loss(g, crest, z_toe))               # m
        # No lower than the toe; a crest already below its toe stays where it is
        new_crest = np.maximum(crest - lowering, np.minimum(z_toe, crest))
        new_top = np.maximum((new_crest - MHW_NAVD88_M) / 10 - b._BermEl, b._DuneRestart)
        safe_top = np.where(top > 0, top, 1.0)
        new_rows = np.where(top[:, None] > 0,
                            np.maximum(rows * (new_top / safe_top)[:, None], b._DuneRestart), rows)
        new_rows = np.minimum(new_rows, rows)                         # erosion only
        removed = rows - new_rows                                     # dam, per cell and row
        b._DuneDomain[t] = new_rows
        # m^3 per m of shore: each row is 10 m cross-shore, averaged over the domain's cells
        vol = float(removed.sum(axis=1).mean()) * 10 * 10
        dx = 0.0
        if to_beach and vol > 0:
            # CASCADE's shoreface_nourishment, in metres: dx_s = 2 V / (2 h_b + D_sf)
            dx = 2 * vol / (2 * b.h_b_TS[-1] * 10 + b._DShoreface * 10)
            b.x_s = b.x_s - dx / 10
            b.x_s_TS[-1] = b.x_s
        out[g] = (float(np.mean((top - new_rows.max(axis=1)) * 10)), vol, dx)
    return out


# Child process: collision erosion after each model year, then the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    rule, cs, to_beach = VARIANTS[variant]
    st = survey_table()

    # Per GIS, 50 cells from a (gis, cell) series, the domain median where a cell has none
    def per_gis(sr):
        dom_median = sr.groupby(level="gis").median()
        out = {}
        for g in range(1, 91):
            v = np.full(50, dom_median.get(g, np.nan))
            if g in sr.index.get_level_values("gis"):
                s = sr.loc[g]
                v[s.index.to_numpy()] = s.to_numpy()
            out[g] = v
        return out

    toe = per_gis(cell_toes(st))
    hours = pd.read_csv(EXP_DIR / "twl_hourly_1996_2009.csv", index_col="time", parse_dates=True)
    if rule == "width":
        fit = pd.read_csv(EXP_DIR / "collision_empirical_fit.csv").iloc[0]
        flown = pd.read_csv(SURVEY_CRESTS).groupby("survey").flown.first()
        widths = {y: per_gis(st[st.survey == y].set_index(["gis", "cell"]).beach_width_m)
                  for y in sorted(st.survey.unique())}

        # The beach widths of the last survey flown before the year's highest water (the first survey if none)
        def width_for(yh):
            peak = yh.TWL.idxmax()
            before = [y for y in widths if flown_span(flown[y])[1] < peak]
            return widths[before[-1] if before else min(widths)]

        def loss_for(yh):
            return width_loss(yh, width_for(yh), fit.a, fit.w0_m, cs)
    else:
        def loss_for(yh):
            return larson_loss(yh, cs)
    log = []
    orig_build = H.build_cascade

    def build_cascade(*a, **k):
        c = orig_build(*a, **k)
        cls = type(c)
        orig_update = cls.update
        calls = [0]

        # Patched on the class, so the saved model still pickles
        def update(self):
            orig_update(self)
            calls[0] += 1
            year = RUN_START + calls[0]                     # update k runs calendar year 1996 + k
            dom = HATTERAS_DOMAINS.for_total(len(self.barrier3d))
            res = collide(self, dom, toe, loss_for(hours.loc[str(year)]), to_beach)
            if to_beach and self._alongshore_transport_module:
                b3 = self.barrier3d
                n = len(b3)
                self._brie_coupler.update_brie_for_human_modifications(
                    [b3[i].x_t_TS[-1] for i in range(n)], [b3[i].x_s_TS[-1] for i in range(n)],
                    [b3[i].x_b_TS[-1] for i in range(n)], [b3[i].h_b_TS[-1] for i in range(n)],
                    [b3[i].s_sf_TS[-1] for i in range(n)])
            yr = pd.DataFrame([dict(year=year, gis=g, crest_lowering_m=dh, volume_m3_m=vol, shoreline_seaward_m=dx)
                               for g, (dh, vol, dx) in res.items()])
            log.append(yr)
            top = yr.loc[yr.crest_lowering_m.idxmax()]
            print(f"\nEXPERIMENT collision {year}: mean crest lowering {yr.crest_lowering_m.mean():.3f} m, "
                  f"max {top.crest_lowering_m:.2f} m (GIS {int(top.gis)}), sand {yr.volume_m3_m.mean():.2f} m3/m"
                  + (f", shoreline {yr.shoreline_seaward_m.mean():.2f} m seaward" if to_beach else ", lost"),
                  flush=True)

        cls.update = update
        return c

    H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    try:
        runpy.run_path(str(HINDCAST), run_name="__main__")
    finally:
        if log:
            pd.concat(log).round(4).to_csv(EXP_DIR / "logs" / f"collision_{variant}.csv", index=False)


def run_dir(variant):
    root = CONTROL_DIR if variant == "control" else EXP_DIR / "runs" / variant
    base = root / f"{RUN_START}_{RUN_END}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# One run in its own process, logged, with its launch record
def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(RUN_START),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO, "HAT_GROIN_ENABLED": "1",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    with open(logs / f"{variant}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rule, cs, to_beach = VARIANTS[variant]
    rec = dict(variant=variant, rule=rule, **({"cs": cs} if rule == "larson" else {"min_twl_m": cs}),
               sand_to_beach=to_beach,
               returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{variant:14s} exit {p.returncode}  {rec['minutes']} min", flush=True)


# Island skill, per-GIS LRR against CoastSat, and the modelled crests against each survey's
def score():
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.saved_runs import load_saved_cascade
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(_HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW

    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=RUN_START,
                               lowess_config=LOWESS_CONFIG, domains=DOM)
    target = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, TARGET_WINDOW)
    target = target.set_index("gis_domain")["target_lrr_m_yr"]
    sv = pd.read_csv(SURVEY_CRESTS)
    sv = sv[sv.status == "ok"].copy()
    sv["row"] = pd.to_datetime(sv.model_state).dt.year - (RUN_START + 1)   # row 0: 1 Jan of the first model year

    frames, summary, crests = {"coastsat": pd.DataFrame({"lrr": target})}, [], []
    for v in ["control", *VARIANTS]:
        rd = run_dir(v)
        if rd is None:
            print(f"  missing {v}")
            continue
        c = load_saved_cascade(rd / f"{rd.name}.npz")
        dom = DOM.for_total(len(c.barrier3d))
        rate = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")
        ow = np.load(next(rd.glob("*_overwash.npy")))
        meta = json.loads((rd / f"{rd.name}_run_metadata.json").read_text(encoding="utf-8"))
        model = {}
        for g in range(1, 91):
            b = c.barrier3d[dom.gis_to_pad(g)]
            model[g] = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + MHW_NAVD88_M
        frames[v] = pd.DataFrame([dict(gis=g, lrr=rate.lrr_m_yr[g], overwash=ow[1:, dom.gis_to_pad(g)].sum())
                                  for g in range(1, 91)]).set_index("gis")
        # Crests in each survey's model state, cell for cell
        for (survey, row), grp in sv.groupby(["survey", "row"]):
            if row < 0 or row >= len(model[1]):
                continue
            m = np.array([model[g][row][cell] for g, cell in zip(grp.gis, grp.cell)])
            o = grp.crest_m_navd88.to_numpy()
            mg = pd.Series(m).groupby(grp.gis.to_numpy()).median()
            og = pd.Series(o).groupby(grp.gis.to_numpy()).median()
            crests.append(dict(run=v, survey=survey, model_median=np.median(m), lidar_median=np.median(o),
                               model_p10=np.percentile(m, 10), lidar_p10=np.percentile(o, 10),
                               model_below_4m=np.mean(m < 4), lidar_below_4m=np.mean(o < 4),
                               r_domains=np.corrcoef(mg, og)[0, 1]))
        lrr = frames[v].lrr
        sk = meta["skill"]
        summary.append(dict(run=v, bias=float(sk["mean_bias_interior_m_yr"]), rmse=float(sk["rmse_interior_m_yr"]),
                            rmse_gis_8_15=float(np.sqrt(((lrr.loc[8:15] - target.loc[8:15]) ** 2).mean())),
                            rmse_gis_80_89=float(np.sqrt(((lrr.loc[80:89] - target.loc[80:89]) ** 2).mean())),
                            island_overwash=frames[v].overwash.loc[2:89].sum()))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    summary = pd.DataFrame(summary).set_index("run")
    crests = pd.DataFrame(crests)
    out = pd.concat(frames, axis=1)
    print(summary.round(3))
    print(out.xs("lrr", axis=1, level=1).loc[list(REPORT_GIS) + list(range(80, 90))].round(2))
    print(crests.pivot(index="survey", columns="run", values="r_domains").round(2))
    print(crests.pivot(index="survey", columns="run", values="model_p10").round(2).assign(
        lidar=crests.groupby("survey").lidar_p10.first().round(2)))
    summary.to_csv(EXP_DIR / "collision_runs_summary.csv")
    out.to_csv(EXP_DIR / "collision_runs_by_gis.csv")
    crests.round(4).to_csv(EXP_DIR / "collision_runs_crests.csv", index=False)
    print(f"\nwrote {EXP_DIR / 'collision_runs_summary.csv'}")


# Figures

# Label, colour and line width per run, all solid: hue for the rule (purple Larson, orange strong, green beach width; blue is
# CoastSat), light for the sand
# lost and dark for the sand on the beach; drawn in this order, so later runs sit on top where they coincide
RUN_STYLES = {
    "control": ("No collision (control)", "#8C8C8C", 1.6),
    "larson_lost": ("Cs 1.7e-4, sand lost", "#BCBDDC", 1.6),
    "larson_beach": ("Cs 1.7e-4, sand to beach", "#6A51A3", 1.8),
    "strong_lost": ("Cs 4e-4, sand lost", "#FDAE6B", 1.6),
    "strong_beach": ("Cs 4e-4, sand to beach", "#D94801", 1.8),
    "width_lost": ("Beach width, sand lost", "#A1D99B", 1.6),
    "width_beach": ("Beach width, sand to beach", "#238B45", 1.8),
    "width_storms_beach": ("Beach width, storm years, sand to beach", "#00441B", 1.8),
}
# (file tag, GIS range or None for the whole island, y half-range in m)
VIEWS = (("island", None, 100.0), ("buxton", (1, 25), 60.0), ("rodanthe", (70, 90), 60.0))
# Map features with nothing to point at in a zoomed view
ZOOM_HIDDEN_FEATURES = {"Shoal influence zone", "Village center", "Pier"}


# Every run's end-minus-start position change against the CoastSat net change, in the layout of the run's
# shoreline_position_change figure (net_change_comparison.py), as the forced-dune-crests plot does
def plot(variants=None):
    import dataclasses
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.lines import Line2D
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, build_coastsat_series
    from cascade_pipeline.plotting.comparison_style import restyle
    from cascade_pipeline.plotting.net_change_comparison import LOWESS_CONFIG, smooth_model
    from cascade_pipeline.plotting.rate_comparison import (
        DEFAULT_RATE_COMPARISON, plot_annotated_rate_comparison)
    from cascade_pipeline.run_info import RunInfo
    from cascade_pipeline.run_layout import resolve
    from cascade_pipeline.shoreline import compute_change_rate
    from site_layer.hat_figure_style import mark_offaxis, offaxis_clause, record_caption
    from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS, WINDOW_ROLE, net_change_dir
    from site_layer.hatteras_site_config import HATTERAS_ANNOTATIONS, HATTERAS_DOMAINS

    # End minus start (m, seaward positive), smoothed as the run figure smooths it
    runs = {}
    for v in ["control", *(variants or [v for v in RUN_STYLES if v != "control"])]:
        rd = run_dir(v)
        if rd is None:
            print(f"  missing {v}")
            continue
        shoreline_m = np.load(resolve(rd, "matrix", rd.name))
        runs[v] = (rd, smooth_model(compute_change_rate(shoreline_m, span_years=1, flip_sign=True)))
    if "control" not in runs:
        print("  the control run is needed to draw the figures")
        return

    tag = f"{RUN_START}_{RUN_END}"
    cs_series = build_coastsat_series(
        [CoastSatDataset(label=f"CoastSat net change ({RUN_START}–{RUN_END})", period_start=RUN_START,
                         csv_path=str(net_change_dir(RUN_START, RUN_END) / f"transect_net_change_{tag}.csv"),
                         rate_col="net_change_m")],
        active_period_start=RUN_START, lowess_config=LOWESS_CONFIG, domains=HATTERAS_DOMAINS)

    control_dir, control = runs["control"]
    label_c, color_c, _ = RUN_STYLES["control"]
    gis = np.arange(HATTERAS_DOMAINS.first_gis_id, HATTERAS_DOMAINS.last_gis_id + 1)
    geometry = HATTERAS_DOMAINS.for_total(len(control))     # the runs' own buffers, not the config's
    real = slice(geometry.start_real_index, geometry.end_real_index)
    role = WINDOW_ROLE.get((RUN_START, RUN_END), "").replace(" (held out)", "")

    for view, gis_range, y_half in VIEWS:
        where = f"GIS {gis_range[0]}–{gis_range[1]}" if gis_range else "whole island"
        subtitle = " · ".join(p for p in (role, where, "full management + groin, edgeBE",
                                          "collision-regime dune erosion") if p)
        config = dataclasses.replace(
            DEFAULT_RATE_COMPARISON, quantity="position", publication_text=True,
            show_features=True, plot_domain_means=False, raw_lrr_southern_only=True,
            smoothing_marker=False, observed_label="CoastSat (7-domain LOWESS)", model_label=label_c,
            title=f"Shoreline position change with dune collision erosion, {RUN_START}–{RUN_END}",
            ylim=(-y_half, y_half), subtitle=subtitle)
        run = RunInfo(run_name=control_dir.name, run_dir=str(control_dir), start_year=RUN_START,
                      end_year=RUN_END, Hs=2.0, flip_sign_model=True)
        # The control is drawn by the run figure's own call, in the control's colour
        fig, ax = plot_annotated_rate_comparison(
            control, cs_series, run, estimator=None, domains=HATTERAS_DOMAINS,
            annotations=dataclasses.replace(HATTERAS_ANNOTATIONS, model_color=color_c),
            lowess_config=LOWESS_CONFIG, config=config)

        in_view = np.ones_like(gis, bool) if gis_range is None else (
            (gis >= gis_range[0]) & (gis <= gis_range[1]))
        named = [(label_c, mark_offaxis(ax, gis[in_view], control[real][in_view], y_half, color=color_c))]
        handles = []
        for v, (rd, change) in runs.items():
            label, color, lw = RUN_STYLES[v]
            handles.append(Line2D([0], [0], color=color, lw=lw + 0.4, label=label))
            if v == "control":
                # Drawn by the run figure's call; given its width here, under the collision runs
                for line in ax.get_lines():
                    if line.get_label() == label_c:
                        line.set_linewidth(lw)
                continue
            ax.plot(gis, change[real], color=color, lw=lw, zorder=6)
            named.append((label, mark_offaxis(ax, gis[in_view], change[real][in_view], y_half, color=color)))

        if gis_range is not None:
            # Zoom: hide the labels of features outside the view; centre the transects on their domain
            ax.set_xlim(gis_range[0] - 0.5, gis_range[1] + 0.5)
            ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
            ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
            for txt in ax.texts:
                x_is_data = txt.get_transform() not in (ax.transAxes, fig.transFigure)
                if x_is_data and not gis_range[0] - 0.5 <= txt.get_position()[0] <= gis_range[1] + 0.5:
                    txt.set_visible(False)
            for coll in ax.collections:
                if "transect" in str(coll.get_label()):
                    coll.set_offsets(coll.get_offsets() - [0.5, 0.0])

        # The figure's legend with every run in place of the one model line, then the house finish
        old = fig.legends.pop()
        hidden = {label_c} | (set() if gis_range is None else ZOOM_HIDDEN_FEATURES)
        entries = [(h, t.get_text()) for h, t in zip(old.legend_handles, old.get_texts())
                   if t.get_text() not in hidden]
        fig.legend(handles + [h for h, _ in entries],
                   [h.get_label() for h in handles] + [t for _, t in entries],
                   loc="outside lower center")
        restyle(fig, ax, subtitle)

        out = EXP_DIR / f"collision_runs_position_{view}.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close("all")
        (s0, s1), (e0, e1) = NET_CHANGE_WINDOWS[(RUN_START, RUN_END)]
        record_caption(out, (
            f"Shoreline position change, {RUN_START}–{RUN_END}, {where}; positive is seaward. Modelled "
            f"(CASCADE): the shoreline at the end of each run minus the start, smoothed with a 7-domain "
            f"LOWESS north of domain 10 as in the run figure, for the full-management edgeBE run with the "
            f"groin without collision erosion (grey) and with every dune cell lowered after each model year "
            f"by Cs·4Σ(TWL − toe)²Δt/T / 20 m over the hours the total water level stood between its lidar toe "
            f"and its crest (Larson et al. 2004): Cs 1.7e-4, the Dennis fit (purple), or 4e-4, what the Dennis "
            f"hotspots need (orange); or, green, in each year whose highest water stood between the toe and the "
            f"crest, by a·max(0, 1 − W/W0) of the crest elevation, W the last survey's beach width (the fit on "
            f"Dennis and 2004–05; darkest, only in years whose highest water reached 3.1 m NAVD88); light, the "
            f"eroded sand lost; dark, placed on the beach as CASCADE's shoreface nourishment places fill. Runs "
            f"drawn: {', '.join(RUN_STYLES[v][0] for v in runs)}. Where runs coincide the last drawn sits on top. "
            f"Observed "
            f"(CoastSat, blue dots and line): per transect, the mean position over {e0} to {e1} minus the "
            f"mean over {s0} to {s1}, smoothed with a 7-domain LOWESS north of domain 10; dots south of it. "
            f"The y axis is fixed at ±{y_half:g} m."
            + offaxis_clause([(k, n) for k, n in named if n], y_half)
            + f" Experiment {TAG}."))
        print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("action", choices=["toes", "twl", "calibrate", "empirical", "run", "score", "plot", "_launch"])
    ap.add_argument("variant", nargs="*")
    args = ap.parse_args()
    if args.action == "_launch":
        _launch(args.variant[0])
    elif args.action == "run":
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, args.variant or list(VARIANTS)))
    elif args.action == "plot":
        plot(args.variant or None)
    else:
        {"toes": toes, "twl": twl, "calibrate": calibrate, "empirical": empirical, "score": score}[args.action]()


if __name__ == "__main__":
    main()
