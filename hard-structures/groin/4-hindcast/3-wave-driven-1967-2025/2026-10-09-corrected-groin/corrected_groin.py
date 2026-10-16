#!/usr/bin/env python3
"""
Fit the corrected groin (the drift groin plus a leak correction at the groin's own bends) on 1996-2009.

    python corrected_groin.py check                       # b 0: must reproduce the no-groin run
    python corrected_groin.py smoke  [--b 0.1 --f 0.6]    # one calibration cell
    python corrected_groin.py grid   [--streams 4]        # the calibration grid
    python corrected_groin.py edge   [--streams 4]        # past the grid's f 0.8 edge, and b 0.25
    python corrected_groin.py score                       # -> grid_scores.csv
    python corrected_groin.py schedules [--streams 4]    # schedules A (instant1996) and C (ramp1996) on 1996-2009
    python corrected_groin.py tests [--streams 4]        # each schedule's best cell on 2009-2025, with and without D
    python corrected_groin.py test --b 0.1 --f 0.6 [--schedule ramp1996] [--off 2017]   # one cell on 2009-2025
    python corrected_groin.py periods --period 2009 --spec "1969:0.2,2004:0.16,2017:0"  # strength per period, set directly
    python corrected_groin.py tracked [--streams 4]      # the tracked reference against the no-groin run, same setup
    python corrected_groin.py refresh [--streams 4]      # every run the README reports, then the score
    python corrected_groin.py --reference baseline score # the 10-09 runs, with the no-groin run as reference

Runs are filed by setup: the calendar-year setup (default since the 2026-10-12 evening) under
groin/2026-10-12-corrected-groin-calendar/, the earlier ones under groin/2026-10-09-corrected-groin/
(--setup six-month). The earlier scores are kept as written (grid_scores_2026-10-12_six-month.csv,
grid_scores_2026-10-09_baseline.csv); rescoring them now would not reproduce them, since the
config's test-period length and the no-groin matrix runs have since changed.

The drift groin of ../../2-blocking-1996-2025/2026-10-08-drift-groin/ (b of BRIE's wave-climate net drift at the updrift
approach face, moved from GIS 5 to GIS 6 each year), plus straight-coast module 6's leak
correction within HALF_WIDTH cells of the groin face, applied to the groin's bends only (option
B): each year the correction is computed for this run's BRIE shoreline and for the groin-free
shoreline in the same year, and only the difference is added to x_s_dt,

    L(x) = 2 (r_i - r_{i-1}) (x_i - x_{i-1}),   x_s_dt += L(x_groin) - L(x_none)   near the groin

The groin-free shoreline x_none is, by default ("tracked", since 2026-10-12), this run's shoreline
minus the groin's own footprint: each year the callback repeats BRIE's alongshore step for last
year with and without its own contribution and keeps the difference, so fills and Barrier3D's
change cancel and no no-groin run is needed; with no groin the footprint is exactly zero. With
--reference baseline (the 10-09 runs) x_none is the no-groin run's shoreline instead. BRIE and
Barrier3D are not changed: the class replaces BlockingGroinCallback in a child process, and the
unchanged runner is driven as the 10-05 setup.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-12
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
REPO = next(p for p in HERE.parents if (p / "pyproject.toml").exists())
GROIN = REPO / "hard-structures" / "groin"
sys.path[:0] = [str(HERE.parent.parent / "2-blocking-1996-2025" / "2026-10-08-drift-groin"), str(HERE.parent.parent / "2-blocking-1996-2025" / "2026-10-08-schedule-refit"),
                str(GROIN / "1-observations" / "analysis" / "2-coastsat_gap"),
                str(GROIN / "3-module-tests" / "8-real-planform" / "6-instant-2004-grid"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]
import drift_groin as dg  # noqa: E402
import schedule_refit as sr  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RUNNER = dg.RUNNER
# Which setup's runs: "calendar" (default since the 2026-10-12 evening: calendar-year start and end lines, storms
# from the year after the DEM, row k of a run = calendar year period + k) or "six-month" (the 10-12 afternoon
# refresh: DEM-date +/-6-month start lines, storms from the DEM year, row k = 1 Jan of period + k). --setup sets it
SETUP = "calendar"
STUDIES = {"calendar": "groin/2026-10-12-corrected-groin-calendar", "six-month": "groin/2026-10-09-corrected-groin"}
STUDY = STUDIES[SETUP]
RAW = dg.RAW
LOGS = HERE / "logs" / "calendar"
HALF_WIDTH = 10                           # cells either side of the groin face; 5 leaks on a grid-aligned coast
B_VALUES = (0.02, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3)   # the drift groin's best sat at b 0.08-0.10
F_VALUES = (0.2, 0.4, 0.6, 0.8)
B_EDGE = (0.15, 0.2, 0.25, 0.3)              # the first grid's best sat at f 0.8, its edge
F_EDGE = (0.9, 1.0)                          # f 1.0: no failure in 2004 (a diagnostic, not a schedule)
SMOKE = (0.1, 0.6)
CALIBRATION, TEST = dg.CALIBRATION, dg.TEST
SCHEDULES = sr.SCHEDULES                     # instant2004 (B, the runner's), instant1996 (A), ramp1996 (C)
A_STRENGTHS = (0.02, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.25, 0.3)   # A: one strength b f all window; run as b, f 1.0
C_B = (0.1, 0.15, 0.2, 0.25, 0.3, 0.4)       # C: full b in 1995, ramping to b f by 2003
C_F = (0.2, 0.4, 0.6, 0.8)
OFF_YEAR = 2017                              # D: the groin passes everything from the 2017 Buxton fill on
DECLINES = ("2015-2021r0", "2015-2021r0.33")  # D': linear decline 2015-2021 to nothing, or to the southernmost of three groins
TWO_FINE = (0.11, 0.13, 0.14, 0.16, 0.17, 0.18, 0.19)   # two-number strengths between schedule A's 0.1 / 0.12 / 0.15 / 0.2
# How the groin-free shoreline for the leak correction is found: "tracked" (default since 2026-10-12) subtracts the
# groin's own tracked footprint, so no no-groin run is needed; "baseline" uses the no-groin run (the 10-09 runs);
# "now" is "baseline" rerun under the 10-12 setup, filed apart for the comparison. --reference sets it
REFERENCE = "tracked"
# -----------------------------------------------------------------------------


# The share of the scheduled strength left in `year` under an end spec: None, a switch-off year ("2017"),
# or a linear decline "Y0-Y1rR" from full at Y0 to R at Y1, held after
def end_factor(spec, year):
    if not spec:
        return 1.0
    spec = str(spec)
    if "-" not in spec:
        return 0.0 if year >= int(spec) else 1.0
    y0, rest = spec.split("-")
    y1, r = rest.split("r")
    y0, y1, r = int(y0), int(y1), float(r)
    if year <= y0:
        return 1.0
    return r if year >= y1 else 1.0 - (1.0 - r) * (year - y0) / (y1 - y0)


# Strength per period, set directly: "1969:0.2,2004:0.16,2017:0" holds 0.2 from 1969, 0.16 from 2004, 0 from 2017
def parse_periods(spec):
    pairs = sorted((int(y), float(v)) for y, v in (item.split(":") for item in str(spec).split(",")))
    return pairs


def strength_for(spec, year):
    value = 0.0
    for y, v in parse_periods(spec):
        if year >= y:
            value = v
    return value


def periods_label(spec):
    return "_".join(f"{y}-{v:g}" for y, v in parse_periods(spec))


# Where a schedule's runs are filed; the 2004 schedule without D keeps the first grid's place.
# Runs with the tracked reference (no no-groin run needed) are filed apart, prefixed "tracked"; "now" is the
# no-groin-run reference rerun under the 2026-10-12 setup (six-month start lines, calendar end windows),
# so the two can be compared on the same setup
PREFIX = {"baseline": "", "tracked": "tracked_", "now": "now_"}


def scores_file():
    if REFERENCE == "baseline":
        return HERE / "grid_scores_2026-10-09_baseline.csv"
    return HERE / ("grid_scores.csv" if SETUP == "calendar" else "grid_scores_2026-10-12_six-month.csv")


def tag(schedule="instant2004", off=None, reference=None):
    pre = PREFIX[reference or REFERENCE]
    if schedule == "periods":
        return f"{STUDY}/{pre}periods_{periods_label(off)}"     # the runner allows four path components
    if schedule == "instant2004" and off is None:
        return f"{STUDY}/{pre[:-1]}" if pre else STUDY
    return f"{STUDY}/{pre}{schedule}" + (f"_off{off}" if off else "")


def run_dir(study, period, b, f):
    hits = sorted((RAW / "experiments" / study / dg.member(b, f)).glob(f"{dg.window(period)}/edgeBE/*/"))
    return hits[-1] if hits else None


def shoreline_file(study, period, b, f):
    d = run_dir(study, period, b, f)
    hits = sorted(d.glob("*_shoreline_matrix.npy")) if d else []
    return hits[-1] if hits else None


# Module 6's leak term for a BRIE shoreline x (landward +), from BRIE's own diffusivity table
def leak(brie, x):
    ny = len(x)
    step = np.roll(x, -1) - x                                   # x[i+1] - x[i]
    theta = 180.0 * np.arctan2(step, brie._dy) / np.pi
    k = np.clip(np.round(90 - theta).astype(int), 1, brie._wave_climl - 1)
    r = np.maximum(0.0, brie._coast_diff[k] * brie._dt / 2 / brie._dy ** 2)
    assert r.shape == (ny,)
    return 2 * (r - np.roll(r, 1)) * np.roll(step, 1)


# BRIE's implicit alongshore step (brie.py:1293-1326) for a shoreline x and shoreline change dx, m
def brie_step(brie, x, dx):
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import spsolve
    theta = 180.0 * np.arctan2(np.roll(x, -1) - x, brie._dy) / np.pi
    k = np.maximum(1, np.minimum(brie._wave_climl, np.round(90 - theta).astype(int)))
    r = np.maximum(0.0, brie._coast_diff[k] * brie._dt / 2 / brie._dy ** 2)
    dv = np.r_[-r[-1], -r[1:], 1 + 2 * r, -r[0:-1], -r[0]]
    rhs = x + r * (np.roll(x, -1) - 2 * x + np.roll(x, 1)) + dx
    return spsolve(csr_matrix((dv, (brie._di, brie._dj))), rhs)


# Child process: the corrected groin in BlockingGroinCallback's place, then the unchanged runner
def _launch():
    import runpy
    import cascade.groin as cg
    reference = os.environ.get("CORRECTED_GROIN_REFERENCE", "baseline")
    base_path = os.environ.get("CORRECTED_GROIN_BASELINE", "")
    # the no-groin run, [year x domain], m, landward +: the reference for "baseline", a check only for "tracked"
    base_m = np.load(base_path) if base_path and Path(base_path).is_file() else None
    schedule = os.environ.get("CORRECTED_GROIN_SCHEDULE", "instant2004")
    delay, mode, ramp = SCHEDULES.get(schedule, SCHEDULES["instant2004"])
    off = os.environ.get("CORRECTED_GROIN_OFF") or None

    class CorrectedGroin(cg.BlockingGroinCallback):
        kind = "blocking"

        def __init__(self, *a, **kw):
            kw.update(deterioration_delay_years=delay, deterioration_mode=mode, deterioration_ramp_years=ramp)
            super().__init__(*a, **kw)

        def __call__(self, cascade, x_s_dt):
            year = self.start_year + self._call_count
            k = self._call_count
            self._call_count += 1
            active = year >= self.install_year
            if schedule == "periods":                    # the strength set directly per period; b and f unused
                b = strength_for(off, year) if active else 0.0
            else:
                b = self._effective_trapping_rate(year) * end_factor(off, year) if active else 0.0
            brie = cascade._brie_coupler._brie
            x_s = np.asarray(brie.x_s, dtype=float)
            ny, lo, hi = len(x_s), self._lo, self._hi
            if not hasattr(self, "_q_net"):
                self._q_net = dg.drift_table(brie)
                self._x0 = x_s.copy()
                self._a = np.zeros(ny)                    # tracked: the groin's own footprint on the shoreline
                self._prev = None
                self.drift_m3_yr_TS, self.approach_deg_TS, self.correction_m_TS, self.baseline_gap_m_TS = [], [], [], []
                self.footprint_m_TS, self.tracked_vs_baseline_m_TS = [], []
            if reference == "tracked" and self._prev is not None:
                # last year's step with and without the groin's own contribution; fills and Barrier3D's change
                # enter both alike and cancel, so with no groin the footprint stays exactly zero
                xp, dxp, gp, ap = self._prev
                self._a = brie_step(brie, xp, dxp) - brie_step(brie, xp - ap, dxp - gp)

            def angle(i):
                return 180.0 * np.arctan2(x_s[(i + 1) % ny] - x_s[i], brie._dy) / np.pi

            def q_net(theta):
                return float(self._q_net[int(np.clip(np.round(90 - theta), 1, brie._wave_climl - 1))])

            # The drift groin, as ../../2-blocking-1996-2025/2026-10-08-drift-groin/drift_groin.py
            toward_lower = q_net(angle(hi)) < 0
            approach = angle(hi) if toward_lower else angle(lo - 1)
            qv = q_net(approach)
            q = b * qv * brie._dt / (brie._dy * (brie._h_b_crit + brie._d_sf))
            dx_lo, dx_hi = (-q, q) if (active and b > 0.0) else (0.0, 0.0)
            x_s_dt[lo] += dx_lo
            x_s_dt[hi] += dx_hi

            # The leak correction, on the groin's bends only: this run minus the groin-free shoreline, taken from
            # the no-groin run ("baseline") or as this run minus the groin's tracked footprint ("tracked")
            x_base = None if base_m is None else self._x0 + (base_m[min(k, len(base_m) - 1)] - base_m[0])
            x_none = x_s - self._a if reference == "tracked" else x_base
            near = np.abs(np.arange(ny) - (lo + 0.5)) <= HALF_WIDTH
            corr = np.where(near, leak(brie, x_s) - leak(brie, x_none), 0.0)
            x_s_dt += corr
            g = corr.copy()
            g[lo] += dx_lo
            g[hi] += dx_hi
            self._prev = (x_s.copy(), np.asarray(x_s_dt, dtype=float).copy(), g, self._a.copy())
            self.footprint_m_TS.append(float(np.abs(self._a).max()))
            self.tracked_vs_baseline_m_TS.append(
                np.nan if x_base is None else float(np.abs((x_s - self._a) - x_base)[near].max()))

            dx_up, dx_down = (dx_lo, dx_hi) if self.updrift_pad == lo else (dx_hi, dx_lo)
            self.year_TS.append(year)
            self.active_TS.append(bool(active))
            self.blocking_applied_TS.append(b)
            self.applied_dx_updrift_TS.append(dx_up)
            self.applied_dx_downdrift_TS.append(dx_down)
            self.trapping_rate_applied_TS.append(abs(dx_up))
            self.face_offset_m_TS.append(float(x_s[hi] - x_s[lo]))
            self.r_ipl_updrift_TS.append(self._r_ipl(brie, x_s, self.updrift_pad))
            self.r_ipl_downdrift_TS.append(self._r_ipl(brie, x_s, self.downdrift_pad))
            self.drift_m3_yr_TS.append(qv)
            self.approach_deg_TS.append(approach)
            self.correction_m_TS.append(float(corr.sum()))
            self.baseline_gap_m_TS.append(float(np.abs(x_s - x_none).max()))
            return x_s_dt

        def diagnostics_frame(self):
            rows = super().diagnostics_frame()
            extra = zip(getattr(self, "drift_m3_yr_TS", []), getattr(self, "approach_deg_TS", []),
                        getattr(self, "correction_m_TS", []), getattr(self, "baseline_gap_m_TS", []),
                        getattr(self, "footprint_m_TS", []), getattr(self, "tracked_vs_baseline_m_TS", []))
            for row, (qv, th, c, g, fp, tv) in zip(rows, extra):
                row.update(net_drift_m3_yr=qv, approach_angle_deg=th, leak_correction_m=c,
                           max_diff_from_no_groin_m=g, footprint_max_m=fp, tracked_vs_baseline_m=tv)
            return rows

    cg.BlockingGroinCallback = CorrectedGroin
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(period, b, f, schedule="instant2004", off=None, reference=None):
    from schedule_refit import baseline_file
    reference = reference or REFERENCE
    study = tag(schedule, off, reference)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "CORRECTED_GROIN_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(b), "HAT_GROIN_DETERIORATION_FRACTION": str(f),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{study}/{dg.member(b, f)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
        "CORRECTED_GROIN_BASELINE": str(baseline_file(period)),
        "CORRECTED_GROIN_SCHEDULE": schedule, "CORRECTED_GROIN_OFF": str(off or ""),
        "CORRECTED_GROIN_REFERENCE": reference,
    })
    LOGS.mkdir(exist_ok=True)
    label = study.replace(STUDY, "").strip("/").replace("/", "_")
    log = LOGS / f"{label + '_' if label else ''}{dg.window(period)}_{dg.member(b, f)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO,
                               env=env, stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(study, period, b, f)
    if code == 0 and d is not None:
        delay, mode, ramp = SCHEDULES.get(schedule, (None, None, None))
        (d / "corrected_groin.txt").write_text(
            "patched in-process by corrected_groin.py: the drift groin (b of BRIE's net drift at the updrift "
            f"approach face, GIS 5 -> GIS 6) plus module 6's leak correction within {HALF_WIDTH} cells, applied "
            + ("as this run's minus that of this run less the groin's tracked footprint (no no-groin run used); "
               if reference == "tracked" else "as this run's minus the no-groin run's; ")
            + "the runner's report and groin labels do not show it. "
            f"Schedule {schedule}: deterioration {delay} yr after the 1969 install, {mode}, ramp {ramp} yr"
            + (f"; end of the groin's working life: {off} (see end_factor)" if off and schedule != "periods" else "")
            + (f"; SCHEDULE IGNORED, strength set per period: {off}" if schedule == "periods" else "") + ".")
    print(f"{time.strftime('%H:%M:%S')}  {label or 'instant2004'} {dg.window(period)} {dg.member(b, f)}: "
          f"exit {code} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


# cells: (period, b, f), (period, b, f, schedule, off) or (period, b, f, schedule, off, reference)
def run_all(cells, streams):
    cells = [tuple(c) + ("instant2004", None, REFERENCE)[len(c) - 3:] for c in cells]
    todo = [c for c in cells if shoreline_file(tag(*c[3:]), *c[:3]) is None]
    print(f"{len(cells) - len(todo)} of {len(cells)} done; running {len(todo)} in {streams} streams", flush=True)
    jobs, failed = queue.Queue(), []
    for c in todo:
        jobs.put(c)

    def worker():
        while True:
            try:
                c = jobs.get_nowait()
            except queue.Empty:
                return
            if run(*c):
                failed.append(c)

    threads = [threading.Thread(target=worker) for _ in range(streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}", flush=True)
    return not failed


# b 0 must reproduce the no-groin run: the correction cancels and the groin moves nothing
def cmd_check(_a):
    from schedule_refit import baseline_file
    ok = run_all([(CALIBRATION, 0.0, 0.6)], 1)
    p = shoreline_file(tag(), CALIBRATION, 0.0, 0.6)
    if not ok or p is None:
        sys.exit(1)
    x, base = np.load(p), np.load(baseline_file(CALIBRATION))
    diag = pd.read_csv(p.parent / "tables" / "groin_diagnostics.csv")
    print(f"max |shoreline - no-groin| over all years and domains: {np.abs(x - base).max():.3e} m")
    print(f"max |leak correction| applied in any year: {diag.leak_correction_m.abs().max():.3e} m")
    print(f"max |rebuilt no-groin BRIE shoreline - this run's|: {diag.max_diff_from_no_groin_m.max():.3e} m")


# The schedule grids: A as one strength (b, f 1.0), C as b x f
def schedule_cells():
    return ([(CALIBRATION, s, 1.0, "instant1996", None) for s in A_STRENGTHS]
            + [(CALIBRATION, b, f, "ramp1996", None) for b in C_B for f in C_F])


# Scoring. "six-month" runs are scored as schedule_refit.py scores them (row k = 1 Jan of period + k, so year y is
# the mean of rows y - period and y - period + 1). "calendar" runs: row k is calendar year period + k (row 0 the
# start line), compared with CoastSat's calendar-year mean of that year minus its calendar-`period` mean, for the
# years after the start through the last modelled year (the handoff of 2026-10-12)
def observed_change(period):
    if SETUP == "calendar":
        from site_layer.hat_observed_rates import NET_CHANGE_WINDOWS
        from site_layer.hatteras_site_config import HATTERAS_PERIODS
        lo, hi = NET_CHANGE_WINDOWS[(period, HATTERAS_PERIODS[period]["end_year"])][0]
        sr.START_WINDOW = {**sr.START_WINDOW, period: (lo, f"{hi} 23:59:59")}
    return sr.coastsat_change(period)


def model_change(path):
    if SETUP != "calendar":
        return sr.model_change(path)
    x = np.load(path)
    g = x[:, sr.DOWN] - x[:, sr.UP]
    return g - g[0]


def score_annual(path, period, obs):
    if SETUP != "calendar":
        return sr.score_annual(path, period, obs)
    from site_layer.hatteras_site_config import last_model_year
    m = model_change(path)
    yrs = [y for y in range(period + 1, last_model_year(period) + 1) if y in obs.index]
    mv, ov = np.array([m[y - period] for y in yrs]), obs.loc[yrs].to_numpy()
    return dict(n_years=len(yrs), rmse_annual=float(np.sqrt(np.mean((mv - ov) ** 2))),
                bias_annual=float(np.mean(mv - ov)), end_model=float(mv[-1]), end_obs=float(ov[-1]))


def score_photos(path, period, photos):
    if SETUP != "calendar":
        return sr.score_photos(path, period, photos)
    from site_layer.hatteras_site_config import last_model_year
    g = model_change(path)
    yrs = [y for y in photos.index if period < y <= last_model_year(period)]
    o = np.array([photos[y] for y in yrs]) - np.interp(period, photos.index, photos.values)
    m = np.array([g[y - period] for y in yrs])
    return dict(rmse_photos=float(np.sqrt(np.mean((m - o) ** 2))), photo_dates=" ".join(map(str, yrs)))


# Annual RMSE before and after the 2017 fill, test period only. Calendar runs score 2010-2016 and 2017-2025
# (2009 is the start line, zero in both); six-month runs 2009-2016 and 2017-2024, as before
SPLITS = {"calendar": (("rmse_2010_2016", (2010, 2016)), ("rmse_2017_2025", (2017, 2025))),
          "six-month": (("rmse_2009_2016", (2009, 2016)), ("rmse_2017_2024", (2017, 2024)))}


def split_scores(path, period, obs):
    if period != TEST:
        return {}
    m = model_change(path)
    out = {}
    for name, (y0, y1) in SPLITS[SETUP]:
        yrs = [y for y in range(y0, y1 + 1) if y in obs.index and y - period < len(m)]
        d = np.array([m[y - period] for y in yrs]) - obs.loc[yrs].to_numpy()
        out[name] = float(np.sqrt(np.mean(d ** 2)))
    return out


def cmd_score(_a):
    from schedule_refit import baseline_file
    from score_instant_grid import observed_series
    photos = observed_series()
    arms = [("corrected", "instant2004", None, tag())] + [
        ("corrected", s, off, tag(s, off)) for s in SCHEDULES for off in (None, str(OFF_YEAR)) + DECLINES if (s, off) != ("instant2004", None)]
    bs = sorted(set(B_VALUES) | set(B_EDGE) | set(A_STRENGTHS) | set(C_B) | set(dg.B_VALUES) | set(dg.B_FINE))
    fs = sorted(set(F_VALUES) | set(F_EDGE) | set(C_F) | set(dg.F_VALUES))
    rows = []
    for period in (CALIBRATION, TEST):
        obs = observed_change(period)
        base = baseline_file(period)
        if base.is_file():
            rows.append(dict(arm="none", schedule="", off=np.nan, period=period, b=0.0, f=np.nan,
                             **score_annual(base, period, obs), **score_photos(base, period, photos),
                             **dg.near_change(base), **split_scores(base, period, obs)))
        finders = [(arm, s, off, study, shoreline_file) for arm, s, off, study in arms]
        pre = PREFIX[REFERENCE] + "periods_"
        for d in sorted((RAW / "experiments" / STUDY).glob(pre + "*")):
            spec = ",".join(item.replace("-", ":", 1) for item in d.name[len(pre):].split("_"))
            finders.append(("corrected", "periods", spec, tag("periods", spec), shoreline_file))
        if REFERENCE == "baseline":                 # the drift groin was only run on the 10-09 setup
            finders.append(("drift", "instant2004", None, dg.STUDY, dg.shoreline_file))
        for arm, s, off, study, finder in finders:
            # a periods run is filed under its first strength with f 1.0
            cells = [(parse_periods(off)[0][1], 1.0)] if s == "periods" else [(b, f) for b in bs for f in fs]
            for b, f in cells:
                if True:
                    p = finder(study, period, b, f)
                    if p is not None:
                        rows.append(dict(arm=arm, schedule=s, off=off if off else "", period=period, b=b, f=f,
                                         **score_annual(p, period, obs), **score_photos(p, period, photos),
                                         **dg.budget(p), **dg.near_change(p), **split_scores(p, period, obs)))
    t = pd.DataFrame(rows)
    for period in (CALIBRATION, TEST):
        none = t[(t.arm == "none") & (t.period == period)]
        if len(none):
            sel = t.period == period
            t.loc[sel, "island_vs_none_m"] = t.loc[sel, "island_seaward_sum_m"] - float(none.island_seaward_sum_m.iloc[0])
    t["joint_rmse"] = np.sqrt((t.rmse_annual ** 2 + t.rmse_photos ** 2) / 2)
    t["strength_after_failure"] = t.b * t.f
    t.to_csv(scores_file(), index=False, float_format="%.3f")
    pd.set_option("display.width", 240)
    cal = t[(t.period == CALIBRATION) & (t.arm == "corrected")]
    print("=== calibration, best joint per schedule (corrected groin) ===")
    print(cal.sort_values("joint_rmse").groupby("schedule").head(1)[
        ["schedule", "b", "f", "rmse_annual", "rmse_photos", "joint_rmse", "island_vs_none_m"]].round(1).to_string(index=False))
    print("=== calibration, best annual and best photos per schedule ===")
    for sc_ in ("rmse_annual", "rmse_photos"):
        print(cal.sort_values(sc_).groupby("schedule").head(1)[
            ["schedule", "b", "f", "rmse_annual", "rmse_photos", "joint_rmse"]].round(1).to_string(index=False))
    test = t[t.period == TEST]
    print("=== test 2009-2025 ===")
    print(test[["arm", "schedule", "off", "b", "f", "strength_after_failure", "rmse_annual", "rmse_photos", "joint_rmse",
                *[n for n, _ in SPLITS[SETUP]], "end_model", "end_obs"]].round(1).to_string(index=False))


# The tracked reference against the no-groin-run reference: b 0 (must reproduce no groin), and the two-number
# groin at the strengths the annual record (0.12) and the joint score (0.19) pick, in both periods
TRACKED_STRENGTHS = (0.12, 0.19)


def cmd_tracked(a):
    from schedule_refit import baseline_file
    from score_instant_grid import observed_series
    spec = {s: f"1969:{s},2017:0" for s in TRACKED_STRENGTHS}
    cells = [(CALIBRATION, 0.0, 0.6, "instant2004", None, "tracked")] + [
        (p, s, 1.0, "periods", spec[s], ref) for ref in ("tracked", "now")
        for p in (CALIBRATION, TEST) for s in TRACKED_STRENGTHS]
    if not run_all(cells, a.streams):
        sys.exit(1)
    photos, rows = observed_series(), []
    x0, base = np.load(shoreline_file(tag(reference="tracked"), CALIBRATION, 0.0, 0.6)), np.load(baseline_file(CALIBRATION))
    print(f"b 0, tracked: max |shoreline - no-groin run| over all years and domains: {np.abs(x0 - base).max():.3e} m")
    for p in (CALIBRATION, TEST):
        obs = observed_change(p)
        for s in TRACKED_STRENGTHS:
            pt = shoreline_file(tag("periods", spec[s], "tracked"), p, s, 1.0)
            pb = shoreline_file(tag("periods", spec[s], "now"), p, s, 1.0)          # same setup, no-groin reference
            xt, xb = np.load(pt), np.load(pb)
            diag = pd.read_csv(pt.parent / "tables" / "groin_diagnostics.csv")
            for ref, path in (("baseline", pb), ("tracked", pt)):
                sa, sp = score_annual(path, p, obs), score_photos(path, p, photos)
                rows.append(dict(period=p, strength=s, reference=ref, rmse_annual=sa["rmse_annual"],
                                 rmse_photos=sp["rmse_photos"], end_model=sa["end_model"], end_obs=sa["end_obs"]))
            rows[-1].update(max_diff_vs_baseline_m=float(np.abs(xt - xb).max()),
                            gap_diff_vs_baseline_m=float(np.abs(model_change(pt) - model_change(pb)).max()),
                            island_diff_vs_baseline_m=float(np.abs((xt - xb).sum(axis=1)).max()),
                            tracked_vs_baseline_ref_m=float(diag.tracked_vs_baseline_m.max()))
    t = pd.DataFrame(rows)
    t["joint_rmse"] = np.sqrt((t.rmse_annual ** 2 + t.rmse_photos ** 2) / 2)
    t.to_csv(HERE / "tracked_vs_baseline.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 240)
    print(t.round(2).to_string(index=False))


# Every run the README reports, under the current setup and reference: the calibration grid, edge and schedules,
# the two-number calibration strengths, the README's test cells, then the score and each schedule's new best cells
def cmd_refresh(a):
    two = lambda s: f"1969:{s},2017:0"                                   # noqa: E731
    cal = ([(CALIBRATION, b, f) for b in B_VALUES for f in F_VALUES]
           + [(CALIBRATION, b, f) for b in B_EDGE for f in F_EDGE] + [(CALIBRATION, 0.25, f) for f in F_VALUES[1:]]
           + schedule_cells() + [(CALIBRATION, s, 1.0, "periods", two(s)) for s in TWO_FINE])
    test = ([(TEST, s, 1.0, "periods", two(s)) for s in (0.12, 0.18, 0.19, 0.25)]
            + [(TEST, 0.2, 1.0, "periods", "1969:0.2,2004:0.16,2017:0")]
            + [(TEST, 0.2, f, "instant2004", off) for f in (0.8, 0.9) for off in (None, str(OFF_YEAR))]
            + [(TEST, 0.2, 1.0, "instant1996", off) for off in (None, str(OFF_YEAR))]
            + [(TEST, 0.2, 0.8, "ramp1996", off) for off in (None, str(OFF_YEAR))]
            + [(TEST, 0.2, 0.8, "instant2004", d) for d in DECLINES])
    ok = run_all(cal + test, a.streams)
    cmd_score(a)
    t = pd.read_csv(scores_file())
    c = t[(t.period == CALIBRATION) & (t.arm == "corrected") & (t.off.fillna("") == "") & (t.schedule != "periods")]
    best = c.sort_values("joint_rmse").groupby("schedule").head(1)
    extra = [(TEST, float(r.b), float(r.f), r.schedule, off) for r in best.itertuples() for off in (None, str(OFF_YEAR))]
    print("best calibration cells:", [(e[1], e[2], e[3]) for e in extra[::2]], flush=True)
    ok = run_all(extra, a.streams) and ok
    cmd_score(a)
    sys.exit(0 if ok else 1)


# The test cells: each schedule's best calibration cell, with and without D
def cmd_tests(a):
    t = pd.read_csv(scores_file())
    cal = t[(t.period == CALIBRATION) & (t.arm == "corrected") & (t.off.fillna("") == "")]
    best = cal.sort_values("joint_rmse").groupby("schedule").head(1)
    cells = [(TEST, float(r.b), float(r.f), r.schedule, off) for r in best.itertuples() for off in (None, OFF_YEAR)]
    for c in cells:
        print("test cell:", c)
    sys.exit(0 if run_all(cells, a.streams) else 1)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reference", choices=("tracked", "baseline"), default=REFERENCE,
                    help="groin-free shoreline: the groin's tracked footprint (default) or the no-groin run (10-09)")
    ap.add_argument("--setup", choices=tuple(STUDIES), default=SETUP,
                    help="calendar-year setup (default) or the 10-12 afternoon six-month runs")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    s = sub.add_parser("smoke")
    s.add_argument("--b", type=float, default=SMOKE[0])
    s.add_argument("--f", type=float, default=SMOKE[1])
    for name in ("grid", "edge", "schedules", "tests", "tracked", "refresh"):
        p = sub.add_parser(name)
        p.add_argument("--streams", type=int, default=4)
    pr = sub.add_parser("periods")
    pr.add_argument("--period", type=int, choices=(CALIBRATION, TEST), required=True)
    pr.add_argument("--spec", required=True, help='strength from each year on, e.g. "1969:0.2,2004:0.16,2017:0"')
    sub.add_parser("score")
    t = sub.add_parser("test")
    t.add_argument("--b", type=float, required=True)
    t.add_argument("--f", type=float, required=True)
    t.add_argument("--schedule", default="instant2004", choices=sorted(SCHEDULES))
    t.add_argument("--off", default=None, help="a switch-off year (2017) or a decline Y0-Y1rR (2015-2021r0)")
    a = ap.parse_args()
    globals()["REFERENCE"] = a.reference
    # the 10-09 no-groin-reference runs and the 10-12 afternoon runs are filed under the 10-09 study
    setup = "six-month" if a.reference == "baseline" else a.setup
    globals().update(SETUP=setup, STUDY=STUDIES[setup],
                     LOGS=HERE / "logs" / ("calendar" if setup == "calendar" else ""))
    if a.cmd == "check":
        return cmd_check(a)
    if a.cmd == "smoke":
        sys.exit(0 if run_all([(CALIBRATION, a.b, a.f)], 1) else 1)
    if a.cmd == "grid":
        sys.exit(0 if run_all([(CALIBRATION, b, f) for b in B_VALUES for f in F_VALUES], a.streams) else 1)
    if a.cmd == "edge":
        cells = [(CALIBRATION, b, f) for b in B_EDGE for f in F_EDGE] + [(CALIBRATION, 0.25, f) for f in F_VALUES[1:]]
        sys.exit(0 if run_all(cells, a.streams) else 1)
    if a.cmd == "schedules":
        sys.exit(0 if run_all(schedule_cells(), a.streams) else 1)
    if a.cmd == "periods":
        first = parse_periods(a.spec)[0][1]
        sys.exit(0 if run_all([(a.period, first, 1.0, "periods", a.spec)], 1) else 1)
    if a.cmd == "tests":
        return cmd_tests(a)
    if a.cmd == "tracked":
        return cmd_tracked(a)
    if a.cmd == "refresh":
        return cmd_refresh(a)
    if a.cmd == "test":
        sys.exit(0 if run_all([(TEST, a.b, a.f, a.schedule, a.off)], 1) else 1)
    cmd_score(a)


if __name__ == "__main__":
    main()
