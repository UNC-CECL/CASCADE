#!/usr/bin/env python3
"""
Fit the corrected groin (the drift groin plus a leak correction at the groin's own bends) on 1996-2009.

    python corrected_groin.py check                       # b 0: must reproduce the no-groin run
    python corrected_groin.py smoke  [--b 0.1 --f 0.6]    # one calibration cell
    python corrected_groin.py grid   [--streams 4]        # the calibration grid
    python corrected_groin.py edge   [--streams 4]        # past the grid's f 0.8 edge, and b 0.25
    python corrected_groin.py score                       # -> grid_scores.csv
    python corrected_groin.py test --b 0.1 --f 0.6        # one cell on 2009-2025

The drift groin of ../2026-10-08-drift-groin/ (b of BRIE's wave-climate net drift at the updrift
approach face, moved from GIS 5 to GIS 6 each year), plus straight-coast module 6's leak
correction within HALF_WIDTH cells of the groin face, applied to the groin's bends only (option
B): each year the correction is computed for this run's BRIE shoreline and for the no-groin run's
in the same year, and only the difference is added to x_s_dt,

    L(x) = 2 (r_i - r_{i-1}) (x_i - x_{i-1}),   x_s_dt += L(x_groin) - L(x_none)   near the groin

so with no groin the two cancel and the no-groin baseline is untouched. The no-groin BRIE
shoreline is rebuilt each year as this run's year-0 BRIE shoreline plus the no-groin run's
shoreline change since year 0. BRIE and Barrier3D are not changed: the class replaces
BlockingGroinCallback in a child process, and the unchanged runner is driven as the 10-05 setup.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-09
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
sys.path[:0] = [str(HERE.parent / "2026-10-08-drift-groin"), str(HERE.parent / "2026-10-08-schedule-refit"),
                str(GROIN / "1-observations" / "gap_across_groins"),
                str(GROIN / "2-module-tests" / "3-real-planform"),
                str(REPO / "scripts" / "hatteras_ms" / "groin-sweep"),
                str(REPO / "scripts" / "hatteras_ms"), str(REPO / "scripts")]
import drift_groin as dg  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
RUNNER = dg.RUNNER
STUDY = "groin/2026-10-09-corrected-groin"
RAW = dg.RAW
LOGS = HERE / "logs"
HALF_WIDTH = 10                           # cells either side of the groin face; 5 leaks on a grid-aligned coast
B_VALUES = (0.02, 0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3)   # the drift groin's best sat at b 0.08-0.10
F_VALUES = (0.2, 0.4, 0.6, 0.8)
B_EDGE = (0.15, 0.2, 0.25, 0.3)              # the first grid's best sat at f 0.8, its edge
F_EDGE = (0.9, 1.0)                          # f 1.0: no failure in 2004 (a diagnostic, not a schedule)
SMOKE = (0.1, 0.6)
CALIBRATION, TEST = dg.CALIBRATION, dg.TEST
# -----------------------------------------------------------------------------


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


# Child process: the corrected groin in BlockingGroinCallback's place, then the unchanged runner
def _launch():
    import runpy
    import cascade.groin as cg
    base_m = np.load(os.environ["CORRECTED_GROIN_BASELINE"])     # no-groin run, [year x domain], m, landward +

    class CorrectedGroin(cg.BlockingGroinCallback):
        kind = "blocking"

        def __call__(self, cascade, x_s_dt):
            year = self.start_year + self._call_count
            k = self._call_count
            self._call_count += 1
            active = year >= self.install_year
            b = self._effective_trapping_rate(year) if active else 0.0
            brie = cascade._brie_coupler._brie
            x_s = np.asarray(brie.x_s, dtype=float)
            ny, lo, hi = len(x_s), self._lo, self._hi
            if not hasattr(self, "_q_net"):
                self._q_net = dg.drift_table(brie)
                self._x0 = x_s.copy()
                self.drift_m3_yr_TS, self.approach_deg_TS, self.correction_m_TS, self.baseline_gap_m_TS = [], [], [], []

            def angle(i):
                return 180.0 * np.arctan2(x_s[(i + 1) % ny] - x_s[i], brie._dy) / np.pi

            def q_net(theta):
                return float(self._q_net[int(np.clip(np.round(90 - theta), 1, brie._wave_climl - 1))])

            # The drift groin, as ../2026-10-08-drift-groin/drift_groin.py
            toward_lower = q_net(angle(hi)) < 0
            approach = angle(hi) if toward_lower else angle(lo - 1)
            qv = q_net(approach)
            q = b * qv * brie._dt / (brie._dy * (brie._h_b_crit + brie._d_sf))
            dx_lo, dx_hi = (-q, q) if (active and b > 0.0) else (0.0, 0.0)
            x_s_dt[lo] += dx_lo
            x_s_dt[hi] += dx_hi

            # The leak correction, on the groin's bends only: this run minus the no-groin run
            x_none = self._x0 + (base_m[min(k, len(base_m) - 1)] - base_m[0])
            near = np.abs(np.arange(ny) - (lo + 0.5)) <= HALF_WIDTH
            corr = np.where(near, leak(brie, x_s) - leak(brie, x_none), 0.0)
            x_s_dt += corr

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
                        getattr(self, "correction_m_TS", []), getattr(self, "baseline_gap_m_TS", []))
            for row, (qv, th, c, g) in zip(rows, extra):
                row.update(net_drift_m3_yr=qv, approach_angle_deg=th, leak_correction_m=c,
                           max_diff_from_no_groin_m=g)
            return rows

    cg.BlockingGroinCallback = CorrectedGroin
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(period, b, f):
    from schedule_refit import baseline_file
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(b), "HAT_GROIN_DETERIORATION_FRACTION": str(f),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{dg.member(b, f)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
        "CORRECTED_GROIN_BASELINE": str(baseline_file(period)),
    })
    LOGS.mkdir(exist_ok=True)
    log = LOGS / f"{dg.window(period)}_{dg.member(b, f)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO,
                               env=env, stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(STUDY, period, b, f)
    if code == 0 and d is not None:
        (d / "corrected_groin.txt").write_text(
            "patched in-process by corrected_groin.py: the drift groin (b of BRIE's net drift at the updrift "
            f"approach face, GIS 5 -> GIS 6) plus module 6's leak correction within {HALF_WIDTH} cells, applied "
            "as this run's minus the no-groin run's; the runner's report and groin labels do not show it\n")
    print(f"{time.strftime('%H:%M:%S')}  {dg.window(period)} {dg.member(b, f)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def run_all(cells, streams):
    todo = [c for c in cells if shoreline_file(STUDY, *c) is None]
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
    p = shoreline_file(STUDY, CALIBRATION, 0.0, 0.6)
    if not ok or p is None:
        sys.exit(1)
    x, base = np.load(p), np.load(baseline_file(CALIBRATION))
    diag = pd.read_csv(p.parent / "tables" / "groin_diagnostics.csv")
    print(f"max |shoreline - no-groin| over all years and domains: {np.abs(x - base).max():.3e} m")
    print(f"max |leak correction| applied in any year: {diag.leak_correction_m.abs().max():.3e} m")
    print(f"max |rebuilt no-groin BRIE shoreline - this run's|: {diag.max_diff_from_no_groin_m.max():.3e} m")


def cmd_score(_a):
    from schedule_refit import baseline_file, coastsat_change, score_annual, score_photos
    from score_instant_grid import observed_series
    photos = observed_series()
    rows = []
    for period in (CALIBRATION, TEST):
        obs = coastsat_change(period)
        base = baseline_file(period)
        if base.is_file():
            rows.append(dict(arm="none", period=period, b=0.0, f=np.nan, **score_annual(base, period, obs),
                             **score_photos(base, period, photos), **dg.near_change(base)))
        for arm, study, finder in (("corrected", STUDY, shoreline_file), ("drift", dg.STUDY, dg.shoreline_file)):
            for b in sorted(set(B_VALUES) | set(B_EDGE) | set(dg.B_VALUES) | set(dg.B_FINE)):
                for f in sorted(set(F_VALUES) | set(F_EDGE) | set(dg.F_VALUES)):
                    p = finder(study, period, b, f)
                    if p is not None:
                        rows.append(dict(arm=arm, period=period, b=b, f=f, **score_annual(p, period, obs),
                                         **score_photos(p, period, photos), **dg.budget(p), **dg.near_change(p)))
    t = pd.DataFrame(rows)
    for period in (CALIBRATION, TEST):
        none = t[(t.arm == "none") & (t.period == period)]
        if len(none):
            sel = t.period == period
            t.loc[sel, "island_vs_none_m"] = t.loc[sel, "island_seaward_sum_m"] - float(none.island_seaward_sum_m.iloc[0])
    t["joint_rmse"] = np.sqrt((t.rmse_annual ** 2 + t.rmse_photos ** 2) / 2)
    t.to_csv(HERE / "grid_scores.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    cal = t[(t.period == CALIBRATION) & (t.arm == "corrected")]
    for score in ("rmse_photos", "rmse_annual", "joint_rmse", "island_vs_none_m"):
        if len(cal):
            print(f"\n=== corrected groin, calibration: {score} (m) ===")
            print(cal.pivot(index="b", columns="f", values=score).to_string(float_format="%.1f"))
    best = t[t.period == CALIBRATION].sort_values("joint_rmse").groupby("arm").head(1)
    print(best[["arm", "b", "f", "rmse_annual", "rmse_photos", "joint_rmse", "island_vs_none_m"]].round(1).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    s = sub.add_parser("smoke")
    s.add_argument("--b", type=float, default=SMOKE[0])
    s.add_argument("--f", type=float, default=SMOKE[1])
    g = sub.add_parser("grid")
    g.add_argument("--streams", type=int, default=4)
    e = sub.add_parser("edge")
    e.add_argument("--streams", type=int, default=4)
    sub.add_parser("score")
    t = sub.add_parser("test")
    t.add_argument("--b", type=float, required=True)
    t.add_argument("--f", type=float, required=True)
    a = ap.parse_args()
    if a.cmd == "check":
        return cmd_check(a)
    if a.cmd == "smoke":
        sys.exit(0 if run_all([(CALIBRATION, a.b, a.f)], 1) else 1)
    if a.cmd == "grid":
        sys.exit(0 if run_all([(CALIBRATION, b, f) for b in B_VALUES for f in F_VALUES], a.streams) else 1)
    if a.cmd == "edge":
        cells = [(CALIBRATION, b, f) for b in B_EDGE for f in F_EDGE] + [(CALIBRATION, 0.25, f) for f in F_VALUES[1:]]
        sys.exit(0 if run_all(cells, a.streams) else 1)
    if a.cmd == "test":
        sys.exit(0 if run_all([(TEST, a.b, a.f)], 1) else 1)
    cmd_score(a)


if __name__ == "__main__":
    main()
