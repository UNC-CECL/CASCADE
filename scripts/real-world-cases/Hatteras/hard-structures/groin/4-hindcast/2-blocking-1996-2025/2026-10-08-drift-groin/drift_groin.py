#!/usr/bin/env python3
"""
Fit the drift-blocking groin (option 1) on 1996-2009, through the unchanged runner.

    python drift_groin.py smoke  [--b 0.1 --f 0.6]        # one calibration cell
    python drift_groin.py grid   [--streams 4]            # the calibration grid
    python drift_groin.py fine   [--streams 4]            # the refinement at b 0.06-0.14
    python drift_groin.py score                           # -> grid_scores.csv
    python drift_groin.py test --b 0.1 --f 0.6            # one cell on 2009-2025

Each year the groin reads the shoreline angle of the face just updrift of the GIS 5|6 face,
looks up the wave climate's net longshore drift there from BRIE's own transport (the wave-angle
pdf convolved with _coast_qs, as straight_coast_test.py builds it), and moves a fraction b of
it from the downdrift cell to the updrift cell through x_s_dt:

    q = b * Q_net(theta_approach) * dt / (dy * (h_b_crit + d_sf))      metres of shoreline
    x_s_dt[GIS 5] += q,  x_s_dt[GIS 6] -= q                            (x_s landward +)

so the module's own books balance every year. BRIE and Barrier3D are not changed: the class
replaces BlockingGroinCallback in a child process and the unchanged runner is driven as the
10-05 setup (full management, the solved edgeBE ends, relocations off, instant failure at the
2004 step, f the deteriorated share of b). Scores as the conserving-groin study.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-08
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

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in HERE.parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
GROIN = REPO / "hard-structures" / "groin"
sys.path[:0] = [str(HERE.parent / "2026-10-08-schedule-refit"),
                str(GROIN / "1-observations" / "analysis" / "2-coastsat_gap"),
                str(GROIN / "3-module-tests" / "8-real-planform" / "6-instant-2004-grid"),
                str(_HAT_DIR / "hatteras_ms" / "groin-sweep"),
                str(_HAT_DIR / "hatteras_ms"), str(_HAT_DIR)]

# --- CONFIG ------------------------------------------------------------------
RUNNER = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
STUDY = "groin-blocking/2026-10-08-drift-groin"
PIN_STUDY = "groin-blocking/2026-10-05-blocking-fit-dem-to-dem"
RAW = REPO / "output" / "raw_runs"
LOGS = workspace_twin(HERE) / "logs"
B_VALUES = (0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6)   # b 0.6 moves ~56 m/yr at Buxton, ~10x the pin
F_VALUES = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
B_FINE = (0.06, 0.07, 0.08, 0.09, 0.11, 0.12, 0.13, 0.14)   # the steep band of the first grid
F_FINE = (0.2, 0.4, 0.6, 0.8)
SMOKE = (0.1, 0.6)
CALIBRATION, TEST = 1996, 2009
UP, DOWN = 15 + 6 - 1, 15 + 5 - 1      # GIS 6 / GIS 5 in the padded array
NEAR = {"GIS 3": 15 + 3 - 1, "GIS 4": 15 + 4 - 1, "GIS 7": 15 + 7 - 1, "GIS 8": 15 + 8 - 1}
# -----------------------------------------------------------------------------


def member(b, f):
    return f"b{b:.2f}_f{f:.1f}"


def window(period):
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    return f"{period}_{HATTERAS_PERIODS[period]['end_year']}"


def run_dir(study, period, b, f):
    hits = sorted((RAW / "experiments" / study / member(b, f)).glob(f"{window(period)}/edgeBE/*/"))
    return hits[-1] if hits else None


def shoreline_file(study, period, b, f):
    d = run_dir(study, period, b, f)
    hits = sorted(d.glob("*_shoreline_matrix.npy")) if d else []
    return hits[-1] if hits else None


# Net drift through each angle bin toward higher domain index, m3/yr, from BRIE's own tables
def drift_table(brie):
    step = 180.0 / (brie._wave_climl - 1)
    pdf = brie._angles.pdf(brie._angle_array) * np.deg2rad(step)
    conv = np.convolve(pdf, brie._coast_qs, mode="full")
    npad = len(brie._coast_qs) - 1
    first = npad - npad // 2
    return -conv[first:first + len(pdf)]


# Child process: the drift-blocking groin in BlockingGroinCallback's place, then the unchanged runner
def _launch():
    import runpy
    import cascade.groin as cg

    class DriftGroin(cg.BlockingGroinCallback):
        kind = "blocking"

        def __call__(self, cascade, x_s_dt):
            year = self.start_year + self._call_count
            self._call_count += 1
            active = year >= self.install_year
            b = self._effective_trapping_rate(year) if active else 0.0
            brie = cascade._brie_coupler._brie
            if not hasattr(self, "_q_net"):
                self._q_net = drift_table(brie)
                self.drift_m3_yr_TS, self.approach_deg_TS = [], []
            x_s = np.asarray(brie.x_s, dtype=float)
            ny, lo, hi = len(x_s), self._lo, self._hi

            def angle(i):                    # the face i -> i+1, as brie.py computes it
                return 180.0 * np.arctan2(x_s[(i + 1) % ny] - x_s[i], brie._dy) / np.pi

            def q_net(theta):
                k = int(np.clip(np.round(90 - theta), 1, brie._wave_climl - 1))
                return float(self._q_net[k])

            toward_lower = q_net(angle(hi)) < 0          # drift toward GIS 5: the face above GIS 6 feeds it
            approach = angle(hi) if toward_lower else angle(lo - 1)
            qv = q_net(approach)
            q = b * qv * brie._dt / (brie._dy * (brie._h_b_crit + brie._d_sf))
            dx_lo, dx_hi = (-q, q) if (active and b > 0.0) else (0.0, 0.0)
            x_s_dt[lo] += dx_lo
            x_s_dt[hi] += dx_hi

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
            return x_s_dt

        def diagnostics_frame(self):
            rows = super().diagnostics_frame()
            for row, qv, th in zip(rows, getattr(self, "drift_m3_yr_TS", []), getattr(self, "approach_deg_TS", [])):
                row.update(net_drift_m3_yr=qv, approach_angle_deg=th)
            return rows

    cg.BlockingGroinCallback = DriftGroin
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(period, b, f):
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_RELOCATIONS": "false", "HAT_GROIN_ENABLED": "true", "HAT_GROIN_KIND": "blocking",
        "HAT_GROIN_BLOCKING_FRACTION": str(b), "HAT_GROIN_DETERIORATION_FRACTION": str(f),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{member(b, f)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8",
    })
    LOGS.mkdir(exist_ok=True)
    log = workspace_twin(LOGS) / f"{window(period)}_{member(b, f)}.log"
    t0 = time.time()
    with open(log, "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO,
                               env=env, stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(STUDY, period, b, f)
    if code == 0 and d is not None:
        (d / "drift_groin.txt").write_text(
            "patched in-process by drift_groin.py: the groin moves b of BRIE's wave-climate net drift at the "
            "updrift approach face across GIS 5|6 each year (not BlockingGroinCallback's exchange); the "
            "runner's report and groin labels do not show it\n")
    print(f"{time.strftime('%H:%M:%S')}  {window(period)} {member(b, f)}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def run_all(cells, streams):
    todo = [c for c in cells if shoreline_file(STUDY, *c) is None]
    print(f"{len(cells) - len(todo)} of {len(cells)} done; running {len(todo)} in {streams} streams",
          flush=True)
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


# The groin's applied shoreline changes, from the run's own diagnostics
def budget(path):
    d = pd.read_csv(path.parent / "tables" / "groin_diagnostics.csv")
    a = d[d.groin_active]
    up, down = a.applied_dx_updrift_m.sum(), a.applied_dx_downdrift_m.sum()
    return dict(applied_up_m=up, applied_down_m=down, applied_net_m=up + down)


# Seaward change (x_s landward +, so start minus end) near the groin, and summed over the island
def near_change(path):
    x = np.load(path)
    out = {f"{k}_seaward_m": float(x[0, i] - x[-1, i]) for k, i in NEAR.items()}
    out["island_seaward_sum_m"] = float((x[0] - x[-1]).sum())
    return out


def cmd_score(_a):
    from schedule_refit import baseline_file, coastsat_change, score_annual, score_photos
    from score_instant_grid import observed_series
    photos = observed_series()
    rows = []
    for period in (CALIBRATION, TEST):
        obs = coastsat_change(period)
        base = baseline_file(period)
        if base.is_file():
            rows.append(dict(arm="none", period=period, b=0.0, f=np.nan,
                             **score_annual(base, period, obs), **score_photos(base, period, photos),
                             **near_change(base)))
        for arm, study in (("drift", STUDY), ("pinned", PIN_STUDY)):
            for b in sorted(set(B_VALUES) | set(B_FINE) | {0.6}):
                for f in F_VALUES:
                    p = shoreline_file(study, period, b, f)
                    if p is not None:
                        rows.append(dict(arm=arm, period=period, b=b, f=f,
                                         **score_annual(p, period, obs),
                                         **score_photos(p, period, photos),
                                         **budget(p), **near_change(p)))
    t = pd.DataFrame(rows)
    t.to_csv(workspace_twin(HERE) / "grid_scores.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    none = t[(t.arm == "none") & (t.period == CALIBRATION)]
    if len(none):
        t["island_vs_none_m"] = t.island_seaward_sum_m - float(none.island_seaward_sum_m.iloc[0])
    cal = t[(t.period == CALIBRATION) & (t.arm == "drift")]
    for score in ("rmse_photos", "rmse_annual"):
        if len(cal):
            print(f"\n=== drift groin, calibration: {score} (m) ===")
            print(cal.pivot(index="b", columns="f", values=score).to_string(float_format="%.1f"))
    print(t[t.period == CALIBRATION].drop(columns=[c for c in t.columns if c.startswith("GIS")])
          .round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("smoke")
    s.add_argument("--b", type=float, default=SMOKE[0])
    s.add_argument("--f", type=float, default=SMOKE[1])
    g = sub.add_parser("grid")
    g.add_argument("--streams", type=int, default=4)
    fn = sub.add_parser("fine")
    fn.add_argument("--streams", type=int, default=4)
    sub.add_parser("score")
    t = sub.add_parser("test")
    t.add_argument("--b", type=float, required=True)
    t.add_argument("--f", type=float, required=True)
    a = ap.parse_args()
    if a.cmd == "smoke":
        sys.exit(0 if run_all([(CALIBRATION, a.b, a.f)], 1) else 1)
    if a.cmd == "grid":
        cells = [(CALIBRATION, b, f) for b in B_VALUES for f in F_VALUES]
        sys.exit(0 if run_all(cells, a.streams) else 1)
    if a.cmd == "fine":
        cells = [(CALIBRATION, b, f) for b in B_FINE for f in F_FINE]
        sys.exit(0 if run_all(cells, a.streams) else 1)
    if a.cmd == "test":
        sys.exit(0 if run_all([(TEST, a.b, a.f)], 1) else 1)
    cmd_score(a)


if __name__ == "__main__":
    main()
