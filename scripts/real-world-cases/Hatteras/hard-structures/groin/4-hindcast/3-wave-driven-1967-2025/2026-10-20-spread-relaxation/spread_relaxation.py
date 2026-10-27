#!/usr/bin/env python3
"""
The groin's step relaxes over three domains: what GIS 6 loses across the groin face each year is shared over
GIS 6, 7 and 8 (1/2, 1/3, 1/6), calibration window 1996-2009 only.

    python spread_relaxation.py run [--streams 4]   # the spread groin over the strength grid
    python spread_relaxation.py score               # -> spread_relaxation_scores.csv, spread_relaxation_domains.csv


Without the groin the model's 1996 step at the groin (GIS 6 about 260 m seaward of GIS 5) relaxes into GIS 6 alone:
GIS 6 loses 71 m where the record has 9 m, and GIS 7-9 lose 0-12 m where it has 14-26 m. The advisor's three-domain
reach is applied to that relaxation. Each year, before BRIE's step, the loss GIS 6 is about to take across the face
is estimated from BRIE's own diffusion term, f = 2 r (x_GIS5 - x_GIS6) with r BRIE's diffusion number at GIS 6, and
GIS 6 is given back half of it while GIS 7 and GIS 8 give up a third and a sixth. The groin's trapping, face and
leak correction are unchanged (`cascade/groin.py` and the runner are not edited); the spread is added to the
groin's own contribution, so the footprint and leak correction include it. Sand is conserved.

Same setup as ../2026-10-20-tapered-fillet/ (runner defaults after ecd6f0ff, dune growth pinned 0.005-0.01), whose
no-groin and runner's-groin runs are the controls here; scoring refuses runs on different settings.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
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
sys.path.insert(0, str(HERE.parent / "2026-10-20-tapered-fillet"))
import tapered_fillet as tf  # noqa: E402  (paths, settings pin, scoring helpers, the control runs)

# --- CONFIG ------------------------------------------------------------------
STUDY = "groin-wave-driven/2026-10-20-spread-relaxation"
RAW = tf.REPO / "output" / "raw_runs" / "experiments" / STUDY
LOGS = workspace_twin(HERE) / "logs"
PERIOD = 1996
STRENGTHS = (0.06, 0.10, 0.13, 0.16, 0.19)
SPREAD = (1 / 2, 1 / 3, 1 / 6)       # share of GIS 6's loss across the face borne by GIS 6, 7, 8 (the advisor's)
FIT_STRENGTHS = (0.10, 0.13, 0.16, 0.19)
# name: (label, shares on GIS 6, 7, 8, mirrored on GIS 5, 4, 3, strengths). "spread2" mirrors the spread onto the
# gain GIS 5 takes across the face; the w* sets are the weight fit of 2026-10-20 (updrift side only)
VARIANTS = {
    "spread": ("spread 1/2, 1/3, 1/6 (GIS 6-8)", SPREAD, False, STRENGTHS),
    "spread2": ("spread both sides (GIS 6-8 and 5-3)", SPREAD, True, STRENGTHS),
    "wB": ("spread 1/3, 1/3, 1/3", (1 / 3, 1 / 3, 1 / 3), False, FIT_STRENGTHS),
    "wC": ("spread 1/3, 1/2, 1/6", (1 / 3, 1 / 2, 1 / 6), False, FIT_STRENGTHS),
    "wD": ("spread 1/4, 1/2, 1/4", (1 / 4, 1 / 2, 1 / 4), False, FIT_STRENGTHS),
    "wE": ("spread 0, 2/3, 1/3", (0.0, 2 / 3, 1 / 3), False, FIT_STRENGTHS),
}
# The test window 2009-2025: the spread groin (0, 2/3, 1/3; run folder wE) only, around its calibration strength (0.10), groin off from 2017 as the runner.
# Its controls are run here ("nogroin", "face" = the runner's groin), not taken from the taper study: the end rates
# were re-solved at 22:39 (c8092db1), after the taper study's 2009 runs
TEST_VARIANTS = {"wE": (0.06, 0.08, 0.10, 0.12, 0.16), "face": (0.08, 0.10, 0.12, 0.16), "nogroin": (0.0,)}
CONTROL_LABELS = {"face": "runner (all in GIS 6)", "nogroin": "no groin"}
# --current: the calibration rerun of 2026-10-21 on the runner exactly as committed (no dune-growth pin; end rates
# and 1996 offset of c8092db1), in its own folders (current-<variant>/) with its own controls
CURRENT_PLAN = [(0.0, "nogroin"), (0.13, "face"), (0.16, "face"), (0.19, "face"),
                (0.08, "wE"), (0.10, "wE"), (0.13, "wE"), (0.16, "wE")]
CURRENT = False                               # set by --current
# -----------------------------------------------------------------------------


# The spread, patched into the GroinManager Cascade builds before the runner runs
def patch_groin(weights=SPREAD, both_sides=False):
    import cascade.cascade as cc
    import cascade.groin as cg
    base = cg.GroinManager

    class SpreadGroin(base):
        def update(self, brie, x_s_dt):
            x = np.asarray(brie.x_s, dtype=float).copy()
            year = self.start_year + self._call_count
            x_s_dt = super().update(brie, x_s_dt)
            lo, hi = self._lo_pad, self._hi_pad
            # GIS 6's coming loss across the face, landward +; only while the groin works and GIS 6 is seaward
            f = 2 * self._diffusion_number(brie, x)[hi] * (x[lo] - x[hi])
            f = f if (self.strength(year) > 0 and f > 0) else 0.0
            # The mirror: GIS 5's coming gain across the face (seaward), from BRIE's diffusion number at GIS 5
            h = 2 * self._diffusion_number(brie, x)[lo] * (x[lo] - x[hi]) if both_sides else 0.0
            h = h if (self.strength(year) > 0 and h > 0) else 0.0
            last = self._last_year
            for arr in (x_s_dt, last.dx_m, last.groin_dx_m):
                arr[hi] -= (1 - weights[0]) * f       # GIS 6 keeps only its share of the loss
                arr[hi + 1] += weights[1] * f         # GIS increases with pad: GIS 7 = hi + 1, GIS 8 = hi + 2
                arr[hi + 2] += weights[2] * f
                arr[lo] += (1 - weights[0]) * h       # GIS 5 keeps only its share of the gain
                arr[lo - 1] -= weights[1] * h         # GIS 4 = lo - 1, GIS 3 = lo - 2
                arr[lo - 2] -= weights[2] * h
            self.record[-1]["spread_loss_m"] = f
            self.record[-1]["spread_gain_m"] = h
            return x_s_dt

    cc.GroinManager = cg.GroinManager = SpreadGroin


def run_dir(s, variant="spread", period=PERIOD):
    folder = f"current-{variant}" if CURRENT else variant     # the run tag allows four path parts
    hits = sorted((RAW / f"{folder}/s{s:.2f}" / tf.WINDOW[period]).glob("*/*/*_shoreline_matrix.npy"))
    return hits[0].parent if hits else None


def _launch():
    import runpy
    import cascade_pipeline.hindcast as hc
    build = hc.build_cascade

    def pinned_build(*a, **kw):
        n = len(kw["rmin"])
        kw["rmin"], kw["rmax"] = [tf.DUNE_GROWTH[0]] * n, [tf.DUNE_GROWTH[1]] * n
        return build(*a, **kw)

    if os.environ.get("SPREAD_PIN_GROWTH", "1") == "1":
        hc.build_cascade = pinned_build
    variant = os.environ["SPREAD_VARIANT"]
    if variant in VARIANTS:                   # the controls run the runner's groin as it is, or none
        _, weights, both_sides, _ = VARIANTS[variant]
        patch_groin(weights, both_sides)
    sys.path.insert(0, str(tf.RUNNER.parent))
    sys.argv = [str(tf.RUNNER)]
    runpy.run_path(str(tf.RUNNER), run_name="__main__")


def run(s, variant, period):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "TAPER", "SPREAD_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
        "HAT_GROIN_ENABLED": "false" if variant == "nogroin" else "true", "HAT_GROIN_STRENGTH": str(s),
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{'current-' if CURRENT else ''}{variant}/s{s:.2f}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8", "SPREAD_VARIANT": variant,
        "SPREAD_PIN_GROWTH": "0" if CURRENT else "1",
    })
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(workspace_twin(LOGS) / f"{tf.WINDOW[period]}_{variant}_s{s:.2f}.log", "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=tf.REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    d = run_dir(s, variant, period)
    if d is not None and variant in VARIANTS:
        (d / "spread_relaxation.txt").write_text(
            f"patched in-process by spread_relaxation.py ({variant}): GIS 6's loss across the groin face shared over "
            f"GIS 6, 7, 8 as {', '.join(f'{w:.3f}' for w in VARIANTS[variant][1])}"
            + (", and GIS 5's gain over GIS 5, 4, 3 the same way" if VARIANTS[variant][2] else "")
            + f", at strength {s:g}. The runner's report shows the plain groin.")
    print(f"{time.strftime('%H:%M:%S')}  {tf.WINDOW[period]} {variant} s{s:.2f}: exit {code} ({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def plan(period):
    """(strength, variant) pairs run in a window."""
    if CURRENT:
        return list(CURRENT_PLAN)
    if period == PERIOD:
        return [(s, v) for v, spec in VARIANTS.items() for s in spec[3]]
    return [(s, v) for v, ss in TEST_VARIANTS.items() for s in ss]


def cmd_run(a):
    every = [(s, v, a.period) for s, v in plan(a.period)]
    todo = [c for c in every if run_dir(*c) is None]
    print(f"{len(every) - len(todo)} of {len(every)} done; running {len(todo)} in {a.streams} streams",
          flush=True)
    jobs, failed = queue.Queue(), []
    for s in todo:
        jobs.put(s)

    def worker():
        while True:
            try:
                s = jobs.get_nowait()
            except queue.Empty:
                return
            if run(*s):
                failed.append(s)

    threads = [threading.Thread(target=worker) for _ in range(a.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}", flush=True)


# Every run scored here, the controls from the taper study first: (label, strength, run folder)
def cells(period=PERIOD):
    if period != PERIOD or CURRENT:           # the test window and the rerun carry their own controls
        out = [(CONTROL_LABELS.get(v) or VARIANTS[v][0], None if v == "nogroin" else s, run_dir(s, v, period))
               for s, v in sorted(plan(period), key=lambda c: (c[1] != "nogroin", c[1] != "face", c[1], c[0]))]
        return [c for c in out if c[2] is not None]
    out = [("no groin", None, tf.run_dir(period, None, None))]
    out += [("runner (all in GIS 6)", s, tf.run_dir(period, "face", s)) for s in tf.STRENGTHS[period]]
    out += [(VARIANTS[v][0], s, run_dir(s, v, period)) for s, v in plan(period)]
    return [c for c in out if c[2] is not None]


# What a run must share with the others to be compared: the Barrier3D settings, the end rates (source/sink digest)
# and the starting shoreline
def run_inputs(rd):
    import hashlib
    import json
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text())
    x0 = np.load(next(rd.glob("*_shoreline_matrix.npy")))[0]
    return tf.run_settings(rd) + (("end rates", meta["source/sink"]["values_digest"]),
                                  ("start", hashlib.md5(x0.round(3).tobytes()).hexdigest()[:8]))


def cmd_score(a):
    period = a.period
    import corrected_groin as cgm
    import site_layer.hat_scoring as common
    from matrix_vs_observed import model_end_change
    from site_layer.hat_observed_rates import net_change_domain_csv
    gap_obs = cgm.observed_change(period)
    d = pd.read_csv(net_change_domain_csv(*tf.NET_KEY[period])).set_index("gis_domain")
    net_obs = d["net_change_lowess7_m"]
    rows, settings = [], set()
    doms = [dict(layout="observed", strength=None, **{f"gis{g}": d.loc[g, "net_change_m"] for g in range(1, 15)})]
    for layout, s, rd in cells(period):
        settings.add(run_inputs(rd))
        path = next(rd.glob("*_shoreline_matrix.npy"))
        raw = model_end_change(rd)
        net = common.smooth_like_target(raw)
        row = dict(layout=layout, strength=s, **cgm.score_annual(path, period, gap_obs))
        row.update(cgm.split_scores(path, period, gap_obs))
        for name, (lo, hi) in (("net_rmse_2_89", (2, 89)), ("net_rmse_2_12", (2, 12)), ("net_rmse_4_9", (4, 9))):
            r = (net - net_obs).loc[lo:hi].dropna()
            row[name] = float(np.sqrt(np.mean(r ** 2)))
        rows.append(row)
        doms.append(dict(layout=layout, strength=s, **{f"gis{g}": raw.loc[g] for g in range(1, 15)}))
    if len(settings) > 1:
        raise RuntimeError(f"the runs used different Barrier3D settings {sorted(settings)}")
    t = pd.DataFrame(rows)
    tag = ("" if period == PERIOD else f"_{tf.WINDOW[period]}") + ("_current" if CURRENT else "")
    t.to_csv(workspace_twin(HERE) / f"spread_relaxation_scores{tag}.csv", index=False, float_format="%.3f")
    pd.DataFrame(doms).to_csv(workspace_twin(HERE) / f"spread_relaxation_domains{tag}.csv", index=False, float_format="%.2f")
    pd.set_option("display.width", 220)
    cols = ["layout", "strength", "rmse_annual", "bias_annual", "rmse_2010_2016", "rmse_2017_2025",
            "net_rmse_2_89", "net_rmse_2_12", "net_rmse_4_9"]
    print(t[[c for c in cols if c in t]].round(2).to_string(index=False))
    print(pd.DataFrame(doms).round(0).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("run")
    sp.add_argument("--streams", type=int, default=4)
    sc = sub.add_parser("score")
    for parser in (sp, sc):
        parser.add_argument("--period", type=int, choices=(1996, 2009), default=PERIOD)
        parser.add_argument("--current", action="store_true", help="the 2026-10-21 rerun on today's runner")
    a = ap.parse_args()
    global CURRENT
    CURRENT = a.current
    {"run": cmd_run, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
