#!/usr/bin/env python3
"""
The field layout in the hindcast windows: refit its strength from each observed start, and score it against
the runner's single-face groin.

    python field_windows.py run --share 0.5 [--streams 4]   # the strength grid in both windows
    python field_windows.py score --share 0.5                # -> field_windows_scores.csv

The field layout takes the groin's sand from GIS 5 (south of the field) and gives a share `share` of it to
GIS 6 (the field) and the rest to GIS 7 (north of it); see face_experiment.py. Setup as
../2026-10-14-window-refit/: edgeBE, full management, no relocation, off from 2017, each window from its
observed shoreline. Scored on the gap across the groins (GIS 5 minus GIS 6) against annual CoastSat, as the
window refit, and on net change against the CoastSat net-change target over GIS 2-89 and GIS 2-12 (7-domain
LOWESS, GIS 1-10 raw), as ../2026-10-14-window-refit/matrix_groin_effect_figure.py. The runner's groin is
scored the same way from the window refit's runs.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
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
for p in (_HAT_DIR, _HAT_DIR / "hatteras_ms", _HAT_DIR / "hatteras_ms" / "experiments",
          _HAT_DIR / "analyze_output" / "compare_runs" / "matrix_vs_observed",
          workspace_twin(HERE.parent) / "2026-10-09-corrected-groin"):
    sys.path.insert(0, str(p))

# --- CONFIG ------------------------------------------------------------------
RUNNER = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
STUDY = "groin-wave-driven/2026-10-14-field-windows"
RAW = REPO / "output" / "raw_runs" / "experiments" / STUDY
REFIT = REPO / "output" / "raw_runs" / "experiments" / "groin-wave-driven" / "2026-10-14-window-refit"
LOGS = workspace_twin(HERE) / "logs" / "windows"
WINDOW = {1996: "1996_2009", 2009: "2009_2025"}
STRENGTHS = {1996: (0.10, 0.13, 0.16, 0.19, 0.22), 2009: (0.06, 0.08, 0.10, 0.12, 0.16)}
REFIT_STRENGTHS = {1996: (0.10, 0.13, 0.16, 0.19, 0.22, 0.25), 2009: (0.06, 0.08, 0.10, 0.12, 0.16, 0.25)}
NET_KEY = {1996: (1996, 2009), 2009: (2009, 2025)}
# -----------------------------------------------------------------------------


# The groin's face moved, patched into cascade.groin before the runner imports it. variant: "face67" (one
# domain north), or "field" (GIS 5 to GIS 7 across GIS 6) with `share` of the gain kept in GIS 6
def patch_groin(variant, share=0.0):
    import cascade.groin as cg
    base = cg.GroinCallback

    class MovedGroin(base):
        def __init__(self, updrift_pad, downdrift_pad, *a, **kw):
            if variant == "face67":                      # both pads one domain north (GIS increases with pad)
                updrift_pad, downdrift_pad = updrift_pad + 1, downdrift_pad + 1
            super().__init__(updrift_pad, downdrift_pad, *a, **kw)
            if variant == "field":                       # take from GIS 5, give to GIS 7; GIS 6 between
                self.updrift_pad = self._hi_pad = updrift_pad + 1
                self._mid = updrift_pad

        def __call__(self, cascade, x_s_dt):
            x_s_dt = super().__call__(cascade, x_s_dt)
            if variant == "field" and share:
                # Move a share of this year's updrift gain from GIS 7 into the field, in the step and in the
                # groin's own contribution the tracked footprint is built from
                dx = share * self.record[-1]["applied_dx_updrift_m"]
                last = self._last_year
                for arr in (x_s_dt, last.dx_m, last.groin_dx_m):
                    arr[self._hi_pad] -= dx
                    arr[self._mid] += dx
            return x_s_dt

    cg.GroinCallback = MovedGroin


def member(share, s):
    return f"share{share:.2f}/s{s:.2f}"


def shoreline_file(period, share, s):
    hits = sorted((RAW / member(share, s) / WINDOW[period]).glob("*/*/*_shoreline_matrix.npy"))
    return hits[0] if hits else None


def refit_file(period, s):
    name = "nogroin" if s is None else f"s{s:.2f}"
    hits = sorted((REFIT / name / WINDOW[period]).glob("*/*/*_shoreline_matrix.npy"))
    return hits[0] if hits else None


def _launch():
    import runpy
    patch_groin("field", float(os.environ["FIELD_SHARE"]))
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(period, share, s):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "CORRECTED_GROIN_", "FIELD_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(period),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management", "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "true", "HAT_GROIN_STRENGTH": str(s), "HAT_GROIN_END_YEAR": "2017",
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{member(share, s)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8", "FIELD_SHARE": str(share),
    })
    LOGS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    with open(workspace_twin(LOGS) / f"{WINDOW[period]}_share{share:.2f}_s{s:.2f}.log", "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    f = shoreline_file(period, share, s)
    if f is not None:
        (f.parent / "field_windows.txt").write_text(
            f"patched in-process by field_windows.py: the field layout (sand from GIS 5; {share:g} of it to GIS 6, "
            f"the rest to GIS 7) at strength {s:g}, off from 2017. The runner's report shows the GIS 5|6 groin.")
    print(f"{time.strftime('%H:%M:%S')}  {WINDOW[period]} share{share:.2f} s{s:.2f}: exit {code} "
          f"({(time.time() - t0) / 60:.1f} min)", flush=True)
    return code


def cmd_run(a):
    cells = [(p, a.share, s) for p, ss in STRENGTHS.items() for s in ss]
    todo = [c for c in cells if shoreline_file(*c) is None]
    print(f"{len(cells) - len(todo)} of {len(cells)} done; running {len(todo)} in {a.streams} streams", flush=True)
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

    threads = [threading.Thread(target=worker) for _ in range(a.streams)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    print(f"DONE: {len(todo) - len(failed)}/{len(todo)} clean; failed {failed}", flush=True)


def cmd_score(a):
    import corrected_groin as cgm
    import HAT_metres_1_offset_units as common
    from matrix_vs_observed import model_end_change
    from site_layer.hat_observed_rates import net_change_domain_csv
    rows = []
    for period in WINDOW:
        gap_obs = cgm.observed_change(period)
        d = pd.read_csv(net_change_domain_csv(*NET_KEY[period]))
        net_obs = d.set_index(d.columns[0])["net_change_lowess7_m"]
        cells = [("no groin", None, refit_file(period, None))]
        cells += [("face 5|6", s, refit_file(period, s)) for s in REFIT_STRENGTHS[period]]
        cells += [(f"field, share {a.share:g}", s, shoreline_file(period, a.share, s)) for s in STRENGTHS[period]]
        for layout, s, path in cells:
            if path is None:
                continue
            net = common.smooth_like_target(model_end_change(path.parent))
            row = dict(period=period, layout=layout, strength=s, **cgm.score_annual(path, period, gap_obs))
            row.update(cgm.split_scores(path, period, gap_obs))
            for name, (lo, hi) in (("net_rmse_2_89", (2, 89)), ("net_rmse_2_12", (2, 12))):
                r = (net - net_obs).loc[lo:hi].dropna()
                row[name] = float(np.sqrt(np.mean(r ** 2)))
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(workspace_twin(HERE) / f"field_windows_scores_share{a.share:.2f}.csv", index=False, float_format="%.3f")
    pd.set_option("display.width", 220)
    cols = ["period", "layout", "strength", "rmse_annual", "bias_annual", "rmse_2010_2016", "net_rmse_2_89",
            "net_rmse_2_12"]
    print(t[[c for c in cols if c in t]].round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "score"):
        sp = sub.add_parser(name)
        sp.add_argument("--share", type=float, required=True)
        if name == "run":
            sp.add_argument("--streams", type=int, default=4)
    a = ap.parse_args()
    {"run": cmd_run, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()
