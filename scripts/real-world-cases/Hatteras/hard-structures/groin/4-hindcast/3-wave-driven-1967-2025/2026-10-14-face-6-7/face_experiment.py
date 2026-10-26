#!/usr/bin/env python3
"""
Where should the groin's face be? The 1967-1995 build-up with the face at 5|6, at 6|7, or across domain 6.

    python face_experiment.py run   [--streams 4]   # every variant and strength, from the 1967 shoreline
    python face_experiment.py check                 # the 5|6 control against the 10-13 full-life runs
    python face_experiment.py score                 # -> face_scores.csv, face_shape.csv
    python face_experiment.py figure                # -> figures/face_shape_1978_1993.png

All three groins sit inside GIS 6 (transects 37, 39 and 41: its south edge, middle and north edge), and
the photos put the updrift build-out at GIS 7 (../2026-10-14-fillet-shape/). The runner's groin is one
face between adjacent domains, GIS 5|6. Variants, patched in-process (cascade/groin.py is not edited):
  face56  the runner's groin as it is (the control)
  face67  the same groin one domain north: from GIS 6 to GIS 7
  field   domain 6 is the field: from GIS 5 (south of it) to GIS 7 (north of it), GIS 6 left to BRIE
  split50, split33  the field, with half or a third of the sand held in GIS 6 itself and the rest in GIS 7
Setup as ../2026-10-13-pregroin-planform/ "full life": the 1967 wet/dry shoreline at GIS 2-12, 1968-1995,
storms 1968-1983 replaying 1984-1999, RSLR 0.004, no fills, the 1996 topography; edgeBE, full management.

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
PREGROIN = HERE.parent / "2026-10-13-pregroin-planform"
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(PREGROIN))
sys.path.insert(0, str(HERE.parent / "2026-10-09-corrected-groin"))

# --- CONFIG ------------------------------------------------------------------
RUNNER = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"
STUDY = "groin-wave-driven/2026-10-14-face-6-7"
RAW = REPO / "output" / "raw_runs" / "experiments" / STUDY
OLD_FULL = REPO / "output" / "raw_runs" / "experiments" / "groin-wave-driven" / "2026-10-13-pregroin-full-life"
LOGS = workspace_twin(HERE) / "logs"
PERIOD = 1996                  # the calibration run's setup, its years swapped as in pregroin_planform
FIRST, LAST = 1968, 1995
START = 1967                   # row k of a run's shoreline matrix is the end of 1967 + k
RSLR = 0.004
CELLS = {"face56": (0.0, 0.16, 0.20, 0.25), "face67": (0.16, 0.25, 0.35), "field": (0.10, 0.16, 0.20, 0.25, 0.35),
         "split50": (0.16, 0.20, 0.25, 0.30, 0.35), "split33": (0.20, 0.25, 0.30)}
SPLIT = {"split50": 0.5, "split33": 1 / 3}   # share of the trapped sand kept in the field (GIS 6)
OBSERVED = (REPO / "hard-structures" / "groin" / "1-observations" / "analysis" / "1-aerial_photo_gap"
            / "Change_from_wetdry_1967_D2_D12.csv")      # landward +, from the 1967 wet/dry line
GIS = np.arange(2, 13)
DOWN, UP = (3, 4, 5), (6, 7, 8, 9)
FIG_YEARS = (1978, 1993)
# -----------------------------------------------------------------------------


def member(variant, s):
    return f"{variant}/s{s:.2f}"


def shoreline_file(variant, s):
    hits = sorted((RAW / member(variant, s)).glob("*/*/*/*_shoreline_matrix.npy"))
    return hits[0] if hits else None


# Child: the period's offset, years, storms, sea level and fills swapped (as pregroin_planform's full life),
# the groin's face moved for the variant, then the unchanged runner
def _launch():
    import runpy
    import pregroin_planform as pp
    from site_layer import hatteras_site_config as hsc
    p = hsc.HATTERAS_PERIODS[PERIOD]
    storms = os.path.relpath(pp.STORMS_OUT, REPO / "data" / "hatteras_init")
    p.update(island_offset_file=str(pp.OFFSET_OUT), first_model_year=FIRST, last_model_year=LAST,
             storm_file=storms, sea_level_rise_rate=RSLR, enable_nourishment=False)
    strict = hsc.last_model_year
    hsc.last_model_year = lambda y: LAST if int(y) == PERIOD else strict(y)

    from field_windows import patch_groin
    variant = os.environ["FACE_VARIANT"]
    if variant in SPLIT:
        patch_groin("field", SPLIT[variant])
    elif variant != "face56":
        patch_groin(variant)
    sys.path.insert(0, str(RUNNER.parent))
    sys.argv = [str(RUNNER)]
    runpy.run_path(str(RUNNER), run_name="__main__")


def run(variant, s):
    env = {k: v for k, v in os.environ.items() if not k.startswith(("HAT_", "CORRECTED_GROIN_", "FACE_"))}
    env.update({
        "HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(PERIOD),
        "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management", "HAT_RELOCATIONS": "false",
        "HAT_GROIN_ENABLED": "false" if s == 0 else "true",
        "HAT_GROIN_STRENGTH": str(s), "HAT_GROIN_END_YEAR": "2017",
        "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{STUDY}/{member(variant, s)}",
        "HAT_SAVE_MODEL_STATE": "false", "HAT_SHOW_FIGURES": "false", "HAT_MAKE_GIFS": "false",
        "HAT_OVERWRITE": "false", "MPLBACKEND": "Agg", "PYTHONIOENCODING": "utf-8", "FACE_VARIANT": variant,
    })
    LOGS.mkdir(exist_ok=True)
    t0 = time.time()
    with open(workspace_twin(LOGS) / f"{variant}_s{s:.2f}.log", "w", encoding="utf-8") as fh:
        code = subprocess.call([sys.executable, str(Path(__file__)), "_launch"], cwd=REPO, env=env,
                               stdout=fh, stderr=subprocess.STDOUT)
    f = shoreline_file(variant, s)
    if f is not None:
        (f.parent / "face_experiment.txt").write_text(
            f"patched in-process by face_experiment.py: variant {variant} at strength {s:g} from 1969; years "
            f"{FIRST}-{LAST} from the 1967 wet/dry shoreline at GIS 2-12, storms 1968-1983 replaying 1984-1999, "
            f"RSLR {RSLR}, no fills, the 1996 topography. The runner's report shows the 1996 period and the "
            f"GIS 5|6 groin.")
    print(f"{time.strftime('%H:%M:%S')}  {variant} s{s:.2f}: exit {code} ({(time.time() - t0) / 60:.1f} min)",
          flush=True)
    return code


def cmd_run(a):
    cells = [(v, s) for v, ss in CELLS.items() for s in ss]
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


# The control against the 10-13 full-life runs (made through the old in-process patch, commit 4d937b79)
def cmd_check(_a):
    for s in CELLS["face56"]:
        old = next(OLD_FULL.glob(f"periods_1969-{s:g}/*/*/*/*/*_shoreline_matrix.npy"))
        new = shoreline_file("face56", s)
        print(f"face56 s{s:.2f} vs 10-13 run: max |diff| {np.abs(np.load(new) - np.load(old)).max():.3f} m")


# Seaward + change since 1967 at GIS 2-12, the straight line between the end covered domains taken off
def shape(landward):
    v = -pd.Series(landward, index=GIS, dtype=float)
    ok = v.dropna()
    return v - np.interp(GIS, [ok.index[0], ok.index[-1]], [ok.iloc[0], ok.iloc[-1]])


def observed():
    t = pd.read_csv(OBSERVED).set_index("Domain_ID")
    out = {}
    for c in t.columns:
        if not c.startswith("change_from_wetdry_1967_wetdry_") or c.endswith("1967_m"):
            continue
        y = int(c.split("_")[5][:4])
        if START < y <= LAST and c.count("_") == 6 and t[c].reindex(GIS).notna().sum() >= 8:
            out[y] = shape(t[c].reindex(GIS).to_numpy())
    return out


def modelled(variant, s, year):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS
    x = np.load(shoreline_file(variant, s))
    pads = HATTERAS_DOMAINS.gis_to_pad(GIS)
    return shape(x[year - START, pads] - x[0, pads])


def cmd_score(_a):
    obs = observed()
    rows, shapes = [], []
    for v, ss in CELLS.items():
        for s in ss:
            if shoreline_file(v, s) is None:
                continue
            d = {g: [] for g in ("all", "down", "up")}
            peaks = []
            for y, o in obs.items():
                m = modelled(v, s, y)
                r = (m - o).dropna()
                d["all"] += list(r.values)
                d["down"] += list(r.loc[list(DOWN)].dropna().values)
                d["up"] += list(r.loc[list(UP)].dropna().values)
                peaks.append(int(m.loc[list(UP)].idxmax()))
                shapes.append(dict(variant=v, strength=s, year=y, **{f"gis{g}": m[g] for g in GIS}))
            rows.append(dict(variant=v, strength=s, n_dates=len(obs),
                             rmse_shape=np.sqrt(np.mean(np.square(d["all"]))),
                             rmse_downdrift_3_5=np.sqrt(np.mean(np.square(d["down"]))),
                             rmse_updrift_6_9=np.sqrt(np.mean(np.square(d["up"]))),
                             updrift_peak_gis=pd.Series(peaks).mode().iloc[0]))
    for y, o in obs.items():
        shapes.append(dict(variant="observed", strength=np.nan, year=y, **{f"gis{g}": o[g] for g in GIS}))
    t = pd.DataFrame(rows)
    t.round(2).to_csv(workspace_twin(HERE) / "face_scores.csv", index=False)
    pd.DataFrame(shapes).round({f"gis{g}": 1 for g in GIS}).to_csv(workspace_twin(HERE) / "face_shape.csv", index=False)
    print("photo dates:", " ".join(map(str, obs)))
    print("observed updrift peak (mode over dates):",
          pd.Series([int(o.loc[list(UP)].idxmax()) for o in obs.values()]).mode().iloc[0])
    print(t.round(1).to_string(index=False))


def cmd_figure(_a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from site_layer import hat_figure_style as fs
    t = pd.read_csv(workspace_twin(HERE) / "face_scores.csv")
    best = {v: float(t[t.variant == v].sort_values("rmse_shape").strength.iloc[0])
            for v in ("face67", "field", "split50")}
    styles = {("face56", 0.0): ("no groin", dict(color="#7f7f7f", lw=1.3)),
              ("face56", 0.25): ("face 5|6 (runner), 0.25", dict(color="#7b3294", lw=1.5)),
              ("face67", best["face67"]): (f"face 6|7, {best['face67']:.2f}", dict(color="#1b7f8c", lw=1.5)),
              ("field", best["field"]): (f"field 5 to 7, {best['field']:.2f}", dict(color="#c8880f", lw=1.5)),
              ("split50", best["split50"]): (f"field, half kept in GIS 6, {best['split50']:.2f}",
                                             dict(color="#2c6e49", lw=1.5))}
    obs = observed()
    fs.apply_style()
    fig, axes = plt.subplots(1, len(FIG_YEARS), figsize=fs.figsize("double", height=3.4), sharey=True,
                             constrained_layout=True)
    for ax, y in zip(axes, FIG_YEARS):
        ax.axhline(0, color=fs.INK_MUTED, lw=0.6)
        ax.axvspan(5.5, 6.5, color="0.92", zorder=0)
        ax.text(6, 0.97, "groin\nfield", transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6.5, color=fs.C["GROIN"])
        for (v, s), (label, st) in styles.items():
            if shoreline_file(v, s) is not None:
                ax.plot(GIS, modelled(v, s, y).values, "o-", ms=2.5, label=label, **st)
        ax.plot(GIS, obs[y].values, "o-", ms=3, color="black", lw=1.8, label="observed (wet/dry photo)", zorder=5)
        ax.set_title(f"1967 to {y}", fontsize=8.5, loc="left")
        ax.set_xlim(1.5, 12.5)
        ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    axes[0].set_ylabel("Shoreline change since 1967 (m), seaward +\n(end-to-end line removed)")
    axes[0].legend(loc="lower right", frameon=False, fontsize=6.3)
    fig.suptitle("Where the groin's face sits: the 1967-1995 build-up, observed against each placement",
                 fontsize=9)
    print(*fs.save(fig, workspace_twin(HERE) / "figures" / "face_shape_1978_1993", dpi=300, close=True), sep="\n")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--streams", type=int, default=4)
    for name in ("check", "score", "figure"):
        sub.add_parser(name)
    a = ap.parse_args()
    {"run": cmd_run, "check": cmd_check, "score": cmd_score, "figure": cmd_figure}[a.cmd](a)


if __name__ == "__main__":
    main()
