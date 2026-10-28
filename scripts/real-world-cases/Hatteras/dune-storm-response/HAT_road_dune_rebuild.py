"""
Does a gentler, more realistic post-storm dune repair along NC-12 change the shoreline?

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_road_dune_rebuild.py run
    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/dunes/HAT_road_dune_rebuild.py score

RoadwayManager rebuilds a road domain's dunes only when a front-row cell falls 0.3 m above the berm
(about 2.0 m NAVD88), and then resets every cell of the domain to the design crest, 4.3 m NAVD88, so
6 m dunes come down with it. It fired twice in the 1996-2009 managed run, both at GIS 11. The lidar on
the road domains has damaged cells pushed back up only 0.3-0.5 m between surveys, to about 2.3-3.4 m.
Here the rebuild raises cells and never lowers them, to a floor that it also triggers at: every cell
below the floor after the year's storms goes up to it, the rest are untouched. In-process only, nothing
forced from the surveys; managed edgeBE 1996-2009 runs with the groin, as the headline matrix run:

    floor_2p7          floor 2.7 m NAVD88, the lowest design CASCADE allows (berm + 1 m)
    floor_3p0          floor 3.0 m NAVD88
    floor_3p0_anyroad  floor 3.0 m, also where the road sits 20 m or more behind the dune, where
                       RoadwayManager otherwise skips the rebuild

The control is the forced-dune-crests experiment's control run (same code, nothing forced). The run
metadata still prints the setup's 3.94 / 0.0 m MHW rebuild heights; the EXPERIMENT line in each log
gives the ones used.

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
from cascade_pipeline.saved_runs import load_saved_cascade  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-26-road-dune-rebuild"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
CONTROL_DIR = (PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes" / "2026-10-26-forced-dune-crests"
               / "runs" / "control")
HINDCAST = _HAT_DIR / "hatteras_ms" / "HAT_hindcast.py"

START, END, PRESET, SCENARIO = 1996, 2009, "edgeBE", "full_management"
MHW_NAVD88_M = 0.36
# Rebuild floor (m NAVD88) and the road setback beyond which no rebuild happens (m; obx default 20)
VARIANTS = {"floor_2p7": (2.7, 20.0),
            "floor_3p0": (3.0, 20.0),
            "floor_3p0_anyroad": (3.0, 1e6)}
REPORT_GIS = range(5, 21)
WORKERS = 3
# -----------------------------------------------------------------------------


# Child process: raise-only rebuild at the variant's floor, then the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade.roadway_manager as RM
    import cascade_pipeline.hindcast as H
    floor_navd, setback_trigger = VARIANTS[variant]
    floor_mhw = floor_navd - MHW_NAVD88_M
    orig_rebuild = RM.rebuild_dunes

    # Cells below the design come up to it; cells above keep their height
    def rebuild_dunes(yxz_dune_grid, *a, **k):
        new, _ = orig_rebuild(yxz_dune_grid, *a, **k)
        new = np.maximum(yxz_dune_grid, new)
        return new, np.sum(new - yxz_dune_grid)

    RM.rebuild_dunes = rebuild_dunes
    orig_build = H.build_cascade

    def build_cascade(*a, **k):
        n = len(k["dune_design_elevation"])
        k.update(dune_design_elevation=[floor_mhw] * n, dune_minimum_elevation=[floor_mhw] * n,
                 road_dune_rebuild_trigger_m=setback_trigger)
        print(f"\nEXPERIMENT raise-only road dune rebuild: trigger and design {floor_navd:.2f} m NAVD88 "
              f"({floor_mhw:.2f} m MHW), setback trigger {setback_trigger:g} m", flush=True)
        return orig_build(*a, **k)

    H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(variant):
    root = CONTROL_DIR if variant == "control" else EXP_DIR / "runs" / variant
    base = root / f"{START}_{END}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


# One run in its own process, logged, with its launch record
def launch(variant):
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(START),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": SCENARIO, "HAT_GROIN_ENABLED": "1",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{variant}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(variant=variant, floor_m_navd88=VARIANTS[variant][0], setback_trigger_m=VARIANTS[variant][1],
               returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{variant:18s} exit {p.returncode}  {rec['minutes']} min", flush=True)


# Island skill, rebuild counts, and per GIS: LRR against CoastSat, overwash and median crest
def score():
    from cascade_pipeline.coastsat_lowess import build_coastsat_series
    from cascade_pipeline.hindcast import build_target_table
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    sys.path.insert(0, str(_HAT_DIR / "figure_making" / "model_output"))
    from rerender_run_figures import COASTSAT_DATASETS, LOWESS_CONFIG
    from run_figures import TARGET_WINDOW

    cs = build_coastsat_series(COASTSAT_DATASETS, active_period_start=START,
                               lowess_config=LOWESS_CONFIG, domains=DOM)
    target = build_target_table(next(c for c in cs if c["active"]), LOWESS_CONFIG, DOM, TARGET_WINDOW)
    target = target.set_index("gis_domain")["target_lrr_m_yr"]

    frames, summary = {"coastsat": pd.DataFrame({"lrr": target})}, []
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
        rows, years, volume, disabled = [], 0, 0.0, 0
        for g in range(1, 91):
            p = dom.gis_to_pad(g)
            b = c.barrier3d[p]
            crest = (np.asarray(b.DuneDomain, float).max(axis=2) + b._BermEl) * 10 + MHW_NAVD88_M
            rw = c.roadways[p] if p < len(c.roadways) else None
            n_rebuilt = int(np.nansum(rw._dunes_rebuilt_TS)) if rw is not None else 0
            if rw is not None:
                years += n_rebuilt
                volume += float(np.nansum(rw._rebuild_dune_volume_TS))
                disabled += int(np.any(rw._road_dune_rebuild_disabled_TS))
            rows.append(dict(gis=g, lrr=rate.lrr_m_yr[g], rebuilds=n_rebuilt, overwash=ow[1:, p].sum(),
                             crest_p10_2009=np.percentile(crest[-1], 10),
                             crest_median_2009=np.median(crest[-1])))
        frames[v] = pd.DataFrame(rows).set_index("gis")
        sk = meta["skill"]
        lrr = frames[v].lrr.loc[2:89]
        summary.append(dict(run=v, bias=float(sk["mean_bias_interior_m_yr"]), rmse=float(sk["rmse_interior_m_yr"]),
                            rmse_gis_8_15=float(np.sqrt(((lrr.loc[8:15] - target.loc[8:15]) ** 2).mean())),
                            domain_rebuild_years=years, rebuild_volume_m3=round(volume),
                            road_domains_rebuild_skipped=disabled,
                            island_overwash=frames[v].overwash.loc[2:89].sum()))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    summary = pd.DataFrame(summary).set_index("run")
    print(summary.round(3))
    out = pd.concat(frames, axis=1)
    print(out.xs("lrr", axis=1, level=1).loc[list(REPORT_GIS)].round(2))
    print(out.xs("rebuilds", axis=1, level=1).loc[list(REPORT_GIS)])
    out.to_csv(EXP_DIR / "road_dune_rebuild_by_gis.csv")
    summary.to_csv(EXP_DIR / "road_dune_rebuild_summary.csv")
    print(f"\nwrote {EXP_DIR / 'road_dune_rebuild_by_gis.csv'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("action", choices=["run", "score", "_launch"])
    ap.add_argument("variant", nargs="*")
    args = ap.parse_args()
    if args.action == "_launch":
        _launch(args.variant[0])
    elif args.action == "run":
        variants = args.variant or list(VARIANTS)
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, variants))
    else:
        score()


if __name__ == "__main__":
    main()
