"""
Does a faster dune growth rate at Cape Point alone stop the cape eroding where it accretes?

    python scripts/hatteras_ms/experiments/HAT_cape_growth_rate.py run
    python scripts/hatteras_ms/experiments/HAT_cape_growth_rate.py score

Why (2026-10-21): at the calibrated rates (0.01-0.03; 0.02-0.04 for 2009-2025)
the Cape Point dunes (GIS 2-4) are flattened early and, restarting at 7.5 cm
under logistic growth, never rebuild; every later storm overwashes them
(900-1000 m3/m over 2009-2025) and GIS 1-3 retreat 20-60 m where CoastSat shows
the cape accreting 3-9 m/yr. Not the groin (with and without it GIS 1-3 lose the
same) and not the buffers (the 30-buffer closure leaves GIS 1-3 unchanged). The
surveyed cape dunes rebuild fast (GIS 1-5 median crest 5.0 m in 2014 to 5.6 m
in 2019).

The change: rmin-rmax at GIS CAPE only, per domain, which Barrier3D takes as an
input already (no new process); every other domain keeps the runner's rate. Set
in-process on the build_cascade call. zeroBE (no source/sink, so nothing hides
the cape), full management, groin on, the main code's 30-buffer quintic closure,
both windows. Reference: cape_calibrated, the same setup with the runner's own
rates at the cape (run C of source-sink/2026-10-21-buffer-closure, re-run with
its model state saved). Scores: GIS 1-5 median dune crest against
each survey, overwash at GIS 1-5 (observed cases reproduced, model-only cases),
GIS 1-5 LRR against CoastSat, interior shoreline RMSE and bias.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-21
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(p for p in _HERE.parents if (p / "pyproject.toml").exists())
HINDCAST = PROJECT_ROOT / "scripts" / "hatteras_ms" / "HAT_hindcast.py"

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-21-cape-growth-rate"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
CAPE = (1, 2, 3, 4, 5)
VARIANTS = {"cape_calibrated": None, "cape_r010_020": (0.10, 0.20), "cape_r025_045": (0.25, 0.45), "cape_r055_095": (0.55, 0.95)}
STARTS = (1996, 2009)
PRESET = "zeroBE"
WORKERS = 4
# -----------------------------------------------------------------------------


# Child process: the cape's rates on the build_cascade call, then the unchanged hindcast
def _launch(variant):
    import runpy
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    import cascade_pipeline.hindcast as H
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    if VARIANTS[variant] is None:
        print("EXPERIMENT cape growth: the runner's own rates everywhere (reference)", flush=True)
    else:
        lo, hi = VARIANTS[variant]
        orig = H.build_cascade

        def build_cascade(*a, **k):
            rmin, rmax = list(k["rmin"]), list(k["rmax"])
            for g in CAPE:
                rmin[DOM.gis_to_pad(g)], rmax[DOM.gis_to_pad(g)] = lo, hi
            k["rmin"], k["rmax"] = rmin, rmax
            print(f"EXPERIMENT cape growth: GIS {CAPE[0]}-{CAPE[-1]} rmin/rmax {lo}/{hi}; elsewhere "
                  f"{sorted(set(k['rmin']) - {lo})}/{sorted(set(k['rmax']) - {hi})}", flush=True)
            return orig(*a, **k)

        H.build_cascade = build_cascade
    sys.path.insert(0, str(HINDCAST.parent))
    sys.argv = [str(HINDCAST)]
    runpy.run_path(str(HINDCAST), run_name="__main__")


def run_dir(variant, start, root=None):
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    from site_layer.hatteras_site_config import HATTERAS_PERIODS
    base = (root or EXP_DIR / "runs" / variant) / f"{start}_{HATTERAS_PERIODS[start]['end_year']}" / PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    variant, start = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": PRESET, "HAT_SCENARIO": "full_management", "HAT_GROIN_ENABLED": "True",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{variant}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    t0 = time.time()
    with open(logs / f"{variant}_{start}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", variant],
                           env=env, cwd=HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    print(f"{start} {variant} exit {p.returncode}  {round((time.time() - t0) / 60, 1)} min", flush=True)


# Survey crests (GIS CAPE) against the model, overwash at the cape, cape LRR and interior skill
def score(label, rd, start):
    import numpy as np
    import pandas as pd
    sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
    sys.path.insert(0, str(_HERE.parent))
    sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate"))
    import HAT_dune_ceiling_uniform as U
    import HAT_dune_crest_vs_surveys as V
    import HAT_low_growth_test_period as T
    import be_edge_domain_solve as B
    G = T.G
    G.run_dir = lambda s, v: rd                                      # T.cells and V read the run through G
    c = U.load(rd)
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    first = V.PERIODS[start][1]
    surv = V.surveys(start)
    errs = []
    for year, g in surv[surv.gis.isin(CAPE)].groupby("survey"):
        k = int(g.model_state.iloc[0][:4]) - first
        mod = np.concatenate([np.array(c.barrier3d[DOM.gis_to_pad(gg)].DuneDomain)[k].max(axis=1) * 10 + 1.7
                              for gg in sorted(g.gis.unique())])
        obs = g.dropna(subset=["crest"]).crest
        errs.append(f"{int(year)}:{np.median(mod) - np.median(obs):+.1f}")
    d = T.cells(start, "x")
    d = d[d.gis.isin(CAPE)]
    o, mo = d.observed.astype(int), d.model.astype(int)
    rate = pd.read_csv(next(rd.glob("tables/shoreline_change_rate.csv"))).set_index("gis_domain").lrr_m_yr
    target = B.load_target(start, U.PERIODS[start][0])
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    return dict(window=start, cape_rate=label,
                cape_crest_err_by_survey=" ".join(errs),
                cape_overwash_hits=f"{int(((o == 1) & (mo == 1)).sum())}/{int(o.sum())}",
                cape_model_only=int(((o == 0) & (mo == 1)).sum()),
                cape_lrr=" ".join(f"{rate.loc[g]:+.1f}" for g in CAPE),
                cape_lrr_obs=" ".join(f"{target.get(g, np.nan):+.1f}" for g in CAPE),
                shoreline_rmse=round(float(meta["skill"]["rmse_interior_m_yr"]), 3),
                shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch(sys.argv[2])
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        todo = [(v, s) for v in VARIANTS for s in STARTS if run_dir(v, s) is None]
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, todo))
    import pandas as pd
    rows = []
    for s in STARTS:
        for v, r in VARIANTS.items():
            rd = run_dir(v, s)
            if rd is not None:
                rows.append(score("calibrated" if r is None else f"{r[0]:.2f}-{r[1]:.2f}", rd, s))
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "cape_growth_scores.csv", index=False)
    with pd.option_context("display.width", 320, "display.max_columns", 30, "display.max_colwidth", 80):
        print(t.to_string(index=False))


if __name__ == "__main__":
    main()
