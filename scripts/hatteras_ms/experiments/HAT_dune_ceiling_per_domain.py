r"""
HAT_dune_ceiling_per_domain.py -- does a dune ceiling taken from each domain's own dunes keep the low spots storms break through?
==============================================================================
WHY (storms-and-overwash/2026-09-28-dune-ceiling-and-rebuild): one island-
wide Dmaxel of 5.5 m NAVD88 matches the 2009 lidar on average and fixes
1996-2010 (PSS 0.60), but every dune grows toward it within a few years. That
erases the real low spots. In 2010-2024 the hit rate falls to 0.23, and Irene
overtops 8 of the 60 higher-dune domains against 47 observed.

THE CEILINGS (Hannah, 2026-09-28: "run the per-domain ceiling test"), each
built from the run's OWN starting dunes (DuneDomain[0]: the 1996-survey mosaic
for 1996-2010, the 2009 lidar for 2010-2024), each held at least FLOOR_M above
the berm, since a ceiling at the berm divides by zero in DuneGrowth:
    dom_median   each domain's Dmaxel = its median starting crest
    dom_p25      each domain's Dmaxel = its 25th-percentile starting crest
    cell         each dune CELL's ceiling = that cell's own starting crest
                 (Barrier3d.DuneGrowth wrapped in-process to take an array;
                 the scalar Dmax it returns, used by the flux limiter and
                 CASCADE's growth-rate reset, is the domain median)
    rebuild rule as now (it made no difference at realistic ceilings);
    storms trim24 and drop72; both windows; managed and natural.
Controls: the current model (3.4 m everywhere) and the uniform 5.5 m ceiling,
both from the earlier experiments.

NOTHING IN THE MAIN CODE CHANGES: in each run's process,
cascade_pipeline.hindcast.build_cascade is wrapped to set the ceilings on the
constructed model before the first step; the storm file is swapped as before.
Scoring applies the same DuneGrowth wrapper, so its crest reconstruction
matches the run.

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-dune-ceiling-per-domain/

USAGE
    python HAT_dune_ceiling_per_domain.py run [--workers 6]
    python HAT_dune_ceiling_per_domain.py score
==============================================================================
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
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_length_selection as S  # noqa: E402
import HAT_dune_ceiling_rebuild as E  # noqa: E402

TAG = "storms-and-overwash/2026-09-28-dune-ceiling-per-domain"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
CEILINGS = ("dom_median", "dom_p25", "cell")
STORMS = ("trim24", "drop72")
SCENARIOS = ("full_management", "natural")
FLOOR_M = 0.5                    # m above the berm


def group(storm, ceiling, scenario):
    return f"{storm}_{ceiling}_{scenario}"


def run_dir(w, st, ceiling, sc):
    if ceiling == "uniform3p4":
        return E.run_dir(w, st, 3.4, "current", sc)
    if ceiling == "uniform5p5":
        return E.run_dir(w, st, 5.5, "current", sc)
    base = EXP_DIR / "runs" / group(st, ceiling, sc) / S.wtag(w) / "edgeBE"
    hits = sorted(base.glob("HAT_*")) if base.exists() else []
    return hits[-1] if hits else None


# --- the ceilings --------------------------------------------------------------------

def install_cell_growth():
    """Barrier3d.DuneGrowth taking a per-cell ceiling where the object carries
    one (_hat_cell_dmax, dam above the berm, shape (BarrierLength,)); every
    other object runs the original. Same arithmetic as barrier3d.py."""
    from barrier3d import Barrier3d
    if getattr(Barrier3d, "_hat_cell_growth", False):
        return
    orig = Barrier3d.DuneGrowth

    def DuneGrowth(self, DuneDomain, t):
        cellmax = getattr(self, "_hat_cell_dmax", None)
        if cellmax is None:
            return orig(self, DuneDomain, t)
        Cf, Qdg = 3, 0
        for q in range(self._DuneWidth):
            reduc = 1 / (Cf ** q)
            G = self._growthparam * DuneDomain[t - 1, :, q] * (1 - DuneDomain[t - 1, :, q] / cellmax) * reduc
            DuneDomain[t, :, q] = G + DuneDomain[t - 1, :, q]
            Qdg = Qdg + (np.sum(G) / self._BarrierLength)
        return DuneDomain, float(np.median(cellmax)), Qdg

    Barrier3d.DuneGrowth = DuneGrowth
    Barrier3d._hat_cell_growth = True


def set_ceilings(cascade, ceiling):
    floor = FLOOR_M / 10.0
    for b in cascade.barrier3d:
        crest_h = np.asarray(b.DuneDomain[0]).max(axis=1)            # dam above the berm, per cell
        if ceiling == "cell":
            b._hat_cell_dmax = np.maximum(crest_h, floor)
            h = float(np.median(b._hat_cell_dmax))
        else:
            q = 50 if ceiling == "dom_median" else 25
            h = max(float(np.percentile(crest_h, q)), floor)
        b._Dmaxel = h + b._BermEl                                     # dam MHW, as load_input stores it
        b._Dmax = h


def _launch(start, storm_path, ceiling):
    import cascade_pipeline.hindcast as H
    install_cell_growth()
    orig = H.build_cascade

    def build_cascade(*a, **k):
        c = orig(*a, **k)
        set_ceilings(c, ceiling)
        return c

    H.build_cascade = build_cascade
    sys.argv = [sys.argv[0]]
    S.MD._launch(start, storm_path)


def launch(job):
    w, st, ceiling, sc = job
    g = group(st, ceiling, sc)
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": sc,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{g}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    with open(logs / f"{g}_{S.wtag(w)}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), str(S.storm_file(w, st)), ceiling],
                           env=env, cwd=S.MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    rec = dict(window=S.wtag(w), storm=st, ceiling=ceiling, scenario=sc, returncode=p.returncode,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{S.wtag(w)} {g:40s} exit {p.returncode} {rec['minutes']} min", flush=True)
    return rec


def jobs():
    return [(w, st, c, sc) for w in S.WINDOWS for st in STORMS for c in CEILINGS for sc in SCENARIOS]


def run(workers):
    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(launch, [j for j in jobs() if run_dir(*j) is None]))

# The matrix controls these experiments compared against were archived on
# 2026-09-28 (archive/2026-09-28-loess10-ends/, when the runner's target moved
# to LOWESS-7 and the ends were re-solved in a parallel session). Read from there.
ARCHIVED_MATRIX = PROJECT_ROOT / "output" / "raw_runs" / "archive" / "2026-09-28-loess10-ends" / "matrix"
TARGET_WINDOW = 7            # the runner's target since 2026-09-28


def shoreline_lowess7(rd, w):
    """Interior RMSE and bias against ONE target for every run (CoastSat LRR,
    LOWESS-7, raw for GIS 1-10, as the runner builds it since 2026-09-28), so
    runs made before and after the target change compare on equal terms."""
    from cascade_pipeline.hindcast import build_target_table
    from cascade_pipeline.coastsat_lowess import CoastSatDataset, LowessConfig, build_coastsat_series
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS, SCORE_INTERIOR_GIS
    from site_layer.hat_observed_rates import lrr_csv
    from cascade_pipeline.run_registry import skill_vs_target
    key = tuple(w)
    if key not in _TARGETS:
        cfg = LowessConfig(window_domains=(TARGET_WINDOW,), skip_southern_domains=10)
        series = build_coastsat_series([CoastSatDataset(label=f"CoastSat {w[0]}", period_start=w[0],
                                                        csv_path=str(lrr_csv(*w)))], w[0], cfg,
                                       domains=HATTERAS_DOMAINS)
        _TARGETS[key] = build_target_table(series[0], cfg, HATTERAS_DOMAINS, TARGET_WINDOW)
    lrr = pd.read_csv(rd / "tables" / "shoreline_change_rate.csv").set_index("gis_domain")["lrr_m_yr"]
    model = np.full(HATTERAS_DOMAINS.total_domains, np.nan)
    for g, v in lrr.items():
        model[HATTERAS_DOMAINS.gis_to_pad(int(g))] = v
    sk = skill_vs_target(model, _TARGETS[key], HATTERAS_DOMAINS, interior_gis=SCORE_INTERIOR_GIS)
    return dict(rmse_interior_m_yr=float(sk["rmse_interior_m_yr"]),
                bias_interior_m_yr=float(sk["mean_bias_interior_m_yr"]))


_TARGETS = {}


# --- scores --------------------------------------------------------------------------

def score():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    install_cell_growth()
    S.MATRIX = ARCHIVED_MATRIX                  # the controls, where they now live
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    lid_c = S.load_state(S.MATRIX / "2010_2024/edgeBE" / S.CONTROLS[(2010, "full_management")])
    lidar = E.crest_series(lid_c, pads)[0]
    low = lidar < np.percentile(lidar, 33)
    rows, cells_all = [], []
    for w in S.WINDOWS:
        for st in STORMS:
            for ceiling in ("uniform3p4", "uniform5p5") + CEILINGS:
                for sc in SCENARIOS:
                    rd = run_dir(w, st, ceiling, sc)
                    if rd is None:
                        print(f"  missing {S.wtag(w)} {st} {ceiling} {sc}")
                        continue
                    c = S.load_state(rd)
                    summ = pd.read_csv(S.summary_csv(w, st))
                    cells = S.overwash_cells(c, w, summ, obs, ovm)
                    cells_all.append(cells.assign(window=S.wtag(w), storm=st, ceiling=ceiling, scenario=sc))
                    cs = E.crest_series(c, pads)
                    row = dict(window=S.wtag(w), storm=st, ceiling=ceiling, scenario=sc,
                               crest_start_m=float(np.median(cs[0])), crest_end_m=float(np.median(cs[-1])),
                               be_gis90=float(json.loads(next(rd.glob("*_run_metadata.json")).read_text())
                                              ["index row"]["be_rate_gis90_m_yr"]),
                               **S.overwash_scores(cells, 0.0), **shoreline_lowess7(rd, w))
                    if w[0] == 1996:
                        d = cs[-1] - lidar
                        row.update(crest_2010_minus_lidar_m=float(np.median(d)),
                                   crest_2010_vs_lidar_r=float(np.corrcoef(cs[-1], lidar)[0, 1]))
                    else:
                        ir = cells[cells.obs_id == "OBS-018"].set_index("gis").reindex(range(1, 91))
                        m = (ir.model_m3_per_m > 0).values
                        o = ir.observed.values
                        for name, sel in (("low", low), ("rest", ~low)):
                            k = sel & ~np.isnan(o)
                            row[f"irene_{name}_obs"] = int(np.nansum(o[k]))
                            row[f"irene_{name}_model"] = int(m[k].sum())
                        row["crest_2011_low_third_m"] = float(np.median(cs[1][low]))
                    rows.append(row)
                    print(f"  scored {S.wtag(w)} {st} {ceiling} {sc}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    pd.concat(cells_all).to_csv(EXP_DIR / "tables" / "cells.csv", index=False)
    cols = ["window", "scenario", "storm", "ceiling", "be_gis90", "crest_end_m", "crest_2010_minus_lidar_m",
            "crest_2010_vs_lidar_r", "POD", "POFD", "PSS", "timing_r", "space_r", "rmse_interior_m_yr",
            "irene_low_obs", "irene_low_model", "irene_rest_obs", "irene_rest_model", "crest_2011_low_third_m"]
    with pd.option_context("display.width", 280, "display.max_columns", 30):
        print(t[[c for c in cols if c in t.columns]].round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:5])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "score"])
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    run(a.workers) if a.action == "run" else score()


if __name__ == "__main__":
    main()
