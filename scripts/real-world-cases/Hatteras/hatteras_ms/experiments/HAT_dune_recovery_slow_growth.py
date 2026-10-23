"""
Does a recovery term for flattened dunes remove the chronic overwash that slow dune growth leaves at Cape Point?

    python scripts/hatteras_ms/experiments/HAT_dune_recovery_slow_growth.py run
    python scripts/hatteras_ms/experiments/HAT_dune_recovery_slow_growth.py score

Why (2026-10-20): with Dmax 10.5 m above the berm and rmin-rmax 0.01-0.03, a dune
a storm flattens restarts at 7.5 cm and logistic growth, proportional to its
height, takes decades to rebuild it. GIS 2-4 stay flat through 2009-2025, are
overwashed by every storm (900-1000 m3/m against 0-374 with the previous dunes)
and retreat where CoastSat shows Cape Point accreting
(dunes/2026-10-20-shoreline-rescore-new-ends).

The change, in-process only (Barrier3d.DuneGrowth replaced in the run's own
process; nothing in Barrier3D or CASCADE is edited): the current DuneGrowth plus
A * max(0, 1 - D / D_REC) per year, A to the front dune row and A/3 to the back
(Barrier3D's Cf = 3 decay). It tapers to zero at the recovery height D_REC, the
NC-12 dune design crest (4.3 m NAVD88 = 2.6 m above the 1.7 m berm), not at
Dmax: a taper to Dmax 10.5 m would add nearly A to every dune and undo the
stable-dune calibration, whose crests (6-7 m NAVD88) sit above D_REC and so
keep the calibrated growth alone. The added sand counts in Qdg, which Barrier3D
takes from the beach and shoreface, as dune growth does. (The 2026-10-01
recovery test tapered to the ceiling, under the old fast growth.)

Runs: the adopted dunes (Dmaxel 12.2 m, 0.01-0.03) with A = 0.15, 0.30, 0.50
m/yr, groin on, both windows, today's setup and ends, through the grid study's
build; the A = 0 reference is the same setup's run in
dunes/2026-10-20-shoreline-rescore-new-ends. Scores: interior shoreline RMSE
(model unsmoothed and LOWESS-7) and bias, overwash hits / model-only cases, GIS
2-4 LRR against CoastSat, flattened dune cell-years, and, 1996 only, the
stable-dune 2010 crest error (GIS 59-62, 20-22).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-20
"""
from __future__ import annotations

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
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402
import HAT_dune_ceiling_growth_grid as G  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "dunes/2026-10-20-dune-recovery-slow-growth"
STUDY = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
REFERENCE = U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "dunes" / "2026-10-20-shoreline-rescore-new-ends"
BERM_NAVD88 = 1.7
D_REC_M = 4.3 - BERM_NAVD88          # recovery height above the berm: the NC-12 design crest, m
RECOVERY = {"rec015": 0.15, "rec030": 0.30, "rec050": 0.50}   # A, m/yr to the front dune row
VARIANT = G.name(12.2, (0.01, 0.03), groin=True)               # the adopted dunes
STARTS = (1996, 2009)
FLAT_M = 0.5                         # a cell this low above the berm counts as flattened
CAPE = (2, 3, 4)
WORKERS = 4
# -----------------------------------------------------------------------------


# Barrier3d.DuneGrowth as it stands, plus the recovery term (dam)
def install_recovery(a_m, d_rec_m):
    from barrier3d import Barrier3d

    def DuneGrowth(self, DuneDomain, t):
        Dmax = self._Dmaxel - self._BermEl
        if Dmax < 0:
            Dmax = 0
        Cf, Qdg = 3, 0
        for q in range(self._DuneWidth):
            reduc = 1 / (Cf ** q)
            d = DuneDomain[t - 1, :, q]
            G_ = self._growthparam * d * (1 - d / Dmax) * reduc
            G_ = G_ + (a_m / 10.0) * np.clip(1 - d / (d_rec_m / 10.0), 0, None) * reduc
            DuneDomain[t, :, q] = G_ + d
            Qdg = Qdg + (np.sum(G_) / self._BarrierLength)
        return DuneDomain, Dmax, Qdg

    Barrier3d.DuneGrowth = DuneGrowth
    print(f"EXPERIMENT dune recovery A {a_m} m/yr, tapering to zero at {d_rec_m:.2f} m above the berm", flush=True)


def _launch(rec):
    install_recovery(RECOVERY[rec], D_REC_M)
    G._launch(VARIANT)                    # sets the ceiling and rates, then runs the unchanged hindcast


# Runs sit at <study>/<rec>/<variant>/<window>/<preset>: the runner allows four tag components, so the
# grid layout's runs/ level is dropped; run_dir is swapped for that layout while a recovery run is read
_GRID_RUN_DIR = G.run_dir


def _rec_run_dir(start, variant):
    base = G.EXP_DIR / variant / f"{start}_{U.PERIODS[start][0]}" / U.PRESET
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def point_at(rec):
    G.TAG, G.EXP_DIR = f"{TAG}/{rec}", STUDY / rec
    G.run_dir = _rec_run_dir


def point_at_reference():
    G.TAG, G.EXP_DIR = "dunes/2026-10-20-shoreline-rescore-new-ends", REFERENCE
    G.run_dir = _GRID_RUN_DIR


def launch(job):
    start, rec = job
    logs = STUDY / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("HAT_")}
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(start),
                "HAT_SOURCE_SINK_PRESET": U.PRESET, "HAT_SCENARIO": U.SCENARIO, "HAT_GROIN_ENABLED": "True",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/{rec}/{VARIANT}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8", "MPLBACKEND": "Agg"})
    t0 = time.time()
    with open(logs / f"{rec}_{start}.log", "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", rec],
                           env=env, cwd=U.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    r = dict(start=start, recovery=rec, returncode=p.returncode, minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(r) + "\n")
    print(f"{start} {rec} exit {p.returncode}  {r['minutes']} min", flush=True)


def jobs():
    todo = []
    for rec in RECOVERY:
        point_at(rec)
        todo += [(s, rec) for s in STARTS if G.run_dir(s, VARIANT) is None]
    return todo


# One run's scores; G must point at the run's folder
def score(start):
    tag, exp_dir, rd_fn = G.TAG, G.EXP_DIR, G.run_dir   # before the imports: each re-points the launcher
    import HAT_low_growth_test_period as T
    import HAT_stable_dune_growth_rate as R
    from statsmodels.nonparametric.smoothers_lowess import lowess
    G.TAG, G.EXP_DIR, G.run_dir = tag, exp_dir, rd_fn
    rd = G.run_dir(start, VARIANT)
    meta = json.loads(next(rd.glob("*_run_metadata.json")).read_text(encoding="utf-8"))
    rate = pd.read_csv(next(rd.glob("tables/shoreline_change_rate.csv"))).set_index("gis_domain").lrr_m_yr
    sys.path.insert(0, str(U.PROJECT_ROOT / "scripts" / "input_prep" / "7-source-sink" / "2-calibrate"))
    import be_edge_domain_solve as B
    target = B.load_target(start, U.PERIODS[start][0])
    gis = rate.index.to_numpy()
    obs = np.array([target.get(int(g), np.nan) for g in gis])
    m = rate.to_numpy()
    sm = lowess(m, gis.astype(float), frac=7 / len(gis), return_sorted=False)
    keep = (gis <= 10) | (gis == gis.max())
    sm[keep] = m[keep]
    inner = (gis >= 2) & (gis <= 89) & np.isfinite(obs)
    d = T.cells(start, VARIANT)
    o, mo = d.observed.astype(int), d.model.astype(int)
    c = U.load(rd)
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    flat = np.mean([np.mean(np.array(c.barrier3d[DOM.gis_to_pad(g)].DuneDomain)[1:].max(axis=2) * 10 < FLAT_M)
                    for g in range(1, 91)])
    row = dict(window=f"{start}_{U.PERIODS[start][0]}",
               shoreline_rmse=round(float(np.sqrt(np.mean((m[inner] - obs[inner]) ** 2))), 3),
               shoreline_rmse_smoothed=round(float(np.sqrt(np.mean((sm[inner] - obs[inner]) ** 2))), 3),
               shoreline_bias=round(float(meta["skill"]["mean_bias_interior_m_yr"]), 3),
               hits=int(((o == 1) & (mo == 1)).sum()), observed=int(o.sum()),
               model_only=int(((o == 0) & (mo == 1)).sum()),
               cape_lrr=" / ".join(f"{rate.loc[g]:+.1f}" for g in CAPE),
               cape_obs=" / ".join(f"{target.get(g, np.nan):+.1f}" for g in CAPE),
               flat_cell_years=round(float(flat), 3))
    if start == 1996:
        mod = R.V.model_crests(start, VARIANT)
        surv = R.V.surveys(start)
        g09 = surv[surv.survey == 2009]
        k = int(g09.model_state.iloc[0][:4]) - R.V.PERIODS[start][1]
        errs = []
        for doms in R.GROUPS.values():
            s = g09[g09.gis.isin(doms)].dropna(subset=["crest"])
            idx = ((s.gis - 1) * R.CELLS + s.cell).to_numpy()
            errs.append(float(np.median(mod[k][idx]) - np.median(s.crest)))
        row["stable_dune_error_20_22_59_62"] = " / ".join(f"{e:+.2f}" for e in errs[::-1])
    G.TAG, G.EXP_DIR, G.run_dir = tag, exp_dir, rd_fn
    return row


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        return _launch(sys.argv[2])
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        todo = jobs()
        with ThreadPoolExecutor(WORKERS) as ex:
            list(ex.map(launch, todo))
    rows = []
    for rec, a in [("rec000", 0.0)] + list(RECOVERY.items()):
        for start in STARTS:
            if rec == "rec000":
                point_at_reference()
            else:
                point_at(rec)
            if G.run_dir(start, VARIANT) is None:
                print(f"missing {rec} {start}")
                continue
            rows.append(dict(recovery_m_yr=a, **score(start)))
    t = pd.DataFrame(rows)
    (STUDY / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(STUDY / "tables" / "recovery_scores.csv", index=False)
    with pd.option_context("display.width", 300, "display.max_columns", 30):
        for w, g in t.groupby("window"):
            print(g.dropna(axis=1, how="all").to_string(index=False), "\n")


if __name__ == "__main__":
    main()
