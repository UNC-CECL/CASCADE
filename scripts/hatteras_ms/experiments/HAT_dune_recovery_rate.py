"""
Does a dune that recovers from flat, instead of regrowing in proportion to its height, improve the overwash skill?

    python scripts/hatteras_ms/experiments/HAT_dune_recovery_rate.py run
    python scripts/hatteras_ms/experiments/HAT_dune_recovery_rate.py score

Adds a linear recovery term to Barrier3D's dune growth, in-process only, at
0.15, 0.3 and 0.5 m/yr; managed and natural, both windows; scored on overwash
skill, flattened dunes and shoreline. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
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
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "8-overwash-analysis" / "4-vs-model"))
import HAT_storm_max_duration as MD  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-and-overwash/2026-10-01-dune-recovery-rate"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
RATES_M = {"control": 0.0, "rec015": 0.15, "rec030": 0.30, "rec050": 0.50}   # m/yr to the front dune row
SCENARIOS = ("full_management", "natural")
ARM = {"full_management": "managed", "natural": "natural"}
AT_BERM_M = 0.5
# -----------------------------------------------------------------------------


def wtag(w):
    return f"{w[0]}_{w[1]}"


# The adopted storm file, as a path the runner can record
def storm_rel(w):
    from site_layer import hat_env_forcings as env
    f = env.storm_series_file(*w, variant=env.DEFAULT_STORM_VARIANT)
    return os.path.relpath(f, PROJECT_ROOT / "data" / "hatteras_init")


# Barrier3d.DuneGrowth as adopted, plus A * (1 - D/ceiling) where D is below its ceiling
def install_recovery(rate_m):
    from barrier3d import Barrier3d

    def DuneGrowth(self, DuneDomain, t):
        Dmax = max(self._Dmaxel - self._BermEl, 0)
        cell = getattr(self, "_DuneCeiling", None)
        if cell is not None:
            Dmax = cell
        Cf, Qdg = 3, 0
        for q in range(self._DuneWidth):
            reduc = 1 / (Cf ** q)
            d = DuneDomain[t - 1, :, q]
            room = 1 - d / Dmax
            G = self._growthparam * d * room * reduc
            G = G + (rate_m / 10.0) * np.clip(room, 0, None) * reduc       # the recovery term, dam
            DuneDomain[t, :, q] = G + d
            Qdg = Qdg + (np.sum(G) / self._BarrierLength)
        if cell is not None:
            Dmax = float(np.median(cell))
        return DuneDomain, Dmax, Qdg

    Barrier3d.DuneGrowth = DuneGrowth


# Child process: patch the growth, then run the unchanged hindcast
def _launch(start, storm_path, variant):
    install_recovery(RATES_M[variant])
    sys.argv = [sys.argv[0]]
    MD._launch(start, storm_path)


def run_dir(w, v, sc):
    base = EXP_DIR / "runs" / f"{v}_{sc}" / wtag(w) / "edgeBE"
    hits = sorted(d for d in base.glob("HAT_*") if any(d.glob("*_run_metadata.json"))) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    w, v, sc = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": sc,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{v}_{sc}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8"})
    t0 = time.time()
    log = logs / f"{v}_{sc}_{wtag(w)}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), storm_rel(w), v],
                           env=env, cwd=MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    b3d = next((ln.strip() for ln in log.read_text(encoding="utf-8", errors="replace").splitlines()
                if ln.startswith("BARRIER3D =")), "")
    rec = dict(window=wtag(w), variant=v, scenario=sc, returncode=p.returncode, barrier3d=b3d,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{wtag(w)} {v:8s} {sc:16s} exit {p.returncode} {rec['minutes']} min  {b3d}", flush=True)


# Overwash outcomes of one run, by overwash_vs_model's rule, with the run's own overwash swapped in
def outcomes(rd, w, sc, obs):
    import storms_vs_overwash as SV
    c = np.load(rd / f"{rd.name}.npz", allow_pickle=True)["cascade"][0]
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    orig = SV.ovm.load_qow
    SV.ovm.load_qow = lambda window, arm: (q, rd.name)
    try:
        d = SV.cells(w, obs)
    finally:
        SV.ovm.load_qow = orig
    return c, d


def score():
    import HAT_dune_ceiling_per_domain as P
    import storms_vs_overwash as SV
    from site_layer import hat_overwash as ow
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    zone = lambda g: "cape_point" if g <= 6 else ("community" if (7 <= g <= 8 or 21 <= g <= 31 or 68 <= g <= 83)  # noqa: E731
                                                  else "road")
    rows = []
    for sc in SCENARIOS:
        for w in MD.WINDOWS:
            for v in RATES_M:
                rd = run_dir(w, v, sc)
                if rd is None:
                    print(f"  missing {wtag(w)} {v} {sc}")
                    continue
                c, d = outcomes(rd, w, sc, obs)
                sk = SV.skill(d)
                flat = []
                for gis in range(1, 91):
                    b = c.barrier3d[SV.DOM.gis_to_pad(gis)]
                    cr = np.asarray(b.DuneDomain, float).max(axis=2)[1:] * 10
                    flat.append((cr < AT_BERM_M).mean())
                fa = d[d.outcome == "false_alarm"].assign(z=lambda x: x.gis.map(zone)).groupby("z").size()
                rows.append(dict(scenario=sc, window=wtag(w), variant=v, rate_m_yr=RATES_M[v],
                                 hit_rate=sk["hit_rate"], false_alarm_rate=sk["false_alarm_rate"], skill=sk["skill"],
                                 hits=sk["hits"], misses=sk["misses"], false_alarms=sk["false_alarms"],
                                 fa_cape_point=int(fa.get("cape_point", 0)), fa_community=int(fa.get("community", 0)),
                                 fa_road=int(fa.get("road", 0)), share_cell_years_flat=float(np.mean(flat)),
                                 **P.shoreline_lowess7(rd, w)))
                print(f"  scored {wtag(w)} {v} {sc}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    with pd.option_context("display.width", 260, "display.max_columns", 30):
        print(t.round(3).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:5])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "score"])
    ap.add_argument("--workers", type=int, default=5)
    a = ap.parse_args()
    if a.action == "run":
        jobs = [(w, v, sc) for sc in SCENARIOS for w in MD.WINDOWS for v in RATES_M if run_dir(w, v, sc) is None]
        with ThreadPoolExecutor(a.workers) as ex:
            list(ex.map(launch, jobs))
    else:
        score()


if __name__ == "__main__":
    main()
