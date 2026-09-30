r"""
HAT_trim_length_adopted.py -- does the storm trim length still matter once the dunes are realistic?
==============================================================================
WHY (Hannah, 2026-09-28: "run the trim-length check first"). trim24 was chosen
in storms-and-overwash/2026-09-28-storm-length-selection on the OLD dunes (held
near 3 m MHW by the default Dmaxel), and 24 h was the shortest length tried.
Barrier3D applies a storm's peak Rhigh for its whole duration, so the length is
a lever on how much sand overwash moves. This re-asks the question on the
setup being adopted:
    Barrier3D  branch hatteras/adopted (worktree ../Barrier3D-adopted): the
               three overwash fixes + per-cell dune ceilings, switched on
               (DuneCeilingFromStart true, floor 0.5 m) in each run's process
    storms     every event kept, trimmed to 12, 24, 48, 72 h, or full length
               (12: the builder, --long-events trim --max-duration 12, written
               here; 24: the adopted hindcast_storms v3_trim24 files; 48, 72,
               full: the storm-length selection's verified files)
    runs       managed (full_management), both windows, the site config's
               current end rates (the LOWESS-7 solve)
SCORES: overwash against the imagery (as before), the 2010 dune crest against
the lidar, interior RMSE/bias against LOWESS-7 (run_registry.skill_vs_target).
Scoring runs under the same Barrier3D, so the storm sharing uses the fixed
DuneGaps and the per-cell DuneGrowth.

NOTHING IN THE MAIN CODE CHANGES.

WHERE: output/raw_runs/experiments/storms-and-overwash/2026-09-28-trim-length-adopted/
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

TAG = "storms-and-overwash/2026-09-28-trim-length-adopted"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
# The worktree these runs used was removed on 2026-09-28 when hatteras/adopted
# was checked out in ../Barrier3D itself (the editable install); same branch.
WORKTREE = next((p for p in (PROJECT_ROOT.parent / "Barrier3D-adopted", PROJECT_ROOT.parent / "Barrier3D")
                 if (p / "barrier3d").exists()), PROJECT_ROOT.parent / "Barrier3D")
BUILDER = S.MD.BUILDER
TRIMS = ("trim12", "trim24", "trim48", "trim72", "full")


def storm_paths(w, v):
    """(npy, summary csv) for a variant."""
    from site_layer import hat_env_forcings as env
    if v == "trim12":
        d = EXP_DIR / "storms"
        return d / f"{S.wtag(w)}_storms_v3_trim12.npy", d / f"{S.wtag(w)}_storms_v3_trim12_summary.csv"
    if v == "trim24":
        f = env.storm_series_file(*w, variant="v3_trim24")
        return f, f.with_name(f.stem + "_summary.csv")
    return S.storm_file(w, v), S.summary_csv(w, v)


def build():
    (EXP_DIR / "storms").mkdir(parents=True, exist_ok=True)
    for w in S.WINDOWS:
        subprocess.run([sys.executable, str(BUILDER), "--start-year", str(w[0]), "--end-year", str(w[1]),
                        "--long-events", "trim", "--max-duration", "12", "--save-dir", str(EXP_DIR / "storms")],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        for v in TRIMS:
            a = np.load(storm_paths(w, v)[0])
            print(f"{S.wtag(w)} {v:7s} {len(a):4d} storms  max {a[:, 4].max():4.0f} h  storm-hours {a[:, 4].sum():6.0f}")


def _launch(start, storm_path):
    import cascade.brie_coupler as bc
    orig = bc.set_yaml

    def set_yaml(var_name, new_vals, file_name):
        orig(var_name, new_vals, file_name)
        if var_name not in ("DuneCeilingFromStart", "DuneCeilingFloor"):
            orig("DuneCeilingFromStart", True, file_name)
            orig("DuneCeilingFloor", 0.5, file_name)

    bc.set_yaml = set_yaml
    sys.argv = [sys.argv[0]]
    S.MD._launch(start, storm_path)


def run_dir(w, v):
    base = EXP_DIR / "runs" / f"{v}_full_management" / S.wtag(w) / "edgeBE"
    hits = sorted(base.glob("HAT_*")) if base.exists() else []
    return hits[-1] if hits else None


def launch(job):
    w, v = job
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": "full_management",
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/{v}_full_management",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8",
                "PYTHONPATH": str(WORKTREE) + os.pathsep + os.environ.get("PYTHONPATH", "")})
    t0 = time.time()
    log = logs / f"{v}_{S.wtag(w)}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), str(storm_paths(w, v)[0])],
                           env=env, cwd=S.MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    b3d = next((ln.strip() for ln in log.read_text(encoding="utf-8", errors="replace").splitlines()
                if ln.startswith("BARRIER3D =")), "")
    rec = dict(window=S.wtag(w), variant=v, returncode=p.returncode, barrier3d=b3d,
               minutes=round((time.time() - t0) / 60, 1))
    with open(logs / "launches.jsonl", "a", encoding="utf-8") as fh:
        fh.write(json.dumps(rec) + "\n")
    print(f"{S.wtag(w)} {v:7s} exit {p.returncode} {rec['minutes']} min  {b3d}", flush=True)


def score():
    import barrier3d
    if not Path(barrier3d.__file__).resolve().is_relative_to(WORKTREE.resolve()):
        raise SystemExit(f"score under the adopted Barrier3D: PYTHONPATH={WORKTREE} (got {barrier3d.__file__})")
    import HAT_dune_ceiling_per_domain as P
    import HAT_dune_ceiling_rebuild as E
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    from site_layer import hat_overwash as ow
    ovm = S.overwash_module()
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    pads = [DOM.gis_to_pad(g) for g in range(1, 91)]
    lidar = E.crest_series(S.load_state(P.ARCHIVED_MATRIX / "2010_2024/edgeBE" / S.CONTROLS[(2010, "full_management")]), pads)[0]
    low = lidar < np.percentile(lidar, 33)
    rows = []
    for w in S.WINDOWS:
        for v in TRIMS:
            rd = run_dir(w, v)
            if rd is None:
                print(f"  missing {S.wtag(w)} {v}")
                continue
            c = S.load_state(rd)
            summ = pd.read_csv(storm_paths(w, v)[1])
            cells = S.overwash_cells(c, w, summ, obs, ovm)
            cs = E.crest_series(c, pads)
            npy = np.load(storm_paths(w, v)[0])
            row = dict(window=S.wtag(w), storms=v, storm_hours=int(npy[:, 4].sum()),
                       **S.overwash_scores(cells, 0.0), **P.shoreline_lowess7(rd, w),
                       overwash_total_m3_per_m=float(np.array([np.asarray(b.QowTS)[1:].sum() for b in c.barrier3d])[pads].sum()))
            if w[0] == 1996:
                row["crest_2010_minus_lidar_m"] = float(np.median(cs[-1] - lidar))
            else:
                ir = cells[cells.obs_id == "OBS-018"].set_index("gis").reindex(range(1, 91))
                m, o = (ir.model_m3_per_m > 0).values, ir.observed.values
                for name, sel in (("low", low), ("rest", ~low)):
                    k = sel & ~np.isnan(o)
                    row[f"irene_{name}"] = f"{int(m[k].sum())}/{int(np.nansum(o[k]))}"
            rows.append(row)
            print(f"  scored {S.wtag(w)} {v}", flush=True)
    t = pd.DataFrame(rows)
    (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
    t.to_csv(EXP_DIR / "tables" / "scores.csv", index=False)
    cols = ["window", "storms", "storm_hours", "POD", "POFD", "PSS", "timing_r", "space_r", "rmse_interior_m_yr",
            "bias_interior_m_yr", "overwash_total_m3_per_m", "crest_2010_minus_lidar_m", "irene_low", "irene_rest"]
    with pd.option_context("display.width", 240, "display.max_columns", 20):
        print(t[[x for x in cols if x in t.columns]].round(2).to_string(index=False))


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:4])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["build", "run", "score"])
    a = ap.parse_args()
    if a.action == "build":
        build()
    elif a.action == "run":
        with ThreadPoolExecutor(5) as ex:
            list(ex.map(launch, [(w, v) for w in S.WINDOWS for v in TRIMS if run_dir(w, v) is None]))
    else:
        score()


if __name__ == "__main__":
    main()
