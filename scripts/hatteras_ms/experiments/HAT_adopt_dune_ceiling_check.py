r"""
HAT_adopt_dune_ceiling_check.py -- does the committed per-cell ceiling reproduce the experiment?
==============================================================================
Step 4 of adopting the per-cell dune ceiling (Hannah, 2026-09-28).

    reproduce  the four trim24 per-cell cases of
               storms-and-overwash/2026-09-28-dune-ceiling-per-domain
               (1996/2010 x full_management/natural) run on Barrier3D branch
               feature/per-cell-dune-ceiling (worktree ../Barrier3D-dune-ceiling,
               PYTHONPATH) with DuneCeilingFromStart switched on, against the
               experiment's runs, which set the same ceilings by wrapping the
               model in-process. The shoreline matrix and every domain's dune
               domain must be identical.
    off        with the setting off, the branch must equal the code every run
               uses now (49fd069): the natural 1996 run against the matrix run's
               shoreline matrix is not usable (the matrix was re-run on other
               ends), so the unit-level check is in the Barrier3D tests.

Nothing in the main code changes: the setting reaches the model through
cascade.brie_coupler.set_yaml in the run's own process, and the storm file is
swapped as in the experiments.

WHERE: output/raw_runs/experiments/code-checks/2026-09-28-per-cell-dune-ceiling-reproduces/
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

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(_HERE.parent))
import HAT_storm_length_selection as S  # noqa: E402
import HAT_dune_ceiling_per_domain as P  # noqa: E402

TAG = "code-checks/2026-09-28-per-cell-dune-ceiling-reproduces"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
WORKTREE = PROJECT_ROOT.parent / "Barrier3D-dune-ceiling"
CASES = [((1996, 2010), "full_management"), ((1996, 2010), "natural"),
         ((2010, 2024), "full_management"), ((2010, 2024), "natural")]


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


def run_dir(w, sc):
    base = EXP_DIR / "runs" / f"trim24_{sc}" / S.wtag(w) / "edgeBE"
    hits = sorted(base.glob("HAT_*")) if base.exists() else []
    return hits[-1] if hits else None


def launch(case):
    w, sc = case
    logs = EXP_DIR / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"HAT_IGNORE_SETTINGS": "1", "HAT_START_YEAR": str(w[0]),
                "HAT_SOURCE_SINK_PRESET": "edgeBE", "HAT_SCENARIO": sc,
                "HAT_RUN_KIND": "experiment", "HAT_RUN_TAG": f"{TAG}/runs/trim24_{sc}",
                "HAT_SAVE_MODEL_STATE": "1", "HAT_MAKE_GIFS": "0", "HAT_SHOW_FIGURES": "0",
                "HAT_OVERWRITE": "1", "PYTHONIOENCODING": "utf-8",
                "PYTHONPATH": str(WORKTREE) + os.pathsep + os.environ.get("PYTHONPATH", "")})
    t0 = time.time()
    log = logs / f"trim24_{sc}_{S.wtag(w)}.log"
    with open(log, "w", encoding="utf-8") as fh:
        p = subprocess.run([sys.executable, str(_HERE), "_launch", str(w[0]), str(S.storm_file(w, "trim24"))],
                           env=env, cwd=S.MD.HINDCAST.parent, stdout=fh, stderr=subprocess.STDOUT)
    b3d = next((ln.strip() for ln in log.read_text(encoding="utf-8", errors="replace").splitlines()
                if ln.startswith("BARRIER3D =")), "")
    print(f"{S.wtag(w)} {sc:16s} exit {p.returncode} {round((time.time() - t0) / 60, 1)} min  {b3d}", flush=True)


def compare():
    P.install_cell_growth()           # to unpickle nothing special; the experiment objects carry _hat_cell_dmax
    ok_all = True
    for w, sc in CASES:
        new, old = run_dir(w, sc), P.run_dir(w, "trim24", "cell", sc)
        mn = json.loads(next(new.glob("*_run_metadata.json")).read_text())
        mo = json.loads(next(old.glob("*_run_metadata.json")).read_text())
        ends_same = all(mn["index row"][k] == mo["index row"][k] for k in ("be_rate_gis1_m_yr", "be_rate_gis90_m_yr"))
        a = np.load(new / f"{new.name}_shoreline_matrix.npy")
        b = np.load(old / f"{old.name}_shoreline_matrix.npy")
        cn, co = S.load_state(new), S.load_state(old)
        dd = max(float(np.abs(np.asarray(x.DuneDomain) - np.asarray(y.DuneDomain)).max())
                 for x, y in zip(cn.barrier3d, co.barrier3d))
        ceil = max(float(np.abs(x._DuneCeiling - y._hat_cell_dmax).max()) for x, y in zip(cn.barrier3d, co.barrier3d))
        ok = ends_same and np.array_equal(a, b) and dd == 0.0 and ceil == 0.0
        ok_all &= ok
        print(f"{S.wtag(w)} {sc:16s} ends same {ends_same}  shoreline max|diff| {np.abs(a - b).max():.3g} m  "
              f"dunes max|diff| {dd:.3g}  ceilings max|diff| {ceil:.3g}  -> {'IDENTICAL' if ok else 'DIFFERENT'}  "
              f"[{mn['identity']['barrier3d_branch']}@{str(mn['identity']['barrier3d_commit'])[:7]}]")
    print("ALL IDENTICAL" if ok_all else "NOT IDENTICAL -- stop")


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(*sys.argv[2:4])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["run", "compare"])
    a = ap.parse_args()
    if a.action == "run":
        with ThreadPoolExecutor(4) as ex:
            list(ex.map(launch, CASES))
    else:
        compare()


if __name__ == "__main__":
    main()
