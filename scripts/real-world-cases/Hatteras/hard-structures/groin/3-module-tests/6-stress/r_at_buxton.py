#!/usr/bin/env python3
"""
BRIE's diffusion number r and the groin face angle at Buxton, in every groin-on run on disk.

    python r_at_buxton.py   ->  tables/r_at_buxton.csv

For each matrix and wave-sweep run with the groin on, r = D dt / (2 dy^2) at the groin face (pads 19|20) and
the largest r in the 10-cell correction zone, every year, from the run's BRIE shoreline and its own wave
climate (read from the run name). Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-15
"""
# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import stress_test as st  # noqa: E402
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
ROOT = str(st.REPO / "output" / "raw_runs")
LO = 19
near = np.abs(np.arange(120) - (LO + 0.5)) <= 10
def r_of(b, x):
    th = np.degrees(np.arctan2(np.roll(x, -1) - x, 500.0))
    k = np.clip(np.round(90 - th).astype(int), 1, b._wave_climl - 1)
    return np.maximum(0.0, b._coast_diff[k] * 1.0 / 2 / 500.0 ** 2), th
def climate(name):
    c = dict(st.OPTION_A)
    m = re.search(r"wave(Hs|Tp|ahf|asym)(\d+(?:p\d+)?)", name)
    if m:
        c[m.group(1)] = float(m.group(2).replace("p", "."))
    return c, (m.group(0) if m else "option A")
rows = []
files = glob.glob(f"{ROOT}/matrix/*/*/*_groin*/*_shoreline_matrix.npy") + glob.glob(f"{ROOT}/sensitivity/*/*/*/*_groin*/*_shoreline_matrix.npy")
for f in files:
    name = os.path.basename(os.path.dirname(f)); c, lab = climate(name)
    xs = np.load(f); b = st.brie_like(c)
    rf, rz, ang = [], [], []
    for x in xs:
        r, th = r_of(b, x)
        rf.append(r[LO]); rz.append(r[near].max()); ang.append(th[LO])
    rows.append(dict(set=("sweep " if "sensitivity" in f else "matrix ") + lab, run=name,
                     r_face_min=min(rf), r_face_max=max(rf), r_zone_max=max(rz), ang_min=min(ang), ang_max=max(ang)))
t = pd.DataFrame(rows)
print(len(files), "runs")
g = t.groupby("set").agg(n=("run","count"), r_face_min=("r_face_min","min"), r_face_max=("r_face_max","max"),
                        r_zone_max=("r_zone_max","max"), ang_min=("ang_min","min"), ang_max=("ang_max","max"))
print(g.round(3).to_string())
t.to_csv(workspace_twin(HERE) / "tables" / "r_at_buxton.csv", index=False, float_format="%.4f")
