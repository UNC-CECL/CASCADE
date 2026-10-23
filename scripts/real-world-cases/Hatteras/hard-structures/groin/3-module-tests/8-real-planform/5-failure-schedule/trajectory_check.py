#!/usr/bin/env python3
"""
Trajectory check: score the instant-failure candidates against the observed gap at its own dates.

    python trajectory_check.py

An OLS trend can be matched by the wrong shape, so each candidate's D5-D6 gap
change is sampled at the wet/dry dates inside the window (full-model no-groin
offset spread linearly) and compared with the observed change. Candidates as in
failure_schedule_test.py. Writes trajectory_check.csv.
Details: README.md beside this script.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
# Non-code files beside this script live in the workspace on D: (hatteras_paths.workspace_twin)
import sys as _hat_sys  # noqa: E402
from pathlib import Path as _HatPath  # noqa: E402
_HAT_DIR = next(_p for _p in _HatPath(__file__).resolve().parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import workspace_twin  # noqa: E402
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents  if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path[:0] = [str(HERE), str(HERE.parent / "4-blocking-emulator"), str(HERE.parent / "1-dipole-stability"),
                str(workspace_twin(HERE.parents[1]) / "1-solver-audit"),
                str(_HAT_DIR / "hatteras_ms" / "groin-sweep"),
                str(_HAT_DIR / "hatteras_ms"),
                str(_HAT_DIR / "site_layer"), str(_HAT_DIR)]
import blocking_groin_emulator as bg  # noqa: E402
import failure_schedule_test as fs  # noqa: E402
from HAT_groin_sweep_config import WETDRY_CHANGE_TABLE  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
FAIL = 2004                                # the instant failure year
# -----------------------------------------------------------------------------

# Observed D5-D6 gap per wet/dry date, and its value at each window start
tab = pd.read_csv(WETDRY_CHANGE_TABLE).set_index("Domain_ID")
obs = {}
for c in tab.columns:
    m = re.match(r"change_from_wetdry_1967_wetdry_(\d{4})", c)
    if m:
        obs.setdefault(int(m.group(1)), []).append(tab.loc[5, c] - tab.loc[6, c])
obs = pd.Series({y: np.mean(v) for y, v in obs.items()}).sort_index()
start_gap = {1996: obs[1996], 2010: np.interp(2010, obs.index, obs.values)}

# Candidates to score, (kind, strength, f); the instant schedule swapped into the emulator
cands = ([("none", 0, 0)]
         + [("dipole", M, f) for M in (4, 6, 7, 8, 9, 10, 11, 12, 15) for f in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)]
         + [("blocking", round(b, 2), f) for b in np.arange(0.2, 1.01, 0.05)
            for f in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)])
bg.b_of_year = lambda b0, f, year: b0 if year < FAIL else b0 * f

# Score every candidate at the observed dates in both windows
rows = []
for start in (1996, 2010):
    x0 = bg.load_planform(start)
    base = bg.solve(x0, start, 0.0, 1.0, "exact")
    off = bg.FULL_NOGROIN[start] - bg.gap_change(base)
    yrs = [y for y in obs.index if start < y <= start + bg.YEARS]
    for kind, s, f in cands:
        st = (base if kind == "none" else
              bg.solve(x0, start, s, f, "callback") if kind == "blocking" else
              fs.dipole(x0, start, s, f, FAIL))
        g = st[:, bg.D] - st[:, bg.U]
        g = g - g[0] + off * np.arange(g.size) / (g.size - 1)
        model = np.interp(np.array(yrs) - start, np.arange(g.size), g)
        o = np.array([obs[y] for y in yrs]) - start_gap[start]
        rows.append(dict(start=start, kind=kind, strength=s, f=f,
                         ols=bg.gap_change(st) + off,
                         rmse_dates=float(np.sqrt(np.mean((model - o) ** 2))),
                         peak=float(g.max()), model=" ".join(f"{v:+.0f}" for v in model),
                         observed=" ".join(f"{v:+.0f}" for v in o)))
# Write the table, print each window, then the joint ranking
t = pd.DataFrame(rows)
(workspace_twin(HERE) / "tables").mkdir(exist_ok=True)
t.to_csv(workspace_twin(HERE) / "tables" / "trajectory_check.csv", index=False)
pd.set_option("display.width", 250, "display.max_colwidth", 60)
for start, sub in t.groupby("start"):
    print(f"\n=== {start}: dates {[y for y in obs.index if start < y <= start + 14]}, observed {sub.observed.iloc[0]}")
    print(sub.drop(columns=["start", "observed"]).round(1).to_string(index=False))
j = t.pivot_table(index=["kind", "strength", "f"], columns="start", values="rmse_dates")
j["joint"] = np.sqrt((j[1996] ** 2 + j[2010] ** 2) / 2)
print("\nJOINT date-RMSE (m):"); print(j.sort_values("joint").head(15).round(1).to_string())
