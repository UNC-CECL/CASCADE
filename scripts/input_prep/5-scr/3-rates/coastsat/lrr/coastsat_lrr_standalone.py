"""
Shoreline change rate (LRR) for every CoastSat transect in a folder.

Input:  CoastSat time-series CSVs, one per transect (columns "dates UTC",
        "chainage (m)"), anywhere under FOLDER.
Output: CSV of transect_id, lrr_m_yr, unc_m_yr (95% CI half-width), n_obs.

Rate = ordinary least-squares slope of position against time, using every
position from 1 Jan START to 31 Dec END. No filtering, no weighting; at
least 3 positions. Positive = seaward (accretion), negative = erosion.

    python coastsat_lrr_standalone.py FOLDER START END OUT.csv
    python coastsat_lrr_standalone.py coastsat_timeseries 1984 2025 rates.csv

START and END are whole calendar years, both included -- the window
convention of every 5-scr script and of 5-scr/template/, whose
shoreline_rates_template.py is the fuller version of this file.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
"""
import sys
from pathlib import Path

import pandas as pd
from scipy import stats

folder, start, end, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]

rows = []
for f in sorted(Path(folder).rglob("*.csv")):
    d = pd.read_csv(f)
    date = pd.to_datetime(d.iloc[:, 0], utc=True)
    pos = pd.to_numeric(d.iloc[:, 1], errors="coerce")
    keep = date.dt.year.between(start, end) & pos.notna()
    date, pos = date[keep], pos[keep]
    if len(pos) < 3:
        continue
    years = (date - date.min()).dt.total_seconds() / (365.25 * 86400)
    fit = stats.linregress(years, pos)
    ci = stats.t.ppf(0.975, len(pos) - 2) * fit.stderr
    rows.append((f.stem, fit.slope, ci, len(pos)))

table = pd.DataFrame(rows, columns=["transect_id", "lrr_m_yr", "unc_m_yr", "n_obs"])
table.round(4).to_csv(out, index=False)
print(f"{len(table)} transects -> {out}")
