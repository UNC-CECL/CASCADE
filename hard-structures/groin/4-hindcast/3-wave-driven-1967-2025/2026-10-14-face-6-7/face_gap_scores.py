#!/usr/bin/env python3
"""
Score the face experiment's build-up runs on the gap across each candidate boundary, GIS 5/6 and 6/7.

    python face_gap_scores.py   ->  face_gap_scores.csv

face_experiment.py scores each run on the shape of its line along GIS 2-12. This scores the gap across
each boundary instead: the north domain's change since 1967 minus the south domain's, seaward positive,
against the wet/dry photos at the same 10 dates (1970-1995), so a run with the groin at 6/7 can be scored
on its own boundary as well as on 5/6. The observed gaps across both boundaries are in
../../../1-observations/analysis/1-aerial_photo_gap/tables/photo_boundary_gaps.csv. Reads the runs and the
photos with face_experiment.py's own settings and file lookup.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-14
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import face_experiment as fe  # noqa: E402  (also puts scripts/ on the path)
from site_layer.hatteras_site_config import HATTERAS_DOMAINS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
BOUNDARIES = ((5, 6), (6, 7))        # (south, north) GIS domain either side of each candidate face
REPORT_YEARS = (1978, 1993, 1995)    # model and observed gap written out at these dates
# -----------------------------------------------------------------------------


# Photo change since 1967 at each date the experiment scores, seaward positive, indexed by GIS domain
def observed():
    t = pd.read_csv(fe.OBSERVED).set_index("Domain_ID")
    out = {}
    for c in t.columns:
        if not c.startswith("change_from_wetdry_1967_wetdry_") or c.endswith("1967_m"):
            continue
        y = int(c.split("_")[5][:4])
        if fe.START < y <= fe.LAST and c.count("_") == 6 and t[c].reindex(fe.GIS).notna().sum() >= 8:
            out[y] = -t[c]
    return dict(sorted(out.items()))


def main():
    obs = observed()
    years = list(obs)
    gaps = {f"{s}_{n}": pd.Series({y: obs[y][n] - obs[y][s] for y in years}) for s, n in BOUNDARIES}
    rows = [dict(variant="observed", strength=np.nan,
                 **{f"gap_{k}_{y}_m": g[y] for k, g in gaps.items() for y in REPORT_YEARS})]
    for v, ss in fe.CELLS.items():
        for s in ss:
            f = fe.shoreline_file(v, s)
            if f is None:
                continue
            x = np.load(f)
            row = dict(variant=v, strength=s)
            for (south, north), (k, g) in zip(BOUNDARIES, gaps.items()):
                ps, pn = HATTERAS_DOMAINS.gis_to_pad(np.array([south, north]))
                # The runner's matrix is landward +: seaward change is minus the change in the matrix
                m = pd.Series({y: (x[y - fe.START, ps] - x[0, ps]) - (x[y - fe.START, pn] - x[0, pn]) for y in years})
                d = m - g
                row[f"rmse_gap_{k}_m"] = float(np.sqrt((d ** 2).mean()))
                row[f"bias_gap_{k}_m"] = float(d.mean())
                row.update({f"gap_{k}_{y}_m": float(m[y]) for y in REPORT_YEARS})
            rows.append(row)
    t = pd.DataFrame(rows)
    t.round(2).to_csv(HERE / "face_gap_scores.csv", index=False)
    print("photo dates:", " ".join(map(str, years)))
    pd.set_option("display.width", 250)
    print(t.round(1).to_string(index=False))


if __name__ == "__main__":
    main()
