"""
With the dune ceiling at the measured dune height, can a slower growth rate keep the low dunes low?

    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_dem_growth.py run
    python scripts/hatteras_ms/experiments/HAT_dune_ceiling_dem_growth.py score

Sets Dmaxel in-process to each period's own starting-DEM crest P99 (9.03 m
NAVD88 in 1996, 9.50 m in 2009) and rmin/rmax to 0.05-0.15, 0.10-0.20 and
0.15-0.30; and, keeping the adopted 5.5 m ceiling, rmin/rmax 0.05-0.15,
0.10-0.20, 0.15-0.30 and 0.25-0.45 (the c55_ variants). The control is the
parameter file as adopted (5.5 m, 0.55-0.95).
Managed edgeBE runs of both periods, scored as HAT_dune_ceiling_uniform.py
plus the low dunes at model year 3. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-19
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parent))
import HAT_dune_ceiling_uniform as U  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-overwash-and-dunes/2026-10-19-dune-ceiling-dem-growth"
CEILING_PCT = 99                     # percentile of the starting-DEM crests used as the ceiling
RANGES = {"control": None,           # parameter file: Dmaxel 5.5 m NAVD88, rmin/rmax 0.55/0.95
          "p99_r005_015": (0.05, 0.15), "p99_r010_020": (0.10, 0.20), "p99_r015_030": (0.15, 0.30),
          # the adopted 5.5 m ceiling kept, only the growth rate slowed
          "c55_r005_015": (0.05, 0.15), "c55_r010_020": (0.10, 0.20), "c55_r015_030": (0.15, 0.30),
          "c55_r025_045": (0.25, 0.45)}
EARLY_YEAR = 3                       # model year at which the low dunes are compared with the DEM
LOW_PCT = 10                         # percentile of the crests reported as "the low dunes"
# -----------------------------------------------------------------------------


# The period's ceiling from its own starting dunes (GIS 1-90), m NAVD88
def dem_ceiling(cascade):
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    crest = np.concatenate([U.cell_crests(cascade, DOM)[g - 1][0] for g in range(1, 91)])
    return float(U.to_navd88(cascade, DOM, np.percentile(crest, CEILING_PCT)))


# Child process: override rmin/rmax and set the ceiling, then run the unchanged hindcast
def _launch(variant):
    import runpy
    import cascade_pipeline.hindcast as H
    if RANGES[variant] is not None:
        lo, hi = RANGES[variant]
        orig = H.build_cascade

        def build_cascade(*a, **k):
            k["rmin"] = [lo] * len(k["rmin"])
            k["rmax"] = [hi] * len(k["rmax"])
            c = orig(*a, **k)
            if variant.startswith("c55_"):
                print(f"EXPERIMENT rmin/rmax = {lo}/{hi} (parameter file's 0.55/0.95 overridden; "
                      f"ceiling 5.5 m NAVD88 as adopted)", flush=True)
                return c
            value = dem_ceiling(c)
            U.set_uniform(c, value)
            print(f"EXPERIMENT dune ceiling = {value:.2f} m NAVD88 (starting-DEM P{CEILING_PCT}), "
                  f"rmin/rmax = {lo}/{hi} (parameter file's 5.5 m and 0.55/0.95 overridden)", flush=True)
            return c

        H.build_cascade = build_cascade
    sys.path.insert(0, str(U.HINDCAST.parent))
    sys.argv = [str(U.HINDCAST)]
    runpy.run_path(str(U.HINDCAST), run_name="__main__")


# Re-point the uniform study at this one; its lidar reference (the 2009 start crest) is this control's
def adopt():
    U.TAG, U.EXP_DIR, U._HERE = TAG, U.PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG, _HERE
    U.VARIANTS = {v: None for v in RANGES}
    own = U.run_dir
    U.run_dir = lambda start, variant: own(start, "control" if variant == "percell" else variant)
    U.overwash_skill = overwash_skill


# The current overwash scorer (overwash_vs_model.compare, keyed by start year since 10-12) on an experiment run
def overwash_skill(c, obs, dom, storm_years):
    import overwash_vs_model as ovm
    import storms_vs_overwash as SV
    start = next(s for s, (_, yrs) in U.PERIODS.items() if yrs == storm_years)
    window = (start, U.PERIODS[start][0])
    q = np.array([np.asarray(b.QowTS) for b in c.barrier3d]).T
    orig = ovm.load_qow
    ovm.load_qow = lambda w, arm: (q, "experiment")
    try:
        d = ovm.compare(window, "managed", obs, ovm.HEADLINE_THR).dropna(subset=["observed"])
    finally:
        ovm.load_qow = orig
    o, m = d.observed.astype(int), d.model.astype(int)
    d["outcome"] = np.select([(o == 1) & (m == 1), (o == 1) & (m == 0), (o == 0) & (m == 1)],
                             ["hit", "miss", "false_alarm"], "neither")
    return SV.skill(d)


# The low dunes and the ceiling at model year EARLY_YEAR, added to the uniform study's scores
def early_table():
    from site_layer.hatteras_site_config import HATTERAS_DOMAINS as DOM
    rows = []
    for start, (end, _) in U.PERIODS.items():
        for v in RANGES:
            rd = U.run_dir(start, v)
            if rd is None:
                continue
            c = U.load(rd)
            cc = U.cell_crests(c, DOM)
            first = np.concatenate([x[0] for x in cc])
            early = np.concatenate([x[EARLY_YEAR] for x in cc])
            b = c.barrier3d[DOM.gis_to_pad(1)]
            rows.append(dict(period=f"{start}_{end}", variant=v,
                             ceiling_m_navd88=round(b._Dmaxel * 10 + U.MHW_NAVD88_M, 2),
                             rmin=float(np.min(b._growthparam)), rmax=float(np.max(b._growthparam)),
                             low_start_m_navd88=round(float(U.to_navd88(c, DOM, np.percentile(first, LOW_PCT))), 2),
                             low_year3_m_navd88=round(float(U.to_navd88(c, DOM, np.percentile(early, LOW_PCT))), 2),
                             median_year3_m_navd88=round(float(U.to_navd88(c, DOM, np.median(early))), 2)))
    return pd.DataFrame(rows)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "_launch":
        _launch(sys.argv[2])
        return
    adopt()
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        from concurrent.futures import ThreadPoolExecutor
        jobs = [(s, v) for s in U.PERIODS for v in RANGES if U.run_dir(s, v) is None]
        with ThreadPoolExecutor(U.WORKERS) as ex:
            list(ex.map(U.launch, jobs))
    else:
        U.score()
        t = early_table()
        t.to_csv(U.EXP_DIR / "tables" / "low_dunes_year3.csv", index=False)
        with pd.option_context("display.width", 300, "display.max_columns", 40):
            print(t.to_string(index=False))


if __name__ == "__main__":
    main()
