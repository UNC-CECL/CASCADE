"""
Are the beach nourishment volumes applied as intended, for every fill in both windows?

    python scripts/hatteras_ms/experiments/HAT_nourishment_volume_check.py

Follows each fill from the nourishment sheet to the model shoreline: source vs
config, cubic yards to m^3/m, what each BeachDuneManager recorded, the shoreline
step it caused against the Ashton and Lorenzo-Trueba formula, and the total
volume placed. Reads saved full_management runs; runs nothing.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-03
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from cascade_pipeline import nourishment as N  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS, HATTERAS_NOURISHMENT_PROJECTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / "management" / "2026-10-03-nourishment-volume-check"
RAW = PROJECT_ROOT / "output" / "raw_runs" / "experiments"
RUNS = {
    (1996, 2015): EXP_DIR / "runs" / "full_management" / "1996_2015" / "edgeBE",
    (2010, 2026): RAW / "management" / "2026-10-03-buxton-2017-fill" / "runs" / "with2017" / "2010_2026" / "edgeBE",
}
SHEET = PROJECT_ROOT / "data" / "hatteras_init" / "4-mgmt-forcing" / "nourishment" / "1-sources" / "Hatteras_BN_data.xlsx"
# Sheet row for each model project: (location, yearCompleted); Buxton 2017 is completed 2018
SHEET_ROW = {("Rodanthe emergency fill", 2014): ("Pea Island/Rodanthe", 2014),
             ("Buxton beach nourishment", 2017): ("Buxton", 2018),
             ("Buxton shore protection", 2022): ("Buxton", 2022),
             ("Avon shore protection", 2022): ("Avon", 2022)}
FT_TO_M = 0.3048
TOL_M = 0.01                  # shoreline step agreement, m
# -----------------------------------------------------------------------------


def source_vs_config():
    sh = pd.read_excel(SHEET)
    rows = []
    for p in HATTERAS_NOURISHMENT_PROJECTS:
        loc, yr = SHEET_ROW[(p.name, p.year)]
        s = sh[(sh.location == loc) & (sh.yearCompleted == yr)].iloc[0]
        model_len = len(p.gis_domains) * HATTERAS_DOMAINS.domain_spacing_m
        rows.append(dict(project=p.name, model_year=p.year, sheet_year=yr,
                         sheet_cy=int(s.volume), model_cy=int(p.volume_cubic_yards),
                         cy_diff_pct=100 * (p.volume_cubic_yards - s.volume) / s.volume,
                         sheet_length_m=round(s.length * FT_TO_M), model_length_m=model_len,
                         gis=f"{min(p.gis_domains)}-{max(p.gis_domains)}",
                         model_m3_per_m=p.volume_m3_per_m(HATTERAS_DOMAINS.domain_spacing_m),
                         hand_m3_per_m=p.volume_cubic_yards * 0.764555 / model_len))
    return pd.DataFrame(rows)


def load(window):
    base = RUNS[window]
    rd = sorted(d for d in base.glob("HAT_*_road_bdm_nourish_nogroin") if any(d.glob("*.npz")))[-1]
    return rd, np.load(rd / f"{rd.name}.npz", allow_pickle=True)["cascade"][0]


def in_model(window):
    rd, c = load(window)
    sched = N.build_schedule(HATTERAS_NOURISHMENT_PROJECTS, HATTERAS_DOMAINS, *window)
    want = {(r["gis"], r["year"]): r["volume_m3_per_m"] for r in sched.events()}
    rows, seen = [], set()
    for pad, mgr in enumerate(c.nourishments):
        if mgr is None:
            continue
        b = c.barrier3d[pad]
        gis = int(HATTERAS_DOMAINS.pad_to_gis(pad))
        for i in np.flatnonzero(np.asarray(mgr._nourishment_TS) == 1):
            year = window[0] + int(i) - 1
            v = float(mgr._nourishment_volume_TS[i])
            pre, post = float(mgr._post_storm_x_s[i]), float(b.x_s_TS[i])
            h_b, dsf = float(b.h_b_TS[i]), float(b.DShoreface)
            # The formula, in dam as the manager calls it: dx = 2V / (2 h_b + D_sf)
            expect_m = 10 * 2 * (v / 100) / (2 * h_b + dsf)
            seen.add((gis, year))
            rows.append(dict(window=f"{window[0]}-{window[1]}", year=year, gis=gis,
                             scheduled_m3_per_m=want.get((gis, year), np.nan), applied_m3_per_m=v,
                             step_measured_m=(pre - post) * 10, step_formula_m=expect_m,
                             h_b_m=h_b * 10, shoreface_depth_m=dsf * 10))
    missing = sorted(k for k in want if k not in seen)
    return pd.DataFrame(rows), missing, rd.name


def main():
    pd.set_option("display.width", 220)
    src = source_vs_config()
    print("\n1. Source sheet vs model config")
    print(src.round(2).to_string(index=False))

    out = []
    for w in RUNS:
        try:
            d, missing, name = in_model(w)
        except (IndexError, FileNotFoundError):
            print(f"\n  no saved run for {w}")
            continue
        print(f"\n2. {w[0]}-{w[1]}: {name}")
        if d.empty:
            print("  no fills fired")
        else:
            d["volume_ok"] = (d.applied_m3_per_m - d.scheduled_m3_per_m).abs() < 1e-6
            d["step_ok"] = (d.step_measured_m - d.step_formula_m).abs() < TOL_M
            g = d.groupby("year").agg(domains=("gis", "size"), gis_lo=("gis", "min"), gis_hi=("gis", "max"),
                                      applied_m3_per_m=("applied_m3_per_m", "mean"),
                                      volume_ok=("volume_ok", "all"), step_m_min=("step_measured_m", "min"),
                                      step_m_max=("step_measured_m", "max"), step_ok=("step_ok", "all"),
                                      unscheduled=("scheduled_m3_per_m", lambda s: int(s.isna().sum())))
            g["total_m3_placed"] = d.groupby("year").applied_m3_per_m.sum() * HATTERAS_DOMAINS.domain_spacing_m
            print(g.round(2).to_string())
        print(f"  scheduled but never applied: {missing or 'none'}")
        out.append(d)
    if out:
        (EXP_DIR / "tables").mkdir(parents=True, exist_ok=True)
        pd.concat(out).to_csv(EXP_DIR / "tables" / "fills_by_domain.csv", index=False)
        src.to_csv(EXP_DIR / "tables" / "source_vs_config.csv", index=False)


if __name__ == "__main__":
    main()
