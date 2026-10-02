"""
Why do breached dunes never rebuild in the model?

    python scripts/hatteras_ms/experiments/HAT_dune_recovery_diagnosis.py

Reads the four matrix runs (1996 and 2010, managed and natural) and, per domain,
sets the dune cells left at the berm against the three things that can keep them
there: a ceiling below the storms, regrowth too slow, and the shoreline moving
the dune rows. Details: scripts/hatteras_ms/experiments/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
PROJECT_ROOT = next(_p for _p in _HERE.parents if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
sys.path.insert(0, str(PROJECT_ROOT / "scripts" / "input_prep" / "8-overwash-analysis" / "4-vs-model"))
import storms_vs_overwash as SV  # noqa: E402
from site_layer import hat_overwash as ow  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TAG = "storms-and-overwash/2026-10-01-dune-recovery-diagnosis"
EXP_DIR = PROJECT_ROOT / "output" / "raw_runs" / "experiments" / TAG
AT_BERM_M = 0.5          # a dune cell within this of the berm counts as flattened
MAX_YEARS = 50           # regrowth longer than this is reported as never
# -----------------------------------------------------------------------------


# Years for one cell to regrow from the restart height to clear `target` (m above berm), logistic as Barrier3D
def regrowth_years(r, ceiling, restart, target):
    if r <= 0 or ceiling <= target:
        return np.inf
    d, n = restart, 0
    while d < target and n < MAX_YEARS:
        d += r * d * (1 - d / ceiling)
        n += 1
    return n if d >= target else np.inf


# One run's per-domain table
def diagnose(window, arm, storms, false_alarms):
    name = SV.ovm.WINDOWS[window][arm]
    c = np.load(SV.ovm.MATRIX / f"{window[0]}_{window[1]}" / "edgeBE" / name / f"{name}.npz",
                allow_pickle=True)["cascade"][0]
    s = storms[storms.pi == list(SV.ovm.WINDOWS).index(window)]
    typical = s.groupby("calendar_year").rhigh_m.max().median()          # m MHW
    rows = []
    for gis in range(1, 91):
        p = SV.DOM.gis_to_pad(gis)
        b = c.barrier3d[p]
        berm = b._BermEl * 10                                          # m MHW
        dd = np.asarray(b.DuneDomain, float) * 10                      # m above berm
        crest = dd.max(axis=2)                                         # (years, cells)
        ceil = np.asarray(b._DuneCeiling, float) * 10                  # m above berm
        r = np.asarray(b._growthparam, float).ravel()
        r0 = np.asarray(c.roadways[p]._original_growth_param, float).ravel() if c.roadway_management_module[p] \
            and getattr(c.roadways[p], "_original_growth_param", None) is not None else r
        flat = crest[1:] < AT_BERM_M                                   # (run years, cells)
        sc = np.asarray(b.ShorelineChangeTS)[1:].astype(int)
        managed = bool(c.roadway_management_module[p])
        rebuilt = int(np.nansum(c.roadways[p]._dunes_rebuilt_TS)) if managed else 0
        target = typical - berm                                        # m above berm
        regrow = np.array([regrowth_years(ri, ci, b._DuneRestart * 10, target) for ri, ci in zip(r0, ceil)])
        ever = flat.any(axis=0)
        rows.append(dict(
            window=f"{window[0]}_{window[1]}", arm=arm, gis=gis,
            typical_storm_m_mhw=round(typical, 2),
            start_low_crest_m_mhw=crest[0].min() + berm, start_median_crest_m_mhw=np.median(crest[0]) + berm,
            ceiling_low_m_mhw=ceil.min() + berm, ceiling_median_m_mhw=np.median(ceil) + berm,
            share_cells_ceiling_below_storm=float((ceil + berm < typical).mean()),
            share_cell_years_flat=float(flat.mean()), cells_ever_flat=int(ever.sum()),
            share_end_flat=float(flat[-1].mean()),
            median_regrowth_years_flat_cells=float(np.median(regrow[ever])) if ever.any() else np.nan,
            share_flat_cells_never_clear=float(np.isinf(regrow[ever]).mean()) if ever.any() else np.nan,
            accretion_cells=int(sc[sc > 0].sum()), retreat_cells=int(-sc[sc < 0].sum()),
            road_managed=managed, dune_rebuilds=rebuilt,
            false_alarms=int(false_alarms.get((window[0], gis), 0)) if arm == "managed" else np.nan))
    return pd.DataFrame(rows)


# The dominant reason a domain's flattened cells stay flat
def cause(row):
    if row.cells_ever_flat == 0:
        return "no flattened cells"
    if row.accretion_cells >= 3 and row.share_end_flat > 0.9:
        return "progradation resets the dune rows"
    if row.share_flat_cells_never_clear >= 0.5:
        return "ceiling below a typical year's storm"
    if row.retreat_cells >= 3 and row.share_end_flat > 0.5:
        return "retreat brings low interior into the dune rows"
    return "regrowth slower than the storms"


def main():
    forcing = SV.sf.load_forcing()
    storms = SV.sf.classify(SV.sf.load_record(SV.env.DEFAULT_STORM_VARIANT, forcing), forcing)
    obs = pd.read_csv(ow.OBSERVATIONS / "overwash_observations.csv")
    fa = {}
    for w in SV.ovm.WINDOWS:
        d = SV.cells(w, obs)
        for gis, n in d[d.outcome == "false_alarm"].groupby("gis").size().items():
            fa[(w[0], gis)] = n
    t = pd.concat([diagnose(w, arm, storms, fa) for w in SV.ovm.WINDOWS for arm in ("managed", "natural")],
                  ignore_index=True)
    t["cause"] = t.apply(cause, axis=1)
    EXP_DIR.mkdir(parents=True, exist_ok=True)
    (EXP_DIR / "tables").mkdir(exist_ok=True)
    t.round(3).to_csv(EXP_DIR / "tables" / "dune_recovery_by_domain.csv", index=False)

    with pd.option_context("display.width", 250, "display.max_columns", 40, "display.max_rows", 400):
        print("Domains by cause (domains with any flattened cell):")
        print(t[t.cells_ever_flat > 0].groupby(["window", "arm", "cause"]).size().unstack(fill_value=0))
        hot = t[(t.arm == "managed") & (t.gis.isin([1, 2, 3, 4, 5, 6, 77, 78, 79, 80, 81, 82, 83, 84]))]
        print("\nHotspots, managed runs:")
        print(hot[["window", "gis", "start_low_crest_m_mhw", "ceiling_median_m_mhw", "typical_storm_m_mhw",
                   "share_cells_ceiling_below_storm", "share_end_flat", "median_regrowth_years_flat_cells",
                   "accretion_cells", "retreat_cells", "road_managed", "dune_rebuilds", "false_alarms",
                   "cause"]].round(2).to_string(index=False))
        m = t[t.arm == "managed"]
        print("\nManaged, island-wide: false alarms per domain against the share of cells whose ceiling "
              "sits below a typical year's storm")
        for w, g in m.groupby("window"):
            r = g[["false_alarms", "share_cells_ceiling_below_storm"]].corr(method="spearman").iloc[0, 1]
            r2 = g[["false_alarms", "share_cell_years_flat"]].corr(method="spearman").iloc[0, 1]
            print(f"  {w}: Spearman r (ceiling below storm) {r:.2f}; (cell-years flat) {r2:.2f}; "
                  f"domains road-managed {int(g.road_managed.sum())}, with a dune rebuild {int((g.dune_rebuilds > 0).sum())}")


if __name__ == "__main__":
    main()
