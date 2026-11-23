"""
Berm elevation per domain from the lidar, to check the model's berm (1.7 m NAVD88; source not recorded) (2026-11-23).

    python scripts/real-world-cases/Hatteras/input_prep/11-datums/HAT_berm_elevation.py [survey ...]

On each 100 m transect, sampled at 1 m from the ocean end and median-filtered over 5 m as the beach slope is
(input_prep/10-beach-slope/HAT_beach_slope_mhw_berm.py): from the MHW crossing (0.36 m NAVD88) landward, the berm
crest is the first point at least MIN_ABOVE_MHW above MHW where the slope over the next 5 m is FLAT_SLOPE or less (the
beach turns flat). A transect that climbs past NO_BERM_ABOVE without flattening runs straight from the foreshore into
the dune: no berm, counted. Each domain's berm is the median of its transects with a berm. Surveys default to 1996.

Writes (D:, hatteras_init/11-datums/1-berm/), keeping surveys measured in earlier runs:
    berm_transects.csv   one row per transect and survey: berm crest elevation (NaN where none)
    berm_domains.csv     per domain and survey: median, count with a berm, count without

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_HAT_DIR = next(p for p in _HERE.parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "10-beach-slope"))

import numpy as np
import pandas as pd
import HAT_beach_slope_mhw_berm as M
from hatteras_paths import WORKSPACE
from site_layer import hat_lidar_strips as L

# --- CONFIG ------------------------------------------------------------------
MHW_NAVD88_M = 0.36
MIN_ABOVE_MHW = 0.5       # m: the berm crest stands at least this far above MHW
FLAT_SLOPE = 0.02         # rise over run across the next 5 m at which the beach counts as flat
NO_BERM_ABOVE = 4.0       # m NAVD88: climbing past this without flattening = straight into the dune
SURVEYS = tuple(int(a) for a in sys.argv[1:]) or (1996,)
OUT = WORKSPACE / "data" / "hatteras_init" / "11-datums" / "1-berm"
# -----------------------------------------------------------------------------


def berm_crest(raw):
    """(berm crest elevation or NaN, whether the profile was readable)."""
    from scipy.ndimage import median_filter
    p = median_filter(np.where(np.isfinite(raw), raw, -99.0), size=5)
    p[(p < -5) | (p > 12)] = np.nan
    ok = np.isfinite(p)
    if ok.sum() < 50:
        return np.nan, False
    i_m = int(np.argmax(ok & (p >= MHW_NAVD88_M)))
    if p[i_m] < MHW_NAVD88_M or i_m == 0 or not (ok[i_m - 1] and p[i_m - 1] < MHW_NAVD88_M):
        return np.nan, False
    for i in range(i_m, len(p) - 5):
        if not (ok[i] and ok[i + 5]):
            return np.nan, False
        if p[i] > NO_BERM_ABOVE:
            return np.nan, True                       # into the dune without a flat berm
        if p[i] >= MHW_NAVD88_M + MIN_ABOVE_MHW and (p[i + 5] - p[i]) / 5.0 <= FLAT_SLOPE:
            return float(p[i]), True
    return np.nan, False


def main():
    import geopandas as gpd
    from scipy.ndimage import map_coordinates
    entries, tabs = {}, []
    for start, (_, ent, *_) in L.S.STARTS.items():
        entries.update(ent)
        tabs.append(pd.read_csv(L.S.tv.product_dir(start) / "survey-crests" / "survey_sources.csv"))
    src_tab = pd.concat(tabs).drop_duplicates("survey").set_index("survey")
    rows = []
    for lo in range(1, 91, M.CHUNK):
        gis = range(lo, min(lo + M.CHUNK, 91))
        g = L.grid(L.Region("chunk", "chunk", gis, (lo,)))
        t, _, crs = g
        tr = gpd.read_file(L.TRANSECTS).to_crs(crs)
        dom = gpd.read_file(L.S.GIS_ROOT / "domains.geojson").to_crs(crs)
        for year in SURVEYS:
            z = M.load(year, g, src_tab, entries)
            for d in gis:
                poly = dom[dom.domain_id == d].geometry.iloc[0]
                for _, row in tr[tr.intersects(poly)].iterrows():
                    c = np.asarray(row.geometry.coords, float)
                    u = (c[-1] - c[0]) / np.hypot(*(c[-1] - c[0]))
                    s = np.arange(0, row.geometry.length, 1.0)
                    xy = c[0] + s[:, None] * u
                    col, rr = (xy[:, 0] - t.c) / t.a - 0.5, (xy[:, 1] - t.f) / t.e - 0.5
                    zb, readable = berm_crest(map_coordinates(z, [rr, col], order=1, mode="constant",
                                                              cval=np.nan))
                    if readable:
                        rows.append(dict(survey=year, gis=d, transect=row["Transects_100m.LineID"],
                                         berm_m_navd88=round(zb, 3) if np.isfinite(zb) else np.nan))
        print(f"GIS {gis.start}-{gis.stop - 1} done", flush=True)
    t = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    f = OUT / "berm_transects.csv"
    if f.exists():
        prev = pd.read_csv(f)
        t = pd.concat([prev[~prev.survey.isin(SURVEYS)], t]).sort_values(["survey", "gis", "transect"])
    t.to_csv(f, index=False)
    d = (t.groupby(["survey", "gis"]).berm_m_navd88
         .agg(median="median", with_berm="count", without_berm=lambda x: int(x.isna().sum())).round(3).reset_index())
    d.to_csv(OUT / "berm_domains.csv", index=False)
    for y in sorted(d.survey.unique()):
        x = d[d.survey == y]
        tt = t[t.survey == y].berm_m_navd88
        print(f"{y}: transects with a berm {tt.notna().sum()} of {len(tt)}; island median of domain medians "
              f"{x['median'].median():.2f} m NAVD88 (mean {x['median'].mean():.2f}, 10-90% "
              f"{x['median'].quantile(.1):.2f}-{x['median'].quantile(.9):.2f}); domains with no berm on any "
              f"transect: {int((x.with_berm == 0).sum())}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
