"""
Foreshore slope per domain, MHW (0.36 m NAVD88) to the berm (1.7 m NAVD88), on the 100 m transects (2026-11-23).

    python scripts/real-world-cases/Hatteras/input_prep/10-beach-slope/HAT_beach_slope_mhw_berm.py [survey ...]

The storm series' run-up (3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py) uses one beach slope, 0.06,
from the NC state team. This measures Hatteras's own: on each 100 m transect, the first point above MHW from the
ocean (with real data just seaward of it, so a data gap is not taken for the water's edge) and the first point above
the berm landward of it; slope = 1.34 m over that distance. Profiles sampled at 1 m, a 5 m median filter for spikes,
slopes outside 0.01-0.3 dropped. Domains in chunks of five so each grid stays small. Surveys default to 1996 and
1997, the two before Bonnie (1998); any survey of either start (1996-2009, 2011-2020) can be named. Elevations on the
2009 datum, each with its own start's offsets (survey_sources.csv), as hat_lidar_strips reads them.

Writes (D:, hatteras_init/10-beach-slope/1-foreshore-mhw-to-berm/), keeping surveys measured in earlier runs:
    foreshore_slope_transects.csv   one row per transect and survey
    foreshore_slope_domains.csv     per domain and survey: median, count, 25-75%

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
"""
import sys
from pathlib import Path

_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))

import numpy as np
import pandas as pd
from hatteras_paths import WORKSPACE
from site_layer import hat_lidar_strips as L

# --- CONFIG ------------------------------------------------------------------
MHW_NAVD88_M = 0.36       # the storm script's MHW
BERM_NAVD88_M = 1.7       # the storm script's berm
SLOPE_RANGE = (0.01, 0.3)  # outside it: flat water returns or a wall, not a beach
CHUNK = 5                 # domains per grid
SURVEYS = tuple(int(a) for a in sys.argv[1:]) or (1996, 1997)
OUT = WORKSPACE / "data" / "hatteras_init" / "10-beach-slope" / "1-foreshore-mhw-to-berm"
# -----------------------------------------------------------------------------


def load(year, g, src_tab, entries):
    """One survey on the grid, m NAVD88 on the 2009 datum (as hat_lidar_strips.surveys)."""
    import rasterio
    kind, path, *_ = entries[year]
    files = sorted(path.parent.glob(path.name)) if "*" in path.name else [path]
    srcs = [rasterio.open(f) for f in files]
    try:
        r = src_tab.loc[year]
        unit = "m" if kind == "mosaic" else r.vertical_unit
        scale = L.S.US_FT if unit == "ft" else 1.0
        resid = pd.to_numeric(r.road_offset_vs_2009_m, errors="coerce")
        z = L.S.raster_on_clip(srcs, g) * scale + float(r.vertical_offset_applied_m) - (0 if pd.isna(resid) else resid)
    finally:
        for s in srcs:
            s.close()
    z[(z < -5) | (z > 30)] = np.nan
    return z


def foreshore_points(p):
    """The filtered profile and the indices of the MHW and berm crossings (ocean first), or None."""
    from scipy.ndimage import median_filter
    p = median_filter(np.where(np.isfinite(p), p, -99.0), size=5)
    p[(p < -5) | (p > 12)] = np.nan
    ok = np.isfinite(p)
    if ok.sum() < 50:
        return None
    i_m = int(np.argmax(ok & (p >= MHW_NAVD88_M)))
    i_b = int(np.argmax(ok & (p >= BERM_NAVD88_M) & (np.arange(len(p)) > i_m)))
    if p[i_m] < MHW_NAVD88_M or p[i_b] < BERM_NAVD88_M or i_b <= i_m:
        return None
    if i_m == 0 or not (ok[i_m - 1] and p[i_m - 1] < MHW_NAVD88_M):
        return None
    return p, i_m, i_b


def foreshore(p, s):
    """Slope MHW to berm on one profile (ocean first), or None."""
    pts = foreshore_points(p)
    if pts is None:
        return None
    _, i_m, i_b = pts
    width = s[i_b] - s[i_m]
    slope = (BERM_NAVD88_M - MHW_NAVD88_M) / width
    return (slope, width) if SLOPE_RANGE[0] <= slope <= SLOPE_RANGE[1] else None


def main():
    import geopandas as gpd
    from scipy.ndimage import map_coordinates
    # every survey of both starts, each with its own start's sources table (units and datum offsets onto 2009)
    entries, tabs = {}, []
    for start, (_, ent, *_) in L.S.STARTS.items():
        entries.update(ent)
        tabs.append(pd.read_csv(L.S.tv.product_dir(start) / "survey-crests" / "survey_sources.csv"))
    src_tab = pd.concat(tabs).drop_duplicates("survey").set_index("survey")
    rows = []
    for lo in range(1, 91, CHUNK):
        gis = range(lo, min(lo + CHUNK, 91))
        g = L.grid(L.Region("chunk", "chunk", gis, (lo,)))
        t, _, crs = g
        tr = gpd.read_file(L.TRANSECTS).to_crs(crs)
        dom = gpd.read_file(L.S.GIS_ROOT / "domains.geojson").to_crs(crs)
        for year in SURVEYS:
            z = load(year, g, src_tab, entries)
            for d in gis:
                poly = dom[dom.domain_id == d].geometry.iloc[0]
                for _, row in tr[tr.intersects(poly)].iterrows():
                    c = np.asarray(row.geometry.coords, float)
                    u = (c[-1] - c[0]) / np.hypot(*(c[-1] - c[0]))          # ocean end -> land end
                    s = np.arange(0, row.geometry.length, 1.0)
                    xy = c[0] + s[:, None] * u
                    col, rr = (xy[:, 0] - t.c) / t.a - 0.5, (xy[:, 1] - t.f) / t.e - 0.5
                    res = foreshore(map_coordinates(z, [rr, col], order=1, mode="constant", cval=np.nan), s)
                    if res:
                        rows.append(dict(survey=year, gis=d, transect=row["Transects_100m.LineID"],
                                         slope=round(res[0], 4), mhw_to_berm_m=res[1]))
        print(f"GIS {gis.start}-{gis.stop - 1} done", flush=True)
    t = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    old = OUT / "foreshore_slope_transects.csv"
    if old.exists():                                    # keep the other surveys, replace the ones measured now
        prev = pd.read_csv(old)
        t = pd.concat([prev[~prev.survey.isin(SURVEYS)], t]).sort_values(["survey", "gis", "transect"])
    t.to_csv(old, index=False)
    d = (t.groupby(["survey", "gis"]).slope
         .agg(median="median", n="size", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75))
         .round(4).reset_index())
    d.to_csv(OUT / "foreshore_slope_domains.csv", index=False)
    for y in sorted(d.survey.unique()):
        m = d[d.survey == y]["median"]
        print(f"{y}: {m.size} domains, median {m.median():.3f}, 10-90% {m.quantile(.1):.3f}-{m.quantile(.9):.3f}, "
              f"below 0.06: {(m < 0.06).sum()}")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
