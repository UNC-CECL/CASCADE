"""
Score every candidate DEM year by how much of the 2009 DEM's dry-land gaps it measures.

    python scripts/input_prep/0-elevation/1-source-selection/HAT_survey_dem_coverage.py [--domains 78,79,80]

Dry-land gaps are 2009 nodata cells the candidates' median puts above MHW.
Candidates that stamp flat water surfaces (hydro-flattened, like 2014 NCFMP) are
flagged. Reads the D: GIS drive; writes a coverage table beside this step's
data. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window, from_bounds

POLY_ROOT = Path(r"D:\Hatteras_GIS\Elevation\Polygons")
BASE_DEM = POLY_ROOT / "2009" / "usace2009_nc_dem_Job1076020" / "2009_full.tif"
DOMAIN_FILE = Path(r"D:\Hatteras_GIS\domains.geojson")
DOMAIN_ID_FIELD = "domain_id"
# Walk up until a directory holds data/hatteras_init
def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
import sys as _elsys
from pathlib import Path as _ELP
_elsys.path.insert(0, str(next(_q for _q in _ELP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_elevation_products as _el  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
# Source-selection/, where the committed copies are
OUT_DIR = _el.source_selection_dir()

MHW = 0.36                 # m NAVD88
CLIP = (500, 2000)
FLAT_FRACTION_WARN = 0.10  # top value > this share of coverage -> flag
# -----------------------------------------------------------------------------


# Every candidate DEM folder under POLY_ROOT except the base
def candidates():
    out = []
    for ydir in sorted(POLY_ROOT.iterdir()):
        if not ydir.is_dir():
            continue
        for dsdir in sorted(p for p in ydir.iterdir() if p.is_dir()):
            if dsdir.resolve() == BASE_DEM.parent.resolve():
                continue
            full = list(dsdir.glob("*_full.tif"))
            tifs = full if full else sorted(
                p for p in dsdir.glob("*.tif")
                if not p.name.endswith("_full.tif"))
            if tifs:
                out.append({"year": ydir.name, "name": dsdir.name, "tifs": tifs})
    return out


# A candidate's EPSG, vertical datum and resolution
def describe(tifs):
    with rasterio.open(tifs[0]) as s:
        crs = CRS.from_wkt(s.crs.to_wkt()) if s.crs else None
        hz, vert = crs, "unknown"
        if crs is not None and crs.is_compound:
            hz, vert = crs.sub_crs_list[0], crs.sub_crs_list[1].name
        return {"epsg": hz.to_epsg() if hz is not None else None,
                "vertical": vert, "res": round(s.transform.a, 3),
                "n_files": len(tifs)}


# WarpedVRT rejects boundless reads, so partial overlap is intersected and pasted at the right offset
def read_on_grid(vrts, bounds, shape):
    bx0, by0, bx1, by1 = bounds
    out = np.full(shape, np.nan, np.float32)
    for vrt in vrts:
        vb = vrt.bounds
        ix0, iy0 = max(bx0, vb.left), max(by0, vb.bottom)
        ix1, iy1 = min(bx1, vb.right), min(by1, vb.top)
        if ix1 <= ix0 or iy1 <= iy0:
            continue
        h, w = int(round(iy1 - iy0)), int(round(ix1 - ix0))
        if h <= 0 or w <= 0:
            continue
        d = vrt.read(1, window=from_bounds(ix0, iy0, ix1, iy1,
                                           transform=vrt.transform),
                     out_shape=(h, w), masked=True,
                     resampling=Resampling.nearest)
        d = np.ma.filled(d.astype(np.float32), np.nan)
        r0, c0 = max(int(round(by1 - iy1)), 0), max(int(round(ix0 - bx0)), 0)
        hh, ww = min(h, shape[0] - r0), min(w, shape[1] - c0)
        if hh <= 0 or ww <= 0:
            continue
        sub = out[r0:r0 + hh, c0:c0 + ww]
        dd = d[:hh, :ww]
        take = np.isnan(sub) & ~np.isnan(dd)
        sub[take] = dd[take]
    return out


# Run: score every candidate over the chosen domains, write the table
def main():
    doms_arg = None
    for i, a in enumerate(sys.argv):
        if a == "--domains" and i + 1 < len(sys.argv):
            doms_arg = [int(v) for v in sys.argv[i + 1].split(",")]

    gdf = gpd.read_file(DOMAIN_FILE)
    doms = doms_arg or sorted(gdf[DOMAIN_ID_FIELD].astype(int).tolist())

    base = rasterio.open(BASE_DEM)
    bcrs = CRS.from_wkt(base.crs.to_wkt())
    bhz = bcrs.sub_crs_list[0] if bcrs.is_compound else bcrs
    print(f"base: {BASE_DEM.name}  EPSG:{bhz.to_epsg()}  {base.transform.a} m")

    cands = candidates()
    for c in cands:
        c["meta"] = describe(c["tifs"])
        c["srcs"] = [rasterio.open(t) for t in c["tifs"]]
        c["vrts"] = [WarpedVRT(s, crs=base.crs, resampling=Resampling.nearest)
                     for s in c["srcs"]]
        c["cov"] = 0
        c["vals"] = []
    keys = [f"{c['year']}_{c['name'][:18]}" for c in cands]
    print(f"\n{len(cands)} candidates, {len(doms)} domains\n")
    for c, k in zip(cands, keys):
        m = c["meta"]
        print(f"  {k:<26} EPSG:{m['epsg']}  {m['res']} m  "
              f"{m['n_files']} file(s)  vert={m['vertical']}")

    per_domain, total_target = [], 0
    print(f"\n{'dom':>4} {'dryland gap':>12}  " +
          "  ".join(f"{k[:11]:>11}" for k in keys))
    for did in doms:
        row = gdf[gdf[DOMAIN_ID_FIELD].astype(int) == did].iloc[0]
        bx0, by0, bx1, by1 = row.geometry.bounds
        c0 = int(round(bx0 - base.transform.c))
        r0 = int(round(base.transform.f - by1))
        a = base.read(1, window=Window(c0, r0, CLIP[1], CLIP[0]),
                      boundless=True, fill_value=base.nodata).astype(np.float32)
        a = np.where((a == base.nodata) | ~np.isfinite(a), np.nan, a)
        bnds = (bx0, by1 - CLIP[0], bx0 + CLIP[1], by1)

        grids = [read_on_grid(c["vrts"], bnds, CLIP) for c in cands]
        stack = np.stack(grids)
        with np.errstate(all="ignore"):
            consensus = np.nanmedian(stack, axis=0)
        # dry-land gap: 2009 missing, and the CONSENSUS of candidates says land
        target = np.isnan(a) & ~np.isnan(consensus) & (consensus > MHW)
        n = int(target.sum())
        total_target += n

        cells, this_domain = [], {}
        for c, k, gr in zip(cands, keys, grids):
            hit = int((target & ~np.isnan(gr)).sum())
            c["cov"] += hit
            if len(c["vals"]) < 400_000:
                v = gr[np.isnan(a) & ~np.isnan(gr)]
                if v.size:
                    c["vals"].append(v[::37].astype(np.float32))
            cells.append(f"{100 * hit / max(n, 1):10.1f}%")
            # Per-domain count, not the running total (which made every row climb)
            this_domain[f"{k}_cells"] = hit
            this_domain[f"{k}_pct"] = (round(100 * hit / n, 2)
                                       if n else None)
        # No dry-land gap: percentages undefined (None), not 0%
        per_domain.append({"domain": did, "dryland_gap": n, **this_domain})
        print(f"{did:>4} {n:>12,}  " + "  ".join(cells), flush=True)

    base.close()

    rows = []
    print("\n" + "=" * 78)
    print(f"DRY-LAND GAP COVERAGE  (target {total_target:,} cells over "
          f"{len(doms)} domains)")
    for c, k in zip(cands, keys):
        pct = 100 * c["cov"] / max(total_target, 1)
        v = np.concatenate(c["vals"]) if c["vals"] else np.array([])
        flat, topval = 0.0, None
        if v.size:
            u, cnt = np.unique(np.round(v, 3), return_counts=True)
            j = int(np.argmax(cnt))
            flat, topval = cnt[j] / v.size, float(u[j])
        warn = ("  <-- FLAT/VOID-FILLED, top value "
                f"{topval:.3f} in {100 * flat:.1f}% of cells"
                if flat > FLAT_FRACTION_WARN else "")
        print(f"  {pct:6.2f}%   {k:<26}{warn}")
        rows.append({"year": c["year"], "dataset": c["name"],
                     "epsg": c["meta"]["epsg"], "res_m": c["meta"]["res"],
                     "vertical": c["meta"]["vertical"],
                     "dryland_target": total_target, "covered": c["cov"],
                     "covered_pct": round(pct, 2),
                     "top_value": topval,
                     "top_value_frac": round(float(flat), 4),
                     "flagged_flat": bool(flat > FLAT_FRACTION_WARN)})
        for vt in c["vrts"]:
            vt.close()
        for s in c["srcs"]:
            s.close()

    genuine = [r for r in rows if not r["flagged_flat"]]
    good = sorted((r for r in genuine if r["covered_pct"] >= 95),
                  key=lambda r: (r["year"], r["dataset"]))
    print("\n  EARLIEST GENUINE DATASET COVERING >=95% OF DRY-LAND GAPS:")
    if good:
        b = good[0]
        print(f"     {b['year']}  {b['dataset']}  -> {b['covered_pct']:.2f}%")
    else:
        print("     none reaches 95%; see the table above.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "dem_coverage_survey.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(OUT_DIR / "dem_coverage_by_domain.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_domain[0].keys()))
        w.writeheader(); w.writerows(per_domain)
    print(f"\n  csv: {OUT_DIR / 'dem_coverage_survey.csv'}")
    print(f"       {OUT_DIR / 'dem_coverage_by_domain.csv'}")


if __name__ == "__main__":
    main()
