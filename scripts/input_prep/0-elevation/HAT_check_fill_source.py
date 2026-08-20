"""
HAT_check_fill_source.py

Pre-flight for a candidate fill source: will this point cloud actually cover the
ground the 2009 DEM is missing, BEFORE spending an hour classifying and filling?

    python HAT_check_fill_source.py <path to .las/.laz file or folder>

Runs in two stages so the cheap answer comes first.

STAGE 1 - HEADERS ONLY, ~1 SECOND
    Point count, CRS, classification (if the header carries it), bounding box
    overlap with the domains, and the z range.

    The z range alone usually decides it. A topographic survey stops at its own
    waterline: 2008 bottomed out at -1.33 m NAVD88, 2009 at its waterline too,
    and neither could see the wet ground behind NC-12 - which is why the 2008
    fill added 14, 3 and 0 cells of 50 in domains 78/79/80. A topobathymetric
    survey goes well below that. If stage 1 reports topo-only, the rest is
    largely a formality: it will fail the same way, in the same place.

STAGE 2 - COVERAGE, MINUTES
    Streams the points once, restricted to the domain footprint, and answers
    the only question that matters for the roadway bug:

        of the 50 alongshore cells in the strip immediately landward of the
        island in domains 78/79/80, how many would gain coverage?

    Reported directly against the 2008 baseline, so the comparison is not left
    to interpretation.

Deliberately uses ALL returns, not a ground-classified subset. This is an upper
bound on what any classifier could keep, and SMRF retained 92-95% of the 2008
returns, so raw coverage is a fair proxy and needs no classification pass first.

Requires: laspy, rasterio, geopandas, numpy
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import geopandas as gpd
import laspy
import rasterio
from rasterio.windows import Window

GIS_ROOT = Path(r"D:\Hatteras_GIS")
DEM = (GIS_ROOT / "Elevation" / "Polygons" / "2009"
       / "usace2009_nc_dem_Job1076020" / "2009_full.tif")
DOMAIN_FILE = GIS_ROOT / "domains.geojson"
DOMAIN_ID_FIELD = "domain_id"

ROAD_DOMAINS = [78, 79, 80]     # the roadways that width-drown at t=0
OFFSETS = (1, 2, 3)             # 10 m cells landward of the island edge

# What the 2008 attempt achieved, for side-by-side comparison
BASELINE_2008 = {78: {1: 14, 2: 13, 3: 11},
                 79: {1: 3, 2: 1, 3: 1},
                 80: {1: 0, 2: 0, 3: 0}}
BASELINE_2008_SOUND_PCT = {78: 19.1, 79: 19.6, 80: 17.7}

# A topo survey bottoms out near its waterline; bathy goes well below.
TOPO_ONLY_Z_MIN = -1.5          # m NAVD88
CHUNK = 5_000_000


def _central4(mask, block=10):
    """
    True where all four central cells of each block x block window are True -
    the condition HAT_dem_resample_clip.py imposes for a 10 m cell to carry a
    value (bilinear at 10x reduces to the central 2x2, all four required).
    """
    h, w = mask.shape
    b = mask.reshape(h // block, block, w // block, block)
    lo, hi = block // 2 - 1, block // 2
    return (b[:, lo, :, lo] & b[:, lo, :, hi] &
            b[:, hi, :, lo] & b[:, hi, :, hi])


def find_files(arg):
    p = Path(arg)
    if p.is_dir():
        f = sorted(list(p.glob("*.laz")) + list(p.glob("*.las")))
        if not f:
            sys.exit(f"no .las/.laz in {p}")
        return f
    if not p.exists():
        sys.exit(f"not found: {p}")
    return [p]


def stage1(files, dom_bounds):
    print("=" * 74)
    print("STAGE 1 - headers")
    print("=" * 74)
    zmin, zmax, total = np.inf, -np.inf, 0
    any_overlap = False
    dminx, dminy, dmaxx, dmaxy = dom_bounds
    for f in files:
        with laspy.open(str(f)) as fh:
            h = fh.header
            total += h.point_count
            zmin = min(zmin, h.mins[2]); zmax = max(zmax, h.maxs[2])
            ov = not (h.maxs[0] < dminx or h.mins[0] > dmaxx or
                      h.maxs[1] < dminy or h.mins[1] > dmaxy)
            any_overlap |= ov
            try:
                crs = h.parse_crs()
                crs_s = f"EPSG:{crs.to_epsg()}" if crs and crs.to_epsg() else (
                    crs.name if crs else "none")
            except Exception:
                crs_s = "unparsed"
            print(f"  {f.name}")
            print(f"    {h.point_count:>12,} pts  fmt {h.point_format.id}  {crs_s}")
            print(f"    x {h.mins[0]:.0f}..{h.maxs[0]:.0f}  "
                  f"y {h.mins[1]:.0f}..{h.maxs[1]:.0f}  "
                  f"z {h.mins[2]:.2f}..{h.maxs[2]:.2f}"
                  f"   domains: {'OVERLAP' if ov else 'no overlap'}")

    print(f"\n  total {total:,} pts,  z {zmin:.2f} .. {zmax:.2f} m")
    if not any_overlap:
        print("\n  STOP: no file overlaps the domain footprint. Wrong area or "
              "wrong CRS.")
        return False, zmin
    if zmin > TOPO_ONLY_Z_MIN:
        print(f"\n  >>> TOPOGRAPHIC ONLY (z stops at {zmin:.2f} m NAVD88).")
        print("      This is the same limitation as 2008 (-1.33) and 2009. It")
        print("      will not see the submerged/saturated ground behind NC-12,")
        print("      so expect it to fail there the way 2008 did. Stage 2 will")
        print("      quantify it, but temper expectations.")
    else:
        print(f"\n  >>> Reaches {zmin:.2f} m NAVD88, below the "
              f"{TOPO_ONLY_Z_MIN:.1f} m topo cutoff - consistent with")
        print("      topobathymetric coverage. This is what could actually fix "
              "the roadway bug.")
    return True, zmin


def stage2(files, gdf, dom_bounds):
    print("\n" + "=" * 74)
    print("STAGE 2 - coverage of the strip landward of NC-12")
    print("=" * 74)
    minx, miny, maxx, maxy = dom_bounds

    xs, ys = [], []
    for f in files:
        with laspy.open(str(f)) as fh:
            h = fh.header
            if (h.maxs[0] < minx or h.mins[0] > maxx or
                    h.maxs[1] < miny or h.mins[1] > maxy):
                continue
            print(f"  scanning {f.name} ({h.point_count:,} pts)...", flush=True)
            for ch in fh.chunk_iterator(CHUNK):
                x = np.asarray(ch.x); y = np.asarray(ch.y)
                m = (x >= minx) & (x <= maxx) & (y >= miny) & (y <= maxy)
                if m.any():
                    xs.append(x[m]); ys.append(y[m])
    if not xs:
        print("  no points inside the domain footprint.")
        return
    px = np.concatenate(xs); py = np.concatenate(ys)
    print(f"  {len(px):,} returns inside the domain footprint\n")

    src = rasterio.open(DEM); t = src.transform
    print(f"  {'dom':>4} {'sound nodata':>13} {'w/ returns':>11} "
          f"{'(2008 was)':>11}")
    sound = {}
    for did in ROAD_DOMAINS:
        row = gdf[gdf[DOMAIN_ID_FIELD].astype(int) == did].iloc[0]
        bx0, by0, bx1, by1 = row.geometry.bounds
        c0 = int(round(bx0 - t.c)); r0 = int(round(t.f - by1))
        win = Window(c0, r0, 2000, 500)
        a = src.read(1, window=win, boundless=True,
                     fill_value=src.nodata).astype(float)
        a = np.where((a == src.nodata) | ~np.isfinite(a), np.nan, a)
        tw = src.window_transform(win)

        hit = np.zeros(a.shape, bool)
        col = ((px - tw.c) / 1.0).astype(np.int64)
        rw = ((tw.f - py) / 1.0).astype(np.int64)
        ok = (rw >= 0) & (rw < a.shape[0]) & (col >= 0) & (col < a.shape[1])
        hit[rw[ok], col[ok]] = True

        nodata = np.isnan(a)
        meas = ~nodata
        west = np.zeros_like(nodata)
        for i in range(a.shape[0]):
            m = np.flatnonzero(meas[i])
            if len(m):
                west[i, :m[0]] = True
        n_nd = int((nodata & west).sum())
        n_hit = int((nodata & west & hit).sum())
        pct = 100 * n_hit / max(n_nd, 1)
        sound[did] = (n_nd, n_hit, pct)
        print(f"  {did:>4} {n_nd:>13,} {n_hit:>7,} ({pct:4.1f}%) "
              f"{BASELINE_2008_SOUND_PCT[did]:>10.1f}%")

    # The headline: the 10 m strip the drowning test actually reads
    print(f"\n  10 m cells gained in the strip landward of the island "
          f"(of 50 alongshore):")
    print(f"  {'dom':>4} {'offset':>7} {'this dataset':>13} {'2008':>6}")
    total_new, total_old = 0, 0
    for did in ROAD_DOMAINS:
        row = gdf[gdf[DOMAIN_ID_FIELD].astype(int) == did].iloc[0]
        bx0, by0, bx1, by1 = row.geometry.bounds
        c0 = int(round(bx0 - t.c)); r0 = int(round(t.f - by1))
        win = Window(c0, r0, 2000, 500)
        a = src.read(1, window=win, boundless=True,
                     fill_value=src.nodata).astype(float)
        a = np.where((a == src.nodata) | ~np.isfinite(a), np.nan, a)
        tw = src.window_transform(win)
        hit = np.zeros(a.shape, bool)
        col = ((px - tw.c) / 1.0).astype(np.int64)
        rw = ((tw.f - py) / 1.0).astype(np.int64)
        ok = (rw >= 0) & (rw < a.shape[0]) & (col >= 0) & (col < a.shape[1])
        hit[rw[ok], col[ok]] = True
        meas = ~np.isnan(a)
        # Collapse to the 10 m grid using the SAME rule the pipeline uses.
        # HAT_dem_resample_clip.py resamples by bilinear, which at a 10x
        # reduction is the mean of the CENTRAL 2x2 source cells, and
        # BILINEAR_REQUIRE_ALL_FOUR means all four must be valid or the output
        # is nodata. So a 10 m cell is gained only if all four central 1 m
        # cells have a return - not if any of the 100 does.
        #
        # Validated: the generous "any of 100" rule predicted 50/44/47 for the
        # 2008 cloud where the pipeline actually delivered 14/3/0. Any-of-100
        # is not a usable proxy.
        meas10 = _central4(meas)
        hit10 = _central4(hit)
        for off in OFFSETS:
            gained = 0
            for i in range(50):
                m = np.flatnonzero(meas10[i])
                if not len(m):
                    continue
                c = m[0] - off
                if c >= 0 and hit10[i, c]:
                    gained += 1
            old = BASELINE_2008[did][off]
            if off == 1:
                total_new += gained; total_old += old
            print(f"  {did:>4} {-off:>7} {gained:>13} {old:>6}")
    src.close()

    print("\n" + "=" * 74)
    print(f"  VERDICT at offset -1, summed over domains 78/79/80:")
    print(f"     this dataset: {total_new} of 150 cells")
    print(f"     2008        : {total_old} of 150 cells")
    if total_new > total_old * 2 and total_new >= 30:
        print("\n  Materially better coverage. Worth running the full pipeline.")
    elif total_new > total_old:
        print("\n  Some improvement, but check whether it clears the >20% "
              "threshold\n  in bulldoze before assuming it fixes the roadways.")
    else:
        print("\n  No better than 2008 in the strip that matters. Filling from "
              "this\n  source will not fix the roadway drowning - the fix has "
              "to come from\n  excluding never-surveyed cells in "
              "roadway_manager.bulldoze instead.")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip().splitlines()[4])
    files = find_files(sys.argv[1])
    gdf = gpd.read_file(DOMAIN_FILE)
    sel = gdf[gdf[DOMAIN_ID_FIELD].astype(int).isin(ROAD_DOMAINS)]
    dom_bounds = tuple(sel.total_bounds)
    print(f"\ncandidate fill source: {sys.argv[1]}")
    print(f"checking against domains {ROAD_DOMAINS}, "
          f"bbox x {dom_bounds[0]:.0f}..{dom_bounds[2]:.0f} "
          f"y {dom_bounds[1]:.0f}..{dom_bounds[3]:.0f}\n")
    ok, _ = stage1(files, dom_bounds)
    if ok:
        stage2(files, gdf, dom_bounds)


if __name__ == "__main__":
    main()
