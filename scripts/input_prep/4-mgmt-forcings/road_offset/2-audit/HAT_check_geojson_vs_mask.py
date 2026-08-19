# ==============================================================================
# HAT_check_geojson_vs_mask.py
#
# Does the rasterized road mask actually sit where the source geojson says the
# road is? Read-only QC on the rasterization step alone.
#
# WHY THIS IS SEPARATE FROM THE PLACEMENT FIGURE
#   HAT_method_comparison_figures.py scores modelled-vs-real road position,
#   which folds together three things: rasterization, the alongshore collapse to
#   one scalar, and cell quantization. This script isolates the FIRST one, and it
#   does so in the ORIGINAL raster frame -- no orientation, no alongshore flip,
#   no shear, no trim -- so it is independent of every transform the extractor
#   applies. If this passes, a registration error downstream is not the
#   rasterizer's fault.
#
# WHAT IT CHECKS
#   1. CONTAINMENT  in each profile the geojson crosses, does it fall inside the
#      mask's cell footprint? The mask is a 6 m buffer of that same line with
#      all_touched, so containment should be essentially universal. Anything else
#      means a CRS, snap-raster or extent mismatch.
#   2. OFFSET       signed distance from the centerline to the mask's centre, in
#      metres. Should sit near zero with a spread of well under a cell.
#   3. COVERAGE     profiles where one source has road and the other does not.
#   4. YEAR IDENTITY  do the 1984 and 2004 masks actually differ, and only where
#      NC-12 was relocated? This is the check that the 2004 rasterization used
#      the 2004 line -- a patched driver could silently re-burn 1984.
#
# USAGE
#     python HAT_check_geojson_vs_mask.py
# ==============================================================================

from __future__ import annotations

from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio

PROJECT_ROOT = Path(r"C:\Users\hanna\PycharmProjects\CASCADE")
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
ROADS_ROOT = INIT_ROOT / "4-mgmt-forcing" / "road_offset"

MASK_FMT = ROADS_ROOT / "raster" / "{year}" / "masks" / "domain_{d}_road_{year}.npy"
GEOJSON_FMT = ROADS_ROOT / "raw_offset" / "{year}" / "nc12_{year}.geojson"
TIF_FMT = (INIT_ROOT / "1-barrier3d-domains" / "2009-raw"
           / "2009-domain-clipresample" / "domain_{d}" / "resampled_domain_{d}.tif")

YEARS = [1984, 2004]
DOMAINS = list(range(1, 91))
CELL_SIZE_M = 10.0

# Documented relocation blocks: the 1999 Buxton-Avon work and the 1989 Pea
# Island move. Used only to interpret check 4, never to gate it.
RELOCATION_BLOCKS = [(9, 16), (84, 88)]


def geojson_cols_by_row(line_gdf, tif_path: Path):
    """Centerline position per raster row, in ORIGINAL (row, col) cell units.

    The geojson is EPSG:2264 in US survey feet; the domains are UTM 18N metres,
    so this reprojects first. Vertices are binned to integer rows and averaged,
    which is what a near-north-south line through a north-up grid needs.
    """
    with rasterio.open(tif_path) as src:
        crs, inv, bounds, shape = src.crs, ~src.transform, src.bounds, src.shape

    line = line_gdf.to_crs(crs)
    rows, cols = [], []
    for geom in line.geometry:
        parts = geom.geoms if geom.geom_type.startswith("Multi") else [geom]
        for part in parts:
            xy = np.asarray(part.coords)
            # densify so a long segment crossing the domain still lands in every
            # row it passes through
            for (x0, y0), (x1, y1) in zip(xy[:-1], xy[1:]):
                n = max(2, int(max(abs(x1 - x0), abs(y1 - y0)) / 2.0))
                for x, y in zip(np.linspace(x0, x1, n), np.linspace(y0, y1, n)):
                    if not (bounds.left <= x <= bounds.right
                            and bounds.bottom <= y <= bounds.top):
                        continue
                    c, r = inv * (x, y)
                    if 0 <= r < shape[0] and 0 <= c < shape[1]:
                        rows.append(r)
                        cols.append(c)
    if not rows:
        return {}
    rows = np.floor(np.asarray(rows)).astype(int)
    cols = np.asarray(cols)
    return {int(r): float(cols[rows == r].mean()) for r in np.unique(rows)}


def check_year(year: int) -> dict:
    geo = gpd.read_file(str(GEOJSON_FMT).format(year=year))
    contained = miss = 0
    offsets = []
    geo_no_mask, mask_no_geo = [], []
    per_domain = {}

    for d in DOMAINS:
        tif = Path(str(TIF_FMT).format(d=d))
        mask_path = Path(str(MASK_FMT).format(year=year, d=d))
        if not tif.is_file() or not mask_path.is_file():
            continue
        mask = np.load(mask_path)
        mask = np.isfinite(mask) & (mask > 0)
        line_cols = geojson_cols_by_row(geo, tif)

        d_contained = d_miss = 0
        d_off = []
        for r in range(mask.shape[0]):
            cells = np.flatnonzero(mask[r])
            has_line = r in line_cols
            if not cells.size and not has_line:
                continue
            if has_line and not cells.size:
                geo_no_mask.append((d, r))
                continue
            if cells.size and not has_line:
                mask_no_geo.append((d, r))
                continue
            c = line_cols[r]
            # `c` is a FRACTIONAL column from the inverse affine, where cell k
            # occupies [k, k+1). A cell's centre is therefore k + 0.5, so the
            # mask's centre is cells.mean() + 0.5. Comparing against
            # cells.mean() alone injects a spurious +0.5 cell (+5 m) bias into
            # every profile.
            inside = cells.min() <= c <= cells.max() + 1.0
            d_contained += inside
            d_miss += (not inside)
            d_off.append((c - (cells.mean() + 0.5)) * CELL_SIZE_M)

        contained += d_contained
        miss += d_miss
        offsets.extend(d_off)
        if d_off:
            per_domain[d] = (d_miss, float(np.mean(np.abs(d_off))))

    off = np.asarray(offsets) if offsets else np.array([np.nan])
    return {
        "year": year, "contained": contained, "miss": miss,
        "offset": off, "geo_no_mask": geo_no_mask, "mask_no_geo": mask_no_geo,
        "per_domain": per_domain,
    }


def check_year_identity() -> None:
    """Do the two years' masks differ, and only where the road was relocated?"""
    print("\n" + "=" * 78)
    print("CHECK 4  do the 1984 and 2004 masks actually differ, and where?")
    print("=" * 78)
    differ, same = [], []
    for d in DOMAINS:
        p84 = Path(str(MASK_FMT).format(year=1984, d=d))
        p04 = Path(str(MASK_FMT).format(year=2004, d=d))
        if not (p84.is_file() and p04.is_file()):
            continue
        a = np.load(p84) > 0
        b = np.load(p04) > 0
        if a.shape != b.shape:
            print(f"  [shape mismatch] domain {d}: {a.shape} vs {b.shape}")
            continue
        n_diff = int((a ^ b).sum())
        (differ if n_diff else same).append((d, n_diff))

    print(f"  domains whose masks DIFFER between years : {len(differ)}")
    print(f"  domains whose masks are IDENTICAL        : {len(same)}")
    if not differ:
        print("  *** ALL MASKS IDENTICAL -- the 2004 run almost certainly "
              "re-burned the 1984 geojson. Check ROAD_GEOJSON in the driver.")
        return

    in_block = [(d, n) for d, n in differ
                if any(lo <= d <= hi for lo, hi in RELOCATION_BLOCKS)]
    out_block = [(d, n) for d, n in differ
                 if not any(lo <= d <= hi for lo, hi in RELOCATION_BLOCKS)]
    print(f"  differing INSIDE  known relocation blocks {RELOCATION_BLOCKS}: "
          f"{len(in_block)} -> {[d for d, _ in in_block]}")
    print(f"  differing OUTSIDE them                   : {len(out_block)}"
          + (f" -> {[d for d, _ in out_block]}" if len(out_block) <= 24 else ""))

    # Magnitude is what separates a real relocation from digitizing noise: a
    # moved road changes hundreds of cells, a redrawn line changes a handful.
    for label, group in (("inside blocks", in_block), ("outside blocks", out_block)):
        if not group:
            continue
        n = np.array([c for _, c in group])
        print(f"    {label:14s}: cells changed median {np.median(n):6.0f}, "
              f"min {n.min():4d}, max {n.max():4d}")
    top = sorted(differ, key=lambda t: -t[1])[:8]
    print("  largest differences (domain, cells changed): "
          + ", ".join(f"D{d}:{n}" for d, n in top))


def main() -> None:
    print("=" * 78)
    print("geojson vs rasterized mask -- ORIGINAL raster frame, no transforms")
    print("=" * 78)

    for year in YEARS:
        r = check_year(year)
        total = r["contained"] + r["miss"]
        off = r["offset"]
        print(f"\n--- {year}")
        print(f"  profiles with both line and mask   : {total}")
        if total:
            print(f"  centerline INSIDE mask footprint   : {r['contained']} "
                  f"({r['contained'] / total:.2%})")
            print(f"  centerline OUTSIDE mask footprint  : {r['miss']} "
                  f"({r['miss'] / total:.2%})")
        print(f"  centerline-to-mask-centre offset   : "
              f"median {np.nanmedian(off):+.2f} m, "
              f"p90 |{np.nanpercentile(np.abs(off), 90):.2f}| m, "
              f"max |{np.nanmax(np.abs(off)):.2f}| m")
        print(f"  line present but mask empty        : {len(r['geo_no_mask'])} "
              f"profile(s)")
        print(f"  mask present but line absent       : {len(r['mask_no_geo'])} "
              f"profile(s)")
        worst = sorted(r["per_domain"].items(), key=lambda kv: -kv[1][1])[:5]
        if worst:
            print("  largest mean |offset| by domain    : "
                  + ", ".join(f"D{d}:{m:.1f}m" for d, (_, m) in worst))

    check_year_identity()
    print()


if __name__ == "__main__":
    main()
