"""
How far seaward of the model's interior row 0 does a digitized dune line sit, per domain?

    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/1-measurement/HAT_measure_duneline_shift.py
    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/1-measurement/HAT_measure_duneline_shift.py --year 2004
    python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/1-measurement/HAT_measure_duneline_shift.py --against 1997

Measures a year's dune line against interior row 0 of the extraction, per
profile and per domain (or one line against another, which cancels row 0),
and writes the shift tables under the duneline-shift folder. Details: scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

from __future__ import annotations

import argparse
import csv
import importlib.util as _iu
import json
import sys
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from shapely.geometry import LineString, shape
from shapely.ops import transform as sh_transform, unary_union


# Walk up until a directory holds data/hatteras_init
def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root (no data/hatteras_init above me)")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_topo_version import (duneline_shift_dir,  # noqa: E402
                              product_for_year)

# --- CONFIG ------------------------------------------------------------------
OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                 / "1-produce" / "HAT_road_offset_from_dune_start.py")

# The repository copy of the lines wins; the D: drive is only a fallback
from site_layer import hat_topo_version as _tv  # noqa: E402
DUNELINE_DIRS = (
    _tv.DUNELINE_DIR,
    Path(r"D:\Hatteras_GIS\Dunelines"),
)


# First existing duneline_<year>.geojson across DUNELINE_DIRS
def duneline_path(year: int) -> Path:
    tried = []
    for d in DUNELINE_DIRS:
        cand = d / "duneline_{}.geojson".format(year)
        tried.append(cand)
        if cand.is_file():
            return cand
    raise SystemExit("\nno dune line for {}. Looked in:\n  {}\n".format(
        year, "\n  ".join(str(t) for t in tried)))


# The 10 m resampled raster a year's domain is read from
def resampled_tif(year: int, domain: int) -> Path:
    from site_layer.hat_elevation_products import product
    prod = {1984: "2009-2014-1996", 1997: "2009-2014-1996",
            1967: "2009-2014-1996", 2004: "2009-2014"}[int(year)]
    d = product(prod).resampled_10m
    for name in ("resampled_domain_{}_filled.tif".format(domain),
                 "resampled_domain_{}.tif".format(domain)):
        cand = d / name
        if cand.is_file():
            return cand
    raise SystemExit(
        "\nno resampled raster for domain {} in {}\n".format(domain, d))
CELL_M = 10.0
# -----------------------------------------------------------------------------


# HAT_road_offset.py, loaded as a module for its extractor helpers
def load_offset_module():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    mod = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = mod
    spec.loader.exec_module(mod)
    return mod


# The dune line, dissolved and reprojected into the domain grid's CRS
def load_line(path: Path, dst_crs):
    if not path.is_file():
        raise SystemExit("\nno dune line at {}\n".format(path))
    gj = json.load(open(path))
    src_crs = gj.get("crs", {}).get("properties", {}).get("name", "EPSG:26918")
    line = unary_union([shape(f["geometry"]) for f in gj["features"]])
    tr = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    return sh_transform(lambda x, y, z=None: tr.transform(x, y), line)


# Per profile and per domain: the line's cross-shore cell against interior row 0
def measure(ext, line_for_crs, year: int, domains):
    windows = json.load(open(ext.WINDOW_JSON))
    line = None
    bounds = None
    rows, per_profile = [], []

    for D in domains:
        tif = resampled_tif(year, D)
        npy = ext.LOAD_PATH / "domain_{}.npy".format(D)
        if not (tif.is_file() and npy.is_file()):
            continue
        with rasterio.open(tif) as src:
            T, n_rows, n_cols, crs = src.transform, src.shape[0], src.shape[1], src.crs
        if line is None:
            line = line_for_crs(crs)
            bounds = line.bounds

        dom = ext.load_profiles(npy)
        prof = ext.masked_profiles(dom["z"])
        w = windows.get("domain_{}".format(D))
        i0, i1 = ((int(w["i0"]), int(w["i1"])) if w
                  else ext.default_window(prof, dom["start_beach"]))
        _elev, dune_loc = ext.find_dunes(prof, dom["start_beach"], i0, i1)
        row0, _lead = ext.interior_row0_line(prof, dune_loc)

        n_along = min(ext.ALONG_COLS, prof.shape[0])
        hits = []
        for p in range(n_along):
            r = (n_rows - 1) - p if ext.ALONGSHORE_FLIP else p
            if not (0 <= r < n_rows) or row0[p] < 0:
                continue
            y = (T * (0, r + 0.5))[1]
            # Map x of cross-shore cell 0 on this profile, inverting the extractor's chain
            j0 = int(dom["c0"]) + int(dom["shear"][p])
            c_pix = (n_cols - 1) - j0
            if not (0 <= c_pix < n_cols):
                continue
            x0 = (T * (c_pix + 0.5, r + 0.5))[0]

            cut = line.intersection(
                LineString([(bounds[0] - 50, y), (bounds[2] + 50, y)]))
            if cut.is_empty:
                continue
            pts = [cut] if cut.geom_type == "Point" else list(getattr(cut, "geoms", []))
            xs = [g.x for g in pts if g.geom_type == "Point"]
            if not xs:
                continue
            # the cross-shore index grows eastward-to-westward, so x decreases
            ks = [(x0 - x) / CELL_M for x in xs]
            # A meandering line can cross one profile more than once
            k = min(ks, key=lambda v: abs(v - row0[p]))
            hits.append((p, float(k), int(row0[p])))
            per_profile.append({"year": year, "topo_version": ext.VERSION,
                                "domain": D, "profile": p,
                                "duneline_cell": round(float(k), 2),
                                "interior_row0_cell": int(row0[p]),
                                "shift_m": round((row0[p] - k) * CELL_M, 1)})

        if not hits:
            print("  D{:3d}  no dune-line crossing on any profile".format(D))
            continue

        a = np.array([[h[1], h[2]] for h in hits])
        shift = (a[:, 1] - a[:, 0]) * CELL_M
        rows.append({
            "domain": D,
            # Stamp the topography version: the numbers only hold for that extraction
            "topo_version": ext.VERSION,
            "n_profiles": len(hits),
            "shift_m_median": round(float(np.median(shift)), 1),
            "shift_cells_median": round(float(np.median(shift)) / CELL_M, 2),
            "shift_m_p10": round(float(np.percentile(shift, 10)), 1),
            "shift_m_p90": round(float(np.percentile(shift, 90)), 1),
            "shift_m_min": round(float(shift.min()), 1),
            "shift_m_max": round(float(shift.max()), 1),
            "row0_cell_median": round(float(np.median(a[:, 1])), 1),
            "duneline_cell_median": round(float(np.median(a[:, 0])), 2),
        })
        print("  D{:3d}  n={:2d}  shift {:+7.1f} m  (p10 {:+6.1f}, p90 {:+6.1f})  "
              "row0 {:.0f}  duneline cell {:.1f}".format(
                  D, len(hits), np.median(shift), np.percentile(shift, 10),
                  np.percentile(shift, 90), np.median(a[:, 1]), np.median(a[:, 0])))
    return rows, per_profile


# Run: measure the chosen year (or difference two), write the tables
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--year", type=int, default=1984,
                    choices=(1967, 1984, 1996, 1997, 2004))
    ap.add_argument("--domains", default="all")
    ap.add_argument("--against", type=int, default=None,
                    choices=(1967, 1984, 1996, 1997, 2004),
                    help="difference against ANOTHER dune line instead of "
                         "interior row 0. Both lines are measured in the "
                         "SAME frame, so the feature term cancels exactly "
                         "and what is left is pure date. Use it when a "
                         "same-feature pair exists.")
    args = ap.parse_args()

    # 1996 and 1967 are measured in the 1984-start frame, the one being corrected
    product = (product_for_year(args.year)
               if args.year in (1984, 2004) else "1984-start")
    domains = (list(range(1, 91)) if args.domains == "all"
               else [int(x) for x in args.domains.split(",")])

    off = load_offset_module()
    ext = off.load_extractor(product)

    print("\n{} dune line vs interior row 0  |  product {}/{}".format(
        args.year, product, ext.VERSION))
    print("+ = dune line SEAWARD of row 0 = cells to move the island seaward\n")

    path = duneline_path(args.year)
    rows, per_profile = measure(ext, lambda crs: load_line(path, crs),
                                args.year, domains)

    if args.against is not None:
        # Line minus line, in one frame, so row 0 cancels
        other = duneline_path(args.against)
        print("\n--- differencing against the {} dune line ---\n"
              .format(args.against))
        rows_b, _ = measure(ext, lambda crs: load_line(other, crs),
                            args.against, domains)
        b = {r["domain"]: r for r in rows_b}
        kept = []
        for r in rows:
            o = b.get(r["domain"])
            if o is None:
                continue
            # Negated so the stored number reads as retreat (positive = later line landward)
            r["shift_m_median"] = round(r["shift_m_median"] - o["shift_m_median"], 1)
            r["shift_cells_median"] = round(r["shift_m_median"] / CELL_M, 2)
            r["duneline_cell_median_other"] = o["duneline_cell_median"]
            for k in ("shift_m_p10", "shift_m_p90", "shift_m_min", "shift_m_max"):
                r[k] = ""      # a difference of two medians has no such spread
            kept.append(r)
        rows = kept
        print("\n{} -> {} dune-line retreat, {} domains".format(
            args.year, args.against, len(rows)))
    if not rows:
        raise SystemExit("\nnothing measured.\n")

    # NOT symmetric between products - 1984-start lives under 2-domain-reconstruction-1984/
    out_dir = duneline_shift_dir(product)
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / ("duneline_shift_{}.csv".format(args.year)
                     if args.against is None else
                     "duneline_retreat_{}_{}.csv".format(args.year, args.against))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    out_p = out_dir / "duneline_shift_{}_profiles.csv".format(args.year)
    with open(out_p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_profile[0].keys()))
        w.writeheader()
        w.writerows(per_profile)

    s = np.array([r["shift_m_median"] for r in rows])
    print("\n{} domains | island-wide median {:+.1f} m | IQR {:+.1f} to {:+.1f} "
          "| min {:+.1f} max {:+.1f}".format(
              len(rows), np.median(s), np.percentile(s, 25), np.percentile(s, 75),
              s.min(), s.max()))
    print("wrote {}".format(out))
    print("wrote {}".format(out_p))


if __name__ == "__main__":
    main()
