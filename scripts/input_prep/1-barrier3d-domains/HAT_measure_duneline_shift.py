r"""
HAT_measure_duneline_shift.py
==============================================================================
How far seaward of the model's interior row 0 does a digitized dune line sit?

WHY
---
The 1984-start topography is a 1996 ALACE beach on a 2009 backdune, and the
NC-12 lines it is measured against are 1984-vintage. Where the island migrated
far enough between 1984 and 1996, the 1984 roadbed ends up SEAWARD of the 1996
dune crest and the setback goes negative (GIS 85: -10 m, floored to 0, which
makes roadway_manager relocate the road in year 1). This measures the offset
that gap represents, so HAT_insert_seaward_rows.py can act on a number rather
than on a target.

METHOD -- and why it is exact rather than approximate
-----------------------------------------------------
A domain clip is north-up, so one extractor "profile" is a raster row of
constant y. Intersect the dune-line geometry with that horizontal line, convert
the crossing's easting to a cross-shore cell index, and difference it against
interior row 0 for the same profile.

The frame comes from the extractor itself -- its own c0 and per-profile shear,
through the same inversion cell_to_map documents -- so this is measured in the
frame CASCADE indexes, not in a re-derived one. That distinction is the whole
reason the legacy RoadSetback numbers and the dunestart numbers disagree by a
median 38 m: the legacy pass measured raw and ocean-first, unstraightened.

SIGN: positive = the dune line lies SEAWARD of interior row 0, i.e. the number
of cells the island has to move seaward for row 0 to land on the digitized line.

THE CONTROL THAT MAKES THE 1984 NUMBER READABLE
-----------------------------------------------
Run it on the 2004 pair (2004 dune line against 2004-start, whose DEM is 2009)
and the stable mid-island comes out at +1.4 / -4.1 / -3.6 m on GIS 40/50/60 --
under half a cell. So a digitized dune line and the extractor's interior row 0
are the SAME feature to within the grid, and a large 1984 number is a date
difference rather than a definitional one.

That control matters because it settles a confound the road-offset work had
recorded as unresolvable without a same-year DEM. It is resolvable without one:
the 2004 line is close enough in date to the 2009 surface that its residual IS
the feature term, and the feature term is ~0.

Values well above zero at GIS 10/11/84/85 in the 2004 pass are NOT method error
-- those are the relocation blocks, where five years of hotspot erosion is real.
Read the mid-island domains for the method check.

INPUT   D:\Hatteras_GIS\Dunelines\duneline_<year>.geojson   (EPSG:26918)
        domain-clips-1m/domain_<N>/resampled_domain_<N>.tif (EPSG:3725)
        the extractor's own picks for the resolved version

OUTPUT  hat_topo_version.duneline_shift_dir(<product>)/duneline_shift_<year>.csv
        (1984-start: .../1984-start/row-insert-scope/duneline-shift/)
            one row per domain: median/p10/p90 shift, row 0, dune-line cell

USAGE
-----
    python HAT_measure_duneline_shift.py                  # 1984, all domains
    python HAT_measure_duneline_shift.py --year 2004      # the control
    python HAT_measure_duneline_shift.py --domains 84,85,86
==============================================================================
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


def _find_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit("could not locate the project root (no data/hatteras_init above me)")


REPO = _find_root(Path(__file__).resolve())
sys.path.insert(0, str(REPO / "scripts"))

from hat_topo_version import (duneline_shift_dir,  # noqa: E402
                              product_for_year)

OFFSET_SCRIPT = (REPO / "scripts" / "input_prep" / "4-mgmt-forcings" / "road_offset"
                 / "1-produce" / "HAT_road_offset_from_dune_start.py")

# The REPO copy wins. The lines used to be read straight off D:\Hatteras_GIS,
# which is not version-controlled, not present on another machine, and not
# something a run can record the state of. Hannah placed a curated set inside the
# repo on 2026-09-02; that is now the source, and the external drive is only a
# fallback so older invocations keep working.
DUNELINE_DIRS = (
    REPO / "data" / "hatteras_init" / "1-barrier3d-domains" / "raw-duneline-geojson",
    Path(r"D:\Hatteras_GIS\Dunelines"),
)


def duneline_path(year: int) -> Path:
    """First existing duneline_<year>.geojson across DUNELINE_DIRS."""
    tried = []
    for d in DUNELINE_DIRS:
        cand = d / "duneline_{}.geojson".format(year)
        tried.append(cand)
        if cand.is_file():
            return cand
    raise SystemExit("\nno dune line for {}. Looked in:\n  {}\n".format(
        year, "\n  ".join(str(t) for t in tried)))


# GEOREFERENCE AGAINST THE PRODUCT'S OWN RASTER, not the clip tree.
#
# This was wrong until 2026-09-03. It used
#   1-barrier3d-domains/domain-clips-1m/domain_<N>/resampled_domain_<N>.tif
# which is a DIFFERENT DEM product from the one the 1984 npy arrays were
# exported from: at GIS 85 its origin is offset +0.50 m in easting and its
# elevations differ (max 6.80 m against 6.19 m). Only the transform is used
# here, so the cost was 0.05 cells -- N at GIS 85 went 5.87 instead of 5.82 and
# rounded to 6 either way -- but mixing the extractor's frame with another
# grid's transform is not a thing to leave in place.
#
# Resolved through hat_elevation_products so the period -> product pairing is
# the same single definition every other reader uses.
def resampled_tif(year: int, domain: int) -> Path:
    from hat_elevation_products import product
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


def load_offset_module():
    spec = _iu.spec_from_file_location("hat_off", OFFSET_SCRIPT)
    mod = _iu.module_from_spec(spec)
    sys.modules["hat_off"] = mod
    spec.loader.exec_module(mod)
    return mod


def load_line(path: Path, dst_crs):
    """The dune line, dissolved and reprojected into the domain grid's CRS.

    Reprojected through pyproj from the file's own declared CRS. NAD83 and
    NAD83(NSRS2007) differ by centimetres here, but the transform is done rather
    than assumed away -- the same rule the 1996 aerial chips follow.
    """
    if not path.is_file():
        raise SystemExit("\nno dune line at {}\n".format(path))
    gj = json.load(open(path))
    src_crs = gj.get("crs", {}).get("properties", {}).get("name", "EPSG:26918")
    line = unary_union([shape(f["geometry"]) for f in gj["features"]])
    tr = Transformer.from_crs(src_crs, dst_crs, always_xy=True)
    return sh_transform(lambda x, y, z=None: tr.transform(x, y), line)


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
            # map x of cross-shore cell 0 on this profile, inverting the
            # extractor's chain exactly as cell_to_map does
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
            # a meandering line can cross one profile more than once; take the
            # crossing nearest row 0 rather than the first, which would silently
            # pick a sound-side meander on the wide domains
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
            # STAMP THE TOPOGRAPHY. These numbers are measured against interior
            # row 0, so they are only valid for the extraction that produced it.
            # Without the stamp a v1-era duneline_shift_1984.csv sat unmarked
            # beside v3 files -- 65.9 m against the correct 75.0 m at GIS 85 --
            # and nothing on disk said which was which.
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

    # 1996 and 1967 are not hindcast period starts, so they have no product
    # of their own. They are measured in the frame of the period being
    # corrected -- 1984-start -- which is also the only frame in which a
    # difference against the 1984 line is meaningful.
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
        # LINE MINUS LINE, in one frame. Each line is first measured against the
        # same interior row 0, then the two are differenced -- so row 0 drops out
        # algebraically and no assumption about it survives into the answer. That
        # is the whole point: the row-0 reference is the part that cannot be
        # validated, and differencing removes it rather than bounding it.
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
            # shift = row0 - line, so (row0 - lineEARLY) - (row0 - lineLATE)
            # = lineLATE - lineEARLY. Negated here so the stored number reads as
            # RETREAT: positive = the later line is LANDWARD of the earlier one,
            # i.e. the dune moved landward by that many metres. Storing the raw
            # difference would put a negative sign on ordinary erosion and invert
            # every consumer that reuses this file expecting the shift_m_median
            # convention of the un-differenced output.
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

    # NOT symmetric between products - 1984-start lives under
    # row-insert-scope/. hat_topo_version.duneline_shift_dir owns that.
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
