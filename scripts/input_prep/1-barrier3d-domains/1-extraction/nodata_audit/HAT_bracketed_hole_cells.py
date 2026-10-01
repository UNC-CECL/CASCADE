"""
The unsurveyed holes that stop FindWidths with measured land on both sides, cell by cell.

    python scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/HAT_bracketed_hole_cells.py

Reads a dune-topo version's topography and nodata arrays, its picks and the
product's npy-arrays; writes bracketed_hole_cells.csv (interior row, npy cell,
UTM centre per hole cell) to the version's nodata-audit/. Refuses if the saved
arrays are not a plain extraction. Details: scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-01
"""

import contextlib
import csv
import io
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "scripts" / "input_prep" / "1-barrier3d-domains" / "1-extraction"))
import HAT_dune_topo_extractor as ex  # noqa: E402
from site_layer import hat_topo_version as htv  # noqa: E402
from site_layer.hat_elevation_products import product as elevation_product  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TOPO_PRODUCT = "1984-start"
VERSION_OVERRIDE = None

# The DEM each product's npy-arrays were exported from (HAT_export_to_numpy.py)
ELEVATION_FOR_PRODUCT = {"1984-start": "2009-2014-1996", "2004-start": "2009-2014"}

FIRST_GIS, LAST_GIS = 1, 90
DAM_TO_M = 10.0
CELL_SIZE_M = 10.0
WATER_AT_M = 0.0                  # FindWidths stops at the first cell <= this (MHW)

# The surviving reviewed hole geometry, diffed against so moved holes are re-reviewed
REVIEWED_CELLS = "aerial-review/bracketed_hole_cells_v1.csv"

# All outputs in one nodata-audit/ folder beside the extraction; audit_dir() names it
AUDIT_SUBDIR = "nodata-audit"
OUT_NAME = "bracketed_hole_cells.csv"
# -----------------------------------------------------------------------------


# <product>/dune-topo/<version>/nodata-audit/, created on demand
def audit_dir(topo_dir):
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# {domain: (origin_x, origin_y)}, the top-left corner of each 10 m clip
def load_origins(topo_product):
    name = ELEVATION_FOR_PRODUCT.get(topo_product)
    if name is None:
        raise SystemExit(f"\nno elevation product listed for {topo_product!r} "
                         f"in ELEVATION_FOR_PRODUCT\n")
    with elevation_product(name).audit_10m.open() as f:
        return {int(r["domain"]): (float(r["origin_x"]), float(r["origin_y"]))
                for r in csv.DictReader(f)}


# The extractor's processed profiles, dune row and interior for one domain, rebuilt from npy-arrays
def rebuild_domain(npy_path, window):
    with contextlib.redirect_stdout(io.StringIO()):
        dom = ex.load_profiles(npy_path)
        _, dune_loc = ex.find_dunes(dom["z"], dom["start_beach"],
                                    int(window["i0"]), int(window["i1"]))
        topo, _ = ex.build_interior(dom["z"], dune_loc)
        if ex.TRIM_INTERIOR_ROWS:
            topo = ex.remove_water_rows(topo, ex.SENTINEL_WATER_M)
        row0, _ = ex.interior_row0_line(dom["z"], dune_loc)
    return dom, topo, row0


# True if the rebuilt interior is the saved one, values and nodata mask both
def matches_saved(topo_rebuilt, topo_m, nodata):
    nod = topo_rebuilt <= ex.NODATA_SENTINEL_M + 1e-9
    vals = np.where(nod, ex.SENTINEL_WATER_M, topo_rebuilt)
    return (vals.shape == topo_m.shape and np.allclose(vals, topo_m, atol=1e-5)
            and np.array_equal(nod, nodata))


# Per profile: the bracketed hole (rows) and land hidden behind it, or the reason it is not one
def classify_profile(col, nod):
    n = col.size
    water = np.nonzero(col <= WATER_AT_M)[0]
    if not water.size or not nod[water[0]]:
        return "not truncated", None, 0
    r = int(water[0])
    e = r
    while e + 1 < n and nod[e + 1]:
        e += 1
    hidden = int(((col[r + 1:] > WATER_AT_M) & ~nod[r + 1:]).sum())
    if r >= 1 and e + 1 < n and col[e + 1] > WATER_AT_M and not nod[e + 1]:
        return "bracketed", list(range(r, e + 1)), hidden
    return ("hides land, not bracketed" if hidden else "nothing behind"), None, hidden


# Saved interior (row, profile) -> npy-arrays (row, col): undo trim, shear, ocean-first and the flip
def npy_cell(dom, row0, r, p, npy_shape):
    n_along, n_cross = npy_shape
    j = int(row0[p]) + r + int(dom["c0"]) + int(dom["shear"][p])
    if row0[p] < 0 or j >= n_cross or ex.OCEAN_LOC != "right" or not ex.ALONGSHORE_FLIP:
        raise SystemExit(f"\n{dom['name']} profile {p} row {r}: no npy-arrays cell "
                         f"(row0 {row0[p]}, raw column {j}, OCEAN_LOC "
                         f"{ex.OCEAN_LOC!r}); the mapping only covers the "
                         f"right/flipped layout\n")
    return n_along - 1 - p, n_cross - 1 - j


# Which reviewed holes moved, appeared or went: the keys to re-review before trusting old verdicts
def diff_reviewed(rows, reviewed_path):
    if not reviewed_path.is_file():
        print(f"  [review] no {reviewed_path.name}; nothing to diff against")
        return
    old, new = defaultdict(set), defaultdict(set)
    with reviewed_path.open() as f:
        for r in csv.DictReader(f):
            old[(int(r["domain"]), int(r["profile"]))].add(
                (float(r["utm_x"]), float(r["utm_y"])))
    for r in rows:
        new[(r["domain"], r["profile"])].add((r["utm_x"], r["utm_y"]))
    moved = sorted(k for k in old.keys() & new.keys() if old[k] != new[k])
    added = sorted(new.keys() - old.keys())
    gone = sorted(old.keys() - new.keys())
    print(f"  [review] against {reviewed_path.name}: "
          f"{len(old.keys() & new.keys()) - len(moved)} holes unchanged, "
          f"{len(moved)} moved, {len(added)} new, {len(gone)} gone")
    for label, keys in (("moved", moved), ("new", added), ("gone", gone)):
        if keys:
            print(f"    {label}: " + ", ".join(f"D{d}/p{p}" for d, p in keys))


# Run: rebuild and check each domain, collect the bracketed hole cells, write the CSV
def main():
    topo_dir, _, version = htv.topo_dirs(TOPO_PRODUCT, VERSION_OVERRIDE)
    npy_dir = htv.npy_dirs(TOPO_PRODUCT)[0]
    picks_path = htv.picks_dir(TOPO_PRODUCT) / f"HAT_dune_search_windows_{version}.json"
    if not picks_path.is_file():
        raise SystemExit(f"\nno picks for {version}: {picks_path}\n")
    picks = json.loads(picks_path.read_text(encoding="utf-8"))
    origins = load_origins(TOPO_PRODUCT)
    print(f"product {TOPO_PRODUCT}, version {version}, picks {picks_path.name}")

    # Rebuild each domain and refuse any whose saved arrays are not a plain extraction
    rows, counts, hidden = [], defaultdict(int), defaultdict(int)
    for gis in range(FIRST_GIS, LAST_GIS + 1):
        stem = f"domain_{gis}"
        npy_path = npy_dir / f"{stem}.npy"
        dom, topo_rebuilt, row0 = rebuild_domain(npy_path, picks[stem])
        topo_m = np.load(topo_dir / htv.array_name("topography", gis)) * DAM_TO_M
        nodata = np.load(topo_dir / htv.array_name("nodata", gis))
        if not matches_saved(topo_rebuilt, topo_m, nodata):
            raise SystemExit(f"\nD{gis}: the saved {version} topography is not "
                             f"the extractor's output from {npy_path.name} with "
                             f"{picks_path.name}; the npy mapping would be wrong. "
                             f"See the README.\n")
        npy_shape = np.load(npy_path, mmap_mode="r").shape
        ox, oy = origins[gis]

        # Classify every profile; keep the cells of each bracketed hole
        for p in range(topo_m.shape[1]):
            kind, hole, hid = classify_profile(topo_m[:, p], nodata[:, p])
            counts[kind] += 1
            hidden[kind] += hid
            for r in hole or []:
                nr, nc = npy_cell(dom, row0, r, p, npy_shape)
                rows.append(dict(domain=gis, profile=p, interior_row=r,
                                 npy_row=nr, npy_col=nc,
                                 utm_x=ox + CELL_SIZE_M * (nc + 0.5),
                                 utm_y=oy - CELL_SIZE_M * (nr + 0.5)))

    # Write the cells
    out = audit_dir(topo_dir) / OUT_NAME
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["domain", "profile", "interior_row",
                                          "npy_row", "npy_col", "utm_x", "utm_y"])
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}\n")

    # Report the counts and the diff against the reviewed holes
    n_holes = len({(r["domain"], r["profile"]) for r in rows})
    n_trunc = sum(v for k, v in counts.items() if k != "not truncated")
    print(f"  profiles truncated by an unsurveyed cell : {n_trunc:,} of "
          f"{sum(counts.values()):,}")
    print(f"  bracketed by measured land (written)     : {n_holes} holes, "
          f"{len(rows)} cells, {hidden['bracketed'] * CELL_SIZE_M:,.0f} m hidden")
    print(f"  hide land but not bracketed (not written): "
          f"{counts['hides land, not bracketed']} holes, "
          f"{hidden['hides land, not bracketed'] * CELL_SIZE_M:,.0f} m hidden")
    diff_reviewed(rows, htv.extraction_dir(TOPO_PRODUCT) / REVIEWED_CELLS)


if __name__ == "__main__":
    main()
