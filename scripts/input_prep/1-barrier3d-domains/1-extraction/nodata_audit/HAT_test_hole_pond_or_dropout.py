"""
Are the unsurveyed holes that truncate the model island real ponds, or lidar dropouts?

    python scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/HAT_test_hole_pond_or_dropout.py

Votes each hole with reference A (the NCFMP stamp test) and reference B (the
aerial review), and writes hole_verdicts.csv and per-domain dropout masks. Details: scripts/input_prep/1-barrier3d-domains/1-extraction/nodata_audit/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from scipy import ndimage

REPO = next(
    _p for _p in Path(__file__).resolve().parents
    if (_p / "pyproject.toml").exists())   # 1-extraction/nodata_audit/ since 2026-09-09
sys.path.insert(0, str(REPO / "scripts"))
from site_layer import hat_topo_version as htv  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
TOPO_PRODUCT = "1984-start"
NCFMP_DIR = Path(r"D:\Hatteras_GIS\Elevation\Polygons\2014"
                 r"\2014_NCFMP_NC_DEM_P1_J1437737")
DOMAIN_EPSG = 3725                 # the resampled 10 m rasters
GRID_M = 10.0

# Reference a thresholds
STAMP_MODAL_POND = 0.60            # this share of one repeated value -> stamped
STAMP_MODAL_DROPOUT = 0.25         # below this -> genuine varying ground
STAMP_DEPRESS_M = 0.10             # kept for reporting; not part of the rule
RING_CELLS = 3                     # 30 m ring of measured ground around the hole

# Reference c thresholds
ELONG_MAX = 3.0                    # bbox aspect above this reads as a flight line
FILL_MIN = 0.45                    # blob area / bbox area below this reads sparse

POND, DROPOUT, UNKNOWN = "POND", "DROPOUT", "UNKNOWN"

# All outputs in one nodata-audit/ folder beside the extraction; audit_dir() names it
AUDIT_SUBDIR = "nodata-audit"
# -----------------------------------------------------------------------------


# <product>/dune-topo/<version>/nodata-audit/, created on demand
def audit_dir(topo_dir):
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d


# {(domain, profile)
def load_holes(version_dir):
    p = version_dir / "bracketed_hole_cells.csv"
    if not p.is_file():
        raise SystemExit(f"\nmissing {p}\nRun HAT_plot_island_nodata.py first.\n")
    holes = defaultdict(list)
    with p.open() as f:
        for r in csv.DictReader(f):
            holes[(int(r["domain"]), int(r["profile"]))].append(
                (int(r["npy_row"]), int(r["npy_col"]),
                 float(r["utm_x"]), float(r["utm_y"])))
    return holes


# The NCFMP tiles as one sampler
class NCFMP:

    def __init__(self, folder):
        self.src = [rasterio.open(p) for p in sorted(folder.glob("*C0.tif"))]
        if not self.src:
            raise SystemExit(f"\nno NCFMP tiles under {folder}\n")
        self.tf = Transformer.from_crs(DOMAIN_EPSG, self.src[0].crs.to_epsg(),
                                       always_xy=True)
        print(f"[NCFMP] {len(self.src)} tiles, "
              f"EPSG {self.src[0].crs.to_epsg()}, "
              f"{self.src[0].res[0]:.3f} m")

    # Every NCFMP pixel inside the 10 m cell footprint around each point
    def sample(self, xs, ys, half=GRID_M / 2.0):
        if not len(xs):
            return np.array([])
        tx, ty = self.tf.transform(np.asarray(xs), np.asarray(ys))
        vals = []
        for s in self.src:
            for x, y in zip(tx, ty):
                try:
                    win = rasterio.windows.from_bounds(
                        x - half, y - half, x + half, y + half, s.transform)
                    a = s.read(1, window=win, boundless=False).astype(float)
                except (ValueError, rasterio.errors.WindowError):
                    continue
                if a.size == 0:
                    continue
                if s.nodata is not None:
                    a = a[a != s.nodata]
                a = a[np.isfinite(a) & (a > -1e5)]
                if a.size:
                    vals.append(a)
        return np.concatenate(vals) if vals else np.array([])


# Measured-ground cell centres in a ring around the hole footprint
def ring_points(cells, nodata_raw, origin):
    ox, oy = origin
    have = {(r, c) for r, c, _, _ in cells}
    rows = [r for r, _, _, _ in cells]
    cols = [c for _, c, _, _ in cells]
    xs, ys = [], []
    for r in range(min(rows) - RING_CELLS, max(rows) + RING_CELLS + 1):
        for c in range(min(cols) - RING_CELLS, max(cols) + RING_CELLS + 1):
            if (r, c) in have:
                continue
            if not (0 <= r < nodata_raw.shape[0] and 0 <= c < nodata_raw.shape[1]):
                continue
            if nodata_raw[r, c]:
                continue
            xs.append(ox + GRID_M * (c + 0.5))
            ys.append(oy - GRID_M * (r + 0.5))
    return xs, ys


# NCFMP stamp test -> (verdict, modal fraction, depression)
def verdict_a(hole_v, ring_v):
    h = hole_v[np.isfinite(hole_v)]
    r = ring_v[np.isfinite(ring_v)]
    if h.size < 4:
        return UNKNOWN, np.nan, np.nan
    vals, counts = np.unique(np.round(h, 3), return_counts=True)
    top = int(counts.max())
    frac = float(top / h.size)
    modal = float(vals[np.argmax(counts)])
    dep = float(np.median(r) - modal) if r.size else np.nan
    if frac >= STAMP_MODAL_POND:
        return POND, frac, dep
    if frac <= STAMP_MODAL_DROPOUT:
        return DROPOUT, frac, dep
    return UNKNOWN, frac, dep


# Reference B
def read_aerial(vdir):
    p = vdir / "aerial_1996_conflicts" / "aerial_review.csv"
    if not p.is_file():
        print("[B] no aerial_review.csv - running on A and C only")
        return {}
    got = {}
    for r in csv.DictReader(p.open()):
        col = next((c for c in r if c.startswith("aerial_verdict")), None)
        v = (r[col] or "").strip().upper() if col else ""
        if v in (POND, DROPOUT, "UNCLEAR"):
            got[(int(r["domain"]), int(r["profile"]))] = v
    n = {k: sum(1 for x in got.values() if x == k)
         for k in (POND, DROPOUT, "UNCLEAR")}
    print(f"[B] aerial review: {len(got)} holes judged   "
          f"POND {n[POND]}  DROPOUT {n[DROPOUT]}  UNCLEAR {n['UNCLEAR']}")
    return got


# Shape of the connected nodata blob this hole belongs to
def blob_shape(nodata_raw, cells):
    lab, _ = ndimage.label(nodata_raw, structure=np.ones((3, 3), int))
    r0, c0 = cells[0][0], cells[0][1]
    k = lab[r0, c0]
    if k == 0:
        return UNKNOWN, np.nan, np.nan, 0
    ys, xs = np.nonzero(lab == k)
    h = ys.max() - ys.min() + 1
    w = xs.max() - xs.min() + 1
    area = int(ys.size)
    elong = float(max(h, w) / max(min(h, w), 1))
    fill = float(area / (h * w))
    if elong <= ELONG_MAX and fill >= FILL_MIN:
        return POND, elong, fill, area
    return DROPOUT, elong, fill, area


# Run: both votes per hole, then the verdicts and the dropout masks
def main():
    topo_dir, _, version = htv.topo_dirs(TOPO_PRODUCT)
    vdir = audit_dir(topo_dir)
    arr_dir, _ = htv.npy_dirs(TOPO_PRODUCT)
    holes = load_holes(vdir)
    print(f"product {TOPO_PRODUCT}, version {version}: {len(holes)} holes, "
          f"{sum(len(v) for v in holes.values())} cells")

    origins = {}
    from site_layer.hat_elevation_products import product as _elprod
    with _elprod("2009-2014-1996").audit_10m.open() as f:
        for r in csv.DictReader(f):
            origins[int(r["domain"])] = (float(r["origin_x"]),
                                         float(r["origin_y"]))

    aerial = read_aerial(vdir)
    ncfmp = NCFMP(NCFMP_DIR)
    nodata_by_domain, rows = {}, []
    for (dom, prof), cells in sorted(holes.items()):
        if dom not in nodata_by_domain:
            raw = np.load(arr_dir / f"domain_{dom}.npy").astype(float)
            nodata_by_domain[dom] = raw <= -9.0
        nod = nodata_by_domain[dom]

        hv = ncfmp.sample([c[2] for c in cells], [c[3] for c in cells])
        rx, ry = ring_points(cells, nod, origins[dom])
        rv = ncfmp.sample(rx, ry)
        va, rng, dep = verdict_a(hv, rv)
        vc, elong, fill, area = blob_shape(nod, cells)

        # THE RULE, with reference B folded in
        vb = aerial.get((dom, prof), "")
        if va == vc and va != UNKNOWN:
            final, why = va, "A+C agree"
        elif vb in (POND, DROPOUT):
            final, why = vb, "B breaks tie"
        else:
            final, why = POND, ("B unclear" if vb else "no tiebreak")
        rows.append(dict(domain=dom, profile=prof, cells=len(cells),
                         utm_x=round(np.mean([c[2] for c in cells]), 1),
                         utm_y=round(np.mean([c[3] for c in cells]), 1),
                         ncfmp_verdict=va,
                         ncfmp_modal_frac=round(rng, 3) if np.isfinite(rng) else "",
                         ncfmp_depression_m=round(dep, 3) if np.isfinite(dep) else "",
                         shape_verdict=vc,
                         blob_elongation=round(elong, 2) if np.isfinite(elong) else "",
                         blob_fill=round(fill, 3) if np.isfinite(fill) else "",
                         blob_cells=area,
                         aerial_verdict=vb,
                         decided_by=why,
                         verdict=final))

    out = vdir / "hole_verdicts.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {out}")

    # Dropout mask, cleared cells only
    mdir = vdir / "dropout_mask"
    mdir.mkdir(exist_ok=True)
    cleared = defaultdict(list)
    for r in rows:
        if r["verdict"] == DROPOUT:
            cleared[r["domain"]] += holes[(r["domain"], r["profile"])]
    for dom in sorted(nodata_by_domain):
        m = np.zeros_like(nodata_by_domain[dom])
        for rr, cc, _, _ in cleared.get(dom, []):
            m[rr, cc] = True
        np.save(mdir / f"domain_{dom}.npy", m)
    n_cleared = sum(len(v) for v in cleared.values())
    print(f"wrote {mdir}  ({n_cleared} cells cleared as dropouts, "
          f"in {len(cleared)} domains)")

    # Report
    def tally(key):
        d = defaultdict(int)
        for r in rows:
            d[r[key]] += 1
        return dict(d)

    print(f"\n  NCFMP stamp (A)   : {tally('ncfmp_verdict')}")
    print(f"  blob shape  (C)   : {tally('shape_verdict')}")
    print(f"  agreed verdict    : {tally('verdict')}")
    agree = sum(1 for r in rows
                if r["ncfmp_verdict"] == r["shape_verdict"] != UNKNOWN)
    conflict = sum(1 for r in rows
                   if UNKNOWN not in (r["ncfmp_verdict"], r["shape_verdict"])
                   and r["ncfmp_verdict"] != r["shape_verdict"])
    print(f"  A and C agree     : {agree} of {len(rows)}")
    print(f"  A and C conflict  : {conflict} of {len(rows)}")
    print(f"  aerial (B)        : {tally('aerial_verdict')}")
    print(f"  decided by        : {tally('decided_by')}")

    print("\n  per domain (holes / cleared as dropout):")
    per = defaultdict(lambda: [0, 0])
    for r in rows:
        per[r["domain"]][0] += 1
        per[r["domain"]][1] += (r["verdict"] == DROPOUT)
    for d in sorted(per, key=lambda d: -per[d][0]):
        print(f"    D{d:<3} {per[d][0]:>3} holes   {per[d][1]:>3} dropout")

    fin = [r["ncfmp_modal_frac"] for r in rows if r["ncfmp_modal_frac"] != ""]
    if fin:
        print(f"\n  NCFMP modal fraction under a hole: median "
              f"{np.median(fin):.2f}  (stamped >= {STAMP_MODAL_POND}, "
              f"varying <= {STAMP_MODAL_DROPOUT})")
    el = [r["blob_elongation"] for r in rows if r["blob_elongation"] != ""]
    fl = [r["blob_fill"] for r in rows if r["blob_fill"] != ""]
    if el:
        print(f"  blob elongation: median {np.median(el):.2f} "
              f"(pond-like <= {ELONG_MAX})")
        print(f"  blob fill      : median {np.median(fl):.2f} "
              f"(pond-like >= {FILL_MIN})")


if __name__ == "__main__":
    main()
