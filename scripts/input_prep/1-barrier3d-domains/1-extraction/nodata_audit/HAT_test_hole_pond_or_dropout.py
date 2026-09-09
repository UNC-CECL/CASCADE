"""
HAT_test_hole_pond_or_dropout.py

Are the unsurveyed holes that truncate the model island real ponds, or lidar
dropouts?

THE QUESTION AND WHY IT IS NOT OBVIOUS
---------------------------------------
99 profiles in 1984-start have their island truncated by an unsurveyed cell:
barrier3d.FindWidths stops at the first cell at or below sea level, and an
unsurveyed cell is written to the -3.0 m water sentinel, so the scan stops there
and every measured cell behind it is discarded. 731 cells are involved.

Whether that is wrong depends on what those cells are:

    a POND      the survey was right to return nothing, water is water, and the
                truncation is the island's real shape. Change nothing.
    a DROPOUT   the survey failed over dry ground, and the model is running on
                an island several hundred metres narrower than the data.

An earlier diagnostic claimed the holes are dropouts because they are bracketed
by measured land at ~0.9 m MHW on both sides. That argument is wrong and is
recorded here so it is not made again: a pond is BY DEFINITION surrounded by dry
ground, so "dry on both sides" is equally consistent with either. The elevation
test discriminates a survey stopping at a false shoreline from a real bay
margin, which is a different question about a different set of cells.

Evidence pointing the other way, which any result has to beat: a cell is
unsurveyed in this product only if 1996 ALACE, 2009 USACE AND 2014 NOAA
Post-Sandy all failed at it. Three independent surveys over 18 years failing at
one spot is what persistent water looks like. "Mostly ponds, change nothing" is
a live outcome, not a failed test.

TWO INDEPENDENT REFERENCES, BOTH MUST AGREE TO OVERTURN
--------------------------------------------------------
A  2014 NCFMP hydro-flattening stamp.
   NCFMP is DISQUALIFIED as an elevation source - 94.6% of its coverage in the
   gap is two stamped constants, -0.762 and -0.914 m, i.e. -2.5 and -3.0 ft.
   That is exactly what makes it authoritative here. A hydro-flattening
   compiler delineates water-body polygons and stamps a flat surface inside
   them, so the constant IS a water classification, made independently of
   whether the lidar returned anything. This reads its water mask, never its
   elevation, so its "unknown" vertical datum and foot-derived values do not
   matter: the test is whether the value is CONSTANT, not what it means.

C  Blob shape.
   Marsh ponds are compact and roughly convex. Lidar dropouts follow flight
   lines and scan geometry, so they come out elongated and sparse in their
   bounding box. Computed on the raw nodata mask, which no other reference
   touches, and independent of time - which matters because A is 2014 standing
   in for 1996.

A third reference, the 1996 aerial frames, was considered and dropped: they are
scanned historical imagery with per-frame exposure variation, so classifying
them needs either a manual pass over 99 chips or a threshold that would not
survive review. Dropping it leaves the temporal assumption resting on A alone.
That is a real weakness of this design and belongs in the methods paragraph.

THE DECISION RULE IS DELIBERATELY ASYMMETRIC
---------------------------------------------
    A and C agree             ->  their verdict
    A and C conflict, B votes ->  B decides
    otherwise                 ->  WATER, and the DEM is left alone

B is the 1996 aerial review from HAT_hole_aerial_picker.py. Chips were rendered
only where A and C conflict, so B exists exactly where the automated pair
cancels out - which makes "2 of 3 agree" and "B breaks the tie" the same rule.
B is also the only CONTEMPORANEOUS reference, so where it has an opinion it
deserves to carry the decision rather than be outvoted by two proxies.

Unknown, absent or conflicting evidence defaults to water, which is the
current behaviour. The DEM changes only where there is affirmative evidence it
is wrong. The cost of that choice is that the island stays too narrow wherever
the test cannot resolve a hole - a known, stated, one-directional bias.

INPUT   dune-topo/<version>/bracketed_hole_cells.csv   from HAT_plot_island_nodata
        <product>/npy-arrays/domain_<N>.npy            raw nodata, for shape
        D:\\Hatteras_GIS\\Elevation\\Polygons\\2014\\2014_NCFMP_*\\*.tif

OUTPUT  dune-topo/<version>/hole_verdicts.csv          per hole, both votes
        dune-topo/<version>/dropout_mask/domain_<N>.npy  bool, cleared cells only
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import rasterio
from pyproj import Transformer
from scipy import ndimage

REPO = Path(__file__).resolve().parents[5]   # 1-extraction/nodata_audit/ since 2026-09-09
sys.path.insert(0, str(REPO / "scripts"))
import hat_topo_version as htv  # noqa: E402

TOPO_PRODUCT = "1984-start"
NCFMP_DIR = Path(r"D:\Hatteras_GIS\Elevation\Polygons\2014"
                 r"\2014_NCFMP_NC_DEM_P1_J1437737")
DOMAIN_EPSG = 3725                 # the resampled 10 m rasters
GRID_M = 10.0

# --- reference A thresholds ---------------------------------------------
STAMP_MODAL_POND = 0.60            # this share of one repeated value -> stamped
STAMP_MODAL_DROPOUT = 0.25         # below this -> genuine varying ground
STAMP_DEPRESS_M = 0.10             # kept for reporting; not part of the rule
RING_CELLS = 3                     # 30 m ring of measured ground around the hole

# --- reference C thresholds ---------------------------------------------
ELONG_MAX = 3.0                    # bbox aspect above this reads as a flight line
FILL_MIN = 0.45                    # blob area / bbox area below this reads sparse

POND, DROPOUT, UNKNOWN = "POND", "DROPOUT", "UNKNOWN"

# Every output of this folder lands under one directory beside the extraction it
# describes, rather than being scattered through the run folder it did not
# produce. audit_dir() is the only place that name is spelled.
AUDIT_SUBDIR = "nodata-audit"


def audit_dir(topo_dir):
    """<product>/dune-topo/<version>/nodata-audit/, created on demand."""
    d = topo_dir.parent / AUDIT_SUBDIR
    d.mkdir(parents=True, exist_ok=True)
    return d



def load_holes(version_dir):
    """{(domain, profile): [(npy_row, npy_col, x, y), ...]} from the cell CSV."""
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


class NCFMP:
    """The NCFMP tiles as one sampler. Returns nan where no tile covers a point."""

    def __init__(self, folder):
        self.src = [rasterio.open(p) for p in sorted(folder.glob("*C0.tif"))]
        if not self.src:
            raise SystemExit(f"\nno NCFMP tiles under {folder}\n")
        self.tf = Transformer.from_crs(DOMAIN_EPSG, self.src[0].crs.to_epsg(),
                                       always_xy=True)
        print(f"[NCFMP] {len(self.src)} tiles, "
              f"EPSG {self.src[0].crs.to_epsg()}, "
              f"{self.src[0].res[0]:.3f} m")

    def sample(self, xs, ys, half=GRID_M / 2.0):
        """Every NCFMP pixel inside the 10 m cell footprint around each point.

        Point-sampling one NCFMP pixel per 10 m cell was the first version and
        it broke the flatness test: NCFMP is 1.52 m, so a 10 m cell contains
        ~43 of its pixels, and taking one made a single-cell hole trivially
        "flat" with a range of exactly zero. 45 of 99 holes returned UNKNOWN
        for that reason alone. Reading the footprint measures flatness WITHIN
        one 10 m cell, which is what the stamp test actually needs.
        """
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


def ring_points(cells, nodata_raw, origin):
    """Measured-ground cell centres in a ring around the hole footprint.

    Excludes the hole itself and every other unsurveyed cell, so the ring is
    ground some survey actually saw. Without that the ring could be more
    no-data, and the depression test would compare two stamps.
    """
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


def verdict_a(hole_v, ring_v):
    """NCFMP stamp test -> (verdict, modal fraction, depression).

    The statistic is the MODAL FRACTION - what share of the NCFMP pixels under
    the hole carry a single repeated value. That is deliberately the same
    metric HAT_survey_dem_coverage.py uses to disqualify a DEM as void-filled
    ("genuine surveys show their most common value in ~0.2% of cells"), applied
    here at one hole instead of island-wide.

    It replaced a max-minus-min range test, which failed in both directions and
    is worth recording:
      * one NCFMP pixel per 10 m cell made a single-cell hole trivially flat,
        range exactly 0, and returned UNKNOWN for 45 of 99 holes;
      * the whole footprint fixed that but picked up the EDGE of a stamped
        water polygon, where values vary by construction, so a pond-edge cell
        read as varying and therefore as a dropout.
    A modal fraction survives both: a pond cell is mostly one stamped value
    even when its footprint clips the polygon edge, and dry ground is not.
    """
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


def read_aerial(vdir):
    """Reference B: {(domain, profile): POND|DROPOUT|UNCLEAR} from the review.

    Absent file, or an unreviewed row, simply yields no vote - which the
    decision rule below treats as no evidence, not as evidence of water.
    """
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


def blob_shape(nodata_raw, cells):
    """Shape of the connected nodata blob this hole belongs to.

    8-connectivity, per domain. A blob straddling a domain seam is measured
    only within its own 500 m tile; that under-measures a few blobs and is
    noted rather than corrected, because the domains are the unit everything
    else here is expressed in.
    """
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


def main():
    topo_dir, _, version = htv.topo_dirs(TOPO_PRODUCT)
    vdir = audit_dir(topo_dir)
    arr_dir, _ = htv.npy_dirs(TOPO_PRODUCT)
    holes = load_holes(vdir)
    print(f"product {TOPO_PRODUCT}, version {version}: {len(holes)} holes, "
          f"{sum(len(v) for v in holes.values())} cells")

    origins = {}
    with (REPO / "data" / "hatteras_init" / "0-elevation" / "2009-2014-1996"
          / "2-resampled-10m" / "resample_audit.csv").open() as f:
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

        # THE RULE, with reference B folded in.
        #
        # Chips were rendered only where A and C conflict, so B exists exactly
        # where the automated pair cancels out. That makes "2 of 3 agree" and
        # "B breaks the tie" the same rule, not two:
        #
        #     A and C agree            -> their verdict
        #     A and C conflict, B votes -> B decides (B + one of A/C = 2 of 3)
        #     otherwise                 -> POND, the conservative default
        #
        # UNCLEAR is not a vote for water. It is the absence of a vote, and it
        # falls through to the default for the same reason anything unresolved
        # does: the DEM changes only on affirmative evidence that it is wrong.
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

    # --- dropout mask, cleared cells only ------------------------------------
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

    # --- report --------------------------------------------------------------
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
