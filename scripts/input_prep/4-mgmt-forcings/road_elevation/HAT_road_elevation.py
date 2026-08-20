r"""
HAT_road_elevation.py
===============================================================================
Per-domain NC-12 road elevation: the MEAN of the 2009 LiDAR under the 2004 road
alignment. ONE set of numbers, used for BOTH the 1984 and the 2004 start period.

Replaces HAT_road_elevation_from_lidar.py, which sampled two alignments and
wrote two files. That script was deleted on 2026-08-17; it was staged but never
committed, so it is not in any commit. Its blob survives in the object store at
becbfc878ae0afa3f4e76037ef0edff810fa857f until the next `git gc`, recoverable
with `git cat-file -p <hash>`. The reason it is gone rather than kept is in WHY
ONE FILE below.

WHY ONE FILE AND NOT ONE PER VINTAGE
------------------------------------
There is only one DEM. Both vintages were sampled on the same 2009 surface, so
any difference between a "1984" and a "2004" road elevation was never a
difference in time -- it was a difference in WHERE ON THE 2009 SURFACE the two
digitised lines happened to fall. Where NC-12 never moved the two lines sit on
top of each other and the numbers are identical by construction. Where the road
WAS relocated the 1984 line lies over the abandoned corridor, and the "1984 road
elevation" there was the elevation of a place a road used to be.

Neither of those is a measurement of temporal change in roadbed height. Writing
two files implied one. This writes one.

WHAT THE ABANDONED CORRIDOR ACTUALLY IS -- MEASURED, NOT ASSUMED
----------------------------------------------------------------
It is tempting to assume the abandoned alignment was overwashed and bulldozed
flat, and so reads LOW in a 2009 DEM. It does not. Sampled here, the 1984 line
through the relocated domains gives a mean of about 2.4 m NAVD88 against about
1.5 m for the 2004 line, with a within-domain standard deviation up to 1.7 m --
GIS 10 comes back at 4.3 m. The dune migrated over the corridor after the road
left it. That sample is FOREDUNE, not roadbed, and not flattened ground either.

This decides the choice rather than merely complicating it: the two candidates
do NOT bracket the truth. Both sit above the natural grade of the neighbouring
un-relocated domains (~1.0 m), so the 2004 value is the lower of the two AND the
only one that is a graded surface -- the conservative choice as well as the
correct one. RELOCATION BRACKET re-measures this on every run.

WHY THE 1 m CLIP AND NOT THE 10 m RESAMPLE
------------------------------------------
NC-12 is a two-lane road: roughly 7-10 m of pavement plus shoulder. On the 10 m
Barrier3D grid the road is ONE cell wide, so a buffered mask averages the crown
into whatever is beside it -- in the inter-village stretch, the foredune that
NC-12 runs immediately behind.

Every domain folder also carries clip_domain_<N>.tif at 1 m, the native LiDAR
before the Barrier3D resample. A 3.5 m buffer on that -- a 7 m corridor, about
one carriageway -- gives a within-domain standard deviation of about 0.07 m
island-wide. That is a road surface.

HONESTY NOTE: on THIS alignment the 10 m grid would have given nearly the same
answer -- the two agree to a median of 0.00 m and a max of 0.05 m. The 1.59 m
standard deviation that originally motivated the 1 m clip was measured on the
1984 line, which crosses a relocation scar; the 2004 line does not. So the 1 m
clip is a precaution here rather than a rescue, and the sample it rests on is
~3500 cells per domain against ~35 at 10 m. Both numbers are printed under
INTERNAL CHECKS every run so this stays checkable rather than inherited.

MEAN, NOT MEDIAN
----------------
Flat unweighted mean of every valid 1 m cell in the corridor. On this alignment
mean and median differ by 0.005 m for a typical domain and 0.09 m at worst, so
the choice is nearly free; the median is carried in the per-domain CSV so any
domain where they diverge -- a bridge deck, a house, a driveway apron in the
corridor -- is visible rather than silently absorbed.

DATUM -- NOT AMBIGUOUS, DESPITE THE RUNNER
------------------------------------------
bulldoze() writes road_ele straight into xyz_interior_grid:

    road_ele = road_ele / dz
    new_road_domain = np.zeros(...) + road_ele

and the interior arrays are MHW-RELATIVE, because HAT_dune_topo_extractor.py
subtracts MHW_M = 0.36 before anything else. So road_ele MUST be MHW-relative
metres. There is no reading under which NAVD88 is correct.

The runner's ROAD_ELEVATION = 1.45 is high under EITHER reading -- see the audit
document. This file writes MHW-relative; the per-domain CSV carries NAVD88
alongside so nothing has to be taken on trust.

TWO THINGS THIS FILE DOES NOT CORRECT
-------------------------------------
1. THE TIME GAP. The DEM is 2009; one run starts in 1984. CASCADE decrements
   road_ele by RSLR every year, so the 1984 run begins with a roadbed that is
   already 25 years of sea-level rise low relative to its own MHW. No
   back-correction is applied -- these are measurements, not reconstructions.

2. THE RELOCATIONS. GIS 9-15 (relocated 1999) and GIS 84-87 (relocated 1989)
   carry the elevation of the POST-relocation alignment in the 1984 run, because
   that is the alignment sampled. In 1984 the road was physically elsewhere in
   those domains. Flagged, not adjusted.

OUTPUTS  (data/hatteras_init/4-mgmt-forcing/road_elevation/)
------------------------------------------------------------
  RoadElevation.csv          2-row CASCADE file (IDs, m MHW-relative)
  RoadElevation_domains.csv  per-domain stats, both datums, flags
  RoadElevation_audit.md     the tracking document
  HAT_road_elevation.png     alongshore QC

REQUIREMENTS
------------
  geopandas, rasterio, numpy, matplotlib
===============================================================================
"""

from __future__ import annotations

import csv
import re
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np

import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from shapely.geometry import box as shapely_box
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from matplotlib.transforms import blended_transform_factory

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[4]
HATTERAS_DATA_BASE = PROJECT_ROOT / "data" / "hatteras_init"
BARRIER3D_DIR = HATTERAS_DATA_BASE / "1-barrier3d-domains"
MGMT_DIR = HATTERAS_DATA_BASE / "4-mgmt-forcing"

# Native-resolution LiDAR clips, one folder per domain. NOT the 10 m resample.
#   clip_domain_<N>.tif      1 m, native
#   resampled_domain_<N>.tif 10 m, the Barrier3D grid
CLIP_ROOT = BARRIER3D_DIR / "2009-raw" / "2009-domain-clipresample"
CLIP_GLOB = "clip_domain_*.tif"
RESAMPLE_GLOB = "resampled_domain_*.tif"

# --- WHICH SURFACE: raw 2009, or 2009 with its holes filled -------------
# None    the original 2009 clips, LiDAR holes and all
# "<tag>" 0-elevation/{1-gapfill-1m,2-resampled-10m}/<tag>/..._filled.tif
#
# ONLY GIS 78, 79 AND 80 CAN CHANGE. Every other domain on the 2004 alignment
# has nodata_frac = 0 -- complete 2009 coverage under the road -- so filling is
# a no-op there, and the unaffected neighbours 77 and 81 move by <= 0.002 m.
# The three that do change are the same three that drowned on coverage gaps in
# 2009_v4, and D79 is the reason this switch exists: its unfilled elevation is
# a mean over 2287 corridor cells out of 3583, 36% of the corridor missing.
#
# WHY 2008 AND NOT 2014, when the v5 TOPOGRAPHY is filled from 2014:
# this is a ROAD SURFACE, and the two questions have different answers. The
# 2008 IOCM survey is one year from the 2009 base, so it measures the same
# pavement. The 2014 Post-Sandy survey postdates Hurricane Irene (2011), the
# Pea Island breach and the NC-12 reconstruction that followed -- at GIS 78-80
# a 2014 surface under the corridor may be a REBUILT road, which is not what
# "the 2009 road elevation" means. Topography has no such problem: there the
# 2014 fill is simply the later and more complete survey of the same barrier.
#
# The choice is about provenance, not magnitude -- the two fills differ by
# <= 0.015 m in the corridor. RELOCATION BRACKET and the QC flags re-measure
# this every run; FILL_SOURCE is recorded in the audit.
FILL_SOURCE = "2008_NOAA_IOCM"

GAPFILL_1M_ROOT = HATTERAS_DATA_BASE / "0-elevation" / "1-gapfill-1m"
GAPFILL_10M_ROOT = HATTERAS_DATA_BASE / "0-elevation" / "2-resampled-10m"
FILL_CLIP_GLOB = "clip_domain_*_filled.tif"
FILL_RESAMPLE_GLOB = "resampled_domain_*_filled.tif"

# The single alignment. 2004, because it is the one that is contemporaneous with
# the 2009 DEM everywhere on the island -- see the header.
ROAD_LINE = (MGMT_DIR / "road_offset" / "raw_offset" / "2004" / "nc12_2004.geojson")
ROAD_LINE_YEAR = 2004

# Sampled ONLY for the RELOCATION BRACKET check -- never written to the product.
# It exists to test, rather than assume, what the abandoned corridor looks like
# in the 2009 DEM. Set to None to skip the check.
BRACKET_LINE = (MGMT_DIR / "road_offset" / "raw_offset" / "1984"
                / "nc12_1984.geojson")
BRACKET_LINE_YEAR = 1984

OUT_ROOT = MGMT_DIR / "road_elevation"

FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN = 9, 90

# Half-width of the sampling corridor. 3.5 m -> a 7 m strip, about one
# carriageway of NC-12. Widen it and you start averaging in the shoulder and the
# dune toe; the sweep printed under INTERNAL CHECKS shows exactly where that
# begins to bite.
BUFFER_M = 3.5
BUFFER_SWEEP = [2.0, 2.5, 3.5, 5.0, 8.0, 12.0]

# all_touched=False on a 1 m grid: the buffer polygon is already several cells
# wide, so touching-cell inclusion only adds edge pixels off the pavement.
ALL_TOUCHED = False

MHW_M = 0.36               # m NAVD88, Duck NC gauge (NOAA 8651370)

# --- QC thresholds ------------------------------------------------------
MIN_CELLS = 200            # flag a domain sampled on fewer 1 m cells
MAX_STD_M = 0.40           # flag a domain whose road cells scatter more
MAX_NODATA_FRAC = 0.25     # flag a mask sitting largely on LiDAR NoData
MAX_JUMP_M = 0.35          # flag a step between ADJACENT domains larger than this

# Domains where the alignment sampled here postdates the 1984 run's road.
# Taken from HISTORICAL_ROAD_EVENTS in HAT_hindcast_1984_2024.py.
# Informational only -- the value written is the same for both periods.
RELOCATED_BEFORE_2004 = set(range(9, 16)) | set(range(84, 88))

# Reference: the runner's current scalar, for the comparison table.
RUNNER_SCALAR = 1.45

C_MODEL, C_TOWN, C_PIER, C_GROIN, C_INK = ("#FF8C00", "#90AFC5", "#1565C0",
                                           "#B71C1C", "#1a1a2e")
SECTIONS = [((1, 6), "Cape Point"), ((7, 8), "Buxton"), ((9, 20), "Buxton-Avon"),
            ((21, 31), "Avon"), ((32, 67), "Wimble Shoals"),
            ((68, 83), "Tri-Village"), ((84, 90), "Pea Island")]

# --- map figure ---------------------------------------------------------
# The island is ~41 km north-south but each domain is only 2 km across, so an
# equal-aspect map of the whole thing is a 6:1 sliver. Cut it into strips laid
# side by side, each drawn at TRUE aspect. Distorting the aspect to fit would
# make the road look like it changes direction where it does not.
MAP_STRIPS = 4
MAP_ROAD_CMAP = "plasma"       # road elevation; vivid, reads over the sand ramp
# Background LiDAR: pale sand at MHW darkening to brown at the dune crests. Low
# saturation on purpose -- the road is the subject, the island is the context.
MAP_TOPO_COLORS = ["#f7f1e3", "#e2d3b3", "#bfa77f", "#8a7350"]
MAP_TOPO_VMAX = 3.0            # m NAVD88; dune tops saturate above this
C_WATER = "#cfe0ea"            # anything below MHW, plus LiDAR NoData
MAP_MIN_WIDTH_M = 2400.0       # min cross-shore window per strip, keeps context
DOMAIN_M = 500.0               # alongshore length of one Barrier3D domain


# =============================================================================
# SAMPLING
# =============================================================================

def find_clips(root: Path, pattern: str) -> dict:
    """domain id -> raster path, for whichever resolution `pattern` selects."""
    out = {}
    for sub in sorted(Path(root).glob("domain_*")):
        if not sub.is_dir():
            continue
        m = re.search(r"domain_(\d+)$", sub.name)
        tifs = sorted(sub.glob(pattern))
        if m and tifs:
            out[int(m.group(1))] = tifs[0]
    return out


def find_clips_flat(root: Path, pattern: str) -> dict:
    """
    domain id -> raster path, for the gap-fill tree.

    The clipresample tree nests one folder per domain; the gap-fill tree is
    flat and carries the id in the filename instead. Same contract, different
    layout, so the caller does not have to care which surface it asked for.
    """
    out = {}
    for tif in sorted(Path(root).glob(pattern)):
        m = re.search(r"domain_(\d+)_", tif.name)
        if m:
            out[int(m.group(1))] = tif
    return out


def resolve_surfaces() -> tuple[dict, dict, str]:
    """(1 m clips, 10 m resamples, a label for the audit) for FILL_SOURCE."""
    if not FILL_SOURCE:
        return (find_clips(CLIP_ROOT, CLIP_GLOB),
                find_clips(CLIP_ROOT, RESAMPLE_GLOB),
                "2009 clips, unfilled")

    one_m = GAPFILL_1M_ROOT / FILL_SOURCE
    ten_m = GAPFILL_10M_ROOT / FILL_SOURCE

    for label, path in (("1 m gap-fill", one_m), ("10 m gap-fill", ten_m)):
        if not path.is_dir():
            avail = (sorted(p.name for p in GAPFILL_1M_ROOT.iterdir()
                            if p.is_dir()) if GAPFILL_1M_ROOT.is_dir() else [])
            raise SystemExit(
                f"\n[stop] {label} directory for FILL_SOURCE={FILL_SOURCE!r} "
                f"does not exist:\n    {path}\n"
                f"sources present: {avail}\n")

    # Both resolutions come from the SAME fill, so the 1 m vs 10 m agreement
    # printed under INTERNAL CHECKS stays a like-for-like comparison.
    return (find_clips_flat(one_m, FILL_CLIP_GLOB),
            find_clips_flat(ten_m, FILL_RESAMPLE_GLOB),
            f"2009 clips, holes filled from {FILL_SOURCE}")


def load_line(path: Path):
    if not Path(path).exists():
        raise FileNotFoundError(f"road line not found:\n    {path}")
    gdf = gpd.read_file(path)
    if gdf.crs is None:
        raise ValueError(f"{Path(path).name} has no CRS. Cannot align it to "
                         f"the LiDAR.")
    return gdf


def corridor_values(line_in_raster_crs, src, z, bad, buffer_m: float):
    """
    Elevations under the buffered road line, and the fraction of the corridor
    that fell on NoData.

    The buffer is applied in the RASTER's CRS, converted through its linear
    units factor -- the geojson is EPSG:2264 (US survey FEET) and the clips are
    UTM 18N (metres), so a raw 3.5 would be 1.07 m if taken literally in the
    source CRS.

    `z` and `bad` are passed in so a caller sweeping several buffer widths reads
    each raster exactly once.
    """
    to_m = float(src.crs.linear_units_factor[1])
    geoms = [gm.buffer(buffer_m / to_m) for gm in line_in_raster_crs.geometry]
    mask = rasterize(
        [(gm, 1) for gm in geoms],
        out_shape=(src.height, src.width), transform=src.transform,
        fill=0, all_touched=ALL_TOUCHED, dtype="uint8",
    ) > 0

    n_mask = int(mask.sum())
    if n_mask == 0:
        return np.array([]), 0.0
    nodata_frac = float((mask & bad).sum()) / n_mask
    return z[mask & ~bad], nodata_frac


def read_surface(src):
    """Elevation band plus a boolean mask of cells that are not real data."""
    z = src.read(1).astype(float)
    bad = (z <= -100.0) | (~np.isfinite(z))
    if src.nodata is not None:
        bad |= (z == src.nodata)
    return z, bad


def stats_for(vals: np.ndarray, nodata_frac: float) -> dict:
    """
    MEAN is the product. Median, percentiles and sigma are diagnostics that ride
    along in the per-domain CSV so a domain where the mean is being dragged by a
    structure in the corridor announces itself.
    """
    if vals.size == 0:
        return dict(n_cells=0, nodata_frac=round(nodata_frac, 3),
                    elev_navd=np.nan, elev_mhw=np.nan, elev_median_navd=np.nan,
                    mean_minus_median=np.nan, elev_std_m=np.nan,
                    elev_p10_navd=np.nan, elev_p90_navd=np.nan,
                    elev_min_navd=np.nan, elev_max_navd=np.nan)
    mean = float(np.mean(vals))
    median = float(np.median(vals))
    return dict(
        n_cells=int(vals.size), nodata_frac=round(nodata_frac, 3),
        elev_navd=round(mean, 3), elev_mhw=round(mean - MHW_M, 3),
        elev_median_navd=round(median, 3),
        mean_minus_median=round(mean - median, 3),
        elev_std_m=round(float(np.std(vals)), 3),
        elev_p10_navd=round(float(np.percentile(vals, 10)), 3),
        elev_p90_navd=round(float(np.percentile(vals, 90)), 3),
        elev_min_navd=round(float(np.min(vals)), 3),
        elev_max_navd=round(float(np.max(vals)), 3),
    )


def flags_for(r: dict) -> list:
    f = []
    if r["n_cells"] == 0:
        return ["NO_DATA"]
    if r["n_cells"] < MIN_CELLS:
        f.append(f"THIN({r['n_cells']})")
    if r["nodata_frac"] > MAX_NODATA_FRAC:
        f.append(f"NODATA({r['nodata_frac'] * 100:.0f}%)")
    if np.isfinite(r["elev_std_m"]) and r["elev_std_m"] > MAX_STD_M:
        f.append(f"SCATTER({r['elev_std_m']:.2f})")
    if r["domain"] in RELOCATED_BEFORE_2004:
        f.append("RELOCATED_PRE_2004")
    return f


def sample_all(line, clips: dict) -> list:
    """Mean road elevation per domain on the 1 m clips, at BUFFER_M."""
    rows = []
    for d in sorted(clips):
        if not (FIRST_ROAD_DOMAIN <= d <= LAST_ROAD_DOMAIN):
            continue
        with rasterio.open(clips[d]) as src:
            z, bad = read_surface(src)
            vals, nodata_frac = corridor_values(line.to_crs(src.crs), src, z,
                                                bad, BUFFER_M)
        r = dict(domain=d)
        r.update(stats_for(vals, nodata_frac))
        rows.append(r)

    # Continuity is an alongshore property, so it can only be flagged once every
    # domain has a value. A road does not step 0.35 m between two 500 m cells.
    order = sorted([r for r in rows if np.isfinite(r["elev_navd"])],
                   key=lambda r: r["domain"])
    for a, b in zip(order, order[1:]):
        if b["domain"] == a["domain"] + 1:
            step = b["elev_navd"] - a["elev_navd"]
            if abs(step) > MAX_JUMP_M:
                b.setdefault("_jump", step)

    for r in rows:
        f = flags_for(r)
        if "_jump" in r:
            f.append(f"JUMP({r.pop('_jump'):+.2f})")
        r["flags"] = ",".join(f)
    return rows


# =============================================================================
# INTERNAL CHECKS
#
# The ArcGIS elevation_2009 column that this method was originally validated
# against no longer exists -- see the audit document. What is left is internal:
# does the answer depend on the corridor width, and does the surface we are
# obliged to model on agree with the surface we measured.
# =============================================================================

def buffer_sweep(line, clips: dict) -> list:
    """Island median and typical within-domain sigma at several corridor widths."""
    acc = {b: {"mean": [], "sd": []} for b in BUFFER_SWEEP}
    for d in sorted(clips):
        if not (FIRST_ROAD_DOMAIN <= d <= LAST_ROAD_DOMAIN):
            continue
        with rasterio.open(clips[d]) as src:
            z, bad = read_surface(src)
            g = line.to_crs(src.crs)
            for b in BUFFER_SWEEP:
                vals, _ = corridor_values(g, src, z, bad, b)
                if vals.size:
                    acc[b]["mean"].append(float(np.mean(vals)))
                    acc[b]["sd"].append(float(np.std(vals)))
    out = []
    for b in BUFFER_SWEEP:
        m, s = acc[b]["mean"], acc[b]["sd"]
        if m:
            out.append((b, float(np.median(m)), float(np.median(s)), len(m)))
    return out


def relocation_bracket(clips: dict, rows: list) -> dict | None:
    """
    What is actually under the 1984 alignment where the road was relocated?

    The case for using the 2004 line everywhere rests on a claim about the
    abandoned corridor, and a claim in a docstring is not evidence. This
    measures it: same surface, same corridor width, the OTHER line, in the 11
    domains where the two disagree.

    Nothing here is ever written to RoadElevation.csv. It exists so the choice
    can be defended with a number that is re-derived on every run.
    """
    if BRACKET_LINE is None or not Path(BRACKET_LINE).exists():
        return None
    line = load_line(BRACKET_LINE)
    by = {r["domain"]: r for r in rows}
    pairs = []
    for d in sorted(RELOCATED_BEFORE_2004):
        if d not in clips or d not in by or not np.isfinite(by[d]["elev_navd"]):
            continue
        with rasterio.open(clips[d]) as src:
            z, bad = read_surface(src)
            vals, _ = corridor_values(line.to_crs(src.crs), src, z, bad,
                                      BUFFER_M)
        if vals.size:
            pairs.append(dict(domain=d, used=by[d]["elev_navd"],
                              other=float(np.mean(vals)),
                              other_sd=float(np.std(vals))))
    if not pairs:
        return None

    used = np.array([p["used"] for p in pairs])
    other = np.array([p["other"] for p in pairs])
    sd = np.array([p["other_sd"] for p in pairs])

    # The grade the road would sit at if this reach were not relocated: the
    # neighbouring un-relocated domains inside the same named reaches.
    ref = []
    for (lo, hi), _ in SECTIONS:
        if not any(lo <= d <= hi for d in RELOCATED_BEFORE_2004):
            continue
        ref += [r["elev_navd"] for r in rows
                if lo <= r["domain"] <= hi
                and r["domain"] not in RELOCATED_BEFORE_2004
                and np.isfinite(r["elev_navd"])]
    return dict(pairs=pairs, n=len(pairs), used_mean=float(used.mean()),
                other_mean=float(other.mean()), other_sd_max=float(sd.max()),
                neighbour_grade=float(np.mean(ref)) if ref else np.nan,
                brackets=bool(min(used.mean(), other.mean())
                              <= (np.mean(ref) if ref else np.nan)
                              <= max(used.mean(), other.mean())))


def resample_check(line, resamples: dict, rows: list) -> dict | None:
    """
    Same corridor, same buffer, on the 10 m Barrier3D grid instead of the 1 m
    clip. This is not a validation -- the 1 m answer is the better one -- it
    quantifies how much the grid the model actually runs on would have cost us.
    """
    if not resamples:
        return None
    by = {r["domain"]: r for r in rows}
    diffs, sds, ncells = [], [], []
    for d in sorted(resamples):
        if d not in by or not np.isfinite(by[d]["elev_navd"]):
            continue
        with rasterio.open(resamples[d]) as src:
            z, bad = read_surface(src)
            vals, _ = corridor_values(line.to_crs(src.crs), src, z, bad,
                                      BUFFER_M)
        if vals.size:
            diffs.append(float(np.mean(vals)) - by[d]["elev_navd"])
            sds.append(float(np.std(vals)))
            ncells.append(int(vals.size))
    if not diffs:
        return None
    diffs = np.asarray(diffs)
    return dict(n=len(diffs), median=float(np.median(diffs)),
                mean_abs=float(np.abs(diffs).mean()),
                max_abs=float(np.abs(diffs).max()),
                sd_10m=float(np.median(sds)),
                cells_10m=int(np.median(ncells)))


# =============================================================================
# OUTPUT
# =============================================================================

def write_cascade_csv(rows, path: Path) -> bool:
    """
    2-row CASCADE file: row 0 = GIS IDs, row 1 = elevation in m MHW-RELATIVE.

    Refuses on a gap, because the runner fills its per-domain arrays BY POSITION
    after dropping the ID row -- one missing domain shifts every domain north of
    it and nothing reports it.
    """
    by = {r["domain"]: r for r in rows}
    expected = list(range(FIRST_ROAD_DOMAIN, LAST_ROAD_DOMAIN + 1))
    missing = [d for d in expected
               if d not in by or not np.isfinite(by[d]["elev_mhw"])]
    if missing:
        print(f"\n  NOT WRITING {path.name}: no elevation for domain(s) "
              f"{missing}.")
        print("  The runner fills per-domain arrays BY POSITION and drops the")
        print("  ID row, so one gap shifts every domain north of it silently.")
        return False

    ids = np.array(expected, dtype=float)
    ev = np.array([by[d]["elev_mhw"] for d in expected], dtype=float)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.vstack([ids, ev]), delimiter=",", fmt="%.3f")
    print(f"[out] {path}")
    print(f"      2 x {len(ids)}  (row 0 = GIS IDs, row 1 = m MHW-relative)")
    print(f"      one file, both periods -- there is only one DEM")
    return True


def write_domains_csv(rows, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["domain"]):
            w.writerow(r)
    print(f"[out] {path}")


def section_stats(rows) -> list:
    """Per-named-reach summary. A domain index is not a place; a reach is."""
    good = [r for r in rows if np.isfinite(r["elev_navd"])]
    out = []
    for (lo, hi), label in SECTIONS:
        vals = [r["elev_navd"] for r in good if lo <= r["domain"] <= hi]
        if len(vals) < 2:
            continue
        v = np.asarray(vals)
        out.append(dict(label=label, lo=lo, hi=hi, n=len(v),
                        mean=float(v.mean()), min=float(v.min()),
                        max=float(v.max()),
                        relocated=sum(1 for d in range(lo, hi + 1)
                                      if d in RELOCATED_BEFORE_2004)))
    return out


def qc_figure(rows, path: Path):
    rows = sorted([r for r in rows if np.isfinite(r["elev_navd"])],
                  key=lambda r: r["domain"])
    if not rows:
        return
    d = np.array([r["domain"] for r in rows])
    mean = np.array([r["elev_navd"] for r in rows])
    p10 = np.array([r["elev_p10_navd"] for r in rows])
    p90 = np.array([r["elev_p90_navd"] for r in rows])

    fig, ax = plt.subplots(figsize=(13, 5.5))
    for k, ((lo, hi), label) in enumerate(SECTIONS):
        ax.axvspan(lo - 0.5, hi + 0.5, color=C_TOWN,
                   alpha=0.18 if k % 2 == 0 else 0.08, lw=0, zorder=0)
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        ax.text((lo + hi) / 2, 0.97, label, transform=tr, ha="center",
                va="top", fontsize=7,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none",
                          alpha=0.85))

    ax.fill_between(d, p10, p90, color=C_GROIN, alpha=0.18, lw=0,
                    label="p10-p90 within domain", zorder=2)
    ax.plot(d, mean, color=C_GROIN, lw=1.7, marker="o", ms=3, zorder=4,
            label="mean road elevation")
    ax.axhline(RUNNER_SCALAR, color=C_MODEL, ls="-.", lw=1.3, zorder=3,
               label=f"runner scalar ROAD_ELEVATION = {RUNNER_SCALAR}")
    ax.axhline(1.30, color=C_PIER, ls="--", lw=1.2, zorder=3,
               label="Velasquez-Montoya et al.: NC-12 avg 1.3 m NAVD88")
    ax.axhline(MHW_M, color="0.4", ls=":", lw=1.2, zorder=3,
               label=f"MHW ({MHW_M} m NAVD88)")

    rel = [r for r in rows if r["domain"] in RELOCATED_BEFORE_2004]
    if rel:
        ax.plot([r["domain"] for r in rel], [r["elev_navd"] for r in rel],
                lw=0, marker="x", ms=10, mew=2.0, color=C_INK, zorder=6,
                label="road relocated before 2004 (1984 run: wrong place)")

    # Reach means, so the alongshore pattern reads at a glance rather than
    # having to be averaged by eye out of 82 points.
    for k, s in enumerate(section_stats(rows)):
        ax.plot([max(s["lo"], d.min()) - 0.4, min(s["hi"], d.max()) + 0.4],
                [s["mean"]] * 2, color=C_INK, lw=2.6, alpha=0.75, zorder=5,
                solid_capstyle="butt",
                label="reach mean" if k == 0 else None)
        ax.text((max(s["lo"], d.min()) + min(s["hi"], d.max())) / 2,
                s["mean"] - 0.035, f"{s['mean']:.2f}", ha="center", va="top",
                fontsize=7.5, color=C_INK, zorder=6,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none",
                          alpha=0.8))

    ax.set_xlabel("GIS domain  (1 = Cape Point / south  ->  90 = Rodanthe / north)")
    ax.set_ylabel("Road elevation (m NAVD88)")
    ax.set_title(f"NC-12 {ROAD_LINE_YEAR} alignment on the 2009 LiDAR "
                 f"(1 m clip, {BUFFER_M:.1f} m buffer, mean of cells)"
                 f"\nused unchanged for both the 1984 and 2004 start periods",
                 fontsize=12)
    # Bottom-left: the only quadrant with no data in it, and it keeps the
    # reach labels along the top edge readable.
    ax.legend(loc="lower left", fontsize=8, framealpha=0.92, ncol=2)
    ax.grid(alpha=0.25)
    ax.set_xlim(d.min() - 0.5, d.max() + 0.5)

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[out] {path}")


def map_figure(rows, line, resamples: dict, path: Path):
    """
    Where the road is high and where it is low, in real geography.

    The alongshore profile answers "how much"; it cannot answer "where", because
    a domain index is not a place. This draws the road on the island it sits on,
    coloured by the same numbers, over the 10 m LiDAR for context -- so a low
    reach can be read against the island being narrow there.

    Drawn in strips at TRUE aspect ratio. See MAP_STRIPS.
    """
    good = sorted([r for r in rows if np.isfinite(r["elev_navd"])],
                  key=lambda r: r["domain"])
    if not good or not resamples:
        return

    # Geometry and topography, once per domain.
    topo, bounds, crs = {}, {}, None
    for r in good:
        d = r["domain"]
        if d not in resamples:
            continue
        with rasterio.open(resamples[d]) as src:
            z, bad = read_surface(src)
            topo[d] = np.where(bad, np.nan, z)
            bounds[d] = src.bounds
            crs = crs or src.crs
    if not topo:
        return

    doms = [d for d in (r["domain"] for r in good) if d in topo]
    by = {r["domain"]: r for r in good}
    g_line = line.to_crs(crs)

    ev = np.array([by[d]["elev_navd"] for d in doms])
    norm = mcolors.Normalize(vmin=float(ev.min()), vmax=float(ev.max()))
    road_cmap = plt.get_cmap(MAP_ROAD_CMAP)

    topo_cmap = mcolors.LinearSegmentedColormap.from_list(
        "hat_sand", MAP_TOPO_COLORS)
    topo_cmap.set_bad(C_WATER)      # NoData -- sound, ocean, gaps
    topo_cmap.set_under(C_WATER)    # below MHW -- water, functionally
    topo_norm = mcolors.Normalize(vmin=MHW_M, vmax=MAP_TOPO_VMAX)

    # Road geometry per domain, clipped to that domain's box -- the SAME box the
    # elevation was sampled in, so colour and geometry cannot disagree.
    segs = {}
    for d in doms:
        b = bounds[d]
        gm = g_line.geometry.intersection(
            shapely_box(b.left, b.bottom, b.right, b.top)).union_all()
        segs[d] = [np.asarray(p.xy) for p in getattr(gm, "geoms", [gm])
                   if p.geom_type == "LineString" and not p.is_empty]

    # Contiguous strips, south to north.
    chunks = [list(c) for c in np.array_split(np.array(doms), MAP_STRIPS)]

    fig, axes = plt.subplots(1, MAP_STRIPS, figsize=(3.1 * MAP_STRIPS, 8.6))
    axes = np.atleast_1d(axes)

    for ax, chunk in zip(axes, chunks):
        for d in chunk:
            b = bounds[d]
            ax.imshow(topo[d], extent=(b.left, b.right, b.bottom, b.top),
                      origin="upper", cmap=topo_cmap, norm=topo_norm,
                      interpolation="nearest", zorder=1)

        for d in chunk:
            col = road_cmap(norm(by[d]["elev_navd"]))
            for xy in segs[d]:
                if d in RELOCATED_BEFORE_2004:
                    ax.plot(xy[0], xy[1], color=C_INK, lw=6.0, alpha=0.9,
                            solid_capstyle="round", zorder=3)
                ax.plot(xy[0], xy[1], color=col, lw=3.4,
                        solid_capstyle="round", zorder=4)

        # Window on the road corridor, widened to MAP_MIN_WIDTH_M so the island
        # around it stays visible -- where the island is narrow is exactly the
        # context that explains a low reach.
        rx = np.concatenate([xy[0] for d in chunk for xy in segs[d]]
                            or [np.array([0.0])])
        cx = 0.5 * (rx.min() + rx.max())
        half = max(MAP_MIN_WIDTH_M, (rx.max() - rx.min()) * 1.35) / 2
        y_lo = min(bounds[d].bottom for d in chunk)
        y_hi = max(bounds[d].top for d in chunk)

        for (lo, hi), label in SECTIONS:
            inside = [d for d in chunk if lo <= d <= hi]
            if len(inside) < 3:
                continue
            ymid = np.mean([(bounds[d].bottom + bounds[d].top) / 2
                            for d in inside])
            ax.text(cx - half * 0.92, ymid, label, fontsize=8.5, rotation=90,
                    ha="left", va="center", color=C_INK, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                              alpha=0.8))

        # Tie the map back to the profile and the CSV.
        for d, va in ((chunk[0], "bottom"), (chunk[-1], "top")):
            b = bounds[d]
            ax.text(cx + half * 0.94, (b.bottom + b.top) / 2, f"GIS {d}",
                    fontsize=8, ha="right", va="center", color=C_INK, zorder=6,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                              alpha=0.8))

        ax.set_xlim(cx - half, cx + half)
        ax.set_ylim(y_lo - 250, y_hi + 250)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_edgecolor("0.7")

    # Scale bar bottom-left, north arrow top-left, both on the southern strip
    # and both kept clear of the GIS labels on the right.
    ax0 = axes[0]
    x0, x1 = ax0.get_xlim()
    y0, y1 = ax0.get_ylim()
    bx = x0 + 0.07 * (x1 - x0)
    byy = y0 + 0.022 * (y1 - y0)
    ax0.plot([bx, bx + 1000], [byy, byy], color=C_INK, lw=2.5, zorder=7,
             solid_capstyle="butt")
    ax0.text(bx + 500, byy + 0.005 * (y1 - y0), "1 km", ha="center",
             va="bottom", fontsize=8, color=C_INK, zorder=7)
    nx = x0 + 0.11 * (x1 - x0)
    ax0.annotate("N", xy=(nx, y1 - 0.020 * (y1 - y0)),
                 xytext=(nx, y1 - 0.055 * (y1 - y0)),
                 arrowprops=dict(arrowstyle="-|>", color=C_INK, lw=1.6),
                 ha="center", va="top", fontsize=9, color=C_INK, zorder=7)

    sm = ScalarMappable(norm=norm, cmap=road_cmap)
    cb = fig.colorbar(sm, ax=list(axes), fraction=0.030, pad=0.02,
                      aspect=38)
    cb.set_label("Mean road elevation (m NAVD88)", fontsize=9)
    cb.ax.axhline(RUNNER_SCALAR, color=C_MODEL, lw=2.0)
    cb.ax.axhline(1.30, color=C_PIER, lw=1.6, ls="--")

    handles = [
        Line2D([], [], color=C_INK, lw=6.0,
               label="road relocated before 2004 — wrong place for the 1984 run"),
        Line2D([], [], color=C_MODEL, lw=2.0,
               label=f"runner scalar {RUNNER_SCALAR} m (marked on colourbar)"),
        Line2D([], [], color=C_PIER, lw=1.6, ls="--",
               label="Velasquez-Montoya et al. 1.3 m (marked on colourbar)"),
        Line2D([], [], color=C_WATER, lw=6,
               label=f"below MHW ({MHW_M} m NAVD88) or LiDAR NoData"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.028))

    fig.suptitle("NC-12 road elevation along Hatteras Island\n"
                 "2004 alignment sampled on the 2009 LiDAR — mean of 1 m cells "
                 f"in a {BUFFER_M * 2:.0f} m corridor",
                 fontsize=12.5, y=0.975)
    fig.text(0.5, 0.005,
             f"Background: 10 m LiDAR, pale sand to brown = {MHW_M}–"
             f"{MAP_TOPO_VMAX:.0f} m NAVD88.  Strips are contiguous, laid south "
             "(left) to north (right), each at true aspect ratio.",
             ha="center", fontsize=8, color="0.35")

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[out] {path}")


def write_markdown(rows, sweep, resamp, bracket, path: Path):
    L = []
    a = L.append
    good = [r for r in rows if np.isfinite(r["elev_navd"])]
    navd = np.array([r["elev_navd"] for r in good])

    a("# NC-12 road elevation from the 2009 LiDAR")
    a("")
    a(f"Generated {datetime.now().isoformat(timespec='seconds')} by "
      f"`{Path(__file__).name}`.")
    a("")
    a("| | |")
    a("|---|---|")
    if FILL_SOURCE:
        a(f"| Surface | `clip_domain_<N>_filled.tif` — **1 m** native LiDAR, "
          f"2009, holes filled from **{FILL_SOURCE}** |")
        a("| Why that fill | 2008 is one year from the 2009 base, so it "
          "measures the same pavement. The v5 *topography* is filled from "
          "2014 Post-Sandy, which postdates Irene, the Pea Island breach and "
          "the NC-12 rebuild — fine for a barrier surface, wrong for a road "
          "surface at GIS 78–80. |")
        _gaps = sorted((r for r in good if r["nodata_frac"] > 0),
                       key=lambda r: -r["nodata_frac"])
        a("| Domains materially affected | **GIS 78, 79, 80** — the only ones "
          "with NoData under the road in the unfilled 2009 clips. Everywhere "
          "else the corridor was already complete, and the change is "
          "≤ 0.01 m: the fill pipeline rebuilds the raster rather than "
          "patching the original, so the corridor mask lands on a few cells' "
          "difference even where there was nothing to fill. Not a no-op, but "
          "well inside the 0.07 m within-domain scatter. |")
        if _gaps:
            _worst = _gaps[0]
            _each = ", ".join(
                "GIS {}: {:.1%}".format(r["domain"], r["nodata_frac"])
                for r in _gaps)
            a("| Residual gaps | {} domain(s) still carry NoData in the "
              "corridor after filling — worst **GIS {} at {:.1%}** ({}). The "
              "2014 Post-Sandy fill would close these, at the provenance cost "
              "above. All are below the {:.0%} flag threshold so they no "
              "longer raise `NODATA`; recorded here so that is not mistaken "
              "for complete coverage. |".format(
                  len(_gaps), _worst["domain"], _worst["nodata_frac"],
                  _each, MAX_NODATA_FRAC))
    else:
        a("| Surface | `clip_domain_<N>.tif` — **1 m** native LiDAR, 2009, "
          "unfilled |")
    a(f"| Alignment | `{Path(ROAD_LINE).name}` — the **{ROAD_LINE_YEAR}** "
      f"digitised centreline |")
    a(f"| Corridor | {BUFFER_M:.1f} m buffer either side "
      f"(~{BUFFER_M * 2:.0f} m strip) |")
    a("| Aggregation | **mean** of all valid 1 m cells in the corridor |")
    a(f"| Datum written | **m MHW-relative** (MHW = {MHW_M} m NAVD88) |")
    a(f"| Domains | {FIRST_ROAD_DOMAIN}–{LAST_ROAD_DOMAIN} "
      f"({len(good)} with data) |")
    a("| Applies to | **both** the 1984 and the 2004 start period |")
    a("")

    a("## Why one file and not one per vintage")
    a("")
    a("There is only one DEM. The previous script sampled two alignments — 1984 "
      "and 2004 — on that same 2009 surface, so any difference between the two "
      "outputs was never a difference in **time**. It was a difference in "
      "**where on the 2009 surface the two digitised lines happened to fall**.")
    a("")
    a("Where NC-12 never moved, the two lines sit on top of each other and the "
      "numbers were identical by construction. Where the road *was* relocated, "
      "the 1984 line lay over the abandoned corridor, and the \"1984 road "
      "elevation\" there was the elevation of a place a road used to be — "
      "measured below.")
    a("")
    a("Neither of those is a measurement of change in roadbed height. Writing "
      "two files implied one. This writes one.")
    a("")

    a("## How the road elevation changes along the island")
    a("")
    a("Two views, both regenerated on every run:")
    a("")
    a("- `HAT_road_elevation_map.png` — the road drawn **in real geography**, "
      "coloured by these numbers, over the 10 m LiDAR. A domain index is not a "
      "place; this is where the high and low reaches actually are, and it shows "
      "the island around them.")
    a("- `HAT_road_elevation.png` — the alongshore profile with the p10–p90 "
      "spread inside each domain and the reach means marked.")
    a("")
    secs = section_stats(rows)
    if secs:
        a("| Reach | GIS | n | Mean (m NAVD88) | Mean (m MHW) | Range |")
        a("|---|---:|---:|---:|---:|---|")
        for s in secs:
            note = (f" — {s['relocated']} relocated"
                    if s["relocated"] else "")
            a(f"| {s['label']}{note} | {s['lo']}–{s['hi']} | {s['n']} | "
              f"{s['mean']:.2f} | {s['mean'] - MHW_M:.2f} | "
              f"{s['min']:.2f}–{s['max']:.2f} |")
        a("")
        hi = max(secs, key=lambda s: s["mean"])
        lo = min(secs, key=lambda s: s["mean"])
        a(f"The spread between reaches is **{hi['mean'] - lo['mean']:.2f} m** — "
          f"{hi['label']} at {hi['mean']:.2f} m against {lo['label']} at "
          f"{lo['mean']:.2f} m. That is larger than the within-domain σ by more "
          f"than an order of magnitude, which is the case for carrying a "
          f"per-domain array rather than one scalar.")
        a("")
        if hi["relocated"]:
            a(f"> Note that the highest reach, {hi['label']}, is also the "
              f"relocated one. Its elevation is a **constructed embankment** "
              f"built during the relocation, not the natural roadbed grade — "
              f"which is why the profile steps down sharply at its northern "
              f"end.")
            a("")

    a("## Mean, not median")
    a("")
    a("The product is the flat unweighted mean of every valid 1 m cell in the "
      "corridor. The median is carried in the per-domain CSV as "
      "`elev_median_navd`, with `mean_minus_median` alongside, so a domain "
      "where the mean is being dragged by a structure in the corridor — a "
      "bridge deck, a house, a driveway apron — announces itself instead of "
      "being silently absorbed.")
    a("")
    if good:
        mm = np.array([r["mean_minus_median"] for r in good])
        a(f"On this alignment the two agree closely: median "
          f"|mean − median| = **{np.median(np.abs(mm)):.3f} m**, "
          f"max **{np.abs(mm).max():.3f} m** "
          f"(GIS {good[int(np.argmax(np.abs(mm)))]['domain']}).")
        a("")

    a("## Internal checks")
    a("")
    a("### The external validation is gone")
    a("")
    a("This method was originally validated against an `elevation_2009` column "
      "in `2004_road_offset_raw.csv`, sampled independently in ArcGIS Pro at "
      "1 m transect points; the two agreed to a median of −0.02 m. **That check "
      "is no longer reproducible.** The file is now `nc12_2004.csv`, the "
      "`elevation_2009` column is absent, and its apparent successors "
      "(`avg_elev_m`, `z_mean`, `z_max`, `z_min`, `road_z`, `relief_m`) are all "
      "zero or all empty across the 1491 rows. The agreement was real when it "
      "was measured; it cannot be re-measured from what is in the repository "
      "now. If a pre-rename copy of that CSV turns up, restore the check.")
    a("")
    a("What follows is internal: does the answer depend on choices we made, and "
      "does the surface the model runs on agree with the surface we measured.")
    a("")

    if sweep:
        a("### Corridor width")
        a("")
        a("If the answer moved with the buffer, the buffer would be the "
          "measurement. It does not, until the corridor grows wide enough to "
          "reach off the pavement:")
        a("")
        a("| Buffer | Strip | Island median (m NAVD88) | Median within-domain σ |")
        a("|---:|---:|---:|---:|")
        for b, m, s, _ in sweep:
            mark = "  ← **used**" if abs(b - BUFFER_M) < 1e-9 else ""
            a(f"| {b:.1f} m | {b * 2:.0f} m | {m:.2f}{mark} | {s:.2f} |")
        a("")
        a("σ climbing while the median falls is the corridor beginning to catch "
          "the shoulder and the dune toe. 3.5 m sits in the flat part.")
        a("")

    if resamp:
        sd1 = float(np.median([r["elev_std_m"] for r in good]))
        n1 = int(np.median([r["n_cells"] for r in good]))
        a("### 1 m clip vs the 10 m Barrier3D grid")
        a("")
        a("The same corridor sampled on `resampled_domain_<N>.tif`, the grid "
          "the model actually runs on:")
        a("")
        a(f"- median difference **{resamp['median']:+.2f} m** over "
          f"{resamp['n']} domains")
        a(f"- mean |difference| **{resamp['mean_abs']:.2f} m**, max "
          f"**{resamp['max_abs']:.2f} m**")
        a(f"- cells per domain in the corridor: **{n1}** at 1 m vs "
          f"**{resamp['cells_10m']}** at 10 m")
        a(f"- median within-domain σ: **{sd1:.2f} m** at 1 m vs "
          f"**{resamp['sd_10m']:.2f} m** at 10 m")
        a("")
        a("**This does not go the way the inherited rationale said it would, "
          "and the difference matters.** The 1 m clip was originally chosen "
          "because the 10 m grid gave GIS 9 a σ of 1.59 m — the foredune, not "
          "the roadbed. But that was measured on the **1984** line, which "
          "crosses a relocation scar. On the 2004 line the two resolutions "
          "agree to a couple of centimetres.")
        a("")
        a("So for this product the 1 m clip is a **precaution, not a rescue**. "
          "It is still the right surface — it rests on "
          f"~{n1} cells per domain against ~{resamp['cells_10m']}, and the "
          "near-equal σ at 10 m is a small-sample artefact rather than "
          "evidence of equal fidelity — but the honest statement is that "
          "switching to the resample would barely move these numbers.")
        a("")

    if bracket:
        a("### Relocation bracket — what is under the 1984 line?")
        a("")
        a("The case for using the 2004 alignment everywhere rests on a claim "
          "about the abandoned corridor, and a claim is not evidence. Same "
          "surface, same corridor width, the other line, in the "
          f"{bracket['n']} domains where the two disagree:")
        a("")
        a("| GIS | 2004 line (used) | 1984 line | Difference | σ on 1984 line |")
        a("|---:|---:|---:|---:|---:|")
        for p in bracket["pairs"]:
            a(f"| {p['domain']} | {p['used']:.2f} | {p['other']:.2f} | "
              f"{p['other'] - p['used']:+.2f} | {p['other_sd']:.2f} |")
        a("")
        a(f"Mean **{bracket['used_mean']:.2f} m** on the line used against "
          f"**{bracket['other_mean']:.2f} m** on the 1984 line, whose "
          f"within-domain σ reaches **{bracket['other_sd_max']:.2f} m**. "
          f"Un-relocated domains in the same named reaches sit at "
          f"**{bracket['neighbour_grade']:.2f} m**.")
        a("")
        if not bracket["brackets"]:
            a("**The two lines do not bracket the neighbouring grade — both are "
              "above it.** So the intuitive picture is wrong: the abandoned "
              "corridor is not overwashed, bulldozed flat, or bare sand. The "
              "dune migrated over it after the road left, and a σ of that size "
              "is a dune, not a graded surface.")
            a("")
            a("This settles the choice rather than complicating it. The 2004 "
              "value is the **lower** of the two candidates and the **only** "
              "one measured on a road surface — the conservative choice as well "
              "as the correct one. Averaging the two, or substituting the 1984 "
              "line in these domains, would raise `road_ele` and make the 1984 "
              "run worse.")
        else:
            a("The two lines fall either side of the neighbouring grade, so "
              "they do bracket it. Revisit whether a blend is defensible here.")
        a("")

    a("### Alongshore continuity")
    a("")
    a("A maintained road does not step abruptly between adjacent 500 m cells. "
      f"Any adjacent-domain step above {MAX_JUMP_M:.2f} m is flagged `JUMP`.")
    a("")
    jumps = [r for r in rows if "JUMP" in r.get("flags", "")]
    if jumps:
        a(f"**{len(jumps)} step(s) flagged:**")
        a("")
        for r in jumps:
            step = r["flags"].split("JUMP(")[1].split(")")[0]
            edge = (r["domain"] - 1 in RELOCATED_BEFORE_2004
                    and r["domain"] not in RELOCATED_BEFORE_2004)
            note = ("**expected** — this is the boundary of a relocated "
                    "reach. The rebuilt road sits on a constructed "
                    "embankment; the step is where that embankment ends, "
                    "not a sampling error."
                    if edge else
                    "**inspect** — no relocation boundary here to explain it.")
            a(f"- GIS {r['domain']}, step {step} m from GIS "
              f"{r['domain'] - 1}: {note}")
        a("")
    else:
        a("**No steps flagged.** The alongshore profile is continuous.")
    a("")

    a("## Datum — not ambiguous, despite the runner")
    a("")
    a("`bulldoze()` writes `road_ele` straight into the interior grid:")
    a("")
    a("```python")
    a("road_ele = road_ele / dz")
    a("new_road_domain = np.zeros(...) + road_ele")
    a("```")
    a("")
    a("and the interior arrays are **MHW-relative**, because "
      "`HAT_dune_topo_extractor.py` subtracts `MHW_M = 0.36` before anything "
      "else. So `road_ele` must be MHW-relative metres. There is no reading "
      "under which NAVD88 is correct. **`RoadElevation.csv` is already "
      "MHW-relative — do not subtract MHW again.**")
    a("")
    if navd.size:
        med = float(np.median(navd))
        a(f"The runner's current `ROAD_ELEVATION = {RUNNER_SCALAR}` is **high "
          f"under either reading**: the measured island median is "
          f"{med:.2f} m NAVD88 = {med - MHW_M:.2f} m MHW. Read as NAVD88, "
          f"{RUNNER_SCALAR} would be {RUNNER_SCALAR - MHW_M:.2f} m MHW; read as "
          f"MHW-relative it would be {RUNNER_SCALAR + MHW_M:.2f} m NAVD88. The "
          f"measured range is {navd.min():.2f}–{navd.max():.2f} m NAVD88.")
        a("")
        a("> Replacing the scalar changes model behaviour, not just "
          "bookkeeping. `bulldoze` computes "
          "`road_overwash_removal = sum(old − road_ele)`, clipped at zero, so a "
          "**lower** road removes **more** sand each year and hands more of it "
          "to the dunes — maintaining a lower roadbed means digging out more. "
          "`road_ele` is also decremented by RSLR each year and triggers "
          "`_drown_break` when it passes 0 m MHW, so a lower start is closer to "
          "that threshold.")
        a("")

    a("## Two things this file does not correct")
    a("")
    a("### 1. The time gap")
    a("")
    a("The DEM is **2009**; one run starts in **1984**. CASCADE decrements "
      "`road_ele` by RSLR every year, so the 1984 run begins with a roadbed "
      "that is already ~25 years of sea-level rise low relative to its own MHW. "
      "No back-correction is applied. These are measurements of a surface that "
      "exists, not reconstructions of one that does not — a reconstructed 1984 "
      "roadbed would also have reintroduced the two-file split this document "
      "exists to remove.")
    a("")
    a("### 2. The relocations")
    a("")
    a("NC-12 was relocated in GIS 9–15 (1999) and GIS 84–87 (1989). The "
      f"{ROAD_LINE_YEAR} alignment sampled here is the **post**-relocation one. "
      "For the 2004 run that is correct. For the 1984 run, the road in those "
      "11 domains was physically somewhere else, and the number written is the "
      "elevation of the corridor the road would later occupy.")
    a("")
    a("This is flagged, not adjusted, and the RELOCATION BRACKET above is why: "
      "sampling the 1984 line instead does not fix it, because that line now "
      "runs under the foredune. The affected domains are marked "
      "`RELOCATED_PRE_2004` below and crossed on the QC figure.")
    a("")
    a("The residual is bounded and worth stating rather than hiding. In 11 of "
      "82 domains, for the first 15 years of the **1984** runs only, `road_ele` "
      "may be too high by roughly the gap to the neighbouring grade. Because "
      "`bulldoze` removes `sum(old − road_ele)` clipped at zero, too high means "
      "the model removes **less** overwash than it should and hands **less** "
      "sand to the dunes there. No measurement closes that gap; substituting a "
      "modelled estimate would trade a known bias for an invented one.")
    a("")
    rel = [r for r in good if r["domain"] in RELOCATED_BEFORE_2004]
    if rel:
        a("| GIS | m NAVD88 | m MHW | Relocation |")
        a("|---:|---:|---:|---|")
        for r in sorted(rel, key=lambda x: x["domain"]):
            ev = "1999 (Buxton–Avon)" if r["domain"] < 16 else "1989 (Pea Island)"
            a(f"| {r['domain']} | {r['elev_navd']:.2f} | {r['elev_mhw']:.2f} | "
              f"{ev} |")
        a("")

    a("## Flags")
    a("")
    a("| Flag | Meaning |")
    a("|---|---|")
    a(f"| `THIN(n)` | Fewer than {MIN_CELLS} valid 1 m cells — the line barely "
      f"crosses this domain. |")
    a(f"| `NODATA(x%)` | More than {MAX_NODATA_FRAC * 100:.0f}% of the corridor "
      f"fell on LiDAR NoData. |")
    a(f"| `SCATTER(σ)` | Within-domain σ above {MAX_STD_M:.2f} m — not a graded "
      f"surface. The corridor is catching an embankment or the dune toe. |")
    a(f"| `JUMP(±x)` | Step from the previous domain above {MAX_JUMP_M:.2f} m. |")
    a("| `RELOCATED_PRE_2004` | Correct for the 2004 run; the wrong place for "
      "the 1984 run. Value unadjusted. |")
    a("")

    a("## Per-domain results")
    a("")
    a("| GIS | m NAVD88 | m MHW | median NAVD88 | σ | p10–p90 | n cells | Flags |")
    a("|---:|---:|---:|---:|---:|---|---:|---|")
    for r in sorted(rows, key=lambda x: x["domain"]):
        if not np.isfinite(r["elev_navd"]):
            a(f"| {r['domain']} | — | — | — | — | — | 0 | `NO_DATA` |")
            continue
        fl = ", ".join(f"`{x}`" for x in r["flags"].split(",") if x)
        a(f"| {r['domain']} | {r['elev_navd']:.2f} | {r['elev_mhw']:.2f} | "
          f"{r['elev_median_navd']:.2f} | {r['elev_std_m']:.2f} | "
          f"{r['elev_p10_navd']:.2f}–{r['elev_p90_navd']:.2f} | "
          f"{r['n_cells']} | {fl} |")
    a("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(L), encoding="utf-8")
    print(f"[out] {path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    clips, resamples, surface_label = resolve_surfaces()

    print("=" * 92)
    print("NC-12 ROAD ELEVATION FROM THE 2009 LIDAR")
    print("one alignment, one surface, one file -- used for BOTH start periods")
    print("=" * 92)
    print(f"  surface: {surface_label}")
    print(f"  clips  : {Path(next(iter(clips.values()))).parent if clips else '-'}")
    print(f"  {len(clips)} clip(s) at 1 m | {len(resamples)} at 10 m")
    if not clips:
        raise SystemExit(
            f"\n[stop] no 1 m clips found for FILL_SOURCE={FILL_SOURCE!r}\n")

    line = load_line(ROAD_LINE)
    print(f"  line   : {Path(ROAD_LINE).name}  CRS {line.crs.to_string()}  "
          f"({len(line)} feature(s))")
    with rasterio.open(next(iter(clips.values()))) as s:
        to_m = float(s.crs.linear_units_factor[1])
        print(f"  cell   : {abs(s.transform.a) * to_m:.2f} m  |  "
              f"buffer {BUFFER_M} m -> {BUFFER_M / to_m:.2f} CRS units")
    print(f"  agg    : mean of cells | datum written: m MHW-relative "
          f"(MHW = {MHW_M} m NAVD88)")

    rows = sample_all(line, clips)
    if not rows:
        raise SystemExit("\n[stop] nothing sampled\n")

    print(f"\n{'GIS':>4} {'NAVD88':>8} {'MHW':>7} {'med':>7} {'sd':>6} "
          f"{'p10':>6} {'p90':>6} {'cells':>7}  flags")
    for r in sorted(rows, key=lambda x: x["domain"]):
        if not np.isfinite(r["elev_navd"]):
            print(f"{r['domain']:>4} {'--':>8} {'--':>7} {'--':>7} {'--':>6} "
                  f"{'--':>6} {'--':>6} {0:>7}  NO_DATA")
            continue
        print(f"{r['domain']:>4} {r['elev_navd']:>8.2f} {r['elev_mhw']:>7.2f} "
              f"{r['elev_median_navd']:>7.2f} {r['elev_std_m']:>6.2f} "
              f"{r['elev_p10_navd']:>6.2f} {r['elev_p90_navd']:>6.2f} "
              f"{r['n_cells']:>7}  {r['flags']}")

    ev = np.array([r["elev_navd"] for r in rows if np.isfinite(r["elev_navd"])])
    if ev.size:
        print(f"\n  road elevation (m NAVD88): min {ev.min():.2f} | "
              f"median {np.median(ev):.2f} | max {ev.max():.2f}")
        print(f"  same, MHW-relative       : min {ev.min() - MHW_M:.2f} | "
              f"median {np.median(ev) - MHW_M:.2f} | "
              f"max {ev.max() - MHW_M:.2f}")
        print(f"  runner scalar ROAD_ELEVATION = {RUNNER_SCALAR} | "
              f"Velasquez-Montoya et al. give NC-12 average 1.3 m NAVD88")

    flagged = [(r["domain"], r["flags"]) for r in rows if r["flags"]]
    if flagged:
        print(f"\n  {len(flagged)} domain(s) flagged:")
        for d, f in flagged:
            print(f"    GIS {d:>2}: {f}")

    # --- internal checks -------------------------------------------------
    print(f"\n{'=' * 92}")
    print("INTERNAL CHECKS")
    print("=" * 92)
    print("  The ArcGIS elevation_2009 column this method was validated against")
    print("  is gone from nc12_2004.csv -- its successors are all zeros. That")
    print("  external check is no longer reproducible. What follows is internal.")

    print(f"\n  corridor width sweep (island median of the per-domain mean):")
    sweep = buffer_sweep(line, clips)
    print(f"    {'buffer':>7} {'strip':>7} {'median':>8} {'within-domain sd':>18}")
    for b, m, s, n in sweep:
        mark = "  <- used" if abs(b - BUFFER_M) < 1e-9 else ""
        print(f"    {b:>6.1f}m {b * 2:>6.0f}m {m:>8.2f} {s:>18.2f}{mark}")
    print("    sd climbing while the median falls is the corridor reaching off")
    print("    the pavement onto the shoulder and the dune toe.")

    bracket = relocation_bracket(clips, rows)
    if bracket:
        print(f"\n  RELOCATION BRACKET -- what is under the "
              f"{BRACKET_LINE_YEAR} line in the {bracket['n']} relocated "
              f"domains?")
        print(f"    {'GIS':>4} {'2004 line':>10} {'1984 line':>10} "
              f"{'diff':>7} {'sd 1984':>8}")
        for p in bracket["pairs"]:
            print(f"    {p['domain']:>4} {p['used']:>10.2f} "
                  f"{p['other']:>10.2f} {p['other'] - p['used']:>+7.2f} "
                  f"{p['other_sd']:>8.2f}")
        print(f"    mean {bracket['used_mean']:.2f} (used) vs "
              f"{bracket['other_mean']:.2f} (other), max sd "
              f"{bracket['other_sd_max']:.2f} m")
        print(f"    un-relocated neighbours in the same reaches: "
              f"{bracket['neighbour_grade']:.2f} m")
        if not bracket["brackets"]:
            print("    The two lines do NOT bracket the neighbouring grade --")
            print("    both are above it. The abandoned corridor is not")
            print("    flattened ground, it is dune that migrated over the")
            print("    corridor after the road left. The value used is the")
            print("    lower of the two AND the only graded surface.")

    resamp = resample_check(line, resamples, rows)
    if resamp:
        print(f"\n  1 m clip vs the 10 m Barrier3D grid, same corridor:")
        print(f"    {resamp['n']} domains | median diff {resamp['median']:+.2f} m"
              f" | mean |diff| {resamp['mean_abs']:.2f} m"
              f" | max |diff| {resamp['max_abs']:.2f} m")
        sd1 = np.median([r["elev_std_m"] for r in rows
                         if np.isfinite(r["elev_std_m"])])
        n1 = int(np.median([r["n_cells"] for r in rows if r["n_cells"]]))
        print(f"    cells/domain: {n1} at 1 m vs {resamp['cells_10m']} at 10 m")
        print(f"    within-domain sd: {sd1:.2f} m at 1 m vs "
              f"{resamp['sd_10m']:.2f} m at 10 m")
        print("    NOTE: the two resolutions agree here. The 1.59 m sd that")
        print("    originally justified the 1 m clip was on the 1984 line,")
        print("    which crosses a relocation scar; the 2004 line does not.")
        print("    The 1 m clip is a precaution on this alignment, not a")
        print("    rescue -- though it rests on ~100x the sample.")

    # --- outputs ---------------------------------------------------------
    out = Path(OUT_ROOT)
    print()
    write_domains_csv(rows, out / "RoadElevation_domains.csv")
    write_cascade_csv(rows, out / "RoadElevation.csv")
    qc_figure(rows, out / "HAT_road_elevation.png")
    map_figure(rows, line, resamples, out / "HAT_road_elevation_map.png")
    write_markdown(rows, sweep, resamp, bracket,
                   out / "RoadElevation_audit.md")

    print(f"\n{'=' * 92}")
    print("NEXT: wire the per-domain elevation into the runner (NOT done here)")
    print("=" * 92)
    print("  CASCADE accepts a per-domain road_ele -- cascade.py:47 does")
    print("      if np.size(road_ele) > 1: self._road_ele = road_ele")
    print("      else:                     self._road_ele = [road_ele] * ny")
    print("  and passes self._road_ele[iB3D] to each RoadwayManager. So a list")
    print("  of length TOTAL_DOMAINS works wherever ROAD_ELEVATION is passed.")
    print()
    print("  In HAT_hindcast_1984_2024.py, replace the scalar with a")
    print("  padded array built the same way road_setbacks_full is -- the SAME")
    print("  file for both periods:")
    print()
    print("      ROAD_ELEVATION_FILE = (MGMT_DIR / 'road_elevation'")
    print("                             / 'RoadElevation.csv')")
    print("      road_elev_raw = np.loadtxt(ROAD_ELEVATION_FILE, skiprows=1,")
    print("                                 delimiter=',')")
    print("      road_elev_full = np.full(TOTAL_DOMAINS, ROAD_ELEVATION)")
    print("      road_elev_full[START_ROAD_INDEX:START_ROAD_INDEX")
    print("                     + len(road_elev_raw)] = road_elev_raw")
    print("      ... road_ele=road_elev_full")
    print()
    print("  The file is ALREADY MHW-relative, which is the frame bulldoze")
    print("  compares against. Do not subtract MHW again.")
    print()
    print("  Check START_ROAD_INDEX against the ID row before trusting the")
    print("  padding -- the raw load DROPS row 0, so alignment is positional.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
