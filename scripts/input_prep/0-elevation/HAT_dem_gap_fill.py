"""
HAT_dem_gap_fill.py

Step 1 of 3: clips the 2009 DEM to each domain and fills its gaps from the 2014
NOAA Post-Sandy DEM. Writes one gap-filled 1 m clip and one survey-year raster
per domain; HAT_dem_resample_clip.py resamples them to 10 m.

WHAT THIS IS ACTUALLY FIXING
----------------------------
Not "voids in a surface". Measured on the real DEM, each domain contains exactly
TWO nodata regions and both touch the domain edge - interior enclosed nodata is
0 cells, 0.00%. With OCEAN_LOC="right" in HAT_dune_topo_extractor.py, the east
region (~480 m) is the Atlantic and the west region (~1045 m) is Pamlico Sound.
The 2009 survey simply stops at the waterline on each side.

The gap that matters is the sound-side margin, and HAT_dune_topo_extractor.py
already documents why (lines 281-297): roadway_manager.bulldoze drowns a roadway
when >20% of the cells BORDERING it sit at or below 0 m MHW, and a no-data cell
passes that test. In GIS 78/79/80 the row landward of NC-12 is 17-25 no-data
cells and ZERO genuinely wet ones, so all three roadways width-drowned at t=0 on
missing coverage alone - while the profiles were still 0.5-0.7 m ABOVE MHW.

So this fills measured ground the 2009 survey missed. It does not invent
elevation anywhere.

THE FOUR RULES THAT BOUND THE FILL
-----------------------------------
1. COVERAGE.   Only cells the 2014 DEM actually has a value for. It covers
               97.34% of the island's DRY-LAND gaps - cells 2009 missed where
               the consensus of candidate DEMs puts the ground above MHW. See
               FILL_DEM_PATH for the full scoring and why the alternatives lost.
2. CONNECTIVITY. Only cells contiguous with the island's valid 2009 surface, so
               detached marsh hummocks and any water returns SMRF kept out in
               the sound cannot become new land in the barrier interior. This
               cannot exclude the target cells: the fringe landward of NC-12 is
               contiguous with the island by definition.
               Computed on a BUFFERED window - on the bare 500 m strip, marsh
               that connects to the island just outside the domain would be
               severed by the crop.
3. ELEVATION.  A guard at the DOWNSTREAM water threshold (-3.0 m MHW, the
               extractor's WATER_CLAMP_M), not at MHW. Flooring at MHW would
               discard real low marsh and would not have protected the drowning
               fix anyway - see the note on FILL_MIN_ELEV_NAVD. Sub-MHW fills
               are counted, not rejected.
4. VERTICAL.   Nothing is applied. Both bias correction and feathering are OFF,
               so a filled cell is the 2014 measurement unchanged. Both bias
               estimates are still computed and written to the audit every run.

Everything each rule rejects is counted per domain in the audit CSV. Nothing is
dropped silently.

THE SEAM IS REAL, KNOWN, AND DELIBERATELY LEFT IN
--------------------------------------------------
Where the fill meets measured 2009 ground there is a step. The numbers below
were measured on the SUPERSEDED 2008 point-cloud attempt; seam_check() re-runs
every time, so the current source's figures are in the audit CSV
(seam_median_abs_m against ctrl_median_abs_m). The reasoning for leaving it
uncorrected carries over:

  measured <-> measured, whole domain      0.028 m   terrain roughness
  measured <-> measured, near the boundary 0.028 m   same - the margin is smooth
  fill <-> measured (what we ship)         0.341 m   ~12x either control
  2008 <-> 2008 across the same boundary   0.030 m   <- the decisive one

The last row grids 2008 on BOTH sides of the boundary, so any inter-survey
offset cancels. It comes out flat, at the roughness floor. The ground therefore
does NOT drop at the marsh edge: terrain accounts for ~9% of the step and the
rest is the two surveys disagreeing where they meet.

The obvious fix - shift 2008 onto 2009 per domain - was rejected because the
disagreement is not a datum offset and does not have a consistent sign. Signed
step by domain: -0.685, +0.136, -0.787, -0.117, +0.277, +0.455, +0.109. The fill
is too low in some domains and too high in others, and the boundary estimate
(+0.25 median) disagrees with the whole-domain overlap (-0.05 median). Any
single-number shift removes the seam at the boundary and introduces a comparable
disagreement across the interior instead - trading a visible artifact for an
invisible one.

Feathering was rejected for the same reason: it hides the step over 5 m without
addressing the disagreement, and smooths measured data to do it.

So the step stays, and clip_domain_<N>_survey.tif marks exactly which cells came
from 2008 so any consumer can find it. At the 10 m Barrier3D grid it collapses
to one cell boundary - 0.3 m over 10 m, ~3% slope, within the range of real
back-barrier relief, but it is fabricated and worth knowing about when reading
overwash behaviour near a fill margin.

INPUTS
    D:/Hatteras_GIS/.../2009_full.tif   base, 1 m, EPSG:3725 + NAVD88
    D:/Hatteras_GIS/.../2014_full.tif   fill, 1 m, EPSG:6347 + NAVD88
                                        (reprojected on read by WarpedVRT)
    D:/Hatteras_GIS/domains.geojson     90 boxes, 2000 x 500 m

OUTPUTS (data/hatteras_init/0-elevation/1-gapfill-1m/)
    clip_domain_<N>_filled.tif      gap-filled clip, m NAVD88
    clip_domain_<N>_survey.tif      which survey each cell came from:
                                        2009 = measured by the 2009 DEM
                                        2014 = filled from the 2014 DEM
                                           0 = neither survey saw it
    gapfill_audit.csv               per-domain counts for every rule above

Filenames follow the legacy 2009-domain-clipresample convention
(clip_domain_N.tif) with _filled marking the new set. Layout is flat rather than
per-domain subfolders: step 1 writes the 1 m clips and step 2 writes the 10 m
ones, so one folder per step means re-running a step is a single delete.

Requires: rasterio, geopandas, numpy, scipy
"""

from __future__ import annotations

import csv
import os
from pathlib import Path

import numpy as np
import geopandas as gpd
import rasterio
from pyproj import CRS
from rasterio.enums import Resampling
from rasterio.vrt import WarpedVRT
from rasterio.windows import Window, from_bounds
from scipy.interpolate import griddata
from scipy.ndimage import binary_dilation, distance_transform_edt, label

# =============================================================================
# CONFIG
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parents[3]
ELEVATION_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "0-elevation"
GIS_ROOT = Path(r"D:\Hatteras_GIS")

BASE_DEM_PATH = (GIS_ROOT / "Elevation" / "Polygons" / "2009"
                 / "usace2009_nc_dem_Job1076020" / "2009_full.tif")

# FILL SOURCE: 2014 NOAA Post-Sandy DEM.
#
# Chosen by measurement, not vintage. Every candidate DEM was scored against the
# 2009 DRY-LAND gaps (cells 2009 missed where the consensus of candidates puts
# the ground above MHW), over all 90 domains:
#
#     2014 NCFMP          100.00%   DISQUALIFIED - hydro-flattened, 70.5% of its
#                                   values in the gap are the constant -0.762 m
#     2014 NOAA Post-Sandy 97.34%   <- earliest genuine, and the best
#     2017 USACE           23.74%
#     2016 post-Matthew    21.98%
#     2019 DUNEX           21.98%
#     2018 post-Florence    9.74%
#
# The 2016/2017/2019 collapse is spatial extent: they score 95-100% on domains
# 78-80 but are localised surveys - 2019 is below 50% in 60 of 90 domains. A
# choice made on the developed reaches alone would have picked a dataset
# covering less than a quarter of the island's gaps.
#
# It is a DEM, not a point cloud, so HAT_laz_ground_classify.py does not run for
# it. Different CRS from the base (EPSG:6347 NAD83(2011) vs EPSG:3725
# NAD83(NSRS2007)); WarpedVRT reprojects on read.
#
# Superseded point-cloud attempts, kept for the record:
#   2008 NOAA IOCM  topo-only, ~19% of the sound-side gap, 17/150 in the NC-12 strip
#   2011 post-Irene topo-only, ~35%, 28/150
FILL_DEM_PATH = (GIS_ROOT / "Elevation" / "Polygons" / "2014"
                 / "2014_NOAA_Post_Sandy_DEM_Job1076021" / "2014_full.tif")
FILL_SOURCE_TAG = "2014_NOAA_PostSandy"
FILL_SOURCE_YEAR = 2014

# "dem" or "point_cloud". A DEM is read through WarpedVRT; a point cloud is
# gridded to the base grid by MEDIAN of the returns in each 1 m cell (minimum
# would take the water surface over wet marsh, and cells below 0 m MHW are
# exactly what the drowning test counts).
FILL_SOURCE_TYPE = "dem"
FILL_POINTS_DIR = GIS_ROOT / "Elevation" / "Points" / "2008" / "ground_smrf"
FILL_POINTS_GLOB = "*_ground_*.laz"

DOMAIN_FILE = GIS_ROOT / "domains.geojson"
DOMAIN_ID_FIELD = "domain_id"

OUTPUT_DIR = ELEVATION_DIR / "1-gapfill-1m" / FILL_SOURCE_TAG
AUDIT_CSV = "gapfill_audit.csv"

GRID_SIZE_M = 10.0    # the eventual Barrier3D cell; the clip must divide by it
EXPECTED_CLIP = (500, 2000)   # rows, cols at 1 m

# --- vertical datum (matches HAT_dune_topo_extractor / hindcast config) ---
MHW_ELEVATION = 0.36          # m NAVD88, Duck NC gauge 8651370

# --- RULE 3: elevation floor ---
# Set to the DOWNSTREAM threshold, not to MHW, and it is a guard rather than a
# filter. Reasoning, because this was reversed once already:
#
#   * HAT_dune_topo_extractor.py does the MHW referencing itself (line 1017,
#     z = raw - 0.36) and applies its own water threshold at WATER_CLAMP_M =
#     -3.0 m MHW = -2.64 m NAVD88. That -3.0 was picked deliberately -
#     "keeps back-barrier marsh cells (Lexi's v3 edit)", up from -1.0. Flooring
#     at MHW here would undo that decision one step upstream, invisibly.
#   * A floor at MHW never protected the road-drowning fix it was added for.
#     roadway_manager.bulldoze drowns when >20% of bordering cells sit at or
#     below 0 m MHW; a filled cell at -0.2 m MHW counts as wet exactly as the
#     -3.0 sentinel did when it was unsurveyed. The bug was cells genuinely
#     ABOVE MHW reading as wet because nobody surveyed them, and filling those
#     with their true elevation fixes it whatever the floor is.
#   * The 2008 cloud bottoms out at -1.33 m NAVD88, so this floor rejects
#     nothing in practice. It exists to catch a future fill source that could.
#
# Sub-MHW fills are still COUNTED per domain (cand_below_mhw in the audit), so
# lowering the floor costs no visibility.
APPLY_ELEV_FLOOR = True
FILL_MIN_ELEV_NAVD = MHW_ELEVATION - 3.0   # -2.64 m NAVD88, the extractor's clamp

# --- RULE 2: connectivity ---
REQUIRE_ISLAND_CONNECTION = True
CONTEXT_BUFFER_M = 200.0      # window padding for connectivity + boundary context

# Strict connectivity severs a 300 m marsh platform if a 20 m tidal creek that
# is water in BOTH surveys separates it from the island. Gaps up to this width
# are bridged before the connectivity test, so creek-separated back barrier is
# kept while genuinely detached patches out in the sound are still rejected.
# 0.0 disables bridging (strict connectivity).
#
# 20 m, measured rather than guessed. Across 12 sampled domains there were 2969
# detached components holding 1,055,741 candidate cells. Recovery vs bridging
# distance, by component count and by CELL count (cells are what matters - the
# goal is captured area, and one 200k-cell platform outweighs 50 specks):
#
#     bridge   % components   % detached cells
#        2 m           7.1%              3.5%
#        5 m          10.3%             16.2%   <- first step
#       15 m          14.6%             18.6%
#       20 m          16.2%             30.4%   <- second step, then a plateau
#       30 m          18.9%             31.6%
#      100 m          31.6%             39.9%
#      200 m          47.7%             47.6%   <- columns converge
#
# 20 m sits on the plateau right after the second step: 20->30 m buys 1.2 more
# points, 30->100 m buys 8.3 for five times the reach. Below 200 m the cell
# column runs at ~2x the component column, i.e. bridging is selectively catching
# large platforms; by 200 m they converge, which means it has stopped gaining
# area preferentially and is just admitting open sound. 20 m is also a credible
# tidal-creek width, which 100 m is not.
#
# In context of ALL candidates in those domains (4,420,718 cells):
#   strict 76.1%  |  5 m 80.0%  |  20 m 83.4%  |  100 m 85.7%
#
# Note this barely moves the domains that motivated the fill. In GIS 78/79/80
# only ~20k of ~300k candidates are detached at all, so ~93-94% is already
# connected and their detached patches sit at 49-174 m - open water, not creeks.
# Bridging is about how much back barrier domains 10/50/90 carry, not the road.
GAP_BRIDGE_M = 20.0

# --- RULE 4: vertical reconciliation ---
#
# BOTH VALUE-MODIFYING STEPS ARE OFF. Of the rules here, only these two change
# what a cell says; coverage, connectivity and the floor merely select which
# measured cells get used. With both off, a filled cell is the 2008 measurement
# and nothing else.
#
# Bias correction was ON and was misfiring. It estimated a single offset from a
# 10 m collar around the fill, which by construction sits on the 2009 waterline
# - the one place the two surveys disagree for reasons that are not a datum
# offset. min-z gridding there picks water-surface returns, so the collar
# measured "min-z reads low in wet cells" and applied it hundreds of metres into
# dry marsh. Across 28 domains it produced +0.824 to -0.104 m (median +0.125),
# the signature of an unstable estimator rather than real per-domain offsets.
# Full-domain overlap says the surveys actually agree to ~0.03 m.
#
# Both estimates are still COMPUTED and written to the audit every run, so the
# decision to not apply them stays visible and checkable.
APPLY_BIAS_CORRECTION = False
APPLY_FEATHER = False

OVERLAP_RING_PX = 10          # collar used for the (reported) ring estimate
FEATHER_WIDTH_PX = 5          # blend distance, only if APPLY_FEATHER
MIN_RING_CELLS = 25           # below this the ring estimate is not trustworthy

# How the 2008 ground returns become a 1 m raster.
#   "median"  median of ground returns in the cell. Closer to how a gridded DEM
#             surface is built, so it is comparable to the 2009 values it sits
#             beside, and one stray low return cannot set the cell.
#   "min"     lowest return. The classic bare-earth proxy, but over wet marsh
#             the lowest return is often the water surface - and cells below
#             0 m MHW are exactly what the drowning test counts, so this is
#             conservative in the wrong direction for this particular bug.
GRID_STAT = "median"

NODATA_OUT = -9999.0

# The _survey raster stores the year each cell's elevation came from, so it
# needs no legend. 0 means neither survey saw the cell.
SURVEY_2009, SURVEY_NONE = 2009, 0
SURVEY_FILL = FILL_SOURCE_YEAR   # the year written for a filled cell
SURVEY_NODATA = 65535   # unused sentinel; 0 is a real value here, not nodata


# =============================================================================
# GEOMETRY - the domain window, snapped to the DEM's own grid
# =============================================================================

def snap_window(bounds, transform, res_x, res_y, block):
    """
    Polygon bounds -> integer window on the source grid, trimmed to whole
    `block`-sized blocks so the later 10 m resample divides evenly.

    Nearest cell edge, not floor/ceil: a polygon 0.4 m off grid should snap to
    the near edge, not grow the window by a whole cell. Returns the adjustment
    so it is reported rather than absorbed.
    """
    minx, miny, maxx, maxy = bounds
    left, top = transform.c, transform.f

    col0 = int(round((minx - left) / res_x))
    row0 = int(round((top - maxy) / res_y))
    col1 = int(round((maxx - left) / res_x))
    row1 = int(round((top - miny) / res_y))

    width, height = col1 - col0, row1 - row0
    trim_w, trim_h = width % block, height % block
    width -= trim_w
    height -= trim_h

    adj = {"snap_x_m": abs((left + col0 * res_x) - minx),
           "snap_y_m": abs((top - row0 * res_y) - maxy),
           "trim_x_m": trim_w * res_x, "trim_y_m": trim_h * res_y}
    return Window(col0, row0, width, height), adj


def pad_window(win, pad_px):
    return Window(win.col_off - pad_px, win.row_off - pad_px,
                  win.width + 2 * pad_px, win.height + 2 * pad_px)


def read_window(src, win, nodata_in):
    """Boundless so a domain hanging off the DEM edge yields nodata, not an
    error - the domains do overshoot the DEM's east edge by ~6 m."""
    fill = nodata_in if nodata_in is not None else NODATA_OUT
    arr = src.read(1, window=win, boundless=True, fill_value=fill).astype(np.float64)
    if nodata_in is not None and not np.isnan(nodata_in):
        arr = np.where(arr == nodata_in, np.nan, arr)
    return arr


# =============================================================================
# =============================================================================
# THE FILL SOURCE DEM
# =============================================================================

class PointCloudSource:
    """
    A ground-classified point cloud as a fill source, gridded to the base grid.

    Kept alongside FillSource so a superseded point-cloud run can be reproduced
    - which is the only way to re-render its figures after the rasters have been
    deleted. Tiles are indexed by header bbox so each domain reads only the
    cells it overlaps, and the last few are cached since consecutive domains
    reuse them.
    """

    def __init__(self, points_dir, pattern, cache_size=4):
        import laspy
        self._laspy = laspy
        self.tiles = []
        for f in sorted(Path(points_dir).glob(pattern)):
            with laspy.open(str(f)) as fh:
                h = fh.header
                if h.point_count == 0:
                    continue
                self.tiles.append((f, (h.mins[0], h.mins[1], h.maxs[0], h.maxs[1])))
        if not self.tiles:
            raise FileNotFoundError(f"no {pattern} in {points_dir}")
        self.epsg, self.res, self.nodata = None, 1.0, None
        self.cache, self.order, self.cache_size = {}, [], cache_size

    def _load(self, f):
        if f not in self.cache:
            las = self._laspy.read(str(f))
            self.cache[f] = (np.asarray(las.x), np.asarray(las.y), np.asarray(las.z))
            self.order.append(f)
            while len(self.order) > self.cache_size:
                self.cache.pop(self.order.pop(0), None)
        return self.cache[f]

    def read_on_grid(self, bounds, shape):
        import pandas as pd
        bx0, by0, bx1, by1 = bounds
        xs, ys, zs = [], [], []
        for f, (tx0, ty0, tx1, ty1) in self.tiles:
            if tx1 < bx0 or tx0 > bx1 or ty1 < by0 or ty0 > by1:
                continue
            x, y, z = self._load(f)
            m = (x >= bx0) & (x < bx1) & (y >= by0) & (y < by1)
            if m.any():
                xs.append(x[m]); ys.append(y[m]); zs.append(z[m])
        out = np.full(shape, np.nan, np.float32)
        if not xs:
            return out
        x = np.concatenate(xs); y = np.concatenate(ys); z = np.concatenate(zs)
        col = ((x - bx0) / self.res).astype(np.int64)
        row = ((by1 - y) / self.res).astype(np.int64)
        ok = (row >= 0) & (row < shape[0]) & (col >= 0) & (col < shape[1])
        if not ok.any():
            return out
        flat = row[ok].astype(np.int64) * shape[1] + col[ok]
        med = pd.Series(z[ok]).groupby(flat).median()
        out.ravel()[med.index.to_numpy()] = med.to_numpy()
        return out

    def close(self):
        self.cache.clear(); self.order.clear()


class FillSource:
    """
    The candidate DEM, read onto the base DEM's grid on demand.

    WarpedVRT handles the CRS difference (EPSG:6347 -> EPSG:3725) and any
    resolution difference on read, so nothing downstream needs to know the
    source is in a different realisation of NAD83.

    Nearest-neighbour resampling deliberately: this fills cells the 2009 survey
    missed, and interpolating would invent values at the very edges where the
    two surveys meet, which is where they are least comparable.
    """

    def __init__(self, path, dst_crs):
        self.src = rasterio.open(path)
        self.vrt = WarpedVRT(self.src, crs=dst_crs,
                             resampling=Resampling.nearest)
        crs = CRS.from_wkt(self.src.crs.to_wkt()) if self.src.crs else None
        hz = crs.sub_crs_list[0] if (crs is not None and crs.is_compound) else crs
        self.epsg = hz.to_epsg() if hz is not None else None
        self.res = self.src.transform.a
        self.nodata = self.src.nodata

    def read_on_grid(self, bounds, shape):
        """WarpedVRT rejects boundless reads, so partial overlap is intersected
        and pasted at the right offset rather than erroring."""
        bx0, by0, bx1, by1 = bounds
        out = np.full(shape, np.nan, np.float32)
        vb = self.vrt.bounds
        ix0, iy0 = max(bx0, vb.left), max(by0, vb.bottom)
        ix1, iy1 = min(bx1, vb.right), min(by1, vb.top)
        if ix1 <= ix0 or iy1 <= iy0:
            return out
        h, w = int(round(iy1 - iy0)), int(round(ix1 - ix0))
        if h <= 0 or w <= 0:
            return out
        d = self.vrt.read(1, window=from_bounds(ix0, iy0, ix1, iy1,
                                                transform=self.vrt.transform),
                          out_shape=(h, w), masked=True,
                          resampling=Resampling.nearest)
        d = np.ma.filled(d.astype(np.float32), np.nan)
        r0, c0 = max(int(round(by1 - iy1)), 0), max(int(round(ix0 - bx0)), 0)
        hh, ww = min(h, shape[0] - r0), min(w, shape[1] - c0)
        if hh > 0 and ww > 0:
            out[r0:r0 + hh, c0:c0 + ww] = d[:hh, :ww]
        return out

    def close(self):
        self.vrt.close()
        self.src.close()


# =============================================================================
# THE FILL
# =============================================================================

def island_connected(valid, candidate, bridge_px=0):
    """
    Candidate cells reachable from the island through valid ground or other
    candidates. The island is the largest connected component of valid 2009
    cells in the (buffered) window.

    bridge_px dilates the land mask before the reachability test, so a gap up
    to 2*bridge_px wide is crossed. Dilation is used ONLY to decide
    reachability - the returned mask is still a subset of the real candidates,
    so no cell is invented by bridging.
    """
    if not valid.any():
        return np.zeros_like(candidate)
    lab_v, n_v = label(valid, structure=np.ones((3, 3)))
    if n_v == 0:
        return np.zeros_like(candidate)
    sizes = np.bincount(lab_v.ravel()); sizes[0] = 0
    island = lab_v == int(np.argmax(sizes))

    land = valid | candidate
    if bridge_px > 0:
        land = binary_dilation(land, iterations=int(bridge_px))
    lab_a, _ = label(land, structure=np.ones((3, 3)))
    keep = set(np.unique(lab_a[island])) - {0}
    return candidate & np.isin(lab_a, list(keep))


def _neighbour_diffs(z, mask_a, mask_b):
    """z[a] - z[b] over every 4-connected pair with a in mask_a, b in mask_b."""
    out = []
    for sl_a, sl_b in ((np.s_[:, :-1], np.s_[:, 1:]),
                       (np.s_[:, 1:], np.s_[:, :-1]),
                       (np.s_[:-1, :], np.s_[1:, :]),
                       (np.s_[1:, :], np.s_[:-1, :])):
        m = mask_a[sl_a] & mask_b[sl_b]
        if m.any():
            out.append(z[sl_a][m] - z[sl_b][m])
    return np.concatenate(out) if out else np.empty(0)


def seam_check(z, measured, filled):
    """
    How big a step does the fill create where it meets measured 2009 ground?

    Reported against a CONTROL: the same statistic between adjacent measured
    cells. Real terrain is not flat, so a seam step only means something
    relative to how much neighbouring cells normally differ. seam ~ control
    means the fill is indistinguishable from the surface it joins, and no
    feathering is warranted. seam >> control is a genuine cliff.
    """
    seam = _neighbour_diffs(z, filled, measured)
    ctrl = _neighbour_diffs(z, measured, measured)
    # Every key present every time, so the audit CSV has stable columns even
    # for a domain with no fill at all.
    res = {"seam_n": int(seam.size), "ctrl_n": int(ctrl.size),
           "seam_median_signed_m": None, "seam_median_abs_m": None,
           "seam_p90_abs_m": None, "ctrl_median_abs_m": None,
           "ctrl_p90_abs_m": None, "seam_over_ctrl": None}
    if seam.size:
        res["seam_median_signed_m"] = float(np.median(seam))
        res["seam_median_abs_m"] = float(np.median(np.abs(seam)))
        res["seam_p90_abs_m"] = float(np.percentile(np.abs(seam), 90))
    if ctrl.size:
        res["ctrl_median_abs_m"] = float(np.median(np.abs(ctrl)))
        res["ctrl_p90_abs_m"] = float(np.percentile(np.abs(ctrl), 90))
    if res["seam_median_abs_m"] is not None and res["ctrl_median_abs_m"]:
        res["seam_over_ctrl"] = round(
            res["seam_median_abs_m"] / res["ctrl_median_abs_m"], 3)
    return res


def estimate_bias(base, fill, fill_mask, ring_px):
    ring = binary_dilation(fill_mask, iterations=ring_px) & ~fill_mask
    both = ring & ~np.isnan(base) & ~np.isnan(fill)
    n = int(both.sum())
    if n < MIN_RING_CELLS:
        return 0.0, n
    return float(np.nanmedian(base[both] - fill[both])), n


def boundary_extrapolation(base, fill_mask, buffer_px):
    """Locally-consistent continuation of the 2009 surface into the fill area,
    used only as the blend target at the seam so the merge has no hard step."""
    dil = binary_dilation(fill_mask, iterations=buffer_px)
    border = dil & ~fill_mask & ~np.isnan(base)
    br, bc = np.where(border)
    fr, fc = np.where(fill_mask)
    out = np.full(base.shape, np.nan)
    if len(br) < 10 or len(fr) == 0:
        return out
    vals = griddata((br, bc), base[br, bc], (fr, fc), method="linear")
    bad = np.isnan(vals)
    if bad.any():
        vals[bad] = griddata((br, bc), base[br, bc], (fr[bad], fc[bad]),
                             method="nearest")
    out[fr, fc] = vals
    return out


def feathered_merge(base, fill_vals, fill_mask, extrap, feather_px):
    dist = distance_transform_edt(fill_mask)
    w = np.clip(dist / max(feather_px, 1e-9), 0, 1)
    blended = np.where(np.isnan(extrap), fill_vals,
                       (1 - w) * extrap + w * fill_vals)
    out = base.copy()
    out[fill_mask] = blended[fill_mask]
    return out


# =============================================================================
# IO
# =============================================================================

def write_raster(arr, transform, crs, path, dtype, nodata):
    profile = {"driver": "GTiff", "height": arr.shape[0], "width": arr.shape[1],
               "count": 1, "dtype": dtype, "crs": crs, "transform": transform,
               "nodata": nodata, "compress": "deflate"}
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype), 1)


def resolve_crs(src, gdf):
    """The DEM is a COMPOUND CRS (EPSG:3725 + NAVD88), whose to_epsg() is None,
    so a plain `gdf.crs != src.crs` reports a reprojection that is a no-op."""
    if src.crs is None:
        print("  WARNING: DEM has no CRS; assuming domains already match.")
        return gdf
    dem_crs = CRS.from_wkt(src.crs.to_wkt())
    hz = dem_crs.sub_crs_list[0] if dem_crs.is_compound else dem_crs
    if dem_crs.is_compound:
        print(f"  DEM CRS compound: horizontal EPSG:{hz.to_epsg()}, "
              f"vertical {dem_crs.sub_crs_list[1].name}")
    if gdf.crs is None:
        print("  WARNING: domains have no CRS; assuming they match the DEM.")
    elif gdf.crs.to_epsg() != hz.to_epsg():
        print(f"  Reprojecting domains {gdf.crs} -> EPSG:{hz.to_epsg()}")
        gdf = gdf.to_crs(hz)
    else:
        print(f"  Domains already EPSG:{hz.to_epsg()} - no reprojection")
    return gdf


# =============================================================================
# MAIN
# =============================================================================

def main():
    _src = FILL_POINTS_DIR if FILL_SOURCE_TYPE == "point_cloud" else FILL_DEM_PATH
    for label_, path in (("base DEM", BASE_DEM_PATH),
                         ("domain file", DOMAIN_FILE),
                         ("fill source", _src)):
        if not path.exists():
            raise FileNotFoundError(
                f"{label_} not found: {path}\n"
                f"  Fix the paths at the top of this script.")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    src = rasterio.open(BASE_DEM_PATH)
    res_x, res_y = src.transform.a, -src.transform.e
    print(f"\n2009 DEM: {src.height} x {src.width}, {res_x} m, nodata={src.nodata}")
    block = int(round(GRID_SIZE_M / res_x))
    pad_px = int(round(CONTEXT_BUFFER_M / res_x))

    gdf = gpd.read_file(DOMAIN_FILE)
    print(f"\n{len(gdf)} domains from {DOMAIN_FILE.name}")
    gdf = resolve_crs(src, gdf)

    if FILL_SOURCE_TYPE == "point_cloud":
        fill = PointCloudSource(FILL_POINTS_DIR, FILL_POINTS_GLOB)
        print(f"\nfill source: {FILL_POINTS_DIR.name}  ({FILL_SOURCE_TAG})")
        print(f"  {len(fill.tiles)} cell file(s), gridded to 1 m by median")
    elif FILL_SOURCE_TYPE == "dem":
        fill = FillSource(FILL_DEM_PATH, src.crs)
        print(f"\nfill source: {FILL_DEM_PATH.name}  ({FILL_SOURCE_TAG})")
        print(f"  EPSG:{fill.epsg}  {fill.res} m  nodata={fill.nodata}"
              f"{'  -> reprojected on read' if fill.epsg else ''}")
    else:
        raise ValueError(f"FILL_SOURCE_TYPE must be 'dem' or 'point_cloud', "
                         f"got {FILL_SOURCE_TYPE!r}")

    print(f"\nRules: coverage={FILL_SOURCE_YEAR} {FILL_SOURCE_TYPE} | "
          f"connectivity={REQUIRE_ISLAND_CONNECTION} (bridge {GAP_BRIDGE_M:g} m)"
          f" | floor={'%.2f m NAVD88' % FILL_MIN_ELEV_NAVD if APPLY_ELEV_FLOOR else 'off'}"
          f" | bias={APPLY_BIAS_CORRECTION} feather={APPLY_FEATHER}")

    audit = []
    for _, row in gdf.iterrows():
        dom = row[DOMAIN_ID_FIELD]
        try:
            dom = int(dom)
        except (TypeError, ValueError):
            pass

        win, adj = snap_window(row.geometry.bounds, src.transform, res_x, res_y, block)
        big = pad_window(win, pad_px)

        base_b = read_window(src, big, src.nodata)
        t_big = src.window_transform(big)

        bx0, by1 = t_big.c, t_big.f
        bx1 = bx0 + big.width * res_x
        by0 = by1 - big.height * res_y
        g2008 = fill.read_on_grid((bx0, by0, bx1, by1), base_b.shape)

        gap = np.isnan(base_b)
        valid = ~gap
        has08 = ~np.isnan(g2008)

        cand = gap & has08
        n_cand = int(cand.sum())

        if REQUIRE_ISLAND_CONNECTION:
            # bridging a gap of GAP_BRIDGE_M needs half that dilation per side
            bridge_px = int(round(GAP_BRIDGE_M / res_x / 2.0))
            conn = island_connected(valid, cand, bridge_px)
        else:
            conn = cand
        n_conn_drop = n_cand - int(conn.sum())

        # Both estimates are computed and reported every run even when nothing
        # is applied, so "we chose not to correct" stays a checkable claim.
        bias_ring, n_ring = estimate_bias(base_b, g2008, conn, OVERLAP_RING_PX)
        ov = valid & has08
        n_ov = int(ov.sum())
        bias_overlap = (float(np.nanmedian(base_b[ov] - g2008[ov]))
                        if n_ov >= MIN_RING_CELLS else 0.0)

        bias = bias_ring if APPLY_BIAS_CORRECTION else 0.0
        fill_vals = g2008 + bias

        n_floor_drop = 0
        # Counted on the CROPPED domain window, not the padded context window,
        # so it is comparable with `filled` in the audit. Counting it on the
        # padded window made this column exceed `filled` and read as >100%.
        _c = np.s_[pad_px:pad_px + win.height, pad_px:pad_px + win.width]
        n_sub_mhw = int((conn[_c] & (fill_vals[_c] < MHW_ELEVATION)).sum())
        if APPLY_ELEV_FLOOR:
            too_low = conn & (fill_vals < FILL_MIN_ELEV_NAVD)
            n_floor_drop = int(too_low.sum())
            conn = conn & ~too_low

        n_fill = int(conn.sum())
        if n_fill and APPLY_FEATHER:
            extrap = boundary_extrapolation(base_b, conn, pad_px)
            merged = feathered_merge(base_b, fill_vals, conn, extrap, FEATHER_WIDTH_PX)
        elif n_fill:
            merged = base_b.copy()
            merged[conn] = fill_vals[conn]
        else:
            merged = base_b

        seam = seam_check(merged, valid, conn)

        # back to the exact domain window
        r0 = pad_px; c0 = pad_px
        out = merged[r0:r0 + win.height, c0:c0 + win.width]
        filled_mask = conn[r0:r0 + win.height, c0:c0 + win.width]
        gap_out = np.isnan(out)

        survey = np.full(out.shape, SURVEY_2009, np.uint16)
        survey[filled_mask] = SURVEY_FILL
        survey[gap_out] = SURVEY_NONE

        if out.shape != EXPECTED_CLIP:
            print(f"  domain {dom}: ERROR clip {out.shape}, expected {EXPECTED_CLIP}")

        t_out = src.window_transform(win)
        write_raster(np.where(np.isnan(out), NODATA_OUT, out), t_out, src.crs,
                     OUTPUT_DIR / f"clip_domain_{dom}_filled.tif",
                     "float32", NODATA_OUT)
        write_raster(survey, t_out, src.crs,
                     OUTPUT_DIR / f"clip_domain_{dom}_survey.tif",
                     "uint16", SURVEY_NODATA)

        gap_before = int(np.isnan(base_b[r0:r0 + win.height, c0:c0 + win.width]).sum())
        n_filled_out = int(filled_mask.sum())
        print(f"  domain {dom:>3}: nodata {gap_before:6d} -> {int(gap_out.sum()):6d}  "
              f"filled {n_filled_out:5d}  "
              f"seam {seam.get('seam_median_abs_m', float('nan')):.3f} vs ctrl "
              f"{seam.get('ctrl_median_abs_m', float('nan')):.3f} m "
              f"(x{seam.get('seam_over_ctrl', float('nan')):.2f})  "
              f"bias ring {bias_ring:+.3f} / overlap {bias_overlap:+.3f}  "
              f"drop: conn {n_conn_drop}, floor {n_floor_drop}", flush=True)

        audit.append({
            "domain": dom, "rows": out.shape[0], "cols": out.shape[1],
            "nodata_before": gap_before, "nodata_after": int(gap_out.sum()),
            "filled": n_filled_out,
            "cand_2008_coverage": n_cand,
            "dropped_connectivity": n_conn_drop,
            "dropped_elev_floor": n_floor_drop,
            "cand_below_mhw": n_sub_mhw,
            "bias_applied_m": round(bias, 4),
            "bias_ring_m": round(bias_ring, 4), "bias_ring_cells": n_ring,
            "bias_overlap_m": round(bias_overlap, 4), "bias_overlap_cells": n_ov,
            **{k: (round(v, 4) if isinstance(v, float) else v)
               for k, v in seam.items()},
            "fill_cells_in_window": int(has08.sum()),
            "snap_x_m": round(adj["snap_x_m"], 4), "snap_y_m": round(adj["snap_y_m"], 4),
            "trim_x_m": adj["trim_x_m"], "trim_y_m": adj["trim_y_m"],
        })

    src.close()
    fill.close()

    path = OUTPUT_DIR / AUDIT_CSV
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(audit[0].keys()))
        w.writeheader(); w.writerows(audit)

    tot_before = sum(a["nodata_before"] for a in audit)
    tot_after = sum(a["nodata_after"] for a in audit)
    tot_fill = sum(a["filled"] for a in audit)
    tot_conn = sum(a["dropped_connectivity"] for a in audit)
    tot_floor = sum(a["dropped_elev_floor"] for a in audit)
    print(f"\n{'=' * 70}\n{len(audit)} domains")
    print(f"  nodata cells   {tot_before:,} -> {tot_after:,}  "
          f"({100 * (tot_before - tot_after) / max(tot_before, 1):.1f}% recovered)")
    print(f"  filled         {tot_fill:,}")
    print(f"  rejected by connectivity {tot_conn:,}   by elevation floor {tot_floor:,}")
    print(f"  audit: {path}")
    print("\nNext: HAT_dem_resample_clip.py resamples these 1 m clips to 10 m.")


if __name__ == "__main__":
    main()
