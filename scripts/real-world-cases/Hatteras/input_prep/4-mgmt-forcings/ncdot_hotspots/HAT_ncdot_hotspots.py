"""
NCDOT's NC 12 hot spots on the model domains: the three in the reach, georeferenced from NCDOT's own maps.

    python scripts/real-world-cases/Hatteras/input_prep/4-mgmt-forcings/ncdot_hotspots/HAT_ncdot_hotspots.py

NCDOT Division 1's "NC 12 Hot Spots" (Midgett, June 2021) numbers seven hot spots from Oregon Inlet to Ocracoke;
three lie in the reach: 3 Rodanthe 'S' Curves, 4 Avon, 5 Buxton. Its overview map (PDF page 2) marks each with
the same 3.8 km oval, so an oval places a hot spot but does not measure it. The chosen ranges take the oval and
trim it where NCDOT's text or detailed maps say the hot spot ends sooner:

    5 Buxton            GIS 8-14   starts at the Cape Hatteras Motel / NC-12 bend at Old Lighthouse Road (the text;
                                   the Alternative 5 nourishment map, PDF p. 18, georeferenced), ends where the oval
                                   ends; the nourishment line runs on to 16, a design length, not the hot spot
    4 Avon              GIS 28-35  the oval: the 2021 report has no detailed Avon map
    3 Rodanthe S-Curves GIS 84-88  the pavement removed when the Jug Handle Bridge opened (the construction plan,
                                   PDF p. 10, vector), inside the oval

Every georeferencing is a similarity transform onto the 2008 NC-12 line and the island outline, its scale fixed
from the map's own scale bar or dimension labels where it has one; the fit statistics go to
georeference_fits.csv. The 2015 Buxton-to-Avon feasibility study (R-4070B, Figure 1) is georeferenced as a
cross-check: it places the Avon end of its one Buxton-Avon hot-spot project in south Avon (Nino Road, GIS ~26),
south of the 2021 oval; the README discusses it.

Reads the two PDFs in data/hatteras_init/4-mgmt-forcing/ncdot_hotspots/sources/, the domain boxes, the island
outline and the 2008 NC-12 line. Writes, in data/hatteras_init/4-mgmt-forcing/ncdot_hotspots/:
    ncdot_hotspot_ranges.csv         the three ranges, the oval spans, the evidence
    ncdot_hotspots_georef.geojson    every georeferenced feature, UTM 18N
    georeference_fits.csv            scale, rotation and residual of each map's fit
    ncdot_hotspots_map.png           the hot spots and ovals on the domains (PDF in supporting/)
    ncdot_hotspots_georeference.png  the evidence behind each range (PDF in supporting/)
    supporting/CAPTIONS.md

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-13
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_HAT_DIR = next(_p for _p in _HERE.parents if (_p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in sys.path:
    sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as PROJECT_ROOT  # noqa: E402
import site_layer.hat_figure_style as fs  # noqa: E402
from site_layer import hat_topo_version as tv, hat_map_layers as ml  # noqa: E402
from site_layer.hat_observed_rates import DOMAIN_BOXES  # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_NCDOT_HOTSPOTS  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
CRS = "EPSG:26918"
OUT_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "4-mgmt-forcing" / "ncdot_hotspots"
SRC = OUT_DIR / "sources"
PDF_2021 = SRC / "NCDOT_2021_NC12_Hot_Spots.pdf"
PDF_2015 = SRC / "NCDOT_2015_R-4070B_Buxton_to_Avon_feasibility_study.pdf"

# The chosen ranges (GIS, inclusive) and the NCDOT number of each, from the site config, south to north
HOTSPOTS = {name: dict(no=no, lo=lo, hi=hi)
            for no, (name, lo, hi) in sorted(HATTERAS_NCDOT_HOTSPOTS.items(), key=lambda kv: kv[1][1])}
FT = 0.3048
MILE = 1609.344

# 2021 report, PDF page 2 (overview map) at 300 dpi: two landmarks on the island and the oval centres
P2_PAGE, P2_DPI = 1, 300
P2_NORTH_TIP_PX = (1967, 229)        # the north tip of Hatteras Island at Oregon Inlet
P2_CAPE_POINT_PX = (2043, 2444)      # Cape Point, the southernmost land east of x = 1950 px
P2_OVAL_CENTRES_PX = {1: (2044, 336), 2: (2123, 558), 3: (2184, 839), 4: (2134, 1820),
                      5: (2080, 2207), 6: (1629, 2463), 7: (1138, 2677)}
P2_OVAL_NAMES = {1: "Canal Zone", 2: "Pea Island Visitor Center", 3: "Rodanthe S-Curves", 4: "Avon",
                 5: "Buxton", 6: "Hatteras Village-Frisco", 7: "Ocracoke"}

# 2021 report, PDF page 18 (Long-Term Alternative 5, Buxton beach nourishment) at 300 dpi, cropped to the map
B_PAGE, B_DPI, B_CROP = 17, 300, (300, 440, 2220, 2960)
B_ROAD_A, B_ROAD_B, B_ROAD_HALF = (819, 2255), (1031, 76), 31   # corridor around NC-12's north-south run
B_EW_BOX = (567, 2249, 781, 2283)                                # NC-12's east-west run into Buxton
B_EW_REF_X = 706
B_BAR_BOX = (1361, 2211, 1764, 2262)                             # the 0-1,250-2,500 ft scale bar
B_BLUE_MAX = (1260, 2394)                                        # the nourishment line lies left of and above this

# 2021 report, PDF page 10 (Jug Handle Bridge construction plan): vector, in PDF points
R_PAGE = 9
R_GREY = (0.75, 0.76, 0.77)          # existing pavement to be removed (and the work-area hatch)
R_GREEN = (0.41, 0.74, 0.27)         # NC-12 that stays
R_RED = (0.93, 0.12, 0.14)           # the bridge's bents
R_YELLOW = (0.97, 0.92, 0.08)        # dimension lines
R_OLD_A, R_OLD_B, R_OLD_HALF = (288, 130), (370, 592), 12        # corridor of the removed pavement
R_DIMS_FT = (1124, 1324, 1372)       # the three horizontal dimension labels, top to bottom
R_ROUNDABOUT_PT, R_BRIDGE_END_PT = (351, 751), (288, 128)

# 2015 study, PDF page 11 (Figure 1, vicinity map) at 200 dpi
F_PAGE, F_DPI = 10, 200
F_ROAD_BOX = (1089, 319, 1472, 1606)                            # NC-12 from the Buxton bend into Avon
F_SHIELDS = ((1367, 349), (1235, 1496))
F_SHIELD_R = 44
F_BAR_BOX = (957, 1914, 1430, 1958)                             # the 0-2 mile scale bar
F_POINTS = {"R-4070B Begin Project (Old Lighthouse Road)": (1155, 1584),
            "R-4070B End Project (Nino Road)": (1351, 402),
            "Park Drive (south edge of Avon)": (1293, 743)}

OVAL_FILL, OVAL_EDGE = "#9ecae1", fs.C_1997
RANGE_C = "#e66101"                  # orange: the domain_framework figures' hot-spot colour
SOURCE_C = {"nourish": "#2166ac", "removed": fs.C["ACCENT"], "bridge": "0.35"}   # bridge grey, clear of the orange
# -----------------------------------------------------------------------------


# Layers

def layers():
    import geopandas as gpd
    d = gpd.read_file(DOMAIN_BOXES).to_crs(CRS)
    d["ID"] = d.domain_id.astype(int)
    d = d[(d.ID >= 1) & (d.ID <= 90)].sort_values("ID").reset_index(drop=True)
    outline = gpd.read_file(ml.ISLAND_OUTLINE).to_crs(CRS)
    road = gpd.read_file(tv.road_line_file(2008)).to_crs(CRS)
    xy = np.vstack([np.asarray(g.coords)[:, :2] for G in road.geometry for g in getattr(G, "geoms", [G])])
    dense = np.vstack([np.linspace(xy[i], xy[i + 1], max(2, int(np.hypot(*(xy[i + 1] - xy[i])) // 10) + 1))
                       for i in range(len(xy) - 1)])
    return d, outline, road, xy, dense


# Fractional GIS of a northing (centroid interpolation) and the box that holds it
def gis_of(d, n):
    cy = d.geometry.centroid.y.values
    b = d.geometry.bounds
    hit = d[(b.miny <= n) & (b.maxy > n)].ID.values
    return float(np.interp(n, cy, d.ID.values)), (int(hit[0]) if len(hit) else None)


def render(pdf, page, dpi):
    import pymupdf
    pix = pymupdf.open(pdf)[page].get_pixmap(dpi=dpi)
    return np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, pix.n)[..., :3].astype(int)


def rigid(X, Y):
    mx, my = X.mean(0), Y.mean(0)
    U, _, Vt = np.linalg.svd((Y - my).T @ (X - mx))
    D = np.eye(2)
    if np.linalg.det(U @ Vt) < 0:
        D[1, 1] = -1
    R = U @ D @ Vt
    return R, my - R @ mx


# A fixed-scale ICP of page points (already in metres, y up) onto a road, from a grid of starts
def icp(src, road_pts, starts, rot0=(-3, 0, 3), keep_pct=90, iters=60, max_rot=8.0):
    from scipy.spatial import cKDTree
    tree = cKDTree(road_pts)
    best = None
    for c0 in starts:
        for a0 in rot0:
            th = np.radians(a0)
            R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
            t = np.asarray(c0) - R @ src.mean(0)
            for _ in range(iters):
                dd, ii = tree.query(src @ R.T + t)
                k = dd < np.percentile(dd, keep_pct)
                R, t = rigid(src[k], road_pts[ii[k]])
            rot = np.degrees(np.arctan2(R[1, 0], R[0, 0]))
            dd, _ = tree.query(src @ R.T + t)
            if abs(rot) <= max_rot and (best is None or np.median(dd) < best[0]):
                best = (float(np.median(dd)), float(np.percentile(dd, 90)), rot, R, t)
    return best


# Georeferencing

# PDF p. 2: a two-point similarity (north tip, Cape Point) and the seven ovals as polygons
def georef_overview(outline):
    from scipy import ndimage
    from shapely.geometry import MultiPoint
    xy = np.vstack([np.asarray(r.coords)[:, :2] for G in outline.geometry for g in getattr(G, "geoms", [G])
                    for r in [g.exterior]])
    north = xy[np.argmax(xy[:, 1])]
    m = xy[:, 0] > 445000
    cape = xy[m][np.argmin(xy[m][:, 1])]
    P = lambda x, y: complex(x, -y)  # noqa: E731
    p1, p2 = P(*P2_NORTH_TIP_PX), P(*P2_CAPE_POINT_PX)
    w1, w2 = complex(*north), complex(*cape)
    s = (w1 - w2) / (p1 - p2)
    t = w1 - s * p1
    a = render(PDF_2021, P2_PAGE, P2_DPI)
    blue = (a[..., 2] > 120) & (a[..., 0] < 60) & (a[..., 1] < 110) & (a[..., 2] - a[..., 0] > 80)
    lab, n = ndimage.label(ndimage.binary_closing(blue, iterations=3))
    ovals = {}
    for i in range(1, n + 1):
        yy, xx = np.where(lab == i)
        if len(xx) < 3000:
            continue
        c = np.array([xx.mean(), yy.mean()])
        no = min(P2_OVAL_CENTRES_PX, key=lambda k: np.hypot(*(c - P2_OVAL_CENTRES_PX[k])))
        if np.hypot(*(c - P2_OVAL_CENTRES_PX[no])) > 60:
            continue                                             # the legend's oval
        z = s * (xx - 1j * yy) + t
        ovals[no] = MultiPoint(np.c_[z.real, z.imag][::5]).convex_hull
    # a third landmark the fit did not use: the island's easternmost point, at the Rodanthe bulge
    east_px = (2221, 1004)
    e = s * P(*east_px) + t
    east = xy[np.argmax(xy[:, 0])]
    fit = dict(map="2021 overview (p. 2)", method="two-point similarity: Oregon Inlet tip, Cape Point",
               scale=abs(s), scale_unit="m/px at 300 dpi", rotation_deg=float(np.degrees(np.angle(s))),
               median_m=np.nan, p90_m=np.nan,
               check=f"easternmost land off by {np.hypot(e.real - east[0], e.imag - east[1]):.0f} m")
    return ovals, fit


# PDF p. 18: the Buxton nourishment map, scale from its bar, rotation from NC-12's north-south run
def georef_buxton(xy_road):
    x0, y0, x1, y1 = B_CROP
    a = render(PDF_2021, B_PAGE, B_DPI)[y0:y1, x0:x1]
    # scale: the bar's 0 and 2,500 ft ticks are the first and last narrow tick clusters
    bx0, by0, bx1, by1 = B_BAR_BOX
    dk = a[by0:by1, bx0:bx1].max(2) < 90
    r = int(np.argmax(dk.sum(1)))
    cols = np.where(dk[max(r - 12, 0):r - 2].sum(0) > 4)[0]
    groups = np.split(cols, np.where(np.diff(cols) > 1)[0] + 1)
    ticks = [g.mean() for g in groups if 1 <= len(g) <= 5]
    mpp = 2500 * FT / (ticks[-1] - ticks[0])
    # NC-12's north-south run on the page and on the ground
    A, B = np.array(B_ROAD_A, float), np.array(B_ROAD_B, float)
    ys, xs = np.where(a.max(2) < 75)
    P = np.c_[xs, ys].astype(float)
    u = (B - A) / np.linalg.norm(B - A)
    nn = np.array([-u[1], u[0]])
    t_ = (P - A) @ u
    Q = P[(np.abs((P - A) @ nn) < B_ROAD_HALF) & (t_ > 50) & (t_ < np.linalg.norm(B - A) - 25)]
    c = Q.mean(0)
    dp = np.linalg.svd(Q - c)[2][0]
    dp = dp if dp[1] < 0 else -dp
    az_page = np.degrees(np.arctan2(dp[0], -dp[1]))
    seg = xy_road[(xy_road[:, 1] > 3903000) & (xy_road[:, 1] < 3909000)]
    cc = seg.mean(0)
    dw = np.linalg.svd(seg - cc)[2][0]
    dw = dw if dw[1] > 0 else -dw
    az_utm = np.degrees(np.arctan2(dw[0], dw[1]))
    th = np.radians(-(az_utm - az_page))
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    lin = lambda px: (np.c_[px[:, 0], -px[:, 1]] * mpp) @ R.T  # noqa: E731
    # translation: the page's north-south run on the road; its east-west run at the road's east-west northing
    ex0, ey0, ex1, ey1 = B_EW_BOX
    ew_y = ey0 + float(np.median(np.where(a[ey0:ey1, ex0:ex1].max(2) < 75)[0]))
    ew_n = float(np.median(xy_road[xy_road[:, 1] < 3902700][:, 1]))
    nv = np.array([dw[1], -dw[0]])
    T = np.linalg.solve(np.array([nv, [0, 1]]),
                        np.array([nv @ (cc - lin(c[None])[0]), ew_n - lin(np.array([[B_EW_REF_X, ew_y]]))[0][1]]))
    W = lambda px: lin(np.atleast_2d(np.asarray(px, float))) + T  # noqa: E731
    # the nourishment line (blue) and a check: the page's road far up the map against the real road
    bl = (a[..., 2] > 150) & (a[..., 0] < 80) & (a[..., 1] < 120) & (a[..., 2] - a[..., 1] > 60)
    bl[B_BLUE_MAX[1]:] = False
    bl[:, B_BLUE_MAX[0]:] = False
    yy, xx = np.where(bl)
    line_px = np.c_[xx, yy].astype(float)
    # keep the long line, not the Haulover Beach Parking dot: points within 15 px of its principal axis
    lc = line_px.mean(0)
    ld = np.linalg.svd(line_px - lc)[2][0]
    line_px = line_px[np.abs((line_px - lc) @ np.array([-ld[1], ld[0]])) < 15]
    top, bot = line_px[np.argmin(line_px[:, 1])], line_px[np.argmax(line_px[:, 1])]
    from scipy.spatial import cKDTree
    dd, _ = cKDTree(xy_road).query(W(Q))
    fit = dict(map="2021 Buxton Alternative 5 (p. 18)",
               method="scale bar (2,500 ft); rotation from NC-12's north-south run; NC-12 bend for position",
               scale=mpp, scale_unit="m/px at 300 dpi", rotation_deg=float(np.degrees(th)),
               median_m=float(np.median(dd)), p90_m=float(np.percentile(dd, 90)),
               check="Haulover Beach Parking falls in GIS 15")
    return dict(nourish=W(np.vstack([bot, top])), nourish_pts=W(line_px[::4])), fit


def _vector_points(page, col):
    out = []
    for x in page.get_drawings():
        if tuple(round(v, 2) for v in (x["color"] or ())) != col:
            continue
        for it in x["items"]:
            if it[0] == "l":
                a, b = it[1], it[2]
                for f in np.linspace(0, 1, 6):
                    out.append((a.x + (b.x - a.x) * f, a.y + (b.y - a.y) * f))
            elif it[0] == "c":
                out += [(it[1].x, it[1].y), (it[4].x, it[4].y)]
    return np.array(out)


def _centreline(P, step=2.0):
    out = []
    for y in np.arange(P[:, 1].min(), P[:, 1].max(), step):
        m = (P[:, 1] >= y) & (P[:, 1] < y + step)
        if m.any():
            out.append((np.median(P[m, 0]), y + step / 2))
    return np.array(out)


# PDF p. 10: the Jug Handle Bridge plan, vector; scale from its dimension lines, fit on the old NC-12 path
def georef_rodanthe(dense_road):
    import pymupdf
    page = pymupdf.open(PDF_2021)[R_PAGE]
    dims = []
    for x in page.get_drawings():
        if tuple(round(v, 2) for v in (x["color"] or ())) != R_YELLOW:
            continue
        for it in x["items"]:
            if it[0] == "l" and abs(it[1].y - it[2].y) < 0.5 and abs(it[1].x - it[2].x) > 40:
                dims.append((it[1].y, abs(it[1].x - it[2].x)))
    dims = sorted(set((round(y, 1), round(l, 2)) for y, l in dims))
    mpt = float(np.mean([ft * FT / l for (_, l), ft in zip(dims, R_DIMS_FT)]))
    green = _vector_points(page, R_GREEN)
    green = green[(green[:, 1] < 130) | ((green[:, 1] > 590) & (green[:, 1] < 735))]
    grey = _vector_points(page, R_GREY)
    A, B = np.array(R_OLD_A, float), np.array(R_OLD_B, float)
    u = (B - A) / np.linalg.norm(B - A)
    nn = np.array([-u[1], u[0]])
    t_ = (grey - A) @ u
    old = _centreline(grey[(t_ > 0) & (t_ < np.linalg.norm(B - A)) & (np.abs((grey - A) @ nn) < R_OLD_HALF)])
    path = np.vstack([_centreline(green), old])
    road = dense_road[dense_road[:, 1] > 3925000]
    src = mpt * np.c_[path[:, 0], -path[:, 1]]
    starts = [(road[:, 0].mean(), n0) for n0 in np.arange(3928000, 3956000, 200)]
    med, p90, rot, R, t = icp(src, road, starts, rot0=(-4, 0, 4), keep_pct=95)
    W = lambda P: mpt * np.c_[np.atleast_2d(P)[:, 0], -np.atleast_2d(P)[:, 1]] @ R.T + t  # noqa: E731
    bridge = _vector_points(page, R_RED)
    fit = dict(map="2021 Jug Handle Bridge plan (p. 10)",
               method=f"dimension lines ({', '.join(f'{v} ft' for v in R_DIMS_FT)}) for scale; "
                      "ICP of the old NC-12 path onto the 2008 line",
               scale=mpt, scale_unit="m/pt", rotation_deg=rot, median_m=med, p90_m=p90,
               check="bridge landings at the Rodanthe roundabout and on Pea Island")
    return dict(removed=W(old), bridge=W(bridge), roundabout=W(R_ROUNDABOUT_PT)[0],
                bridge_end=W(R_BRIDGE_END_PT)[0]), fit


# 2015 study Figure 1: scale bar, ICP of NC-12 from the Buxton bend into Avon
def georef_r4070b(dense_road):
    from scipy import ndimage
    a = render(PDF_2015, F_PAGE, F_DPI)
    bx0, by0, bx1, by1 = F_BAR_BOX
    dk = a[by0:by1, bx0:bx1].max(2) < 80
    r = int(np.argmax(dk.sum(1)))
    cols = np.where(dk[r])[0]
    runs = np.split(cols, np.where(np.diff(cols) > 1)[0] + 1)
    bar = max(runs, key=len)                                     # the bar's outline, not the "2 Miles" type beside it
    mpp = 2 * MILE / (bar.max() - bar.min())
    x0, y0, x1, y1 = F_ROAD_BOX
    blk = np.zeros(a.shape[:2], bool)
    blk[y0:y1, x0:x1] = a[y0:y1, x0:x1].max(2) < 60
    blk = ndimage.binary_opening(blk, iterations=2)
    ys, xs = np.where(blk)
    P = np.c_[xs, ys].astype(float)[::2]
    for sx, sy in F_SHIELDS:
        P = P[np.hypot(P[:, 0] - sx, P[:, 1] - sy) > F_SHIELD_R]
    road = dense_road[dense_road[:, 1] < 3916000]
    src = mpp * np.c_[P[:, 0], -P[:, 1]]
    starts = [(453300 + dx, 3907500 + dy) for dx in np.arange(-1500, 1501, 500) for dy in np.arange(-3000, 3001, 500)]
    med, p90, rot, R, t = icp(src, road, starts, rot0=(-2, 0, 2), keep_pct=85)
    W = lambda px: mpp * np.c_[np.atleast_2d(px)[:, 0], -np.atleast_2d(px)[:, 1]] @ R.T + t  # noqa: E731
    pts = {k: W(v)[0] for k, v in F_POINTS.items()}
    fit = dict(map="2015 R-4070B Figure 1 (p. 11)", method="scale bar (2 mi); ICP of NC-12 onto the 2008 line",
               scale=mpp, scale_unit="m/px at 200 dpi", rotation_deg=rot, median_m=med, p90_m=p90,
               check="Park Drive at the start of the Avon community zone (GIS 21)")
    return pts, fit


# Outputs

def tables(d, ovals, bux, rod, r4070, fits):
    rows = []
    for name, h in HOTSPOTS.items():
        ov = ovals[h["no"]]
        o_lo, o_hi = gis_of(d, ov.bounds[1])[0], gis_of(d, ov.bounds[3])[0]
        o_c = gis_of(d, ov.centroid.y)[0]
        if name == "Buxton":
            ev = (f"text: 'at the northern end of Buxton at the Cape Hatteras Motel'; Alternative 5 nourishment "
                  f"line GIS {gis_of(d, bux['nourish'][0][1])[0]:.1f}-{gis_of(d, bux['nourish'][1][1])[0]:.1f} "
                  "(a design length); range ends at the oval's north end")
        elif name == "Avon":
            ev = (f"oval only (no detailed map in the 2021 report); the 2015 R-4070B study ends its Buxton-Avon "
                  f"project at Nino Road, GIS {gis_of(d, r4070['R-4070B End Project (Nino Road)'][1])[0]:.1f}")
        else:
            lo, hi = rod["removed"][:, 1].min(), rod["removed"][:, 1].max()
            ev = (f"pavement removed for the Jug Handle Bridge, GIS {gis_of(d, lo)[0]:.1f}-{gis_of(d, hi)[0]:.1f}; "
                  f"bridge lands at GIS {gis_of(d, rod['roundabout'][1])[0]:.1f} and "
                  f"{gis_of(d, rod['bridge_end'][1])[0]:.1f}")
        rows.append(dict(hotspot=name, ncdot_no=h["no"], gis_lo=h["lo"], gis_hi=h["hi"],
                         oval_gis_lo=round(o_lo, 1), oval_gis_hi=round(o_hi, 1), oval_gis_centre=round(o_c, 1),
                         evidence=ev))
    rng = pd.DataFrame(rows)
    rng.to_csv(OUT_DIR / "ncdot_hotspot_ranges.csv", index=False)
    pd.DataFrame(fits).round(3).to_csv(OUT_DIR / "georeference_fits.csv", index=False)
    return rng


def geojson(d, ovals, bux, rod, r4070):
    from shapely.geometry import LineString, Point, mapping
    feats = []

    def add(geom, **props):
        feats.append({"type": "Feature", "properties": props, "geometry": mapping(geom)})
    for no, ov in sorted(ovals.items()):
        add(ov, kind="NCDOT 2021 oval (p. 2)", ncdot_no=no, name=P2_OVAL_NAMES[no])
    for name, h in HOTSPOTS.items():
        g = d[(d.ID >= h["lo"]) & (d.ID <= h["hi"])].geometry.union_all()
        add(g, kind="chosen hot-spot range", ncdot_no=h["no"], name=name, gis_lo=h["lo"], gis_hi=h["hi"])
    add(LineString(bux["nourish"]), kind="2021 Buxton Alternative 5 nourishment line (p. 18)", name="Buxton")
    o = rod["removed"][np.argsort(rod["removed"][:, 1])]
    add(LineString(o), kind="2021 pavement removed for the Jug Handle Bridge (p. 10)", name="Rodanthe S-Curves")
    add(Point(rod["roundabout"]), kind="Jug Handle Bridge south landing (Rodanthe roundabout)")
    add(Point(rod["bridge_end"]), kind="Jug Handle Bridge north landing")
    for k, v in r4070.items():
        add(Point(v), kind="2015 R-4070B Figure 1 point", name=k)
    fc = {"type": "FeatureCollection", "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::26918"}},
          "features": feats}
    (OUT_DIR / "ncdot_hotspots_georef.geojson").write_text(json.dumps(fc))


# Figures

def _base(ax, d, outline, road, win, label_every=None, label_side="left"):
    import geopandas as gpd
    from shapely.geometry import box
    x0, x1, y0, y1 = win
    ax.set_facecolor(fs.WATER_FILL if hasattr(fs, "WATER_FILL") else "#e9eff4")
    clip = box(x0 - 2000, y0 - 2000, x1 + 2000, y1 + 2000)
    land = outline.geometry.intersection(clip)
    gpd.GeoSeries(land[~land.is_empty], crs=outline.crs).plot(ax=ax, facecolor="#ede9df", edgecolor="0.55",
                                                              lw=0.35, zorder=1)
    d.plot(ax=ax, facecolor="none", edgecolor="0.45", lw=0.3, zorder=3)
    road.clip(clip).plot(ax=ax, color=fs.INK, lw=0.7, zorder=5)
    if label_every:
        for _, r in d.iterrows():
            b = r.geometry.bounds
            if (r.ID == 1 or r.ID % label_every == 0) and y0 < (b[1] + b[3]) / 2 < y1:
                ax.text(b[0] - 250, (b[1] + b[3]) / 2, str(r.ID), ha="right", va="center", clip_on=True,
                        **fs.MAP_TEXT_DARK)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fs.spines_for_image(ax)


def _range_outline(ax, d, lo, hi, buf=110, lw=1.4):
    g = d[(d.ID >= lo) & (d.ID <= hi)].geometry.union_all().buffer(buf, join_style=2)
    for p in getattr(g, "geoms", [g]):
        x, y = p.exterior.xy
        ax.plot(x, y, color=RANGE_C, lw=lw, ls=(0, (4, 2)), zorder=7, solid_capstyle="butt")
    return g


def _oval(ax, ov, alpha=0.55):
    import geopandas as gpd
    gpd.GeoSeries([ov]).plot(ax=ax, facecolor=OVAL_FILL, alpha=alpha, edgecolor="none", zorder=2)
    gpd.GeoSeries([ov]).boundary.plot(ax=ax, color=OVAL_EDGE, lw=0.7, zorder=2.5)


# A zoom: `height_m` tall around the range (and any `include` northings), its left edge just clear of the
# sound-side domain numbers so the beach and some ocean stay in frame
def _zoom_window(d, h, aspect, height_m, include=(), label_room_m=650):
    sel = d[(d.ID >= h["lo"]) & (d.ID <= h["hi"])]
    x0, y0, _, y1 = sel.total_bounds
    lo, hi = min([y0, *include]), max([y1, *include])
    yc = (lo + hi) / 2
    near = d[(d.geometry.bounds.maxy > yc - height_m / 2) & (d.geometry.bounds.miny < yc + height_m / 2)]
    left = near.total_bounds[0] - label_room_m
    return (left, left + height_m * aspect, yc - height_m / 2, yc + height_m / 2)


def _legend_handles(extra=()):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    return [Patch(facecolor=OVAL_FILL, edgecolor=OVAL_EDGE, lw=0.7, label="NCDOT (2021) symbol"),
            Line2D([], [], color=RANGE_C, lw=1.4, ls=(0, (4, 2)), label="hot-spot range"),
            Patch(facecolor="none", edgecolor="0.45", lw=0.5, label="model domain"),
            Line2D([], [], color=fs.INK, lw=0.8, label="NC-12, 2008"), *extra]


def fig_map(d, outline, road, ovals):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=fs.figsize("double", height=5.6))
    fw, fh = fig.get_size_inches()
    bottom, top = 0.085, 0.975
    ph = (top - bottom) * fh                                       # panel height, inches
    bx0, by0, bx1, by1 = d.total_bounds
    win_a = (bx0 - 3000, bx1 + 6000, by0 - 2000, by1 + 2000)
    wa = ph * (win_a[1] - win_a[0]) / (win_a[3] - win_a[2])
    ax_a = fig.add_axes([0.035, bottom, wa / fw, top - bottom])
    _base(ax_a, d, outline, road, win_a, label_every=10)
    for name, h in HOTSPOTS.items():
        _oval(ax_a, ovals[h["no"]])
        g = _range_outline(ax_a, d, h["lo"], h["hi"], buf=150, lw=1.1)
        b = g.bounds
        short = {"Rodanthe S-Curves": "S-Curves"}.get(name, name)
        fs.place_label(ax_a, b[2] + 700, (b[1] + b[3]) / 2, f"{h['no']}\n{short}", fs.MAP_TEXT_DARK, ha="left",
                       fontstyle="normal")
    # scale bar and north arrow in the open sound, clear of the domains and their numbers
    fs.scale_bar_km(ax_a, length_m=10_000, x=0.08, y=0.62, text_kw=fs.MAP_TEXT_DARK)
    fs.north_dart(ax_a, fs.axes_point(ax_a, 0.22, 0.80), north=(0, 1), arrow_m=2600, text_kw=fs.MAP_TEXT_DARK)
    fs.letter_at(ax_a, 0, 0.03, 0.985)
    # three zooms, south to north
    x_left = 0.035 + wa / fw + 0.03
    gap = 0.012
    zw = (0.995 - x_left - 2 * gap) / 3
    aspect = (zw * fw) / ph
    for i, (name, h) in enumerate(HOTSPOTS.items()):
        ax = fig.add_axes([x_left + i * (zw + gap), bottom, zw, top - bottom])
        win = _zoom_window(d, h, aspect, height_m=4400 / aspect)     # 4.4 km across: sound, domains, beach, the symbol
        _base(ax, d, outline, road, win, label_every=2)
        _oval(ax, ovals[h["no"]])
        _range_outline(ax, d, h["lo"], h["hi"])
        fs.scale_bar_km(ax, length_m=1000, segments=1, x=0.55, y=0.05, text_kw=fs.MAP_TEXT_DARK)
        fs.letter_at(ax, i + 1, 0.04, 0.985)
        ax.text(0.96, 0.985, f"{h['no']} {name}\nGIS {h['lo']}–{h['hi']}", transform=ax.transAxes, ha="right",
                va="top", fontsize=8, color=fs.INK, zorder=20,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
        ax_a.add_patch(plt.Rectangle((win[0], win[2]), win[1] - win[0], win[3] - win[2], facecolor="none",
                                     edgecolor=fs.INK, lw=0.6, zorder=9))
    fig.legend(handles=_legend_handles(), loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0),
               **{**fs.MAP_LEGEND, "framealpha": 0})
    out = fs.save(fig, OUT_DIR / "ncdot_hotspots_map.png", vector=True, dpi=300)
    fs.record_caption(out[0],
        "NCDOT's NC 12 hot spots on the model domains. (a) The reach, north up, every tenth domain numbered; the "
        "three hot spots of NCDOT (2021) inside it (3 Rodanthe 'S' Curves, 4 Avon, 5 Buxton; 1, 2, 6 and 7 lie "
        "outside) with their map symbols georeferenced from the report's overview map (blue) and the ranges taken "
        "from them (orange dashed; GIS 8-14, 28-35, 84-88). The boxes mark panels (b)-(d). (b-d) Each hot spot, every "
        "second domain numbered. NCDOT draws every hot spot with the same 3.8 km symbol, so the symbol places a hot "
        "spot but does not measure it; a range follows the symbol and is trimmed where NCDOT's text or detailed "
        "maps end the hot spot sooner (Buxton's south end at the Cape Hatteras Motel, GIS 8; Rodanthe at the "
        "pavement removed for the Jug Handle Bridge, GIS 84-88). Georeferencing error of the symbols is about "
        "1 km alongshore. NC-12 as digitised on 2008 imagery.")
    plt.close(fig)
    return out[0]


def fig_georeference(d, outline, road, ovals, bux, rod, r4070):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    fig = plt.figure(figsize=fs.figsize("double", height=5.6))
    fw, fh = fig.get_size_inches()
    bottom, top, gap = 0.12, 0.975, 0.012
    ph = (top - bottom) * fh
    zw = (0.99 - 0.01 - 2 * gap) / 3
    aspect = zw * fw / ph
    for i, (name, h) in enumerate(HOTSPOTS.items()):
        ax = fig.add_axes([0.01 + i * (zw + gap), bottom, zw, top - bottom])
        inc = [v[1] for k, v in r4070.items() if not k.startswith("R-4070B Begin")] if name == "Avon" else []
        win = _zoom_window(d, h, aspect, height_m=10500, include=inc)
        _base(ax, d, outline, road, win, label_every=2)
        _oval(ax, ovals[h["no"]], alpha=0.35)
        _range_outline(ax, d, h["lo"], h["hi"])
        if name == "Buxton":
            p = bux["nourish_pts"]
            ax.plot(p[:, 0], p[:, 1], ".", ms=1.2, color=SOURCE_C["nourish"], zorder=6)
            for (x, y), lab in zip(bux["nourish"], ("Begin", "End")):
                ax.plot(x, y, "o", ms=3.5, mfc="white", mec=SOURCE_C["nourish"], mew=1.0, zorder=8)
                ax.text(x + 350, y, lab, ha="left", va="center", **fs.MAP_TEXT_DARK)
        elif name == "Avon":
            for k, v in r4070.items():
                if k.startswith("R-4070B Begin"):
                    continue
                ax.plot(*v, "D", ms=3.8, mfc="white", mec=fs.C["REF"], mew=1.1, zorder=8)
                lab = "R-4070B end\n(Nino Rd)" if "End" in k else "Park Dr"
                ax.text(v[0] + 350, v[1], lab, ha="left", va="center", **fs.MAP_TEXT_DARK)
        else:
            p = rod["bridge"]
            ax.plot(p[:, 0], p[:, 1], ".", ms=0.8, color=SOURCE_C["bridge"], zorder=6)
            p = rod["removed"]
            ax.plot(p[:, 0], p[:, 1], ".", ms=1.4, color=SOURCE_C["removed"], zorder=6)
        fs.scale_bar_km(ax, length_m=1000, segments=1, x=0.55, y=0.05, text_kw=fs.MAP_TEXT_DARK)
        fs.letter_at(ax, i, 0.04, 0.985)
        ax.text(0.96, 0.985, f"{h['no']} {name}\nGIS {h['lo']}–{h['hi']}", transform=ax.transAxes, ha="right",
                va="top", fontsize=8, color=fs.INK, zorder=20,
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", boxstyle="square,pad=0.2"))
    extra = [Line2D([], [], color=SOURCE_C["nourish"], lw=1.6, label="Alt. 5 nourishment (2021)"),
             Line2D([], [], ls="none", marker="D", ms=3.8, mfc="white", mec=fs.C["REF"], mew=1.1,
                    label="R-4070B points (2015)"),
             Line2D([], [], color=SOURCE_C["removed"], lw=1.6, label="NC-12 removed (2021)"),
             Line2D([], [], color=SOURCE_C["bridge"], lw=1.6, label="Jug Handle Bridge (2021)")]
    fig.legend(handles=_legend_handles(extra), loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0),
               **{**fs.MAP_LEGEND, "framealpha": 0})
    out = fs.save(fig, OUT_DIR / "ncdot_hotspots_georeference.png", vector=True, dpi=300)
    rf = rod["removed"][:, 1]
    fs.record_caption(out[0],
        "The evidence behind each hot-spot range (orange dashed), every NCDOT map georeferenced onto the 2008 NC-12 "
        "line. (a) Buxton: the beach nourishment line of NCDOT (2021) Long-Term Alternative 5 (blue), Begin at the "
        "NC-12 bend at Old Lighthouse Road, the hot spot's south end in the report's text; it runs on north of the "
        "hot-spot symbol as a design length. (b) Avon: the 2021 report has no detailed map, so the range is the "
        "symbol; the 2015 Buxton-to-Avon feasibility study (R-4070B) ends its one Buxton-Avon hot-spot project at "
        "Nino Road in south Avon (green diamonds), south of the symbol. (c) Rodanthe 'S' Curves: the pavement "
        f"removed when the Jug Handle Bridge opened (purple, GIS {gis_of(d, rf.min())[0]:.1f}-"
        f"{gis_of(d, rf.max())[0]:.1f}) and the bridge (grey). Scale from each map's own bar or dimension lines; "
        "fit statistics in georeference_fits.csv.")
    plt.close(fig)
    return out[0]


def main():
    import matplotlib
    matplotlib.use("Agg")
    fs.apply_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    d, outline, road, xy, dense = layers()
    ovals, f1 = georef_overview(outline)
    bux, f2 = georef_buxton(xy)
    rod, f3 = georef_rodanthe(dense)
    r4070, f4 = georef_r4070b(dense)
    fits = [f1, f2, f3, f4]
    rng = tables(d, ovals, bux, rod, r4070, fits)
    geojson(d, ovals, bux, rod, r4070)
    print(pd.DataFrame(fits)[["map", "scale", "scale_unit", "rotation_deg", "median_m", "p90_m", "check"]]
          .round(2).to_string(index=False))
    print(rng.drop(columns="evidence").to_string(index=False))
    for k, v in r4070.items():
        print(f"{k}: GIS {gis_of(d, v[1])[0]:.1f}")
    print("wrote", fig_map(d, outline, road, ovals))
    print("wrote", fig_georeference(d, outline, road, ovals, bux, rod, r4070))


if __name__ == "__main__":
    main()
