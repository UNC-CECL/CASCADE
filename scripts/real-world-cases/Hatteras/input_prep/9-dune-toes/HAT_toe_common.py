"""
Shared by the dune toe scripts of input_prep/9-dune-toes (2026-11-05): the settings, the 1996 lidar and the 100 m
transects, a profile along a transect, the model's dune search windows on it, the figure palette, and the files of
your crest picks and toe edits (read by every method's run). Not run on its own.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-11-05
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
sys.path.insert(0, str(_HERE.parents[1] / "1-barrier3d-domains" / "3-dune-toes"))
import HAT_dune_toes as T  # noqa: E402

S = T.S

# --- CONFIG ------------------------------------------------------------------
TRANSECTS = S.REPO / "data" / "hatteras_init" / "2-brie-offset" / "transects" / "transects_100m.geojson"
TRANSECT_ID = "Transects_100m.LineID"
YEAR = int(os.environ.get("HAT_TOE_YEAR", 1996))   # the transect survey: 1996, or 2009 for the 2009 start
START = {1996: "1996-start", 2009: "2009-start"}[YEAR]   # whose dune windows and tables
SURVEY_START = "1996-start"      # S.STARTS and the vertical offsets list both surveys under this start
STEP_M = 1.0                    # profile spacing
MHW = T.MHW_NAVD88_M             # 0.36 m NAVD88
SHORE_RUN_M = 10.0               # the shoreline: the first point above MHW with this much dry profile behind it
CREST_SEARCH_M = 120.0           # the crest is looked for this far landward of the shoreline
FOREDUNE_DROP_M = T.FOREDUNE_DROP_M        # 1.0: the front dune is the first peak with this deep a dip behind it
PEAK_HALF_M = 5.0                # a peak is the highest point within this distance each side
TOE_MIN_M_NAVD88 = T.TOE_MIN_M_NAVD88      # 1.0
TOE_BELOW_CREST_M = T.TOE_BELOW_CREST_M    # 0.5
BERM_M = S.X.BERM_ELEV_NAVD_M    # 1.7 m NAVD88: the model's berm; the dune stands on it, so no toe below it
BEND_SHARE = 0.6                 # a second bend: a point at least this share as far below the line as the toe
BEND_M = 1.0                     # ... whose toe would be at least this far from the one picked
MIN_KEPT = 3                     # transects kept for a domain median
EXAMPLE_GIS = 22
OUT = S.REPO / "data" / "hatteras_init" / "9-dune-toes" / f"1-transects-{YEAR}"   # D:; subfolders by topic:
TOES = OUT / "1-toes"                       # transect and domain toes, the viewer/picker JSON
PICKS = OUT / "2-picks"                     # your crest picks and toe edits, the picker views
METHODS = OUT / "3-methods"                 # one folder per method: its tables and its method figure
M_LINE = METHODS / "1-line"                 # ours (the toes themselves are in 1-toes/)
M_AUTOMORPH = METHODS / "2-automorph"       # Itzkin et al. (2020)
M_CURVATURE = METHODS / "3-curvature"       # Stockdon et al. (2007)
M_PYBEACH = METHODS / "4-pybeach"           # Beuzen (2019)
M_COMPARE = METHODS / "5-comparison"        # all of them side by side
CHECKS = OUT / "4-checks"                   # crest-choice/, dip/
EXAMPLES = (119, 326)                       # the example transects of the method figures: a single foredune (GIS 22),
                                            # a dune field (GIS 63)
# -----------------------------------------------------------------------------


def survey():
    import rasterio
    kind, path, *_ = S.STARTS[SURVEY_START][1][YEAR]
    src = rasterio.open(path)
    scale = S.US_FT if S.stated_unit(src) == "ft" else 1.0
    return src, scale, S.SURVEY_VERTICAL_OFFSET_M.get(SURVEY_START, {}).get(YEAR, 0.0)


def transects():
    import geopandas as gpd
    t = gpd.read_file(TRANSECTS)[[TRANSECT_ID, "geometry"]].rename(columns={TRANSECT_ID: "ObjectID"})
    d =gpd.read_file(S.GIS_ROOT / "domains.geojson").to_crs(t.crs)
    j = gpd.sjoin(t, d[["domain_id", "geometry"]], predicate="intersects", how="inner")
    j = j.rename(columns={"domain_id": "gis"}).sort_values(["gis", "ObjectID"])
    return j[["ObjectID", "gis", "geometry"]].reset_index(drop=True), d


# Sample points along one transect: s = distance from its ocean end (the first vertex, offshore to the east)
def line_points(geom):
    c = np.asarray(geom.coords, float)
    u = c[-1] - c[0]
    u /= np.hypot(*u)                                            # unit vector, ocean -> sound
    s = np.arange(0.0, geom.length + STEP_M / 2, STEP_M)
    return s, c[0][None, :] + s[:, None] * u[None, :], c[0], u


def sample(src, scale, offset, xy, crs):
    from rasterio.warp import transform
    from rasterio.windows import from_bounds
    from scipy.ndimage import map_coordinates
    xs, ys = transform(crs, src.crs, xy[:, 0], xy[:, 1])
    xs, ys = np.asarray(xs), np.asarray(ys)
    pad = 4 * abs(src.res[0])
    win = from_bounds(xs.min() - pad, ys.min() - pad, xs.max() + pad, ys.max() + pad,
                      transform=src.transform).round_offsets().round_lengths()
    a = S.clean(src.read(1, window=win, boundless=True, fill_value=np.nan), src.nodata)
    inv = ~src.window_transform(win)
    cols, rows = inv * (xs, ys)
    z = map_coordinates(a, [rows - 0.5, cols - 0.5], order=1, mode="constant", cval=np.nan)
    return z * scale + offset


def first_dune(seg, drop=None):
    """The front dune's crest (index into seg, from the shoreline): the most seaward peak with a dip of at least
    FOREDUNE_DROP_M behind it before the profile climbs above it again; the highest point where none has one."""
    z = np.where(np.isfinite(seg), seg, -np.inf)
    w = int(PEAK_HALF_M / STEP_M)
    for i in range(w, len(z) - 1):
        if z[i] == -np.inf or z[i] < z[max(0, i - w):i + w + 1].max():
            continue                                             # not a local peak
        behind = z[i + 1:]
        higher = np.nonzero(behind > z[i])[0]
        stop = higher[0] if len(higher) else len(behind)
        if stop and z[i] - behind[:stop].min() >= (FOREDUNE_DROP_M if drop is None else drop):
            return i
    return int(np.nanargmax(seg))


class Windows:
    """The model's dune search windows (HAT_survey_dune_crests: the start's per-domain windows, moved per cell for
    the survey) as distances along a transect. The 1 m clip of each domain is the model array's grid: oriented
    ocean first (columns flipped) and south first (ALONGSHORE_FLIP), each 10 m cell's profile shifted by its shear
    and the water trim c0, so window column i is clip columns 10 * (i + c0 + shear) onward, counted from the east."""

    def __init__(self):
        import json
        self.tv = S.tv
        self.w = json.loads(S.tv.picks_dir(START).joinpath(
            f"HAT_dune_search_windows_{S.tv.topo_dirs(START)[2]}.json").read_text(encoding="utf-8"))
        self.clips = S.product(S.STARTS[START][0]).gapfill_1m
        self.cache = {}

    def domain(self, g):
        if g not in self.cache:
            t, (h, w), _ = S.clip_grid(self.clips / f"clip_domain_{g}_filled.tif")
            shear, c0, i0, i1 = S.frame(START, g, self.w)
            a, b = S.survey_windows(START, YEAR, g, i0, i1, n=h // S.BLOCK)
            self.cache[g] = (t, h, w, np.asarray(shear, int), int(c0), np.asarray(a, int), np.asarray(b, int))
        return self.cache[g]

    def along(self, g, p0, u, y):
        """(cell, s_seaward, s_landward) of the window crossed at northing y, s from the transect's ocean end."""
        t, h, w, shear, c0, a, b = self.domain(g)
        bottom = t.f - h                                         # oriented row 0 is the south edge
        cell = int((y - bottom) // S.BLOCK)
        if not 0 <= cell < len(a):
            return None
        j0 = S.BLOCK * (a[cell] + c0 + shear[cell])              # oriented 1 m columns, ocean first
        j1 = S.BLOCK * (b[cell] + c0 + shear[cell]) - 1
        e_sea, e_land = t.c + w - 0.5 - j0, t.c + w - 0.5 - j1   # eastings
        s = lambda e: float((e - p0[0]) * u[0] + (y - p0[1]) * u[1])
        return cell, s(e_sea), s(e_land)


# --- figures -----------------------------------------------------------------
# Journal style: no titles in the images (CAPTIONS.md under supporting/ carries the text), one size for every
# single panel, muted colours, and the same panels together as one lettered figure (dune_toe_method_1996.png)
PROFILE, FILL, WATER, CREST, TOE, CHORD = "#1A1A1A", "#E3E3E3", "0.30", "#8C510A", "#08519C", "0.40"
CURRENT, DROPPED = "0.45", "#A8A8A8"
WINDOW = "#1F8A84"                               # the model's dune search window method
VIEW_M = (-40.0, CREST_SEARCH_M + 30.0)          # cross-shore view on the profile panels


def _ends(ax, fs, left="Ocean", right="Land"):
    for x, ha, txt in ((0.012, "left", left), (0.988, "right", right)):
        ax.text(x, 0.025, txt, transform=ax.transAxes, ha=ha, va="bottom", fontsize=7, style="italic",
                color=fs.INK_MUTED)


def _legend(ax, **kw):
    kw = dict(loc="upper left", fontsize=7, frameon=False, handlelength=1.6, borderaxespad=0.3) | kw
    ax.legend(**kw)


CREST_FILE = PICKS / f"crest_picks_{YEAR}.json"   # your crests: {transect: {status, crest_x, crest_z, updated}};


def crest_picks():
    import json
    return json.loads(CREST_FILE.read_text(encoding="utf-8")) if CREST_FILE.exists() else {}


QC_FILE = PICKS / f"toe_qc_{YEAR}.json"   # the picker's checks: {transect: {status, toe_x, toe_z, updated}}


def qc_edits():
    import json
    return json.loads(QC_FILE.read_text(encoding="utf-8")) if QC_FILE.exists() else {}


def apply_qc(t):
    """The checked toe per transect (qc_ columns): your moved toe, nothing where rejected, else the automatic pick
    where it passed the screen; the status is 'unchecked' where the picker has not been."""
    e = qc_edits()
    st = t.transect.astype(str).map(lambda i: e.get(i, {}).get("status", "unchecked"))
    z = t.transect.astype(str).map(lambda i: e.get(i, {}).get("toe_z"))
    t["qc_status"] = st
    t["qc_toe_m_navd88"] = np.where(st == "moved", z.astype(float), np.where(t.kept, t.toe_m_navd88, np.nan))
    t["qc_kept"] = (st != "rejected") & np.isfinite(t.qc_toe_m_navd88)
    t.loc[st == "rejected", "qc_toe_m_navd88"] = np.nan
    return t
