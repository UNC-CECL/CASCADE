"""
HAT_plot_gapfill.py

Review figure for the 2009 DEM gap fill: what the 2009 survey alone gives,
what the fill adds, and exactly which cells came from where.

Three panels, all on the same 10 m grid and the same colour scale:
    A  2009 only        cells the 2009 DEM measured; everything else blank
    B  2009 + fill      the product that goes to the dune/topo extractor
    C  survey source    2009 measured / fill-year filled / never surveyed

Panels A and B differ ONLY in the filled cells, so flipping between them shows
the fill directly. Panel C is the same information as a categorical map, which
is easier to read where the fill is thin.

Domain boxes are drawn over every panel: no fill, thin white outline, so they
locate a domain without hiding the data under it.

COLOUR
------
Elevation uses `terrain`, and the reason is specific: it is built for
topography, with a blue water band occupying the first 25% of the ramp. That
only reads correctly if sea level lands exactly on that internal break, hence
vmin is DERIVED - vmin = SEA_LEVEL_M - (vmax - SEA_LEVEL_M) / 3 - rather than
taken from a percentile. Set vmin from a percentile and the blue/green boundary
drifts to an arbitrary elevation, drawing dry ground as water.

The categorical panel's two colours are sampled from `terrain` itself,
terrain(0.05) water blue for 2009 and terrain(0.30) low-land green for the
fill, so it belongs to the same palette as the elevation maps. They carry a
deliberate luminance ladder - 80 / 153 / 228 against the grey, spacing 73 and
75 - because blue-vs-green is a weak colour-vision-deficiency axis and
brightness has to carry what hue cannot.

Nodata is neutral grey in every panel and never a step on the elevation ramp -
"not surveyed" is not a low elevation, and the whole point of this work is that
conflating those two drowned three roadways at t=0.

SOURCE SELECTION
----------------
Tag, fill year and long label live together in SOURCES, so a re-render cannot
put one source's label on another's data - hand-editing the three constants
separately already printed "cells the 2014 DEM measured" on the 2008 figures
once.

    python HAT_plot_gapfill.py                       # 2009-2014
    python HAT_plot_gapfill.py --source <PRODUCT>

INPUT   data/hatteras_init/0-elevation/<SOURCE_TAG>/2-resampled-10m/
            resampled_domain_<N>_filled.tif
            resampled_domain_<N>_survey.tif
OUTPUT  data/hatteras_init/0-elevation/<SOURCE_TAG>/figures/
            HAT_gapfill_<SOURCE_TAG>_island.png         whole island, 3 panels
            HAT_gapfill_<SOURCE_TAG>_domains_78_80.png  zoom on 78-80
            HAT_gapfill_<SOURCE_TAG>_roads_78_80.png    the zoom, + NC-12
        SOURCE_TAG is in the name, so a new source cannot overwrite the
        existing figures.

Requires: rasterio, geopandas, numpy, matplotlib
"""

from pathlib import Path
import re
import sys

import numpy as np
import geopandas as gpd
import rasterio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

def _find_project_root(start: Path) -> Path:
    """
    Walk up until a directory holds data/hatteras_init.

    NOT parents[N]. This file moved into 3-figures/ on 2026-08-25, and the
    old parents[3] then resolved to input_prep/ rather than the project root.
    That raises nothing - it just makes every path below it wrong, silently,
    until some glob comes back empty. Same helper and same reason as
    4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
    """
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
ELEVATION_DIR = PROJECT_ROOT / "data" / "hatteras_init" / "0-elevation"
IN_DIR = None   # set from SOURCE_TAG below
FIG_DIR = ELEVATION_DIR / "figures"
DOMAIN_FILE = Path(r"D:\Hatteras_GIS\domains.geojson")

# NC-12 alignments. These are EPSG:2264 (NC State Plane, US survey FEET) while
# the maps are EPSG:3725 (UTM 18N, metres), so they are reprojected on load -
# plotted raw they would land thousands of km off the map.
ROAD_DIR = (PROJECT_ROOT / "data" / "hatteras_init" / "4-mgmt-forcing"
            / "road_offset" / "raw_offset")
ROAD_FILES = {1984: ROAD_DIR / "1984" / "nc12_1984.geojson",
              2004: ROAD_DIR / "2004" / "nc12_2004.geojson"}
# The two alignments are very nearly coincident through 78-80. Drawn in the same
# colour with the solid one last, 2004 simply paints over 1984 and only ONE road
# appears. So 2004 is solid black underneath and 1984 is WHITE dashed on top:
# where they coincide you see a black line with white dashes (both present),
# and where they diverge each is legible on its own. Each carries a casing in
# the opposite colour so it survives terrain running from dark water to near-
# white dune crest.
ROAD_STYLE = {2004: dict(color="black", linestyle="-", linewidth=1.8),
              1984: dict(color="white", linestyle=(0, (4.5, 3.0)), linewidth=1.5)}
ROAD_CASING = {2004: dict(color="white", linewidth=3.4),
               1984: dict(color="black", linewidth=3.0)}
ROAD_ORDER = [2004, 1984]   # draw order: solid first, dashed on top

# Fill sources this script knows how to plot. Keeping tag, year and label in ONE
# place stops them drifting apart - hand-editing three constants per re-render
# already put "cells the 2014 DEM measured" on the 2008 figures once.
#
#     python HAT_plot_gapfill.py                      # default (below)
#     python HAT_plot_gapfill.py --source <PRODUCT>
#
# Keyed by PRODUCT, not by fill source - the directories were renamed for
# composition on 2026-08-25 (2014_NOAA_PostSandy -> 2009-2014).
#
# 2008_NOAA_IOCM was an entry here until 2026-08-26. Its rasters are gone from
# disk (they were never tracked - *.tif is gitignored), and the point-cloud
# path that produced them has been removed from HAT_dem_gap_fill.py, so
# --source 2008_NOAA_IOCM could not have re-rendered anything. A dead option
# that reads as a live one is worse than no option.
SOURCES = {
    "2009-2014": (
        2014,
        "2014 NOAA Post-Sandy DEM (Job1076021), 1 m, EPSG:6347 + NAVD88"),
}
DEFAULT_SOURCE = "2009-2014"

SOURCE_TAG = DEFAULT_SOURCE
if "--source" in sys.argv:
    SOURCE_TAG = sys.argv[sys.argv.index("--source") + 1]
if SOURCE_TAG not in SOURCES:
    raise SystemExit(f"unknown --source {SOURCE_TAG!r}; "
                     f"known: {', '.join(SOURCES)}")
SURVEY_FILL, SOURCE_LONG = SOURCES[SOURCE_TAG]

# Built from SURVEY_FILL so it cannot disagree with the data being plotted.
SOURCE_NOTE = (f"Fill is limited to cells the {SURVEY_FILL} source measured, "
               f"contiguous with the island (20 m bridging), above -2.64 m NAVD88. "
               f"No bias correction, no feathering - filled cells are the "
               f"{SURVEY_FILL} measurement unchanged.")

# The superseded fallback that used to live here is gone: the resolver knows
# which products are superseded and where they sit, so there is one place that
# has to be right rather than a probe-two-paths-and-hope in every consumer.
# A missing directory is reported there, not as an empty glob later.
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from hat_elevation_products import product as _product  # noqa: E402

_P = _product(SOURCE_TAG)
IN_DIR = _P.resampled_10m
FIG_DIR = _P.figures

GRID = 10.0

# Breathing room around the island-wide mosaic. Without it the northernmost and
# southernmost domains sit flush against the axes frame, which reads as the data
# being cut off rather than ending.
ISLAND_PAD_M = 700.0

# Zoom figure geometry. Height is derived from the data aspect at draw time;
# ZOOM_CHROME_IN is the vertical allowance for the suptitle, colorbar and tick
# labels, which do not scale with the map. NOT the footnote - that sits outside
# the axes and bbox_inches="tight" adds it after layout, so reserving space for
# it just opens a gap under the title. Measured: gap above the axes tracks this
# value almost 1:1, and a two-line suptitle needs ~0.35 in.
# Zoom legends sit upper-left: the island occupies the centre-right of these
# panels, so the top-left corner is the one reliably empty area. The ISLAND
# figure keeps lower-right - there the island runs bottom-left to top-right and
# the bottom-right corner is the clear one.
ZOOM_LEGEND_LOC = "upper left"

# Road-overlay zooms: (domain ids, which NC-12 years, filename slug, title).
#
# 78-80 carries BOTH alignments - those are the domains the extractor names as
# width-drowning at t=0, and seeing 1984 against 2004 there is the point.
#
# 8-15 carries BOTH as well. An earlier version drew only 2004 here, on the
# grounds that it is the alignment contemporaneous with this DEM and that
# putting the 1984 line over a 2009+2014 surface invites comparing a road to a
# DEM holding no information from its era. That was overruled deliberately:
# seeing where the road WAS against where it WENT is the point, and this view
# is meant to be read as a pair with the 1984-start DEM's own 8-15 figure,
# which now draws the same two lines.
ROAD_ZOOMS = [
    ([78, 79, 80], [2004, 1984], "roads_78_80",
     "Domains 78-80 with NC-12 alignments"),
    (list(range(8, 16)), [2004, 1984], "roads_8_15",
     "Domains 8-15 with NC-12 alignments"),
    (list(range(82, 89)), [2004, 1984], "roads_82_88",
     "Domains 82-88 with NC-12 alignments"),
]

ZOOM_FIG_W = 15.0
ZOOM_CHROME_IN = 1.0
ID_RE = re.compile(r"resampled_domain_(\w+)_filled\.tif$")

SURVEY_2009, SURVEY_NONE = 2009, 0

# Categorical colours for the survey-source panel, SAMPLED FROM `terrain` so
# panel C is built from the same palette as the elevation maps:
#
#     C_2009  terrain(0.05)  #2353b9  deep water blue
#     C_FILL  terrain(0.30)  #31d670  the green of terrain's low-land band
#
# Luminance ladder, which is what keeps the three readable:
#
#     2009 measured  #2353b9   luminance  80
#     fill           #31d670   luminance 153   (73 from blue, 75 from grey)
#     never surveyed #E4E4E4   luminance 228
#
# Spacing 73 and 78 is nearly even, so all three separate by brightness alone.
# That matters more here than usual: blue-vs-green is the WEAKEST colour-vision
# -deficiency axis of the pairings tried, so brightness is doing the work that
# hue cannot be relied on for. Without the gap this pair would be a poor choice;
# with it, it holds up in greyscale and under CVD.
#
# The dataviz validator could not be run here (no node on this machine), so the
# ladder was computed directly rather than machine-checked.
C_2009 = "#2353b9"   # terrain's water blue
C_FILL = "#31d670"   # terrain's low-land green
C_NONE = "#E4E4E4"   # neutral grey, off the elevation ramp entirely

ELEV_CMAP = "terrain"
ELEV_PCT = (2, 98)   # clip the ramp to percentiles so a few spikes don't flatten it

# matplotlib's `terrain` is built for topography: its blue water band occupies
# the FIRST 25% of the ramp, then green -> brown -> white for land. That is only
# meaningful if sea level lands exactly on that internal boundary, so vmin is
# derived rather than taken from a percentile:
#
#     0 maps to  |vmin| / (|vmin| + vmax)  ==  0.25   ->   vmin = -vmax / 3
#
# Set from a percentile instead and the blue/green break drifts to some
# arbitrary elevation, so the map would draw dry ground as water or vice versa.
SEA_LEVEL_M = 0.0        # m NAVD88; use MHW (0.36) to key the break to MHW
TERRAIN_WATER_FRAC = 0.25


def load_mosaic():
    """Places every domain on one 10 m grid covering the island."""
    paths = sorted(IN_DIR.glob("resampled_domain_*_filled.tif"))
    if not paths:
        raise FileNotFoundError(f"no domain rasters in {IN_DIR} - run steps 1-2 first")

    boxes = []
    for p in paths:
        with rasterio.open(p) as s:
            t = s.transform
            boxes.append((t.c, t.f - s.height * GRID, t.c + s.width * GRID, t.f))
    minx = min(b[0] for b in boxes); miny = min(b[1] for b in boxes)
    maxx = max(b[2] for b in boxes); maxy = max(b[3] for b in boxes)

    W = int(round((maxx - minx) / GRID))
    H = int(round((maxy - miny) / GRID))
    elev = np.full((H, W), np.nan)
    surv = np.zeros((H, W), np.uint16)

    for p in paths:
        dom = ID_RE.search(p.name).group(1)
        with rasterio.open(p) as s:
            a = s.read(1).astype(float)
            nd = s.nodata
            t = s.transform
        if nd is not None and not np.isnan(nd):
            a = np.where(a == nd, np.nan, a)
        sp = IN_DIR / f"resampled_domain_{dom}_survey.tif"
        if sp.exists():
            with rasterio.open(sp) as s:
                sv = s.read(1).astype(np.uint16)
        else:
            sv = np.where(np.isnan(a), SURVEY_NONE, SURVEY_2009).astype(np.uint16)

        c0 = int(round((t.c - minx) / GRID))
        r0 = int(round((maxy - t.f) / GRID))
        h, w = a.shape
        sub_e = elev[r0:r0 + h, c0:c0 + w]
        sub_s = surv[r0:r0 + h, c0:c0 + w]
        # Domains are placed independently and their boxes are not a perfect
        # tiling (505 m spacing against a 500 m extent), so overlaps exist.
        # Keep whatever is already there rather than letting the later domain
        # silently overwrite its neighbour.
        take = np.isnan(sub_e) & ~np.isnan(a)
        sub_e[take] = a[take]
        sub_s[take] = sv[take]
        fresh = (sub_s == 0) & (sv != 0)
        sub_s[fresh] = sv[fresh]

    extent = [minx, maxx, miny, maxy]
    return elev, surv, extent, len(paths)


def draw_domains(ax, gdf, lw=0.35):
    gdf.boundary.plot(ax=ax, color="white", linewidth=lw, zorder=5)


def km_axes(ax, nx=3, ny=6):
    """
    UTM eastings here are 6-digit metres (450439..458392). Three panels side by
    side cannot fit those without colliding, and matplotlib's shared offset
    label is easy to miss. Kilometres with a small tick count is legible at any
    panel width and needs no offset text.
    """
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nx, prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=ny))

    # Decimal places from the span, not fixed: the island figure covers ~45 km
    # of northing where 0 dp is right, the zoom covers ~1.5 km where 0 dp
    # renders every tick as the same number.
    def dp(span_m):
        return 0 if span_m > 20_000 else (1 if span_m > 2_000 else 2)

    xs = abs(np.diff(ax.get_xlim())[0])
    ys = abs(np.diff(ax.get_ylim())[0])
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _, d=dp(xs): f"{v / 1000:.{d}f}"))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _, d=dp(ys): f"{v / 1000:.{d}f}"))
    ax.tick_params(labelsize=8)


def load_roads(dst_crs, clip_to=None):
    """
    Loads the NC-12 alignments, reprojects them, and CLIPS them to the domain
    footprint.

    The geojsons run the full length of the highway, well beyond the 90 domains
    at both ends. Unclipped, the island figure shows road where there is no
    model domain, which reads as coverage that does not exist.
    """
    out = {}
    for yr, f in ROAD_FILES.items():
        if not f.exists():
            print(f"  WARNING: {f} missing - {yr} road not drawn")
            continue
        g = gpd.read_file(f)
        if g.crs is not None:
            g = g.to_crs(dst_crs)
        if clip_to is not None:
            before = float(g.geometry.length.sum())
            g = gpd.clip(g, clip_to)
            after = float(g.geometry.length.sum()) if len(g) else 0.0
            print(f"  NC-12 {yr}: clipped to domains, "
                  f"{after / 1000:.1f} km of {before / 1000:.1f} km kept")
        if len(g):
            out[yr] = g
    return out


def draw_roads(ax, roads, scale=1.0):
    """
    Casing then line, in ROAD_ORDER so the dashed 1984 lands on top of the solid
    2004 rather than under it (see ROAD_STYLE). `scale` thins the lines for the
    island-wide figure, where the same widths would smother the island.
    """
    for yr in ROAD_ORDER:
        if yr not in roads:
            continue
        cas = dict(ROAD_CASING[yr]); cas["linewidth"] *= scale
        roads[yr].plot(ax=ax, linestyle="-", zorder=6 + ROAD_ORDER.index(yr) * 2,
                       alpha=0.9, **cas)
        st = dict(ROAD_STYLE[yr]); st["linewidth"] *= scale
        roads[yr].plot(ax=ax, zorder=7 + ROAD_ORDER.index(yr) * 2, **st)


def road_legend_handles(roads):
    return [Line2D([], [], label=f"NC-12 {y}",
                   **{**ROAD_STYLE[y],
                      "color": "black" if ROAD_STYLE[y]["color"] == "white"
                      else ROAD_STYLE[y]["color"]})
            for y in ROAD_ORDER if y in roads]


def panel_elev(ax, arr, extent, vmin, vmax, title):
    ax.set_facecolor(C_NONE)
    im = ax.imshow(arr, extent=extent, origin="upper", cmap=ELEV_CMAP,
                   vmin=vmin, vmax=vmax, interpolation="nearest", zorder=1)
    ax.set_title(title, fontsize=11, pad=8)
    return im


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    elev, surv, extent, n = load_mosaic()
    print(f"mosaic {elev.shape} from {n} domains, extent {[round(v) for v in extent]}")

    gdf = gpd.read_file(DOMAIN_FILE)
    # union of the 90 domain boxes - the road is shown only where a domain is
    dom_union = (gdf.union_all() if hasattr(gdf, 'union_all')
                 else gdf.unary_union)
    roads = load_roads(gdf.crs, clip_to=dom_union)

    only09 = np.where(surv == SURVEY_2009, elev, np.nan)
    filled = elev
    n_fill = int((surv == SURVEY_FILL).sum())
    n_meas = int((surv == SURVEY_2009).sum())
    print(f"cells: 2009 measured {n_meas:,}  {SURVEY_FILL} filled {n_fill:,} "
          f"({100 * n_fill / max(n_meas + n_fill, 1):.1f}% of land cells)")

    valid = filled[~np.isnan(filled)]
    _, vmax = np.percentile(valid, ELEV_PCT)
    # place SEA_LEVEL_M exactly on terrain's water/land break (see above)
    f = TERRAIN_WATER_FRAC
    vmin = SEA_LEVEL_M - (vmax - SEA_LEVEL_M) * f / (1 - f)
    n_clip = int((valid < vmin).sum())
    print(f"elev ramp: {vmin:.2f} .. {vmax:.2f} m, sea level {SEA_LEVEL_M:.2f} "
          f"at {100 * f:.0f}% of terrain; {n_clip:,} cells clip low "
          f"({100 * n_clip / valid.size:.2f}%)")

    fig, axes = plt.subplots(1, 3, figsize=(13, 19), constrained_layout=True)

    im = panel_elev(axes[0], only09, extent, vmin, vmax,
                    f"A  2009 DEM only\n{n_meas:,} cells")
    panel_elev(axes[1], filled, extent, vmin, vmax,
               f"B  2009 + {SURVEY_FILL} fill\n{n_meas + n_fill:,} cells "
               f"(+{n_fill:,})")

    # Boundaries must ascend and the fill year (2014) is now GREATER than the
    # measured year (2009), so measured precedes fill in the colour list. With
    # the old 2008 ordering these two were swapped and the map lied.
    lo, hi = sorted([SURVEY_2009, SURVEY_FILL])
    cmap_s = ListedColormap([C_NONE,
                             C_2009 if lo == SURVEY_2009 else C_FILL,
                             C_FILL if hi == SURVEY_FILL else C_2009])
    norm_s = BoundaryNorm([-0.5, 0.5, lo + 0.5, hi + 0.5], cmap_s.N)
    axes[2].set_facecolor(C_NONE)
    axes[2].imshow(surv, extent=extent, origin="upper", cmap=cmap_s, norm=norm_s,
                   interpolation="nearest", zorder=1)
    axes[2].set_title("C  survey source", fontsize=11, pad=8)

    for ax in axes:
        draw_domains(ax, gdf)
        # thinner at island scale: 45 km of line at zoom widths would smother
        # the island. Dashed vs solid does not resolve at this scale - the zoom
        # figure is where that distinction is readable.
        draw_roads(ax, roads, scale=0.45)
        ax.set_xlim(extent[0] - ISLAND_PAD_M, extent[1] + ISLAND_PAD_M)
        ax.set_ylim(extent[2] - ISLAND_PAD_M, extent[3] + ISLAND_PAD_M)
        ax.set_xlabel("Easting (km)", fontsize=9)
        km_axes(ax)
        ax.set_aspect("equal")
    axes[0].set_ylabel("Northing (km)", fontsize=9)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    # ax=all three, not axes[:2] - a colorbar sized against a subset shrinks
    # only those axes and leaves the rest misaligned.
    cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                      fraction=0.035, pad=0.01, aspect=45)
    cb.set_label("Elevation (m NAVD88)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    axes[2].legend(handles=[Patch(facecolor=C_2009, label="2009 measured"),
                            Patch(facecolor=C_FILL, label=f"{SURVEY_FILL} fill"),
                            Patch(facecolor=C_NONE, label="never surveyed")]
                   + road_legend_handles(roads),
                   loc="lower right", fontsize=8, framealpha=0.9)

    fig.suptitle(f"Hatteras 2009 DEM gap fill\nfill source: {SOURCE_LONG}",
                 fontsize=12)
    fig.text(0.5, -0.012, SOURCE_NOTE, ha="center", va="top", fontsize=8,
             wrap=True, color="#444444")
    out = FIG_DIR / f"HAT_gapfill_{SOURCE_TAG}_island.png"
    fig.savefig(out, dpi=170, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # ---- zoom: the domains the extractor names as width-drowning at t=0 ----
    sel = gdf[gdf["domain_id"].astype(int).isin([78, 79, 80])]
    if not sel.empty:
        zminx, zminy, zmaxx, zmaxy = sel.total_bounds
        pad = 150
        # Height from the DATA aspect, not hard-coded. The panels are
        # set_aspect("equal") and the zoom extent is wider than it is tall, so
        # matplotlib shrinks each axes to match and a fixed tall figure leaves
        # slack that constrained_layout splits above and below the panels -
        # which reads as a big empty gap under the title.
        _zw = (zmaxx + pad) - (zminx - pad)
        _zh = (zmaxy + pad) - (zminy - pad)
        _panel_w = ZOOM_FIG_W / 3.0
        _fig_h = _panel_w / (_zw / _zh) + ZOOM_CHROME_IN
        fig, axes = plt.subplots(1, 3, figsize=(ZOOM_FIG_W, _fig_h),
                                 constrained_layout=True)
        im = panel_elev(axes[0], only09, extent, vmin, vmax, "A  2009 DEM only")
        panel_elev(axes[1], filled, extent, vmin, vmax,
                   f"B  2009 + {SURVEY_FILL} fill")
        axes[2].set_facecolor(C_NONE)
        axes[2].imshow(surv, extent=extent, origin="upper", cmap=cmap_s,
                       norm=norm_s, interpolation="nearest", zorder=1)
        axes[2].set_title("C  survey source", fontsize=11, pad=8)
        for ax in axes:
            draw_domains(ax, gdf, lw=0.9)
            ax.set_xlim(zminx - pad, zmaxx + pad)
            ax.set_ylim(zminy - pad, zmaxy + pad)
            ax.set_aspect("equal")
            km_axes(ax, nx=3, ny=5)
            ax.set_xlabel("Easting (km)", fontsize=9)
        axes[0].set_ylabel("Northing (km)", fontsize=9)
        axes[2].legend(handles=[Patch(facecolor=C_2009, label="2009 measured"),
                                Patch(facecolor=C_FILL, label=f"{SURVEY_FILL} fill"),
                                Patch(facecolor=C_NONE, label="never surveyed")],
                       loc=ZOOM_LEGEND_LOC, fontsize=8, framealpha=0.9)
        cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                          fraction=0.05, pad=0.02, aspect=40)
        cb.set_label("Elevation (m NAVD88)", fontsize=9)
        fig.suptitle("Domains 78-80 - the roadways that width-drowned at t=0 "
                     f"on missing survey coverage\nfill source: {SOURCE_LONG}",
                     fontsize=11)
        fig.text(0.5, -0.02, SOURCE_NOTE, ha="center", va="top", fontsize=8,
                 wrap=True, color="#444444")
        out2 = FIG_DIR / f"HAT_gapfill_{SOURCE_TAG}_domains_78_80.png"
        fig.savefig(out2, dpi=170, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out2}")

        # ---- road-overlay zooms, one per entry in ROAD_ZOOMS ----
        # Drawn at zoom rather than island scale on purpose: island-wide the
        # road is a ~1 px line over 45 km, where dashed and solid are
        # indistinguishable and the overlay would carry no information.
        def _road_zoom(dom_ids, years, slug, title):
            """One A/B/C zoom with the named NC-12 alignments drawn over it.

            `years` selects WHICH alignments, and every entry in ROAD_ZOOMS
            currently asks for both. It stays a parameter rather than being
            hardcoded because the restriction was tried and reversed once
            already - see the note on ROAD_ZOOMS - and a future zoom may well
            want one line only.
            """
            zsel = gdf[gdf["domain_id"].astype(int).isin(dom_ids)]
            if zsel.empty:
                print(f"  domains {dom_ids} not in the domain file - "
                      f"{slug} skipped")
                return
            rsub = {y: g for y, g in roads.items() if y in years}
            if not rsub:
                print(f"  no road data for {years} - {slug} skipped")
                return
            zx0, zy0, zx1, zy1 = zsel.total_bounds
            pad_ = 150
            zw, zh = (zx1 + pad_) - (zx0 - pad_), (zy1 + pad_) - (zy0 - pad_)
            fh = (ZOOM_FIG_W / 3.0) / (zw / zh) + ZOOM_CHROME_IN
            fig, axes = plt.subplots(1, 3, figsize=(ZOOM_FIG_W, fh),
                                     constrained_layout=True)
            im = panel_elev(axes[0], only09, extent, vmin, vmax,
                            "A  2009 DEM only")
            panel_elev(axes[1], filled, extent, vmin, vmax,
                       f"B  2009 + {SURVEY_FILL} fill")
            axes[2].set_facecolor(C_NONE)
            axes[2].imshow(surv, extent=extent, origin="upper", cmap=cmap_s,
                           norm=norm_s, interpolation="nearest", zorder=1)
            axes[2].set_title("C  survey source", fontsize=11, pad=8)
            for ax in axes:
                draw_domains(ax, gdf, lw=0.9)
                draw_roads(ax, rsub)
                ax.set_xlim(zx0 - pad_, zx1 + pad_)
                ax.set_ylim(zy0 - pad_, zy1 + pad_)
                ax.set_aspect("equal")
                km_axes(ax, nx=3, ny=5)
                ax.set_xlabel("Easting (km)", fontsize=9)
            axes[0].set_ylabel("Northing (km)", fontsize=9)
            rh = road_legend_handles(rsub)
            axes[0].legend(handles=rh, loc=ZOOM_LEGEND_LOC, fontsize=8,
                           framealpha=0.9)
            axes[2].legend(
                handles=[Patch(facecolor=C_2009, label="2009 measured"),
                         Patch(facecolor=C_FILL, label=f"{SURVEY_FILL} fill"),
                         Patch(facecolor=C_NONE, label="never surveyed")] + rh,
                loc=ZOOM_LEGEND_LOC, fontsize=8, framealpha=0.9)
            cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                              fraction=0.05, pad=0.02, aspect=40)
            cb.set_label("Elevation (m NAVD88)", fontsize=9)
            fig.suptitle(f"{title}\nfill source: {SOURCE_LONG}", fontsize=11)
            fig.text(0.5, -0.02, SOURCE_NOTE, ha="center", va="top",
                     fontsize=8, wrap=True, color="#444444")
            outp = FIG_DIR / f"HAT_gapfill_{SOURCE_TAG}_{slug}.png"
            fig.savefig(outp, dpi=170, bbox_inches="tight")
            plt.close(fig)
            print(f"wrote {outp}")

        if roads:
            for _ids, _yrs, _slug, _title in ROAD_ZOOMS:
                _road_zoom(_ids, _yrs, _slug, _title)


if __name__ == "__main__":
    main()
