"""
HAT_plot_1984_mosaic.py

Review figures for the 1984-start topography: what each of the THREE sources
contributed, and how far the beach start moved as a result.

HAT_plot_gapfill.py draws the two-source 2009-start product and is left alone.
This is a separate script rather than a fourth entry in its SOURCES dict
because the 1984 product is a different thing to look at - it has a third
source, a boundary (the 1984 NC-12 line) that the two-source product has no
concept of, and a downstream consequence (the window origin moving) that is the
main reason to review it at all.

OUTPUTS (data/hatteras_init/0-elevation/2009-2014-1996/figures/)
    HAT_mosaic1984_island.png      3 panels, whole island
    HAT_mosaic1984_zoom_76_81.png  the developed reach, + both NC-12 lines
    HAT_mosaic1984_roads_8_15.png  the southern end, + both NC-12 lines
    HAT_mosaic1984_roads_82_88.png the northern reach, + both NC-12 lines
    HAT_mosaic1984_shift.png       per-domain start_beach shift

    Each product owns its figures/ folder, so these cannot collide with
    the baseline product's HAT_gapfill_* figures.

THE THREE PANELS
----------------
    (a) 2009 survey         what the base survey measured; everything else blank
    (b) 1984-start surface  1996 ocean-side, 2009, 2014 - the product
    (c) survey source       which survey each cell came from

(a) and (b) differ only in filled/overridden cells, so flipping between them
shows the whole intervention at once. (c) is the same information as a
categorical map, which is the only readable form where the 1996 band is thin.

STYLE. Every figure here is drawn under `scripts/hat_figure_style.py` at the
printed width (190 mm), and carries no title, statistics line or footnote on the
canvas: that text is written to CAPTIONS.md beside the PNGs. The terrain ramp is
the house style's one sanctioned exception to drawing elevation in classes.

COLOUR
------
Elevation uses `terrain` with a DERIVED vmin, for the reason HAT_plot_gapfill.py
gives at length: terrain's blue water band occupies the first 25% of the ramp,
so sea level has to land exactly on that internal break or the map draws dry
ground as water. vmin = SEA_LEVEL - (vmax - SEA_LEVEL)/3, never a percentile.

The categorical panel extends that script's existing ladder rather than
inventing a palette, so the two products' figures stay comparable - a reader
who has learned that blue means "2009 measured" and green means "filled" does
not have to relearn it here:

    2009 measured   #2353b9   terrain(0.05), water blue      Y =  26
    1996 ALACE      #e6550d   ColorBrewer Oranges            Y =  60
    2014 fill       #31d670   terrain(0.30), low-land green  Y = 127
    never surveyed  #E4E4E4   neutral grey, off the ramp     Y = 198

The 1996 step is NOT sampled from terrain, and that is deliberate. Every
terrain slot in the usable lightness band is a blue, a green or a tan, and all
of them sit too close to the two colours already spoken for once the palette is
run through a colour-vision-deficiency check. An orange is the nearest hue
family terrain does not use, and it also reads correctly as "the thing this
product added".

VALIDATED, NOT EYEBALLED. The dataviz skill's JS validator needs node, which is
not on this machine - the same wall HAT_plot_gapfill.py hit - so the identical
checks were computed in Python: OKLab dE x100 between every pair, under normal
vision and under Vienot deuteranopia / protanopia / tritanopia simulation.

    worst pair, any deficiency   1996 vs 2014, protanopia   dE 16.6   (target >= 8)
    worst pair, normal vision    1996 vs 2014               dE 34.2   (floor 15)

The luminance ladder 26 / 60 / 127 / 198 is monotonic with gaps of 34, 67 and
71, so all four also separate in greyscale and in print with no hue at all.

ONE PRE-EXISTING WEAKNESS, INHERITED AND NOT INTRODUCED: the 2009-blue against
2014-green pair scores dE 5.6 under TRITANOPIA, below the target. It comes from
HAT_plot_gapfill.py's existing palette, not from anything added here, and
changing those two would break comparability with every figure already
published from the 2009-start product. Tritanopia is far rarer than the
red-green deficiencies (both of which that pair clears comfortably, at 25.5 and
30.3), and the luminance ladder carries the distinction regardless. Recorded so
it is a known trade rather than an unnoticed one.

Nodata is neutral grey in every panel and never a step on the elevation ramp:
"not surveyed" is not a low elevation, and conflating the two is what drowned
three roadways at t=0 in the first place.

    python HAT_plot_1984_mosaic.py

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
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
ELEVATION_DIR = INIT_ROOT / "0-elevation"

SOURCE_TAG = "2009-2014-1996"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from hat_elevation_products import product as _product  # noqa: E402
from hat_figure_style import (  # noqa: E402
    apply_style, figsize, save, caption, C, C_1984, C_1997, INK, INK_MUTED,
    DOMAIN_AXIS_LABEL, town_bands, open_frame, spines_for_image, _title)

apply_style()

_P = _product(SOURCE_TAG)
IN_DIR = _P.resampled_10m
AUDIT_1M = _P.audit_1m
FIG_DIR = _P.figures

DOMAIN_FILE = Path(r"D:\Hatteras_GIS\domains.geojson")
# NC-12 alignments. EPSG:2264 (NC State Plane, US survey FEET) while the maps
# are EPSG:3725 (UTM 18N, metres), so they are reprojected on load - plotted raw
# they would land thousands of km off the map.
#
# BOTH are drawn, and the styling matches HAT_plot_gapfill.py so a reader moving
# between the two products' figures does not have to relearn the key. Two
# vintages of the same line, so they take the house vintage pair: the EARLIER
# alignment (1984) red, the LATER one (2004) blue. The two are very nearly
# coincident over much of the island, so 2004 is solid underneath and 1984 is
# dashed on top - where they coincide you see a blue line with red dashes, and
# where they diverge each is legible alone. Both carry a white casing so they
# survive terrain running from dark water to near-white dune crest.
ROAD_DIR = (INIT_ROOT / "4-mgmt-forcing" / "road_offset" / "raw_offset")
ROAD_FILES = {1984: ROAD_DIR / "1984" / "nc12_1984.geojson",
              2004: ROAD_DIR / "2004" / "nc12_2004.geojson"}
ROAD_STYLE = {2004: dict(color=C_1997, linestyle="-", linewidth=1.5),
              1984: dict(color=C_1984, linestyle=(0, (3.6, 2.4)), linewidth=1.5)}
ROAD_CASING = {2004: dict(color="white", linewidth=3.0),
               1984: dict(color="white", linewidth=3.0)}
ROAD_ORDER = [2004, 1984]   # draw order: solid first, dashed on top

GRID = 10.0
ISLAND_PAD_M = 700.0

# Chrome geometry, matching HAT_plot_gapfill.py exactly. Height is DERIVED from
# the data aspect at draw time: the panels are set_aspect("equal") and a zoom
# extent is usually wider than tall, so matplotlib shrinks each axes to match
# and a fixed tall figure leaves slack that constrained_layout splits above and
# below - which reads as a large empty gap under the title. ZOOM_CHROME_IN is
# the vertical allowance for panel titles, colorbar and tick labels, which do
# not scale with the map.
#
# The width is the house double-column width (190 mm) since 2026-09-10: a figure
# is drawn at the width it is printed, so its 8-9 pt type is 8-9 pt on the page.
# It was 15 in, which reduced to a page turned every label into 4 pt.
ZOOM_FIG_W = figsize("double")[0]
ZOOM_CHROME_IN = 1.05
# Upper-left: on these zooms the island runs up the centre-right of the frame,
# so the top-left corner is the one reliably empty area.
ZOOM_LEGEND_LOC = "upper left"

# (domain ids, filename slug, title). BOTH NC-12 alignments are drawn on every
# zoom. An earlier version restricted each product to the alignment
# contemporaneous with it - 1984 here, 2004 on the baseline DEM - on the
# grounds that comparing a road to a DEM holding no information from its era
# invites a false reading. That was overruled deliberately: seeing where the
# road WAS against where it WENT is the point of the comparison, and the two
# products' 8-15 views are meant to be read as a pair.
ZOOMS = [
    (list(range(76, 82)), "zoom_76_81",
     "Domains 76-81, the developed reach, where the beach start moves furthest "
     "seaward (47-72 m, 5-7 Barrier3D cells). (a) the 2009 survey alone; "
     "(b) the surface the extractor reads; (c) the survey each cell came from. "
     "Both NC-12 alignments are drawn: 1984 dashed red, 2004 solid blue."),
    (list(range(8, 16)), "roads_8_15",
     "Domains 8-15, the southern end, where the NC-12 exports begin; domains "
     "1-7 have no road line and, since the boundary was dropped, take 1996 "
     "anyway. (a) the 2009 survey alone; (b) the surface the extractor reads; "
     "(c) the survey each cell came from. Both NC-12 alignments are drawn: "
     "1984 dashed red, 2004 solid blue."),
    # Subtitle numbers are from mosaic_1984_audit.csv, not from eyeballing the
    # map. The first draft said "mostly NEW land", which the audit contradicts:
    # this reach is 28% new against 72% overwrite. RE-READ AFTER THE 2026-08-26
    # NO-BOUNDARY REBUILD: it used to be 45% against an 17% island-wide share,
    # i.e. nearly 3x. Dropping the boundary admitted a large band of backdune
    # that 2009 HAD seen, so both shares fell and the ratio with them - the
    # reach is now 1.75x the island figure, not 3x. Still the highest on the
    # island, which is the point, but the figure must not keep claiming 3x.
    (list(range(82, 89)), "roads_82_88",
     "Domains 82-88, the northern reach, where 28% of what 1996 writes is land "
     "the 2009 survey never saw (16% island-wide) and the mean beach-start "
     "shift is +29 m. (a) the 2009 survey alone; (b) the surface the extractor "
     "reads; (c) the survey each cell came from. Both NC-12 alignments are "
     "drawn: 1984 dashed red, 2004 solid blue."),
]

ID_RE = re.compile(r"resampled_domain_(\w+)_filled\.tif$")

SURVEY_NONE, SURVEY_1996, SURVEY_2009, SURVEY_2014 = 0, 1996, 2009, 2014

C_1996 = "#e6550d"   # ColorBrewer Oranges - see the module docstring
C_2009 = "#2353b9"   # terrain(0.05), water blue
C_2014 = "#31d670"   # terrain(0.30), low-land green
C_NONE = "#E4E4E4"   # neutral grey, off the elevation ramp entirely

ELEV_CMAP = "terrain"
ELEV_PCT_HI = 98
SEA_LEVEL_M = 0.0
TERRAIN_WATER_FRAC = 0.25

# The Barrier3D cell. The shift figure's whole point is which domains moved by
# at least one of these, so it is drawn rather than left to the reader.
CELL_M = 10.0


# =============================================================================
# LOAD
# =============================================================================

def load_mosaic():
    """Every domain placed on one 10 m grid covering the island."""
    paths = sorted(IN_DIR.glob("resampled_domain_*_filled.tif"))
    if not paths:
        raise SystemExit(
            f"no resampled rasters in {IN_DIR}\n"
            f"  run HAT_dem_1984_mosaic.py, then\n"
            f"  HAT_dem_resample_clip.py --product {SOURCE_TAG}")

    boxes = []
    for p in paths:
        with rasterio.open(p) as s:
            t = s.transform
            boxes.append((t.c, t.f - s.height * GRID,
                          t.c + s.width * GRID, t.f))
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
            nd, t = s.nodata, s.transform
        if nd is not None and not np.isnan(nd):
            a = np.where(a == nd, np.nan, a)
        sp = IN_DIR / f"resampled_domain_{dom}_survey.tif"
        with rasterio.open(sp) as s:
            sv = s.read(1).astype(np.uint16)

        c0 = int(round((t.c - minx) / GRID))
        r0 = int(round((maxy - t.f) / GRID))
        h, w = a.shape
        sub_e, sub_s = elev[r0:r0 + h, c0:c0 + w], surv[r0:r0 + h, c0:c0 + w]
        # Domain boxes are 505 m apart on a 500 m extent, so they overlap.
        # First writer wins, as in HAT_plot_gapfill.load_mosaic - a later
        # domain must not silently repaint its neighbour.
        take = np.isnan(sub_e) & ~np.isnan(a)
        sub_e[take] = a[take]
        sub_s[take] = sv[take]
        fresh = (sub_s == 0) & (sv != 0)
        sub_s[fresh] = sv[fresh]

    return elev, surv, [minx, maxx, miny, maxy], len(paths)


# =============================================================================
# DRAW
# =============================================================================

def km_axes(ax, nx=3, ny=6):
    """UTM eastings are 6-digit metres and collide at panel width; km do not."""
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nx, prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=ny))

    def dp(span):
        return 0 if span > 20_000 else (1 if span > 2_000 else 2)

    xs = abs(np.diff(ax.get_xlim())[0])
    ys = abs(np.diff(ax.get_ylim())[0])
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _, d=dp(xs): f"{v / 1000:.{d}f}"))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _, d=dp(ys): f"{v / 1000:.{d}f}"))
    ax.tick_params(labelsize=8)


def elev_limits(elev):
    """vmax from a percentile, vmin DERIVED so 0 m lands on terrain's internal
    water/land break. Never both from percentiles - see the docstring."""
    v = elev[np.isfinite(elev)]
    vmax = float(np.percentile(v, ELEV_PCT_HI))
    vmin = SEA_LEVEL_M - (vmax - SEA_LEVEL_M) * (
        TERRAIN_WATER_FRAC / (1 - TERRAIN_WATER_FRAC))
    return vmin, vmax


def panel_elev(ax, i, arr, extent, vmin, vmax, title):
    ax.set_facecolor(C_NONE)
    im = ax.imshow(arr, extent=extent, origin="upper", cmap=ELEV_CMAP,
                   vmin=vmin, vmax=vmax, interpolation="nearest", zorder=1)
    _title(ax, i, title)
    return im


def panel_survey(ax, i, surv, extent, title):
    """Categorical provenance. Codes are mapped to contiguous indices so the
    colours cannot slide if a code is absent from a crop."""
    codes = [SURVEY_NONE, SURVEY_1996, SURVEY_2009, SURVEY_2014]
    cols = [C_NONE, C_1996, C_2009, C_2014]
    idx = np.zeros(surv.shape, np.uint8)
    for i_c, c in enumerate(codes):
        idx[surv == c] = i_c
    cmap = ListedColormap(cols)
    norm = BoundaryNorm(np.arange(-0.5, len(codes) + 0.5), cmap.N)
    ax.set_facecolor(C_NONE)
    ax.imshow(idx, extent=extent, origin="upper", cmap=cmap, norm=norm,
              interpolation="nearest", zorder=1)
    _title(ax, i, title)


def survey_legend(concise=False):
    """Handles for the three surveys plus the unsurveyed background.

    concise=True keeps only the year. The island figure needs it: its panels
    are ~2 in wide (a 46 km island in an 8 km window at equal aspect) and the
    full labels are wider than the axes they sit in, so the legend spilled
    across the panel border. The sources are named in this module's docstring,
    so the year is enough to read the colours by.
    """
    if concise:
        return [Patch(facecolor=C_1996, label="1996"),
                Patch(facecolor=C_2009, label="2009"),
                Patch(facecolor=C_2014, label="2014"),
                Patch(facecolor=C_NONE, label="no survey")]
    return [Patch(facecolor=C_1996, label="1996 ALACE (no road boundary)"),
            Patch(facecolor=C_2009, label="2009 USACE, measured"),
            Patch(facecolor=C_2014, label="2014 NOAA Post-Sandy, gap fill"),
            Patch(facecolor=C_NONE, label="never surveyed")]


def counts_note(surv):
    n = {c: int((surv == c).sum())
         for c in (SURVEY_1996, SURVEY_2009, SURVEY_2014)}
    tot = sum(n.values())
    return (f"10 m cells: 1996 {n[SURVEY_1996]:,} ({100 * n[SURVEY_1996] / tot:.1f}%), "
            f"2009 {n[SURVEY_2009]:,} ({100 * n[SURVEY_2009] / tot:.1f}%), "
            f"2014 {n[SURVEY_2014]:,} ({100 * n[SURVEY_2014] / tot:.1f}%)")


# The method paragraph. It used to be printed under every figure as a
# `fig.text` footnote; under the house style nothing on the canvas belongs in a
# caption, so it is now the tail of each CAPTIONS.md entry instead.
METHOD = ("1996 is admitted wherever it has data — there is no road boundary, "
          "and the landward limit is the ALACE swath edge: above -2.64 m "
          "NAVD88 where no other survey saw the cell and above MHW where one "
          "did, below a 12 m ceiling, contiguous with the island. No bias "
          "correction and no feathering: every cell is its own survey's "
          "measurement, unchanged. Axes are UTM eastings and northings in km; "
          "domain outlines are white.")

# The legend swatches carry only the year, because the panels are too narrow for
# the full wording; the surveys are named here instead.
SURVEY_KEY = ("Survey sources: 1996 ALACE lidar (admitted with no road "
              "boundary), 2009 USACE, 2014 NOAA Post-Sandy gap fill; grey is "
              "never surveyed.")


def load_roads(dst_crs, clip_to=None):
    """
    Both NC-12 alignments, reprojected and CLIPPED to the domain footprint.

    The geojsons run the full length of the highway, well beyond the 90
    domains at both ends. Unclipped, a figure shows road where there is no
    model domain, which reads as coverage that does not exist.
    """
    out = {}
    for yr, f in ROAD_FILES.items():
        if not f.exists():
            print(f"  WARNING: {f} missing - {yr} road not drawn")
            continue
        g = gpd.read_file(f)
        if g.crs is not None and dst_crs is not None:
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
    """BOTH casings first, then both lines in ROAD_ORDER so the dashed 1984
    lands on top of the solid 2004 rather than under it.

    Casings-then-lines, not casing-line-casing-line: the two alignments are
    nearly coincident for most of the island, and the second casing then
    painted out the first line, so 2004 disappeared wherever it mattered.
    `scale` thins the lines for the island figure, where the same widths would
    smother the island."""
    for yr in ROAD_ORDER:
        if yr not in roads:
            continue
        cas = dict(ROAD_CASING[yr]); cas["linewidth"] *= scale
        roads[yr].plot(ax=ax, linestyle="-", alpha=0.9, zorder=6, **cas)
    for yr in ROAD_ORDER:
        if yr not in roads:
            continue
        st = dict(ROAD_STYLE[yr]); st["linewidth"] *= scale
        roads[yr].plot(ax=ax, zorder=8 + ROAD_ORDER.index(yr), **st)


def road_legend_handles(roads):
    """Both alignments in their map colours. They are the house vintage pair,
    so neither is white and the swatches carry straight over."""
    return [Line2D([], [], label=f"NC-12 {y}", **ROAD_STYLE[y])
            for y in ROAD_ORDER if y in roads]


# =============================================================================
# FIGURES
# =============================================================================

def fig_island(elev, surv, extent, gdf, roads):
    only09 = np.where(surv == SURVEY_2009, elev, np.nan)
    vmin, vmax = elev_limits(elev)
    n09 = int((surv == SURVEY_2009).sum())
    n96 = int((surv == SURVEY_1996).sum())
    n14 = int((surv == SURVEY_2014).sum())

    # SAME CANVAS AS HAT_plot_gapfill.py's island figure, deliberately: the two
    # products are read side by side, so they have to share panel proportions
    # and gaps. The island spans ~9 km east-west and ~47 km north-south at equal
    # aspect, so panel width follows figure HEIGHT, not the width asked for.
    # At the house double-column width (190 mm) a full page of height gives
    # three ~40 mm panels, which is the whole strip at one look; the figure is
    # no longer 13 x 19 in reduced to a page, where the type became 4 pt.
    # A shade under FIG_H_MAX: the legend sits outside the axes and
    # bbox_inches="tight" adds it after layout, so the saved page is ~0.1 in
    # taller than the canvas asked for.
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=9.15),
                             sharex=True, sharey=True, constrained_layout=True)
    # Wording parallels the gapfill figure's panel titles. Cell counts are in
    # the caption, not on the canvas.
    im = panel_elev(axes[0], 0, only09, extent, vmin, vmax, "2009 survey")
    panel_elev(axes[1], 1, elev, extent, vmin, vmax, "1984-start surface")
    panel_survey(axes[2], 2, surv, extent, "survey source")

    for ax in axes:
        gdf.boundary.plot(ax=ax, color="white", linewidth=0.35, zorder=5)
        draw_roads(ax, roads, scale=0.45)
        ax.set_xlim(extent[0] - ISLAND_PAD_M, extent[1] + ISLAND_PAD_M)
        ax.set_ylim(extent[2] - ISLAND_PAD_M, extent[3] + ISLAND_PAD_M)
        ax.set_aspect("equal")
        km_axes(ax)
        spines_for_image(ax)
        ax.set_xlabel("Easting (km)")
    axes[0].set_ylabel("Northing (km)")

    cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                      fraction=0.028, pad=0.01, aspect=45)
    cb.set_label("Elevation (m NAVD88)")

    rh = road_legend_handles(roads)
    # Concise survey labels: the panels are ~40 mm wide, and the full wording is
    # wider than the axes it would sit in. The sources are named in the caption.
    fig.legend(handles=survey_legend(concise=True) + rh,
               loc="outside lower center", ncol=6, frameon=False)
    caption(fig, "The 1984-start topography, every domain on one 10 m grid. "
                 f"(a) the 2009 survey alone, {n09:,} cells; (b) the surface "
                 f"the extractor reads, {n09 + n96 + n14:,} cells "
                 f"(+{n96 + n14:,}); (c) the survey each cell came from. "
                 + counts_note(surv) + ". Domain 1 is at Cape Point in the "
                 "south, domain 90 at Pea Island in the north. "
                 + SURVEY_KEY + " " + METHOD)
    out = FIG_DIR / "HAT_mosaic1984_island.png"
    save(fig, out, vector=False, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_zoom(elev, surv, extent, gdf, roads, dom_ids, slug, cap):
    sel = gdf[gdf["domain_id"].astype(int).isin(dom_ids)]
    if sel.empty:
        print(f"  no domains {dom_ids} in the domain file - {slug} skipped")
        return None
    zx0, zy0, zx1, zy1 = sel.total_bounds
    pad = 150.0
    zw, zh = (zx1 + pad) - (zx0 - pad), (zy1 + pad) - (zy0 - pad)
    fh = (ZOOM_FIG_W / 3.0) / (zw / zh) + ZOOM_CHROME_IN

    only09 = np.where(surv == SURVEY_2009, elev, np.nan)
    vmin, vmax = elev_limits(elev)

    # sharey: all three panels show the same extent, so repeating the northing
    # labels three times only narrows the maps.
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=fh),
                             sharex=True, sharey=True, constrained_layout=True)
    im = panel_elev(axes[0], 0, only09, extent, vmin, vmax, "2009 survey")
    panel_elev(axes[1], 1, elev, extent, vmin, vmax, "1984-start surface")
    panel_survey(axes[2], 2, surv, extent, "survey source")

    for ax in axes:
        gdf.boundary.plot(ax=ax, color="white", linewidth=0.9, zorder=5)
        draw_roads(ax, roads)
        ax.set_xlim(zx0 - pad, zx1 + pad)
        ax.set_ylim(zy0 - pad, zy1 + pad)
        ax.set_aspect("equal")
        km_axes(ax, nx=3, ny=5)
        spines_for_image(ax)
        ax.set_xlabel("Easting (km)")
    axes[0].set_ylabel("Northing (km)")

    cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                      fraction=0.045, pad=0.02, aspect=40)
    cb.set_label("Elevation (m NAVD88)")

    fig.legend(handles=survey_legend(concise=True) + road_legend_handles(roads),
               loc="outside lower center", ncol=6, frameon=False)

    caption(fig, cap + " " + SURVEY_KEY + " " + METHOD)
    out = FIG_DIR / f"HAT_mosaic1984_{slug}.png"
    save(fig, out, vector=False, bbox_inches="tight")
    plt.close(fig)
    return out


def fig_shift():
    """
    Per-domain movement of the extractor's beach start.

    A diverging encoding, because the quantity has a meaningful zero and a
    meaningful sign: seaward carries the 1996 orange, landward the 2009 blue,
    so the bar colour says which survey won that domain. Zero is a neutral
    rule, not a hue.

    The +/-10 m band is drawn because it is the only threshold that matters -
    the extractor works at 10 m, so anything inside that band cannot move a
    Barrier3D cell by more than one, and a reader should not have to do the
    arithmetic to see which bars clear it.
    """
    if not AUDIT_1M.exists():
        print(f"  {AUDIT_1M} missing - shift figure skipped")
        return None
    import csv
    rows = list(csv.DictReader(open(AUDIT_1M, newline="")))
    dom = [int(r["domain"]) for r in rows]
    shift = [float(r["start_beach_shift_m"]) for r in rows]
    has_road = [r["road_line"] == "True" for r in rows]

    fig, ax = plt.subplots(figsize=figsize("double", aspect=0.44),
                           constrained_layout=True)
    # The +/- one-cell band is a reference threshold, so it takes the house
    # reference green rather than a second grey: the village bands at the top
    # of the panel are grey, and two greys on one panel cannot be told apart.
    ax.axhspan(-CELL_M, CELL_M, color=C["REF"], alpha=0.08, lw=0, zorder=0)

    # Domains 1-7 get no 1996 and so have a shift of exactly zero - a bar of
    # zero height, which is invisible. A legend swatch pointing at nothing
    # reads as "these are missing from the chart", so the span is shaded and
    # labelled in place instead.
    no_road = [d for d, hr in zip(dom, has_road) if not hr]

    cols = [C_1996 if s > 0 else C_2009 for s in shift]
    ax.bar(dom, shift, color=cols, width=0.8, zorder=3)
    ax.axhline(0, color=INK, lw=0.8, zorder=4)
    for y in (-CELL_M, CELL_M):
        ax.axhline(y, color=C["REF"], lw=0.8, ls=(0, (3, 2)), zorder=4)

    n_sea = sum(1 for s, hr in zip(shift, has_road) if hr and s >= CELL_M)
    n_land = sum(1 for s, hr in zip(shift, has_road) if hr and s <= -CELL_M)
    ax.set_xlabel(DOMAIN_AXIS_LABEL)
    ax.set_ylabel("beach start shift (m)\n+ seaward")
    ax.set_xlim(0, 91)
    # Explicit limits rather than margins(): the shifts run -5 to +72, and a
    # symmetric margin then opened a band of empty axes below -20.
    lo, hi = min(min(shift), -CELL_M), max(shift)
    ax.set_ylim(lo - 0.12 * (hi - lo), hi + 0.24 * (hi - lo))
    open_frame(ax)
    ax.grid(axis="y")
    ax.set_axisbelow(True)

    if no_road:
        ax.axvspan(min(no_road) - 0.5, max(no_road) + 0.5,
                   color=C["BASE_FILL"], zorder=1)
        ax.text((min(no_road) + max(no_road)) / 2, 0.78,
                f"{min(no_road)}-{max(no_road)}\nno 1996",
                transform=ax.get_xaxis_transform(), ha="center", va="center",
                fontsize=7.5, color=INK_MUTED, zorder=4)

    # Villages as a strip against the top edge, so they cannot be confused with
    # the shaded domain span below them. After set_xlim, as the helper
    # requires, and after the span so Buxton's name is not painted over.
    town_bands(ax, strip=0.085)

    fig.legend(handles=[
        Patch(facecolor=C_1996, label="seaward: 1996 adds beach above 0.50 m MHW"),
        Patch(facecolor=C_2009, label="landward"),
        Patch(facecolor=C["REF"], alpha=0.25,
              label="within one 10 m cell")],
        loc="outside lower center", ncol=3, frameon=False)
    caption(fig, "Where the 1984-start surface moves each domain's cross-shore "
                 f"window origin. {n_sea} domains move at least one 10 m "
                 f"Barrier3D cell seaward and {n_land} move one landward; "
                 "shifts inside the band are below the model's resolution. "
                 "Bar colour is the survey that won the beach start: 1996 "
                 "orange for seaward, 2009 blue for landward. Domains "
                 f"{min(no_road)}-{max(no_road)} have no 1984 NC-12 line and "
                 "take no 1996, so their shift is exactly zero (shaded). "
                 "Domain 1 is at Cape Point in the south, domain 90 at Pea "
                 "Island in the north; village spans are shaded along the top."
            if no_road else
            "Where the 1984-start surface moves each domain's cross-shore "
            f"window origin. {n_sea} domains move at least one 10 m Barrier3D "
            f"cell seaward and {n_land} move one landward. Domain 1 is at Cape "
            "Point in the south, domain 90 at Pea Island in the north.")
    out = FIG_DIR / "HAT_mosaic1984_shift.png"
    save(fig, out, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    elev, surv, extent, n = load_mosaic()
    print(f"mosaic {elev.shape} from {n} domains")
    print(f"  {counts_note(surv)}")

    gdf = gpd.read_file(DOMAIN_FILE)
    dom_union = (gdf.union_all() if hasattr(gdf, "union_all")
                 else gdf.unary_union)
    roads = load_roads(gdf.crs, clip_to=dom_union)

    outs = [fig_island(elev, surv, extent, gdf, roads)]
    for _ids, _slug, _cap in ZOOMS:
        outs.append(fig_zoom(elev, surv, extent, gdf, roads,
                             _ids, _slug, _cap))
    outs.append(fig_shift())
    for out in outs:
        if out:
            print(f"  wrote {out}")


if __name__ == "__main__":
    main()
