"""
Review figure for a DEM gap fill: 2009 alone, what the fill adds, and which cells came from where.

    python scripts/input_prep/0-elevation/3-figures/HAT_plot_gapfill.py
    python scripts/input_prep/0-elevation/3-figures/HAT_plot_gapfill.py --source <PRODUCT>

Writes the island, domains 78-80 and road-overlay figures to
data/hatteras_init/0-elevation/<product>/figures/. Details: scripts/input_prep/0-elevation/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
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

# Walk up until a directory holds data/hatteras_init
def _find_project_root(start: Path) -> Path:
    for p in [start, *start.parents]:
        if (p / "data" / "hatteras_init").is_dir():
            return p
    raise SystemExit(f"cannot find data/hatteras_init above {start}")


PROJECT_ROOT = _find_project_root(Path(__file__).resolve())

# The house style, before anything is drawn or any colour is named.
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, figsize, save, caption, C_1984, C_1997, spines_for_image,
    _title)

apply_style()

import sys as _elsys
from pathlib import Path as _ELP
_elsys.path.insert(0, str(next(_q for _q in _ELP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_elevation_products as _el  # noqa: E402
ELEVATION_DIR = _el.ELEVATION_ROOT
IN_DIR = None   # set from SOURCE_TAG below
FIG_DIR = None  # set from SOURCE_TAG below (was the pooled 0-elevation/figures/, gone since 08-25)
# The repository copy of D:/Hatteras_GIS/domains.geojson (identical; 2026-09-18).
from site_layer.hat_map_layers import DOMAIN_BOXES as DOMAIN_FILE  # noqa: E402

# NC-12 alignments, reprojected from NC State Plane feet on load
from site_layer.hat_topo_version import ROAD_LINE_ROOT as ROAD_DIR  # noqa: E402
# Keyed by period; the files are the 1978 and 2008 lines those periods read
ROAD_FILES = {1984: ROAD_DIR / "1978" / "nc12_1978.geojson",
              2004: ROAD_DIR / "2008" / "nc12_2008.geojson"}
# Two vintages of the same line, so they take the house vintage pair
ROAD_STYLE = {2004: dict(color=C_1997, linestyle="-", linewidth=1.5),
              1984: dict(color=C_1984, linestyle=(0, (3.6, 2.4)), linewidth=1.5)}
ROAD_CASING = {2004: dict(color="white", linewidth=3.0),
               1984: dict(color="white", linewidth=3.0)}
ROAD_ORDER = [2004, 1984]   # draw order: solid first, dashed on top

# Fill sources this script knows how to plot
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

# Built from SURVEY_FILL so it cannot disagree with the data being plotted
SOURCE_NOTE = (f"Fill is limited to cells the {SURVEY_FILL} source measured, "
               f"contiguous with the island (20 m bridging), above -2.64 m NAVD88. "
               f"No bias correction, no feathering — filled cells are the "
               f"{SURVEY_FILL} measurement unchanged. Fill source: "
               f"{SOURCE_LONG}. Axes are UTM eastings and northings in km; "
               f"domain outlines are white.")

# The superseded fallback that used to live here is gone
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_elevation_products import product as _product  # noqa: E402

_P = _product(SOURCE_TAG)
IN_DIR = _P.resampled_10m
FIG_DIR = _P.figures

# --- CONFIG ------------------------------------------------------------------
GRID = 10.0

# Breathing room around the island-wide mosaic
ISLAND_PAD_M = 700.0

# Zoom figure geometry: height from the data aspect, legends outside the axes

# Road-overlay zooms: (domain ids, NC-12 years, filename slug, caption sentence)
ROAD_ZOOMS = [
    ([78, 79, 80], [2004, 1984], "roads_78_80",
     "Domains 78-80, the roadways the extractor names as width-drowning at "
     "t=0 on missing survey coverage, with both NC-12 alignments drawn: "
     "1984 dashed red, 2004 solid blue."),
    (list(range(8, 16)), [2004, 1984], "roads_8_15",
     "Domains 8-15, the southern end, with both NC-12 alignments drawn: "
     "1984 dashed red, 2004 solid blue."),
    (list(range(82, 89)), [2004, 1984], "roads_82_88",
     "Domains 82-88, the northern reach, with both NC-12 alignments drawn: "
     "1984 dashed red, 2004 solid blue."),
]

# The house double-column width (190 mm) since 2026-09-10
ZOOM_FIG_W = figsize("double")[0]
ZOOM_CHROME_IN = 1.05
ID_RE = re.compile(r"resampled_domain_(\w+)_filled\.tif$")

SURVEY_2009, SURVEY_NONE = 2009, 0

# Survey-source colours sampled from terrain, spaced by luminance so they read in greyscale
C_2009 = "#2353b9"   # terrain's water blue
C_FILL = "#31d670"   # terrain's low-land green
C_NONE = "#E4E4E4"   # neutral grey, off the elevation ramp entirely

ELEV_CMAP = "terrain"
ELEV_PCT = (2, 98)   # clip the ramp to percentiles so a few spikes don't flatten it

# Matplotlib's `terrain` is built for topography
SEA_LEVEL_M = 0.0        # m NAVD88; use MHW (0.36) to key the break to MHW
TERRAIN_WATER_FRAC = 0.25
# -----------------------------------------------------------------------------


# Places every domain on one 10 m grid covering the island
def load_mosaic():
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
        # Overlapping boxes: keep what is already placed rather than overwrite a neighbour
        take = np.isnan(sub_e) & ~np.isnan(a)
        sub_e[take] = a[take]
        sub_s[take] = sv[take]
        fresh = (sub_s == 0) & (sv != 0)
        sub_s[fresh] = sv[fresh]

    extent = [minx, maxx, miny, maxy]
    return elev, surv, extent, len(paths)


# Domain box outlines, in white
def draw_domains(ax, gdf, lw=0.35):
    gdf.boundary.plot(ax=ax, color="white", linewidth=lw, zorder=5)


# UTM eastings here are 6-digit metres (450439..458392)
def km_axes(ax, nx=3, ny=6):
    ax.xaxis.set_major_locator(MaxNLocator(nbins=nx, prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=ny))

    # Decimal places from the span, not fixed
    def dp(span_m):
        return 0 if span_m > 20_000 else (1 if span_m > 2_000 else 2)

    xs = abs(np.diff(ax.get_xlim())[0])
    ys = abs(np.diff(ax.get_ylim())[0])
    ax.xaxis.set_major_formatter(
        FuncFormatter(lambda v, _, d=dp(xs): f"{v / 1000:.{d}f}"))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, _, d=dp(ys): f"{v / 1000:.{d}f}"))
    ax.tick_params(labelsize=8)


# Loads the NC-12 alignments, reprojects them, and CLIPS them to the domain footprint
def load_roads(dst_crs, clip_to=None):
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


# BOTH casings first, then both lines in ROAD_ORDER so the dashed 1984 lands on top of the solid 2004 ...
def draw_roads(ax, roads, scale=1.0):
    for yr in ROAD_ORDER:
        if yr not in roads:
            continue
        cas = dict(ROAD_CASING[yr]); cas["linewidth"] *= scale
        roads[yr].plot(ax=ax, linestyle="-", zorder=6, alpha=0.9, **cas)
    for yr in ROAD_ORDER:
        if yr not in roads:
            continue
        st = dict(ROAD_STYLE[yr]); st["linewidth"] *= scale
        roads[yr].plot(ax=ax, zorder=8 + ROAD_ORDER.index(yr), **st)


# Both alignments in their map colours
def road_legend_handles(roads):
    return [Line2D([], [], label=f"NC-12 {y}", **ROAD_STYLE[y])
            for y in ROAD_ORDER if y in roads]


# One elevation panel
def panel_elev(ax, i, arr, extent, vmin, vmax, title):
    ax.set_facecolor(C_NONE)
    im = ax.imshow(arr, extent=extent, origin="upper", cmap=ELEV_CMAP,
                   vmin=vmin, vmax=vmax, interpolation="nearest", zorder=1)
    _title(ax, i, title)
    return im


# Run: load the mosaic and roads, draw the island, zoom and road figures
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

    # Panel width follows figure height here (a 9 x 47 km island at equal aspect)
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=9.15),
                             sharex=True, sharey=True, constrained_layout=True)

    # Cell counts are in the caption, not on the canvas.
    im = panel_elev(axes[0], 0, only09, extent, vmin, vmax, "2009 survey")
    panel_elev(axes[1], 1, filled, extent, vmin, vmax,
               f"with {SURVEY_FILL} fill")

    # Boundaries ascend: measured 2009 before the 2014 fill
    lo, hi = sorted([SURVEY_2009, SURVEY_FILL])
    cmap_s = ListedColormap([C_NONE,
                             C_2009 if lo == SURVEY_2009 else C_FILL,
                             C_FILL if hi == SURVEY_FILL else C_2009])
    norm_s = BoundaryNorm([-0.5, 0.5, lo + 0.5, hi + 0.5], cmap_s.N)
    axes[2].set_facecolor(C_NONE)
    axes[2].imshow(surv, extent=extent, origin="upper", cmap=cmap_s, norm=norm_s,
                   interpolation="nearest", zorder=1)
    _title(axes[2], 2, "survey source")

    for ax in axes:
        draw_domains(ax, gdf)
        # Thinner roads at island scale
        draw_roads(ax, roads, scale=0.45)
        ax.set_xlim(extent[0] - ISLAND_PAD_M, extent[1] + ISLAND_PAD_M)
        ax.set_ylim(extent[2] - ISLAND_PAD_M, extent[3] + ISLAND_PAD_M)
        ax.set_xlabel("Easting (km)")
        km_axes(ax)
        ax.set_aspect("equal")
        spines_for_image(ax)
    axes[0].set_ylabel("Northing (km)")

    # Colorbar against all three axes, so none is shrunk out of line
    cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                      fraction=0.028, pad=0.01, aspect=45)
    cb.set_label("Elevation (m NAVD88)")

    fig.legend(handles=[Patch(facecolor=C_2009, label="2009 measured"),
                        Patch(facecolor=C_FILL, label=f"{SURVEY_FILL} fill"),
                        Patch(facecolor=C_NONE, label="never surveyed")]
               + road_legend_handles(roads),
               loc="outside lower center", ncol=5, frameon=False)

    caption(fig, "The 2009 DEM gap fill, every domain on one 10 m grid. "
                 f"(a) the 2009 survey alone, {n_meas:,} cells; (b) the "
                 f"surface the extractor reads, {n_meas + n_fill:,} cells "
                 f"(+{n_fill:,}); (c) the survey each cell came from. Domain 1 "
                 "is at Cape Point in the south, domain 90 at Pea Island in "
                 "the north. " + SOURCE_NOTE)
    out = FIG_DIR / f"HAT_gapfill_{SOURCE_TAG}_island.png"
    save(fig, out, vector=False, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")

    # zoom: the domains the extractor names as width-drowning at t=0
    sel = gdf[gdf["domain_id"].astype(int).isin([78, 79, 80])]
    if not sel.empty:
        zminx, zminy, zmaxx, zmaxy = sel.total_bounds
        pad = 150
        # Figure height from the zoom's aspect, so no gap opens under the title
        _zw = (zmaxx + pad) - (zminx - pad)
        _zh = (zmaxy + pad) - (zminy - pad)
        _panel_w = ZOOM_FIG_W / 3.0
        _fig_h = _panel_w / (_zw / _zh) + ZOOM_CHROME_IN
        # Shared y: all three panels show the same extent
        fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=_fig_h),
                                 sharex=True, sharey=True,
                                 constrained_layout=True)
        im = panel_elev(axes[0], 0, only09, extent, vmin, vmax, "2009 survey")
        panel_elev(axes[1], 1, filled, extent, vmin, vmax,
                   f"with {SURVEY_FILL} fill")
        axes[2].set_facecolor(C_NONE)
        axes[2].imshow(surv, extent=extent, origin="upper", cmap=cmap_s,
                       norm=norm_s, interpolation="nearest", zorder=1)
        _title(axes[2], 2, "survey source")
        for ax in axes:
            draw_domains(ax, gdf, lw=0.9)
            ax.set_xlim(zminx - pad, zmaxx + pad)
            ax.set_ylim(zminy - pad, zmaxy + pad)
            ax.set_aspect("equal")
            km_axes(ax, nx=3, ny=5)
            ax.set_xlabel("Easting (km)")
            spines_for_image(ax)
        axes[0].set_ylabel("Northing (km)")
        fig.legend(handles=[Patch(facecolor=C_2009, label="2009 measured"),
                            Patch(facecolor=C_FILL, label=f"{SURVEY_FILL} fill"),
                            Patch(facecolor=C_NONE, label="never surveyed")],
                   loc="outside lower center", ncol=3, frameon=False)
        cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                          fraction=0.045, pad=0.02, aspect=40)
        cb.set_label("Elevation (m NAVD88)")
        caption(fig, "Domains 78-80, the roadways the extractor names as "
                     "width-drowning at t=0 on missing survey coverage. "
                     "(a) the 2009 survey alone; (b) the surface the extractor "
                     "reads; (c) the survey each cell came from. " + SOURCE_NOTE)
        out2 = FIG_DIR / f"HAT_gapfill_{SOURCE_TAG}_domains_78_80.png"
        save(fig, out2, vector=False, bbox_inches="tight")
        plt.close(fig)
        print(f"wrote {out2}")

        # Road-overlay zooms, one per entry in ROAD_ZOOMS, drawn at zoom scale
        def _road_zoom(dom_ids, years, slug, cap):
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
            fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=fh),
                                     sharex=True, sharey=True,
                                     constrained_layout=True)
            im = panel_elev(axes[0], 0, only09, extent, vmin, vmax,
                            "2009 survey")
            panel_elev(axes[1], 1, filled, extent, vmin, vmax,
                       f"with {SURVEY_FILL} fill")
            axes[2].set_facecolor(C_NONE)
            axes[2].imshow(surv, extent=extent, origin="upper", cmap=cmap_s,
                           norm=norm_s, interpolation="nearest", zorder=1)
            _title(axes[2], 2, "survey source")
            for ax in axes:
                draw_domains(ax, gdf, lw=0.9)
                draw_roads(ax, rsub)
                ax.set_xlim(zx0 - pad_, zx1 + pad_)
                ax.set_ylim(zy0 - pad_, zy1 + pad_)
                ax.set_aspect("equal")
                km_axes(ax, nx=3, ny=5)
                ax.set_xlabel("Easting (km)")
                spines_for_image(ax)
            axes[0].set_ylabel("Northing (km)")
            rh = road_legend_handles(rsub)
            fig.legend(
                handles=[Patch(facecolor=C_2009, label="2009 measured"),
                         Patch(facecolor=C_FILL, label=f"{SURVEY_FILL} fill"),
                         Patch(facecolor=C_NONE, label="never surveyed")] + rh,
                loc="outside lower center", ncol=5, frameon=False)
            cb = fig.colorbar(im, ax=list(axes), orientation="horizontal",
                              fraction=0.045, pad=0.02, aspect=40)
            cb.set_label("Elevation (m NAVD88)")
            caption(fig, cap + " (a) the 2009 survey alone; (b) the surface "
                          "the extractor reads; (c) the survey each cell came "
                          "from. " + SOURCE_NOTE)
            outp = FIG_DIR / f"HAT_gapfill_{SOURCE_TAG}_{slug}.png"
            save(fig, outp, vector=False, bbox_inches="tight")
            plt.close(fig)
            print(f"wrote {outp}")

        if roads:
            for _ids, _yrs, _slug, _cap in ROAD_ZOOMS:
                _road_zoom(_ids, _yrs, _slug, _cap)


if __name__ == "__main__":
    main()
