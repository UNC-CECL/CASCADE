"""
Review figures for the 1984-start DEM: what each of the three surveys contributed, and how far beach start moved.

    python scripts/input_prep/0-elevation/3-figures/HAT_plot_1984_mosaic.py

Writes the island, zoom and road figures and the beach-start shift figure to
data/hatteras_init/0-elevation/2009-2014-1996/figures/. Details: scripts/input_prep/0-elevation/README.md.

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
INIT_ROOT = PROJECT_ROOT / "data" / "hatteras_init"
import sys as _elsys
from pathlib import Path as _ELP
_elsys.path.insert(0, str(next(_q for _q in _ELP(__file__).resolve().parents
                               if (_q / "pyproject.toml").exists()) / "scripts"))
from site_layer import hat_elevation_products as _el  # noqa: E402
ELEVATION_DIR = _el.ELEVATION_ROOT

SOURCE_TAG = "2009-2014-1996"
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from site_layer.hat_elevation_products import product as _product  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    apply_style, figsize, save, caption, C, C_1984, C_1997, INK, INK_MUTED,
    DOMAIN_AXIS_LABEL, town_bands, open_frame, spines_for_image, _title)

apply_style()

_P = _product(SOURCE_TAG)
IN_DIR = _P.resampled_10m
AUDIT_1M = _P.audit_1m
FIG_DIR = _P.figures

# The repository copy of D:/Hatteras_GIS/domains.geojson (identical; 2026-09-18).
from site_layer.hat_map_layers import DOMAIN_BOXES as DOMAIN_FILE  # noqa: E402
# NC-12 alignments, reprojected from NC State Plane feet; both drawn, 1984 dashed over 2004
from site_layer.hat_topo_version import ROAD_LINE_ROOT as ROAD_DIR  # noqa: E402
# --- CONFIG ------------------------------------------------------------------
# Keyed by period; the files are the 1978 and 2008 lines those periods read
ROAD_FILES = {1984: ROAD_DIR / "1978" / "nc12_1978.geojson",
              2004: ROAD_DIR / "2008" / "nc12_2008.geojson"}
ROAD_STYLE = {2004: dict(color=C_1997, linestyle="-", linewidth=1.5),
              1984: dict(color=C_1984, linestyle=(0, (3.6, 2.4)), linewidth=1.5)}
ROAD_CASING = {2004: dict(color="white", linewidth=3.0),
               1984: dict(color="white", linewidth=3.0)}
ROAD_ORDER = [2004, 1984]   # draw order: solid first, dashed on top

GRID = 10.0
ISLAND_PAD_M = 700.0

# Chrome geometry, matching HAT_plot_gapfill.py exactly
ZOOM_FIG_W = figsize("double")[0]
ZOOM_CHROME_IN = 1.05
# Upper left: the one reliably empty corner on these zooms
ZOOM_LEGEND_LOC = "upper left"

# (domain ids, filename slug, title); both NC-12 alignments on every zoom
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
    # Subtitle numbers are from mosaic_1984_audit.csv, not from eyeballing the map
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

# One Barrier3D cell, drawn on the shift figure
CELL_M = 10.0
# -----------------------------------------------------------------------------


# Load

# Every domain placed on one 10 m grid covering the island
def load_mosaic():
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
        # Domain boxes are 505 m apart on a 500 m extent, so they overlap
        take = np.isnan(sub_e) & ~np.isnan(a)
        sub_e[take] = a[take]
        sub_s[take] = sv[take]
        fresh = (sub_s == 0) & (sv != 0)
        sub_s[fresh] = sv[fresh]

    return elev, surv, [minx, maxx, miny, maxy], len(paths)


# Draw

# UTM eastings are 6-digit metres and collide at panel width
def km_axes(ax, nx=3, ny=6):
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


# vmax from a percentile, vmin DERIVED so 0 m lands on terrain's internal water/land break
def elev_limits(elev):
    v = elev[np.isfinite(elev)]
    vmax = float(np.percentile(v, ELEV_PCT_HI))
    vmin = SEA_LEVEL_M - (vmax - SEA_LEVEL_M) * (
        TERRAIN_WATER_FRAC / (1 - TERRAIN_WATER_FRAC))
    return vmin, vmax


# One elevation panel
def panel_elev(ax, i, arr, extent, vmin, vmax, title):
    ax.set_facecolor(C_NONE)
    im = ax.imshow(arr, extent=extent, origin="upper", cmap=ELEV_CMAP,
                   vmin=vmin, vmax=vmax, interpolation="nearest", zorder=1)
    _title(ax, i, title)
    return im


# Categorical provenance panel; codes mapped to contiguous indices so colours cannot slide
def panel_survey(ax, i, surv, extent, title):
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


# Handles for the three surveys plus the unsurveyed background
def survey_legend(concise=False):
    if concise:
        return [Patch(facecolor=C_1996, label="1996"),
                Patch(facecolor=C_2009, label="2009"),
                Patch(facecolor=C_2014, label="2014"),
                Patch(facecolor=C_NONE, label="no survey")]
    return [Patch(facecolor=C_1996, label="1996 ALACE (no road boundary)"),
            Patch(facecolor=C_2009, label="2009 USACE, measured"),
            Patch(facecolor=C_2014, label="2014 NOAA Post-Sandy, gap fill"),
            Patch(facecolor=C_NONE, label="never surveyed")]


# The 10 m cell count and share for each survey
def counts_note(surv):
    n = {c: int((surv == c).sum())
         for c in (SURVEY_1996, SURVEY_2009, SURVEY_2014)}
    tot = sum(n.values())
    return (f"10 m cells: 1996 {n[SURVEY_1996]:,} ({100 * n[SURVEY_1996] / tot:.1f}%), "
            f"2009 {n[SURVEY_2009]:,} ({100 * n[SURVEY_2009] / tot:.1f}%), "
            f"2014 {n[SURVEY_2014]:,} ({100 * n[SURVEY_2014] / tot:.1f}%)")


# The method paragraph, appended to each caption (nothing on the canvas)
METHOD = ("1996 is admitted wherever it has data — there is no road boundary, "
          "and the landward limit is the ALACE swath edge: above -2.64 m "
          "NAVD88 where no other survey saw the cell and above MHW where one "
          "did, below a 12 m ceiling, contiguous with the island. No bias "
          "correction and no feathering: every cell is its own survey's "
          "measurement, unchanged. Axes are UTM eastings and northings in km; "
          "domain outlines are white.")

# The legend swatches carry only the year, because the panels are too narrow for the full wording
SURVEY_KEY = ("Survey sources: 1996 ALACE lidar (admitted with no road "
              "boundary), 2009 USACE, 2014 NOAA Post-Sandy gap fill; grey is "
              "never surveyed.")


# Both NC-12 alignments, reprojected and CLIPPED to the domain footprint
def load_roads(dst_crs, clip_to=None):
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


# BOTH casings first, then both lines in ROAD_ORDER so the dashed 1984 lands on top of the solid 2004 ...
def draw_roads(ax, roads, scale=1.0):
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


# Both alignments in their map colours
def road_legend_handles(roads):
    return [Line2D([], [], label=f"NC-12 {y}", **ROAD_STYLE[y])
            for y in ROAD_ORDER if y in roads]


# Figures

# The whole island: 2009 alone, the mosaic, and which survey each cell came from
def fig_island(elev, surv, extent, gdf, roads):
    only09 = np.where(surv == SURVEY_2009, elev, np.nan)
    vmin, vmax = elev_limits(elev)
    n09 = int((surv == SURVEY_2009).sum())
    n96 = int((surv == SURVEY_1996).sum())
    n14 = int((surv == SURVEY_2014).sum())

    # SAME CANVAS AS HAT_plot_gapfill.py's island figure, deliberately
    fig, axes = plt.subplots(1, 3, figsize=figsize("double", height=9.15),
                             sharex=True, sharey=True, constrained_layout=True)
    # Wording parallels the gapfill figure's panel titles
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
    # Concise survey labels: the panels are ~40 mm wide
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


# One zoom on a few domains, with both NC-12 alignments
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

    # Shared y: all three panels show the same extent
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


# Per-domain movement of the extractor's beach start
def fig_shift():
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
    # The +/- one-cell band, in the house reference green
    ax.axhspan(-CELL_M, CELL_M, color=C["REF"], alpha=0.08, lw=0, zorder=0)

    # Domains 1-7 shift exactly zero, so their span is shaded and labelled in place
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
    # Explicit limits rather than margins()
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

    # Villages as a strip along the top edge, drawn after the limits and the span
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


# Run: load the mosaic and roads, draw every figure
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
