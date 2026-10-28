"""
Does the domain setup match the island? Aerial imagery above, the 1 m domain clips below, the run's BRIE shoreline on both.

    python scripts/figure_making/pipeline/2-brie-offset/brie_domains_vs_imagery.py
    python scripts/figure_making/pipeline/2-brie-offset/brie_domains_vs_imagery.py --year 1996 --sections 1

    python scripts/figure_making/pipeline/2-brie-offset/brie_domains_vs_imagery.py --sections 0

One figure per 15 domains (six cover GIS 1-90) for each start year, plus one of the whole island
(--sections 0), its cross-shore axis exaggerated so the 45 km island reads. The panels share one frame:
northing along the island (x) and distance landward of the offshore datum the transects start on (y),
the frame the offset is measured in. (a) is the imagery of the start's own CoastSat window with the
mean CoastSat shoreline and the BRIE shoreline (one step per domain) on it; (b) is each domain's 1 m
clip, before the resample to 10 m, where it was cut, with the BRIE step on top, nothing moved; (c) is
each domain's Barrier3D grid, the arrays the run reads, built on the BRIE step as Barrier3D builds it
(fixed beach, two dune rows, interior). Reality, then the DEM, then the model (Hannah, 2026-10-28).
Writes
output/figures/3-model-inputs/2-brie-offset/domains_vs_imagery/<year>/.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-28
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import geopandas as gpd  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
import yaml  # noqa: E402
from affine import Affine  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.vrt import WarpedVRT  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
import sys as _hat_sys  # noqa: E402
_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
if str(_HAT_DIR) not in _hat_sys.path:
    _hat_sys.path.insert(0, str(_HAT_DIR))
from hatteras_paths import WORKSPACE as REPO  # noqa: E402
sys.path.insert(0, str(_HAT_DIR))
from site_layer import hat_topo_version as tv  # noqa: E402
from site_layer.hat_aerial_imagery import AERIAL_ROOT, Imagery  # noqa: E402
from site_layer.hat_elevation_products import product as elevation_product  # noqa: E402
from site_layer.hat_observed_rates import mean_shoreline_geojson  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    C, DOMAIN_AXIS_LABEL, FIG_W_DOUBLE, _title, apply_style, caption, elevation_cmap,
    figure_dir, figsize, save,
)

# --- CONFIG ------------------------------------------------------------------
OUT = figure_dir("inputs", "2-brie-offset") / "domains_vs_imagery"
YEARS = (1996, 2009)
# One figure per 15 domains, the same sixths as offset_sources_on_domains.py
SECTIONS = tuple((lo, lo + 14) for lo in range(1, 91, 15))
# Plus one figure of the whole island (GIS 1-90); 45 km along by ~1 km across is a sliver at true
# scale, so its cross-shore axis is exaggerated to this panel height, and read coarser to match
FULL = (1, 90)
FULL_PANEL_H_IN = 1.6
FULL_IMAGERY_RES_M = 5.0
FULL_CLIP_STEP = 5
DST_CRS = "EPSG:3725"      # the clips, transects and domain boxes
DATUM_X = 460198.45        # the offshore datum line the transects start on (EPSG:3725)
IMAGERY_RES_M = 2.5        # the imagery is read at this, about one printed pixel at 300 dpi
CLIP_STEP = 2              # the 1 m clips drawn every 2nd cell; the figure cannot show finer
PAD_M = 150                # room seaward and landward of the domains
CELL_M = 10.0              # Barrier3D cell, m
DAM = 10.0                 # decametre -> metre
# The 1 m DEM each topography product was resampled from (as HAT_bracketed_hole_cells.py)
ELEVATION_FOR_PRODUCT = {"1996-start": "2009-2014-1996", "2009-start": "2009-2014"}
# Imagery inside each start's CoastSat window (Hannah, 2026-10-28): the 1996 USGS frames and the 2009
# Google Earth captures; the 2008 NOAA mosaic is cleaner but falls outside 2009
IMAGERY = {
    1996: {"label": "USGS-georeferenced USACE photographs of 14 Oct 1996 (Henderson, in prep.)"},
    2009: {"label": "Google Earth captures of 30 May 2009 (provider not recorded)",
           "captures": AERIAL_ROOT / "2009" / "google_earth_captures" / "raw_ge",
           "crs": "EPSG:26918"},   # NAD 1983 UTM 18N, as the captures' aux.xml names it
}
PARAMETER_FILE = tv.INIT_ROOT / "Hatteras-CASCADE-parameters.yaml"
HALO = [pe.withStroke(linewidth=2.2, foreground="white")]
# -----------------------------------------------------------------------------


# The run's offset build: per-domain BRIE shoreline, metres landward of the datum (GIS 1-90)
def load_offset(year):
    source = tv.RUN_OFFSET_SOURCE
    version = tv.offset_version(year, source)
    build = tv.offset_build_dir(year, version, source)
    raw = pd.read_csv(next(build.glob("*_raw.csv")))
    datum = raw.groupby("domain_id")["ORIG_LEN"].mean().loc[1:90].to_numpy(float)
    unpadded = pd.read_csv(tv.offset_file(year, "unpadded", source=source, version=version))
    # The model input is this minus its smallest value; refuse a build whose files disagree
    if np.abs((datum - datum.min()) - unpadded.iloc[:, 1].to_numpy(float)).max() > 1e-3:
        raise SystemExit(f"{year}: {build} raw and unpadded offsets disagree; rebuild it")
    return source, version, build, datum


# One domain's 1 m clip in m above MHW (NaN where no data) and its bounds
def load_clip(folder, gis, mhw):
    with rasterio.open(folder / f"clip_domain_{gis}_filled.tif") as src:
        z = src.read(1, masked=True).filled(np.nan).astype(np.float32) - mhw
        return z, src.bounds


# The 2009 captures: world files (.jgwx) are their georeferencing, which GDAL does not read for most
def _captures(folder, crs):
    out = []
    for jpg in sorted(glob.glob(str(folder / "*.jpg"))):
        a, d, b, e, c, f = (float(v) for v in Path(jpg).with_suffix(".jgwx").read_text().split())
        # A world file names the centre of the top-left pixel; the affine wants its corner
        t = Affine(a, b, c - a / 2 - b / 2, d, e, f - d / 2 - e / 2)
        out.append((Path(jpg), t, crs))
    return out


# (H, W, 3) uint8 of the captures over the box, north up, 0 where none; first capture wins
def read_captures(caps, x0, x1, y0, y1, res):
    W, H = int(round((x1 - x0) / res)), int(round((y1 - y0) / res))
    dst_t = Affine(res, 0.0, x0, 0.0, -res, y1)
    out = np.zeros((3, H, W), dtype=np.uint8)
    for path, t, crs in caps:
        with rasterio.open(path) as src, WarpedVRT(
                src, src_crs=crs, src_transform=t, crs=DST_CRS, transform=dst_t, width=W, height=H,
                nodata=0, resampling=Resampling.average) as vrt:
            arr = vrt.read(indexes=[1, 2, 3])
        empty = out.max(axis=0) == 0
        out[:, empty] = arr[:, empty]
    return np.moveaxis(out, 0, -1)


# A north-up map array turned into the figure frame: northing along x, landward (west) up
def to_frame(a):
    return np.swapaxes(a[::-1, ::-1], 0, 1)


# GIS numbers at the domain centres, every 5th labelled
def domain_ticks(ax, centres, lo, hi):
    gis = np.arange(lo, hi + 1)
    ax.set_xticks(centres[(gis % 5) == 0])
    ax.set_xticklabels(gis[(gis % 5) == 0])
    ax.set_xticks(centres, minor=True)
    ax.tick_params(axis="x", which="minor", length=2)


def fig_section(year, lo, hi, ctx):
    gis = np.arange(lo, hi + 1)
    clips = [ctx["clips"][g] for g in gis]
    d_brie = ctx["datum"][gis - 1]
    grids = [ctx["grids"][g] for g in gis]
    beach_m = ctx["beach_m"]
    n0, n1 = clips[0]["b"].bottom, clips[-1]["b"].top
    # Every domain box and every Barrier3D grid, with room either side, in metres landward of the datum
    y_lo = min(DATUM_X - c["b"].right for c in clips) - PAD_M
    y_hi = max(max(DATUM_X - c["b"].left for c in clips),
               max(d + beach_m + len(gr) * CELL_M for d, gr in zip(d_brie, grids))) + PAD_M
    x0, x1 = DATUM_X - y_hi, DATUM_X - y_lo
    full = (lo, hi) == FULL
    res = FULL_IMAGERY_RES_M if full else IMAGERY_RES_M
    clip_step = FULL_CLIP_STEP if full else CLIP_STEP

    if "captures" in ctx:
        img = read_captures(ctx["captures"], x0, x1, n0, n1, res)
    else:
        img = ctx["imagery"].read(x0, x1, n0, n1, res, DST_CRS)
    rgba = np.dstack([img, np.where(img.max(axis=2) > 0, 255, 0).astype(np.uint8)])
    coverage = float((img.max(axis=2) > 0).mean())

    # Three panels stacked, one colour bar beside (b) and (c); at true scale unless the whole island
    pw = FIG_W_DOUBLE * 0.84
    ph = pw * (y_hi - y_lo) / (n1 - n0)
    vex = max(1.0, FULL_PANEL_H_IN / ph) if full else 1.0
    ph *= vex
    fig, axes = plt.subplots(3, 1, figsize=figsize("double", height=3 * ph + 1.4), sharex=True,
                             sharey=True, constrained_layout=True)
    ax_a, ax_b, ax_c = axes
    centres = np.array([(c["b"].bottom + c["b"].top) / 2 for c in clips])
    edges = np.r_[[c["b"].bottom for c in clips], clips[-1]["b"].top]

    # (a) the imagery, the domain boxes as clipped, the CoastSat mean shoreline and the BRIE step
    ax_a.set_facecolor("0.85")
    ax_a.imshow(to_frame(rgba), extent=(n0, n1, y_lo, y_hi), origin="lower", interpolation="bilinear",
                zorder=0)
    for c in clips:
        b = c["b"]
        ax_a.add_patch(plt.Rectangle((b.bottom, DATUM_X - b.right), b.top - b.bottom, b.right - b.left,
                                     fill=False, ec="white", lw=0.5, alpha=0.7, zorder=2))
    for line in ctx["shoreline"].geometry:
        for part in getattr(line, "geoms", [line]):
            x, y = part.xy
            # Above the BRIE step (zorder 4) and thinner, with no halo, so both show where they coincide
            ax_a.plot(y, DATUM_X - np.asarray(x), color=C["LATE"], lw=0.5, zorder=5)

    # (b) each clip where it was cut; the cmap's under colour is the water
    cmap, norm, bounds = elevation_cmap()
    ax_b.set_facecolor("0.92")
    for c in clips:
        b = c["b"]
        z = c["z"][::clip_step, ::clip_step]
        ax_b.imshow(to_frame(np.ma.masked_invalid(z)), cmap=cmap, norm=norm, origin="lower",
                    interpolation="nearest", aspect="auto", zorder=0,
                    extent=(b.bottom, b.top, DATUM_X - b.right, DATUM_X - b.left))

    # (c) each domain's Barrier3D grid as the model builds it: BRIE shoreline, the fixed beach at berm
    # height, the two dune rows, then the extracted interior; columns run south to north, row 0 is ocean
    ax_c.set_facecolor("0.92")
    beach_col = cmap(norm(ctx["berm_m"]))
    for c, d, gr in zip(clips, d_brie, grids):
        b = c["b"]
        ax_c.fill_between([b.bottom, b.top], d, d + beach_m, color=beach_col, lw=0, zorder=0)
        y0 = d + beach_m
        ax_c.imshow(gr, cmap=cmap, norm=norm, origin="lower", interpolation="nearest", aspect="auto",
                    zorder=0, extent=(b.bottom, b.top, y0, y0 + len(gr) * CELL_M))

    for ax in (ax_b, ax_c):
        for edge in edges[1:-1]:
            ax.axvline(edge, color="white", lw=0.3 if full else 0.5, ls=(0, (2, 3)), alpha=0.7, zorder=2)

    # The BRIE shoreline, a step at each domain's mean CoastSat station, in every panel
    for ax in axes:
        ax.stairs(d_brie, edges, baseline=None, color=C["ACCENT"], lw=1.0 if full else 1.3, zorder=4,
                  path_effects=HALO)
        ax.set_aspect(vex)
        ax.set_xlim(n0, n1)
        ax.set_ylim(y_lo, y_hi)
        ax.set_ylabel("Landward of datum (m)")
        domain_ticks(ax, centres, lo, hi)
    for ax in (ax_a, ax_b):
        ax.tick_params(axis="x", which="both", labelbottom=False)
    ax_c.set_xlabel(DOMAIN_AXIS_LABEL)
    w0, w1 = (pd.Timestamp(t) for t in ctx["window"])
    _title(ax_a, 0, f"{ctx['imagery_label'].split(' (')[0]}")
    _title(ax_b, 1, "1 m domain clips, before the resample, with the BRIE shoreline")
    _title(ax_c, 2, "Barrier3D grids (10 m) as the model builds them on the BRIE shoreline")
    ax_a.legend(handles=[
        Line2D([], [], color=C["LATE"], lw=0.5, label=f"Mean CoastSat shoreline, {w0:%b %Y}–{w1:%b %Y}"),
        Line2D([], [], color=C["ACCENT"], lw=1.3, label="BRIE shoreline (one step per domain)"),
        plt.Rectangle((0, 0), 1, 1, fill=False, ec="0.4", lw=0.6, label="Domain box, as clipped")],
        loc="upper right", fontsize=6.5, framealpha=0.9)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    # The same bar beside (b) and beside (c) (Hannah, 2026-10-28); a hidden one beside (a) keeps the
    # three panels one width
    for ax in axes:
        cb = fig.colorbar(sm, ax=ax, ticks=bounds[1:-1], extend="min", fraction=0.025, pad=0.01)
        cb.set_label("Elevation (m above MHW)")
        cb.outline.set_linewidth(0.5)
        if ax is ax_a:
            cb.ax.set_visible(False)

    if full:
        scope = (f"The whole island, GIS 1–90 of the {year} start (the {len(SECTIONS)} section figures show "
                 f"it at true scale): does the domain setup stand where the island is? All three panels are "
                 f"cross-shore exaggerated {vex:.0f}× in the frame the BRIE offset is measured in")
    else:
        scope = (f"GIS {lo}–{hi} of the {year} start, one of {len(SECTIONS)} figures covering GIS 1–90: does "
                 f"the domain setup stand where the island is? All three panels are at true scale in the "
                 f"frame the BRIE offset is measured in")
    caption(fig, (
        f"{scope}, northing along the island (shown as GIS domain, south to north) and metres "
        f"landward of the offshore datum line the 100 m transects run due west from; ocean is down. (a) "
        f"{ctx['imagery_label']}, inside the start's CoastSat window"
        f"{'' if coverage > 0.98 else f' (photographs cover {coverage:.0%} of the panel; grey is none)'}; "
        f"white boxes are the 2 km × 500 m domain clips where they were cut. Blue is the mean CoastSat "
        f"shoreline of {w0:%d %b %Y}–{w1:%d %b %Y}; purple is the BRIE shoreline the run receives "
        f"({ctx['source']}/{ctx['version']}), each domain's mean of its ~5 transect stations on that line, "
        f"drawn as a flat step because BRIE treats each 500 m domain as straight. (b) The {ctx['dem']} 1 m "
        f"DEM clip of every domain (the {ctx['product']} topography before its resample to 10 m), drawn "
        f"where it was cut, nothing moved, in m above MHW ({ctx['mhw']:.2f} m NAVD88), with the same BRIE "
        f"step (purple). (c) What the model reads: each domain's Barrier3D grid ({ctx['product']}, "
        f"dune-topo {ctx['dune_version']}) built on the BRIE shoreline as Barrier3D builds it, the fixed "
        f"{beach_m} m initial beach at berm height (int(BermEl / beta) cells, the same for every domain), "
        f"then the two dune rows (dune height + berm) and the extracted 10 m interior, landward up. "
        f"Where the step runs along the beach in (a) and (b), the domains and the BRIE offset stand where "
        f"the island is; an oblique coast would show as a step cutting across the shoreline in (a). "
        f"Differences between (b) and (c) are what the extraction changes: the fixed beach replacing the "
        f"surveyed one, each profile straightened against its beach, and water deeper than 3 m below MHW "
        f"trimmed from the back (shallower sound is kept, so some grids reach well into it). "
        f"The model only receives the offsets up to one constant (zeroed on their smallest); pinning them "
        f"to the datum here is what puts the step on the real coast."))
    # One subfolder per start year (Hannah, 2026-10-28)
    out = save(fig, OUT / str(year) / f"brie_domains_vs_imagery_{year}_GIS{lo:02d}-{hi:02d}", vector=False,
               close=True)
    return out[0]


def load_year(year):
    source, version, build, datum = load_offset(year)
    product = tv.product_for_year(year)
    dem = ELEVATION_FOR_PRODUCT[product]
    folder = elevation_product(dem).gapfill_1m
    with open(PARAMETER_FILE) as f:
        p = yaml.safe_load(f)
    mhw = float(p["MHW"])
    berm_m = float(p["BermEl"]) - mhw
    # Barrier3D's initial beach width: int(BermEl / beta) in dam, the same for every domain
    beach_m = int(berm_m / DAM / float(p["beta"])) * 10
    clips = {}
    for g in range(1, 91):
        z, b = load_clip(folder, g, mhw)
        clips[g] = {"z": z, "b": b}
    # The arrays the run reads, in dam: the two dune rows (above the berm) then the interior, ocean first
    elev, dunes = tv.domain_arrays(product, first_gis=1, last_gis=90)
    grids = {g: np.vstack([np.load(dn) * DAM + berm_m] * 2 + [np.load(e) * DAM])
             for g, e, dn in zip(range(1, 91), elev, dunes)}
    window = tv.shoreline_window_for_year(year)
    shoreline = gpd.read_file(mean_shoreline_geojson(*window)).to_crs(DST_CRS)
    ctx = {"source": source, "version": version, "build": build, "datum": datum, "product": product,
           "dem": dem, "mhw": mhw, "clips": clips, "window": window, "shoreline": shoreline,
           "grids": grids, "beach_m": beach_m, "berm_m": berm_m, "dune_version": tv.topo_dirs(product)[2],
           "imagery_label": IMAGERY[year]["label"]}
    if "captures" in IMAGERY[year]:
        ctx["captures"] = _captures(IMAGERY[year]["captures"], IMAGERY[year]["crs"])
        print(f"  {year}: {len(ctx['captures'])} Google Earth captures")
    else:
        ctx["imagery"] = Imagery(year, DST_CRS)
    return ctx


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--year", type=int, nargs="+", default=list(YEARS))
    ap.add_argument("--sections", type=int, nargs="+", default=None,
                    help="1-6, south to north, 0 the whole island (default: all six and the whole island)")
    args = ap.parse_args(argv)
    apply_style()
    sections = ([FULL if i == 0 else SECTIONS[i - 1] for i in args.sections] if args.sections
                else SECTIONS + (FULL,))
    for year in args.year:
        ctx = load_year(year)
        print(f"{year}: offset {ctx['source']}/{ctx['version']}, clips {ctx['dem']}, MHW {ctx['mhw']} m")
        for lo, hi in sections:
            path = fig_section(year, lo, hi, ctx)
            print(f"  GIS {lo:2d}-{hi:2d}  ->  {path.relative_to(REPO)}")


if __name__ == "__main__":
    main()
