"""
Observed cross-shore elevation profiles from every lidar survey, 1996-2019: one figure per GIS domain, plus named examples.

    python scripts/figure_making/shoreline/transect_profiles.py
    python scripts/figure_making/shoreline/transect_profiles.py --gis 80 83
    python scripts/figure_making/shoreline/transect_profiles.py --examples-only
    python scripts/figure_making/shoreline/transect_profiles.py --all-surveys      # every survey, to check

Each domain's CoastSat transects (about 10 per 500 m domain) are sampled every 1 m on every
survey HAT_survey_dune_crests.py reads, with its files, vertical units and offsets, and put
on x = 0 at that transect's CoastSat shoreline (median, Aug 2008 to Aug 2010). by_domain/ draws the median over the domain's
transects for six surveys spread over the record, over the beach and dune. Needs geopandas, rasterio,
pyproj, scipy.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-10-27
"""

from __future__ import annotations

import argparse
import glob
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pyproj  # noqa: E402
import rasterio  # noqa: E402
import rasterio.windows  # noqa: E402
from matplotlib.colors import LightSource, Normalize  # noqa: E402
from rasterio.windows import Window  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

# Hatteras/ holds the scripts; data/ and output/ live in the workspace on D: (hatteras_paths)
_HAT_DIR = next(p for p in Path(__file__).resolve().parents if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(_HAT_DIR))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "5-scr" / "lib"))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "5-scr" / "1-observations" / "mean_shoreline"))
sys.path.insert(0, str(_HAT_DIR / "input_prep" / "1-barrier3d-domains" / "2-survey-crests"))
import scr_paths  # noqa: E402,F401
import HAT_survey_dune_crests as S  # noqa: E402
from coastsat_lrr import filter_dates, load_timeseries  # noqa: E402
from coastsat_mean_shoreline import TARGET_CRS, timeseries_file, transect_geometry  # noqa: E402
from site_layer.hat_figure_style import (  # noqa: E402
    INK, INK_MUTED, apply_style, caption, elevation_cmap, figsize, figure_dir, north_dart,
    open_frame, save, scale_bar_km, spines_for_image, support_dir, title)
from site_layer.hat_observed_rates import transect_lookup  # noqa: E402

# --- CONFIG ------------------------------------------------------------------
OUT = figure_dir("observations", "transect_profiles")
LANDWARD_M = 350                   # drawn landward of each transect's shoreline
SEAWARD_M = 120                    # drawn seaward of it
ANCHOR_WINDOW = ("2008-08-17", "2010-08-17")   # x = 0: the transect's median CoastSat position, the 2009 start's
MIN_ANCHOR_OBS = 5                 # fewer CoastSat positions in the window: the transect is left out
REF_YEAR = 2009                    # the survey a raster's vertical unit is checked against
BASEMAP_YEAR = 2019                # the plan-view panel's survey
MHW_M = S.X.MHW_M                  # m NAVD88
MIN_TRANSECT_SHARE = 0.5           # the median is drawn where at least this share of transects has data
SPIKE_M = S.SPIKE_M                # above both the survey before and after at the same x: blanked
UNIT_CHECK_GIS = (10, 30, 50, 70)  # where a raster's vertical unit is read against the 2009 lidar
UNIT_SEPARATION = 2.0              # the other unit must fit at least this many times worse
CMAP = plt.cm.rainbow_r            # oldest red through to newest purple
POST_STORM = {1999: "Dennis", 2016: "Matthew", 2018: "Florence"}   # flown within weeks: dashed
SHOW_YEARS = (1997, 2001, 2005, 2009, 2014, 2019)   # the surveys drawn; --all-surveys draws every one
ZOOM_M = (-150, 100)               # the beach and dune: the plotted extent, unless an example is wide
WIDE_EXAMPLES = ("pea_island_overwash_gis87", "road_relocation_1999_gis13")   # the interior is the point: full extent
BAND_PCT = (25, 75)                # spread across transects, shaded for the first and last survey drawn
BARE_EARTH_FROM = 2009             # earlier surveys are first-return grids: buildings and trees stay in
BEHIND_DUNE_M = 120                # landward of this, those surveys are drawn faint (README)
# Datum check, reported and never applied: each survey against REF_YEAR on ground behind the dune that the
# bare-earth surveys agree has not changed
STABLE_RANGE_M = 0.2               # bare-earth surveys within this of each other: unchanged ground
STABLE_ABOVE_MHW_M = 0.3           # and this far above MHW in REF_YEAR, so not water
MIN_STABLE_POINTS = 200            # fewer: no check in that domain
DATUM_FLAG_M = 0.3                 # a survey further than this from REF_YEAR there is starred
# stem -> (GIS domain, panel title, why it is an example: the caption's last sentence)
EXAMPLES = {
    "rodanthe_gis80": (80, "Rodanthe", "The fastest-eroding domain on the island (CoastSat 1996-2025, -2.5 m/yr)."),
    "rodanthe_fill_2014_gis83": (83, "Rodanthe, 2014 fill", "Inside the 2014 Rodanthe emergency fill (GIS 82-88): "
                                 "the 2014 survey is before it, 2016 after."),
    "buxton_nourishment_gis07": (7, "Buxton, north of the groins", "North of the Buxton groin field, inside the "
                                 "2017 nourishment (placed June 2017 to February 2018): 2018 is the first survey after."),
    "avon_coastsat_hotspot_gis25": (25, "Avon", "The erosion hotspot of the CoastSat rates at the south end of "
                                    "Avon (-2.2 m/yr)."),
    "avon_mapped_hotspot_gis34": (34, "North of Avon", "Inside the Avon erosion hotspot of domains.geojson (GIS "
                                  "32-36), where CoastSat has the shoreline advancing (+1.7 m/yr)."),
    "salvo_wimble_accretion_gis67": (67, "Salvo, Wimble Shoals", "The fastest-accreting domain away from Cape "
                                     "Point (CoastSat +2.9 m/yr)."),
    "pea_island_overwash_gis87": (87, "Rodanthe-Pea Island overwash", "Inside the Rodanthe-Pea Island overwash "
                                  "hotspot (GIS 82-90), where NC-12 was moved landward in 1989 and replaced by the "
                                  "Jug Handle Bridge in 2022, after the last survey."),
    "control_gis41": (41, "Undeveloped, slow change", "A slowly changing domain outside every hotspot (CoastSat "
                      "+0.6 m/yr), with the longest survey record."),
    "road_relocation_1999_gis13": (13, "1999 NC-12 relocation", "NC-12 was moved landward here in 1999 (GIS 9-14); "
                                   "the 1996-1999 surveys are before it."),
    "s_curves_gis86": (86, "S-curves", "The Rodanthe S-curves erosion hotspot and the Rodanthe-Pea Island overwash "
                       "hotspot; the 2011 post-Irene survey is not drawn (point cloud only)."),
    "cape_point_gis03": (3, "Cape Point", "Cape Point, accreting fast (CoastSat +2.7 m/yr) and an overwash hotspot."),
}
# -----------------------------------------------------------------------------


# Every survey either start reads, oldest first: year -> file glob
def survey_files():
    out = {}
    for _, table in S.STARTS.values():
        for year, (_, pattern, *_) in table.items():
            out[year] = pattern
    return dict(sorted(out.items()))


# The vertical shift (m) each survey is given in the dune-crest work; none for the rest
def survey_offsets():
    out = {}
    for table in S.SURVEY_VERTICAL_OFFSET_M.values():
        out.update(table)
    return out


# Open each survey's files once, with a transformer from UTM 18N into each file's horizontal CRS
def open_surveys():
    surveys = {}
    for year, pattern in survey_files().items():
        files = sorted(glob.glob(str(pattern)))
        if not files:
            raise SystemExit(f"{year}: no files at {pattern}")
        srcs = []
        for f in files:
            src = rasterio.open(f)
            crs = pyproj.CRS.from_wkt((src.crs or S.ASSUMED_CRS).to_wkt())
            if crs.is_compound:
                crs = crs.sub_crs_list[0]
            srcs.append((src, pyproj.Transformer.from_crs(TARGET_CRS, crs, always_xy=True)))
        surveys[year] = srcs
    return surveys


# Elevation at UTM 18N points in the source's own unit, bilinear; NaN where no file has all four cells
def sample(srcs, x, y):
    z = np.full(x.shape, np.nan)
    for src, to_src in srcs:
        sx, sy = to_src.transform(x, y)
        col, row = ~src.transform * (sx, sy)
        col, row = col - 0.5, row - 0.5
        hit = (col >= 0) & (row >= 0) & (col < src.width - 1) & (row < src.height - 1) & np.isnan(z)
        if not hit.any():
            continue
        c0, r0 = int(col[hit].min()), int(row[hit].min())
        c1, r1 = int(col[hit].max()) + 2, int(row[hit].max()) + 2
        a = S.clean(src.read(1, window=Window(c0, r0, c1 - c0, r1 - r0), boundless=True, fill_value=np.nan),
                    src.nodata)
        z[hit] = map_coordinates(a, [row[hit] - r0, col[hit] - c0], order=1, mode="constant", cval=np.nan)
    return z


# Points every 1 m along each transect, from LANDWARD_M behind its origin to SEAWARD_M past its end
def transect_points(geom):
    out = {}
    for tid, (x0, y0, ux, uy, length) in geom.items():
        s = np.arange(-LANDWARD_M, int(length) + SEAWARD_M + 1, 1.0)
        out[tid] = (s, x0 + s * ux, y0 + s * uy)
    return out


# Every survey along every transect of one domain: year -> transect -> elevation (source unit)
def sample_domain(surveys, points):
    tids = list(points)
    x = np.concatenate([points[t][1] for t in tids])
    y = np.concatenate([points[t][2] for t in tids])
    cuts = np.cumsum([len(points[t][0]) for t in tids])[:-1]
    return {year: dict(zip(tids, np.split(sample(srcs, x, y), cuts))) for year, srcs in surveys.items()}


# Metres per source unit: stated in the CRS, else whichever of m or ft puts its land closer to the 2009 lidar's
def unit_factors(surveys, raw_by_gis):
    factors = {}
    for year, srcs in surveys.items():
        stated = {S.stated_unit(src) for src, _ in srcs} - {None}
        if len(stated) == 1:
            factors[year] = S.US_FT if stated.pop() == "ft" else 1.0
            continue
        a = np.concatenate([np.concatenate(list(r[year].values())) for r in raw_by_gis.values()])
        ref = np.concatenate([np.concatenate(list(r[REF_YEAR].values())) for r in raw_by_gis.values()])
        m = np.isfinite(a) & np.isfinite(ref) & (ref > 1.0) & (a > 0)
        if m.sum() < 200:
            raise SystemExit(f"{year}: under 200 land points shared with the {REF_YEAR} lidar on GIS {UNIT_CHECK_GIS}")
        miss = {f: float(np.median(np.abs(a[m] * f - ref[m]))) for f in (1.0, S.US_FT)}
        best, other = sorted(miss, key=miss.get)
        if miss[other] < UNIT_SEPARATION * miss[best]:
            raise SystemExit(f"{year}: cannot tell the vertical unit, median misfit {miss[1.0]:.2f} m as metres, "
                             f"{miss[S.US_FT]:.2f} m as feet, against {REF_YEAR}")
        factors[year] = best
        print(f"  {year}: vertical unit {'ft' if best != 1.0 else 'm'} (median misfit to {REF_YEAR} "
              f"{miss[best]:.2f} m, {miss[other]:.2f} m the other way)")
    return factors


# n values of arr from index start, NaN where that runs off either end
def window(arr, start, n):
    out = np.full(n, np.nan)
    lo, hi = max(start, 0), min(start + n, len(arr))
    if hi > lo:
        out[lo - start:hi - start] = arr[lo:hi]
    return out


# The transect's median CoastSat chainage (m from its origin) over ANCHOR_WINDOW; None if too few positions
def coastsat_anchor(tid):
    ts = filter_dates(load_timeseries(str(timeseries_file(tid))), *ANCHOR_WINDOW)
    return float(ts.chainage_m.median()) if len(ts) >= MIN_ANCHOR_OBS else None


# Blank points more than SPIKE_M above both the survey before and after at the same x (the dune-crest rule)
def drop_spikes(stack):
    n_year, n_x = stack.shape
    prev, nxt = np.full_like(stack, np.nan), np.full_like(stack, np.nan)
    last = np.full(n_x, np.nan)
    for k in range(n_year):
        prev[k] = last
        last = np.where(np.isfinite(stack[k]), stack[k], last)
    last = np.full(n_x, np.nan)
    for k in range(n_year - 1, -1, -1):
        nxt[k] = last
        last = np.where(np.isfinite(stack[k]), stack[k], last)
    with np.errstate(invalid="ignore"):
        up_prev = (stack - prev > SPIKE_M) | np.isnan(prev)
        up_next = (stack - nxt > SPIKE_M) | np.isnan(nxt)
    spike = up_prev & up_next & (np.isfinite(prev) | np.isfinite(nxt)) & np.isfinite(stack)
    return np.where(spike, np.nan, stack), spike.sum(axis=1)


# One domain's profiles on x: year -> (transects, x) in m NAVD88, the anchors (UTM), and a coverage table
def domain_profiles(gis, raw, points, geom, factors, offsets):
    years = list(raw)
    xs = np.arange(-LANDWARD_M, SEAWARD_M + 1, 1.0)
    rows, anchors, stacks, spike_counts = [], {}, [], []
    for tid, (s, _, _) in points.items():
        chain = coastsat_anchor(tid)
        if chain is None:
            print(f"  GIS {gis}: {tid} has under {MIN_ANCHOR_OBS} CoastSat positions in {ANCHOR_WINDOW}, left out")
            continue
        k = int(round(chain - s[0]))
        stack = np.array([window(raw[y][tid], k - LANDWARD_M, len(xs)) * factors[y] + offsets.get(y, 0.0)
                          for y in years])
        stack, spikes = drop_spikes(stack)
        stacks.append(stack)
        x0, y0, ux, uy, _ = geom[tid]
        anchors[tid] = (x0 + s[k] * ux, y0 + s[k] * uy, ux, uy)
        spike_counts.append(spikes)
    if not stacks:
        raise SystemExit(f"GIS {gis}: no transect has a CoastSat shoreline in {ANCHOR_WINDOW}")
    cube = np.stack(stacks, axis=1)                    # (year, transect, x)
    spikes = np.sum(spike_counts, axis=0)
    datum = interior_offsets(xs, dict(zip(years, cube)))
    for j, y in enumerate(years):
        rows.append(dict(gis=gis, survey=y, transects=cube.shape[1],
                         transects_with_data=int(np.isfinite(cube[j]).any(axis=1).sum()),
                         median_coverage=float(np.isfinite(median_profile(cube[j])).mean()),
                         spike_points=int(spikes[j]), **datum[y]))
    return xs, dict(zip(years, cube)), anchors, rows


# Each survey's median difference from REF_YEAR on unchanged ground behind the dune, with the point count
def interior_offsets(xs, profiles):
    bare = np.stack([b for y, b in profiles.items() if y >= BARE_EARTH_FROM])
    ref = profiles[REF_YEAR]
    with np.errstate(invalid="ignore"):
        stable = (np.isfinite(bare).all(axis=0) & (bare.max(axis=0) - bare.min(axis=0) <= STABLE_RANGE_M)
                  & (ref > MHW_M + STABLE_ABOVE_MHW_M) & (xs < -BEHIND_DUNE_M)[None, :])
    out = {}
    for y, block in profiles.items():
        d = (block - ref)[stable & np.isfinite(block)]
        ok = len(d) >= MIN_STABLE_POINTS and y != REF_YEAR
        out[y] = dict(stable_points=int(len(d)), interior_offset_m=round(float(np.median(d)), 3) if ok else np.nan)
    return out


# Median over transects at each x, where at least MIN_TRANSECT_SHARE of them have data
def median_profile(block):
    n = np.isfinite(block).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(block, axis=0)
    return np.where(n >= np.ceil(MIN_TRANSECT_SHARE * block.shape[0]), med, np.nan)


# The basemap survey around the drawn transects, m MHW, with its extent in UTM 18N
def basemap(srcs, factor, offset, bounds):
    if len(srcs) > 1:
        raise SystemExit(f"basemap survey {BASEMAP_YEAR} must be one file, not {len(srcs)}")
    src, to_src = srcs[0]
    sx, sy = to_src.transform([bounds[0], bounds[2]], [bounds[1], bounds[3]])
    win = rasterio.windows.from_bounds(sx[0], sy[0], sx[1], sy[1], transform=src.transform)
    win = win.round_offsets().round_lengths()
    a = S.clean(src.read(1, window=win, boundless=True, fill_value=np.nan), src.nodata)
    t = src.window_transform(win)
    ex, ey = to_src.transform([t.c, t.c + t.a * a.shape[1]], [t.f + t.e * a.shape[0], t.f], direction="INVERSE")
    return a * factor + offset - MHW_M, (ex[0], ex[1], ey[0], ey[1])


# (a) the basemap with the transects over the plotted extent and the shoreline they are measured from
def draw_map(ax, z, extent, anchors, xlim):
    cmap, norm, _ = elevation_cmap()
    ax.set_facecolor("0.85")
    shade = LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(z, nan=np.nanmin(z)), vert_exag=2)
    ax.imshow(shade, cmap="gray", extent=extent, interpolation="none")
    ax.imshow(np.ma.masked_invalid(z), cmap=cmap, norm=norm, alpha=0.75, extent=extent, interpolation="none")
    for ax_, ay, ux, uy in anchors.values():
        ax.plot([ax_ + xlim[0] * ux, ax_ + xlim[1] * ux], [ay + xlim[0] * uy, ay + xlim[1] * uy],
                color="white", lw=0.9, solid_capstyle="butt",
                path_effects=[pe.withStroke(linewidth=1.8, foreground="0.12")])
    a = np.array([v[:2] for v in anchors.values()])
    ax.plot(a[:, 0], a[:, 1], color=INK, lw=0.9, ls=":", zorder=5)
    ax.set(xlim=extent[:2], ylim=extent[2:], xticks=[], yticks=[], aspect="equal")
    spines_for_image(ax)
    scale_bar_km(ax, length_m=100, segments=1, unit="m")
    w, h = extent[1] - extent[0], extent[3] - extent[2]
    north_dart(ax, (extent[0] + 0.1 * w, extent[3] - 0.13 * h), arrow_m=0.08 * h)


# The BAND_PCT percentiles over transects at each x, where the median is drawn
def spread(block):
    n = np.isfinite(block).sum(axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        lo, hi = np.nanpercentile(block, BAND_PCT, axis=0)
    keep = n >= np.ceil(MIN_TRANSECT_SHARE * block.shape[0])
    return np.where(keep, lo, np.nan), np.where(keep, hi, np.nan)


# (b) the domain's median for each survey drawn, with the spread for the first and last
def draw_profiles(ax, xs, profiles, years, xlim):
    norm = Normalize(min(survey_files()), max(survey_files()))
    view = (xs >= xlim[0]) & (xs <= xlim[1])
    behind = xs < -BEHIND_DUNE_M
    drawn = [y for y in years if y in profiles and np.isfinite(median_profile(profiles[y])[view]).any()]
    zs = []
    for i, year in enumerate(drawn):
        block = profiles[year]
        med = median_profile(block)
        col, ls = CMAP(norm(year)), ("--" if year in POST_STORM else "-")
        faint = behind if year < BARE_EARTH_FROM else np.zeros_like(behind)
        if i in (0, len(drawn) - 1):
            lo, hi = spread(block)
            ax.fill_between(xs, np.where(faint, np.nan, lo), np.where(faint, np.nan, hi), color=col, alpha=0.18,
                            lw=0, zorder=1)
        ax.plot(xs, np.where(faint, np.nan, med), color=col, lw=1.3, ls=ls, zorder=3)
        ax.plot(xs, np.where(faint | np.roll(faint, 1), med, np.nan), color=col, lw=0.8, ls=ls, alpha=0.3,
                zorder=2)
        zs.append(np.where(faint | ~view, np.nan, med))
    ax.axhline(MHW_M, color=INK_MUTED, lw=0.6, ls=":", zorder=0)
    ax.text(xlim[0] + 0.01 * (xlim[1] - xlim[0]), MHW_M + 0.05, "MHW", ha="left", va="bottom", fontsize=7,
            color=INK_MUTED)
    lo, hi = np.nanmin(zs), np.nanmax(zs)
    ax.set(xlim=xlim, ylim=(max(lo, -2.0) - 0.3, hi + 0.5),
           xlabel="Distance seaward of the 2008-2010 shoreline (m)", ylabel="Elevation (m NAVD88)")
    open_frame(ax)
    return matplotlib.cm.ScalarMappable(norm=norm, cmap=CMAP), drawn


# One figure: plan view and profiles, saved with its caption
def draw(gis, xs, profiles, anchors, surveys, factors, offsets, datum, path, heading, why, years, xlim):
    a = np.array([v for v in anchors.values()])
    ends = np.vstack([a[:, :2] + xlim[0] * a[:, 2:], a[:, :2] + xlim[1] * a[:, 2:]])
    pad = 30.0
    bounds = (ends[:, 0].min() - pad, ends[:, 1].min() - pad, ends[:, 0].max() + pad, ends[:, 1].max() + pad)
    z, extent = basemap(surveys[BASEMAP_YEAR], factors[BASEMAP_YEAR], offsets.get(BASEMAP_YEAR, 0.0), bounds)
    fig = plt.figure(figsize=figsize("double", 0.42), layout="constrained")
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 2.1])
    axm, axp = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    draw_map(axm, z, extent, anchors, xlim)
    title(axm, 0, f"GIS {gis}, {BASEMAP_YEAR} lidar")
    sm, drawn = draw_profiles(axp, xs, profiles, years, xlim)
    title(axp, 1, heading)
    cb = fig.colorbar(sm, ax=axp, fraction=0.04, pad=0.01)
    flagged = {y: d for y, d in datum.items() if y in drawn and abs(d) > DATUM_FLAG_M}
    cb.set_ticks(drawn, labels=[f"{y}*" if y in flagged else str(y) for y in drawn])
    cb.ax.tick_params(labelsize=7, length=2)
    cb.set_label("Survey year")
    shifts = ", ".join(f"{y} {offsets[y]:+.2f} m" for y in drawn if y in offsets)
    storms = ", ".join(f"{y} {n}" for y, n in POST_STORM.items() if y in drawn)
    left_out = [y for y in profiles if y not in drawn and np.isfinite(median_profile(profiles[y])).any()]
    faint = xlim[0] < -BEHIND_DUNE_M and any(y < BARE_EARTH_FROM for y in drawn)
    caption(fig, (
        f"Observed cross-shore elevation profiles at GIS {gis}, lidar surveys {', '.join(map(str, drawn))}"
        + (f" (also surveyed, not drawn: {', '.join(map(str, left_out))})" if left_out else "")
        + f". (a) The {BASEMAP_YEAR} lidar (m above MHW) over a hillshade, with the domain's {len(anchors)} CoastSat "
        f"transects (white) over the plotted extent and the CoastSat shoreline they are measured from (dotted: each "
        f"transect's median position, {ANCHOR_WINDOW[0]} to {ANCHOR_WINDOW[1]}, the 2009 model start's shoreline). "
        f"(b) Median elevation over the transects at each distance from that shoreline, drawn where at least "
        f"{MIN_TRANSECT_SHARE:.0%} of them have data, coloured by survey year from red (oldest) to purple (newest); "
        f"shaded: the {BAND_PCT[0]}th-{BAND_PCT[1]}th percentile across transects for {drawn[0]} and {drawn[-1]}"
        + (f"; dashed: flown within weeks of a storm ({storms})" if storms else "")
        + (f"; faint: before {BARE_EARTH_FROM} the surveys are first-return grids that keep buildings and trees, so "
           f"landward of {BEHIND_DUNE_M} m they are drawn faint" if faint else "")
        + f". Surveys sampled bilinearly every 1 m; the 1996-2005 grids are 3-9 m and smooth sharp crests, and 2001 "
        f"is unclassified returns. Vertical shifts to the {REF_YEAR} lidar on NC-12, as in the dune-crest work: "
        f"{shifts or 'none'}. Points more than {SPIKE_M:g} m above both the survey before and after at the same "
        f"distance are lidar spikes and removed. "
        + (f"Starred: on ground behind the dune that the {BARE_EARTH_FROM}-2019 surveys agree has not "
           f"changed (within {STABLE_RANGE_M} m), the survey reads "
           + ", ".join(f"{y} {d:+.2f} m" for y, d in flagged.items())
           + f" from {REF_YEAR}, more than {DATUM_FLAG_M} m; not corrected. " if flagged else "")
        + f"MHW {MHW_M} m NAVD88." + (f" {why}" if why else "")))
    save(fig, path, close=True)


# Run: transects, surveys, units, then each domain's figure and the examples
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--gis", type=int, nargs="+", help="only these domains (default: all 90)")
    ap.add_argument("--examples-only", action="store_true", help="draw only the examples")
    ap.add_argument("--all-surveys", action="store_true",
                    help="draw every survey at full width, into *_all_surveys/, to check")
    args = ap.parse_args()
    apply_style()

    # Transects per domain, and the domains to draw
    lookup = pd.read_csv(transect_lookup())
    by_gis = {int(g): list(t.transect_id) for g, t in lookup.groupby("domain_number")}
    want = set(args.gis or by_gis)
    domains = [] if args.examples_only else sorted(want)
    examples = {k: v for k, v in EXAMPLES.items() if v[0] in want}
    geom = transect_geometry(set(lookup.transect_id))

    # Surveys, and their vertical units from the check domains
    print("opening surveys")
    surveys = open_surveys()
    offsets = survey_offsets()
    raw = {}
    for g in sorted(set(domains) | {v[0] for v in examples.values()} | set(UNIT_CHECK_GIS)):
        pts = transect_points({t: geom[t] for t in by_gis[g] if t in geom})
        raw[g] = (pts, sample_domain(surveys, pts))
        print(f"  sampled GIS {g}", flush=True)
    factors = unit_factors(surveys, {g: raw[g][1] for g in UNIT_CHECK_GIS})

    # Each domain's figure, then the examples
    years = list(survey_files()) if args.all_surveys else list(SHOW_YEARS)
    full = (-LANDWARD_M, SEAWARD_M)
    tag = "_all_surveys" if args.all_surveys else ""
    rows, done = [], {}
    for g in sorted(set(domains) | {v[0] for v in examples.values()}):
        pts, r = raw[g]
        xs, profiles, anchors, cover = domain_profiles(g, r, pts, geom, factors, offsets)
        datum = {c["survey"]: c["interior_offset_m"] for c in cover if np.isfinite(c["interior_offset_m"])}
        done[g] = (xs, profiles, anchors, datum)
        rows += cover
        if g in domains:
            draw(g, xs, profiles, anchors, surveys, factors, offsets, datum,
                 OUT / f"by_domain{tag}" / f"gis{g:02d}_profile", f"Median over {len(anchors)} transects", "",
                 years, full if args.all_surveys else ZOOM_M)
    for stem, (g, heading, why) in examples.items():
        xs, profiles, anchors, datum = done[g]
        wide = args.all_surveys or stem in WIDE_EXAMPLES
        draw(g, xs, profiles, anchors, surveys, factors, offsets, datum, OUT / f"examples{tag}" / stem,
             heading, why, years, full if wide else ZOOM_M)
        print(f"  example {stem}")

    # Coverage table beside the figures, from a full run only
    if args.gis or args.examples_only:
        print(f"wrote {len(domains)} domain figures and {len(examples)} examples; tables left as they were "
              f"(written by a full run)")
        return
    out = support_dir(OUT) / "transect_profile_coverage.csv"
    table = pd.DataFrame(rows)
    table.to_csv(out, index=False)
    wide = table.pivot(index="gis", columns="survey", values="interior_offset_m")
    wide.to_csv(support_dir(OUT) / "interior_offset_vs_2009.csv")
    print(f"surveys more than {DATUM_FLAG_M} m from {REF_YEAR} behind the dune, domains flagged / checked:")
    for y in wide.columns:
        col = wide[y].dropna()
        if y != REF_YEAR and len(col):
            print(f"  {y}: {(col.abs() > DATUM_FLAG_M).sum()} / {len(col)}, median {col.median():+.2f} m")
    print(f"wrote {len(domains)} domain figures, {len(examples)} examples, and {out}")


if __name__ == "__main__":
    main()
