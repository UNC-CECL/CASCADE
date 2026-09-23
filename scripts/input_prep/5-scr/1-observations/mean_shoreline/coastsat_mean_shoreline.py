"""
coastsat_mean_shoreline.py -- an averaging window of CoastSat, as a line
==============================================================================
One averaging window's MEAN satellite shoreline, placed back on the ground:
each CoastSat transect's chainage averaged over the window, geolocated, and
the ~906 mean points strung into a single polyline that the 2-brie-offset
intersection step reads exactly as it reads a digitised dune line.

Built 2026-09-22 (Hannah, by interview) so BRIE's island offset can be derived
from the satellite SHORELINE as well as from the digitised DUNE line.

WHY THE GEOLOCATION STEP IS THE WHOLE JOB
    Every other CoastSat product in 5-scr is a DIFFERENCE of chainage -- an
    LRR slope, an endpoint change -- and in a difference each transect's
    arbitrary origin cancels. A position is not a difference, so here the
    origin does not cancel, and it is not small. Aggregated to the 90 domains:

        raw mean chainage        alongshore range  124 m, median step  8.8 m
        geolocated mean position alongshore range 6222 m, median step 81.7 m
        the transect ORIGINS alone                6169 m

    The origins follow the shore around the cape, so ~98% of a raw-chainage
    "island shape" is origin bookkeeping. Each observation is therefore put
    back in space as

        point = origin + chainage * unit_vector_along_transect

    in EPSG:26918, the CRS the other dune lines declare, before anything is
    averaged alongshore.

WHY A MEAN AND NOT A DATE
    A dune line is digitised from imagery flown on ONE day, so it is a moment.
    A single satellite pass is not: it carries tide, wave setup and cloud-edge
    noise worth metres. The position a period starts from is therefore a mean
    over a window of passes. For 1995-1997 that is a median of 28 positions per
    transect, scatter 9-18 m, so the standard error on each transect mean is
    2-3 m -- inside the 10 m Barrier3D cell.

    The window is the CALENDAR span, not a span centred on the 1997 dune
    survey ([[cascade-period-is-the-calendar-year]]). The two are ~9 months
    apart and that mismatch is REPORTED in PROVENANCE.md, not corrected.

WHAT IS NOT DONE TO THE DATA
    No outlier rejection. The 9-18 m scatter within a window is the beach
    moving, not error to be cleaned, and per-transect sd, n and date span are
    written to the CSV so a reader can judge it. No smoothing of the line: the
    mean points are 50 m apart and the 100 m transect frame samples them.
    A transect with fewer than --min-obs positions is EXCLUDED and listed, not
    silently dropped.

ON THE TIDE
    The chainages come from coastsat.space, whose transect layer carries the
    per-transect beach_slope (and its confidence interval) used for tidal
    correction -- which is good evidence these series are already tidally
    corrected, but it is not a statement from the download, and nothing in
    this repository records one. See PROVENANCE.md. It matters less than it
    looks: island_offset_hybrid.py zeroes each build on its own minimum, so a
    UNIFORM tidal bias cancels entirely and only the alongshore VARIATION in
    beach slope (0.04-0.06 here) survives, worth a few metres.

OUTPUT   data/hatteras_init/5-scr/1-observations/mean_shoreline/<start>_<end>/
    shoreline_mean_<start>_<end>.geojson   ONE LineString, EPSG:26918, with the
                                           metadata properties a dune line
                                           carries so step 1 of the offset
                                           build reads it unchanged
    transect_means_<start>_<end>.csv       per CoastSat transect: n, mean, sd,
                                           se, first/last date, the geolocated
                                           mean point, domain, included/why not
    mean_shoreline_<start>_<end>.png       the diagnostic: where the line is,
                                           and how well sampled it is
    mean_shoreline_<start>_<end>_island_outline.png
                                           its panel (a) alone, over the
                                           island outline
    PROVENANCE.md
    Read through hat_observed_rates.mean_shoreline_{dir,geojson,csv}().

USAGE
    python coastsat_mean_shoreline.py
    python coastsat_mean_shoreline.py --window 1995 1997 --min-obs 10

THEN (the offset build, which this script does not do)
    duneline_to_raw_offsets.py --duneline <the geojson above>
        --out 1995_1997_shoreline_offset_raw.csv
    island_offset_hybrid.py --year 1996 --source shoreline --version v1
        --raw-file <the raw file above>
==============================================================================
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import pyproj

_REPO = next(_p for _p in Path(__file__).resolve().parents
             if (_p / "pyproject.toml").exists())
sys.path.insert(0, str(_REPO / "scripts"))
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)

from coastsat_lrr import filter_dates, load_timeseries  # noqa: E402
from site_layer import hat_figure_style as fs  # noqa: E402
from site_layer.hat_observed_rates import (  # noqa: E402
    COASTSAT_TIMESERIES, TRANSECT_LAYER, mean_shoreline_csv,
    mean_shoreline_dir, mean_shoreline_geojson, transect_lookup,
)

# The CRS the 1984, 2009 and 2023 dune lines declare. The 100 m transects are
# EPSG:3725, NAD83(NSRS2007) / UTM 18N, which is the same grid to within the
# null transform, and duneline_to_raw_offsets.py reprojects anyway.
TARGET_CRS = "EPSG:26918"
LAYER_CRS = "EPSG:4326"

DEFAULT_WINDOW = (1995, 1997)
# The house minimum, from coastsat_lrr.MIN_OBS: fewer than ten positions is
# not a mean of a seasonal cycle.
DEFAULT_MIN_OBS = 10

# Copied into the raw offsets file by duneline_to_raw_offsets.LINE_META, so a
# shoreline-derived raw file carries as full a provenance as a dune-derived one.
SOURCE_TYPE = "Landsat 5/7/8 via CoastSat (coastsat.space)"


# =============================================================================
# the transect layer
# =============================================================================

def transect_geometry(wanted):
    """Origin and seaward unit vector of each wanted transect, in TARGET_CRS.

    The layer is the global CoastSat file (233k transects, 80 MB), so it is
    read once and filtered to the ids the domain lookup names. Chainage is
    measured from the FIRST vertex along the line, so that vertex is the
    origin and the line's own direction is the unit vector.
    """
    to_utm = pyproj.Transformer.from_crs(LAYER_CRS, TARGET_CRS, always_xy=True)
    out = {}
    with open(TRANSECT_LAYER, "r", encoding="utf-8") as fh:
        layer = json.load(fh)
    for feat in layer["features"]:
        # The layer spells an id "usa_NC_0032-0001"; the timeseries file and
        # the domain lookup both spell it with an underscore.
        tid = str(feat["properties"]["id"]).replace("-", "_")
        if tid not in wanted:
            continue
        coords = feat["geometry"]["coordinates"]
        x0, y0 = to_utm.transform(*coords[0])
        x1, y1 = to_utm.transform(*coords[-1])
        length = float(np.hypot(x1 - x0, y1 - y0))
        if length <= 0:
            continue
        out[tid] = (x0, y0, (x1 - x0) / length, (y1 - y0) / length, length)
    return out


def timeseries_file(transect_id):
    """The per-transect CSV, which sits in its site's folder."""
    site = transect_id.rsplit("_", 1)[0]
    return COASTSAT_TIMESERIES / "{0}_timeseries".format(site) / "{0}.csv".format(transect_id)


# =============================================================================
# the window mean
# =============================================================================

def window_means(lookup, geometry, window, min_obs):
    """One row per CoastSat transect: its mean position over the window.

    The window is inclusive of both calendar years. A transect is EXCLUDED,
    with the reason recorded, when its geometry or timeseries is missing or
    when it holds fewer than `min_obs` positions -- never dropped in silence.
    """
    start, end = int(window[0]), int(window[1])
    lo, hi = "{0}-01-01".format(start), "{0}-12-31".format(end)
    rows = []
    for tid, domain in zip(lookup["transect_id"], lookup["domain_number"]):
        row = {"transect_id": tid, "site": tid.rsplit("_", 1)[0],
               "transect_number": int(tid.rsplit("_", 1)[1]),
               "domain_number": int(domain)}
        geom = geometry.get(tid)
        path = timeseries_file(tid)
        if geom is None:
            rows.append(dict(row, n_obs=0, included=False,
                             excluded_because="no geometry in the transect layer"))
            continue
        if not path.is_file():
            rows.append(dict(row, n_obs=0, included=False,
                             excluded_because="no timeseries file"))
            continue
        obs = filter_dates(load_timeseries(str(path)), lo, hi)
        n = len(obs)
        row["n_obs"] = n
        if n:
            ch = obs["chainage_m"].to_numpy(dtype=float)
            x0, y0, ux, uy, _ = geom
            sd = float(ch.std(ddof=1)) if n > 1 else np.nan
            row.update(
                mean_chainage_m=float(ch.mean()),
                sd_chainage_m=sd,
                se_chainage_m=sd / np.sqrt(n) if n > 1 else np.nan,
                min_chainage_m=float(ch.min()), max_chainage_m=float(ch.max()),
                first_date=obs["date"].iloc[0].date().isoformat(),
                last_date=obs["date"].iloc[-1].date().isoformat(),
                n_years=int(obs["date"].dt.year.nunique()),
                x=x0 + ch.mean() * ux, y=y0 + ch.mean() * uy,
            )
        keep = bool(n >= min_obs)
        row["included"] = keep
        row["excluded_because"] = ("" if keep else
                                   "n_obs {0} < min_obs {1}".format(n, min_obs))
        rows.append(row)
    df = pd.DataFrame(rows).sort_values(["site", "transect_number"])
    return df.reset_index(drop=True)


def line_vertices(df):
    """The included mean points in alongshore order.

    Ordering is (site, transect number), which is alongshore here: the six
    Hatteras sites chain south to north with monotonically increasing domain
    spans, and the resulting vertex spacing is ~50 m with no gap over 300 m
    (checked 2026-09-22). Where two sites overlap at a domain boundary the
    line can double back a little; that is why the intersection step records
    n_crossings and takes one crossing per transect.
    """
    keep = df[df["included"] & df["x"].notna()]
    return keep.sort_values(["site", "transect_number"]).reset_index(drop=True)


# =============================================================================
# outputs
# =============================================================================

def write_geojson(vertices, path, window, built_on, n_total):
    """ONE LineString, with the properties a digitised dune line carries.

    duneline_to_raw_offsets.py exits on a file holding more than one feature
    and copies LINE_META across, so this has to be a single feature and it is
    worth filling the metadata in: the raw offsets file is where a reader
    meets this line next.
    """
    start, end = int(window[0]), int(window[1])
    props = {
        "feature_type": "Shoreline (CoastSat window mean)",
        # A window, not a moment -- rule 2. The raw offsets file carries this
        # string in its `year` column, where "1997" would be a lie.
        "year": "{0}-{1}".format(start, end),
        "imagery_date": "{0}-01-01/{1}-12-31".format(start, end),
        "source_type": SOURCE_TYPE,
        "method": ("Mean of the satellite positions in calendar {0}-{1} per "
                   "CoastSat transect ({2} of {3} transects), geolocated as "
                   "origin + chainage * unit vector and strung into one line "
                   "in alongshore order. No outlier rejection, no smoothing."
                   .format(start, end, len(vertices), n_total)),
        "editor": "coastsat_mean_shoreline.py",
        "edit_date": built_on,
        "notes": ("Window mean, not a survey. Median {0} positions per "
                  "transect; see transect_means_{1}_{2}.csv and PROVENANCE.md."
                  .format(int(vertices["n_obs"].median()), start, end)),
    }
    feature = {
        "type": "Feature", "id": 1, "properties": props,
        "geometry": {"type": "LineString",
                     "coordinates": [[round(float(x), 4), round(float(y), 4)]
                                     for x, y in zip(vertices["x"], vertices["y"])]},
    }
    doc = {"type": "FeatureCollection",
           "crs": {"type": "name", "properties": {"name": TARGET_CRS}},
           "features": [feature]}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)


def figure(df, vertices, folder, window):
    """Two panels: where the line is, and how well sampled it is."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fs.apply_style()
    start, end = int(window[0]), int(window[1])
    # Two panels STACKED, and the map drawn with northing across the page.
    # The island is ~45 km north-south by ~6 km east-west, so at equal aspect
    # -- which a map has to keep, or the shape it is showing is not the shape
    # on the ground -- a portrait panel is an unreadable sliver. Turned on its
    # side it is a wide ribbon, and alongshore runs left to right as it does
    # in every other figure here.
    fig, axes = plt.subplots(
        2, 1, figsize=fs.figsize("double", aspect=0.62),
        gridspec_kw={"height_ratios": [1, 2.1], "hspace": 0.55})

    ax = axes[0]
    ax.plot(vertices["y"] / 1000.0, vertices["x"] / 1000.0, "-",
            color=fs.C["LATE"], lw=1.1, zorder=3)
    ax.set_xlabel("Northing (km, EPSG:26918)  —  south → north")
    ax.set_ylabel("Easting (km)")
    ax.set_aspect("equal")
    fs._title(ax, 0, "Mean shoreline, {0}–{1}".format(start, end))

    ax = axes[1]
    by_domain = df[df["included"]].groupby("domain_number")
    ax.bar(by_domain.size().index, by_domain.size().values, width=0.9,
           color=fs.C["BASE_FILL"], edgecolor="none", zorder=2)
    ax.set_ylabel("CoastSat transects used", color=fs.C["BASE"])
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    twin = ax.twinx()
    med = by_domain["n_obs"].median()
    twin.plot(med.index, med.values, "-", color=fs.C["LATE"], lw=1.2, zorder=3)
    twin.set_ylabel("median positions per transect", color=fs.C["LATE"])
    twin.set_ylim(bottom=0)
    fs._title(ax, 1, "Sampling per domain")

    fs.caption(fig, (
        "The CoastSat mean shoreline for calendar {0}–{1}. (a) the {2} "
        "per-transect window means, geolocated and strung into the line the "
        "island-offset build intersects with the 100 m transect frame, drawn "
        "with alongshore across the page and at equal aspect. (b) how "
        "many CoastSat transects each Barrier3D domain contributes (bars) and "
        "the median number of satellite positions behind each of those "
        "transect means (line). Positions are means over the window, not a "
        "survey on a date.".format(start, end, len(vertices))))
    out = fs.save(fig, folder / "mean_shoreline_{0}_{1}.png".format(start, end),
                  close=True)
    print("  figure -> {0}".format(out[0].name))


# The ribbon: panel (a) of the diagnostic on its own, over a map (Hannah,
# 2026-09-23). Northing across the page, easting up, equal aspect. Shared with
# coastsat_mean_shoreline_on_imagery.py, which draws the same ribbon on photos.
RIBBON_PAD_M = (3500.0, 1200.0, 600.0)    # landward, seaward, alongshore


def ribbon_extent(vertices):
    """(n0, n1, e0, e1): the line's northing span, and easting wide enough to
    hold the island landward of it (Buxton Woods is ~3 km across)."""
    land, sea, along = RIBBON_PAD_M
    return (vertices["y"].min() - along, vertices["y"].max() + along,
            vertices["x"].min() - land, vertices["x"].max() + sea)


def ribbon_axes(ax, ext):
    """Km ticks in the same words as the diagnostic's panel (a), but easting
    increasing DOWN (Hannah, 2026-09-23: "flip these vertically"): the ocean is
    at the bottom and the ribbon is a north-up map turned 90 degrees clockwise,
    where panel (a)'s easting-up axes draw its mirror image."""
    n0, n1, e0, e1 = ext
    ax.set_xlim(n0, n1)
    ax.set_ylim(e1, e0)
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(lambda v, _: "{0:g}".format(v / 1000.0))
    ax.yaxis.set_major_formatter(lambda v, _: "{0:g}".format(v / 1000.0))
    ax.set_xlabel("Northing (km, EPSG:26918)  —  south → north")
    ax.set_ylabel("Easting (km)")
    fs.spines_for_image(ax)


def draw_island_outline(ax, ext, edge_on_top=None):
    """Water tint and the island outline, in the ribbon's swapped axes.

    `edge_on_top` (a colour) also draws the outline's edge above everything
    at zorder 3.5, between the photographs and the mean line; the imagery
    ribbon passes "white" (Hannah, 2026-09-23) so the outline reads against
    the photographs."""
    import geopandas as gpd
    from shapely.affinity import affine_transform
    from shapely.geometry import box
    from site_layer.hat_map_layers import ISLAND_OUTLINE

    n0, n1, e0, e1 = ext
    # swap the axes: (easting, northing) -> (northing, easting)
    isl = gpd.read_file(ISLAND_OUTLINE).to_crs(TARGET_CRS)
    isl = isl.clip(box(e0, n0, e1, n1)).geometry.apply(
        lambda g: affine_transform(g, [0, 1, 1, 0, 0, 0]))
    ax.set_facecolor("#eef4f9")                           # water
    gpd.GeoSeries(isl).plot(ax=ax, color="0.88", edgecolor=fs.INK_MUTED, lw=0.4, zorder=1)
    if edge_on_top:
        gpd.GeoSeries(isl).boundary.plot(ax=ax, color=edge_on_top, lw=0.6, zorder=3.5)


def outline_figure(vertices, folder, window):
    """The line over the island outline, alone."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fs.apply_style()
    start, end = int(window[0]), int(window[1])
    ext = ribbon_extent(vertices)
    n0, n1, e0, e1 = ext
    w = fs.FIG_W_DOUBLE
    fig, ax = plt.subplots(figsize=(w, (e1 - e0) / (n1 - n0) * (w - 0.8) + 1.0),
                           layout="constrained")
    draw_island_outline(ax, ext)
    ax.plot(vertices["y"], vertices["x"], "-", color=fs.C["LATE"], lw=1.1, zorder=3)
    ribbon_axes(ax, ext)
    ax.set_title("Mean shoreline, {0}–{1}".format(start, end))    # one panel: no letter
    fs.caption(fig, (
        "The CoastSat mean shoreline for calendar {0}–{1} (blue), the {2} per-transect "
        "window means geolocated and strung into one line, over the Hatteras Island "
        "outline (grey; map_elements/hatteras_outline), drawn with alongshore across the "
        "page and at equal aspect, with easting increasing downward, so the ocean is at the bottom. Panel (a) of "
        "mean_shoreline_{0}_{1}.png on its own.".format(start, end, len(vertices))))
    out = fs.save(fig, folder / "mean_shoreline_{0}_{1}_island_outline.png".format(start, end),
                  close=True)
    print("  figure -> {0}".format(out[0].name))


def _excluded_table(excluded):
    """The excluded transects as a markdown table.

    Spelled out rather than DataFrame.to_markdown(), which wants `tabulate`;
    a provenance note is not worth a dependency (2026-09-22).
    """
    if not len(excluded):
        return "No transect was excluded."
    cols = ["transect_id", "domain_number", "n_obs", "excluded_because"]
    lines = ["### Excluded transects", "",
             "| " + " | ".join(cols) + " |",
             "|" + "---|" * len(cols)]
    for _, r in excluded.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


def write_provenance(df, vertices, folder, window, min_obs, built_on):
    start, end = int(window[0]), int(window[1])
    excluded = df[~df["included"]]
    per_domain = df[df["included"]].groupby("domain_number").size()
    step = np.hypot(np.diff(vertices["x"]), np.diff(vertices["y"]))
    excluded_block = _excluded_table(excluded)

    text = """# mean_shoreline/{start}_{end} -- the CoastSat window mean, as a line

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py`
on {built_on}.

## What this is

The mean satellite shoreline position over **calendar {start}-{end}**, per
CoastSat transect, geolocated and strung into one polyline. It is the
shoreline counterpart of a digitised dune line, and `2-brie-offset` consumes
it the same way.

It is **a window mean, not a survey on a date.** A single satellite pass
carries metres of tide, wave setup and cloud-edge noise; the mean over a
window of passes is the quantity a period can start from.

## The numbers

| | |
|---|---|
| CoastSat transects in the domain lookup | {n_total} |
| used (n_obs >= {min_obs}) | {n_used} |
| excluded | {n_excluded} |
| positions per transect | median {n_med}, min {n_min}, max {n_max} |
| within-window scatter (sd) | median {sd_med:.1f} m |
| standard error of a transect mean | median {se_med:.1f} m |
| domains covered | {d_lo}-{d_hi} ({n_dom} of 90) |
| transects per domain | min {t_min}, max {t_max} |
| vertex spacing along the line | median {s_med:.0f} m, p95 {s_p95:.0f} m, max {s_max:.0f} m |

The standard error is the number that matters for an island offset: at
~{se_med:.0f} m it is inside the 10 m Barrier3D cell, so the alongshore shape of
this line is not sampling noise.

## Four things a reader should know

**1. The geolocation is the method, not a formality.** Aggregated to the 90
domains, raw mean chainage spans 124 m alongshore; the geolocated position
spans 6222 m, and the transect origins alone span 6169 m. The CoastSat
transect origins follow the shore around the cape, so a "shape" read straight
off chainage would be ~98% origin bookkeeping. Every other CoastSat product in
`5-scr` is a *difference* of chainage, where the origin cancels; this one is
not.

**2. The window is the calendar span, and it does not match the dune line.**
The 1996 period's dune line is digitised from imagery flown 1997-10-12,
**15.4 months** after the centre of this window (1996-07-01). Per
[[cascade-period-is-the-calendar-year]] the mismatch is reported here rather
than corrected by re-centring the window on the survey date. It matters when
this line is differenced against the dune line: part of any gap is those
fifteen months of shoreline change, not beach width.

**3. The years inside the window are not evenly sampled.** Landsat 7 does not
launch until 1999, so {end} is thinner than the earlier years -- across a
sample of 80 transects, roughly 735 / 976 / 359 positions for 1995 / 1996 /
1997. The plain mean therefore leans slightly toward the early window. A
year-balanced mean was offered and not taken (Hannah, 2026-09-22); it would be
a small change to `window_means`.

**4. Tidal correction is probable but unrecorded.** The transect layer from
coastsat.space carries `beach_slope`, `cil` and `ciu` -- the per-transect
slope used for tidal correction and its confidence interval -- which is good
evidence these chainages are already tidally corrected. **Nothing in this
repository records a statement from the download, so it is not asserted
here.** The exposure is limited: `island_offset_hybrid.py` zeroes each build
on its own minimum, so a uniform tidal bias cancels completely and only the
alongshore variation in beach slope (0.04-0.06 here) survives, worth a few
metres.

## What was not done to the data

No outlier rejection -- the {sd_lo:.0f}-{sd_hi:.0f} m within-window scatter is the beach
moving, and `sd`, `se`, `n_obs` and the date span are in the CSV so a reader
can judge each transect. No smoothing of the line. No gap filling. A transect
below the minimum is excluded and named, never silently dropped.

{excluded_block}

## Files

| file | what it is |
|---|---|
| `shoreline_mean_{start}_{end}.geojson` | one LineString, EPSG:26918, with the metadata properties a dune line carries |
| `transect_means_{start}_{end}.csv` | per transect: n, mean, sd, se, date span, the geolocated point, domain, included/why not |
| `mean_shoreline_{start}_{end}.png` | the diagnostic figure |
| `mean_shoreline_{start}_{end}_island_outline.png` | panel (a) of the diagnostic alone, the line over the island outline |
| `on_imagery/` | the line and its ±1 sd band on the USGS photographs flown inside the window, at six sites and island-wide (three segments, and one ribbon panel); written by `coastsat_mean_shoreline_on_imagery.py`, which needs the D: drive (see its supporting/ folders) |

Resolved through `hat_observed_rates.mean_shoreline_dir/_geojson/_csv`.
Never type these paths.
""".format(
        start=start, end=end, built_on=built_on,
        n_total=len(df), n_used=len(vertices), n_excluded=len(excluded),
        min_obs=min_obs,
        n_med=int(vertices["n_obs"].median()), n_min=int(vertices["n_obs"].min()),
        n_max=int(vertices["n_obs"].max()),
        sd_med=vertices["sd_chainage_m"].median(),
        se_med=vertices["se_chainage_m"].median(),
        sd_lo=vertices["sd_chainage_m"].quantile(0.05),
        sd_hi=vertices["sd_chainage_m"].quantile(0.95),
        d_lo=int(per_domain.index.min()), d_hi=int(per_domain.index.max()),
        n_dom=len(per_domain), t_min=int(per_domain.min()),
        t_max=int(per_domain.max()),
        s_med=np.median(step), s_p95=np.percentile(step, 95), s_max=step.max(),
        excluded_block=excluded_block,
    )
    (folder / "PROVENANCE.md").write_text(text, encoding="utf-8")


# =============================================================================

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--window", nargs=2, type=int, default=list(DEFAULT_WINDOW),
                    metavar=("START", "END"),
                    help="inclusive calendar years to average (default 1995 1997)")
    ap.add_argument("--min-obs", type=int, default=DEFAULT_MIN_OBS,
                    help="positions a transect needs to contribute (default 10)")
    a = ap.parse_args(argv)
    window = (a.window[0], a.window[1])
    if window[1] < window[0]:
        ap.error("the window ends before it starts")
    built_on = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    lookup = pd.read_csv(transect_lookup()).dropna(subset=["domain_number"])
    print("Window     : calendar {0}-{1}".format(*window))
    print("Transects  : {0} in the domain lookup".format(len(lookup)))
    print("Layer      : {0}".format(TRANSECT_LAYER.name))
    geometry = transect_geometry(set(lookup["transect_id"]))
    print("  geolocated {0} of {1}".format(len(geometry), len(lookup)))

    df = window_means(lookup, geometry, window, a.min_obs)
    vertices = line_vertices(df)
    print("  {0} transects used, {1} excluded".format(
        len(vertices), len(df) - len(vertices)))
    if len(vertices) < 2:
        sys.exit("fewer than two usable transects; nothing to build a line from")
    print("  positions per transect: median {0}, min {1}".format(
        int(vertices["n_obs"].median()), int(vertices["n_obs"].min())))
    print("  standard error of a transect mean: median {0:.1f} m".format(
        vertices["se_chainage_m"].median()))

    folder = mean_shoreline_dir(*window)
    folder.mkdir(parents=True, exist_ok=True)
    csv_path = mean_shoreline_csv(*window)
    df.to_csv(csv_path, index=False)
    print("\nWrote {0}  ({1} rows)".format(csv_path, len(df)))
    geo_path = mean_shoreline_geojson(*window)
    write_geojson(vertices, geo_path, window, built_on, len(df))
    print("Wrote {0}  ({1} vertices)".format(geo_path, len(vertices)))
    figure(df, vertices, folder, window)
    outline_figure(vertices, folder, window)
    write_provenance(df, vertices, folder, window, a.min_obs, built_on)
    print("Wrote {0}".format(folder / "PROVENANCE.md"))


if __name__ == "__main__":
    main()
