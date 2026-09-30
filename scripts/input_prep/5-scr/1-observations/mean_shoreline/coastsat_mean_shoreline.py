"""
One averaging window of CoastSat as a line: each transect's mean position, strung into a mean shoreline.

    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py
    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py --window 1995 1997 --min-obs 10
    python scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py --centred-on alace_1996

Writes the mean shoreline (points, line, CSV), its provenance and a
diagnostic figure under the window's folder. Details: scripts/input_prep/5-scr/1-observations/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-29
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
    mean_shoreline_dir, mean_shoreline_geojson, mean_shoreline_label,
    transect_lookup,
)

# --- CONFIG ------------------------------------------------------------------
# The CRS the 1984, 2009 and 2023 dune lines declare
TARGET_CRS = "EPSG:26918"
LAYER_CRS = "EPSG:4326"

DEFAULT_WINDOW = (1995, 1997)
# The house minimum, from coastsat_lrr.MIN_OBS
DEFAULT_MIN_OBS = 10

# Copied into the raw offsets file by duneline_to_raw_offsets.LINE_META
SOURCE_TYPE = "Landsat 5/7/8 via CoastSat (coastsat.space)"

# The lidar surveys a period's start topography is built on, for --centred-on
SURVEY_ANCHORS = {
    "alace_1996": {
        "period_start": 1996,
        "survey": "1996 fall East Coast NOAA/NASA ALACE lidar",
        "flown": ("1996-10-09", "1996-10-16"),
        "centre": "1996-10-12",
        "source": "NOAA InPort item 48147, https://www.fisheries.noaa.gov/inport/item/48147",
        "dem": ("the beach and foredune of `0-elevation/2009-2014-1996`, which the "
                "1996 start's topography (`1984-start`) is built on"),
    },
    "usace_2009": {
        "period_start": 2010,
        "survey": "2009 USACE NCMP topobathy lidar (CHARTS)",
        "flown": ("2009-08-10", "2009-08-24"),
        "centre": "2009-08-17",
        "source": "NOAA InPort item 54934, https://www.fisheries.noaa.gov/inport/item/54934",
        "dem": ("`0-elevation/2009-2014`, the 2010 start's topography (`2004-start`): "
                "2009 USACE wherever it measured, 2014 Post-Sandy only in its nodata"),
    },
}
DEFAULT_HALF_WIDTH_YEARS = 1
# -----------------------------------------------------------------------------


# An averaging window
class Window:

    def __init__(self, lo, hi, calendar, anchor=None, half_width=None):
        self.lo, self.hi = pd.Timestamp(lo), pd.Timestamp(hi)
        if self.hi < self.lo:
            raise ValueError("the window ends before it starts")
        self.calendar = calendar
        self.anchor = anchor            # a SURVEY_ANCHORS key, or None
        self.half_width = half_width    # years either side of the anchor

    @classmethod
    def from_years(cls, start, end):
        return cls("{0}-01-01".format(int(start)), "{0}-12-31".format(int(end)), True)

    @classmethod
    def from_dates(cls, lo, hi):
        return cls(lo, hi, False)

    @classmethod
    def centred_on(cls, anchor, half_width=DEFAULT_HALF_WIDTH_YEARS):
        c = pd.Timestamp(SURVEY_ANCHORS[anchor]["centre"])
        off = pd.DateOffset(years=half_width)
        return cls(c - off, c + off, False, anchor=anchor, half_width=half_width)

    @property
    def lo_iso(self):
        return self.lo.date().isoformat()

    @property
    def hi_iso(self):
        return self.hi.date().isoformat()

    @property
    def key(self):
        if self.calendar:
            return (self.lo.year, self.hi.year)
        return (self.lo_iso, self.hi_iso)

    @property
    def label(self):
        return mean_shoreline_label(*self.key)

    # For titles: `1995–1997`, or `1995-10-12 – 1997-10-12`
    @property
    def span(self):
        if self.calendar:
            return "{0}–{1}".format(*self.key)
        return "{0} – {1}".format(*self.key)

    # For captions and prose: what the window IS, not only its ends
    @property
    def described(self):
        if self.calendar:
            return "calendar {0}".format(self.span)
        text = "{0} to {1}".format(self.lo_iso, self.hi_iso)
        if self.anchor:
            a = SURVEY_ANCHORS[self.anchor]
            text += ", ±{0} yr of the {1} (flown {2} to {3})".format(
                self.half_width, a["survey"], *a["flown"])
        return text

    @property
    def centre(self):
        if self.anchor:
            return pd.Timestamp(SURVEY_ANCHORS[self.anchor]["centre"])
        # midway between the two end days: 1996-07-01 for calendar 1995-1997
        return (self.lo + (self.hi - self.lo) / 2).normalize()

    # The model period this window starts: the anchor's, else the year of its centre
    @property
    def period_start(self):
        if self.anchor:
            return SURVEY_ANCHORS[self.anchor]["period_start"]
        return int(self.centre.year)

    # The observations inside the window
    def clip(self, obs):
        return filter_dates(obs, self.lo_iso, self.hi_iso + " 23:59:59")


# The transect layer

# Origin and seaward unit vector of each wanted transect, in TARGET_CRS
def transect_geometry(wanted):
    to_utm = pyproj.Transformer.from_crs(LAYER_CRS, TARGET_CRS, always_xy=True)
    out = {}
    with open(TRANSECT_LAYER, "r", encoding="utf-8") as fh:
        layer = json.load(fh)
    for feat in layer["features"]:
        # The layer spells an id "usa_NC_0032-0001"
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


# The per-transect CSV, which sits in its site's folder
def timeseries_file(transect_id):
    site = transect_id.rsplit("_", 1)[0]
    return COASTSAT_TIMESERIES / "{0}_timeseries".format(site) / "{0}.csv".format(transect_id)


# The window mean

# One row per CoastSat transect
def window_means(lookup, geometry, window, min_obs):
    rows = []
    by_year = {}
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
        obs = window.clip(load_timeseries(str(path)))
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
        if keep:
            for yr, k in obs["date"].dt.year.value_counts().items():
                by_year[int(yr)] = by_year.get(int(yr), 0) + int(k)
    df = pd.DataFrame(rows).sort_values(["site", "transect_number"])
    return df.reset_index(drop=True), dict(sorted(by_year.items()))


# The included mean points in alongshore order
def line_vertices(df):
    keep = df[df["included"] & df["x"].notna()]
    return keep.sort_values(["site", "transect_number"]).reset_index(drop=True)


# Outputs

# ONE LineString, with the properties a digitised dune line carries
def write_geojson(vertices, path, window, built_on, n_total):
    # A window, not a moment -- rule 2
    if window.calendar:
        year = "{0}-{1}".format(*window.key)
        within = "calendar {0}-{1}".format(*window.key)
    else:
        year = "{0}/{1}".format(window.lo_iso, window.hi_iso)
        within = window.described
    props = {
        "feature_type": "Shoreline (CoastSat window mean)",
        "year": year,
        "imagery_date": "{0}/{1}".format(window.lo_iso, window.hi_iso),
        "source_type": SOURCE_TYPE,
        "method": ("Mean of the satellite positions in {0} per "
                   "CoastSat transect ({1} of {2} transects), geolocated as "
                   "origin + chainage * unit vector and strung into one line "
                   "in alongshore order. No outlier rejection, no smoothing."
                   .format(within, len(vertices), n_total)),
        "editor": "coastsat_mean_shoreline.py",
        "edit_date": built_on,
        "notes": ("Window mean, not a survey. Median {0} positions per "
                  "transect; see transect_means_{1}.csv and PROVENANCE.md."
                  .format(int(vertices["n_obs"].median()), window.label)),
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


# Two panels
def figure(df, vertices, folder, window):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fs.apply_style()
    # Two panels STACKED, and the map drawn with northing across the page
    fig, axes = plt.subplots(
        2, 1, figsize=fs.figsize("double", aspect=0.46),
        gridspec_kw={"height_ratios": [1, 1.1], "hspace": 0.55})

    ax = axes[0]
    ax.plot(vertices["y"] / 1000.0, vertices["x"] / 1000.0, "-",
            color=fs.C["LATE"], lw=1.1, zorder=3)
    ax.set_xlabel("Northing (km, EPSG:26918)  —  south → north")
    ax.set_ylabel("Easting (km)")
    ax.set_aspect("equal")
    # Easting increases DOWN, as in ribbon_axes
    ax.invert_yaxis()
    fs._title(ax, 0, "Mean shoreline, {0}".format(window.span))

    ax = axes[1]
    by_domain = df[df["included"]].groupby("domain_number")
    ax.bar(by_domain.size().index, by_domain.size().values, width=0.9,
           color=fs.C["BASE_FILL"], edgecolor="none", zorder=2)
    ax.set_ylabel("transects used", color=fs.C["BASE"])
    ax.set_xlabel(fs.DOMAIN_AXIS_LABEL)
    twin = ax.twinx()
    med = by_domain["n_obs"].median()
    twin.plot(med.index, med.values, "-", color=fs.C["LATE"], lw=1.2, zorder=3)
    twin.set_ylabel("median positions", color=fs.C["LATE"])
    twin.set_ylim(bottom=0)
    fs._title(ax, 1, "Sampling per domain")

    fs.caption(fig, (
        "The CoastSat mean shoreline for {0}. (a) the {1} "
        "per-transect window means, geolocated and strung into the line the "
        "island-offset build intersects with the 100 m transect frame, drawn "
        "with alongshore across the page and at equal aspect, with easting "
        "increasing downward, so the ocean is at the bottom. (b) how "
        "many CoastSat transects each Barrier3D domain contributes (bars) and "
        "the median number of satellite positions behind each of those "
        "transect means (line). Positions are means over the window, not a "
        "survey on a date.".format(window.described, len(vertices))))
    out = fs.save(fig, folder / "mean_shoreline_{0}.png".format(window.label),
                  close=True)
    print("  figure -> {0}".format(out[0].name))


# The ribbon: panel (a) on its own over a map, shared with the imagery script
RIBBON_PAD_M = (3500.0, 1200.0, 600.0)    # landward, seaward, alongshore


# (n0, n1, e0, e1)
def ribbon_extent(vertices):
    land, sea, along = RIBBON_PAD_M
    return (vertices["y"].min() - along, vertices["y"].max() + along,
            vertices["x"].min() - land, vertices["x"].max() + sea)


# Km ticks in the same words as the diagnostic's panel (a), but easting increasing DOWN (Hannah, ...
def ribbon_axes(ax, ext):
    n0, n1, e0, e1 = ext
    ax.set_xlim(n0, n1)
    ax.set_ylim(e1, e0)
    ax.set_aspect("equal")
    ax.xaxis.set_major_formatter(lambda v, _: "{0:g}".format(v / 1000.0))
    ax.yaxis.set_major_formatter(lambda v, _: "{0:g}".format(v / 1000.0))
    ax.set_xlabel("Northing (km, EPSG:26918)  —  south → north")
    ax.set_ylabel("Easting (km)")
    fs.spines_for_image(ax)


# Water tint and the island outline, in the ribbon's swapped axes
def draw_island_outline(ax, ext, edge_on_top=None):
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


# The line over the island outline, alone
def outline_figure(vertices, folder, window):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fs.apply_style()
    ext = ribbon_extent(vertices)
    n0, n1, e0, e1 = ext
    w = fs.FIG_W_DOUBLE
    fig, ax = plt.subplots(figsize=(w, (e1 - e0) / (n1 - n0) * (w - 0.8) + 1.0),
                           layout="constrained")
    draw_island_outline(ax, ext)
    ax.plot(vertices["y"], vertices["x"], "-", color=fs.C["LATE"], lw=1.1, zorder=3)
    ribbon_axes(ax, ext)
    ax.set_title("Mean shoreline, {0}".format(window.span))    # one panel: no letter
    fs.caption(fig, (
        "The CoastSat mean shoreline for {0} (blue), the {1} per-transect "
        "window means geolocated and strung into one line, over the Hatteras Island "
        "outline (grey; map_elements/hatteras_outline), drawn with alongshore across the "
        "page and at equal aspect, with easting increasing downward, so the ocean is at the bottom. Panel (a) of "
        "mean_shoreline_{2}.png on its own.".format(window.described, len(vertices), window.label)))
    out = fs.save(fig, folder / "mean_shoreline_{0}_island_outline.png".format(window.label),
                  close=True)
    print("  figure -> {0}".format(out[0].name))


# The excluded transects as a markdown table
def _excluded_table(excluded):
    if not len(excluded):
        return "No transect was excluded."
    cols = ["transect_id", "domain_number", "n_obs", "excluded_because"]
    lines = ["### Excluded transects", "",
             "| " + " | ".join(cols) + " |",
             "|" + "---|" * len(cols)]
    for _, r in excluded.iterrows():
        lines.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(lines)


# (vintage, ISO date) of the dune line a period start reads, or None
def dune_line_date(period_start):
    from site_layer.hat_topo_version import dune_line_for_year
    from coastsat_vs_duneline import KNOWN_SURVEY_DATES
    vintage = dune_line_for_year(period_start, strict=False)
    date = KNOWN_SURVEY_DATES.get(vintage) if vintage else None
    return (vintage, date) if date else None


# A time span in months
def _months(delta):
    return delta.days / 30.4375


# Point 2 of the provenance
def _window_paragraph(window):
    dune = dune_line_date(window.period_start)
    dune_text = None
    if dune:
        vintage, date = dune
        gap = _months(pd.Timestamp(date) - window.centre)
        when = ("{0:.1f} months after".format(gap) if gap >= 0
                else "{0:.1f} months before".format(-gap))
        dune_text = ("The {0} period's dune line is digitised from imagery "
                     "flown {1} (the {2} line), **{3}** the centre of this "
                     "window ({4})".format(window.period_start, date, vintage,
                                           when, window.centre.date().isoformat()))

    if window.calendar:
        head = ("**2. The window is the calendar span, and it does not match "
                "the dune line.**")
        body = [
            (dune_text + ".") if dune_text else
            "No dune-line survey date is recorded for the {0} period.".format(
                window.period_start),
            "Per [[cascade-period-is-the-calendar-year]] the mismatch is "
            "reported here rather than corrected by re-centring the window on "
            "the survey date. It matters when this line is differenced against "
            "the dune line: part of any gap is shoreline change over that "
            "interval, not beach width.",
        ]
        return head + "\n" + " ".join(body)

    if window.anchor:
        a = SURVEY_ANCHORS[window.anchor]
        head = ("**2. The window is centred on the start DEM's lidar survey, "
                "not on the calendar.**")
        body = [
            "The {0} was flown {1} to {2} ({3}); the window is ±{4} yr of the "
            "middle of those flights, {5}. That survey is {6}.".format(
                a["survey"], a["flown"][0], a["flown"][1], a["source"],
                window.half_width, a["centre"], a["dem"]),
            "The line becomes the {0} shoreline island offset, a snapshot the "
            "model starts from beside that topography, so it is dated like "
            "the topography (Hannah, 2026-09-29) -- a deliberate exception to "
            "[[cascade-period-is-the-calendar-year]], which still governs the "
            "rates, total change and the scoring target.".format(
                window.period_start),
        ]
    else:
        head = "**2. The window is given by dates, not calendar years.**"
        body = ["{0} to {1}, centred on {2}.".format(
            window.lo_iso, window.hi_iso, window.centre.date().isoformat())]
    if dune_text:
        tail = dune_text
        if dune[1] == window.hi_iso:
            tail += "; the window ends on the dune-line date"
        elif dune[1] == window.lo_iso:
            tail += "; the window starts on the dune-line date"
        body.append(tail + ".")
    return head + "\n" + " ".join(body)


# Point 3
def _sampling_paragraph(window, by_year):
    counts = " / ".join(str(by_year[y]) for y in by_year)
    years = " / ".join(str(y) for y in by_year)
    partial = ""
    if not window.calendar:
        partial = (" {0} and {1} are partial years: the window runs from "
                   "{2} to {3}.".format(window.lo.year, window.hi.year,
                                        window.lo_iso, window.hi_iso))
    landsat = ""
    if window.hi.year < 1999:
        landsat = (" Landsat 7 does not launch until 1999, so this window rests "
                   "on Landsat 5 alone.")
    return (
        "**3. The years inside the window are not evenly sampled.** Positions "
        "behind the included means, by calendar year: {0} for {1}.{2}{3} The "
        "plain mean leans toward the better-sampled stretch of the window. A "
        "year-balanced mean was offered and not taken (Hannah, 2026-09-22); it "
        "would be a small change to `window_means`."
        .format(counts, years, partial, landsat))


# What went into the line: transects kept, excluded, and why
def write_provenance(df, vertices, folder, window, by_year, min_obs, built_on):
    excluded = df[~df["included"]]
    per_domain = df[df["included"]].groupby("domain_number").size()
    step = np.hypot(np.diff(vertices["x"]), np.diff(vertices["y"]))
    excluded_block = _excluded_table(excluded)
    label = window.label
    imagery_row = ""
    if (folder / "on_imagery").is_dir():
        imagery_row = ("| `on_imagery/` | the line and its ±1 sd band on the USGS "
                       "photographs flown inside the window, at six sites and "
                       "island-wide (three segments, and one ribbon panel); written "
                       "by `coastsat_mean_shoreline_on_imagery.py`, which needs the "
                       "D: drive (see its supporting/ folders) |\n")
    if (folder / "storm_check").is_dir():
        imagery_row += ("| `storm_check/` | were there big storms around this window? The storm "
                        "record 3 yr either side, ranked in 1984-2024, and the mean without "
                        "post-storm passes; written by `coastsat_mean_shoreline_storm_check.py` |\n")

    text = """# mean_shoreline/{label} -- the CoastSat window mean, as a line

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline.py`
on {built_on}.

## What this is

The mean satellite shoreline position over **{described}**, per
CoastSat transect, geolocated and strung into one polyline. It is the
shoreline counterpart of a digitised dune line, and `2-brie-offset` consumes
it the same way.

It is **a window mean, not a survey on a date.** A single satellite pass
carries metres of tide, wave setup and cloud-edge noise; the mean over a
window of passes is the quantity a period can start from.

## The numbers

| | |
|---|---|
| window | {lo} to {hi}, both days included |
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

{window_paragraph}

{sampling_paragraph}

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
| `shoreline_mean_{label}.geojson` | one LineString, EPSG:26918, with the metadata properties a dune line carries |
| `transect_means_{label}.csv` | per transect: n, mean, sd, se, date span, the geolocated point, domain, included/why not |
| `mean_shoreline_{label}.png` | the diagnostic figure |
| `mean_shoreline_{label}_island_outline.png` | panel (a) of the diagnostic alone, the line over the island outline |
{imagery_row}
Resolved through `hat_observed_rates.mean_shoreline_dir/_geojson/_csv`.
Never type these paths.
""".format(
        label=label, built_on=built_on, described=window.described,
        lo=window.lo_iso, hi=window.hi_iso,
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
        window_paragraph=_window_paragraph(window),
        sampling_paragraph=_sampling_paragraph(window, by_year),
        excluded_block=excluded_block, imagery_row=imagery_row,
    )
    (folder / "PROVENANCE.md").write_text(text, encoding="utf-8")


# Run: pick the window, average, geolocate, write the line and figures
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    how = ap.add_mutually_exclusive_group()
    how.add_argument("--window", nargs=2, type=int, metavar=("START", "END"),
                     help="inclusive calendar years to average (the default, 1995 1997)")
    how.add_argument("--window-dates", nargs=2, metavar=("FIRST", "LAST"),
                     help="inclusive ISO dates to average, e.g. 1995-10-12 1997-10-12")
    how.add_argument("--centred-on", choices=sorted(SURVEY_ANCHORS),
                     help="+/- --half-width years of a start DEM's lidar flights")
    ap.add_argument("--half-width", type=int, default=DEFAULT_HALF_WIDTH_YEARS,
                    help="years either side of --centred-on (default 1)")
    ap.add_argument("--min-obs", type=int, default=DEFAULT_MIN_OBS,
                    help="positions a transect needs to contribute (default 10)")
    a = ap.parse_args(argv)
    try:
        if a.centred_on:
            window = Window.centred_on(a.centred_on, a.half_width)
        elif a.window_dates:
            window = Window.from_dates(*a.window_dates)
        else:
            window = Window.from_years(*(a.window or DEFAULT_WINDOW))
    except ValueError as exc:
        ap.error(str(exc))
    built_on = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    lookup = pd.read_csv(transect_lookup()).dropna(subset=["domain_number"])
    print("Window     : {0}".format(window.described))
    print("Transects  : {0} in the domain lookup".format(len(lookup)))
    print("Layer      : {0}".format(TRANSECT_LAYER.name))
    geometry = transect_geometry(set(lookup["transect_id"]))
    print("  geolocated {0} of {1}".format(len(geometry), len(lookup)))

    df, by_year = window_means(lookup, geometry, window, a.min_obs)
    vertices = line_vertices(df)
    print("  {0} transects used, {1} excluded".format(
        len(vertices), len(df) - len(vertices)))
    if len(vertices) < 2:
        sys.exit("fewer than two usable transects; nothing to build a line from")
    print("  positions per transect: median {0}, min {1}".format(
        int(vertices["n_obs"].median()), int(vertices["n_obs"].min())))
    print("  standard error of a transect mean: median {0:.1f} m".format(
        vertices["se_chainage_m"].median()))

    print("  positions by year: {0}".format(by_year))

    folder = mean_shoreline_dir(*window.key)
    folder.mkdir(parents=True, exist_ok=True)
    csv_path = mean_shoreline_csv(*window.key)
    df.to_csv(csv_path, index=False)
    print("\nWrote {0}  ({1} rows)".format(csv_path, len(df)))
    geo_path = mean_shoreline_geojson(*window.key)
    write_geojson(vertices, geo_path, window, built_on, len(df))
    print("Wrote {0}  ({1} vertices)".format(geo_path, len(vertices)))
    figure(df, vertices, folder, window)
    outline_figure(vertices, folder, window)
    write_provenance(df, vertices, folder, window, by_year, a.min_obs, built_on)
    print("Wrote {0}".format(folder / "PROVENANCE.md"))


if __name__ == "__main__":
    main()
