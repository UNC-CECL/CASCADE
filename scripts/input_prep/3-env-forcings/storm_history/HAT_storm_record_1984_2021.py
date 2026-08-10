"""
HAT_storm_record_1984_2021.py

Named-storm record for Hatteras Island, NC, 1984-2021.

Data sources:
  - NOAA Historical Hurricane Tracks search (60 nm buffer of Hatteras, Dare
    County, NC) -> "start"/"end"/"cat" per HISTORICAL_STORMS below, plus
    max sustained wind speed / min pressure pulled from the same search.
  - "Hatteras, NC Hurricane History Since 1985" -> local landfall/impact
    detail (surge, evacuations, damage) for the more notable storms.
  - Bertha, Fran (1996) and Isaias (2020): NOT in the 60 nm buffer search
    (their tracks passed outside it), but documented as local impacts in
    the hurricane-history doc. Wind/pressure for these three were cross-
    checked against NHC/NWS tropical cyclone reports (source="local"
    below). This distinction is kept in the data for provenance but is no
    longer flagged on the figure itself (previously a dagger + footnote).

Known remaining gaps (not yet added -- flagged to Hannah, not included
without confirmation): Earl 2010, Sandy 2012, Maria 2017, Florence 2018,
and Michael 2018 all appear in the hurricane-history doc but, like
Bertha/Fran/Isaias, are absent from the 60 nm buffer search. Unlike
Bertha/Fran/Isaias their NC landfalls (or, for Sandy/Michael, their
tracks) were far enough from Hatteras that local impact was minor/indirect
per the doc's own text, so they were left out pending a decision on
whether to include them.

Nor'easters / non-tropical winter storms are still not included (no
dataset available for these).

Chart labels: storms with a documented local-impact note get a plain
asterisk (*); the writeup is in the "Storm Details" sidebar, matched by
name/year rather than a numbered reference (numbers next to the year
abbreviation, e.g. Gloria '85, read as confusingly similar to the year
itself).
"""

import datetime as dt

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox

# --------------------------------------------------------------------------
# Storm data
# --------------------------------------------------------------------------
# cat        = max Saffir-Simpson category reached over the storm's full
#              lifetime (NHC best track), not necessarily at landfall
# wind       = max sustained wind speed (mph) corresponding to `cat`
# pressure   = min central pressure (mb) corresponding to `cat`
# landfall   = True  -> doc confirms direct Hatteras/Dare Co. impact or evac
#              False -> doc/NHC indicates storm passed by with minor/no
#                       local impact (or made landfall well south)
#              None  -> no local-impact narrative available (NOAA-only)
# source     = "noaa"  -> from the 60 nm buffer search
#              "local" -> from the hurricane-history doc only; wind/
#                         pressure cross-checked against NHC reports
HISTORICAL_STORMS = [
    # -------------------------------------------------------------------------
    # 1984-2004  (Period 1: calibration)
    # -------------------------------------------------------------------------
    {"name": "Diana",     "start": "1984-09-08", "end": "1984-09-16", "cat": "H4", "wind": 115, "pressure": 949,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Gloria",    "start": "1985-09-16", "end": "1985-10-02", "cat": "H4", "wind": 125, "pressure": 920,  "landfall": True,  "note": "Direct hit \u2014 Cat 2 at landfall, 6\u20138 ft surge", "source": "noaa"},
    {"name": "Kate",      "start": "1985-11-15", "end": "1985-11-23", "cat": "H3", "wind": 105, "pressure": 954,  "landfall": False, "note": "",                                                  "source": "noaa"},
    {"name": "Charley",   "start": "1986-08-13", "end": "1986-08-30", "cat": "H1", "wind": 70,  "pressure": 980,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Alberto",   "start": "1988-08-05", "end": "1988-08-08", "cat": "TS", "wind": 35,  "pressure": 1002, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bob",       "start": "1991-08-16", "end": "1991-08-29", "cat": "H3", "wind": 100, "pressure": 950,  "landfall": False, "note": "",                                                  "source": "noaa"},
    {"name": "Danielle",  "start": "1992-09-22", "end": "1992-09-26", "cat": "TS", "wind": 55,  "pressure": 1001, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Emily",     "start": "1993-08-22", "end": "1993-09-06", "cat": "H3", "wind": 100, "pressure": 960,  "landfall": False, "note": "Dare Co. evacuated, $12M damage",                    "source": "noaa"},
    {"name": "Allison",   "start": "1995-06-03", "end": "1995-06-11", "cat": "H1", "wind": 65,  "pressure": 982,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Arthur",    "start": "1996-06-17", "end": "1996-06-23", "cat": "ET", "wind": 45,  "pressure": 992,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bertha",    "start": "1996-07-05", "end": "1996-07-17", "cat": "H3", "wind": 115, "pressure": 960,  "landfall": True,  "note": "Landfall nr. Topsail; Cat 2, 104 mph, 5 ft surge",   "source": "local"},
    {"name": "Fran",      "start": "1996-08-23", "end": "1996-09-10", "cat": "H3", "wind": 120, "pressure": 946,  "landfall": False, "note": "Landfall Cape Fear; no Dare Co. evacuation",         "source": "local"},
    {"name": "Josephine", "start": "1996-10-04", "end": "1996-10-16", "cat": "TS", "wind": 60,  "pressure": 970,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bonnie",    "start": "1998-08-19", "end": "1998-08-31", "cat": "H3", "wind": 100, "pressure": 954,  "landfall": True,  "note": "Dare Co. evacuated, 6\u20138 ft surge",              "source": "noaa"},
    {"name": "Dennis",    "start": "1999-08-24", "end": "1999-09-08", "cat": "H2", "wind": 90,  "pressure": 962,  "landfall": True,  "note": "Dare Co. evacuated, $10M damage",                    "source": "noaa"},
    {"name": "Irene",     "start": "1999-10-12", "end": "1999-10-19", "cat": "H2", "wind": 95,  "pressure": 958,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Arthur",    "start": "2002-07-14", "end": "2002-07-19", "cat": "TS", "wind": 50,  "pressure": 992,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Gustav",    "start": "2002-09-08", "end": "2002-09-15", "cat": "H2", "wind": 85,  "pressure": 960,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Kyle",      "start": "2002-09-20", "end": "2002-10-12", "cat": "H1", "wind": 75,  "pressure": 980,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Isabel",    "start": "2003-09-06", "end": "2003-09-20", "cat": "H5", "wind": 145, "pressure": 915,  "landfall": True,  "note": "Breached island; $167M damage",                      "source": "noaa"},
    # -------------------------------------------------------------------------
    # 2004-2024  (Period 2: validation; data currently through 2021)
    # -------------------------------------------------------------------------
    {"name": "Alex",      "start": "2004-07-31", "end": "2004-08-06", "cat": "H3", "wind": 105, "pressure": 957,  "landfall": True,  "note": "Sound-side flooding, $2.4M damage",                  "source": "noaa"},
    {"name": "Bonnie",    "start": "2004-08-03", "end": "2004-08-14", "cat": "TS", "wind": 55,  "pressure": 1001, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Ophelia",   "start": "2005-09-06", "end": "2005-09-23", "cat": "H1", "wind": 75,  "pressure": 976,  "landfall": True,  "note": "Hatteras Island evacuated",                          "source": "noaa"},
    {"name": "Alberto",   "start": "2006-06-10", "end": "2006-06-19", "cat": "TS", "wind": 60,  "pressure": 969,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Barry",     "start": "2007-05-31", "end": "2007-06-05", "cat": "TS", "wind": 50,  "pressure": 990,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Gabrielle", "start": "2007-09-08", "end": "2007-09-11", "cat": "TS", "wind": 50,  "pressure": 1004, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Cristobal", "start": "2008-07-19", "end": "2008-07-23", "cat": "TS", "wind": 55,  "pressure": 998,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Irene",     "start": "2011-08-21", "end": "2011-08-30", "cat": "H3", "wind": 105, "pressure": 942,  "landfall": True,  "note": "Dare Co. evacuated; significant flooding, $54M damage", "source": "noaa"},
    {"name": "Beryl",     "start": "2012-05-25", "end": "2012-06-02", "cat": "TS", "wind": 60,  "pressure": 992,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Arthur",    "start": "2014-06-28", "end": "2014-07-09", "cat": "H2", "wind": 85,  "pressure": 973,  "landfall": True,  "note": "Earliest hurricane on record; Hatteras evacuated",   "source": "noaa"},
    {"name": "Claudette", "start": "2015-07-12", "end": "2015-07-15", "cat": "TS", "wind": 45,  "pressure": 1003, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Bonnie",    "start": "2016-05-27", "end": "2016-06-09", "cat": "TS", "wind": 40,  "pressure": 1006, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Colin",     "start": "2016-06-05", "end": "2016-06-08", "cat": "ET", "wind": 50,  "pressure": 987,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Hermine",   "start": "2016-08-28", "end": "2016-09-08", "cat": "H1", "wind": 70,  "pressure": 981,  "landfall": True,  "note": "Villages flooded; 4 ft surge, $5.4M damage",         "source": "noaa"},
    {"name": "Julia",     "start": "2016-09-13", "end": "2016-09-21", "cat": "TS", "wind": 45,  "pressure": 1007, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Matthew",   "start": "2016-09-28", "end": "2016-10-10", "cat": "H5", "wind": 145, "pressure": 934,  "landfall": False, "note": "Landfall SC; Dare Co. gusts to 94 mph",              "source": "noaa"},
    {"name": "Unnamed",   "start": "2017-08-27", "end": "2017-08-29", "cat": "ET", "wind": 40,  "pressure": 1004, "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Dorian",    "start": "2019-08-24", "end": "2019-09-09", "cat": "H5", "wind": 160, "pressure": 910,  "landfall": True,  "note": "Cat 1 landfall on Hatteras; 4\u20137 ft surge, $14.8M damage", "source": "noaa"},
    {"name": "Arthur",    "start": "2020-05-16", "end": "2020-05-21", "cat": "ET", "wind": 55,  "pressure": 989,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Fay",       "start": "2020-07-05", "end": "2020-07-11", "cat": "TS", "wind": 50,  "pressure": 998,  "landfall": None,  "note": "",                                                  "source": "noaa"},
    {"name": "Isaias",    "start": "2020-07-30", "end": "2020-08-05", "cat": "H1", "wind": 90,  "pressure": 986,  "landfall": True,  "note": "72 mph in Avon; 1\u20133 ft surge, Zone A evacuated", "source": "local"},
    {"name": "Claudette", "start": "2021-06-17", "end": "2021-06-23", "cat": "TS", "wind": 40,  "pressure": 1003, "landfall": None,  "note": "",                                                  "source": "noaa"},
]

PERIOD_BOUNDARY = 2004  # Period 1 (calibration) / Period 2 (validation)


def decimal_year(date_str):
    y, m, d = (int(x) for x in date_str.split("-"))
    start = dt.date(y, 1, 1)
    this = dt.date(y, m, d)
    end = dt.date(y + 1, 1, 1)
    return y + (this - start).days / (end - start).days


for s in HISTORICAL_STORMS:
    s["x"] = decimal_year(s["start"])

# --------------------------------------------------------------------------
# De-clutter: several years have 2-5 storms only weeks apart (e.g. five in
# 2016). Rather than a running left-to-right push (which can cascade a
# whole cluster's positions into a neighboring year), group consecutive
# storms whose gap is below MIN_SEP and evenly redistribute *within* each
# cluster, centered on the cluster's true mean date. Labels still show each
# storm's own year, so nothing displayed becomes inaccurate.
storms = sorted(HISTORICAL_STORMS, key=lambda s: s["x"])
MIN_SEP = 0.16
clusters, current = [], [storms[0]]
for s in storms[1:]:
    if s["x"] - current[-1]["x"] < MIN_SEP:
        current.append(s)
    else:
        clusters.append(current)
        current = [s]
clusters.append(current)

for cluster in clusters:
    n = len(cluster)
    if n == 1:
        continue
    mean_x = sum(c["x"] for c in cluster) / n
    for i, c in enumerate(cluster):
        c["x"] = mean_x + (i - (n - 1) / 2) * MIN_SEP

# Cluster-level spreading can leave the tail of one cluster close to the
# head of the next (e.g. Colin '16 vs Hermine '16). Clean up any residual
# close pairs with a small, local nudge.
storms.sort(key=lambda s: s["x"])
CLEANUP_SEP = 0.10
for i in range(1, len(storms)):
    if storms[i]["x"] - storms[i - 1]["x"] < CLEANUP_SEP:
        storms[i]["x"] = storms[i - 1]["x"] + CLEANUP_SEP

# --------------------------------------------------------------------------
# Styling
# --------------------------------------------------------------------------
CAT_COLORS = {
    "H5": "#7a1f1f", "H4": "#c0392b", "H3": "#e67e22", "H2": "#f2a541",
    "H1": "#f1c40f", "TS": "#2f7fb8", "ET": "#95a5a6",
}
CAT_ORDER = ["H5", "H4", "H3", "H2", "H1", "TS", "ET"]
CAT_LABELS = {
    "H5": "Cat. 5", "H4": "Cat. 4", "H3": "Cat. 3", "H2": "Cat. 2",
    "H1": "Cat. 1", "TS": "Trop. Storm", "ET": "Extratropical",
}

plt.rcParams["font.family"] = "DejaVu Sans"
fig = plt.figure(figsize=(30, 12.5), dpi=150)
gs = fig.add_gridspec(1, 2, width_ratios=[3.5, 1], wspace=0.02)
ax = fig.add_subplot(gs[0, 0])
ax_side = fig.add_subplot(gs[0, 1])
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

x_min, x_max = 1983.3, 2022.6

# Year background stripes
for yr in range(1984, 2023):
    if yr % 2 == 0:
        ax.axvspan(yr - 0.5, yr + 0.5, color="#f7ead6", alpha=0.5, zorder=0)

# Period 1 / Period 2 boundary
ax.axvline(PERIOD_BOUNDARY, color="#4a4a4a", linestyle=":", linewidth=1.4, zorder=1)
ax.text(PERIOD_BOUNDARY - 0.15, 168, "Period 1 (calibration)", ha="right", va="bottom",
        fontsize=13, color="#4a4a4a", style="italic")
ax.text(PERIOD_BOUNDARY + 0.15, 168, "Period 2 (test)", ha="left", va="bottom",
        fontsize=13, color="#4a4a4a", style="italic")

# Storm markers
for s in HISTORICAL_STORMS:
    size = 170 + s["wind"] * 2.1
    s["_size"] = size
    ax.scatter(s["x"], s["wind"], s=size, marker="D",
               color=CAT_COLORS[s["cat"]], edgecolor="#2b2b2b",
               linewidth=1.6 if s["landfall"] else 0.8, zorder=3)

# Axes formatting (set BEFORE label placement so pixel<->data mapping used
# for collision checks below matches the final rendered chart)
ax.set_xlim(x_min, x_max)
ax.set_ylim(0, 172)
ax.set_xticks(range(1984, 2022, 2))
ax.set_xticklabels(range(1984, 2022, 2), rotation=30, ha="right", fontsize=12)
ax.set_yticks(range(0, 161, 20))
ax.tick_params(axis="y", labelsize=12)
ax.set_xlabel("Year", fontsize=14, labelpad=12)
ax.set_ylabel("Max Sustained Wind Speed  (mph)", fontsize=14, labelpad=12)
ax.set_title("Storm Record  \u00b7  Hatteras Island  \u00b7  1984\u20132021",
             fontsize=24, fontweight="bold", pad=34, color="#1c2b39")

ax.grid(axis="y", color="#e2e2e2", linewidth=0.8, zorder=0)
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)
for spine in ["left", "bottom"]:
    ax.spines[spine].set_color("#888888")

# Legend (added now so tight_layout below accounts for it)
handles = [Line2D([0], [0], marker="D", color="w", markerfacecolor=CAT_COLORS[c],
                   markeredgecolor="#2b2b2b", markersize=20, label=CAT_LABELS[c])
           for c in CAT_ORDER]
legend = ax.legend(handles=handles, title="Storm Category (max intensity, NHC):",
                    loc="upper center", bbox_to_anchor=(0.5, -0.12),
                    ncol=7, frameon=False, fontsize=17,
                    title_fontsize=18, columnspacing=2.2, handletextpad=0.7)

plt.tight_layout()
fig.canvas.draw()  # finalize layout so pixel positions below are accurate
renderer = fig.canvas.get_renderer()

# --------------------------------------------------------------------------
# Collision-aware label placement
#
# Rather than assuming above/below alternation is enough (it isn't once 3+
# storms land within the same year), place each label, measure its actual
# rendered bounding box, and check it against every box placed so far
# (marker diamonds + other labels). If it collides, try the next vertical
# tier out. Names are placed most-intense-storm-first, so e.g. Isabel/
# Dorian keep their close, prominent labels.
#
# Notes no longer live on the chart at all -- they're too long to fit
# without crowding, even with tiering. Instead, storms with a note get a
# plain asterisk, and readers find the matching writeup by name in the
# chronological sidebar list (no number-matching needed).
# --------------------------------------------------------------------------
placed_boxes = []


def px_bbox_from_data(x, y, size_pts2):
    disp = ax.transData.transform((x, y))
    r = (size_pts2 ** 0.5) * fig.dpi / 72.0 / 2.0
    return Bbox.from_bounds(disp[0] - r, disp[1] - r, 2 * r, 2 * r)


def collides(bbox, pad=2.0):
    bbox = bbox.padded(pad)
    return any(bbox.overlaps(b) for b in placed_boxes)


# Register every marker diamond as a fixed obstacle first
for s in HISTORICAL_STORMS:
    placed_boxes.append(px_bbox_from_data(s["x"], s["wind"], s["_size"]))

LABEL_TIERS = [16, -16, 32, -32, 48, -48, 64, -64, 80, -80, 96, -96, 112, -112]

for s in sorted(HISTORICAL_STORMS, key=lambda s: -s["wind"]):
    yr_suffix = "'" + s["start"][2:4]
    star = "*" if s["note"] else ""
    label = f"{s['name']} {yr_suffix}{star}"
    weight = "bold" if s["landfall"] else "normal"
    fontsize = 15 if s["wind"] >= 140 else 11.5

    chosen_dy = None
    for dy in LABEL_TIERS:
        va = "bottom" if dy > 0 else "top"
        txt = ax.annotate(label, (s["x"], s["wind"]), xytext=(0, dy),
                           textcoords="raw_offset points", ha="center", va=va,
                           fontsize=fontsize, fontweight=weight, color="#2b2b2b")
        bbox = txt.get_window_extent(renderer=renderer)
        if not collides(bbox):
            placed_boxes.append(bbox)
            chosen_dy = dy
            break
        txt.remove()
    if chosen_dy is None:
        dy = LABEL_TIERS[0]
        txt = ax.annotate(label, (s["x"], s["wind"]), xytext=(0, dy),
                           textcoords="raw_offset points", ha="center", va="bottom",
                           fontsize=fontsize, fontweight=weight, color="#2b2b2b")
        placed_boxes.append(txt.get_window_extent(renderer=renderer))

# --------------------------------------------------------------------------
# Sidebar: chronological detail list for every storm with a note.
#
# Positions are anchored to axes-fraction (0, 1) -- the top-left of the
# panel -- with a running raw_offset in points, so the physical panel-pixel
# height doesn't depend on any data/ylim choice. Each entry's vertical
# footprint is measured directly from its actual rendered bounding box
# (same renderer-based approach used for the main chart), so header/note
# spacing is exact rather than a guessed line-height.
# --------------------------------------------------------------------------
import textwrap

ax_side.axis("off")
ax_side.set_xlim(0, 1)
ax_side.set_ylim(0, 1)

MONTH_ABBR = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun",
              "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

noted_storms = [s for s in HISTORICAL_STORMS if s["note"]]
WRAP_WIDTH = 60
cursor_pt = 6.0  # running raw_offset, in points, down from the top of the panel


def place_side(text, fontsize, gap_before, weight="normal", style="normal",
                color="#1c2b39", indent=0.0, linespacing=1.3):
    global cursor_pt
    cursor_pt += gap_before
    txt = ax_side.annotate(text, xy=(indent, 1), xycoords="axes fraction",
                            xytext=(0, -cursor_pt), textcoords="raw_offset points",
                            va="top", ha="left", fontsize=fontsize, fontweight=weight,
                            style=style, color=color, linespacing=linespacing)
    bbox = txt.get_window_extent(renderer=renderer)
    cursor_pt += bbox.height * 72.0 / fig.dpi
    return txt


place_side("Storm Details", fontsize=17, gap_before=0, weight="bold")
place_side("(* on the chart)", fontsize=10.5, gap_before=2, style="italic", color="#6a6a6a")

for s in noted_storms:
    y, m, d = (int(v) for v in s["start"].split("-"))
    date_str = f"{MONTH_ABBR[m]} {y}"
    yr_suffix = "'" + s["start"][2:4]
    header = f"\u2022 {s['name']} {yr_suffix} \u2014 {date_str}"
    place_side(header, fontsize=12, gap_before=15, weight="bold")

    wrapped = textwrap.fill(s["note"], width=WRAP_WIDTH)
    place_side(wrapped, fontsize=10.5, gap_before=3, style="italic",
               color="#5a5a5a", indent=0.03, linespacing=1.4)

panel_height_pt = ax_side.get_position().height * fig.get_size_inches()[1] * 72.0
print(f"Sidebar content height: {cursor_pt:.0f}pt / panel height: {panel_height_pt:.0f}pt")

plt.savefig(r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\storm_history\HAT_storm_record_1984_2021.png",
            dpi=200, bbox_inches="tight")
plt.savefig(r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\storm_history\HAT_storm_record_1984_2021.pdf",
            bbox_inches="tight")
print("Saved figure.")
