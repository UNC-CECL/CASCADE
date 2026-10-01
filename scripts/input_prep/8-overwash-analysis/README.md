# 8-overwash-analysis — scripts

Figures of the observed overwash record on Hatteras Island, per image and
per CASCADE domain, for the two model periods. The record itself is the
workbook in `data/hatteras_init/8-overwash-analysis/1-observations/`; every
output lands under `data/hatteras_init/8-overwash-analysis/` (see the README
there).

Split into steps 2026-09-22 to match the data tree, which was already
`1-observations/ 2-record/ 3-vs-footprint/` while the code sat flat.

| file | what it does |
|---|---|
| `1-observations/overwash_data.py` | Reads the workbook (observation matrix, storm reference sheet) and decides, from dates, which image first shows each storm. No plotting. Imported by all three of the others, which is why it sits in the step it builds rather than in a folder of its own. |
| `2-record/overwash_heatmap_multiperiod.py` | The heatmap figure, one per period (`period1`, `period2`, `combined`), plus the two tables and `CAPTIONS.md`. |
| `2-record/overwash_map_periods.py` | The island map, one figure per period, with the domains shaded by images with overwash and the date strip beside it; `--both` adds the side-by-side figure. Needs `D:/Hatteras_GIS` (domain boxes, coastline). |
| `3-vs-footprint/overwash_vs_footprint.py` | Sets the overwash seen between the two dune-line frames (Aug 1985 to Oct 1997, nine images) against the rows the 1984 reconstruction adds and removes, per domain; contingency, lists, flags and three figures (alongshore, summary, three-panel island map). Reads the footprint, the road relocation table and the shoreline rates as inputs. |
| `superseded_20260910/` | The May 2026 single-period script. Kept for the record; its paths are dead. |

## Run

```
python 2-record/overwash_heatmap_multiperiod.py                 # all three periods
python 2-record/overwash_heatmap_multiperiod.py period2         # one period
python 2-record/overwash_heatmap_multiperiod.py combined --uniform-years
python 2-record/overwash_map_periods.py                         # one per period
python 2-record/overwash_map_periods.py --both                  # plus side by side
python 3-vs-footprint/overwash_vs_footprint.py                  # vs the 1984 footprint
```

Each script replaces only its own entries in `CAPTIONS.md` (`upsert_caption`
in `1-observations/overwash_data.py`), so the run order does not matter.

## The storm-to-image rule

An image shows a storm if it is the first image taken on or after the
storm's last listed day (`Storm_Reference`, column `Search_GE_After`), with
a 7-day grace because the sheet's date ranges run to dissipation. The one
hand override is the May 2022 nor'easter, routed to the October 2023 image
on Hannah's own note. Storms not in the sheet (Ida 2009, Debby 2024) are
listed in `EXTRA_STORMS`. The decisions are written out to
`1-observations/storms_by_image.csv`, so a change in the sheet is visible there
rather than only in the figure.

Style: `scripts/site_layer/hat_figure_style.py`. No in-image titles or footnotes; the
words are in `CAPTIONS.md`.

## Naming and the sibling imports

No `HAT_` prefix; `overwash_vs_footprint.py` lost its on 2026-09-22, when it
was the only prefixed file of the four. See `../README.md` for which stages
are bare and which are not.

Three scripts import `overwash_data` by name, and `overwash_vs_footprint.py`
also imports the two in `2-record/` for their shared colours and geometry.
Each puts the folders it needs on `sys.path`, anchored on the `pyproject.toml`
found by searching upward — never counted from the file's own depth, so a
script can change step without breaking (rule 5).

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-observations/overwash_data.py

The overwash observation record, read off Hatteras_Overwash_Data.xlsx, in the shape the figures need.

From the script's original header:

```text
The overwash observation record, read off `Hatteras_Overwash_Data.xlsx` in
data/hatteras_init/8-overwash-analysis/1-observations/, in the shape the two
figure scripts need. Nothing is plotted here.

WHAT THE WORKBOOK HOLDS
    Overwash_Matrix   one row per imagery date assessed (28 as of 2026-09-10),
                      one column per CASCADE domain 1..90: 1 = overwash
                      present, 0 = assessed and absent, blank = the image does
                      not cover that domain. Period 1 rows are Hapke and
                      Henderson's delineations; Period 2 rows are Google Earth.
    Storm_Reference   the named storms the imagery search was organised
                      around, with a date range and a "search after" date.

WHAT IS DERIVED HERE, AND WHY
    Whether an image shows a storm is decided from DATES, not from the flags
    the old script carried by hand: the 1985, 1991 and 1992 flags said the
    storms were visible when the image of that year was taken before them
    (Aug 23 1985 vs Gloria in late Sep; Oct 19 1991 vs the Perfect Storm on
    Oct 28; Oct 2 1992 vs the Dec nor'easter). The rule is: the image that
    shows a storm is the first image taken on or after the storm's last day,
    with a 7-day grace because the reference sheet's ranges run to dissipation
    (Emily's range ends Sep 6 1993 but its closest approach was Aug 31, and
    the image is Sep 2; Irene's image is landfall day). One hand override
    remains, Hannah's own routing of the May 2022 nor'easter to the Fall 2023
    image (her note: "first good imagery since").

USAGE
    from overwash_data import load_observations, load_storms, SECTIONS
```

Notes that were in the code:

```text
Every folder is resolved by site_layer/hat_overwash.py (2026-09-18, when the
data folder was regrouped by job). OUT_DIR stays the name the three sibling
scripts import; it is the root, where CAPTIONS.md and the README live.
```

```text
The workbook is the hand-digitised record, so it lives with the data, not
the code (moved 2026-09-10). Edit it there.
```

```text
Matches ANN_TOWN_SPANS in the CASCADE site configuration. "Tri-Village" is
Rodanthe, Waves and Salvo.
```

```text
Images Hannah flagged as poor or partial, by (year, season). The
Image_Quality column is free text ("Medium and cloudy", "Imagery stopped at
domain 17?"), so the judgement is kept here rather than parsed.
```

```text
Storm id -> (year, season) of the image that shows it, where the date rule
is overruled on purpose. See the module docstring.
```

```text
Storms that are not in the reference sheet but were in the old figure
script, kept so the figure does not lose them. End date = last day near NC.
```

```text
Order of the folders in the file, whatever order the scripts ran in: the
order of hat_overwash.CAPTION_FOLDERS, keyed by the label a heading carries.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_observations()`**

```text
(obs, domains, matrix)

obs      DataFrame, one row per image, sorted by Imagery_Date, with
         Obs_ID, Imagery_Date, Year, Season, Image_Quality, Source,
         Linked_Event_ID, poor (bool).
domains  int array, 1..90 in ascending order.
matrix   float (n_obs, 90): 1, 0 or NaN, columns in `domains` order.
```

**`load_storms()`**

```text
List of dicts, one per storm, sorted by end date:
    id, name, cat ('H1'..'H5', 'NE 3'..'NE 5', 'TS', 'ET'),
    end (Timestamp, the reference sheet's Search_GE_After),
    approx (True when the sheet gives only a month), year, month.
```

**`assign_capture()`**

```text
For every storm, the index into `obs` of the image that shows it, or None.

The first image on or after (end - CAPTURE_GRACE_DAYS), CAPTURE_OVERRIDE
winning where set. Prints the rows whose Linked_Event_ID disagrees with
the rule, so a change in the sheet is noticed rather than silently
absorbed.
```

**`upsert_caption()`**

```text
Replace or add the entry for <name> in CAPTIONS.md.

`folder` is the short name ("heatmaps", "map", "vs-footprint"); the
heading carries its path from hat_overwash.CAPTION_FOLDERS.

Every script owns only its own entries; the others are kept as they are,
and the file is re-sorted by folder so the order never depends on which
script ran last.
```

</details>

### 2-record/overwash_heatmap_multiperiod.py

Observed overwash per image and per domain, with the storm record beside it, one figure per period.

From the script's original header:

```text
Observed overwash on Hatteras Island, per image and per CASCADE domain, with
the storm record beside it. One figure per period; all three are written in
one run.

    python overwash_heatmap_multiperiod.py            # period1, period2, combined
    python overwash_heatmap_multiperiod.py period2    # one of them
    python overwash_heatmap_multiperiod.py combined --uniform-years

PANELS
    (a) the observation matrix: one row per image (labelled year and date),
        one column per domain; red = overwash present, white = assessed and
        absent, dotted = the image does not reach that domain, hatched rows =
        no image that year. Runs of years without an image are collapsed to
        one thin row (--uniform-years keeps a row per year).
    (b) how many domains each image shows overwashed, against how many it
        assessed (grey).
    (c) the named storms from the reference sheet, at their year. Bold with a
        filled dot: the image on that row was taken after the storm, so what
        the row shows includes it. Light italic: the image predates the storm
        or there is none that year; the thin line then leads down to the
        first image taken after it, which is where its effects can appear.
    (d) how many images show each domain overwashed.

STYLE
    hat_figure_style; no in-image title or footnote. The words are in
    data/hatteras_init/8-overwash-analysis/CAPTIONS.md, written by this script.

OUTPUT   data/hatteras_init/8-overwash-analysis/ (site_layer/hat_overwash.py)
    2-record/heatmaps/overwash_heatmap_<period>.png
    1-observations/overwash_observations.csv   one row per image and domain
    1-observations/storms_by_image.csv         which image first shows each storm
    CAPTIONS.md

The 2026-05 version of this script had a comparison mode for a modelled
overwash matrix that was never produced. It is gone; a model comparison
should align to 1-observations/overwash_observations.csv.
```

Notes that were in the code:

```text
overwash_data.py is the stage's shared module -- the observation record,
SECTIONS, PERIODS and the loaders. It sits in 1-observations/ because that
is the step it builds. Anchored on REPO, never counted from HERE (rule 5),
so this survives the file changing depth.
```

```text
Storm intensity: hurricanes on a red ramp, nor'easters on a blue ramp,
tropical/extratropical storms grey. (size, colour)
```

```text
A year with several images: sit on the first image after the
storm, else on the last image of the year.
```

```text
A storm's label position, and a small horizontal offset per connector
so the lines leading from several storms to the same image stay apart.
```

```text
the category is the marker (see the key beneath), so the label
carries name and month only: the column is one text width wide.
```

```text
One image row is ROW_IN tall unless the period has too many rows for a
page, when every row shrinks together and the figure still fits. The
storms named on a gap row each need a line of type whatever the rows
come out at, so the two are solved for together rather than fixed.
```

<details><summary>Function notes (the original docstrings)</summary>

**`build_rows()`**

```text
Display rows, top to bottom: dicts with kind ('obs' | 'gap'), label,
years (list), obs (index into `obs` or None), height (row units, set
later once the storms are placed).
```

**`place_storms()`**

```text
Attach each storm in the period to a display row and decide whether that
row's image shows it. Returns {row_index: [storm, ...]} and sets
storm['matched'], storm['capture_row'].
```

**`set_heights()`**

```text
Row heights in row units. A gap row holds `slot` units per storm named
on it; `slot` is set by make_figure so that one storm gets one line of
type however far the rows have had to be squeezed.
```

**`draw_period_bars()`**

```text
Two thin columns, Period 1 over 1984–2004 and Period 2 over 2004–2024.
They overlap at the 2004 row, which is the last image of one and the
first of the other.
```

</details>

### 2-record/overwash_map_periods.py

Where on Hatteras Island overwash was observed, and when, for the two model periods.

From the script's original header:

```text
Where on Hatteras Island overwash was observed, and when, for the two model
periods, drawn on the island outline.

    python overwash_map_periods.py            # one figure per period
    python overwash_map_periods.py --both     # plus the two-period figure

PANELS
    (a) the island (NC 1:80k coastline) for the period, with the
        90 CASCADE domain boxes, each clipped to land and shaded by the number
        of images in the period that show it overwashed. NC-12 as a thin dark
        line. Domain numbers every ten on the ocean side, reach names on the
        sound side.
    (b) [and (d)]: one column per image assessed in the period, in date order,
        registered to the same northing as the island beside it, so a purple
        cell sits level with the domain it belongs to. Column heads (under
        the strip) carry the image date and the named storm(s) that image is
        the first to show.

INPUTS
    hat_map_layers.DOMAIN_BOXES                     the domain boxes (EPSG:3725)
    hat_map_layers.NC_COAST                         the coastline, NC 1:80k clipped
    (both in the repository since 2026-09-18; they were read off D:/Hatteras_GIS)
    data/hatteras_init/4-mgmt-forcing/road_offset/raw_offset/2008/nc12_2008.geojson
    data/hatteras_init/8-overwash-analysis/1-observations/Hatteras_Overwash_Data.xlsx
    via overwash_data.py

OUTPUT
    data/hatteras_init/8-overwash-analysis/2-record/map/overwash_map_period1.png (+ .pdf)
    data/hatteras_init/8-overwash-analysis/2-record/map/overwash_map_period2.png (+ .pdf)
    (overwash_map_periods.png with --both) and their entries in CAPTIONS.md.

STYLE
    hat_figure_style, drawn at the double-column width (figsize("double")):
    the island panel takes the height the width allows, capped at a page.
    Overwash is C["ACCENT"]; the island's count classes are the matching
    Purples ramp so the two panels read as one colour.

The domain boxes live on the external drive; the script stops with a
message if the drive is not there rather than drawing without them.
```

Notes that were in the code:

```text
overwash_data.py is the stage's shared module -- the observation record,
SECTIONS, PERIODS and the loaders. It sits in 1-observations/ because that
is the step it builds. Anchored on REPO, never counted from HERE (rule 5),
so this survives the file changing depth.
```

```text
Two periods do not leave a column two lines of rotated type wide, so
their heads run as one line and the band below the strips grows instead.
```

```text
The island is as tall as the page allows; if that leaves the image
columns thinner than their own labels, the islands give width back.
```

<details><summary>Function notes (the original docstrings)</summary>

**`count_cmap()`**

```text
Discrete purples for 1..vmax, the family of C["ACCENT"]; 0 is drawn
as bare land.
```

**`draw_island()`**

```text
The island with the domain boxes clipped to land.

fill_of        {domain: colour}; domains not in it are drawn as bare land.
alpha_of       {domain: alpha}, optional, for a paler fill on flagged domains.
pad_*          metres of water to leave around `bounds`; pad_w defaults to
               PAD_W (room for the reach names) or 1500 without them.
label_every    domain numbers at 1 and every N.
reach_rotation 90 writes the reach names along the sound side, which is
               what a zoomed section needs, where the panel is narrow.
letter         the panel letter, "a", "b", ...; drawn bold at the left of
               the title by hat_figure_style._title.
Only the domains in `dom` are drawn and labelled, so pass a subset to
draw a section.
```

**`scalebar_and_north()`**

```text
The house scale bar and north arrow (hat_figure_style), for a map
panel without coordinate ticks.
```

**`strip_heads()`**

```text
The column-head strings of the strip: date (* poor image) and, above
it, the named storm the image is the first to show (+n more). `one_line`
joins the two, for a strip whose columns are narrower than two lines of
rotated type.
```

**`draw_strips()`**

```text
One column per image; rows are the domain boxes' northing spans. The
column heads hang below the strip, reading upward into it.
```

**`render()`**

```text
One figure at the double-column width: for each period tag an island
panel and its image strip, side by side. With two periods the islands
give up width to the strips, and the reach names go with it.
```

**`main()`**

```text
Default: one figure per period (overwash_map_period1.png,
overwash_map_period2.png). `--both` adds the two-period figure
(overwash_map_periods.png), the two island-and-strip pairs side by side.
The shade scale is shared either way.
```

</details>

### 3-vs-footprint/overwash_vs_footprint.py

Does the observed overwash of 1984-1997 line up with the rows the 1984 reconstruction adds and removes?

From the script's original header:

```text
Does the observed overwash of 1984–1997 line up with the rows the 1984
reconstruction adds and removes?

    python overwash_vs_footprint.py

THE TWO RECORDS
    The footprint (2-domain-reconstruction-1984/2-extent/footprint_1984_by_domain.csv)
    gives every domain the median shift between the digitised 1984 and 1997
    dune lines, positive where the 1984 line lay seaward, and the row count the
    10 m rule keeps of it: rows ADDED where the dune retreated over the window,
    rows REMOVED where it advanced.

    The overwash record (1-observations/Hatteras_Overwash_Data.xlsx) gives every
    domain, per image, whether washover was visible.

THE WINDOW (Hannah, 2026-09-10: strictly between the line dates)
    Both dune lines were digitised from the same USGS photographs the overwash
    record uses: the 1984 line from the 19 Sep 1984 frame, the 1997 line from
    the 12 Oct 1997 frame. So the images that can speak to the shift are the
    ones taken AFTER the 1984 frame and UP TO AND INCLUDING the 1997 frame:
    nine images, Aug 1985 to Oct 1997. The Diana image is out (the 1984 line
    was drawn on it, so that overwash predates the line); the Bonnie image is
    out (it postdates the 1997 line).

WHAT COUNTS AS AGREEMENT (Hannah: show both readings, do not blame absence)
    Physically, overwash flattens the dune and pushes the vegetation break
    landward, so overwash should go with rows ADDED. Overwash in a rows-REMOVED
    domain is the disagreement worth a look. A rows-added domain with NO
    overwash in the window is listed as unexplained, not as a mismatch: the
    nine images are two to four years apart and washover fades from imagery
    within a few years, so absence is weak evidence. Both readings are
    reported: "given overwash, which action?" and "given the action, was
    there overwash?".

FLAGS (kept, not dropped)
    SPREAD_STRADDLES_ZERO from the footprint (the p10–p90 of the shift
    crosses zero); the erosion hotspots and jetties from domains.geojson;
    NC-12 relocated 1984–2004 (road_relocation_1978_2008.csv); and shoreline
    erosion faster than ERODE_THRESH by the CoastSat 1984–2004 LRR, with the
    DSAS 1978–1997 mean rate carried beside it because it sits closer to the
    window. Erosion retreats the dune line without any overwash, which is
    why the last flag exists.

OUTPUT   data/hatteras_init/8-overwash-analysis/3-vs-footprint/tables/
    overwash_vs_footprint_by_domain.csv    the joined table, one row per domain
    overwash_vs_footprint_contingency.csv  overwashed x action, all and unflagged
    overwash_vs_footprint_summary.txt      the readings in words, with the lists
    ../  (3-vs-footprint/, the figures)
        overwash_vs_footprint_alongshore.png   images, footprint bars, flags, by domain
        overwash_vs_footprint_summary.png      shift by overwash status; share overwashed per action
        overwash_vs_footprint_map.png          three alongshore sections, zoomed, each with
                                               the three layers: overwash, footprint, reading
        overwash_vs_footprint_map_island.png   the whole island, the same three layers
    The map reads the domain boxes and coastline from the repository, through
    overwash_map_periods.load_geometry (off the D: drive since 2026-09-18).
```

Notes that were in the code:

```text
overwash_data.py is the stage's shared module -- the observation record,
SECTIONS, PERIODS and the loaders. It sits in 1-observations/ because that
is the step it builds. Anchored on REPO, never counted from HERE (rule 5),
so this survives the file changing depth.
This one draws on TWO earlier steps: the record from 1-observations and the
shared colours and geometry from the 2-record figures. Both folders go on
the path, anchored on REPO (rule 5), never counted from this file.
```

```text
The repository copy (identical to D:/Hatteras_GIS/domains.geojson in
geometry and every attribute, hotspot and armor included; 2026-09-18).
```

```text
The footprint bars take the LIGHTER RdBu pair so they cannot be read as the
overwash accent red of panel (a); the overwash marker in (b) is the accent.
```

```text
(c) flags
The straddle is not a row here: it is already the hatching in (b).
```

```text
one key under the figure: the panels are full of data and the village
names have the top of (a)
```

```text
the panel is one column wide: the letter alone goes above it and
the legend beneath names the layer
```

```text
one legend per layer, side by side under the figure: the panels are a
single column wide, too narrow to carry a legend under each of them
```

<details><summary>Function notes (the original docstrings)</summary>

**`fig_alongshore()`**

```text
Three stacked panels on one domain axis, at the double-column width:
the images, the footprint, the flags.
```

</details>

### 4-vs-model/overwash_vs_model.py

Does the model overwash where and when the imagery shows washover?

From the script's original header:

```text
Does the model overwash where and when the imagery shows washover? The
observed presence/absence record (1-observations/overwash_observations.csv)
against the overwash Barrier3D records in the hindcast runs, for the two
current windows, 1996-2010 and 2010-2024.

    python scripts/input_prep/8-overwash-analysis/4-vs-model/overwash_vs_model.py

Writes to data/hatteras_init/8-overwash-analysis/4-vs-model/ (figures) and
tables/ inside it (hat_overwash.VS_MODEL, VS_MODEL_TABLES).

HOW AN IMAGE IS MATCHED TO THE MODEL
    An image shows the washover left by the storms since the previous image.
    So each image is compared with the model's overwash between the previous
    image's date and its own (the window is cut at the run's start, and such
    images are flagged "partial"). Images after the last modelled storm year
    are left out.

    Barrier3D records overwash per domain per MODEL YEAR (QowTS, m3/m), not
    per storm. A model year's overwash is dated by that year's LARGEST storm
    (the storm file's EndTime), because every storm in a year is tested
    against the same dune crest (storm_replay.py): if any storm overwashed a
    domain, the largest one did. Smaller storms that year may also have, so
    `model_any_storm` is the upper bound that dates the year's overwash by
    every storm in it.

    Storm dates take the same 7-day grace the observed record's storm table
    uses (overwash_data.assign_capture): a storm counts toward the first image
    dated on or after its last hour above the berm LESS 7 days, because a
    storm's tail runs to dissipation and image dates are approximate. Without
    it, the model's Irene (ends 2011-08-28) falls one day after the 2011-08-27
    image that shows Irene's washover (storms_by_image.csv maps it there).

    A cell counts as modelled overwash when QowTS > the threshold; the
    headline uses 0 (any overwash), and tables/ carries 1 and 5 m3/m too.

WHAT IS AND IS NOT A MISMATCH
    Washover fades and images are months to years apart, so "observed 0,
    model 1" can be real washover the image no longer shows. Both directions
    are reported; neither is discounted.

THE RUNS
    The managed runs are the comparison (the imagery is the managed island):
    edgeBE, road + beach/dune manager (+ fills in 2010-2024), no groin. The
    natural runs (no road, no manager) are in the tables for contrast.
```
