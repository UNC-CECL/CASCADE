# 4-comparisons — one source against another

Answers, not model inputs. Nothing the model reads is produced here; these
scripts take the products in `../3-rates/` and put two of them side by side.
**They read the stored tables and never recompute a rate** — if a comparison
disagrees with a rate figure, the rate product is what to fix.

One folder per question.

```
shoreline_vs_duneline/   Does the digitized dune line move with the CoastSat
                         shoreline?
dsas_vs_coastsat/        Do the two rate sources agree?
duneline_positions/      Where did the dune lines actually sit? (positions,
                         not change)
```

## shoreline_vs_duneline

Both sides are net change in metres, over the same windows, with the same sign
convention (seaward positive). That symmetry was the point of rebuilding the
dune side as an endpoint product in September 2026.

```
coastsat_vs_duneline.py
    The base comparison: 3-rates/duneline/endpoint/ against
    3-rates/coastsat/endpoint/, per window. It is also the module the other
    three import — load_chainage, the survey-date table, the beach-width
    shading — so it is the one to read first.

total_change_vs_duneline.py
    Total shoreline change (the LRR over the window it was fitted on, turned
    into metres) against the dune line's measured net change, over 1996-2010,
    2010-2024 and 1996-2024. Absorbed projected_vs_duneline.py on 2026-09-21
    when the vocabulary was settled: that script's numbers were already here
    as the *_dune_interval_m columns.

net_change_vs_duneline.py
    The same pairing over 1996-2024 and its two halves.

smoothed_lowess7_vs_duneline.py
    Both curves through the model target's alongshore LOWESS at 7 domains
    (3.5 km), in the two readings of the shoreline side, so the effect of the
    smoother is visible rather than assumed.
```

## dsas_vs_coastsat

Two scripts because the CoastSat side can be anchored two ways, and they
answer different questions:

```
dsas_vs_coastsat_raw.py
    Calendar windows, no smoothing: the raw per-domain mean LRR from each
    source. Separate from 6-scr-smooth/dsas_vs_coastsat/ on purpose — that
    folder exists to argue about the smoothing, and every figure in it draws
    a LOWESS.

dsas_vs_coastsat_datematched.py
    The CoastSat side anchored on the shoreline SURVEY DATES instead: the mean
    position within +/-W days of each date, end minus start, at both +/-30
    days and +/-6 months. Not an OLS over a window.
```

## duneline_positions

```
duneline_positions.py
    Where the 1997, 2009 and 2023 dune lines sat — the lines standing for
    model years 1996, 2010 and 2024. Overview maps in three north-up segments,
    imagery zooms, dune-line-to-NC-12 distance, and beach width.
```

Two traps that folder's outputs record: a transect extended landward can cross
the road, and a seaward-drawn road line makes the dune-to-road distance
negative. Both are described in the data-tree README beside the figures.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### dsas_vs_coastsat/dsas_vs_coastsat_datematched.py

DSAS against CoastSat, with the CoastSat side anchored on the shoreline survey dates.

From the script's original header:

```text
DSAS against CoastSat with the CoastSat side anchored on the SHORELINE SURVEY
DATES rather than on a calendar window (Hannah, 2026-09-22, "do both": at
+/-30 days and at +/-6 months).

    4-comparisons/dsas_vs_coastsat/survey_dates/dsas_vs_coastsat_datematched.png
    4-comparisons/dsas_vs_coastsat/survey_dates/slides/ (--slide)
    4-comparisons/dsas_vs_coastsat/survey_dates/supporting/*.csv

THE METHOD
    Not an OLS over a window. For each CoastSat transect, the mean shoreline
    position within +/-W days of the START survey date is subtracted from the
    mean within +/-W days of the END date and divided by the interval -- the
    same endpoint method coastsat_endpoint.py uses against the dune line, run
    through the same endpoint_by_transect(). That is what "match the imagery
    dates" means: both sources then describe motion between the same two
    moments, not between two calendar years.

THE ARCHIVE DID NOT DO THIS
    5-scr/archive/coastsat_lrr_superseded_20260810/dsas_coastsat_specific_dates/
    is named for this method but was produced with SURVEY_DATES = [] in
    coastsat_domain_lrr_specific_dates.py, so it fell through to the
    continuous-range mode. Its median n_obs (417) matches the plain calendar
    fit (414.5) and its statistics are the calendar comparison's to two
    decimals. This script is the analysis that folder's name promised.

THE TWO ANCHORS, AND ONE DISCREPANCY
    1997-09-27 and 2019-09-07. The 2019 date is what nc_shorelines.geojson
    carries in SHR_DATE (Wet-Dry, 32 features). The 1997 date is Hannah's
    (2026-09-22) and matches the date commented into
    coastsat_domain_lrr_specific_dates.py, so two independent records of the
    survey agree on it.

    THE INVENTORY AGREED ONLY AFTER IT WAS FIXED. nc_shorelines.geojson
    stamped every 1997 feature 1/1/1997, a placeholder Hannah entered when the
    date was not to hand. On 2026-09-22 the two features covering this study
    area -- "Outer Banks - National Seashore" and "Outer Banks - North of
    Oregon Inlet" -- were restamped 9/27/1997; the other 21, elsewhere in the
    state and flown on other days, still carry the placeholder. See
    1-observations/shoreline_inventory/PROVENANCE.md.

    THAT EDIT DOES NOT REACH THIS SCRIPT. The dates below are module
    constants; nothing here opens the geojson, so the restamp changed no
    number in this comparison (verified by re-running: identical to three
    decimals). The two records now agree, which is worth having, but they are
    still two records.

    Both window widths are still drawn, because the two anchors are late
    September and early September -- nearly the same point in the seasonal
    cycle, so a tight window is meaningful, and agreement between the widths
    says the result does not depend on how much of the year is swept in.

    The two DSAS ends are different proxies (MHW in 1997, wet/dry in 2019)
    while CoastSat is one proxy throughout. That is a property of the DSAS
    rate, not of the matching, and no window width fixes it.

USAGE
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_datematched.py
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_datematched.py --slide
```

Notes that were in the code:

```text
The survey date, from Hannah (2026-09-22), matching the date commented into
coastsat_domain_lrr_specific_dates.py. The inventory's 1997-01-01 is her own
placeholder, not an alternative reading; see the header.
```

<details><summary>Function notes (the original docstrings)</summary>

**`coastsat_datematched()`**

```text
Per-domain endpoint rate from CoastSat, anchored on the survey dates.

Returns (per-domain Series, per-transect frame). A transect with no
position inside one of the two windows yields NaN and drops out.
```

</details>

### dsas_vs_coastsat/dsas_vs_coastsat_raw.py

The two shoreline-rate sources against each other with no smoothing: raw per-domain LRR, DSAS and CoastSat.

From the script's original header:

```text
The two shoreline-rate sources against each other with NO SMOOTHING: the raw
per-domain mean LRR from DSAS and from CoastSat, on the two DSAS windows
(Hannah, 2026-09-22: "I want to see it without smoothing").

    4-comparisons/dsas_vs_coastsat/calendar_windows/dsas_vs_coastsat_raw.png
    4-comparisons/dsas_vs_coastsat/calendar_windows/slides/ (--slide)
    4-comparisons/dsas_vs_coastsat/calendar_windows/supporting/*.csv

WHY IT IS SEPARATE FROM 6-scr-smooth/dsas_vs_coastsat/
    That folder exists to argue about the SMOOTHING -- every figure in it
    draws a LOWESS, and its own README calls these windows retired. The
    question here is different and prior to it: before any smoothing, do the
    two sources say the same thing about the same 500 m of beach? So it sits
    with the other comparisons, and the figure carries no LOWESS at all.

WHAT THE CoastSat SIDE IS
    REFIT HERE from the current time series (Hannah, 2026-09-22), not read
    from the archive: every transect in the current transect_domain_lookup.csv
    is fitted over the DSAS window with coastsat_lrr.compute_lrr --
    the same loader, the same date filter, the same OLS, the same 3-position
    minimum as the live windows -- and averaged per domain. Only the window
    differs from 3-rates/coastsat/lrr/.

    The refit is NOT written to 3-rates/coastsat/lrr/. hat_observed_rates
    .windows() enumerates that tree by scanning it, so a 1978_1997 folder
    there would become a window for every caller: rates_figures would draw
    figures for it, and the shared y bound of EVERY window figure is the
    largest domain mean over all windows, so the existing figures would
    change. These years are a comparison, not a rate product, and they stay
    inside this folder.

    The archived fits from 5-scr/archive/coastsat_lrr_superseded_20260810/
    are still read, for one purpose: the table carries them beside the refit
    so the two can be differenced. Nothing in this folder feeds a run.

    READ THE FIRST PANEL WITH CARE. CoastSat imagery begins in 1984, so its
    "1978-1997" rate is fitted from 1984-06-17 -- 13.5 years against DSAS's
    19, and it misses the six years at the start entirely. The second window
    is a fair comparison (CoastSat 1997-01-12 to 2019-12-28); the first is
    two different periods with one label, which is worth more of the
    disagreement than any method difference.

SLIDE VERSION
    --slide draws the same two panels on a 3.4 in canvas, as
    dsas_vs_coastsat_raw_slide.png. The reach is NOT cut down the way
    lrr_transect_zoom's slide is: the whole island IS the comparison here, and
    the thing a viewer reads off it -- the two lines apart in (a), together in
    (b) -- survives the smaller canvas, where a nine-domain crop would not
    show it at all. Only the labels, ticks and line weight come down.

USAGE
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_raw.py
    python scripts/input_prep/5-scr/4-comparisons/dsas_vs_coastsat/dsas_vs_coastsat_raw.py --slide
```

Notes that were in the code:

```text
Village names and structure labels are page furniture: at 3.4 in
they collide with each other and with the data (2026-09-22). The
bands stay, unlabelled, so the villages are still locatable.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_series_cache()`**

```text
{transect_id: chainage frame} for every transect in the current lookup.

Read once and fitted to both windows, rather than once per window: the
time series are the same file either way.
```

**`coastsat_refit()`**

```text
Per-domain mean LRR refitted from the current time series.

Same method as the live windows: every position inside the calendar
window, no outlier filter, no weighting, at least MIN_OBS positions.
Returns (per-domain Series, per-transect frame).
```

**`agreement()`**

```text
n, bias, RMSE and r over the domains where both sources have a value.

bias is CoastSat minus DSAS, so positive means CoastSat reports the more
seaward rate.
```

</details>

### duneline_positions/duneline_positions.py

Where the dune line sat across the island in 1997, 2009 and 2023, the lines for model years 1996, 2010 and 2024.

From the script's original header:

```text
Where did the dune line sit across Hatteras Island in 1997, 2009 and 2023 --
the lines that stand for the model years 1996, 2010 and 2024? Positions, not
change: a set of maps and alongshore profiles, designed by interview with
Hannah on 2026-09-18.

THE FIGURES   data/hatteras_init/5-scr/4-comparisons/duneline_positions/
    overview/duneline_positions_overview.png
        The island in three north-up segments side by side (south: Cape Point
        to Avon, GIS 1-30; central: GIS 31-60; north: the Tri-Village to GIS
        90), the three dune lines over the island outline and the 90 domain
        boxes, villages named, a scale bar and north arrow per segment, and a
        locator map. Orientation: at this scale the lines overlap; the zooms
        show the metres.
    zooms/zoom_<site>.png, zooms/duneline_positions_zooms.png
        Each a window 3 domains (1.5 km) alongshore by 950 m across (650 m
        landward, 300 m seaward of the 2023 line), every panel the same extent
        and scale, north-up over the 2023 NOAA orthomosaic (D:, read through
        its overviews): the three dune lines, edged for contrast, and NC-12
        white with the year by dash pattern. Named sites Buxton (GIS 1-15),
        Avon (21-31), the Tri-Village (68-83) and Mirlo Beach / the S-curves
        (84-90), each window centred on the site's domain with the largest
        |net dune change| 1997-2023; plus the domain outside them with the
        largest change. The rule and the choices are in
        supporting/zoom_sites.csv. The combined page puts every site on one
        sheet.
    context/dune_to_nc12.png
        Distance from the dune line to the NC-12 centreline per domain, one
        line per year: along each 100 m transect, the road's station minus
        the dune's (both measured from the offshore datum by
        duneline_to_raw_offsets.intersect, the function that builds the dune
        stations), positive where the road is landward of the dune. Road
        line per year: 1978 export for 1997, 2008 export for 2009 (the
        model's ROAD_LINE_FOR_YEAR), today's NC-12 for 2023
        (road_offset/raw_offset/current/).
    context/beach_width.png
        Dune line to CoastSat shoreline per domain, one line per year: along
        each CoastSat transect, the shoreline position (the mean within +/-6
        months of the dune-line image date, 3-rates/coastsat/endpoint) minus
        the distance at which the dune line crosses that transect.
    supporting/   PDFs, CAPTIONS.md, and the tables behind every figure.

YEAR COLOURS: one ordered ramp, light grey 1997 -> slate 2009 -> ink 2023, so
the order reads at a glance and red / blue stay free for seaward / landward.
Every figure is also published to output/figures/2-observations/duneline/.

USAGE
    python scripts/input_prep/5-scr/4-comparisons/duneline_positions/duneline_positions.py
    python scripts/input_prep/5-scr/4-comparisons/duneline_positions/duneline_positions.py --no-imagery
```

Notes that were in the code:

```text
CoastSat transects often START seaward of the dune line (256-470 of
906 missed it in the first draw, 2026-09-18), so extend each one
EXT_M landward along its own direction and measure from the ORIGINAL
origin: a dune landward of it gets a negative chainage.
```

```text
ONE extent for every segment (centred on each), so the three panels share
a scale and line up; the ocean side carries the village brackets
```

```text
THE ZOOMS (reworked 2026-09-18 after the first draw). At 5-15 domains a panel
was 5-7 km tall and the three lines sat on top of each other; the point of a
zoom is the tens of metres between them. Each site now shows a 3-domain
(1.5 km) WINDOW centred on the domain of that site with the largest
|net dune change| 1997-2023, every panel at the same extent and scale.
```

<details><summary>Function notes (the original docstrings)</summary>

**`road_stations()`**

```text
Per transect: the station of the SEAWARD-most road crossing (the one
nearest the dune). duneline_to_raw_offsets.intersect takes the landward-
most, right for a dune line but not for a road: at Buxton (GIS 8-9) the
transects also cross the leg of NC-12 that turns west toward Frisco, and
the landward-most rule put the road 1.6-2 km inland there (first draw,
2026-09-18). At Rodanthe the 2022 Jug Handle bridge is the only crossing,
so it is kept.
```

**`beach_width()`**

```text
Per CoastSat transect: shoreline chainage minus the chainage where the
dune line crosses the transect (m), per year.
```

**`pick_reach()`**

```text
The PICK_WIDTH-domain reach outside the named sites with the largest
mean |net dune change| 1997-2023.
```

**`_panel_title()`**

```text
Letter and title left-aligned on one line, so a narrow map panel
cannot overlap them (the house _title centres the text).
```

**`_zoom_extent()`**

```text
A fixed window, 3 domains alongshore by 950 m cross-shore, set on the
2023 dune line: 650 m landward (NC-12 is 270 m behind the dune on
average, so it is usually in view) and 300 m seaward. The first rework
widened each window to take every road in, which at Buxton (the Frisco
turn) and Rodanthe (the Jug Handle) stretched every panel and squeezed
the lines back together.
```

</details>

### mean_shoreline_windows/coastsat_mean_shoreline_windows.py

One period's CoastSat mean shoreline over two windows, differenced transect by transect and domain by domain.

From the script's original header:

```text
coastsat_mean_shoreline_windows.py -- one period's mean shoreline, two windows
The CoastSat mean shoreline a period starts from, averaged over two windows,
differenced transect by transect and domain by domain:

    3-yr   the calendar span the shoreline offset v1 was built from
           (1995_1997 for 1996, 2009_2011 for 2010)
    2-yr   +/-1 yr of the start DEM's lidar flights, the window the offset v2
           is built from (1995-10-12_1997-10-12, 2008-08-17_2010-08-17)

Asked for by Hannah on 2026-09-29 ("a new comparison of the 3 year window vs
the 2"), after the interview that moved the offset onto the DEM-centred
window ([[cascade-shoreline-offset-dem-centred-window]]).

WHAT IS COMPARED
    Both windows' per-transect means as coastsat_mean_shoreline.py stored them
    (mean_shoreline/<label>/transect_means_<label>.csv). Nothing is re-averaged
    here. The difference is 2-yr MINUS 3-yr along the transect, so + is
    SEAWARD: the DEM-centred line sits seaward of the calendar one. A domain's
    value is the mean over its transects included in BOTH windows.

WHAT IT CANNOT SAY
    The two windows share most of their satellite passes (1996 shares
    1995-10-12..1997-10-12, all of the 2-yr window), so their means are not
    independent and no significance is attached to a difference. The standard
    errors are drawn beside the difference so a reader can see how big a
    difference is against the sampling noise of either mean.

    The island offset is zeroed on its own minimum, so a UNIFORM shift between
    the windows does not reach the model; only the alongshore-varying part
    does. Both are reported.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/mean_shoreline_windows/<period>/
    mean_shoreline_windows_<period>.png        (a) the difference, (b) the SE
    supporting/mean_shoreline_windows_<period>.pdf, CAPTIONS.md
    supporting/domain_comparison_<period>.csv  per domain: both means' SE and
                                               positions, the difference
    supporting/transect_comparison_<period>.csv per transect, both windows
    PROVENANCE.md
    Resolved through hat_observed_rates.MEAN_SHORELINE_WINDOWS.

USAGE
    python coastsat_mean_shoreline_windows.py                 # both periods
    python coastsat_mean_shoreline_windows.py --periods 1996
```

Notes that were in the code:

```text
the two lines themselves, on a photograph (Hannah, 2026-09-29: "plot the
shorelines themselves to visualize the difference")
```

```text
The photograph nearest the start DEM's survey: the USGS frame of 1996-10-14
(two days before the ALACE flights ended); for 2010 no georeferenced photo
lies in the window, and the NOAA NGS mosaic of March 2008 is the nearest,
as for the 2009_2011 imagery figures.
```

```text
The middle half of each domain: at 500 m alongshore the oblique shore makes
the panel wide and a few metres between the lines is under a line width.
```

```text
On the photographs the two windows are the house RdBu pair, 3-yr red and
2-yr blue (Hannah, 2026-09-29: "more distinct instead of grey and purple").
Grey and purple were the BASE/ACCENT pair of the other figures here, but on
a photograph grey read as a white line on sand and purple sat too close to
it; red against blue separates on hue and on greyscale. Only these figures
carry no red/blue difference scale, so the pair means nothing else on them.
```

```text
The whole island in six north-up segments of 15 domains (Hannah, 2026-09-29,
after the duneline_positions overview). At 7.5 km a panel, 10 m is ~0.2 pt:
the two lines would print as one. So the geometry drawn is the true 3-yr
line, and the DIFFERENCE is carried by colour, transect by transect.
(Dropped and restored the same day: Hannah kept it beside the datum figure.)
```

```text
ONE extent for every panel, centred on each segment, so all six share
a scale; the ocean side carries the domain numbers
```

```text
a grey edge under the colour, so where the windows agree (white)
the shoreline does not vanish against the white ocean
```

```text
the two lines as the offset build sees them (Hannah, 2026-09-29: "I was
thinking of this figure" -- 2-brie-offset/<year>/comparisons/duneline_vs_
shoreline, compare_offset_sources.py). The same 2 x 3 grid of vertical
strips, the same frame: the distance from the shared offshore datum along
the 100 m model transects, per 500 m domain. That is the quantity the
island offset is made of, so this is the offset v1 vs v2 before the build.
```

```text
SECTIONS and GRID_COLS as compare_offset_sources.py sets them (its 2 x 3
note: rows, not six across, keep width on the axis the gap is measured on).
```

```text
stations grow LANDWARD, so 3-yr minus 2-yr is + where the 2-yr line is
seaward -- the same sign as everywhere else in this folder
```

```text
Each section is a PAIR: the two profiles, and beside them a narrow strip
of the difference itself on the same domain axis (Hannah, 2026-09-29).
The profiles span 1-2 km across and the difference is at most ~12 m, so
in the profile panel the two lines print as one; the strip is where the
difference can be read. One symmetric scale for every strip.
```

<details><summary>Function notes (the original docstrings)</summary>

**`line_points()`**

```text
The window's mean line in alongshore order, as coastsat_mean_shoreline
strung it.
```

**`lines_figure()`**

```text
One panel per domain: both lines on the photograph, north up, at one
scale. At the island scale a 10 m difference is a hair; one domain
(500 m) by ~200 m across makes it visible.
```

**`datum_stations()`**

```text
Per domain, each window's mean station from the offshore datum.

Computed in memory with the offset build's own intersection
(duneline_to_raw_offsets.intersect on its 100 m transects), so the 2-yr
line needs no raw file in 2-brie-offset before step 3 writes one, and both
windows go through the identical code. One value per transect first, then
the domain mean, as compare_offset_sources._raw_domain_means does.
```

**`_town_bands_y()`**

```text
Village spans against a VERTICAL alongshore axis, as
compare_offset_sources._town_bands_alongshore_y draws them.
```

</details>

### shoreline_vs_duneline/coastsat_vs_duneline.py

Does the digitized dune line move with the CoastSat shoreline? Net change on both sides, over the same interval.

From the script's original header:

```text
Does the digitized dune line move with the CoastSat shoreline?

WHAT IS COMPARED  (NET CHANGE ON BOTH SIDES since 2026-09-18; Hannah: "I
wanted the coastsat vs duneline comparison to both be using net position
change")
    Dune side      3-rates/duneline/endpoint/<window>/: the end dune line minus
                   the start line per 100 m transect, domain mean, seaward
                   positive.
    Shoreline side 3-rates/coastsat/endpoint/<window>/: the mean CoastSat
                   position within +/-6 months of each dune-line survey date,
                   end minus start, domain mean. The like-for-like quantity
                   for two surveys.
    Both read from the stored products, not computed here. Shown as the net
    change over the survey interval (m/yr) so the four windows share one axis;
    the metres are in domain_comparison.csv.

    Until 09-18 the shoreline was ALSO drawn as the CoastSat LRR, an OLS
    through ~250 dates, and the correlations reported against both. That is
    not a two-survey quantity; it stays the model's scoring target in
    3-rates/coastsat/lrr/ and in model_vs_observed/vs_shoreline/. The helpers below
    (window_mean, endpoint_by_transect, KNOWN_SURVEY_DATES) are kept: the
    stored CoastSat endpoint product is built with them.

SURVEY DATES
    The dune-line files carry no date. KNOWN_SURVEY_DATES holds them by line
    VINTAGE: 1984-09-19 and 1997-10-12 from the Henderson USGS metadata on
    D:, 2004-05-25 and 2009-05-30 from the Google Earth captures (Hannah,
    2026-09-15), and None for the 2023 NOAA set until its flight date is
    known. A None is centred mid-year of the line's year, PROVENANCE.md says
    so, and a sensitivity block reports how far the endpoint rate moves for
    a +/- 6 month shift of that centre. A period year reaches its vintage
    through hat_topo_version.DUNE_LINE_FOR_YEAR (2010 reads the 2009 line,
    2024 the 2023 one).

METHOD
    Every raw dune file is built by duneline_to_raw_offsets.py since
    2026-09-15 (1984 and 2004 were ArcGIS exports until that afternoon, one
    metre landward of the exact crossing; see raw_offsets/PROVENANCE.md), so
    a change between any two years carries no method term.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/endpoint_net_change/<start>_<end>/
         (was 4-comparisons/coastsat_vs_duneline/ until 2026-09-19; the
         alongshore figures are in METRES with the beach-width gap since then)
             scatter_dune_vs_coastsat.png     shoreline vs dune, 1:1
             alongshore_dune_vs_coastsat.png  the two net changes by domain
             supporting/                      the PDFs, CAPTIONS.md,
                 domain_comparison.csv            one row per GIS domain (m and m/yr)
                 PROVENANCE.md
         data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/endpoint_net_change/
             alongshore_four_windows.png      every window on one y axis,
                                              stacked full width (--grid;
                                              --layout grid for the 2 x 2)

USAGE
    python coastsat_vs_duneline.py --start-year 1984 --end-year 2004
        # the dates come from the stored products
    python coastsat_vs_duneline.py --grid                 # every window, stacked
    python coastsat_vs_duneline.py --grid --layout grid   # the 2 x 2 by period
```

Notes that were in the code:

```text
Colours (Hannah, 2026-09-15, third pass): the two CoastSat estimators are
the SAME feature measured two ways, so they share the blue family, one line
each -- the LRR dark (the RdBu blue pole the house uses for the sea side),
the endpoint a lighter blue and dashed so the pair still separates in
greyscale. The dune line is the house RdBu red against them, the pair the
other figures already use, so the red/blue contrast reads at a glance
(Hannah, 2026-09-15). Here the pair means FEATURE, dune vs shoreline, not
vintage and not sign; the caption says so. Grey, purple, sand brown, orange
and REF green were tried and rejected, the green as too dark to see.
```

```text
Since 2026-09-18 the shoreline is ONE line, its net change at the dune
dates, drawn in the dark blue the LRR had.
```

```text
Keyed by LINE VINTAGE (the year in the geojson name), not by period year: a
period finds its vintage through hat_topo_version.DUNE_LINE_FOR_YEAR.
1984, 1997: the Henderson USGS metadata on D: (Calendar_Date). 2004, 2009:
Google Earth capture dates (Hannah, 2026-09-15); the raw_GE frames carry
none. 2023: NOAA NGS imagery under D:\Hatteras_GIS\Aerial\2023 whose
metadata gives only the 2015-2023 series extent, so None until Hannah
supplies the flight date; a None is centred mid-year and flagged.
```

```text
The window is IN the stem (Hannah, 2026-09-21): five windows wrote this
same basename, so a figure lifted out of its folder could not be told
from the other four.
```

```text
The alongshore figures are in METRES since 2026-09-19 (Hannah: one form for
every shoreline-vs-dune figure, net change with the beach-width gap).
```

```text
The gap between the two lines: widened solid grey, narrowed hatched. Shared
by net_change_vs_duneline.py and total_change_vs_duneline.py.
```

```text
the offshore shoals, and the fills placed inside the window, as the
3-rates figures mark them (2026-09-19)
```

```text
Both sides are OBSERVED here -- two snapshots differenced, no rate
anywhere -- which is what separates this folder from total_change/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`dune_position_by_domain()`**

```text
Mean ORIG_LEN per GIS domain, first row per transect, as the hindcast
loader (`load_absolute_dune_distance`) reads it. Grows LANDWARD.
```

**`shade_beach_width()`**

```text
The space between the shoreline and dune-line changes: solid grey where
the beach WIDENED (shoreline change > dune-line change), hatched where it
narrowed.
```

**`draw_alongshore()`**

```text
One alongshore panel in the house form: village bands, groin and piers,
open frame, y grid, symmetric limits. TWO lines since 2026-09-18, both net
change over the same survey interval: the CoastSat shoreline blue, the
dune line red; in METRES with the beach-width gap shaded since 2026-09-19.
Call it after the legend is placed (structures() tests its labels against
the layout).
```

**`_chains()`**

```text
Windows linked end-to-start: [(1984,2004),(2004,2024)],
[(1996,2010),(2010,2024)] -- the same rule as coastsat_lrr_windows.py.
```

**`four_windows_figure()`**

```text
Every model window on ONE y axis, one full-width panel per window in
chain order (`layout="column"`), or the 2 x 2 by period (`"grid"`). Reads
each window's supporting/domain_comparison.csv; run the windows first.
The context window 1996_2024 is left out: it is not a model window.
```

</details>

### shoreline_vs_duneline/net_change_vs_duneline.py

Net shoreline change (CoastSat) against net dune-line change per GIS domain, 1996-2024 and its halves, in metres.

From the script's original header:

```text
Net shoreline change (CoastSat) against net dune-line change, per GIS domain,
over 1996-2024 and its two halves, in METRES. Built 2026-09-18 (Hannah, by
interview: "the total long term shoreline position change and then compare
this to the total change in duneline position").

WHAT IS COMPARED
    Both sides are read from the stored endpoint products, never recomputed:
        shoreline   3-rates/coastsat/endpoint/<window>/   mean CoastSat position
                    within +/-6 months of each dune-line survey date, end minus
                    start
        dune line   3-rates/duneline/endpoint/<window>/   end line minus start
                    line
    Same windows, same survey dates (1997-10-12, 2009-05-30, 2023-07-01
    ASSUMED), seaward positive on both, so the gap is
        beach-width change = shoreline change - dune-line change
    positive where the beach widened (the waterline gained on the dune).
    Both products share the 2009 date, so for each of them the two halves add
    up to the whole exactly; the script checks it.

WHAT IS DRAWN
    Three stacked panels, one per window (1997-2023, 1997-2009, 2009-2023),
    one y axis in metres: the shoreline blue, the dune line red (the house
    pair for FEATURE in coastsat_vs_duneline), the gap between them shaded
    grey. The village spans are a strip along the top of each panel, not the
    usual full-height wash, because the grey gap already shades. Groin and
    piers as hairlines, the offshore shoals as faint hatched amber boxes, the
    model-input fills as bars above the top panel -- the same marks as the
    two halves figures.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/endpoint_net_change/chains/
    (was 4-comparisons/net_change_1996_2024/ until 2026-09-19)
    net_change_chain_1996_2010_2024.png   also published to
                                          output/figures/2-observations/shoreline_vs_duneline/
    supporting/
        domain_comparison.csv   one row per window x domain: shoreline, dune,
                                beach-width change, whether they agree in sign
        island_summary.csv      per window: means, r, slope, RMSE, sign
                                agreement
        net_change_shoreline_vs_dune.pdf, CAPTIONS.md
    The per-transect CoastSat change (with positions per end window) is in
    3-rates/coastsat/endpoint/<window>/transect_endpoint.csv.

USAGE
    python scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/net_change_vs_duneline.py
```

Notes that were in the code:

```text
The chain figure of shoreline_vs_duneline/endpoint_net_change/ since 2026-09-19
(was 4-comparisons/net_change_1996_2024/net_change_shoreline_vs_dune).
```

### shoreline_vs_duneline/smoothed_lowess7_vs_duneline.py

The halves-overlay sheets with both curves through the model target's 7-domain LOWESS.

From the script's original header:

```text
The two halves-overlay sheets again, with BOTH curves passed through the
model target's alongshore LOWESS at a 7-domain (3.5 km) window. Built
2026-09-22 (Hannah, by interview) as a place to see what the smoother does to
the shoreline-vs-dune-line comparison, in the two readings of the shoreline
side side by side.

WHAT IS DRAWN, one sheet per shoreline reading, 1996-2010 above 2010-2024:

    projected     shoreline = the 1996-2024 LRR x 14 yr, the SAME in both
                  panels (the long-term trend carried onto each half)
    total_change  shoreline = each half's OWN LRR x its own 14 yr

    The dune side is the same either way and always follows the sub-period:
    the measured net change between the two digitized lines bounding that
    half. Only the shoreline reading differs between the two sheets.

THE SMOOTHING (Hannah's choices, 2026-09-22)

    window        7 domains = 3.5 km.
    both sides    BOTH curves get the same pass at the same window, at
                  TRANSECT resolution, then average to domains. Smoothing one
                  side and not the other would make the gap between them an
                  artefact of the treatment rather than a beach-width change
                  -- the same trap the 3-rates smoothed panels avoid.
                  CoastSat has ~10 transects per domain and the dune line
                  exactly 5; `rates_figures._along` gives both an even spread
                  inside their domain, so the two are handled identically.
    GIS 1-10      kept at their RAW domain means, the Oregon Inlet boundary
                  treatment (`coastsat_lowess.LowessConfig.skip_southern_domains`).
                  Hannah chose to keep it so the figure shows the target the
                  way the model actually sees it. The cost, stated here so it
                  is not read as a result: those ten domains are IDENTICAL to
                  the unsmoothed sheet by construction, and any difference
                  there is not the smoother.
    raw kept      the raw domain means stay on the figure as faint dots
                  behind each curve, so what the smoother removed is visible
                  without leaving the sheet.

WHAT THIS IS FOR.  It is a test folder, not a product: the question is how
much of the shoreline-vs-dune-line disagreement survives smoothing at the
scale the model resolves. The unsmoothed sheets it is paired with are
`coastsat_{projected,total_change}_vs_duneline_endpoint/all_windows_stacked/
*_halves_overlay.png`.

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_lowess7/
    lowess7_projected_vs_duneline_1996_2010_2024_halves_overlay.png
    lowess7_total_change_vs_duneline_1996_2010_2024_halves_overlay.png
    domain_smoothed.csv     per domain per window per product: both sides raw
                            and smoothed, and the beach-width gap of each
    PROVENANCE.md           what changed, per window and per product
    README.md, supporting/  (PDFs, CAPTIONS.md)

USAGE
    python scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/smoothed_lowess7_vs_duneline.py
    python ... --window 7          # LOWESS width in domain units
```

<details><summary>Function notes (the original docstrings)</summary>

**`smooth_side()`**

```text
One alongshore LOWESS pass at TRANSECT resolution, averaged to domains.

The same two steps the scoring target is built through, so the curve here
is the quantity the model is graded against rather than a different
smoother that happens to look similar.
```

</details>

### shoreline_vs_duneline/total_change_vs_duneline.py

Total shoreline change against the dune line's measured net change per GIS domain, over 1996-2024 and its halves.

From the script's original header:

```text
TOTAL shoreline change against the dune line's MEASURED net change, per GIS
domain, over 1996-2010, 2010-2024 and the whole 1996-2024. Built 2026-09-21
(Hannah, by interview: "the difference between the measured dune line change
between 1996-2010 and 2010-2024, and compare that with the LRR net position
change over those same periods"); renamed from lrr_net_change.py and merged
with projected_vs_duneline.py the same day, after a second interview settled
the vocabulary.

WHY "TOTAL" AND NOT "PROJECTED".  A rate turned into a distance is named by
the window it was FITTED on, never by the arithmetic:

    TOTAL      the rate is evaluated over the SAME window it was fitted on.
               Every window here is total: LRR(1996-2010) x 14 yr,
               LRR(2010-2024) x 14 yr, LRR(1996-2024) x 28 yr. Nothing is
               extrapolated.
    PROJECTED  the rate is carried onto a window it was NOT fitted on. None
               of that happens here; it lives in
               3-rates/coastsat/projected/ and, against the model, in
               output/comparisons/target_comparison/projected/.

THE MERGE (2026-09-21, Hannah, by interview).  This folder absorbed
`shoreline_vs_duneline/projected/1996_2024`, which was the 1996-2024 rate x
the 25.72 yr DUNE-LINE INTERVAL rather than the 28 calendar years. That is
the same fit window with a shorter span, so it was never a projection either.
Its three numbers were already carried here as the `*_dune_interval_m`
columns -- checked identical to 0 m before the merge -- so the headline
figure is the 28 yr CALENDAR span and the dune interval is a column and a
caption line. Hannah chose the calendar span so every window in the tree is
read the same way. The old folder is in superseded_20260921/.

    Hannah, 2026-09-21: "Dont fit on the exact dune dates, the year is most
    important" -- so the fit window stays the stored CALENDAR window,
    1 January to 31 December.

WHAT IS COMPARED, per 500 m GIS domain, SEAWARD POSITIVE, in metres
    shoreline   3-rates/coastsat/lrr/<window>/transect_lrr_full.csv,
                lrr_m_yr x (end_year - start_year), mean over the ~10 CoastSat
                transects of the domain. 14 yr in each half, 28 yr over the
                whole -- the CALENDAR interval (Hannah's choice, 2026-09-21),
                matching the year labels and 3-rates/coastsat/total_change.
    dune line   3-rates/duneline/endpoint/<window>/ as stored, end line minus
                start line, mean over the ~5 dune transects. Read, never
                recomputed.
    beach width shoreline change - dune-line change; positive = the beach
                widened (the waterline gained on the dune).

THE INTERVAL MISMATCH, reported and not corrected.  The dune line is two
photographs, 11.63 yr apart in the first half (1997-10-12 to 2009-05-30),
14.09 yr in the second (to 2023-07-01, assumed) and 25.72 yr over the whole,
while the shoreline is carried over the full calendar span in each. So in the
first half the gap holds about 2.4 yr of shoreline drift on top of beach-width
change. Every table carries the interval-matched alternative beside the
headline value (`*_dune_interval_m`, the same rate x the dune line's own
span) and every PROVENANCE.md reports both, so the size of the term is
visible rather than argued about.

FIGURE TITLES carry quantity, window and method (Hannah, 2026-09-21), e.g.
"Total shoreline change vs dune line, 1996-2010 (LRR x 14 yr)".

WHAT IS DRAWN
    <window>/    two_panel, shaded_gap, overlay -- the three presentations,
                 one folder per window.
    chains/      1996-2024 stacked over its two halves on one y axis, the
                 shape of endpoint_net_change/chains.
    difference/  second half minus first half on both sides: where the
                 shoreline trend sped up or slowed, did the dune line follow?

OUTPUT   data/hatteras_init/5-scr/4-comparisons/shoreline_vs_duneline/total_change/
    <window>/domain_comparison.csv, the three PNGs, PROVENANCE.md,
             supporting/ (PDFs, CAPTIONS.md)
    chains/total_change_chain_1996_2010_2024.png,
             supporting/{domain_comparison,island_summary}.csv
    difference/total_change_difference.png, supporting/domain_difference.csv

USAGE
    python scripts/input_prep/5-scr/4-comparisons/shoreline_vs_duneline/total_change_vs_duneline.py
    python ... --start-year 1996 --end-year 2010
```

Notes that were in the code:

```text
How the SHORELINE side is read. The dune side never varies -- it is the
measured endpoint between the two digitized lines bounding the window -- so
this is the only choice the script makes, and it is what the output folder
is named for. Mirrors 3-rates/coastsat/{total_change,projected}.
```

```text
no 1996_2024: there the rate window IS the change window, so the
answer is the TOTAL product, not a second copy under another name.
```

```text
BOTH halves carry the same 1996-2024 rate over the same 14 yr, so
the shoreline side is IDENTICAL in them. That makes
change_between_periods zero by construction on that side -- not
drawn -- but it is exactly what makes the STACK worth drawing: one
prediction against two different dune-line outcomes, which is the
question this product exists to ask (Hannah, 2026-09-22).
```

```text
The shared fixed metre axis (Hannah, 2026-09-22); see
coastsat_total_change.Y_HALF_M for why it is fixed and not a floor.
```

```text
`w` indexes the DUNE side and the change window; the shoreline rate may
come from a different window (see PRODUCTS).
```

```text
The axis is fixed at +/-100 m, so a domain-mean line can leave it. Mark
it at the edge and hand the values back for the caption; nothing that
walks off the panel should do so silently.
```

```text
Said once above the panels, because it is what the two panels have in
common -- and for `projected` it is also the whole point of the figure.
```

```text
The fit window is IN the label, not implied by the folder (Hannah,
2026-09-21): "CoastSat LRR 1996–2010 × 14 yr" says where the rate came
from, and the reader can see it matches the window in the title -- which
is the whole difference between total change and a projection.
```

```text
Quantity, window, method (Hannah, 2026-09-21). "Total" because the rate
is fitted on the window it is evaluated over; see the module docstring.
```

```text
NOT shade_beach_width: both lines are differences, so the gap between
them is disagreement, not a width. One neutral fill, no hatch, so the
figure does not borrow a vocabulary that would read as beach width.
```

```text
Each half is fitted on itself, so there is no single fit window to name
in the title; the legend carries it per line.
```

```text
The two-panel overlay of the halves ALONGSIDE the three-panel chain, not
instead of it (Hannah, 2026-09-22). The chain answers "how does the whole
period relate to its halves"; this one is the direct counterpart of the
projected product's sheet, so the two can be laid side by side to see what
changes when the rate is fitted per half instead of over the full record.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load()`**

```text
Per-domain change from both sides, plus the per-transect shoreline.

The shoreline is the stored calendar-window LRR carried over the CALENDAR
interval; `*_dune_interval_m` is the same rate over the dune line's own
span, kept beside it as the diagnostic for the mismatch.
```

**`_header()`**

```text
The two sides named on the canvas (Hannah, 2026-09-22).

The interval mismatch is the reason this is worth the space: the shoreline
is carried over the full CALENDAR span and the dune line spans whatever
its two photographs do, so the gap between them holds that difference as
well as beach-width change. Stated here, it cannot be missed.
```

**`_pad_title()`**

```text
Lift the centred title clear of the fill bars.

draw_fills puts its bars at 1.025 in axes fractions with the year label
above them, so a title at the default pad lands on top of "2022 fill".
_title() has already set the bold letter at the left; re-setting only the
centred string keeps it and moves both (pad is per-axes in matplotlib).
```

**`halves_overlay_figure()`**

```text
The two halves' overlay panels on one sheet, 1996-2010 above 2010-2024.

TITLE SPACE (Hannah, 2026-09-22: "be strategic ... concise yet organized").
The method is the SAME in both panels -- that is the whole point of the
projected product -- so it is stated ONCE in the header and never repeated.
Each panel title then carries only what actually differs between them: the
window, and the two dune-line dates with their real interval. Nothing is
said twice, and nothing a reader needs is only in the caption.
```

</details>
