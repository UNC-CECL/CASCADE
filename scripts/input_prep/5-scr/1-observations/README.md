# 1-observations — what was measured, not what we fitted

The shoreline record as it arrives: digitized lines, the NC Coastal Management
historical set, and the CoastSat satellite chainage. Nothing here fits a rate.
That is `../3-rates/`, and it cannot run until `../2-transect-frame/` has tied
these observations to domains.

```
shoreline_inventory/
    shoreline_inventory.py
        Cross-source inventory for the whole study area: what exists from each
        source, how the sources overlap in time, and where the gaps are. The
        study area comes from a spatial filter file, so the same script can be
        pointed at a narrower reach (the Buxton groin transects, say) without
        being edited.
        Writes 1-observations/shoreline_inventory/shoreline_position_output/.

mean_shoreline/
    coastsat_mean_shoreline.py
        One averaging window's MEAN satellite shoreline, as a line on the
        ground: each CoastSat transect's chainage averaged over the window,
        geolocated, and the mean points strung into a single polyline that
        2-brie-offset intersects with the 100 m transect frame exactly as it
        intersects a digitized dune line.

        It is the only script in 5-scr that needs a POSITION rather than a
        difference, and that is the whole job. Every other CoastSat product
        differences chainage, so each transect's arbitrary origin cancels;
        here it does not. Aggregated to the 90 domains, raw chainage spans
        124 m alongshore while the geolocated position spans 6222 m -- the
        origins follow the shore around the cape.
        Writes 1-observations/mean_shoreline/<start>_<end>/.
    coastsat_mean_shoreline_storm_check.py
        Was a window mean shaped by a storm? The storms around each window,
        ranked against 1984-2024, and the mean without post-storm passes.
        Writes <start>_<end>/storm_check/.
    coastsat_mean_shoreline_compared.py
        The three mean shorelines the model starts from and is graded against
        (1996, 2009, 2025), drawn as the island: distance from the offshore
        datum per domain, six sections in one row.
        Writes 1-observations/mean_shoreline/compared/ and a copy to
        output/figures/2-observations/mean_shoreline/.

shoreline_patterns/
    shoreline_trajectory_classification.py
        Is a domain eroding steadily, stable, or reversing? Classifies each
        domain's trajectory from the CoastSat time series over 1984-2004,
        2004-2024 and the full record.
    shoreline_trajectory_map.py
        The same classification, and LRR magnitude, drawn on the island.
        USE_SATELLITE = True wants an Esri basemap and therefore the internet;
        False gives a plain ocean background and works offline.
```

Both trajectory scripts write under `1-observations/shoreline_patterns/`. Their
earlier outputs were deleted as stale on 2026-09-18 — the folder being empty
does not mean the scripts are broken, only that nobody has rerun them since.

Data paths are resolved through `scripts/site_layer/hat_observed_rates.py`;
none is typed here.

## detrended_position

Why do the rate products disagree between windows? Not because the fits are
noisy: there is one signal the whole island shares, on top of every transect's
own trend.

```
coastsat_detrended_position.py
    Detrends every CoastSat transect against its OWN 1996-2024 fit, reduces it
    to an annual median and averages all 906, so anything local cancels. What
    survives is a landward sag through 2005-2020 and a +16.8 m step in the
    single year to 2021. That excursion alone biases the fitted rate by
    -0.37 m/yr over 1996-2010 and +0.88 m/yr over 2010-2024, which is the
    answer to why no short window recovers the long-term rate. Writes the
    detrended matrix the next script reads, so the detrending happens once.

coastsat_position_attribution.py
    What the step IS, in five tests, each stored with its verdict including
    the two that failed: per-transect noise (rejected), nourishment (rejected
    as the cause, confirmed as a visible signal), seasonal sampling (rejected),
    the independent dune line (corroborates), spatial structure (corroborates).
    The step is real, so the grading target's sensitivity to 2021-2024 is a
    question about which period the model should represent, not about which
    data to trust. Reads stored endpoint tables for the dune-line test and
    refits nothing.
```

Built 2026-09-23, out of `3-rates/coastsat/window_convergence_1996_2024/`.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### detrended_position/coastsat_detrended_position.py

Is there one signal the whole island shares, on top of each CoastSat transect's own trend?

From the script's original header:

```text
IS THERE ONE SIGNAL THE WHOLE ISLAND SHARES, ON TOP OF EACH TRANSECT'S TREND?

Detrend every CoastSat transect against its OWN 1996-2024 fit, reduce each to
one value a year (the annual median position), and average across all 906. Any
signal that survives that averaging is common to the island, because anything
local is incoherent between transects and cancels.

WHY IT WAS BUILT (2026-09-23). `3-rates/coastsat/window_convergence_1996_2024/` found
that no window shorter than about 25 years recovers the long-term rate, and
that the answer barely varies between transects -- a transect with a fast clean
trend needs as long as a slow noisy one. A per-transect explanation cannot
produce a per-transect-invariant answer, so the cause had to be shared. This is
the search for it.

WHAT IT FINDS. One coherent excursion: roughly flat 1996-2004, a landward sag
through 2005-2020 bottoming at -7.4 m, then +16.8 m in a single year into 2021,
held through 2024. It is only 17% of the mean transect variance, but it is the
COHERENT part, and coherent is what moves an OLS slope systematically:

    window        bias the departure alone puts on the fitted rate
    1996-2010     -0.367 m/yr      (the graded window, ~37% of a median rate)
    2010-2024     +0.875 m/yr      (the chain's second leg, opposite sign)
    1996-2024      0.000 m/yr      (zero by construction -- it is detrended
                                    against this window)

That is the answer to "why do the windows disagree": not noise, a shared
multi-decadal excursion that an OLS slope cannot separate from the trend until
the window spans the whole of it. `coastsat_position_attribution.py` beside this file
then asks what the excursion IS.

WHAT THIS IS NOT. Not a rate product and not a model input. Nothing is graded
against it; it exists to explain a property of the rate products, which is why
it lives under 1-observations rather than 3-rates.

Inputs
    transect_domain_lookup.csv      2-transect-frame/, via hat_observed_rates
    CoastSat time-series CSVs       1-observations/coastsat_timeseries/

Outputs  (hat_observed_rates.DETRENDED_POSITION)
    annual_medians_detrended.csv    year x transect, m. The matrix everything
                                    else here reads, so the detrending is done
                                    once and cannot drift between scripts
    detrended_position_by_year.csv        a row per year: the index, the counts, and
                                    the nourished / untouched split
    detrended_position_by_domain.csv    a row per domain: its share of the step,
                                    and how well it tracks the index
    detrended_position.png              the index, the alongshore step, the record

Usage
    python .../coastsat_detrended_position.py
```

Notes that were in the code:

```text
The step the index turns out to contain, as two periods to difference. Named
here rather than buried so the number quoted in the README and the figure is
the same one.
```

```text
NOURISHMENT, as the hindcast receives it (hatteras_site_config
.HATTERAS_NOURISHMENT_PROJECTS, drawn in 4-mgmt-forcing/nourishment/). Held
here so the index can be split by it; `coastsat_position_attribution.py` does the test.
```

```text
One label per YEAR, not per project: Buxton and Avon are both 2022 and
two rotated labels on the same bar overprint into a smear.
```

```text
.loc throughout: a bare [2021:2024] on an integer index is POSITIONAL
in pandas and silently returns nothing, which printed the step as NaN.
```

<details><summary>Function notes (the original docstrings)</summary>

**`build_matrix()`**

```text
Detrended annual median position, year x transect, in metres.

Each transect is detrended against ITS OWN 1996-2024 OLS, so what is left
is departure from that transect's long-term behaviour and nothing else. A
transect eroding at 3 m/yr and one accreting at 1 m/yr both come out
centred on zero, which is what lets them be averaged.
```

</details>

### detrended_position/coastsat_position_attribution.py

What is the 2021 step? Five tests of the island-wide shift, including the two that failed.

From the script's original header:

```text
WHAT IS THE 2021 STEP? FIVE TESTS, INCLUDING THE TWO THAT FAILED.

`coastsat_detrended_position.py` found one signal the whole island shares, and a +16.8 m
seaward step in the single year to 2021 that holds through 2024. This script
asks what it is. Each test is here with its verdict, including the hypotheses
that were wrong, because a rejected explanation is evidence and the next person
to look at this will otherwise re-run them (Hannah's standing preference:
diagnostics and honest reporting over a tidy story).

WHY IT MATTERS MORE THAN IT LOOKS. The 2021-2024 years sit at the
highest-leverage end of the 1996-2024 fit, and that fit is the model's grading
target. Refitting on 1996-2020 instead shifts the fitted rate a near-uniform
+0.46 m/yr -- enough to flip the island median from -0.349 m/yr (eroding) to
+0.172 m/yr (accreting) and to flip the sign at 18% of transects. So whether
the step is the shoreline or the satellites decides whether the target is
credible.

THE TESTS

  1. PER-TRANSECT NOISE -- REJECTED. If each transect's own wobble set how long
     its window takes to settle, convergence time would track that transect's
     noise-to-trend ratio. Correlation is 0.01 (forward) and 0.16 (backward);
     the calmest quartile needs 27 years and the swingiest 26. A per-transect
     cause cannot give a per-transect-invariant answer, which is what sent the
     search to a shared signal in the first place.

  2. NOURISHMENT -- REJECTED as the cause, CONFIRMED as a real signal. The
     fills are Rodanthe 2014 (GIS 84-89), Buxton 2022 (6-15) and Avon 2022
     (21-28), 24 of 90 domains. The 2020->2021 jump is +16.5 m in the 672
     NEVER-NOURISHED transects against +17.4 m in the nourished ones, and it
     lands a year BEFORE the 2022 fills. The fills are nonetheless plainly
     visible against the rest of the island in the year after placement, in
     the order their fill densities predict.

  3. SAMPLING -- NOT THE CAUSE, but it sharpens the step. Observations per
     transect per year go 12.8 (2020) to 30.6 (2021), a 2.4x jump in the same
     year. It is NOT a seasonal-mix effect: comparing like quarters the step
     is present in all four. 2020 is both the most landward year and the most
     thinly sampled, so the sparse pre-2021 half is the LESS reliable half,
     not the more -- the denser record is the better one. Dropping 2020
     entirely still leaves +14.6 m between 2019 and 2021.

  4. THE DUNE LINE -- CORROBORATES. The digitized lines are 1997, 2009 and
     2023, so 2023 sits inside the stepped period and the other two before it.
     Over 1997->2009 the dune line retreats a median -11.7 m; over 2009->2023
     it ADVANCES +3.7 m, and CoastSat over the same pairs gives -2.4 m then
     +3.8 m -- the later interval matching to 0.1 m, correlated alongshore at
     r = 0.71. An artefact would have made the stepped interval carry about
     +13.6 m MORE CoastSat-minus-dune; it carries 9.7 m LESS.
     Limit: three snapshots cannot date the advance within 2009-2023, so this
     confirms direction and magnitude, not the 2021 timing.

  5. SPATIAL STRUCTURE -- CORROBORATES. A sensor or waterline bias applies
     nearly the same offset everywhere. This does not: domain means run -11.9
     to +71.5 m (sd 11.5 m between domains) while transects INSIDE a domain
     agree to 2.8 m, and eight domains step landward. The largest, GIS 1 at
     +71.5 m, is the Cape Point shoal attachment already documented in
     3-rates/coastsat/5yr_bins/README.md.

VERDICT. The step is real. The target's sensitivity to including 2021-2024 is
therefore a physical question -- which period should the model represent? --
and not a data-quality one.

Inputs
    annual_medians_detrended.csv    from coastsat_detrended_position.py, via the resolver
    CoastSat time-series CSVs       for the seasonal and sampling tests
    duneline / coastsat endpoint    3-rates/*/endpoint/<window>/, stored
                                    tables only -- nothing is refitted here

Outputs  (hat_observed_rates.DETRENDED_POSITION)
    attribution_tests.csv           one row per test: what it measured, what
                                    it returned, and the verdict
    attribution_nourishment.csv     each fill's local step against the island
    attribution_duneline.csv        both sources over each survey-date pair
    detrended_position_attribution.png                 the three tests that carry a figure

Usage
    python .../coastsat_position_attribution.py        (run coastsat_detrended_position.py first)
```

Notes that were in the code:

```text
The dune-line survey pairs, by the rate-window folder that holds them.
The labels stay ASCII: this script prints them, and a Windows console on
cp1252 raises UnicodeEncodeError on an arrow. Arrows belong in the figure.
```

<details><summary>Function notes (the original docstrings)</summary>

**`test_duneline()`**

```text
Does the independent dune line show the same reversal?

Stored tables only -- both endpoint products are read as written, nothing
is refitted, so this obeys the rule 4-comparisons states for source-against
-source work even though it is filed here with the rest of the attribution.
```

</details>

### mean_shoreline/coastsat_mean_shoreline.py

One averaging window of CoastSat as a line: each transect's mean position, strung into a mean shoreline.

From the script's original header:

```text
coastsat_mean_shoreline.py -- an averaging window of CoastSat, as a line
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

    By default the window is the CALENDAR span, not a span centred on the
    1997 dune survey ([[cascade-period-is-the-calendar-year]]); the mismatch
    is REPORTED in PROVENANCE.md, not corrected.

    Since 2026-09-29 a window can also be given by dates. --centred-on takes
    +/-1 yr of the flights of the lidar a period's start DEM is built on
    (SURVEY_ANCHORS), because the line that becomes the shoreline island
    offset is a snapshot the model starts from beside that DEM, and should
    be dated like it (Hannah, by interview). That exception is for the
    offset only; rates and scoring stay calendar.

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

OUTPUT   data/hatteras_init/5-scr/1-observations/mean_shoreline/<label>/
    <label> is `<start>_<end>` in years for a calendar window (1995_1997),
    in ISO dates otherwise (1995-10-12_1997-10-12).
    shoreline_mean_<label>.geojson         ONE LineString, EPSG:26918, with the
                                           metadata properties a dune line
                                           carries so step 1 of the offset
                                           build reads it unchanged
    transect_means_<label>.csv             per CoastSat transect: n, mean, sd,
                                           se, first/last date, the geolocated
                                           mean point, domain, included/why not
    mean_shoreline_<label>.png             the diagnostic: where the line is,
                                           and how well sampled it is
    mean_shoreline_<label>_island_outline.png
                                           its panel (a) alone, over the
                                           island outline
    PROVENANCE.md
    Read through hat_observed_rates.mean_shoreline_{dir,geojson,csv}().

USAGE
    python coastsat_mean_shoreline.py
    python coastsat_mean_shoreline.py --window 1995 1997 --min-obs 10
    python coastsat_mean_shoreline.py --centred-on alace_1996     # 1995-10-12 .. 1997-10-12
    python coastsat_mean_shoreline.py --centred-on usace_2009     # 2008-08-17 .. 2010-08-17
    python coastsat_mean_shoreline.py --window-dates 1995-10-12 1997-10-12

THEN (the offset build, which this script does not do)
    duneline_to_raw_offsets.py --duneline <the geojson above>
        --out 1995_1997_shoreline_offset_raw.csv
    island_offset_hybrid.py --year 1996 --source shoreline --version v1
        --raw-file <the raw file above>
```

Notes that were in the code:

```text
The CRS the 1984, 2009 and 2023 dune lines declare. The 100 m transects are
EPSG:3725, NAD83(NSRS2007) / UTM 18N, which is the same grid to within the
null transform, and duneline_to_raw_offsets.py reprojects anyway.
```

```text
The house minimum, from coastsat_lrr.MIN_OBS: fewer than ten positions is
not a mean of a seasonal cycle.
```

```text
Copied into the raw offsets file by duneline_to_raw_offsets.LINE_META, so a
shoreline-derived raw file carries as full a provenance as a dune-derived one.
```

```text
The lidar surveys a period's start topography is built on, for --centred-on.
Decided 2026-09-29 (Hannah, by interview): the line that becomes the
shoreline island offset is a SNAPSHOT paired with the start DEM, so it is
averaged over +/-1 yr of the survey's flights rather than over the calendar
span. This is a deliberate exception to [[cascade-period-is-the-calendar-year]]
for the offset only; rates and scoring stay on calendar years. The centre is
the middle of the flights, written out rather than computed so it matches
the dates agreed in the interview.
```

```text
The layer spells an id "usa_NC_0032-0001"; the timeseries file and
the domain lookup both spell it with an underscore.
```

```text
A window, not a moment -- rule 2. The raw offsets file carries `year`
in its `year` column, where "1997" would be a lie. A calendar window
keeps the strings it has always written, so its raw file is unchanged.
```

```text
Two panels STACKED, and the map drawn with northing across the page.
The island is ~45 km north-south by ~6 km east-west, so at equal aspect
-- which a map has to keep, or the shape it is showing is not the shape
on the ground -- a portrait panel is an unreadable sliver. Turned on its
side it is a wide ribbon, and alongshore runs left to right as it does
in every other figure here.
Panel (b) is kept short (Hannah, 2026-09-29: "make panel b thinner");
the figure is shortened with it so panel (a) keeps its size.
```

```text
Easting increases DOWN, as in ribbon_axes: with north to the right, an
easting-up axis draws the mirror image of the map (Hannah, 2026-09-29,
the same flip the outline figure got on 09-23). Ocean at the bottom.
```

```text
The ribbon: panel (a) of the diagnostic on its own, over a map (Hannah,
2026-09-23). Northing across the page, easting up, equal aspect. Shared with
coastsat_mean_shoreline_on_imagery.py, which draws the same ribbon on photos.
```

<details><summary>Function notes (the original docstrings)</summary>

**`Window()`**

```text
An averaging window: first and last day, both inclusive.

Built from calendar years (1 Jan of the first to 31 Dec of the last) or
from two ISO dates. `key` is what the hat_observed_rates resolvers take --
years for a calendar window, so its folder stays `1995_1997`, dates
otherwise, so its folder is `1995-10-12_1997-10-12`.
```

**`transect_geometry()`**

```text
Origin and seaward unit vector of each wanted transect, in TARGET_CRS.

The layer is the global CoastSat file (233k transects, 80 MB), so it is
read once and filtered to the ids the domain lookup names. Chainage is
measured from the FIRST vertex along the line, so that vertex is the
origin and the line's own direction is the unit vector.
```

**`window_means()`**

```text
One row per CoastSat transect: its mean position over the window, and
the positions behind the INCLUDED means counted by calendar year.

`window` is a Window, inclusive of its first and last day. A transect is
EXCLUDED, with the reason recorded, when its geometry or timeseries is
missing or when it holds fewer than `min_obs` positions -- never dropped
in silence.
```

**`line_vertices()`**

```text
The included mean points in alongshore order.

Ordering is (site, transect number), which is alongshore here: the six
Hatteras sites chain south to north with monotonically increasing domain
spans, and the resulting vertex spacing is ~50 m with no gap over 300 m
(checked 2026-09-22). Where two sites overlap at a domain boundary the
line can double back a little; that is why the intersection step records
n_crossings and takes one crossing per transect.
```

**`write_geojson()`**

```text
ONE LineString, with the properties a digitised dune line carries.

duneline_to_raw_offsets.py exits on a file holding more than one feature
and copies LINE_META across, so this has to be a single feature and it is
worth filling the metadata in: the raw offsets file is where a reader
meets this line next.
```

**`ribbon_extent()`**

```text
(n0, n1, e0, e1): the line's northing span, and easting wide enough to
hold the island landward of it (Buxton Woods is ~3 km across).
```

**`ribbon_axes()`**

```text
Km ticks in the same words as the diagnostic's panel (a), but easting
increasing DOWN (Hannah, 2026-09-23: "flip these vertically"): the ocean is
at the bottom and the ribbon is a north-up map turned 90 degrees clockwise,
where panel (a)'s easting-up axes draw its mirror image.
```

**`draw_island_outline()`**

```text
Water tint and the island outline, in the ribbon's swapped axes.

`edge_on_top` (a colour) also draws the outline's edge above everything
at zorder 3.5, between the photographs and the mean line; the imagery
ribbon passes "white" (Hannah, 2026-09-23) so the outline reads against
the photographs.
```

**`_excluded_table()`**

```text
The excluded transects as a markdown table.

Spelled out rather than DataFrame.to_markdown(), which wants `tabulate`;
a provenance note is not worth a dependency (2026-09-22).
```

**`dune_line_date()`**

```text
(vintage, ISO date) of the dune line a period start reads, or None.

The vintage from hat_topo_version.DUNE_LINE_FOR_YEAR, its flight date
from coastsat_vs_duneline.KNOWN_SURVEY_DATES -- neither is typed here.
```

**`_window_paragraph()`**

```text
Point 2 of the provenance: what the window is centred on, and how far
that sits from the period's dune line. Computed, so it is right for every
window (until 2026-09-29 it was 1996's text, printed into 2009_2011 too).
```

**`_sampling_paragraph()`**

```text
Point 3: the positions behind the included means, per calendar year,
counted over every included transect (the 09-22 text quoted a sample of
80).
```

**`period_start()`**

```text
The model period this window starts: the anchor's, else the
calendar year holding the centre (1996 for 1995-1997).
```

**`clip()`**

```text
The observations inside the window. The last day is included WHOLE
(to 23:59:59); until 2026-09-29 the calendar filter stopped at
midnight on 31 Dec, which dropped nothing in the existing windows --
no pass falls on 1997-12-31 or 2011-12-31 -- but would have dropped a
pass on the last day of a date window.
```

</details>

### mean_shoreline/coastsat_mean_shoreline_on_imagery.py

A window's CoastSat mean shoreline on that window's aerial photographs, one panel per flight year.

From the script's original header:

```text
coastsat_mean_shoreline_on_imagery.py -- the window mean, on that window's photographs
The CoastSat mean shoreline for a window (coastsat_mean_shoreline.py) drawn
on the USGS aerial photographs flown inside that window, one panel per
flight year, at a handful of sites along the island. Asked for by Hannah on
2026-09-23 ("the shoreline position imposed on the aerial imagery from that
year").

WHAT IS ON EACH PANEL
    the photograph      the USGS Henderson release (doi 10.5066/P1CXBCDW),
                        read frame by frame from D:\Hatteras_GIS\Aerial
                        through the 1984 imagery review's Imagery class (the
                        seamline rule and the film-fringe handling live there,
                        not here); stated accuracy 1.2 m
    the mean shoreline  the window-mean line, ink with a white halo
    +/-1 sd             the within-window scatter of each transect, placed
                        along the transect's own direction and joined
                        alongshore, as a translucent white band, dashed edges
    the positions       (second version only) every satellite position behind
                        the mean, geolocated the same way (origin + chainage *
                        direction), coloured by date on one scale for the
                        whole window; the scale marks each flight date
                        (added 2026-09-23, Hannah: "with and without the dots",
                        "a gradient to show throughout time")
    the domains         (third version only) the Barrier3D domain boxes
                        (transect_domains/HAT_domains.json), white with a grey
                        edge, each labelled with its GIS number; on the site
                        zooms and the island overview (added 2026-09-28)

WHAT IT CAN AND CANNOT SHOW
    A photograph is ONE October day; the line is a three-year mean of ~28
    satellite passes. The wet/dry line in a photo is not expected to sit on
    the mean -- it is expected to sit, most days, inside the band. A photo
    edge far outside the band at one site is worth looking at; nothing is
    measured from the photographs here.

PHOTOGRAPHS FROM OUTSIDE THE WINDOW (--photo-years, added 2026-09-28)
    By default the photographs are the window's own years. A window with none
    (2009-2011: the 2009 folder is raw Google Earth captures, not georeferenced)
    takes the nearest year instead, e.g. --photo-years 2008 for the NOAA NGS
    mosaic of 26-27 March 2008; the captions then say the photograph is from
    outside the window, and the date scale widens to reach its flight date.
    Sources other than the USGS release are described in PHOTO_SOURCES.

SITES (supporting/sites.csv)
    Each window is three domains (1.5 km) alongshore, centred on the site's
    domain, and every panel of a site shares one extent. Cross-shore it runs
    LAND_M landward and SEA_M seaward of the line. The sites are spread along
    the island and include the two piers, which are fixed in all three photos.

OUTPUT   <mean_shoreline_dir(window)>/on_imagery/
    line_and_band/mean_shoreline_<window>_on_imagery_GIS<NN>_<site>.png
    line_and_band/mean_shoreline_<window>_on_imagery_island_1996.png
        the whole island in three north-up segments (GIS 1-30, 31-60, 61-90)
        at one scale on the 1996 photographs (--island-year; the first photo
        year when 1996 is not among them), the site windows outlined
    line_and_band/mean_shoreline_<window>_on_imagery_ribbon_1996.png
        panel (a) of mean_shoreline_<window>.png alone on the 1996 photographs,
        the island outline beneath where no frame covers
    with_positions/mean_shoreline_<window>_on_imagery_with_positions_GIS<NN>_<site>.png
    with_domains/mean_shoreline_<window>_on_imagery_with_domains_GIS<NN>_<site>.png
    with_domains/mean_shoreline_<window>_on_imagery_with_domains_island_1996.png
        each subfolder with supporting/CAPTIONS.md (no PDFs: raster panels)
    supporting/sites.csv   the windows, both files and the position counts per site
    Also published to output/figures/2-observations/mean_shoreline/<the same subfolders>.

USAGE
    python coastsat_mean_shoreline_on_imagery.py
    python coastsat_mean_shoreline_on_imagery.py --window 1995 1997 --sites 26 79
    python coastsat_mean_shoreline_on_imagery.py --only island   (or sites, ribbon)
    python coastsat_mean_shoreline_on_imagery.py --window 2009 2011 --photo-years 2008
    python coastsat_mean_shoreline_on_imagery.py --centred-on alace_1996
    python coastsat_mean_shoreline_on_imagery.py --centred-on usace_2009 --photo-years 2008
    A window by dates (2026-09-29) is the one coastsat_mean_shoreline.py built
    with the same flags; a photograph is inside it by its flight DATE, and the
    date scale runs over the window itself rather than whole years.
    Needs the D: drive; the .venv Python (rasterio).
```

Notes that were in the code:

```text
The +/-1 sd band is neutral -- translucent white with a dashed ink edge -- so
the only colour on a panel is the positions' date scale. It was yellow until
2026-09-23, when the dates moved to viridis, whose light end is yellow.
```

```text
The positions' date scale (Hannah, 2026-09-23: "academic and professional"):
viridis is perceptually uniform, reads in greyscale and to colour-blind
readers, and is the scale reviewers expect for an ordered variable. The thin
ink edge keeps the light end visible on the brightest beach.
```

```text
One subfolder per version, each with its own supporting/CAPTIONS.md
(Hannah, 2026-09-23); sites.csv covers all and stays in on_imagery/supporting.
The key is the version, the value its subfolder; with_domains added 2026-09-28.
```

```text
The domain boxes follow the island outline's convention in the ribbon figure
(white, grey edge), so they cannot be mistaken for the ink shoreline or the
dashed band edges.
```

```text
Photographs that are not the USGS Henderson release. The imagery reader dates
a year from Henderson frame names only, so the flight date is kept here, read
from the source's own metadata.
```

```text
D:\Hatteras_GIS\Aerial\2008, 2008_IOCM_NaturalColorImagery_J1129187_metadata.xml:
"2008 NOAA NGS Ortho-rectified Color Mosaic from Ocracoke, NC to Virginia
Beach, VA", beginPosition 2008-03-26, endPosition 2008-03-27
```

```text
the island overview (Hannah, 2026-09-23: "across the island in 3 vertical
panels with the 1996 imagery")
```

<details><summary>Function notes (the original docstrings)</summary>

**`_inside()`**

```text
Was this photograph flown inside the window? By date since 2026-09-29,
when windows stopped being whole calendar years; a source whose reader
knows only the year (PHOTO_SOURCES) is dated from its table entry.
```

**`positions()`**

```text
Every satellite position behind the included means, on the ground.

The same window filter as coastsat_mean_shoreline.window_means, so these
are exactly the positions each mean was taken over.
```

**`draw_line()`**

```text
The mean line and its band. At island scale (edges=False) the band's
dashed edges would merge with the line, so only its fill is drawn.
```

**`draw_domains()`**

```text
The Barrier3D domain boxes crossing the window, each labelled at its
landward (west) side, clear of the shoreline.
```

**`site_figure()`**

```text
Three versions of one site: the line and band, the same with the
positions, and the same with the domain boxes.
```

**`island_figure()`**

```text
Three north-up segments side by side at one scale, on one year's photos,
with the site windows of the zoom figures outlined; drawn twice, the
second time with every domain box (the photographs are read once).
```

**`ribbon_figure()`**

```text
Panel (a) of the diagnostic, alone, over one year's photographs: the
same extent and axes as coastsat_mean_shoreline.outline_figure.
```

</details>

### mean_shoreline/coastsat_mean_shoreline_storm_check.py

Was a window mean shaped by a storm? A check on the +/-1 yr means that become the shoreline offset.

From the script's original header:

```text
coastsat_mean_shoreline_storm_check.py -- was a window mean shaped by a storm?
A check on the +/-1 yr mean shorelines that become the shoreline island
offset (1995-10-12 .. 1997-10-12 for 1996, 2008-08-17 .. 2010-08-17 for
2010). A two-year mean still leans on whatever the beach was doing in those
two years, and a big storm inside the window, or just before it, can pull a
share of the passes landward. This puts each window in its longer storm
record and asks three questions:

    1. WHICH STORMS. Every high-water event from 3 yr before the window to
       3 yr after it, ranked against the whole 1984-2024 record.
    2. WAS THE WINDOW STORMY. Storm-hours above the berm in the window
       against every other 2-yr span of 1984-2024.
    3. DID IT MOVE THE MEAN. Each transect's window mean recomputed without
       the passes that fall within --recovery-days after a major storm; the
       shift is what that storm's aftermath contributed to the line.

Built 2026-09-29 (Hannah asked for it as a check inside the mean-shoreline
folders).

THE STORMS ARE THE MODEL'S STORMS
    The events, their peak total water level (Rhigh, m above MHW) and their
    hours above the berm are the hindcast series (hat_env_forcings'
    DEFAULT_STORM_VARIANT): Duck gauge 8651370 + Stockdon (2006) R2% from WIS
    63228 waves, berm 1.7 m NAVD88. The full 1984-2024 record is years
    <= 2003 from the 1984_2004 file and >= 2004 from the 2004_2024 file.
    Peak hours and the tropical/other class come from storm_figures.py, so
    the names and colours match 3-storms/figures/. Water levels are an
    estimate from the gauge and hindcast waves, not observations at
    Hatteras.

"MAJOR" IS A RETURN LEVEL, NOT A PERCENTILE OF EVENTS
    A storm is major when its Rhigh OR its hours above the berm reach the
    MEDIAN ANNUAL MAXIMUM of 1984-2024 -- a level the record reaches in half
    its years, the ~2-yr event. Height alone missed Nor'Ida (Nov 2009,
    2.85 m but 103 h, the 13th longest event of the record), and a long
    nor'easter can move more sand than a higher, shorter storm. A percentile
    of all events would move with how many small events the berm threshold
    admits; the annual maximum does not. --major-rhigh / --major-hours
    override the two levels.

THE SHORELINE SERIES
    Each included transect of the window's transect_means CSV, its CoastSat
    positions over the context span, minus that transect's window mean (so
    zero IS the mean line, positive is seaward). One point per image date:
    the median across transects, drawn only where at least --min-coverage
    of the transects have a position that day (a cloud-clipped scene
    samples one end of the island).

WHAT IS NOT DONE
    No pass is removed from the mean line. The sensitivity in question 3 is
    reported, never applied. Same rule as the producer: the scatter is the
    beach moving.

OUTPUT   <mean_shoreline window folder>/storm_check/
    storm_check_<label>.png                  the storms around the window
    storm_check_<label>_events.csv           every event in the context span,
                                             ranked in 1984-2024
    storm_check_<label>_mean_shift.csv       per transect: the window mean
                                             with and without post-storm passes
    README.md                                the findings, computed
    supporting/  PDF, CAPTIONS.md, the island-median series

USAGE
    python coastsat_mean_shoreline_storm_check.py                 # both DEM-centred windows
    python coastsat_mean_shoreline_storm_check.py --centred-on alace_1996
    python coastsat_mean_shoreline_storm_check.py --window-dates 1995-10-12 1997-10-12
    python coastsat_mean_shoreline_storm_check.py --context-years 5 --recovery-days 60
```

Notes that were in the code:

```text
Beach recovery after a storm runs weeks to months; 90 days is the middle
of that, and it is a parameter because nothing here measures it.
```

```text
Named storms carry their HURDAT2 name; an unnamed storm that is major by
length only carries its hours, since its dot sits below the dashed line.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_events()`**

```text
Every event of 1984-2024, one row each, with its peak hour, Rhigh (m
MHW), hours above the berm before trimming, and storm type.
```

**`storm_hours_2yr()`**

```text
Storm-hours above the berm in every 2-yr span of the record, stepped
monthly: (span start, span end, hours, events, max Rhigh).
```

**`load_positions()`**

```text
Every CoastSat position of the window's included transects over the
context span, as the anomaly from that transect's window mean.
```

**`mean_shift()`**

```text
Each transect's window mean without the passes inside `recovery_days`
after any of `storms`; one column per storm and one for all of them.
```

**`figure()`**

```text
The storms around the window, one panel. The shoreline series and the
2-yr storminess are in the README and the CSVs, not drawn (Hannah,
2026-09-29: keep only the storm panel).
```

**Events after 2024 (2026-10-06)**

```text
The 2025 window needs storms through 2025, so load_events also reads the
2009_2025 series for years after 2024, with the forcing swapped to the
extended Duck and WIS records (identical to the 1984-2024 files where both
have data). The reference is still 1984-2024: major thresholds, ranks and
the storminess spans come from it alone, and a later event is ranked as if
added to it, so the 1996 and 2009 checks are unchanged. The storminess
comparison uses spans as long as the window (1 yr for 2025, 2 yr for the
others). The storm record ends 2025-12-31; the README says so when a
window runs past it.

    python coastsat_mean_shoreline_storm_check.py --window-dates 2025-02-17 2026-02-17
```

**`link_from_provenance()`**

```text
One row in the window's PROVENANCE.md Files table, if it is missing.
coastsat_mean_shoreline.py writes the same row when the folder exists.
```

</details>

### mean_shoreline/coastsat_mean_shoreline_compared.py

The three mean shorelines the model starts from and is graded against, drawn as the island.

Built 2026-10-06. The three windows come from `hat_observed_rates.NET_CHANGE_WINDOWS`,
so the figure always shows the lines the net-change targets difference: the
1996 start (±1 yr of the ALACE flights), the 2009 calibration end and test
start (±1 yr of the USACE flights), and the 2025 test end (±6 months of
2025-08-17; no DEM). Each line is placed in the offshore-datum frame with the
offset build's own intersection (`duneline_to_raw_offsets.intersect`), averaged
per domain, so the 1996 and 2009 stations equal the stored shoreline-offset
raw files (checked to 0.005 m). Colours: 1996 red, 2009 blue, 2025 amber
(purple is the dune line's in the duneline_vs_shoreline figure). No spread
band and no smoothing, by request. One row of six panels, wider than the
house double column, for slides and posters.

    python coastsat_mean_shoreline_compared.py

Close the image in any viewer before re-running: an open file blocks the save.

### shoreline_inventory/shoreline_inventory.py

Every shoreline observation available across the Hatteras study area, by source and year.

From the script's original header:

```text
Shoreline Data Inventory — Hatteras Island study area
Cross-source inventory of shoreline observations available across the
full Hatteras Island study area (all CASCADE domains, Cape Point to
north of Rodanthe). Reports what data you have from each source, how
they overlap in time, and where the gaps are.

The study area is defined by a spatial filter file (see CONFIG) so this
tool can be reused for narrower analyses later (e.g., point it at the
Buxton transects shapefile to focus on the groin study area).

Sources handled
1. User-digitized wet-dry lines           (GeoJSON/shapefile with date attr)
2. NC Coastal Management historical lines (GeoJSON/shapefile with DATE_ attr)
3. CoastSat satellite-derived shorelines  (folder of per-transect CSVs)

Outputs
  shoreline_inventory_by_year.csv    — one row per year, count per source
  shoreline_inventory_observations.csv — one row per (source, date) observation
  shoreline_inventory_timeline.png   — three-panel visual timeline

Usage
Edit the CONFIG section below, then run:
    python shoreline_inventory.py
```

Notes that were in the code:

```text
--- Study area spatial filter ---
Path to a shapefile or GeoJSON defining the study area. Can be either:
• Polygon features (e.g., bounding box, CASCADE domains): merged and
optionally buffered.
• Line features (e.g., transects): merged and buffered by
STUDY_AREA_BUFFER_M to create a filter polygon.
The script auto-detects which case you're using.
Anchored on this file 2026-09-12. The literals here were
drive-rooted and had never resolved; the data they name also
moved out of the scripts tree on that date.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
Buffer distance (m) applied around the filter geometry.
For a bounding box: 500 m to catch drifted historic shorelines.
For CASCADE domain polygons: 500 m.
For transect lines: 1000+ m for a corridor filter.
```

```text
--- Source 3: CoastSat time-series CSVs ---
GeoJSON with CoastSat transect geometry (used to spatially filter)
```

```text
--- Visualization ---
Whether to produce spatial and alongshore visualizations
```

```text
Colormap for year-based shoreline coloring.
'viridis' is perceptually uniform, colorblind-safe, and reads intuitively
as "dark = original, bright = new". Other good options: 'plasma', 'cividis'.
```

```text
Number of representative CoastSat snapshot dates to reconstruct.
The script picks these to span the full temporal record (evenly-spaced bins)
and picks the most-complete date within each bin (highest transect coverage).
```

```text
Basemap for the spatial map view. Requires the 'contextily' library:
pip install contextily
Options:
"none"           — plain background (default, no dependencies)
"esri_satellite" — ESRI World Imagery (aerial photography)
"esri_topo"      — ESRI Topographic (labeled road_offset, terrain)
"carto_light"    — Carto positron (clean, high-contrast for shorelines)
If contextily is unavailable, the script falls back to "none" with a warning.
```

```text
Zoom regions used for the multi-panel detailed sub-region map. Each tuple
is (region_name, y_min_northing, y_max_northing). The four regions cover
the island south-to-north with small overlaps for visual continuity. To
change: edit the values or add/remove regions and the function auto-scales.
```

```text
Geographic labels drawn on the right edge of each map panel.
Values are (name, northing_meters_UTM18N). Approximate positions of major
Hatteras Island communities. TO FIX A LABEL POSITION:
1. Open your data in ArcGIS Pro or QGIS
2. Click on the actual location of the town on the map
3. Read the y-coordinate (northing) from the status bar (should be a
value between ~3,890,000 and ~3,960,000 for Hatteras)
4. Update the tuple below with that value
```

```text
Zoom sub-regions for the "detailed regions" plot. Each region becomes one
panel in the multi-panel figure, at much higher zoom than the full-island
view. Small overlap between adjacent regions is intentional so shorelines
spanning a boundary appear in both.
Format: (region_name, y_min_meters_UTM18N, y_max_meters_UTM18N)
```

```text
Compute the length of each shoreline feature within the study area
(proxy for alongshore extent of that observation)
```

```text
Aggregate by unique date — many shoreline shapefiles store single
shorelines as multiple line segments (topologically clipped, split at
grid edges, etc.). Counting each segment as a separate observation
would inflate high-density years like 1980 (~140 segments for one
shoreline). We collapse to one row per date and sum the segment lengths.
```

```text
Floor to day using UTC to collapse same-day satellite passes, then strip
the timezone so these dates can be concat'd/sorted alongside tz-naive
dates from the shapefile sources.
```

```text
─── Figure and axes ──────────────────────────────────────────────────
Independent y-axes on each panel (sharey with equal-aspect + geopandas
is a known incompatibility — the axes have to be independent).
```

```text
─── Use the study area polygon bounds so every panel shows the same
extent as the input bounding box. Extra right padding gives labels
room outside the island's east edge, and importantly makes the plot
aspect wide enough that panels fill their grid cells rather than
shrinking to a thin strip (which would create empty gaps between them).
```

```text
No polygon backdrop is drawn — the plot extent alone reflects the
bounding-box limits set by the user's cascade_area.geojson.
```

```text
Location labels on the right side of every panel, in the reserved
right padding area (outside the island's footprint). Subtle white
background with rounded corners keeps text readable if any part
happens to fall over a shoreline feature.
```

```text
Helper to batch-plot a GeoDataFrame by unique year in one geopandas
call per year (much faster and avoids many aspect resets).
```

```text
─── Continuous year colorbar with tick marks at observed years ───────
Uses the rank-based color mapping (each unique year gets equal color
space) but presented as a smooth continuous bar rather than swatches.
Ticks are placed at every N-th observed year to avoid overlap.
```

```text
─── Plot each source with distinct style ─────────────────────────────
NC state — medium dashed (drawn oldest first so newer on top)
```

```text
─── Figure ───────────────────────────────────────────────────────────
Wider than tall since panels are laid out horizontally. Height chosen
so each panel gets enough vertical resolution for detail.
```

```text
Compute x-extent from actual shoreline data within this y-range
so the panel is tightly zoomed on the shoreline itself.
```

```text
─── Plot each source in the region ───────────────────────────────
NC state (bottom layer, dashed)
```

```text
Reload shapefile geometries with parsed year (used by the map view).
Cheap since these files are small.
```

```text
Multi-region zoomed detail view — the island split into sub-regions
so individual shorelines are actually visible per region.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_study_area_filter()`**

```text
Load a shapefile/GeoJSON defining the study area and return a
single buffered polygon used as the spatial filter for every source.

Auto-detects whether the input contains polygon features (CASCADE
domains) or line features (transects), and handles each accordingly:
  - Polygons: unioned, then buffered by `buffer_m`.
  - Lines:    buffered by `buffer_m`, then unioned.
  - Points:   buffered by `buffer_m`, then unioned.
```

**`parse_date_column()`**

```text
Parse a date column that could be in one of several formats:

- String dates: '1972-07-30', '7/22/1998', 'July 30, 1972', etc.
- Numeric milliseconds since Unix epoch: 81302400000 → 1972-07-30
  (This format is produced by some ArcGIS / QGIS exports.)
- Numeric year-only integers: 1972 → 1972-01-01

Returns tz-naive datetime series with NaT for unparseable values.
```

**`load_shapefile_dates()`**

```text
Load a shoreline shapefile/GeoJSON, spatially filter to the study area,
and extract observation dates.

Returns a DataFrame with columns:
  source, date, year, extent_m
where `extent_m` is the length of the shoreline within the study area
(in meters), giving a proxy for alongshore coverage per feature.
```

**`collect_csv_map()`**

```text
Walk one level of subfolders under root_dir and return
{csv_stem: full_filepath} for every CSV found.
```

**`load_coastsat_dates()`**

```text
Load CoastSat observations, filter to study area, and aggregate.

Returns:
    per_date_summary : DataFrame with columns source/date/year/n_transects_covered
    raw_observations : DataFrame with columns transect_id/date/chainage_m
                       (used later for shoreline reconstruction; empty if none)
    transect_geoms   : GeoDataFrame of study-area CoastSat transects
                       with columns transect_id/geometry (LineString in PROJECTED_CRS)
                       (used later for reconstructing shoreline points)
```

**`pick_representative_coastsat_dates()`**

```text
Pick N CoastSat dates that (a) span the full temporal record and
(b) each have high transect coverage.

Method: divide the time span into N evenly-spaced bins. In each bin,
pick the date with the highest n_transects_covered.
```

**`get_transect_direction()`**

```text
Return unit direction vector (dx, dy) from a transect LineString's
origin (first vertex) to its endpoint (last vertex).
```

**`reconstruct_coastsat_shoreline_points()`**

```text
For each (date, transect_id) in target_dates, compute the shoreline
point (x, y) in the projected CRS by walking along the transect from
its origin for `chainage_m` meters.

Returns DataFrame with columns:
  date, transect_id, x, y, chainage_m, origin_x, origin_y
```

**`compute_alongshore_position()`**

```text
Assign each CoastSat transect an alongshore position (m) — cumulative
distance from the southernmost transect origin, walking north via
nearest-neighbor.

Returns DataFrame: transect_id, origin_x, origin_y, alongshore_m
```

**`_compute_data_bounds()`**

```text
Return zoom bounds from the union of all shoreline data, with a small
buffer. Falls back to the filter polygon bounds if no data provided.
```

**`plot_shoreline_map()`**

```text
Four-panel horizontal map view of shoreline coverage over time.

Design decisions:
  - Panels share y-axis (only leftmost shows northing labels)
  - Named locations labeled as secondary y-tick labels on leftmost panel
  - Axes in km (raw UTM /1000) rather than meters with 1e6 raw_offset
  - Rank-based year coloring so historic (1849) and modern (2024) years are
    both distinguishable rather than compressed into narrow color bands
  - Zoom to actual data extent + horizontal buffer rather than the full
    study area polygon (island is much narrower than the polygon)
  - Optional satellite / topo / positron basemap via contextily
```

**`plot_shoreline_map_detailed()`**

```text
Standalone large-format all-sources overlaid map for detailed
inspection. Same data as the 4-panel view but in one big panel where
the shorelines are individually visible.

Sources distinguished by line style + width:
  wet-dry   → thick solid   (highest precision, user-digitized)
  NC state  → medium dashed (historical survey data)
  CoastSat  → thin dotted   (satellite-derived, densely sampled)

Year encoded by rank-based viridis coloring.
```

**`plot_shoreline_map_zoomed_regions()`**

```text
Multi-panel high-detail map. The island is split into ZOOM_REGIONS
(config), each rendered as its own panel with its own tight x/y limits.
Because each panel is zoomed in on a small y-range AND to the actual
shoreline data on x, individual shorelines are visually separable.
```

**`plot_alongshore_profile()`**

```text
CoastSat shoreline chainage vs. alongshore position, colored by year.

X-axis: alongshore position (km, south to north)
Y-axis: chainage (m from transect origin, positive = seaward)
```

**`load_shapefile_geoms_for_plot()`**

```text
Reload a shapefile and return the filtered geometries + parsed year.
Used by the map-view plotting step (separate from the inventory step).
```

</details>

### shoreline_patterns/shoreline_trajectory_classification.py

Classify each domain's shoreline trajectory from the CoastSat transects, over two periods and the full record.

From the script's original header:

```text
Classifies shoreline trajectory stability across Hatteras Island domains
using CoastSat transect time-series over three periods:
    - Period 1 (calibration): 1984–2004
    - Period 2 (validation):  2004–2024
    - Full record:            1984–2024

For each domain and period, computes:
    - Linear rate of change (LRR) via OLS on annual median chainage
    - Sign consistency (fraction of year-on-year steps in dominant direction)
    - Interannual variability (std dev of annual positions around trend)

Classification scheme (applied per domain, per period pair):
    Persistent Erosion      — both periods erosional  (LRR < -THRESHOLD)
    Persistent Accretion    — both periods accretional (LRR > +THRESHOLD)
    Persistently Stable     — both periods within ±THRESHOLD
    Switching: Acc→Ero      — full reversal: P1 accretional, P2 erosional
    Switching: Ero→Acc      — full reversal: P1 erosional,   P2 accretional
    Decelerating Erosion    — P1 erosional, P2 stable (erosion slowing)
    Accelerating Erosion    — P1 stable, P2 erosional (recently destabilised)
    Decelerating Accretion  — P1 accretional, P2 stable (accretion pulse fading)
    Accelerating Accretion  — P1 stable, P2 accretional (recently gaining)

Outputs
1.  Along-island classification bar chart (domain × period)
2.  Period 1 vs Period 2 LRR scatter plot (per domain, coloured by class)
3.  Hovmöller heatmap (domain × year annual deviation) with classification overlay
4.  CSV summary table of all metrics

Usage
    python shoreline_trajectory_classification.py

Dependencies
    pip install pandas numpy matplotlib scipy tqdm
```

Notes that were in the code:

```text
Anchored 2026-09-14: this named a home directory, or a tree renamed since.
Rule 5 of ORGANIZATION.md.
```

```text
The three paths below were driveless (str(_PATH_REPO / "scripts" / "..."), resolving to
C:\scripts) and one named the pre-2026 "input_preperation" folder, so
this script could not run. Anchored on the repo root (2026-09-10).
```

<details><summary>Function notes (the original docstrings)</summary>

**`compute_lrr()`**

```text
OLS linear rate of change (m/yr) from annual median series.
Returns (lrr, r2, n_years) or (nan, nan, 0) if insufficient data.
```

**`sign_consistency()`**

```text
Fraction of year-on-year steps in the dominant direction (0.5–1.0).
1.0 = monotonic, 0.5 = random walk.
```

**`classify()`**

```text
Classify trajectory from Period 1 and Period 2 LRRs.

Uses 8 classes capturing both direction and acceleration/deceleration:
  - Persistent:    same side of threshold both periods
  - Switching:     full reversal across threshold
  - Decelerating:  was outside threshold, now inside (slowing down)
  - Accelerating:  was inside threshold, now outside (speeding up)
```

**`add_geo_annotations()`**

```text
Add community spans, Wimble Shoals, piers, groin, village lines.
orientation: 'horizontal' (domain on x-axis) or 'vertical' (domain on y-axis).
```

**`compute_domain_metrics()`**

```text
For every domain, compute LRR, sign consistency, and variability
for Period 1, Period 2, and the full record.
Returns a DataFrame indexed by domain number.
```

</details>

### shoreline_patterns/shoreline_trajectory_map.py

Maps of the island with the domains coloured by trajectory class or LRR.

From the script's original header:

```text
Geographic maps of Hatteras Island with CASCADE domains coloured by
shoreline trajectory classification or LRR magnitude.

CONFIG options:
  USE_SATELLITE = True   → Esri WorldImagery satellite basemap (requires internet)
  USE_SATELLITE = False  → plain ocean-blue background (no internet needed)

Dependencies
  pip install geopandas contextily shapely matplotlib pyproj
  pip install matplotlib-scalebar   (optional — for accurate scale bar)
```

Notes that were in the code:

```text
Resolved through hat_observed_rates.py (2026-09-18); the typed paths named
coastsat_lrr/ and used _PATH_REPO before it was defined.
```

```text
Anchored on this file 2026-09-12. The literals here were
drive-rooted and had never resolved; the data they name also
moved out of the scripts tree on that date.
```

```text
Anchored 2026-09-14: absolute into a home directory, or into a tree
renamed since. Rule 5 of ORGANIZATION.md.
```

```text
Web Mercator (EPSG:3857) bounds — used to set axis limits for satellite tiles
(computed from WGS84 bounds above via pyproj)
```

<details><summary>Function notes (the original docstrings)</summary>

**`plot_delta_lrr_map()`**

```text
Single-panel map showing the change in LRR between Period 1 and Period 2.
ΔLRR = P2 LRR − P1 LRR
  Positive (blue) = beach accreted more / eroded less in P2 than P1
  Negative (red)  = beach eroded more / accreted less in P2 than P1
```

</details>
