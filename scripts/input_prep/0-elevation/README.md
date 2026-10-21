# 0-elevation - the initial topography Barrier3D starts from

Turns the 2009 USACE DEM into the `domain_<N>.npy` arrays
`HAT_dune_topo_extractor.py` reads, filling the gaps the 2009 survey left with
measured ground from a second survey.

**Current fill source: 2014 NOAA Post-Sandy.**

## What this is actually fixing

Not "voids in a surface". Each domain contains exactly two nodata regions and
both touch the domain edge - interior enclosed nodata is 0 cells. The 2009
survey simply stops at the waterline on each side, and the side that matters is
the sound.

`roadway_manager.bulldoze` drowns a roadway when >20% of the cells bordering it
sit at or below 0 m MHW, **and a no-data cell passes that test**. In GIS
78/79/80 the row landward of NC-12 was 17-25 no-data cells and zero genuinely
wet ones, so all three roadways width-drowned at t=0 on missing coverage alone -
while the profiles were still 0.5-0.7 m *above* MHW.

So this fills measured ground the 2009 survey missed. It does not invent
elevation anywhere.

## Layout

```
1-source-selection/   which survey to fill from - run before committing
2-produce/            the chain that makes the model input
3-figures/            review figures
```

Every script, by folder:

```
1-source-selection/HAT_survey_dem_coverage.py   score candidate DEMs against the 2009 dry-land gaps
2-produce/HAT_dem_gap_fill.py                   step 1, 2009-start: clip 2009, fill gaps from 2014
2-produce/HAT_dem_1984_mosaic.py                step 1, 1984-start: 2009 + 2014 + the 1996 ALACE overlay
2-produce/HAT_dem_resample_clip.py              step 2: 1 m clips -> the 10 m Barrier3D grid
2-produce/HAT_export_to_numpy.py                step 3: 10 m rasters -> the extractor's .npy arrays
3-figures/HAT_plot_gapfill.py                   review: the 2009-start fill
3-figures/HAT_plot_1984_mosaic.py               review: the 1984-start mosaic and beach-start shift
3-figures/HAT_plot_dem_holes.py                 review: nodata and sub-MHW holes in a start DEM
```

## Datasets

| role | dataset |
|---|---|
| base DEM | `elevation\dem\2009\noaa_usace2009_nc_dem_Job1076020\2009_full.tif` - 1 m, EPSG:3725 + NAVD88 (m) |
| fill source | `elevation\dem\2014\noaa_2014_NOAA_Post_Sandy_DEM_Job1076021\2014_full.tif` - 1 m, EPSG:6347 + NAVD88 (m), reprojected on read |
| domains | `D:\Hatteras_GIS\domains.geojson` - 90 boxes, 2000 x 500 m |
| NC-12 lines | `4-mgmt-forcing/road_offset/raw_offset/{1984,2004}/nc12_*.geojson` - EPSG:2264 (NC State Plane, US survey FEET), reprojected and clipped on load |

## The chain, in order

Each script's output is the next one's input. Run from anywhere - they locate
the project root by walking up for `data/hatteras_init`, not by counting
parents.

```
2-produce/HAT_dem_gap_fill.py       clip + fill   -> 0-elevation/<PRODUCT>/1-gapfill-1m/
2-produce/HAT_dem_resample_clip.py  1 m -> 10 m   -> 0-elevation/<PRODUCT>/2-resampled-10m/
2-produce/HAT_export_to_numpy.py    -> 1-barrier3d-domains/2009-raw/2009-npy-arrays/
                                       2009_pea_hatteras_filled{,_survey}/
3-figures/HAT_plot_gapfill.py       -> 0-elevation/<PRODUCT>/figures/
```

Only the first three produce model input. `HAT_plot_gapfill.py` is review only -
nothing downstream reads it, but it is the one check that the fill landed where
you think it did.

## Which scripts are conditional, and why

Nothing in `1-source-selection/` runs in a normal rebuild. The one script there
exists to choose and justify the fill source.

| script | when |
|---|---|
| `HAT_survey_dem_coverage.py` | choosing among **gridded DEMs**. Scores every candidate against the 2009 dry-land gaps across all 90 domains, with the hydro-flattening check. This is the script that picked 2014, and the script to re-run when a new candidate appears. |

### The point-cloud path was removed on 2026-08-26

Two scripts and one branch went together, because they only ever worked
together:

| removed | was |
|---|---|
| `1-source-selection/HAT_check_fill_source.py` | pre-flight for a point-cloud candidate - header-only stage 1, coverage stage 2 |
| `1-source-selection/HAT_laz_ground_classify.py` | SMRF ground classification, for a cloud that ships unclassified as 2008 did |
| `HAT_dem_gap_fill.py`: `FILL_SOURCE_TYPE`, `FILL_POINTS_*`, `PointCloudSource` | ~75 lines that consumed the classifier's output |

They served the 2008 NOAA IOCM attempt, which lost to 2014 and whose product
folder is no longer on disk. Keeping the classifier without the branch, or the
branch without the classifier, would have left a path that reads as live and
cannot run - so both went, and with them the `2008_NOAA_IOCM` entries in
`scripts/site_layer/hat_elevation_products.py` and in `HAT_plot_gapfill.py`'s `SOURCES`.

**What this costs.** The 2008 comparison figures can no longer be regenerated,
and `HAT_road_elevation.py` - which deliberately samples the road surface from
the 2008 fill rather than the 2014 one - can no longer re-run. Both were
already blocked by the missing rasters; this makes the block explicit rather
than a path that resolves to nothing. See the note at `FILL_SOURCE` in that
script: the RoadElevation CSVs it produced are committed and unaffected, and
what to regenerate them from is an open provenance decision.

**A point-cloud candidate is still allowed.** It has to arrive as a DEM. Grid
it outside this repo, then register it as a product like any other - which is
what `HAT_dem_gap_fill.py` already assumes about every source it reads.

The reasoning the pre-flight script encoded is worth keeping even though the
code is gone: a topographic-only survey stops at its own waterline, so it
cannot see the wet ground landward of NC-12 that this whole chain exists to
fill. 2008 bottomed out at -1.33 m NAVD88 and added 17 cells of 150 in the
strip that mattered. The z range in the header answers that in about a second,
before anything is classified or filled.

## Why 2014 NOAA

Every candidate scored against the 2009 **dry-land** gaps - cells 2009 missed
where the consensus of candidates puts the ground above MHW - across all 90
domains:

| dataset | coverage | |
|---|---|---|
| 2014 NCFMP | 100.00% | **disqualified** - hydro-flattened, 70.5% of its gap values are the single constant -0.762 m |
| **2014 NOAA Post-Sandy** | **97.34%** | earliest genuine, and the best |
| 2017 USACE | 23.74% | |
| 2016 post-Matthew | 21.98% | |
| 2019 DUNEX | 21.98% | |
| 2018 post-Florence | 9.74% | |

2016/2017/2019 score 95-100% on domains 78-80 but are localised surveys - 2019
is below 50% in 60 of 90 domains. Choosing on the developed reach alone would
have picked a dataset covering under a quarter of the island's gaps. That is why
the scoring runs island-wide.

The consensus is the **median** of every candidate covering the cell,
deliberately not any single dataset: using one year as the "is this land?"
reference makes that year score 100% by construction.

Superseded point-cloud attempts, both topographic-only, in the strip landward of
NC-12 (of 150 cells): 2008 NOAA IOCM 17, 2011 post-Irene 28. Their output is
under `1-gapfill-1m/superseded/`.

## Fill rules - selection only, no value is modified

1. **Coverage** - the fill DEM has a value
2. **Connectivity** - contiguous with the island's valid 2009 surface, bridging
   gaps <= 20 m, computed on a *buffered* window so marsh that connects just
   outside the domain is not severed by the crop
3. **Elevation** - >= -2.64 m NAVD88, the extractor's own `WATER_CLAMP_M`, as a
   guard rather than a filter. Flooring at MHW would undo the extractor's
   "keeps back-barrier marsh" decision one step upstream, invisibly.
4. **Vertical** - nothing applied. Bias correction and feathering are **off**, so
   a filled cell is the fill-source measurement unchanged. Both bias estimates
   are still computed and written to the audit every run.

Everything each rule rejects is counted per domain in the audit CSV. Nothing is
dropped silently.

## Result

| | 2008 attempt | 2014 NOAA |
|---|---|---|
| nodata recovered | 23.6% | **43.7%** |
| filled cells, 1 m | 13,823,194 | **25,591,292** |
| filled cells, 10 m | 44,179 | **254,760** |

82.1% of the fill is below MHW - the 2014 DEM covers the sound, and the rules
admit anything above -2.64 m. Set `FILL_MIN_ELEV_NAVD = MHW_ELEVATION` for a
dry-land-only fill (~5M cells).

## Units - verified

All stages metre / NAVD88. The empirical check: median(2009 - 2014) on 3,053,999
co-measured cells = **-0.062 m**. A unit mismatch would show metres of
disagreement.

The seam where fill meets measured 2009 ground is real, measured, and
deliberately left uncorrected - the numbers and the reasoning are in
[`0-elevation/FIGURES.md`](../../../data/hatteras_init/0-elevation/FIGURES.md).

## To change fill source

The tag and fill year appear in three scripts and must agree, or the survey
rasters and the figure legend will disagree about what a filled cell is:

| file | set |
|---|---|
| `scripts/site_layer/hat_elevation_products.py` | **add the product first** - a `Product` entry and its `FILL_CODES` |
| `2-produce/HAT_dem_gap_fill.py` | `FILL_DEM_PATH`, `FILL_SOURCE_TAG`, `FILL_SOURCE_YEAR`, `PRODUCT_TAG` |
| `2-produce/HAT_dem_resample_clip.py` | nothing - pass `--product <NAME>` |
| `2-produce/HAT_export_to_numpy.py` | `SURVEY_FILL`, `DEM_NAME`; pass `--product <NAME>` |
| `3-figures/HAT_plot_gapfill.py` | add one entry to `SOURCES`, then `--source <NAME>` |

The product name is the TOP-LEVEL folder and appears in the figure filenames,
so a new product cannot overwrite an existing one or its figures. Paths are
resolved by `scripts/site_layer/hat_elevation_products.py`; do not rebuild them by joining
strings, which is how `HAT_road_elevation.py` stopped resolving when 2008 moved
under `superseded/`.

Two things are **not** tagged: `DEM_NAME`, the extractor's input folder - point
the extractor at the same name - and the audit CSVs, which are written to fixed
names (`gapfill_audit.csv`, `resample_audit.csv`, `export_audit.csv`) inside
each tag's own directory.

## Before committing to a new source

Run the relevant selector first. On domains 78-80 alone, 2016/2017/2019 all
looked like fine choices at 95-100%, and island-wide they cover ~22%.

## The 1984-start product

A second, parallel product for the 1984 hindcast start. It does not touch the
2009-start chain above: different tag, different output tree, and
`HAT_dem_gap_fill.py` is unmodified.

```
2-produce/HAT_dem_1984_mosaic.py     -> 0-elevation/2009-2014-1996/1-gapfill-1m/
2-produce/HAT_dem_resample_clip.py --product 2009-2014-1996
3-figures/HAT_plot_1984_mosaic.py    -> 0-elevation/2009-2014-1996/figures/
```

`HAT_dem_resample_clip.py` now takes `--product` instead of a hand-edited
constant, and its survey downsampler handles more than one fill code. Run with
no arguments it still produces the baseline product exactly as before -
verified byte-identical on 300 synthetic survey rasters.

### What it adds

A **fifth stage** that OVERWRITES measured 2009 ground with the 1996
NOAA/NASA ALACE survey, wherever ALACE has data. The four rules above only ever
write into nodata; this one replaces measurements, which is a different
operation and carries its own justification in the script's docstring.

    wherever 1996 has data :   1996  >  2009  >  2014
    everywhere else        :           2009  >  2014

**There is no road boundary (since 2026-08-26).** The override used to be
confined to the ocean side of `nc12_1978.geojson`. Measurement showed that
boundary was buying almost nothing: switching it to the 2004 alignment would
have recovered just 2,050 cells of new land island-wide, while removing it
recovers 235,563 and gives domains 1-7 their first 1996 data. ALACE stops
itself 429-979 m from the ocean edge against an island 1274-1999 m wide, so
the survey - not the road - was always the binding landward limit. Both
alignments are still rasterized and reported per domain; neither gates
anything. Full numbers in the product README.

### The boundary is not where the request says it is

ALACE surveyed "from the low water line to the landward base of the sand
dunes". Measured: it covers 30-84% of the ocean-side band, **0-16% landward of
the road**, and reaches 70-133 m further seaward than the 2009 waterline. So
this is a 1996 beach and foredune grafted onto a 2009 backdune and interior,
and **the seam lands at the dune toe, not at the road** - which is where
`HAT_dune_topo_extractor.py` picks its dune windows. Panel C of
`HAT_mosaic1984_zoom_76_81.png` shows it directly: a blue 2009 strip sits
between the road and the orange 1996 band.

That "0-16% landward of the road" is now **admitted rather than discarded** -
2,723,389 cells, 18.5% of everything 1996 writes. It is backdune, not sound:
the swath ends on its own far short of the bay shore in every domain measured.
The audit reports it per domain as `n1996_landward_of_1984_road` and
`n1996_landward_of_2004_road`.

### Three rules ALACE needed that the 2009 chain does not

| rule | value | why |
|---|---|---|
| split floor, into a gap | -2.64 m NAVD88 | no other survey saw the cell; same admission 2014 gets |
| split floor, to replace | MHW, 0.36 m | 2009 **or 2014** saw it. One survey may replace another's measurement, but not with a wet return |
| ceiling | 12 m NAVD88 | rejects the uncorrected-return tail |

Both were found by failure, not foresight, and both guards are still in the
script:

* **The wet edge.** ALACE runs to the low water line, so it carries swash
  returns between -2.64 m and MHW. At a single -2.64 m floor those overwrote
  dry beach and moved 33 of 83 domains' beach start *landward*, up to 27 m -
  while 1996 measured *higher* than 2009 in the overlap. The first fix tested
  only `~valid09`, which still left 16 movers (D73 -21 m, D90 -16 m) because
  the 2009 holes are not empty in the product: 2014 fills them. Testing
  coverage by **any** other survey takes every whole-cell landward mover to
  zero. The 11 that remain are all sub-cell, at most 6 m of a 10 m cell.
* **The ceiling.** The 1996 grid runs to **+256.65 m**. With no ceiling the
  first build put 250 m spikes into 66 of 90 domains against a 2009+2014
  island-wide max of 10.18 m. The floors only ever guarded the low side.

Nothing is clamped or shifted by either - a rejected cell falls through to
2009, then 2014, and every rejection is counted per domain in the audit.

### Vertical datum is inferred, not declared

The .prj carries the horizontal CRS only; the vertical rests on `geoid18` in
the NOAA Digital Coast path the metadata cites. What actually supports it:
median(2009 - 1996) over co-measured cells is -0.11 to -0.36 m through the
developed reach, i.e. 1996 sits decimetres *higher*, which is the sign 12 years
of erosion predicts. A NGVD29 or ellipsoid mistake would show metres. Domains
11, 85 and 86 report ~-1.6 to -1.9 m on slivers of overlap - not a datum
finding. Bias is computed and written to the audit every run; none is applied.

### Domains 1-7 - resolved 2026-08-26

Both NC-12 exports end at domain 8, and 1-7 sit 265 m to 3.3 km beyond the
terminus, so there was no line to define an ocean side for them. Under the old
boundary they got no 1996 at all and carried a 2009 beach while 8-90 carried a
1996 one - an alongshore discontinuity in the initial condition at the 7/8
seam. The stubbed alternative, `NO_ROAD_POLICY = "extend"`, would have invented
a boundary by carrying domain 8's alignment 3.3 km south across Hatteras
village; it was never implemented.

Dropping the boundary dissolves the problem rather than trading it. 1-7 now
receive 1996 on the same terms as everywhere else - 77k to 175k cells each,
6k to 60k of it ground no other survey saw. `NO_ROAD_POLICY` is gone with the
boundary it configured. `road_line=False` survives as a purely descriptive
flag: these domains have no NC-12 export, not a different treatment.

### What moves downstream, and what does not

**Moves: the cross-shore window origin.** `start_beach` is the first cell above
0.50 m MHW, so 1996 beach above that threshold slides each domain's window
seaward. 36 of 90 domains move at least one 10 m Barrier3D cell seaward, none
move a whole cell landward, and the developed reach 76-81 moves 47-72 m - 5 to
7 cells. Per-domain before/after is in `mosaic_1984_audit.csv` as
`start_beach_before_m` / `_after_m` / `_shift_m`, and drawn in
`HAT_mosaic1984_shift.png`. **A new pick set is therefore required**: the
2009_v5 dune windows were picked against a different origin.

**Does not move: `shoreline_offset`.** That is
`2-brie-offset/1984/Island_Dune_Offsets_1984_CASCADE_Input.csv`,
measured independently of any DEM, already existing for 1984, passed straight
to `Cascade()`. Nothing here touches it - a deliberate choice to leave the two
independent rather than derive one from the other.

### Downstream state

`HAT_export_to_numpy.py` has been run against this tag: 90 domains sit in
`1-barrier3d-domains/1984-start/1-extraction/npy-arrays/`, rebuilt with the product itself
after the road boundary came out. `1984-start/dune-topo/v1` exists. What the
"Still to do" note here used to say - that nothing downstream consumed the 1984
product - is no longer true.

## Still open

The roadway drowning in domains 78-80 is a separate matter. The strip landward
of NC-12 is genuinely below MHW in 79 and 80 (2017 and 2019 both agree), so
those roads border real water and drowning there may be correct rather than an
artifact. Only domain 78 showed dry marsh behind the road.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-source-selection/HAT_survey_dem_coverage.py

Score every candidate DEM year by how much of the 2009 DEM's dry-land gaps it measures.

From the script's original header:

```text
Which DEM year can fill the 2009 base DEM's DRY-LAND gaps, and what is the
earliest one that does?

WHAT THE TARGET IS, AND WHY IT CHANGED
The 2009 DEM's nodata is not one thing. Two very different populations sit
inside it, and conflating them produced a wrong answer once already:

  SOUND MARGIN     west of the community, genuinely below MHW. 2017 and 2019
                   agree it is 0% above MHW out to ~600 m from the island. No
                   survey should "fill" this - it is water, and the model
                   should see water.
  DEVELOPED GAPS   holes inside the community itself, where the 2009 survey is
                   only 64-80% complete across parts of the island. Both 2017
                   and 2019 show these at ~1.2-1.5 m NAVD88, 100% above MHW.
                   This is real dry land the 2009 survey missed, and it is what
                   is worth filling.

So the target here is DRY-LAND GAPS only: cells where 2009 has no value and the
consensus of the candidate DEMs says the ground is above MHW.

The consensus is the MEDIAN of every candidate covering the cell, deliberately
not any single dataset. Using one year as the "is this land?" reference makes
that year score 100% by construction, because cells it does not cover are
excluded from the target it is then measured against.

AUTHENTICITY CHECK - COVERAGE IS NOT THE SAME AS MEASUREMENT
A gridded DEM product may be hydro-flattened and void-filled, in which case it
reports coverage everywhere without having measured anything. 2014 NCFMP scored
100% here and is disqualified for exactly this: 85.9% of its values in the gap
are the single constant -0.762 m, another 8.7% are -0.914 m, so 94.6% of its
"coverage" is two stamped water surfaces. Genuine surveys show their most common
value in ~0.2% of cells.

Every candidate is therefore scored for value repetition, and any dataset whose
top value exceeds FLAT_FRACTION_WARN of its gap coverage is flagged.

    python HAT_survey_dem_coverage.py [--domains 78,79,80]

Requires: rasterio, geopandas, numpy, pyproj
```

Notes that were in the code:

```text
source-selection/, where the committed copies are. This wrote to a pooled
0-elevation/figures/ that the 2026-08-25 product-first inversion removed,
so a re-run recreated that folder beside the real one (fixed 2026-09-18).
```

```text
per-domain, NOT the running total c["cov"] - writing the
cumulative sum here made every row look like a monotonic climb
```

```text
n == 0 means the domain has no dry-land gap at all; percentages are
undefined there, not 0%, or the summary reads as 26 failing domains
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data/hatteras_init.

NOT parents[N]. This file moved into 1-source-selection/ on 2026-08-25, and the
old parents[3] then resolved to input_prep/ rather than the project root.
That raises nothing - it just makes every path below it wrong, silently,
until some glob comes back empty. Same helper and same reason as
4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
```

**`read_on_grid()`**

```text
WarpedVRT rejects boundless reads, so partial overlap is intersected and
pasted at the right offset. Errors are not swallowed - an earlier version
caught everything and silently scored every dataset 0%.
```

</details>

### 2-produce/HAT_dem_1984_mosaic.py

Build the 1984-start DEM: 2009 + 2014, with the 1996 ALACE beach and foredune laid over them.

From the script's original header:

```text
Step 1 of 3 for the 1984-START topography. The 2009-start product is made by
HAT_dem_gap_fill.py and is NOT touched by this script; the two write to
different tags and neither can overwrite the other.

WHAT THIS DOES THAT HAT_dem_gap_fill.py DOES NOT
HAT_dem_gap_fill.py only ever writes into NODATA. Its four rules SELECT which
measured cells are allowed into a hole; a cell the 2009 survey saw is never
changed.

This script adds a fifth stage that OVERWRITES measured 2009 ground with the
1996 NOAA/NASA ALACE survey. That is a different operation and it needs its own
justification, which is below.

    wherever 1996 has data :   1996  >  2009  >  2014
    everywhere else        :           2009  >  2014

There is NO road boundary. Until 2026-08-26 the override was confined to the
ocean side of the 1984 NC-12 alignment; it no longer is. The landward limit is
the ALACE swath's own edge. Why, and what it cost, is the next section.

WHY THERE IS NO ROAD BOUNDARY
The original design used "ocean side of NC-12" as the landward limit, on the
reasoning that a rule you can state is better than a survey artefact. Choosing
WHICH NC-12 - the 1984 or the 2004 alignment - forced the question of what the
boundary was actually doing, and it turned out to be doing very little.

All three options were measured across all 90 domains by running stage 1 three
times with only the mask changed, every other rule identical:

    boundary        1996 cells written      new land      overwrote measured
    1984 line             11,983,935       2,084,834               9,899,101
    2004 line             12,251,459       2,086,884              10,164,575
    none                  14,707,324       2,320,397              12,386,927

Read the middle column. Moving 1984 -> 2004 buys 267,524 cells of which only
2,050 are ground no other survey saw. 99.2% of that gain is 1996 replacing a
2009 or 2014 measurement in the backdune - the part of the ALACE swath with
the WORST coverage. It re-vintages a strip; it recovers almost nothing.

Dropping the boundary buys 2,723,389 cells, 235,563 of them new land - 115x
the recovered ground of the 2004 line - and it is the only option that gives
domains 1-7 any 1996 at all.

The two objections to dropping it both failed on measurement:

  * IT DOES NOT REACH THE SOUND. Cross-shore reach of what 1996 writes, metres
    from the ocean edge, median over profiles, against the island's own extent:

        domain      island      road 84 / 04      1996 reach, no boundary
            11        1999          416 / 499                          560
            60        1274          621 / 621                          979
            77        1627          717 / 717                          429
            85        1782          185 / 305                          432
            90        1999          382 / 382                          441

    ALACE stops itself at 429-979 m against an island 1274-1999 m wide. The
    road was never what held it back. Domain 77 is the proof: the reach is
    429 m under every option while the road sits at 717 m, so there the
    boundary was not even binding.

  * IT ADMITS NO EXTRA JUNK. Island-wide maximum stays 12.00 m under all three
    options - the ceiling, not the mask, is what holds the return tail.

And the window origin is untouched by the choice. start_beach moves in ZERO
domains between the 1984 and 2004 masks (max |delta| 0.00 m) and in 6 domains
between 1984 and none, max 26 m. All three agree at the seaward edge; every
difference is landward, behind the beach start.

WHAT WAS GIVEN UP. A stated rule was traded for a survey footprint. "Ocean
side of NC-12" is defensible in a methods section; "wherever ALACE flew" is a
survey artefact standing in for a landform boundary. The honest description of
this product is that its 1996 extent is the ALACE swath, and the audit reports
n1996_landward_of_1984_road / _2004_road per domain so the discarded band is a
number rather than a claim.

Both alignments are still read and still rasterized every run. They gate
nothing. They are reported.

WHAT THE 1996 SOURCE ACTUALLY IS
1996 Fall East Coast NOAA/NASA Airborne LiDAR Assessment of Coastal Erosion
(ALACE), tiles J1441002 R0_C0 and R1_C0. 3 m grid, EPSG:6347 (NAD83(2011) /
UTM 18N), nodata -999999, +/-15 cm stated vertical accuracy.

Its abstract is the whole methodological problem in one sentence: the aircraft
surveyed "from the low water line to the landward base of the sand dunes."

Measured consequences, over the ocean-side band in 13 sampled domains:

    coverage of the ocean-side band   30-84%, median ~53%
    coverage LANDWARD of the road     0-16%   - it stops before NC-12
                                      ^ this is now ADMITTED, not discarded,
                                        and is 17.2% of what 1996 writes
    seaward reach vs 2009 waterline   +70 to +133 m
    median(1996 - 2009) on overlap    +0.11 to +0.36 m in the developed reach

So this is NOT "the ocean side becomes 1996". It is a 1996 beach and foredune
grafted onto a 2009 backdune and interior, and the seam lands at the ALACE
swath's landward edge - near the dune toe - not at the road. Which is the same
place HAT_dune_topo_extractor.py picks its dune windows. Read a dune-crest
number near that line knowing it sits beside a survey boundary.

VERTICAL DATUM IS AN INFERENCE, NOT A DECLARATION
The .prj carries the horizontal CRS only. The sole evidence for the vertical is
`geoid18` in the NOAA Digital Coast S3 path the metadata cites, i.e. NAVD88 via
GEOID18, metres. The empirical check is what actually supports it: median
(1996 - 2009) over co-measured cells is +0.11 to +0.36 m in the developed
reach. A NGVD29 or ellipsoid-height mistake would show METRES of disagreement,
not decimetres, and the sign is the one 12 years of erosion predicts. Both
per-domain bias estimates are written to the audit every run so this stays
checkable rather than asserted.

Domains 11 and 86 report +1.64 and +1.66 m. Both have 4% and 22% 2009 coverage
in the band, so the overlap is a sliver at the waterline. Do not read those two
as a datum finding.

THE SPLIT FLOOR, AND WHY IT IS NOT THE SAME NUMBER TWICE
The first build used the existing -2.64 m NAVD88 floor for both cases. It moved
33 of 83 domains' beach start LANDWARD, by up to 27 m, while 1996 measured
HIGHER than 2009 in the overlap - the opposite of what erosion predicts.

The cause is the ALACE swath running to the low water line. It carries wet
swash and beach-face returns between -2.64 m and MHW. Ocean-side of the road
those cells overwrote 2009's dry beach and pushed the first cell above
BEACH_START_THR_M (0.50 m MHW) landward. Measured on the worst offenders,
median shift of start_beach in metres, positive = seaward:

    domain   floor -2.64, overwrite   floor MHW, overwrite   1996 into gaps only
      73             -27.5                    0.0                    0.0
      42             -23.0                    0.0                    0.0
      58             -22.0                    0.0                    0.5
      90             -16.0                   +1.0                   +1.0
      79             +52.5                  +52.5                  +51.5
      80             +71.0                  +71.0                  +71.0
      77             +59.0                  +59.0                  +59.0
      11             +57.0                  +57.0                  +57.0

Two things fall out of that table and both matter:

  * Raising the floor to MHW *for the overwrite case only* takes every landward
    mover to zero and leaves every seaward mover untouched.
  * The seaward signal is identical when 1996 is only allowed into nodata. The
    real gain is where the 2009 survey saw NOTHING, not where the two surveys
    disagree. The overwrite is still kept, because it is what carries the 1996
    beach and foredune ELEVATIONS - it just does not move the window.

Hence two floors, for two different questions:

    NO other survey    ->  admit 1996 at -2.64 m, the extractor's WATER_CLAMP_M,
    saw the cell           exactly as 2014 is admitted. Nothing is lost, and the
                           sub-MHW seaward gain is kept.
    2009 OR 2014 saw   ->  admit 1996 only above MHW. One survey may replace
    the cell               another's measurement, but not with a wet return.

    The first build tested only `~valid09` for the first case and shipped 16
    landward movers, D73 at -21 m and D90 at -16 m. The 2009 holes are not
    empty in the product - 2014 fills them - so that let a -2 m 1996 swash
    return displace a +1 m 2014 measurement. The guard at the end of main()
    is what caught it, and it stays there for the next time.

No value is modified by either. Only the selection differs, and every cell each
floor rejects is counted per domain in the audit.

Resampling 3 m -> 1 m was tested and is irrelevant here: nearest and bilinear
agree to within 1 m of start_beach in every domain checked. Nearest is used, as
in HAT_dem_gap_fill.FillSource, because it invents no values at survey edges.

DOMAINS 1-7 - RESOLVED 2026-08-26
Both NC-12 exports end at domain 8. Domains 1-7 sit 265 m to 3.3 km beyond the
terminus, so there was no line to define an ocean side for them, and under the
old boundary they received no 1996 at all. They carried a 2009 beach while
8-90 carried a 1996 one, leaving an alongshore discontinuity in the initial
condition at the 7/8 seam. The documented alternative was to invent a boundary
by extending domain 8's alignment south across Hatteras village.

Removing the boundary dissolves the problem instead of trading it. 1-7 now
receive 1996 on exactly the same terms as every other domain - 77k to 175k
cells each, 6k to 60k of it ground no other survey saw - and nothing is
invented. The seam is gone.

They are still flagged `road_line=False` in the audit. That flag is now purely
descriptive: it says these domains have no NC-12 export, not that they were
treated differently.

WHAT MOVES DOWNSTREAM, AND WHAT DOES NOT
MOVES.  Each domain's cross-shore window origin. HAT_dune_topo_extractor.py
        sets start_beach to the first cell above 0.50 m MHW, so adding 1996
        land above that threshold slides the window seaward. Measured with the
        split floor: 35 of 83 domains move >= 1 model cell seaward, 6 are
        unchanged, none move landward. The developed reach 76-81 moves 47-72 m,
        i.e. 5-7 Barrier3D cells. Every domain's before/after is in the audit as
        start_beach_before_m / start_beach_after_m / start_beach_shift_m, so
        this is verifiable per domain rather than asserted island-wide.
        A NEW PICK SET IS THEREFORE REQUIRED - the 2009_v5 dune windows were
        picked against a different origin.

DOES NOT.  `shoreline_offset`. That comes from
        2-brie-offset/1984/Island_Dune_Offsets_1984_CASCADE_Input.csv,
        is measured independently of any DEM, already exists for 1984, and is
        passed straight to Cascade(). Nothing here touches it. It was a
        deliberate choice to leave the two independent rather than re-derive one
        from the other.

INPUTS
    D:/Hatteras_GIS/.../2009_full.tif                base, 1 m, EPSG:3725
    D:/Hatteras_GIS/.../2014_full.tif                gap fill, 1 m, EPSG:6347
    D:/Hatteras_GIS/elevation/dem/1996/noaa_1996_FallEC_J1441002/*.tif   override, 3 m, EPSG:6347
    D:/Hatteras_GIS/domains.geojson                  90 boxes, 2000 x 500 m
    4-mgmt-forcing/.../1978/nc12_1978.geojson        EPSG:2264, US survey FEET

OUTPUTS (data/hatteras_init/0-elevation/2009-2014-1996/1-gapfill-1m/)
    clip_domain_<N>_filled.tif   the mosaic, m NAVD88
    clip_domain_<N>_survey.tif   provenance per cell:
                                     1996 = ALACE, wherever it has data
                                     2009 = measured by the 2009 DEM
                                     2014 = filled from 2014 NOAA Post-Sandy
                                        0 = no survey saw it
    mosaic_1984_audit.csv        per-domain counts for every rule above

Filenames match HAT_dem_gap_fill.py's so HAT_dem_resample_clip.py needs only
its SOURCE_TAG changed. The survey raster now carries FOUR codes rather than
three, which its downsample_survey() has been extended to handle.

Requires: rasterio, geopandas, numpy, scipy

    python HAT_dem_1984_mosaic.py
```

Notes that were in the code:

```text
Shared helpers come from the 2009 script rather than being copied. Same
directory, and importing it runs no IO - its module level is paths and
constants only, and main() is guarded.
```

```text
Named for its COMPOSITION rather than the start period - see the note on
PRODUCT_TAG in HAT_dem_gap_fill.py. Resolved through
scripts/site_layer/hat_elevation_products.py so the layout lives in one place.
```

```text
--- THE ROAD LINES ARE DIAGNOSTIC ONLY (2026-08-26) ------------------------

They no longer gate anything. 1996 is admitted wherever it has data that
clears the floors, the ceiling and connectivity, and the landward limit is
the ALACE swath's own edge. See "WHY THERE IS NO ROAD BOUNDARY" above.

Both vintages are still read and still reported per domain, because "how far
landward of NC-12 did 1996 actually write?" is the first question a reader
will put to this product, and the audit should answer it rather than leave
it to be re-derived. NO_ROAD_POLICY is gone with the boundary it configured.
Keyed by PERIOD; each period's LINE comes from ROAD_LINE_FOR_YEAR.
```

```text
--- THE SPLIT FLOOR. See the module docstring; these are not the same number
twice by accident, and collapsing them re-introduces the 33 landward
movers. ---
```

```text
--- THE CEILING. 2009 and 2014 are clean gridded DEMs and need none; this is
specific to ALACE and it is not optional. ---

The 1996 grid runs from -16.56 to +256.65 m NAVD88. Its 99th percentile is
7.61 m, so everything above roughly 12 m is the uncorrected-return tail that
ALACE-era ATM data is known for - cloud, bird and aircraft returns the
vendor's "visual inspection" filter did not catch.

Shipped without a ceiling, the first build put 250 m spikes into 66 of 90
domains, against a 2009+2014 island-wide max of 10.18 m. The floors only ever
guarded the low side; nothing guarded this one.

12 m, from the distribution of the 12,021,270 cells 1996 actually wrote,
1 m bins:

5 m 459288   6 m 222927   7 m 86532   8 m 34755   9 m 18645
10 m  15561  11 m  13362  12 m 10458  13 m  4881  14 m  2328
15 m   1923  16 m   1752  17 m   885  18 m  1224  19 m   675
20 m    594  ... a flat few-hundred-per-bin tail to 58 m, then to 256 m

The terrain population decays steeply to 12 m and HALVES at 13 m; past ~15 m
the histogram is flat, which is an artifact population, not a landform. 12 m
also sits 18% above the 10.18 m island-wide max of the 2009+2014 product -
enough headroom for a 1996 dune genuinely taller than anything the 2009
survey measured, which is the erosion signal this product exists to capture,
without admitting the tail.

Rejects 37,335 cells, 0.31% of what 1996 wrote. A rejected cell is NOT lost:
it falls through to 2009, then to 2014, so the ceiling costs coverage only
where no other survey saw the cell at all. Counted per domain in the audit as
dropped_1996_ceiling.

NOTHING IS CLAMPED. A cell over the ceiling is rejected, not pulled down to
it - the same posture as every other rule here.
```

```text
Reproduces HAT_dune_topo_extractor.py so the audit's shift column means the
same thing the extractor will do. Changing either without the other makes the
diagnostic quietly wrong.
```

```text
Built for the AUDIT ONLY - neither is applied. They answer "how far
landward of each alignment did 1996 write?", which is what a reader
needs in order to judge this product now that the survey's own swath
edge is the only landward limit.
```

```text
--- STAGE 1: the 1996 override, ocean side only -------------------
"Gap" means NO OTHER SURVEY SAW THIS CELL - 2009 or 2014. Not just
2009. The first build tested `~valid09` and the guard below caught
it: 16 domains still moved landward, D73 by 21 m and D90 by 16 m.
The 2009 holes are not empty in the product, 2014 fills them, so
`~valid09` let 1996 put a -2 m swash return over a +1 m 2014
measurement - precisely the replacement the MHW floor exists to
stop, just wearing a different survey's name. The stated principle
was always "one survey may replace another's measurement, but not
with a wet return"; ANOTHER means any, so the test is coverage by
any other survey.

has14 is used raw, before 2014's own floor and connectivity. That is
deliberate and conservative in the right direction: where 2014 has
any value at all, 1996 must clear MHW to displace it.
```

```text
NO ROAD BOUNDARY - the one substantive change of 2026-08-26. `ocean`
is built above and deliberately NOT applied here. The landward limit
is the ALACE swath edge, measured at 429-979 m from the ocean edge
against an island extent of 1274-1999 m: the survey stops far short
of the sound unaided, so the road was never what held it back. In
domain 77 the 1996 reach is 429 m while the road sits at 717 m - the
boundary was not even binding there. Every rule below is unchanged.
```

```text
Reported every run whether or not anything is applied, so "we chose
not to correct" stays a checkable claim rather than a comment.
base - fill, the sign convention of gf.estimate_bias. NEGATIVE
means 1996 sits ABOVE 2009, which is what erosion predicts and
what the run reports (about -0.25 m through the developed reach).
```

```text
The ocean-side-of-1984 figure the docstring's datum argument was
built on. Kept so that argument stays checkable on its ORIGINAL
footprint now that the written population is a fifth larger, rather
than being quietly restated against a different set of cells.
```

```text
--- STAGE 2: the 2014 gap fill ------------------------------------
Run against the POST-1996 surface: 1996 cells are measured ground and
are legitimate connectivity anchors. Landward of the road this is
identical to HAT_dem_gap_fill.py, because nothing changed there.
```

```text
--- the before/after the window origin depends on ------------------
"Before" is rebuilt here rather than read from the 2014 product, so
the two sides of the comparison come from one code path and a stale
product on disk cannot silently change the answer.
```

```text
DIAGNOSTIC ONLY from 2026-08-26 - neither alignment gates the
override any more. road_line is retained because readers of this
audit key off it, and because "does this domain have a road line
at all" remains a real property of the domain.
```

```text
EVERY domain now receives 1996, so the shift statistics run over all of
them. Before 2026-08-26 they ran over `with_road` only, because domains
with no road line got no override at all and their shift was
structurally zero - averaging those in would have diluted the number.
That no longer applies. `with_road` survives to report the split.
```

```text
The threshold that matters is one Barrier3D cell, not one metre. This
diagnostic runs at 1 m; the extractor works at 10 m, so a shift under
10 m cannot move the window by more than a single cell and usually moves
it by none. Both bands are printed - the sub-cell one is residual, the
whole-cell one is the wet-edge failure the split floor exists to prevent
and means the rule needs re-checking before this product is used.
```

<details><summary>Function notes (the original docstrings)</summary>

**`TiledSource()`**

```text
Several single-band tiles read as one raster on the base grid.

ALACE ships J1441002 as two row-tiles that abut at y=3947502 rather than as
one mosaic. They are read independently and pasted in order; a later tile
never overwrites a cell an earlier one already filled, so the shared edge
resolves deterministically instead of by read order.

Each tile goes through gf.FillSource, so the CRS difference (EPSG:6347 ->
EPSG:3725) and the 3 m -> 1 m resolution difference are handled exactly as
they are for the 2014 fill, by the same code, with the same nearest-
neighbour choice.
```

**`ocean_side_mask()`**

```text
True ocean-side (east) of the road, per raster row.

OCEAN_LOC is "right" throughout this pipeline - the domains are 2000 m
cross-shore in x by 500 m alongshore in y, and the Atlantic is at increasing
easting. So the boundary is one column index per row.

The MOST SEAWARD road cell in a row is used, not the first. Where the
alignment runs diagonally through a row it occupies several columns, and
taking the max yields the SMALLEST ocean region - the choice that cannot
accidentally place a road cell on the ocean side of itself.

Rows the road does not reach are interpolated from the rows it does, then
held flat past the ends. Without that, the 200 m context pad and any row
where the line steps outside the window would punch holes in the mask, and
the override would be ragged along the domain edges for no physical reason.

Returns (mask, n_rows_with_road). An all-False mask with 0 rows means this
window has no road at all - the caller decides what that means.
```

**`start_beach_median_m()`**

```text
Median over the alongshore profiles of the extractor's start_beach, in
metres from the ocean edge of the window. Smaller = further seaward.

Mirrors HAT_dune_topo_extractor.load_domain exactly:
    ocean-first (arr[:, ::-1], since OCEAN_LOC="right")
    z = raw - MHW ; z below WATER_CLAMP_M pinned to the sentinel
    first index where z > BEACH_START_THR_M
then default_window()'s median over profiles. Nodata is treated as the water
sentinel here, which is what the clamp does to it downstream.

This is computed at 1 m. The extractor works at 10 m, so divide by 10 for
Barrier3D cells - the audit carries both.
```

**`select()`**

```text
Connectivity, shared by both stages. `cand` is already coverage- and
floor-filtered; this only drops what the island cannot reach.
```

</details>

### 2-produce/HAT_dem_duneline_coverage.py - deleted 2026-10-19

Retired with the 1984 start: it measured the 1984 dune line against the 1996-start DEM. Its outputs (`0-elevation/2009-2014-1996-duneline/`) are on `D:\archive\CASCADE_offload\data\hatteras_init\0-elevation\`. Recover the script with `git show 2046d67b:scripts/input_prep/0-elevation/2-produce/HAT_dem_duneline_coverage.py`.

### 2-produce/HAT_dem_gap_fill.py

Step 1 of 3 (2009-start DEM): clip the 2009 DEM per domain and fill its gaps from the 2014 NOAA Post-Sandy DEM.

From the script's original header:

```text
Step 1 of 3: clips the 2009 DEM to each domain and fills its gaps from the 2014
NOAA Post-Sandy DEM. Writes one gap-filled 1 m clip and one survey-year raster
per domain; HAT_dem_resample_clip.py resamples them to 10 m.

WHAT THIS IS ACTUALLY FIXING
Not "voids in a surface". Measured on the real DEM, each domain contains exactly
TWO nodata regions and both touch the domain edge - interior enclosed nodata is
0 cells, 0.00%. With OCEAN_LOC="right" in HAT_dune_topo_extractor.py, the east
region (~480 m) is the Atlantic and the west region (~1045 m) is Pamlico Sound.
The 2009 survey simply stops at the waterline on each side.

The gap that matters is the sound-side margin, and HAT_dune_topo_extractor.py
already documents why (lines 281-297): roadway_manager.bulldoze drowns a roadway
when >20% of the cells BORDERING it sit at or below 0 m MHW, and a no-data cell
passes that test. In GIS 78/79/80 the row landward of NC-12 is 17-25 no-data
cells and ZERO genuinely wet ones, so all three roadways width-drowned at t=0 on
missing coverage alone - while the profiles were still 0.5-0.7 m ABOVE MHW.

So this fills measured ground the 2009 survey missed. It does not invent
elevation anywhere.

THE FOUR RULES THAT BOUND THE FILL
1. COVERAGE.   Only cells the 2014 DEM actually has a value for. It covers
               97.34% of the island's DRY-LAND gaps - cells 2009 missed where
               the consensus of candidate DEMs puts the ground above MHW. See
               FILL_DEM_PATH for the full scoring and why the alternatives lost.
2. CONNECTIVITY. Only cells contiguous with the island's valid 2009 surface, so
               detached marsh hummocks and any water returns SMRF kept out in
               the sound cannot become new land in the barrier interior. This
               cannot exclude the target cells: the fringe landward of NC-12 is
               contiguous with the island by definition.
               Computed on a BUFFERED window - on the bare 500 m strip, marsh
               that connects to the island just outside the domain would be
               severed by the crop.
3. ELEVATION.  A guard at the DOWNSTREAM water threshold (-3.0 m MHW, the
               extractor's WATER_CLAMP_M), not at MHW. Flooring at MHW would
               discard real low marsh and would not have protected the drowning
               fix anyway - see the note on FILL_MIN_ELEV_NAVD. Sub-MHW fills
               are counted, not rejected.
4. VERTICAL.   Nothing is applied. Both bias correction and feathering are OFF,
               so a filled cell is the 2014 measurement unchanged. Both bias
               estimates are still computed and written to the audit every run.

Everything each rule rejects is counted per domain in the audit CSV. Nothing is
dropped silently.

THE SEAM IS REAL, KNOWN, AND DELIBERATELY LEFT IN
Where the fill meets measured 2009 ground there is a step. The numbers below
were measured on the SUPERSEDED 2008 point-cloud attempt; seam_check() re-runs
every time, so the current source's figures are in the audit CSV
(seam_median_abs_m against ctrl_median_abs_m). The reasoning for leaving it
uncorrected carries over:

  measured <-> measured, whole domain      0.028 m   terrain roughness
  measured <-> measured, near the boundary 0.028 m   same - the margin is smooth
  fill <-> measured (what we ship)         0.341 m   ~12x either control
  2008 <-> 2008 across the same boundary   0.030 m   <- the decisive one

The last row grids 2008 on BOTH sides of the boundary, so any inter-survey
offset cancels. It comes out flat, at the roughness floor. The ground therefore
does NOT drop at the marsh edge: terrain accounts for ~9% of the step and the
rest is the two surveys disagreeing where they meet.

The obvious fix - shift 2008 onto 2009 per domain - was rejected because the
disagreement is not a datum offset and does not have a consistent sign. Signed
step by domain: -0.685, +0.136, -0.787, -0.117, +0.277, +0.455, +0.109. The fill
is too low in some domains and too high in others, and the boundary estimate
(+0.25 median) disagrees with the whole-domain overlap (-0.05 median). Any
single-number shift removes the seam at the boundary and introduces a comparable
disagreement across the interior instead - trading a visible artifact for an
invisible one.

Feathering was rejected for the same reason: it hides the step over 5 m without
addressing the disagreement, and smooths measured data to do it.

So the step stays, and clip_domain_<N>_survey.tif marks exactly which cells came
from 2008 so any consumer can find it. At the 10 m Barrier3D grid it collapses
to one cell boundary - 0.3 m over 10 m, ~3% slope, within the range of real
back-barrier relief, but it is fabricated and worth knowing about when reading
overwash behaviour near a fill margin.

INPUTS
    D:/Hatteras_GIS/.../2009_full.tif   base, 1 m, EPSG:3725 + NAVD88
    D:/Hatteras_GIS/.../2014_full.tif   fill, 1 m, EPSG:6347 + NAVD88
                                        (reprojected on read by WarpedVRT)
    D:/Hatteras_GIS/domains.geojson     90 boxes, 2000 x 500 m

OUTPUTS (data/hatteras_init/0-elevation/1-gapfill-1m/)
    clip_domain_<N>_filled.tif      gap-filled clip, m NAVD88
    clip_domain_<N>_survey.tif      which survey each cell came from:
                                        2009 = measured by the 2009 DEM
                                        2014 = filled from the 2014 DEM
                                           0 = neither survey saw it
    gapfill_audit.csv               per-domain counts for every rule above

Filenames follow the legacy 2009-domain-clipresample convention
(clip_domain_N.tif) with _filled marking the new set. Layout is flat rather than
per-domain subfolders: step 1 writes the 1 m clips and step 2 writes the 10 m
ones, so one folder per step means re-running a step is a single delete.

Requires: rasterio, geopandas, numpy, scipy
```

Notes that were in the code:

```text
FILL SOURCE: 2014 NOAA Post-Sandy DEM.

Chosen by measurement, not vintage. Every candidate DEM was scored against the
2009 DRY-LAND gaps (cells 2009 missed where the consensus of candidates puts
the ground above MHW), over all 90 domains:

2014 NCFMP          100.00%   DISQUALIFIED - hydro-flattened, 70.5% of its
values in the gap are the constant -0.762 m
2014 NOAA Post-Sandy 97.34%   <- earliest genuine, and the best
2017 USACE           23.74%
2016 post-Matthew    21.98%
2019 DUNEX           21.98%
2018 post-Florence    9.74%

The 2016/2017/2019 collapse is spatial extent: they score 95-100% on domains
78-80 but are localised surveys - 2019 is below 50% in 60 of 90 domains. A
choice made on the developed reaches alone would have picked a dataset
covering less than a quarter of the island's gaps.

Different CRS from the base (EPSG:6347 NAD83(2011) vs EPSG:3725
NAD83(NSRS2007)); WarpedVRT reprojects on read.

Superseded point-cloud attempts, kept for the record:
2008 NOAA IOCM  topo-only, ~19% of the sound-side gap, 17/150 in the NC-12 strip
2011 post-Irene topo-only, ~35%, 28/150
Neither is reproducible from this script any more - the point-cloud path was
removed on 2026-08-26 along with the classifier that fed it. A point-cloud
candidate now has to be gridded to a DEM before it gets here.
```

```text
The product this script builds. Named for its COMPOSITION, not for the fill
source and not for a hindcast period - this DEM currently serves both the
1984 and the 2004 period, so a period in the name would be a false claim.
FILL_SOURCE_TAG above still names the SOURCE, and appears in the console
output and the figure labels.
```

```text
Paths come from scripts/site_layer/hat_elevation_products.py, not from string
concatenation here. Six scripts used to build them by hand and that is how
HAT_road_elevation.py silently stopped finding its rasters - see the note at
the top of that module.
```

```text
--- RULE 3: elevation floor ---
Set to the DOWNSTREAM threshold, not to MHW, and it is a guard rather than a
filter. Reasoning, because this was reversed once already:

* HAT_dune_topo_extractor.py does the MHW referencing itself (line 1017,
z = raw - 0.36) and applies its own water threshold at WATER_CLAMP_M =
-3.0 m MHW = -2.64 m NAVD88. That -3.0 was picked deliberately -
"keeps back-barrier marsh cells (Lexi's v3 edit)", up from -1.0. Flooring
at MHW here would undo that decision one step upstream, invisibly.
* A floor at MHW never protected the road-drowning fix it was added for.
roadway_manager.bulldoze drowns when >20% of bordering cells sit at or
below 0 m MHW; a filled cell at -0.2 m MHW counts as wet exactly as the
-3.0 sentinel did when it was unsurveyed. The bug was cells genuinely
ABOVE MHW reading as wet because nobody surveyed them, and filling those
with their true elevation fixes it whatever the floor is.
* The 2008 cloud bottoms out at -1.33 m NAVD88, so this floor rejects
nothing in practice. It exists to catch a future fill source that could.

Sub-MHW fills are still COUNTED per domain (cand_below_mhw in the audit), so
lowering the floor costs no visibility.
```

```text
Strict connectivity severs a 300 m marsh platform if a 20 m tidal creek that
is water in BOTH surveys separates it from the island. Gaps up to this width
are bridged before the connectivity test, so creek-separated back barrier is
kept while genuinely detached patches out in the sound are still rejected.
0.0 disables bridging (strict connectivity).

20 m, measured rather than guessed. Across 12 sampled domains there were 2969
detached components holding 1,055,741 candidate cells. Recovery vs bridging
distance, by component count and by CELL count (cells are what matters - the
goal is captured area, and one 200k-cell platform outweighs 50 specks):

bridge   % components   % detached cells
2 m           7.1%              3.5%
5 m          10.3%             16.2%   <- first step
15 m          14.6%             18.6%
20 m          16.2%             30.4%   <- second step, then a plateau
30 m          18.9%             31.6%
100 m          31.6%             39.9%
200 m          47.7%             47.6%   <- columns converge

20 m sits on the plateau right after the second step: 20->30 m buys 1.2 more
points, 30->100 m buys 8.3 for five times the reach. Below 200 m the cell
column runs at ~2x the component column, i.e. bridging is selectively catching
large platforms; by 200 m they converge, which means it has stopped gaining
area preferentially and is just admitting open sound. 20 m is also a credible
tidal-creek width, which 100 m is not.

In context of ALL candidates in those domains (4,420,718 cells):
strict 76.1%  |  5 m 80.0%  |  20 m 83.4%  |  100 m 85.7%

Note this barely moves the domains that motivated the fill. In GIS 78/79/80
only ~20k of ~300k candidates are detached at all, so ~93-94% is already
connected and their detached patches sit at 49-174 m - open water, not creeks.
Bridging is about how much back barrier domains 10/50/90 carry, not the road.
```

```text
--- RULE 4: vertical reconciliation ---

BOTH VALUE-MODIFYING STEPS ARE OFF. Of the rules here, only these two change
what a cell says; coverage, connectivity and the floor merely select which
measured cells get used. With both off, a filled cell is the 2008 measurement
and nothing else.

Bias correction was ON and was misfiring. It estimated a single offset from a
10 m collar around the fill, which by construction sits on the 2009 waterline
- the one place the two surveys disagree for reasons that are not a datum
offset. min-z gridding there picks water-surface returns, so the collar
measured "min-z reads low in wet cells" and applied it hundreds of metres into
dry marsh. Across 28 domains it produced +0.824 to -0.104 m (median +0.125),
the signature of an unstable estimator rather than real per-domain offsets.
Full-domain overlap says the surveys actually agree to ~0.03 m.

Both estimates are still COMPUTED and written to the audit every run, so the
decision to not apply them stays visible and checkable.
```

```text
How the 2008 ground returns become a 1 m raster.
"median"  median of ground returns in the cell. Closer to how a gridded DEM
surface is built, so it is comparable to the 2009 values it sits
beside, and one stray low return cannot set the cell.
"min"     lowest return. The classic bare-earth proxy, but over wet marsh
the lowest return is often the water surface - and cells below
0 m MHW are exactly what the drowning test counts, so this is
conservative in the wrong direction for this particular bug.
```

```text
The _survey raster stores the year each cell's elevation came from, so it
needs no legend. 0 means neither survey saw the cell.
```

```text
Every key present every time, so the audit CSV has stable columns even
for a domain with no fill at all.
```

```text
Both estimates are computed and reported every run even when nothing
is applied, so "we chose not to correct" stays a checkable claim.
```

```text
Counted on the CROPPED domain window, not the padded context window,
so it is comparable with `filled` in the audit. Counting it on the
padded window made this column exceed `filled` and read as >100%.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data/hatteras_init.

NOT parents[N]. This file moved into 2-produce/ on 2026-08-25, and the
old parents[3] then resolved to input_prep/ rather than the project root.
That raises nothing - it just makes every path below it wrong, silently,
until some glob comes back empty. Same helper and same reason as
4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
```

**`snap_window()`**

```text
Polygon bounds -> integer window on the source grid, trimmed to whole
`block`-sized blocks so the later 10 m resample divides evenly.

Nearest cell edge, not floor/ceil: a polygon 0.4 m off grid should snap to
the near edge, not grow the window by a whole cell. Returns the adjustment
so it is reported rather than absorbed.
```

**`read_window()`**

```text
Boundless so a domain hanging off the DEM edge yields nodata, not an
error - the domains do overshoot the DEM's east edge by ~6 m.
```

**`FillSource()`**

```text
The candidate DEM, read onto the base DEM's grid on demand.

WarpedVRT handles the CRS difference (EPSG:6347 -> EPSG:3725) and any
resolution difference on read, so nothing downstream needs to know the
source is in a different realisation of NAD83.

Nearest-neighbour resampling deliberately: this fills cells the 2009 survey
missed, and interpolating would invent values at the very edges where the
two surveys meet, which is where they are least comparable.
```

**`island_connected()`**

```text
Candidate cells reachable from the island through valid ground or other
candidates. The island is the largest connected component of valid 2009
cells in the (buffered) window.

bridge_px dilates the land mask before the reachability test, so a gap up
to 2*bridge_px wide is crossed. Dilation is used ONLY to decide
reachability - the returned mask is still a subset of the real candidates,
so no cell is invented by bridging.
```

**`seam_check()`**

```text
How big a step does the fill create where it meets measured 2009 ground?

Reported against a CONTROL: the same statistic between adjacent measured
cells. Real terrain is not flat, so a seam step only means something
relative to how much neighbouring cells normally differ. seam ~ control
means the fill is indistinguishable from the surface it joins, and no
feathering is warranted. seam >> control is a genuine cliff.
```

**`boundary_extrapolation()`**

```text
Locally-consistent continuation of the 2009 surface into the fill area,
used only as the blend target at the seam so the merge has no hard step.
```

**`resolve_crs()`**

```text
The DEM is a COMPOUND CRS (EPSG:3725 + NAVD88), whose to_epsg() is None,
so a plain `gdf.crs != src.crs` reports a reprojection that is a no-op.
```

**`read_on_grid()`**

```text
WarpedVRT rejects boundless reads, so partial overlap is intersected
and pasted at the right offset rather than erroring.
```

</details>

### 2-produce/HAT_dem_resample_clip.py

Step 2 of 3: resample each 1 m domain clip to the 10 m Barrier3D grid (50 x 200 per domain).

From the script's original header:

```text
Step 2 of 3: resamples each gap-filled 1 m domain clip (output of
HAT_dem_gap_fill.py) to the 10 m Barrier3D grid, 50 x 200 per domain.
HAT_export_to_numpy.py converts these to .npy in the final step.

CLIPPING NOW HAPPENS IN STEP 1 - AND THE ORDER MATTERS
Barrier3D needs each domain to be an exact 50 x 200 array of true 10 x 10 m
cells. That only works if the 10 m grid is built from the domain's own corner,
so each 10 m cell is an exact 10 x 10 block of 1 m source cells.

Resampling the whole island first and clipping second cannot deliver that. The
domain boxes do not sit on a global 10 m grid - their origins land at arbitrary
sub-10 m offsets (450439.120, 454507.796, ...). Cutting them from a grid whose
cell edges fall on multiples of 10 snaps outward to the enclosing cells, giving
51 x 201 for most domains and shifting every domain up to half a cell off its
polygon.

Confirmed against the existing ArcGIS outputs: clip_domain_N.tif and
resampled_domain_N.tif share a byte-identical origin for every domain checked,
so the 10 m grid was built inside the clip. There is also no global grid worth
preserving - domains 50 and 51 share an x origin but their y origins differ by
505.126 m against a 500 m extent, so the boxes were never a tiling of one grid.

Clipping therefore lives in step 1, which needs the 1 m window anyway to do the
fill. This step only reduces 1 m -> 10 m, which keeps the window definition in
exactly one place.

RESAMPLE METHOD - SETTLED EMPIRICALLY
Reconstructed resampled_domain_N.tif from clip_domain_N.tif for domains 1, 2, 3,
25, 50, 75, 100 and 120. ArcGIS used BILINEAR, not the Nearest Neighbor default
the old docstring assumed:

    method                          agreement with existing files
    nearest (any of the 100 cells)  0.1% of cells (chance)
    block mean (all 100 cells)      max diff 1.16 m
    bilinear                        exact, max diff 2.4e-07 (float32 epsilon)

At an exact 10x reduction the 10 m cell center lands on the boundary between
source cells 4 and 5 on both axes, so bilinear collapses to a tie: the plain
mean of the central 2 x 2 source cells. It uses 4 of the 100 cells under each
output cell and ignores the other 96. That reproduces your existing domains, so
it is the default - see AGGREGATION for the alternative.

INPUTS  (data/hatteras_init/0-elevation/1-gapfill-1m/)
    clip_domain_<N>_filled.tif
    clip_domain_<N>_survey.tif

OUTPUTS (data/hatteras_init/0-elevation/2-resampled-10m/)
    resampled_domain_<N>_filled.tif     the 50 x 200 Barrier3D grid
    resampled_domain_<N>_survey.tif     2009 / fill year / 0 per cell
    resample_audit.csv

Requires: rasterio, numpy
```

Notes that were in the code:

```text
Must match FILL_SOURCE_TAG in HAT_dem_gap_fill.py, or PRODUCT_TAG in
HAT_dem_1984_mosaic.py - each source keeps its own subfolder so a re-run
cannot clobber another source's rasters or its audit CSV.

Selected on the command line rather than by editing, because there are now
two live products and hand-editing a constant per run is how the figures
ended up labelled with the wrong source once already:

python HAT_dem_resample_clip.py                              # baseline
python HAT_dem_resample_clip.py --product 2009-2014-1996     # 1984 start
```

```text
AGGREGATION
"arcgis_bilinear"  mean of the central 2 x 2 source cells - reproduces your
existing domains exactly. Uses 4 of 100 cells.
"mean"             mean of all 100 cells. More defensible for elevation, but
will NOT reproduce the existing files (up to ~1.2 m at
dune crests, where sampling 4 cells is least
representative).
"nearest"          single source cell. NOT what produced the existing files.
```

```text
ArcGIS emitted a partial-weight value at some nodata edges and nodata at
others. Strict (all 4 required) matched it exactly in the interior and
differed only at <= 13 cells per domain, all on nodata margins. Strict does
not invent elevation at the water edge, so it is the default.
```

```text
Every non-base code a survey raster may carry, MOST SPECIFIC FIRST. This is
the precedence downsample_survey resolves a mixed 2 x 2 block with, so the
order is a decision, not a list:

1996 first  it is the override that defines the 1984-start product. A 10 m
cell that drew any of its four read cells from ALACE should say
so - that is the flag a reader uses to find the graft.
2014 next   the gap fill, and the only code the 2009-start product has.

SURVEY_NONE outranks both: an unsurveyed cell in the core is also when the
elevation output is nodata under BILINEAR_REQUIRE_ALL_FOUR, so the two agree
by construction. Blocks that mix two fill codes are counted per domain in the
audit as `mixed_source_cells`, so the precedence never hides how often it had
to choose.
```

```text
Kept because the audit and the console line report "filled cells" against one
code. For the 1984 product that is the 1996 count; the 2014 count is reported
beside it.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data/hatteras_init.

NOT parents[N]. This file moved into 2-produce/ on 2026-08-25, and the
old parents[3] then resolved to input_prep/ rather than the project root.
That raises nothing - it just makes every path below it wrong, silently,
until some glob comes back empty. Same helper and same reason as
4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
```

**`downsample()`**

```text
Reduces a (block*R, block*C) array to (R, C). Every output cell is an exact
block x block window of source cells - that is what makes the 10 m cells
true 10 x 10 m cells rather than resampled approximations.
```

**`downsample_survey()`**

```text
Survey year for the SAME four cells bilinear actually reads, so the flag
describes the value that was written rather than the whole block. A 10 m
cell reads a fill code if any of the central 2 x 2 came from that fill, and
0 if any of them was unsurveyed - which is also when the elevation output is
nodata under BILINEAR_REQUIRE_ALL_FOUR, so the two agree by construction.

Returns (survey_10m, n_mixed) where n_mixed counts blocks whose four read
cells carried more than one fill code and SURVEY_FILL_CODES had to break
the tie.
```

</details>

### 2-produce/HAT_export_to_numpy.py

Step 3 of 3: convert each 10 m domain raster to the .npy array HAT_dune_topo_extractor.py reads.

From the script's original header:

```text
Step 3 of 3: converts each 10 m domain raster into the .npy array that
HAT_dune_topo_extractor.py reads, matching the documented ArcGIS
RasterToNumPyArray convention:
    - nodata cells filled with -10 (not NaN)
    - no unit conversion - stays in METRES NAVD88
    - no axis transpose - rasterio's row-major read uses the same
      north-at-top / west-at-left convention as arcpy.RasterToNumPyArray

THE CONTRACT THIS HAS TO SATISFY
Read out of HAT_dune_topo_extractor.py rather than assumed:

  LOAD_PATH              INIT_ROOT/1-barrier3d-domains/{TOPO_PRODUCT}/
                         npy-arrays
  filenames (line 2557)  startswith("domain_") and endswith(".npy")
  load      (line 994)   np.load(...).astype(float), must be 2D
  nodata    (line 1015)  raw <= RAW_NODATA_MAX_NAVD (-9.0); raw nodata is
                         exactly -10.0 m NAVD88
  units     (line 1017)  z = raw - MHW_M, so raw must be m NAVD88
  shape                  ALONG_COLS=50 alongshore, TOPO_ROWS=200 cross-shore,
                         OCEAN_LOC="right" -> ocean at the HIGH column index

Our 10 m rasters are 50 rows (alongshore, 500 m) x 200 cols (cross-shore,
2000 m) with east at the high column index, and east is the ocean side. So the
arrays go through as-is: no transpose, no flip, no unit change.

THE FILENAME IS NOT domain_N_topography_2009.npy
That is what the extractor WRITES. What it READS is domain_<N>.npy. Getting
this backwards produces a folder the extractor silently finds zero domains in.

ONE FOLDER PER START PERIOD
Output goes to 1-barrier3d-domains/<TOPO_TARGET>/npy-arrays/, where TOPO_TARGET
is the period the arrays are for - "1984-start" or "2004-start". The two periods
start from different DEMs:

    --product 2009-2014-1996  ->  1984-start
    --product 2009-2014       ->  2004-start

Before 2026-08-25 the tree was keyed on the DEM year and both periods read one
set of arrays. Dune picks are keyed per version WITHIN a product, so a new
version starts from defaults rather than inheriting another version's windows.

THE SURVEY ARRAYS GO IN A SIBLING FOLDER, DELIBERATELY
The extractor globs domain_*.npy. A survey array named domain_5_survey.npy would
match that glob and be loaded as if it were a domain. So the survey arrays go to
a separate directory that the extractor never looks at, under the same
domain_<N>.npy name.

INPUTS  (data/hatteras_init/0-elevation/<PRODUCT>/2-resampled-10m/)
    resampled_domain_<N>_filled.tif
    resampled_domain_<N>_survey.tif

OUTPUTS
    data/hatteras_init/1-barrier3d-domains/<TOPO_TARGET>/
        npy-arrays/domain_<N>.npy          m NAVD88, -10 nodata
        npy-arrays_survey/domain_<N>.npy   provenance codes

    python HAT_export_to_numpy.py --product 2009-2014-1996   # -> 1984-start
    python HAT_export_to_numpy.py                            # -> 2004-start

Requires: rasterio, numpy
```

Notes that were in the code:

```text
Which elevation product to export. "2009-2014" is the baseline;
"2009-2014-1996" is the 1984-start DEM. Resolved through
scripts/site_layer/hat_elevation_products.py so a layout change cannot leave this
pointing at a directory that is no longer there.
```

```text
WHICH PERIOD PRODUCT these arrays are for. Must match TOPO_PRODUCT in
HAT_dune_topo_extractor.py - that is the script that reads them.
--product 2009-2014-1996  ->  TOPO_TARGET "1984-start"
--product 2009-2014       ->  TOPO_TARGET "2004-start"
```

```text
Every non-base code this product's survey rasters may carry. Taken from the
resolver rather than hardcoded to 2014: the 1984 product also carries 1996,
and a hardcoded 2014 reported its fill count as if the 1996 graft were not
there.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data/hatteras_init.

NOT parents[N]. This file moved into 2-produce/ on 2026-08-25, and the
old parents[3] then resolved to input_prep/ rather than the project root.
That raises nothing - it just makes every path below it wrong, silently,
until some glob comes back empty. Same helper and same reason as
4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
```

</details>

### 3-figures/HAT_plot_1984_mosaic.py

Review figures for the 1984-start DEM: what each of the three surveys contributed, and how far beach start moved.

From the script's original header:

```text
Review figures for the 1984-start topography: what each of the THREE sources
contributed, and how far the beach start moved as a result.

HAT_plot_gapfill.py draws the two-source 2009-start product and is left alone.
This is a separate script rather than a fourth entry in its SOURCES dict
because the 1984 product is a different thing to look at - it has a third
source, a boundary (the 1984 NC-12 line) that the two-source product has no
concept of, and a downstream consequence (the window origin moving) that is the
main reason to review it at all.

OUTPUTS (data/hatteras_init/0-elevation/2009-2014-1996/figures/)
    HAT_mosaic1984_island.png      3 panels, whole island
    HAT_mosaic1984_zoom_76_81.png  the developed reach, + both NC-12 lines
    HAT_mosaic1984_roads_8_15.png  the southern end, + both NC-12 lines
    HAT_mosaic1984_roads_82_88.png the northern reach, + both NC-12 lines
    HAT_mosaic1984_shift.png       per-domain start_beach shift

    Each product owns its figures/ folder, so these cannot collide with
    the baseline product's HAT_gapfill_* figures.

THE THREE PANELS
    (a) 2009 survey         what the base survey measured; everything else blank
    (b) 1984-start surface  1996 ocean-side, 2009, 2014 - the product
    (c) survey source       which survey each cell came from

(a) and (b) differ only in filled/overridden cells, so flipping between them
shows the whole intervention at once. (c) is the same information as a
categorical map, which is the only readable form where the 1996 band is thin.

STYLE. Every figure here is drawn under `scripts/site_layer/hat_figure_style.py` at the
printed width (190 mm), and carries no title, statistics line or footnote on the
canvas: that text is written to CAPTIONS.md beside the PNGs. The terrain ramp is
the house style's one sanctioned exception to drawing elevation in classes.

COLOUR
Elevation uses `terrain` with a DERIVED vmin, for the reason HAT_plot_gapfill.py
gives at length: terrain's blue water band occupies the first 25% of the ramp,
so sea level has to land exactly on that internal break or the map draws dry
ground as water. vmin = SEA_LEVEL - (vmax - SEA_LEVEL)/3, never a percentile.

The categorical panel extends that script's existing ladder rather than
inventing a palette, so the two products' figures stay comparable - a reader
who has learned that blue means "2009 measured" and green means "filled" does
not have to relearn it here:

    2009 measured   #2353b9   terrain(0.05), water blue      Y =  26
    1996 ALACE      #e6550d   ColorBrewer Oranges            Y =  60
    2014 fill       #31d670   terrain(0.30), low-land green  Y = 127
    never surveyed  #E4E4E4   neutral grey, off the ramp     Y = 198

The 1996 step is NOT sampled from terrain, and that is deliberate. Every
terrain slot in the usable lightness band is a blue, a green or a tan, and all
of them sit too close to the two colours already spoken for once the palette is
run through a colour-vision-deficiency check. An orange is the nearest hue
family terrain does not use, and it also reads correctly as "the thing this
product added".

VALIDATED, NOT EYEBALLED. The dataviz skill's JS validator needs node, which is
not on this machine - the same wall HAT_plot_gapfill.py hit - so the identical
checks were computed in Python: OKLab dE x100 between every pair, under normal
vision and under Vienot deuteranopia / protanopia / tritanopia simulation.

    worst pair, any deficiency   1996 vs 2014, protanopia   dE 16.6   (target >= 8)
    worst pair, normal vision    1996 vs 2014               dE 34.2   (floor 15)

The luminance ladder 26 / 60 / 127 / 198 is monotonic with gaps of 34, 67 and
71, so all four also separate in greyscale and in print with no hue at all.

ONE PRE-EXISTING WEAKNESS, INHERITED AND NOT INTRODUCED: the 2009-blue against
2014-green pair scores dE 5.6 under TRITANOPIA, below the target. It comes from
HAT_plot_gapfill.py's existing palette, not from anything added here, and
changing those two would break comparability with every figure already
published from the 2009-start product. Tritanopia is far rarer than the
red-green deficiencies (both of which that pair clears comfortably, at 25.5 and
30.3), and the luminance ladder carries the distinction regardless. Recorded so
it is a known trade rather than an unnoticed one.

Nodata is neutral grey in every panel and never a step on the elevation ramp:
"not surveyed" is not a low elevation, and conflating the two is what drowned
three roadways at t=0 in the first place.

    python HAT_plot_1984_mosaic.py

Requires: rasterio, geopandas, numpy, matplotlib
```

Notes that were in the code:

```text
NC-12 alignments. EPSG:2264 (NC State Plane, US survey FEET) while the maps
are EPSG:3725 (UTM 18N, metres), so they are reprojected on load - plotted raw
they would land thousands of km off the map.

BOTH are drawn, and the styling matches HAT_plot_gapfill.py so a reader moving
between the two products' figures does not have to relearn the key. Two
vintages of the same line, so they take the house vintage pair: the EARLIER
alignment (1984) red, the LATER one (2004) blue. The two are very nearly
coincident over much of the island, so 2004 is solid underneath and 1984 is
dashed on top - where they coincide you see a blue line with red dashes, and
where they diverge each is legible alone. Both carry a white casing so they
survive terrain running from dark water to near-white dune crest.
```

```text
Keyed by PERIOD; the files are the 1978 and 2008 LINES those periods read
(hat_topo_version.ROAD_LINE_FOR_YEAR), filed by vintage since 2026-09-15.
```

```text
Chrome geometry, matching HAT_plot_gapfill.py exactly. Height is DERIVED from
the data aspect at draw time: the panels are set_aspect("equal") and a zoom
extent is usually wider than tall, so matplotlib shrinks each axes to match
and a fixed tall figure leaves slack that constrained_layout splits above and
below - which reads as a large empty gap under the title. ZOOM_CHROME_IN is
the vertical allowance for panel titles, colorbar and tick labels, which do
not scale with the map.

The width is the house double-column width (190 mm) since 2026-09-10: a figure
is drawn at the width it is printed, so its 8-9 pt type is 8-9 pt on the page.
It was 15 in, which reduced to a page turned every label into 4 pt.
```

```text
Upper-left: on these zooms the island runs up the centre-right of the frame,
so the top-left corner is the one reliably empty area.
```

```text
(domain ids, filename slug, title). BOTH NC-12 alignments are drawn on every
zoom. An earlier version restricted each product to the alignment
contemporaneous with it - 1984 here, 2004 on the baseline DEM - on the
grounds that comparing a road to a DEM holding no information from its era
invites a false reading. That was overruled deliberately: seeing where the
road WAS against where it WENT is the point of the comparison, and the two
products' 8-15 views are meant to be read as a pair.
```

```text
Subtitle numbers are from mosaic_1984_audit.csv, not from eyeballing the
map. The first draft said "mostly NEW land", which the audit contradicts:
this reach is 28% new against 72% overwrite. RE-READ AFTER THE 2026-08-26
NO-BOUNDARY REBUILD: it used to be 45% against an 17% island-wide share,
i.e. nearly 3x. Dropping the boundary admitted a large band of backdune
that 2009 HAD seen, so both shares fell and the ratio with them - the
reach is now 1.75x the island figure, not 3x. Still the highest on the
island, which is the point, but the figure must not keep claiming 3x.
```

```text
The Barrier3D cell. The shift figure's whole point is which domains moved by
at least one of these, so it is drawn rather than left to the reader.
```

```text
Domain boxes are 505 m apart on a 500 m extent, so they overlap.
First writer wins, as in HAT_plot_gapfill.load_mosaic - a later
domain must not silently repaint its neighbour.
```

```text
The method paragraph. It used to be printed under every figure as a
`fig.text` footnote; under the house style nothing on the canvas belongs in a
caption, so it is now the tail of each CAPTIONS.md entry instead.
```

```text
The legend swatches carry only the year, because the panels are too narrow for
the full wording; the surveys are named here instead.
```

```text
SAME CANVAS AS HAT_plot_gapfill.py's island figure, deliberately: the two
products are read side by side, so they have to share panel proportions
and gaps. The island spans ~9 km east-west and ~47 km north-south at equal
aspect, so panel width follows figure HEIGHT, not the width asked for.
At the house double-column width (190 mm) a full page of height gives
three ~40 mm panels, which is the whole strip at one look; the figure is
no longer 13 x 19 in reduced to a page, where the type became 4 pt.
A shade under FIG_H_MAX: the legend sits outside the axes and
bbox_inches="tight" adds it after layout, so the saved page is ~0.1 in
taller than the canvas asked for.
```

```text
Wording parallels the gapfill figure's panel titles. Cell counts are in
the caption, not on the canvas.
```

```text
Concise survey labels: the panels are ~40 mm wide, and the full wording is
wider than the axes it would sit in. The sources are named in the caption.
```

```text
sharey: all three panels show the same extent, so repeating the northing
labels three times only narrows the maps.
```

```text
The +/- one-cell band is a reference threshold, so it takes the house
reference green rather than a second grey: the village bands at the top
of the panel are grey, and two greys on one panel cannot be told apart.
```

```text
Domains 1-7 get no 1996 and so have a shift of exactly zero - a bar of
zero height, which is invisible. A legend swatch pointing at nothing
reads as "these are missing from the chart", so the span is shaded and
labelled in place instead.
```

```text
Explicit limits rather than margins(): the shifts run -5 to +72, and a
symmetric margin then opened a band of empty axes below -20.
```

```text
Villages as a strip against the top edge, so they cannot be confused with
the shaded domain span below them. After set_xlim, as the helper
requires, and after the span so Buxton's name is not painted over.
```

<details><summary>Function notes (the original docstrings)</summary>

**`elev_limits()`**

```text
vmax from a percentile, vmin DERIVED so 0 m lands on terrain's internal
water/land break. Never both from percentiles - see the docstring.
```

**`panel_survey()`**

```text
Categorical provenance. Codes are mapped to contiguous indices so the
colours cannot slide if a code is absent from a crop.
```

**`survey_legend()`**

```text
Handles for the three surveys plus the unsurveyed background.

concise=True keeps only the year. The island figure needs it: its panels
are ~2 in wide (a 46 km island in an 8 km window at equal aspect) and the
full labels are wider than the axes they sit in, so the legend spilled
across the panel border. The sources are named in this module's docstring,
so the year is enough to read the colours by.
```

**`load_roads()`**

```text
Both NC-12 alignments, reprojected and CLIPPED to the domain footprint.

The geojsons run the full length of the highway, well beyond the 90
domains at both ends. Unclipped, a figure shows road where there is no
model domain, which reads as coverage that does not exist.
```

**`draw_roads()`**

```text
BOTH casings first, then both lines in ROAD_ORDER so the dashed 1984
lands on top of the solid 2004 rather than under it.

Casings-then-lines, not casing-line-casing-line: the two alignments are
nearly coincident for most of the island, and the second casing then
painted out the first line, so 2004 disappeared wherever it mattered.
`scale` thins the lines for the island figure, where the same widths would
smother the island.
```

**`road_legend_handles()`**

```text
Both alignments in their map colours. They are the house vintage pair,
so neither is white and the swatches carry straight over.
```

**`fig_shift()`**

```text
Per-domain movement of the extractor's beach start.

A diverging encoding, because the quantity has a meaningful zero and a
meaningful sign: seaward carries the 1996 orange, landward the 2009 blue,
so the bar colour says which survey won that domain. Zero is a neutral
rule, not a hue.

The +/-10 m band is drawn because it is the only threshold that matters -
the extractor works at 10 m, so anything inside that band cannot move a
Barrier3D cell by more than one, and a reader should not have to do the
arithmetic to see which bars clear it.
```

</details>

### 3-figures/HAT_plot_dem_holes.py

Where the nodata and sub-MHW holes are in a Barrier3D start DEM.

From the script's original header:

```text
Where the nodata and sub-MHW holes are in a Barrier3D start DEM.

WHAT THIS ANSWERS
"Do I have a lot of cells by my dunes that are under MHW, or data gaps?"
Counting them per domain (the audit CSVs) says how many. This says WHERE, and
separates the two things that both look like "low" on an elevation ramp:

    a cell that is water because it is the ocean or the sound   - expected
    a cell that is water INSIDE the island                      - a hole

Only the second is coloured loudly. Everything else is deliberately pale, so a
figure of a clean DEM is a quiet figure.

THE CATEGORIES, per profile, in the extractor's own frame
Profiles are read ocean-first (arr[:, ::-1], OCEAN_LOC="right"), z = raw - MHW,
exactly as HAT_dune_topo_extractor.load_domain does.

    beach_start  first cell with z > 0.50 m       (BEACH_START_THR_M)
    last_land    last cell with z > 0

Cells seaward of beach_start or landward of last_land are open water. Cells
BETWEEN them that are not land are holes:

    land          z > 0                             pale sand
    open water    outside [beach_start, last_land]  pale blue
    wet hole      z <= 0, valid data, inside        strong blue
    gap           nodata (raw <= -9), inside        strong red

Nodata is never a step on an elevation ramp here, for the reason FIGURES.md
gives: "not surveyed" is not a low elevation, and conflating the two is what
drowned three roadways at t=0.

The two affected categories are #1f6fb4 blue against #d7191c red. The red is
saturated on purpose. Unsurveyed cells are the category a reader must not skim
past - they are the ones that become a fictitious elevation downstream - and
they are also the rarer of the two at 0.32% of cells, so at this scale they
have to hold their own against the blue at single-pixel widths. A muted
#d6604d, sampled off the same RdBu ramp as the blue, was tried and is the
better choice for a figure meant to sit quietly in a page of body text; it lost
too much at one-cell width here.

Under a red-green deficiency the red darkens towards olive while the blue
holds, so the pair still parts. Magenta against this blue, the first draft, is
the pair to avoid: it holds neither the hue nor the luminance gap.

ORIENTATION
Ocean is at the BOTTOM of panels A and B, landward upward. That is the
convention HAT_dune_topo_extractor.pick_window draws for picking a dune search
window, so this figure and the picker read the same way round.

THE THREE PANELS
A  the island unrolled. Alongshore runs left-right; cross-shore is UTM easting
   with the island trend removed by a cubic fit through the 90 domain origins,
   so a 45 km arc lies flat instead of drifting 5 km across the panel. The
   detrend is a rigid per-domain shift - no cell is resampled, and cross-shore
   distances within a domain are untouched. This is the locator panel: it keeps
   each domain's own shoreline shape, which panel B removes. The cross-shore
   axis is metres landward of the seaward edge of the detrended strip, so its
   zero is a drawing origin and not a landform.

B  the same cells straightened: every profile shifted so its own beach_start
   sits at cross-shore 0. The dune band becomes a horizontal stripe instead of
   following the shoreline curve, so a hole IN THE DUNES is separable from a
   hole 500 m behind them. Panel A cannot show that; the shoreline moves.

C  per-domain percentages, dune band vs interior.

ONE ALONGSHORE AXIS
All three panels share x, in kilometres of UTM northing measured from the
southern edge of domain 1. Every domain is painted at its true northing, so a
vertical line means the same place in all three panels. Domains are spaced
~504 m and carry 500 m of data, so there are ~4 m unpainted seams between them;
that is real, not a plotting artefact.

Left to right is south -> north, which is the model alongshore direction
(ALONGSHORE_FLIP = True flips the raster north-at-top rows so profile index
increases northward). Within a domain, profile p is raster row 49 - p.

INPUT   data/hatteras_init/1-barrier3d-domains/<PRODUCT>/npy-arrays/domain_<N>.npy
        the arrays CASCADE reads - m NAVD88, -10 nodata - not the .tif, so what
        is drawn is what the model ingests.
        Georeferencing comes from the elevation product resample_audit.csv.

OUTPUT  <elevation product>/figures/HAT_<slug>_holes.png
```

Notes that were in the code:

```text
--- the shared alongshore axis -----------------------------------------
One 10 m grid in UTM northing. Domain n occupies columns [x0, x0+50),
counting from the SOUTH, so column index increases northward like the
model profile index does.
```

```text
--- panel A: unroll the island -----------------------------------------
Remove the island trend from easting so the strip lies flat. The fit is
evaluated once per domain, so each block is shifted rigidly - no cell is
resampled and no cross-shore distance changes.
```

```text
cats[n] is (profile S->N, cross ocean-first). Flip cross so index
increases eastward, then transpose to (easting, alongshore).
```

```text
Row 0 is the WEST (sound) edge. Flip so row 0 is the ocean edge, then
origin="lower" puts the ocean at the bottom with landward upward, the
same way round as the extractor's pick_window.
```

<details><summary>Function notes (the original docstrings)</summary>

**`classify()`**

```text
(n_along, n_cross) ocean-first raw m NAVD88 -> codes, beach_start, last_land.

beach_start / last_land are -1 on a profile with no land at all. None exist
in this product, but a forecast domain could have one.
```

**`load_domain()`**

```text
Raw array, ocean-first cross-shore, profile index increasing NORTHWARD.

The .npy is raster order: row 0 north, column 199 east = ocean. [:, ::-1]
puts the ocean first, which is what the extractor does. [::-1] on the rows
is ALONGSHORE_FLIP, so profile 0 is the southern edge of the domain.
```

</details>

### 3-figures/HAT_plot_duneline_offset.py - deleted 2026-10-19

Retired with the 1984 start: it measured the 1984 dune line against the 1996-start DEM. Its outputs (`0-elevation/2009-2014-1996-duneline/`) are on `D:\archive\CASCADE_offload\data\hatteras_init\0-elevation\`. Recover the script with `git show 2046d67b:scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py`.

### 3-figures/HAT_plot_gapfill.py

Review figure for a DEM gap fill: 2009 alone, what the fill adds, and which cells came from where.

From the script's original header:

```text
Review figure for the 2009 DEM gap fill: what the 2009 survey alone gives,
what the fill adds, and exactly which cells came from where.

Three panels, all on the same 10 m grid and the same colour scale:
    (a) 2009 survey      cells the 2009 DEM measured; everything else blank
    (b) with fill        the product that goes to the dune/topo extractor
    (c) survey source    2009 measured / fill-year filled / never surveyed

Panels (a) and (b) differ ONLY in the filled cells, so flipping between them
shows the fill directly. Panel (c) is the same information as a categorical map,
which is easier to read where the fill is thin.

STYLE. Every figure here is drawn under `scripts/site_layer/hat_figure_style.py` at the
printed width (190 mm), and carries no title, statistics line or footnote on the
canvas: that text is written to CAPTIONS.md beside the PNGs. The terrain ramp is
the house style's one sanctioned exception to drawing elevation in classes.

Domain boxes are drawn over every panel: no fill, thin white outline, so they
locate a domain without hiding the data under it.

COLOUR
Elevation uses `terrain`, and the reason is specific: it is built for
topography, with a blue water band occupying the first 25% of the ramp. That
only reads correctly if sea level lands exactly on that internal break, hence
vmin is DERIVED - vmin = SEA_LEVEL_M - (vmax - SEA_LEVEL_M) / 3 - rather than
taken from a percentile. Set vmin from a percentile and the blue/green boundary
drifts to an arbitrary elevation, drawing dry ground as water.

The categorical panel's two colours are sampled from `terrain` itself,
terrain(0.05) water blue for 2009 and terrain(0.30) low-land green for the
fill, so it belongs to the same palette as the elevation maps. They carry a
deliberate luminance ladder - 80 / 153 / 228 against the grey, spacing 73 and
75 - because blue-vs-green is a weak colour-vision-deficiency axis and
brightness has to carry what hue cannot.

Nodata is neutral grey in every panel and never a step on the elevation ramp -
"not surveyed" is not a low elevation, and the whole point of this work is that
conflating those two drowned three roadways at t=0.

SOURCE SELECTION
Tag, fill year and long label live together in SOURCES, so a re-render cannot
put one source's label on another's data - hand-editing the three constants
separately already printed "cells the 2014 DEM measured" on the 2008 figures
once.

    python HAT_plot_gapfill.py                       # 2009-2014
    python HAT_plot_gapfill.py --source <PRODUCT>

INPUT   data/hatteras_init/0-elevation/<SOURCE_TAG>/2-resampled-10m/
            resampled_domain_<N>_filled.tif
            resampled_domain_<N>_survey.tif
OUTPUT  data/hatteras_init/0-elevation/<SOURCE_TAG>/figures/
            HAT_gapfill_<SOURCE_TAG>_island.png         whole island, 3 panels
            HAT_gapfill_<SOURCE_TAG>_domains_78_80.png  zoom on 78-80
            HAT_gapfill_<SOURCE_TAG>_roads_78_80.png    the zoom, + NC-12
        SOURCE_TAG is in the name, so a new source cannot overwrite the
        existing figures.

Requires: rasterio, geopandas, numpy, matplotlib
```

Notes that were in the code:

```text
NC-12 alignments. These are EPSG:2264 (NC State Plane, US survey FEET) while
the maps are EPSG:3725 (UTM 18N, metres), so they are reprojected on load -
plotted raw they would land thousands of km off the map.
```

```text
Keyed by PERIOD; the files are the 1978 and 2008 LINES those periods read
(hat_topo_version.ROAD_LINE_FOR_YEAR), filed by vintage since 2026-09-15.
```

```text
Two vintages of the same line, so they take the house vintage pair: the
EARLIER alignment (1984) red, the LATER one (2004) blue. They are very nearly
coincident through 78-80, so 2004 is solid underneath and 1984 dashed on top -
where they coincide you see a blue line with red dashes, and where they
diverge each is legible on its own. Both carry a white casing so they survive
terrain running from dark water to near-white dune crest.
```

```text
Fill sources this script knows how to plot. Keeping tag, year and label in ONE
place stops them drifting apart - hand-editing three constants per re-render
already put "cells the 2014 DEM measured" on the 2008 figures once.

python HAT_plot_gapfill.py                      # default (below)
python HAT_plot_gapfill.py --source <PRODUCT>

Keyed by PRODUCT, not by fill source - the directories were renamed for
composition on 2026-08-25 (2014_NOAA_PostSandy -> 2009-2014).

2008_NOAA_IOCM was an entry here until 2026-08-26. Its rasters are gone from
disk (they were never tracked - *.tif is gitignored), and the point-cloud
path that produced them has been removed from HAT_dem_gap_fill.py, so
--source 2008_NOAA_IOCM could not have re-rendered anything. A dead option
that reads as a live one is worse than no option.
```

```text
Built from SURVEY_FILL so it cannot disagree with the data being plotted.
It used to be printed under every figure as a `fig.text` footnote; under the
house style nothing on the canvas belongs in a caption, so it is now the tail
of each CAPTIONS.md entry instead.
```

```text
The superseded fallback that used to live here is gone: the resolver knows
which products are superseded and where they sit, so there is one place that
has to be right rather than a probe-two-paths-and-hope in every consumer.
A missing directory is reported there, not as an empty glob later.
```

```text
Breathing room around the island-wide mosaic. Without it the northernmost and
southernmost domains sit flush against the axes frame, which reads as the data
being cut off rather than ending.
```

```text
Zoom figure geometry. Height is derived from the data aspect at draw time;
ZOOM_CHROME_IN is the vertical allowance for the suptitle, colorbar and tick
labels, which do not scale with the map. NOT the footnote - that sits outside
the axes and bbox_inches="tight" adds it after layout, so reserving space for
it just opens a gap under the title. Measured: gap above the axes tracks this
value almost 1:1, and a two-line suptitle needs ~0.35 in.
Legends sit OUTSIDE the axes, under the panels, since 2026-09-10: an inside
legend on a map this narrow either covers the island or sits in the nodata
grey, and the house rule is frameless and outside wherever the layout allows.
```

```text
Road-overlay zooms: (domain ids, which NC-12 years, filename slug, title).

78-80 carries BOTH alignments - those are the domains the extractor names as
width-drowning at t=0, and seeing 1984 against 2004 there is the point.

8-15 carries BOTH as well. An earlier version drew only 2004 here, on the
grounds that it is the alignment contemporaneous with this DEM and that
putting the 1984 line over a 2009+2014 surface invites comparing a road to a
DEM holding no information from its era. That was overruled deliberately:
seeing where the road WAS against where it WENT is the point, and this view
is meant to be read as a pair with the 1984-start DEM's own 8-15 figure,
which now draws the same two lines.
The third element is the figure's own caption sentence; the common method
paragraph (SOURCE_NOTE) and panel key are appended when it is written.
```

```text
The house double-column width (190 mm) since 2026-09-10: a figure is drawn at
the width it is printed, so its 8-9 pt type is 8-9 pt on the page. It was
15 in, which reduced to a page turned every label into 4 pt.
```

```text
Categorical colours for the survey-source panel, SAMPLED FROM `terrain` so
panel C is built from the same palette as the elevation maps:

C_2009  terrain(0.05)  #2353b9  deep water blue
C_FILL  terrain(0.30)  #31d670  the green of terrain's low-land band

Luminance ladder, which is what keeps the three readable:

2009 measured  #2353b9   luminance  80
fill           #31d670   luminance 153   (73 from blue, 75 from grey)
never surveyed #E4E4E4   luminance 228

Spacing 73 and 78 is nearly even, so all three separate by brightness alone.
That matters more here than usual: blue-vs-green is the WEAKEST colour-vision
-deficiency axis of the pairings tried, so brightness is doing the work that
hue cannot be relied on for. Without the gap this pair would be a poor choice;
with it, it holds up in greyscale and under CVD.

The dataviz validator could not be run here (no node on this machine), so the
ladder was computed directly rather than machine-checked.
```

```text
matplotlib's `terrain` is built for topography: its blue water band occupies
the FIRST 25% of the ramp, then green -> brown -> white for land. That is only
meaningful if sea level lands exactly on that internal boundary, so vmin is
derived rather than taken from a percentile:

0 maps to  |vmin| / (|vmin| + vmax)  ==  0.25   ->   vmin = -vmax / 3

Set from a percentile instead and the blue/green break drifts to some
arbitrary elevation, so the map would draw dry ground as water or vice versa.
```

```text
Domains are placed independently and their boxes are not a perfect
tiling (505 m spacing against a 500 m extent), so overlaps exist.
Keep whatever is already there rather than letting the later domain
silently overwrite its neighbour.
```

```text
Decimal places from the span, not fixed: the island figure covers ~45 km
of northing where 0 dp is right, the zoom covers ~1.5 km where 0 dp
renders every tick as the same number.
```

```text
The island spans ~9 km east-west and ~47 km north-south at equal aspect,
so panel width follows figure HEIGHT, not the width asked for. At the
house double-column width (190 mm) a full page of height gives three
~40 mm panels, which is the whole strip at one look; the figure is no
longer 13 x 19 in reduced to a page, where the type became 4 pt. A shade
under FIG_H_MAX because the legend sits outside the axes and
bbox_inches="tight" adds it after layout.
```

```text
Boundaries must ascend and the fill year (2014) is now GREATER than the
measured year (2009), so measured precedes fill in the colour list. With
the old 2008 ordering these two were swapped and the map lied.
```

```text
thinner at island scale: 45 km of line at zoom widths would smother
the island. Dashed vs solid does not resolve at this scale - the zoom
figure is where that distinction is readable.
```

```text
ax=all three, not axes[:2] - a colorbar sized against a subset shrinks
only those axes and leaves the rest misaligned.
```

```text
Height from the DATA aspect, not hard-coded. The panels are
set_aspect("equal") and the zoom extent is wider than it is tall, so
matplotlib shrinks each axes to match and a fixed tall figure leaves
slack that constrained_layout splits above and below the panels -
which reads as a big empty gap under the title.
```

```text
sharey: all three panels show the same extent, so repeating the
northing labels three times only narrows the maps.
```

```text
---- road-overlay zooms, one per entry in ROAD_ZOOMS ----
Drawn at zoom rather than island scale on purpose: island-wide the
road is a ~1 px line over 45 km, where dashed and solid are
indistinguishable and the overlay would carry no information.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_find_project_root()`**

```text
Walk up until a directory holds data/hatteras_init.

NOT parents[N]. This file moved into 3-figures/ on 2026-08-25, and the
old parents[3] then resolved to input_prep/ rather than the project root.
That raises nothing - it just makes every path below it wrong, silently,
until some glob comes back empty. Same helper and same reason as
4-mgmt-forcings/road_offset/2-audit/HAT_road_setback_audit.py.
```

**`km_axes()`**

```text
UTM eastings here are 6-digit metres (450439..458392). Three panels side by
side cannot fit those without colliding, and matplotlib's shared offset
label is easy to miss. Kilometres with a small tick count is legible at any
panel width and needs no offset text.
```

**`load_roads()`**

```text
Loads the NC-12 alignments, reprojects them, and CLIPS them to the domain
footprint.

The geojsons run the full length of the highway, well beyond the 90 domains
at both ends. Unclipped, the island figure shows road where there is no
model domain, which reads as coverage that does not exist.
```

**`draw_roads()`**

```text
BOTH casings first, then both lines in ROAD_ORDER so the dashed 1984 lands
on top of the solid 2004 rather than under it (see ROAD_STYLE).

Casings-then-lines, not casing-line-casing-line: the two alignments are
nearly coincident through the reaches these figures zoom on, and the second
casing then painted out the first line, so 2004 disappeared wherever it
mattered. `scale` thins the lines for the island-wide figure, where the same
widths would smother the island.
```

**`road_legend_handles()`**

```text
Both alignments in their map colours. They are the house vintage pair,
so neither is white and the swatches carry straight over.
```

**`_road_zoom()`**

```text
One A/B/C zoom with the named NC-12 alignments drawn over it.

`years` selects WHICH alignments, and every entry in ROAD_ZOOMS
currently asks for both. It stays a parameter rather than being
hardcoded because the restriction was tried and reversed once
already - see the note on ROAD_ZOOMS - and a future zoom may well
want one line only.
```

</details>
