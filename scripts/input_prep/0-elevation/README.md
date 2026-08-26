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

## Datasets

| role | dataset |
|---|---|
| base DEM | `usace2009_nc_dem_Job1076020 / 2009_full.tif` - 1 m, EPSG:3725 + NAVD88 (m) |
| fill source | `2014_NOAA_Post_Sandy_DEM_Job1076021 / 2014_full.tif` - 1 m, EPSG:6347 + NAVD88 (m), reprojected on read |
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
`scripts/hat_elevation_products.py` and in `HAT_plot_gapfill.py`'s `SOURCES`.

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
| `scripts/hat_elevation_products.py` | **add the product first** - a `Product` entry and its `FILL_CODES` |
| `2-produce/HAT_dem_gap_fill.py` | `FILL_DEM_PATH`, `FILL_SOURCE_TAG`, `FILL_SOURCE_YEAR`, `PRODUCT_TAG` |
| `2-produce/HAT_dem_resample_clip.py` | nothing - pass `--product <NAME>` |
| `2-produce/HAT_export_to_numpy.py` | `SURVEY_FILL`, `DEM_NAME`; pass `--product <NAME>` |
| `3-figures/HAT_plot_gapfill.py` | add one entry to `SOURCES`, then `--source <NAME>` |

The product name is the TOP-LEVEL folder and appears in the figure filenames,
so a new product cannot overwrite an existing one or its figures. Paths are
resolved by `scripts/hat_elevation_products.py`; do not rebuild them by joining
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
confined to the ocean side of `nc12_1984.geojson`. Measurement showed that
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
`2-brie-offset/hindcast_1984/Island_Dune_Offsets_1984_CASCADE_Input.csv`,
measured independently of any DEM, already existing for 1984, passed straight
to `Cascade()`. Nothing here touches it - a deliberate choice to leave the two
independent rather than derive one from the other.

### Downstream state

`HAT_export_to_numpy.py` has been run against this tag: 90 domains sit in
`1-barrier3d-domains/1984-start/npy-arrays/`, rebuilt with the product itself
after the road boundary came out. `1984-start/dune-topo/v1` exists. What the
"Still to do" note here used to say - that nothing downstream consumed the 1984
product - is no longer true.

## Still open

The roadway drowning in domains 78-80 is a separate matter. The strip landward
of NC-12 is genuinely below MHW in 79 and 80 (2017 and 2019 both agree), so
those roads border real water and drowning there may be correct rather than an
artifact. Only domain 78 showed dry marsh behind the road.
