# 2009 DEM gap fill - current source: 2014 NOAA Post-Sandy

## Datasets

| role | dataset |
|---|---|
| base DEM | `usace2009_nc_dem_Job1076020 / 2009_full.tif` - 1 m, EPSG:3725 + NAVD88 (m) |
| fill source | `2014_NOAA_Post_Sandy_DEM_Job1076021 / 2014_full.tif` - 1 m, EPSG:6347 + NAVD88 (m), reprojected on read |
| domains | `D:\Hatteras_GIS\domains.geojson` - 90 boxes, 2000 x 500 m |
| NC-12 lines | `4-mgmt-forcing/road_offset/raw_offset/{1984,2004}/nc12_*.geojson` - EPSG:2264 (NC State Plane, US survey FEET), reprojected and clipped on load |

## Figures

| file | shows |
|---|---|
| `HAT_gapfill_2014_NOAA_PostSandy_island.png` | whole island: 2009 only / 2009+2014 / survey source, with NC-12 |
| `HAT_gapfill_2014_NOAA_PostSandy_domains_78_80.png` | same panels, zoomed on the road-drowning domains |
| `HAT_gapfill_2014_NOAA_PostSandy_roads_78_80.png` | the zoom with both NC-12 alignments overlaid |
| `HAT_gapfill_2008_NOAA_IOCM_*.png` | **superseded** first attempt, kept for comparison |

All three share panels A (2009 only), B (2009 + fill) and C (survey source). A and
B differ *only* in filled cells, so flipping between them shows the fill directly.

## Figure design decisions

These are choices with reasons, not defaults - changing them without the reason
will quietly break something:

**Elevation uses `terrain`, with `vmin` derived rather than taken from a
percentile.** `terrain` puts its blue water band in the first 25% of the ramp,
so sea level has to land exactly on that internal break or the map draws dry
ground as water. Hence `vmin = SEA_LEVEL_M - (vmax - SEA_LEVEL_M) / 3`. Set
`vmin` from a percentile and the blue/green boundary drifts to an arbitrary
elevation. `SEA_LEVEL_M` is 0.0 m NAVD88; set it to 0.36 to key the break to MHW.

**Panel C's two colours are sampled from `terrain` itself** - `terrain(0.05)`
water blue for 2009, `terrain(0.30)` low-land green for the fill - so the
categorical panel belongs to the same palette as the elevation maps. They carry
a deliberate luminance ladder, 80 / 153 / 228 against the grey, spacing 73 and
75. That near-even spacing is doing real work: blue-vs-green is a weak
colour-vision-deficiency axis, so brightness carries what hue cannot. Pick a
prettier pair at similar lightness and the panel stops working in greyscale and
under CVD.

**Nodata is neutral grey in every panel and never a step on the elevation
ramp.** "Not surveyed" is not a low elevation, and conflating those two is what
drowned three roadways at t=0 in the first place.

**NC-12: 2004 solid black underneath, 1984 white-dashed on top, each with a
casing in the opposite colour.** The two alignments are nearly coincident
through 78-80, so drawing both in black with the solid one last makes 2004 paint
straight over 1984 and only ONE road appears. The casings are needed because the
lines cross terrain running from dark blue water to near-white dune crest.

**Roads are clipped to the union of the domain boxes** (42.8 of 58.4 km kept for
1984, 42.7 of 58.4 for 2004). The geojsons run the full length of the highway,
well past the modelled reach at both ends; unclipped they imply model coverage
that does not exist. Clipping to the union rather than a bounding box also drops
the line in the gaps between domains.

**At island scale the road is drawn at 45% width and dashed-vs-solid does not
resolve** - 45 km of line at zoom widths smothers the island. The island view
shows *where* NC-12 runs; the roads zoom is where the vintages are separable.

**Axes are kilometres with few ticks, decimals set from the span.** Six-digit UTM
metres collide across three narrow panels. The island spans ~45 km of northing
where 0 dp is right; the zoom spans ~1.5 km where 0 dp renders every tick
identically.

## Why 2014 NOAA

Every candidate DEM scored against the 2009 **dry-land** gaps (cells 2009 missed
where the consensus of candidates puts the ground above MHW), all 90 domains:

| dataset | coverage | |
|---|---|---|
| 2014 NCFMP | 100.00% | **disqualified** - hydro-flattened, 70.5% of its gap values are the single constant -0.762 m |
| **2014 NOAA Post-Sandy** | **97.34%** | earliest genuine, and the best |
| 2017 USACE | 23.74% | |
| 2016 post-Matthew | 21.98% | |
| 2019 DUNEX | 21.98% | |
| 2018 post-Florence | 9.74% | |

2016/2017/2019 score 95-100% on domains 78-80 but are localised surveys - 2019 is
below 50% in 60 of 90 domains. Choosing on the developed reaches alone would have
picked a dataset covering under a quarter of the island's gaps.

Superseded point-cloud attempts, both topographic-only, in the strip landward of
NC-12 (of 150 cells): 2008 NOAA IOCM 17, 2011 post-Irene 28, combined 42.

## Fill rules - selection only, no value is modified

1. **Coverage** - the 2014 DEM has a value
2. **Connectivity** - contiguous with the island, bridging gaps <= 20 m (chosen
   from measurement: 20 m recovers 30.4% of detached *cells* against 16.2% of
   components, on the plateau after that step)
3. **Elevation** - >= -2.64 m NAVD88, the extractor's own `WATER_CLAMP_M`

Bias correction and feathering are **off**. A filled cell is the 2014
measurement unchanged.

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

## The seam, measured

| | n | median abs | p90 | p99 | signed |
|---|---|---|---|---|---|
| seam, dry land | 4,507 | **0.158** | 0.519 | 2.156 | +0.092 |
| seam, wet (<= MHW) | 3,494 | 0.477 | 1.115 | 1.640 | -0.477 |
| control, measured-to-measured | 12,583,950 | 0.027 | 0.109 | 0.297 | +0.000 |

Where it matters - the dry developed interior - the step is 0.158 m, 6x terrain
roughness, a **1.6% slope** at the 10 m grid. The wet seam's -0.477 m is largely
real: 2009's last cell is the water's edge, the 2014 cell beyond it is water.

Offset on co-measured cells by distance from the boundary: +0.053, +0.035,
+0.009, -0.014, -0.014, -0.026, -0.032, -0.087 m at 1, 2, 3, 5, 10, 20, 50 and
100 m. Flat within +/-0.09 m - no edge effect, no datum offset, hence nothing to
correct.

The 2.156 m p99 is not buildings: the 2014 fill has **zero** cells above 4 m in
the developed domains (max 2.88-3.47 m), while 2009 reaches 6-7 m. It is low fill
abutting dune ridges - real terrain contrast.

## Regenerating

Figures only: `python scripts/input_prep/0-elevation/HAT_plot_gapfill.py`
(reads `2-resampled-10m/`, writes all three PNGs).

Full pipeline, in order:

```
HAT_dem_gap_fill.py       clip + fill  -> 0-elevation/1-gapfill-1m/
HAT_dem_resample_clip.py  1 m -> 10 m  -> 0-elevation/2-resampled-10m/
HAT_export_to_numpy.py    -> 2009-npy-arrays/2009_pea_hatteras_filled{,_survey}/
HAT_plot_gapfill.py       -> 0-elevation/figures/
```

`HAT_laz_ground_classify.py` is **not** in this chain - it is only needed when
the fill source is an unclassified point cloud, as 2008 was. A DEM source skips
it.

### To change fill source

The fill year appears in five places and they must agree, or the survey rasters
and the figure legend will disagree about what a filled cell is:

| file | set |
|---|---|
| `HAT_dem_gap_fill.py` | `FILL_DEM_PATH`, `FILL_SOURCE_YEAR` |
| `HAT_dem_resample_clip.py` | `SURVEY_FILL` |
| `HAT_export_to_numpy.py` | `SURVEY_FILL` |
| `HAT_plot_gapfill.py` | `SURVEY_FILL`, `SOURCE_TAG`, `SOURCE_LONG`, `SOURCE_NOTE` |

`SOURCE_TAG` also names the output PNGs, so a new source will not overwrite the
current figures. Consider tagging the audit CSVs the same way - the pipeline
writes `gapfill_audit.csv` etc. to fixed names and a re-run overwrites them
(the 2008 set here was copied out with a `_2008_NOAA_IOCM` suffix for this
reason).

### Before committing to a new source

`HAT_check_fill_source.py <path>` scores a candidate point cloud in two stages,
the first taking about a second from headers alone. `HAT_survey_dem_coverage.py`
does the same for the DEM collection across all 90 domains, with the
hydro-flattening check that disqualified 2014 NCFMP.

Both are worth running first: on domains 78-80 alone, 2016/2017/2019 all looked
like fine choices at 95-100%, and island-wide they cover ~22%.

## Still open

The roadway drowning in domains 78-80 is a separate matter. The strip landward
of NC-12 is genuinely below MHW in 79 and 80 (2017 and 2019 both agree), so those
roads border real water and drowning there may be correct rather than an
artifact. Only domain 78 showed dry marsh behind the road.
