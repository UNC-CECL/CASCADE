# 4-mgmt-forcings - the management the hindcast applies

NC-12 (where it sits, how high, how far it moved) and beach nourishment,
as inputs the runner reads.

```
beach_nourishment.py                 when and where the beach was nourished (figures)
road_elevation/HAT_road_elevation.py the per-domain road elevation, one set for both starts
road_offset/                         the road setback per domain: produce, audit, figures, compare (own README)
road_relocation/                     how far the road moved between vintages (own README)
```


## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### beach_nourishment.py

When and where the beach was nourished: figures of the fill projects the hindcast fires.

From the script's original header:

```text
When and where the beach was nourished: figures for
data/hatteras_init/4-mgmt-forcing/nourishment/.

THE QUESTION
    The hindcast fires the projects in HATTERAS_NOURISHMENT_PROJECTS
    (scripts/site_layer/hatteras_site_config.py) when they fall inside a run window.
    That list is three entries and the record behind it is a spreadsheet
    (Hatteras_Management_Timelines.xlsx, sheet Nourishment_Timeline), and
    nothing in the data tree shows the two side by side, or shows the reader
    which stretch of the island was filled in which year. These figures do.

TWO SOURCES, DRAWN TOGETHER, NEVER MERGED
    * the MODEL INPUT: the site-config projects, their extents and their
      volumes spread to m^3/m over 500 m domains. This is what a run receives.
    * the RECORD: the spreadsheet. Its domain flags for the same three
      projects are narrower than the site config's (Rodanthe 85-88 against
      84-89; Avon 23-26 against 21-28; Buxton the same 6-15) because the site
      config re-derived the footprints from the project descriptions -- the
      reasons are in the comments beside each entry. The record also carries
      fourteen Pea Island / Oregon Inlet navigation fills (1990-2004, 2013)
      that lie NORTH of GIS 90, off the modelled reach, and that no run sees.
    The figures draw the model extent as the fill and the record's flags as
    a darker inner bar, so a difference is visible rather than reconciled.

OUTPUTS   data/hatteras_init/4-mgmt-forcing/nourishment/
    nourishment_when_where.png/.pdf       year x domain event chart
    nourishment_volume_alongshore.png/.pdf m^3/m per domain, as delivered
    nourishment_domain_map.png/.pdf       the filled domains on the island
    nourishment_projects.csv              every row drawn, both sources
    CAPTIONS.md                           the text that is not on the canvas

RUN
    python scripts/input_prep/4-mgmt-forcings/beach_nourishment.py
```

Notes that were in the code:

```text
Vintage colours: the earlier fill is red, the later one blue, as everywhere
two vintages share a figure (hat_figure_style).
```

```text
One year can flag two separate stretches (2022: Buxton AND Avon),
so each contiguous run of flags is its own row.
```

```text
Left margin: the hindcast windows as vertical brackets. Right margin: the
off-reach fills. Both are OFF the domain axis, which runs 1-90.
```

```text
Rotate 90 deg clockwise so south is left, north right, ocean below;
rotation preserves distance, so the scale bar holds.
```

```text
Every tenth domain, unless a footprint end sits within one domain of it
(20 beside 21, 90 beside 89), plus the footprint ends themselves.
```

```text
Scale bar in data units, 5 km, and a north arrow that points along the
strip since the map is rotated off north.
```

<details><summary>Function notes (the original docstrings)</summary>

**`record_projects()`**

```text
The spreadsheet: one row per year that carries a note or a flag.
Domain flags (1) give the footprint; a note with no flags is a fill the
record places outside the 90 domains (Pea Island / Oregon Inlet).
```

**`_windows_sentence()`**

```text
'1984–2004 and 1996–2010 carry no fill, 2004–2024 and 2010–2024 carry
all three', computed from the config so the caption cannot go stale.
```

</details>

### road_elevation/HAT_road_elevation.py

Per-domain NC-12 road elevation: the mean of the 2009 LiDAR under the 2008 road line, for both start periods.

From the script's original header:

```text
Per-domain NC-12 road elevation: the MEAN of the 2009 LiDAR under the 2004 road
alignment. ONE set of numbers, used for BOTH the 1984 and the 2004 start period.

Replaces HAT_road_elevation_from_lidar.py, which sampled two alignments and
wrote two files. That script was deleted on 2026-08-17; it was staged but never
committed, so it is not in any commit. Its blob survives in the object store at
becbfc878ae0afa3f4e76037ef0edff810fa857f until the next `git gc`, recoverable
with `git cat-file -p <hash>`. The reason it is gone rather than kept is in WHY
ONE FILE below.

WHY ONE FILE AND NOT ONE PER VINTAGE
There is only one DEM. Both vintages were sampled on the same 2009 surface, so
any difference between a "1984" and a "2004" road elevation was never a
difference in time -- it was a difference in WHERE ON THE 2009 SURFACE the two
digitised lines happened to fall. Where NC-12 never moved the two lines sit on
top of each other and the numbers are identical by construction. Where the road
WAS relocated the 1984 line lies over the abandoned corridor, and the "1984 road
elevation" there was the elevation of a place a road used to be.

Neither of those is a measurement of temporal change in roadbed height. Writing
two files implied one. This writes one.

WHAT THE ABANDONED CORRIDOR ACTUALLY IS -- MEASURED, NOT ASSUMED
It is tempting to assume the abandoned alignment was overwashed and bulldozed
flat, and so reads LOW in a 2009 DEM. It does not. Sampled here, the 1984 line
through the relocated domains gives a mean of about 2.4 m NAVD88 against about
1.5 m for the 2004 line, with a within-domain standard deviation up to 1.7 m --
GIS 10 comes back at 4.3 m. The dune migrated over the corridor after the road
left it. That sample is FOREDUNE, not roadbed, and not flattened ground either.

This decides the choice rather than merely complicating it: the two candidates
do NOT bracket the truth. Both sit above the natural grade of the neighbouring
un-relocated domains (~1.0 m), so the 2004 value is the lower of the two AND the
only one that is a graded surface -- the conservative choice as well as the
correct one. RELOCATION BRACKET re-measures this on every run.

WHY THE 1 m CLIP AND NOT THE 10 m RESAMPLE
NC-12 is a two-lane road: roughly 7-10 m of pavement plus shoulder. On the 10 m
Barrier3D grid the road is ONE cell wide, so a buffered mask averages the crown
into whatever is beside it -- in the inter-village stretch, the foredune that
NC-12 runs immediately behind.

Every domain folder also carries clip_domain_<N>.tif at 1 m, the native LiDAR
before the Barrier3D resample. A 3.5 m buffer on that -- a 7 m corridor, about
one carriageway -- gives a within-domain standard deviation of about 0.07 m
island-wide. That is a road surface.

HONESTY NOTE: on THIS alignment the 10 m grid would have given nearly the same
answer -- the two agree to a median of 0.00 m and a max of 0.05 m. The 1.59 m
standard deviation that originally motivated the 1 m clip was measured on the
1984 line, which crosses a relocation scar; the 2004 line does not. So the 1 m
clip is a precaution here rather than a rescue, and the sample it rests on is
~3500 cells per domain against ~35 at 10 m. Both numbers are printed under
INTERNAL CHECKS every run so this stays checkable rather than inherited.

MEAN, NOT MEDIAN
Flat unweighted mean of every valid 1 m cell in the corridor. On this alignment
mean and median differ by 0.005 m for a typical domain and 0.09 m at worst, so
the choice is nearly free; the median is carried in the per-domain CSV so any
domain where they diverge -- a bridge deck, a house, a driveway apron in the
corridor -- is visible rather than silently absorbed.

DATUM -- NOT AMBIGUOUS, DESPITE THE RUNNER
bulldoze() writes road_ele straight into xyz_interior_grid:

    road_ele = road_ele / dz
    new_road_domain = np.zeros(...) + road_ele

and the interior arrays are MHW-RELATIVE, because HAT_dune_topo_extractor.py
subtracts MHW_M = 0.36 before anything else. So road_ele MUST be MHW-relative
metres. There is no reading under which NAVD88 is correct.

The runner's ROAD_ELEVATION = 1.45 is high under EITHER reading -- see the audit
document. This file writes MHW-relative; the per-domain CSV carries NAVD88
alongside so nothing has to be taken on trust.

TWO THINGS THIS FILE DOES NOT CORRECT
1. THE TIME GAP. The DEM is 2009; one run starts in 1984. CASCADE decrements
   road_ele by RSLR every year, so the 1984 run begins with a roadbed that is
   already 25 years of sea-level rise low relative to its own MHW. No
   back-correction is applied -- these are measurements, not reconstructions.

2. THE RELOCATIONS. GIS 9-15 (relocated 1999) and GIS 84-87 (relocated 1989)
   carry the elevation of the POST-relocation alignment in the 1984 run, because
   that is the alignment sampled. In 1984 the road was physically elsewhere in
   those domains. Flagged, not adjusted.

OUTPUTS  (data/hatteras_init/4-mgmt-forcing/road_elevation/)
  RoadElevation.csv          2-row CASCADE file (IDs, m MHW-relative)
  RoadElevation_domains.csv  per-domain stats, both datums, flags
  RoadElevation_audit.md     the tracking document
  HAT_road_elevation.png     alongshore QC

REQUIREMENTS
  geopandas, rasterio, numpy, matplotlib
```

Notes that were in the code:

```text
Native-resolution LiDAR clips, one folder per domain. NOT the 10 m resample.
clip_domain_<N>.tif      1 m, native
resampled_domain_<N>.tif 10 m, the Barrier3D grid
MOVED. Was BARRIER3D_DIR/"2009-raw"/"2009-domain-clipresample", a path the
2026-08-25 restructure removed; the clips then sat under superseded/ and
were lifted out on 2026-08-26 because four scripts read them. Same files,
same per-domain layout - see 1-barrier3d-domains/LINEAGE.md.
```

```text
--- WHICH SURFACE: raw 2009, or 2009 with its holes filled -------------
None    the original 2009 clips, LiDAR holes and all
"<tag>" 0-elevation/{1-gapfill-1m,2-resampled-10m}/<tag>/..._filled.tif

ONLY GIS 78, 79 AND 80 CAN CHANGE. Every other domain on the 2004 alignment
has nodata_frac = 0 -- complete 2009 coverage under the road -- so filling is
a no-op there, and the unaffected neighbours 77 and 81 move by <= 0.002 m.
The three that do change are the same three that drowned on coverage gaps in
2009_v4, and D79 is the reason this switch exists: its unfilled elevation is
a mean over 2287 corridor cells out of 3583, 36% of the corridor missing.

WHY 2008 AND NOT 2014 -- THE ARGUMENT, WHICH NO LONGER DECIDES THE SETTING.
Kept because it is still the right way to think about a road surface, and
because it is the cost of the 2026-08-26 change recorded below: 2008 is not
on disk any more, so this is now an argument about three domains and
centimetres rather than a live choice. Read it, then read that note.

this is a ROAD SURFACE, and the two questions have different answers. The
2008 IOCM survey is one year from the 2009 base, so it measures the same
pavement. The 2014 Post-Sandy survey postdates Hurricane Irene (2011), the
Pea Island breach and the NC-12 reconstruction that followed -- at GIS 78-80
a 2014 surface under the corridor may be a REBUILT road, which is not what
"the 2009 road elevation" means. Topography has no such problem: there the
2014 fill is simply the later and more complete survey of the same barrier.

The choice is about provenance, not magnitude -- the two fills differ by
<= 0.015 m in the corridor. RELOCATION BRACKET and the QC flags re-measure
this every run; FILL_SOURCE is recorded in the audit.

2008 was SUPERSEDED as a TOPOGRAPHY fill (2014 replaced it) while remaining
the right answer for a ROAD SURFACE, for the reason above. It therefore lived
under 0-elevation/superseded/, and this script reaches its product through
scripts/site_layer/hat_elevation_products.py rather than by joining strings - which is
exactly what broke on 2026-08-25, when 2008 was moved there and the
hand-built path stopped resolving. It is now deleted rather than superseded;
resolving through the registry is why that reads as an error instead of as
"no fill available" in 90 domains.

THE PROVENANCE CALL, TAKEN 2026-08-26: FILL_SOURCE = "2009-2014".

This script could not re-run for one day. FILL_SOURCE was "2008_NOAA_IOCM"
and that product is gone: its rasters were never tracked (*.tif is
gitignored) and are not on disk, its registry entry was removed, and the
point-cloud path that built it was removed from HAT_dem_gap_fill.py along
with HAT_laz_ground_classify.py. _elev_product() raised immediately - the
intended failure, not a silent stale read, but a forcing you cannot
regenerate is a forcing you cannot check.

WHY THE BASELINE AND NOT THE 1984 PRODUCT. There are two live products, and
under the road they are NOT interchangeable. Sampling this script's own 2004
alignment in the same 3.5 m corridor on both:

2009-2014-1996 minus 2009-2014, corridor mean:  median +0.222 m
54 of 82 domains move more than 0.05 m
cell counts IDENTICAL in every domain

Identical counts means ALACE REPLACED measured 2009 pavement rather than
filling holes in it, and +0.222 m is not a roadbed: it is the island-wide
1996-vs-2009 survey offset, which mosaic_1984_audit.csv reports per domain at
median +0.255 m (p10 +0.14, p90 +0.33) and HAT_dem_1984_mosaic.py leaves
UNCORRECTED by design ("bias correction OFF, feathering OFF"). Building a
1984 road elevation from the 1984 DEM would push road_ele up ~0.22 m
island-wide and that increment would be the offset, not the road. A higher
road is buried by overwash less often, so it would reach the model.

So: ONE elevation set, on the baseline, for both periods - see the note at
HATTERAS_ROAD_ELEVATION_FILE in hatteras_site_config.py.

WHAT IS LOST BY NOT USING 2008. The 2008 IOCM survey was one year from the
2009 base, so it measured the same pavement, and the header above argues a
2014 surface under GIS 78-80 may be a REBUILT road (post-Irene, post-breach).
That argument still stands and is not resolved by this change - it is
bounded. Only GIS 78, 79 and 80 have any nodata under the 2004 alignment, so
only those three can move at all, and the two fills were measured to differ
by <= 0.015 m in the corridor. Three domains, centimetres, against a forcing
nobody could rebuild. The RELOCATION BRACKET check and the QC flags
re-measure it every run, and FILL_SOURCE is recorded in the audit.

Override from the shell to compare products without editing this file:
HAT_ROAD_ELEV_FILL=2009-2014-1996 python HAT_road_elevation.py
HAT_ROAD_ELEV_FILL=none           python HAT_road_elevation.py
"none" samples the raw 2009 clips, holes and all. Whatever is used is written
into RoadElevation_audit.md, so a comparison run cannot be mistaken for the
shipped one afterwards.
```

```text
The single alignment: the 2008 line, which is the 2004 period's road
(hat_topo_version.ROAD_LINE_FOR_YEAR) and the one contemporaneous with the
2009 DEM everywhere on the island -- see the header. Filed under its true
vintage since 2026-09-15; it was raw_offset/2004/nc12_2004.geojson before.
```

```text
Sampled ONLY for the RELOCATION BRACKET check -- never written to the product.
It exists to test, rather than assume, what the abandoned corridor looks like
in the 2009 DEM. Set to None to skip the check.
```

```text
Half-width of the sampling corridor. 3.5 m -> a 7 m strip, about one
carriageway of NC-12. Widen it and you start averaging in the shoulder and the
dune toe; the sweep printed under INTERNAL CHECKS shows exactly where that
begins to bite.
```

```text
all_touched=False on a 1 m grid: the buffer polygon is already several cells
wide, so touching-cell inclusion only adds edge pixels off the pavement.
```

```text
Domains where the alignment sampled here postdates the 1984 run's road.
Taken from HISTORICAL_ROAD_EVENTS in HAT_hindcast_1984_2024.py.
Informational only -- the value written is the same for both periods.
```

```text
--- map figure ---------------------------------------------------------
The island is ~41 km north-south but each domain is only 2 km across, so an
equal-aspect map of the whole thing is a 6:1 sliver. Cut it into strips laid
side by side, each drawn at TRUE aspect. Distorting the aspect to fit would
make the road look like it changes direction where it does not.
```

```text
Background LiDAR: pale sand at MHW darkening to brown at the dune crests. Low
saturation on purpose -- the road is the subject, the island is the context.
```

```text
_elev_product raises with the known product names if FILL_SOURCE is not
one of them, and again if it is known but not on disk - so the probes
that used to be inlined here cannot go stale when the layout moves.
```

```text
Both resolutions come from the SAME fill, so the 1 m vs 10 m agreement
printed under INTERNAL CHECKS stays a like-for-like comparison.
```

```text
Continuity is an alongshore property, so it can only be flagged once every
domain has a value. A road does not step 0.35 m between two 500 m cells.
```

```text
INTERNAL CHECKS

The ArcGIS elevation_2009 column that this method was originally validated
against no longer exists -- see the audit document. What is left is internal:
does the answer depend on the corridor width, and does the surface we are
obliged to model on agree with the surface we measured.
```

```text
The grade the road would sit at if this reach were not relocated: the
neighbouring un-relocated domains inside the same named reaches.
```

```text
Reach means, so the alongshore pattern reads at a glance rather than
having to be averaged by eye out of 82 points.
```

```text
Bottom-left: the only quadrant with no data in it, and it keeps the
reach labels along the top edge readable.
```

```text
Road geometry per domain, clipped to that domain's box -- the SAME box the
elevation was sampled in, so colour and geometry cannot disagree.
```

```text
Window on the road corridor, widened to MAP_MIN_WIDTH_M so the island
around it stays visible -- where the island is narrow is exactly the
context that explains a low reach.
```

```text
Scale bar bottom-left, north arrow top-left, both on the southern strip
and both kept clear of the GIS labels on the right.
```

<details><summary>Function notes (the original docstrings)</summary>

**`find_clips_flat()`**

```text
domain id -> raster path, for the gap-fill tree.

The clipresample tree nests one folder per domain; the gap-fill tree is
flat and carries the id in the filename instead. Same contract, different
layout, so the caller does not have to care which surface it asked for.
```

**`corridor_values()`**

```text
Elevations under the buffered road line, and the fraction of the corridor
that fell on NoData.

The buffer is applied in the RASTER's CRS, converted through its linear
units factor -- the geojson is EPSG:2264 (US survey FEET) and the clips are
UTM 18N (metres), so a raw 3.5 would be 1.07 m if taken literally in the
source CRS.

`z` and `bad` are passed in so a caller sweeping several buffer widths reads
each raster exactly once.
```

**`stats_for()`**

```text
MEAN is the product. Median, percentiles and sigma are diagnostics that ride
along in the per-domain CSV so a domain where the mean is being dragged by a
structure in the corridor announces itself.
```

**`relocation_bracket()`**

```text
What is actually under the 1984 alignment where the road was relocated?

The case for using the 2004 line everywhere rests on a claim about the
abandoned corridor, and a claim in a docstring is not evidence. This
measures it: same surface, same corridor width, the OTHER line, in the 11
domains where the two disagree.

Nothing here is ever written to RoadElevation.csv. It exists so the choice
can be defended with a number that is re-derived on every run.
```

**`resample_check()`**

```text
Same corridor, same buffer, on the 10 m Barrier3D grid instead of the 1 m
clip. This is not a validation -- the 1 m answer is the better one -- it
quantifies how much the grid the model actually runs on would have cost us.
```

**`write_cascade_csv()`**

```text
2-row CASCADE file: row 0 = GIS IDs, row 1 = elevation in m MHW-RELATIVE.

Refuses on a gap, because the runner fills its per-domain arrays BY POSITION
after dropping the ID row -- one missing domain shifts every domain north of
it and nothing reports it.
```

**`map_figure()`**

```text
Where the road is high and where it is low, in real geography.

The alongshore profile answers "how much"; it cannot answer "where", because
a domain index is not a place. This draws the road on the island it sits on,
coloured by the same numbers, over the 10 m LiDAR for context -- so a low
reach can be read against the island being narrow there.

Drawn in strips at TRUE aspect ratio. See MAP_STRIPS.
```

</details>

