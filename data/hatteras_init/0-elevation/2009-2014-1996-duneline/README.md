# 2009-2014-1996-duneline — where the two dune lines sit

The 1984-start DEM with the **1984 and 1997 digitized dune lines** drawn on it,
and the cross-shore distance between those two lines, per domain.

**This folder holds no elevation raster and forks nothing.** The 178 MB of 1 m
tifs in `../2009-2014-1996/` are read in place, so this and the product it
describes cannot drift apart on disk. It is deliberately not registered in
`scripts/site_layer/hat_elevation_products.py` — that resolves products with
`1-gapfill-1m` and `2-resampled-10m` stages, and this has neither.

Built by `scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py`.

## What is measured

One thing: **how far apart the 1984 and 1997 lines are.** Nothing here says
whether the 1996 ALACE swath reaches either of them — see *The question this
does not answer* below.

Every domain box is axis-aligned, 2000 m in easting by 500 m in northing, and
this pipeline's `OCEAN_LOC` is `"right"` — the Atlantic is at increasing
easting. So easting is cross-shore, northing is alongshore, and the separation
between two shore-parallel lines is a difference in easting at a shared
northing. The script **checks** the 2000 × 500 box shape rather than assuming
it, and raises if it ever stops holding.

    offset_m = x_1984 − x_1997     at the same northing

    POSITIVE  =  the 1984 line lies SEAWARD of the 1997 line
                 — the sign 13 years of erosion predicts.

Sampled every 1 m along each box's northing axis, so each domain contributes
500 independent measurements and the reported number is their median with the
quartiles beside it. **All 90 domains returned 500/500 samples for both lines** —
neither line has a gap anywhere in the modelled reach.

## The result

| | m | 10 m cells |
|---|---:|---:|
| median over 90 domains | **+1.2** | +0.12 |
| p25 / p75 | −10.0 / +14.3 | −1.0 / +1.4 |
| min (domain 64) | **−58.9** | −5.9 |
| max (domain 80) | **+70.2** | +7.0 |
| mean absolute | 18.2 | 1.8 |
| RMS | 24.8 | 2.5 |

**47 of 90 domains positive, 43 negative.** The island-wide median of +1.2 m
corroborates the +0.8 m "date" term in
`../../2-brie-offset/dunelines/README.md`, which was derived a
different way — but that agreement is the least informative thing in the table.
The two lines disagree by **at least one Barrier3D cell in 52 of 90 domains**
and by three or more cells in 17:

    3, 5, 16, 49, 50, 63, 64, 65, 66, 67, 73, 79, 80, 81, 82, 84, 85

and the sign is not consistent alongshore. Domains 63–67 have the 1984 line
40–59 m **landward** of the 1997 line — the wrong direction for erosion, across
five neighbouring domains, with a within-domain standard deviation of only
5–10 m, so it is not noise. Domains 79–85 run the other way, up to +70 m.

An island-wide mean near zero built out of ±60 m local disagreement is not
evidence that the two lines describe the same feature. It is consistent with
that, and equally consistent with digitizing differences that happen to cancel.
Nothing here can separate the two, because the 1984 file carries **no metadata
at all** — no `feature_type`, no `method`, no editor — while 1997 records
*"digitized from light/dark elevation break (no DEM available)"*. That caveat is
in the raw-duneline README and it is the limiting factor on reading these
numbers, not the measurement.

## What the zooms show — an observation, not a measurement

In `HAT_duneline_offset_zooms.png` the two excursion reaches do not look like
one line displaced from the other. They look like the two lines **straddling
the same dune ridge from opposite sides**, and swapping which side they take:

* **62–68**, where the offset is −40 to −59 m, 1984 runs along the landward
  flank of the bright ridge and 1997 along its seaward flank.
* **78–85**, where the offset is +23 to +70 m, the roles are reversed.
* **17–21**, the control, has both lines on the crest and nearly coincident.

Real cross-shore dune movement between 1984 and 1997 would be systematic in one
direction. Swapping sides is what a **definitional** difference looks like —
which edge of the light/dark break a digitizer picked — and 1997's recorded
method is exactly *"digitized from light/dark elevation break (no DEM
available)"*.

This is read off a rendered figure against a 2009-backdune DEM, not measured,
and it is written down here as the thing to test next rather than as a result.
Testing it properly means comparing each line to the dune crest in a surface of
its own era, which for 1984 does not exist.

## The orientation check

`nearest_med_m` is the plain nearest-point distance from each 1984 sample to
the 1997 line — no axis, no sign, no assumption about which way the ocean is.
It is reported next to the easting difference as a check on the frame.
`offset_over_nearest` is their ratio; where the island runs obliquely to the
grid the two must diverge, by one over the cosine of the obliquity.

**Worst ratio island-wide is 1.11.** The easting-based number is within 11% of
the orientation-free one everywhere, so the axis-aligned frame is sound and the
signed offset can be read directly.

## Files

Every figure is drawn to one house style (`apply_style()` in the script): Arial, thin dark-grey axes, ColorBrewer red/blue for the two lines and for the sign of the offset, panel letters, a north arrow and scale bar on maps without coordinate ticks, and no in-image titles or footnote paragraphs. The words that used to be on the figures are in `figures/CAPTIONS.md`.

`figures/` is sorted by what a figure IS (2026-09-08): `island/` holds the
whole-island maps, `detail/` the true-scale crops, `offset/` the two readings
that are not maps. `fig_path()` in the script is the only way a figure name
becomes a path, so nothing can land at the folder root. Every map in the folder
is drawn in one style — grey relief from the 1 m gap-filled DEM, the 1984 line
solid red and the 1997 line solid blue, a scale bar in place of coordinate
ticks; the two terrain-coloured figures the folder used to carry (the locator
and the reach zooms) are retired and redrawn respectively.

| file | what |
|---|---|
| `duneline_offset_by_domain.csv` | 90 rows. Median, quartiles, min/max, mean, sd, cells, the nearest-point check, and the median easting of each line |
| `figures/detail/HAT_duneline_offset_simple.png` | **the one to look at.** The two lines on grey relief, four two-domain pairs, tight crop, no elevation values |
| `figures/island/HAT_duneline_offset_simple_island.png` | the same two lines over the whole island, with each domain's measured offset as a bar aligned to the map |
| `figures/island/HAT_duneline_offset_simple_island_mean.png` | the same map-and-bar figure with the per-domain MEAN on the bars instead of the median; `--simple --simple-stat mean` |
| `figures/CAPTIONS.md` | a caption per figure with the numbers filled from the table; the figures carry no title sentences or footnotes, so use these under them. `--captions` rewrites it alone |
| `figures/island/HAT_duneline_offset_lines_island.png` | the whole island as maps only, no bar strips: nine ~5 km panels, each cropped at equal aspect to the strip the two lines occupy, so the offset is visible on the map itself. `--lines-island` renders it alone |
| `figures/island/HAT_duneline_offset_lines_island_3panel.png` | the same, three panels of 30 domains to match the simple_island layout; drawn by the main run since 2026-09-08 (or `--lines-island --lines-per-panel 30 --lines-island-out <this path>`) |
| `figures/offset/HAT_duneline_offset_ribbon.png` | **the one to read.** Both lines against a smoothed midline, band filled by sign, at full 1 m alongshore resolution |
| `figures/detail/HAT_duneline_offset_zooms.png` | the same two lines at true scale on three reaches of five to eight domains, the wider ±300 m crop; since 2026-09-08 drawn in the simple style (grey relief, both lines solid) through `fig_zooms_simple`, so nothing in the folder carries the terrain colour ramp any more |
| `figures/detail/HAT_duneline_offset_zoom_83_87.png` | one extra reach, GIS 83–87 — the five domains around 85, each labelled with the rows the 1984 footprint adds; also copied to `2-domain-reconstruction-1984/figures/1-measurement/rows-added/` |
| ~~`figures/HAT_duneline_offset_island.png`~~ | retired 2026-09-08: the terrain-coloured locator with the domain boxes, superseded by `island/HAT_duneline_offset_lines_island_3panel.png`, which shows the same boxes and lines on grey relief. `fig_island` stays in the script, uncalled |
| `figures/offset/HAT_duneline_offset_bydomain.png` | the offset per domain, with each domain's IQR |

### The line key

Every figure here draws the same two lines the same way: **1984 solid red
(`#d62728`), 1997 solid blue (`#1f77b4`)**, both with a white casing, 1984 on
top. Colour is the only thing that separates them. Until 2026-09-03 the 1997
line was dashed on the DEM figures and solid on the simple ones, so one folder
carried two keys for one pair of lines; `LINE_STYLE` in the script is now the
single key and `SIMPLE_LINE_STYLE` is a copy of it.

What still varies between figures is line WIDTH, through
`draw_lines(scale=...)`: `LINE_SCALE_ISLAND` for the two whole-island figures,
`LINE_SCALE_DETAIL` for the two detail figures and the zooms,
`LINE_SCALE_RIBBON` for the trace. They are paired constants, not per-figure
numbers, so a change cannot land on one of a pair and not the other. A 46 km
locator and a 300 m crop cannot carry the same line weight.

**On the two whole-island maps this means the 1997 line is mostly invisible** —
solid red is drawn over solid blue and at that scale they coincide. That is the
same scale limit the locator already had, and it is why
`HAT_duneline_offset_simple_island.png` carries a bar strip. Read the offset
off the bar or off `HAT_duneline_offset_ribbon.png`, never off either island
map.

### The simple figure

`HAT_duneline_offset_simple.png` is the zooms with everything that is not the
two lines taken out: no terrain colour ramp, no colourbar, no elevation values
and no coordinate ticks. What is behind them is a greyscale hillshade, so the
dune ridge is still there to see which side of it each line takes, but nothing
about its height can be read off, and nothing is labelled with a number.

It is built by the same function family as the zooms and shares their line
key, draw order and equal aspect, so the two cannot disagree about which line
is which or which is on top.

Three differences from `HAT_duneline_offset_zooms.png` are worth knowing:

**Two domains per panel, not five to eight.** The island runs about 7° oblique
to the UTM grid, so a dune line drifts ~130 m in easting for every km of
alongshore. A crop tight enough to make a 50 m offset obvious therefore cannot
hold a whole reach — over 1.5 km the two lines sweep 190–330 m in easting and
leave the frame. Two neighbours sweep 150–223 m, which fits inside ±150 m with
margin. So each panel is the **pair that carries its reach's extreme**, not the
reach. Four of them, south to north:

| pair | what | offset |
|---|---|---:|
| 3–4 | the far south | +40, +29 m |
| 19–20 | the quietest pair on the island | −2, −7 m |
| 63–64 | 1984 landward of 1997 | −51, −59 m |
| 79–80 | 1984 seaward of 1997 | +46, +70 m |

The same half-width, 150 m, is used for all four, because the control only
works if it is drawn at exactly the scale of the others. Each pair was checked
against that before being chosen — 3–4 needs 115 m, 19–20 needs 88, 63–64
needs 112, 79–80 needs 75. **4–5 carries the south's single largest offset**
(+62 m at domain 5) and is deliberately not used: its lines need 196 m and
would have forced a wider crop on all four panels.

**The backdrop is the 1 m gapfilled tile, not the 10 m mosaic.** Across a 300 m
frame the 10 m product is 30 cells wide and draws the dune as a staircase. The
1 m tiles carry **no CRS tag**, so the script does not trust them: each tile's
bounds are checked against its domain box and a disagreement over 1.5 m is
fatal. That check is the only thing tying the raster to the dune lines.

**The shading is vertically exaggerated ×2.2 and smoothed over 3 m; the map is
not exaggerated at all.** These are two different things and the figure says so
in its footer. The dune is ~5 m of relief over ~50 m cross-shore, which at 1:1
shading is nearly flat grey, and 1 m lidar over a vegetated backdune is speckly
enough at the cell scale that exaggerating the slope buries the ridge in noise.
Both parameters touch a smoothed **copy** used for shading only. The map plane
is equal aspect and 1:1, the elevation array is untouched, and nothing on this
figure is measured off the backdrop — the numbers all come from
`duneline_offset_by_domain.csv`.

A scale bar replaces the coordinate ticks: 50 m, which is five Barrier3D cells.

### The panel titles are derived, not written

Each detail panel is titled in two lines:

```
Domains 63-64 . Wimble Shoals
1984 landward by 51-59 m (5-6 cells)
```

Only the *role* of a pair is written by hand, in the third field of
`SIMPLE_REACHES` - and only 19-20 has one, `"control"`. Everything else is
generated:

* **the place** from `HATTERAS_ANNOTATIONS`, most specific first - a community
  containing the pair (naming the village centre it sits on, where there is
  one), else a named shoal zone, else the gap between the two nearest
  communities. So 79-80 is "Tri-Village, at Rodanthe", 63-64 is "Wimble
  Shoals", 19-20 is "Buxton-Avon" and 3-4 is "south of Buxton". No place name
  in this figure was invented for it.
* **the direction** from the SIGN of the measured medians, not from a word
  typed next to them. A panel cannot read SEAWARD over numbers that are
  negative.
* **the magnitude and the cell count** from the same table the in-panel labels
  read, with `round(offset / 10 m)` - the same rounding the row-insert scope
  uses.

The titles used to carry hand-written descriptions such as
`"1984 line LANDWARD of 1997"`, which restated in words what the numbers
underneath already said and could drift from them if the table were remeasured.

### The whole island

`HAT_duneline_offset_simple_island.png` is the same two solid lines and the
same grey relief over all 90 domains, in three columns of a third of the island
each. Every column is **two axes sharing a northing axis**:

* **left, a true map.** Equal aspect, both lines, the four detail pairs
  labelled, and the communities bracketed on the seaward margin. It is a
  *locator*. It does not show the offset and it cannot: 46 km of island against
  a 70 m offset is under half a line width, so on the map the two lines lie on
  top of each other nearly everywhere. That is a property of the scale, not of
  the lines.
* **right, the measured offset as a bar**, one bar per domain, aligned to the
  map row for row, red where 1984 lies seaward and blue where it lies landward,
  with the ±10 m Barrier3D cell marked. Same medians the detail panels print,
  off the same CSV.

Read together they answer the two halves of the question — the bar says where
along the island the lines disagree and by how much, the map says what that
part of the island looks like and where each detail panel is cut from. The bar
exists **because** the map cannot carry the offset, and the figure's footer
says so rather than leaving the coincident lines to be read as agreement.

**The detail pairs are no longer boxed.** At this scale a two-domain box is
2 km of a 2 km-wide island, so the rectangle enclosed the whole width and read
as a feature of the island rather than a crop mark. The bold pair label on the
map and the grey band on the bar carry the same information without drawing a
rectangle over the only two lines the map has.

**The place names are not defined here.** Communities, village centres and the
two end labels are read from `HATTERAS_ANNOTATIONS` in
`scripts/site_layer/hatteras_site_config.py` — the same object the shoreline-rate figures
annotate from, and the same spans the model consumes as
`HATTERAS_COMMUNITY_ZONES`. Buxton is GIS 7–8, Avon 21–31, Tri-Village 68–83
(Salvo 69, Waves 74, Rodanthe 80), all inclusive GIS ids in the 1 = south frame
this figure already uses. A town that moves in the config moves here, and this
figure cannot disagree with the rest of the repo about where Avon is.

Two placement rules are worth knowing. Avon straddles the 1–30/31–60 panel
break, so its **bracket is drawn in both panels** — the community really does
continue past the break — but its **name is drawn only on the panel holding
most of it**, or the one-domain sliver at the foot of panel 2 reads as a second
Avon. And the communities are bracketed in the ocean margin rather than washed
across the panel: a translucent band would sit on the island and on both dune
lines, and its colour (`#90AFC5`) is a blue close enough to the 1997 line's to
be read as belonging to it.

**Structures and the alongshore ruler.** Both piers (Avon, GIS 26; Rodanthe,
GIS 79) and the Buxton groin (GIS 5.5, the boundary between domains 5 and 6)
are drawn as short marks running seaward off the 1984 line, in the config's own
`color_pier` and `color_groin`. Two things about them:

* **The length is a drawing constant** (`STRUCTURE_LEN_M`), not a measurement.
  No structure in this repo has a surveyed length, and a 46 km panel could not
  resolve the difference between 200 m and 400 m of pier anyway. Read them as
  positions, not extents.
* **They are named in the legend, not on the map.** At this scale a label
  beside the Rodanthe pier lands on the Rodanthe village tick and on the 79-80
  detail label -- three labels inside one 500 m domain, which no amount of
  nudging fixes.

They are drawn perpendicular to a pair of shore-parallel lines, which is what
keeps `color_pier` (`#1565C0`) from being read as the 1997 line despite the two
being close blues.

The ruler on the left of each bar is **km north of the south end of domain 1**,
the same origin and direction `HAT_duneline_offset_ribbon.png` uses on its
x-axis, so the two figures can be read against each other. It sits on the bar
and not on the map deliberately: the map is pinned to equal aspect and its
column width is set to its own data aspect, so hanging tick labels off it would
letterbox it and break the row-for-row alignment that is the entire reason the
bar sits beside it.

The piers, the groin and the ruler are on the locator only. The ribbon does
not carry them yet -- it has a domain x-axis and could take
`cascade_pipeline.annotations.add_geographic_annotations` directly.

One implementation note that is not cosmetic: **the map column widths are
computed, not chosen.** Every axes in a one-row figure gets the same height, so
a map pinned to equal aspect in a column wider than its own data aspect is
letterboxed — it shrinks vertically and stops lining up with the bar beside it,
which destroys the only thing the pairing is for. Each map column is therefore
set to exactly its own x-span/y-span. The island really is 5.7 km wide at the
south end and 4.1 km at the north, so the three columns are genuinely different
widths.

### Rendering either one

```
python scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py \
    --simple --simple-span 84-85 --simple-out <path>
```

`--simple` renders both figures. It reads `duneline_offset_by_domain.csv`
rather than re-measuring, so it takes seconds. `--simple-halfwidth` widens the
crop if a chosen pair's lines run outside ±150 m, `--no-island` skips the
locator, and `--simple-island-out` sets the locator's path.

### Seeing a 20 m difference on a 46 km island

The island map cannot show this and no styling fixes that. At equal aspect the
island is 46 km long and about 2 km wide; a panel 3 km across renders two lines
18 m apart as one line. Splitting the map into thirds helps a little and is
worth having as a locator, but the separation is resolved two other ways:

**The ribbon** drops the map. Plotted as raw easting the two lines sweep 6.5 km
across the island's curve, which is what drowns the signal — so a **smoothed
midline of the two, boxcar over 2 km alongshore**, is subtracted and what
differs between them is left. The baseline is a drawing device carrying no
claim: it is symmetric in the two lines, so it cannot move one relative to the
other, and the filled band's width is the offset exactly. What the window
length does control is how much of each line's own sinuosity survives — shorter
flattens both toward the axis, longer lets shared meanders back in.

**The zooms** stay on the DEM at equal aspect with nothing exaggerated. What
makes the offset visible there is the **cross-shore crop**: 300 m either side of
the local line position instead of the full 2000 m domain box, so 50 m of offset
is about a twelfth of the frame rather than a fortieth. Three reaches, picked
from the table rather than by eye — 62–68 (largest sustained negative run),
78–85 (largest positive), and 17–21, **the quietest five domains on the island**.
The control is not optional: without it every figure of this kind reads as a
discrepancy, and there is no way to see what agreement looks like at the same
scale.

Both lines are **clipped to the domain footprint for drawing only** — each runs
~16 km past domain 90 unclipped, which would draw dune line where there is no
model domain. The measurement uses the unclipped lines, so a 1984 sample near a
box edge can still find its true nearest 1997 point.

Any other reach can be rendered without re-measuring:

```
python scripts/input_prep/0-elevation/3-figures/HAT_plot_duneline_offset.py \
    --zoom 83-87 --zoom-out <path>
```

`--zoom` reads `duneline_offset_by_domain.csv` rather than re-running the 1 m
sampling, so it takes seconds, and it goes through the same code path as the
three standard reaches — the crop rule, the line styling and the equal aspect
cannot quietly differ. `HAT_duneline_offset_zoom_83_87.png` is that command.
It overlaps the 78–85 panel at domains 83–85 and extends it to 86–87.

Add `--no-row-labels` to suppress the Barrier3D row-count annotation, which is
on by default when `1-barrier3d-domains/1984-start/2-domain-reconstruction-1984/` is on
disk. The copies in THIS folder are un-annotated: row counts belong to the
insert, not to the dune-line measurement.

## The question this does NOT answer

Whether the 1996 ALACE swath actually reaches the 1984 dune. That was the
original question and it was set aside deliberately, because it does not have a
clean answer:

`scripts/input_prep/0-elevation/2-produce/HAT_dem_duneline_coverage.py` was
written for it and works — its recomputation of the mosaic's stage 1 matches
the shipped `clip_domain_*_survey.tif` **cell for cell, 0 mismatches** — but the
verdict it produces is entirely a function of one free parameter. Holes inside
the ALACE swath have a median of 12 m and a p90 of 51 m, so "how far landward
does 1996 reach" depends on how large a hole you are willing to walk across:

| hole tolerance | domains with the 1984 dune inside the swath | median gap |
|---:|---:|---:|
| 0 m | 14 / 90 | +103 m |
| 9 m | 23 / 90 | +88 m |
| 15 m | 33 / 90 | +73 m |
| 30 m | 63 / 90 | −185 m |
| 60 m | 85 / 90 | −237 m |

There is no plateau. Between 15 and 30 m the walk starts bridging the swath's
own landward edge into detached backdune patches and the island-wide sign
flips. Any single tolerance is a chosen answer, not a measured one.

**That script has not been run over all 90 domains and writes nothing into this
folder.** Its `--domains` outputs were deleted rather than left to be mistaken
for a result. Resolving it needs a decision about what the flag rests on; until
then, this folder contains the line-to-line distance and nothing more.

## Status

Written 2026-09-03. Current against `2009-2014-1996` as rebuilt 2026-08-26.
`HAT_duneline_offset_simple.png` and `HAT_duneline_offset_simple_island.png`
added 2026-09-03; both are redrawings of the existing table, and no measurement
changed when they were added.

All seven figures redrawn 2026-09-03 on the single line key described above;
the island locator gained the community brackets and the two end labels, and
the detail panels gained derived titles. The locator also gained the two
piers, the Buxton groin and an alongshore km ruler. Drawing only — every number in this
README still holds: median +1.2 m, range −58.9 to +70.2 m, 47 of 90 positive.

`duneline_offset_by_domain.csv` was **deliberately left at its committed
version**. Re-running `main()` remeasures, and on a machine where
`D:\Hatteras_GIS\domains.geojson` is not mounted the domain boxes are rebuilt
from the 90 resampled rasters instead — which shifts each box by centimetres
and rewrites all 90 rows. The rewrite is noise, not a result: max change in
`offset_med_m` 0.12 m, mean 0.014 m, and every summary statistic identical to
four significant figures. If you remeasure on a machine with the GIS drive
mounted, commit that CSV; a rewrite produced without it should be reverted.
