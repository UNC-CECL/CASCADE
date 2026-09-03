# 2009-2014-1996-duneline — where the two dune lines sit

The 1984-start DEM with the **1984 and 1997 digitized dune lines** drawn on it,
and the cross-shore distance between those two lines, per domain.

**This folder holds no elevation raster and forks nothing.** The 178 MB of 1 m
tifs in `../2009-2014-1996/` are read in place, so this and the product it
describes cannot drift apart on disk. It is deliberately not registered in
`scripts/hat_elevation_products.py` — that resolves products with
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
`../../1-barrier3d-domains/raw-duneline-geojson/README.md`, which was derived a
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

| file | what |
|---|---|
| `duneline_offset_by_domain.csv` | 90 rows. Median, quartiles, min/max, mean, sd, cells, the nearest-point check, and the median easting of each line |
| `figures/HAT_duneline_offset_ribbon.png` | **the one to read.** Both lines against a smoothed midline, band filled by sign, at full 1 m alongshore resolution |
| `figures/HAT_duneline_offset_zooms.png` | the same two lines on the DEM at true scale, three reaches |
| `figures/HAT_duneline_offset_zoom_83_87.png` | one extra reach, GIS 83–87 — the five domains around 85 |
| `figures/HAT_duneline_offset_island.png` | the whole island in three panels — a locator, not a measurement |
| `figures/HAT_duneline_offset_bydomain.png` | the offset per domain, with each domain's IQR |

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
on by default when `1-barrier3d-domains/1984-start/row-insert-scope/` is on
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
