# Elevation figures - design decisions

Shared by every product under `0-elevation/`. Each product keeps its own
rendered figures in its own `figures/` folder.

Review figures for the gap fill, and the design decisions behind them.

The pipeline that produces them - datasets, the chain, fill rules, how to change
fill source - is documented next to the scripts, in
[`scripts/input_prep/0-elevation/README.md`](../../../scripts/input_prep/0-elevation/README.md).
Kept there rather than duplicated here so the two cannot drift apart.

## Figures

| file | shows |
|---|---|
| `2009-2014/figures/HAT_gapfill_2009-2014_island.png` | whole island: 2009 only / 2009+2014 / survey source, with NC-12 |
| `2009-2014/figures/HAT_gapfill_2009-2014_domains_78_80.png` | same panels, zoomed on the road-drowning domains |
| `2009-2014/figures/HAT_gapfill_2009-2014_roads_78_80.png` | the zoom with both NC-12 alignments overlaid |
| `2009-2014-1996/figures/HAT_mosaic1984_*.png` | the 1984-start DEM: three sources, and the window shift |

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

```
python scripts/input_prep/0-elevation/3-figures/HAT_plot_gapfill.py
python scripts/input_prep/0-elevation/3-figures/HAT_plot_1984_mosaic.py
```

`HAT_plot_gapfill.py` reads `<PRODUCT>/2-resampled-10m/` and writes all three
PNGs into `<PRODUCT>/figures/`. Paths come from
`scripts/hat_elevation_products.py`, so a product that is not on disk is a
loud error rather than an empty glob.

**The 2008 comparison figures are not reproducible (since 2026-08-26).** They
were regenerated with `--source 2008_NOAA_IOCM`; that option is gone, because
everything behind it is: the 2008 rasters are off disk (`*.tif` is gitignored,
so they were never in the repo), the product entry has been removed from the
resolver, and the point-cloud path that built them has been removed from
`HAT_dem_gap_fill.py` together with `HAT_laz_ground_classify.py`. The numbers
that mattered survive in this file and in the 0-elevation README; the PNGs do
not. Re-running that comparison now means re-gridding the point cloud to a DEM
outside this repo and registering it as a new product.

Tag, fill year and long label live together in the script's `SOURCES` dict, so a
re-render cannot put one source's label on another's data - hand-editing the
three constants separately already printed "cells the 2014 DEM measured" on the
2008 figures once.

To regenerate the full chain, not just the figures, see the
[pipeline README](../../../../scripts/input_prep/0-elevation/README.md).
