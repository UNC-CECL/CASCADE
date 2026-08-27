# 2009-2014 - the baseline DEM

**2009 USACE**, with its nodata filled from the **2014 NOAA Post-Sandy** survey.

Built by `scripts/input_prep/0-elevation/2-produce/HAT_dem_gap_fill.py`, then
`HAT_dem_resample_clip.py` (no arguments - this is the default product).

## What it is for

Everything except the 1984 start. It is the DEM behind the `2009_v5`
dune/topography run, and it currently feeds **both** hindcast periods - which
is why this folder is not named for one of them.

## The rule that matters

The fill **only ever writes into nodata**. A cell the 2009 survey measured is
never changed. The four selection rules - coverage, connectivity, the -2.64 m
floor, and no vertical reconciliation at all - are documented with their
reasoning in
[`scripts/input_prep/0-elevation/README.md`](../../../../scripts/input_prep/0-elevation/README.md).

Result: 43.7% of nodata recovered, 25,591,292 cells at 1 m, 254,760 at 10 m.

## Survey codes in `*_survey.tif`

| code | meaning |
|---|---|
| 2009 | measured by the 2009 USACE DEM |
| 2014 | filled from 2014 NOAA Post-Sandy |
| 0 | no survey saw this cell |

## The seam is real and deliberately uncorrected

Where fill meets measured 2009 ground there is a step. Bias correction and
feathering are both **off**, so a filled cell is the 2014 measurement unchanged.
Both bias estimates are computed and written to `gapfill_audit.csv` every run,
so "we chose not to correct" stays checkable. The reasoning is in
[`../FIGURES.md`](../FIGURES.md).

## Figures

| file | shows |
|---|---|
| `HAT_gapfill_2009-2014_island.png` | whole island: 2009 only / 2009+2014 / survey source |
| `HAT_gapfill_2009-2014_domains_78_80.png` | zoom on the road-drowning domains |
| `HAT_gapfill_2009-2014_roads_78_80.png` | that zoom with BOTH NC-12 alignments |
| `HAT_gapfill_2009-2014_roads_8_15.png` | domains 8-15 with BOTH alignments |
| `HAT_gapfill_2009-2014_roads_82_88.png` | domains 82-88 with BOTH alignments |

Both alignments are drawn: **2004 solid black, 1984 white-dashed on top**, so
where they coincide you see a black line with white dashes and where they
diverge each is legible alone. `2009-2014-1996` carries the same 8-15 view over
the 1984-start DEM - the two are meant to be read side by side, with the roads
held constant and the DEM the only thing that changes.

Add or change a zoom in `ROAD_ZOOMS` at the top of
`3-figures/HAT_plot_gapfill.py`: `(domain ids, which NC-12 years, filename
slug, title)`.
