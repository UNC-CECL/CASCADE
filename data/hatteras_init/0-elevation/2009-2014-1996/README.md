# 2009-2014-1996 - the 1984-start DEM

The `2009-2014` product, **plus the 1996 NOAA/NASA ALACE survey overwriting
measured ground wherever ALACE has data**.

    wherever 1996 has data :   1996  >  2009  >  2014
    everywhere else        :           2009  >  2014

Built by `HAT_dem_1984_mosaic.py`, then
`HAT_dem_resample_clip.py --product 2009-2014-1996`.

## Read this before using it

**There is no road boundary.** Until 2026-08-26 the 1996 override was confined
to the ocean side of the 1984 NC-12 alignment. It no longer is. The landward
limit is the ALACE swath's own edge, which is roughly the landward base of the
dunes. Both NC-12 alignments are still rasterized every run and reported per
domain, but they gate nothing.

**The graft boundary is the dune toe, not the road.** ALACE surveyed "from the
low water line to the landward base of the sand dunes." So this is a 1996 beach
and foredune on a 2009 backdune, and the seam lands where
`HAT_dune_topo_extractor.py` picks its dune windows. Read a dune-crest number
near that line knowing it sits beside a survey boundary.

**Nothing in it is 1984.** 1996 beach, 2009 interior, 2014 sound.

**The window origin moves** relative to the plain 2009-2014 product. 36 of 90
domains move at least one 10 m Barrier3D cell seaward; domains 76-81 move
47-72 m. Per-domain numbers are in `1-gapfill-1m/mosaic_1984_audit.csv`
(`start_beach_before_m` / `_after_m` / `_shift_m`).
**A new dune pick set is required** - the `2009_v5` windows were picked against
a different origin.

**`shoreline_offset` is NOT affected.** That comes from
`2-brie-offset/hindcast_1984/`, is measured independently of any DEM, and is
passed straight to `Cascade()`.

## Why the road boundary was dropped, 2026-08-26

The original design used "ocean side of NC-12" as the landward limit, on the
reasoning that a rule you can state beats a survey artefact. Deciding WHICH
NC-12 - the 1984 or the 2004 alignment - forced the question of what the
boundary was doing, and the answer was: very little.

All three options were measured across all 90 domains by running stage 1 three
times with only the mask changed, every other rule identical:

| boundary | 1996 cells written | new land | overwrote measured |
|---|---:|---:|---:|
| 1984 line | 11,983,935 | 2,084,834 | 9,899,101 |
| 2004 line | 12,251,459 | 2,086,884 | 10,164,575 |
| **none** (shipped) | **14,707,324** | **2,320,397** | **12,386,927** |

Read the middle column. Moving 1984 to 2004 buys 267,524 cells of which only
**2,050** are ground no other survey saw - 99.2% of that gain is 1996
replacing a 2009 or 2014 measurement in the backdune, the part of the ALACE
swath with the worst coverage. Dropping the boundary buys 2,723,389 cells,
235,563 of them new land: **115x the recovered ground of the 2004 line**.

Both objections to dropping it failed on measurement:

**It does not reach the sound.** Cross-shore reach of what 1996 writes, metres
from the ocean edge, median over profiles, against the island's own extent:

| domain | island extent | road 84 / 04 | 1996 reach, no boundary |
|---:|---:|---:|---:|
| 11 | 1999 | 416 / 499 | 560 |
| 60 | 1274 | 621 / 621 | 979 |
| 77 | 1627 | 717 / 717 | 429 |
| 85 | 1782 | 185 / 305 | 432 |
| 90 | 1999 | 382 / 382 | 441 |

ALACE stops itself at 429-979 m against an island 1274-1999 m wide. The road
was never what held it back. Domain 77 is the proof: reach is 429 m under every
option while the road sits at 717 m - there the boundary was not even binding.

**It admits no extra junk.** Island-wide maximum stays 12.00 m under all three
options. The ceiling, not the mask, is what holds the ALACE return tail.

**And the choice barely touches the window origin.** `start_beach` moves in
ZERO domains between the 1984 and 2004 masks, and in only 6 between 1984 and
none - all of them domains 1-7, which previously received no 1996 at all.
Domains 8-90 have the same window origin they had under the old boundary.

### What was given up

A stated rule was traded for a survey footprint. "Ocean side of NC-12" is
defensible in a methods section; "wherever ALACE flew" is a survey artefact
standing in for a landform boundary. The honest description of this product is
that its 1996 extent is the ALACE swath. The audit reports
`n1996_landward_of_1984_road` and `n1996_landward_of_2004_road` per domain so
the band the old boundary discarded is a number rather than a claim:
**2,723,389 cells, 18.5% of everything 1996 wrote.**

## Domains 1-7 - resolved

Both NC-12 exports end at domain 8, so under the old boundary domains 1-7 got
no 1996 and carried a 2009 beach while 8-90 carried a 1996 one - an alongshore
discontinuity at the 7/8 seam. The documented alternative was to invent a
boundary by extending domain 8's alignment 3.3 km south across Hatteras
village.

Removing the boundary dissolves it instead of trading it. 1-7 now receive 1996
on the same terms as everywhere else - 77k to 175k cells each, 6k to 60k of it
ground no other survey saw - and nothing is invented. The seam is gone.

They are still flagged `road_line=False` in the audit. That flag is now purely
descriptive: these domains have no NC-12 export, not a different treatment.

## Survey codes in `*_survey.tif`

| code | meaning |
|---|---|
| 1996 | ALACE |
| 2009 | measured by the 2009 USACE DEM |
| 2014 | filled from 2014 NOAA Post-Sandy |
| 0 | no survey saw this cell |

## Two guards ALACE needed that the baseline does not

Both were found by failure. Both are still in the script, and both count what
they reject per domain in the audit. No value is clamped or shifted - a
rejected cell falls through to 2009, then 2014.

| guard | value | what it stops |
|---|---|---|
| split floor, into a gap | -2.64 m NAVD88 | nothing - this is the normal admission, where **no other survey** saw the cell |
| split floor, to replace | MHW, 0.36 m | a wet swash return displacing dry measured beach |
| ceiling | 12 m NAVD88 | the uncorrected-return tail, which reaches **+256.65 m** |

Without the split floor, 33 of 83 domains' beach start moved *landward*, up to
27 m. Testing only "2009 saw it" was not enough and still left 16 - the 2009
holes are not empty, 2014 fills them. Without the ceiling, 250 m spikes landed
in 66 of 90 domains against a baseline island-wide max of 10.18 m.

The guards matter more without a road boundary, not less: 918,762 cells are now
rejected by the floor-to-replace test and 72,948 by the ceiling. The run's
landward-mover alarm stays silent - 12 domains move landward, all by less than
one 10 m model cell, worst -6.0 m.

## Vertical datum is inferred, not declared

The .prj carries the horizontal CRS only; the vertical rests on `geoid18` in
the NOAA Digital Coast path the metadata cites. What supports it empirically:
median(2009 - 1996) over co-measured cells is **-0.255 m** island-wide -
decimetres, and correctly signed for 12 years of erosion. A NGVD29 or ellipsoid
mistake would show metres.

Dropping the boundary enlarged the co-measured population by a fifth, so the
audit now carries the check on **both** footprints:
`bias_2009_minus_1996_m` over all co-measured cells and
`bias_2009_minus_1996_ocean1984_m` over the old ocean-side-of-1984 band. Their
medians are -0.255 and -0.244 m. The datum argument survives the change on its
original footprint rather than being quietly restated against a different set
of cells.

Domains 11, 85 and 86 report ~-1.6 to -1.9 m on slivers of overlap; do not read
those as a datum finding.

## Figures

Regenerated 2026-08-26 against the no-boundary product. At 10 m the mosaic is
**1996 25.0% / 2009 32.3% / 2014 42.7%** of cells.

| file | shows |
|---|---|
| `HAT_mosaic1984_island.png` | whole island: 2009 only / the mosaic / survey source |
| `HAT_mosaic1984_zoom_76_81.png` | the developed reach, biggest seaward shift |
| `HAT_mosaic1984_roads_8_15.png` | domains 8-15 with BOTH alignments |
| `HAT_mosaic1984_roads_82_88.png` | domains 82-88 with BOTH alignments |
| `HAT_mosaic1984_shift.png` | per-domain beach-start shift |
| `HAT_1984dem_holes.png` | unsurveyed and sub-MHW cells, dune search band vs interior |

`HAT_1984dem_holes.png` is drawn by `3-figures/HAT_plot_dem_holes.py` and reads
the **exported `.npy` arrays**, not the rasters, so it shows what CASCADE
ingests rather than what the fill chain produced.

Each cell is classified within a per-profile **island envelope** running from
beach start (first cell above 0.50 m MHW) to the last cell above 0 m MHW. Only
cells inside that envelope count as affected, so ocean and sound are excluded
by construction rather than by a threshold, and a clean DEM draws quiet.
Unsurveyed cells stay a separate category from measured-but-low ones for the
reason in FIGURES.md.

The **dune search band** - 0-80 m landward of each profile's own beach start,
the extractor's default window - carries **0.54%** affected cells island-wide,
worst domain 17 at 3.25%. Everything else sits in the back-barrier: 2.25% of
all cells are at or below MHW inside the envelope and 0.32% are unsurveyed,
concentrated in domains 10 (63%), 11 (58%), 9 (30%) and 66 (26%). That is low
marsh and the ALACE swath's landward edge, not a fill failure - but it is why
those domains are the fragile ones, and worth having in view before the v1 dune
pick pass picks windows that could reach into it.

Both alignments are drawn: **2004 solid black, 1984 white-dashed on top**, the
same key `HAT_plot_gapfill.py` uses. Now that neither gates the product, they
read as annotation rather than as the boundary.

Add or change a zoom in `ZOOMS` at the top of
`3-figures/HAT_plot_1984_mosaic.py`: `(domain ids, filename slug, title)`.

## Status

Rebuilt 2026-08-26 with no road boundary. `1-gapfill-1m/`, `2-resampled-10m/`
and the `1984-start/npy-arrays/` export are all current. **No dune/topo
extraction exists yet** - that pick pass is the next step.
