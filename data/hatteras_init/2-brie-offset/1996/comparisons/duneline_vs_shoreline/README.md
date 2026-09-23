# 1996 island offset: duneline vs shoreline

Two builds of the **same** 1996 island offset, from two **different
features** on the island:

| source | what it is | build |
|---|---|---|
| `duneline` | dune line (1997 imagery) | `../../v1/` |
| `shoreline` | CoastSat shoreline (1995–1997 mean) | `../../shoreline/v1/` |

Written by `scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py`.
This is a comparison, not a build: nothing here is read by a model run.

## What it found

| | |
|---|---|
| shoreline seaward of duneline (fixed datum) | mean +37.4 m, median +31.7, sd 17.2, range +16.7 to +96.6 |
| domains with shoreline seaward | 90 of 90 |
| gap between the two zeroing baselines | +45.6 m |
| difference in the model frame | mean +8.2 m, sd 17.2, range -51.0 to +28.9 |

## Read the correlation, not the way round you expect

The two profiles correlate at **r = 1.0000**. That is not the dune line and
the shoreline agreeing about the beach — it is both of them being dominated by
the same ~6.2 km of cape curvature, against which the 17 m sd of their
difference is 0.3%. **Difference these profiles; never correlate them.**

## The figure

Six consecutive sections of the island on a **2 × 3 grid**, each a **vertical
strip**: alongshore up the page, cross-shore across it, south at the bottom,
ocean on the right. The panels read as the island.

Both profiles are drawn as the **distance from the shared offshore datum**,
not as the min-zeroed offset the model is handed — see the next section for
why that matters. They are the **absolute** profiles, the shape each source
gives, not a residual.

## The dune line is behind the shoreline, and the figure has to say so

On the ground it always is: the shoreline is seaward of the dune line in
**90 of 90 domains**, by 17–97 m. The band between the two lines is that beach.

This figure was drawn in the **model frame** until 2026-09-22, and in that
frame it was not true. Each build is zeroed on its own most seaward domain,
and those two minima are **45.6 m** apart, so

```
model_diff = -(beach width) + 45.6 m
```

and the dune line comes out *apparently seaward* wherever the beach is
narrower than 45.6 m — **69 of 90 domains**, which the old version shaded as
"dune line seaward". That was an artefact of the zeroing, not anything that
happens on this island.

The datum frame carries no such constant and is the **same shape** (a model
offset is this minus the build's own minimum), so nothing about the profiles
was lost by moving to it. Both frames are in the CSV.

**The lesson generalises:** never difference two min-zeroed offset files and
read the sign as physical. `2-brie-offset/README.md` has warned about this
since before the shoreline source existed; this is what it looks like when it
bites.

## Why a grid, and not just more panels

These panels are columns, so each extra one **alongside** takes width from the
offset axis — the axis the gap is measured on — faster than the tighter zoom
gives back. Measured as the widest gap actually rendered, at 300 dpi:

| layout | panel width | gap on paper |
|---|---|---|
| 3 across | 2.07 in | 13 px |
| 4 across | 1.52 in | 10 px |
| **6 across** | 0.96 in | **7 px** — more panels, *less* gap |
| 10 across | 0.52 in | 6 px |

Splitting into **rows** gives the width back, so the zoom is kept and the
offset axis is not squeezed:

| layout | panel width | gap on paper |
|---|---|---|
| **6 as 2 × 3** (this figure) | 2.07 in | **15 px** |
| 9 as 3 × 3 | 2.07 in | 23 px |

Removing a smooth trend from both sources would do better still (14–43% of the
panel), and was built and then taken back out — the absolute shape is the
point. So the gap is carried by the filled band and, where the band is a
sliver, by the **widest-gap number written into each panel**. The per-domain
numbers, in both frames, are in the CSV.

## Two frames, and the sign flip between them

`ORIG_LEN` grows **landward** from the shared offshore datum, so in the
fixed-datum frame `duneline − shoreline` is positive where shoreline is the more seaward
feature — a beach width. The model frame is not that: each build is zeroed on
its own most seaward domain, so differencing the two subtracts the +45.6 m
gap between those baselines and flips the sign,

```
model_diff = -(seaward gap) + +45.6 m
```

which is why a mean beach width of +37.4 m appears in the model frame as
+8.2 m. **The band in the figure is the model-frame gap, so it is not a
beach width.** The `seaward_gap_m` column of the CSV is.

## Files

| file | what it is |
|---|---|
| `offset_1996_duneline_vs_shoreline.csv` | per domain: both sources in both frames, columns named for the source |
| `offset_1996_duneline_vs_shoreline.png` | the three-panel figure (PDF and caption under `supporting/`) |

## Rebuild

```
python scripts/input_prep/2-brie-offset/2-figures/compare_offset_sources.py --year 1996
```

Both sources resolve through their `CURRENT`, so this re-reads whatever each
source currently points at rather than the builds that were current on the day
it was written.
