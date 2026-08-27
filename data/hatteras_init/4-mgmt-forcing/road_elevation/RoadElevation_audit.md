# NC-12 road elevation from the 2009 LiDAR

Generated 2026-08-26T17:00:04 by `HAT_road_elevation.py`.

| | |
|---|---|
| Surface | `clip_domain_<N>_filled.tif` — **1 m** native LiDAR, 2009, holes filled from **2009-2014** |
| Why that fill | The 2008 NOAA IOCM survey was preferred here until 2026-08-26 — one year from the 2009 base, so the same pavement, where a 2014 Post-Sandy surface at GIS 78–80 may be a REBUILT road (post-Irene, post-breach). That product is deleted and not reproducible from this repo, so the baseline is used instead. The cost is bounded: only GIS 78–80 have NoData under the 2004 alignment, and the two fills differ by ≤ 0.015 m in the corridor. |
| Why NOT the 1984 product | `2009-2014-1996` sits **+0.222 m** higher through this corridor (median; 54 of 82 domains beyond 0.05 m, cell counts identical, so ALACE REPLACED measured 2009 pavement). That is the uncorrected island-wide 1996-vs-2009 survey offset — `mosaic_1984_audit.csv` puts it at median +0.255 m — not a roadbed, and it is kept out of the forcing. One elevation set, on the baseline, for both periods. |
| Domains materially affected | **GIS 78, 79, 80** — the only ones with NoData under the road in the unfilled 2009 clips. Everywhere else the corridor was already complete, and the change is ≤ 0.01 m: the fill pipeline rebuilds the raster rather than patching the original, so the corridor mask lands on a few cells' difference even where there was nothing to fill. Not a no-op, but well inside the 0.07 m within-domain scatter. |
| Alignment | `nc12_2004.geojson` — the **2004** digitised centreline |
| Corridor | 3.5 m buffer either side (~7 m strip) |
| Aggregation | **mean** of all valid 1 m cells in the corridor |
| Datum written | **m MHW-relative** (MHW = 0.36 m NAVD88) |
| Domains | 9–90 (82 with data) |
| Applies to | **both** the 1984 and the 2004 start period |

## Why one file and not one per vintage

There is only one DEM. The previous script sampled two alignments — 1984 and 2004 — on that same 2009 surface, so any difference between the two outputs was never a difference in **time**. It was a difference in **where on the 2009 surface the two digitised lines happened to fall**.

Where NC-12 never moved, the two lines sit on top of each other and the numbers were identical by construction. Where the road *was* relocated, the 1984 line lay over the abandoned corridor, and the "1984 road elevation" there was the elevation of a place a road used to be — measured below.

Neither of those is a measurement of change in roadbed height. Writing two files implied one. This writes one.

## How the road elevation changes along the island

Two views, both regenerated on every run:

- `HAT_road_elevation_map.png` — the road drawn **in real geography**, coloured by these numbers, over the 10 m LiDAR. A domain index is not a place; this is where the high and low reaches actually are, and it shows the island around them.
- `HAT_road_elevation.png` — the alongshore profile with the p10–p90 spread inside each domain and the reach means marked.

| Reach | GIS | n | Mean (m NAVD88) | Mean (m MHW) | Range |
|---|---:|---:|---:|---:|---|
| Buxton-Avon — 7 relocated | 9–20 | 12 | 1.57 | 1.21 | 1.24–1.94 |
| Avon | 21–31 | 11 | 0.98 | 0.62 | 0.76–1.25 |
| Wimble Shoals | 32–67 | 36 | 1.01 | 0.65 | 0.76–1.30 |
| Tri-Village | 68–83 | 16 | 0.96 | 0.60 | 0.66–1.22 |
| Pea Island — 4 relocated | 84–90 | 7 | 1.17 | 0.81 | 0.92–1.44 |

The spread between reaches is **0.61 m** — Buxton-Avon at 1.57 m against Tri-Village at 0.96 m. That is larger than the within-domain σ by more than an order of magnitude, which is the case for carrying a per-domain array rather than one scalar.

> Note that the highest reach, Buxton-Avon, is also the relocated one. Its elevation is a **constructed embankment** built during the relocation, not the natural roadbed grade — which is why the profile steps down sharply at its northern end.

## Mean, not median

The product is the flat unweighted mean of every valid 1 m cell in the corridor. The median is carried in the per-domain CSV as `elev_median_navd`, with `mean_minus_median` alongside, so a domain where the mean is being dragged by a structure in the corridor — a bridge deck, a house, a driveway apron — announces itself instead of being silently absorbed.

On this alignment the two agree closely: median |mean − median| = **0.005 m**, max **0.094 m** (GIS 16).

## Internal checks

### The external validation is gone

This method was originally validated against an `elevation_2009` column in `2004_road_offset_raw.csv`, sampled independently in ArcGIS Pro at 1 m transect points; the two agreed to a median of −0.02 m. **That check is no longer reproducible.** The file is now `nc12_2004.csv`, the `elevation_2009` column is absent, and its apparent successors (`avg_elev_m`, `z_mean`, `z_max`, `z_min`, `road_z`, `relief_m`) are all zero or all empty across the 1491 rows. The agreement was real when it was measured; it cannot be re-measured from what is in the repository now. If a pre-rename copy of that CSV turns up, restore the check.

What follows is internal: does the answer depend on choices we made, and does the surface the model runs on agree with the surface we measured.

### Corridor width

If the answer moved with the buffer, the buffer would be the measurement. It does not, until the corridor grows wide enough to reach off the pavement:

| Buffer | Strip | Island median (m NAVD88) | Median within-domain σ |
|---:|---:|---:|---:|
| 2.0 m | 4 m | 1.04 | 0.06 |
| 2.5 m | 5 m | 1.04 | 0.07 |
| 3.5 m | 7 m | 1.03  ← **used** | 0.07 |
| 5.0 m | 10 m | 1.01 | 0.08 |
| 8.0 m | 16 m | 0.97 | 0.11 |
| 12.0 m | 24 m | 1.00 | 0.18 |

σ climbing while the median falls is the corridor beginning to catch the shoulder and the dune toe. 3.5 m sits in the flat part.

### 1 m clip vs the 10 m Barrier3D grid

The same corridor sampled on `resampled_domain_<N>.tif`, the grid the model actually runs on:

- median difference **-0.00 m** over 82 domains
- mean |difference| **0.01 m**, max **0.06 m**
- cells per domain in the corridor: **3534** at 1 m vs **35** at 10 m
- median within-domain σ: **0.07 m** at 1 m vs **0.07 m** at 10 m

**This does not go the way the inherited rationale said it would, and the difference matters.** The 1 m clip was originally chosen because the 10 m grid gave GIS 9 a σ of 1.59 m — the foredune, not the roadbed. But that was measured on the **1984** line, which crosses a relocation scar. On the 2004 line the two resolutions agree to a couple of centimetres.

So for this product the 1 m clip is a **precaution, not a rescue**. It is still the right surface — it rests on ~3534 cells per domain against ~35, and the near-equal σ at 10 m is a small-sample artefact rather than evidence of equal fidelity — but the honest statement is that switching to the resample would barely move these numbers.

### Relocation bracket — what is under the 1984 line?

The case for using the 2004 alignment everywhere rests on a claim about the abandoned corridor, and a claim is not evidence. Same surface, same corridor width, the other line, in the 11 domains where the two disagree:

| GIS | 2004 line (used) | 1984 line | Difference | σ on 1984 line |
|---:|---:|---:|---:|---:|
| 9 | 1.69 | 3.10 | +1.41 | 1.70 |
| 10 | 1.74 | 4.26 | +2.52 | 0.82 |
| 11 | 1.93 | 2.06 | +0.14 | 0.64 |
| 12 | 1.73 | 2.52 | +0.79 | 0.16 |
| 13 | 1.66 | 2.63 | +0.97 | 0.59 |
| 14 | 1.68 | 1.98 | +0.30 | 0.14 |
| 15 | 1.94 | 1.93 | -0.00 | 0.15 |
| 84 | 1.24 | 1.56 | +0.32 | 0.87 |
| 85 | 1.17 | 2.18 | +1.01 | 0.95 |
| 86 | 1.05 | 1.72 | +0.67 | 0.62 |
| 87 | 1.16 | 2.05 | +0.90 | 1.26 |

Mean **1.54 m** on the line used against **2.36 m** on the 1984 line, whose within-domain σ reaches **1.70 m**. Un-relocated domains in the same named reaches sit at **1.26 m**.

**The two lines do not bracket the neighbouring grade — both are above it.** So the intuitive picture is wrong: the abandoned corridor is not overwashed, bulldozed flat, or bare sand. The dune migrated over it after the road left, and a σ of that size is a dune, not a graded surface.

This settles the choice rather than complicating it. The 2004 value is the **lower** of the two candidates and the **only** one measured on a road surface — the conservative choice as well as the correct one. Averaging the two, or substituting the 1984 line in these domains, would raise `road_ele` and make the 1984 run worse.

### Alongshore continuity

A maintained road does not step abruptly between adjacent 500 m cells. Any adjacent-domain step above 0.35 m is flagged `JUMP`.

**1 step(s) flagged:**

- GIS 16, step -0.57 m from GIS 15: **expected** — this is the boundary of a relocated reach. The rebuilt road sits on a constructed embankment; the step is where that embankment ends, not a sampling error.


## Datum — not ambiguous, despite the runner

`bulldoze()` writes `road_ele` straight into the interior grid:

```python
road_ele = road_ele / dz
new_road_domain = np.zeros(...) + road_ele
```

and the interior arrays are **MHW-relative**, because `HAT_dune_topo_extractor.py` subtracts `MHW_M = 0.36` before anything else. So `road_ele` must be MHW-relative metres. There is no reading under which NAVD88 is correct. **`RoadElevation.csv` is already MHW-relative — do not subtract MHW again.**

The runner's current `ROAD_ELEVATION = 1.45` is **high under either reading**: the measured island median is 1.03 m NAVD88 = 0.67 m MHW. Read as NAVD88, 1.45 would be 1.09 m MHW; read as MHW-relative it would be 1.81 m NAVD88. The measured range is 0.66–1.94 m NAVD88.

> Replacing the scalar changes model behaviour, not just bookkeeping. `bulldoze` computes `road_overwash_removal = sum(old − road_ele)`, clipped at zero, so a **lower** road removes **more** sand each year and hands more of it to the dunes — maintaining a lower roadbed means digging out more. `road_ele` is also decremented by RSLR each year and triggers `_drown_break` when it passes 0 m MHW, so a lower start is closer to that threshold.

## Two things this file does not correct

### 1. The time gap

The DEM is **2009**; one run starts in **1984**. CASCADE decrements `road_ele` by RSLR every year, so the 1984 run begins with a roadbed that is already ~25 years of sea-level rise low relative to its own MHW. No back-correction is applied. These are measurements of a surface that exists, not reconstructions of one that does not — a reconstructed 1984 roadbed would also have reintroduced the two-file split this document exists to remove.

### 2. The relocations

NC-12 was relocated in GIS 9–15 (1999) and GIS 84–87 (1989). The 2004 alignment sampled here is the **post**-relocation one. For the 2004 run that is correct. For the 1984 run, the road in those 11 domains was physically somewhere else, and the number written is the elevation of the corridor the road would later occupy.

This is flagged, not adjusted, and the RELOCATION BRACKET above is why: sampling the 1984 line instead does not fix it, because that line now runs under the foredune. The affected domains are marked `RELOCATED_PRE_2004` below and crossed on the QC figure.

The residual is bounded and worth stating rather than hiding. In 11 of 82 domains, for the first 15 years of the **1984** runs only, `road_ele` may be too high by roughly the gap to the neighbouring grade. Because `bulldoze` removes `sum(old − road_ele)` clipped at zero, too high means the model removes **less** overwash than it should and hands **less** sand to the dunes there. No measurement closes that gap; substituting a modelled estimate would trade a known bias for an invented one.

| GIS | m NAVD88 | m MHW | Relocation |
|---:|---:|---:|---|
| 9 | 1.69 | 1.33 | 1999 (Buxton–Avon) |
| 10 | 1.74 | 1.38 | 1999 (Buxton–Avon) |
| 11 | 1.93 | 1.56 | 1999 (Buxton–Avon) |
| 12 | 1.73 | 1.37 | 1999 (Buxton–Avon) |
| 13 | 1.66 | 1.30 | 1999 (Buxton–Avon) |
| 14 | 1.68 | 1.32 | 1999 (Buxton–Avon) |
| 15 | 1.94 | 1.58 | 1999 (Buxton–Avon) |
| 84 | 1.24 | 0.88 | 1989 (Pea Island) |
| 85 | 1.17 | 0.81 | 1989 (Pea Island) |
| 86 | 1.05 | 0.69 | 1989 (Pea Island) |
| 87 | 1.16 | 0.80 | 1989 (Pea Island) |

## Flags

| Flag | Meaning |
|---|---|
| `THIN(n)` | Fewer than 200 valid 1 m cells — the line barely crosses this domain. |
| `NODATA(x%)` | More than 25% of the corridor fell on LiDAR NoData. |
| `SCATTER(σ)` | Within-domain σ above 0.40 m — not a graded surface. The corridor is catching an embankment or the dune toe. |
| `JUMP(±x)` | Step from the previous domain above 0.35 m. |
| `RELOCATED_PRE_2004` | Correct for the 2004 run; the wrong place for the 1984 run. Value unadjusted. |

## Per-domain results

| GIS | m NAVD88 | m MHW | median NAVD88 | σ | p10–p90 | n cells | Flags |
|---:|---:|---:|---:|---:|---|---:|---|
| 9 | 1.69 | 1.33 | 1.70 | 0.06 | 1.61–1.77 | 3526 | `RELOCATED_PRE_2004` |
| 10 | 1.74 | 1.38 | 1.72 | 0.09 | 1.65–1.88 | 3529 | `RELOCATED_PRE_2004` |
| 11 | 1.93 | 1.56 | 1.91 | 0.08 | 1.82–2.05 | 3538 | `RELOCATED_PRE_2004` |
| 12 | 1.73 | 1.37 | 1.70 | 0.09 | 1.64–1.89 | 3537 | `RELOCATED_PRE_2004` |
| 13 | 1.66 | 1.30 | 1.66 | 0.05 | 1.60–1.72 | 3551 | `RELOCATED_PRE_2004` |
| 14 | 1.68 | 1.32 | 1.68 | 0.05 | 1.62–1.74 | 3555 | `RELOCATED_PRE_2004` |
| 15 | 1.94 | 1.58 | 1.97 | 0.12 | 1.74–2.07 | 3571 | `RELOCATED_PRE_2004` |
| 16 | 1.37 | 1.01 | 1.28 | 0.20 | 1.20–1.73 | 3535 | `JUMP(-0.57)` |
| 17 | 1.24 | 0.88 | 1.24 | 0.05 | 1.18–1.30 | 3522 |  |
| 18 | 1.26 | 0.90 | 1.26 | 0.05 | 1.20–1.33 | 3527 |  |
| 19 | 1.30 | 0.94 | 1.29 | 0.07 | 1.22–1.40 | 3597 |  |
| 20 | 1.34 | 0.98 | 1.36 | 0.10 | 1.20–1.45 | 3589 |  |
| 21 | 1.25 | 0.89 | 1.24 | 0.12 | 1.09–1.42 | 3516 |  |
| 22 | 1.02 | 0.66 | 1.02 | 0.06 | 0.94–1.11 | 3531 |  |
| 23 | 1.10 | 0.74 | 1.10 | 0.05 | 1.04–1.16 | 3584 |  |
| 24 | 1.04 | 0.68 | 1.04 | 0.05 | 0.98–1.10 | 3590 |  |
| 25 | 1.05 | 0.69 | 1.05 | 0.06 | 0.97–1.13 | 3569 |  |
| 26 | 0.97 | 0.61 | 0.97 | 0.05 | 0.91–1.03 | 3514 |  |
| 27 | 1.06 | 0.70 | 1.06 | 0.06 | 0.97–1.13 | 3521 |  |
| 28 | 0.91 | 0.55 | 0.93 | 0.08 | 0.80–1.00 | 3595 |  |
| 29 | 0.76 | 0.40 | 0.76 | 0.04 | 0.70–0.81 | 3596 |  |
| 30 | 0.78 | 0.42 | 0.79 | 0.04 | 0.73–0.84 | 3592 |  |
| 31 | 0.82 | 0.46 | 0.82 | 0.06 | 0.76–0.90 | 3646 |  |
| 32 | 0.81 | 0.45 | 0.81 | 0.04 | 0.75–0.86 | 3609 |  |
| 33 | 0.79 | 0.42 | 0.78 | 0.06 | 0.71–0.88 | 3552 |  |
| 34 | 0.77 | 0.41 | 0.77 | 0.04 | 0.71–0.83 | 3563 |  |
| 35 | 0.76 | 0.40 | 0.76 | 0.05 | 0.70–0.82 | 3602 |  |
| 36 | 0.77 | 0.41 | 0.77 | 0.04 | 0.71–0.82 | 3589 |  |
| 37 | 0.86 | 0.50 | 0.85 | 0.08 | 0.77–0.99 | 3593 |  |
| 38 | 0.98 | 0.62 | 0.99 | 0.07 | 0.88–1.07 | 3591 |  |
| 39 | 1.14 | 0.78 | 1.12 | 0.10 | 1.01–1.28 | 3549 |  |
| 40 | 1.26 | 0.90 | 1.25 | 0.11 | 1.13–1.42 | 3513 |  |
| 41 | 1.18 | 0.82 | 1.19 | 0.09 | 1.04–1.29 | 3513 |  |
| 42 | 1.05 | 0.69 | 0.99 | 0.16 | 0.86–1.26 | 3505 |  |
| 43 | 0.93 | 0.57 | 0.93 | 0.06 | 0.86–1.00 | 3504 |  |
| 44 | 0.94 | 0.58 | 0.94 | 0.07 | 0.86–1.03 | 3503 |  |
| 45 | 1.12 | 0.76 | 1.11 | 0.08 | 1.02–1.24 | 3508 |  |
| 46 | 1.26 | 0.91 | 1.26 | 0.05 | 1.20–1.34 | 3508 |  |
| 47 | 1.23 | 0.87 | 1.23 | 0.06 | 1.16–1.30 | 3509 |  |
| 48 | 1.30 | 0.94 | 1.30 | 0.07 | 1.20–1.39 | 3506 |  |
| 49 | 1.30 | 0.94 | 1.29 | 0.08 | 1.21–1.43 | 3504 |  |
| 50 | 1.20 | 0.83 | 1.20 | 0.05 | 1.13–1.27 | 3509 |  |
| 51 | 1.18 | 0.82 | 1.19 | 0.08 | 1.05–1.28 | 3508 |  |
| 52 | 0.94 | 0.57 | 0.93 | 0.07 | 0.85–1.04 | 3513 |  |
| 53 | 0.95 | 0.59 | 0.96 | 0.12 | 0.78–1.09 | 3521 |  |
| 54 | 1.17 | 0.81 | 1.18 | 0.06 | 1.09–1.25 | 3518 |  |
| 55 | 1.16 | 0.80 | 1.17 | 0.08 | 1.06–1.25 | 3523 |  |
| 56 | 1.01 | 0.65 | 1.01 | 0.07 | 0.91–1.10 | 3520 |  |
| 57 | 0.94 | 0.58 | 0.94 | 0.05 | 0.87–1.01 | 3516 |  |
| 58 | 0.94 | 0.58 | 0.94 | 0.05 | 0.88–1.00 | 3521 |  |
| 59 | 0.92 | 0.56 | 0.91 | 0.04 | 0.86–0.97 | 3517 |  |
| 60 | 0.81 | 0.45 | 0.81 | 0.06 | 0.73–0.88 | 3528 |  |
| 61 | 0.94 | 0.58 | 0.94 | 0.08 | 0.83–1.04 | 3551 |  |
| 62 | 0.91 | 0.55 | 0.92 | 0.10 | 0.78–1.03 | 3553 |  |
| 63 | 0.99 | 0.63 | 0.98 | 0.06 | 0.91–1.08 | 3550 |  |
| 64 | 1.05 | 0.69 | 1.05 | 0.07 | 0.96–1.15 | 3548 |  |
| 65 | 1.00 | 0.64 | 0.98 | 0.08 | 0.90–1.12 | 3549 |  |
| 66 | 0.90 | 0.54 | 0.90 | 0.07 | 0.81–0.98 | 3555 |  |
| 67 | 0.91 | 0.55 | 0.89 | 0.08 | 0.80–1.02 | 3532 |  |
| 68 | 1.00 | 0.64 | 0.98 | 0.10 | 0.88–1.15 | 3533 |  |
| 69 | 1.22 | 0.86 | 1.22 | 0.07 | 1.12–1.31 | 3664 |  |
| 70 | 0.92 | 0.56 | 0.95 | 0.14 | 0.74–1.10 | 3809 |  |
| 71 | 0.66 | 0.30 | 0.66 | 0.07 | 0.57–0.76 | 3541 |  |
| 72 | 0.83 | 0.47 | 0.81 | 0.08 | 0.73–0.95 | 3500 |  |
| 73 | 0.95 | 0.59 | 0.94 | 0.10 | 0.82–1.09 | 3500 |  |
| 74 | 0.84 | 0.48 | 0.84 | 0.11 | 0.70–0.99 | 3528 |  |
| 75 | 1.03 | 0.67 | 1.02 | 0.07 | 0.95–1.12 | 3580 |  |
| 76 | 0.99 | 0.63 | 0.99 | 0.07 | 0.91–1.09 | 3509 |  |
| 77 | 0.84 | 0.48 | 0.83 | 0.07 | 0.75–0.94 | 3500 |  |
| 78 | 0.91 | 0.55 | 0.91 | 0.11 | 0.78–1.05 | 3501 |  |
| 79 | 1.01 | 0.66 | 1.03 | 0.10 | 0.86–1.13 | 3583 |  |
| 80 | 0.98 | 0.62 | 0.97 | 0.07 | 0.90–1.09 | 3512 |  |
| 81 | 1.08 | 0.72 | 1.08 | 0.07 | 0.99–1.17 | 3662 |  |
| 82 | 1.14 | 0.78 | 1.15 | 0.08 | 1.03–1.25 | 3572 |  |
| 83 | 0.99 | 0.63 | 0.98 | 0.09 | 0.87–1.10 | 3500 |  |
| 84 | 1.24 | 0.88 | 1.25 | 0.09 | 1.11–1.34 | 3599 | `RELOCATED_PRE_2004` |
| 85 | 1.17 | 0.81 | 1.20 | 0.10 | 1.02–1.28 | 3609 | `RELOCATED_PRE_2004` |
| 86 | 1.05 | 0.69 | 1.04 | 0.06 | 0.98–1.13 | 3504 | `RELOCATED_PRE_2004` |
| 87 | 1.16 | 0.80 | 1.16 | 0.08 | 1.06–1.26 | 3620 | `RELOCATED_PRE_2004` |
| 88 | 1.44 | 1.08 | 1.38 | 0.31 | 1.09–1.93 | 3639 |  |
| 89 | 1.21 | 0.85 | 1.22 | 0.08 | 1.09–1.31 | 3529 |  |
| 90 | 0.92 | 0.56 | 0.90 | 0.12 | 0.78–1.08 | 3536 |  |
