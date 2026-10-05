# 1-observations/detrended_position — the signal the whole island shares

Every CoastSat transect detrended against its **own** 1996–2024 fit,
reduced to one annual median, and averaged over all 906 of them. Anything local
is incoherent between transects and cancels; what survives is common to the
island.

## Why it exists

`3-rates/coastsat/window_convergence_1996_2024/` found that no window shorter than about
25 years recovers the long-term rate, and that the answer **barely varies
between transects** — a fast, clean transect needs as long as a slow, noisy
one (correlation between a transect's noise-to-trend ratio and its convergence
time: 0.01). A per-transect cause cannot produce a per-transect-invariant
answer, so the cause had to be shared. This is the search for it.

## What it found

One coherent excursion: roughly flat 1996–2004, a landward sag through
2005–2020 bottoming at −7.4 m, then **+16.7 m in the single year
to 2021**, held through 2024. It is only 17% of the mean transect variance
— but it is the *coherent* part, and coherence is what moves an OLS slope
systematically while incoherent scatter averages out inside the window.

The bias it alone puts on a fitted rate:

| window | bias from the anomaly | against a median rate of ~1 m/yr |
|---|---|---|
| 1996–2010 | −0.367 m/yr | 37% of it, erosional |
| 2010–2024 | +0.875 m/yr | 87% of it, the other way |
| 1996–2024 | 0.000 m/yr | zero by construction |

That is the answer to *why the windows disagree*: not noise, but a shared
multi-decadal excursion an OLS slope cannot separate from the trend until the
window spans the whole of it.

## The 2021 step: what is known (2026-10-01)

**Start with `shoreline_position_since_1996.png`**: the same jump in plain
shoreline positions, no detrending (each transect's yearly median minus its
1996 median, never-nourished transects only). The median transect stays
within about 10 m of 1996 for 24 years, moves +18.9 m seaward from 2020 to
2021, the largest year-to-year move in the record by double, and eases back
to +4.4 m by 2024. Table: `shoreline_position_since_1996.csv`; producer
`coastsat_position_since_1996.py`.

Detrended figure: `detrended_position_2021_step.png`.

**When.** The island mean sits at −7.4 m in 2020, its most landward year, and
+9.4 m in 2021: **+16.7 m in one year**. It holds and drifts slowly back
(+9.2, +8.2, +6.4 m through 2024). The never-nourished transects (672) and
the later-nourished ones (234) jump together, +16.5 and +17.4 m.

**Not the nourishments.** The fills in the modelled reach are Rodanthe 2014
and Buxton + Avon 2022. The jump lands a year before the 2022 fills and is
the same on ground never filled. The fills are visible separately, as a
local step the year after placement: Rodanthe +15.4, Avon +11.3,
Buxton +5.0 m.

**Where, measured before any fill.** `step_2021_by_domain_prefill.csv`, 2021
minus 2020 per domain (and 2021 minus 2019, which skips the thin 2020):

| reach | 2020 → 2021 | 2019 → 2021 |
|---|---|---|
| Cape Point–Buxton (1–15) | +15.9 m | +11.7 m |
| Avon (16–31) | +23.0 m | +20.0 m |
| open reach (32–66) | +16.9 m | +15.4 m |
| Tri-Village–Rodanthe (67–90) | +12.6 m | +9.7 m |

Island median +16.3 m; 87 of 90 domains seaward (landward: 15, 35, 36). The
step is seaward almost everywhere and only mildly larger in the south. The
**2021–2024 step in `detrended_position.png` (b) overstates the south**:
its window includes the 2022 Buxton and Avon fills (Avon GIS 24 +41.7 m) and
the Cape Point shoal attachment after 2021 (GIS 1 +71.5 m against +23.6 m
from 2020 to 2021).

**Real, not a satellite artefact.**
- The dune line agrees: 2009 → 2023, CoastSat +3.8 m, dune line +3.7 m
  (alongshore r = 0.71). Three snapshots cannot date it to 2021.
- The step appears inside every calendar quarter on its own (+10 to +13 m),
  so it is not a shift in which season was imaged.
- Neighbouring transects agree to 2.8 m while domains differ by 9–11 m; a
  sensor or waterline bias would be one offset everywhere.

**Sharpened by sampling.** 2020 has 12.8 images per transect, 2021 has 30.6.
2020 is both the thinnest and the most landward year, so the single-year jump
is somewhat exaggerated; 2019 → 2021 is still +14.6 m.

**Cause: not tested.** Recovery after Dorian (2019) and Isaias (2020), or a
sediment pulse, would fit the pattern. Nothing here tests either.

**Why it matters.** The step alone puts −0.37 m/yr on the 1996–2010 rate and
+0.88 m/yr on 2010–2024, and drives most of the window_convergence bias
(`3-rates/coastsat/window_convergence_1996_2024/`).

```
annual_medians_detrended.csv   year x transect, m. The matrix every other
                               script here reads, so the detrending happens
                               once and cannot drift
detrended_position_by_year.csv       a row per year: the index, the nourished /
                               untouched split, the sampling density
detrended_position_by_domain.csv   a row per domain: its share of the step, and
                               how well it tracks the index
detrended_position.png             the index, the alongshore step, the record
detrended_position_2021_step.png   the step alone: when (never-nourished vs
                               nourished), and where, from years before the fills
step_2021_by_domain_prefill.csv    a row per domain: 2020->2021 and 2019->2021
attribution_*                  written by coastsat_position_attribution.py — what
                               the step IS
```

## What this is not

Not a rate product and not a model input. Nothing is graded against it. It
lives under `1-observations/` because its subject is the observed record
itself, not a rate fitted from it.

Producers: `scripts/input_prep/5-scr/1-observations/detrended_position/`
(`coastsat_detrended_position.py` builds it, `coastsat_position_attribution.py` tests it). Built
2026-09-23 by interview (Hannah).
