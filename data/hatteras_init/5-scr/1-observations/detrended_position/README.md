# 1-observations/detrended_position — the signal the whole island shares

Every CoastSat transect detrended against its **own** 1996–2024 fit,
reduced to one annual median, and averaged over all 906 of them. Anything local
is incoherent between transects and cancels; what survives is common to the
island.

## Why it exists

`3-rates/coastsat/window_convergence/` found that no window shorter than about
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

```
annual_medians_detrended.csv   year x transect, m. The matrix every other
                               script here reads, so the detrending happens
                               once and cannot drift
detrended_position_by_year.csv       a row per year: the index, the nourished /
                               untouched split, the sampling density
detrended_position_by_domain.csv   a row per domain: its share of the step, and
                               how well it tracks the index
detrended_position.png             the index, the alongshore step, the record
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
