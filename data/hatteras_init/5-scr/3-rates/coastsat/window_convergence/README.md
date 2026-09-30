# window_convergence — how long a window does the shoreline rate need?

**Start here.** The model is graded on a CoastSat shoreline change rate over
1996–2010 (and runs 2010–2024 as the second leg). Is a 15-year window long
enough to give the same rate as the full 1996–2024 record, and if not, how
long does it have to be?

Every product here fits the target's own OLS rate (`coastsat_lrr.compute_lrr`)
over a family of **nested** windows that grow one year at a time until they
reach 1996–2024:

- **forward_from_1996**: the start is pinned at 1996 and the end walks out
  (1996–1997, 1996–1998 … 1996–2024). How much record do you need from 1996?
- **backward_from_2024**: the end is pinned at 2024 and the start walks back
  (2023–2024, 2022–2024 … 1996–2024). How late can a window start?

Because every window sits inside 1996–2024, every curve reaches the reference
**by construction**. Read *how soon* it gets there, never *whether*.

## The two questions

```
1-rate_profiles/       Does the WHOLE alongshore profile look like 1996–2024?
    forward_from_1996/     one line per window over the reference, all 906
    backward_from_2024/    transects; scored by alongshore Pearson r (shape)

2-settling_window/     How many years does each PLACE need before its rate
                       stays within ±0.25 / ±0.5 / ±1.0 m/yr of 1996–2024?
    years_needed_alongshore.png   THE figure: both directions, 906 transects
    forward_from_1996/     a-eight_sites/     positions through time at 8
    backward_from_2024/                       transects (supporting figure)
                           b-every_transect/  tables behind the figure
                           c-domain_means/    tables, averaged to 90 domains

```

## The answer

| | forward from 1996 | backward from 2024 |
|---|---|---|
| **1. Whole profile.** r of the 15-yr model window | 0.69 (1996–2010) | 0.80 (2010–2024) |
| Shortest window with r ≥ 0.9 | 23 yr (1996–2018) | 20 yr (2005–2024) |
| **2. Each place.** Transects where 15 yr is already within ±1.0 / ±0.5 / ±0.25 m/yr | 36% / 10% / 0% | 46% / 19% / 5% |
| Median domain settling window (CI overlap, in the tables) | 1996–2021 (26 yr) | 2003–2024 (22 yr) |

Both halves of the canonical chain are too short to recover the long-term
rate, whether you ask for the profile's shape or each location's value.

**Deleted 2026-09-29 (Hannah):** `experiments/record_cut_2020/`, the settling
sweep on a record cut at 2020, before the 2021 CoastSat step. Its median
settling windows were 1996–2015 (20 yr) forward and 2002–2020 (19 yr) backward,
3–6 years shorter than on the full record, but still not 15. To regenerate it:
`coastsat_window_convergence.py --ref-end 2020`.

## Figures to open, in order

1. [`1-rate_profiles/forward_from_1996/window_profiles_forward_from_1996.png`](1-rate_profiles/forward_from_1996/window_profiles_forward_from_1996.png)
   shows every window over the reference, then r against window length. The
   `_panels_` twin has one panel per window.
2. [`2-settling_window/years_needed_alongshore.png`](2-settling_window/years_needed_alongshore.png)
   gives the years each transect needs at three plain thresholds, with the
   model's 15 years marked.
3. [`2-settling_window/forward_from_1996/a-eight_sites/shoreline_position_window_fits_forward_from_1996.png`](2-settling_window/forward_from_1996/a-eight_sites/shoreline_position_window_fits_forward_from_1996.png)
   shows why, at eight places: the raw positions with the fits drawn on them.

1-rate_profiles and a-eight_sites each have a `backward` twin. Captions are in
each folder's `supporting/CAPTIONS.md`.

## Producers

```
scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/
    coastsat_window_profiles.py      -> 1-rate_profiles/
    coastsat_window_convergence.py   -> 2-settling_window/
```

The layout is resolved in `scripts/site_layer/hat_observed_rates.py`
(`window_profiles_dir`, `window_convergence_dir`, `SETTLING_SCALE_DIRS`).
Reorganised 2026-09-29 from `record_<start>_<end>/<direction>/{sites,
all_transects,domain_means,alongshore_profiles}/`. The data did not change.
