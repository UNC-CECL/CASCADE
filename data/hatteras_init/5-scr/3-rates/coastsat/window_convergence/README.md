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

## The three questions

```
1-rate_profiles/       What does each window's alongshore profile LOOK like?
    forward_from_1996/     overlay (every window over the reference) and
    backward_from_2024/    panels (one per window, the reference last);
                           window_profiles_transects.csv = every fit, the
                           input to 2-

2-r_bias_rmse/         How CLOSE is each window's profile to 1996–2024?
    window_profiles_r.png            r against window length, both directions
    window_profiles_bias_rmse.png    bias and RMSE, both directions
    window_profiles_r_bias_rmse.csv  the values, with 95% intervals
    forward_from_1996/     the same two figures, one direction, with the
    backward_from_2024/    calendar window on the top axis

3-settling_window/     How many years does each PLACE need before its rate
                       stays within ±0.25 / ±0.5 / ±1.0 m/yr of 1996–2024?
    years_needed_alongshore.png   both directions, 906 transects
    forward_from_1996/     a-eight_sites/     positions through time at 8
    backward_from_2024/                       transects (supporting figure)
                           b-every_transect/  tables behind the figure
                           c-domain_means/    tables, averaged to 90 domains

4-split_windows/       What do the two HALVES look like at one transect?
                       positions through time at eight transects picked by
                       behaviour, the 1996–2024 line and the two windows
                       either side of a cutoff (2010), and both window rates
                       against every cutoff 2000–2020;
                       interactive/ = all 906, with a cutoff slider
```

**The three scores in 2-r_bias_rmse**, each window against 1996–2024 over all
906 transects:

- **r** (alongshore Pearson correlation): are the hotspots in the same
  places? Shape only; ignores the overall level and the size of the swings.
- **bias** (mean of window rate minus 1996–2024 rate): is the overall level
  right? Negative means the window is more erosional than the long-term
  record. Misses of opposite sign cancel.
- **RMSE** (root-mean-square of the same difference): how far off is a
  typical transect, sign ignored? RMSE² = bias² + scatter².

The 95% intervals come from 1000 bootstrap resamples of the 90 domains
(domains, not transects, because neighbouring transects move together; a
transect resample makes the interval about three times too narrow).

## The answer

| | forward from 1996 | backward from 2024 |
|---|---|---|
| **2. Profile, at the 15-yr model window** | **1996–2010** | **2010–2024** |
| r | 0.69 (0.59–0.78) | 0.80 (0.74–0.86) |
| bias | −0.49 m/yr (−0.78 to −0.21) | +1.00 m/yr (+0.75 to +1.23) |
| RMSE | 1.51 m/yr (1.31–1.72) | 1.55 m/yr (1.36–1.75) |
| Years until r stays ≥ 0.5 / 0.75 / 0.9 | 10 / 18 / 24 | 6 / 13 / 21 |
| **3. Each place.** Transects where 15 yr is already within ±1.0 / ±0.5 / ±0.25 m/yr | 36% / 10% / 0% | 46% / 19% / 5% |
| Median domain settling window (CI overlap, in the tables) | 1996–2021 (26 yr) | 2003–2024 (22 yr) |

Both halves of the canonical chain are too short to recover the long-term
rate, whether you ask for the profile's shape, its level or each location's
value. 2010–2024 has the better shape (higher r) but the worse level (+1 m/yr
too accretional). The forward bias holds near −0.5 m/yr from 15 to 25 years
and reaches zero only once 2021 is in the window: that is the CoastSat 2021
step, so a longer 1996-start window would not have fixed the level until the
record ran past 2021.

The r ≥ 0.9 forward window was given as 23 yr before 2026-10-01: r there is
0.8996, which rounds to 0.90 but is below it. The rule now is the first window
from which r *stays* at or above the level.

**Deleted 2026-09-29 (Hannah):** `experiments/record_cut_2020/`, the settling
sweep on a record cut at 2020, before the 2021 CoastSat step. Its median
settling windows were 1996–2015 (20 yr) forward and 2002–2020 (19 yr) backward,
3–6 years shorter than on the full record, but still not 15. To regenerate it:
`coastsat_window_convergence.py --ref-end 2020`.

## Figures to open, in order

1. [`1-rate_profiles/forward_from_1996/window_profiles_overlay_forward_from_1996.png`](1-rate_profiles/forward_from_1996/window_profiles_overlay_forward_from_1996.png)
   shows every window over the reference (y-axis zoomed to ±10 m/yr). The
   `_panels_` twin has one panel per window, with the reference alone in the
   last.
2. [`2-r_bias_rmse/window_profiles_r.png`](2-r_bias_rmse/window_profiles_r.png)
   and [`2-r_bias_rmse/window_profiles_bias_rmse.png`](2-r_bias_rmse/window_profiles_bias_rmse.png)
   score each window against 1996–2024 by window length, with the 15-year
   model window and r = 0.5 / 0.75 / 0.9 marked and every marked point
   labelled with its value.
3. [`3-settling_window/years_needed_alongshore.png`](3-settling_window/years_needed_alongshore.png)
   gives the years each transect needs at three plain thresholds, with the
   model's 15 years marked.
4. [`3-settling_window/forward_from_1996/a-eight_sites/shoreline_position_window_fits_forward_from_1996.png`](3-settling_window/forward_from_1996/a-eight_sites/shoreline_position_window_fits_forward_from_1996.png)
   shows why, at eight places: the raw positions with the fits drawn on them.
5. [`4-split_windows/split_windows_2010.png`](4-split_windows/split_windows_2010.png)
   shows the model's two halves at eight transects picked by behaviour (agree,
   disagree, sign flip, 2021 step) against the 1996–2024 line.

Each direction folder has a `backward` twin. Captions are in each folder's
`supporting/CAPTIONS.md`.

## Producers

```
scripts/input_prep/5-scr/3-rates/coastsat/window_convergence/
    coastsat_window_profiles.py      -> 1-rate_profiles/   (fits every window; slow)
    coastsat_window_r_bias_rmse.py   -> 2-r_bias_rmse/     (reads 1-'s transects CSV; no refit)
    coastsat_window_convergence.py   -> 3-settling_window/
    coastsat_split_windows.py        -> 4-split_windows/   (reads 1-'s transects CSV)
```

Run 1- before 2- after any change to the fits. The layout is resolved in
`scripts/site_layer/hat_observed_rates.py` (`window_profiles_dir`,
`window_scores_dir`, `window_convergence_dir`, `SETTLING_SCALE_DIRS`).

Reorganised 2026-09-29 from `record_<start>_<end>/<direction>/{sites,
all_transects,domain_means,alongshore_profiles}/`, and again 2026-10-01: the
scores moved out of `1-rate_profiles/` into their own `2-r_bias_rmse/`, the
settling sweep was renumbered `2-` → `3-`, and the per-direction
`window_profiles_correlation.csv` (r without intervals) was deleted in favour
of `2-r_bias_rmse/window_profiles_r_bias_rmse.csv`. The fits did not change
(`window_profiles_transects.csv` byte-identical on rerun).
