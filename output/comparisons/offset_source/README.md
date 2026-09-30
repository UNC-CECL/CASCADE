# offset_source - dune line or shoreline as BRIE's island offset: how much does it change the model?

Hannah, 2026-09-28: a simple comparison of the model started from the
dune-line offset against the model started from the shoreline offset
(`data/hatteras_init/2-brie-offset/<start>/{duneline,shoreline}/v1`), for
1996-2010 and 2010-2024, mainly to see how much the island's orientation
affects the outcome.

| | |
|---|---|
| script | `scripts/analyze_output/compare_runs/offset_source_comparison.py` |
| runs | the full-management pair of each period in `raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/` (listed in `tables/summary.csv`); no new runs |
| settings | option A waves, zeroBE ends, relocations and groins off; within a pair only the offset differs |
| model change | each run's LRR x 14 yr (m, seaward positive) |
| orientation | atan of the alongshore slope of the seaward position (degrees); **turning** is its rate of change alongshore (degrees per km, positive at an embayment) |
| statistics | interior GIS 2-89 |

## Results (2026-09-28)

| | 1996-2010 | 2010-2024 |
|---|---|---|
| r between the two runs | 0.88 | 0.95 |
| mean \|difference\| | 3.3 m | 3.1 m |
| largest difference | 18.2 m at GIS 80 | 14.8 m at GIS 23 |
| sd of difference / sd of the model profile | 0.48 | 0.32 |
| island-wide mean difference | +0.3 m | +0.1 m |
| orientation difference, sd (max) | 1.1° (3.2°) | 1.1° (2.9°) |
| r, orientation difference vs change difference | -0.08 | -0.11 |
| r, turning difference vs change difference | **+0.81** | **+0.64** |

- The source moves the model locally, by about 3 m on average and up to
  15-18 m, against a profile whose own spread is 10-14 m. It does not move
  the island-wide mean (under 0.5 m) and barely moves the pattern (r 0.88
  and 0.95).
- The two offsets differ in orientation by about 1° (at most 3°). That
  angle difference predicts almost none of the change difference.
- What predicts it is the difference in **turning**: where one offset has a
  local embayment the other does not, that run fills it in, and where it
  has a bulge, that run cuts it back. BRIE's alongshore transport diffuses
  the planform, so the domain-scale wiggles in each line are what the model
  responds to, not the overall lean.
- The largest differences are in Tri-Village (GIS 77-82), Avon Shoals
  (GIS 28-35) and the north end.

## Against the projected target (added 2026-09-29)

Hannah: the change figure again with ONE target on top, projected shoreline
change (CoastSat LRR 1996-2024, LOWESS 7 domains, x 14 yr), the same profile
in both periods; model unsmoothed. Same script, `projected_target()`.

| model minus target, interior GIS 2-89 | bias | RMS residual | explained | r |
|---|---|---|---|---|
| 1996-2010 dune line (1997) | -7.0 m | 16.7 m | -7% | 0.40 |
| 1996-2010 shoreline (1995-1997 mean) | -6.6 m | 16.3 m | -2% | 0.42 |
| 2010-2024 dune line (2009) | -11.0 m | 22.1 m | -87% | 0.20 |
| 2010-2024 shoreline (2009-2011 mean) | -10.9 m | 22.4 m | -92% | 0.13 |

- Against this target the source barely matters in either period: bias and
  RMS differ by under 0.5 m, and both runs miss the same things -- the
  accreting Avon Shoals (GIS 28-35) and Wimble Shoals (GIS 64-72) bulges and
  the Cape Point accretion; in 2010-2024 both also over-erode Tri-Village.
- The 1996-2010-only two-panel version drawn earlier the same day
  (`HAT_offset_source_comparison_option_a.py plot-projected`) is deleted;
  this figure replaces it.

```
offset_source_model_change_full_management.png          (a, b) only: the two runs' change, both periods
offset_source_model_change_vs_projected_full_management.png   the same with the projected target on top
offset_source_difference_full_management.png            profiles: change, orientation difference, change difference
offset_source_orientation_vs_model_full_management.png   scatter: angle vs turning as the predictor
tables/summary.csv, tables/per_domain.csv, tables/vs_projected.csv
supporting/CAPTIONS.md, supporting/*.pdf
```
