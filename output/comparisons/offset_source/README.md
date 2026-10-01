# offset_source - dune line or shoreline as BRIE's island offset: how much does it change the model?

Hannah, 2026-09-28: a simple comparison of the model started from the
dune-line offset against the model started from the shoreline offset
(`data/hatteras_init/2-brie-offset/<start>/{duneline,shoreline}/v1`), for
1996-2010 and 2010-2024, mainly to see how much the island's orientation
affects the outcome.

| | |
|---|---|
| script | `scripts/analyze_output/compare_runs/offset_source/offset_source_comparison.py` |
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

## shoreline_v1_vs_v2/ - does the shoreline offset's averaging window matter? (added 2026-10-01)

Hannah, 2026-10-01: the model-output comparison of shoreline offset **v1**
(CoastSat mean over calendar 1995-1997 / 2009-2011) against **v2** (mean over
+/-1 yr of the start DEM's lidar flights, CURRENT since 09-29).

| | |
|---|---|
| script | `scripts/analyze_output/compare_runs/offset_source/shoreline_v1_vs_v2_comparison.py` |
| runs | the full-management `shoreline_v1` / `shoreline_v2` pair of each period in `raw_runs/experiments/island-offset/2026-09-29-shoreline-offset-v1-vs-v2-adopted-setup/` (listed in `tables/summary.csv`); no new runs |
| settings | adopted setup, option A waves, zeroBE ends, relocations and groins off; within a pair only the offset's averaging window differs |
| input-side comparison | the offsets themselves: `data/hatteras_init/2-brie-offset/<year>/shoreline/v2/offset_<year>_v1_vs_v2.png` |

| interior GIS 2-89, v2 minus v1 | 1996-2010 | 2010-2024 |
|---|---|---|
| r between the two runs | 0.99 | 0.99 |
| mean \|difference\| | 0.8 m | 1.4 m |
| largest difference | 3.7 m at GIS 5 | 5.7 m at GIS 79 |
| sd of difference / sd of the model profile | 0.10 | 0.13 |
| island-wide mean difference | -0.1 m | +0.1 m |
| offset difference, sd (max) | 2.8 m (8.4 m) | 4.8 m (12.9 m) |
| r, turning difference vs change difference | +0.74 | +0.74 |
| bias / RMS vs projected target, v1 then v2 | -7.4 / 16.8, -7.4 / 16.8 m | -6.9 / 19.7, -6.9 / 19.7 m |

- The averaging window moves the offset by a few metres alongshore (up to
  8-13 m in one domain). The model passes on about a third of that, at most
  5.7 m. It does not move the mean or the pattern (r 0.99), and the
  projected-target scores agree to 0.1 m.
- As with the dune-line/shoreline pair, what the model responds to is the
  change in **turning**: a new local embayment fills in, and a new local bulge
  is cut back. The change difference is therefore anti-correlated with the
  offset difference itself (r -0.58, -0.59).

```
shoreline_v1_vs_v2/shoreline_v1_vs_v2_model_change_vs_projected_full_management.png   the two runs' change, target on top
shoreline_v1_vs_v2/shoreline_v1_vs_v2_difference_full_management.png        profiles: offset difference, change difference
shoreline_v1_vs_v2/shoreline_v1_vs_v2_offset_vs_model_full_management.png    scatter: offset vs turning as the predictor
shoreline_v1_vs_v2/tables/summary.csv, per_domain.csv, vs_projected.csv
shoreline_v1_vs_v2/supporting/CAPTIONS.md, supporting/*.pdf
```
