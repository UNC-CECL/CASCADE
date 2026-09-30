# 2026-09-28 — dune line vs shoreline as the island offset, option A waves, 1996–2010

Hannah, 2026-09-28: does the 09-25 result
(`../2026-09-25-metres-offset-duneline-vs-shoreline-waves-hs1-tp8/`) hold at the waves
adopted on 09-27?

| | |
|---|---|
| offsets | `2-brie-offset/1996/duneline/v1` and `1996/shoreline/v1`, metres. No 2010 shoreline offset exists, so 1996–2010 only |
| waves | option A: Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle 0.5. One setting, no sweep (Hannah's choice) |
| ends | zeroBE, as on 09-25. Option A's edgeBE ends were solved on the dune-line offset and would favour it |
| scope | natural and full management; 4 runs, relocations and groins off, Barrier3D with the route_overwash fix |
| score | **raw** share of the alongshore variation explained, interior GIS 2–89 (the score option A was chosen on); smoothed, bias and r beside it |

## Layout

```
tables/all_runs.csv            every run: offset, scenario, scores, status, run folder
figures/rates_duneline_vs_shoreline_1996_2010_option_a.png
                               modelled rate, each offset, against CoastSat
figures/total_change_difference_shoreline_minus_duneline_full_management.png
                               where the two offsets disagree, per domain, by period (m)
logs/<offset>_<scenario>/<settings>.log, logs/driver.log
runs/<offset>_<scenario>/1996_2010/zeroBE/<run_name>/     on disk only
```

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py run --jobs 4
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py plot
```

## Results

Run 2026-09-28; all 4 scored, about 6 min each.

| | raw explained | smoothed | bias (m/yr) | raw r |
|---|---|---|---|---|
| dune line, natural | +19.8% | +24.8% | −0.24 | 0.51 |
| shoreline, natural | **+22.6%** | +25.9% | −0.21 | 0.51 |
| dune line, managed | +18.2% | +22.4% | +0.03 | 0.45 |
| shoreline, managed | **+22.4%** | +22.1% | +0.05 | 0.48 |

- **On the raw score the shoreline offset is ahead in both scenarios**:
  +3 points natural, +4 managed. On the smoothed score they are level
  (±1 point), so the shoreline's gain is at the domain scale, not in the
  broad pattern.
- The rates differ by 0.24–0.28 m/yr on average and by up to 1.3–1.4 m/yr,
  largest at GIS 80 (Tri-Village), then GIS 32–35 and the north end.
- Bias is the same to 0.04 m/yr.
- Neither offset reproduces the accreting peaks at GIS 18, 29 and 42.
- Against 09-25 (Hs 1.0, Tp 8, asym 0.8): the managed tie holds or tips to
  the shoreline; the natural gap narrows (+8 points smoothed then, +1 now).

## Each offset graded on its own feature (added 2026-09-28)

Hannah: a run started from the dune line should be graded on the dune line's
change, a run started from the shoreline on the shoreline's. Full management,
both periods; `run-2010` added the 2010–2024 pair on the new
`2-brie-offset/2010/shoreline/v1` (mean CoastSat shoreline 2009–2011).
Dune line graded on net change in metres (observed 1997→2009 = 11.6 yr and
2009→2023 = 14.1 yr, against the model's 14 calendar years, not corrected);
shoreline on total shoreline change in metres: each period's OWN CoastSat LRR x 14 yr
(1996-2010 LRR for 1996-2010, 2010-2024 LRR for 2010-2024; not the 1996-2024 rate
projected), against the model's own LRR x 14 yr (Hannah's choice). Scaling by 14
leaves the score as it was on the rates. The two are on
different targets and estimators: each ranks an offset against its own
feature, not against the other.

```
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py run-2010
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py grade
python scripts/hatteras_ms/experiments/HAT_offset_source_comparison_option_a.py plot-own
```

`tables/graded_on_own_feature.csv`;
`figures/duneline_offset_vs_duneline_change_full_management.png`,
`figures/shoreline_offset_vs_coastsat_total_change_full_management.png`
`figures/shoreline_offset_vs_coastsat_projected_change_full_management.png`
(third comparison, 2026-09-28: the same shoreline runs against PROJECTED
shoreline change, the CoastSat LRR fitted on 1996-2024 x 14 yr, LOWESS 7 --
one observed profile in both panels; scored in `graded_on_own_feature.csv`)
`figures/total_change_difference_shoreline_minus_duneline_full_management.png`
(they replace the 1996 natural/managed rates and difference figures).
The two starts' total change differs by 3.4 m on average in 1996–2010 (18 m at
GIS 80) and 3.1 m in 2010–2024 (15 m at GIS 23).

Both observations smoothed at **7 domains** (LOWESS, southern 10 raw; the
group's smoothing range, Hannah 2026-09-28); the model unsmoothed. This
study sets 7 itself; the runner's own target is still 10.

| offset, graded on | 1996–2010 explained | bias | 2010–2024 explained | bias |
|---|---|---|---|---|
| dune line, dune-line net change (LOWESS 7) | −38% | +9.6 m | −86% | −10.5 m |
| shoreline, CoastSat total change (own LRR × 14 yr, LOWESS 7) | +21% | +0.7 m | −137% | −24.5 m |

Smoothing the dune observation LOWERED its scores (−11% → −38%, −46% → −86%
unsmoothed): the smoothed line has less variance to explain, while the
unsmoothed model's domain-scale swings still count as error.

- The dune-line runs explain none of the dune line's own pattern in either
  period, smoothed or not.
- 1996–2010: the dune line retreats island-wide while the model's shoreline
  mostly holds, hence the +10.7 m bias. That island-wide retreat may be an
  imagery artefact (see output/comparisons, net change shoreline vs dune line).
- 2010–2024 shoreline: the known failure. CoastSat's target is dominated by
  the 2021 step; the model erodes.

## How this relates to the other study (added 2026-09-28)

`output/comparisons/target_comparison/` and `output/raw_runs/experiments/island-offset/2026-09-28-metres-offset-duneline-vs-shoreline-waves-option-a/` ask different questions of the same option A model.

- **target_comparison** holds the model fixed (dune-line offset) and changes the ruler: which observation should CASCADE be graded against, CoastSat or the dune line? Every run is graded against both targets. There are three sets of end rates (solved on CoastSat, solved on the dune line, unsolved). Scores are bias, RMSE and r in metres over 14 yr.
- **the offset-source study** changes the model's starting island: is the BRIE offset built from the dune line or from the CoastSat shoreline? Every run uses unsolved (zeroBE) ends, and each run is graded on its own feature: the dune-line offset against the dune line, the shoreline offset against CoastSat. Scores are the share of variation explained, then bias and r. Observations are smoothed at 7 domains and the model is left unsmoothed.

| | target_comparison | offset-source study |
|---|---|---|
| varies | target, and which target the ends were solved on | offset source (model input) |
| offset | dune line in every run | dune line vs shoreline (1996/ and 2010/shoreline/v1) |
| ends | CoastSat-solved, dune-solved, unsolved | unsolved (zeroBE) only |
| graded against | both targets | each run's own feature |

**Where they overlap.** The offset study's dune-line, full-management run is effectively the same model as target_comparison's `ends_unsolved` set (option A, zeroBE, dune-line offset). Graded against the dune line in metres, the two agree up to one scaling choice:

| model minus dune line | offset-source study | target_comparison (smoothed_lowess7_with_cascade) |
|---|---|---|
| 1996–2010 | +9.6 m | +10.7 m |
| 2010–2024 | −10.5 m | −10.7 m |

- **Why 1996–2010 differs.** The 1997→2009 dune lines span 11.6 yr. target_comparison scales that change to 14 yr; the offset study compares it unscaled. 11.6 → 14 is a factor of 1.2, which is the 1 m gap. 2009→2023 spans 14.1 yr, so 2010–2024 nearly agrees.
- **The shoreline gradings are not directly comparable.** The offset study compares the model's LRR × 14 with CoastSat LRR × 14; target_comparison uses the model's end-minus-start change.

**Consistent findings:**
- **1996–2010:** the model sits roughly on CoastSat but about 10 m seaward of the dune line. The dune line retreats island-wide while the model holds, and that retreat may be an imagery artefact.
- **2010–2024:** the model erodes past CoastSat (about −24 m against the window's own rate, the 2021 step).
- **Only the offset study adds:** a shoreline-built offset fits slightly better at the domain scale (+3–4 points raw, level when smoothed). So the offset source is a small lever and not the cause of the 2010–2024 failure.
