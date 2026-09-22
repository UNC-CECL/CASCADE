# hindcast_calibrated - the hindcast against the CoastSat target

## Since 2026-09-18: the 1996 -> 2010 -> 2024 chain, edgeBE and zeroBE

The figure moved onto the canonical chain the same day it was redrawn for the
house-style column. What changed:

* **periods** 1996-2010 and 2010-2024, read from HATTERAS_PERIODS; the target
  is `hat_observed_rates.lrr_csv(start, end)`
* **presets** edgeBE (the default now) and zeroBE. calibBE is not solved on
  this chain, so there is no calibrated figure, and the folder name is now
  historical
* **groin off.** The matrix on this chain has only `nogroin` arms, so the
  figure draws `road_bdm_nogroin` (1996, no fill scheduled) and
  `road_bdm_nourish_nogroin` (2010). Everything below this section describes
  groin-on runs
* **layout.** No title or note on the canvas; that text is the caption. The
  same figure is also written to `output/figures/shoreline/hindcast_<preset>.png`
  with its caption in `supporting/CAPTIONS.md`

```
hindcast_edgeBE_loess_reference.png   edgeBE, both periods, both scoring windows
hindcast_zeroBE_loess_reference.png   the same, no source/sink field at all
```

The last calibBE render (1984/2004) was retired to
`output/archive/2026-09-18_hindcast-calibrated/`; its WHY.md says why.

| | 1996-2010 | 2010-2024 |
|---|---|---|
| edgeBE, LOESS D11-D89 | RMSE 1.12, bias -0.09, r 0.34 | RMSE 2.34, bias -1.58, r 0.08 |
| edgeBE, D2-D89 | RMSE 1.14, bias +0.03, r 0.44 | RMSE 2.31, bias -1.33, r 0.28 |
| zeroBE, LOESS D11-D89 | RMSE 1.04, bias -0.19, r 0.46 | RMSE 2.44, bias -1.86, r 0.10 |
| zeroBE, D2-D89 | RMSE 1.09, bias -0.22, r 0.47 | RMSE 2.61, bias -2.02, r 0.08 |

zeroBE scores better than edgeBE over the interior in 1996-2010, the reverse
of what the edge solve is for. That was not investigated on 09-18.

Runs, all `matrix` rows in `output/raw_runs/run_index.csv`:
`HAT_1996_2010_{edgeBE,zeroBE}_road_bdm_nogroin` and
`HAT_2010_2024_{edgeBE,zeroBE}_road_bdm_nourish_nogroin`.

---

## History: the 1984/2004 figure, groin on (to 2026-09-17)


The headline result, both periods, and its uncalibrated companion. Everything
here is full management (roadway + beach/dune, with nourishment in 2004-2024)
with the groin on at the fitted M = 60 m/yr, f = 0.6; the only thing that
varies between the calibrated pair and the edgeBE pair is the source/sink
field. `../scenario_grid/scenario_grid_by_preset.png` is where presets and management
scenarios are compared properly; these four are a single deliberate contrast.

```
hindcast_calibrated_loess_reference.png   THE figure. Modelled rate against the
                                          CoastSat LOESS reference, with the
                                          geographic annotation layer and BOTH
                                          scoring windows printed per panel
hindcast_edgeBE_loess_reference.png       the same, calibration removed
```

## There used to be two, and why there is now one

A second figure, `HAT_hindcast_final_figure.py`, drew the same two runs scored
over D2-D89 instead of D11-D89. The gap between those two numbers is the D2-D10
strip alone, not the smoothing: north of D10 the LOESS curve and the
calibration target are identical numbers.

Keeping that difference in two files was not free. Each figure quoted the
other's number in its caption, and the headline's value sat hardcoded here for
twelve days after the 1984 run was remade on topography v2 -- the caption said
0.526 / 0.558 while the headline printed 0.547 / 0.580.

So on 2026-09-14 both windows were put on each panel of this figure, computed
rather than quoted, and the headline was retired to
`scripts/figure_making/model_output/superseded_20260914/` (see the WHY.md there). Its
last two PNGs were deleted the same day rather than kept beside their
replacements -- the retired script still runs and regenerates them, so keeping
a copy bought nothing but a second set of numbers to go stale.

The headline's three shaded bands were NOT ported, and that was not an
oversight: they had already been tried on this canvas and removed, because with
the geographic annotation layer in place the Buxton area carried four
overlapping fills and read as clutter. D5-D7 reserved and D1/D90 locked are
stated in the caption instead, and the groin line still marks D5.5. The grey
frozen-zone band is the one fact now drawn in no live figure -- judged the
least load-bearing of the three, and recoverable from the retired script.

Written by `scripts/figure_making/model_output/hindcast_final_figure_loess.py`,
`--preset calibBE` (default) or `--preset edgeBE`. `zeroBE` is wired up but not
built here. PNG only - it does not write a vector copy.

Until 2026-09-18 these carried their title sentence ON the canvas rather than
in a CAPTIONS.md, because they were opened alone, months later, with nothing
around them. The 09-17 move to the 190 mm column broke that layout, and the
title is the caption now.

## What edgeBE is, and why it is the right "uncalibrated"

edgeBE keeps the D1 and D90 edge values and removes every interior correction.
The edges stay because they are solved by buffer-cell reproduction rather than
fitted to a residual - they are boundary absorbers, not a sediment budget - so
removing them would take out something that was never calibrated in the first
place. What edgeBE removes is exactly what the calibration put in, which makes
the gap between the two pairs a fair picture of the job the source/sink field
does. `zeroBE` removes the edges too and is the rawer baseline.

Two things follow in the figures themselves. The grey "outside the frozen zone
set" band is not drawn under edgeBE, because no interior domain is corrected
and the zone set therefore divides nothing; the footnote says so in place of
the band. D5-D7 stay marked, because the groin reservation is independent of
the source/sink field.

| | 1984-2004, calibration | 2004-2024, test |
|---|---|---|
| calibBE, D2-D89 | RMSE 0.547, bias +0.008, r 0.93 | RMSE 0.580, bias +0.104, r 0.87 |
| edgeBE, D2-D89 | RMSE 1.235, bias +0.149, r 0.58 | RMSE 1.779, bias -0.994, r 0.23 |
| calibBE, LOESS D11-D89 | RMSE 0.448, bias +0.077, r 0.95 | RMSE 0.467, bias +0.046, r 0.91 |
| edgeBE, LOESS D11-D89 | RMSE 1.254, bias +0.221, r 0.43 | RMSE 1.843, bias -1.120, r 0.08 |

## The runs behind them

| period | preset | run | topography |
|---|---|---|---|
| 1984-2004 | calibBE | `HAT_1984_2004_calibBE_road_bdm_groin` | 1984-start v2 |
| 2004-2024 | calibBE | `HAT_2004_2024_calibBE_road_bdm_nourish_groin` | 2004-start v1 |
| 1984-2004 | edgeBE | `HAT_1984_2004_edgeBE_road_bdm_groin` | 1984-start v2 |
| 2004-2024 | edgeBE | `HAT_2004_2024_edgeBE_road_bdm_nourish_groin` | 2004-start v1 |

All four are the `arm='calibration'` rows in `output/raw_runs/run_index.csv`.
That qualifier matters: the 1984 calibBE name carries FIVE index rows and only
one is the run in `output/raw_runs/`. The others (`behindroad-copy`,
`version-pair/v2`, `version-pair/v3`, `pea1989basenoreloc`) live in arm
subtrees and differ - the v3 pair gives RMSE 0.544 / bias +0.013 rather than
0.547 / +0.008.

All four are on their product's CURRENT topography, so the calibBE-to-edgeBE
gap is a preset difference and nothing else. It was not always: the edgeBE pair
was built 2026-09-01 on 1984-start v1, and both cells were re-run on
2026-09-14 to close that - the 1984 cell onto v2 (CURRENT), the 2004 cell under
current code, its 2004-start v1 already being the only and current version of
that product. Two of the 93-run v1 backlog, neither on its exclusion list.

Moving the 1984 cell from v1 to v2 changed it about as much as the topography
version usually does: RMSE 1.228 -> 1.235, bias +0.168 -> +0.149 over D2-D89.
The 2004 cell moved by 0.001 m/yr RMSE, a code-only refresh.

## Currency, checked 2026-09-14

The calibrated pair sits on its product's CURRENT topography and the current
calibBE source/sink field, scored with the LRR estimator. Neither run has
relocations enabled, so the cell-rounding change that moved the relocation arm
does not touch them.

Both were made from a dirty tree one to two weeks before that check, so the
code axis was measured rather than assumed - re-run under current code into arm
`currency-20260914` and differenced per domain:

- **calibration period reproduces to exactly zero.** All 90 domains,
  `lrr_m_yr`, `change_rate_m_yr` and `lrr_r2` bit-identical.
- **test period moves 56 of 90 domains by at most 0.0084 m/yr** (largest at
  GIS 15). Skill goes RMSE 0.580353 -> 0.579999, bias +0.103500 -> +0.104062.
  That is ~70x below the RMSE; the only thing it changes on the page is the
  test-period bias in the third decimal, +0.103 -> +0.104.

What produces that residual was not isolated. The calibration run is six days
newer and reproduces exactly, so it is somewhere in the commits between them.

The edgeBE pair does not need that check: both cells were rebuilt from scratch
on 2026-09-14 under current code and CURRENT topography, so they carry no
drift.

A figure here is only as current as the runs behind it, and the run index
records the topography version and git commit of every run.
