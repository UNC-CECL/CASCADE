# 2026-09-25 — four-parameter wave grid, scored on the smoothed model

> **Record (zeroBE, smoothed score).** The settings to use: [`../2026-09-27-wave-recommendation/README.md`](../2026-09-27-wave-recommendation/README.md).

Follow-up to `wave-climate/2026-09-24-metres-2-wave-sensitivity/` (see `2026-09-24-metres-INDEX.md`):
that study searched one parameter at a time plus two 2-D grids, never all four
wave parameters together, and never a grid under full management. Designed
with Hannah on 2026-09-25.

## Design

| | |
|---|---|
| score | share of the alongshore variation explained, 1 − SSE/SST, with the **model smoothed like the CoastSat target** (LOWESS over 10 domains, the southern 10 raw), interior GIS 2–89, against each window's CoastSat LRR target. Bias, RMSE, correlation and the raw (unsmoothed) score beside it |
| coarse grid | Hs 0.75, 1, 1.5, 2 × Tp 7, 8, 10 s × asymmetry 0.5, 0.7, 0.9 × high-angle 0.3, 0.45, 0.55 = 108 per window × scenario, 432 runs (Hs 0.65 and Tp 12 left out: they drowned the barrier in step 2) |
| scope | natural and full management; 1996–2010 first, then 2010–2024 |
| refine | per window × scenario, a 3×3×3×3 grid around the best setting at half the coarse step (the midpoints to its coarse neighbours), launched automatically |
| cross | the top 5 per window × scenario, run in the other window, so a shared setting has candidates |
| shared pick | lowest mean of smoothed RMSE ÷ that window's flat-line RMSE, among settings run in both windows |
| fixed | metres offset (dune line v1), zeroBE, no groin, no relocations; Barrier3D with the route_overwash fix (the driver refuses to run without it) |

The step-2 runs are **not** reused as grid cells: most predate the Barrier3D
fix, and keeping every cell on one model keeps the grid clean. They are
rescored the same way in `tables/step2_rescored_smoothed.csv`.

## Layout

```
README.md
tables/
  all_runs.csv                    every cell: settings, smoothed and raw scores, status, run folder, log
  best_settings.csv               best per window x scenario, and the shared pick per scenario
  observed_targets.csv            each window's observed mean and flat-line RMSE
  step2_rescored_smoothed.csv     the 2026-09-24 step-2 runs, scored the same way
figures/                          (after the runs)
logs/<phase>_<scenario>/<period>/<settings>.log, logs/drivers/
runs/<phase>_<scenario>/<period>/zeroBE/<run_name>/     on disk only; phase = coarse, refine, cross
```

## Reproduce

```
python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py run all
python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py score
python scripts/hatteras_ms/experiments/HAT_wave_grid_smoothed_score.py rescore-step2
```

## Status

Launched 2026-09-25 11:10 (8 at a time from 11:22). **Paused at 14:01 after the
1996–2010 coarse grid** (Hannah), before 2010–2024: 216 cells, 198 scored, 18
drowned in year 4 (every Hs 0.75 / Tp 10 cell, both scenarios). Resume with
`run all --jobs 8`: it skips what has run.

## 1996–2010 coarse grid

| | best setting | smoothed | raw | bias (m/yr) |
|---|---|---|---|---|
| natural | Hs 1.5, Tp 10, asym 0.7, high-angle 0.45 | +23% | +8% | −0.51 |
| natural, #2 | Hs 1.0, Tp 7, asym 0.9, high-angle 0.45 | +23% | +14% | −0.48 |
| full management | Hs 1.5, Tp 10, asym 0.7, high-angle 0.45 | +23% | +12% | −0.19 |
| full management, #2 | Hs 1.0, Tp 8, asym 0.9, high-angle 0.45 | +22% | +16% | −0.05 |

The four-parameter grid reaches the same ceiling as step 2 (rescored: +24%
natural, +23% managed, both at asymmetry 0.8, which the coarse grid does not
contain; the refine step adds it). High-angle 0.45 is in nearly every top
setting, and larger Hs pairs with longer Tp. The top profiles all have one
shape (`figures/best/grid/top5_smoothed_profiles_grid.png`): near −0.5 m/yr
through the middle of the island and an erosion trough at Tri-Village, and
none of the observed accreting peaks at GIS 18, 29 and 42, so no combination
of the four wave parameters produces them.

## Final (sweep finished 2026-09-25, ~21:00)

578 runs: coarse 432, refine 136, cross 10; 18 drowned (every Hs 0.75 / Tp 10
cell), no crashes. `tables/best_settings.csv`; figures in `figures/best/grid/`.

| | best setting (smoothed score) | smoothed | bias (m/yr) |
|---|---|---|---|
| natural 1996–2010 | Hs 1.25, Tp 10, asym 0.8, high-angle 0.45 (refine) | +24% | −0.39 |
| managed 1996–2010 | Hs 1.5, Tp 10, asym 0.7, high-angle 0.45 | +23% | −0.19 |
| natural 2010–2024 | Hs 2.0, Tp 10, asym 0.5, high-angle 0.50 | −530% | −3.42 |
| managed 2010–2024 | Hs 2.0, Tp 7, asym 0.5, high-angle 0.55 | −124% | −1.54 |

- **1996–2010 is a plateau.** The top five in each scenario lie within about
  1 point (+22–24%), along a ridge from Hs 1 / Tp 7 to Hs 1.75 / Tp 10 with
  high-angle 0.45; the full search confirms step 2's best (natural Hs 1.25 /
  Tp 10 / asym 0.8) rather than finding a better one.
- **2010–2024 runs to the edges of the grid** (Hs 2, the largest; asymmetry
  0.5, the smallest; high-angle 0.5–0.55, the largest) and still explains
  far less than a flat line. The search is pushing the waves to reduce the
  model's erosion, which the waves cannot fix (see the 2010–2024 diagnosis:
  the 2021 CoastSat step, dune gaps, storms). The best managed run is better
  than step 2's (−124% against −175%) but is not a fit.
- **The shared rule picks a 2010–2024 setting.** Natural Hs 2 / Tp 10 /
  asym 0.5 / high-angle 0.5 (1996–2010 +16%); managed Hs 2 / Tp 7 / asym 0.5 /
  high-angle 0.55 (1996–2010 +8%, bias +0.28). The mean of RMSE ÷ flat line
  is dominated by 2010–2024, so the pick gives up most of the 1996–2010 fit
  for a window no setting can match: not a usable compromise as it stands.
- No setting in either window makes the accreting peaks at GIS 18, 29, 42.
