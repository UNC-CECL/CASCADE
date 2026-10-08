# 2026-10-08-conserving-groin

**Question.** The pinned blocking groin (b 0.6, f 0.6) does not conserve sand at its face. Does a groin that does still fit the 1996-2009 gap, and with what b and f?

**Why.** `BlockingGroinCallback` cancels a fraction b of each cell's coupling to the face, and each cell uses its own BRIE diffusion number. At Buxton GIS 5's is about 4 times GIS 6's, so each year GIS 5 loses about 22 m while GIS 6 gains about 5 m. Over the calibration run the pin applies +60 m seaward at GIS 6 and 246 m landward at GIS 5, a net loss of 186 m of shoreline (123 m over 2009-2025). That loss comes from BRIE's row-scaled diffusion, which is itself non-conservative across the face; the groin cancels b of an exchange that was already uneven.

**Setup.** `conserving_groin.py` drives the unchanged runner with `_r_ipl` patched in a child process, so both sides use one coefficient, r_face: BRIE's diffusion number for the shoreline angle across the GIS 5|6 face (the lower cell's forward difference). Each year `dx_GIS6 = -dx_GIS5`. Everything else is the 10-05 setup: full management, the solved edgeBE ends, relocations off, instant failure at the 2004 step. Grid: b 0.05-0.6 x f 0.1-0.8 (63 runs). Runs are under `output/raw_runs/experiments/groin/2026-10-08-conserving-groin/`, and each run dir has a `conserving.txt` because the runner's report does not show the patch.

**Scores.** `grid_scores.csv` holds the 10-05 photo-date RMSE and the 10-08 annual CoastSat RMSE of the GIS 5|6 gap change, the groin's applied budget (`applied_net_m`, which must be 0 here), and the seaward change at GIS 3, 4, 7, 8. It holds the conserving arm, the pinned arm (10-05 runs) and the no-groin baseline.

## Smoke cell, b 0.6 f 0.6 (calibration)

| run | photo RMSE | annual RMSE | gap 1997 / 2004 / 2008 | GIS 6 / GIS 5 seaward | GIS 3 / 4 seaward | net applied |
|---|---|---|---|---|---|---|
| no groin | 69.2 | 41.6 | −20 / −74 / −86 | −53 / +36 | +9.8 / +26.0 | — |
| pinned | 4.0 | 29.2 | +4 / +13 / −16 | −21 / −2 | +4.6 / +12.1 | −186 m |
| conserving | 77.0 | 111.3 | +19 / +120 / +72 | +112 / +45 | +6.8 / +24.4 | 0.0 m |

Observed gap change from the 1996 start: +2 / +16 / −10. With the exchange balanced, the same b traps about four times more at GIS 6 (21-42 m/yr), and r_face rises from 0.067 to 0.091 as the step grows, so b 0.6 overshoots by about 100 m at 2004. The pin's fit relied partly on the downdrift sink. The grid was moved down to b 0.05-0.6.

## The three modules side by side

`compare_groin_modules.py` animates the dipole (M 12, f 0.3), the pinned blocking groin and the conserving blocking groin year by year on 1996-2009: shoreline change along GIS 1-12, what each module applied at GIS 5 and 6 that year with the running net, and the gap against the photos. Output: `figures/groin_module_comparison_1996_2009.gif` and the final frame as `.png`. It draws the conserving column at b 0.6 f 0.6 for now; after the grid, re-draw it at the best cell with `--conserving b<b>_f<f>`.

## Results

**Grid stopped before it ran (2026-10-08).** The premise was wrong. The pinned groin's imbalance offsets BRIE's own non-conservative solve, so at system level it is the groin that conserves sand: over 1996-2009 the island's total change against no groin is −22 m for the pinned groin, +186 m for this conserving variant and +138 m for the dipole (all 120 domains). The "−186 m" above is the module's own books, not the island's. See `../../../2-module-tests/1-straight-coast/README.md`, findings 3-4. Only the b 0.6 f 0.6 smoke run exists; relaunch with `python conserving_groin.py grid` if it is ever wanted.
