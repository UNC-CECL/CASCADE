# end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-split12

**Question.** After the storms changed to `v3_split12_trim24` (2026-09-29), what end rates at GIS 1 and 90 make the model match the full-period 1996–2024 CoastSat LRR in each window? This solve supplies the `ends_solved_on_coastsat` set of `target_comparison/projected/`, through `FULL_SOLVE_DIR`. (Hannah, 2026-09-29: "re-solve the 1996-2024 LRR ends and redraw target_comparison".)

**Method.** Unchanged from `../2026-09-29-ends-solved-on-lrr-1996-2024-dunecap/`, except that it is capped at 6 steps. The last solve stalled at the 0.1 m/yr probe resolution and had to be stopped by hand.

    be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-split12
        --windows 1996 2010 --target coastsat --coastsat-window 1996_2024 --tol 0.02 --max-steps 6

It starts from the re-run matrix (split12 storms).

**Answer** (`loop_log.csv`; `target_comparison.py` reads the last step, 6):

| window | ends, GIS 1 / 90 (m/yr) | residuals (m/yr) | trim24 storms (09-29 dunecap solve) |
|---|---|---|---|
| 1996–2010 | +4.0 / +27.7 | +0.022 / 0.000 | +3.9 / +27.7 |
| 2010–2024 | +3.6 / +17.5 | −0.030 / −0.008 | +3.6 / +17.5 |

Both chains finished just outside 0.02 at GIS 1, which is the 0.1 m/yr probe resolution again. 1996 repeated the same probe from step 3. 2010 overshot GIS 90 to +18.7 at step 4 (+0.187) and came back to +17.5.

**Runs.** `coastsat/step<k>/<window>/edgeBE/<run>/` stay on disk only (`.gitignore`). The driver output is `output/logs/driver/lrr1996_2024_endsolve_split12_20260929.log`.
