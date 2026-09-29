# end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-dunecap

**Question.** After the dune-cap fix (2026-09-28, the other session), what end rates at GIS 1 and 90 make the model match the full-period 1996–2024 CoastSat LRR in each window?

**Why again.** It is the same reason as the dune-line re-solve beside it (`../2026-09-29-ends-solved-on-duneline-dunecap/NOTE.md`). The fix reran every beach/dune-managed run. That changes the interior at Buxton, Avon and Tri-Village, not the end domains. (Hannah, 2026-09-29: "Re-solve both, then redraw".) It supplies the `ends_solved_on_coastsat` set of `target_comparison/projected/`, through `FULL_SOLVE_DIR`.

**Method.** Unchanged from 09-28:

- `be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-29-ends-solved-on-lrr-1996-2024-dunecap --windows 1996 2010 --target coastsat --coastsat-window 1996_2024 --tol 0.02 --max-steps 8`.
- The target is the 1996–2024 table, with GIS 1 against the raw mean and GIS 90 against the LOESS-7 value.
- It starts from the fixed matrix runs.

**Answer** (`loop_log.csv`; `target_comparison.py` reads the last step, 6):

| window | ends, GIS 1 / 90 (m/yr) | residuals (m/yr) | 09-28, before the fix |
|---|---|---|---|
| 1996–2010 | +3.9 / +27.7 | −0.047 / +0.001 | +3.9 / +27.7 |
| 2010–2024 | +3.6 / +17.5 | +0.038 / +0.008 | +3.6 / +18.6 |

**Both chains stalled just outside 0.02.**

- The cause is the 0.1 m/yr probe resolution, as on 09-27 and 09-28.
- From step 4 (1996) and step 5 (2010) the solver printed no further step, so steps 5–6 re-ran the same probe.
- The loop was stopped by hand during step 7 instead of repeating that probe to step 8. The two unfinished step-7 probes were deleted.

**Runs.** `coastsat/step<k>/<window>/edgeBE/<run>/` stay on disk only (`.gitignore`). Logs are in `logs/`.
