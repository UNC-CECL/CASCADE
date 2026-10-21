# 2026-10-08-schedule-refit

**Question.** Does the blocking groin fit 1996–2009 better if it starts failing in 1996, after the last repair in 1995, rather than at the 2004 step?

**Answer.** On the annual CoastSat gap the three schedules score about the same (9.8–10.7 m RMSE). Every best cell is a weak groin of about the same strength, so CoastSat can't tell when the groin failed. The photos still prefer the 2004 step at b 0.6 / f 0.6. The current pin and the CoastSat best share b × f = 0.36, so their 2009–2025 runs are identical. The choice between them is Hannah's and is still open.

**Write-up.** The scripts, scores (`grid_scores.csv`), figures and the options for the decision are in `hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-08-schedule-refit/README.md`. This folder holds only the runs.

**Folder names.**
- `instant2004/`: full strength b to 2003, b × f from 2004 (the runner's own schedule).
- `instant1996/`: b × f from 1996.
- `ramp1996/`: b in 1995, falling linearly to b × f in 2003.
- `<schedule>/b<b>_f<f>/`: one cell of the b 0.3–0.9 × f 0.1–0.8 grid, 49 per schedule.

Each run's `schedule.json` records its schedule. The runner's own report in the run metadata still prints the default (2004) schedule, because the schedule is patched in a child process.

**Status.** Current, decision pending.
