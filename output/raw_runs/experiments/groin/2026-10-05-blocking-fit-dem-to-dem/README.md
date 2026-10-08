# 2026-10-05-blocking-fit-dem-to-dem

**Question.** Which groin strength fits the GIS 5|6 gap through the 1996–2009 calibration period?

**Answer.** Blocking b 0.6 / f 0.6: photo-date RMSE 4.0 m against 69.2 m with no groin. f is loose (0.5–0.8). Pinned in `output/calibration/groin/joint_fit.json`.

**Write-up.** The scripts, scores (`grid_scores.csv`), figures and logs are in `hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-05-blocking-fit-calibration/`. This folder holds only the runs.

**Folder names.**
- `b<b>_f<f>/`: blocking groin. b is the fraction of the GIS 5|6 exchange blocked while the groin is intact; f is the fraction of b left after it fails at the 2004 step. The grid is b 0.3–0.9 × f 0.1–0.8.
- `dipole_M<M>_f<f>/`: the dipole groin, for comparison.
- Inside each cell, `1996_2009/` is the calibration run. Only b 0.6 with f 0.5, 0.6 and 0.8, and the two dipole cells, also have a `2009_2025/` test run.

**Status.** Current. The 2026-10-08 schedule refit (`../2026-10-08-schedule-refit/`) re-scores this grid against the annual CoastSat gap; the choice between the two is open.
