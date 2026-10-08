# 2026-10-08-schedule-refit

**Question.** Does the blocking groin fit the 1996–2009 calibration period better if it starts failing in 1996 rather than at the 2004 step? The condition analysis (`hard-structures/groin/1-observations/coastsat_groin_condition/`) found that the CoastSat gap across the groin stopped widening in 1995, at the last repair, rather than at Isabel. This study refits strength b and post-failure fraction f under three failure schedules and scores them on that annual CoastSat gap.

**Setup.** `schedule_refit.py` drives the unchanged runner. A child process patches `cascade.groin.BlockingGroinCallback` to the schedule, and each run records its schedule in `schedule.json`. The runner's own report still prints its default schedule.

| schedule | delay from 1969 install | mode | strength each year |
|---|---|---|---|
| `instant2004` | 35 | instant | b to 2003, b × f from 2004 (the runner's schedule, pinned 10-05) |
| `instant1996` | 27 | instant | b × f from 1996, the first year after the 1995 repair |
| `ramp1996` | 26 | linear ramp, 8 yr | b in 1995, falling linearly to b × f in 2003 |

Common to every run:
- full management, the solved edgeBE ends, relocations off;
- the grid is b 0.3–0.9 × f 0.1–0.8, 49 cells per schedule, 147 runs, all clean;
- runs are filed under `output/raw_runs/experiments/groin/2026-10-08-schedule-refit/<schedule>/<b_f>/`.

**Checks.**
- The 2004 schedule at b 0.6/f 0.6 is bit-identical to the 10-05 run.
- The yearly applied strength in `groin_diagnostics.csv` follows each schedule.

**Scores.**
- **Primary: RMSE against the annual CoastSat gap, 1996–2008.** The modelled GIS 5|6 gap change is taken at mid-year (the mean of 1 Jan of that year and the next) and compared with the CoastSat calendar-year mean. Both are relative to the start: the model's t = 0, and the CoastSat mean over the DEM-centred start window (Oct 1995–Oct 1997).
- **Secondary: the 10-05 photo-date RMSE** at 1997, 2004 and 2008.

## Calibration 1996–2009

No groin scores 41.6 m on the annual CoastSat gap and 69.2 m on the photo dates.

| schedule | best on CoastSat (b, f) | annual RMSE | its photo RMSE | best on photos (b, f) | photo RMSE | its annual RMSE |
|---|---|---|---|---|---|---|
| failure at 2004 | 0.4, 0.8 | 10.7 | 26.9 | **0.6, 0.6 (the pin)** | **4.0** | 29.2 |
| failure from 1996 | b × f = 0.36 | **9.8** | 29.7 | 0.9, 0.6 | 9.8 | 27.5 |
| ramp 1996 to 2003 | 0.4, 0.8 | 10.4 | 31.8 | 0.8, 0.6 | 7.4 | 35.9 |

**Reading.**
1. **On CoastSat, the three schedules score about the same (9.8–10.7 m).** Every best cell is a weak groin, at roughly 0.32–0.4 of full blocking, held nearly constant through the window. "Failure from 1996" is that constant-strength family exactly, and it scores best. The 2004 and ramp bests both sit on the f = 0.8 grid edge, which is barely a failure.
2. **CoastSat can't tell when the groin failed; it can tell how strong the groin was.** The annual series falls steadily from 1996, so no fit wants a strong groin before 2004.
3. **The photos want the opposite:** a strong groin (b 0.6) until a sharp failure at 2004. That rests on the 2004 survey, the one date where the photos and CoastSat disagree by 42 m (condition analysis, figure 7). Every CoastSat-best cell misses the photos by 27–32 m, and every photo-best cell misses CoastSat by 28–36 m.
4. **Failure from 1996 can only constrain the product b × f.** Equal products give identical runs (b 0.6/f 0.6 = b 0.9/f 0.4).

## Test 2009–2025 (not fitted)

Every schedule has failed by 2009, so a test run depends only on b × f.

| b × f | stands for | RMSE 2009–2016 | RMSE 2017–2024 | all years |
|---|---|---|---|---|
| none (no groin) | | 31.7 | 9.8 | 23.5 |
| 0.36 | failure-from-1996 best, **and the current pin** | **7.5** | 43.1 | 30.9 |
| 0.32 | 2004-step best, ramp best | 9.1 | 36.9 | 26.9 |

Up to the 2017 Buxton fill, both groins track the test period (7.5–9.1 m against 31.7 m with no groin). After the fill, every groin run is off by 37–43 m. The cause is the fill placement already diagnosed on 10-05: the observed sand moved south past the groin, while the model keeps it at GIS 6. It is not caused by the groin schedule.

**The current pin and the CoastSat best have the same post-failure strength (b × f = 0.36).** Their test runs, and any forward scenario from 2009 on, are identical. Switching to "failure from 1996" changes only how the groin behaves in the calibration years 1996–2003.

## For the decision (Hannah's call)

- **A. Failure from 1996, b × f = 0.36.** This fits CoastSat best and matches the maintenance record (no repair after 1995). In the 2009 test and forward runs it behaves exactly like the current pin. In code, the runner's schedule changes from delay 35 to delay 27, with b and f kept (or any pair with the same product).
- **B. Keep the 2004 step, b 0.6/f 0.6.** This fits the photos best, but CoastSat scores it 29 m. It rests on the single 2004 survey.
- **C. Ramp, b 0.4/f 0.8.** This sits between A and B on both scores. It is on the f grid edge.

**Knock-on effects of A:**
- The per-domain BE set 1 (`domainBE`) was solved with the groin as pinned, so it would need re-solving. The calibration change is local to GIS 2–8.
- The edgeBE ends should be re-checked, as they were on 10-05 after the pin. The GIS 1 residual was −0.04 m then.

## Files

| file | contents |
|---|---|
| `schedule_refit.py` | `grid`, `score`, `test`; the in-process schedule patch |
| `refit_figures.py` | the three figures |
| `grid_scores.csv` | every cell: annual and photo RMSE, bias, end values; both periods |
| `coastsat_gap_change.csv` | the observed annual gap change relative to each start window |
| `grid_run.log`, `score.log`, `logs/` | run and score logs |
| `figures/schedule_refit_1_misfit_maps.png` | annual CoastSat (top) and photo-date (bottom) RMSE over b and f, per schedule |
| `figures/schedule_refit_2_calibration_gap.png` | the gap through the calibration run: best cell per schedule, the pin, no groin, CoastSat and photos |
| `figures/schedule_refit_3_test_gap.png` | the same on the test period, one line per b × f |

Captions are in `figures/supporting/CAPTIONS.md`.
