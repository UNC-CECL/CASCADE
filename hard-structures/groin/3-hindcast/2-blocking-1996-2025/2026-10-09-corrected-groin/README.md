# 2026-10-09-corrected-groin

**Question.** Does the corrected groin (straight-coast module 6: the drift groin plus a correction for BRIE's sand leak at the groin's bends) fit the 1996-2009 gap across the groins in the full model, and how does it do on 2009-2025?

**Why.** The drift groin (`../2026-10-08-drift-groin/`) fitted better than the pin, but leaked about 107 m of summed shoreline onto the island through BRIE's solve. Module 6's correction removes that leak on a straight coast (`../../../2-module-tests/`). As built, though, the correction also acted on Buxton's natural bend, changing the no-groin baseline by about 100 m (`../../../2-module-tests/3-real-planform/corrected_groin/`).

**Setup.** `corrected_groin.py` replaces `BlockingGroinCallback` with the corrected groin in a child process and drives the unchanged runner, as the drift-groin study does: full management, the solved edgeBE ends, relocations off, instant failure at the 2004 step (b falls to b·f). Each year:
- **The groin** moves b of BRIE's wave-climate net drift at the updrift approach face from GIS 5 to GIS 6, as `drift_groin.py` does.
- **The correction (option B, chosen 2026-10-09)** is applied to the groin's own bends only. Within 10 cells of the groin face, it is module 6's leak term for this run's BRIE shoreline minus the same term for the no-groin run's in that year. The no-groin BRIE shoreline is rebuilt as this run's year-0 BRIE shoreline plus the no-groin run's change since year 0. The zone is 10 cells, not module 6's 5, because 5 leaks on a grid-aligned coast (`../../../2-module-tests/6-long-term/`).
- **Check (`check`):** at b 0 the run reproduces the no-groin run to 5 × 10⁻¹³ m in every year and domain. The correction applied was at most 2 × 10⁻¹⁴ m, and the rebuilt no-groin shoreline matched the run's own to 9 × 10⁻¹³ m. So option B leaves the baseline untouched.

Each run dir has `corrected_groin.txt`, because the runner's report and run name (`groinblock`) do not show the patch. The diagnostics CSV gains `net_drift_m3_yr`, `approach_angle_deg`, `leak_correction_m` and `max_diff_from_no_groin_m`. Runs are under `output/raw_runs/experiments/groin/2026-10-09-corrected-groin/`.

    python corrected_groin.py check
    python corrected_groin.py grid [--streams 4]     # b 0.02-0.3 x f 0.2-0.8, 32 runs
    python corrected_groin.py edge [--streams 4]     # f 0.9-1.0 at b 0.15-0.3, and b 0.25; 11 runs
    python corrected_groin.py score                  # -> grid_scores.csv (corrected, drift and no groin; joint = RMS of the two scores)
    python corrected_groin.py test --b 0.2 --f 0.8

## Calibration, 1996-2009

All 43 runs clean (`grid_run.log`, `edge_run.log`). Scores are as the drift-groin study: annual CoastSat RMSE and photo-date RMSE of the GIS 5|6 gap change.

| run | annual RMSE | photo RMSE | joint | island vs no groin |
|---|---|---|---|---|
| no groin | 41.6 | 69.2 | 57.1 | 0 |
| pinned b 0.6 f 0.6 (the 10-05 pin) | 29.2 | 4.0 | 20.9 | −22 m |
| drift b 0.1 f 0.4 (10-08 best) | 19.0 | 11.7 | 15.8 | +107 m |
| **corrected b 0.2 f 0.8** (best with the 2004 failure in the grid) | 19.6 | 21.1 | **20.4** | +19 m |
| corrected b 0.2 f 0.9 / 1.0 | 20.6 / 21.6 | 18.7 / 17.6 | 19.7 / 19.7 | +20 m |

1. **It fits about as well as the pin, and conserves sand.** Its best joint score is 20.4 m, against the pin's 20.9 m. Across the whole grid the island gains only 3 to 21 m of summed shoreline against no groin, against the drift groin's 19 to about 200 m. Part of the drift groin's better score (15.8 m) came from that leaked sand piling up updrift.
2. **It needs about twice the drift groin's b** (0.2 against 0.1), because without the leak the same b holds a smaller gap. At b 0.05 f 0.8 the corrected groin scores 28.0 / 55.2 m, the drift groin 13.7 / 37.6 m.
3. **The two scores still want different groins:** the annual CoastSat score is best at b 0.15 with little failure (11.6 m at f 1.0), the photo score at b 0.25 f 0.9 (14.0 m).
4. **The calibration barely sees the 2004 failure.** The best cells sit at f 0.8 to 1.0, and removing the failure altogether (f 1.0) gains 0.7 m.

## Test, 2009-2025

| run | annual RMSE | photo RMSE | joint | gap change at the end, model / CoastSat |
|---|---|---|---|---|
| no groin | 23.5 | 31.1 | 27.6 | −57.8 / −50.4 |
| drift b 0.1 f 0.4 | 26.2 | 31.7 | 29.0 | −17.5 / −50.4 |
| **corrected b 0.2 f 0.8** | 38.7 | 43.3 | **41.1** | **−2.8** / −50.4 |
| corrected b 0.2 f 1.0 | 48.5 | 53.4 | 51.0 | +9.1 / −50.4 |

5. **Out of sample it does worse than no groin and than the drift groin.** It holds the gap almost unchanged after 2009 (−2.8 m) while CoastSat's closes by 50 m. After the 2004 step it still blocks b·f = 0.16 of the drift, four times the drift groin's 0.04. The calibration chose a groin that hardly weakened in 2004, and the later record needs one that largely gave way. This is the failure-schedule conflict open in `../2026-10-08-schedule-refit/`, sharper here because the corrected groin holds more sand. The 2017 fill placement also costs every groin run 37-43 m in this period (10-05 README).

**What this means.** The corrected groin's mechanism carries over to Buxton: with the correction confined to its own bends, it traps on the right side, conserves sand and leaves the no-groin baseline unchanged. It fits 1996-2009 as well as the pin. It does not fix the test period, because what the groin does after 2004 is set by the failure schedule, not by how it traps. The next step is the schedule decision (options A-C in `../2026-10-08-schedule-refit/`), then refitting this module under the chosen schedule.
