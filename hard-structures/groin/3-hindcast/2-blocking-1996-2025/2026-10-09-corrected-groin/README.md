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

## The failure schedules, A-D (2026-10-09)

Asked the same day: refit under every schedule option of `../2026-10-08-schedule-refit/`, plus a fourth. The schedules are taken from `schedule_refit.SCHEDULES`; D is new here.
- **A, failure from 1996:** one strength all window, run as b with f 1.0 (`instant1996`, 9 runs).
- **B, the 2004 step:** the grid above (`instant2004`).
- **C, ramp:** full in 1995, falling to b·f by 2003 (`ramp1996`, b 0.1-0.4 × f 0.2-0.8, 24 runs).
- **D:** any of these, with the groin passing everything from the 2017 Buxton fill on (`OFF_YEAR`). It is the hypothesis that a 2 million m³ fill buried or overtopped the groins. It changes only the test period. Each schedule's best calibration cell is run on 2009-2025 with and without it.

All 39 runs clean (`schedules_run.log`, `tests_run.log`). The test RMSE is also split at the fill: 2009-2016 and 2017-2024.

**Calibration, best joint cell per schedule:**

| schedule | b, f | annual | photos | joint |
|---|---|---|---|---|
| B, 2004 step | 0.2, 0.9 | 20.6 | 18.7 | 19.7 |
| A, from 1996 | 0.2 (constant) | 21.6 | 17.6 | 19.7 |
| C, ramp | 0.2, 0.8 | 16.1 | 25.1 | 21.1 |

6. **The calibration cannot tell the schedules apart.** All three best cells score 19.7-21.1 m, all at b 0.2 with little or no weakening, as the pin study found.

**Test, 2009-2025:**

| run | strength after failure | annual | photos | joint | 2009-2016 | 2017-2024 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|---|
| no groin | 0 | 23.5 | 31.1 | 27.6 | 31.7 | 9.8 | −57.8 / −50.4 |
| drift b 0.1 f 0.4 | 0.04 | 26.2 | 31.7 | 29.0 | 11.9 | 35.0 | −17.5 / −50.4 |
| corrected B or C (b 0.2 f 0.8) | 0.16 | 38.7 | 43.3 | 41.1 | 18.4 | 51.5 | −2.8 / −50.4 |
| corrected A (b 0.2) | 0.2 | 48.5 | 53.4 | 51.0 | 27.2 | 63.0 | +9.1 / −50.4 |
| **corrected B or C + D** | 0.16, off from 2017 | **16.8** | **20.1** | **18.5** | 18.4 | **15.0** | **−53.3** / −50.4 |
| corrected B (f 0.9) + D | 0.18, off from 2017 | 19.5 | 22.1 | 20.8 | 22.2 | 16.3 | −52.7 / −50.4 |
| corrected A + D | 0.2, off from 2017 | 23.2 | 24.5 | 23.8 | 27.2 | 18.2 | −51.7 / −50.4 |

7. **Turning the groin off at the 2017 fill fixes most of the test period.** With D, the corrected groin beats no groin (18.5 m against 27.6 m joint) and the drift groin (29.0 m), and ends within 3 m of CoastSat's gap. Its 2017-2024 error falls from 51.5 to 15.0 m (no groin: 9.8 m), and 2009-2016 is unchanged.
8. **B and C are the same run on the test period.** Both have weakened to b·f 0.16 by 2009, so their test runs are identical. The schedules differ only in 1996-2003, where the calibration cannot separate them (finding 6).
9. **2009-2016 favours a weaker groin than the calibration chose.** 18.4 m at strength 0.16, against the drift groin's 11.9 m at 0.04 and the pin's 7.5 m at 0.36 (schedule-refit). Choosing the strength from the test period would be fitting to it, so this is not done here.

**Caveat on D.** The 2017 switch-off was suggested by the test-period misfit, not by evidence about the structure.

**The record, checked 2026-10-09** (`../../../1-observations/structure_history.md`, "After 1995"):
- *Supports a switch-off around 2015-2017, for a different reason.* The groins fell apart rather than being buried. Holes and structural losses were reported from 2015-2016, large portions were gone by 2021, and the northern and middle groins are now submerged most days.
- *Placement.* The 2017-18 fill (Jun 2017 to Feb 2018) was placed across the whole groin field, from the Haulover to the southernmost groin.
- *No source says the fill buried or overtopped the groins.* Burial stays unconfirmed; structural loss is documented, but only from press accounts.
- *Not a clean step.* The southernmost groin stayed partly intact until its 2026 rebuild, and the county manager linked faster erosion after the 2022 fill to its decay.

What the record supports is a groin that lost most of its function over 2015-2021, which an instant switch-off at 2017 approximates. A ramp from 2015 to 2021 would match the record more closely; it has not been run.

## D as the record describes it: a decline 2015-2021 (2026-10-09)

The record (above) describes structural loss over 2015-2021, not a switch-off at the fill. Two test runs at b 0.2 f 0.8 hold the post-2004 strength to 2015, then decline linearly to a remnant by 2021, held after (`end_factor`, off spec `2015-2021rR`). The remnant is 0, or 1/3 for the southernmost of three groins still more than half intact. The 1/3 is a stand-in, not a measured value.

| end of the groin's working life | annual | photos | joint | 2009-2016 | 2017-2024 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|
| none (the 2004 step only) | 38.7 | 43.3 | 41.1 | 18.4 | 51.5 | −2.8 / −50.4 |
| off at the 2017 fill | 16.8 | 20.1 | **18.5** | 18.4 | **15.0** | −53.3 / −50.4 |
| decline 2015-2021 to 0 | 22.2 | 28.9 | 25.7 | 18.5 | 25.4 | −50.7 / −50.4 |
| decline 2015-2021 to 1/3 | 26.0 | 32.1 | 29.2 | 18.4 | 31.8 | −35.5 / −50.4 |
| no groin | 23.5 | 31.1 | 27.6 | 31.7 | 9.8 | −57.8 / −50.4 |

10. **The decline the record describes fits worse than switching off at the fill.** It ends at the right gap (−50.7 m against −50.4 m), but its 2017-2024 error is 25.4 m, against 15.0 m for the switch-off. While it is still declining, in 2017-2020, the groin holds the fill's sand at GIS 6. A remnant of 1/3 holds too much (29.2 m joint).
11. **So the 2017-2024 misfit is mostly the fill, not the groin's decay.** The data want the fill's sand to pass the groin field at once. That fits a fill placed across the whole field, filling the compartments so the groins bypassed everything, more than a gradual structural loss. No source confirms the bypass (see the record above); it is what the shoreline data imply.
