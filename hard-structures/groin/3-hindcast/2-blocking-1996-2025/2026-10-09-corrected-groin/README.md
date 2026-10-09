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

## Strength per period, set directly (2026-10-09)

b, f and the schedules are one way of saying how strong the groin is in each stretch of its life. `periods` sets those strengths directly, as the share of the arriving drift the groin holds back from each year on:

    python corrected_groin.py periods --period 2009 --spec "1969:0.2,2004:0.16,2017:0"

That spec reads: 0.2 from installation, 0.16 from the 2004 failure, 0 from the 2017 fill. b, f, the schedule and the switch-off are not used. Runs are filed under `periods_<spec>`, and the run's `corrected_groin.txt` records the spec. **Check:** that spec on 2009-2025 reproduces the b 0.2 / f 0.8 / off-at-2017 run to 8.5 × 10⁻¹⁴ m (0.2 × 0.8 is 0.16000000000000003 in floating point), with the same scores: joint 18.5 m, 2009-2016 18.4 m, 2017-2024 15.0 m.

In this form the corrected groin's best description so far is three numbers: **0.2 to 2003, 0.16 from 2004, 0 from 2017.** The calibration barely separates the first two (findings 4 and 6). The third rests on the test period and the 2017 placement record (findings 10-11).

## The two-number groin: one strength, one end year (2026-10-09)

Because the calibration cannot separate the strength before and after 2004, the groin is described by two numbers: a constant strength s from installation, and 0 from an end year. The end year is **not fitted**: 1996-2009 never reaches it, and fitting it on 2009-2025 would be fitting to the test. It is fixed at **2017** from the placement record (the fill was placed across the groin field) and the documented decay from 2015 (`../../../1-observations/structure_history.md`). The test's 2017-2024 half still informed that choice (findings 7, 10-11), so only 2009-2016 is fully independent.

**Calibration, 1996-2009.** The schedule-A runs (constant strength) plus `periods` runs at 0.16-0.19:

| s | 0.10 | 0.12 | 0.15 | 0.16 | 0.17 | 0.18 | **0.19** | 0.20 | 0.25 | 0.30 |
|---|---|---|---|---|---|---|---|---|---|---|
| annual | 14.7 | **11.3** | 11.6 | 13.2 | 14.6 | 16.6 | 18.9 | 21.6 | 34.7 | 47.6 |
| photos | 41.5 | 34.6 | 27.5 | 25.4 | 23.5 | 21.5 | 19.4 | 17.6 | **14.5** | 21.8 |
| joint | 31.1 | 25.8 | 21.1 | 20.2 | 19.6 | 19.2 | **19.2** | 19.7 | 26.6 | 37.0 |

The joint score is flat at 0.18-0.19; at 0.19 the two scores are nearly equal. The annual record alone wants 0.12, the photos alone 0.25.

**Test, 2009-2025, end year 2017:**

| s (chosen by) | annual | photos | joint | 2009-2016 | 2017-2024 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|
| no groin | 23.5 | 31.1 | 27.6 | 31.7 | 9.8 | −57.8 / −50.4 |
| **0.12 (annual CoastSat)** | **11.7** | **18.6** | **15.5** | **11.5** | 12.0 | −54.5 / −50.4 |
| 0.19 (joint) | 21.3 | 23.3 | 22.3 | 24.6 | 17.4 | −52.2 / −50.4 |
| 0.25 (photos) | 32.0 | 31.4 | 31.7 | 39.5 | 22.0 | −51.4 / −50.4 |
| three numbers: 0.2, 0.16 from 2004, 0 from 2017 | 16.8 | 20.1 | 18.5 | 18.4 | 15.0 | −53.3 / −50.4 |

12. **The strength fitted to the annual CoastSat record predicts the test period best**, on both scores and on the independent 2009-2016 half (11.5 m against 24.6 m for the joint choice and 31.7 m for no groin). The photo-fitted strength does worse than no groin. The photos' pull toward a stronger groin rests on the 2004 survey, the one date where the photos and CoastSat disagree by 42 m (`../../../1-observations/gap_across_groins/`, figure 7).
13. **This bears on the calibration target (decision 2):** the annual CoastSat gap generalises; the three photo dates do not.
14. **The two-number groin with s 0.12 is the best-predicting groin so far**, ahead of the three-number one (18.5 m). It uses two numbers, s = 0.12 (holding back about 11 m of shoreline a year at Buxton) and an end year of 2017. Neither was fitted to 2009-2016.

## Figures (2026-10-09)

`python results_figures.py` draws four figures into `figures/`, with captions in `figures/supporting/CAPTIONS.md`:

| figure | question |
|---|---|
| `corrected_groin_1_strength_over_time` | How strong is the groin each year in each option, against the structure's history? |
| `corrected_groin_2_gap_over_time` | Does the corrected groin reproduce the observed gap, in calibration and test? (no groin, the pinned groin, the two-number groin at 0.12 and 0.19) |
| `corrected_groin_3_choosing_strength` | Which strength each observed record picks, and how each does on the test |
| `corrected_groin_4_test_summary` | Every groin version on 2009-2025, before and after the 2017 fill |
