# 2026-10-09-corrected-groin

**Question.** Does the corrected groin (straight-coast module 6: the drift groin plus a correction for BRIE's sand leak at the groin's bends) fit the 1996-2009 gap across the groins in the full model, and how does it do on 2009-2025?

**Strength superseded, 2026-10-13: see `../2026-10-13-pregroin-planform/`.** Run as one continuous history from the 1967 shoreline, the groin needs about 0.25 to 2009 and 0.20 to 2016, not the 0.16 and 0.10 fitted below. Each window here starts from the observed shoreline, which already holds the fillet, so its fitted strength is the rate that fillet relaxes at, not the groin's strength.

**Refit on the calendar-year setup, 2026-10-12 night: read that section first.** The calibration and test now run on calendar-year shorelines (offset v4, storms 1997-2009 and 2010-2025, RSLR 0.002 and 0.006; commit `e9c5e663`). Only the corrected groin was refit; its results are in "The calendar-year refit" below, in `grid_scores.csv`, and in figures 1-4. Every section between here and that one reports the six-month setup and is kept as written (`grid_scores_2026-10-12_six-month.csv`). Its numbers are superseded, and rescoring those runs today would not reproduce them, because the no-groin runs and run lengths they were scored against have been replaced.

**Numbers refreshed 2026-10-12 (six-month setup, superseded).** Every corrected-groin and no-groin number below was rerun (`python corrected_groin.py refresh`, 92 runs) under the 10-12 setup (commit `6219a55b`: six-month start lines, calendar end windows, the test period now 2009-2025 with CoastSat's gap closing by 62.0 m, not 50.4 m) and with the tracked reference (below), now the default. The conclusions held; the numbers moved by a few metres. The pinned groin (b 0.6, f 0.6) and the drift groin (b 0.1, f 0.4) were rerun on the same setup (`same_setup_runs.py`, under `groin/2026-10-12-same-setup/`); only the figures marked † come from the drift groin's own 10-08 grid, which was not rerun. The 10-09 values are in `grid_scores_2026-10-09_baseline.csv` (`python corrected_groin.py --reference baseline score`).

**Why.** The drift groin (`../2026-10-08-drift-groin/`) fitted better than the pin, but leaked about 107 m of summed shoreline onto the island through BRIE's solve. Module 6's correction removes that leak on a straight coast (`../../../2-module-tests/`). As built, though, the correction also acted on Buxton's natural bend, changing the no-groin baseline by about 100 m (`../../../2-module-tests/3-real-planform/corrected_groin/`).

**Setup.** `corrected_groin.py` replaces `BlockingGroinCallback` with the corrected groin in a child process and drives the unchanged runner, as the drift-groin study does: full management, the solved edgeBE ends, relocations off, instant failure at the 2004 step (b falls to b·f). Each year:
- **The groin** moves b of BRIE's wave-climate net drift at the updrift approach face from GIS 5 to GIS 6, as `drift_groin.py` does.
- **The correction (option B, chosen 2026-10-09)** is applied to the groin's own bends only. Within 10 cells of the groin face, it is module 6's leak term for this run's BRIE shoreline minus the same term for the groin-free shoreline in that year. Since 10-12 the groin-free shoreline is this run's minus the groin's own tracked footprint, so no no-groin run is needed (section below); the 10-09 runs rebuilt it from the no-groin run instead. The zone is 10 cells, not module 6's 5, because 5 leaks on a grid-aligned coast (`../../../2-module-tests/6-long-term/`).
- **Check (`check`):** at b 0 the run reproduces the no-groin run exactly (0.000 m in every year and domain with the tracked reference; 5 × 10⁻¹³ m with the 10-09 one). So the correction leaves the baseline untouched.

Each run dir has `corrected_groin.txt`, because the runner's report and run name (`groinblock`) do not show the patch. The diagnostics CSV gains `net_drift_m3_yr`, `approach_angle_deg`, `leak_correction_m`, `max_diff_from_no_groin_m`, `footprint_max_m` and `tracked_vs_baseline_m`. Calendar-year runs are under `output/raw_runs/experiments/groin/2026-10-12-corrected-groin-calendar/` (`--setup calendar`, the default). Six-month runs are under `output/raw_runs/experiments/groin/2026-10-09-corrected-groin/` (`--setup six-month`): `tracked*/` the tracked ones, the unprefixed folders the 10-09 runs, `now_*/` the comparison runs.

    python corrected_groin.py check
    python corrected_groin.py grid [--streams 4]     # b 0.02-0.3 x f 0.2-0.8, 32 runs
    python corrected_groin.py edge [--streams 4]     # f 0.9-1.0 at b 0.15-0.3, and b 0.25; 11 runs
    python corrected_groin.py score                  # -> grid_scores.csv (corrected and no groin; joint = RMS of the two scores)
    python corrected_groin.py test --b 0.2 --f 0.8
    python corrected_groin.py refresh [--streams 4]  # every run this README reports, then the score

## Calibration, 1996-2009

All runs clean (`refresh_run.log`). Scores are as the drift-groin study: annual CoastSat RMSE and photo-date RMSE of the GIS 5|6 gap change.

| run | annual RMSE | photo RMSE | joint | island vs no groin |
|---|---|---|---|---|
| no groin | 37.6 | 64.8 | 53.0 | 0 |
| pinned b 0.6 f 0.6 (the 10-05 pin) | 32.3 | 1.9 | 22.9 | −28 m |
| drift b 0.1 f 0.4 (10-08 best) | 23.1 | 7.6 | 17.2 | +104 m |
| **corrected b 0.2 f 0.8** | 20.9 | 19.8 | **20.4** | +12 m |
| corrected b 0.2 f 0.9 / 1.0 | 21.8 / 23.3 | 18.8 / 17.4 | 20.3 / 20.6 | +12 m |

1. **It fits a little better than the pin, and conserves sand.** Its best joint score is 20.3-20.4 m, against the pin's 22.9 m. Across the whole grid the island gains only 0.3 to 12.5 m of summed shoreline against no groin, against the drift groin's +104 m at its best cell (19 to about 200 m across its grid †). Part of the drift groin's better score (17.2 m) came from that leaked sand piling up updrift.
2. **It needs about twice the drift groin's b** (0.2 against 0.1), because without the leak the same b holds a smaller gap. At b 0.05 f 0.8 the corrected groin scores 24.8 / 52.1 m, the drift groin 13.7 / 37.6 m †.
3. **The two scores still want different groins:** the annual CoastSat score is best at b 0.12 (12.7 m at f 0.8 in this grid; 11.1 m as one constant strength, below), the photo score at b 0.25 f 0.9 (14.8 m).
4. **The calibration barely sees the 2004 failure.** The best cells sit at f 0.8 to 1.0, within 0.3 m of each other; removing the failure altogether (f 1.0) costs 0.2 m.

## Test, 2009-2025

| run | annual RMSE | photo RMSE | joint | gap change at the end, model / CoastSat |
|---|---|---|---|---|
| no groin | 24.8 | 33.1 | 29.2 | −65.4 / −62.0 |
| pinned b 0.6 f 0.6 | 29.4 | 33.6 | 31.6 | −18.2 / −62.0 |
| drift b 0.1 f 0.4 | 25.7 | 30.3 | 28.1 | −24.2 / −62.0 |
| **corrected b 0.2 f 0.8** | 33.1 | 36.7 | **35.0** | **−15.9** / −62.0 |
| corrected b 0.2 f 1.0 | 42.6 | 45.5 | 44.0 | −4.6 / −62.0 |

5. **Out of sample it does worse than no groin and than the drift groin.** It holds the gap after 2009 (−15.9 m) while CoastSat's closes by 62 m. After the 2004 step it still blocks b·f = 0.16 of the drift, four times the drift groin's 0.04. The calibration chose a groin that hardly weakened in 2004, and the later record needs one that largely gave way. This is the failure-schedule conflict open in `../2026-10-08-schedule-refit/`, sharper here because the corrected groin holds more sand. The 2017 fill placement also costs every groin run 37-43 m in this period (10-05 README).

**What this means.** The corrected groin's mechanism carries over to Buxton: with the correction confined to its own bends, it traps on the right side, conserves sand and leaves the no-groin baseline unchanged. It fits 1996-2009 as well as the pin. It does not fix the test period, because what the groin does after 2004 is set by the failure schedule, not by how it traps. The next step is the schedule decision (options A-C in `../2026-10-08-schedule-refit/`), then refitting this module under the chosen schedule.

## The failure schedules, A-D (2026-10-09)

Asked the same day: refit under every schedule option of `../2026-10-08-schedule-refit/`, plus a fourth. The schedules are taken from `schedule_refit.SCHEDULES`; D is new here.
- **A, failure from 1996:** one strength all window, run as b with f 1.0 (`instant1996`, 9 runs).
- **B, the 2004 step:** the grid above (`instant2004`).
- **C, ramp:** full in 1995, falling to b·f by 2003 (`ramp1996`, b 0.1-0.4 × f 0.2-0.8, 24 runs).
- **D:** any of these, with the groin passing everything from the 2017 Buxton fill on (`OFF_YEAR`). It is the hypothesis that a 2 million m³ fill buried or overtopped the groins. It changes only the test period. Each schedule's best calibration cell is run on 2009-2025 with and without it.

The test RMSE is also split at the fill: 2009-2016 and 2017-2024.

**Calibration, best joint cell per schedule:**

| schedule | b, f | annual | photos | joint |
|---|---|---|---|---|
| B, 2004 step | 0.2, 0.9 | 21.8 | 18.8 | 20.3 |
| A, from 1996 | 0.2 (constant) | 23.3 | 17.4 | 20.6 |
| C, ramp | 0.2, 0.8 | 18.0 | 24.2 | 21.3 |

6. **The calibration cannot tell the schedules apart.** All three best cells score 20.3-21.3 m, all at b 0.2 with little or no weakening, as the pin study found.

**Test, 2009-2025:**

| run | strength after failure | annual | photos | joint | 2009-2016 | 2017-2024 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|---|
| no groin | 0 | 24.8 | 33.1 | 29.2 | 34.3 | 11.2 | −65.4 / −62.0 |
| pinned b 0.6 f 0.6 | 0.36 | 29.4 | 33.6 | 31.6 | 8.7 | 38.9 | −18.2 / −62.0 |
| drift b 0.1 f 0.4 | 0.04 | 25.7 | 30.3 | 28.1 | 14.4 | 31.9 | −24.2 / −62.0 |
| corrected B or C (b 0.2 f 0.8) | 0.16 | 33.1 | 36.7 | 35.0 | 15.3 | 42.8 | −15.9 / −62.0 |
| corrected B (f 0.9) | 0.18 | 38.2 | 41.8 | 40.0 | 18.4 | 49.3 | −10.7 / −62.0 |
| corrected A (b 0.2) | 0.2 | 42.6 | 45.5 | 44.0 | 23.0 | 54.0 | −4.6 / −62.0 |
| **corrected B or C + D** | 0.16, off from 2017 | **13.5** | **18.4** | **16.1** | 15.3 | **12.3** | **−65.0** / −62.0 |
| corrected B (f 0.9) + D | 0.18, off from 2017 | 15.6 | 19.5 | 17.7 | 18.4 | 13.4 | −65.0 / −62.0 |
| corrected A + D | 0.2, off from 2017 | 18.5 | 20.6 | 19.5 | 23.0 | 14.0 | −64.7 / −62.0 |

7. **Turning the groin off at the 2017 fill fixes most of the test period.** With D, the corrected groin beats no groin (16.1 m against 29.2 m joint), the drift groin (28.1 m) and the pin (31.6 m), and ends within 3 m of CoastSat's gap. Its 2017-2024 error falls from 42.8 to 12.3 m (no groin: 11.2 m), and 2009-2016 is unchanged.
8. **B and C are the same run on the test period.** Both have weakened to b·f 0.16 by 2009, so their test runs are identical. The schedules differ only in 1996-2003, where the calibration cannot separate them (finding 6).
9. **2009-2016 favours a weaker groin than the calibration chose.** 15.3 m at strength 0.16, against 14.4 m for the drift groin at 0.04 and 10.3 m for the two-number groin at 0.12 (below). The pin does best on this half (8.7 m), then misses 2017-2024 by 38.9 m. Choosing the strength from the test period would be fitting to it, so this is not done here.

**Caveat on D.** The 2017 switch-off was suggested by the test-period misfit, not by evidence about the structure.

**The record, checked 2026-10-09** (`../../../1-observations/data/structure_history.md`, "After 1995"):
- *Supports a switch-off around 2015-2017, for a different reason.* The groins fell apart rather than being buried. Holes and structural losses were reported from 2015-2016, large portions were gone by 2021, and the northern and middle groins are now submerged most days.
- *Placement.* The 2017-18 fill (Jun 2017 to Feb 2018) was placed across the whole groin field, from the Haulover to the southernmost groin.
- *No source says the fill buried or overtopped the groins.* Burial stays unconfirmed; structural loss is documented, but only from press accounts.
- *Not a clean step.* The southernmost groin stayed partly intact until its 2026 rebuild, and the county manager linked faster erosion after the 2022 fill to its decay.

What the record supports is a groin that lost most of its function over 2015-2021, which an instant switch-off at 2017 approximates. A ramp from 2015 to 2021 matches the record more closely; it is run in the next section.

## D as the record describes it: a decline 2015-2021 (2026-10-09)

The record (above) describes structural loss over 2015-2021, not a switch-off at the fill. Two test runs at b 0.2 f 0.8 hold the post-2004 strength to 2015, then decline linearly to a remnant by 2021, held after (`end_factor`, off spec `2015-2021rR`). The remnant is 0, or 1/3 for the southernmost of three groins still more than half intact. The 1/3 is a stand-in, not a measured value.

| end of the groin's working life | annual | photos | joint | 2009-2016 | 2017-2024 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|
| none (the 2004 step only) | 33.1 | 36.7 | 35.0 | 15.3 | 42.8 | −15.9 / −62.0 |
| off at the 2017 fill | 13.5 | 18.4 | **16.1** | 15.3 | **12.3** | −65.0 / −62.0 |
| decline 2015-2021 to 0 | 17.3 | 24.8 | 21.4 | 15.4 | 19.9 | −62.6 / −62.0 |
| decline 2015-2021 to 1/3 | 21.0 | 27.4 | 24.4 | 15.4 | 25.9 | −46.2 / −62.0 |
| no groin | 24.8 | 33.1 | 29.2 | 34.3 | 11.2 | −65.4 / −62.0 |

10. **The decline the record describes fits worse than switching off at the fill.** It ends at the right gap (−62.6 m against −62.0 m), but its 2017-2024 error is 19.9 m, against 12.3 m for the switch-off. While it is still declining, in 2017-2020, the groin holds the fill's sand at GIS 6. A remnant of 1/3 holds too much (24.4 m joint).
11. **So the 2017-2024 misfit is mostly the fill, not the groin's decay.** The data want the fill's sand to pass the groin field at once. That fits a fill placed across the whole field, filling the compartments so the groins bypassed everything, more than a gradual structural loss. No source confirms the bypass (see the record above); it is what the shoreline data imply.

## Strength per period, set directly (2026-10-09)

b, f and the schedules are one way of saying how strong the groin is in each stretch of its life. `periods` sets those strengths directly, as the share of the arriving drift the groin holds back from each year on:

    python corrected_groin.py periods --period 2009 --spec "1969:0.2,2004:0.16,2017:0"

That spec reads: 0.2 from installation, 0.16 from the 2004 failure, 0 from the 2017 fill. b, f, the schedule and the switch-off are not used. Runs are filed under `periods_<spec>` (`tracked_periods_<spec>` since 10-12), and the run's `corrected_groin.txt` records the spec. **Check:** that spec on 2009-2025 reproduces the b 0.2 / f 0.8 / off-at-2017 run to 1.1 × 10⁻¹³ m (0.2 × 0.8 is 0.16000000000000003 in floating point), with the same scores: joint 16.1 m, 2009-2016 15.3 m, 2017-2024 12.3 m.

In this form the corrected groin's best description so far is three numbers: **0.2 to 2003, 0.16 from 2004, 0 from 2017.** The calibration barely separates the first two (findings 4 and 6). The third rests on the test period and the 2017 placement record (findings 10-11).

## The two-number groin: one strength, one end year (2026-10-09)

Because the calibration cannot separate the strength before and after 2004, the groin is described by two numbers: a constant strength s from installation, and 0 from an end year. The end year is **not fitted**: 1996-2009 never reaches it, and fitting it on 2009-2025 would be fitting to the test. It is fixed at **2017** from the placement record (the fill was placed across the groin field) and the documented decay from 2015 (`../../../1-observations/data/structure_history.md`). The test's 2017-2024 half still informed that choice (findings 7, 10-11), so only 2009-2016 is fully independent.

**Calibration, 1996-2009.** The schedule-A runs (constant strength) plus `periods` runs at 0.16-0.19:

| s | 0.10 | 0.12 | 0.15 | 0.16 | 0.17 | **0.18** | 0.19 | 0.20 | 0.25 | 0.30 |
|---|---|---|---|---|---|---|---|---|---|---|
| annual | 13.0 | **11.1** | 13.2 | 15.2 | 16.7 | 19.3 | 21.5 | 23.3 | 35.2 | 49.1 |
| photos | 37.5 | 33.0 | 26.7 | 23.4 | 22.2 | 19.8 | 18.7 | 17.4 | **15.5** | 23.2 |
| joint | 28.0 | 24.6 | 21.0 | 19.7 | 19.6 | **19.6** | 20.1 | 20.6 | 27.2 | 38.4 |

The joint score is flat at 0.16-0.18, lowest at 0.18, where the two scores are nearly equal (the 10-09 runs put it at 0.18-0.19). The annual record alone wants 0.12, the photos alone 0.25.

**Test, 2009-2025, end year 2017:**

| s (chosen by) | annual | photos | joint | 2009-2016 | 2017-2024 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|
| no groin | 24.8 | 33.1 | 29.2 | 34.3 | 11.2 | −65.4 / −62.0 |
| **0.12 (annual CoastSat)** | **9.7** | **18.1** | **14.6** | **10.3** | **9.7** | −63.5 / −62.0 |
| 0.18 (joint) | 15.6 | 19.5 | 17.7 | 18.4 | 13.4 | −65.0 / −62.0 |
| 0.19 (the 10-09 joint choice) | 17.1 | 20.0 | 18.6 | 20.8 | 13.7 | −64.6 / −62.0 |
| 0.25 (photos) | 26.3 | 26.0 | 26.1 | 34.3 | 17.2 | −62.3 / −62.0 |
| three numbers: 0.2, 0.16 from 2004, 0 from 2017 | 13.5 | 18.4 | 16.1 | 15.3 | 12.3 | −65.0 / −62.0 |

12. **The strength fitted to the annual CoastSat record predicts the test period best**, on both scores, on both halves, and on the independent 2009-2016 half (10.3 m against 18.4 m for the joint choice and 34.3 m for no groin). The photo-fitted strength is the worst of the groins: it beats no groin overall (26.1 against 29.2 m), but not on 2009-2016 (34.3 m each). The photos' pull toward a stronger groin rests on the 2004 survey, the one date where the photos and CoastSat disagree by 42 m (`../../../1-observations/analysis/2-coastsat_gap/`, figure 7).
13. **This bears on the calibration target (decision 2):** the annual CoastSat gap generalises; the three photo dates do not.
14. **The two-number groin with s 0.12 is the best-predicting groin so far**, ahead of the three-number one (16.1 m). It uses two numbers, s = 0.12 (holding back about 11.5 m of shoreline a year at Buxton) and an end year of 2017. Neither was fitted to 2009-2016.

## No no-groin run needed: the tracked reference (2026-10-12)

**Why.** Option B needs the no-groin run's shoreline every year, so every groin run depended on a matching
no-groin run. Two ways to drop that were tried.

1. **Freeze the natural bend at the start** (correction = L(x) − L(x₀)). Rejected without a run: from the
   no-groin runs, the natural bend's leak changes so much that this would put up to 175-220 m of false
   correction into the cells beside the groin, more than the groin's own effect (about 80 m by 2009).
2. **Track the groin's own footprint** (`tracked`, adopted). Each year the callback repeats BRIE's alongshore
   step for last year (`brie_step`, a copy of `brie.py:1293-1326`) twice, with and without its own
   contribution (the drift it moved and the correction it added), and keeps the difference as the groin's
   footprint a. The groin-free shoreline is then x − a. Fills, Barrier3D's shoreline change and anything else
   enter both steps alike and cancel. With no groin, a is exactly zero, so the no-groin run is reproduced by
   construction. What it leaves out: Barrier3D responding differently to the groin's sand (overwash,
   shoreface), since a is spread only by BRIE's step. It adds 15 lines to the callback.

**Test (`python corrected_groin.py tracked`, `tracked_vs_baseline.csv`).** Both references on the same
setup (10-12; the option-B runs rerun as `now_*`), the two-number groin, off from 2017:

| period | strength | reference | annual | photos | joint | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|
| 1996-2009 | 0.12 | no-groin run | 10.7 | 30.3 | 22.7 | −38.4 / −28.1 |
| | | tracked | 11.1 | 33.0 | 24.6 | −42.6 / −28.1 |
| | 0.19 | no-groin run | 23.5 | 15.6 | 19.9 | −15.3 / −28.1 |
| | | tracked | 21.5 | 18.7 | 20.1 | −19.7 / −28.1 |
| 2009-2025 | 0.12 | no-groin run | 10.4 | 18.0 | 14.7 | −62.9 / −62.0 |
| | | tracked | 9.7 | 18.1 | 14.6 | −63.5 / −62.0 |
| | 0.19 | no-groin run | 18.6 | 21.3 | 20.0 | −61.9 / −62.0 |
| | | tracked | 17.1 | 20.0 | 18.6 | −64.6 / −62.0 |

15. **The tracked reference reproduces no groin exactly** (0.000 m at b 0, every year and domain) and stays
    within 4-5 m of option B's shoreline, 4.7 m on the gap and 6 m on the island total. Its groin-free
    shoreline is within 5-8 m of the no-groin run's near the groin.
16. **It changes no conclusion.** The annual record still picks about 0.12, and 0.12 still predicts the test
    period best (joint 14.6, option B 14.7). The one consistent difference: in calibration the tracked groin
    holds about 4 m less gap, the Barrier3D response it leaves out.
17. **Adopted as the default** (`REFERENCE = "tracked"`); `--reference baseline` reaches the 10-09 runs and
    their scores (`grid_scores_2026-10-09_baseline.csv`). Tracked runs are filed under `tracked*/`.

## The calendar-year refit (2026-10-12 night)

**Setup.** The calibration and test were moved to calendar-year shorelines (start and end lines are calendar-year means; offset v4; storms 1997-2009 and 2010-2025; RSLR 0.002 and 0.006; the solved edgeBE ends for each). Row k of a run is calendar year start + k and the start line is row 0, so both periods start from a gap change of exactly zero. The runner passes the first model year to the callback, so the 2004 failure and the 2017 switch-off land in the right calendar years. `corrected_groin.py` scores this with its own `observed_change`, `model_change`, `score_annual`, `score_photos` and `split_scores`, which fall back to `schedule_refit.py`'s under `--setup six-month`. The test splits are now **2010-2016** (fully independent) and **2017-2025**. CoastSat's gap changes by −18.6 m over 1996-2009 and −60.2 m over 2009-2025; the photo dates are 1997, 2004 and 2008 in calibration, 2014-2023 in the test.

**Runs.** `python corrected_groin.py check` (b 0 reproduces no groin to 0.000 m, as before), then `refresh` (98 runs plus 2 re-tests, all clean), then by hand with `periods`: calibration at s 0.20, 0.22 and 0.25, test at s 0.16 and 0.20; then `score`. Logs are in `logs/calendar/`. The pinned and drift groins were not rerun, so they have no calendar-year numbers.

**Calibration, 1996-2009, the two-number groin (constant strength s; 0 from 2017, which 1996-2009 never reaches):**

| s | 0.10 | 0.12 | 0.14 | 0.15 | **0.16** | 0.17 | 0.18 | **0.20** | 0.22 | **0.25** | 0.30 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| annual | 18.9 | 15.3 | 12.9 | 13.0 | **12.6** | 13.9 | 14.7 | 19.1 | 23.5 | 30.3 | 42.7 |
| photos | 41.2 | 36.0 | 31.4 | 29.3 | 26.1 | 24.5 | 23.6 | 18.3 | 16.4 | **15.1** | 20.6 |
| joint | 32.0 | 27.6 | 24.0 | 22.6 | 20.5 | 19.9 | 19.7 | **18.7** | 20.3 | 23.9 | 33.5 |

No groin: 46.9 / 69.8 / 59.4 m (the gap closes by 91 m against 19 m observed). Across the b × f grid the best joint cell of each schedule scores 18.7-21.4 m, all at b 0.2 with f 0.8-1.0: the calibration still cannot see the 2004 failure.

**Test, 2009-2025:**

| run (s chosen by) | annual | photos | joint | 2010-2016 | 2017-2025 | gap change at end, model / CoastSat |
|---|---|---|---|---|---|---|
| no groin | 21.8 | 23.5 | 22.7 | 32.4 | **5.2** | −61.0 / −60.2 |
| 0.12, off from 2017 (weaker than any target) | **10.1** | **9.2** | **9.7** | **12.0** | 8.2 | −60.4 / −60.2 |
| **0.16, off from 2017 (annual CoastSat)** | 14.4 | 11.2 | 12.9 | 18.5 | 10.3 | −59.6 / −60.2 |
| 0.18, off from 2017 | 17.6 | 13.1 | 15.5 | 23.3 | 11.5 | −58.5 / −60.2 |
| 0.20, off from 2017 (both records) | 20.9 | 15.7 | 18.5 | 28.4 | 12.3 | −58.4 / −60.2 |
| 0.25, off from 2017 (photos) | 28.3 | 21.3 | 25.0 | 39.6 | 14.4 | −58.0 / −60.2 |
| three numbers: 0.2, 0.16 from 2004, off from 2017 | 14.4 | 11.2 | 12.9 | 18.5 | 10.3 | −59.6 / −60.2 |
| 0.2, 0.16 from 2004, never off | 37.3 | 38.4 | 37.9 | 18.5 | 47.0 | −12.2 / −60.2 |
| 0.2, 0.16 from 2004, declining 2015-2021 to 0 | 19.6 | 20.0 | 19.8 | 18.6 | 20.4 | −58.4 / −60.2 |

18. **Every target now picks a stronger groin than on the six-month setup:** annual CoastSat 0.16 (was 0.12), both records 0.20 (was 0.18), the photos still 0.25. The annual minimum is flat from 0.14 to 0.16 (12.6-13.0 m).
19. **The ranking of the targets holds.** On the test the annual-fitted strength beats the joint and photo ones on every score, and beats no groin by about 10 m (joint 12.9 against 22.7; 2010-2016 18.5 against 32.4). The photo-fitted strength does worse than no groin on 2010-2016 (39.6 m). The case for annual CoastSat as the target (findings 12-13) stands.
20. **The test prefers a weaker groin than any target picks.** The test error rises steadily with strength, and 0.12 is the best run on every score (joint 9.7, 2010-2016 12.0 m). At 0.16 the model widens the gap by about 12 m in 2010-2011 while CoastSat's closes by 15 m. Choosing the strength from the test would be fitting to it, so 0.12 is reported, not adopted.
21. **The groin still has to stop at the 2017 fill.** Left on, it misses 2017-2025 by 47 m; declining over 2015-2021, by 20 m; off from 2017, by 8-14 m. No groin is still best on 2017-2025 (5.2 m): every groin run lags the observed drop in 2017 by about 20 m (figure 2b), holding the sand it trapped before the fill for a year.
22. **The three-number groin and the 0.16 two-number groin are the same run on the test**, since both hold 0.16 from 2009 to 2016. With the 2004 failure invisible in calibration, the two-number groin at 0.16, off from 2017, is the one to carry forward if the target is annual CoastSat.

## Figures (2026-10-09, figures 1-4 redrawn 2026-10-12 night on the calendar-year setup)

`python results_figures.py` draws figures 1-4 into `figures/`, with captions in `figures/supporting/CAPTIONS.md`, from the calendar-year runs. The pinned and drift groins have no calendar-year runs and are no longer drawn. Figure 5 and the GIF have not been redrawn and still show the six-month runs (`same_setup_runs.py` made their pinned-groin runs):

| figure | question |
|---|---|
| `corrected_groin_1_strength_over_time` | How strong is the groin each year in each option, against the structure's history? |
| `corrected_groin_2_gap_over_time` | Does the corrected groin reproduce the observed gap, in calibration and test? (no groin, the two-number groin at 0.16 and 0.20) |
| `corrected_groin_3_choosing_strength` | Which strength each observed record picks, and how each (and 0.12) does on the test |
| `corrected_groin_4_test_summary` | Every corrected-groin version on 2009-2025, before and after the 2017 fill |
| `corrected_groin_5_plan_view` | The shorelines on a map (1.5 km windows, true scale) in 2008, 2016 and 2024: observed, no groin, pinned and corrected (`plan_view_figure.py`) |
| `corrected_groin_on_off_<calibration|test>.gif` | The groin on against off, year by year, on the 1.5 km map with the alongshore change beside it (`groin_on_off_gif.py`; git-ignored) |
| `corrected_groin_panels_<calibration|test>.gif` | Stacked plan views along GIS 2-12, land up: no groin, fills, the groin; accretion yellow, erosion red, cross-shore exaggerated about 3 times (`groin_panels_gif.py`, whose `render` the 10-13 study reuses; git-ignored) |
| `corrected_groin_effect.gif` | The shoreline with and without the corrected groin, year by year on the 1.5 km map, with the gap alongside (`groin_effect_gif.py`; git-ignored, rerun to make it) |
