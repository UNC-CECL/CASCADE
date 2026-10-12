# DEM-to-DEM calibration and test: where it stands

Started 2026-10-05 from the advisor's plan. Last updated **2026-10-12**. Start here to pick the work back up.

## The plan (advisor, 2026-10-05)

1. **Periods.** Calibration from the 1996 DEM to the 2009 DEM; test from the 2009 DEM to 2025. The calibration end shoreline and the test start shoreline are the same.
2. **Shorelines.** Each start is the CoastSat mean over ±1 yr of its DEM date. Each target is the mean over ±1 yr of its end date, except 2025, which uses ±6 months because CoastSat stops on 2026-01-13.
3. **Calibrate in layers.** Base, then end effects, then nourishment, then the wave-diffusivity change (Aidan).
4. **Sources and sinks.** Tune them until the modelled start-to-end change matches the observed change. Do not tune to the long-term rate, which also carries storms, fills and the groin.
5. **Test.** Run the test period with its fills and the calibration sources and sinks.
6. **Second set.** The test residual becomes a second set of sources and sinks for forward scenarios.

## Settled choices

| | |
|---|---|
| Period keys | `HATTERAS_PERIODS` 1996 (end 2009) and 2009 (end 2025); the old 2010 key was renamed 2009 |
| Model years | 13 (1996–2008, ends 1 Jan 2009) and **17 (2009–2025, ends 1 Jan 2026; was 16 until 2026-10-09)**. Runs start 1 Jan of the DEM year; the ~7.5-month offset from the DEM date is reported, not corrected |
| Windows | the same rules in both periods since 2026-10-12. **Start** = the island offset's window, the DEM midpoint ±6 months: 1996-04-12→1997-04-12 and 2009-02-17→2010-02-17. **End** = the last modelled year ±6 months (1 Jul→30 Jun): 2007-07-01→2009-06-30 and 2024-07-01→2026-06-30 (CoastSat ends 2026-01-13). No end is tied to a DEM date (`hat_observed_rates.NET_CHANGE_WINDOWS`). Study: `experiments/target-end-window/2026-10-10-end-options-calendar/` |
| Start shoreline | shoreline offset **v3** (DEM midpoint ±6 months, min 5 positions) for both periods, CURRENT since 2026-10-12 (v2, ±1 yr, before) |
| Storms | all storms kept for now. Open with the advisor: drop the storms before each DEM flight (2009: 3, Feb–Apr; 1996: 9 of 11, incl. the 31 Aug and 5 Sep events) or keep all |
| Target | net change in metres, 7-domain LOWESS (GIS 1–10 raw), scored over the interior GIS 2–89 |
| Groin | included; **blocking, b 0.6, f 0.6**, pinned from the calibration fit |
| Relocations | off ("let the model behave naturally") |
| 2017 Buxton fill | footprint kept as reported. The observed sand moved south past the groin, which the model cannot do; reported, not corrected (`data/hatteras_init/4-mgmt-forcing/README.md`) |
| Runs | new matrix folders `output/raw_runs/matrix/1996_2009/` and `2009_2025/` |

## Done

| step | result | where | commit |
|---|---|---|---|
| Periods, storms, RSLR, CoastSat LRR | storms v3_split12_trim24 per window; RSLR 0.003 (1996–2009, fit 0.0025 ± 0.0025) and 0.005 | `hatteras_site_config.py`, `3-env-forcings/` | 4832ae6b |
| Net-change targets | island mean −9.2 m (12.85 yr) and +14.1 m (16 yr) | `5-scr/3-rates/coastsat/net_change/` | 0ae1316e |
| End rates, calibration | GIS 1 **+1.4981**, GIS 90 **+10.5659** m/yr | `experiments/end-domain-boundaries/2026-10-05-ends-solved-on-net-change-1996_2009/` | 12411289 |
| End rates, test (edge-only base run) | GIS 1 **+32.7049**, GIS 90 **+21.0679** m/yr (GIS 1 noisy, ±1–2 m) | `…-ends-solved-on-net-change-2009_2025/` | 1c43b4d5 |
| Blocking groin fit + pin | calibration date RMSE 4.0 m (no groin 69.2); test fits until the 2017 fill | `hard-structures/groin/3-hindcast/2-blocking-1996-2025/2026-10-05-blocking-fit-calibration/` | bf498225, 89a54887 |
| Dipole vs blocking comparison | blocking keeps the real step at the groin; its two flanks are unbalanced (inherits BRIE's non-conservation) | same folder, `figures/` | e55b2bf1 |
| Per-domain source/sink, calibration (**BE set 1**) | interior RMSE 0.49 m after 10 passes; GIS 81–85 alternate on purpose; smoothing the field raises RMSE to 4.3 m | `experiments/source-sink/2026-10-05-be-domain-solve-1996_2009/`; preset **`domainBE`** | ca2940f4, 7fd83b5c, 232a5fa1 |
| Test on BE set 1 | interior bias −18.2 m, RMSE 32.1 m, r 0.36 (worse bias than edge only) | `matrix/2009_2025/domainBE/` | a14e37dd |
| 2021 CoastSat step check | the test bias **is** the step: −18.2 → −1.1 m with the step removed; RMSE unchanged, so the alongshore pattern does not carry across periods | `experiments/source-sink/2026-10-06-test-target-and-the-2021-step/` | db7cd3b9 |
| Net-change run figure | every run on these windows draws `figures/shoreline_position_change_with_buffers.png` (both sides LOWESS-7, ±100 m axis) | `cascade_pipeline/plotting/net_change_comparison.py` | dbfeb437 |
| Groin condition from CoastSat | the gap across the groin stopped widening in **1995** (90% 1992–98), at the last repair, then declined slowly; no separate drop at Isabel. The photos agree except the single 2004 survey (+42 m off), which set the 2004 failure step. Last repair corrected 1996 → 1995 (CSE 2013) | `hard-structures/groin/1-observations/gap_across_groins/` | 69ff268c, dc16bab4 |
| End rates re-solved (2026-10-09) | calibration GIS 1 **+1.4041**, GIS 90 **+10.3439** (calendar 2008 ±6-month end); test **+32.2422** / **+21.5341** (17-yr run). Both moved < 0.5 m/yr; interior scores unchanged | `experiments/end-domain-boundaries/2026-10-09-ends-solved-on-net-change-*/` | not yet |
| Clean re-run (2026-10-12) | v3 offsets + new windows; every old matrix run deleted; targets rebuilt (island mean −4.7 m calibration, +12.6 m test, 905 transects); ends re-solved: calibration **+1.6765 / +10.8869** (5 steps), test **+29.1320 / +21.7947** (8 steps), full window **+139.19 / +68.1128** (on the LRR, 3 steps); zeroBE + edgeBE re-run for all three windows (edgeBE bit-identical to the last solve step). **domainBE not re-run** | `experiments/end-domain-boundaries/2026-10-12-*` | not yet |
| Groin schedule refit | 3 failure schedules × the b-f grid on the annual CoastSat gap: all reach ~10 m (no groin 41.6); best is **failure from 1996, b × f 0.36** (9.8 m), i.e. a weak, constant groin. The pin scores 29.2 m on CoastSat but 4.0 m on the photos. Test runs depend only on b × f, so the pin and the CoastSat best are identical from 2009 | `hard-structures/groin/3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/` | not yet |

### Current scores (interior GIS 2–89, both sides LOWESS-7)

Since 2026-10-12 (v3 offsets, the 10-12 windows, all storms):

| run | bias | RMSE | r |
|---|---|---|---|
| 1996–2009 zeroBE, no groin | +0.6 m | 16.7 m | 0.55 |
| 1996–2009 edgeBE, no groin | +1.3 | 16.9 | 0.54 |
| 2009–2025 zeroBE, no groin | −12.1 | 27.2 | 0.33 |
| 2009–2025 edgeBE (own ends), no groin | −10.4 | 24.4 | 0.46 |
| 1996–2025 zeroBE (net change to the test end window) | −12.9 | 36.8 | 0.38 |
| 1996–2025 edgeBE | −9.1 | 33.1 | 0.50 |

domainBE is not re-run (BE set 1 and the groin pin were fitted on the old targets).

2026-10-09 to 10-12 (v2 offsets): calibration edgeBE +2.4 / 17.0 / 0.53, domainBE −3.3 / 5.6 / 0.98; test edgeBE −11.6 / 25.1 / 0.42, domainBE −20.2 / 34.0 / 0.36.

## Where we left off

### 2026-10-12: rules decided, clean re-run

**Done (uncommitted):** the rules above (v3 offset CURRENT, `SHORELINE_WINDOW_FOR_YEAR` and `NET_CHANGE_WINDOWS` updated, `coastsat_net_change.py` labels follow the window lengths); targets rebuilt; all old matrix runs deleted at Hannah's request; ends re-solved and zeroBE + edgeBE re-run for 1996_2009, 2009_2025 and 1996_2025; full-window figures redrawn.

**Next:**
1. **Advisor's storm answer** (keep all vs drop pre-DEM). Dropping them means new storm files for both periods and another full re-run and end solve.
2. **Groin decision**, then the **per-domain solve** (BE set 1) on the new calibration target, then domainBE runs in both periods.
3. **Commit** by path.


### 2026-10-09: periods, targets and end rates

**Done (all uncommitted):**
- **Periods** match Hannah's definition: calibration 1996 through end of 2008 (13 yr, unchanged); test 2009 through end of **2025** (17 yr, ending 1 Jan 2026; `last_model_year` 2025 in `hatteras_site_config.py`, was 2024).
- **Calibration target end** is calendar 2008 ±6 months, 2007-07-01 → 2009-06-30 (the advisor's pick; leaves out Nor'Ida and 2010). `NET_CHANGE_WINDOWS` updated and targets rebuilt. The calibration end is no longer the test start. Options and scores: `experiments/target-end-window/2026-10-07-end-window-options/`.
- **±6-month start lines tested** (Hannah's definition of the start): offset `shoreline/v3` for 1996 and 2009, built with a 5-position minimum (the 6-month windows are too thin at 10), **not CURRENT**. Scores move about 1 m, so v2 (±1 yr) stays the default: `experiments/target-end-window/2026-10-09-6mo-starts/`.
- **Matrix re-run** (`matrix/1996_2009`, `2009_2025`); the old runs are in `output/raw_runs/archive/2026-10-09-pre-test-2025/` (calibration bit-identical; deletable once checked).
- **End rates re-solved** on the new targets (Done table above): 1996 +1.4041 / +10.3439, 2009 +32.2422 / +21.5341; matrix edgeBE re-run on them.

**Next:**
1. **Decide the start window.** Keep ±1 yr (v2), or adopt ±6 months (v3: set `CURRENT`, `SHORELINE_WINDOW_FOR_YEAR` and the starts in `NET_CHANGE_WINDOWS`, then re-run the matrix and re-solve the ends).
2. **Groin check** on the new calibration target, then the **per-domain solve** (BE set 1, `domainBE`). Both were fitted to the old targets; the runs use them as they stand. Check the 2026-10-09 groin commits (corrected groin, strength + end year) against the open groin decision below first.
3. **Commit** the 10-09 work by path (config, targets, offsets v3 PROVENANCE, experiments, matrix runs).

**Loose ends:** `HAT_compare_offset_versions.py` fails for 2009 (v2's file names its column `2010`); `target_end_window.py` (10-07 study) assumes a 1 Jan 2025 test end, so re-running it would mis-date the 17-yr runs.

**Open decision: the groin's failure schedule** (from the 2026-10-08 refit; README in that folder).
- (A) **Failure from 1996, b 0.6 f 0.6 kept (recommended).** It fits CoastSat best and matches the repair record. The runner's delay goes 35 → 27. It changes only 1996–2003; test and forward runs are unchanged. It needs `domainBE` (BE set 1) re-solved and the end-domain check repeated.
- (B) Keep the 2004 step (b 0.6 f 0.6). It fits the photos best, but rests on the single 2004 survey that CoastSat contradicts.
- (C) Ramp 1996 → 2003, b 0.4 f 0.8. It falls between A and B and sits on the f grid edge.

**Open decision: what BE set 2 is built from.** Worth raising with the advisor.
- (a) The full test residual, as the plan states. This bakes the one-time 2021 step (about +17 m) into the forward rates as roughly +1 m/yr everywhere.
- (b) **The residual with the 2021 step removed (recommended).** Set 2 corrects the alongshore pattern only, and the step is reported as an event.
- (c) Hold set 2 until Aidan's wave change is in.

**Waiting on:** Aidan's wave-diffusivity change between the periods (plan step 3, last layer). When it arrives, rerun the calibration layers on it: the end solve, the groin check and the per-domain solve.

**Nourishment layer:** nothing to add for calibration (no fill falls in 1996–2008).

## Things to say in the write-up

- The groin is fitted on three photo dates (1997, 2004, 2008); f rests on 2008 alone. The 2004 failure timing came from the full photo record.
- CoastSat puts the end of the groin's effect at the 1995 last repair, not Isabel. On the annual CoastSat gap the failure timing cannot be identified, only the strength (about 0.3–0.4 of full blocking). The photos and CoastSat disagree only at the 2004 survey.
- The blocking groin partly stands in for the cape, which BRIE cannot hold. Its two flanks are unbalanced.
- The 1996–2009 source/sink pattern does not persist into 2009–2025 (the two targets correlate at r = 0.28 alongshore).
- The test-period cape advance (GIS 1–3, up to +160 m) is not reproduced by any version of the model.
- The Buxton fill sand at the groin end moved south past the groin; the model keeps it at GIS 6.

## How to resume

```
# Score / redraw: any run on these windows
python scripts/figure_making/model_output/rerender_run_figures.py --arm matrix/2009_2025

# End rates for a period, on its net-change target
python scripts/hatteras_ms/HAT_end_solve_net_change.py solve --period 1996

# Per-domain source/sink (calibration), and the smoothing test
python scripts/hatteras_ms/HAT_be_domain_solve_net_change.py solve --period 1996
python scripts/hatteras_ms/HAT_be_domain_solve_net_change.py smooth --period 1996 --step 10

# Groin condition (CoastSat gap, breakpoint) and the failure-schedule refit
python hard-structures/groin/1-observations/gap_across_groins/coastsat_gap_across_groins.py
python hard-structures/groin/3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/schedule_refit.py score

# Targets
python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py
```

Runs take about 2–5 minutes each. Logs go to `output/logs/driver/dem_to_dem/`.
