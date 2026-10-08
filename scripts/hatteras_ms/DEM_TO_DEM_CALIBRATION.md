# DEM-to-DEM calibration and test: where it stands

Started 2026-10-05 from the advisor's plan. Last updated **2026-10-06**. Start here to pick the work back up.

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
| Model years | 13 (1996–2008, ends 1 Jan 2009) and 16 (2009–2024, ends 1 Jan 2025). Runs start 1 Jan of the DEM year; the ~7.5-month offset from the DEM date is reported, not corrected |
| Windows | calibration 1995-10-12→1997-10-12 to 2008-08-17→2010-08-17; test 2008-08-17→2010-08-17 to 2025-02-17→2026-02-17 (`hat_observed_rates.NET_CHANGE_WINDOWS`) |
| Start shoreline | shoreline offset v2 (DEM-centred CoastSat mean) for both periods |
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

### Current scores (interior GIS 2–89, both sides LOWESS-7)

| run | bias | RMSE | r |
|---|---|---|---|
| 1996–2009 edgeBE, no groin | +6.0 m | 19.0 m | 0.51 |
| 1996–2009 edgeBE + groin | +5.5 | 18.1 | 0.67 |
| 1996–2009 domainBE + groin | ≈0 | 0.5 (raw model vs target) | — |
| 2009–2025 edgeBE (own ends), no groin | −10.4 | 25.0 | 0.39 |
| 2009–2025 domainBE + groin | −18.2 | 32.1 | 0.36 |

## Where we left off

**Open decision: what BE set 2 is built from.** Worth raising with the advisor.
- (a) The full test residual, as the plan states. This bakes the one-time 2021 step (about +17 m) into the forward rates as roughly +1 m/yr everywhere.
- (b) **The residual with the 2021 step removed (recommended).** Set 2 corrects the alongshore pattern only, and the step is reported as an event.
- (c) Hold set 2 until Aidan's wave change is in.

**Waiting on:** Aidan's wave-diffusivity change between the periods (plan step 3, last layer). When it arrives, rerun the calibration layers on it: the end solve, the groin check and the per-domain solve.

**Nourishment layer:** nothing to add for calibration (no fill falls in 1996–2008).

## Things to say in the write-up

- The groin is fitted on three photo dates (1997, 2004, 2008); f rests on 2008 alone. The 2004 failure timing came from the full photo record.
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

# Targets
python scripts/input_prep/5-scr/3-rates/coastsat/net_change/coastsat_net_change.py
```

Runs take about 2–5 minutes each. Logs go to `output/logs/driver/dem_to_dem/`.
