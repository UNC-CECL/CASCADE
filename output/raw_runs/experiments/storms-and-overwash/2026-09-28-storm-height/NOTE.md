# Are the storm water levels too low for Hatteras? (2026-09-28)

**Question.** With dunes that match the 2009 lidar (`../2026-09-28-dune-ceiling-per-domain/`), Irene still overwashes only 8–10 of the 60 domains with dunes of 4.3 m and up, where the imagery shows 47. The storm series takes water level from Duck (8651370), about 80 km north, adds Stockdon R2% run-up on WIS ST63228 waves, and uses one beach slope, 0.06. Is that too low for Hatteras?

Driver: `scripts/hatteras_ms/experiments/HAT_storm_height_test.py` (`gauges`, `build`, `run`, `score`). The slope check was a one-off script, recorded below.

## Part 1: observations

**Gauges** (peak water level, m above MHW, NOAA CO-OPS; `tables/gauge_peaks_m_above_mhw.csv`, raw series in `data/`):

- **Ocean side: Duck is representative.** The Cape Hatteras Fishing Pier (8654400) is the only ocean-side gauge in the reach. Against Duck it gives a median of +0.05 m over 5 storms (Fran, Bonnie, Dennis, Floyd and Isabel; range −0.25 to +0.14).
- **Irene was a sound-side event.**

  | gauge | Irene peak, m above MHW |
  |---|---|
  | Duck | 0.56 |
  | Oregon Inlet Marina (sound side) | **2.02** |
  | USCG Station Hatteras (sound side) | **1.09** |

  Matthew 2016 (USCG Hatteras 1.78 against Duck 0.61), Floyd 1999 and Arthur 2014 look similar. Barrier3D models ocean-side overwash only.

**Foreshore slope.** Measured on the 2009–2014 1 m lidar (`0-elevation/2009-2014/1-gapfill-1m`), along each east–west row from the ocean-side MHW contour to the first cell at the berm (1.7 m NAVD88), in 81 east-facing domains (GIS 10–90): median 0.067, P10–P90 0.050–0.101. The builder's 0.06 is close.

## Part 2: sensitivity (12 runs, all ran to the end)

- **Storm variants of trim24:**
  - `slope0p10`: run-up recomputed with slope 0.10 and the events re-found. That gives about twice as many events, with peaks up to 6.9 m in 1996–2010.
  - `plus0p25` and `plus0p50`: Rhigh and Rlow raised 0.25 and 0.5 m, on the same events.
- **Dunes:** both realistic ceilings, uniform 5.5 m NAVD88 and per-cell.
- **Runs:** managed, both windows, plus the trim24 controls. RMSE is against LOESS-7 for every run (see the per-domain NOTE).

| window | ceiling | storms | POD | POFD | PSS | timing r | RMSE | Irene low third (obs 25) | Irene rest (obs 47) |
|---|---|---|---|---|---|---|---|---|---|
| 1996–2010 | uniform 5.5 | trim24 | 0.67 | 0.06 | 0.60 | 0.97 | 1.20 | | |
| | | slope 0.10 | 0.75 | 0.12 | 0.64 | 0.96 | 1.25 | | |
| | | +0.25 / +0.50 m | 0.70 / 0.71 | 0.08 / 0.09 | 0.61 / 0.62 | 0.98 | 1.20 / 1.19 | | |
| | per-cell | trim24 | 0.72 | 0.15 | 0.57 | 0.85 | 1.17 | | |
| | | slope 0.10 | 0.87 | 0.42 | 0.45 | 0.60 | 1.43 | | |
| | | +0.25 / +0.50 m | 0.77 / 0.81 | 0.21 / 0.27 | 0.56 / 0.55 | 0.82 / 0.78 | 1.17 / 1.19 | | |
| 2010–2024 | uniform 5.5 | trim24 | 0.23 | 0.09 | 0.15 | 0.78 | 1.91 | 15 | 8 |
| | | slope 0.10 | 0.83 | 0.56 | 0.27 | 0.51 | 2.52 | 30 | 46 |
| | | +0.25 / +0.50 m | 0.29 / 0.32 | 0.12 / 0.16 | 0.17 / 0.16 | 0.77 / 0.74 | 1.92 / 1.93 | 19 / 23 | 10 / 11 |
| | per-cell | trim24 | 0.49 | 0.34 | 0.14 | 0.35 | 2.09 | 26 | 10 |
| | | slope 0.10 | 0.85 | 0.69 | 0.16 | 0.36 | 2.87 | 30 | 42 |
| | | +0.25 / +0.50 m | 0.57 / 0.63 | 0.42 / 0.51 | 0.15 / 0.12 | 0.36 / 0.31 | 2.05 / 2.11 | 28 / 28 | 16 / 21 |

## Answer

**The ocean-side storm heights are not too low.** The ocean gauge inside the reach matches Duck, and the measured beach slope matches the builder's.

- **A uniform surge of 0.25–0.5 m does not recover Irene.** Irene's higher-dune domains go from 8–10 to only 10–21 of 47, and false alarms rise.
- **A steeper run-up slope does recover it (42–46 of 47), but only by making every storm overtop more.** False alarms rise to 56–69% in 2010–2024, the timing correlation falls, and RMSE worsens to 2.5–2.9 m/yr. Slope 0.10 is also the steep end (P90) of the measured beaches, not typical.
- **Why the model misses Irene.** Irene's high water on Hatteras came from Pamlico Sound (Oregon Inlet 2.02 m against Duck 0.56 m). Much of the washover and breaching in the 2011 image is sound-side surge, a process Barrier3D does not have. The 2010–2024 overwash shortfall on high-dune domains is best read as that missing process, not as a storm-series error.

**Reading.** Keep the storm series at the observed ocean levels and slope 0.06 (trim24), with a realistic dune ceiling. Report the 2010–2024 overwash limitation as sound-side surge the model does not represent.

**Status: record.** Nothing adopted.
