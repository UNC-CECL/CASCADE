# Dune growth rate, 2009-2025 (test period)

Where the runner's 2009 dune growth rate comes from: `DUNE_GROWTH_RATES[2009] = (0.02, 0.04)` in
`HAT_hindcast_1984_2024.ipynb` / `.py` (per period since 2026-10-21, Hannah). The calibration rate,
0.01-0.03, is in [`../../1996-start/dune-growth-rate/`](../../1996-start/dune-growth-rate/README.md).

**Question (Hannah, 2026-10-20):** did the dunes build faster after 2009? The island-wide median crest rose
~0.06 m/yr over 2009-2019, against ~0 over 1996-2009. Base Barrier3D growth; only rmin-rmax changes. The 2009
start, Dmaxel 12.2 m, groin on, the ends re-solved on 2026-10-20 (for 0.01-0.03). Scored against the 2014,
2016, 2017, 2018 and 2019 lidar (corrected to the 2009 lidar on NC-12, `../survey-crests/`) on the 2009 dune
windows: island-wide median and low-dune 10th-percentile crest per survey, overwash (Irene, after Sandy,
March 2024) and the interior shoreline.

| rate | median crest error, 2014 / 16 / 17 / 18 / 19 (m) | RMSE | P10 RMSE | Irene | shoreline RMSE |
|---|---|---|---|---|---|
| 0.01-0.03 (calibration rate) | -0.14 / -0.25 / -0.31 / -0.26 / -0.39 | 0.28 | 0.36 | 30/72 | 2.072 |
| **0.02-0.04 (adopted)** | -0.00 / -0.09 / -0.16 / -0.05 / -0.10 | 0.10 | 0.45 | 31/72 | 2.046 |
| 0.03-0.07 | +0.13 / +0.24 / +0.18 / +0.34 / +0.31 | 0.25 | 0.76 | 29/72 | 2.055 |
| 0.05-0.15 | +0.62 / +1.09 / +1.02 / +1.45 / +1.45 | 1.17 | 1.23 | 29/72 | 2.042 |

Full table (also 0.10-0.20 and 0.25-0.45, hits and false alarms per image): `tables/test_period_growth_scores.csv`.

Reading: 0.02-0.04 tracks the median crest within 0.16 m. The low dunes are slightly worse (P10 RMSE 0.45 vs
0.36), and overwash and shoreline barely change. Faster rates overshoot. **The rate is fitted to the test
window's own surveys, so 2009-2025 no longer tests dune growth**, only overwash and shoreline. The 2009 ends
were solved on 0.01-0.03.

Script: `scripts/hatteras_ms/experiments/HAT_test_period_growth_rate.py` (`run`, `score`). The model runs and
launch logs stay in `output/raw_runs/experiments/dunes/2026-10-20-test-period-growth-rate/`, because the run
registry files every run under `output/raw_runs/`; the 0.01-0.03 reference run is in
`dunes/2026-10-20-shoreline-rescore-new-ends/`. Moved here from the experiment folder on 2026-10-21.
