# 2026-09-18-dune-edgesolve - results

Written 2026-09-18 by scripts/input_prep/7-source-sink/2-calibrate/HAT_be_dune_edgesolve_results.py. See NOTE.md for the question and the layout.

## The solved pairs, GIS 1 / GIS 90, m/yr

| window | reading | step | dune-line solve | target at the ends | CoastSat solve |
|---|---|---|---|---|---|
| 1984-2004 | mean3 | 2 | -1.6 / +6.2 | +0.32 / -0.94 | -42.6 / +13.0 |
| 1984-2004 | raw | 3 | +31.8 / +9.5 | +4.05 / -0.50 | -42.6 / +13.0 |
| 1996-2010 | mean3 | 3 | -16.0 / -3.3 | -2.61 / -2.04 | +32.2 / +10.0 |
| 1996-2010 | raw | 4 | +0.3 / +1.9 | -0.01 / -1.12 | +32.2 / +10.0 |
| 2004-2024 | mean3 | 3 | -0.4 / +23.6 | -1.26 / +0.62 | +50.3 / +46.7 |
| 2004-2024 | raw | 2 | -2.3 / +24.2 | -1.53 / +0.71 | +50.3 / +46.7 |
| 2010-2024 | mean3 | 3 | +20.1 / +19.2 | +1.32 / +1.18 | +72.6 / +31.3 |
| 2010-2024 | raw | 3 | +18.2 / +18.7 | +1.04 / +1.13 | +72.6 / +31.3 |

## Interior skill, GIS 2-89, model minus observation, m/yr

Each run scored against both observations: the CoastSat LRR target (model OLS slope, as run_index.csv) and the dune-line endpoint rate (model endpoint rate, as rate_windows/duneline/endpoint). The CoastSat row is the run the dune solve started from.

| window | solve | ends | vs CoastSat bias | RMSE | vs dune line bias | RMSE |
|---|---|---|---|---|---|---|
| 1984-2004 | coastsat | -42.6 / +13.0 | +0.16 | 1.22 | -0.09 | 2.21 |
| 1984-2004 | dune-mean3 | -1.6 / +6.2 | +0.42 | 1.44 | +0.16 | 2.02 |
| 1984-2004 | dune-raw | +31.8 / +9.5 | +0.62 | 1.83 | +0.37 | 2.15 |
| 1996-2010 | coastsat | +32.2 / +10.0 | +0.02 | 1.13 | +0.99 | 2.32 |
| 1996-2010 | dune-mean3 | -16.0 / -3.3 | -0.37 | 1.22 | +0.60 | 1.87 |
| 1996-2010 | dune-raw | +0.3 / +1.9 | -0.21 | 1.08 | +0.76 | 2.00 |
| 2004-2024 | coastsat | +50.3 / +46.7 | -1.00 | 1.79 | -0.43 | 2.30 |
| 2004-2024 | dune-mean3 | -0.4 / +23.6 | -1.47 | 2.03 | -0.90 | 1.93 |
| 2004-2024 | dune-raw | -2.3 / +24.2 | -1.48 | 2.05 | -0.91 | 1.93 |
| 2010-2024 | coastsat | +72.6 / +31.3 | -1.34 | 2.31 | -0.34 | 2.19 |
| 2010-2024 | dune-mean3 | +20.1 / +19.2 | -1.77 | 2.41 | -0.75 | 1.86 |
| 2010-2024 | dune-raw | +18.2 / +18.7 | -1.78 | 2.42 | -0.77 | 1.86 |
