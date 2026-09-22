# 2026-09-16-dune-edgesolve - results

Written 2026-09-16 by scripts/input_prep/7-source-sink/2-calibrate/be_dune_edgesolve_results.py. See NOTE.md for the question and the layout.

## The solved pairs, GIS 1 / GIS 90, m/yr

| window | reading | step | dune-line solve | target at the ends | CoastSat solve |
|---|---|---|---|---|---|
| 1984-2004 | mean3 | 2 | -1.6 / +6.2 | +0.32 / -0.94 | -42.6 / +13.0 |
| 1984-2004 | raw | 3 | +31.8 / +9.5 | +4.05 / -0.50 | -42.6 / +13.0 |
| 1996-2010 | mean3 | 4 | -28.6 / -3.1 | -4.67 / -2.04 | +32.2 / +10.0 |
| 1996-2010 | raw | 2 | -2.7 / +2.1 | -0.40 / -1.12 | +32.2 / +10.0 |
| 2004-2024 | mean3 | 3 | -6.2 / +15.9 | -2.02 / -0.36 | +50.3 / +46.7 |
| 2004-2024 | raw | 2 | -2.7 / +18.1 | -1.58 / -0.07 | +50.3 / +46.7 |
| 2010-2024 | mean3 | 2 | +27.2 / +10.8 | +2.39 / -0.16 | +72.6 / +31.3 |
| 2010-2024 | raw | 2 | +27.8 / +12.3 | +2.48 / +0.08 | +72.6 / +31.3 |

## Interior skill, GIS 2-89, model minus observation, m/yr

Each run scored against both observations: the CoastSat LRR target (model OLS slope, as run_index.csv) and the dune-line endpoint rate (model endpoint rate, as duneline_vs_modeled_windows). The CoastSat row is the run the dune solve started from.

| window | solve | ends | vs CoastSat bias | RMSE | vs dune line bias | RMSE |
|---|---|---|---|---|---|---|
| 1984-2004 | coastsat | -42.6 / +13.0 | +0.16 | 1.22 | -0.09 | 2.21 |
| 1984-2004 | dune-mean3 | -1.6 / +6.2 | +0.42 | 1.44 | +0.16 | 2.02 |
| 1984-2004 | dune-raw | +31.8 / +9.5 | +0.62 | 1.83 | +0.37 | 2.15 |
| 1996-2010 | coastsat | +32.2 / +10.0 | +0.03 | 1.14 | +0.90 | 2.95 |
| 1996-2010 | dune-mean3 | -28.6 / -3.1 | -0.47 | 1.45 | +0.41 | 2.35 |
| 1996-2010 | dune-raw | -2.7 / +2.1 | -0.23 | 1.11 | +0.64 | 2.56 |
| 2004-2024 | coastsat | +50.3 / +46.7 | -1.00 | 1.79 | -0.07 | 2.29 |
| 2004-2024 | dune-mean3 | -6.2 / +15.9 | -1.57 | 2.11 | -0.64 | 1.73 |
| 2004-2024 | dune-raw | -2.7 / +18.1 | -1.53 | 2.07 | -0.60 | 1.74 |
| 2010-2024 | coastsat | +72.6 / +31.3 | -1.33 | 2.31 | -0.31 | 2.01 |
| 2010-2024 | dune-mean3 | +27.2 / +10.8 | -1.77 | 2.37 | -0.74 | 1.69 |
| 2010-2024 | dune-raw | +27.8 / +12.3 | -1.75 | 2.36 | -0.72 | 1.69 |
