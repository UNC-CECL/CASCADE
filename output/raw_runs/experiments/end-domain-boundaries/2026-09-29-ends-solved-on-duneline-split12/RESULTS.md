# end-domain-boundaries/2026-09-29-ends-solved-on-duneline-split12 - results

Written 2026-09-29 by scripts/input_prep/7-source-sink/2-calibrate/be_dune_edgesolve_results.py. See NOTE.md for the question and the layout.

## The solved pairs, GIS 1 / GIS 90, m/yr

| window | reading | step | dune-line solve | target at the ends | CoastSat solve |
|---|---|---|---|---|---|
| 1996-2010 | mean3 | 3 | -2.7 / +8.0 | -2.61 / -2.04 | +4.4 / +19.1 |
| 1996-2010 | raw | 3 | +0.6 / +10.7 | -0.01 / -1.12 | +4.4 / +19.1 |
| 2010-2024 | mean3 | 5 | +2.2 / +13.8 | +1.32 / +1.18 | +8.0 / +21.3 |
| 2010-2024 | raw | 4 | +2.0 / +13.7 | +1.04 / +1.13 | +8.0 / +21.3 |

## Interior skill, GIS 2-89, model minus observation, m/yr

Each run scored against both observations: the CoastSat LRR target (model OLS slope, as run_index.csv) and the dune-line endpoint rate (model endpoint rate, as model_vs_observed/vs_duneline/endpoint_net_change). The CoastSat row is the run the dune solve started from.

| window | solve | ends | vs CoastSat bias | RMSE | vs dune line bias | RMSE |
|---|---|---|---|---|---|---|
| 1984-2004 | coastsat | -42.6 / +13.0 | +0.17 | 1.34 | -0.09 | 2.21 |
| 1996-2010 | coastsat | +4.4 / +19.1 | +0.06 | 1.17 | +1.08 | 2.03 |
| 1996-2010 | dune-mean3 | -2.7 / +8.0 | +0.01 | 1.20 | +1.04 | 1.96 |
| 1996-2010 | dune-raw | +0.6 / +10.7 | +0.03 | 1.18 | +1.05 | 1.98 |
| 2004-2024 | coastsat | +50.3 / +46.7 | -1.04 | 1.90 | -0.43 | 2.30 |
| 2010-2024 | coastsat | +8.0 / +21.3 | -1.37 | 2.07 | -0.28 | 1.58 |
| 2010-2024 | dune-mean3 | +2.2 / +13.8 | -1.40 | 2.09 | -0.31 | 1.58 |
| 2010-2024 | dune-raw | +2.0 / +13.7 | -1.40 | 2.10 | -0.31 | 1.58 |
