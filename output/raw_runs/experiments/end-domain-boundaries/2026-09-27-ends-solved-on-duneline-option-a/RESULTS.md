# end-domain-boundaries/2026-09-27-ends-solved-on-duneline-option-a - results

Written 2026-09-27 by scripts/input_prep/7-source-sink/2-calibrate/be_dune_edgesolve_results.py. See NOTE.md for the question and the layout.

## The solved pairs, GIS 1 / GIS 90, m/yr

| window | reading | step | dune-line solve | target at the ends | CoastSat solve |
|---|---|---|---|---|---|
| 1996-2010 | mean3 | 4 | -3.0 / +7.6 | -2.61 / -2.04 | +4.8 / +17.5 |
| 1996-2010 | raw | 2 | +0.3 / +10.6 | -0.01 / -1.12 | +4.8 / +17.5 |
| 2010-2024 | mean3 | 5 | +3.4 / +15.0 | +1.32 / +1.18 | +18.8 / +24.5 |
| 2010-2024 | raw | 6 | +3.1 / +14.8 | +1.04 / +1.13 | +18.8 / +24.5 |

## Interior skill, GIS 2-89, model minus observation, m/yr

Each run scored against both observations: the CoastSat LRR target (model OLS slope, as run_index.csv) and the dune-line endpoint rate (model endpoint rate, as model_vs_observed/vs_duneline/endpoint_net_change). The CoastSat row is the run the dune solve started from.

| window | solve | ends | vs CoastSat bias | RMSE | vs dune line bias | RMSE |
|---|---|---|---|---|---|---|
| 1984-2004 | coastsat | -42.6 / +13.0 | +0.16 | 1.22 | -0.09 | 2.21 |
| 1996-2010 | coastsat | +4.8 / +17.5 | +0.09 | 1.05 | +1.07 | 2.04 |
| 1996-2010 | dune-mean3 | -3.0 / +7.6 | +0.05 | 1.08 | +1.02 | 1.98 |
| 1996-2010 | dune-raw | +0.3 / +10.6 | +0.07 | 1.05 | +1.05 | 2.01 |
| 2004-2024 | coastsat | +50.3 / +46.7 | -1.00 | 1.79 | -0.43 | 2.30 |
| 2010-2024 | coastsat | +18.8 / +24.5 | -1.64 | 2.25 | -0.63 | 1.65 |
| 2010-2024 | dune-mean3 | +3.4 / +15.0 | -1.67 | 2.28 | -0.66 | 1.67 |
| 2010-2024 | dune-raw | +3.1 / +14.8 | -1.68 | 2.28 | -0.66 | 1.67 |
