# 2026-09-16-peaisland-ext: results

Interior RMSE is GIS 2-89 against the surveyed CoastSat target in every geometry. End values are what the run imposed (m/yr). "near-end" is the RMSE over GIS 80-90 against the same target. Stage 0 is zeroBE; the solved rows are the last Newton probe.

| member | geometry | mode | preset | step | interior RMSE | interior bias | reach RMSE | near-end RMSE | south end | north end |
|---|---|---|---|---|---|---|---|---|---|---|
| base-asrun | base | asrun | zeroBE | 0 | 1.0852 | -0.2226 | nan | 0.9446 | +0.0 | +0.0 |
| base-asrun | base | asrun | edgeBE | 0 | 1.1407 | +0.0292 | nan | 1.3417 | +32.2 | +10.0 |
| base-check | base | asrun | edgeBE | 0 | 1.1407 | +0.0292 | 1.1407 | 1.3417 | +32.2 | +10.0 |
| base-detrended | base | detrended | zeroBE | 0 | 8.9786 | -3.8218 | 8.9786 | 8.8298 | +0.0 | +0.0 |
| base-detrended | base | detrended | edgeBE | 4 | 9.0251 | -3.9322 | 9.0251 | 9.6716 | +90.0 | -84.4 |
| n115-asrun | n115 | asrun | zeroBE | 0 | 1.1121 | -0.2204 | 1.2224 | 0.9374 | +0.0 | +0.0 |
| n115-asrun | n115 | asrun | edgeBE | 3 | 1.0981 | -0.0180 | 1.0233 | 0.9374 | +37.1 | +42.0 |
| n115-detrended | n115 | detrended | zeroBE | 0 | 8.9531 | -4.4274 | 8.9099 | 9.7619 | +0.0 | +0.0 |
| n115-detrended | n115 | detrended | edgeBE | 4 | 8.6121 | -3.7588 | 9.0464 | 9.7639 | +109.4 | -117.6 |

## Rates on GIS 80-90, m/yr (model LRR; target is the surveyed LOESS)

| GIS | target | base-asrun zeroBE | base-asrun edgeBE | base-check edgeBE | base-detrended zeroBE | base-detrended edgeBE | n115-asrun zeroBE | n115-asrun edgeBE | n115-detrended zeroBE | n115-detrended edgeBE |
|---|---|---|---|---|---|---|---|---|---|---|
| 80 | -2.417 | -1.909 | -1.740 | -1.740 | -17.641 | -19.497 | -1.901 | -1.901 | -18.550 | -18.550 |
| 81 | -2.425 | -1.843 | -1.616 | -1.616 | -15.977 | -18.359 | -1.829 | -1.829 | -17.491 | -17.491 |
| 82 | -2.493 | -1.629 | -1.330 | -1.330 | -12.959 | -15.962 | -1.606 | -1.606 | -15.215 | -15.216 |
| 83 | -2.540 | -1.449 | -1.061 | -1.061 | -10.014 | -13.771 | -1.411 | -1.411 | -13.207 | -13.207 |
| 84 | -2.528 | -1.276 | -0.781 | -0.781 | -6.922 | -11.550 | -1.214 | -1.214 | -11.251 | -11.251 |
| 85 | -2.444 | -1.097 | -0.474 | -0.474 | -3.506 | -9.146 | -0.999 | -0.999 | -9.229 | -9.230 |
| 86 | -2.144 | -0.981 | -0.215 | -0.215 | -0.260 | -7.025 | -0.831 | -0.831 | -7.674 | -7.676 |
| 87 | -1.663 | -0.976 | -0.054 | -0.054 | +2.415 | -5.552 | -0.750 | -0.749 | -7.020 | -7.023 |
| 88 | -1.172 | -1.048 | +0.031 | +0.031 | +5.185 | -3.981 | -0.712 | -0.712 | -6.597 | -6.602 |
| 89 | -0.662 | -1.208 | +0.012 | +0.012 | +7.953 | -2.292 | -0.719 | -0.718 | -6.496 | -6.504 |
| 90 | -0.136 | -1.462 | -0.143 | -0.143 | +10.889 | -0.141 | -0.764 | -0.762 | -6.503 | -6.515 |
