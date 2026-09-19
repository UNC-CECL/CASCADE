# 2026-09-16-peaisland-ext: results

Interior RMSE is GIS 2-89 against the surveyed CoastSat target in every geometry. End values are what the run imposed (m/yr). "near-end" is the RMSE over GIS 80-90 against the same target. Stage 0 is zeroBE; the solved rows are the last Newton probe.

| member | geometry | mode | preset | step | interior RMSE | interior bias | reach RMSE | near-end RMSE | south end | north end |
|---|---|---|---|---|---|---|---|---|---|---|
| base-asrun | base | asrun | zeroBE | 0 | 1.0705 | -0.2278 | 1.0705 | 0.9189 | +0.0 | +0.0 |
| base-asrun | base | asrun | edgeBE | 0 | 1.1288 | +0.0232 | 1.1288 | 1.3224 | +32.2 | +10.0 |
| base-check | base | asrun | edgeBE | 0 | 1.1288 | +0.0232 | 1.1288 | 1.3224 | +32.2 | +10.0 |
| base-detrended | base | detrended | zeroBE | 0 | 8.9122 | -3.6803 | 8.9122 | 9.1462 | +0.0 | +0.0 |
| base-detrended | base | detrended | edgeBE | 4 | 9.0296 | -3.8898 | 9.0296 | 9.8270 | +81.7 | -90.9 |
| n115-asrun | n115 | asrun | zeroBE | 0 | 1.0968 | -0.2252 | 1.2085 | 0.9151 | +0.0 | +0.0 |
| n115-asrun | n115 | asrun | edgeBE | 3 | 1.0869 | -0.0240 | 1.0140 | 0.9151 | +36.9 | +42.0 |
| n115-detrended | n115 | detrended | zeroBE | 0 | 8.8655 | -4.3457 | 8.8334 | 9.8040 | +0.0 | +0.0 |
| n115-detrended | n115 | detrended | edgeBE | 5 | 8.5888 | -3.7122 | 9.0215 | 9.8053 | +100.0 | -122.7 |

## Rates on GIS 80-90, m/yr (model LRR; target is the surveyed LOESS)

| GIS | target | base-asrun zeroBE | base-asrun edgeBE | base-check edgeBE | base-detrended zeroBE | base-detrended edgeBE | n115-asrun zeroBE | n115-asrun edgeBE | n115-detrended zeroBE | n115-detrended edgeBE |
|---|---|---|---|---|---|---|---|---|---|---|
| 80 | -2.417 | -2.072 | -1.903 | -1.903 | -18.364 | -20.345 | -2.064 | -2.064 | -19.356 | -19.356 |
| 81 | -2.425 | -1.927 | -1.701 | -1.701 | -16.086 | -18.611 | -1.913 | -1.913 | -17.693 | -17.693 |
| 82 | -2.493 | -1.682 | -1.384 | -1.384 | -12.748 | -15.934 | -1.659 | -1.658 | -15.121 | -15.122 |
| 83 | -2.540 | -1.475 | -1.088 | -1.088 | -9.657 | -13.661 | -1.436 | -1.436 | -13.005 | -13.005 |
| 84 | -2.528 | -1.285 | -0.791 | -0.791 | -6.488 | -11.430 | -1.222 | -1.222 | -11.017 | -11.018 |
| 85 | -2.444 | -1.098 | -0.477 | -0.477 | -3.028 | -9.070 | -1.000 | -1.000 | -9.023 | -9.024 |
| 86 | -2.144 | -1.003 | -0.238 | -0.238 | +0.031 | -7.238 | -0.853 | -0.853 | -7.728 | -7.730 |
| 87 | -1.663 | -0.967 | -0.045 | -0.045 | +3.032 | -5.540 | -0.742 | -0.741 | -6.813 | -6.816 |
| 88 | -1.172 | -1.032 | +0.047 | +0.047 | +5.905 | -3.966 | -0.700 | -0.699 | -6.357 | -6.360 |
| 89 | -0.662 | -1.187 | +0.033 | +0.033 | +8.762 | -2.284 | -0.704 | -0.703 | -6.239 | -6.244 |
| 90 | -0.136 | -1.436 | -0.117 | -0.117 | +11.770 | -0.137 | -0.750 | -0.748 | -6.240 | -6.246 |
