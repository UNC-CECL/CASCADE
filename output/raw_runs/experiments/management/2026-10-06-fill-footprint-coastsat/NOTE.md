# Fill footprints: reported vs CoastSat-observed (2026-10-06)

Question: does placing each fill on the range where CoastSat saw the shoreline move after it
(`nourishment/4-extent-checks/coastsat/nourishment_extent_coastsat_summary.csv`) improve the 2009-2025 test run?

Setup: 2009-2025 test run, domainBE (set 1), blocking groin, full management, relocations off, run twice by
`scripts/hatteras_ms/experiments/HAT_fill_footprint_coastsat.py`. The coastsat member swaps each fill's
domains in-process; volumes are unchanged, so the narrower footprints carry more sand per metre. The main
code is unchanged. The reported member reproduces the matrix run exactly (max |net change difference| 0).

| fill | reported | CoastSat |
|---|---|---|
| Rodanthe 2014 | GIS 82-88 | 83-87 |
| Buxton 2017 | 6-16 | 8-13 |
| Buxton 2022 | 6-16 | 7-9 |
| Avon 2022 | 21-28 | 22-26 |

| | interior bias | interior RMSE | r | GIS 6-16 RMSE | GIS 21-28 RMSE | GIS 82-88 RMSE |
|---|---|---|---|---|---|---|
| reported | -18.2 | 32.1 | 0.36 | 36.4 | 14.1 | 38.9 |
| CoastSat | -18.6 | 35.9 | 0.24 | 54.7 | 12.8 | 43.5 |

Net change in metres against the test target. Interior is GIS 2-89, both sides 7-domain LOWESS (GIS 1-10 raw);
the footprint columns are raw domain values over the union of both footprints.

Answer: no. The CoastSat footprints make the island-wide fit worse (RMSE 32.1 -> 35.9 m, r 0.36 -> 0.24).
At Buxton the two fills stacked on GIS 7-9 build a +56 to +61 m spike where CoastSat shows erosion at GIS 6-8,
and GIS 14-16 lose 30-37 m. Rodanthe loses more at its ends than it gains in the middle. Avon is slightly
better (14.1 -> 12.8 m), the only footprint that improves. The observed ranges are where the change cleared
the noise, not the placement footprint (Buxton 2022 was 7-9 because its 9 m signal barely beat a 7.8 m
background). Not adopted; the reported footprints stay.

Figures (captions in `figures/supporting/CAPTIONS.md`), drawn by `HAT_fill_footprint_coastsat_plot.py`:
- `figures/fill_footprint_reported_vs_coastsat_net_change_2009_2025.png`: net change, island-wide, and the difference between runs.
- `figures/fill_footprint_reported_vs_coastsat_positions_2009_2025.png` (`--positions`): shoreline positions in the three fill reaches (2009 start, observed 2025, both runs), each from a straight baseline fitted to the reach's start so the island's curvature doesn't hide the differences.
Tables: `tables/scores.csv`, `tables/per_domain_net_change.csv`.
