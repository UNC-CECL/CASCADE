# Buxton 2017 fill: with vs without (2026-10-03)

Question: does adding the 2017-18 Buxton fill (2.6 M cy, GIS 6-15, fired 2017) improve the 2010-2026 hindcast?

Setup: 2010-2026 edgeBE full_management, nogroin, run twice by
`scripts/hatteras_ms/experiments/HAT_buxton_2017_fill.py`. The without2017 member drops the
project from HATTERAS_NOURISHMENT_PROJECTS in-process; the main code is unchanged. The with2017 member
reproduces the matrix run exactly (max |LRR difference| 0).

| | interior RMSE | interior bias | GIS 6-15 RMSE | GIS 6-15 bias |
|---|---|---|---|---|
| with 2017 | 1.836 | -0.970 | 1.580 | +0.782 |
| without 2017 | 1.923 | -1.272 | 2.215 | -1.754 |

m/yr, vs the CoastSat 2010-2026 LOWESS-7 target (GIS 1-10 are raw domain means in that target).

Answer: yes. The fill adds ~2.6 m/yr across GIS 6-15. That fixes GIS 9-15, where the model now sits on
the target. It overshoots GIS 6-8, where CoastSat shows erosion to slight accretion (-1.6 to +0.5) and the
model shows +2.7. Little sand spreads out of the footprint: +0.7 at GIS 5, +0.3 at GIS 16.

Tables: `tables/scores.csv`, `tables/per_domain_lrr.csv`.
