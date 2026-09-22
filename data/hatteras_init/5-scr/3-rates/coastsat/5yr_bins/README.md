# 3-rates/coastsat/5yr_bins — the CoastSat rate in 5-year bins

When, inside a window, did the change happen? The same per-transect OLS as
`../lrr/`, fitted separately inside successive 5-year bins of each window.
Rebuilt 2026-09-18 on the canonical chain.

```
1996_2010/lrr_bins_5yr.csv    1996–2000, 2001–2005, 2006–2010
2010_2024/lrr_bins_5yr.csv    2010–2014, 2015–2019, 2020–2024
1996_2024/lrr_bins_5yr.csv    1996–2000 … 2016–2020, 2021–2024
```

Rows are bins and columns are GIS domains; values are the domain MEAN of the
transect rates, in m/yr, seaward positive.

**Rules:** calendar years inclusive, with whole-year bins; a bin shorter than
3.75 years is dropped. A transect needs at least 3 positions in a bin, and a
transect rate beyond ±50 m/yr is dropped. `compute_lrr` is the same function
as `../lrr/` uses.

**A caution.** A 5-year fit goes through about 80–130 noisy positions, so a
single-domain spike is weak evidence. The largest value is Cape Point (GIS 1)
at +34 to +36 m/yr in the last bin, the 2021–24 shoal attachment. The
figures share a ±18 m/yr axis taken over GIS 2–90 (GIS 2 reaches about +16
in the 2020s), and GIS 1 is clipped at the top edge, marked with a triangle
and labelled with its value (Hannah, 2026-09-18).

Producer: `scripts/input_prep/5-scr/3-rates/coastsat/5yr_bins/coastsat_5yr_bins.py`
(tables). Figure: `<window>/lrr_5yr_bins_<window>.png` beside the table, drawn
by `coastsat_5yr_bins_figure.py`, which `rates_figures.py` calls. The 1984 / 2004 runs of 2026-05/06 are in
`5-scr/archive/coastsat_5yr_bins/`.
