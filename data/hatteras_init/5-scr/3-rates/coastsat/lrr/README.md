# 3-rates/coastsat/lrr — the CoastSat rate fits (the model's scoring target)

One folder per window, `<start>_<end>`: 1984_2004, 1996_2010, 2004_2024,
2010_2024, and 1996_2024 (context only; no run is graded against it).

```
<window>/
    transect_lrr_full.csv    per CoastSat transect: lrr_m_yr (OLS slope of
                             position on date), r_squared, p_value, unc_m_yr
                             (95 % CI half-width), n_obs, first and last date
    domain_lrr_summary.csv   per GIS domain: n_valid, mean_lrr (the value every
                             figure draws), median, std, min, max, pct_eroding,
                             n_transects
    ext/                     1996_2010 and 1996_2024 only: the Pea Island
                             extension (GIS 91-115), from coastsat_extension_lrr.py
```

**The window** runs from 1 January of the start year to 31 December of the
end year, inclusive. A transect needs at least 3 positions. The fit uses
every position with no outlier filter and no weighting. Recent years carry
more weight because they have more positions: about 27 a year per transect
since 2017, about 16 before.

**Seaward positive.** CoastSat chainage grows seaward.

**Read by** the hindcast runner (section 8 builds the scoring target from
`transect_lrr_full.csv`), the edge solve, and every rate comparison, all
through `hat_observed_rates.lrr_csv()` / `COASTSAT_LRR_ROOT`.

**Rebuild a window.** The fit script writes tables only:
```
python scripts/input_prep/5-scr/CoastSat/coastsat_domain_lrr_fixed.py --start-year 1996 --end-year 2010
```
then redraw the figure beside it, `<window>/lrr_<window>.png`, with
`python scripts/input_prep/5-scr/rates_figures.py`. The 2 x 2 of the four model windows is `lrr_four_windows.png` here (moved from `4-comparisons/coastsat_windows/` on 2026-09-19).
