# duneline_lrr/2004_2024 - provenance

Written 2026-09-16 by scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py.

## Surveys in the fit

| vintage | file | survey date | |
|---|---|---|---|
| 2004 | 2004_duneline_offset_raw.csv | 2004-05-25 |  |
| 2009 | 2009_duneline_offset_raw.csv | 2009-05-30 |  |
| 2023 | 2023_duneline_offset_raw.csv | 2023-07-01 | ASSUMED mid-year, flight date unknown |

3 island-wide surveys on 450 transects; ordinary least squares of seaward position (-ORIG_LEN) against decimal survey date per transect. unc_m_yr is the slope's standard error with 1 degree(s) of freedom.

Buxton-only clips (1967, 2017) and the 1978 line are not in the fit; see the script header.

## Island summary

mean of transect rates -0.332 m/yr, median -0.120, 54% of transects landward.
