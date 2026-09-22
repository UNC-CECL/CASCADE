# 3-rates/duneline/endpoint — net dune-line change per window

The stored dune-line observation (2026-09-18): the **net change** between the
two dune lines that bound each window, per 100 m transect and per GIS domain.
Written by `scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py`,
read through `hat_observed_rates.dune_endpoint_csv(start, end, level)`.
`rate_windows.py` draws and scores the dune line from here. It replaced
`duneline_lrr/`, an OLS through every line in the window; see
`../../../archive/duneline_lrr_retired_20260918/WHY.md`.

```
<start>_<end>/
    transect_endpoint.csv         per transect: both line positions (m from the
                                  offshore datum, growing landward), change_m,
                                  rate_m_yr, the line vintages, survey dates,
                                  whether each date is assumed, interval_yr
    domain_endpoint_summary.csv   per domain: n_transects, mean/std/min/max
                                  change_m, mean/std rate_m_yr, pct_landward
    PROVENANCE.md                 the raw files and dates it was built from
```

**Seaward is positive.** `change_m` = start position minus end position.
`rate_m_yr` = `change_m` / the interval between the two survey dates. The 2023
flight date is not known and is assumed to be 1 July. The `rate_m_yr` column
inherits that assumption; `change_m` does not.

Windows: 1984_2004, 1996_2010, 2004_2024, 2010_2024 (the model windows) and
1996_2024 (context). A period year reads its line through
`hat_topo_version.DUNE_LINE_FOR_YEAR` (1996 → 1997, 2010 → 2009,
2024 → 2023). Built from the lines as re-digitized on 2026-09-18.

| window | lines | interval | mean change | domains landward |
|---|---|---|---|---|
| 1984_2004 | 1984 → 2004 | 19.68 yr | −15.4 m (−0.78 m/yr) | 66 / 90 |
| 1996_2010 | 1997 → 2009 | 11.63 yr | −16.2 m (−1.39 m/yr) | 73 / 90 |
| 2004_2024 | 2004 → 2023 | 19.10 yr | −1.7 m (−0.09 m/yr) | 42 / 90 |
| 2010_2024 | 2009 → 2023 | 14.09 yr | +1.4 m (+0.10 m/yr) | 37 / 90 |
| 1996_2024 | 1997 → 2023 | 25.72 yr | −14.8 m (−0.58 m/yr) | 61 / 90 |

Rebuild after any change to a dune line:
`python scripts/input_prep/5-scr/3-rates/duneline/duneline_endpoint.py`

**Figure:** `<window>/duneline_endpoint_<window>.png`, drawn by `scripts/input_prep/5-scr/3-rates/rates_figures.py`, on the metre axis shared by both endpoint products.
