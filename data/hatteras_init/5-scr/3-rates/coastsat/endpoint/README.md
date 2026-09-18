# 3-rates/coastsat/endpoint — net CoastSat shoreline change per window

The CoastSat counterpart of `../../duneline/endpoint/` (2026-09-18). For each
CoastSat transect, the mean shoreline position within **±6 months of each of
the window's dune-line survey dates**, end minus start, **seaward positive**,
in metres and m/yr. The domain value is the mean over the domain's transects.
Written by `scripts/input_prep/5-scr/coastsat_endpoint/coastsat_endpoint.py`,
read through `hat_observed_rates.coastsat_endpoint_csv(start, end, level)`.

```
<start>_<end>/
    transect_endpoint.csv         per transect: both window means, positions in
                                  each window, first/last date inside each,
                                  change_m, rate_m_yr, the survey dates
    domain_endpoint_summary.csv   per domain: n, mean/std/min/max change_m,
                                  mean rate, pct_landward, median positions per
                                  end window, transects with an empty window
    PROVENANCE.md
```

**Why the dune dates.** The windows are centred on the same moments as the
dune lines, so this product and the dune product difference like for like.
A gap between them is then beach-width change, not a date mismatch. The
window means come from the functions `duneline_vs_coastsat.py` uses (imported,
not copied), and its CoastSat endpoint agrees with this to 0.0000 m/yr.

The dates are 1984-09-19, 1997-10-12, 2004-05-25, 2009-05-30 and 2023-07-01.
The 2023 date is **assumed**, so the end window is centred on it. The 1984
start window is one-sided, because the CoastSat record begins 1984-09-21.

| window | lines | positions per end window | mean change | domains landward |
|---|---|---|---|---|
| 1984_2004 | 1984 → 2004 | 7 / 20 | −19.4 m (−0.99 m/yr) | 62 / 90 |
| 1996_2010 | 1997 → 2009 | 6 / 9 | −2.0 m (−0.17 m/yr) | 47 / 90 |
| 2004_2024 | 2004 → 2023 | 20 / 39 | +7.8 m (+0.41 m/yr) | 38 / 90 |
| 2010_2024 | 2009 → 2023 | 9 / 39 | +7.0 m (+0.50 m/yr) | 42 / 90 |
| 1996_2024 | 1997 → 2023 | 6 / 39 | +5.0 m (+0.20 m/yr) | 40 / 90 |

No transect had an empty end window in any window (906 transects each). The
1997 and 2009 ends rest on only 6 and 9 positions, so a storm inside either
window moves that end more than it moves the 39-position 2023 end.

**Figure:** `<window>/coastsat_endpoint_<window>.png`, drawn by `scripts/input_prep/5-scr/rates_figures.py`, on the metre axis shared by both endpoint products.
