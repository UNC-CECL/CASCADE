# 3-rates — the observed change products

**Tables, plus one house-style figure per window** beside them. The figures
are drawn by `scripts/input_prep/5-scr/rates_figures.py` after the tables are built. They show the domain value
filled blue for seaward and red for landward, the single transects as dots
coloured by their own sign (the same blue / red), and the village, groin,
pier, shoal and fill marks. The PDF and caption
sit under each window's `supporting/`. Figures that compare products or runs
are in `../4-comparisons/`.

**Per model chain** there is one more figure per product: the chain's two
windows stacked, earlier above, on the same axis as the product's window
figures. The chains are 1984 → 2004 → 2024 and 1996 → 2010 → 2024:

```
coastsat/lrr/chains/lrr_chain_1984_2004_2024.png                     lrr_chain_1996_2010_2024.png
coastsat/endpoint/chains/coastsat_endpoint_chain_1984_2004_2024.png  ..._1996_2010_2024.png
duneline/endpoint/chains/duneline_endpoint_chain_1984_2004_2024.png  ..._1996_2010_2024.png
```

`5yr_bins` has no chain figure: it covers only the 1996 chain, and its
`1996_2024` figure already shows it bin by bin.

Y axes:
- `lrr`: ±8 m/yr, the bound of every window figure.
- The two `endpoint` products: ONE shared metre axis (±120 m, set by
  1984–2004 at Rodanthe), so a shoreline figure reads against its dune-line
  figure.
- `5yr_bins`: one ±18 m/yr axis for all three windows, taken over GIS 2–90.
  Cape Point (GIS 1, about +35 m/yr in the 2020s bins) is clipped at the
  edge, marked with a triangle and labelled with its value.

The folders are grouped **by source**
(2026-09-18, Hannah). Every path is resolved through
`scripts/site_layer/hat_observed_rates.py`, so no script should type one.

```
coastsat/                      the satellite waterline
    lrr/<window>/              OLS rate per transect and per domain       MODEL TARGET
    endpoint/<window>/         net change between ±6-month means at the dune-line dates
    5yr_bins/<window>/         OLS in successive 5-year bins
duneline/                      the digitized dune line
    endpoint/<window>/         net change between the two lines bounding the window
```

| product | what | units | windows | producer (`scripts/input_prep/5-scr/`) | read by |
|---|---|---|---|---|---|
| `coastsat/lrr` | per transect, the OLS slope of CoastSat position against date over the calendar window (1 Jan start to 31 Dec end, ≥ 3 positions); per domain the mean | m/yr | 1984_2004, 1996_2010, 2004_2024, 2010_2024; 1996_2024 (context only); `ext/` for the Pea Island reach | `CoastSat/coastsat_domain_lrr_fixed.py`, `CoastSat/coastsat_extension_lrr.py` | **the hindcast runner** (the scoring target), the edge solve, every rate comparison |
| `coastsat/endpoint` | per transect, the mean position within ±6 months of the END dune-line date minus the same at the START date; per domain the mean | m and m/yr | the four model windows + 1996_2024 | `coastsat_endpoint/coastsat_endpoint.py` | `net_change_1996_2024.py` |
| `coastsat/5yr_bins` | the lrr fit, repeated inside successive 5-year bins (bins under 3.75 yr dropped, \|rate\| > 50 m/yr dropped) | m/yr | 1996_2010, 2010_2024, 1996_2024 | `CoastSat_timeseries/coastsat_lrr_5year_bins.py` | `coastsat_5yr_bins_figure.py` |
| `duneline/endpoint` | per 100 m transect, the start line's distance from the offshore datum minus the end line's; per domain the mean | m and m/yr | the four model windows + 1996_2024 | `duneline_endpoint/duneline_endpoint.py` | `HAT_rate_windows.py`, the dune edge solve, `net_change_1996_2024.py` |

**Sign:** seaward positive in every product.

**Dates:** the dune lines are 1984-09-19, 1997-10-12, 2004-05-25, 2009-05-30
and 2023-07-01. The 2023 date is **assumed**, so every m/yr column that ends
in 2023 inherits it. The metres columns don't depend on it.

**Why two estimators.** The model is graded on the LRR (`coastsat/lrr`), a
trend through every satellite date. The dune line is two surveys, so it can
only be a net change, and it is compared with the CoastSat net change taken
at the same dates (`coastsat/endpoint`). Hannah, 2026-09-18: "we are tracking
net change".

## Changed on 2026-09-18

- The folders were regrouped. Before, they were flat: `coastsat_lrr/`,
  `coastsat_endpoint/`, `coastsat_5yr_bins/`, `duneline_endpoint/`.
- The quick-look PNGs that sat beside each LRR fit
  (`domain_lrr_bar.png`, `transect_lrr_scatter.png`) were archived to
  `../archive/coastsat_lrr_quicklooks_20260918/`. They were autoscaled,
  titled, and coloured by magnitude, so they clashed with the house figures.
  House-style figures replaced them the same evening (`rates_figures.py`),
  one per window for EVERY product; Hannah wanted all of them to have
  figures.
- Two stale summaries were deleted: `domain_lrr_1984_2004_summary.csv` and
  `domain_lrr_2004_2024_summary.csv`, from an older transect mapping, which
  put 29 transects in GIS 1 against 9 now.
- The 5-year bins were rebuilt on 1996 → 2024. The 06-02 run on 1984/2004 is
  in `../archive/coastsat_5yr_bins/`.
- The dune-line OLS (`duneline_lrr/`) was retired for the endpoint product;
  see `../archive/duneline_lrr_retired_20260918/`.
