# duneline_lrr - a linear regression rate through the dune lines

One folder per window, the `coastsat_lrr/` layout (`transect_lrr_full.csv`,
`domain_lrr_summary.csv`, `PROVENANCE.md`), so `CoastSatDataset`,
`build_coastsat_series` and `build_target_table` read it unchanged. Written
2026-09-16 by `scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py`
(Hannah: "score it using the other dune line geojsons I have measured to
make LRR rates").

Per 100 m transect, an ordinary-least-squares slope of the dune-line station
(seaward positive) against survey date over every ISLAND-WIDE line inside the
window, start and end vintages inclusive through `DUNE_LINE_FOR_YEAR`:

| window | surveys | n |
|---|---|---|
| 1984-2004 | 1984, 1997, 2004 | 3 |
| 1996-2010 | 1997, 2004, 2009 | 3 |
| 2004-2024 | 2004, 2009, 2023 | 3 |
| 2010-2024 | 2009, 2023 | 2 (the endpoint rate) |

The Buxton-only clips (1967, 2017, GIS 2-12, ArcGIS exports about a metre
landward of the shapely build) and the island-wide 1978 line (before every
window) are not in the fit, so every domain in a window has the same design
and every line the same method. Dates: `duneline_vs_coastsat.KNOWN_SURVEY_DATES`;
the 2023 flight date is unknown and assumed 2023-07-01, flagged in each
window's PROVENANCE.md.

Read through `hat_observed_rates.dune_lrr_csv(start, end)`. Drawn against
the model in `output/comparisons/rate_windows/duneline/lrr/`.
