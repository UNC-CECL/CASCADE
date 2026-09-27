# CoastSat shoreline change rates, Cape Point to the Virginia line, 1984-2025

One rate per CoastSat transect over the whole satellite record, for the
2,984 transects from the south end of the Hatteras model domain
(`usa_NC_0032_0021`, Cape Point) north to the North Carolina / Virginia line
(`usa_NC_0049_0230`, 36.550 N). Prepared 2026-09-27
for the Murray lab.

**File:** `coastsat_lrr_obx_1984_2025.csv`, one row per transect, south to north: the
transect ID, where it is, its rate and the rate's uncertainty. Every other
fit statistic and CoastSat field is in `supporting/coastsat_lrr_obx_1984_2025_full.csv`.

## How the rate is measured

- **Data:** CoastSat shoreline time series (Vos et al.), downloaded per
  transect from the CoastSat portal (coastsat.space), sites `usa_NC_0032` to
  `usa_NC_0049`: one cross-shore position (chainage, metres along the
  transect) per usable satellite image. CoastSat tidally corrects these
  positions itself, with the FES2022 tide model and a satellite-derived beach
  slope for each transect. The positions used run 1984-04-14 to
  2025-12-28; the downloaded files continue into January 2026, so they come
  from the continuously updated portal and extend past the archived US East
  Coast release (Zenodo v1.0, 9 June 2025,
  doi:10.5281/zenodo.15626280, CC-BY-4.0, whose description gives the
  processing above). Please cite CoastSat if you use these numbers.
- **Window:** 1 January 1984 to 31 December 2025, the whole
  record through the last full calendar year. The first images are from 1984;
  the few January 2026 positions are left out.
- **Rate:** linear regression rate (LRR), the ordinary least-squares slope of
  shoreline position against time, in m/yr. Every position in the window is
  used: **we apply no outlier filter and no weighting** on top of CoastSat's
  own processing. A transect needs at least
  3 positions to get a rate.
- **Sign:** CoastSat chainage increases seaward, so **positive = seaward
  movement (accretion), negative = landward (erosion).**
- This is the same fit, on the same software, that produces the rate the
  Hatteras CASCADE model is scored against, so inside the model domain these
  numbers are directly comparable to that work (over a different window).

Things worth knowing about the record:

- Sampling is uneven in time. Averaged over these transects: about 8
  positions a year in 1984-1998, 18 a year in 1999-2020, and 29 a year in
  2021-2025. An unweighted OLS fit therefore leans toward recent years.
- `unc_m_yr` and `p_value` assume independent residuals. Shoreline position
  has seasonal and storm-driven memory, so consecutive positions are not
  independent and the true uncertainty is wider than stated; treat
  `unc_m_yr` as a lower bound and `p_value` as optimistic.
- Along Hatteras Island there is a real, abrupt seaward shift of about 17 m
  in 2021 (confirmed against the dune line, not a nourishment artefact). Over
  a 42-year fit it moves the rate much less than it moves
  a 10-15 year window; if you cut the record into sub-periods, it matters.

## Columns

`coastsat_lrr_obx_1984_2025.csv`:

| column | meaning |
|---|---|
| `transect_id` | CoastSat transect ID, `usa_NC_<site>_<transect>` |
| `alongshore_km` | distance north along the coast from the origin of `usa_NC_0032_0021` (km) |
| `origin_lon`, `origin_lat` | landward end of the transect (WGS84, decimal degrees) |
| `lrr_m_yr` | shoreline change rate (m/yr; positive = seaward/accretion, negative = landward/erosion) |
| `unc_m_yr` | 95 % confidence half-width on the rate (m/yr); see the caveat above |
| `flag` | reasons to look twice, `;`-separated; empty when none (see Flags) |

`supporting/coastsat_lrr_obx_1984_2025_full.csv` has those columns plus:

| column | meaning |
|---|---|
| `site_id` | CoastSat site |
| `seaward_lon`, `seaward_lat` | seaward end of the transect (WGS84) |
| `r_squared` | R² of the fit (dimensionless) |
| `p_value` | p-value of the slope (two-sided, null: zero trend), 3 significant figures |
| `n_obs` | positions used in the fit |
| `start_date`, `end_date` | first and last position used (`YYYY-MM-DD`, UTC) |
| `span_yr` | years between them |
| `beach_slope` | CoastSat's satellite-derived beach-face slope (tan β, dimensionless), the value CoastSat used to tidally correct this transect |
| `beach_slope_ci_lower`, `beach_slope_ci_upper` | the confidence interval CoastSat reports on that slope (its `cil` / `ciu`) |

The three beach-slope columns are copied unchanged from the CoastSat
transect layer; nothing here estimates them. CoastSat reports the slope in
discrete steps (median 0.07, range 0.015-0.18 over these transects), and
neighbouring transects often share a value.

## Flags

Flagged transects are **kept, with their rate**; the flag says why a reader
might treat them with care. 103 of 2,984 transects carry one or more.

| flag | transects |
|---|---|
| `within_2_km_of_oregon_inlet` | 60 |
| `within_2_km_of_cape_point` | 43 |
| `fewer_than_50_positions` | 0 |
| `span_under_20_yr` | 0 |
| `no_fit` | 0 |

- `fewer_than_50_positions`: too few positions for a stable slope.
- `span_under_20_yr`: the positions cover less than
  20 years of the window, so the rate is not the full-record rate.
- `within_2_km_of_oregon_inlet`: inlet-flank shoreline (Bodie Island
  spit and the north end of Pea Island), where inlet migration and the
  terminal groin, not open-coast processes, set the rate.
- `within_2_km_of_cape_point`: the first 2 km north of the
  tip of Cape Hatteras, where the shoreline turns and the cape shoals shelter it.
- `no_fit`: fewer than 3 positions; no rate.

## Alongshore distance

Each step between neighbouring transect origins is projected onto the local
shore-parallel direction (perpendicular to the transects' bearing) and
summed. That keeps the landward origins' cross-shore wander (up to ~300 m on
the Currituck Banks) from being counted as length. Oregon Inlet sits at
62.0 km, the one real gap in the transects (~1.1 km). The transect
spacing is nominally 50 m; around `usa_NC_0049_0113`-`0116`, in the
northernmost site, CoastSat's transects are 300-570 m apart.

## Summary

| | south of Oregon Inlet | north of Oregon Inlet | all |
|---|---|---|---|
| transects | 1231 | 1753 | 2984 |
| km | 0-61.4 | 62.5-153.5 | 153.5 |
| median LRR (m/yr) | -0.33 | +0.29 | +0.04 |
| eroding (%) | 62 | 39 | 49 |

2,984 of 2,984 transects have a rate.

## Figures

- `lrr_obx_1984_2025.png`: every transect's rate against alongshore distance,
  coloured red (erosion) to blue (accretion), with place names along the top.
  The black line ("1 km median") is drawn over the points as a guide; it
  is not applied to the data. At each transect it is the median rate of all
  transects within 0.5 km either side (about 21), computed separately on each
  side of Oregon Inlet. It is left blank where fewer than 5 transects fall in
  that window, which happens only at about 144-148 km (see Alongshore
  distance). Twelve transects just north of Oregon Inlet are faster than
  -8 m/yr and sit off the axis, marked with a triangle; their values are in the
  CSV.
- `lrr_obx_1984_2025_map_overview.png`: map of the whole reach,
  transects coloured by rate, beside the rates plotted against northing;
  boxes A-D mark the regional maps.
- `..._map_A_cape_point_to_salvo.png`, `..._map_B_rodanthe_to_south_nags_head.png`,
  `..._map_C_nags_head_to_duck.png`, `..._map_D_corolla_to_virginia.png`:
  the same two panels zoomed to ~35-40 km each (1 km overlap). The colour
  scale is fixed at +/-3 m/yr in every map, so colours compare across them.
- The figures draw every transect the same way; the flags are in the table
  only.

Captions in `supporting/CAPTIONS.md`, PDFs in `supporting/`. The maps are
drawn by `coastsat_obx_lrr_maps.py` (basemap: Esri World Shaded Relief,
which carries no labels, shown in greyscale; every name on a map is one
placed here). Latitude is marked on each map's left edge; it is exact at
the coast and the (b) panels share that axis.

## Reproduce

```
python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr.py        # table, README, profile figure
python scripts/input_prep/5-scr/3-rates/coastsat/lrr/coastsat_obx_lrr_maps.py   # the five maps (needs internet for the basemap)
```

in the CASCADE repository (branch `hannahaline/hatteras-cascade`). The fit
is `coastsat_lrr.compute_lrr` (`scripts/input_prep/5-scr/lib/`); the
transect geometry is the CoastSat global transect layer.
