# mean_shoreline_windows/1996 -- the 1996 mean shoreline, 3-yr calendar vs 2-yr DEM-centred

Written by `scripts/input_prep/5-scr/4-comparisons/mean_shoreline_windows/coastsat_mean_shoreline_windows.py`
on 2026-09-29.

## The two windows

| | window | folder | positions per transect (median) |
|---|---|---|---|
| 3-yr | calendar 1995–1997 | `1-observations/mean_shoreline/1995_1997/` | 28 |
| 2-yr | 1995-10-12 to 1997-10-12, ±1 yr of the 1996 fall East Coast NOAA/NASA ALACE lidar (flown 1996-10-09 to 1996-10-16; NOAA InPort item 48147, https://www.fisheries.noaa.gov/inport/item/48147) | `1-observations/mean_shoreline/1995-10-12_1997-10-12/` | 20 |

The 3-yr line built the shoreline island offset v1; the 2-yr line builds v2
(Hannah, 2026-09-29). Both are read as stored; nothing is re-averaged here.

## The difference (2-yr minus 3-yr, + seaward)

| | |
|---|---|
| transects in both windows | 905 (not in both: none) |
| per transect | mean +0.62 m, sd 3.09 m, p5 -4.5 m, p95 +6.3 m, largest 10.7 m |
| per domain | mean +0.62 m, range -6.3 to +8.6 m |
| per domain, island mean removed | sd 2.74 m |
| domains beyond one 10 m cell | 0 of 90 |
| median standard error, 3-yr / 2-yr | 2.44 / 2.37 m |

The five domains that move most once the island mean is removed:

| domain | difference (m) | minus island mean (m) | SE 3-yr / 2-yr (m) |
|---|---|---|---|
| GIS 30 | +8.6 | +8.0 | 3.5 / 2.8 |
| GIS 76 | -6.3 | -6.9 | 5.0 / 5.0 |
| GIS 20 | +7.2 | +6.6 | 3.4 / 2.9 |
| GIS 57 | +6.2 | +5.6 | 3.0 / 2.5 |
| GIS 29 | +6.2 | +5.6 | 2.8 / 2.3 |

## How to read it

- **Only the alongshore-varying part reaches the model.** The island offset
  is zeroed on its own minimum, so the island-mean shift (+0.6 m) cancels;
  the 2.7 m sd after removing it is what can change BRIE's shoreline shape.
- **No significance is attached.** The windows overlap and share passes, so
  their means are not independent samples; the standard errors are shown so a
  difference can be read against the sampling noise of either mean.
- **Positive is seaward.** CoastSat chainage grows offshore along every
  transect here, so a positive difference puts the DEM-centred line seaward
  of the calendar one.

## Files

| file | what it is |
|---|---|
| `mean_shoreline_windows_1996.png` | (a) the difference per transect and domain, (b) the standard error of each window's mean |
| `mean_shoreline_windows_1996_datum.png` | both lines as the offset build sees them: distance from the shared offshore datum along the 100 m model transects, per domain, six sections on a 2 x 3 grid (the layout of 2-brie-offset's duneline_vs_shoreline figure), each with a strip of the difference beside it |
| `mean_shoreline_windows_1996_island.png` | the whole island in six north-up segments of 15 domains, the 3-yr line coloured by the difference (blue seaward, red landward) |
| `supporting/datum_stations_1996.csv` | per domain: each window's station from the offshore datum and the seaward shift of the 2-yr line |
| `mean_shoreline_windows_1996_lines_largest.png` | both lines on the photograph nearest the DEM survey, the six domains that move most once the island mean is removed (middle 250 m of each) |
| `mean_shoreline_windows_1996_lines_sites.png` | the same at the centre domains of the six mean_shoreline imagery sites |
| `supporting/domain_comparison_1996.csv` | per domain: n transects, mean / min / max difference, difference minus the island mean, median SE and positions for each window |
| `supporting/transect_comparison_1996.csv` | per transect: both means, n, SE, and the difference |
| `supporting/CAPTIONS.md`, `supporting/*.pdf` | the caption and the vector copy |
