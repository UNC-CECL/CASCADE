# mean_shoreline_windows/2010 -- the 2010 mean shoreline, 3-yr calendar vs 2-yr DEM-centred

Written by `scripts/input_prep/5-scr/4-comparisons/mean_shoreline_windows/coastsat_mean_shoreline_windows.py`
on 2026-09-29.

## The two windows

| | window | folder | positions per transect (median) |
|---|---|---|---|
| 3-yr | calendar 2009–2011 | `1-observations/mean_shoreline/2009_2011/` | 48 |
| 2-yr | 2008-08-17 to 2010-08-17, ±1 yr of the 2009 USACE NCMP topobathy lidar (CHARTS) (flown 2009-08-10 to 2009-08-24; NOAA InPort item 54934, https://www.fisheries.noaa.gov/inport/item/54934) | `1-observations/mean_shoreline/2008-08-17_2010-08-17/` | 28 |

The 3-yr line built the shoreline island offset v1; the 2-yr line builds v2
(Hannah, 2026-09-29). Both are read as stored; nothing is re-averaged here.

## The difference (2-yr minus 3-yr, + seaward)

| | |
|---|---|
| transects in both windows | 906 (not in both: none) |
| per transect | mean +1.67 m, sd 5.22 m, p5 -7.4 m, p95 +10.5 m, largest 15.7 m |
| per domain | mean +1.77 m, range -11.1 to +11.2 m |
| per domain, island mean removed | sd 4.69 m |
| domains beyond one 10 m cell | 5 of 90 |
| median standard error, 3-yr / 2-yr | 2.20 / 2.82 m |

The five domains that move most once the island mean is removed:

| domain | difference (m) | minus island mean (m) | SE 3-yr / 2-yr (m) |
|---|---|---|---|
| GIS 51 | -11.1 | -12.9 | 2.5 / 2.7 |
| GIS 37 | -8.8 | -10.6 | 2.3 / 2.9 |
| GIS 21 | +11.2 | +9.4 | 2.6 / 3.1 |
| GIS 50 | -7.2 | -9.0 | 2.3 / 2.6 |
| GIS 48 | +10.5 | +8.7 | 2.2 / 2.5 |

## How to read it

- **Only the alongshore-varying part reaches the model.** The island offset
  is zeroed on its own minimum, so the island-mean shift (+1.8 m) cancels;
  the 4.7 m sd after removing it is what can change BRIE's shoreline shape.
- **No significance is attached.** The windows overlap and share passes, so
  their means are not independent samples; the standard errors are shown so a
  difference can be read against the sampling noise of either mean.
- **Positive is seaward.** CoastSat chainage grows offshore along every
  transect here, so a positive difference puts the DEM-centred line seaward
  of the calendar one.

## Files

| file | what it is |
|---|---|
| `mean_shoreline_windows_2010.png` | (a) the difference per transect and domain, (b) the standard error of each window's mean |
| `mean_shoreline_windows_2010_datum.png` | both lines as the offset build sees them: distance from the shared offshore datum along the 100 m model transects, per domain, six sections on a 2 x 3 grid (the layout of 2-brie-offset's duneline_vs_shoreline figure), each with a strip of the difference beside it |
| `mean_shoreline_windows_2010_island.png` | the whole island in six north-up segments of 15 domains, the 3-yr line coloured by the difference (blue seaward, red landward) |
| `supporting/datum_stations_2010.csv` | per domain: each window's station from the offshore datum and the seaward shift of the 2-yr line |
| `mean_shoreline_windows_2010_lines_largest.png` | both lines on the photograph nearest the DEM survey, the six domains that move most once the island mean is removed (middle 250 m of each) |
| `mean_shoreline_windows_2010_lines_sites.png` | the same at the centre domains of the six mean_shoreline imagery sites |
| `supporting/domain_comparison_2010.csv` | per domain: n transects, mean / min / max difference, difference minus the island mean, median SE and positions for each window |
| `supporting/transect_comparison_2010.csv` | per transect: both means, n, SE, and the difference |
| `supporting/CAPTIONS.md`, `supporting/*.pdf` | the caption and the vector copy |
