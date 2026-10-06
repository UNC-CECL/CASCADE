# mean_shoreline/compared -- the three model shorelines, side by side

Written by `scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_compared.py`
on 2026-10-06.

## What this is

The three CoastSat mean shorelines the hindcast starts from and is graded
against, drawn as the island. They are the windows of
`hat_observed_rates.NET_CHANGE_WINDOWS`:

| shoreline | role | centred on | window | folder |
|---|---|---|---|---|
| 1996 | calibration start | 1996-10-12, middle of the ALACE lidar flights (1996-10-09 to 10-16) | ±1 yr | `../1995-10-12_1997-10-12/` |
| 2009 | calibration end, test start | 2009-08-17, middle of the USACE NCMP lidar flights (2009-08-10 to 08-24) | ±1 yr | `../2008-08-17_2010-08-17/` |
| 2025 | test end | 2025-08-17, 16 yr after the 2009 centre; no DEM | ±6 months | `../2025-02-17_2026-02-17/` |

CoastSat ends 2026-01-13, so the 2025 mean covers about 11 months. It still
has the most positions per transect (median 35, against 28 for 2009 and 20
for 1996), because more satellites were flying. Its storm check
(`../2025-02-17_2026-02-17/storm_check/`) finds that Hurricane Erin and two
long nor'easters in fall 2025 pull it landward (median 5.7 m).

Each line is the distance from the shared offshore datum along the 100 m model
transects, averaged per 500 m domain, computed with the offset build's own
intersection (`duneline_to_raw_offsets.intersect`). The 1996 and 2009 stations
equal the stored shoreline-offset raw files (`2-brie-offset/<year>/shoreline/v2/`)
to 0.005 m. No smoothing, no spread band.

## Files

| file | what it is |
|---|---|
| `mean_shoreline_compared_1996_2009_2025.png` | six sections in one row, south to north; 1996 red, 2009 blue, 2025 amber; the range of change per section printed in each panel |
| `supporting/mean_shoreline_compared_1996_2009_2025.csv` | per GIS domain: each line's station from the datum (m, landward +) and the change over each period (m, seaward +) |
| `supporting/mean_shoreline_compared_1996_2009_2025.pdf` | vector copy |
| `supporting/CAPTIONS.md` | caption |

A copy of the PNG is published to `output/figures/2-observations/mean_shoreline/`.
