# survey-crests: dune crests of the lidar surveys inside the 2009-2025 window

Written by `scripts/input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py`.

The dune crest (highest point across the dune) at every 10 m cell, GIS 1-90, from five later surveys, taken
exactly as the 2009 start's dunes were: the 1 m mosaic cut on the 2009-2014 product's domain clips (same pixel
grid), resampled to 10 m by the same rule, straightened in the 2009 frame, and the maximum inside the 2009
picked dune window (`1-extraction/picks/HAT_dune_search_windows_v1.json`), clipped to the survey's own beach.
The script first runs the 2009 arrays through the same code and stops unless it reproduces the 2009 dune
files (it does, to 1e-15 m).

| survey | mosaic (D:\hatteras_ch1\gis_sources\elevation\dem, folder `noaa_<job>`) | flown | compared with the model state of |
|---|---|---|---|
| 2014 | 2014_NOAA_Post_Sandy_DEM_Job1076021 | 2013-11 to 2014-06 | 1 Jan 2014 |
| 2016 | 2016_USACE_post_matthew_Job1076022 | 2016-10 to 2016-12 | 1 Jan 2017 |
| 2017 | usace2017_east_coast_dem_Job1076023 | 2017-06 to 2017-08 | 1 Jan 2017 |
| 2018 | usace2018_post_florence_dem_Job1076024 | 2018-10 | 1 Jan 2019 |
| 2019 | nc2019_dunex_Job1076025 | 2019-09 to 2019-10 | 1 Jan 2020 |

Files:
- `survey_dune_crests.csv`: one row per survey, domain and cell. `crest_m_navd88` is the crest;
  `status` is `ok`, or `beach_landward_of_window` where the 2009 dune position is now beach (the dune is gone;
  no crest written). `crest_widened_window_m_navd88` widens the window 20 m each side, a sensitivity only:
  on the 2009 surface itself widening raises 10% of crests by 1.3 m or more (taller dunes just landward).
- `survey_coverage.csv`: per survey and domain, cells with a crest and why the others have none.

Every cell is surveyed in every survey; the only missing crests are dunes that became beach (115, 44, 20, 88
and 108 of 4,500 cells in 2014-2019; by 2019 all of GIS 83-84).

What the surveys show (2009 -> 2019): the median crest rose from 5.39 to 5.93 m NAVD88, the low dunes did not
(10th percentile 3.05 -> 2.8-3.2 m), about a quarter of cells stayed below 4.3 m throughout, and the
alongshore pattern held (domain medians 2009 vs 2019, r 0.82).
