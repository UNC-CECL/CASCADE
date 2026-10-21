# survey-crests: dune crests of the lidar surveys inside the 1996-2009 window

Written by `scripts/input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py --start 1996-start`.
Same method as `../../2009-start/survey-crests/README.md`, on the 1996 start's clips, frame and picked dune
windows (`1-extraction/picks/HAT_dune_search_windows_v2.json`). The script first reproduces the 1996-start
dune files from its own arrays (to 1e-15 m).

The 1996-2005 surveys are the NOAA Data Access Viewer grids as downloaded (`D:\Hatteras_GIS\Elevation\Points`),
3-9 m or ft, reprojected bilinearly onto the 1 m clips; the vertical unit is read from the CRS or from the
land-elevation ratio to the 1996 arrays (`survey_sources.csv`) and converted to metres. Coarse grids blur
individual crests: the 1996 grid against the 1 m 1996 survey matches in median (+0.01 m), 10th percentile
(3.97 vs 3.93 m), share below 4.3 m (15% vs 16%) and domain medians (r 0.99), but single cells only at r 0.52.
Use the distribution and domain values; trust single cells less. Replace with 1 m gridded point clouds when
downloaded. 2001 is unclassified (all returns): buildings or vegetation can raise a crest. 2009 is the 1 m
USACE DEM, on these 1996 windows.

| survey | flown | model state | 10th pct (m) | median (m) | cells < 4.3 m |
|---|---|---|---|---|---|
| 1996 (1 m, the start) | 1996-08 to 10 | 1 Jan 1997 | 3.93 | 5.73 | 16% |
| 1997 | 1997-09 to 10 | 1 Jan 1998 | 3.88 | 5.89 | 16% |
| 1998 | 1998-09 | 1 Jan 1999 | 3.69 | 5.77 | 19% |
| 1999 post-Dennis | 1999-09 | 1 Jan 2000 | 3.36 | 5.83 | 22% |
| 2000 | 2000-02 to 03 | 1 Jan 2000 | 3.13 | 5.58 | 26% |
| 2001 | 2001-01 to 03 | 1 Jan 2001 | 2.91 | 5.40 | 27% |
| 2004 pre-Ivan | 2004-04 to 09 | 1 Jan 2004 | 2.47 | 5.38 | 31% |
| 2005 | 2005-08 to 11 | 1 Jan 2006 | 2.34 | 5.46 | 30% |
| 2009 | 2009-08 | 1 Jan 2010 | 2.35 | 5.39 | 31% |

The low dunes eroded steadily from 1996 to 2004, before and through Isabel, and did not recover by 2009.
