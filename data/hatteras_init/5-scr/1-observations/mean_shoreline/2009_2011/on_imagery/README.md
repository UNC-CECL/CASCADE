# on_imagery (2009-2011): one panel, from outside the window

The 1995-1997 figures have one panel per flight year. These have one panel,
the NOAA NGS colour orthomosaic flown 2008-03-26/27 (InPort 48695), because
no georeferenced photographs from 2009, 2010 or 2011 are on the imagery
drive (checked 2026-09-28):

- **2009**: `D:\Hatteras_GIS\Aerial\2009_Henry\raw_GE\HAT_2009_*.jpg` are
  Google Earth screenshots with no georeferencing, so the line cannot be
  placed on them.
- **2010, 2011**: no folders. The nearest years on the drive are 2008 and 2014.

The 2008 photograph was flown nine months before the window opens, so its
waterline can sit off the band through real change over that gap, not only
through the day's water level and waves. The captions say so.

## To get the three panels

1. Georeference the 2009 captures in ArcGIS against the 2008 mosaic and
   export GeoTIFFs (worse and undocumented accuracy; covers 2009 only).
2. Download orthophotos for 2010 and 2011. Candidates, not yet checked for
   coverage or resolution: USDA NAIP North Carolina 2010 (~1 m, summer), and
   NOAA NGS post-Irene imagery (late August 2011).

The imagery reader only finds a non-USGS year as
`D:\Hatteras_GIS\Aerial\<year>*\<year>_full_aerial.tif`, so save each
mosaic under that name. Then add an entry to `PHOTO_SOURCES` in
`scripts/input_prep/5-scr/1-observations/mean_shoreline/coastsat_mean_shoreline_on_imagery.py`
(flight date, label, source reference), then:

    python coastsat_mean_shoreline_on_imagery.py --window 2009 2011 --photo-years 2009 2010 2011
