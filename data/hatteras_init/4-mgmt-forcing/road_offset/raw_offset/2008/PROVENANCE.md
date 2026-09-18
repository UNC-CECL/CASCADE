# The 2008 road line, as digitised -- the 2004 period's road

`nc12_2008.geojson` is the source alignment everything else in `road_offset/` is derived from.

**Filed under its imagery year since 2026-09-15.** This folder was `raw_offset/2004/` and the files
`nc12_2004.{csv,geojson}` -- the period they stand in for, not the year they were digitised off.
The rename makes `2008` mean a LINE everywhere in `4-mgmt-forcing` and `2004` mean a PERIOD; which
period reads this line is `scripts/site_layer/hat_topo_version.py:ROAD_LINE_FOR_YEAR`. The masks burned from
it are under `raster/2008/`, and the ArcGIS origin (`\\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\nc12_2004.csv`)
keeps its old name.

## The line behind it is not from 2004

The two NC-12 alignments in this repository were digitised off **1978 and 2008
imagery** and stand in for 1984 and 2004. That is deliberate and recorded in
`hatteras_site_config.py` and in `road_offset/README.md`; it is not a mistake
to be corrected.

What follows from it, and matters when reading any road number:

* about 70% of the two lines is the same vertices, because the later one was
  digitised by editing a copy of the earlier. Most of the "no movement" between
  them is therefore an editing artefact rather than a measurement.
* the 2004 period inherits everything the 2008 imagery gets wrong about
  2004, which is four years of it.
