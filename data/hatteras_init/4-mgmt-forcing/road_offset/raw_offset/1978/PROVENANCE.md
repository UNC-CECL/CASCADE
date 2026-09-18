# The 1978 road line, as digitised -- the 1984 period's road

`nc12_1978.geojson` is the source alignment everything else in `road_offset/` is derived from.

**Filed under its imagery year since 2026-09-15.** This folder was `raw_offset/1984/` and the files
`nc12_1984.{csv,geojson}` -- the period they stand in for, not the year they were digitised off.
The rename makes `1978` mean a LINE everywhere in `4-mgmt-forcing` and `1984` mean a PERIOD; which
period reads this line is `scripts/site_layer/hat_topo_version.py:ROAD_LINE_FOR_YEAR`. The masks burned from
it are under `raster/1978/`, and the ArcGIS origin (`\\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\nc12_1984.csv`)
keeps its old name.

## The line behind it is not from 1984

The two NC-12 alignments in this repository were digitised off **1978 and 2008
imagery** and stand in for 1984 and 2004. That is deliberate and recorded in
`hatteras_site_config.py` and in `road_offset/README.md`; it is not a mistake
to be corrected.

What follows from it, and matters when reading any road number:

* about 70% of the two lines is the same vertices, because the later one was
  digitised by editing a copy of the earlier. Most of the "no movement" between
  them is therefore an editing artefact rather than a measurement.
* the 1984 period inherits everything the 1978 imagery gets wrong about
  1984, which is six years of it.
