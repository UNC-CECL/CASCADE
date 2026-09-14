# 1984 road line, as digitised

`nc12_1984.geojson` is the source alignment everything else in `road_offset/` is derived from.

## The line behind it is not from 1984

The two NC-12 alignments in this repository were digitised off **1978 and 2008
imagery** and stand in for 1984 and 2004. That is deliberate and recorded in
`hatteras_site_config.py` and in `road_offset/README.md`; it is not a mistake
to be corrected.

What follows from it, and matters when reading any road number:

* about 70% of the two lines is the same vertices, because the later one was
  digitised by editing a copy of the earlier. Most of the "no movement" between
  them is therefore an editing artefact rather than a measurement.
* the 1984 file inherits everything the 1978 imagery gets wrong about
  1984, which is six years of it.
