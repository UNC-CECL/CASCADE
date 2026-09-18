# 2004 road setback, dune-start method

`RoadSetback_2004_dunestart.csv` is the model-facing forcing: metres landward of **interior row 0** of the 2004-start extraction, measured on the **2008 NC-12 line** (`raw_offset/2008/`, masks `raster/2008/`).

Under `dunestart_offset/measured/` since 2026-09-15 because it IS a measurement; the `derived/` sibling holds the files built from it.

## The line behind it is not from 2004

The two NC-12 alignments in this repository were digitised off **1978 and 2008
imagery** and stand in for 1984 and 2004. That is deliberate and recorded in
`hatteras_site_config.py` and in `road_offset/README.md`; it is not a mistake
to be corrected.

What follows from it, and matters when reading any road number:

* about 70% of the two lines is the same vertices, because the later one was
  digitised by editing a copy of the earlier. Most of the "no movement" between
  them is therefore an editing artefact rather than a measurement.
* the 2004 file inherits everything the 2008 imagery gets wrong about
  2004, which is four years of it.

## It belongs to one extraction

A setback measured from interior row 0 is only valid against the arrays it was measured on. Spending it on another topography version measures from a row that moved. Resolve the version through `scripts/site_layer/hat_topo_version.py`; never pin one by hand.
