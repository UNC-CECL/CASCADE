# RoadSetback_2010_dunestart.csv

A COPY of `RoadSetback_2004_dunestart.csv`, written by
`scripts/input_prep/4-mgmt-forcings/road_offset/1-produce/HAT_road_setback_derived_vintages.py`
on 2026-09-11.

Under `dunestart_offset/derived/` since 2026-09-15: the address says this is not a
measurement. The source is `dunestart_offset/measured/2004/RoadSetback_2004_dunestart.csv`,
reached through `hat_topo_version.road_setback_file(2004)`.

The 2010 period reads the same topography product and the same road line as
the 2004 period, and no relocation in the record falls between the two dates --
the next event is the 2022 bridge, which is inside the run rather than at year
zero. So the measurement is identical.

It is written as its own file rather than referenced across periods so that
every period names its own forcing. The two are byte-identical by
construction: re-run the producer rather than editing either one.
