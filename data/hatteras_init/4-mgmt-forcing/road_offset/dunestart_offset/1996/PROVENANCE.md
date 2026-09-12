# RoadSetback_1996_dunestart.csv

DERIVED, not measured. Written by
`scripts/input_prep/4-mgmt-forcings/road_offset/1-produce/HAT_road_setback_derived_vintages.py`
on 2026-09-11.

    1996 = RoadSetback_1984_dunestart.csv + the 1989 Pea Island relocation

The 1989 relocation has already happened by 1996 and the 1999 one has not, so
the 1996 road is the 1984 road with exactly one event applied. The
displacements are read from `HATTERAS_ROAD_EVENTS`, which is the same source
the model uses when the event fires, so the two cannot disagree.

Measured against interior row 0 of the **1984-start** extraction, inherited
from the 1984 file. A 1996 run must read that topography product.

No NC-12 line of 1996 vintage exists in the repo. Digitising one off the 1997
imagery and re-measuring is the thing that would replace this file.
