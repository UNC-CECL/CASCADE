# ==============================================================================
# HAT_road_setback_derived_vintages.py
#
# A road setback file for each of the two hindcast periods added 2026-09-11,
# 1996-2010 and 2010-2024.
#
# NEITHER IS A MEASUREMENT. No NC-12 line was digitised for either vintage --
# the repo holds two road lines, exported for 1978 and 2008 and labelled 1984
# and 2004 (see "Road line vintages" in hatteras_site_config.py). So each new
# file is DERIVED from a measured one, and the derivation lives here rather
# than in a hand-edited CSV, because a hand-edit leaves no record of what was
# added to what.
#
#   1996  =  RoadSetback_1984_dunestart.csv  +  the 1989 Pea Island relocation
#
#       The 1989 relocation (GIS 84-87) has ALREADY HAPPENED by 1996 and the
#       1999 one has not, so the 1996 road is the 1984 road with exactly one
#       event applied. The displacements are not retyped here: they are read
#       from HATTERAS_ROAD_EVENTS, so this file and the event the model fires
#       in a 1984-start run cannot disagree about how far the road moved.
#
#       WHAT THIS INHERITS. The 1984 setbacks are measured against interior
#       row 0 of the 1984-start extraction, so these are too, and the 1996
#       period must read that same product. Everything the 1978-vintage line
#       gets wrong about 1984 it also gets wrong about 1996, eighteen years
#       later rather than six.
#
#   2010  =  RoadSetback_2004_dunestart.csv, unchanged
#
#       The 2010 period reads the same topography product as the 2004 period
#       and the same road line, and no relocation in the record falls between
#       the two dates -- the next event is the 2022 bridge, which is inside the
#       run rather than at year zero. So the measurement is identical and this
#       is a copy under the period's own name, written rather than referenced
#       so every period names its own file (Hannah, 2026-09-11).
#
# OUTPUTS
#   dunestart_offset/1996/RoadSetback_1996_dunestart.csv
#   dunestart_offset/2010/RoadSetback_2010_dunestart.csv
#   dunestart_offset/<year>/PROVENANCE.md   what was derived from what
#
# Nothing existing is read for writing, and both outputs are new files.
#
# Author: Hannah A. Henry, UNC CECL
# ==============================================================================

from __future__ import annotations

from pathlib import Path
import sys as _sys

_HERE = Path(__file__).resolve()
_sys.path.insert(0, str(_HERE.parents[4]))          # scripts/
_sys.path.insert(0, str(_HERE.parent))              # this producer folder

from hatteras_site_config import HATTERAS_ROAD_EVENTS   # noqa: E402
from HAT_road_offset_from_dune_start import (           # noqa: E402
    read_two_row_csv, write_two_row_csv, OUT_ROOT)

RELOCATION_YEAR = 1989      # the one event between the 1984 line and 1996
SOURCE_YEAR = {1996: 1984, 2010: 2004}


def setback_path(year: int) -> Path:
    return OUT_ROOT / str(year) / "RoadSetback_{0}_dunestart.csv".format(year)


def event_displacements(year: int) -> dict:
    """The displacement dict of the relocation event dated `year`.

    Read from HATTERAS_ROAD_EVENTS rather than from the measurement CSV
    directly, so this applies the SAME rounded, whole-cell displacement the
    model applies when the event fires.
    """
    events = [e for e in HATTERAS_ROAD_EVENTS
              if getattr(e, "year", None) == year
              and getattr(e, "displacement_m", None)]
    if len(events) != 1:
        raise ValueError(
            "expected exactly one relocation event dated {0}, found {1}"
            .format(year, len(events)))
    return dict(events[0].displacement_m)


def build_1996() -> dict:
    """The 1984 setbacks with the 1989 relocation applied."""
    base = read_two_row_csv(setback_path(SOURCE_YEAR[1996]))
    if not base:
        raise FileNotFoundError(setback_path(SOURCE_YEAR[1996]))

    moved = event_displacements(RELOCATION_YEAR)
    absent = sorted(gis for gis in moved if gis not in base)
    if absent:
        raise ValueError(
            "the 1989 event moves GIS {0}, which the 1984 setback file does "
            "not carry -- adding a displacement to a setback that was never "
            "measured would invent a road position".format(absent))

    derived = dict(base)
    for gis, displacement_m in moved.items():
        derived[gis] = base[gis] + displacement_m

    print("\n1996 = 1984 + the 1989 Pea Island relocation")
    print("  {0:>5}  {1:>10}  {2:>12}  {3:>10}".format(
        "GIS", "1984 m", "moved m", "1996 m"))
    for gis in sorted(moved):
        print("  {0:>5}  {1:>10.1f}  {2:>+12.1f}  {3:>10.1f}".format(
            gis, base[gis], moved[gis], derived[gis]))
    print("  {0} of {1} domains changed".format(len(moved), len(base)))
    return derived


def build_2010() -> dict:
    """The 2004 setbacks, unchanged."""
    base = read_two_row_csv(setback_path(SOURCE_YEAR[2010]))
    if not base:
        raise FileNotFoundError(setback_path(SOURCE_YEAR[2010]))
    print("\n2010 = 2004, unchanged ({0} domains)".format(len(base)))
    return dict(base)


PROVENANCE = {
    1996: """# RoadSetback_1996_dunestart.csv

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
""",
    2010: """# RoadSetback_2010_dunestart.csv

A COPY of `RoadSetback_2004_dunestart.csv`, written by
`scripts/input_prep/4-mgmt-forcings/road_offset/1-produce/HAT_road_setback_derived_vintages.py`
on 2026-09-11.

The 2010 period reads the same topography product and the same road line as
the 2004 period, and no relocation in the record falls between the two dates --
the next event is the 2022 bridge, which is inside the run rather than at year
zero. So the measurement is identical.

It is written as its own file rather than referenced across periods so that
every period names its own forcing. The two are byte-identical by
construction: re-run the producer rather than editing either one.
""",
}


def main() -> None:
    built = {1996: build_1996(), 2010: build_2010()}
    for year, values in built.items():
        path = setback_path(year)
        if path.exists():
            raise FileExistsError(
                "{0} already exists. Delete it deliberately if you mean to "
                "rebuild -- a silent overwrite of a forcing file is how a run "
                "ends up on numbers nobody chose.".format(path))
        write_two_row_csv(path, values)
        (path.parent / "PROVENANCE.md").write_text(
            PROVENANCE[year], encoding="utf-8")
        print("\nwrote {0}".format(path))
        print("wrote {0}".format(path.parent / "PROVENANCE.md"))


if __name__ == "__main__":
    main()
