# 1996 island offsets

**Derived, not surveyed.** Built from the **1997** dune line, the nearest
island-wide survey, because no 1996 line exists.

Produced by `scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1996`
reading `raw_offsets/1996_duneline_offset_raw.csv`, which is itself a
byte-identical copy of the 1997 file. See `raw_offsets/PROVENANCE.md`.

What follows from the one-year gap:

* the 1996 period starts from the island as surveyed a year later;
* a 1996-to-2010 shoreline change computed from these files spans thirteen
  years of survey while the run spans fourteen.

Sanity check that held: the 1996 profile sits between the 1984 and 2004
profiles, which is what chronology requires.

Digitising a 1996-vintage line and re-running the intersection is what would
replace this. Delete the copy in the same commit if that happens, so two files
never claim to be the same survey.
