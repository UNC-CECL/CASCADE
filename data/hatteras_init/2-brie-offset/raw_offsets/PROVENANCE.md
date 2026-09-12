# raw_offsets — which survey is behind each file

One file per dune-line survey, holding the per-transect distance from the
shared offshore baseline (`ORIG_LEN`). Everything in `hindcast_<year>/` is
built from these, and the end-year validation target differences two of them,
so the baseline must be the same layer in every file here.

## Files named for a HINDCAST PERIOD rather than a survey

The loaders resolve a file by the period's **start year**
(`<year>_duneline_offset_raw.csv`), the same convention the road setback files
follow. Where no survey exists in that year, the nearest one is copied under
the period's name and recorded here.

| File | Actual survey | Note |
|---|---|---|
| `1996_duneline_offset_raw.csv` | **1997** | byte-identical copy of `1997_duneline_offset_raw.csv`, added 2026-09-11 for the 1996-2010 period |

Consequences of that one copy, stated rather than left to be discovered:

* the 1996 period starts from the island as surveyed a year later;
* a 1996-to-2010 shoreline change computed from these files spans thirteen
  years of survey while the run spans fourteen.

Re-digitising a 1996-vintage line and re-running the intersection is what would
replace the copy. Delete it in the same commit if that happens, so two files
never claim to be the same survey.

## Coverage

`2017_duneline_offset_raw.csv` is a **Buxton-only clip**, eleven domains, not
an island-wide survey. It cannot stand in for a hindcast start or an end-year
target. `1978` and `1997` carry domains past GIS 90; everything downstream
keeps 1 to 90 and drops the rest.
