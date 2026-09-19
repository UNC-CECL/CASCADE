# 1996 island offsets, v1

**Superseded by `../v2/` on 2026-09-15** (`../CURRENT` says `v2`). Kept
whole because it is reproducible from the repo and is the reference the v2
change is measured against. Moved from `1996/` into `1996/v1/` that day; the
files are unchanged.

**Derived, not surveyed.** Built from the **1997** dune line
(`2-brie-offset/dunelines/duneline_1997.geojson`, digitised
2026-09-02), the nearest island-wide survey, because no 1996 line exists.

Produced 2026-09-11 by `island_offset_hybrid.py --year 1996` reading
`raw_offsets/1996_duneline_offset_raw.csv`, which was then a byte-identical
copy of the ArcGIS export `1997_duneline_offset_raw.csv`. That copy now holds
the v2 survey; to rebuild THIS version, name the v1 export explicitly:

```
python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1996 --version v1 \
    --raw-file data/hatteras_init/2-brie-offset/raw_offsets/1997_duneline_offset_raw.csv
```

Done 2026-09-15: byte-identical to the files here.

What follows from the one-year gap:

* the 1996 period starts from the island as surveyed a year later;
* a 1996-to-2010 shoreline change computed from these files spans thirteen
  years of survey while the run spans fourteen.

Sanity check that held: the 1996 profile sits between the 1984 and 2004
profiles, which is what chronology requires.

The three 1996-2010 runs made 2026-09-11 to 09-13 (`HAT_1996_2010_zeroBE_road_bdm_nogroin`,
`HAT_1996_2010_edgeBE_road_bdm_nogroin`, `HAT_1996_2010_zeroBE_noroad_nobdm_nogroin`)
read THIS version; see `../PROVENANCE.md`.
