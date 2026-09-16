# 1996 island offsets — version index

Two builds of the 1996 start, side by side. `CURRENT` names the one every
reader takes; `hatteras_site_config._island_offset_file(1996)` resolves it
(env `HAT_OFFSET_VERSION_1996` outranks the file, for a per-run selection that
does not change the shared default). Both are **derived, not surveyed**: no
1996 dune line exists, and each is built from a digitisation of the **1997**
line, the nearest island-wide survey.

```
v1   from duneline_1997.geojson     ArcGIS 1.5 m buffer / 1 m point intersection   2026-09-11
v2   from duneline_1997_v2.geojson  shapely intersection (duneline_to_raw_offsets.py)   2026-09-15   CURRENT
```

| | dune line | raw file | intersection | what changed |
|---|---|---|---|---|
| `v1` | `2-brie-offset/dunelines/duneline_1997.geojson` (digitised 2026-09-02) | `raw_offsets/1997_duneline_offset_raw.csv` | ArcGIS export | — |
| `v2` | `.../duneline_1997_v2.geojson` (local corrections, 2026-09-15) | `raw_offsets/1997_duneline_offset_raw.csv` (named `1997_v2_...` for a few hours; the vintage's current build) | `duneline_to_raw_offsets.py`, validated against the v1 export | 22 of 90 domains, all moved landward, up to 63 m — see `v2/PROVENANCE.md` |

Each version folder holds the three `Island_Dune_Offsets_1996_*` files the
model and the poster read, the buffer diagnostic, a copy of the raw file it
was built from (`1997_duneline_offset_raw.csv`: the GIS export in v1, the
shapely build in v2), and its own PROVENANCE.

The raw file a 1996 build reads is resolved through
`hat_topo_version.DUNE_LINE_FOR_YEAR[1996] == 1997` since the evening of
2026-09-15. Before that, `raw_offsets/1996_duneline_offset_raw.csv` was a
byte-identical copy of the 1997 file; the copy is gone.

## Rebuilding either

```
python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1996 --version v2 \
    --raw-file data/hatteras_init/2-brie-offset/1996/v2/1997_duneline_offset_raw.csv
python scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1996 --version v1 \
    --raw-file data/hatteras_init/2-brie-offset/1996/v1/1997_duneline_offset_raw.csv
```

Without `--raw-file` the script reads the vintage's current raw, which is the
v2 build today. v1 was rebuilt from its own raw on 2026-09-15 and came out
byte-identical. A v3 would come from `build_island_offset.py --duneline
<a new 1997 or a true 1996 line> --year 1996` (a 1996 line needs
`DUNE_LINE_FOR_YEAR[1996]` changed to 1996 first).

## The one-year gap, unchanged by versioning

* the 1996 period starts from the island as surveyed a year later;
* a 1996-to-2010 shoreline change computed from the raw files spans thirteen
  years of survey while the run spans fourteen.

## Runs

Three 1996-2010 runs exist and **all of them read v1** (they predate v2):

```
HAT_1996_2010_zeroBE_road_bdm_nogroin        2026-09-11   calibration
HAT_1996_2010_edgeBE_road_bdm_nogroin        2026-09-12   calibration, and a version-check/v1 arm
HAT_1996_2010_zeroBE_noroad_nobdm_nogroin    2026-09-13   calibration
```

under `output/raw_runs/1996_2010/`. They are not wrong, they are on the
superseded start: a re-run on v2 replaces them, and until then a 1996 figure
should say which offset build it was drawn from. Any 1996 run from here reads
v2 unless `HAT_OFFSET_VERSION_1996` says otherwise.
