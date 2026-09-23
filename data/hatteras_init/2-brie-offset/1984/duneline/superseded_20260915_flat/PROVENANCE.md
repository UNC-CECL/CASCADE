# 1984 island offsets — the flat build, superseded 2026-09-15

This was `2-brie-offset/1984/duneline/` itself until 2026-09-15: built from the ArcGIS export of the 1984 line (now `raw_offsets/superseded_20260915_gis_exports/`). Replaced by `../v1/`, built by `build_island_offset.py` from the same geojson through the shapely intersection; `../CURRENT` names v1. The model-frame difference is under 0.6 m at every domain (`../v1/offset_1984_superseded_20260915_flat_vs_v1.csv`). Every run before that date read this build; its metadata says `island_offset_version` is unset, which is how to tell.

---

The note as it stood:

# 1984 island offsets

Built from `raw_offsets/1984_duneline_offset_raw.csv`, which is a **genuine
1984 dune-line survey**. No stand-in.

Produced by `scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1984`.

The padded file is zeroed on its own most seaward domain, so it cannot be
differenced against another year's padded file. Difference the raw files
instead: they share a fixed offshore datum.
