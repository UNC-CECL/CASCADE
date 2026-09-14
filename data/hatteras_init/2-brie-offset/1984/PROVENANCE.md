# 1984 island offsets

Built from `raw_offsets/1984_duneline_offset_raw.csv`, which is a **genuine
1984 dune-line survey**. No stand-in.

Produced by `scripts/input_prep/2-brie-offset/island_offset_hybrid.py --year 1984`.

The padded file is zeroed on its own most seaward domain, so it cannot be
differenced against another year's padded file. Difference the raw files
instead: they share a fixed offshore datum.
