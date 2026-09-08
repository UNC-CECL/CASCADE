# `v3` — the 1984 footprint behind the road, copy fill

Built 2026-09-08 08:36 by `HAT_build_footprint_version.py` from `v2` and `row-insert-scope/footprint_1984_by_domain.csv`. **Nothing is re-measured here**: the footprint (`HAT_footprint_1984.py`) and the fill (`HAT_fill_copy_scope.py`) carry the argument; this folder applies them.

| | |
|---|---|
| what | `v2` + the symmetric 1984 footprint placed behind NC-12, filled by copying the N rows that follow the insert point |
| changed | 52 of 90 domains: +73 rows at 29, −47 rows at 23 |
| dune array | unchanged |
| setback CSV | the 1984 setbacks (`setback_new_m`), no floor: 82 road domains change, GIS 85/86 go from 0 to 44/12 m; the road sits on measured cells and the block is directly behind it |
| audit | `HAT_footprint_audit.csv` |

Verified on write: unchanged domains byte-identical to `v2`; changed domains have rows_before + N rows, identical rows before the insert point, and a block equal to the rows that follow it.

To run on it for one run: `HAT_TOPO_VERSION_1984_START=v3` AND copy this folder's `RoadSetback_1984_dunestart.csv` over the forcing-tree one for the run (restore after) - `hatteras_site_config.py` hardcodes that path. `CURRENT` is not changed by building.
