# ARCHIVE — dune-start road offset on the **2009_v3** topography

**Frozen 2026-08-19. Do not regenerate into this folder, and do not feed it to a
model run.** The live outputs are in `../dunestart_offset/`, which is built on
`2009_v4`.

This is the complete `dunestart_offset/` tree exactly as it stood before the v4
dune windows were re-picked — every file, including the figures and the two
audit markdowns.

## Why it exists

The dune search windows were re-picked in `2009_v4` (2026-08-19) with NC-12 drawn
on the picker, which moved 21 of the 90 windows and changed 18 domains' interior
arrays. The road setback is measured from interior row 0, so those setbacks moved
with it. This folder keeps the pre-re-pick answer so the two can be compared.

## Provenance — this is a CONSISTENT v3 set

Everything here was produced against `2009_v3`, both the setbacks and the
downstream audits. That consistency is worth stating because the intermediate
state was *not* consistent: for a short window the produce step had moved to v4
while four downstream scripts still read v3 interiors from a hardcoded path
(since fixed — see `road_offset/hat_topo_version.py`). No file in this folder
comes from that mixed state.

Identifying marks, if you ever need to tell a v3 file from a v4 one:

| | this archive (v3) | live (v4) |
|---|---|---|
| 1984 negative setbacks | **7** — GIS 10, 11, 12, 13, 84, 85, 86 | **6** — GIS 11, 12, 13, 84, 85, 86 |
| `RoadOffset_*_profiles.csv` | no coordinate columns | has `interior_x/y`, `road_x/y` |
| `RoadOffset_dunestart_audit.md` | says `2009_v3` | says `2009_v4` |
| `RoadSetback_audit.csv` | 166 FITS / 9 DROWNS | 167 FITS / 8 DROWNS |

## What changed between this and v4

**12 domains moved in each year**, all by 5–60 m:

| GIS | 1984 v3 → v4 | 2004 v3 → v4 |
|---|---|---|
| 9 | +20 → +40 | +40 → +50 |
| **10** | **−40 → +20** | +10 → +70 |
| 11 | −70 → −60 | +10 → +20 |
| 14 | +20 → +30 | +50 → +55 |
| 31 | +230 → +250 | +230 → +250 |
| 32 | +220 → +230 | +220 → +230 |
| 33 | +250 → +260 | +250 → +260 |
| 62 | +340 → +350 | +340 → +350 |
| 66 | +340 → +370 | +340 → +370 |
| 79 | +510 → +500 | +510 → +500 |
| 85 | −70 → −60 | +30 → +50 |
| 86 | −70 → −60 | +20 → +30 |

The one qualitative change is **GIS 10**, which stopped being negative in 1984
(−40 → +20), dropping the floored count from 7 to 6. It also flipped from DROWNS
to FITS in the 1999 inter-village relocation scenario — though note that verdict
flip only became visible once the hardcoded-v3 bug was fixed, so the DROWNS
verdict recorded in this archive's `RoadSetback_audit.csv` was itself computed on
v3 interiors and is internally correct for this folder.

Unchanged between the two: road span GIS 9–90 (82 domains), median setback
195 m, the three drowning relocations (GIS 78/79/80 moved 20/100/50 m seaward),
D8 still `EXCLUDED_FROM_SPAN`, and 2004 still has zero negatives.

## Orphans — this is now their ONLY copy

These four have no producing script any more. They are dated 2026-08-17 and
predate the `1-produce / 2-audit / 3-figures / 4-compare` reorg; nothing in
`scripts/` writes those names today:

- `HAT_road_method_on_island.png`
- `HAT_road_on_b3d_domains.png`
- `HAT_road_placement_accuracy.csv`
- `HAT_road_placement_accuracy.png`

The nearest current equivalents are `HAT_dunestart_road_on_domains.png` and the
placement statistics inside `4-compare/HAT_method_comparison_figures.py`.

**They were deleted from the live `../dunestart_offset/` on 2026-08-19**, after
each was verified byte-identical to the copy here, because sitting in the live
folder they read as current output when they are two topography versions stale
and unregenerable. This archive is the only copy that remains.

### `figures_1984_v3_leftovers/`

Same story, different folder. `HAT_road_domain_views.py` used to write into
`road_offset/figures/<year>/`; it now writes flat into `road_offset/figures/`.
Three files were stranded under the old layout:

- `road_init_GIS052_1984.png`
- `road_init_overview_1984.png`
- `barrier3d_grid/b3d_grid_GIS052_1984.png`

All three have current v4 equivalents of the same name in `road_offset/figures/`,
so nothing was lost by removing them from the live tree — but they are v3-era
output, so they are kept here. Note these are the ONLY v3-era files from
`figures/` that survive; the rest of that tree was overwritten before any
snapshot was taken (see "Not archived" below).

## Not archived

Only `dunestart_offset/` was captured (plus the three stranded figures above).
The v3-era versions of these were overwritten in place and are **not
recoverable** (untracked in git, no snapshot):

- `road_offset/HAT_method_comparison_on_domains.png`, `HAT_method_vs_actual_road.png`
- `road_offset/figures/` (per-domain views, apart from the three in
  `figures_1984_v3_leftovers/`)
- `road_offset/old_method_offset/HAT_old_method_road_on_domains.png`
- `road_offset/raster/HAT_road_geojson_on_2009_dem.png`

The road masks (`raster/*/masks/`) and the source geojson (`raw_offset/`) are
**inputs**, unaffected by any dune/topo version, and were never regenerated.
