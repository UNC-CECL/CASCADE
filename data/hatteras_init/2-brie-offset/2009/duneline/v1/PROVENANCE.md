# 2010 island offsets (duneline), v1

Built 2026-09-24 by `scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 2010 --version v1 --source duneline --raw-file data/hatteras_init/2-brie-offset/2010/duneline/v1/2009_duneline_offset_raw.csv`, from the SAME raw file as the build it replaces, `../superseded_20260924_pre-metres/v1/` (copied here). Built as `v2` and renumbered `v1` the same day (Hannah: "restart all the years as v1" once the offset went in as metres); the step output below keeps the name it was built under.

## What changed from the superseded build: the buffer domains only

| | `superseded_20260924_pre-metres/v1` | this build |
|---|---|---|
| real domains GIS 1-90 | metres | **identical** (checked value for value) |
| unpadded files | | **identical** |
| 15 + 15 buffer domains | local slope segment + linear bridge, clipped at 0 | the smooth wrap-around the model uses: `cascade_pipeline.hindcast.pad_offset_ring` (cubic Hermite from GIS 90 back round to GIS 1, matched to the island's end slopes) |

Why: BRIE's domain is periodic, so the buffers must carry the shoreline from the last real domain back to the first. In `metres` mode the runner never used the old buffers -- `build_island_offset` threw them away and re-closed the ring with this same curve -- so the old file, and its buffer diagnostic, showed a buffer the model did not see (they differ by up to about 2 km). This build writes the closure itself (Hannah, 2026-09-24, "option (a)"): **the padded file is now exactly what offset_mode `metres` hands Cascade.** Checked: `build_island_offset(<this file>, "metres")` returns it unchanged, and equals `build_island_offset(<superseded file>, "metres")`.

## Units

Metres throughout: the raw file's `ORIG_LEN` is measured in EPSG:3725 (UTM 18N, metres), and nothing in the build converts it. Real span 0.0 - 6194.4 m.

## What this does and does not change for runs

- **metres and detrended runs: nothing.** They receive the same array from this build or the superseded one (checked).
- **asrun (offset / 10) runs: the buffers.** asrun divides the whole file, buffers included, so a run made before 2026-09-24 reproduces only with the build it was made from: `HAT_OFFSET_VERSION_2010=superseded_20260924_pre-metres/v1`.

## Shoreline angles BRIE reads

Largest angle between neighbouring domains: real reach 26.1 deg; buffers 27.8 deg (BRIE goes anti-diffusive past ~42). Both under the ~42 degree limit where `(1.2 sin^2 - cos^2)` changes sign and the shoreline goes anti-diffusive. `Island_Dune_Offsets_2010_buffer_diagnostic.png` draws the padded profile and these angles (caption in `supporting/CAPTIONS.md`).

## CURRENT

`../CURRENT` = `v1` since 2026-09-24. `HAT_OFFSET_VERSION_2010` in the environment outranks it for one run.

## Step output

```
--- Processing 2010 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\duneline\v1\2009_duneline_offset_raw.csv
  90 domains processed.
  Baseline distance = 2004.070 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\duneline\v2\Island_Dune_Offsets_2010_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\duneline\v2\Island_Dune_Offsets_2010_CASCADE_Input.csv

Padding summary (smooth wrap-around, pad_offset_ring):
  Buffer domains per side : 15
  Real span               : 0.0 - 6194.4 m
  Buffer span             : 893.0 - 6253.9 m
  Largest angle, real     : 26.1 deg
  Largest angle, buffer   : 27.8 deg (BRIE goes anti-diffusive past ~42)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\duneline\v2\Island_Dune_Offsets_2010_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2010\duneline\v2\Island_Dune_Offsets_2010_buffer_diagnostic.png
```
