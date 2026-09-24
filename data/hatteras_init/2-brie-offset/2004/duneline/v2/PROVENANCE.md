# 2004 island offsets (duneline), v2

Built 2026-09-24 by `scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py --year 2004 --version v2 --source duneline --raw-file data/hatteras_init/2-brie-offset/2004/duneline/v1/2004_duneline_offset_raw.csv`, from the SAME raw file as `../v1/` (copied here). A fix to how the file is padded, not new source data, so the count continues (v1 -> v2).

## What changed from v1: the buffer domains only

| | v1 | v2 |
|---|---|---|
| real domains GIS 1-90 | metres | **identical** (checked value for value) |
| unpadded files | | **identical** |
| 15 + 15 buffer domains | local slope segment + linear bridge, clipped at 0 | the smooth wrap-around the model uses: `cascade_pipeline.hindcast.pad_offset_ring` (cubic Hermite from GIS 90 back round to GIS 1, matched to the island's end slopes) |

Why: BRIE's domain is periodic, so the buffers must carry the shoreline from the last real domain back to the first. In `metres` mode the runner never used the v1 buffers -- `build_island_offset` threw them away and re-closed the ring with this same curve -- so the v1 file, and its buffer diagnostic, showed a buffer the model did not see (they differ by up to about 2 km). v2 writes the closure itself (Hannah, 2026-09-24, "option (a)"): **the padded file is now exactly what offset_mode `metres` hands Cascade.** Checked: `build_island_offset(v2, "metres")` returns the v2 file unchanged, and equals `build_island_offset(v1, "metres")`.

## Units

Metres throughout: the raw file's `ORIG_LEN` is measured in EPSG:3725 (UTM 18N, metres), and nothing in the build converts it. Real span 0.0 - 6163.8 m.

## What this does and does not change for runs

- **metres and detrended runs: nothing.** They receive the same array from v1 or v2 (checked).
- **asrun (offset / 10) runs: the buffers.** asrun divides the whole file, buffers included, so a run made before 2026-09-24 reproduces only with the build it was made from: `HAT_OFFSET_VERSION_2004=v1`.

## Shoreline angles BRIE reads

Largest angle between neighbouring domains: real reach 26.7 deg; buffers 27.7 deg (BRIE goes anti-diffusive past ~42). Both under the ~42 degree limit where `(1.2 sin^2 - cos^2)` changes sign and the shoreline goes anti-diffusive. `Island_Dune_Offsets_2004_buffer_diagnostic.png` draws the padded profile and these angles (caption in `supporting/CAPTIONS.md`).

## CURRENT

`../CURRENT` = `v2` since 2026-09-24. `HAT_OFFSET_VERSION_2004` in the environment outranks it for one run.

## Step output

```
--- Processing 2004 ---
Input file: C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\duneline\v1\2004_duneline_offset_raw.csv
  90 domains processed.
  Baseline distance = 1990.853 m (min mean).

Unpadded offsets saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\duneline\v2\Island_Dune_Offsets_2004_CASCADE_Input_unpadded.csv
Unpadded CASCADE-format file saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\duneline\v2\Island_Dune_Offsets_2004_CASCADE_Input.csv

Padding summary (smooth wrap-around, pad_offset_ring):
  Buffer domains per side : 15
  Real span               : 0.0 - 6163.8 m
  Buffer span             : 903.6 - 6228.0 m
  Largest angle, real     : 26.7 deg
  Largest angle, buffer   : 27.7 deg (BRIE goes anti-diffusive past ~42)

SUCCESS: Padded CASCADE input saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\duneline\v2\Island_Dune_Offsets_2004_PADDED_120.csv

Diagnostic figure saved to:
  C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\2-brie-offset\2004\duneline\v2\Island_Dune_Offsets_2004_buffer_diagnostic.png
```
