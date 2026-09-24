# 2004 island offsets — version index

`CURRENT` names the build every reader takes; `hatteras_site_config._island_offset_file(2004)` resolves it (env `HAT_OFFSET_VERSION_2004` outranks the file).

## Builds

Written by `build_island_offset.py`, one row per build; each version's own `PROVENANCE.md` has the detail.

| version | built | line | vintage | zero domain | compared with | |
|---|---|---|---|---|---|---|
| `superseded_20260924_pre-metres/v1` | 2026-09-15 | `duneline_2004.geojson` | 2004 | GIS 78 | superseded_20260915_flat | |
| `v1` | 2026-09-24 | same raw as the superseded build (`2004_duneline_offset_raw.csv`) | 2004 | GIS 78 | superseded_20260924_pre-metres/v1 (real domains identical) | CURRENT |

`v1` (2026-09-24) is the build now in `superseded_20260924_pre-metres/v1/`, re-padded with the model's smooth wrap-around when the offset went in as metres; the real domains are identical. Built as `v2`, renumbered `v1` the same day (numbering restarted). See `v1/PROVENANCE.md`.
