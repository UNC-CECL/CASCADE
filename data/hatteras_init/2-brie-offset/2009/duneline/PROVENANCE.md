# 2010 island offsets — version index

`CURRENT` names the build every reader takes; `hatteras_site_config._island_offset_file(2010)` resolves it (env `HAT_OFFSET_VERSION_2010` outranks the file).

**Numbering restarted 2026-09-19** (Hannah): `v1` is the first build from the
re-digitized 2009 line (2026-09-18). The builds from the earlier 2009 line
are in `superseded_20260919_pre-redigitized/`, under their old numbers.

## Builds

Written by `build_island_offset.py`, one row per build; each version's own `PROVENANCE.md` has the detail.

| version | built | line | vintage | zero domain | compared with | |
|---|---|---|---|---|---|---|
| `superseded_20260924_pre-metres/v1` | 2026-09-18 | `duneline_2009.geojson` (re-digitized 2026-09-18) | 2009 | GIS 76 | superseded v1 | |
| `v1` | 2026-09-24 | same raw as the superseded build (`2009_duneline_offset_raw.csv`) | 2009 | GIS 76 | superseded_20260924_pre-metres/v1 (real domains identical) | CURRENT |

`v1` (2026-09-24) is the build now in `superseded_20260924_pre-metres/v1/`, re-padded with the model's smooth wrap-around when the offset went in as metres; the real domains are identical. Built as `v2`, renumbered `v1` the same day (numbering restarted). See `v1/PROVENANCE.md`.

## Superseded (`superseded_20260919_pre-redigitized/`)

| old name | built | line |
|---|---|---|
| `v1` | 2026-09-15 | `duneline_2009.geojson` as digitized before 2026-09-18 |

Runs made on it are in `output/raw_runs/archive/2026-09-18-pre-redigitized-dunelines/`
and `.../2026-09-19-pre-redigitized-sens-exp/`; their metadata names it
`superseded_20260919_pre-redigitized/v1`.
