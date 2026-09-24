# 1996 island offsets — version index

`CURRENT` names the build every reader takes; `hatteras_site_config._island_offset_file(1996)`
resolves it (env `HAT_OFFSET_VERSION_1996` outranks the file, for a per-run
selection that does not change the shared default). Every build is **derived,
not surveyed**: no 1996 dune line exists, and each is built from a
digitisation of the **1997** line, the nearest island-wide survey
(`hat_topo_version.DUNE_LINE_FOR_YEAR[1996] == 1997`).

**Numbering restarted 2026-09-19** (Hannah): `v1` is the first build from the
re-digitized 1997 line (2026-09-18). The builds from earlier digitizations
are in `superseded_20260919_pre-redigitized/`, under their old numbers.

## Builds

Written by `build_island_offset.py`, one row per build; each version's own `PROVENANCE.md` has the detail.

| version | built | line | vintage | zero domain | compared with | |
|---|---|---|---|---|---|---|
| `v1` | 2026-09-18 (built as v3) | `duneline_1997.geojson` (re-digitized 2026-09-18) | 1997 | GIS 77 | superseded v2 | |
| `v2` | 2026-09-24 | same raw as v1 (`1997_duneline_offset_raw.csv`) | 1997 | GIS 77 | v1 (real domains identical) | CURRENT |

v2 (2026-09-24) is v1 re-padded with the model's smooth wrap-around; the real domains are identical. See `v2/PROVENANCE.md`.

A new build is `build_island_offset.py --duneline <line> --year 1996`; it
takes the next free number (v3), skipping the superseded folder.

## Superseded (`superseded_20260919_pre-redigitized/`)

| old name | built | line | intersection |
|---|---|---|---|
| `v1` | 2026-09-11 | `duneline_1997.geojson` as digitised 2026-09-02 | ArcGIS 1.5 m buffer / 1 m point export |
| `v2` | 2026-09-15 | `duneline_1997_v2.geojson` (local corrections) | shapely (`duneline_to_raw_offsets.py`); 22 domains moved landward, up to 63 m |

Runs made on them are in `output/raw_runs/archive/2026-09-18-pre-redigitized-dunelines/`
and `.../2026-09-19-pre-redigitized-sens-exp/`; their metadata names them
`superseded_20260919_pre-redigitized/v1` and `/v2`.

## The one-year gap, unchanged by versioning

* the 1996 period starts from the island as surveyed a year later;
* a 1996-to-2010 shoreline change computed from the raw files spans thirteen
  years of survey while the run spans fourteen.

## Extended reaches (experiment inputs, not versions)

Written by `build_island_offset.py --geometry`; each reach's own `ext/<geometry>/PROVENANCE.md` has the detail. GIS 1-90 is identical to the version named under *built as*; `CURRENT` is untouched.

| reach | built | line | vintage | built as | | |
|---|---|---|---|---|---|---|
| `ext/n115` | 2026-09-18 | `duneline_1997.geojson` (re-digitized) | 1997 | v1 (then named v3) | polygon join onto `1-barrier3d-domains/domain-geojson/domains_pea_hatteras_120.geojson` | GIS 1 to 115 |
