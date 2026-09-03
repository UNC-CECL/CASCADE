# buffer/ - the padding domains

Three arrays that pad the alongshore ends of the model grid. They are **not**
produced by any script in the current chain and carry no product and no
version: the same three files serve every period. See `../LINEAGE.md`.

## What is here

| file | shape | units | sentinel |
|---|---|---|---|
| `domain_buffer.npy` | (50, 200) | m NAVD88, raw DEM frame (alongshore, cross-shore) | `-10` nodata |
| `sample_1_topography.npy` | (200, 50) | dam, processed interior | `-0.3` dam (= -3.0 m MHW-relative) |
| `sample_1_dune.npy` | (50,) | dam above berm | - |

`domain_buffer.npy` is a raw export - 1 m-era orientation, `-10` nodata, the
`nodata_to_value=-10` convention of the retired ArcGIS exporters. The other two
are already through the extractor and are in Barrier3D's native decameters.

`sample_1_dune.npy` is **(50,)**, one row. Barrier3D runs `DuneWidth = 2` and
builds row 1 by copying row 0 - a packed `(100,)` file is silently truncated,
so one row is correct here, not a bug.

## Provenance

All three were copied on 2025-11-07 from **GIS domain 111** of the pre-reorg
tree, by a one-shot script (`buffers/buffer_creation.py`) that was deleted on
2026-08-26:

| source, in the tree as it stood then | copied to |
|---|---|
| `hatteras_init/elevations/2009_pea_hatteras/domain_111_resampled.npy` | `domain_buffer.npy` |
| `hatteras_init/topography/2009/domain_111_resampled_topography_2009.npy` | `sample_1_topography.npy` |
| `hatteras_init/dunes/2009/domain_111_resampled_dune_2009.npy` | `sample_1_dune.npy` |

None of those three source directories exist any more - `elevations/`,
`topography/` and `dunes/` were folded into the numbered tree - so the copy
is **not reproducible from the current repo**.

**These three files are also not in git.** `.gitignore:164` excludes all of
`data/hatteras_init/1-barrier3d-domains/`, so what is in this folder is the
*only* instance anywhere: not recoverable from history, not regenerable by
any script. Back them up. Do not delete them expecting to rebuild them.

The deleted script is recoverable at `7cd5af0`:

    git show 7cd5af0:scripts/input_prep/1-barrier3d-domains/buffers/buffer_creation.py

It could not be run in any case: its `DESTINATION_DIR` was a root-absolute
`/data/hatteras_init/buffer`, which on Windows resolves to `C:\data\...`, not
into this repo.

## Why domain 111

The model reach is GIS domains **1-90** (`num_real_domains=90` in
`scripts/hatteras_site_config.py`), padded by 15 buffer domains a side. The
clipped-DEM tree was built out to **1-136** - `domains.geojson` covers more of
the island than the model uses. Domain 111 is one of those extra boxes: a real
surveyed 2000 x 500 m section (elevations up to 7.1 m NAVD88), north of the
modelled reach.

So the padding is a clone of genuine island, from outside the reach being
modelled, rather than a synthetic profile or a repeat of an end domain. The
same clip is still on disk at `domain-clips-1m/domain_111/`, which is what
makes the choice checkable even though the copy itself is not reproducible.

The out-of-reach clips **91-136 were purged on 2026-08-26** as unread by any
script - except `domain_111`, held back precisely so the check above still
works. See `../archive_purge_20260826.csv`.
