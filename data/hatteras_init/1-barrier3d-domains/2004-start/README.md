# 2004-start

Barrier3D domains for the **2004-2024** hindcast period.

| | |
|---|---|
| DEM | `0-elevation/2009-2014/` - 2009 USACE, gaps filled from 2014 NOAA Post-Sandy |
| current version | `v5` (see `dune-topo/CURRENT`) |
| picks | `picks/HAT_dune_search_windows_v5.json` |

## This is the former `2009-dune-topo/2009_v5`, moved and not rebuilt

Nothing was re-extracted. v5 was always built from the 2009+2014 DEM, which is
what the 2004 start uses, so the migration was a rename. The arrays are
byte-identical to what every run before 2026-08-25 used.

Two things changed name and are worth knowing:

* the run folder is `v5`, not `2009_v5` - the period folder now carries what
  the `2009` prefix used to. **Files inside keep their original names**
  (`HAT_dune_topo_settings_2009_v5.csv`, and `RUN_MANIFEST.txt` still records
  the old paths). Those are outputs of a run that happened; renaming them would
  rewrite history. A re-run would write them as `v5`.
* the pick file is `HAT_dune_search_windows_v5.json`, renamed from
  `..._2009_v5.json` to match `PICK_SET = RUN_NAME`. **This rename was
  required, not cosmetic**: `save_windows()` writes back after every domain, so
  had the extractor been run against a pick file it could not find, it would
  have started from defaults and overwritten v5's windows - which the manifest
  records as not regenerable.

## Caveat on provenance

Runs published before 2026-08-25 used this topography for the **1984** period
too, because the runner had no per-period concept. Their provenance label now
reads `2004-start`, which is where the arrays live but not what those runs
were. The arrays themselves have not changed.

## Superseded versions

`v3` and `v4` are under `../superseded/`. Both were picked against the
*unfilled* 2009 DEM. To regenerate an old figure:
`topo_dirs("2004-start", override=...)` after moving the version back.
