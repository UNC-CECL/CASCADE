# SUPERSEDED — 1984-2004 runs made BEFORE the 1996-mosaic topography

Archived 2026-08-27. **Do not use these in analysis.** They are kept for
method comparison only.

## Why they are superseded

Every run in here was started from `1-barrier3d-domains/2004-start/dune-topo/
2009_v5` — the gap-filled 2009+2014 DEM. That is the correct topography for the
**2004-2024** period. It is the wrong one for **1984-2004**, which now starts
from the `1984-start` product: the 2009+2014 DEM with the **1996 ALACE survey
grafted over it**, so the barrier begins the run with its 1996 planform rather
than its 2009 one.

The two products are not close. All 90 domains differ and **65 have a different
interior shape** — GIS 11 is 165 rows on `1984-start` and 157 on `2004-start`.
A 1984 run on `2004-start` is a different island, not a rounding difference.

The period/product pairing is now defined once, in `scripts/hat_topo_version.py`
(`YEAR_PRODUCT`), and consumed by `hatteras_site_config.HATTERAS_PERIODS`.

## What is in here

36 run directories: 30 that were live in `run_index.csv` (calibBE, edgeBE and
the zeroBE groin arms, 12 scenario cells each across three presets), plus the
`calibBE_pregroin_20260824/` set that had already been set aside for a
different reason (pre-groin-fit) and shares this topography era.

`run_index_snapshot.csv` holds the 30 index rows exactly as they stood before
archiving. They have been removed from the live `output/raw_runs/run_index.csv`,
because that index is meant to describe what is currently on disk and in use.

## The .npz model states are GONE

All of these had their ~240 MB model pickles deleted on 2026-08-26 — see
`output/raw_runs/npz_purge_20260826.csv`, reason "1984 run on pre-restructure
shared topography".

Consequence for comparison: you can still read the **shoreline change rate
CSVs, the road management summaries, the run metadata and the shoreline
matrices**. You **cannot** re-run `HAT_relocation_comparison.py` against these,
because it pulls `_road_setback_TS` / `_road_relocated_TS` out of the pickled
Cascade. Comparing relocation timing against this era means re-running them.

## Also stale in here, independently of the topography

- **Road setbacks.** These predate the 2026-08-26 16:54 regeneration of
  `RoadSetback_1984_dunestart.csv`, which was the first one measured against
  `1984-start` row 0. See the sibling `superseded_presetbackfix_20260827/`
  README for the size of that change on the relocation domains.
- **Source/sink presets.** `edgeBE`'s GIS 90 value and the whole `calibBE`
  table were solved against base runs on THIS topography. They have not yet
  been re-solved on `1984-start`, so an edgeBE or calibBE run today is carrying
  numbers fit to the island in this folder.

## Known consumers that will notice these moved

- `HAT_run_all.py` resumes from `output/driver/driver_manifest.jsonl`, which
  still records these as complete. A future matrix re-run will SKIP them rather
  than rebuild them. Migrate or archive that manifest before re-running the
  matrix.
- The groin sweep's drift guard differences a rate curve against the published
  matrix run `HAT_1984_2004_edgeBE_road_bdm_groin`, which now lives in here.
