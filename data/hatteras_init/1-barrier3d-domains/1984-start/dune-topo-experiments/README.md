# dune-topo-experiments

**Empty of topography as of 2026-09-03.** All five remaining variants were
deleted; only the two comparison figures and this note are left.

## What was here, and where it went

| variant | arm | removed because |
|---|---|---|
| `v1_blocks_measured_duneline` | `blocksduneline` | the N-SOURCE comparison, settled by the 1997 dune line |
| `v1_blocks_measured_dsas` | `blocksdsas` | same, and the `dsas` source itself was removed |
| `v1_blocks_minimum` | `blocksminimum` | same |
| `v1_pea1989_keepcrest` | `pea1989keep` | the crest-shaving pair |
| `v1_pea1989_lowercrest` | `pea1989lower` | the crest-shaving pair |

Three earlier variants (`v1_D85_translate`, `v1_blocks_none_dsas`,
`v1_blocks_none_duneline`) went the same day, having never been referenced.

**The run outputs were NOT touched.** `output/raw_runs/blocksduneline`,
`blocksdsas`, `blocksminimum`, `pea1989keep`, `pea1989lower` and their
`*noreloc` twins are all still there, and
`HAT_plot_crest_experiment.py`, `HAT_score_relocation_timing.py` and
`HAT_score_road_position.py` read run outputs rather than topography — so every
comparison built on these arms still works. What is gone is the ability to
**re-run** them from their inputs.

Sizes and reasons: `../../archive_purge_20260903.csv`. Lineage: `../../LINEAGE.md`.

## The two figures

| file | what |
|---|---|
| `HAT_pea1989_crest_compare.png` | keep-crest vs shave-crest at GIS 84–86 |
| `HAT_seaward_insert_compare.png` | island width and profile under the insert variants |

Both are the only surviving picture of comparisons whose topography is gone.

## Why this folder existed outside `dune-topo/`

`hat_topo_version.versions()` lists the directories under `dune-topo/`, so
leaving thirteen variants there made every resolution message unreadable and
left the "only one version present" fallback permanently dead. If variants are
ever built again, put them here for the same reason.
