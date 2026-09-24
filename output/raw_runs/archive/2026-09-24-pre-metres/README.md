# 2026-09-24-pre-metres

Every production matrix run (88) and every sensitivity-sweep cell (103), moved
here intact (`.npz` included, 21 GB) on 2026-09-24, when the island offset
started going into the model in metres. Hannah: "all of these runs before the
meter offset correction [should] be moved to the archive since they aren't
active in my analysis anymore".

## What these runs used

`offset_mode: asrun`: the island offset divided by ten. BRIE's shoreline and
the offset file are both metres, so every run here carries one tenth of the
measured planform. See
`experiments/2026-09-24-island-offset-scale-wave-tuning/README.md` for the
study that led to the switch, and `2-brie-offset/<start>/<source>/v1/PROVENANCE.md`
for the builds the runner reads now.

| folder | runs | periods | what |
|---|---|---|---|
| `matrix/` | 88 | 1984-2004, 1996-2010, 2004-2024, 2010-2024 | the production scenario matrix, calibBE / edgeBE / zeroBE |
| `sensitivity/` | 103 | 1984, 1996, 2004 starts | the rset, waveHs, waveTp, waveahf and waveasym sweeps |

Each run records the offset build it used in `island_offset_version`:

- `superseded_20260924_pre-metres/v1` (relabelled 2026-09-24): the 1996,
  2004 and 2010 builds with the old slope-and-bridge buffers, now in
  `2-brie-offset/<start>/duneline/superseded_20260924_pre-metres/v1/`.
- empty: the flat pre-versioning builds of the 1984 and 2004 starts, now
  `2-brie-offset/<start>/duneline/superseded_20260915_flat/`.

Many 1984-start runs were already flagged in `../../SUPERSEDED_CANDIDATES.md`
as built on dune-topo v1 while v2 is CURRENT.

## Reproducing one

```
HAT_OFFSET_MODE=asrun
HAT_OFFSET_VERSION_<start>=superseded_20260924_pre-metres/v1   # or superseded_20260915_flat
```

plus the run's own settings from its `_run_metadata.json`. asrun divides the
whole file by ten, buffers included, so only the build a run was made from
reproduces it.

## Also here

- `run_index_snapshot.csv`: the index as it stood before the move.
- `manifests/`: moved out of the live tree so a re-run is not skipped or
  mixed with these cells:
  - `driver_manifest.jsonl` (was `output/logs/driver/`). `HAT_run_all.py`
    decides a matrix job is done by `stage|period|preset|scenario|groin|reloc|be_digest`,
    with no offset mode, so left in place it would have skipped every metres
    re-run as already done.
  - `sensitivity_1984.jsonl`, `sensitivity_1996.jsonl`, `sensitivity_2004.jsonl`,
    `wave_sensitivity_1984.jsonl` (were `output/calibration/sensitivity/`),
    the sweep manifests `plot_sensitivity.py` reads.

## A departure from the usual rule

`../../README.md` says superseded runs are archived at the moment of the
re-run that replaces them. These were archived before any metres re-run
exists, at Hannah's request, because they are no longer part of the
analysis. Until the matrix and sweeps are re-run in metres, anything that reads
production runs from `matrix/` or `sensitivity/` finds none:

- the figure scripts that draw the hindcast and model-vs-observed figures
  (the figures already made are unchanged; redrawing them needs the re-runs);
- `be_edge_domain_solve.py` (edge solves read matrix runs);
- `plot_sensitivity.py`;
- the groin sweep's validation lookup (it validates against a matrix run).

The calibration these runs carry (edgeBE end rates, groin M and f, Hs 2.5) was
fitted under asrun and has to be re-solved in metres before the re-run.
