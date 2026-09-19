# NC-12 relocation comparison — 1984-2004

One directory per source/sink preset, written by
`scripts/hatteras_ms/HAT_relocation_comparison.py --preset <name>`.
`dune_position_check/` holds the dune-line-vs-road figures from
`HAT_relocation_dune_position_check.py`.

## report.txt

Every run now writes `report.txt` — the console output of the run that produced
the CSVs and GIFs beside it, with a provenance header naming both arms, their
topography product/version, and their git commit.

Two guarantees, both added 2026-08-27:

- **It is regenerated on every run.** The report is a tee of the run's own
  stdout, so it cannot disagree with the CSVs next to it.
- **A crashed run leaves no report rather than the previous one.** Any existing
  `report.txt` is deleted before the work starts; a run that dies partway
  writes a partial report ending in a `RUN FAILED` banner.

## WARNING: the reports in `superseded_*/` are STALE

The report-writing step was absent from the script between roughly 2026-08-22
and 2026-08-27, so a re-run on 2026-08-25 rewrote every CSV and GIF in those
folders and left the 2026-08-22 report untouched:

| folder | report.txt | its CSVs |
|---|---|---|
| `superseded_calibBE_20260826` | 2026-08-22 01:51 | 2026-08-25 10:16 |
| `superseded_edgeBE_20260826`  | 2026-08-22 01:19 | 2026-08-25 10:18 |
| `superseded_zeroBE_20260826`  | 2026-08-22 01:19 | 2026-08-25 10:20 |
| `superseded_edgeBE_20260821`  | (none)           | 2026-08-21 23:27 |

**In those three folders the CSVs are authoritative and the report is not.**
The reports describe an earlier pair of runs — they name the pre-restructure
flat run paths (`raw_runs/1984_2004/HAT_...` rather than
`raw_runs/1984_2004/<preset>/HAT_...`), which is the giveaway.

All of these predate the 1996-mosaic topography and the road-setback
re-measurement, so they are kept for method comparison only. See
`output/raw_runs/1984_2004/superseded_*/README.md` for what makes the
underlying runs superseded, and note that the pre-mosaic runs no longer have
the `.npz` model states these comparisons need — they cannot be regenerated.
