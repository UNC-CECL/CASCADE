# 2026-09-18-pre-redigitized-dunelines

The 24 production matrix runs of the 1996 and 2010 starts, moved here intact
(.npz included) on 2026-09-18, before being re-run on the island offsets
rebuilt from the re-digitized 1997 and 2009 dune lines.

| start | island offset these runs used | the re-run uses |
|---|---|---|
| 1996 | `2-brie-offset/1996/v2` (duneline_1997_v2) | `1996/v3` |
| 2010 | `2-brie-offset/2010/v1` (the old duneline_2009) | `2010/v2` |

Why: Hannah re-digitized the 1997, 2009 and 2023 lines on 2026-09-18. They
moved 32 and 50 domains respectively, by up to about 66 m, every one of them
seaward (`2-brie-offset/dunelines/README.md`). These are the "before" runs for
that change.

Also here:
- `run_index_snapshot.csv`, the index as it stood before the move.
- `driver_manifest.jsonl`, moved out of `output/driver/` (now `output/logs/driver/`) so that
  HAT_run_all.py does not skip the re-run as already done (trap 2 in the
  run-archive notes).

Reproducible from tracked config plus git history: the old lines are in git
before 2026-09-18, and each old offset build keeps its raw file. They were
still moved, not deleted.
