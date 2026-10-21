# Folder renames

Old folder names still appear in run metadata (`*_run_metadata.json` / `.txt`, `run_tag`), logs and per-run captions written before the move. Those files were left as written. `run_index.csv`, the scripts, the READMEs and the NOTEs were updated to the new paths.

## 2026-10-20: the two largest topics split, loose checks filed

| old | new |
|---|---|
| `storms-overwash-and-dunes/2026-10-13-dune-*`, `2026-10-19-dune-*`, `2026-10-20-dune-findings-figures`, `2026-10-20-stable-dune-growth-rate` (13 studies) | `dunes/<same study>` |
| `storms-overwash-and-dunes/2026-10-13-overwash-scoring-tolerance`, `2026-10-13-soundside-overwash`, `2026-10-19-false-alarm-locations`, `2026-10-19-isabel-pattern-scores`, `2026-10-20-storm-notebook-rule` | `storms-overwash/<same study>` |
| `groin/2026-10-05-blocking-fit-dem-to-dem`, `2026-10-08-conserving-groin`, `2026-10-08-drift-groin`, `2026-10-08-schedule-refit` | `groin-blocking/<same study>` |
| `groin/<every study from 2026-10-09 on>` (18 studies) | `groin-wave-driven/<same study>` |
| `2026-10-15-groin-guard-check`, `2026-10-20-groin-cleanup-check` (at the top level) | `code-checks/<same study>` |

The split follows the write-up tree: `hard-structures/groin/4-hindcast/2-blocking-1996-2025/` and `3-wave-driven-1967-2025/`. The topics are siblings, not `groin/blocking/...`, because a run tag can have at most four parts (`run_registry.check_tag`) and the groin tags already use four.

The offloaded studies on `D:\archive\CASCADE_offload\output\raw_runs\experiments\` keep their old topics (`groin/`, `storms-overwash-and-dunes/`), and scripts for those studies still name them. `storms-overwash/` is a new name. The older `storms-and-overwash/` was renamed `storms-overwash-and-dunes/` on 2026-10-05.

Renames before 2026-10-07 are in `chains/RENAMES.md` on D:.
