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

## 2026-10-20, later: the end solves split by window

| old | new |
|---|---|
| `end-domain-boundaries/<study>_1996_2009` (5 studies) | `end-domain-boundaries-1996-2009/<same study>` |
| `end-domain-boundaries/<study>_2009_2025` (5 studies) | `end-domain-boundaries-2009-2025/<same study>` |
| `end-domain-boundaries/<study>_1996_2025` (3 studies) | `end-domain-boundaries-1996-2025/<same study>` |

`HAT_end_solve_net_change.py` picks the topic from the solve's window. The topic names use hyphens because the run registry reads any `<4 digits>_<4 digits>` folder as a period folder. The offloaded solves on D: keep `end-domain-boundaries/`.

Renames before 2026-10-07 are in `chains/RENAMES.md` on D:.

## 2026-10-21: the dune growth rate fits moved to data/

The two studies that set the runner's `DUNE_GROWTH_RATES` are model inputs, not experiments. Their write-up,
tables and figures moved; the runs and logs stay in `dunes/` because the run registry files every run under
`output/raw_runs/`, so the run tags are unchanged.

| old | new |
|---|---|
| `dunes/2026-10-20-stable-dune-growth-rate-topo-v3/{NOTE.md, tables/, figures/}` | `data/hatteras_init/1-barrier3d-domains/1996-start/dune-growth-rate/{README.md, tables/, figures/}` |
| `dunes/2026-10-20-test-period-growth-rate/tables/` | `data/hatteras_init/1-barrier3d-domains/2009-start/dune-growth-rate/tables/` (README.md new) |
