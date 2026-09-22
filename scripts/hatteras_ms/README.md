# hatteras_ms — the hindcast, and everything that reads its output

Four kinds of thing live here, and until 2026-09-13 they were one flat list of
31 items sharing a `HAT_` prefix that sorted them without separating them.

```
HAT_hindcast_1984_2024.ipynb   THE RUN. The notebook is authoritative;
HAT_hindcast_1984_2024.py      the .py is its headless mirror, with no
                               features of its own. A change goes in BOTH.
HAT_hindcast_config.py         which run happens, and where the value came
hat_run.yaml                   from: env > yaml > the default in the module
HAT_run_all.py                 the batch driver: the matrix, then the sweep
HAT_hindcast_methods.md        the written method
HINDCAST_PLAN.md               planning notes: the order the notebook builds
                               the run in. Was `HAT_hindcast_plan`, a file
                               with no extension that this line described as
                               a folder (2026-09-22)

tools/        read or repair the run record; none of them run the model
experiments/  one-off studies, each a run-it then plot-it pair
figures/      figures built from finished runs
groin-sweep/  the M and f fit, its own config and worker
old_drafts/   superseded
old_versions/ superseded runners, and the inherited driver that predates them
```

## tools/

| Script | What it answers |
|---|---|
| `HAT_period_input_check.py` | is this period runnable, and what is missing |
| `HAT_index_runs.py` | rebuild `run_index.csv` from the runs on disk |
| `HAT_list_runs.py` | what is on disk, grouped |
| `HAT_run_supersession_report.py` | which runs a later one has superseded |
| `HAT_migrate_run_layout.py` | move runs to the current directory layout |

Start with the period check before a run and `HAT_list_runs.py` after one.

## experiments/

Two live studies, plus one retired.

* **the crest experiment** -- `HAT_run_crest_experiment.py` and its plotter.
* **the relocation set** -- `HAT_relocation_comparison.py` (takes `--period`
  since 2026-09-15: 1984 or 1996, one output root per window),
  `HAT_relocation_period_compare.py` (the 1999 event read across both
  windows, from the per-period tables),
  `HAT_relocation_dune_position_check.py`, `HAT_score_relocation_timing.py`
  and `HAT_score_road_position.py`. `RELOCATION_COMPARISON_RESULTS.md` is
  what it concluded.
* `superseded_20260907/` -- the 1984 seaward row-insert set, which cannot be
  re-run: the topography layers it studied were deleted, its output folders
  are empty, and its arms are not in the run tree. Kept because the four
  scripts are the only record of how it was driven. Its `WHY.md` has the
  evidence.

A driver spawns the runner as a subprocess with `HAT_IGNORE_SETTINGS=1`, so
whatever is sitting in `hat_run.yaml` cannot reach an experiment.

## figures/ — moved

The five figure scripts that were here (`hindcast_final_figure_loess`,
`scenario_grid`, `rerender_run_figures`, `planview_evolution_gif`,
`gis11_relocation_drown_figure`) and their `superseded_20260914/` moved to
`scripts/figure_making/model_output/` on 2026-09-18, so every figure script is
in one tree.

## Paths: search upward, do not count

Every script here finds the project root by walking up until it sees a marker,
never by counting parent directories:

```python
REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
```

Six files already did this; the other fourteen counted depth and were
converted on 2026-09-13, before the move, so that the move itself could not
break them. Keep it that way — a counted depth is correct only until the file
is filed somewhere better.

The runner stays at the top level because the config module, the settings yaml
and the batch driver all reach it as a sibling. A script that needs it names
`REPO / "scripts" / "hatteras_ms"` rather than its own folder.
