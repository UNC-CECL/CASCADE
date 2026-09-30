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

## Which Barrier3D the run uses

Barrier3D is a separate repository, installed editable, so the branch checked out at `../Barrier3D` is the model. Every run records the branch and commit it imported (`run_registry.barrier3d_provenance`). None of the fixes below has been pushed upstream (Hannah, 2026-09-28).

**To do: push the Barrier3D fixes eventually** (Hannah, 2026-09-29). They stay local for now. Until they are pushed, this branch runs only on a machine that has `hatteras/adopted` checked out at `../Barrier3D`: the runner refuses a Barrier3D without per-cell ceilings. The full record, with evidence paths, is `HATTERAS_FIXES.md` on `hatteras/adopted`. Changes to CASCADE's own `cascade/` package are recorded in `HATTERAS_CASCADE_CHANGES.md` at the repository root.

| branch / tag | what it has | in use? |
|---|---|---|
| **`hatteras/adopted`** (2b8f167 code; 8a588ea records it) | everything below merged: the route_overwash fix, the three overwash fixes, and per-cell dune ceilings (`DuneCeilingFromStart`, switched on in `data/hatteras_init/Hatteras-CASCADE-parameters.yaml`) | **yes, since 2026-09-28**: `../Barrier3D` is on it. The runner refuses a Barrier3D without per-cell ceilings. The matrix before this is in `output/raw_runs/archive/2026-09-28-pre-ceiling/`. |
| `feature/per-cell-dune-ceiling` (d343461) | per-cell dune ceilings alone, on 49fd069; worktree `../Barrier3D-dune-ceiling` | merged into `hatteras/adopted` |
| `fix/overwash-gaps-momentum` (db0ba30, tag `hat-fix-overwash-gaps-momentum`; worktree `../Barrier3D-overwashfix`) | `DuneGaps` dropping cells; the gap discharge slice; the inundation momentum constant reset to 0 | merged into `hatteras/adopted` |
| `fix/route-overwash-axis-swap` (49fd069, tag `hat-fix-route-overwash`) | the `route_overwash` index swap: wrong cells read, and a silent crash on long storms | merged; every run from 2026-09-24 to 09-28 used it alone |

The storm series moved to `v3_trim24` on the same day (`hat_env_forcings.DEFAULT_STORM_VARIANT`). Every run records its storm file and dune-ceiling mode in its metadata. To reproduce a run made before 2026-09-28, check out `fix/route-overwash-axis-swap` and use that date's parameter template and `v3_72` storms.

`../Barrier3D-prefix-ce36866` is a detached worktree of the code **before** the route_overwash fix. It is kept only for the storm-duration cause test (`output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-max-duration/`) and can be removed with `git -C ../Barrier3D worktree remove ../Barrier3D-prefix-ce36866`.

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

The five figure scripts that were here (`hindcast_final_figure_lowess`,
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
