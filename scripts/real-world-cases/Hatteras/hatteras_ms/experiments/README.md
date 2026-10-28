# experiments

One-off studies, each asking one question of the hindcast. A study drives the
unchanged runner (`../HAT_hindcast.py`) as a subprocess with
`HAT_IGNORE_SETTINGS=1`, or in-process with one input swapped, and files its
runs under `output/raw_runs/experiments/<theme>/<date>-<name>/`. Each script's
header says what it asks and how to run it; `run` / `score` / `plot` actions
are separate so a finished sweep can be re-scored without re-running.

The studies are filed by theme. On 2026-10-23, when the hindcast moved to CASCADE-obx, only the
calibrations that may need re-running on the new model were kept; the settled studies were retired
(`RETIRED_STUDIES.md` lists them, keeps their READMEs, and says how to recover any script from git).
Their runs and write-ups stay on D: under `output/raw_runs/experiments/<theme>/`.

| Folder | Scripts | What it calibrates |
|---|---|---|
| `waves/` | `HAT_wave_grid.py` | the wave climate for each window |
| `dunes/` | `HAT_stable_dune_growth_rate.py`, `HAT_test_period_growth_rate.py` and their figures, plus the seven modules they import | the dune growth rate for each window |
| `ends/` | `HAT_buffer_closure_test.py`, `HAT_cape_source_shape.py` | the buffer ring closure; the Cape Point source shape |

The helpers the matrix analysis and figures import (`coastsat_target`, `interior`,
`smooth_like_target`) are in `site_layer/hat_scoring.py`.

## Importing across folders

Many studies import another's module by bare name and re-point its folders
(the `_div10` / `_option_a` variants, the storm studies' shared `MD` / `S`,
the dune studies' `U`), often across themes. `experiment_paths.py` puts every
theme folder on `sys.path`, so a script reaches it by searching upward and
importing it:

```python
sys.path.insert(0, str(next(_p for _p in _HERE.parents if (_p / "experiment_paths.py").exists())))
import experiment_paths  # noqa: E402,F401  (every experiments/<theme>/ folder onto sys.path)
```

A script outside `experiments/` puts `hatteras_ms/experiments` on the path and
imports it the same way. Moving a study to another theme folder needs no
import changes; giving two studies the same file name would, and
`python experiment_paths.py` reports that.

`superseded_20260907/`, the retired 1984 seaward row-insert set that could
not be re-run, was deleted on 2026-10-01 (recover from git; see
`../README.md`, "Deleted 2026-10-01"). The crest experiment and the two
relocation scoring scripts were retired with the 1984 and 2004 starts
(`fe452b83`, 2026-10-19); see `../README.md`, "Deleted 2026-10-23".
