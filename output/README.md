# `output/` - what is in here, and where to look

Everything under here is **produced**. The one exception is
`calibration/groin/joint_fit.json`, an input the pipeline reads back. Rewritten
2026-09-18, after superseded material, logs and calibration were each given one home.

## Where do I find...

| I want... | look in |
|---|---|
| a model run | `raw_runs/`. Resolve it with `cascade_pipeline.run_registry`; never build the path by hand |
| a figure for the manuscript or a talk | `figures/<subject>/` (captions in `supporting/CAPTIONS.md`), talk versions in `figures/talk/` |
| a figure of the observed record, not of a run | `observations/` |
| a cross-run comparison (rates, relocation, scenario grid) | `comparisons/<question>/` |
| why a calibration value is what it is | the three decision records: `calibration/groin/SELECTED_M60_f0.60/README.md` (M and f), `calibration/hs/DECISION.md` (Hs = 2.5), `calibration/sensitivity/figures/README.md` |
| a one-off study | `raw_runs/experiments/<date>-<tag>/` (runs + `NOTE.md`), or `experiments/` for older studies with their own drivers |
| what a batch did overnight | `logs/driver/driver_manifest.jsonl` |
| anything retired | `archive/`, or `raw_runs/archive/` for retired runs |

## The map

| directory | written by | holds |
|---|---|---|
| `raw_runs/` | `scripts/hatteras_ms/HAT_hindcast_1984_2024.ipynb` (and its headless twin) | every model run, filed by purpose (see below), and `run_index.csv` |
| `comparisons/` | `HAT_relocation_comparison.py`, `HAT_scenario_grid.py`, `HAT_rate_windows.py`, the hindcast figure | cross-run figures and tables, one folder per question. Its README has the tracking rule |
| `figures/` | `scripts/figure_making/` via `site_layer/hat_figure_style.FIGURES_ROOT` | finished figures, by subject: `site`, `forcing`, `initialization`, `management`, `shoreline`, `style`, `talk` |
| `observations/` | `scripts/figure_making/shoreline/chainage/` | figures of the observed record. ~1200 files, mostly gif frames |
| `calibration/groin/` | `scripts/hatteras_ms/groin-sweep/` (paths from `HAT_groin_sweep_config.GROIN_SWEEP_ROOT`) | the (M, f) calibration. Holds `joint_fit.json`, the input exception. The study's code is in `hard-structures/groin/` |
| `calibration/hs/` | `scripts/input_prep/7-source-sink/2-calibrate/HAT_be_zone_residual_fit.py`, `scripts/sensitivity_analysis/HAT_plot_hs_experiment.py` | the Hs 3.0 test, its `DECISION.md`, and its arms in `runs/` |
| `calibration/sensitivity/` | `scripts/sensitivity_analysis/` | the parameter sweep's manifests, logs and figures. The sweep's runs are in `raw_runs/sensitivity/` |
| `calibration/groin_rig/` | `hard-structures/groin/.../HAT_groin_hindcast_1967_2017.py`, read by `scripts/hatteras_ms/groin-sweep/` | the 1967-2018 groin rig, the only window spanning the deterioration ramp |
| `experiments/` | `scripts/hatteras_ms/experiments/` | older one-off studies that have their own drivers and layouts |
| `logs/driver/` | `scripts/hatteras_ms/HAT_run_all.py`, `tools/HAT_rerun_arm.py` | the unattended driver's manifest, stdout and per-job logs |
| `logs/scratch/` | by hand | terminal captures from one-off comparisons. Untracked, safe to delete |
| `archive/` | nothing (moved aside by hand) | retired material, `YYYY-MM-DD_<what>/`. **Do not use for analysis** |

### Where things go

- **Logs.** A run's log goes in its run directory. A study's logs go in that
  study's `logs/`. Driver batches go in `logs/driver/`. Hand captures go in
  `logs/scratch/`. Nothing else gets loose `.log` files.
- **Retired material.** It goes to `archive/YYYY-MM-DD_<what>/` with a `WHY.md`.
  Retired *runs* are the exception: they go to `raw_runs/archive/`, because
  `archive` is a registry kind and they stay indexed there.
- **Dates.** ISO `YYYY-MM-DD` in every new folder name.

## How a run is addressed

    raw_runs/matrix/<start>_<end>/<preset>/<run_name>/
    raw_runs/sensitivity/<axis>/<start>_<end>/<preset>/<run_name>_<token>/
    raw_runs/experiments/<tag>/<start>_<end>/<preset>/<run_name>/
    raw_runs/versions/<tag>/<start>_<end>/<preset>/<run_name>/
    raw_runs/archive/<date>-<reason>/<start>_<end>/<preset>/<run_name>/

The **name describes the scenario**. The **path describes the purpose**. The run name
is derived from the runner's management switches and is never typed. A
sensitivity cell is its baseline's name plus one trailing token. Every run has a
`kind` and, except in the matrix, a `tag`, set through `HAT_RUN_KIND` /
`HAT_RUN_TAG`. `raw_runs/README.md` has the whole design.

**Never build one of these paths by hand.** `cascade_pipeline.run_registry`
(`preset_dir_for`, `run_dir_for`, `find_run_dir`, `run_dir_for_index_row`)
resolves them. When a run isn't where you asked, it raises and names the places
the run *is*, rather than returning a path that doesn't exist.

`run_index.csv` is keyed on **(run_name, kind, tag)** and is DERIVED. It is rebuilt
from disk by `tools/HAT_index_runs.py`. On 2026-09-18 it held 314 runs: 88 matrix,
103 sensitivity, 87 experiment, 12 version and 24 archive. `retired_runs.csv`
records runs whose directories are gone.

## Three number spellings, all deliberate

Do not "fix" these. Figure scripts hardcode the paths.

| tree | spelling | why |
|---|---|---|
| `raw_runs/` | `waveHs1p2`, `rset40` | `.` becomes `p` and `-` becomes `m`, so a name can be split on `_` and read by eye (`cascade_pipeline.hindcast._number_token`) |
| `calibration/groin/` | `M110_be-10_f0.60` | predates the rule. Its README records that renaming would break the figure scripts |
| `calibration/hs/` | `02_zones_Hs2p5` | ordinal prefixes, because the stages are read in order |

## What is tracked, and what is not

`output/` is ignored wholesale, except for these re-includes:

- `raw_runs/`: the small text products (rate CSVs, road summary, run metadata,
  `run_index.csv`) and the ~20 KB `*_shoreline_matrix.npy`. The `.npz`, GIFs and
  PNGs are regenerable and are ignored.
- `comparisons/`: every README, CAPTIONS.md, `.csv` and `.txt`. The images are
  ignored.
- This README, and the READMEs and WHY.md files under `archive/`.
- `calibration/`: every README and DECISION.md, `groin/joint_fit.json`, and
  `groin/SELECTED_M60_f0.60/README.md`. That covers all three decision records.

**Consequence:** the sweeps, runs and figures under `calibration/`, `figures/` and `logs/` exist only on the machine that made them.
If a decision record needs to survive a fresh clone, re-include it in
`.gitignore`.

## Size

About 30 GB, nearly all of it `.npz` model state in `raw_runs/`: 21 GB in `matrix/`
and 4.9 GB in `archive/`. A run's `.npz` is ~250 MB and is written only when
`save_model_state` is on. It is read by the planview GIF, the relocation
comparison and the source/sink calibration, so it isn't dead weight. But a run
takes ~1.5 minutes, so regenerating one is cheap if disk gets tight. The obvious
place to reclaim space is the `.npz` files under `raw_runs/archive/`.
