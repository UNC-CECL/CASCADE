# `output/` - what is in here, and where to look

Everything under here is **produced**, and everything here is **current**. The one exception
to "produced" is `calibration/groin/joint_fit.json`, an input the pipeline reads back.

Superseded material is not kept here. On 2026-10-07 it all moved to the Seagate drive,
`D:\CASCADE_offload\output\`, under the same relative paths (see the bottom).

## Where do I find...

| I want... | look in |
|---|---|
| the current model runs | `raw_runs/matrix/`: `1996_2009` (calibration), `2009_2025` (test), `1996_2025` (full window). `raw_runs/matrix/README.md` says which preset is which. Resolve a run with `cascade_pipeline.run_registry`; never build the path by hand |
| how a run compares to CoastSat | the run's own `figures/shoreline_position_change_with_buffers.png` (calibration and test), or `comparisons/full_window/` (1996-2025) |
| the plan, scores and open decisions | `scripts/hatteras_ms/DEM_TO_DEM_CALIBRATION.md` |
| a study under the current plan (end solves, BE set 1, groin fit, 2021 step...) | `raw_runs/experiments/<topic>/<date>-<what>/`; `raw_runs/experiments/README.md` is the map |
| the pinned groin | `calibration/groin/joint_fit.json` and its README |
| a figure for the manuscript or a talk | `figures/<n>-<subject>/` (captions in `supporting/CAPTIONS.md`), talk versions in `figures/talk/` |
| what a batch did | `logs/driver/driver_manifest.jsonl`; the DEM-to-DEM and full-window logs are in `logs/driver/dem_to_dem/` and `logs/driver/full_window_1996_2025/` |
| anything older | `D:\CASCADE_offload\output\`, same path as it had here |

## The map

| directory | written by | holds |
|---|---|---|
| `raw_runs/` | `scripts/hatteras_ms/HAT_hindcast_1984_2024.ipynb` (and its headless twin) | the current runs (`matrix/`), the current studies (`experiments/`), and `run_index.csv` |
| `comparisons/` | `scripts/analyze_output/compare_runs/` | cross-run figures and tables on the current windows, one folder per question |
| `figures/` | `scripts/figure_making/` via `site_layer/hat_figure_style.FIGURES_ROOT`, `regenerate_all_figures.py` | finished figures by subject: `1-site` to `5-results`, `style`, `talk`. Not yet redrawn for the current windows (some model-input figures still show 1996_2010 / 2010_2024) |
| `calibration/groin/` | the 2026-10-05 pin, by hand | `joint_fit.json` (pipeline input) and its README |
| `logs/driver/` | `scripts/hatteras_ms/HAT_run_all.py` and the DEM-to-DEM drivers | the driver's manifest and stdout, and one folder of logs per current batch |

### Where things go

- **Logs.** A run's log goes in its run directory. A study's logs go in that
  study's folder. Driver batches go in `logs/driver/<batch>/`.
- **Superseded material.** Copy it to `D:\CASCADE_offload\output\` under the same
  path, check every file arrived at full size, then delete it here. Git-tracked
  text is removed in a commit, so history keeps it.
- **Dates.** ISO `YYYY-MM-DD` in every new folder name.

## How a run is addressed

    raw_runs/matrix/<start>_<end>/<preset>/<run_name>/
    raw_runs/experiments/<tag>/<start>_<end>/<preset>/<run_name>/

(`sensitivity/`, `versions/` and `archive/` are also registry kinds; none is on C: now.)

The **name describes the scenario**. The **path describes the purpose**. The run name
is derived from the runner's management switches and is never typed. `raw_runs/README.md`
has the whole design. `run_index.csv` is keyed on **(run_name, kind, tag)** and is
DERIVED: `tools/HAT_index_runs.py` rebuilds it from disk, and `retired_runs.csv`
records runs whose directories are gone (including everything moved to D:).

## Number spellings

`raw_runs/` names spell numbers with `p` for `.` and `m` for `-` (`waveHs1p2`), so a
name can be split on `_` and read by eye (`cascade_pipeline.hindcast._number_token`).
Groin cells (`b0.60_f0.6`) predate the rule. Do not "fix" either: scripts hardcode them.

## What is tracked, and what is not

`output/` is ignored wholesale, except for re-includes in `.gitignore`: the small
text products of `raw_runs/` (rate CSVs, run metadata, `run_index.csv`, the ~20 KB
`*_shoreline_matrix.npy`), every README / CAPTIONS.md / `.csv` / `.txt` under
`comparisons/`, this README, and `calibration/**/README.md` plus `groin/joint_fit.json`.
Figures and logs exist only on this machine.

## Moved to D: on 2026-10-07

| what | why |
|---|---|
| `raw_runs/matrix/1996_2010`, `2010_2024`, `1996_2015`, `2010_2026`, `figures/` | windows superseded on 10-03 and by the 10-05 calibration/test plan |
| `raw_runs/sensitivity/`, `raw_runs/versions/`, `raw_runs/archive/`, `SUPERSEDED_CANDIDATES.md` | the 09-28 wave sweep, version checks and retired runs, all on old windows |
| `raw_runs/experiments/`: every study dated before 2026-10-05 | see `raw_runs/experiments/OFFLOADED.md` |
| `comparisons/` except `full_window/` | all on superseded windows |
| `calibration/hs`, `sensitivity`, `groin_rig`, and the dipole groin sweep | superseded calibrations (÷10 offset, dipole groin) |
| `archive/`, `observations/`, `logs/scratch/`, per-job driver logs before 10-05 | retired material, the 09-10 gif frames, old captures |

The earlier offload the same morning (`raw_runs/archive` and every experiment `.npz`)
is described in `raw_runs/experiments/OFFLOADED.md`. Figure producer
`scripts/figure_making/pipeline/7-source-sink/be_method_figures.py` reads `ends.json`
from three September end-solve studies that are now on D:; it runs only with the
drive plugged in.
