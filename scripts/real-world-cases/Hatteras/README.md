# Hatteras Island hindcast (Hannah Henry)

The scripts for the CASCADE hindcast of Hatteras Island, NC (dissertation
chapter 1, Coastal Environmental Change Lab, UNC Chapel Hill). This folder
holds **scripts only**. Data, model runs, figures and every other output live
outside the repository, in a workspace on the D: drive.

Moved here from the branch `hannahaline/hatteras-cascade` on 2026-10-23. That
branch keeps the full history of everything below, including the data and
outputs that were tracked there.

## Where things are

| what | where |
|---|---|
| scripts (this folder) | `scripts/real-world-cases/Hatteras/` on `obx-cascade` |
| data, runs, figures, study results | the workspace, `D:\hatteras_ch1\model_runs\` |
| raw GIS sources, ArcGIS project | `D:\hatteras_ch1\gis_sources\`, `D:\hatteras_ch1\gis_project\` |
| the model | this repository's `cascade/`, plus [`cascade_patch.py`](cascade_patch.py) |
| Barrier3D | published Barrier3D plus [`barrier3d-patch/`](barrier3d-patch/README.md) |

The workspace keeps the old repository's layout, so a path there reads the
same as it did in the repository: `data/hatteras_init/...`, `output/raw_runs/...`,
`hard-structures/groin/...`, and `scripts/...` for the non-code files that
used to sit beside scripts.

## How a script finds its files

Every script reads its locations from [`hatteras_paths.py`](hatteras_paths.py):

- `HATTERAS_DIR`: this folder, for imports and for running other scripts.
- `WORKSPACE`: the D: workspace, for every data and output path. Set the
  environment variable `HATTERAS_WORKSPACE` to use another folder (a second
  drive under another letter, say).
- `workspace_twin(folder)`: where a script folder keeps its non-code files.
  A study that wrote `figures/` beside its script writes them to the same
  place in the workspace instead.

So nothing a script writes lands in the repository. A script that does is a
bug.

## Layout

The folders keep the names they had in the old repository's `scripts/`:

| folder | what |
|---|---|
| `hatteras_ms/` | the hindcast runner (`HAT_hindcast.py`, with its setup steps in `hindcast/HAT_hindcast_setup.py`; `HAT_hindcast.ipynb` walks through one run), its settings (`hat_run.yaml`), the matrix driver, calibration solves, experiments (by theme) and tools |
| `cascade_pipeline/` | the runner's library: domains, nourishment, roadway, shoreline, run registry, plotting |
| `site_layer/` | Hatteras site facts and resolvers: topography versions, forcings, observed rates, figure style |
| `input_prep/` | the scripts that built the model inputs, by stage (0-elevation ... 8-overwash-analysis) |
| `figure_making/`, `analyze_output/`, `sensitivity_analysis/` | figures and run comparisons |
| `hard-structures/` | the groin and revetment studies (their results are in the workspace) |
| `repo_tools/`, `other_ms/` | style checks and templates; the Pea Island manuscript scripts |
| `tests/` | `test_groin.py`, `test_hatteras_config.py` |
| `docs/` | project notes carried from the old repository (`UNITS.md`, `ORGANIZATION.md`, ...). Paths inside them describe the old layout |
| `barrier3d-patch/`, `cascade-patch/`, `cascade_patch.py` | the model fixes the runs need beyond this repository's code |

## Running the hindcast

1. Use this clone's own environment (`.venv` at the repository root). It
   needs CASCADE installed from this clone and Barrier3D with the patch
   applied ([`barrier3d-patch/README.md`](barrier3d-patch/README.md)).
2. The D: drive must be attached; `hatteras_paths.py` stops with a message
   if the workspace is missing.
3. Choose the run in `hatteras_ms/hat_run.yaml`, then run
   `hatteras_ms/HAT_hindcast.ipynb` (or `HAT_hindcast.py`). Results go to
   `D:\hatteras_ch1\model_runs\output\raw_runs\`.

## How the runs differ from the old branch

The move to obx-cascade changed the model as well as the paths (decided
2026-10-23):

- Fills on road domains go through `RoadwayManager`. obx-cascade allows one
  manager per domain. A domain whose road drowns is handed to
  `BeachDuneManager` the year after, so later fills are still placed.
- Road dunes follow obx-cascade's rules, which rebuild only within 10 m of the
  road. With NC-12 mostly ~200 m back, no road dunes are rebuilt on Hatteras.
- The dune cap clips the whole dune cell (obx-cascade).
- obx-cascade's Outwasher is switched off.
- One crash fix, `resize_interior_domain`, is patched in at run time
  ([`cascade-patch/`](cascade-patch/README.md)).

Same settings, 2009 start: island-wide RMSE 2.12 m/yr on the old branch,
2.47 on obx-cascade; overwash skill 0.12 against 0.23.

## Working on a shared branch

`obx-cascade` is shared with other lab members. Commit only paths under this
folder, pull before pushing, and never force-push.

## The scripts

How a script is laid out, commented and signed is [`STYLE.md`](STYLE.md).

The sections below come from the old repository's `scripts/README.md`, updated for this folder.

### `site_layer/`: the seam

Not stray one-offs. The seam between the site-agnostic `cascade_pipeline/`
package and the task folders, which all consume Hatteras' own paths, place
names and house style. Each module answers one question, in one place.

| Module | Answers | Imported by |
|---|---|---|
| `hatteras_site_config.py` | *What is this place* — domain geometry, town spans, periods, BE presets, road events, nourishment projects. | 87 |
| `hat_figure_style.py` | *What does a Hatteras figure look like* — one type scale, one palette, one elevation ramp, one way to letter a panel. | 102 |
| `hat_topo_version.py` | *Which* Barrier3D domains does this script read — which product (`1984-start` / `2004-start` / `forecast`), and which version inside it? | 77 |
| `hat_observed_rates.py` | *Where is the observed shoreline* — the CoastSat chainage, the per-window rate fits the model is graded against, and the transect-to-domain lookup. | 22 |
| `hat_elevation_products.py` | *Which* elevation product and stage — `2009-2014` or `2009-2014-1996`, gapfilled 1 m or resampled 10 m? | 13 |
| `hat_extension_domains.py` | *Which alongshore reach does a run model* — the named span, and the 500 m bins that number the coast beyond the 90 surveyed domains. | 11 |

**These resolvers exist because a hand-built path fails silently.** Four road
scripts once hardcoded `2009-dune-topo/2009_v3`; when the dune windows were
re-picked into `v4` they kept reading `v3` interiors while consuming `v4`
setbacks, and 18 domains — including two of the three managed roadways — had
their drown verdicts computed on the wrong grid, with nothing raised. The same
shape of bug hit `HAT_road_elevation.py` when a fill source moved under
`superseded/`. Each location is resolved **once**, there, and a name that is not
on disk is an immediate, loud error listing what is.

So: **never rebuild one of these paths by hand.** Call `topo_dirs()`,
`domain_arrays()`, `array_path()`, or `product()` and let it raise.

[`site_layer/README.md`](site_layer/README.md) covers the rest: why `__init__.py` re-exports nothing,
and the `__package__` guard the two runnable modules need.

They must **not** move into `cascade_pipeline/`. That package deliberately ships
no site content — a different study site writes its own sibling of `site_layer/`
and never touches the package. One leak already exists
(`cascade_pipeline/hindcast.py:61` imports `hat_topo_version`); folding the site
config in would make the package permanently Hatteras-only.

### `repo_tools/`: tools that act on the tree

`hat_layout_check.py` audits the repo against the seven rules in
`ORGANIZATION.md`. Nothing imports it; it is run
(`python scripts/real-world-cases/Hatteras/repo_tools/hat_layout_check.py`). It was written for the old repository's layout and has not been adapted and it always exits zero.

## The folders

```
site_layer/            what and where Hatteras is: the six site modules the
                       task folders all consume. Import the module, never the
                       package -- site_layer/README.md says why.

repo_tools/            tools that act on the repository itself rather than on
                       the science. hat_layout_check.py, run and never
                       imported.

cascade_pipeline/      the library: post-run analysis and figures.
                       Site-agnostic by contract — geometry, shoreline
                       extraction, CoastSat/LOWESS, plotting/. Consumes a
                       finished run; does not drive the simulation.

input_prep/            builds the model's inputs. Stage folders 0-7 mirror
                       data/hatteras_init/ (in the workspace) one-for-one, so a script sits
                       under the same number as the product it writes.
                       (8-overwash-analysis/ has no data counterpart, and
                       note the folder here is 4-mgmt-forcing*s*, plural,
                       against the singular one in data/.)

hatteras_ms/           the manuscript runs. HAT_hindcast.py runs one; its
                       setup steps are in hindcast/HAT_hindcast_setup.py and its figures
                       are drawn after, from the run folder
                       (figure_making/model_output/run_figures.py).
                       HAT_hindcast.ipynb walks through one run step by step,
                       with its QC plots. Settings are typed
                       in hat_run.yaml and read by hindcast/HAT_hindcast_config.py
                       (env var > yaml > default); HAT_run_matrix.py drives the
                       matrix, without and with the groin, unattended.

sensitivity_analysis/  standalone sweep driver and plotter, working off the
                       run registry rather than off a live run.

analyze_output/        read-only comparisons across finished runs —
                       compare_runs/, overwash/.
                       Scripts only: their products go to
                       output/comparisons/<topic>/ in the workspace, never beside the script.

figure_making/         the older figure tree, largely predating
                       cascade_pipeline/. New figures belong with the run
                       that produces them, or in cascade_pipeline/plotting/
                       if they generalise.

other_ms/              other manuscripts, not Hatteras: pea_island_ms/.
```

Several of these folders carry their own README with the decisions behind
them — `input_prep/0-elevation/`, `input_prep/1-barrier3d-domains/`,
`input_prep/4-mgmt-forcings/road_offset/` and `.../road_relocation/` are the
substantial ones. Read those before changing a forcing.

## Retired code is deleted

A script that is retired is deleted, not parked in a folder (ORGANIZATION.md
rule 4, since 2026-10-01). Git keeps it. The folder's `README.md` gets a short
"Deleted <date>" record naming each file, why it went, and the recovery
commands:

```
git log --diff-filter=D --oneline -- <path>
git show <commit>^:<path>
```

Before deleting, check that the file is committed (`git cat-file -e
HEAD:<path>`). A file that is only staged, or never committed, is not
protected. Check also that nothing live imports it or loads it by path, and
that any data folder beside it is either tracked or no longer wanted, because
git-ignored data deleted with the code is gone for good.

History: until 2026-09-02 retired code sat in `old/`-style folders (six
spellings across fourteen folders), then in `old_<what>/`, then in
`superseded_<date>/`. On 2026-10-01 the last fourteen `superseded_<date>/`
folders under `scripts/` were deleted. Each one is recorded in the README of
the folder that held it.

## The import idiom

A script finds this folder by walking up to `hatteras_paths.py`, puts it on
`sys.path`, and takes its data and output locations from there:

```python
import sys
from pathlib import Path

HATTERAS = next(p for p in Path(__file__).resolve().parents
                if (p / "hatteras_paths.py").exists())
sys.path.insert(0, str(HATTERAS))

from hatteras_paths import WORKSPACE, workspace_twin          # noqa: E402
from site_layer.hat_topo_version import topo_dirs, array_name  # noqa: E402
```

No folder counting (`parents[N]`): the upward search does not depend on how
deep a script sits, and a script moved between folders keeps working. In the
old repository the same search looked for `pyproject.toml`; here
`pyproject.toml` belongs to the shared repository and the marker is
`hatteras_paths.py`. Data paths start from `WORKSPACE`, never from the script's
own folder; a file that used to sit beside a script is reached through
`workspace_twin(<script folder>)`. The same anchor is what makes
`from cascade_pipeline.domains import DomainGeometry` resolve; the package is
imported by path, not installed.
