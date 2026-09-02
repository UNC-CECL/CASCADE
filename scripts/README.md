# scripts

Everything that prepares CASCADE's Hatteras inputs, drives its runs, and turns
the output into figures.

The model itself is **not** here — it is `cascade/` at the repo root, installed
as a package (`setup.py`, name `cascade`). Nothing under `scripts/` is packaged.
These scripts reach each other by putting **`scripts/` on `sys.path`**, and that
one fact explains the whole layout below, including why three modules sit loose
at this level.

## The three modules at the root

They are not stray one-offs. They are the **site layer** — the seam between the
site-agnostic `cascade_pipeline/` package and the task folders that all consume
Hatteras' own paths and place names.

| Module | Answers |
|---|---|
| `hat_topo_version.py` | *Which* Barrier3D domains does this script read — which product (`1984-start` / `2004-start` / `forecast`), and which version inside it? |
| `hat_elevation_products.py` | *Which* elevation product and stage — `2009-2014` or `2009-2014-1996`, gapfilled 1 m or resampled 10 m? |
| `hatteras_site_config.py` | *What is this place* — domain geometry, town spans, periods, BE presets, road events, nourishment projects. |

**Both resolvers exist because a hand-built path fails silently.** Four road
scripts once hardcoded `2009-dune-topo/2009_v3`; when the dune windows were
re-picked into `v4` they kept reading `v3` interiors while consuming `v4`
setbacks, and 18 domains — including two of the three managed roadways — had
their drown verdicts computed on the wrong grid, with nothing raised. The same
shape of bug hit `HAT_road_elevation.py` when a fill source moved under
`superseded/`. Each location is now resolved **once**, in these modules, and a
name that is not on disk is an immediate, loud error listing what is.

So: **never rebuild one of these paths by hand.** Call `topo_dirs()`,
`domain_arrays()`, `array_path()`, or `product()` and let it raise.

### Why they live at the root and not in a folder

Two reasons, both mechanical:

1. **`scripts/` is the `sys.path` anchor.** Some 58 files insert it and then
   import by bare name — `from hat_topo_version import topo_dirs`. Root is what
   makes that spelling work from every task folder.
2. **They pin the repo root.** Both resolvers compute
   `PROJECT_ROOT = Path(__file__).resolve().parents[1]`, i.e. `scripts/ → repo
   root`, and hang `INIT_ROOT = data/hatteras_init` off it. Move either module
   one level deeper and `parents[1]` silently becomes `scripts/` — exactly the
   quiet-wrong-path failure they were written to end.

Moving them into a package is defensible but is not a file move: it is ~76
import lines across ~58 files (the hindcast notebook included) plus a fix to
`parents[1]` in both resolvers.

They must **not** move into `cascade_pipeline/`. That package deliberately ships
no site content — a different study site writes its own sibling of
`hatteras_site_config.py` and never touches the package. One leak already exists
(`cascade_pipeline/hindcast.py:61` imports `hat_topo_version`); folding the site
config in would make the package permanently Hatteras-only.

## The folders

```
cascade_pipeline/      the library: post-run analysis and figures.
                       Site-agnostic by contract — geometry, shoreline
                       extraction, CoastSat/LOESS, plotting/. Consumes a
                       finished run; does not drive the simulation.

input_prep/            builds the model's inputs. Stage folders 0-7 mirror
                       data/hatteras_init/ one-for-one, so a script sits
                       under the same number as the product it writes.
                       (8-overwash-analysis/ has no data counterpart, and
                       note the folder here is 4-mgmt-forcing*s*, plural,
                       against the singular one in data/.)

hatteras_ms/           the manuscript runs. HAT_hindcast_1984_2024.ipynb is
                       the source of truth; the .py beside it is a headless
                       mirror with no features of its own. Settings are typed
                       in hat_run.yaml and read by HAT_hindcast_config.py
                       (env var > yaml > default); HAT_run_all.py drives the
                       matrix and the sweep unattended. groin-sweep/ holds
                       the groin calibration and its figures.

sensitivity_analysis/  standalone sweep driver and plotter, working off the
                       run registry rather than off a live run.

analyze_output/        read-only comparisons across finished runs —
                       compare_runs/, overwash/, smoothing_vs_cascade/.

figure_making/         the older figure tree, largely predating
                       cascade_pipeline/. New figures belong with the run
                       that produces them, or in cascade_pipeline/plotting/
                       if they generalise.

other_ms/              other manuscripts, not Hatteras: chom_ms/,
                       ocracoke_ms/, pathways_ms/.
```

Several of these folders carry their own README with the decisions behind
them — `input_prep/0-elevation/`, `input_prep/1-barrier3d-domains/`,
`input_prep/4-mgmt-forcings/road_offset/` and `.../road_relocation/` are the
substantial ones. Read those before changing a forcing.

## Superseded code: `old_<what it holds>/`

A folder holding retired scripts is named **`old_` plus what it holds** — one
rule, everywhere under `scripts/`:

```
figure_making/old_dsas_scripts/            input_prep/old_source_sink_search/
figure_making/old_plot_tests/              input_prep/1-barrier3d-domains/old_extractors/
figure_making/shoreline_change/old_rate_analysis/
hatteras_ms/old_drafts/                    input_prep/5-scr/CoastSat/old_dsas_comparisons/
hatteras_ms/old_versions/                  input_prep/5-scr/CoastSat/old_time_periods/
sensitivity_analysis/old_guides/           input_prep/6-scr-smooth/.../old_smoothing/
input_prep/4-mgmt-forcings/road_offset/1-produce/old_method/
```

Before 2026-09-02 there were six spellings across fourteen folders — `old/`,
`drafts/`, `old_versions/`, `old_method/`, `old_plot_tests/`,
`old_time_periods/` — and seven of them were bare `old/`, including an `old/`
nested inside another `old/`. Nothing referenced them by path, so the renames
broke nothing.

**The descriptor is not decoration.** `figure_making/` holds two retired
folders side by side, so a single uniform name — every one of them
`old_method/` — would collide. Keeping what each holds is also what makes a
bare listing legible: `old_extractors` and `old_smoothing` say which stage they
were retired from, where two `old/` entries at different depths do not.

**`old_method/` under `road_offset/1-produce/` is the one to leave alone.** It
is not dormant: `road_offset/README.md` cites it three times as the legacy
method kept for comparison against `HAT_road_offset_from_dune_start.py`, and
that comparison is the argument for the current setbacks.

Retired code is kept, not deleted, and all of it is tracked — so anything here
is recoverable from history if a folder is ever removed. What is NOT uniformly
tracked is the data beside it: several of these folders hold figures and CSVs
that exist nowhere else and are not in git. Check before deleting one.

## The import idiom

A task script reaches the site layer by walking up to `scripts/`:

```python
import sys
from pathlib import Path

# parents[N] IS scripts/ -- count the folders between this file and it.
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from hat_topo_version import topo_dirs, array_name   # noqa: E402
from hatteras_site_config import HATTERAS_DOMAINS    # noqa: E402
```

`N` is positional. **Moving a script between folders changes it**, and getting
it wrong yields an `ImportError` at best or a wrong `PROJECT_ROOT` at worst —
so re-count it whenever a script changes depth, and keep the comment naming
which parent is `scripts/`.

The same anchor is what makes `from cascade_pipeline.domains import
DomainGeometry` resolve; the package is imported by path, not installed.
