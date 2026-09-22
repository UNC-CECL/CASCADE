# scripts

Everything that prepares CASCADE's Hatteras inputs, drives its runs, and turns
the output into figures.

The model itself is **not** here — it is `cascade/` at the repo root, installed
as a package (`setup.py`, name `cascade`). Nothing under `scripts/` is packaged.
These scripts reach each other by putting **`scripts/` on `sys.path`**, and
that one fact explains the layout below.

## The root is folders and this file

Nothing else. On 2026-09-18 it held eleven loose files — six site modules, the
layout checker, and four `hatteras_site_config_prebe_<stamp>.py` snapshots of
the solved BE field. They are in folders now:

```
site_layer/     what and where Hatteras is -- the six site modules
repo_tools/     tools that act on the repository itself
```

The snapshots were never code at all: the calibrate step writes one before each
pass, beside the config it copies. They read as stray duplicates, which is how
the equivalent 2026-08-24 snapshot came to be discarded, taking the only record
of the one-shot solve with it and leaving `plot_be_zones.py` undrawable for
three weeks. They now live in
`data/hatteras_init/7-source-sink/2-calibrate/prebe/`, and the apply step writes
new ones there.

`hat_layout_check.py` holds the root to this shape: Rule 1's second case reports
any file here but `README.md`, saying which folder it should have gone to.

### `site_layer/` — the seam

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

`site_layer/README.md` covers the rest: why `__init__.py` re-exports nothing,
and the `__package__` guard the two runnable modules need.

They must **not** move into `cascade_pipeline/`. That package deliberately ships
no site content — a different study site writes its own sibling of `site_layer/`
and never touches the package. One leak already exists
(`cascade_pipeline/hindcast.py:61` imports `hat_topo_version`); folding the site
config in would make the package permanently Hatteras-only.

### `repo_tools/` — tools that act on the tree

`hat_layout_check.py` audits the repo against the seven rules in
`ORGANIZATION.md`. Nothing imports it; it is run
(`python scripts/repo_tools/hat_layout_check.py`) and it always exits zero.

### Why the move was safe

The old warning here was that both resolvers computed
`PROJECT_ROOT = Path(__file__).resolve().parents[1]`, so dropping a module one
level deeper would make `parents[1]` silently mean `scripts/` — the quiet
wrong-path failure they were written to end. **That is no longer how they
resolve.** Every one of them searches upward for `pyproject.toml`, which is
Rule 5 of `ORGANIZATION.md`, so depth does not enter into it.

What remained was import spelling, and it was mechanical: 302 lines in 170
files, `from hat_topo_version import x` becoming
`from site_layer.hat_topo_version import x`. `scripts/` is still the `sys.path`
anchor and every consumer already inserts it, so no `sys.path` line changed.
Three scripts that read `hatteras_site_config.py` as a **file** rather than
importing it had to be repointed by hand — they are in
`input_prep/7-source-sink/`, and they are the reason a grep for the module
names is not the same as a grep for the import.

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
                       Scripts only: their products go to
                       output/comparisons/<topic>/, never beside the script.

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
figure_making/old_plot_tests/              input_prep/1-barrier3d-domains/1-extraction/old_extractors/
figure_making/shoreline/old_rate_analysis/
hatteras_ms/old_drafts/                    input_prep/5-scr/CoastSat/old_dsas_comparisons/
hatteras_ms/old_versions/                  input_prep/5-scr/CoastSat/old_time_periods/
sensitivity_analysis/old_guides/           input_prep/6-scr-smooth/old_smoothing/
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

from site_layer.hat_topo_version import topo_dirs, array_name   # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS    # noqa: E402
```

`N` is positional. **Moving a script between folders changes it**, and getting
it wrong yields an `ImportError` — so re-count it whenever a script changes
depth, and keep the comment naming which parent is `scripts/`.

It no longer yields a wrong `PROJECT_ROOT`, which it used to. Every module in
`site_layer/` finds the repo root by searching upward for `pyproject.toml`
(Rule 5 of `ORGANIZATION.md`), so a miscounted `N` now fails loudly at import
instead of quietly resolving paths against the wrong root. Better still, use
the upward search here too and skip the counting:

```python
REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))
```

The same anchor is what makes `from cascade_pipeline.domains import
DomainGeometry` resolve; the package is imported by path, not installed.
