# 5-scr — the observed shoreline record, and the rates fitted to it

This is the code. Everything it reads and writes is in the data tree at the
matching path, `data/hatteras_init/5-scr/`, and the two trees are numbered the
same way on purpose: a folder here has a twin there.

| here | there | what it is |
|---|---|---|
| `1-observations/` | `1-observations/` | what was measured or digitized, not fitted by us |
| `2-transect-frame/` | `2-transect-frame/` | which CoastSat transect belongs to which Barrier3D domain |
| `3-rates/` | `3-rates/` | **the model targets** — the rate fits, and one figure per window |
| `4-comparisons/` | `4-comparisons/` | one source against another; answers, not inputs |
| `lib/` | — | modules the scripts share; produces nothing |
| `tools/` | — | checks and indexes; produces no product of its own |
| `template/` | — | a bare, standalone version of this process, to hand to colleagues |

The numbers are the order the work runs. An observation is collected, tied to a
domain, fitted into a rate, and only then compared against something else.
Nothing in `3-rates/` can run before `2-transect-frame/` has built the lookup.

---

## The scripts, in run order

### 1-observations — the record

```
shoreline_inventory/shoreline_inventory.py
    What shoreline data exists across the study area, from all three sources
    (digitized wet-dry lines, NC Coastal Management, CoastSat), how they
    overlap in time, and where the gaps are.
shoreline_patterns/shoreline_trajectory_classification.py
    Classifies each domain's trajectory as stable, eroding or reversing.
shoreline_patterns/shoreline_trajectory_map.py
    The same classification drawn on the island.
```

### 2-transect-frame — the join

```
coastsat_domain_mapping.py
    Spatially joins the CoastSat transects to the 90 domain polygons and
    writes transect_domain_lookup.csv. Point-in-polygon only; nearest-feature
    snapping is deliberately off. Every rate fit downstream reads this file,
    so it runs first and rarely.
```

### 3-rates — the model targets

`rates_figures.py` sits at the top of this folder rather than inside any one
product, because it draws the house-style figure for **all** of them, beside
their tables. Run it after whichever product you rebuilt.

```
coastsat/
    lrr/              the OLS rate fit — what the model is graded against
        coastsat_domain_lrr.py             the fit itself, one window per run
        coastsat_lrr_windows.py            every window on ONE y axis, so a
                                           2 m/yr swing is not drawn as tall
                                           as a 7 m/yr one
        coastsat_lrr_smoothing_windows.py  the 3 / 5 / 10-domain LOESS on one
                                           field
        coastsat_lrr_transect_zoom.py      one window at transect resolution,
                                           over a short reach
    endpoint/         net change between the +/-6-month means at the dune-line
        coastsat_endpoint.py               survey dates, so the shoreline and
                                           the dune line differ like for like
    5yr_bins/         WHEN inside a window the change happened
        coastsat_5yr_bins.py               the table
        coastsat_5yr_bins_figure.py        the figure
    total_change/     the rate as a DISTANCE. --product total is named for the
        coastsat_total_change.py           window it was FITTED on; --product
                                           projected carries a rate onto a
                                           window it was not fitted on.
    extension/        the same fit beyond the 90 surveyed domains, for the
        coastsat_extension_lrr.py          Pea Island extension experiment
duneline/
    duneline_endpoint.py      End dune line minus start line, per 100 m
                              transect and per domain. Endpoint, not LRR: this
                              is net change, and an OLS through the
                              intermediate lines would not be.
rates_figures.py              One figure per window, for every product above.
```

### 4-comparisons — one source against another

```
shoreline_vs_duneline/   does the dune line move with the shoreline?
    coastsat_vs_duneline.py         the base comparison, net change on both
                                    sides; also the module the others import
                                    for the chainage loader
    total_change_vs_duneline.py     total shoreline change against the dune
                                    line's measured change, per window
    net_change_vs_duneline.py       the same over 1996-2024 and its halves
    smoothed_loess7_vs_duneline.py  both curves through a 7-domain LOESS
dsas_vs_coastsat/        the two rate sources against each other
    dsas_vs_coastsat_raw.py         no smoothing, on calendar windows
    dsas_vs_coastsat_datematched.py CoastSat anchored on the survey dates
                                    instead, at +/-30 d and +/-6 months
duneline_positions/      where the 1997 / 2009 / 2023 lines actually sat
    duneline_positions.py           maps, zooms, dune-to-NC-12, beach width
```

### lib and tools

```
lib/coastsat_lrr.py            load_timeseries, compute_lrr, filter_dates.
                               The one OLS every CoastSat script uses, so
                               they cannot drift apart.
lib/scr_paths.py               where the shared modules live (see below)
tools/windows_index.py         regenerates data/.../5-scr/WINDOWS.md from
                               hat_observed_rates.WINDOW_ROLE
tools/coastsat_rates_check.py  consistency checks on the rate outputs:
                               NaN audit, domain means, transect counts
```

### template — for people outside this project

`template/shoreline_rates_template.py` is the whole process — load, window,
fit, screen, group, write — in one standalone file that imports nothing from
this repository. It is for handing to a colleague starting the same work at
another site, and it is deliberately not wired into anything here. See
`template/README.md`, which includes a synthetic dataset with known rates so
the script can be checked before it is trusted.

---

## How a script finds its data, and its siblings

**Data** is resolved through `scripts/site_layer/hat_observed_rates.py`. No
script in this folder types a data path. A window that is not on disk is a
loud error naming the ones that are.

**Sibling modules** are resolved through `lib/scr_paths.py`. Twelve scripts
here import a module from another folder — `rates_figures` for the drawing,
`coastsat_vs_duneline` for the chainage loader, `coastsat_lrr` for the fit.
Each one carries these two lines:

```python
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
```

Importing `scr_paths` is what does the work: it puts every module-bearing
folder on `sys.path`, so a plain `import rates_figures` then resolves from
anywhere in the tree. **To move a module, edit its row in
`scr_paths.MODULE_DIRS` and nothing else.** Running `python lib/scr_paths.py`
checks that table against disk.

Every script finds the repository by searching upward for `pyproject.toml`
(rule 5), so depth is never counted and a file can change folder without
breaking.

---

## How scripts here are named

    <subject>_<product>[_<variant>].py

A **noun phrase, subject first, no prefix and no verb.** The subject is the
data source (`coastsat_`, `duneline_`, `dsas_`, `shoreline_`) or, for a
comparison, both sides in the order the folder reads them
(`total_change_vs_duneline.py` inside `shoreline_vs_duneline/`).

Three rules fall out of that, and one hard constraint:

1. **A name must stand alone, not lean on its folder.** `coastsat/endpoint/`
   holds `coastsat_endpoint.py`, not `endpoint.py`. The repetition looks
   redundant in a path and is not: see the constraint below.
2. **No dates or years in a file name.** A window belongs in an argument, not
   a filename — `net_change_1996_2024.py` had to be renamed the moment it grew
   a second window.
3. **Name the product, not the method.** `smoothed_loess7_vs_duneline.py` says
   what it compares; the `loess7` is the variant, and it survives only because
   it matches the data folder it writes to.

**The constraint: every file name here is also a Python module name, and they
share one flat namespace.** `scr_paths` puts five folders on `sys.path` at
once, so two files called `endpoint.py` in different folders would shadow each
other, and the one that won would depend on import order. File names under
5-scr must therefore be unique across the whole tree, not merely within a
folder.

**No `HAT_` prefix**, unlike much of `input_prep/`. Seven files here are
imported as modules, where `import HAT_rates_figures` reads badly, and the
repo's importable trees — `site_layer/`, `cascade_pipeline/`, `repo_tools/` —
carry no prefix either. Settled 2026-09-22, when 23 of 26 files were already
bare.

`lib/` and `tools/` follow the same rule; `tools/` holds the only names that
describe a job rather than a product (`windows_index.py`,
`coastsat_rates_check.py`), because neither produces one.

---

## 2026-09-22 — this folder was reorganised

Until this date 5-scr was the one stage under `input_prep/` that did not
mirror its data tree. It held twenty-one folders named for their producers —
`CoastSat/`, `CoastSat_timeseries/`, `coastsat_endpoint/`,
`total_change_vs_duneline/` — with four separate `superseded_*` folders, a
folder with an ampersand in its name, and `coastsat_lrr_analysis.py` present
twice as byte-identical copies. `input_prep/README.md` recorded the reason it
had been left alone: the `CoastSat/` scripts imported that module as a
sibling, so splitting them by job meant moving the module first. That is what
`lib/` and `scr_paths.py` are.

The twelve hand-built `sys.path` inserts went at the same time. Each spelled
the folder layout out again, in five different spellings, and a stale one does
not fail where it is written — the insert succeeds against a directory that no
longer exists, and the error surfaces forty lines later naming the module
rather than the path that is wrong. Rule 6 says a location is decided once;
now it is.

**Retired scripts were deleted, not parked.** That departs from rule 4 of
`ORGANIZATION.md` (nothing is deleted for being superseded), and it was a
deliberate call by Hannah on the day, to get the clutter out of the working
tree. Everything was committed before deletion and is recoverable:

```
git log --diff-filter=D --oneline -- scripts/input_prep/5-scr/<old path>
git show <commit>^:scripts/input_prep/5-scr/<old path>
```

What went, and what it had been:

| deleted | what it was |
|---|---|
| `CoastSat/superseded_20260810/` | the two DSAS-vs-CoastSat comparison scripts |
| `superseded_20260918/duneline_lrr.py` | the dune-line LRR, replaced by endpoint |
| `superseded_20260919/duneline_windows/` | dune-line windows, a 3-rates duplicate |
| `superseded_20260921/projected_vs_duneline/` | folded into `total_change_vs_duneline.py` |
| `shoreline_inventory/version_control/` | an older copy of the inventory script |
| `CoastSat/custom_range&dates/` | one-off Buxton press plots |
| `CoastSat/coastsat_domain_lrr_specific_dates.py` | wrote only to a superseded product |
| `CoastSat/coastsat_two_period_comparison.py` | its output had been deleted as stale |
| `CoastSat_timeseries/coastsat_lrr_analysis.py` | byte-identical to the `CoastSat/` copy |

The data products those scripts made are untouched. Where one of them is named
as a producer in a `PROVENANCE.md` or `WHY.md`, that file now carries a note
saying the path has gone and how to read the script again.

`coastsat_domain_lrr_fixed.py` lost its `_fixed` — there is no unfixed sibling
any more — and is now `3-rates/coastsat/lrr/coastsat_domain_lrr.py`.
