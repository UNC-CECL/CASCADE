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

The numbers are the order the work runs. An observation is collected, tied to a
domain, fitted into a rate, and only then compared against something else.
Nothing in `3-rates/` can run before `2-transect-frame/` has built the lookup.

---

## The scripts, in run order

### 1-observations — the record

```
shoreline_inventory/HAT_shoreline_inventory.py
    What shoreline data exists across the study area, from all three sources
    (digitized wet-dry lines, NC Coastal Management, CoastSat), how they
    overlap in time, and where the gaps are.
shoreline_patterns/HAT_shoreline_trajectory_classification.py
    Classifies each domain's trajectory as stable, eroding or reversing.
shoreline_patterns/HAT_trajectory_map.py
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
    lrr/            the OLS rate fit — what the model is graded against
        coastsat_domain_lrr.py       the fit itself, one window per run
        coastsat_lrr_windows.py      every window on ONE y axis, so a 2 m/yr
                                     swing is not drawn as tall as a 7 m/yr one
        lrr_smoothing_windows.py     the 3 / 5 / 10-domain LOESS on one field
        lrr_transect_zoom.py         one window at transect resolution, over a
                                     short reach
    endpoint/       net change between the +/-6-month means at the dune-line
        coastsat_endpoint.py         survey dates, so the shoreline and the
                                     dune line difference like for like
    5yr_bins/       WHEN inside a window the change happened
        coastsat_lrr_5year_bins.py   the table
        coastsat_5yr_bins_figure.py  the figure
    total_change/   the rate as a DISTANCE. --product total is named for the
        coastsat_total_change.py     window it was FITTED on; --product
                                     projected carries a rate onto a window it
                                     was not fitted on. Writes both.
    extension/      the same fit beyond the 90 surveyed domains, for the
        coastsat_extension_lrr.py    Pea Island extension experiment
duneline/
    duneline_endpoint.py    End dune line minus start line, per 100 m transect
                            and per domain. Endpoint, not LRR: this is net
                            change, and an OLS through the intermediate lines
                            would not be.
rates_figures.py            One figure per window, for every product above.
```

### 4-comparisons — one source against another

```
shoreline_vs_duneline/   does the dune line move with the shoreline?
    duneline_vs_coastsat.py        the base comparison, net change on both
                                   sides; also the module the others import
                                   for the chainage loader
    total_change_vs_duneline.py    total shoreline change against the dune
                                   line's measured change, per window
    net_change_1996_2024.py        the same over 1996-2024 and its two halves
    smoothed_loess7.py             both curves through a 7-domain LOESS
dsas_vs_coastsat/        the two rate sources against each other
    dsas_vs_coastsat_raw.py          no smoothing, on calendar windows
    dsas_vs_coastsat_datematched.py  CoastSat anchored on the survey dates
                                     instead, at +/-30 d and +/-6 months
duneline_positions/      where the 1997 / 2009 / 2023 lines actually sat
    duneline_positions.py            maps, zooms, dune-to-NC-12, beach width
```

### lib and tools

```
lib/coastsat_lrr.py     load_timeseries, compute_lrr, filter_dates. The one
                        OLS every CoastSat script uses, so they cannot drift.
lib/scr_paths.py        where the shared modules live (see below)
tools/write_windows_md.py         Regenerates data/.../5-scr/WINDOWS.md from
                                  hat_observed_rates.WINDOW_ROLE.
tools/coastsat_verify_pipeline.py Consistency checks on the rate outputs:
                                  NaN audit, domain means, transect counts.
```

---

## How a script finds its data, and its siblings

**Data** is resolved through `scripts/site_layer/hat_observed_rates.py`. No
script in this folder types a data path. A window that is not on disk is a
loud error naming the ones that are.

**Sibling modules** are resolved through `lib/scr_paths.py`. Twelve scripts
here import a module from another folder — `rates_figures` for the drawing,
`duneline_vs_coastsat` for the chainage loader, `coastsat_lrr` for the fit.
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
