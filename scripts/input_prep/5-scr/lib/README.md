# lib — what the 5-scr scripts share

Two modules. Neither produces a product, and neither has a data twin, which is
why this folder is not numbered.

```
coastsat_lrr.py   The OLS. load_timeseries, filter_dates, compute_lrr,
                  _empty_lrr — imported by the LRR fit, the 5-year bins, the
                  extension fit and the raw DSAS comparison. One copy so the
                  four cannot drift apart; there were two byte-identical
                  copies of it until 2026-09-22.
                  (Was coastsat_lrr_analysis.py, under CoastSat/.)

scr_paths.py      Where the shared modules live. Import it and every
                  module-bearing folder in 5-scr goes onto sys.path.
```

## Using scr_paths

A script that imports a sibling module carries these two lines, and no others:

```python
sys.path.insert(0, str(_REPO / "scripts" / "input_prep" / "5-scr" / "lib"))
import scr_paths  # noqa: E402,F401  (5-scr sibling modules onto sys.path)
```

After that, `import rates_figures` or `from coastsat_lrr import compute_lrr`
resolves from anywhere in the tree.

**To move a module, edit its row in `MODULE_DIRS` and nothing else.**
`python lib/scr_paths.py` checks every row against disk and names the ones that
are wrong — which is the thing the twelve hand-built `sys.path` inserts this
replaced could never do. A stale insert does not fail where it is written: it
succeeds against a directory that no longer exists, and the error arrives forty
lines later naming the module instead of the path.

Six modules are registered. `coastsat_lrr` and `scr_paths` live here; the other
four live beside the product they belong to, because they are working scripts
that other scripts happen to import:

| module | lives in | because |
|---|---|---|
| `rates_figures` | `3-rates/` | it draws every 3-rates product |
| `coastsat_lrr_windows` | `3-rates/coastsat/lrr/` | it is the multi-window LRR figure |
| `duneline_endpoint` | `3-rates/duneline/` | it is the dune-line product |
| `coastsat_vs_duneline` | `4-comparisons/shoreline_vs_duneline/` | it is the base comparison |
| `total_change_vs_duneline` | `4-comparisons/shoreline_vs_duneline/` | same folder, same question |

Pulling those four in here would file them away from the product they make,
which is the larger harm. Only code that exists *solely* to be imported
belongs in `lib/`.
