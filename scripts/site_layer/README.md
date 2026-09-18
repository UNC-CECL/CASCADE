# `site_layer/` — what and where Hatteras is

The seam between the site-agnostic `cascade_pipeline/` package and the task
folders, which all consume Hatteras' own paths, place names and house style.
Six modules, each answering one question, each answering it in one place.

| Module | Answers | Imported by |
|---|---|---|
| `hatteras_site_config.py` | *What is this place* — domain geometry, town spans, periods, BE presets, road events, nourishment projects. | 87 |
| `hat_figure_style.py` | *What does a Hatteras figure look like* — one type scale, one palette, one elevation ramp, one way to letter a panel. | 102 |
| `hat_topo_version.py` | *Which* Barrier3D domains — which product (`1984-start` / `2004-start` / `forecast`) and which version inside it. | 77 |
| `hat_observed_rates.py` | *Where is the observed shoreline* — the CoastSat chainage, the per-window rate fits the model is graded against, the transect-to-domain lookup. | 22 |
| `hat_elevation_products.py` | *Which* elevation product and stage — `2009-2014` or `2009-2014-1996`, gapfilled 1 m or resampled 10 m. | 13 |
| `hat_extension_domains.py` | *Which alongshore reach* a run models, and the 500 m bins that number the coast beyond the 90 surveyed domains. | 11 |

## Why the indirection exists

**A hand-built path fails silently.** Four road scripts once hardcoded
`2009-dune-topo/2009_v3`; when the dune windows were re-picked into `v4` they
kept reading `v3` interiors while consuming `v4` setbacks, and 18 domains —
two of the three managed roadways among them — had their drown verdicts
computed on the wrong grid, with nothing raised. The same shape of bug hit
`HAT_road_elevation.py` when a fill source moved under `superseded/`.

Each location is resolved **once**, here, and a name that is not on disk is an
immediate, loud error listing what is. So never rebuild one of these paths by
hand: call `topo_dirs()`, `domain_arrays()`, `array_path()` or `product()` and
let it raise.

## How to import from it

```python
import sys
from pathlib import Path

REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
sys.path.insert(0, str(REPO / "scripts"))

from site_layer.hat_topo_version import topo_dirs, array_name   # noqa: E402
from site_layer.hatteras_site_config import HATTERAS_DOMAINS    # noqa: E402
```

Import the **module**, not the package. `__init__.py` deliberately re-exports
nothing: pulling the six in would make every consumer pay for all of them —
`hatteras_site_config` alone reads several CSVs at import — and would turn the
one real cycle in here (`hatteras_site_config` imports `hat_topo_version`;
`hat_figure_style` lazily imports `hatteras_site_config`) into an import-order
problem.

**Depth no longer matters.** Every module finds the repo root by searching
upward for `pyproject.toml`, not by counting `parents[N]`, which is Rule 5 of
`ORGANIZATION.md`. That is what made moving them off the `scripts/` root on
2026-09-18 a safe operation — the older `parents[1]` spelling would have
silently resolved to `scripts/` one level down.

## The two that also run

`hat_figure_style.py` writes the style sheet (`python scripts/site_layer/hat_figure_style.py`
→ `data/hatteras_init/9-figures/`) and `hat_observed_rates.py` prints a
locations diagnostic. Running a file directly puts *its own folder* on
`sys.path` rather than `scripts/`, so both open with a
`if __package__ in (None, ""):` guard that adds `scripts/` back. Importing
them normally never takes that branch.

## Not for another site

`cascade_pipeline/` ships no site content by contract — a different study site
writes its own sibling of this package and never touches it. One leak already
exists (`cascade_pipeline/hindcast.py` imports `hat_topo_version`); folding
site content into the package would make it permanently Hatteras-only.
