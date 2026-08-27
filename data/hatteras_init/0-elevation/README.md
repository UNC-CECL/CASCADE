# 0-elevation - the DEMs Barrier3D starts from

One folder per **product**. A product is a composition of surveys clipped to
the 90 GIS domains, carried through the same two stages, with its own figures.

```
2009-2014/          the baseline DEM      2009 USACE + 2014 NOAA Post-Sandy
2009-2014-1996/     the 1984-start DEM    the above + 1996 ALACE, no boundary
superseded/         earlier attempts, same shape, not for use
source-selection/   island-wide DEM-candidate scoring - belongs to no product
FIGURES.md          figure design decisions, shared by every product
```

Each product holds:

```
<product>/
    1-gapfill-1m/      clip_domain_<N>_filled.tif   + _survey.tif  + audit CSV
    2-resampled-10m/   resampled_domain_<N>_*.tif   + resample_audit.csv
    figures/
    README.md          what this DEM is, and what it is for
```

## Naming

Folders are named for their **composition**, not for a hindcast period and not
for the fill source alone.

* Not the period: `2009-2014` currently serves **both** the 1984 and the 2004
  runs, so calling it `2004-start` would assert something that is not true.
* Not the fill source: the old name `2014_NOAA_PostSandy` said what was *added*
  but not what the product *contained*, which is the question a reader
  actually has.

## Do not hardcode these paths

Resolve them through **`scripts/hat_elevation_products.py`**:

```python
import sys; sys.path.insert(0, str(REPO / "scripts"))
from hat_elevation_products import product

p = product("2009-2014")
p.gapfill_1m, p.resampled_10m, p.figures, p.audit_1m
```

It raises with the list of known products if the name is wrong, and again if
the name is known but the directory is not there. That guard exists because
`HAT_road_elevation.py` built its path by joining strings, and stopped
resolving the moment 2008 moved under `superseded/`.

## The .tif files are git-ignored

`.gitignore` line 148 is `*.tif`, so the 720 rasters here (~370 MB) live on
your machine only. **The audit CSVs and READMEs are tracked** - those carry the
per-domain numbers, so a reviewer can check the work without the bulk. Rebuild
the rasters from the D: drive sources with the chain in
[`scripts/input_prep/0-elevation/README.md`](../../../scripts/input_prep/0-elevation/README.md).

## History

Before 2026-08-25 this was stage-first: `1-gapfill-1m/<TAG>/` and
`2-resampled-10m/<TAG>/`, with every product's figures pooled in one
`figures/`. That worked while there was one product. With two it scattered each
DEM across three top-level directories and left the figures ambiguous, so the
tree was inverted to product-first.
