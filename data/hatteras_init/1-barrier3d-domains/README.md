> **Version numbers restart at `v1` per product** - `1984-start/v1` and
> `2004-start/v1` are DIFFERENT surfaces. See [LINEAGE.md](LINEAGE.md) for
> what each version was, what it is now, and which DEM built it.

# 1-barrier3d-domains - the arrays Barrier3D starts from

One folder per **start period**. Each period owns its arrays, its dune/topo
versions and its dune-window picks, because the periods start from different
DEMs.

```
1984-start/     from DEM 2009-2014-1996   (1996 ALACE grafted ocean-side)
2004-start/     from DEM 2009-2014        (the baseline gap fill)
forecast/       from a 2025 DEM, not built yet
buffer/         shared padding domains - not a survey of anything
```

Shared, not owned by a period:

```
domain-clips-1m/         per-domain 1 m clips + the 10 m resampled tifs every
                         measurement georeferences against
control-picks/           window sets kept across a version clear, so a deleted
                         extraction stays reproducible from its picks
raw-duneline-geojson/    digitized dune lines, 1967 / 1984 / 1997 / 2004
npy-arrays_2009_unfilled/  legacy pre-gap-fill arrays; nothing live reads them
```

Each period holds:

```
<period>/
    npy-arrays/          domain_<N>.npy   m NAVD88, -10 nodata   <- extractor INPUT
    npy-arrays_survey/   domain_<N>.npy   provenance codes
    dune-topo/
        CURRENT          which version to use, when it is ambiguous
        <version>/       topography/  dunes/  figures/  RUN_MANIFEST.txt
    picks/               HAT_dune_search_windows_<version>.json
    README.md
```

1984-start also carries, and 2004-start may:

```
    duneline-shift/          dune-line vs interior-row-0 measurements
    dune-topo-experiments/   superseded variants, OUTSIDE dune-topo/ so
                             versions() does not list them
    aerial-review/           manual imagery verdicts, survives a re-pick
```

## Which period uses which DEM

Set in `HATTERAS_PERIODS[start]["topo_product"]` in
`scripts/hatteras_site_config.py`, beside `storm_file`, `island_offset_file`
and `road_setback_file` - it is the same kind of thing, a per-period input.

**Before 2026-08-25 there was no such key.** The runner hardcoded one
topography and *both* periods read it, so a 1984 run and a 2004 run started
from the same barrier. They no longer do.

## Do not hardcode these paths

Resolve through **`scripts/hat_topo_version.py`**:

```python
from hat_topo_version import topo_dirs, BUFFER_DIR
TOPO, DUNE, VERSION = topo_dirs("1984-start")
TOPO, DUNE, VERSION = topo_dirs()          # 2004-start, the default
```

`topo_dirs()` with no argument still resolves what it resolved before the
restructure, so the road tree, the groin sweep and the poster script were not
moved by this change.

**The version is resolved, never pinned.** Order:

1. an explicit `override=` argument
2. `HAT_TOPO_VERSION_<PRODUCT>` in the environment, e.g.
   `HAT_TOPO_VERSION_1984_START=v2` - product-scoped, because both products have
   a `v1` and a global override would resolve to a real but wrong directory
3. the extractor's `VERSION`, if it is pointed at this product
4. `CURRENT`
5. the only version present

Ambiguity raises rather than guesses - pinning by hand is what let the runner sit
on v3 while the setbacks were rebuilt on v4.

**The env rule outranks the extractor deliberately** (added 2026-09-02, after it
cost a run). A batch that selected its arm by writing `CURRENT` was ignored,
because rule 3 fires first when the extractor happens to sit on the same product.
Two arms of a three-arm experiment silently duplicated the control, exit code 0.
`CURRENT` is a persistent shared DEFAULT; per-run selection needs something that
does not mutate state every other reader sees.

## Filenames carry NO year

`domain_<N>_topography.npy`, `_dune.npy`, `_nodata.npy`. There is no year in the
name and there should not be one - see the long note in
`scripts/hat_topo_version.py`, which records both the `2009` literal that was
false (neither live product is a 2009 DEM) and the period retag that was tried
and reverted the same day, because the tag reached twelve scripts four different
ways and no single grep could audit it.

The period lives in the DIRECTORY and in each run's `RUN_MANIFEST.txt`. Build
these names with `array_name()` / `array_path()` / `domain_arrays()` rather than
by hand.

> **Stale caller:** `scripts/figure_making/road_check/diagnose_road_drowning.py`
> still builds `domain_{n}_topography_{TOPO_DUNE_INIT_YEAR}.npy` and will not
> find its files.

## Everything here is git-ignored

`.gitignore` line 164 ignores this whole directory. Nothing in it is tracked -
the arrays are rebuilt from `0-elevation/` products, and each run records what
made it in its own `RUN_MANIFEST.txt`.
