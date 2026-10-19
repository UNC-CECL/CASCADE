> **Version numbers restart at `v1` per product** - `1996-start/v1` (offloaded)
> and `2009-start/v1` are DIFFERENT surfaces. See [LINEAGE.md](LINEAGE.md) for
> what each version was, what it is now, and which DEM built it.

# 1-barrier3d-domains - the arrays Barrier3D starts from

One folder per **start period**. Each period owns its arrays, its dune/topo
versions and its dune-window picks, because the periods start from different
DEMs. Only the 1996 and 2009 starts are run (Hannah, 2026-10-19).

```
1996-start/     the 1996 run   DEM 2009-2014-1996 (1996 ALACE grafted ocean-side)   dune-topo v2
2009-start/     the 2009 run   DEM 2009-2014 (the baseline gap fill)                dune-topo v1
forecast/       from a 2025 DEM, not built yet
shared/         inputs every product uses (below)
figures/        figures that compare the two starts, so belong to neither
README.md  LINEAGE.md
```

**Renamed 2026-10-19.** `1996-start/` was `1984-start/` and `2009-start/` was
`2004-start/`: the same DEMs and the same versions, named for the run that reads
them once the 1984 and 2004 starts were retired. Run metadata and RUN_MANIFESTs
written before then say the old names; `hat_topo_version.canonical_product()`
translates them, and the old `HAT_TOPO_VERSION_1984_START` /
`_2004_START` variables are still honoured.

## shared/

```
buffer/                    padding domains either side of the 90 real ones - not a survey of anything
domain-clips-1m/           per-domain 1 m clips + the 10 m resampled tifs every
                           measurement georeferences against
control-picks/             window sets kept across a version clear, so a deleted
                           extraction stays reproducible from its picks
npy-arrays_2009_unfilled/  the pre-gap-fill arrays; still read by
                           HAT_rasterize_road_to_domains.py as a grid check, so LIVE
domain-geojson/            the 120-domain Pea-Hatteras polygons (extension work)
```

(The digitized dune lines moved to `../2-brie-offset/dunelines/` on 2026-09-15:
they are the island-offset input, not topography.)

## figures/

```
dune_crest_heights/{p10_p90,min_max}/   dune crest height by domain, 1996 vs 2009
```

Written by `scripts/figure_making/pipeline/1-barrier3d-domains/dune_crest_figures.py`,
which also publishes the same files to
`output/figures/3-model-inputs/1-domains/dune_crest_heights/`. Figures of ONE
version stay in that version's `dune-topo/<version>/figures/`.

## Each period holds

```
<period>/
    1-extraction/
        npy-arrays/          domain_<N>.npy   m NAVD88, -10 nodata   <- extractor INPUT
        npy-arrays_survey/   domain_<N>.npy   provenance codes
        picks/               HAT_dune_search_windows_<version>.json   (not regenerable)
        aerial-review/       1996-start only: manual imagery verdicts on the 1996
                             survey holes, survives a re-pick - LIVE
    dune-topo/
        CURRENT          which version to use
        <version>/       topography/  dunes/  figures/  RUN_MANIFEST.txt
    README.md
```

## Offloaded 2026-10-19

Used only by the retired 1984 and 2004 starts, moved (SHA-256 verified) to
`D:\archive\CASCADE_offload\data\hatteras_init\1-barrier3d-domains\`
under their OLD paths:

```
1984-start/2-domain-reconstruction-1984/   the 1984 footprint work (509 MB)
1984-start/dune-topo/v1/                   original pick set; the 09-01 calibration tree and pea1989base ran on it
1984-start/dune-topo/v3/                   v2 + the 1984 footprint behind the road
2004-start/duneline-shift/                 dune-line vs interior-row-0 for the 2004 start
archive_purge_20260826/0903/0907.csv       records of the August/September deletions
```

Their tracked READMEs and CSVs are also in git history. The scripts that
made and read them were deleted the same day and are in git at 2046d67b
(`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/`,
the crest experiment and its scorers, `overwash_vs_footprint.py`, the 1984
dune-line checks in `0-elevation/`). `HAT_rerun_arm.py` can still pin a v1 run
once v1 is copied back.

## Which period uses which DEM

Set in `HATTERAS_PERIODS[start]["topo_product"]` in
`scripts/site_layer/hatteras_site_config.py`, from `YEAR_PRODUCT` in
`hat_topo_version.py`, beside `storm_file`, `island_offset_file` and
`road_setback_file` - it is the same kind of thing, a per-period input.

## Do not hardcode these paths

Resolve through **`scripts/site_layer/hat_topo_version.py`**:

```python
from site_layer.hat_topo_version import topo_dirs, product_for_year, BUFFER_DIR
TOPO, DUNE, VERSION = topo_dirs(product_for_year(1996))   # 1996-start
TOPO, DUNE, VERSION = topo_dirs()                         # 2009-start, the default
```

The rest has a name too: `product_dir()`, `npy_dirs()`, `picks_dir()`,
`dune_topo_root()` per period, and `DOMAIN_ROOT`, `SHARED_DIR`, `BUFFER_DIR`,
`DOMAIN_CLIPS_DIR` / `domain_clip_file(n)`, `CONTROL_PICKS_DIR`,
`UNFILLED_2009_DIR`, `DOMAIN_GEOJSON_DIR`, `DOMAIN_FIGURES_DIR` for what is
shared. Moving `shared/` on 2026-10-19 touched two scripts besides this one.

**The version is resolved, never pinned.** Order:

1. an explicit `override=` argument
2. `HAT_TOPO_VERSION_<PRODUCT>` in the environment, e.g.
   `HAT_TOPO_VERSION_1996_START=v2` - product-scoped, because both products have
   had a `v1` and a global override would resolve to a real but wrong directory
3. `CURRENT`
4. the extractor's `VERSION`, if it is pointed at this product
5. the only version present

Ambiguity raises rather than guesses - pinning by hand is what let the runner sit
on v3 while the setbacks were rebuilt on v4.

**The env rule outranks the rest deliberately** (added 2026-09-02, after it
cost a run). A batch that selected its arm by writing `CURRENT` was ignored.
Two arms of a three-arm experiment silently duplicated the control, exit code 0.
`CURRENT` is a persistent shared DEFAULT; per-run selection needs something that
does not mutate state every other reader sees.

## Filenames carry NO year

`domain_<N>_topography.npy`, `_dune.npy`, `_nodata.npy`. There is no year in the
name and there should not be one - see the long note in
`scripts/site_layer/hat_topo_version.py`. The period lives in the DIRECTORY and
in each run's `RUN_MANIFEST.txt`. Build these names with `array_name()` /
`array_path()` / `domain_arrays()` rather than by hand.

## The arrays are git-ignored; the record of them is not

`.gitignore` ignores this directory's CONTENTS (`1-barrier3d-domains/**`) and
re-includes the written record: every `README.md` and `LINEAGE.md`. The arrays
are rebuilt from `0-elevation/` products, and each run records what made it in
its own `RUN_MANIFEST.txt`.
