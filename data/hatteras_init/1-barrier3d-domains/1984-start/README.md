# 1984-start

Barrier3D domains for the **1984–2004** hindcast period.

| | |
|---|---|
| DEM | `0-elevation/2009-2014-1996/` — the baseline plus 1996 ALACE grafted wherever ALACE has data (no road boundary; rebuilt 2026-08-26) |
| npy arrays | 90 domains, exported 2026-08-26 from the no-boundary DEM |
| **what loads** | **`v2`** — `CURRENT` says so, and since 2026-09-04 `CURRENT` decides |
| picks | `picks/HAT_dune_search_windows_v2.json` (`_v1.json` for `v1`) |

## Layout

```
1-extraction/                  THE FIRST HALF: what the extractor reads and records
    npy-arrays/                domain_<N>.npy   m NAVD88, -10 nodata   <- extractor INPUT
    npy-arrays_survey/         domain_<N>.npy   provenance codes (0/1996/2009/2014)
    picks/                     the dune search windows that define each version
    aerial-review/             58 holes of 1996 imagery review, keyed on (domain, profile)
    dune-topo-experiments/     emptied 2026-09-03; two comparison figures + a record
2-domain-reconstruction-1984/  THE SECOND HALF: v2 -> v3, the 1984 domains reconstructed,
                               in six steps (1-measurement ... 6-result) + figures/;
                               was row-insert-scope/ until 2026-09-09
dune-topo/                     THE PRODUCT, written by both halves: v1, v2 (extraction),
                               v3 (reconstruction) + CURRENT. Stays at the root: what loads.
```

Two halves, mirroring `scripts/input_prep/1-barrier3d-domains/1-extraction/`
and `2-domain-reconstruction-1984/` (2026-09-09, Hannah). Paths resolve through
`hat_topo_version` (`extraction_dir`, `npy_dirs`, `picks_dir`,
`insert_scope_dir`, `insert_scope_step`); the extractor derives its own from
the product folder. 2004-start has the same `1-extraction/` and no
reconstruction; its `duneline-shift/` stays at its root.

## Which topography loads

**`CURRENT` says `v2` and, since 2026-09-04, `CURRENT` decides** (until then the
extractor's `VERSION` literal outranked it and the file was inert). The full
explanation, the resolution order, and the two steps that changing the default
takes are in **`dune-topo/README.md`**; the deleted layers are described in
`2-domain-reconstruction-1984/DUNE_TOPO_VERSION_GUIDE.md` — read those before pointing
anything at a version.

Short form:

| | what it is |
|---|---|
| `v1` | extraction on the ORIGINAL pick set (2026-08-27). The 2026-09-01 calibration tree and the `pea1989base` arm ran on it |
| **`v2`** | extraction on the 2026-09-02 re-pick, NC-12 visible. **What `CURRENT` names and the road tree measures against** |

Both are unmodified extractions: no rows added, no cell edited. **The layers
`v3`–`v8`** (`v2` + rows inserted seaward where the 1984 dune line stood
seaward of the 1996 one; one per scope/fill rule) **were deleted 2026-09-07**
together with every run made on modified topography — Hannah's decision to
keep only unmodified topography. `dune-topo/README.md` has the list; the
guide has what each was. The pre-re-pick `v2` had gone on 2026-09-03.

## duneline-shift moved into 2-domain-reconstruction-1984

It is a measurement and its interpretation, and they now sit together.
`2-domain-reconstruction-1984/1-measurement/duneline-shift/` holds `duneline_retreat_1984_1997.csv`,
which is *where N comes from* — `HAT_insert_seaward_rows.py` builds the layers
from it. Delete it and no future insert version can be built. It is 301 KB.

The path is resolved by **`hat_topo_version.duneline_shift_dir(product)`**, not
built by hand. 2004-start keeps the plain `2004-start/duneline-shift/` layout;
the asymmetry is contained in that one function.

## History

Cleared 2026-08-27 back to a blank slate: the earlier `v1` (2026-08-26) and its
bridged `v2` were deleted along with their window sets, so numbering restarted.
Sizes and reasons are in `../archive_purge_20260826.csv`; the lineage entry is
in `../LINEAGE.md`. The window sets from before the clear survive at
`../control-picks/HAT_dune_search_windows_1984-start_v{1,2}.json`, and the
aerial review survived in `aerial-review/` because a cross-shore re-pick cannot
move a `(domain, profile)` key.

Tidied again 2026-09-03: `dune-topo/v2` and three unreferenced experiment
variants removed (19.3 MB), logged in `../archive_purge_20260903.csv`.

Reduced to the two extractions 2026-09-07: `dune-topo/v3`–`v8` (29 MB) and
every run on inserted or edited topography (8.7 GB under `output/`) removed;
the calibBE `road_bdm_groin` calibration-tree run re-run on `v2`. Logged in
`../archive_purge_20260907.csv`; lineage entry in `../LINEAGE.md`.
