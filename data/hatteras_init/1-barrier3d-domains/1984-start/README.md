# 1984-start

Barrier3D domains for the **1984–2004** hindcast period.

| | |
|---|---|
| DEM | `0-elevation/2009-2014-1996/` — the baseline plus 1996 ALACE grafted wherever ALACE has data (no road boundary; rebuilt 2026-08-26) |
| npy arrays | 90 domains, exported 2026-08-26 from the no-boundary DEM |
| **what actually loads** | **`v3`** — see the warning below; it is *not* what `CURRENT` says |
| picks | `picks/HAT_dune_search_windows_v3.json` |

## Layout

```
npy-arrays/            domain_<N>.npy   m NAVD88, -10 nodata   <- extractor INPUT
npy-arrays_survey/     domain_<N>.npy   provenance codes (0/1996/2009/2014)
picks/                 the dune search windows that define each version
dune-topo/             v1, v3, v4, v5 + CURRENT + a version index README
dune-topo-experiments/ emptied 2026-09-03; two comparison figures + a record
row-insert-scope/      EVERYTHING for the seaward-row insert:
                         duneline-shift/  the measurement of N (moved here)
                         the scope report, the fill comparison, five figures
aerial-review/         58 holes of 1996 imagery review, keyed on (domain, profile)
```

## Which topography loads

**`CURRENT` says `v5`. `resolve_version()` returns `v3`.** The extractor's
`VERSION` literal outranks the file, so `CURRENT` is currently inert. The full
explanation, the version-by-version table, and the two steps needed to actually
run on `v5` are in **`dune-topo/README.md`** — read that before pointing
anything at a version.

Short form:

| | what it is |
|---|---|
| `v1` | extraction on the OLD pick set. Superseded; kept for the `pea1989base` run arm |
| **`v3`** | extraction on the 2026-09-02 re-pick. **The base everything resolves to** |
| `v4` | `v3` + rows at the two relocation blocks (8 domains) |
| `v5` | `v3` + rows at all 38 domains the measurement selects. **The one to take forward** |

`v2` was deleted 2026-09-03 — see `dune-topo/README.md`.

## duneline-shift moved into row-insert-scope

It is a measurement and its interpretation, and they now sit together.
`row-insert-scope/duneline-shift/` holds `duneline_retreat_1984_1997.csv`,
which is *where N comes from* — `HAT_insert_seaward_rows.py` builds v4 and v5
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
