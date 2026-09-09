# 1984-start dune-topo — version guide

**What each version under `1984-start/dune-topo/` is, how it was made, what is
measured in it and what is not, and what has been run on it.** Written
2026-09-04, the day the versions were renumbered. This is the guide;
`../dune-topo/README.md` is the index and points here.

> **2026-09-07 — `v3`–`v8` no longer exist.** Only the two extractions `v1`
> and `v2` are kept (Hannah's decision: unmodified topography only). Every
> run listed below as made on a layer was deleted with it, as were
> `output/experiments/row_insert_set/` and the relocation-comparison
> reports; the calibration-tree `HAT_1984_2004_calibBE_road_bdm_groin` was
> re-run on `v2`. This guide is kept as the description of what the layers
> were, how they were built (the recipe in *What every layer shares* still
> works against `v2`) and what they showed. `../../archive_purge_20260907.csv`.

## The numbering

A plain sequence, one number per version, in build order:

| version | what | built |
|---|---|---|
| `v1` | extraction, **the original pick set** | 2026-08-27 |
| `v2` | extraction, re-pick with NC-12 visible — **the base every layer is built on, and what `CURRENT` says** | 2026-09-02 |
| `v3` | v2 + rows at the ten relocation-block domains, measured + floor | 2026-09-02/03 |
| `v4` | v2 + rows at all 38 measured domains, measured + floor — the rule that shipped 09-03 | 2026-09-03 |
| `v5` | same rows, measured + median — the fill decided 09-04 | 2026-09-04 |
| `v6` | same rows, flat backdune platform — reference | 2026-09-04 |
| `v7` | same rows, matched backdune, crest kept | 2026-09-04 |
| `v8` | same rows, matched backdune, crest skipped | 2026-09-04 |

`v1` and `v2` are **extractions** (all 90 domains from the DEM). `v3`–`v8` are
**layers**: rows prepended to `v2`, bit-identical to it behind the added rows,
dune arrays copied unchanged, no cell of the extracted interior modified. That
was a constraint set in the 2026-09-04 interview and every layer obeys it.

### Rename map — read this before trusting a version name in old output

The layers were numbered differently when built, and for about an hour on
2026-09-04 carried lineage names. **Run metadata, `run_index.csv` rows, dated
reports, figure file names and arm tags written before the renumber use the
as-built names**, and several of those numbers now mean something else.

| as built | interim (09-04) | **final** | what it is |
|---|---|---|---|
| `v1` | `v1` | **`v1`** | original extraction |
| `v3` | `v3` | **`v2`** | re-pick base |
| `v4` | `v3-blocks-floor` | **`v3`** | blocks, measured + floor |
| `v5` | `v3-island-floor` | **`v4`** | island, measured + floor |
| `v6` | `v3-island-median` | **`v5`** | island, measured + median |
| `v7` / `test_backdune` | `v3-island-platform` | **`v6`** | island, flat platform |
| `v8` | `v3-island-matchedcrest` | **`v7`** | island, matched, crest kept |
| `v9` | `v3-island-matchednocrest` | **`v8`** | island, matched, crest skipped |

So: a run whose metadata says `v3` ran on today's `v2`; `v5` means today's
`v4`; `v6` means today's `v5`. The arm tags `islandv5` and `blocksv4` in `output/raw_runs/` were named for
as-built versions and ran on today's `v4` and `v3`. (The four-arm `*noreloc`
test of 2026-09-04 was deleted once the set below had run.) The
calibration-tree `HAT_1984_2004_calibBE_road_bdm_groin` of 2026-09-04 says `v6`
in its metadata and ran on today's `v5`. Each renamed folder's
`RUN_MANIFEST.txt` keeps its build-time header with the rename notes appended.
`v2` (as built: `v1` + block-scope rows on the old picks) was deleted
2026-09-03, before any of this.

**The picture of all this**: `figures/superseded-layers/HAT_insert_explainer_grid_GIS85.png` (plan-view grids) and `figures/superseded-layers/HAT_insert_explainer_GIS85.png` (one cross-shore line)
(`scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/4-fill/HAT_plot_insert_explainer.py --domain N`
redraws it for any domain).

## What every layer shares

| | |
|---|---|
| **Source DEM** | 2009-2014-1996 mosaic: 1996 ALACE beach and foredune grafted at the dune toe onto the 2009+2014 backdune (`0-elevation/2009-2014-1996-duneline/`) |
| **What the interior is made of** | measured, from the per-cell survey year in `npy-arrays_survey/` mapped through the extractor's orient/shear (2026-09-04): interior row 0 is 1996 ALACE in 88 of 90 domains; the rows behind it stay mostly 1996 to a median depth of 20 rows (IQR 15–24; GIS 85: 15, all of rows 0–12), then 2009. The pick decides only where row 0 falls; nothing landward of it is changed by any version. The inserted rows land on cells the DEM holds as 1996 too, but there they are the 1996 beach and dune face, a later landform than the 1984 ground they stand for |
| **Pick set** | the 2026-09-02 re-pick with NC-12 drawn on the picker: `1984-start/1-extraction/picks/HAT_dune_search_windows_v2.json` (v1 has its own) |
| **N, the number of rows** | `round((line_1997 − line_1984) / 10 m)`, floored at 0, per domain, from `2-domain-reconstruction-1984/1-measurement/duneline-shift/duneline_retreat_1984_1997.csv`. Measured, not tuned. Island scope: 38 domains, 98 rows, largest 7 (GIS 80) |
| **Where the rows go** | prepended at the seaward end of the interior array, so interior row 0 moves N cells seaward and the road's setback becomes `(road − row0) + 10·N` |
| **Setback CSV** | each layer carries `RoadSetback_1984_dunestart.csv` matched to its arrays. `v4`–`v8` are identical (same N). GIS 85: −15 → +45 m; GIS 86: −10 → +20 m |
| **Known caveats** | the dune lines are 1984 and 1997 but the surface is 1996 (N overstates 1984→1996 retreat by close to a year; moves N at GIS 12, 84, 85 by one); progradation is not removal (N ≥ 0); the dune arrays assume the 1984 crest height equalled the 1996 one; no survey covers the added ground |
| **Build** | `python scripts/input_prep/1-barrier3d-domains/2-domain-reconstruction-1984/4-fill/HAT_insert_seaward_rows.py --variant pad --fill <rule> --shift-source date --n-rule measured --domains measured --dst-version vN` (block scope: `--domains 9,10,...`) |
| **Audit** | `HAT_seaward_row_insert_audit.csv` in each folder: N, setback before/after, land rows, mean land elevation, cells inserted, cells from the DEM |

## The versions

### `v1` — extraction, the original pick set

| | |
|---|---|
| built | 2026-08-27, `HAT_dune_topo_extractor.py` |
| what | all 90 domains from the 2009-2014-1996 DEM with the pre-re-pick dune windows (`picks/HAT_dune_search_windows_v1.json`) |
| status | superseded by `v2` for new work; kept as the original, and the `pea1989*` experiment arms ran on it |
| setback CSV | `v1/RoadSetback_1984_dunestart.csv`, copied 2026-09-04 from `4-mgmt-forcing/.../dunestart_offset_ARCHIVE_1984start_v1/` — the road tree's measurement on these interiors, frozen before the re-pick (differs from v2's at 27 domains) |
| runs | the published relocation comparison of 2026-08-31 and every 1984 calibration run before 2026-09-02; set arm **`original`** (reference, 2026-09-04) |

### `v2` — extraction, re-pick with the road visible. **The base.**

| | |
|---|---|
| built | 2026-09-02 (as `v3`); the extractor's `VERSION` literal is `"v2"`, still what it writes |
| what | all 90 domains; every dune window re-picked with NC-12 drawn on the picker, so a crest argmax could not lock onto the road embankment. 39 domains got a taller crest; GIS 85 +1.52 m |
| rows added | none. GIS 85 setback −15 m → 0 in the model CSV, GIS 86 −10 → 0 |
| setback CSV | `v2/RoadSetback_1984_dunestart.csv`, the road tree's measurement on these interiors, saved 2026-09-04 (also the live forcing-tree file) |
| role | the extraction every layer is built on; the `none` control arm of the fill set |
| default | **what `CURRENT` says today** |
| runs | set arm **`none`** (2026-09-04); every calibration-tree 1984 run made between 2026-09-02 and 09-04 |

### `v3` — rows at the ten relocation-block domains, measured + floor

| | |
|---|---|
| built | 2026-09-02/03 (as `v4`), `--fill measured --domains 9,10,11,12,13,14,84,85,86,87` |
| scope | block; rows at 8 of the 10 (GIS 9 and 14 round to 0) |
| fill | **measured + floor**: a dry DEM cell at the row's coordinates is kept where it is at or above the backdune platform (interior rows 1–3 median) and floored to the platform otherwise; water cells get the platform |
| status | superseded by `v4` (identical N at all ten block domains; `v3` is a strict subset). Kept for the `blocksv4` arm and the scope-change figure. **Not in the fill set** |

### `v4` — the shipped rule, island scope

| | |
|---|---|
| built | 2026-09-03 (as `v5`), `--fill measured --domains measured` |
| fill | **measured + floor**, as above |
| measured share | 3732 of 4900 inserted cells are DEM cells at their own coordinates (76%); min 42% (GIS 3, 79), 15 domains at 100% |
| what it asserts | 1996 is a lower bound on 1984 — and then raises 44% of the GIS 85 block to the platform, which asserts more than a lower bound |
| GIS 85 block | 1.65 1.65 1.70 1.92 3.17 4.96 m (row medians, seaward → landward); road rows 4–5 on 3.17 / 4.96 m |
| runs | `islandv5` arm (2026-09-03); set arm **`measured-floor`** |

### `v5` — measured + median, the decided rule

| | |
|---|---|
| built | 2026-09-04 (as `v6`), `--fill median --domains measured` |
| fill | **measured + median**: every dry DEM cell at the row's coordinates kept as measured; only cells at or below MHW get a value, the median of the block's own dry cells. One guard, one constant, no measurement raised |
| measured share | 4765 of 4900 (97%); min 79% (GIS 79), 32 domains round to 100% |
| what it asserts | 1996 is a lower bound on 1984, and nothing further |
| GIS 85 block | 1.35 0.78 1.20 1.84 3.17 4.96 m; 28 water cells took 1.74 m; road rows 4–5 on 3.17 / 4.96 m |
| known cost | the road sits inside the measured 1996 crest at GIS 85 and 86. `bulldoze()` scrapes it in year 1 and splits it over the dune: +3.3 m and +2.5 m per dune cell (32,700 and 25,300 m³), dune to ~7 m at GIS 85. Accepted and recorded, not corrected |
| vs `v4` | setbacks and land rows identical in every domain; mean land elevation differs by 0.001 m median, −0.043 m at GIS 85; indistinguishable in the model except one year-1 dune height at GIS 81 |
| runs | the calibration-tree `HAT_1984_2004_calibBE_road_bdm_groin` of 2026-09-04 (metadata says `v6`); set arm **`median`** |
| status | decided in the 2026-09-04 interview; was the default for a few hours; default reverted to `v2` while the fill set is compared |

### `v6` — flat backdune platform, reference

| | |
|---|---|
| built | 2026-09-04 (as `v7`, earlier as `test_backdune`), `--fill backdune --domains measured` |
| fill | **flat platform**: every inserted cell at the per-column median of interior rows 1–3. Entirely fabricated (0 of 4900 from the DEM), flat across the domain by construction |
| GIS 85 block | 1.65 m in every row; road on 1.65 m; year-1 scrape 0.86 m per dune cell |
| role | dropped as a candidate 2026-09-03; in the set as the cleanest "road on backdune" **reference**, because it removes the crest cost with no other change |
| runs | set arm **`platform`** |

### `v7` — matched backdune, crest kept

| | |
|---|---|
| built | 2026-09-04 (as `v8`), `--fill matched-crest --domains measured` (rule added that day) |
| fill | **profile copy**: block row k = interior row k, k = 0 … N−1, per column. Panel (b) of the fill figure exactly. Because v2 row 0 IS the 1996 crest, the crest appears at the new seaward edge and again N cells landward |
| what it asserts | the 1984 backdune looked like the present one, further seaward — form preserved, position shifted — with the crest duplicated |
| measured share | 0 of 4900 at their own coordinates; every value a real cell of the same domain |
| GIS 85 block | 3.15 2.22 1.63 1.51 1.38 1.21 m; road on 1.38 / 1.21 m |
| runs | set arm **`matched-crest`** (2026-09-04) |

### `v8` — matched backdune, crest skipped

| | |
|---|---|
| built | 2026-09-04 (as `v9`), `--fill matched-nocrest --domains measured` |
| fill | **profile copy from row 1**: block row k = interior row k+1. Same assertion as `v7` without the second crest: one dune array, backdune, then the measured 1996 crest as a relict ridge |
| measured share | 0 of 4900 at their own coordinates |
| GIS 85 block | 2.22 1.63 1.51 1.38 1.21 1.12 m; road on 1.21 / 1.12 m; the lowest mean land elevation of the five island layers |
| runs | set arm **`matched-nocrest`** (2026-09-04) |

Both matched builds have one cell at GIS 81 (block row 1, column 3) that is a
−3.0 m water sentinel inside the near-dune interior, re-floored at the sentinel
and differing from its source by float rounding only.

## The fill comparison set

Seven runs — `v1`, `v2` and the five island layers — all 1984–2004 calibBE full
management, groin on (M 60, f 0.6), prescribed relocations off: the
calibration-tree run's settings exactly.

| arm | version |
|---|---|
| `original` | `v1` (reference: the original picks with their own setbacks; added for the relocation comparison) |
| `none` | `v2` |
| `measured-floor` | `v4` |
| `median` | `v5` |
| `platform` | `v6` |
| `matched-crest` | `v7` |
| `matched-nocrest` | `v8` |

Driver `scripts/hatteras_ms/HAT_run_row_insert_set.py` → runs under
`output/raw_runs/row-insert/<arm>/1984_2004/calibBE/` (a two-level arm, allowed
by `run_registry.arm_component` since 2026-09-04). Comparison
`HAT_plot_row_insert_set.py` → `output/experiments/row_insert_set/`.
**Run 2026-09-04**; results in `output/experiments/row_insert_set/` (README there).
Each arm also has a relocation-ON partner (run name with the `reloc` token, same
folder), the arm B of `HAT_relocation_comparison.py`; those reports are under
`output/comparisons/relocation_1984_2004/row-insert/<arm>/`.

## Which version a script reads

`scripts/hat_topo_version.py`, in order: explicit `override=`;
`HAT_TOPO_VERSION_1984_START`; the `CURRENT` file in `dune-topo/`; the
extractor's `VERSION` literal (what it *writes*); the only version present.
`CURRENT` outranks the extractor since 2026-09-04 so a layer can be adopted
without editing the extractor to a name it would overwrite. **`CURRENT` says
`v2`.** Adopting a layer is two steps: write its name into `CURRENT`, and copy
its `RoadSetback_1984_dunestart.csv` over
`4-mgmt-forcing/road_offset/dunestart_offset/1984/RoadSetback_1984_dunestart.csv`,
which `hatteras_site_config.py` hardcodes. A version's arrays paired with
another version's CSV restores the off-by-N error the layers exist to remove.

## Adding a version

1. Add the rule to `fabricate_rows()` in `HAT_insert_seaward_rows.py` and to
   the `--fill` choices; give it one token.
2. Build with `--dst-version v9` (the next number).
3. Verify: `v9[n:] == v2` in all 90 domains, dune arrays equal, setback CSV
   identical to `v4`'s (island scope).
4. Write the folder README (what it asserts, GIS 85 block, measured share),
   add it here and to `../dune-topo/README.md`, and a line to `../../LINEAGE.md`.
5. Add it as an arm to `HAT_run_row_insert_set.py` and `HAT_plot_row_insert_set.py`.
