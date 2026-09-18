# road_offset — what each file is, and what wrote it

Where NC-12 sits per Barrier3D domain. The **method** documentation lives with
the code, at `scripts/input_prep/4-mgmt-forcings/road_offset/README.md` — that
is where the two methods are explained and compared. This file is the inventory:
for every file here, what produced it, whether the model reads it, and whether
losing it costs anything.

Two measured vintages, `1984` and `2004`, plus two derived from them, `1996`
and `2010`. Each measurement is made against **its own period's extraction** —
`1984-start` and `2004-start` — because the two are different islands, and on
**its own period's line** — the 1978 and 2008 exports. Both pairings are
defined once, in `scripts/site_layer/hat_topo_version.py` (`YEAR_PRODUCT`,
`ROAD_LINE_FOR_YEAR`).

**Two kinds of integer in this tree (since 2026-09-15).** A folder named
`1984`, `1996`, `2004` or `2010` is a PERIOD START and lives under
`dunestart_offset/`. A folder named `1978` or `2008` is a LINE VINTAGE and
lives under `raw_offset/` and `raster/`. Before 2026-09-15 the lines and masks
were filed under the start years they stand in for, so the same integer meant
two things; the rename is recorded in each `raster/<vintage>/RUN_MANIFEST.txt`
footer and `raw_offset/<vintage>/PROVENANCE.md`.

## Only one file here is a model input

```
dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv
```

`hatteras_site_config.py` resolves it as `PERIOD["road_setback_file"]` and
`HAT_hindcast_1984_2024.py:531` loads it by path. Two rows: GIS ids, then
setback in metres landward of interior row 0. Everything else in this tree is
provenance, diagnostics, or figures.

`raster/<vintage>/masks/` is an input too, but to the *pipeline* rather than
the model — both `HAT_road_offset_from_dune_start.py` and the dune-topo
extractor read it, through `hat_topo_version.road_mask_file()`. See the
warning under `raster/` below.

`dunestart_offset/derived/<year>/` holds the other two model inputs, the 1996
and 2010 setbacks. They are built from the measured files by
`1-produce/HAT_road_setback_derived_vintages.py` and never measured; the
folder name is what says so.

## Inventory

| path | written by | notes |
|---|---|---|
| `dunestart_offset/measured/<year>/RoadSetback_<year>_dunestart.csv` | `1-produce/HAT_road_offset_from_dune_start.py` | **the forcing** for 1984 and 2004. Floored at 0 |
| `dunestart_offset/derived/<year>/RoadSetback_<year>_dunestart.csv` | `1-produce/HAT_road_setback_derived_vintages.py` | **the forcing** for 1996 (= 1984 + the 1989 event) and 2010 (= 2004). Refuses to overwrite |
| `dunestart_offset/derived/<year>/PROVENANCE.md` | same | what each was derived from |
| `dunestart_offset/measured/<year>/RoadOffset_<year>_domains.csv` | same | per-domain detail; carries the **signed** `setback_dunestart_m` |
| `dunestart_offset/measured/<year>/RoadOffset_<year>_profiles.csv` | same | per-profile detail |
| `dunestart_offset/measured/<year>/RoadOffset_<year>_domains_rawframe.csv` | same | the **unstraightened control** pass, not the forcing |
| `dunestart_offset/measured/<year>/RoadElevation_<year>_dunestart.csv` | same | road elevation over the dune-start road profiles. Not the model's elevation file — that is `../road_elevation/RoadElevation.csv` |
| `dunestart_offset/RoadOffset_dunestart_audit.md` | same | **one document for both vintages.** A half run publishes a write-up missing a period |
| `dunestart_offset/RoadSetback_audit.{csv,md}` | `2-audit/HAT_road_setback_audit.py` | the only thing that tells you a domain would spend the hindcast as an unmanaged barrier wearing a road label |
| `dunestart_offset/HAT_dunestart_road_on_domains.png` | `1-produce/HAT_road_placement_on_domains.py` | |
| `dunestart_offset/HAT_road_island_planview_<year>.{png,pdf}` | `3-figures/HAT_road_island_planview.py` | the whole island in one view — the road as `roadway_manager.py` builds it (a flat 2-cell band per domain), coloured by `RoadElevation.csv`. Styled and framed as the extractor's plan view, which it overlays. 300 dpi raster + vector |
| `dunestart_offset/HAT_road_island_planview_<year>_caption.txt` | same | the figure caption. The figure carries **no title** — this is it, regenerated from the same constants the figure is drawn with so the two cannot drift |
| `dunestart_offset/modifications/*` | `3-figures/HAT_dunestart_modification_stages.py`, `3-figures/HAT_oceanfloor_offset_check.py` | what the flooring and the seaward relocation actually moved |
| `archive/superseded_20260911/<year>/RoadSetback_<year>.csv` (was `old_method_offset/`) | `1-produce/old_method/road_offset_pipeline.py` | the legacy 5-transect method. **Not stale, not a model input** — see below |
| `archive/superseded_20260911/RoadSetback_oldmethod_{audit.md,domains.csv}` | `1-produce/old_method/HAT_old_method_figures.py` | how the legacy number is built, and where it strains |
| `archive/superseded_20260911/HAT_old_method_*.png` | same, and `HAT_road_placement_on_domains.py` | |
| `raster/<vintage>/masks/domain_<N>_road_<vintage>.npy` | `1-produce/HAT_rasterize_road_to_domains.py` | 131 per LINE vintage (1978, 2008). **The only script that masks the road** |
| `raster/<vintage>/{HAT_road_mask_diagnostics,HAT_road_mask_summary}*` | same | |
| `raster/<vintage>/figures/` | same | 14 QC domains only, by design — `QC_DOMAINS` at that script's line 168 |
| `raster/<vintage>/RUN_MANIFEST.txt` | same | the run record as written; see the 2026-08-28 footer (two paths renamed by the 08-25 restructure) and the 2026-09-15 footer (the folder was `raster/<period>/`) |
| `raster/HAT_road_geojson_on_2009_dem.png` | `3-figures/HAT_road_geojson_map.py` | |
| `raw_offset/<vintage>/nc12_<vintage>.{csv,geojson}` | **ArcGIS, not this repo** | source data, filed by the imagery year (1978, 2008). See provenance below |
| `method_comparison/HAT_method_comparison_on_domains.png`, `.../HAT_method_vs_actual_road.png` | `4-compare/HAT_method_comparison_figures.py` | |
| `method_comparison/HAT_road_method_diagnostic.{csv,png}` | `4-compare/HAT_road_method_diagnostic.py` | |

**4-compare output lives in `method_comparison/`, never inside a method's
folder** — a legacy-vs-dune-start result belongs to neither method, and it is
not part of the forcing either. `HAT_road_method_diagnostic` wrote into
`dunestart_offset/` until 2026-08-28; that was a drift, now fixed. Both
4-compare scripts create `method_comparison/` on demand, so it can be deleted
whole and regenerated.

## The shape of this folder

```
road_offset/
    README.md              this file
    dunestart_offset/      the method that produces the forcing
        measured/1984/     measured: 1978 line vs 1984-start row 0
        measured/2004/     measured: 2008 line vs 2004-start row 0
        derived/1996/      1984 + the 1989 relocation (no line of its own)
        derived/2010/      2004, unchanged (no relocation 2004-2010)
        modifications/     what the flooring and the seaward move changed
        RoadSetback_audit.{csv,md}, RoadOffset_dunestart_audit.md, figures
    raster/1978/, 2008/    road masks per LINE vintage -- a pipeline input
    raw_offset/1978/, 2008/  the digitised NC-12 lines, from ArcGIS
    method_comparison/     old-vs-new. Read when choosing a method, not when
                           running the model
    archive/               (since 2026-09-18, as in 5-scr and 7-source-sink)
        superseded_20260907/   the v1-era 1984 measurement (see WHY.md there)
        superseded_20260911/   the legacy 5-transect method, kept for comparison
                               only (was old_method_offset/; see WHY.md there)
```

The top level holds the forcing, its source, and its inputs. Anything that
exists to *compare* the two methods sits in `method_comparison/` so it is not
mistaken for part of the product. The legacy method sat beside
`dunestart_offset/` as `old_method_offset/` until 2026-09-11, when it was
retired to a dated folder; it is under `archive/` now. Every script reaches it
through `hat_topo_version.LEGACY_SETBACK_ROOT` / `legacy_setback_file()`, and
every other path in this tree through the same module (2026-09-18).

## Source provenance for raw_offset/

The two NC-12 line files were digitised in ArcGIS and imported. Their ArcGIS
metadata sidecars (`nc12_<year>.csv.xml`, ~360 KB each) held the only record of
where they came from, and were removed on **2026-08-28** after the fact was
rescued here — the same treatment the `.tif.xml` sidecars got in the
2026-08-26 purge, and for the same reason: the sidecar is a display byproduct,
the origin is not.

```
raw_offset/1978/nc12_1978.csv   from \\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\nc12_1984.csv
raw_offset/2008/nc12_2008.csv   from \\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\nc12_2004.csv
created                         2024-10-10 15:56:55, both files
```

The origin files on the laptop are still named `nc12_1984` and `nc12_2004`;
they were imported under those names and renamed to their imagery vintage in
this repo on 2026-09-15 (`git log --follow` sees through the rename).

These two line vintages are **1978 and 2008 exports standing in for 1984 and
2004 on purpose** — not a bug. About 70% of their vertices are shared, so most
apparent "no movement" between them is an editing artefact rather than a
finding. Anything reading relocation out of these two lines needs to know that.

## What is deliberately NOT kept

`figures/` — 330 per-domain PNGs, 70 MB, deleted 2026-08-28. Nothing read them
and `.gitignore` (line 158, `*.png`) never tracked them. Regenerate what you
need:

```
python 3-figures/HAT_road_domain_views.py --domains 52 --year 1984
python 3-figures/HAT_road_domain_views.py --domains all --year 1984
```

`old_method_offset/RoadSetback_audit.{csv,md}` — deleted 2026-08-28. Orphans
from before `HAT_road_setback_audit.py` was repointed at `dunestart_offset/`
(it hardcodes that destination at line 152, so nothing had written these since
2026-08-17). They also collided by name with the live audit one folder over,
which is the actual reason they had to go. The legacy method's real audit is
`RoadSetback_oldmethod_audit.md`, beside them and current.

## Why the legacy setbacks are still here (archive/superseded_20260911/)

It is 1.9 MB and no model reads it, but three live scripts do:
`HAT_road_offset_from_dune_start.py:920` for `setback_legacy_m` /
`delta_vs_legacy_m`, `HAT_road_placement_on_domains.py` for its `old` entry,
and both `4-compare/` scripts.

Deleting it would **not raise**. `read_two_row_csv()` returns `{}` for a missing
path, so the producer would keep running and publish an audit whose migration
diagnostic is silently all-NaN, with the prose explaining that diagnostic still
in place. Removing the legacy method is a code change, not a folder delete.

## Last full rebuild

**2026-08-28.** 1984 re-measured on the live `1984-start/v1` after the 08-27
clear; 15 of 83 road domains moved, max 25 m; 2004 verified byte-identical;
all 83 road domains agree to 0.00 m with the extractor's independently computed
road column. The changelog in the scripts-side README has the detail.
