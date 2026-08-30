# road_offset — what each file is, and what wrote it

Where NC-12 sits per Barrier3D domain. The **method** documentation lives with
the code, at `scripts/input_prep/4-mgmt-forcings/road_offset/README.md` — that
is where the two methods are explained and compared. This file is the inventory:
for every file here, what produced it, whether the model reads it, and whether
losing it costs anything.

Two vintages throughout, `1984` and `2004`. Each is measured against **its own
period's extraction** — `1984-start` and `2004-start` — because the two are
different islands. That pairing is defined once, in
`scripts/hat_topo_version.py:YEAR_PRODUCT`.

## Only one file here is a model input

```
dunestart_offset/<year>/RoadSetback_<year>_dunestart.csv
```

`hatteras_site_config.py` resolves it as `PERIOD["road_setback_file"]` and
`HAT_hindcast_1984_2024.py:531` loads it by path. Two rows: GIS ids, then
setback in metres landward of interior row 0. Everything else in this tree is
provenance, diagnostics, or figures.

`raster/<year>/masks/` is an input too, but to the *pipeline* rather than the
model — both `HAT_road_offset_from_dune_start.py` and the dune-topo extractor
read it. See the warning under `raster/` below.

## Inventory

| path | written by | notes |
|---|---|---|
| `dunestart_offset/<year>/RoadSetback_<year>_dunestart.csv` | `1-produce/HAT_road_offset_from_dune_start.py` | **the forcing.** Floored at 0 |
| `dunestart_offset/<year>/RoadOffset_<year>_domains.csv` | same | per-domain detail; carries the **signed** `setback_dunestart_m` |
| `dunestart_offset/<year>/RoadOffset_<year>_profiles.csv` | same | per-profile detail |
| `dunestart_offset/<year>/RoadOffset_<year>_domains_rawframe.csv` | same | the **unstraightened control** pass, not the forcing |
| `dunestart_offset/<year>/RoadElevation_<year>_dunestart.csv` | same | road elevation over the dune-start road profiles. Not the model's elevation file — that is `../road_elevation/RoadElevation.csv` |
| `dunestart_offset/RoadOffset_dunestart_audit.md` | same | **one document for both vintages.** A half run publishes a write-up missing a period |
| `dunestart_offset/RoadSetback_audit.{csv,md}` | `2-audit/HAT_road_setback_audit.py` | the only thing that tells you a domain would spend the hindcast as an unmanaged barrier wearing a road label |
| `dunestart_offset/HAT_dunestart_road_on_domains.png` | `1-produce/HAT_road_placement_on_domains.py` | |
| `dunestart_offset/HAT_road_island_planview_<year>.{png,pdf}` | `3-figures/HAT_road_island_planview.py` | the whole island in one view — the road as `roadway_manager.py` builds it (a flat 2-cell band per domain), coloured by `RoadElevation.csv`. Styled and framed as the extractor's plan view, which it overlays. 300 dpi raster + vector |
| `dunestart_offset/HAT_road_island_planview_<year>_caption.txt` | same | the figure caption. The figure carries **no title** — this is it, regenerated from the same constants the figure is drawn with so the two cannot drift |
| `dunestart_offset/modifications/*` | `3-figures/HAT_dunestart_modification_stages.py`, `3-figures/HAT_oceanfloor_offset_check.py` | what the flooring and the seaward relocation actually moved |
| `old_method_offset/<year>/RoadSetback_<year>.csv` | `1-produce/old_method/road_offset_pipeline.py` | the legacy 5-transect method. **Not stale, not a model input** — see below |
| `old_method_offset/RoadSetback_oldmethod_{audit.md,domains.csv}` | `1-produce/old_method/HAT_old_method_figures.py` | how the legacy number is built, and where it strains |
| `old_method_offset/HAT_old_method_*.png` | same, and `HAT_road_placement_on_domains.py` | |
| `raster/<year>/masks/domain_<N>_road_<year>.npy` | `1-produce/HAT_rasterize_road_to_domains.py` | 131 per year. **The only script that masks the road** |
| `raster/<year>/{HAT_road_mask_diagnostics,HAT_road_mask_summary}*` | same | |
| `raster/<year>/figures/` | same | 14 QC domains only, by design — `QC_DOMAINS` at that script's line 168 |
| `raster/<year>/RUN_MANIFEST.txt` | same | see the 2026-08-28 footer: two paths in it were renamed by the 08-25 restructure |
| `raster/HAT_road_geojson_on_2009_dem.png` | `3-figures/HAT_road_geojson_map.py` | |
| `raw_offset/<year>/nc12_<year>.{csv,geojson}` | **ArcGIS, not this repo** | source data. See provenance below |
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
    old_method_offset/     the legacy method, kept for comparison only
    raster/                road masks -- an input to the pipeline
    raw_offset/            the digitised NC-12 lines, from ArcGIS
    method_comparison/     old-vs-new. Read when choosing a method, not when
                           running the model
```

The top level holds the forcing, its source, and its inputs. Anything that
exists to *compare* the two methods sits in `method_comparison/` so it is not
mistaken for part of the product. `old_method_offset/` stays a sibling of
`dunestart_offset/` rather than moving under it: they are two methods, and
that parallel is a fact about the tree. Which one the model reads is stated
above, in words, rather than encoded in the nesting.

## Source provenance for raw_offset/

The two NC-12 line files were digitised in ArcGIS and imported. Their ArcGIS
metadata sidecars (`nc12_<year>.csv.xml`, ~360 KB each) held the only record of
where they came from, and were removed on **2026-08-28** after the fact was
rescued here — the same treatment the `.tif.xml` sidecars got in the
2026-08-26 purge, and for the same reason: the sidecar is a display byproduct,
the origin is not.

```
nc12_1984.csv   \\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\nc12_1984.csv
nc12_2004.csv   \\HANNAHS-LAPTOP\D$\Hatteras_GIS\Roads\nc12_2004.csv
created         2024-10-10 15:56:55, both files
```

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

## Why old_method_offset/ is still here

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
