# 1984-start `dune-topo` — version index

Two versions, both **extractions** (all 90 domains cut from the DEM by
`HAT_dune_topo_extractor.py`), no rows added, no cell modified:

```
v1   extraction, the ORIGINAL pick set (2026-08-27)        67 MB   the 2026-09-01 calibration tree and the pea1989base arm ran on it
v2   extraction, re-pick with NC-12 visible (2026-09-02)   67 MB   THE BASE: what CURRENT says, what the road tree measures against
v3   v2 + the 1984 footprint BEHIND THE ROAD, copy fill (2026-09-08)   a layer on v2, built by HAT_build_footprint_version.py; test arm `behindroad-copy`
```

`v3` is the first version built on the new approach: the symmetric 1984
footprint (rows added where the 1984 dune line lay seaward of the 1997 line,
removed where it lay landward; 52 domains, +73/−47 rows) placed directly behind
the model's two roadway rows and filled by copying the N rows that follow the
insert point. Row 0 and the dune array do not move; **its setback CSV carries
the 1984 setbacks** (`setback_new_m`: the road against the 1984 dune line,
row-0 convention, no floor - GIS 85/86 go from 0 to 44/12 m); the block sits
directly behind the road as placed, so the road is on measured cells and the
added width is behind it; removed rows come out of the interior directly in
front of today's road (2026-09-08), so the road's cells and everything behind
them are kept. Placement checked against the 1984 lines on the lidar:
`2-domain-reconstruction-1984/3-placement/HAT_road_placement_check_1984.txt`. `CURRENT` still says
`v2`. To run on it, set `HAT_TOPO_VERSION_1984_START=v3` AND copy its CSV over
the forcing-tree one for the run. Its README and
`HAT_footprint_audit.csv` say what was done to every domain; the argument is in
`../2-domain-reconstruction-1984/`.

| | full extraction? | pick set | rows added | own setback CSV |
|---|---|---|---|---|
| `v1` | yes, 90 domains | original (08-27), `../picks/HAT_dune_search_windows_v1.json` | no | yes — the v1-era measurement (GIS 85 −10 m floored to 0), saved 09-04 from the `dunestart_offset_ARCHIVE_1984start_v1` archive |
| `v2` | yes, 90 domains | **re-pick (09-02)**, `../picks/HAT_dune_search_windows_v2.json` | no | yes — the road tree's measurement on v2, saved 09-04; byte-identical to the forcing-tree file |
| `v3` | no — a layer on `v2` | v2's | **yes, both signs**: +73 rows at 29 domains behind NC-12 as placed, −47 rows at 23 from directly in front of it | yes — the 1984 setbacks from the footprint table, no floor |

Each version holds `topography/` (90 `domain_<N>_topography.npy` + nodata
masks), `dunes/`, its settings CSV/XLSX, plan-view and offset figures, a
`RUN_MANIFEST.txt`, and `RoadSetback_1984_dunestart.csv` matched to its own
arrays. **A setback CSV read against the other version's arrays is off by the
pick difference** — keep them paired.

## Which one loads — `CURRENT` decides (since 2026-09-04), and it says `v2`

`scripts/site_layer/hat_topo_version.py` resolves in this order:

1. an explicit `override=` argument
2. the `HAT_TOPO_VERSION_1984_START` environment variable
3. **the `CURRENT` file in this folder**
4. the extractor's own `VERSION` literal (`v2`, what it *writes*)
5. the only version present, if there is exactly one

Until 2026-09-04 rules 3 and 4 were the other way round, so `CURRENT` recorded
intent and did nothing. They were swapped so a version the extractor did not
write could be made the default. What survives of that: a fresh extraction is
**not** adopted until `CURRENT` says so, and the road tree
(`HAT_road_offset_from_dune_start.py`) measures against `CURRENT` too.

### The forcing-tree setback CSV is the other half

`hatteras_site_config.py:142` hardcodes
`4-mgmt-forcing/road_offset/dunestart_offset/measured/1984/RoadSetback_1984_dunestart.csv`.
It is the **v2-measured** file, a copy of which is saved in `v2/`. Changing the
default version means both steps:

```
echo v1 > dune-topo\CURRENT
copy dune-topo\v1\RoadSetback_1984_dunestart.csv ^
     ..\..\4-mgmt-forcing\road_offset\dunestart_offset\1984\
```

For one run without touching either, set `HAT_TOPO_VERSION_1984_START=<version>`
and copy that version's CSV over the live one inside a `finally`, as
`HAT_run_row_insert_set.py` does.

## Which calibration-tree runs sit on which version

The 1984–2004 calibration tree (`output/raw_runs/1984_2004/<preset>/`) was run
2026-09-01 on **v1**. One run, `calibBE/HAT_1984_2004_calibBE_road_bdm_groin`,
was re-run 2026-09-04 on the median-fill layer and again on **2026-09-07 on
v2**, so it is the only tree run on v2 and the tree is not on one version.
Read a run's `*_run_metadata.txt` (`topo_dune_version`), not this file, for
what it ran on; `run_index.csv` carries the same column.

## What was removed, and when

**The layers `v3`–`v8` — deleted 2026-09-07.** Six versions built on `v2` by
`HAT_insert_seaward_rows.py` for the 1984 seaward-row insert: rows prepended at
the seaward edge where the 1984 dune line stood seaward of the 1996 one (`v3`
at the 10 relocation-block domains; `v4`–`v8` at all 38 measured domains, one
per fill rule — measured + floor, measured + median, flat platform, matched
backdune with and without the crest). Hannah's decision: keep only unmodified
topography. Every run made on modified topography went too — the six-fill set
under `output/raw_runs/row-insert/` (controls included), `blocksv4*`,
`islandv5`, and the 2026-09-03 arms whose topography was already gone — with
their `run_index.csv` rows, `output/experiments/row_insert_set/` and
`output/comparisons/relocation_1984_2004/row-insert/`. The 09-04
calibration-tree run on the median fill was replaced by a v2 run. Sizes and
reasons: `../../archive_purge_20260907.csv`; lineage in `../../LINEAGE.md`.

What the layers established is recorded, not lost: the fill rule does not
matter island-wide (interior RMSE 0.540–0.547 across all seven arms), the
insert itself makes every emergent relocation 3–10 years late, and the
no-insert interiors reproduce the published 0.30 / 0.40 recall. The method,
the measurement of N and the figures stay in `../2-domain-reconstruction-1984/`, whose
`DUNE_TOPO_VERSION_GUIDE.md` describes each deleted version and the recipe
that built it. The scripts that ran or drew the layers are kept with those
arms retired; `hat_topo_version.require_version()` fails loudly if one is
pointed at a version that is not here.

**`v2` (as built: `v1` + block-scope rows on the old picks) — deleted 2026-09-03.**
Superseded by the re-pick layers before any of the above; its run output
(`blocksdate*`) survived until 2026-09-07. `../../archive_purge_20260903.csv`.

**The 2026-08-27 clear.** The earlier `v1`/`v2` of 08-26 were deleted and
numbering restarted; `../../archive_purge_20260826.csv`.

Retired experimental variants lived in `../dune-topo-experiments/`, outside
this folder because `hat_topo_version.versions()` lists the directories here.
It has held no topography since 2026-09-03.

## No figures in this folder

Insert figures live in `../2-domain-reconstruction-1984/figures/`; each version keeps its
own `figures/` subfolder. Plotters default there via
`hat_topo_version.insert_figures_dir`.

Written 2026-09-03; layers added 2026-09-04; rewritten 2026-09-07 when the
layers were removed.
