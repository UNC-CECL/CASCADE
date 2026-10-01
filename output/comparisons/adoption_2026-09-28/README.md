# The matrix before and after the 2026-09-28 adoption

- **Before:** `raw_runs/archive/2026-09-28-pre-ceiling/matrix`. Barrier3D 49fd069, Dmaxel default, `v3_72` storms, LOWESS-7 option A ends.
- **After:** `raw_runs/matrix`. Barrier3D `hatteras/adopted` (overwash fixes and per-cell dune ceilings), `v3_trim24` storms, the beach/dune manager's 4 m cap limited to the sand it adds (see below), re-solved ends (1996 +4.3509/+19.0935, 2010 +8.0/+21.2582).

Produced by `scripts/analyze_output/compare_runs/adoption_2026-09-28/adoption_before_after.py`. Each side was scored under the Barrier3D it ran on.

| File | Contents |
|---|---|
| `scores_<side>.csv` | one row per run |
| `cells_<side>.csv` | every image × domain: observed and model overwash |
| `crest_<side>.csv` | end-of-run dune crest per domain |
| `crest_lidar_2009.csv` | the 2010-start crest |
| `before_after.csv` | the two sides paired |
| `supporting_precap/` | the after side as scored before the dune-cap fix |

The shoreline score is the interior (GIS 2-89) RMSE and bias against CoastSat LOWESS-7. The overwash score compares each run with the imagery (the date rule in `8-overwash-analysis`).

## Result (edgeBE; zeroBE is the same within 0.1)

| | 1996-2010 managed | 1996-2010 natural | 2010-2024 managed | 2010-2024 natural |
|---|---|---|---|---|
| RMSE, m/yr | 1.19 → 1.17 | 1.14 → 1.14 | 2.31 → 2.07 | 4.18 → 2.72 |
| bias, m/yr | +0.10 → +0.06 | −0.18 → −0.09 | −1.66 → −1.36 | −3.87 → −2.19 |
| overwash PSS | −0.07 → 0.60 | −0.14 → 0.56 | 0.36 → 0.17 | 0.03 → 0.12 |
| POD / POFD | 0.17/0.24 → 0.78/0.19 | 0.28/0.41 → 0.81/0.25 | 0.81/0.44 → 0.57/0.40 | 0.95/0.92 → 0.61/0.49 |
| end crest, median m MHW | 3.04 → 5.48 | 1.42 → 5.51 | 2.41 → 4.95 | 1.55 → 4.88 |

The 2010 crest from the 1996 runs minus the 2009 lidar went from −2.7 m (managed) and −3.4 m (natural) to −0.17 m and −0.27 m.

The 2010-2024 managed PSS fell because the old model overwashed about half the island in every image window. Its hits came from blanket coverage: POFD was 0.44. It now misses cells the imagery shows. Irene (2011) is sound-side, which Barrier3D cannot produce, so part of that loss is outside the model.

## The dune-cap fix (2026-09-28, late)

Scored first with the adopted Barrier3D, the managed runs showed village crests flat at 5.34 m MHW, 4.0 m above the berm. The cause was CASCADE's `beach_dune_manager`. `filter_overwash` clipped every dune cell to `_artificial_maximum_dune_height = 4` m above the berm ("parameterized for Nags Head, NC"), every year, in every domain the manager covers. The cap exists to stop the bulldozed overwash sand building 10 m dunes. Written that way, it also cut natural dunes: 24 domains in 1996-2010 and 17 in 2010-2024, about 150,000 m³ of starting dune each time. It had never applied before, because the Dmaxel default kept dunes near 3 m.

Hannah chose "option 1, clip only the bulldozed sand". The cap now limits only what the manager adds: a cell already above 4 m keeps its height and takes no sand (`DUNE_CAP_APPLIES_TO = "added sand only"`, recorded in run metadata). The 2010 GIS 90 end moved +22.4937 → +21.2582. The superseded runs are in `raw_runs/archive/2026-09-28-pre-dunecap/`.

Effect on the managed runs, after the adoption, cap on whole cells → cap on added sand:

| | 1996-2010 full mgmt | 1996-2010 beach/dune only | 2010-2024 full mgmt | 2010-2024 beach/dune only |
|---|---|---|---|---|
| RMSE, m/yr | 1.17 → 1.17 | 1.16 → 1.16 | 2.12 → 2.07 | 2.26 → 2.20 |
| bias, m/yr | −0.03 → +0.06 | −0.11 → −0.02 | −1.45 → −1.36 | −1.74 → −1.65 |
| 2010 crest − lidar, m | −0.31 → −0.17 | −0.31 → −0.15 | | |

Overwash scores did not change at two decimals. The natural and roadway-only runs do not use the manager.

## Figures (`figures/`, made by `scripts/analyze_output/compare_runs/adoption_2026-09-28/adoption_before_after_figures.py`)

The captions are in `figures/supporting/CAPTIONS.md`.

| Figure | What it shows |
|---|---|
| `adoption_scorecard.png` | every score, before → after, one row per scenario |
| `adoption_shoreline_alongshore.png` | model LRR vs CoastSat, managed and natural, both windows |
| `adoption_overwash_map_full_management.png`, `_natural.png` | image × domain: hit, miss, false alarm, before and after |
| `adoption_overwash_by_image.png` | when: domains overwashed per image, observed vs before vs after |
| `adoption_dune_crest_2010.png` | the 1996-2010 runs' 2010 crest vs the 2009 lidar |
