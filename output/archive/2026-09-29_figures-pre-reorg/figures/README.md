# figures — by subject

Every figure this project draws for a manuscript, a poster or a talk.
**Generated** by `scripts/figure_making/tools/figure_index.py`; re-run it
after adding a figure rather than editing this file.

One folder per subject, not per script (ORGANIZATION.md rule 1). A figure's
PNG sits at the top of its subject folder; its PDF and its caption go under
`supporting/` (scripts/figure_making/STYLE.md). `talk/` mirrors the subjects with
projector versions: water almost white, a point more type, heavier lines.
A retired figure goes to a dated `superseded_*` folder with a note, never
to the bin (rule 4).

## site

Where the reach is, how the 90 domains tile it, and what one domain is

| figure | shows | drawn by |
|---|---|---|
| `domain_framework.png` | The domain framework. | `scripts/figure_making/island/study_area_figures.py` |
| `domain_framework_vertical.png` | The domain framework, north up. | `scripts/figure_making/island/study_area_figures.py` |
| `domain_grid.png` | One domain as the model reads it (GIS 45, north of Avon, 2004-start product). | `scripts/figure_making/island/study_area_figures.py` |
| `domain_metrics.png` | What the 2004-start extraction gives each domain. | `scripts/figure_making/island/study_area_figures.py` |
| `domain_schematic.png` | How one domain is built (GIS 45, 2004-start dune-topo v1). | `scripts/figure_making/island/study_area_figures.py` |
| `reach_elevation.png` | The reach at the model's 10 m resolution, before the dune search: the per-domain elevation arrays of the 1984-start extraction (the 2009–2014 lidar with the 1996 ALACE survey grafted on the ocean side, upper strip of each pair) and the 2004-start extraction (the 2009–2014 lidar alone, lower strip),  | `scripts/figure_making/island/study_area_figures.py` |
| `site_overview.png` | The study reach. | `scripts/figure_making/island/study_area_figures.py` |
| `study_area.png` | Study area. | `scripts/figure_making/island/study_area_figures.py` |

## forcing

The record the model eats: storms, sea level, and the management timeline

| figure | shows | drawn by |
|---|---|---|
| `forcing_timeline_1984.png` | The forcing record over the 1984 chain, 1984–2024. | `scripts/figure_making/island/study_area_figures.py` |
| `forcing_timeline_1996.png` | The forcing record over the 1996 chain, 1996–2024. | `scripts/figure_making/island/study_area_figures.py` |

## management

NC-12 and beach nourishment: where, when, and under what rules

| figure | shows | drawn by |
|---|---|---|
| `gis11_relocation_drown.png` | GIS 11: a standard relocation setback moves where a prescribed relocation lands. | `scripts/figure_making/model_output/gis11_relocation_drown_figure.py` |
| `management_footprint.png` | Where the reach has been managed, as the model prescribes it, on one reach in two panels. | `scripts/figure_making/island/study_area_figures.py` |
| `rules_table.png` | The management rules the hindcast applies. | `—` |
| `rules_table_slide.png` | The presentation rendering of rules_table.png: the same rows, set larger, carrying the title and the note on the canvas because a slide has no caption beside it. | `—` |
| `timeline_1984_2024.png` | The management record, 1984–2024, by domain and year: the beach-nourishment projects, the NC-12 relocations and the Rodanthe bridge, against the community zones and the run periods. | `—` |
| `timeline_1996_2024.png` | The management the model applies, 1996–2024, by domain and year: the beach-nourishment projects, the NC-12 relocations and the Rodanthe bridge, against the community zones and the run periods. | `—` |

## shoreline

Shoreline change, observed and modelled

| figure | shows | drawn by |
|---|---|---|
| `coastsat_calibration_periods.png` | Observed shoreline change rate per model domain, CoastSat, one curve per run period: 1996–2010 in red and 2010–2024 in blue, positive seaward. | `—` |
| `coastsat_endpoint_vs_duneline_1996_2010_2024_stacked.png` | Net change in shoreline and dune-line position by GIS domain (1 at Cape Point, 90 at Pea Island), seaward positive, in metres: (a) 1997–2023, (b) 1997–2009, (c) 2009–2023, the digitized dune lines standing in for the model years 1996, 2010 and 2024. | `—` |
| `dune_lines.png` | The digitised dune line through the record. | `scripts/figure_making/island/study_area_figures.py` |
| `hindcast_edgeBE.png` | CASCADE hindcast with the source/sink calibration removed, against the CoastSat LOESS reference, Cape Hatteras. | `—` |
| `hindcast_zeroBE.png` | CASCADE hindcast with no source/sink field, against the CoastSat LOESS reference, Cape Hatteras. | `—` |
| `observed_rates.png` | Observed shoreline change, the calibration target. | `scripts/figure_making/island/study_area_figures.py` |
| `scenario_grid.png` | Modelled shoreline change rate by domain for every management scenario, 1996–2010 and 2010–2024, under each source/sink preset. | `scripts/figure_making/model_output/scenario_grid.py` |
| `two_periods_7_domains.png` | The LOESS-smoothed companion to coastsat_calibration_periods.png: the same CoastSat domain-mean LRR rates over the same two run periods, 1996–2010 in red and 2010–2024 in blue, positive seaward, smoothed with a LOESS window of 7 domains (3.5 km, frac=0.078). | `—` |
| `duneline_positions/beach_width.png` | Beach width, the distance from the dune line to the CoastSat shoreline, by GIS domain (1 at Cape Point, 90 at Pea Island) in 1997, 2009 and 2023, the dune lines that stand for the model years 1996, 2010 and 2024: along each CoastSat transect, the shoreline position (the mean satellite position withi | `—` |
| `duneline_positions/dune_to_nc12.png` | Distance from the dune line to the NC-12 centreline by GIS domain (1 at Cape Point, 90 at Pea Island) in 1997, 2009 and 2023, the dune lines that stand for the model years 1996, 2010 and 2024: along each 100 m transect, the road's distance from the fixed offshore datum minus the dune line's (both fo | `—` |
| `duneline_positions/duneline_positions_overview.png` | Where the digitized dune line sat along Hatteras Island in 1997, 2009 and 2023, the lines that stand for the model years 1996, 2010 and 2024, north-up in three segments at one common scale: (a) Cape Point to Avon, GIS 1–30; (b) Avon to the Tri-Village, GIS 31–60; (c) the Tri-Village to GIS 90. | `—` |
| `duneline_positions/duneline_positions_zooms.png` | The dune line in 1997, 2009 and 2023 at five sites, each a window three model domains (1.5 km) alongshore by 950 m across, 650 m landward and 300 m seaward of the 2023 line, at one scale: (a) Buxton, centred on GIS 7; (b) Avon, centred on GIS 24; (c) Tri-Village, centred on GIS 78; (d) Mirlo Beach S | `—` |
| `duneline_positions/zoom_avon.png` | The dune line at Avon, GIS 23–25 (named site GIS 21-31; centred on its largest |net dune change| 1997-2023, GIS 24 (-73.3 m)), in 1997, 2009 and 2023. | `—` |
| `duneline_positions/zoom_buxton.png` | The dune line at Buxton, GIS 6–8 (named site GIS 1-15; centred on its largest |net dune change| 1997-2023, GIS 7 (-65.9 m)), in 1997, 2009 and 2023. | `—` |
| `duneline_positions/zoom_mirlo.png` | The dune line at Mirlo Beach S-curves, GIS 83–85 (named site GIS 84-90; centred on its largest |net dune change| 1997-2023, GIS 84 (-79.3 m)), in 1997, 2009 and 2023. | `—` |
| `duneline_positions/zoom_picked.png` | The dune line at Largest change elsewhere, GIS 32–34 (the domain outside every named site with the largest |net dune change| 1997-2023, GIS 33 (+43.0 m)), in 1997, 2009 and 2023. | `—` |
| `duneline_positions/zoom_trivillage.png` | The dune line at Tri-Village, GIS 77–79 (named site GIS 68-83; centred on its largest |net dune change| 1997-2023, GIS 78 (-90.4 m)), in 1997, 2009 and 2023. | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_GIS04_buxton.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_GIS26_avon_pier.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_GIS40_central.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_GIS55_central_north.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_GIS79_rodanthe_pier.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_GIS88_mirlo_beach.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_island_1996.png` | — | `—` |
| `mean_shoreline/line_and_band/mean_shoreline_1995_1997_on_imagery_ribbon_1996.png` | — | `—` |
| `mean_shoreline/with_positions/mean_shoreline_1995_1997_on_imagery_with_positions_GIS04_buxton.png` | — | `—` |
| `mean_shoreline/with_positions/mean_shoreline_1995_1997_on_imagery_with_positions_GIS26_avon_pier.png` | — | `—` |
| `mean_shoreline/with_positions/mean_shoreline_1995_1997_on_imagery_with_positions_GIS40_central.png` | — | `—` |
| `mean_shoreline/with_positions/mean_shoreline_1995_1997_on_imagery_with_positions_GIS55_central_north.png` | — | `—` |
| `mean_shoreline/with_positions/mean_shoreline_1995_1997_on_imagery_with_positions_GIS79_rodanthe_pier.png` | — | `—` |
| `mean_shoreline/with_positions/mean_shoreline_1995_1997_on_imagery_with_positions_GIS88_mirlo_beach.png` | — | `—` |

## initialization

The island as the model starts it

| figure | shows | drawn by |
|---|---|---|
| `1984/classes/island_1984.png` | The island as the model starts it in 1984, from the 1984-start v2 extraction: every domain's elevation array at the model's 10 m resolution, in the elevation classes of the house style, m MHW, water (below 0 m) one colour. | `scripts/figure_making/island/initialization_figures.py` |
| `1984/classes/island_1984_absolute.png` | The island as the model starts it in 1984, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `1984/classes/island_1984_absolute_buffers.png` | The island as the model starts it in 1984, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `1984/terrain/island_1984.png` | The island as the model starts it in 1984, from the 1984-start v2 extraction: every domain's elevation array at the model's 10 m resolution, on a continuous terrain ramp over the land range 0-4 m MHW, with water (below 0 m) a single deeper blue under the same rule the class scheme uses. | `scripts/figure_making/island/initialization_figures.py` |
| `1984/terrain/island_1984_absolute.png` | The island as the model starts it in 1984, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `1984/terrain/island_1984_absolute_buffers.png` | The island as the model starts it in 1984, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `1996/classes/island_1996.png` | The island as the model starts it in 1996, from the 1984-start v2 extraction: every domain's elevation array at the model's 10 m resolution, in the elevation classes of the house style, m MHW, water (below 0 m) one colour. | `scripts/figure_making/island/initialization_figures.py` |
| `1996/classes/island_1996_absolute.png` | The island as the model starts it in 1996, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `1996/classes/island_1996_absolute_buffers.png` | The island as the model starts it in 1996, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `1996/terrain/island_1996.png` | The island as the model starts it in 1996, from the 1984-start v2 extraction: every domain's elevation array at the model's 10 m resolution, on a continuous terrain ramp over the land range 0-4 m MHW, with water (below 0 m) a single deeper blue under the same rule the class scheme uses. | `scripts/figure_making/island/initialization_figures.py` |
| `1996/terrain/island_1996_absolute.png` | The island as the model starts it in 1996, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `1996/terrain/island_1996_absolute_buffers.png` | The island as the model starts it in 1996, from the 1984-start v2 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `2004/classes/island_2004.png` | The island as the model starts it in 2004, from the 2004-start v1 extraction: every domain's elevation array at the model's 10 m resolution, in the elevation classes of the house style, m MHW, water (below 0 m) one colour. | `scripts/figure_making/island/initialization_figures.py` |
| `2004/classes/island_2004_absolute.png` | The island as the model starts it in 2004, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `2004/classes/island_2004_absolute_buffers.png` | The island as the model starts it in 2004, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `2004/terrain/island_2004.png` | The island as the model starts it in 2004, from the 2004-start v1 extraction: every domain's elevation array at the model's 10 m resolution, on a continuous terrain ramp over the land range 0-4 m MHW, with water (below 0 m) a single deeper blue under the same rule the class scheme uses. | `scripts/figure_making/island/initialization_figures.py` |
| `2004/terrain/island_2004_absolute.png` | The island as the model starts it in 2004, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `2004/terrain/island_2004_absolute_buffers.png` | The island as the model starts it in 2004, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `2010/classes/island_2010.png` | The island as the model starts it in 2010, from the 2004-start v1 extraction: every domain's elevation array at the model's 10 m resolution, in the elevation classes of the house style, m MHW, water (below 0 m) one colour. | `scripts/figure_making/island/initialization_figures.py` |
| `2010/classes/island_2010_absolute.png` | The island as the model starts it in 2010, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `2010/classes/island_2010_absolute_buffers.png` | The island as the model starts it in 2010, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |
| `2010/terrain/island_2010.png` | The island as the model starts it in 2010, from the 2004-start v1 extraction: every domain's elevation array at the model's 10 m resolution, on a continuous terrain ramp over the land range 0-4 m MHW, with water (below 0 m) a single deeper blue under the same rule the class scheme uses. | `scripts/figure_making/island/initialization_figures.py` |
| `2010/terrain/island_2010_absolute.png` | The island as the model starts it in 2010, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains only. | `scripts/figure_making/island/initialization_figures.py` |
| `2010/terrain/island_2010_absolute_buffers.png` | The island as the model starts it in 2010, from the 2004-start v1 extraction, with every domain at its TRUE cross-shore position: each domain's elevation array placed at its own BRIE shoreline offset, so this canvas is the initial condition the run begins from, the 90 real domains and the 15 buffer  | `scripts/figure_making/island/initialization_figures.py` |

## model

How Barrier3D, BRIE and CASCADE work, drawn from a finished run: the grid through a storm, alongshore diffusion, the coupling loop, overwash routing

| figure | shows | drawn by |
|---|---|---|
| `barrier3d_domain_budget.png` | One Barrier3D domain over a hindcast window: GIS 45, 1996-2010, the natural run (edgeBE, no management). | `scripts/figure_making/model/model_mechanics_figures.py` |
| `barrier3d_storm_year.png` | Barrier3D through one storm year: GIS 6, model year 2006 (the natural 1996-2010 run, edgeBE, no road, no beach/dune management). | `scripts/figure_making/model/model_mechanics_figures.py` |
| `brie_diffusion.png` | How BRIE moves sand alongshore in CASCADE. | `scripts/figure_making/model/model_mechanics_figures.py` |
| `cascade_coupling_loop.png` | CASCADE's annual loop (cascade/cascade_groin.py, Cascade.update) and the unit each exchange is made in. | `scripts/figure_making/model/model_mechanics_figures.py` |
| `cascade_island_grids.png` | The island as CASCADE holds it at the start of the 1996 window. | `scripts/figure_making/model/model_mechanics_figures.py` |
| `cascade_shoreline_split.png` | How CASCADE divides a shoreline's change between its two models: the natural run, 1996-2010 (edgeBE, no management), positive seaward. | `scripts/figure_making/model/model_mechanics_figures.py` |
| `management_modules.png` | What the two CASCADE management modules do to the grid. | `scripts/figure_making/model/model_mechanics_figures.py` |
| `storm_routing_check.png` | The storm replay is the model: GIS 6, the 2006 storms of the natural 1996-2010 run. | `scripts/figure_making/model/overwash_routing_figures.py` |
| `storm_routing_hours.png` | One storm hour by hour: GIS 6 entering 2006, a 24 h storm with Rhigh 3.4 m MHW (over the dune, run-up routing), replayed through Barrier3D. | `scripts/figure_making/model/overwash_routing_figures.py` |
| `storm_routing_ladder.png` | The same Barrier3D grid (GIS 6, entering 2006, natural 1996-2010 run) hit by four 24 h storms of rising strength, each replayed through the model's own update from that grid. | `scripts/figure_making/model/overwash_routing_figures.py` |
| `storm_routing_response.png` | Does overwash grow with the storm? GIS 6 entering 2006 (natural 1996-2010 run), one 24 h storm at a time replayed through Barrier3D, Rhigh stepped from below the lowest dune crest to above the highest (crest after the year's dune growth: lowest 2.17 m dotted, mean 2.95 m solid; berm 1.34 m). | `scripts/figure_making/model/overwash_routing_figures.py` |

## pipeline

How each model input is built from its source data, one folder per input_prep step

| figure | shows | drawn by |
|---|---|---|
| `0-elevation/dem_resample_one_domain.png` | From surveys to the 10 m grid, one domain (GIS 45, the 2009-2014-1996 product; the 2009-2014 product is the same without the 1996 layer). | `scripts/figure_making/pipeline/0-elevation/dem_composition_figures.py` |
| `0-elevation/dem_sources_alongshore.png` | Which survey supplies the model's topography, domain by domain: the share of each domain's 1 m land cells (above MHW) taken from each survey. | `scripts/figure_making/pipeline/0-elevation/dem_composition_figures.py` |
| `1-barrier3d-domains/domain_extraction_gis45.png` | How one domain's Barrier3D arrays are made: GIS 45, the 2004-start product, dune-topo v1, drawn by the extractor's own functions (HAT_dune_topo_extractor.py) on its saved window pick, and checked to reproduce the saved arrays exactly. | `—` |
| `2-brie-offset/offset_build_1996.png` | How the BRIE shoreline offset for the 1996 start is built from the 1997 dune line (build data/hatteras_init/2-brie-offset/1996/duneline/v1). | `—` |
| `2-brie-offset/offset_build_2010.png` | How the BRIE shoreline offset for the 2010 start is built from the 2009 dune line (build data/hatteras_init/2-brie-offset/2010/duneline/v1). | `—` |
| `3-storms/storm_construction_steps.png` | How the model's storm series is built, on September 2003. | `scripts/figure_making/pipeline/3-storms/storm_construction_figures.py` |
| `3-storms/storm_events_by_duration.png` | Every storm event the generator finds at Duck, 1996-2024, and which ones reach the model at a 72 h maximum duration. | `scripts/figure_making/pipeline/3-storms/storm_construction_figures.py` |
| `4-mgmt-forcings/road_setback_inputs.png` | The NC-12 setbacks the two runs read, metres landward of interior row 0, and where they come from. | `scripts/figure_making/pipeline/4-mgmt-forcings/road_setback_figures.py` |
| `4-mgmt-forcings/road_setback_measurement.png` | How the NC-12 setback is measured on one domain: GIS 31, the 2004 measurement (the 2008-digitised centreline on the 2004-start extraction), which the 2010 run reads unchanged. | `scripts/figure_making/pipeline/4-mgmt-forcings/road_setback_figures.py` |
| `5-scr/observed_target_1996_2010.png` | How the observed shoreline-change target for 1996-2010 is built from CoastSat. | `—` |
| `5-scr/observed_target_2010_2024.png` | How the observed shoreline-change target for 2010-2024 is built from CoastSat. | `—` |
| `7-source-sink/be_end_solve.png` | How the two end-domain source/sink values that the current (edgeBE) runs carry were solved, 2026-09-27, at wave option A (Hs 2.0 m, Tp 7.5 s, asymmetry 0.6, high-angle fraction 0.5) on the metres island offset, full-management run, against each window's CoastSat LRR (GIS 1 the raw domain mean, GIS 9 | `scripts/figure_making/pipeline/7-source-sink/be_method_figures.py` |
| `7-source-sink/be_zone_field.png` | The zone-by-zone source/sink calibration (calibBE) for the 1996-2010 / 2010-2024 pair, as made on 2026-09-18 (7-source-sink/2-calibrate/1996_2010__2010_2024/): a single pass, BEFORE the island offset went to metres and before wave option A. | `scripts/figure_making/pipeline/7-source-sink/be_method_figures.py` |

## talk

The same 13 figures for a projector, under `talk/<subject>/`. Drawn by the
same scripts with `--talk`.
