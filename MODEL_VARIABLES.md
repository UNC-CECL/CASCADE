# CASCADE model variables

Every input CASCADE takes, its default, and the value the Hatteras hindcast uses.

- **Default**: the upstream value (constructor default, BRIE coupler constant, or
  the CASCADE template `data/default_cascade_variables/barrier3d-default-parameters.yaml`).
- **Hatteras**: what `scripts/hatteras_ms/HAT_hindcast_1984_2024.py` passes through
  `build_cascade()` (`scripts/cascade_pipeline/hindcast.py`), with wave and run
  settings from `hat_run.yaml` / `HAT_*` variables, and Barrier3D settings from
  `data/hatteras_init/Hatteras-CASCADE-parameters.yaml`.
- **Units**: as you enter them. Barrier3D works internally in decametres relative
  to MHW and converts the scalar parameters itself (`/10`, and `/10 - MHW` for
  elevations); the `.npy` elevation and dune grids must already be in dam. Full
  datum and conversion record: `UNITS.md`.
- **Changed**: `yes` = differs from default; `no` = same; `period` = depends on the
  calibration/test period; `n/a` = has no effect in the Hatteras setup.

Compiled 2026-10-13 from the code, not from a run's saved output.

## 1. `Cascade(...)` constructor (`cascade/cascade.py`)

| Variable | Group | Default | Hatteras | Changed | Units | Meaning |
|---|---|---|---|---|---|---|
| `datadir` | Files | — | `data/hatteras_init` | yes | — | Folder the inputs resolve against |
| `name` | Files | "default" | run name | yes | — | Run name |
| `elevation_file` | Files | barrier3d-default-elevation.npy | one per domain (topo version resolved) | yes | — | Initial interior topography |
| `dune_file` | Files | barrier3d-default-dunes.npy | one per domain | yes | — | Initial dune heights |
| `parameter_file` | Files | barrier3d-default-parameters.yaml | Hatteras-CASCADE-parameters.yaml (copied per run) | yes | — | Barrier3D parameter file |
| `storm_file` | Files | cascade-default-storms.npy | 1997–2009 or 2010–2025 series | period | — | Storm series |
| `num_cores` | Run | 1 | 1 | no | — | Parallel cores |
| `alongshore_section_count` | Run | 6 | 120 (90 real + 15 buffer each side) | yes | — | Number of 500 m domains |
| `time_step_count` | Run | 200 | run years + 1 (14 for 1996, 17 for 2009) | period | yr | Number of states |
| `alongshore_transport_module` | Modules | True | True | no | — | BRIE alongshore transport |
| `beach_nourishment_module` | Modules | True | per-domain list (community and nourished domains) | yes | — | Nourishment and dune management |
| `roadway_management_module` | Modules | False | per-domain list along NC-12 | yes | — | Road management |
| `community_economics_module` | Modules | False | False | no | — | CHOM |
| `sandbag_management_on` | Modules | False | False | no | — | Sandbag height floor |
| `enable_shoreline_offset` | Modules | False | True | yes | — | Apply a starting shoreline offset |
| `shoreline_offset` | Modules | [] | measured offset per domain | yes | m | Cross-shore start position |
| `wave_height` | Waves | 1 | 2.0 | yes | m | Deepwater Hs |
| `wave_period` | Waves | 7 | 7.5 | yes | s | Peak period |
| `wave_asymmetry` | Waves | 0.8 | 0.6 | yes | — | Fraction of waves from the left, looking onshore |
| `wave_angle_high_fraction` | Waves | 0.2 | 0.5 | yes | — | Fraction of waves at more than 45° |
| `bay_depth` | Environment | 3.0 | 3.0 | no | m | Back-barrier bay depth |
| `s_background` | Environment | 0.001 | 0.001 | no | — | Background slope |
| `berm_elevation` | Environment | 1.9 | 1.7 | yes | m NAVD88 | Berm elevation |
| `MHW` | Environment | 0.46 | 0.36 | yes | m NAVD88 | Mean high water (Duck gauge) |
| `beta` | Environment | 0.04 | 0.04 | no | — | Beach slope for run-up |
| `sea_level_rise_rate` | Environment | 0.004 | 0.002 (1996) / 0.006 (2009) | period | m/yr | RSLR |
| `sea_level_rise_constant` | Environment | True | True | no | — | Linear RSLR |
| `background_erosion` | Environment | 0.0 | per-domain list from the source/sink preset (all 0 under zeroBE) | yes | m/yr | Imposed shoreline change |
| `min_dune_growth_rate` | Dunes | 0.25 | 0.55 | yes | — | Minimum dune growth rate |
| `max_dune_growth_rate` | Dunes | 0.65 | 0.95 | yes | — | Maximum dune growth rate |
| `road_ele` | Road | 1.7 | measured per domain | yes | m MHW | Road elevation |
| `road_width` | Road | 30 | 20 | yes | m | Road width |
| `road_setback` | Road | 30 | measured per domain (setback file) | yes | m | Road setback from the dune line |
| `road_relocation_setback` | Road | None | 20 | yes | m | Setback of a relocated road |
| `dune_design_elevation` | Road | 3.7 | 3.94 (= 4.3 m NAVD88, Velasquez-Montoya et al. 2021; 3.0 until 2026-10-13) | yes | m MHW | Dune rebuild target |
| `dune_minimum_elevation` | Road | 2.2 | 0.0 (raised to berm + 0.3) | yes | m MHW | Dune rebuild trigger |
| `trigger_dune_knockdown` | Road | False | False | no | — | Reset dunes after the road is abandoned |
| `group_roadway_abandonment` | Road | None | None | no | — | Abandon roads in groups |
| `nourishment_interval` | Nourishment | None | None (fills injected by year from the schedule) | no | yr | Years between fills |
| `nourishment_volume` | Nourishment | 300 | 0 at start; historical volumes injected per year | yes | m³/m | Fill volume |
| `overwash_filter` | Nourishment | 40 | 40 in community zones, 0 elsewhere | yes | % | Overwash removed from the interior |
| `overwash_to_dune` | Nourishment | 10 | 9 | yes | % | Overwash moved to the dunes |
| `beach_full_cross_shore` | Nourishment | 70 | 70 | no | m | Target beach width |
| `sandbag_elevation` | Nourishment | 1.5 | 0 | n/a | m above the berm | Dune height sandbags restore (sandbags off) |
| `number_of_communities` | CHOM | 1 | 1 | n/a | — | CHOM communities |
| `sand_cost` | CHOM | 10 | 10 | n/a | $/m³ | Unit cost of sand |
| `taxratio_oceanfront` | CHOM | 1 | 1 | n/a | — | Oceanfront share of the tax burden |
| `external_housing_market_value_oceanfront` | CHOM | 6e5 | 6e5 | n/a | $ | Comparable house outside the system |
| `external_housing_market_value_nonoceanfront` | CHOM | 4e5 | 4e5 | n/a | $ | Same, non-oceanfront |
| `fixed_cost_beach_nourishment` | CHOM | 2e6 | 2e6 | n/a | $ | Fixed cost per fill |
| `fixed_cost_dune_nourishment` | CHOM | 2e5 | 2e5 | n/a | $ | Fixed cost per dune build |
| `nourishment_cost_subsidy` | CHOM | 10e6 | 10e6 | n/a | $ | Subsidy on the nourishment plan |
| `house_footprint_x` | CHOM | 15 | 15 | n/a | m | House footprint, cross-shore |
| `house_footprint_y` | CHOM | 20 | 20 | n/a | m | House footprint, alongshore |

## 2. Fixed inside the BRIE coupler (`cascade/brie_coupler.py`)

| Variable | Default | Hatteras | Changed | Units | Meaning |
|---|---|---|---|---|---|
| `ast_model` | True | True | no | — | BRIE alongshore transport on |
| `barrier_model` | False | False | no | — | BRIE's own overwash off |
| `inlet_model` | False | False | no | — | Inlets off |
| `b3d` | True | True | no | — | Barrier3D overwash on |
| `z` | 10.0 | 10.0 | no | m | Initial sea-level reference |
| `Jmin` | 10000 | 10000 | no | m | Minimum inlet spacing |
| `a0` | 0.5 | 0.5 | no | m | Tide amplitude |
| `marsh_cover` | 0.5 | 0.5 | no | — | Back-barrier marsh fraction |
| `dy` | 500 | 500 | no | m | Domain length |
| `dt` | 1 | 1 | no | yr | Time step |
| `dtsave` | 1 | 1 | no | yr | Save spacing |
| `h_b_crit` | 1.44 (= 1.9 − 0.46) | 1.34 (= 1.7 − 0.36) | yes | m above MHW | Critical barrier height = berm − MHW; sets shoreline diffusivity |

## 3. Barrier3D parameters

"Set by CASCADE" = written into the parameter file at start-up
(`brie_coupler.py:117-185`), so editing the file has no effect. Keys missing from
the Hatteras file take Barrier3D's built-in default
(`Barrier3D/barrier3d/configuration.py`), which is not always the CASCADE
template's value.

| Variable | Group | Default (CASCADE template) | Hatteras | Changed | Set by CASCADE | Meaning |
|---|---|---|---|---|---|---|
| `BarrierLength` | Geometry | 500 | 500 | no | ✓ | Domain length [m] |
| `BayDepth` | Geometry | 3.0 | 3.0 | no | ✓ | Bay depth [m] |
| `DShoreface` | Shoreface | 8.9 | 17.8 (= 8.9 × Hs) | yes | ✓ | Shoreface depth [m] |
| `LShoreface` | Shoreface | 524.4 | computed by BRIE | yes | ✓ | Shoreface length [m] |
| `ShorefaceToe` | Shoreface | 1100 | computed by BRIE | yes | ✓ | Shoreface toe position [m] |
| `s_sf_eq` | Shoreface | 0.0170 | computed by BRIE | yes | ✓ | Equilibrium shoreface slope [–] |
| `k_sf` | Shoreface | 19325 | computed by BRIE | yes | ✓ | Shoreface flux constant [m³/m/yr] |
| `MHW` | Water level | 0.46 | 0.36 | yes | ✓ | Mean high water [m NAVD88] |
| `TideAmp` | Water level | 1.2 | — | n/a | | Not read by Barrier3D |
| `RSLR_Constant` | Water level | true | true | no | ✓ | Linear RSLR |
| `RSLR_const` | Water level | 0.004 | 0.002 / 0.006 | period | ✓ | RSLR [m/yr] |
| `beta` | Storms | 0.04 | 0.04 | no | ✓ | Beach slope for run-up |
| `storm_file` | Storms | default series | period storm series | period | ✓ | Storm series file |
| `StormSeries` | Storms | [] | [] | no | | Storm series held in memory |
| `StormStart` | Storms | 2 | 2 (built-in) | no | | Year storms start |
| `BermEl` | Dunes | 1.9 | 1.7 | yes | ✓ | Berm elevation [m NAVD88] |
| `Dmaxel` | Dunes | 3.4 | 5.5, one value island-wide, both periods (matches the surveyed crest height; top of the 1996-2009 overwash plateau 4.5-5.5, collapse above ~5.6) | yes | | Maximum dune elevation [m NAVD88]; the height every dune grows toward |
| `DuneCeilingFromStart` | Dunes | not in base Barrier3D (local fork only; built-in false) | removed 2026-10-13 (was true 09-28 to 10-13) | no | | Per-cell ceiling = starting crest; dropped to stay on the base model |
| `DuneCeilingFloor` | Dunes | not in base Barrier3D (local fork only; built-in 0.5) | removed 2026-10-13 | no | | Lowest per-cell ceiling, above the berm [m]; unused without per-cell ceilings |
| `Dstart` | Dunes | 0.5 | 0.5 (built-in) | no | | Initial dune height above the berm [m] |
| `DuneWidth` | Dunes | 20 | 20 (built-in) | no | | Dune width [m] |
| `DuneRestart` | Dunes | 0.075 | 0.075 (built-in) | no | | Restart height of a destroyed dune [m] |
| `DuneParamStart` | Dunes | true | true | no | ✓ | Read dune heights from file |
| `DuneParamMultipleRows` | Dunes | not in template (built-in false) | false (built-in) | no | | Multi-row dune file |
| `GrowthParamStart` | Dunes | false | false | no | ✓ | Read growth rates from file |
| `growth_param_file` | Dunes | barrier3d-default-growthparam.npy | built-in name | n/a | | Growth-rate file (unused) |
| `rmin` | Dunes | 0.25 | 0.55 | yes | ✓ | Minimum dune growth rate |
| `rmax` | Dunes | 0.65 | 0.95 | yes | ✓ | Maximum dune growth rate |
| `C1` | Dunes | 8.8 | 8.8 (built-in) | no | | Empirical dune erosion parameter [m] |
| `C2` | Dunes | 4.6 | 4.6 (built-in) | no | | Empirical dune erosion parameter [m] |
| `HdDiffu` | Dunes | 0.75 | 0.75 (built-in) | no | | Largest height step between neighbouring dune cells [m] |
| `dune_file` | Files | domain_1_resampled_dune.npy | per domain | yes | ✓ | Dune file |
| `elevation_file` | Files | domain_1_resampled_topography.npy | per domain | yes | ✓ | Topography file |
| `Rat` | Shoreline | -1.0 | per-domain background erosion | yes | ✓ | Background erosion rate [m/yr] |
| `Cx` | Overwash | 10 | 10 (built-in) | no | | Slope multiplier in flow routing |
| `nn` | Overwash | 0.5 | 0.5 (built-in) | no | | Flow-routing exponent |
| `mm` | Overwash | 2.0 | 2.0 (built-in) | no | | Sediment-transport exponent |
| `MaxUpSlope` | Overwash | 0.25 | 0.25 (built-in) | no | | Steepest upslope flow can climb [m/m] |
| `threshold_in` | Overwash | 0.25 | 0.25 (built-in) | no | | Inundation vs run-up threshold [m³/hr] |
| `Ki` | Overwash | 7.5e-6 | 7.5e-6 (built-in) | no | | Transport coefficient, inundation |
| `Kr` | Overwash | 7.5e-5 | 7.5e-5 (built-in) | no | | Transport coefficient, run-up |
| `Rin_i` | Overwash | 0.1 | **0.25** (built-in) | yes | | Infiltration, inundation |
| `Rin_r` | Overwash | 2.0 | 2.0 (built-in) | no | | Infiltration, run-up |
| `OWss_i` | Overwash | 2 | 2 (built-in) | no | | Overwash substeps, inundation |
| `OWss_r` | Overwash | 1 | 1 (built-in) | no | | Overwash substeps, run-up |
| `Qs_min` | Overwash | 1.0 | 1.0 (built-in) | no | | Minimum discharge for transport [m³/hr] |
| `Qs_bb_min` | Overwash | 1 | 1.0 (built-in) | no | | Minimum sediment flux into the bay [m³/hr] |
| `Cbb_i` | Overwash | 0.85 | **0.8** (built-in) | yes | | Bay sediment-load decay, inundation |
| `Cbb_r` | Overwash | 0.7 | **0.5** (built-in) | yes | | Bay sediment-load decay, run-up |
| `Shrub_ON` | Shrubs | 0 | 0 | no | ✓ | Shrubs off (forced) |
| `Dshrub` | Shrubs | 2.75 | 2.75 (built-in) | n/a | | Minimum fronting-dune elevation for shrubs [m] |
| `MaxShrubHeight` | Shrubs | 5.3 | 5.2 (built-in) | n/a | | Maximum shrub height [m] |
| `ShrubEl_min` / `ShrubEl_max` | Shrubs | 1.2 / 2.3 | 1.2 / 2.3 (built-in) | n/a | | Shrub elevation range [m NAVD88] |
| `Qshrub_max` | Shrubs | 0.15 | 0.15 (built-in) | n/a | | Largest overwash reduction by a shrub cell |
| `BurialLimit` / `UprootLimit` | Shrubs | 0.75 / -0.2 | 0.75 / -0.2 (built-in) | n/a | | Shrub burial limit (Barrier3D labels it m but describes a percentage of height) and erosion limit [m] |
| `SalineLimit` | Shrubs | 5 | 0.05 (built-in) | n/a | | Saline-flooding discharge limit for shrubs [m³/hr] |
| `SprayDist` | Shrubs | 170 | 170 (built-in) | n/a | | Salt-spray distance [m] |
| `TimeFruit` | Shrubs | 5 | 5 (built-in) | n/a | | Age shrubs start fruiting [yr] |
| `Seedmin` / `Seedmax` | Shrubs | 100 / 1000 | 1000 / 10000 (built-in) | n/a | | Seeds per shrub per year |
| `GermRate` / `Female` | Shrubs | 0.6 / 0.5 | 0.6 / 0.5 (built-in) | n/a | | Germination rate, female fraction |
| `disp_mu` / `disp_sigma` | Shrubs | -0.72 / 1.5 | -0.72 / 1.5 (built-in) | n/a | | Seed dispersal distribution |
| `TMAX` | Run | 19 (left over) | run years + 1 | period | ✓ | Run length |
| `SeededRNG` | Run | true | true (built-in) | no | | Fixed random seed |
| `enable_sandbags` | Sandbags | false | false | no | ✓ | Sandbag switch |
| `Sandbag_elevation` | Sandbags | 1.9 | 0 | n/a | ✓ | Sandbag dune height above the berm [m] |

## 4. Hatteras additions (not CASCADE inputs)

Run settings from `scripts/hatteras_ms/hat_run.yaml`, overridable by `HAT_*` variables
(`HAT_hindcast_config.py`).

| Variable | Code default | Meaning |
|---|---|---|
| `start_year` | 1996 | Period key (DEM year) |
| `source_sink_preset` | zeroBE | Background erosion preset |
| `scenario` | full_management | Management scenario |
| `relocations` | None | Historical road relocations on/off |
| `offset_mode` | metres | How the island offset is built |
| `geometry` | base | Reach (90 real domains) |
| `groin_enabled` | False | GroinCallback on/off |
| `groin_strength_1996` | 0.16 | Share of the arriving drift the groin holds back, 1996 start |
| `groin_strength_2009` | 0.10 | Share of the arriving drift the groin holds back, 2009 start |
| `groin_strength` | None | Overrides both when set; needed for any other start year (0.25 = whole-life strength) |
| `groin_end_year` | 2017 | First year the groin no longer works (the 2017 Buxton fill) |
| `relocation_setback_m` | 20.0 | Relocated-road setback [m] |
| `sandbags` | False | Sandbag placement |
