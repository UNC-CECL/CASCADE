# 3-on-vs-off — shoreline position change with the collision module on and off

The scripts stay here; their figures are in the workspace's `dune-collision/4-hindcast/` (D:), one folder per
question since 2026-11-12, each with its own `supporting/` (CSV, PDF, CAPTIONS.md). The full index (folder, question,
testing-log section, script) is that folder's `README.md`. Each script names its folder in `OUT` or the figure name:

| Folder | Script |
|---|---|
| `1-module-on-off/` | `HAT_collision_on_off_position.py onoff`, `HAT_collision_crest_compare.py`, `HAT_collision_missed_drops.py`, `HAT_collision_toe_cs_figures.py` |
| `2-model-steps/` | `HAT_collision_model_steps.py`, `HAT_collision_toe_cs_figures.py` |
| `3-collision-design/` (`dune-rows/`, `dune-search-window/`, `logs/`) | `HAT_collision_row_mode_sweep.py`, `HAT_dune_rows_test.py`, `HAT_overwash_compare.py`, `HAT_collision_on_off_position.py maxmed / maxmedlat / win8` |
| `4-cs-and-toes/` | `HAT_collision_cs_sweep.py` |
| `5-road/` (`road-rule/`, `road-relocation-1999/`, `bulldoze-fix/`) | `HAT_road_rule_*.py`, `HAT_road_relocation_test.py` |
| `6-missed-drops/` (`1-zones/`, `2-by-design/`, `3-compared/`) | `HAT_be_hotspot_zones.py`, `HAT_collision_on_off_position.py behot` |
| `7-storm-examples/` | `HAT_collision_buxton_dennis.py`, `HAT_collision_isabel.py` |
| `8-other-forcings/` | `HAT_collision_on_off_position.py` (other modes) |
| `_superseded/` | Automorph-toe, before-road-fix and beach-factor figures; no script writes them any more (`HAT_collision_toe_cs_figures.py` still names the old all-toes and domain-toes figures there) |

What was tested and found, with the figure names: `../TESTING_2026-11-09_to_11.md`.

Most of these were drawn by session scripts from runs already made (nothing in them is refit); the crest maps go
through `HAT_collision_crest_compare.py` and the position plots through `HAT_dune_collision.py plot`.

Since 2026-11-04 the main figures here are the model-change steps, `HAT_collision_model_steps.py` (Larson alone,
Cs 1.9e-5, since the evening of 2026-11-04; the beach-factor version is in the workspace's
`dune-collision/4-hindcast/_superseded/beach-factor-steps-until-2026-11-04/`):
`model_steps_position_1996_2009.png`, `model_steps_position_2009_2025.png` (each step's shoreline change against
CoastSat) and `model_steps_fit.png` (shoreline and crest RMSE per step, island and zones; table
`supporting/model_steps_fit.csv`, the earlier toe and width versions beside it). The steps and what they show:
`../../README.md`, Settings. `HAT_collision_on_off_position.py` also draws, from runs already made: `village`
(villages 40% against 90%), `toes` and `toesnow` (start against median toes), `widths`
(`beach_widths_vs_domain_factor_position_<window>.png`: per-cell beach widths against one beach factor per
domain), `upkeep`, `sand`, `hs`. The tables below are from 2026-11-03, before Rodanthe 90% and the domain toes.

`HAT_collision_cs_sweep.py` (2026-11-04): Larson alone at Cs 1e-5, 1.9e-5, 3e-5, 5e-5 and 8e-5, both windows
(runs `cal_la_cs*`, `cal_step6_toeavg_la`; `la_cs*`, `nocape_rd198_dry_toeavg_la`), against the base and the beach
factor: `cs_sweep_fit.png` (island RMSE against Cs), `cs_sweep_position_<window>.png` (all five) and
`cs_pair_position_<window>.png` (1e-5 against 1.9e-5); table `supporting/cs_sweep_fit.csv`. More Cs is worse in both
windows; 1.9e-5 kept (`../../README.md`, Settings). Island RMSE, shoreline (m/yr) / crest (m):

| Cs | Calibration | Test |
|---|---|---|
| 1e-5 | 1.18 / 0.89 | 2.08 / 1.56 |
| 1.9e-5 | 1.21 / 0.90 | 2.16 / 1.61 |
| 3e-5 | 1.23 / 0.96 | 2.29 / 1.78 |
| 5e-5 | 1.30 / 1.14 | 2.51 / 2.05 |
| 8e-5 | 1.37 / 1.50 | 2.77 / 2.42 |
| module off (base) | 1.14 / 0.91 | 2.24 / 1.41 |
| beach factor, Cs 1.02e-4 | 1.19 / 0.87 | 2.10 / 1.50 |

`HAT_collision_on_off_position.py` (2026-11-03). One figure per window, drawn as the run figures draw position change,
from runs already made (nothing rerun here):

| Window | Module off | Module on (dry beach, 25 m, Cs 1.02e-4, first-dune toes, off at GIS 1-5, 198 m road rule) |
|---|---|---|
| 1996-2009 calibration | `cal_control_g16` | `cal_g16_nocape_rd198_dry` |
| 2009-2025 test | `control_rd198` | `nocape_rd198_dry` |

The 2009-2025 module-off run has the 198 m road-dune rule too (`control_rd198`, added 2026-11-03), so there the
difference is the module alone; the 1996-2009 module-off run does not have it.

Why the module makes the test window worse (2026-11-03, `control_rd198` against `nocape_rd198_dry`):

| Zone | Dune minus lidar, off / on (m) | Net overwash, off / on (m3/m) | Shoreline LRR RMSE, off / on (m/yr) |
|---|---|---|---|
| Buxton-Avon road | -0.38 / -0.52 | 36 / 140 | 1.86 / 2.06 |
| Avon | +0.43 / -0.05 | 38 / 263 | 1.62 / 2.19 |
| Avon-Salvo road | -0.10 / -0.42 | 122 / 162 | 2.31 / 2.35 |
| Tri-Village | +0.23 / -0.08 | 236 / 317 | 2.38 / 2.69 |
| Rodanthe north | +0.87 / +0.14 | 156 / 416 | 0.80 / 0.53 |

Island: 2.19 off, 2.36 on. The road rule alone improves on the plain control (2.28). The module only moves the
shoreline landward (mean -6.7 m, nowhere seaward): it lowers the dunes, toward the lidar at Avon, the Tri-Village and
Rodanthe, and the lower dunes then overwash 1.3-7 times more, which rolls the shoreline back. Where the module-off
dunes were too tall (Avon +0.43, Rodanthe +0.87 m) they were holding the shoreline in place; with the dunes right,
the missing sand (CoastSat accretes at Avon and the Tri-Village, which no Hs reproduces) shows. The sand the module
erodes is also deleted, not put on the beach (open to-do), which adds to the loss.

Figures, in the workspace's `dune-collision/4-hindcast/1-module-on-off/`: `collision_on_off_position_1996_2009.png`
and `collision_on_off_position_2009_2025.png` (captions under `supporting/`).

What they show (2026-11-03):

- 1996-2009: the module adds about 10 m of retreat at Buxton (GIS 8-14), toward CoastSat, and 5-10 m at Rodanthe
  (GIS 76-88), toward CoastSat north of GIS 83 but past it at GIS 79-82. Elsewhere on and off overlap.
- 2009-2025: better at Buxton (GIS 7-14: more retreat south of 10, more accretion north of it, both toward
  CoastSat). Worse at Avon (GIS 22-31, the shoreline cut back where CoastSat advanced) and at Rodanthe (GIS 77-82,
  about 10 m more retreat where the model already over-erodes by 30-40 m). These two stretches are why the module
  scores worse island-wide in the test.

Eroded sand onto the beach (2026-11-03, `nocape_rd198_dry_sand` in HAT_collision_test_2009.py, an experiment patch
run on joblib threads; figure `collision_sand_to_beach_position_2009_2025.png`): shoreline LRR RMSE 2.32, against
2.36 with the sand deleted and 2.19 with the module off. Putting the sand back recovers about a fifth of the
module's penalty (mean on-minus-off -5.4 m against -6.7 m); the rest is the extra overwash over the lowered dunes.
