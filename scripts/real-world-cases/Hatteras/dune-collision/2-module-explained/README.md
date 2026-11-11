# 2-module-explained — how DuneCollisionManager works

Figures that explain the dune collision module (`cascade/dune_collision.py`) step by step, as
`../../hard-structures/groin/2-module-explained/` does for the groin. Nothing here is a new test or a fit:
`collision_explained_figures.py` puts one example dune cell (a Buxton dune before Dennis: toe 2.4 m, crest
7.3 m NAVD88) through the module's own functions (`beach_factor`, `collision_cells`, `larson_volume`,
`lower_dunes`) at the hindcast's settings (Cs 1.02e-4, W0 25 m on the dry beach, beach start to toe; since 2026-11-03,
when the widths moved from the crest to the toe and Cs from 1.24e-4 and W0 from 45 m; toes to the first dune
since 15:50 that day, the figures reading `dune_toe_cells.csv` redrawn on them), with storms from the 1997-2009 hindcast
storm series. Runs in seconds.

    python collision_explained_figures.py

The figures are written to the workspace twin of this folder,
`D:\research\UNC\hatteras_ch1\model_runs\dune-collision\2-module-explained\figures\`; each figure's
Tests / How to read / Shows caption is in `figures/supporting/CAPTIONS.md` there, with the PDFs.

| figure | question |
|---|---|
| `collision_1_regimes.png` | Which storms does the module act on, and how far does it lower a dune? (a) the three regimes on a dune profile: below the toe nothing, toe to crest the module, above the crest Barrier3D's own overwash; the sand one Dennis-sized storm removes. (b) the crest after one storm against its water level: zero below the toe, growing with the square of the water above it, floored at the toe (Dennis, no beach: -1.02 m) |
| `collision_2_beach.png` | What does the beach in front of the dune do? (a) the beach factor max(0, 1 - W / 25 m), W the dry beach; (b) the same storm behind 0, 8, 16 and 25 m dry beaches (Dennis: -1.02, -0.69, -0.37, 0 m): the beach changes how much a storm takes, not which storms reach the dune, and nothing is lost behind 25 m or more |
| `collision_3_in_the_model_year.png` | Where does the module act in Barrier3D's year? (a) after dune growth it keeps the pre-storm crest; after each storm's overwash lowering, before dune diffusion, it lowers the cells that storm left in collision. (b) the 1999 storms on the example cell behind a 10 m dry beach: 6 of 12 reach the toe, Dennis takes 0.61 m, the crest goes 7.30 -> 6.61 m over the year, later storms see the lowered dune |
| `collision_4_year_of_storms.gif` | The 1999 storms on the real Buxton dunes, storm by storm (`collision_year_gif.py`): every 10 m cell of GIS 8-15 starts at its 1998-09 lidar crest, with the toes and beach widths the hindcast reads; each storm shows the water arriving, the cells it puts in collision (blue ticks) or overtops (grey ticks, Barrier3D's), then the lowered dunes. Dennis puts 286 of 400 cells in collision; by December 313 cells are lower, 0.64 m on average (up to 4.6 m, dunes taken to their toes), most behind narrow beaches and none behind 25 m or more of dry beach. Its last frame is `collision_4_year_of_storms.png`. Dune growth is left out |
| `collision_5_one_domain.gif` | The simple version: the 1999 storms on the 50 dunes of GIS 12 alone, as bars from the berm to the crest (`collision_domain_gif.py`; `--gis` and `--year` for another). GIS 12 has a full lidar record, a 0-20 m dry beach everywhere and no 1999 overtopping, so all it shows is collision. Dennis takes 0.89 m on average, the year 1.07 m (the lidar's 1998-1999 drop on the same cells: 1.36 m), most from the low dunes with low toes. Last frame `collision_5_one_domain.png` |
| `collision_6_full_model.gif` | GIS 12 through the full model, 1997-2009, step by step (`collision_full_model_gif.py record`, one ~3 min hindcast run as `hat_run.yaml` sets it, recording the domain's 50 dunes at every step; then `gif`). Each year: Barrier3D's dune growth (green), each storm that changes the dunes with what overwash (grey, Barrier3D) and collision (purple, the module) took, Barrier3D's dune migration when the shoreline retreats a full cell (teal, end of 2003 and 2008), and the road manager (orange: overwash cleared off NC 12 onto the dune, and rebuilds of every dune to 4.35 m at the end of 2000, 2004, 2008 and 2009). On the final calibration (Cs 1.02e-4, W0 25 m dry beach; re-recorded 2026-11-03), per cell over the run: growth +0.66 m, overwash -2.61, collision -3.86, migration +1.40, management +2.42; median crest 6.54 -> 4.34 m. Dennis is all collision (0.81 m), Isabel (5.58 m) all overwash (1.17 m). The legend keeps every entry in place and fades those not acting; the lower panel colours each step by what moved the median; steps that move a dune 1 m or more are held. 2 frames a second, about 90 s. Recorded states: `data/gis12_dune_states.csv` in the workspace twin |
| `collision_7_one_storm_profile.gif` | One storm on one dune in cross-section (`collision_storm_profile_gif.py`; `--width` for another beach): the example cell behind a 10 m dry beach through Dennis. The water rises past the toe, holds at 3.65 m for 116 h, the module takes Larson's volume times the beach factor (0.60) off both dune rows, 12.3 m³ per m of shore, crest 7.30 -> 6.69 m; the sand leaves the model. Last frame `collision_7_one_storm_profile.png` |
| `collision_8_dennis_island.png` | One storm along the whole island (`collision_dennis_island.py`; `--year`): every dune cell of GIS 1-90 at its 1998 lidar crest against Dennis. (a) each domain's share of cells overtopped (Barrier3D's overwash), in collision and lowered by the module, in collision behind a dry beach of 25 m or more (nothing taken), or below the toe; (b) the module's mean crest lowering per domain from this storm alone. 33% of cells are lowered, 55% are in collision behind wide beaches, 10% overtopped, 2% below the toe; the lowering sits at Buxton (GIS 8-13, up to 1.1 m), GIS 51-55 (up to 1.1 m) and Rodanthe (GIS 77-89, up to 1.7 m), the three places the lidar shows Dennis lowering the dunes |

**Since these figures (2026-11-04 to 11), they are out of date in three ways:** the beach factor they show was
removed (2026-11-04); the toe is now the berm, 1.7 m NAVD88, not the lidar toe, with Cs about 1e-6 to 2e-6
(2026-11-10); and the module has options it did not: `rows` ("front", "equal", "lateral": the face retreating through
the front row, which keeps its height until its sand is gone) and per-cell tracking of the eroded volume and of
dune cell collapse (`../4-hindcast/TESTING_2026-11-09_to_11.md`, section 10). Also worth knowing when reading the
dune rows here: Barrier3D's two dune rows start equal but drift apart, because it grows the back row at a third of
the front row's rate (front minus back, median +0.59 m by 2009 with the module off); the module and the scores use
the higher row.

Run `python collision_year_gif.py --year 2005 --gis 80 89` for another year or stretch (the year before needs a
lidar survey: 1996-2001, 2004, 2005; storms run to 2009). Still to add: the same on the model's own crests from a hindcast run.
