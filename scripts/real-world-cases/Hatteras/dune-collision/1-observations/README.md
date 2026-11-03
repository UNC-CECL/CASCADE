# 1-observations — the dune loss the model misses

True whatever the model does. Measured from the 1996-2009 lidar surveys (`survey_dune_crests.csv`) and
the storm series, by the measurement actions of
`../4-hindcast/1-calibration-1996-2009/HAT_dune_collision.py` (to move here when it is split):

| Action | What it measures | Writes (under `output/raw_runs/dune-storm-response/2026-10-26-dune-collision/`) |
|---|---|---|
| `toes` | the dune toe per survey and 10 m cell (now also `../../input_prep/1-barrier3d-domains/3-dune-toes/HAT_dune_toes.py`, which the runner reads) | `dune_toes.csv` |
| `twl` | the hourly total water level the storm series is built from | `twl_hourly_1996_2009.csv` |
| `calibrate` | crest change against Larson exposure on the cells never overtopped | `collision_calibration*.csv/png` |
| `empirical` | crest loss over pre-storm crest against R_high/D_high and against beach width | `collision_empirical*.csv/png` |

The motivation figures (2026-11-03) go to
`D:\research\UNC\hatteras_ch1\model_runs\dune-collision\1-observations\collision-motivation\`, one folder per topic,
each with its PNGs and a `supporting/` (PDFs, `CAPTIONS.md`, the numbers behind each panel):

| Folder | Figures | Script |
|---|---|---|
| `1-mechanism` | `collision_motivation.png`, `two_dunes.png` | `HAT_collision_motivation.py profiles`, `plot`, `two_dunes` |
| `2-mismatch-along-island` | `collision_offset_heatmap.png`, `collision_alongshore_crests.png` | `HAT_collision_motivation_alongshore.py` |
| `3-hotspots` | `collision_hotspot_crests.png` (and `collision_hotspots.csv`) | same |
| `4-drivers` | `collision_exposure_offset.png`, `collision_beach_offset.png` | same |

`collision_motivation.png` has four panels,
against the calibration control (`cal_control_g16`, module off): (a) the 1998 and 1999 lidar at GIS 12 cell 22,
Dennis's peak water level on the dune face, the crest down 2.6 m while the model's does not move; (b) GIS 10-12
crests, 1996-2009; (c) the share of domains whose modelled crest is more than 1 m above the lidar's:
0% in 1996, 3% in 1997, 8% after Bonnie, 19% after Dennis, 27% by 2009; (d) crest loss in collision
against the module's beach width (per-cell median, 1998 left out: the 1998 grid's own widths put most
Buxton cells behind beaches over 60 m and reverse the pattern there).

Where along the island (2026-11-03): `HAT_collision_motivation_alongshore.py`, five more figures (folders 2-4),
each with a calibration (1996-2009) and a test (2009-2025) part against the module-off controls. They are
the offset heatmap (domain x survey), the crests at the four hotspots (GIS 7-13, 22-25, 75-78, 84-88), the crests
along the island in three surveys per period, and collision-storm exposure and dry beach width each above the
crest offset. Domains whose dunes saw more collision storms are too high in the model (Spearman r 0.51 / 0.44);
domains behind narrower beaches are too (r -0.68 / -0.44). These are with the dune toes regenerated on
2026-11-03 15:50 (the 2009-start toes now over all eight surveys); on the earlier toes they were 0.41 / 0.47
and -0.64 / -0.40. The 2009-start lidar now includes 2011, 2012 and
2020 (see that start's `survey-crests/README.md`).

Whether the lidar toes sit on the dune traced from the aerial imagery (2026-11-03): `2-toes-vs-dunelines/`.
The traced line is on the dune face, a median 7-17 m landward of the lidar toe; GIS ~32-36 and ~63-67 disagree.

What they show (from the 2026-10-27/28 record, `../4-hindcast/1-calibration-1996-2009/README.md`):

- The lidar dunes drop at GIS 2-3, 10-12 and 84-86 in storms that reach the dune face but not the crest.
- Long, de Bakker & Plant (2014, GRL 41, 2899-2906) give no collision equation: the crest changes by
  less than ~10%, up to ~50% just below overtopping.
- In the Hatteras lidar, crest loss in collision is set by beach width, not R_high/D_high:
  dz/D_high = 0.50 max(0, 1 - W / 45 m), fit on Dennis and 2004-05 (r 0.48 cells, 0.63 domain medians;
  R_high/D_high alone 0.29). Figure `collision_empirical.png`.
- Forcing the lidar crests into the model makes the missing Buxton shoreline erosion appear
  (`1-forced-dune-crests/`).
