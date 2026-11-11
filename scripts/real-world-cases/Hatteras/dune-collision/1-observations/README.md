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

The motivation figures (2026-11-03) go to the first four topic folders of
`D:\research\UNC\hatteras_ch1\model_runs\dune-collision\1-observations\`, in reading order (the problem, where it is worst, what happens there, what sets it;
then `5-forced-dune-crests`, why it matters, and `6-toes-vs-dunelines`, a check on the toe input),
each with its PNGs and a `supporting/` (PDFs, `CAPTIONS.md`, the numbers behind each panel):

| Folder | Figures | Script |
|---|---|---|
| `1-mismatch-along-island` | `collision_offset_heatmap.png`, `collision_alongshore_crests.png` | `HAT_collision_motivation_alongshore.py` |
| `2-hotspots` | `collision_hotspot_crests.png` (and `collision_hotspots.csv`); `collision_hotspot_classification.png`, `collision_hotspot_storms.png`, `collision_hotspot_profiles.png` | same; `HAT_collision_hotspot_diagnostics.py` |
| `3-mechanism` | `collision_motivation.png`, `two_dunes.png` | `HAT_collision_motivation.py profiles`, `plot`, `two_dunes` |
| `4-storms-and-beach-width` | `collision_drivers_scatter.png` | `HAT_collision_motivation_alongshore.py` |

`collision_motivation.png` has four panels,
against the calibration control (`cal_control_g16`, module off): (a) the 1998 and 1999 lidar at GIS 12 cell 22,
Dennis's peak water level on the dune face, the crest down 2.6 m while the model's does not move; (b) GIS 10-12
crests, 1996-2009; (c) the share of domains whose modelled crest is more than 1 m above the lidar's:
0% in 1996, 3% in 1997, 8% after Bonnie, 19% after Dennis, 27% by 2009; (d) crest loss in collision
against the module's beach width (per-cell median, 1998 left out: the 1998 grid's own widths put most
Buxton cells behind beaches over 60 m and reverse the pattern there).

Where along the island (2026-11-03): `HAT_collision_motivation_alongshore.py`, five more figures (folders 1, 2 and 4),
each with a calibration (1996-2009) and a test (2009-2025) part against the module-off controls. They are
the offset heatmap (domain x survey, the hotspots boxed), the crests at the collision hotspots, the crests
along the island at each period's end (and in every survey, the `_all` supplements), and one scatter of each
domain's crest offset against its collision-storm count and its dry beach width (was two bar figures in `4-drivers/`
until 2026-11-03).

Collision hotspots (since 2026-11-03), found in each period on its own: stretches of at least two domains where the
model is more than 1 m too high in at least half the surveys after the start (1 m on average) AND the lidar crest
fell at least 1 m from the period's first survey to its last while the model's fell at least 1 m less, the crest
lowering the module is there to supply. Lidar crests in the dune window moved landward by the observed shoreline
retreat, never seaward (the crest tables of 2026-11-03 22:09; `HAT_survey_dune_crests.py`). Calibration: GIS 12-13
(Buxton); test: GIS 23-25 (Avon); missed lowering 3.2 and 2.2 m. Domains meeting the rule alone (no neighbour): see
`2-hotspots/supporting/collision_hotspots_by_domain.csv`. Too high without that missed loss, so not the module's to
fix (dune growth or rebuilding too fast, or a start mismatch): calibration GIS 3, 10-11, 55; test GIS 7-9, 22, 39,
47, 60, 84. Both hotspots hold on the 1 m lidar: the 10 m block-centre crests read 0.4-0.6 m (median) below the
1 m crests in the same window, but the 1 m crest still falls 1.1-1.8 m at GIS 12-13 (1996-2009) and 1.5-1.9 m at
GIS 23-25 (2009-2020).

How the list changed with the dune window (2026-11-03): fixed start windows gave calibration GIS 10-13, 51-53, 84-88
and test GIS 22-25, 75-76, but counted dune retreat as lowering (in the last survey 67-91% of those cells had the
crest on the window's landward edge). A window that followed the crest landward found back dunes already standing
at the start (60% of moved cells; 98% at Avon). The shoreline window, adopted, reads the crest where Barrier3D's
dune row sits, a fixed distance behind its shoreline: at GIS 51-53, 84-88 and 75-76 the dune moved landward with
the shoreline and kept its height, so no missed lowering remains there.

`HAT_collision_hotspot_diagnostics.py` adds three figures to `2-hotspots`: each domain's mean offset against its
missed lowering (the rule at work); lidar and model crest change in each survey interval with the largest storm
between (by R_high in the model's storm series); and a lidar profile at a representative cell of each hotspot, start
against last survey, in the crest script's straightened frame with the start's and the moved dune window shaded.
Domains whose dunes saw more collision storms are too high in the model (Spearman r 0.40 / 0.40); domains behind
narrower beaches are too (r -0.56 / -0.43), with the shoreline-window crests and the toes of 22:16-22:18. The
2009-start lidar now includes 2011, 2012 and
2020 (see that start's `survey-crests/README.md`).

Whether the lidar toes sit on the dune traced from the aerial imagery (2026-11-03): `6-toes-vs-dunelines/`
(figures in the workspace's `6-toes-vs-dunelines/`).
The traced line is on the dune face, a median 7-17 m landward of the lidar toe; GIS ~32-36 and ~63-67 disagree.

Rodanthe village, 2009-2014 (2026-11-03): `7-rodanthe-dunes-2009-2014/` (figure in the workspace's
`7-rodanthe-dunes-2009-2014/`). Asked while finding why the 2009-2025 model over-erodes the village shoreline
(GIS 77-81): the 2009, 2011 and 2014 lidar on one datum beside the 2008, 2009 and 2014 imagery. The crests barely
change (no rebuild after 2009; the model starts from the right dunes); the beach at GIS 79-81 is cut back 20-40 m
in the post-Irene 2011 survey and partly back by 2014. It also showed that the 1996-start table had read the 2009
lidar in the 1996 dune windows (crests at the berm at GIS 78-80): fixed the same day with SURVEY_WINDOWS_FROM in
`input_prep/1-barrier3d-domains/2-survey-crests/HAT_survey_dune_crests.py`.

What they show (from the 2026-10-27/28 record, `../4-hindcast/1-calibration-1996-2009/README.md`):

- The lidar dunes drop at GIS 2-3, 10-12 and 84-86 in storms that reach the dune face but not the crest.
- Long, de Bakker & Plant (2014, GRL 41, 2899-2906) give no collision equation: the crest changes by
  less than ~10%, up to ~50% just below overtopping.
- In the Hatteras lidar, crest loss in collision is set by beach width, not R_high/D_high:
  dz/D_high = 0.50 max(0, 1 - W / 45 m), fit on Dennis and 2004-05 (r 0.48 cells, 0.63 domain medians;
  R_high/D_high alone 0.29). Figure `collision_empirical.png`.
- Forcing the lidar crests into the model makes the missing Buxton shoreline erosion appear
  (`5-forced-dune-crests/`; figures in the workspace's `5-forced-dune-crests/`).

## 8-buxton-dune-change (2026-11-11)

`8-buxton-dune-change/HAT_buxton_dune_change.py`: every 1996-2009 lidar survey of GIS 10-15 on the 2009 datum (the
model's vertical offsets, then each survey's remaining NC-12 difference from 2009; on the road at GIS 15 every survey is
within 0.05 m of 2009, `supporting/buxton_datum_check.csv`), on the model clips' 1 m grid: maps by survey, change between
surveys, and one cross-shore profile per domain. Asked because the collision module cannot lower GIS 10-12 in time
(`../4-hindcast/TESTING_2026-11-09_to_11.md`, section 4).

Found: at GIS 10-12 the original foredune was cut down and pushed landward in 1997-1999 (GIS 11 already by 1998;
1998-1999 the strongest loss, along the dune front), ending at about the Dennis peak water level (3.6-4 m); then a new,
taller dune (6-7 m) was built about 50-80 m landward, beside the relocated NC-12, from 2001 to 2004. GIS 13 kept its
front dune (about 6 m) and gained a second, lower ridge behind it; GIS 14-15 held. So the lidar "crest" at GIS 10-12
is one dune until 1999 and a different, built dune after 2001: the model, with one dune line and no road relocation,
cannot follow either the loss or the rebuild. The 2005 survey does not cover GIS 10-11.
