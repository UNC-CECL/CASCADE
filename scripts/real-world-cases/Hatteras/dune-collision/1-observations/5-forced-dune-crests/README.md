# 5-forced-dune-crests — do the lidar's dunes make the shoreline erosion appear?

`HAT_forced_dune_crests.py`. A diagnosis, not a model change: the evidence that the missing dune lowering matters for the shoreline. Runs and logs: `output/raw_runs/dune-storm-response/2026-10-26-forced-dune-crests/`, which also holds the control run the other collision studies compare against. Figures, with the tables and captions under `supporting/`: `D:\research\UNC\hatteras_ch1\model_runs\dune-collision\1-observations\5-forced-dune-crests\` (moved there 2026-11-03; the script writes them there).

Asked (Hannah, 2026-10-26) because the managed edgeBE 1996-2009 run misses the
CoastSat erosion at GIS ~8-15. Dennis (1999) runs up to 3.65 m NAVD88, the
modelled dunes there stand 5-7 m, and Barrier3D lowers a dune only where a
storm overtops it, so Dennis leaves them untouched while the lidar shows
GIS 10-12 dropping to ~1.7 m. Datums (MHW 0.36, berm 1.7 m NAVD88) and the
Duck gauge were checked and are not the cause: the Cape Hatteras Pier gauge
peaked lower than Duck during Dennis.

In-process, before the given storm year, each cell's two dune rows are scaled
so the crest equals the survey crest (`survey_dune_crests.csv`; beach landward
of the window counts as the berm). Managed edgeBE, groin, 1996-2009.
Variants: control, buxton (GIS 8-12, 1999), buxton_pre (8-12, 1998),
buxton_both (8-12, 1998 and 1999), buxton_8_15 (8-15, 1999), island (1-90, 1999),
events (1-90; 1998, 2000, 2004), buxton_events (8-15; 1998, 2000, 2004),
all_surveys (1-90; every survey).

    python HAT_forced_dune_crests.py run [variants...]
    python HAT_forced_dune_crests.py score
    python HAT_forced_dune_crests.py plot     # position change vs CoastSat
    python HAT_forced_dune_crests.py crests   # control vs lidar, GIS 8-15

Results 2026-10-26 (figures in the workspace's `dune-collision/1-observations/5-forced-dune-crests/`):
the control reproduces the matrix run (interior bias -0.445 m/yr). Forcing the
1999 crests at GIS 8-12 or 8-15 makes the erosion appear: GIS 11-12 match
CoastSat, GIS 8-10 overshoot (GIS 10 -4.0 against ~-0.9 m/yr). 1998 forcing
alone barely changes anything. GIS 13-15 do not move: their lidar crests
(5.7-6.9 m) stand above every storm's run-up, so the erosion there needs
another mechanism (background or alongshore loss). Interior bias: buxton
-0.538, buxton_8_15 -0.538, island -0.689 m/yr. Forcing every event makes it
worse (events -0.865, all_surveys -0.931; GIS 9-12 overshoot 2-4x).

The modelled and lidar crests agree through 1997 and part in 1998-1999,
when Bonnie and Dennis lowered the real dunes and not the modelled ones.
Island-wide, per run: `cascade_pipeline/plotting/dune_crest_comparison.py`
(drawn by `run_figures.py` into each run's `figures/dune_crests/`); in the
managed 1996-2009 run 17-21 domains end more than 1 m above the lidar after
1999. Not adopted: a diagnosis, not a model change.
