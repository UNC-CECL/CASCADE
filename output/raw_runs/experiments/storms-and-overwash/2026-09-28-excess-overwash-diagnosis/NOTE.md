# Why the model overwashes more than the imagery shows (2026-09-28)

**Question (Hannah).** In the storm-length selection (stage 1), every storm series overwashed far more domains than the observed record in most image windows. Why?

It was investigated read-only, on existing runs and inputs (driver `scripts/hatteras_ms/experiments/HAT_excess_overwash_diagnosis.py`, plus the checks recorded below). The managed runs used are `drop72` (the matrix) and `trim24`, for both windows.

## Answer: the model's dunes are held at about 3 m, but Hatteras's are about 5 m

**1. The dune height ceiling is Barrier3D's Virginia default.**

- `Dmaxel`, the maximum dune elevation, is set in no Hatteras parameter file. Every run therefore uses Barrier3D's default of 3.4 m NAVD88 (`configuration.py`), which is 3.04 m MHW.
- Barrier3D grows dunes logistically towards that ceiling: G = r·D·(1 − D/Dmax). A dune taller than Dmax gets negative growth, so it **shrinks every year, storms or not**.
- The 2009 lidar puts the island-median foredune crest at 4.9 m MHW (10th–90th percentile 2.8–7.1 m).

**2. The NC-12 dune rebuild cuts tall dunes down to the design height.** When any front-row dune cell falls below the rebuild trigger (1.64 m MHW), `roadway_manager.rebuild_dunes` replaces the **whole** dune field with the design height (`DUNE_DESIGN_ELEVATION_M` = 3.0 m MHW). That includes cells that were 5–7 m. Its own docstring cites 4.3 m NAVD88 as the crest NC-12 needs to be safe from overwash.

**Evidence: island-median dune crest (m MHW), 1996–2010 runs, by year**

| run | 1996 | 1997 | 1998 | … | 2006 | 2007 | 2010 |
|---|---|---|---|---|---|---|---|
| natural | 5.5 | 3.9 | 3.6 | … | 3.1 | 2.3 | 1.4 |
| managed | 5.5 | 3.1 | 3.1 | … | 3.0 | 2.4 | 3.0 |

- **The first-year drop is not the storms.** 1996's largest storm reached Rhigh 2.35 m.
- **The managed run sits at exactly the 3.0 m design height** for the rest of the window.
- **Against the lidar at the same place and time:** the modelled 2010 crest is lower than the 2009 lidar crest (the 2010–2024 run's starting dunes) by a median 2.7 m (managed) and 3.4 m (natural). The model is more than 1 m lower in 75 of 90 domains (managed) and 86 of 90 (natural).

**The false alarms follow from this.** When the model overwashes where the imagery shows none, the storm's pre-storm lowest crest is about 1.45 m MHW (median), barely above the 1.34 m berm, and the storm clears even the average crest by 0.6–1.0 m. Storms of Rhigh 2.5–3.8 m, which a 5 m dune would stop, overtop the model's 3 m dune line nearly everywhere (`tables/cells.csv`).

## What was ruled out, or is secondary

- **Clipped dune search windows (2004-start v1).** These are not the cause. Domains whose true crest sits more than 1 m behind the dune row have *fewer* false alarms in 2010–2024 (24%) than the rest (46%). The clipping is still an input defect worth fixing (`tables/hidden_crest.csv`: 12 domains with more than 1 m hidden in 2010–2024, 1 in 1996–2010).
- **Overwash too small to see.** This is part of the story in 1996–2010: hits carry about 37 m³/m (trim24) against about 10 m³/m for false alarms, and requiring more than 20 m³/m lifts the Peirce skill from 0.37 to 0.61. It is not the story in 2010–2024, where hits and false alarms have similar volumes and no threshold helps.
- **Storm levels.** These were not tested directly. The margins above are measured against the model's lowered crests; against the lidar crests most of these storms would not overtop.

**Status: diagnosis complete. No change made.** The next step is Hannah's choice: an experiment that sets `Dmaxel` for Hatteras, and changes the rebuild so it does not lower dunes taller than the design height. That changes what the model ingests.

Files: `tables/hidden_crest.csv`, `tables/cells.csv`.
