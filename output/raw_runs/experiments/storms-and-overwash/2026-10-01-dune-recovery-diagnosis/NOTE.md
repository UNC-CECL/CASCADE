# Why breached dunes never rebuild (2026-10-01)

**Question (Hannah: "look into why the breached dunes never rebuild").** The storms-vs-overwash figure (`data/hatteras_init/8-overwash-analysis/4-vs-model/storms_vs_overwash_1996_2024.png`) showed false-alarm hotspots at GIS 1–6 and GIS 78–84. In those domains the model's lowest dune crest sits at the berm for most of each run, although it started 1.4–3.3 m MHW.

**Setup.** Read-only, on the four adopted matrix runs (1996/2010, managed/natural). There are no new runs. Driver: `scripts/hatteras_ms/experiments/HAT_dune_recovery_diagnosis.py`. Table: `tables/dune_recovery_by_domain.csv`.

- A cell counts as **flattened** when its crest is within 0.5 m of the berm.
- The **typical storm** is the median of each year's largest storm: 2.75 m MHW (1996) and 3.09 m MHW (2010).
- **Regrowth time** is Barrier3D's logistic growth, iterated from the restart height to that typical storm, using the cell's own growth rate and ceiling.

## What keeps a flattened cell flat

1. **Regrowth is slower than the storms.** A cell lowered to the berm is reset to `DuneRestart` = 7.5 cm (`barrier3d.py:1478`). Growth is logistic, G = r·D·(1 − D/Dmax), so it scales with the height already there. From 7.5 cm it takes a median of **8–9 years** to clear a typical year's storm, and that storm comes every year. Once flattened, a cell is overtopped and reset every year. This is the main cause in 1996–2010 (13 of 17 managed domains) and for 12 of 31 in 2010–2024.
2. **The ceiling is below a typical year's storm.** The per-cell ceiling is the starting crest (`DuneCeilingFromStart`, floored at 0.5 m above the berm). In the 2010 start (2009 lidar), GIS 2, 6 and 78–84 have starting crests of 1.4–2.5 m MHW. Most of their cells can never grow above 3.09 m, even fully regrown (64–100% of cells at 78, 80, 81, 83 and 84). One lidar snapshot sets the ceiling, so a breach captured in that flight stays a breach for the whole run. This is the cause in 18 domains in 2010 and 2 in 1996. Island-wide, the share of a domain's cells with a ceiling below the storm tracks its false alarms: Spearman r = 0.81 in 2010 and 0.54 in 1996.
3. **Nothing rebuilds them.** The hotspots sit where the managed run doesn't rebuild dunes:
   - GIS 1–6 is Cape Point, with no NC-12.
   - GIS 68–83 is the Tri-Village community zone, where road management is off. There, the beach/dune manager only returns 9% of the overwash to the dune, spread along it.
   
   Only the 50–55 road-managed domains get the NC-12 rebuild, and 16–21 of them used it.
4. **Progradation resets the dune rows (GIS 1 and 5).** When the shoreline advances a cell, Barrier3D moves the dune into the interior and inserts a near-zero (5 cm) dune row at the front (`migrate_dunes`). Two cells of advance replace the dune entirely. The south end advances under the solved edge rate, so its dune is wiped every year or two.

## Size

| | flattened-cell domains | of which road-managed | their share of false alarms | false alarms, Cape Point GIS 1–6 | community zones | road-managed |
|---|---|---|---|---|---|---|
| 1996–2010 | 17 | 4 | 57 of 123 | 25 | 44 | 54 |
| 2010–2024 | 31 | 10 | 82 of 132 | 16 | 37 | 79 |

The road-managed domains carry their own false alarms but are rarely flattened (≤1% of cell-years). Those come from storms overtopping intact or rebuilt dunes: a different mechanism, not covered here.

## Options (not tried; they change model inputs or code, so they are Hannah's decision)

- **Faster recovery from flat:** a higher `DuneRestart`, or a linear aeolian growth term, so a flattened cell regains decimetres a year instead of growing in proportion to its height.
- **Ceilings that don't lock in breaches:** for example, floor each cell's ceiling at its alongshore neighbours' median, or take the starting crest from more than one lidar flight.
- **Dune rebuilding in the community zones** (sand pushed back after storms in the Tri-Village), through the beach/dune manager.
- **Progradation reset:** upstream Barrier3D behaviour. It matters only at the accreting south end.
