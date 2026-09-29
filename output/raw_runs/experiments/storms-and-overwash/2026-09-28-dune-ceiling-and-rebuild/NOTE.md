# Dune height ceiling (Dmaxel) and the NC-12 rebuild rule (2026-09-28)

**Question.** The model's dunes are held at about 3 m MHW, while the 2009 lidar puts Hatteras foredunes at about 4.9 m (`../2026-09-28-excess-overwash-diagnosis/`). Does a Hatteras ceiling, or a rebuild that never lowers dunes, fix the excess overwash?

**Design** (56 runs, all ran to the end; driver `scripts/hatteras_ms/experiments/HAT_dune_ceiling_rebuild.py`):

| factor | levels |
|---|---|
| `Dmaxel` | 3.4 m NAVD88 (current: Barrier3D's default), 5.5, 7.5, 9.0 |
| rebuild rule (managed runs only) | `current`: reset the whole dune field to 3.0 m MHW · `nolower`: same trigger, never lower a cell · `nolower43`: nolower, design height 4.3 m NAVD88 |
| storms | `drop72` (committed) and `trim24` (from the length selection) |
| windows and scenarios | 1996–2010 and 2010–2024; managed for every cell, natural for `Dmaxel` only |

- **Nothing in the main code changed.** Inside each run's process, `set_yaml` writes `Dmaxel` into the run's parameter copy, and `rebuild_dunes` is wrapped in `roadway_manager` and `beach_dune_manager`. Each run's saved model confirms the `Dmaxel` it loaded.
- **Scoring** is as in the storm-length selection: overwash against the imagery (dated by storm), the 2010 crest against the 2009 lidar, and interior RMSE. RMSE is at the matrix edge rates, which were solved on the 3.4 m ceiling and `drop72`, so it is not yet a fair shoreline comparison.

## Results (managed, current rebuild rule)

| window, storms | Dmaxel | 2010 crest − lidar | POD | POFD | PSS | timing r | RMSE |
|---|---|---|---|---|---|---|---|
| 1996–2010, drop72 | 3.4 (now) | −2.70 m | 0.17 | 0.24 | −0.07 | −0.18 | 1.05 |
| 1996–2010, drop72 | 5.5 | +0.23 m | 0.04 | 0.02 | 0.03 | 0.28 | 1.20 |
| 1996–2010, trim24 | 3.4 | −2.68 m | 0.77 | 0.41 | 0.36 | 0.42 | 1.11 |
| **1996–2010, trim24** | **5.5** | **+0.18 m** | **0.67** | **0.06** | **0.60** | **0.97** | 1.20 |
| 1996–2010, trim24 | 7.5 / 9.0 | +1.3 / +2.3 m | 0.06 / 0.04 | 0.02 | 0.04 / 0.02 | | 1.24 / 1.33 |
| 2010–2024, drop72 | 3.4 (now) | — | 0.81 | 0.44 | 0.36 | 0.66 | 2.25 |
| 2010–2024, drop72 | 5.5 | — | 0.23 | 0.07 | 0.15 | 0.79 | 1.91 |
| 2010–2024, trim24 | 3.4 | — | 0.80 | 0.75 | 0.05 | 0.33 | 2.26 |
| 2010–2024, trim24 | 5.5 | — | 0.23 | 0.09 | 0.15 | 0.78 | 1.91 |
| 2010–2024, either | 7.5 / 9.0 | — | 0.13–0.16 | 0.04–0.05 | 0.09–0.11 | 0.91–0.95 | 2.21 / 2.44 |

Full table (including natural runs and every rebuild rule): `tables/scores.csv`. Figure: `figures/dmaxel_scores.png`.

**What it shows.**

1. **The ceiling is the dominant control.** Dunes grow to it within about five years, so `Dmaxel` acts as the island's equilibrium dune height.
   - At 5.5 m NAVD88 (5.14 m MHW) the modelled 2010 crest is within 0.2 m of the lidar median, against −2.7 m now.
   - At 7.5 and 9 m the dunes grow too tall and overwash nearly vanishes.
2. **1996–2010 is fixed by 5.5 m together with keeping the long storms.** With trim24 at 5.5 m, PSS is 0.60 (from −0.07), false alarms are 6%, and the per-image timing correlation is 0.97. With `drop72` at 5.5 m it collapses to 0.03, because the storms that can overtop realistic dunes, Isabel above all, are the ones the 72 h rule deletes.
3. **2010–2024 improves in timing, false alarms and shoreline error, but now underpredicts.** At 5.5 m, false alarms fall from 44% to about 8% and timing r rises from 0.66 to 0.79. Interior RMSE improves from 2.25 to 1.91 even at the old ends. But the hit rate falls from 0.81 to 0.23, so PSS is 0.15.
4. **Why 2010–2024 underpredicts: one ceiling erases the low spots.**
   - Within a year the model grows the lowest third of domains from their lidar crest of 3.56 m to 4.33 m MHW.
   - Irene (Rhigh 3.80 m MHW) then overtops 15 of 30 low-dune domains against 25 observed, and 8 of 60 of the rest against 47 observed.
   - Real washover goes through local low points, which a single island-wide ceiling that every dune grows toward removes.
5. **The rebuild rule hardly matters once the ceiling is realistic.** Tall dunes rarely fall below the trigger, and `nolower` and `nolower43` give nearly identical results. It matters only at the current 3.4 m ceiling, where `nolower43` lifts the managed crest to 3.94 m.

**Reading.** `Dmaxel` must be set for Hatteras. About 5.5 m NAVD88, the lidar median crest, is right on average, and it must go with a storm series that keeps the long events (trim24). What limits 2010–2024 is now spatial: the model has no persistent low spots. The next experiment is a per-domain ceiling taken from each domain's lidar crest. `Dmaxel` is a per-Barrier3D-object value, so each domain can carry its own, set in-process without a code change. How fast dunes grow into gaps is the other lever.

**Status: record, decision pending (Hannah).** Nothing adopted. Shoreline skill needs the edge rates re-solved on whatever is chosen.
