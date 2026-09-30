# Dune ceiling from each domain's own dunes (2026-09-28)

**Question.** One island-wide `Dmaxel` of 5.5 m NAVD88 matches the lidar on average and fixes 1996–2010, but it erases the real low spots: in 2010–2024 the hit rate drops to 0.23. Does a ceiling taken from each domain's own dunes keep those low spots?

**Ceilings** (24 runs, all ran to the end; driver `scripts/hatteras_ms/experiments/HAT_dune_ceiling_per_domain.py`). Each is built from the run's own starting dunes, held at least 0.5 m above the berm:

| ceiling | rule |
|---|---|
| `dom_median` | each domain's `Dmaxel` = its median starting crest |
| `dom_p25` | each domain's `Dmaxel` = its 25th-percentile starting crest |
| `cell` | each dune cell's ceiling = its own starting crest. `DuneGrowth` is wrapped in-process to take an array; the scalar it returns is the domain median. |

- **Other settings.** Current rebuild rule; trim24 and drop72 storms; both windows; managed and natural.
- **Controls.** `uniform3p4` (the current model) and `uniform5p5` (`../2026-09-28-dune-ceiling-and-rebuild/`).
- **Nothing in the main code changed.** `build_cascade` is wrapped in each run's process to set the ceilings after construction.

**A change made during these runs, and how it is handled.** A parallel session switched the runner's CoastSat target from LOWESS-10 to LOWESS-7 and re-solved the end rates on 2026-09-28. The matrix was archived to `archive/2026-09-28-loess10-ends/` and is being re-run.

| runs | end rates (GIS 1 / 90) |
|---|---|
| these runs | 1996 4.8394 / 18.2545 · 2010 18.8657 / 24.2358 |
| the controls | 1996 4.8394 / 17.545 · 2010 18.8 / 24.535 |
| 3 runs of the uniform experiment (1996) | 1996 4.8394 / 18.2545 |

- **Overwash scores are unaffected.** They use the imagery, not CoastSat, and the ends only change the two end domains.
- **RMSE and bias are recomputed for EVERY run against one target**, LOWESS-7, with `run_registry.skill_vs_target` over GIS 2–89. This reproduces the runner's own values to 1e-4. The controls are read from the archive.

## Results (managed; overwash threshold 0; RMSE vs LOWESS-7)

| storms | ceiling | PSS 96 | timing r 96 | space r 96 | crest vs lidar r (2010) | RMSE 96 | PSS 10 | POD 10 | POFD 10 | RMSE 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| trim24 | uniform 3.4 (now) | 0.36 | 0.42 | 0.05 | 0.10 | 1.24 | 0.05 | 0.80 | 0.75 | 2.32 |
| trim24 | uniform 5.5 | **0.60** | **0.97** | 0.39 | 0.21 | 1.20 | **0.15** | 0.23 | 0.09 | **1.91** |
| trim24 | dom_median | 0.44 | 0.86 | 0.43 | 0.53 | 1.16 | 0.11 | 0.34 | 0.23 | 2.04 |
| trim24 | dom_p25 | 0.52 | 0.86 | 0.43 | 0.58 | 1.14 | 0.14 | 0.43 | 0.30 | 2.05 |
| trim24 | **cell** | 0.57 | 0.85 | **0.63** | 0.52 | 1.17 | 0.14 | 0.49 | 0.34 | 2.09 |
| drop72 | uniform 3.4 (now) | −0.07 | −0.18 | 0.12 | 0.09 | 1.19 | **0.36** | 0.81 | 0.44 | 2.31 |
| drop72 | any realistic ceiling | 0.03–0.06 | ≤ 0.28 | 0.32–0.42 | 0.21–0.57 | 1.16–1.20 | 0.12–0.17 | 0.23–0.37 | 0.07–0.25 | 1.91–2.03 |

**Irene (the 2011-08-27 image), managed trim24:**

| ceiling | low-dune third (lidar crest < 4.33 m) | the other 60 domains |
|---|---|---|
| observed | 25 / 30 | 47 / 60 |
| uniform 3.4 | 29 | 49 |
| uniform 5.5 | 15 | 8 |
| dom_median | 22 | 6 |
| dom_p25 | 27 | 8 |
| cell | 26 | 10 |

Full table: `tables/scores.csv`.

**What it shows.**

1. **Per-domain and per-cell ceilings keep the real dune pattern.** The modelled 2010 crests correlate with the lidar at r 0.52–0.58, against 0.21 for the uniform ceiling. They also put the overwash in the right places: the 1996–2010 space r is 0.63 with the per-cell ceiling, against 0.39.
2. **They fix the low spots.** Irene now overwashes 26–27 of the 30 low-dune domains, against 25 observed; the uniform 5.5 m ceiling gave 15.
3. **They do not fix the rest of Irene.** In reality Irene overwashed 47 of the 60 domains whose dunes are 4.3 m or higher. No ceiling that matches the lidar lets a Rhigh of 3.80 m MHW over those dunes; only the current, too-low 3 m dunes did. So what remains in 2010–2024 is storm height, not dunes. The Duck gauge is about 80 km north of the reach; Irene's surge on Hatteras, run-up on a 0.06 slope, and inlet breaching at Pea Island and Rodanthe are not in the storm file. This is test 3 of the diagnosis, which has not yet been run.
4. **The storms still need Isabel.** With `drop72`, every realistic ceiling scores near zero in 1996–2010.
5. **Trade-off between the uniform and per-cell ceilings:**
   - The uniform 5.5 m ceiling has the best overall skill (PSS 0.60 and 0.15) and the best 2010–2024 RMSE (1.91).
   - The per-cell ceiling is nearly as skilful in 1996–2010 (0.57) and has the best spatial pattern (space r 0.63, crest r 0.52). It catches twice the observed washover in 2010–2024 (POD 0.49 against 0.23), at the cost of more false alarms and RMSE 2.09.

**Reading.** Both are large improvements on the current model.

- **For the physical representation of the dunes: the per-cell ceiling.** It is the observed dune line itself, and it keeps the low spots through which washover actually occurs.
- **For bulk skill with the fewest parameters: the uniform 5.5 m ceiling.**

Either way the storms should keep the long events (trim24). The remaining 2010–2024 misses point at storm heights (diagnosis test 3), not at the dunes.

**Status: record, decision pending (Hannah).** Nothing adopted.
