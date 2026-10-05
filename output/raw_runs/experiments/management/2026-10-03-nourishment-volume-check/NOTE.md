# Nourishment volume check (2026-10-03)

Question: are the beach nourishment volumes applied as intended for every fill, in both model windows?

Script: `scripts/hatteras_ms/experiments/HAT_nourishment_volume_check.py` (reads runs, runs nothing).
Runs: 1996-2015 edgeBE full_management under `runs/full_management/` (made for this check, unchanged code);
2010-2026 is the with2017 member of `../2026-10-03-buxton-2017-fill/`, which matches the matrix run exactly.

Answer: yes. Every scheduled fill fired in its year, in its domains, at its volume. Nothing extra fired.

- Sheet vs config: Buxton 2017, Buxton 2022 and Avon 2022 match the sheet volumes. Rodanthe 2014 is
  1.60 M cy in the config vs 1.62 M cy in the sheet (-1.2%), over 3.0 km vs the sheet's 3.43 km.
  Lengths differ by up to 0.3 km because footprints are whole 500 m domains; the total volume is kept.
- yd3 -> m3/m: the config values match a hand calculation (407.76, 397.57, 183.49, 191.14).
- Applied: each BeachDuneManager's own record equals the schedule in all 40 domain-fills (6 + 34).
  Total placed equals the project total.
- Shoreline step: equals 2V/(2 h_b + D_sf) to <0.001 m in 37 of 40. The three exceptions (GIS 8 in 2017
  and 2022, GIS 26 in 2022) are community domains where the overwash filter put sand back on the
  shoreface in the same step. The extra step is a constant 82% of the overwash removed, so it is
  filter sand, not a fill error.
- Persistence: with-minus-without 2017 shoreline is 36-40 m at the fill and ~26 m at the end of the run.

Tables: `tables/source_vs_config.csv`, `tables/fills_by_domain.csv`.

**Follow-up, same day:** Hannah set Rodanthe 2014 to the sheet's 1.62 M cy in HATTERAS_NOURISHMENT_PROJECTS
(footprint unchanged, GIS 84-89). Re-ran 1996-2015: 412.86 m3/m applied in all six domains, step 39.2-43.2 m,
total 1,238,579 m3, all checks pass. The 2010-2026 rows in `tables/fills_by_domain.csv` still come from the
1.60 M cy run.
