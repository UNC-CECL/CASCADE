# 5-long-term — how the groin modules hold up over 200 years

**Question.** Over a long run, does each version of the groin module hold a steady gap, keep conserving sand and stay numerically stable, and does it keep tracking a real-sized groin field?

**Setup.** One groin on the straight-coast rig: 61 domains, option A waves, BRIE's alongshore solve alone (no storms, no Barrier3D, no sea-level rise). The five modules run 200 years at −30° to +30°: dipole, pinned, drift, corrected, and corrected with length (module 7). The 25 m reference field with three 110 m groins (`../4-fine-grid-reference/`) runs alongside. 200 years is the longest run in `../TEST_PLAN.md`. The three scripts imported are used unchanged. `../3-groin-field/groin_field_test.py`'s `run_field()` gained an optional `record_years` argument on 2026-10-09 to return yearly shorelines; its default output is byte-identical.

```
long_term_test.py      the runs and the tables (about 12 minutes with the reference; --skip-reference for seconds)
long_term_figures.py   both figures, from the tables
tables/                long_term_by_year.csv (gap, both sides, net change, steepest angle, every year),
                       long_term_summary.csv (gap at 50, 100, 200 yr, growth over the last 50, net change, steepest angle)
figures/               long_term_1_gap_over_time, long_term_2_end_state; PDFs and CAPTIONS.md in supporting/
logs/run_log.txt       the 2026-10-09 run's printed summary
```

    python long_term_test.py
    python long_term_figures.py

## Findings

1. **No module runs away.** Every run stays finite for 200 years. The steepest shoreline angle at any face in any year is 56°, against the 90° where BRIE's diffusivity table ends. That is the drift module at −30°; the wave-driven groin's steepest is 46°, also at −30°.
2. **The reference field levels off and stays there.** Its gap is 221 to 229 m at year 200 at every orientation, about twice the 110 m groin length. At +30°, where the drift is weakest, it takes about 80 years to fill.
3. **The wave-driven groin holds steady and stays close where it was close.** Its gap changes by at most 1.3 m over years 150-200 at every orientation but +10°. At −30° and −20° it ends at 226 and 219 m, against the reference's 222 and 223 m. Where it was off at 50 years it stays off: 257 m against 229 m at −10°, 335 against 223 m at 0°, 46 against 221 m at +20° and 3 against 222 m at +30°. At +10° it is still creeping (288, 297 and 318 m at years 50, 100 and 200) against a reference of 222 m.
4. **The wave-driven groin leaks sand on a grid-aligned coast.** At 0° it loses 12 m of summed shoreline by year 50, then about 1.5 m a year, 232 m by year 200. Everywhere else it stays within 19 m. The cause is the correction's 5-cell zone. At 0°, BRIE's diffusivity sits at its clip at zero, so any turn of the shoreline switches it on or off between neighbouring cells, and BRIE's row-scaled solve leaks along the whole eroding tail. That tail reaches 8.5 cells downdrift by year 200. Rerun with a 10-cell zone, the leak at 0° is 1 m at year 200, and the gaps change by under 4 m (−20°: 219 m; 0°: 338 m; +10°: 320 m). This is a recommendation, not a change: the module's zone is still 5 cells.
5. **Module 7 is steady and conserves sand** (within 37 m; the same 0° leak, smaller, because its fillet is smaller). But its gap stays at 126 to 185 m at −30° to +10°, as at 50 years.
6. **The others drift further.** The pinned version is still growing at −30° to −10° (+5 to +9 m over the last 50 years) and on the wrong side from +10° on. The fixed-shift groin and drift versions keep creating or destroying sand year after year: by year 200 the fixed-shift groin is at +4,117 m at −20° and −5,082 m at +10°, the drift version at −18,758 m at +20°. The drift version's gap also zigzags year to year at −30° and −20°.

**What this means for the module.** Over 200 years the wave-driven groin is stable and holds its gap. Where it matches the real-sized field at 50 years, it still matches at 200. Its one long-run defect is the slow sand leak on a grid-aligned coast, fixed by a 10-cell correction zone. Its orientation dependence (finding 3) does not change with time.
