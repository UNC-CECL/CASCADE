# 2-module-tests — what the groin module does on its own

Tests of the code, not of Buxton: does each version of the module behave like a groin before it is fitted to anything? Numbered from simplest to most complex. Each folder has its own README.

```
TEST_PLAN.md        the 2026-09-11 design: tier 0 offline (built), tiers 1-2 full-model
                    rigs on fake barriers (not built yet)
1-straight-coast/   2026-10-08. Version 1 (the source/sink dipole) against version 6 (the
                    corrected drift groin) on straight coasts at -40 to +40 deg in BRIE's
                    alongshore solve alone; trimmed to these two on 2026-10-12
2-solver-audit/     2026-09-11. The dipole in BRIE's solve alone: diffusivity, runaway,
                    closure, domain count, groin fields, the chosen rig
3-real-planform/    2026-09-29. The dipole and the first blocking emulator on the real
                    Hatteras planform under option A waves, and the full-model grids
                    that led to the blocking groin
4-groin-field/      2026-10-09. Four groins 1 to 10 km apart on one straight coast: the
                    dipole, pinned, drift and corrected versions scored as a groin field
5-fine-grid-reference/  2026-10-09. A Buxton-sized field (three groins, 200 m apart) on a 25 m
                    one-line model with GENESIS-style bypassing, as the yardstick for the modules
6-long-term/        2026-10-09. The modules and the reference over 200 years: stability,
                    sand conservation, and whether the gap levels off
figures/            groin_module_logic: the dipole's arithmetic as a schematic
archive/            the retired versions (2 pinned, 3 conserving, 4 drift, 5 hybrid,
                    7 length): their straight-coast results, the all-version figures,
                    and why each was dropped. Version 2 is still the groin in the runner
```

**The versions.** Seven were tested; two are kept here. Version 1, the dipole, moves a fixed amount of
sand across the groin whatever the waves do; it was the first module and is the baseline. Version 6, the
corrected groin, holds back part of the real drift and gives back the sand BRIE's solve creates or destroys
at the groin; it is the current candidate. Tests 4-6 still score some retired versions next to version 6,
because those comparisons are the evidence for choosing it.

**What these found, in one line each.**
1. BRIE's row-scaled solve creates or destroys sand at the bend a groin makes; the dipole's results are mostly that error (+752 to −1,481 m of summed shoreline in 50 years), while the corrected groin traps the real drift and conserves sand within 27 m at every orientation.
2. The dipole's fillet is bought by M / diffusion number, the solve is not volume-conserving, and five domains are too few.
3. The old M = 60 is not runnable under option A; a blocking groin, bounded by the transport that arrives, is.
4. With four groins on one coast, the corrected version behaves most like a groin field (105 of 112 checks against 60 to 75 for the others), but close groins share their sand only weakly, since no version has a groin length.
5. Against a real-sized field on a 25 m grid, the corrected version is closest (gap error 71 m against 103 to 241 m), but the field's gap is set by groin length (about twice it), which no module has.
6. Over 200 years nothing runs away and the corrected version holds its gap, but on a grid-aligned coast it leaks sand slowly (232 m by year 200) unless its correction zone is widened from 5 to 10 cells.

**Next step up (not built).** The same straight-coast test inside full CASCADE on a fake barrier, storms off then on (`TEST_PLAN.md` tiers 1-2), once a module version passes step 1.

Two scripts in `3-real-planform/` (`diagnose_no_groin_relaxation.py`, `trajectory_check.py`) read the 1996-2010 matrix runs, which moved to `D:\archive\CASCADE_offload\` on 2026-10-07; they run only with that drive attached.
