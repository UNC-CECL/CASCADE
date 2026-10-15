# 3-module-tests — what the groin module does on its own

**Code these studies ran against.** The fixed-shift (`GroinCallback`, the fixed-shift groin) and blocking (`BlockingGroinCallback`) groins that the scripts here use or patch were removed from `cascade/groin.py` on 2026-10-14, when the wave-driven groin became the runner's only groin. The scripts ran against the code at commit `4d937b79`; check that commit out to rerun one.

Tests of the code, not of Buxton: does each version of the module behave like a groin before it is fitted to anything? Numbered from simplest to most complex (since 2026-10-15): tests 1 to 7 run BRIE's alongshore step alone on idealized coasts, and test 8 runs the full model on the real Hatteras planform. The figures that explain how the fixed-shift groin and the wave-driven groin work are in `../2-module-explained/`.

```
TEST_PLAN.md            the 2026-09-11 design: tier 0 offline (built), tiers 1-2 full-model
                        rigs on fake barriers (not built yet)
1-solver-audit/         2026-09-11. The fixed-shift groin in BRIE's solve alone: diffusivity, runaway,
                        closure, domain count, groin fields, the chosen rig
2-straight-coast/       2026-10-08. The fixed-shift groin (the source/sink dipole) against the
                        wave-driven groin (the corrected drift groin) on straight coasts at -40 to +40 deg in BRIE's
                        alongshore solve alone; trimmed to these two on 2026-10-12. Its
                        1-run/straight_coast_test.py is the rig tests 3 to 7 import
3-groin-field/          2026-10-09. Four groins 1 to 10 km apart on one straight coast: the
                        fixed-shift, pinned, drift and wave-driven versions scored as a groin field
4-fine-grid-reference/  2026-10-09. A Buxton-sized field (three groins, 200 m apart) on a 25 m
                        one-line model with GENESIS-style bypassing, as the yardstick for the modules
5-long-term/            2026-10-09. The modules and the reference over 200 years: stability,
                        sand conservation, and whether the gap levels off
6-stress/               2026-10-15. The runner's own groin (CorrectedGroinCallback, unchanged) pushed
                        past where Buxton uses it: wave climates, strength, a sharp start, switch-off,
                        cycling, year-to-year waves, Barrier3D-like input, placement, 500 years
7-mirror/               2026-10-15. Drift from either direction: the groin's response mirrors exactly;
                        whole runs mirror only to within BRIE's own directional bias
8-real-planform/        2026-09-29. The fixed-shift groin and the first blocking emulator on the real
                        Hatteras planform under option A waves, the full-model grids that led to
                        the blocking groin, and the wave-driven groin on the real planform
archive/                the retired versions (2 pinned, 3 conserving, 4 drift, 5 hybrid,
                        7 length): their straight-coast results, the all-version figures,
                        and why each was dropped. Version 2 is still the groin in the runner
```

**Inside every test.** The same layout throughout (since 2026-10-15):

```
README.md               the question, how to run it, what it found
<script>.py             the test and its figure scripts, at the top of the folder
tables/                 what the scripts write: CSVs, and .txt copies of printed results
figures/                PNGs (git-ignored; rerun to make them) and GIFs
    supporting/         each figure's PDF and CAPTIONS.md
logs/                   run logs (printed output of a run, grid launch logs)
```

A test made of several steps (`2-straight-coast/`, `8-real-planform/`) has one numbered subfolder per
step, each with that same layout, and one README at the test's top level. A folder only has the parts it
uses: `7-mirror/` writes no figures.

**The versions.** Seven were tested; two are kept here. The **fixed-shift groin** (version 1, once "the
dipole") moves a fixed amount of sand across the groin whatever the waves do; it was the first module and
is the baseline. The **wave-driven groin** (version 6, once "the wave-driven groin") holds back part of the real drift and gives back the sand BRIE's solve creates or destroys
at the groin; it is the current candidate. Tests 3-5 still score some retired versions next to the wave-driven groin,
because those comparisons are the evidence for choosing it.

**What these found, in one line each.**
1. The fixed-shift groin's fillet is bought by M / diffusion number, the solve is not volume-conserving, and five domains are too few.
2. BRIE's row-scaled solve creates or destroys sand at the bend a groin makes; the fixed-shift groin's results are mostly that error (+752 to −1,481 m of summed shoreline in 50 years), while the wave-driven groin traps the real drift and conserves sand within 27 m at every orientation.
3. With four groins on one coast, the wave-driven groin behaves most like a groin field (105 of 112 checks against 60 to 75 for the others), but close groins share their sand only weakly, since no version has a groin length.
4. Against a real-sized field on a 25 m grid, the wave-driven groin is closest (gap error 71 m against 103 to 241 m), but the field's gap is set by groin length (about twice it), which no module has.
5. Over 200 years nothing runs away and the wave-driven groin holds its gap, but on a grid-aligned coast it leaks sand slowly (232 m by year 200) unless its correction zone is widened from 5 to 10 cells.
6. At option A the runner's groin holds, but it breaks where BRIE's diffusion at the groin is too strong (its explicit leak correction runs away: high-angle 0.3 and below, Hs 2.5+ on coasts turned with the drift) or clipped to zero (the gap grows without limit: asymmetry 0.7+), and under year-to-year waves. Since 2026-10-15 the module stops a run that leaves its tested range (`GroinOutsideTestedRange`: one wave climate, the groin's own sand balance within 10 m or 5%, the groin face at most 45°, the correction zone inside the domains); both Buxton windows run unaffected.
7. The groin treats drift from either direction identically (one call mirrors to 0.0 m); whole runs do not mirror exactly because BRIE's alongshore step and angle tables are themselves slightly direction-dependent, and the groin adds at most 10 m to that at Buxton-like settings.
8. On the real planform the old M = 60 is not runnable under option A; a blocking groin, bounded by the transport that arrives, is.

**Next step up (not built).** The same straight-coast test inside full CASCADE on a fake barrier, storms off then on (`TEST_PLAN.md` tiers 1-2), once a module version passes step 1.

Two scripts in `8-real-planform/` (`3-no-groin-relaxation/diagnose_no_groin_relaxation.py`,
`5-failure-schedule/trajectory_check.py`) read the 1996-2010 matrix runs, which moved to `D:\archive\CASCADE_offload\` on 2026-10-07; they run only with that drive attached.
