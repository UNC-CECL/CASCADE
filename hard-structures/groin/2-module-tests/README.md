# 2-module-tests — what the groin module does on its own

Tests of the code, not of Buxton: does each version of the module behave like a groin before it is fitted to anything? Numbered from simplest to most complex. Each folder has its own README.

```
TEST_PLAN.md        the 2026-09-11 design: tier 0 offline (built), tiers 1-2 full-model
                    rigs on fake barriers (not built yet)
1-straight-coast/   2026-10-08. Four versions (source/sink dipole, trapping pinned,
                    trapping conserving, trapping drift), one subfolder each, on straight
                    coasts at -40 to +40 deg in BRIE's alongshore solve alone
2-solver-audit/     2026-09-11. The dipole in BRIE's solve alone: diffusivity, runaway,
                    closure, domain count, groin fields, the chosen rig
3-real-planform/    2026-09-29. The dipole and the first blocking emulator on the real
                    Hatteras planform under option A waves, and the full-model grids
                    that led to the blocking groin
figures/            groin_module_logic: the dipole's arithmetic as a schematic
```

**What these found, in one line each.**
1. The pinned groin conserves sand but blocks only the tilt-driven part of the transport, so at 0° it does nothing while 0.91 M m³/yr of net drift passes; a drift-blocking groin traps correctly, but BRIE's row-scaled solve creates or destroys sand around it.
2. The dipole's fillet is bought by M / diffusion number, the solve is not volume-conserving, and five domains are too few.
3. The old M = 60 is not runnable under option A; a blocking groin, bounded by the transport that arrives, is.

**Next step up (not built).** The same straight-coast test inside full CASCADE on a fake barrier, storms off then on (`TEST_PLAN.md` tiers 1-2), once a module version passes step 1.

Two scripts in `3-real-planform/` (`diagnose_no_groin_relaxation.py`, `trajectory_check.py`) read the 1996-2010 matrix runs, which moved to `D:\CASCADE_offload\` on 2026-10-07; they run only with that drive attached.
