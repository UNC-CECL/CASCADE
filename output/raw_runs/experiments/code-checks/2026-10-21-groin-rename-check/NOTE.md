# 2026-10-21 groin rename check

Does renaming the variables in `cascade/groin.py` (single letters and unitless names given full names with units,
per scripts/STYLE.md: `s` to `strength`, `q` to `trapped_m`, `_a` to `_footprint_m`, `_prev` to `_last_year` with
named fields, `_lo`/`_hi` to `_lo_pad`/`_hi_pad`, `_q_net` to `_drift_table_m3_yr`, and the helpers' locals)
change any result?

Two runs, 1996-2009, edgeBE, full management, groin on, everything else at the runner's defaults as of 23:30 on
2026-10-20 (after c8092db1): `old/` with the module from before the rename loaded in its place, `new/` with the
renamed one.

| check | result |
|---|---|
| shoreline matrix | bit-identical |
| groin_diagnostics.csv | byte-identical |
| unit tests (tests/test_groin.py) | 11 of 11 pass |
| stress-test rig, old against new module | 219 runs bit-identical; 78 stopped by the guard in both, same year and message |
| set E spread patch (spread_relaxation.py), old against new | 1,650 model years bit-identical over 45 cases; 6 stopped the same way |

The four experiment scripts that read the module's private names were updated with it: stress_test.py,
tapered_fillet.py, spread_relaxation.py, field_windows.py; and tests/test_groin.py.
