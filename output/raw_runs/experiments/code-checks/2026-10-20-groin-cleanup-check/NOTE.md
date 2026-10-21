# 2026-10-20 groin cleanup check

Does the 2026-10-20 cleanup of `cascade/groin.py` (renamed to `GroinCallback`; `install_year` taken from
the strength history; one `record` instead of 11 per-year lists; shared angle and table-index helpers)
change any result?

Four runs, edgeBE, full_management, groin on, everything else at the code defaults (storms `v3_240`), on
the same runner: `old/` with the pre-cleanup module loaded in its place, `new/` with the cleaned one.

| window | shoreline matrix | groin_diagnostics.csv |
|---|---|---|
| 1996-2009 | bit-identical | byte-identical |
| 2009-2025 | bit-identical | byte-identical |

A first pair of `new/` runs used the `v3_split12_trim24` storms (the default changed at 13:33 during the
check) and was rerun; the logs `run_*.log` are from the final runs. The same old/new comparison on the
stress-test rig (219 runs, real BRIE tables) also matched bit for bit.
