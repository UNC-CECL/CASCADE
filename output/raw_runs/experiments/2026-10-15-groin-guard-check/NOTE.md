# 2026-10-15 groin guard check

Does the tested-range guard added to `cascade.groin.CorrectedGroinCallback` on 2026-10-15 leave the Buxton
hindcast alone, and does it stop the wave settings the stress test found broken?

Four runs, edgeBE, full_management, groin on, everything else at the code defaults:

| run | result |
|---|---|
| 1996-2009, option A | ran; guard silent; shoreline matrix identical to the matrix run (max difference 0.0 m) |
| 2009-2025, option A | ran; guard silent; identical to the matrix run (0.0 m) |
| 1996-2009, wave asymmetry 0.8 | stopped in 2000: groin face 46.3 deg, past 45 (the sweep cell ran on to 59.7 deg) |
| 1996-2009, high-angle fraction 0.3 | stopped in year 3: the groin lost 10.8 m of summed shoreline (gross 111.2 m) |

Logs: `run_*.log`. Results doc: `hard-structures/groin/2-module-tests/7-stress/README.md`.
