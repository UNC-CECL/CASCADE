# 2026-10-08-conserving-groin

**Question.** Does a blocking groin that conserves sand across the GIS 5|6 face still fit the 1996–2009 gap?

**Answer.** Stopped after the smoke cell, because the premise was wrong. The pinned groin's imbalance offsets BRIE's own non-conservative solve, so across the whole island the pinned groin is the one that conserves sand (−22 m against no groin, versus +186 m for this variant). At b 0.6 / f 0.6 the conserving groin overshoots the 2004 gap by about 100 m.

**Write-up.** `hard-structures/groin/groin-module-test/1-dem-to-dem/2026-10-08-conserving-groin/README.md`, and the solver audit it points to, `hard-structures/groin/groin-module-test/0-solver-audit/2026-10-08-straight-coast/README.md`.

**Folder names.** `b0.60_f0.6/1996_2009/`: the one smoke run. Its `conserving.txt` records the patch, because the runner's report does not show it.

**Status.** Record. The grid was never run; relaunch with `python conserving_groin.py grid` from the write-up folder if wanted.
