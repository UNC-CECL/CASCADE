# Why this is archived

`duneline_lrr/` was a per-transect OLS slope through every island-wide dune
line inside a window (2026-09-16). It was retired on 2026-09-18 (Hannah:
"these should not be lrr, they would just be endpoint, we are tracking net
change"). The dune line is now scored on the net change between the two lines
that bound a window: `5-scr/3-rates/duneline_endpoint/`.

The numbers here are also stale for a second reason: they were built from the
1997, 2009 and 2023 lines as they were BEFORE the 2026-09-18 re-digitization.

Nothing reads this folder, and `hat_observed_rates` no longer resolves it. The
producer is `scripts/input_prep/5-scr/superseded_20260918/duneline_lrr.py`.
The figures drawn from it (`output/comparisons/rate_windows/**/duneline/lrr/`)
were removed at the same time.

---

**2026-09-22 — the producer named above is no longer on disk.** The 5-scr
scripts tree was reorganised to mirror this data tree that day, and its
retired scripts were deleted rather than parked in a dated folder (a departure
from rule 4 of `ORGANIZATION.md`, taken deliberately; see
`scripts/input_prep/5-scr/README.md`). The product here is unaffected — only
the path that made it has gone.

To read the script again — note the path below is the one it was last
COMMITTED under, which is not the path named above:

```
git log --diff-filter=D --oneline -- scripts/input_prep/5-scr/superseded_20260918/duneline_lrr.py
git show <commit>^:scripts/input_prep/5-scr/superseded_20260918/duneline_lrr.py
```
