# 4-comparisons/shoreline_vs_duneline/projected/1996_2024 - provenance

Written 2026-09-22 10:11 by scripts/input_prep/5-scr/projected_vs_duneline/projected_vs_duneline.py.

- Shoreline: `3-rates/coastsat/lrr/1996_2024/transect_lrr_full.csv`, lrr_m_yr x 25.7166 yr (the dune-line interval, so both cover 1997-10-12 to 2023-07-01; the end date is ASSUMED). Not the 28-yr `3-rates/coastsat/lrr_projected`.
- Dune line: `3-rates/duneline/endpoint/1996_2024/` as stored (1997 and 2023 lines).
- Beach-width change = shoreline change - dune-line change; positive = widened.
- Domain means only: the two use different transects.

## Island summary

Domain mean shoreline +3.9 m, dune line -14.8 m, beach width +18.7 m (range -28.9 to +101.8); beach narrowed in 17 of 90 domains. r(shoreline, dune line) = 0.77. y axis ±110 m.

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
git log --diff-filter=D --oneline -- scripts/input_prep/5-scr/superseded_20260921/projected_vs_duneline/projected_vs_duneline.py
git show <commit>^:scripts/input_prep/5-scr/superseded_20260921/projected_vs_duneline/projected_vs_duneline.py
```
