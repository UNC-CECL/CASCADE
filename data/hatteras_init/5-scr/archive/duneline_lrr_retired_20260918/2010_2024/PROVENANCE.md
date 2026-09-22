# duneline_lrr/2010_2024 - provenance

Written 2026-09-16 by scripts/input_prep/5-scr/duneline_lrr/duneline_lrr.py.

## Surveys in the fit

| vintage | file | survey date | |
|---|---|---|---|
| 2009 | 2009_duneline_offset_raw.csv | 2009-05-30 |  |
| 2023 | 2023_duneline_offset_raw.csv | 2023-07-01 | ASSUMED mid-year, flight date unknown |

2 island-wide surveys on 450 transects; ordinary least squares of seaward position (-ORIG_LEN) against decimal survey date per transect. With two surveys the slope is the endpoint rate and r_squared, p_value and unc_m_yr are NaN.

Buxton-only clips (1967, 2017) and the 1978 line are not in the fit; see the script header.

## Island summary

mean of transect rates +0.087 m/yr, median +0.107, 34% of transects landward.

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
