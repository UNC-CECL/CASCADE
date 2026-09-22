# smoothing_test — what does the smoother do to this comparison?

**A test folder, not a product.** Built 2026-09-22 (Hannah, by interview) to
see how much of the shoreline-vs-dune-line disagreement survives smoothing at
the scale the model resolves.

Two sheets, both in the same form as the unsmoothed
`*_halves_overlay.png` they are paired with — 1996–2010 above 2010–2024, one
fixed ±100 m axis:

| sheet | shoreline side | dune side |
|---|---|---|
| `loess7_projected_vs_duneline_...` | the **same** 1996–2024 LRR × 14 yr in both panels | that half's own measured change |
| `loess7_total_change_vs_duneline_...` | **each panel's own** LRR × its own 14 yr | that half's own measured change |

The dune line always follows the sub-period. Only the shoreline reading
differs between the two sheets, which is the whole point of having both.

## The smoothing

- **7 domains = 3.5 km**, the width Hannah asked for.
- **Both curves get the same pass at the same width, at transect resolution**,
  then average to domains. Smoothing one side and not the other would make
  the gap between them an artefact of the treatment rather than a beach-width
  change. CoastSat has ~10 transects per domain and the dune line exactly 5;
  both are spread evenly inside their domain by `rates_figures._along`, so the
  two are handled identically.
- **GIS 1–10 keep their raw domain means** — the Oregon Inlet boundary
  treatment the scoring target uses (`coastsat_loess.skip_southern_domains`).
  Hannah chose to keep it so the figure shows the target the way the model
  actually sees it. **Those ten domains are therefore identical to the
  unsmoothed sheet by construction**, and any difference there is not the
  smoother. It is visible on the figure: the two curves stay jagged over GIS
  1–10 and only settle north of it.
- **The raw domain means stay on the figure** as faint dots behind each curve,
  so what the smoother removed is visible without leaving the sheet.

## What it shows

| shoreline reading | window | r raw | r smoothed | beach width raw | smoothed |
|---|---|---|---|---|---|
| projected | 1996–2010 | 0.37 | 0.41 | +18.3 m | +17.8 m |
| projected | 2010–2024 | 0.85 | 0.87 | +0.8 m | +0.6 m |
| total | 1996–2010 | 0.59 | 0.58 | +11.6 m | +10.4 m |
| total | 2010–2024 | 0.66 | 0.69 | +14.7 m | +14.7 m |

**Read the beach width, not r.** A symmetric smoother strips high-frequency
variance that is uncorrelated between the two series, so r drifts up almost
regardless of whether the two features genuinely agree better at 3.5 km — the
same circularity that `target_comparison/smoothing_scale/` found and that the
3-rates smoothed panels warn about. The beach width barely moves (at most
1.2 m, and 0.0 m in one case), which is the honest summary: **smoothing at
3.5 km does not change what this comparison says.**

The substantive result is the one already visible unsmoothed — the long-term
rate tracks the dune line closely in 2010–2024 (r 0.85) and poorly in
1996–2010 (r 0.37), while each half's own rate sits in between at both.

## Files

```
loess7_projected_vs_duneline_1996_2010_2024_halves_overlay.png
loess7_total_change_vs_duneline_1996_2010_2024_halves_overlay.png
domain_smoothed.csv    per domain per window per product: both sides raw and
                       smoothed, and the beach-width gap of each
PROVENANCE.md          the table above, with the reading rule
supporting/            PDFs, CAPTIONS.md, island_summary.csv
```

Rebuild: `python scripts/input_prep/5-scr/smoothing_test/smoothing_test.py`
(`--window N` for a different width). The unsmoothed pair it is read against
is in `../coastsat_projected_vs_duneline_endpoint/all_windows_stacked/` and
`../coastsat_total_change_vs_duneline_endpoint/all_windows_stacked/`.

---

**2026-09-22 — the producer named above is no longer on disk.** The 5-scr
scripts tree was reorganised to mirror this data tree that day, and its
retired scripts were deleted rather than parked in a dated folder (a departure
from rule 4 of `ORGANIZATION.md`, taken deliberately; see
`scripts/input_prep/5-scr/README.md`). The product here is unaffected — only
the path that made it has gone.

**This one is not recoverable.** `smoothing_test/smoothing_test.py` was never
committed — `git log -- scripts/input_prep/5-scr/smoothing_test/` returns
nothing on any branch — so it was already absent before the 2026-09-22
reorganisation, and that reorganisation is not what removed it. The figures it
produced are the only record of it. What replaced it,
`4-comparisons/shoreline_vs_duneline/smoothed_loess7.py`, is committed.
