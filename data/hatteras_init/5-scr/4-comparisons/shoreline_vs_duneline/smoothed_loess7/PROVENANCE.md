# smoothing_test — provenance

Written 2026-09-22 09:53 by `scripts/input_prep/5-scr/smoothing_test/smoothing_test.py`.

Both curves LOESS-smoothed at **7 domains (3.5 km)**, at transect resolution, with GIS 1–10 kept at their raw domain means (the scoring target's Oregon Inlet treatment, Hannah's choice 2026-09-22). The dune line always follows the sub-period.

| shoreline reading | window | r raw | r smoothed | beach width raw (m) | beach width smoothed (m) |
|---|---|---|---|---|---|
| projected | 1996–2010 | 0.37 | 0.41 | +18.3 | +17.8 |
| projected | 2010–2024 | 0.85 | 0.87 | +0.8 | +0.6 |
| total | 1996–2010 | 0.59 | 0.58 | +11.6 | +10.4 |
| total | 2010–2024 | 0.66 | 0.69 | +14.7 | +14.7 |

**r rises with smoothing on both sides.** That is what a symmetric smoother does — it strips high-frequency variance that is uncorrelated between the two series — so it is NOT evidence that the two features agree better at 3.5 km. Read the beach width, which is close to smoothing-invariant, and compare the SHAPE against the unsmoothed sheets in `coastsat_{projected,total_change}_vs_duneline_endpoint/all_windows_stacked/`.

GIS 1–10 are unsmoothed by construction, so no difference there is the smoother's.
