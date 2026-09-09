# aerial-review — LIVE INPUT, not an archive

The 1996 aerial adjudication of unsurveyed holes. **58 holes reviewed by eye**,
one at a time, in `HAT_hole_aerial_picker.py`.

## Why it lives here and not in a version directory

It is keyed on `(domain, profile)`. A dune re-pick moves the CROSS-SHORE window
origin, which changes `interior_row` — it does not change which alongshore
profile a hole sits on. So these verdicts outlive any number of re-extractions,
and they are the single most expensive artifact in the audit chain to reproduce.

It used to live at `dune-topo/v1/nodata-audit/aerial_1996_conflicts/`, and the
2026-08-27 clear of `v1`/`v2` would have taken it. Hence one level up, beside
`npy-arrays/`, out of reach of version churn.

## Contents

| file | what |
|---|---|
| `aerial_review.csv` | the reviewed verdicts — **this is the irreplaceable one** |
| `aerial_review.*.bak.csv` | the picker's own timestamped backups |
| `hole_verdicts_v1.csv` | what the cleared v1 pass concluded, all three references |
| `bracketed_hole_cells_v1.csv` | the cleared v1 hole geometry, for the moved-hole diff |

The `chip_cache/` (97 MB of image chips) was **not** kept — it regenerates from
the 1996 aerial tifs by re-running `HAT_hole_aerial_chips.py`. Slow, not lost.

## How to use it on the next pass

Copy `aerial_review.csv` into the new version's
`nodata-audit/aerial_1996_conflicts/` **before** running
`HAT_test_hole_pond_or_dropout.py` — that script reads it from the resolved
version's own directory (`vdir / "aerial_1996_conflicts" / "aerial_review.csv"`)
and silently proceeds on references A and C alone if it is missing, printing
`[B] no aerial_review.csv`. Since B is what decided 52 of the 99 holes last
time, missing it changes the answer without erroring.

Then diff `utm_x`/`utm_y` per `(domain, profile)` against the new
`bracketed_hole_cells.csv` before trusting it. Each key holds only the FIRST
bracketed hole on its profile; a window that moved far enough could put a
different hole under an old verdict. Re-review anything that moved or is new —
the picker resumes on blank rows and refuses to blank a filled file.
