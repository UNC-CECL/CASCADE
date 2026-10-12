# archive — groin versions no longer carried forward

Moved here on 2026-10-12, when the module tests were trimmed to the two versions that matter: the first
(version 1, the source/sink dipole) and the current candidate (version 6, the corrected drift groin).
Kept in the repo, not deleted, because the hindcast studies under `../../3-hindcast/` cite these results,
and version 2 is still the groin in the runner.

**The code was not removed.** `../1-straight-coast/straight_coast_test.py` still defines and runs all seven
versions (tests 4-6 and the hindcast studies import it). Versions 2-5 and 7 carry `archived=True`, so
rerunning the test writes their tables and figures here, along with the all-version comparison figures.
The tables every version shares (`straight_coast_summary_year50.csv`, `pelnard_considere_check.csv`,
`conservative_solve_check.csv`) stay in `../1-straight-coast/comparison/tables/`.

```
1-straight-coast/
    2-trapping-pinned/          tables/ and figures/ for each retired version
    3-trapping-conserving/
    4-trapping-drift/
    5-trapping-hybrid/
    7-corrected-length/
    comparison/figures/         all seven versions together: rules, mechanics, profiles at year 50,
                                orientation summary, and the year-by-year GIF (git-ignored)
    straight_coast_module_rules_figure.py        the rule of each version, sketched
    straight_coast_trapping_mechanics_figure.py  how the trapping versions work inside BRIE
    straight_coast_gif.py                        all versions year by year
```

The three scripts run from here: `python straight_coast_module_rules_figure.py` and so on, after
`../1-straight-coast/straight_coast_test.py`.

## The retired versions, and why

| # | version | what it does each year | why it was dropped |
|---|---|---|---|
| 2 | **pinned** (`BlockingGroinCallback`, b 0.6, f 0.6) | cancels 60% of BRIE's own diffusive exchange across the groin face | Conserves sand, but blocks the wrong thing: that exchange follows the coast's tilt, not the drift, so on a straight coast it does nothing while 0.91 M m³/yr passes, and on a coast turned with the drift it builds the "fillet" on the downdrift side. **Still the groin in the runner** until the corrected groin replaces it. |
| 3 | **conserving** | as 2, with one face coefficient on both sides | Balances its own books but leaves BRIE's solve error in place (+186 m at Buxton, 1996-2009). |
| 4 | **drift** | holds back a fraction of the wave climate's net drift arriving at the groin | Traps the right quantity, but BRIE's row-scaled solve creates or destroys sand around the bend it makes (+107 m at Buxton). Version 6 is this plus the fix. |
| 5 | **hybrid** | 4, plus 2's row correction at the groin face | Does not close the balance: the solve's error is spread along the whole fillet, not only at the face. |
| 7 | **corrected + length** | 6, with blocking fading as the fillet reaches the groin tip (L 110 m) | Caps the fillet too low: it reads the whole 500 m cell's build-out, which BRIE's diffusion keeps far short of the tip. |

The full straight-coast findings that led here (findings 1-9, and what they meant for the module) are in
`1-straight-coast/README.md`, the straight-coast README as it stood before the trim.

**Tests 4-6 still compare these versions.** `../4-groin-field/`, `../5-fine-grid-reference/` and
`../6-long-term/` score the pinned, drift and length versions alongside the corrected one; their results
were left in place because they are the evidence for choosing version 6.
