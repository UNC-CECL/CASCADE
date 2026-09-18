# Retired — the calibration as it stood before the frozen-set correction

The converged source/sink field of 2026-09-12, and the config that carried it,
kept because the field below it was re-derived on 2026-09-14 and every run built
before that date was built against THESE numbers.

## What was wrong with it

`FROZEN_ZONE_DOMAINS` held seven domains that the zone rules do not put there.
Each one is warranted by the significance and width rules — in the **other**
period:

| domain | was in the set for | the rules warrant it in |
|---|---|---|
| D8  | 2004-2024 | 1984-2004 |
| D9  | 1984-2004 | 2004-2024 |
| D22 | 1984-2004 | 2004-2024 |
| D44 | 1984-2004 | 2004-2024 |
| D48 | 2004-2024 | 1984-2004 |
| D57 | 2004-2024 | 1984-2004 |
| D62 | 1984-2004 | 2004-2024 |

Seven for seven, and each is already a legitimate member of the other period's
tuple — D22 sits inside period 2's D8-D22 run, D48 inside period 1's D48-D57.
They were written into both. That is a transposition when the two tuples were
assembled, not drift and not a hand-edit chasing a residual.

It surfaced as five corrections **one domain wide** (D22, D44, D62 in period 1;
D48, D57 in period 2), which `MIN_ZONE_WIDTH = 3` exists to make impossible. D8
and D9 are the same error but land beside legitimate members, so they never
showed as singletons.

## How it was established

Pass 0 was re-derived against the current edgeBE base runs with
`HAT_BE_OUTPUT_DIR` redirected to a scratch directory, so nothing in the
production calibration was touched. The width rule was then re-applied to the
resulting residuals directly, because `compute_be_rates` overwrites the
significance verdict at D5-D7 (the groin reservation) after the rule has run,
which makes D8 look isolated in the metrics CSV when it is not.

The re-derived set contains no isolated members in either period and is a strict
subset of the committed one: nothing the rules produce was ever missing, only
these seven added.

## What it was worth

Six domains carried a value they should not have had:

```
1984-2004   D9 +0.00   D22 +0.30   D44 +0.30   D62 +0.30      sum +0.90
2004-2024   D48 -0.40  D57 -0.30                              sum -0.70
2004-2024   D8  -1.80                             <- the one that is not small
```

against field totals of -23.8 (period 1) and +64.5 (period 2). D9 had resolved
to exactly 0.0, so its membership was spurious but its value was not.

Nothing here is maintained or expected to run. Rule 4 of `ORGANIZATION.md`:
retirement is a dated folder with a note, and deletes nothing.
