# 2026-10-21 runner spread switch check

The runner was switched to the spread groin (0, 2/3, 1/3) at strength 0.10 in both windows (`groin.spread_weights`
and `groin.strength_1996` in `scripts/hatteras_ms/hat_run.yaml`, the same code defaults in `HAT_hindcast_config.py`).
Does the runner, reading the yaml as a matrix run does, build and apply the spread groin?

One run, 1996-2009, edgeBE, full management, groin on, `hat_run.yaml` read (only the window, groin and output
switches set by environment). The runner carried another session's uncommitted changes (30 buffer domains a side
instead of 15, end rates re-solved for that closure), so it cannot be bit-identical with the spread study's
`current-wE/s0.10` run.

| check | result |
|---|---|
| setup report | `attached, wave-driven, strength history [(1969, 0.1), (2017, 0.0)], ... spread 0, 0.667, 0.333` |
| run metadata | groin `kind` wave-driven, `spread_weights` [0, 2/3, 1/3] |
| groin_diagnostics.csv | `spread_loss_m` 4.6-8.7 m every year 1997-2009 |
| end-of-window change, GIS 3-12 (m, seaward +) | -4 +15 +19 -20 -18 -18 -4 -3 -11 -3; the study run -4 +15 +19 -20 -18 -18 -4 -3 -12 -3 |

A first attempt stopped writing the metadata: the spread weights went in as a 3-tuple, which the metadata writer
reads as (value, comment); they are now passed as a list.
