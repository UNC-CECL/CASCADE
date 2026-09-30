# 3-storms — what is here, and the max-duration choice

The storm series the model runs on, and the two validators' output.

```
hindcast_storms/   the model input, one folder per window
validation/        the validators' tables and figures
figures/           the storm record figures
```

## The max-duration test

**An event longer than `max_storm_dur` hours is not included as a storm.** The
limit was chosen by testing, with Lexi's script, which values the generator
would actually run to completion:

```
Storm testing with Lexis script

Mx hr limits
- 36 worked
- 60 worked
- 72 worked (moving forward with)
- 96 did not
- 120 did not
- 240 did not
```

That is the original note, verbatim. It was a file called `Notes`, with no
extension, until 2026-09-22 — the only record of the test, and invisible to
every tool that looks for documentation by suffix.

**72 is what the generator uses**, and it is not a free parameter any more:

```python
# scripts/input_prep/3-env-forcings/3-storms/historical_storm_creation_v3_HAT.py
max_storm_dur = 72      # maximum duration to include in storm events [hrs]
save_name = "{0}_storms_v3_72".format(PERIOD_TAG)
```

The value is in the **output filename** (`<window>_storms_v3_72`), so a series
built at a different limit cannot silently overwrite this one. If you change
`max_storm_dur`, the name changes with it — that is deliberate, and it is why
the existing files can be trusted to be the 72 h series.

**What "did not work" means is not recorded**, and the note does not say. It
reads as the generator failing or not terminating at 96 h and above rather
than as a judgement about the physics, but that is an inference from the
wording, not something the note states. Treat the boundary between 72 and 96
as an observed limit of the script, not a coastal one, unless it is re-tested.

## 2026-09-28: the 72 h limit was a crash workaround; `v3_trim24` keeps every event

**Why the longer limits "did not work".** Barrier3D before the route_overwash axis-swap fix (commit 49fd069, 2026-09-24) read out of bounds during overwash routing. A long storm routes through enough steps to reach that read, and the process died without a message; Isabel 2003 did it in model year 8 of a 240 h series. On the fixed code every limit runs (`output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-max-duration/`).

**What the 72 h limit cost.** The builder DROPS an event longer than the limit rather than shortening it. At 72 h that removed 29 events from 1996–2024, among them Isabel 2003 (Rhigh 5.22 m MHW, 133 h after grouping), March 2018 (4.30 m), Florence 2018, November 2010, Dennis 1999 and Nor'Ida 2009.

**The adopted rule (Hannah, 2026-09-28): `--long-events trim --max-duration 24`, files `<window>_storms_v3_trim24`.** Every event is kept. One longer than 24 h is cut to the 24 h above the berm centred on its peak total water level, with Rhigh, Rlow and the period taken from what is kept. It was chosen on overwash against the imagery:

- `experiments/storms-and-overwash/2026-09-28-storm-length-selection/`: trim length barely changes where and when overwash happens, but longer trims add retreat.
- `.../2026-09-28-dune-ceiling-per-domain/`: with realistic dunes, the long storms are what make 1996–2010 match the imagery.

**Check.** The builder's defaults (`drop`, 72) rebuild the `v3_72` files exactly, and the `v3_trim24` files equal the ones those experiments ran. Built for all four windows. The `v3_72` files are kept, so earlier runs stay reproducible. Which file a run reads is set by `hat_env_forcings.storm_series_file`'s default variant.

## 2026-09-29: split back-to-back storms, `v3_split12_trim24`

**Why.** The 24 h grouping chains storms a week apart, because the berm is overtopped at most high tides in an active spell. Edouard + Fran 1996 became one 119 h event and Jose + Maria 2017 one of 186 h. `v3_trim24` then kept only the 24 h around each merged event's highest peak, so **Fran (2.62 m MHW) and Jose (2.81 m) were not in the model**. Across 1996–2024, 56 spells of ≥8 h above the berm were cut out of merged events.

**The adopted rule (Hannah, 2026-09-29): `--long-events trim --max-duration 24 --split-gap 12`, files `<window>_storms_v3_split12_trim24`.**

1. Group as before (spells <24 h apart).
2. Split each group wherever the water stays below the berm for ≥12 h. A piece shorter than the 8 h minimum is folded into the piece before it (after, for the first), so no hour the group counted is lost.
3. Each piece is an event of its own, dated by its start and trimmed to the 24 h around its peak.

**Evidence** (`output/raw_runs/experiments/storms-and-overwash/2026-09-29-event-splitting/`): on the adopted model, managed, both windows, the overwash scores (PSS, POD, POFD, timing and space r), interior RMSE and bias are unchanged to two decimals. Overwash volume rises 2–3% in 2010–2024. It is adopted as the more faithful record, not for skill.

**Check.** Without `--split-gap` the builder rebuilds `v3_72` and `v3_trim24` exactly, and with it the files equal the experiment's `split12`. Built for all four windows:

| window | events | storm-hours (split12_trim24) | events | storm-hours (trim24) |
|---|---|---|---|---|
| 1984–2004 | 162 | 3,023 | 155 | 2,936 |
| 1996–2010 | 159 | 2,981 | 150 | 2,839 |
| 2004–2024 | 261 | 4,889 | 245 | 4,649 |
| 2010–2024 | 191 | 3,542 | 178 | 3,345 |

`hat_env_forcings.DEFAULT_STORM_VARIANT` is `v3_split12_trim24`. The `v3_trim24` and `v3_72` files are kept so earlier runs stay reproducible.
