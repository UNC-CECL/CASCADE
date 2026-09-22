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
