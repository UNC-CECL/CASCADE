# `calibration/hs/runs` — the Hs 3.0 arms, moved out of `raw_runs/`

These seven runs are the evidence behind `../DECISION.md` ("keep Hs = 2.5",
concluded 2026-09-01). They were moved here from `output/raw_runs/` on
2026-09-02. Nothing about them changed: the directories, their contents and
their index rows are exactly as the runner wrote them.

## Why they are not in `raw_runs/`

A **forcing arm** is a parallel tree of runs made under forcing that differs
from the calibration, so a run name describes the SCENARIO and the arm
describes the FORCING. That means an arm run shares its name with the
production run it is compared against — and it must, since it is the same
scenario. While these lived in `raw_runs/`, three names existed there four
times each, and `raw_runs/` no longer had one directory per run name.

That is not hypothetical damage. `HAT_plot_sensitivity.load_index()` collapsed
the repeats with `drop_duplicates(keep="last")`, which does not pick the
calibration run — it picks whichever row sorts last. Every edgeBE wave cell was
drawn against `waveHs3_probe`, a Newton probe at Hs 3.0, instead of against its
Hs 2.5 baseline. Six cells were plotted as beating a baseline they do not beat.

So: **a live calibration probe is still written into `raw_runs/` under
`HAT_ARM_TAG`** — that is what stops it overwriting the base run it probes —
**and is moved here when the question it was answering is settled.** `raw_runs/`
holds the production matrix and its sensitivity cells, and one run name means
one directory there.

## What is here

The layout is identical to `raw_runs/`, so the same resolver reads both:

    <arm>/<start>_<end>/<preset>/<run_name>/

| arm | Hs | be1 / be90 | what it is |
|---|---|---|---|
| `waveHs3` | 3.0 | −42.60 / 13.00 | the test arm, edge values frozen from the 2.5 m calibration |
| `waveHs3_edge1` | 3.0 | −51.99 / 27.66 | a Newton probe re-solving the two locked end domains |
| `waveHs3_probe` | 3.0 | −39.60 / 16.00 | a second probe in that solve |

The probes are **steps in a calibration, not results.** Do not read a probe's
skill number as the model's skill at Hs 3.0; the two `waveHs3` rows are the
runs that answer that.

`run_index.csv` here carries these seven rows and only these. They were split
out of `output/raw_runs/run_index.csv` in the same move, which is archived
alongside it as `run_index_archive_20260902_002606.csv`.

## Reproducing them

`be_zone_residual_fit.py` resolves this root automatically: with
`HAT_BE_HS` at the calibration value it reads `raw_runs/`, and off it, this
directory. Nothing needs a path typed in.

## What git tracks here

Since 2026-09-18 this README and `../DECISION.md` are tracked (the
`output/calibration/` block in `.gitignore`), so the Hs = 2.5 decision survives a
fresh clone. The arm runs themselves are still untracked and exist only on this
machine.
