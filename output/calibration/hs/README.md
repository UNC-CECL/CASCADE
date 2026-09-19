# Does a higher Hs need less source/sink correction?

**The question.** The calibrated source/sink field exists to supply what the
physics does not. If raising the wave height makes the model reproduce the
observed shoreline change on its own, the field it needs should be *smaller*,
*narrower*, or *in fewer places*. This measures whether it is.

**Status:** CONCLUDED 2026-09-01 -- keep Hs = 2.5. See `DECISION.md`,
including the edge-re-solve addendum, which closed the frozen-edge
confound and turned period 1 from a marginal gain into a clear loss.

---

## What is being compared

Two calibrations of the SAME model, differing only in wave height, each derived
by the project's own pass-0 method:

| arm | base run | Hs |
|---|---|---|
| control | `edgeBE road_bdm groin-on`, M = 60, f = 0.6 | 2.5 m |
| test | the same run, same groin | 3.0 m |

`edgeBE` because it imposes background erosion only at the two end domains, so
the interior residual **is** the correction the calibration would have to
invent. Groin held at the production (M, f) in both arms, so any change is the
wave height and not the structure.

## Why pass-0 only, and no iteration

`convergence_history.json` records that zone membership is identified ONCE, from
the pass-0 residual, and frozen. Re-deriving it each pass was tried and
abandoned: it let incoherent features cross the 0.5 m/yr threshold, including
the groin's own footprint at D5-D7 in period 2, and it does not terminate.

So the honest comparison is between the two arms' **pass-0 zone
identification** -- "if we calibrated from scratch at this Hs, what would pass 0
select?" -- and nothing further. No iteration, no `--add`, no re-derivation.

## What is NOT being done

**The production calibration is not touched.** It is converged (2026-08-24) and
still inside its own 5% stopping rule against the current runs -- drift is +0.7%
(edgeBE 1984), 0.0% (edgeBE 2004), +3.2% and +4.1% (calibBE). Every run here
sets `HAT_BE_OUTPUT_DIR` to this directory, so `be_zone_metrics.csv`,
`DOMAIN_BE_RATES*.txt` and `convergence_history.json` under
`data/hatteras_init/7-source-sink/2-calibrate/` are never written.

**The groin is not refitted.** M = 60 was fitted at Hs 2.5, and raising Hs
raises alongshore diffusivity, so the fitted M would move. Holding it fixed
confounds nothing *here* -- it isolates the wave height -- but it means this
experiment cannot say what the groin should be at Hs 3.0. That is a separate
question, and `GROIN_PLAN.md` already records a partial answer: at Hs 3.5 the
fillet decay improves from -15.3 m to -20.6 m against -24.4 m observed, at the
cost of the reach RMSE going 15.97 -> 36.58.

## Metrics

Agreed before the runs, so the answer cannot be chosen after seeing them:

1. **RMS of the required correction** -- how large.
2. **Split by reach** -- the totals hide that the ends and the middle move in
   opposite directions, which is the finding, not a detail.
3. **Domains selected by pass 0** -- how many zones must be invented, under the
   analysis's own rule (`|smoothed residual| > 0.5 m/yr`, contiguous runs of at
   least 3 domains).

## Layout

```
logs/             the two Hs 3.0 base runs' logs, and one per pass-0 arm
02_zones_Hs2p5/    pass-0 calibration, control arm
03_zones_Hs3/      pass-0 calibration, test arm
comparison/        the three metrics, side by side
DECISION.md        what we concluded, what would overturn it, and the
                   2026-09-01 edge-re-solve addendum
```

Model output itself stays in the run registry, not here. One tree per ARM, the
component being `HAT_ARM_TAG` in the runner (and the `arm` column in
`run_index.csv`, which carries the same string):

```
output/raw_runs/<period>/edgeBE/                  control, Hs 2.5
output/raw_runs/waveHs3/<period>/edgeBE/          test, edges from the 2.5 solve
output/raw_runs/waveHs3_probe/<period>/edgeBE/    +3.0 m/yr gain probes
output/raw_runs/waveHs3_edge1/<period>/edgeBE/    test, edges re-solved
```

The last two are the addendum's. They exist as separate arms rather than as
overwrites because a probe shares its name, preset and wave climate with the
run it probes -- so without an arm of its own it would derive that run's
directory, and `guard_run_dir` would refuse (or, worse, overwrite the baseline
the probe is measured against).

## Known limits, recorded before the result

- One wave height against one control. Two points do not make a curve, and the
  interior optimum in the managed sweep sat near Hs 1.2, not above 2.5.
- The zone rule has a hard threshold at 0.5 m/yr, so a domain sitting at 0.49
  and one at 0.51 are counted very differently. Small changes in the residual
  can move the count without meaning much.
- Period 2's residual is dominated by a near-uniform +0.99 m/yr offset that the
  wave height barely moves. Expect its zone count to be insensitive.
