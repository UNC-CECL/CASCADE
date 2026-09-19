# groin - the Buxton groin field

A study of its own, with its own conventions, kept apart from the hindcast
trees on purpose.

```
GROIN_PLAN.md               the authority on the M and f fit
HAT-groin-buxton-input/     the structure, the surveys, the fillet record
HAT-groin-buxton-output/    the 1967-2017 runs behind the fit
HAT-groin-gis-analysis/     the GIS work behind the extents
HAT-groin-figures/          the figures
HAT-buxton-hindcast-groin-test/  the groin inside the hindcast
groin-module-test/          the solver audit
```

**Its model output is not here.** The (M, f) sweep, `joint_fit.json` (read by
`HAT_run_all.py` stage 6) and the `SELECTED_M60_f0.60/` record are in
`output/calibration/groin/`. The 1967-2018 rig runs are in `output/calibration/groin_rig/`.

## Where it departs from the rest of the project

**It keeps its data beside its code**, unlike `data/hatteras_init/`, and its
folders are hyphenated with a `HAT-groin-` prefix that matches nothing else.
Both follow from this being a separate study, and neither is written down
anywhere else, which is why it is written down here.

## What to be careful of

The trapping rate and the deterioration fraction are **fitted from different
evidence** and must be quoted as a pair, with `GROIN_PLAN.md` cited for why.
The two routes to them share topography, wave climate and physics, so they are
independent WINDOWS, not independent evidence.
