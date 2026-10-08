# 1-dipole-1967-2017 — the dipole groin fit (SUPERSEDED)

> **Superseded 2026-10-05** by the blocking groin (`../2-blocking-1996-2025/`). Kept as the record of how the dipole (`cascade.groin.GroinCallback`, a fixed −M / +M pair) was fitted. Not maintained: some scripts read runs that are now on the offload drive.

`dipole_fit_notes.md` is the reasoning: how M 60 and f 0.6 were chosen, what the dipole cannot do, and the corrections along the way.

```
dipole_fit_notes.md   the fit, its limits and corrections (was GROIN_PLAN.md, sections 2-5)
inputs/               the 1967 rig's inputs: island offset for GIS 2-12, storms, the dune-line
                      target, and the scripts that make them
runs/                 the rig runner, the (M, f) sweep and its worker, the edge solve, the
                      fillet-trajectory target, plots; the six 1967-1997 PNGs left on disk
results/              the rig's sweep results and the three-run comparison
figures/              groin_timeline_and_hindcast: the fillet against the dipole's M_eff schedule
```

The rig's runs are in `D:\CASCADE_offload\output\calibration\groin_rig\` (moved 2026-10-07). The 1984-2024 hindcast sweeps that pinned M 60 are in `D:\CASCADE_offload\output\calibration\groin\`, and their scripts are in `scripts/hatteras_ms/groin-sweep/`.
