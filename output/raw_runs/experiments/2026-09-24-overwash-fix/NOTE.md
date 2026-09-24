# 2026-09-24 — how much does the Barrier3D route_overwash fix move the results?

## The bug and the fix

`Barrier3D/barrier3d/barrier3d.py`, `route_overwash`: the test deciding whether
overwash flux keeps moving landward or decays into the bay read
`Elevation[TS, i, d+1:d+10]` (row `i`, columns `d+1…d+9`) where the nine cells
landward of the flow, `Elevation[TS, d+1:d+10, i]`, are meant. Wrong cells
whenever `i` < rows; out of bounds whenever `i` ≥ rows (silent in the jitted
code, and the cause of every crash in `sensitivity/2026-09-24-natural-waves/`
and of the 2026-09-11 groin crash). Found and explained in that study's README.

The fix is one line, commit `49fd069` on the local Barrier3D branch
`fix/route-overwash-axis-swap` (not pushed). `master` (`ce36866`) is unchanged
and **is what is checked out now**: every run made from here uses the
unpatched model until the fix is adopted.

## What was run

Script: `scripts/hatteras_ms/experiments/HAT_overwash_fix_check.py`. Every
launch recorded the Barrier3D branch and commit it ran on (`launches.jsonl`).

| member | twin (unpatched) |
|---|---|
| natural and full-management baselines, 1996–2010 and 2010–2024 (metres, zeroBE, no groin; Hs 1.0, Tp 8, asymmetry 0.8, high-angle 0.45) | the natural-waves study's runs, made today on today's code |
| natural 1996–2010, Hs 1.25, high-angle 0.5 (the study's best ridge point) | the study's grid run |
| full management 2010–2024, high-angle 0.4 | none: it crashed in year 13 unpatched |
| ÷10 full management 1996–2010 and 2010–2024 (the old matrix settings: Hs 2.5, pre-metres offset build) | re-run on `master` today (`unpatched_*`) |

Plus the natural 1996–2010 and full-management 2010–2024 baselines on the fix
with `NUMBA_BOUNDSCHECK=1`.

The ÷10 re-runs on `master` reproduce the archived matrix runs
(`archive/2026-09-24-pre-metres/matrix/…/HAT_1996_2010_zeroBE_road_bdm_nogroin`,
`…/HAT_2010_2024_zeroBE_road_bdm_nourish_nogroin`) exactly, domain for domain:
nothing else in the code has moved results since, and the unpatched model is
deterministic despite its out-of-bounds reads.

## Results (`comparison.csv`)

| run | RMSE unpatched → patched | bias | share explained | per-domain rate change (RMS / max), domains moved > 0.1 m/yr |
|---|---|---|---|---|
| natural 1996–2010 baseline | 1.141 → 1.143 | −0.444 → −0.448 | 5.7% → 5.3% | 0.007 / 0.025, 0 |
| managed 1996–2010 baseline | 1.070 → 1.070 | −0.044 → −0.045 | 16.9% → 16.9% | 0.003 / 0.018, 0 |
| natural 1996–2010 ridge | 1.099 → 1.100 | −0.361 → −0.362 | 12.4% → 12.3% | 0.005 / 0.020, 0 |
| natural 2010–2024 baseline | 4.874 → 4.870 | −4.468 → −4.465 | −1008% → −1006% | **0.249 / 2.160, 14** |
| managed 2010–2024 baseline | 2.588 → 2.584 | −1.870 → −1.863 | −212% → −211% | 0.026 / 0.105, 1 |
| managed 2010–2024, high-angle 0.4 | crashed (year 13) → **2.678** | → −1.971 | → −235% | — |
| ÷10 managed 1996–2010 | 1.071 → 1.071 | −0.228 → −0.228 | 16.9% → 16.9% | 0.000 / 0.002, 0 |
| ÷10 managed 2010–2024 | 2.629 → 2.627 | −2.030 → −2.028 | −222% → −222% | 0.008 / 0.028, 0 |
| bounds check on the fix | — | — | — | natural 1996 and managed 2010 baselines: **all 14 years, no IndexError** |

- **The fix removes the crashes and the out-of-bounds reads.** The crashed
  cell now runs 14 years; the bounds-checked runs find nothing else out of
  bounds.
- **Its effect on scores is negligible** in every run compared: RMSE changes
  in the third decimal, share explained by under half a point.
- **The one visible change is local**: natural 2010–2024, where overwash is
  heaviest, moves 14 domains by more than 0.1 m/yr and one by 2.2 m/yr, with
  the island-wide scores unchanged. The wrong-cell reads mostly returned the
  same subaerial answer the right cells would have.
- **So no earlier conclusion changes**, ÷10 or metres. What the bug did cost is
  the crashes, which removed runs from the sweeps (10 in the natural-waves
  study) and were misread on 2026-09-11 as a groin progradation ceiling.

## Adopted (2026-09-24)

Hannah: use the fixed version from now on, keep the original untouched, and
leave the branch local for now.

- Barrier3D now has **`fix/route-overwash-axis-swap` checked out** (`49fd069`);
  `master` (`ce36866`) is the untouched original. Barrier3D is installed
  editable, so the checked-out branch is the model. The branch is **local
  only**: not pushed anywhere.
- Every run now records which Barrier3D it used
  (`run_registry.barrier3d_provenance()`): in its metadata identity
  (`barrier3d_branch`, `barrier3d_commit`, `barrier3d_dirty`,
  `barrier3d_route_overwash_fix`, the last read from the source of the module
  actually imported) and in `run_index.csv` (`barrier3d_commit`,
  `barrier3d_route_overwash_fix`). The runner prints the Barrier3D branch and
  the fix flag at start and warns, without refusing, if the fix is missing, so
  an old run can still be reproduced on purpose. `provenance_check/` is the run
  that confirmed it.
- **A blank in those columns means the run predates the recording**: made on
  unfixed Barrier3D, except the ten `patched*` runs here, which ran on the fix
  before the field existed (their Barrier3D commit is in `launches.jsonl`).

## Open
- The upstream issue draft:
  `sensitivity/2026-09-24-natural-waves/barrier3d_route_overwash_issue_draft.md`
  (now with these numbers).
- The 2026-09-11 groin progradation ceiling should be re-tested on the fix.
