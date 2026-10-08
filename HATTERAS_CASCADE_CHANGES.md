# Changes to the CASCADE model code for the Hatteras hindcast

This is the companion to `HATTERAS_FIXES.md` on the local Barrier3D branch `hatteras/adopted`, which records the Barrier3D side. This file covers the `cascade/` package in this repository. For each change: what was wrong or missing, how it showed up, the change, its measured effect, and the evidence. It then lists:
- accidental edits still in the code,
- behaviour of unchanged CASCADE code that the hindcast works around from outside,
- upstream work this branch does not have.

Written 2026-09-29. The pipeline (`scripts/`), data and figures are not covered here; see `ORGANIZATION.md` and `UNITS.md`.

## What the changes are measured against

This branch (`hannahaline/hatteras-cascade`) began as an orphan commit, f409654c "Fresh start: scripts only", so git has no merge base with upstream. Two reference points are used:

| reference | commit | why |
|---|---|---|
| **published** | 8549847a (2024-02-09, "Updates to version references for publication of Parts I & II") | the upstream commit closest to where this branch began; the hatteras-nps, PeaIsland and benton branches all start from it |
| **upstream now** | `origin/main` 6e896375 (2026-09-11, PR #58 beach_width_mods) | what upstream has today |

Against *published*, eight files differ. Relative to that version, `check_sandbag_need`, the sandbag arguments and `BrieCoupler.offset_shoreline` look new, but they are **upstream code**: upstream `main`, ocracoke-nps and peaisland-hindcast all have them. They are not listed as Hatteras changes.

The fresh-start commit already carried some local edits, made before this branch's history began. Those are listed as found: sections 5 and 6.

## Summary

| # | change | files | kind | commit(s) | used by the hindcast |
|---|---|---|---|---|---|
| 1 | Groin representation: an alongshore source/sink dipole injected before BRIE's transport solve | `groin.py` (new), `cascade_groin.py` (new: `cascade.py` + a hook) | new capability | 058cc821, 5b26cb54 (+ pre-history) | yes, `cascade_groin.Cascade` is the class the runner builds; the groin itself is off in the current matrix (`nogroin`) |
| 2 | Where a relocated road is rebuilt is its own input (`road_relocation_setback`) | `cascade.py`, `cascade_groin.py` | interface fix | ead66e13 | yes |
| 3 | `resize_interior_domain` aligns pre/post-storm grids by anchoring the edge that did not move | `beach_dune_manager.py` | crash fix | 9b89dfee | yes, every managed run |
| 4 | The bulldozer dune cap limits only the sand added, not the whole dune cell | `beach_dune_manager.py` | bug fix | 4d9c3acd | yes, every managed run |
| 5 | Find-and-replace accidents in names, labels and a default file name | `cascade.py`, `cascade_groin.py`, `beach_dune_manager.py`, `roadway_manager.py`, `brie_coupler.py`, `tools/plotters.py` | defect, **all reverted 2026-09-29** | pre-history | no effect on Hatteras runs (see 5) |
| 6 | Empty `res_manager.py`; debug `print`s; formatting drift | various | cruft, **removed 2026-09-29** | pre-history | harmless |

---

## 1. Groins (`groin.py`, `cascade_groin.py`)

**What was missing.** CASCADE has no hard structures. The Buxton groins trap alongshore sand, building an updrift fillet and a downdrift notch, and nothing in BRIE or Barrier3D can represent that.

**The change.**
- `cascade/cascade_groin.py` is a copy of `cascade.py` with one addition: just before BRIE's alongshore-transport solve, `update()` calls `cascade._groin_callback(self, x_s_dt)` if one is attached (`cascade_groin.py:608-618`). With no callback, the model is bit-for-bit the same as `cascade.py`.
- `cascade/groin.py` provides `GroinCallback`: each year it adds −M (m) to the updrift domain and +M to the downdrift one. BRIE's diffusion then spreads that dipole into an emergent fillet and notch. The only knob is the trapping rate M, plus the fraction f in the runner. It also provides `predict_fillet`: the steady-state amplitude A ≈ M / (4 r_ipl), and the extent L ≈ dy·√(2 r_ipl t), which M does not affect. So fitting M to the observed amplitude leaves the extent as an independent test.

**Effect and status.**
- Implemented and fitted (joint two-period fit, M = 60 m/yr, f = 0.6 in the driver's `joint_fit.json`).
- **Not used in the current matrix**, which is all `nogroin`.
- Open problems:
  - M at that size moves about half the reach's annual sediment budget past one groin, so it fails a budget check.
  - The fillet depends on M / r_ipl, and r_ipl scales with Hs, so a fitted M only holds at the Hs it was fitted at.
  - BRIE runs away past about −23° of shoreline angle.
  - BRIE's alongshore solve does not conserve volume. The shoreline score is demeaned, which is what makes that drift harmless.
- The "groin progradation ceiling" of 2026-09-11 was the Barrier3D `route_overwash` crash (`HATTERAS_FIXES.md` §1), not a groin limit.

**Evidence.** Module docstrings; `scripts/hatteras_ms/groin-sweep/`; the solver audit `hard-structures/groin/2-module-tests/2-solver-audit/HAT_groin_solver_audit.py`; `output/calibration/groin*/`.

**Upstreaming note.** `cascade_groin.py` duplicates about 850 lines of `cascade.py`, and every model change has to go into both (change 2 did). Upstream, the hook should go into `cascade.py` itself; `GroinCallback` then needs no separate class.

## 2. `road_relocation_setback` (ead66e13)

**What was wrong.** One `road_setback` did two jobs. It placed the road at t = 0, and the yearly update re-read it as the distance a relocated road is rebuilt at. A caller wanting a standard rebuild clearance had to overwrite `cascade._road_setback` after construction. That only worked because the constructor had already used the value.

**The change.** `Cascade(..., road_relocation_setback=None)` in both classes. The yearly update reads it (`cascade_groin.py:713`). `None` keeps the old behaviour exactly: each domain relocates to its own measured offset.

**Effect.** Proven a no-op for existing runs: same arm and topography, every domain's rate matches to zero and the road-management table is identical row for row. It removes the post-construction overwrite from the pipeline.

**Evidence.** Commit message of ead66e13.

## 3. `resize_interior_domain` (9b89dfee)

**What was wrong.** `BeachDuneManager` differences the pre- and post-storm interior grids. When their cross-shore lengths differed, the old code tried to rebuild the seaward shift from `dune_migration` (`ShorelineChangeTS`), then trimmed trailing all-bay rows. Both halves were wrong:
- The shoreline change is not the number of interior rows lost. In one Hatteras run it read −1.0 at three steps where the front had lost 2, 3 and 2 rows.
- The trim required the excess rows to be entirely at bay depth, which a partly filled back-barrier row is not.

**How it showed up.** When neither branch closed the gap, the function returned mismatched arrays. `filter_overwash` then died about 140 lines later on a numpy broadcast error that named neither the function nor the cause.

**The change.** Anchor on the edge that did not move:
- **Seaward rows lost:** keep the last `n_post` rows of the pre-storm grid.
- **Bay rows gained:** pad the pre-storm grid with bay cells.

A mismatch in alongshore length now raises a `CascadeError` that names itself. `dune_migration` is kept in the signature and no longer read.

**Effect.** On runs that already completed, all 45 `pre > post` calls give byte-identical arrays, so no existing result moves. It resolves only the cases that used to crash.

**Known limit.** A step that loses seaward rows *and* gains bay rows at once has no single anchor, and will be mis-aligned by the rows gained. This was not seen in any Hatteras run (48 of 48 length changes were one-ended), and it cannot be detected from the two grids alone.

**Evidence.** The function's docstring.

## 4. The bulldozer dune cap limits only added sand (4d9c3acd)

**What was wrong.** In managed domains, `filter_overwash` puts `overwash_to_dune` (9% for Hatteras) of the year's overwash onto the dunes. It then clipped the **whole dune cell** to `_artificial_maximum_dune_height = 4` m above the berm. That value is fixed in the code as a Nags Head value and is not a constructor argument. The clip runs every year in every domain the manager covers, because `_overwash_removal` is always True. The cap exists to stop the bulldozed sand building 10 m dunes; clipping the whole cell also cut natural dunes.

**How it showed up.**
- It never bound while Barrier3D's `Dmaxel` default held dunes near 3 m MHW.
- After the per-cell ceilings (`HATTERAS_FIXES.md` §6) it cut 24 domains in 1996-2010 and 17 in 2010-2024, about 150,000 m³ of starting dune each time.
- Avon and Tri-Village crests sat flat at 5.34 m MHW, below the lidar.

**The change.** The new height is `max(post, min(post + added, cap))`: sand is added up to the cap, and a cell already above it keeps its height and takes no sand. As before, any sand beyond the cap is lost. `DUNE_CAP_APPLIES_TO = "added sand only"` is recorded in every run's metadata (`bdm_dune_cap_applies_to`).

**Effect** (managed runs only):
- 2010-2024 RMSE: 2.12 → 2.07 (full management), 2.26 → 2.20 (beach/dune only).
- 1996-2010 runs' 2010 crest against the lidar: −0.31 → −0.17 m.
- Overwash scores unchanged.
- The 2010 GIS 90 end was re-solved: +22.4937 → +21.2582.

**Evidence.** `output/comparisons/adoption_2026-09-28/README.md` ("The dune-cap fix"); `output/raw_runs/experiments/end-domain-boundaries/2026-09-28-ends-resolved-dunecap/README.md`; superseded runs in `output/raw_runs/archive/2026-09-28-pre-dunecap/`.

## 5. Find-and-replace accidents (defect, reverted 2026-09-29)

A project-wide rename (most likely an IDE refactor of folder names) rewrote words inside the model code. Found 2026-09-29 by diffing against the published version:

| original word | became | where |
|---|---|---|
| `storms` | `original` | the **default storm file name** `cascade-default-storms.npy` → `cascade-default-original.npy` (`cascade.py`) and `cascade-default-old.npy` (`cascade_groin.py`); the matching error message; comments and docstrings in `roadway_manager.py`, `brie_coupler.py`, `tools/plotters.py` |
| `topography` | `topography_dunes` | docstrings in `beach_dune_manager.py`, `roadway_manager.py` |
| `Dunes` / `Topo` | `PEA_2011` | **plot labels** "PEA_2011 Rebuilt" and title "Cross-shore PEA_2011 Transects" in `tools/plotters.py` |
| `output` | `comparison` | comments and a docstring in `cascade.py`, `tools/plotters.py` |
| `offset` | `raw_offset` | one comment in `cascade.py` |

**Effect on the hindcast: none.** The runner always passes `storm_file` explicitly and does not use those plotters. But a caller relying on the default storm file would get a file-not-found error, and the default-storms guard now tests a name no file has.

**Fixed 2026-09-29: the default storm file** (Hannah: "fix 1a and 2 now"). The default, the default-storms guard and its message say `cascade-default-storms.npy` again in both classes. The guard had never been able to fire. Hatteras results are bit-identical (`output/raw_runs/experiments/code-checks/2026-09-29-default-storms-sandbag-fix/NOTE.md`).

**The rest reverted 2026-09-29** (Hannah: "revert the rest of the find-and-replace accidents"): 24 lines in all, 21 in this pass:
- the two plot labels and the title in `plotters.py`
- ten docstrings (`topography_dunes`)
- the `storms`/`output`/`offset` comments

A line was reverted only if its sole difference from the published version (8549847a) was one of the known substitution pairs; lines with real edits were left alone. The legitimate uses of "original" (`original_growth_param`, "back to original") and "the old behaviour" (change 2) are untouched. Comments, docstrings and plot labels only, so no model behaviour changes. `grep` finds no remaining `topography_dunes`, `PEA_2011` or `raw_offset` in `cascade/`.

## 6. Cruft (removed 2026-09-29)

Hannah: "remove the leftover cruft too". Hatteras results bit-identical (`output/raw_runs/experiments/code-checks/2026-09-29-cruft-cleanup/NOTE.md`).

- `cascade/res_manager.py`, an empty file referenced by nothing and absent upstream: **deleted**.
- `print`s in `roadway_manager`: most of them ("Roadway relocated", "Elevation is low enough for sandbags", "Sandbags would be added…") turned out to be **upstream code**, so they are kept. Removing them would widen the gap with upstream. The two local-only prints ("Road close enough for sandbags", "Road far enough away") and two commented-out debug prints are **removed**, and the sandbag message is back to upstream's shorter form.
- Formatting drift: **black** (upstream's pre-commit formatter) run on the drifted modules. `groin.py` was left alone because another session is editing it.
- **Fixed earlier the same day:** `sandbag_management_on` was not broadcast from a single value, so the default `False` crashed in year 1. It is now broadcast to every domain, as upstream f7ad676b does. Not hit by the hindcast, which passes a list; results bit-identical (`code-checks/2026-09-29-default-storms-sandbag-fix/`).
- Still different from upstream, and **not cruft**: `check_sandbag_need`'s logic (threshold 0.08 vs upstream 0.101, all dune rows vs row 0). Sandbags are off in the Hatteras runs.

---

## Behaviour of unchanged CASCADE code that the hindcast works around

None of these were changed in `cascade/`. The pipeline or the reporting handles each one. Listed so that nobody "fixes" the pipeline side without knowing why it is there.

| behaviour | where | consequence | how the hindcast handles it |
|---|---|---|---|
| `update()` copies its own `_nourishment_volume[iB3D]` and `_dune_design_elevation[iB3D]` into each `BeachDuneManager` just before calling it | `cascade_groin.py:802-809` | a volume written onto `cascade.nourishments[i]` is overwritten, and every fill silently used the 100 m³/m init default (Avon 2022: 100 instead of 841 m³/m) | volumes are written to the Cascade-level list by `cascade_pipeline/nourishment.py` (`NourishmentSchedule.apply_to_cascade`); `verify_nourishment()` reads the manager's own `_nourishment_volume_TS` to confirm |
| the beach/dune manager runs on the union of community zones and nourished domains, all period | runner footprint | where it overlaps the roadway footprint, both managers run: overwash removed twice and `ShorelineChangeTS` pinned at 0 | reported, not corrected; asserted against the finished run |
| a fill is an instantaneous step in `x_s`, and BRIE's Crank-Nicolson solve rings at the grid scale (r ≈ 1.05) | `shoreface_nourishment`; `brie.py` ~1294 | alternating ±m wiggles next to recent fills, decaying ~40%/yr | model rates are an OLS slope (LRR), which averages the ringing |
| BRIE's alongshore solve is not volume-conserving | `brie.py` | the reach mean drifts; matters for the groin dipole | the shoreline score is demeaned |
| `drown_threshold = 0` is commented "m MSL" but compared against MHW-relative elevations | `roadway_manager.py:766` | effectively 0 m MHW | documented in `UNITS.md` |
| the roadway dune rebuild trigger is floored at `BermEl*10 + 0.3` m and the design height at `BermEl*10 + 1.0` m | `roadway_manager.py:530, 719-722` | a passed `dune_minimum_elevation` of 0 does nothing | documented in `UNITS.md` |
| `beta` (0.04 default) is not the storm run-up slope (0.06) | `Cascade(beta=...)` not passed | affects only the initial beach width in the beach/dune manager | documented in `UNITS.md`; left as is |

## Upstream work this branch does not have

Upstream `main` (6e896375) has moved on since this branch's code was taken. Present upstream, absent here:
- `outwasher.py` and the `outwash_module` (Benton), with its storm and beach files and `percent_washout_to_shoreface`.
- The user-defined and shared beach-width work (PR #58 `beach_width_mods`: `use_defined_beach_width`, `user_inputed_beach_width`, `beach_width_threshold`, and one `beach_width` shared across modules).
- `allow_causeway` in `roadway_manager.bulldoze` (a road not drowned when surrounded by water).
- `build_interior_dunes` (NCDOT right-of-way dune building, from Ocracoke).
- `set_specified_variable_RSLR` in `brie_coupler`.
- `road_relocation_setback` upstream is a plain number (default 30). This branch's version (change 2) defaults to `None`, meaning "each domain's own setback".
- `sandbag_management_on` broadcast from a single value (f7ad676b). **Equivalent fix applied here 2026-09-29.**

Merging upstream would bring these in and would need the Hatteras changes (1-4) re-applied on top. Upstream's `road_relocation_setback` and this branch's version would collide.

## When upstreaming (not yet: Hannah's call)

1. Revert the find-and-replace accidents (5) and the cruft (6).
2. Offer 3 (the `resize_interior_domain` crash fix) and 4 (cap on added sand) as bug fixes. Both are small and self-contained.
3. Offer the groin hook (1) inside `cascade.py`, not as a second class, with `groin.py` as an optional module.
4. Reconcile 2 with upstream's own `road_relocation_setback`.
5. Do this together with the Barrier3D side (`HATTERAS_FIXES.md` on `hatteras/adopted`, still local): the hindcast needs both.
