# Buxton groin — the dipole fit (SUPERSEDED)

> **Superseded.** This fits the DIPOLE groin (M = 60, f = 0.6), the reference note of 2026-08-24 that was `GROIN_PLAN.md`. Since 2026-10-05 the pinned groin is the BLOCKING groin, b = 0.6, f = 0.6 (`../2-blocking-1996-2025/2026-10-05-blocking-fit-calibration/README.md`), and the straight-coast test (`../../2-module-tests/1-straight-coast/`) found that neither module traps the real drift yet. Kept as the record of how the dipole was fitted. The observed history this note opened with moved to `../../1-observations/structure_history.md` on 2026-10-08. Paths below were updated to the 2026-10-08 layout; the dipole sweeps under `output/calibration/groin/` moved to `D:\CASCADE_offload\output\calibration\groin\` on 2026-10-07.

**Figures.** `figures/groin_timeline_and_hindcast.png` (the fillet against the module's `M_eff` schedule) is beside this note; `../../2-module-tests/figures/groin_module_logic.png` (what the dipole adds each year, and why it cannot close the gap) and `../../1-observations/figures/groin_two_shorelines.png` are with their categories. All three came under the project house style on 2026-09-11. Diffusion alone closes the gap at 1.5% of the observed post-2004 rate, not the "about a tenth" the module-logic footnote claimed; that tenth is roughly right for the NET change across the window, which is a different comparison.

---

## 1. The plan (dipole)

**Turn the groin on for both periods with the schedule it already has.**
`GroinCallback` carries an absolute calendar timeline — install 1969,
deterioration onset 1996, linear ramp to 2003, then hold at `M × f`. That
already reproduces the measured history. **No period-specific configuration.**

| | value | where it comes from |
|---|---|---|
| `trapping_rate_m_yr` | **M = 60** | **fitted on PERIOD 1, window D4–D8, production geometry, be1 pinned at the production value, demeaned score.** RMSE 11.69 vs 15.58 with no groin — the groin closes **25%** of the shape misfit, and this is within 0.10 m of the global best. **Not** corroborated by the 1967 rig: the rig improves monotonically to M = 60 and then blows up (M = 70 → RMSE 320–378, M ≥ 100 crashes), so its M = 60 is the largest value it can hold, not an optimum (*corrected 2026-08-30*). Independently reproduced instead by **D3–D9**, the one other window symmetric about the structure, at the same gain. Intercepts ~719,000 m³/yr, marginally above the 5–7×10⁵ drift band (a literature range, not a hard limit) |
| `deterioration_fraction` | **f = 0.6** | same fit; f is only weakly constrained by period 1 (which mostly precedes the 1996–2003 ramp), so it leans on the rig and on period 2 showing trapping ceased. **Updated 2026-08-30:** the rig, re-run on `1984-start/v1`, now returns **f = 0.6** itself (RMSE 23.78, against 27.24 and f = 0.5 on the pre-fix topography) — and 0.6 is bracketed on both sides there (0.5 → 24.21, 0.6 → 23.78, 0.7 → 25.06). **f is the parameter the rig actually resolves** — it rails in M but not in f |
| `install_year` | 1969 | documented |
| `deterioration_delay_years` | 27 (→ 1996) | first year after the 1995 last repair |
| `deterioration_ramp_years` | 7 (→ 2003) | storm damage |
| `updrift / downdrift` | GIS 6 / 5 | field occupies D6 |

Written to `output/calibration/groin/joint_fit.json`, which stage 6 and
`be_zone_residual_fit.py` both read.

### The fit that supports these values

**PERIOD 1 (1984–2004), window D4–D8, be1 = −40, demeaned score:**

| cell | RMSE | interception | note |
|---|---|---|---|
| no groin | 15.58 | — | |
| **M = 60, f = 0.6** | **11.69** | 719k m³/yr | **chosen** — 25% better than no groin, 0.10 m off best |
| M = 70, f = 0.6 | 11.64 | 838k m³/yr | ~1.3× the drift |
| M = 50, f = 1.0 | 11.59 | 599k m³/yr | best score, but **f = 1.0 = no deterioration** — contradicts the 2003 damage and period 2 |
| M = 50, f = 0.6 | 12.03 | 599k m³/yr | *superseded* — M is low for f = 0.6 |

**M and f trade off along a ridge of roughly constant period-1 trapping**: at
f = 1.0 the best M is 50, at f = 0.6 it is 70. So a low score at M = 50, f = 0.6
is not evidence against f = 0.6 — it means M was set too low for that f.

**Period 1 does respond to f** (3.1 m of RMSE across the f range at M = 50), but
it reads high f as "more trapping" because it mostly PRECEDES the 1996–2003
ramp: period-1 cumulative trapping is `M(15.5 + 4.5f)`, which f moves by only
29% across its whole range. Period 2 is `20·M·f`, where f = 0 gives zero — total
leverage. **Set f from the rig and period 2; set M from period 1.**

**Two things had to be right for the groin to show at all:**

1. **Exclude D1.** The cape's error over 1984–2004 is 81–104 m and swamps a
   ~17 m groin signal. On the full D1–D12 window with a raw score, no-groin
   wins by 0.18 m; on D4–D8 the groin wins by 5.06 m.
2. **Demean, or narrow.** A uniform level offset is absorbed by the source/sink
   calibration downstream, so it is not the groin's job. Demeaning makes the
   groin win even on the full D1–D12 window (+2.17 m).

### What the groin explains, and what it doesn't

| period | observed | groin supplies | to source/sink |
|---|---|---|---|
| 1 | +52 m | **≈ +17.2 m (33%)** | ~67% |
| 2 | −76 m | ≈ 0 | ~all |

**Period 1 is the only window where the groin does something the module can
reproduce.** Including it in period 2 is right for consistency of the
structure's timeline, not because it explains that period's shoreline.

---

## 2. What the module cannot do

**Trapping is bounded at ≥ 0.** The groin can stop adding sand; it cannot
actively drain the fillet. Period 2's −76 m release is outside the
parameterisation at any (M, f), and is carried by the source/sink calibration
together with the Cape Point dynamics the dipole does not represent.

Three further limits, all measured:

- **Volume-neutral dipole is wrong at this site.** Observed downdrift extent is
  **0 m**; the model's is **2,500 m**. The real structure accretes updrift
  without a measurable downdrift deficit — the sand comes from the cape, not
  from D5. Tested: removing the sink halves the fillet and quadruples reach
  bias, so "delete the sink" is the wrong fix.
- **Sub-grid.** The real fillet is ~190 m wide; one model domain is 500 m. M is
  an **effective, grid-specific, field-aggregate** rate — not a sediment flux,
  not divisible by four for a per-structure value.
- **Only ~3% of what the dipole injects is retained** *(measured 2026-08-30,
  `output/calibration/groin/figures/sediment_budget.png`)*. Over the rig's 50 years at
  M = 60, f = 0.6 the module applies **2,400 m** of cumulative one-sided
  displacement — ±28.7 million m³ — and holds a fillet of **69 m** (peak 129 m).
  BRIE's alongshore diffusion removes the rest. So M is the rate needed to
  **sustain** a fillet against diffusion, not the rate sand is impounded, and
  the 719,000 m³/yr affordability figure below is a **gross restoring rate set
  against a net transport budget** — not like for like. It remains a fair reason
  to prefer M = 60 over M = 95; it is not a statement that the groin impounds
  more sand than the coast carries.
- **Four groins, one dipole, deliberately.** The field spans northing
  3901373–3901789, entirely inside D6, so four dipoles in one cell would sum to
  the one the cell can express.

---

## 3. Corrections worth remembering

Things that looked like problems and were not, or vice versa:

- **f → 0 in the period-2 sweeps was the right answer, not a rail artefact.**
  The observations show the fillet declining after 2004, i.e. trapping ceased.
  Considerable time was lost re-defining targets to "fix" a result that was
  correct.
- **The stability ceiling is rig-specific.** M ≥ 70 went unstable and M ≥ 100
  drowned the barrier on the 41-domain rig. On the production 120-domain grid,
  all 36 cells including M=70 and M=80 ran clean. **Do not quote that ceiling
  for production runs.**
- **The hindcast windows cannot fit this groin, and nothing is wrong with
  them.** They begin 15 years after installation, so they record the fillet's
  decay rather than its creation. Fit on the 1967 window; apply in the hindcast.
- **The continuous 1984–2024 window is the WRONG instrument for fitting, and
  "no groin fits best" was an artefact of it.** That window NETS period 1's
  build (+52 m) against period 2's collapse (−76 m), so a module that can only
  widen the gap can never win. Its 43-cell sweep is still useful — it bounds M
  from above — but it cannot fit a groin, and it should not be read as evidence
  against one. **Fit on period 1.**
- **Including D1 with a raw score hid the groin entirely.** The cape error
  (81–104 m over period 1) is five times the groin signal (~17 m). Excluding D1,
  or removing the mean, flips the result in 7 of 8 window/score combinations.
- **The model's fillet relaxes too slowly** — −15.3 m over 1984–2024 with no
  groin at all, against −24.4 m observed. Adding a groin makes it more
  persistent still.
- **Raising Hs fixes the decay but is not worth it — tested and rejected.**
  Wave height drives BRIE's alongshore diffusivity, so it was the obvious
  candidate. At Hs = 3.5 the fillet decay improves to −20.6 m (closing ~58% of
  the gap to observed), confirming the mechanism — **but the reach RMSE goes
  from 15.97 to 36.58**, a 2.3× degradation of the exact quantity the
  source/sink calibration is built on. **Keep Hs = 2.5.** Cost ~10 minutes to
  settle instead of a full re-calibration.
  (Hs = 3.0 crashed on a pre-existing CASCADE bug —
  `beach_dune_manager.py:254 filter_overwash`, interior grid 136 vs 137 rows —
  unrelated to the groin, but that configuration is known to fail.)

---

## 4. Where things live

| | |
|---|---|
| fitted values | `output/calibration/groin/joint_fit.json` |
| two-period sweeps | `output/calibration/groin/<period>_<preset>/` |
| continuous 1984–2024 sweep | `output/calibration/groin/fullperiod_1984_2024/` |
| 1967 rig | `results/sensitivity_sweep/` (runs: `runs/1967_2017_run/`) |
| GIS analysis | `../../1-observations/coastsat_shoreline/` |
| observed fillet table | `../../1-observations/wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv` |

**One operational rule:** never run two sweep orchestrators at once. Every
CASCADE construction writes the shared
`data/hatteras_init/Hatteras-CASCADE-parameters.yaml`; concurrent writers
corrupt it, which cost 25 cells of a running sweep on 2026-08-24.
