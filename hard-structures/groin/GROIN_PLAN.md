# Buxton groin — observed history, chosen parameters, and the hindcast plan

Reference note, 2026-08-24.

**Figures** (regenerate with the scripts beside them):

| file | shows | script |
|---|---|---|
| `groin_two_shorelines.png` | **start here** -- the two shorelines and the gap between them; needs no definition | `HAT_groin_two_shorelines_figure.py` |
| `groin_module_logic.png` | **the mechanics** -- what the module adds each year, and why it cannot close the gap | `HAT_groin_module_logic_figure.py` |
| `groin_timeline_and_hindcast.png` | the fillet timeline against the module's `M_eff` schedule | `HAT_groin_timeline_figure.py` |

**The one thing to keep straight:** the "fillet" is the GAP between the updrift
and downdrift shorelines, not a volume of new beach. Both sides eroded; the
sheltered side eroded less. By 2004 the downdrift side had retreated 175 m
since 1967 and the updrift side only 25 m -- that 150 m difference is the
fillet. The post-2004 collapse is ~85% the updrift side eroding once the
structure failed, not impounded sand draining downdrift.

---

## 1. What the GIS analysis shows

Fillet = shoreline offset `D5 − D6` against a fixed **1967 datum**, from 24 dated
wet/dry surveys (`HAT-groin-test-output/shoreline_position_output/Change_from_wetdry_1967_D2_D12.csv`).
Landward-positive: a rising curve means the updrift side is holding while the
downdrift side retreats — what a groin builds.

| phase | window | fillet | rate |
|---|---|---|---|
| **BUILD** | 1967 → 1978 | 0 → 117 m | **+10.7 m/yr** |
| slow growth | 1978 → 2004 | 117 → 150 m | +1.3 m/yr |
| **RELEASE** | 2004 → 2023 | 150 → 74 m | **−4.0 m/yr** |

**Structure history:** installed **1969** · last repaired **1996** · storm damage **2003**.

The fillet peaks in **2004** and declines from there — the turning point matches
the storm, not anything in the model.

### The two hindcast periods carry opposite signals

| period | window | fillet change | rate | groin is |
|---|---|---|---|---|
| **1** | 1985 → 2004 | **+52.0 m** | +2.74 m/yr | **still trapping** |
| **2** | 2004 → 2023 | **−76.4 m** | −4.02 m/yr | **releasing** |

---

## 2. The plan

**Turn the groin on for both periods with the schedule it already has.**
`GroinCallback` carries an absolute calendar timeline — install 1969,
deterioration onset 1996, linear ramp to 2003, then hold at `M × f`. That
already reproduces the measured history. **No period-specific configuration.**

| | value | where it comes from |
|---|---|---|
| `trapping_rate_m_yr` | **M = 60** | **fitted on PERIOD 1, window D4–D8, production geometry, be1 pinned at the production value, demeaned score.** RMSE 11.69 vs 15.58 with no groin — the groin closes **25%** of the shape misfit, and this is within 0.10 m of the global best. Corroborated by the 1967 rig (M=60). Intercepts ~719,000 m³/yr, marginally above the 5–7×10⁵ drift band (a literature range, not a hard limit) |
| `deterioration_fraction` | **f = 0.6** | same fit; f is only weakly constrained by period 1 (which mostly precedes the 1996–2003 ramp), so it leans on the rig's f≈0.5 and on period 2 showing trapping ceased |
| `install_year` | 1969 | documented |
| `deterioration_delay_years` | 27 (→ 1996) | last repair |
| `deterioration_ramp_years` | 7 (→ 2003) | storm damage |
| `updrift / downdrift` | GIS 6 / 5 | field occupies D6 |

Written to `output/groin_sweep/joint_fit.json`, which stage 6 and
`HAT_be_zone_LOESS_analysis.py` both read.

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

## 3. What the module cannot do

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
- **Four groins, one dipole, deliberately.** The field spans northing
  3901373–3901789, entirely inside D6, so four dipoles in one cell would sum to
  the one the cell can express.

---

## 4. Corrections worth remembering

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

## 5. Where things live

| | |
|---|---|
| fitted values | `output/groin_sweep/joint_fit.json` |
| two-period sweeps | `output/groin_sweep/<period>_<preset>/` |
| continuous 1984–2024 sweep | `output/groin_sweep/fullperiod_1984_2024/` |
| 1967 rig | `hard_structures/groin/HAT-hindcast-groin-test/sensitivity_sweep/` |
| GIS analysis | `hard_structures/groin/HAT-groin-gis-analysis/` |
| observed fillet table | `HAT-groin-test-output/shoreline_position_output/Change_from_wetdry_1967_D2_D12.csv` |

**One operational rule:** never run two sweep orchestrators at once. Every
CASCADE construction writes the shared
`data/hatteras_init/Hatteras-CASCADE-parameters.yaml`; concurrent writers
corrupt it, which cost 25 cells of a running sweep on 2026-08-24.
