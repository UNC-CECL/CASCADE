<!-- TRACKED MIRROR. The live copy is
     output/groin_sweep/SELECTED_M60_f0.60/README.md, beside the figures it
     describes -- but .gitignore excludes output/*, so that copy does not
     survive a clone. This mirror exists for the same reason
     CALIBRATION_FIGURES.md does: the reasoning has to outlive the PNGs.
     Edit the live copy, then re-copy it here. -->

# M = 60, f = 0.6 — the selected pair, and the case for it

Decided 2026-08-30. Single value, conditions stated, no ensemble over M, no
further runs. This folder holds the evidence in the order it should be read.

```
trapping_rate_m_yr        M = 60      m/yr
deterioration_fraction    f = 0.6     fraction of M held after 2003
install_year              1969        documented
deterioration onset       1996        last repair
deterioration ramp ends   2003        storm damage
updrift / downdrift       GIS 6 / 5   the four-groin field sits inside D6
```

`joint_fit.json` in this folder is a copy of the pinned file the pipeline
reads. **It is pinned, not fitted** — see "the trap" at the foot of this page.

---

## The one-paragraph version

The Buxton groin's fillet is the **gap** between the updrift and downdrift
shorelines, not a volume of new beach: both sides eroded, the sheltered side
eroded less. It builds from 0 to 117 m over 1967–1978, creeps to 150 m by 2004,
then collapses to 74 m by 2023 — and ~85% of that collapse is the *updrift*
side eroding once the structure failed, not impounded sand draining downdrift.
The module can widen that gap or stop widening it; **it cannot close it**,
because trapping is bounded at ≥ 0. So only the build-and-hold part of the
record is fittable at all. **M = 60** comes from the one hindcast window where
the gap still widens (period 1, +52 m). **f = 0.6** comes from the one window
that contains the deterioration ramp (the 1967–2018 rig). Neither window can
set both.

---

## Read in this order

| # | figure | what it settles |
|---|---|---|
| 01 | `01_why_M60_f0.6.png` | **the choice.** Period-1 D4–D8 fit: no groin 15.68 m → chosen 11.65 m, a 26% improvement. Best cell is 11.54 m, so the chosen pair is 0.11 m off optimal. Draws the tied band rather than a star, because 13 cells sit within 0.5 m |
| 02 | `02_period1_fit_profiles.png` | the top cells and the no-groin baseline against the observed change profile, fit window marked |
| 03 | `03_full_life_1967_2017.png` | **the module against every survey over the structure's whole life**, with the applied trapping-rate schedule underneath. The only window containing the build phase — and the strongest single piece of evidence here. See the note below |
| 04 | `04_fillet_vs_surveys_hindcast.png` | the same comparison inside the two hindcast windows, both presets — and the residual each leaves for the source/sink calibration |
| 05 | `05_rig_f_bracket.png` | **why f = 0.6**: the rig's RMSE bracketed on both sides, and the stability wall that limits what the rig can say about M |
| 06 | `06_three_targets.png` | the same 61 runs scored three ways give M = 95, M = 0 and M = 60. **The fitted groin is set by the target, not by the data** |
| 07 | `07_M_f_identifiability.png` | the (M, f) surface with iso-M·f contours — refutes "only the product is identified" |
| 08 | `08_window_sensitivity_D4D7.png` | D4–D7 gives the *largest* RMSE gain of any window, and why a larger gain is not a better fit |
| 09 | `09_profiles_by_M.png` | every M's end-state profile against observed, in one panel |
| 10 | `10_sediment_budget.png` | **what the module moves, and what it keeps** — annual and cumulative volume against the littoral drift band, the volume-neutrality assumption, and the ~3% retention that says M is not a rate of impoundment. Best read straight after 05 |

### What figure 03 actually shows

Run 2026-08-30 on the corrected topography at the decided pair:

| | observed | model |
|---|---|---|
| fillet peak | 155 m (1995) | **129 m** |
| end 2017 | +112 m | +67 m |
| groin OFF, whole run | — | **flat at ~0 m** |

Three things to take from it:

1. **The groin is doing essentially all of the fillet work.** The groin-OFF
   baseline never leaves zero. Whatever the module's limits, the fillet in this
   model is the module, not a background gradient.
2. **The build phase is reproduced in sign, timing and roughly in amplitude** —
   129 m against 155 m observed, from a standing start at install. This is the
   part of the record neither hindcast window contains.
3. **The model over-declines after 2003.** Observed recovers to ~112 m by 2017;
   the model runs down to 67 m. So `f = 0.6` is, if anything, slightly too
   aggressive a deterioration over the long tail — the opposite of the concern
   that f is a free knob tuned to make the groin look good.

Compare this against `04_fillet_vs_surveys_hindcast.png`, where the same pair
*undershoots* inside period 1 (+14 m modelled against +52 m observed). Both are
true: the module builds a credible fillet from scratch over 50 years, and
cannot reproduce the *rate of widening* inside a 20-year window that starts
15 years in with the fillet already inherited in its initial shoreline.

### The gifs

**`gifs/00_full_life_vs_observations_1967_2017.gif` is the one to watch.**
51 frames, 1967→2017, D2–D12, shoreline change since 1967 against the survey
record — and here **the observations animate too**. The hindcast gifs can only
draw a fixed 1984→2004 endpoint, because inside that window the observation IS
a single endpoint. Over the full life the wet/dry record carries 19 dated
surveys, so each frame shows the most recent survey at or before that model
year, with the last five left behind as fading ghosts. Phase banner, applied
M_eff and the running D6−D5 fillet (model vs survey) are on every frame.

**All seven gifs are landward-positive — erosion moves UP**, so each panel
reads as a plan view of the island with the ocean below the axis and the island
above it. Gifs 01–06 were seaward-positive until 2026-08-30 and were flipped to
match; the fillet is reported throughout as **D5 − D6**, matching
`GROIN_PLAN.md`.

**The static profile figures were flipped to match on the same day** — 01, 02,
08, 09 and `fig_top_profiles.png`. Figure 06 needed no flip: it plots error
against M, not a profile, so it has no direction. Figures 03, 04, 05 and 10 are
time series or score curves and are likewise unaffected.

The flip was not a single sign change — the two data sources have opposite
conventions. `shoreline_matrix.npy` (Barrier3D's `x_s`) and the 1967 wet/dry
table are both **landward**-positive and are now used as they come;
`observed_change_profile()` reads CoastSat chainage, which is **seaward**-positive
and is negated. So the minus sign moved from the model series to the observed
one.

**Every flip is at the plotting layer only.** `HAT_groin_sweep_config.FLIP_SIGN_MODEL`,
`HAT_fullperiod_target` and the sweep worker are untouched, so the scoring
pipeline, the joint fit and every reported bias keep their original convention.
Because both series flip together, the squared difference is invariant and
**every RMSE in these figures is numerically unchanged** — verified against the
documented values: 11.44 best / 11.58 chosen on D4–D8, and 10.09 at M = 60.

Not demeaned, unlike the other gifs: the question here is how well the module
reproduces the *measured position*, so the level is left in. Read a whole-curve
offset as the source/sink term's business and the shape around D5/D6 as the
groin's. Survey coverage is labelled per frame — several surveys cover only
8–10 of the 11 domains and the 2017 one covers 6, so those lines stop early
rather than interpolating across a gap.

A timeline strip beneath the panel carries the structure's four phases (before
the groin / sound / deteriorating / failed) with the installation, last repair
and storm-damage years marked, and a cursor tracking the frame — so the state
of the structure is legible without reading the M_eff number.

What it shows: the D5 trough and the D6 step appear together as soon as the
groin switches on in 1970 and are absent from the groin-OFF curve entirely;
through the 1970s–90s the model tracks the survey shape closely (1978: fillet
model +113 m vs observed +117 m; 1997: +122 vs +136); after the deterioration
ramp the model keeps eroding D6–D12 past the observations and ends 2017 at a
+67 m fillet against +112 m measured.

`gifs/00b_full_life_planform_1967_2017.gif` is the same run's own D1–D15
planform evolution. `gifs/01_CHOSEN_M60_f0.60.gif` is the hindcast-window run. The other five are the
comparison that makes the choice legible — each one is a thing that goes wrong:

| gif | what it shows going wrong |
|---|---|
| `02_compare_M40_f0.60` | too weak — the fillet barely forms |
| `03_compare_M95_f0.60` | too strong — and ~1.1M m³/yr, roughly 2× the littoral drift |
| `04_compare_M160_railed_joint_fit` | what the **joint fit's own ranking** would have run. This is why that ranking is not used |
| `05_compare_f0.00` | no trapping after 2003 at all |
| `06_compare_f1.00` | the structure never deteriorates — contradicts the 2003 damage outright, and it is the top-scoring cell on period 1 |

The last two are the point of f: period 1 alone would happily pick f = 1.0.

---

## Why M = 60

**From period 1, D4–D8, demeaned, be1 = −42.6.** `hatteras_site_config.py`
records the canonical scoring as 15.20 m no-groin → 11.58 m chosen against a
best of 11.44 m; the figure in this folder recomputes marginally different
absolutes (15.68 → 11.65, best 11.54) because it scores off a different
intermediate. **The ranking is the same and the conclusion is the same:** the
production pair is statistically indistinguishable from optimal.

Three things had to be right for the groin to be visible at all:

1. **Exclude D1.** Cape Point's shoreline change over period 1 is 81–104 m,
   about **five times** the groin's ~17 m signal. On the full D1–D12 window
   with a raw score, no-groin wins by 0.18 m; on D4–D8 the groin wins by 5.06 m.
2. **Demean.** A uniform level offset is absorbed by the source/sink
   calibration downstream, so removing it is not the groin's job. Demeaning
   keeps every gradient, unlike a linear detrend, which would partly absorb the
   dipole's own gradient and hide a working groin.
3. **Fit period 1, not the full window.** See the root README's trap 1.

**Corroboration:** D3–D9, the one other window symmetric about the structure,
independently returns M = 60 at the same gain (2.05 vs 2.10).

**Affordability** prefers 60 over 70: M = 60 intercepts ~719,000 m³/yr against
a 5–7 × 10⁵ m³/yr littoral drift, M = 70 about 838,000 (~1.3×). This is a
**literature range, not a hard limit** — a reason to prefer, not a prohibition.

---

## Why f = 0.6

**The hindcast windows cannot see f and should not be asked.** They begin 15
years after installation and contain almost none of the 1996–2003 ramp. Period-1
cumulative trapping is `M(15.5 + 4.5f)` — f moves it by only 29% across its
entire range, so period 1 reads high f as simply "more trapping" and rails at
f = 1.0. Period 2 is `20·M·f` and prefers f = 0.

**The 1967–2018 rig is the only window that spans the ramp.** At M = 60:

```
f      0.1     0.3     0.4     0.5     0.6     0.7     0.9
RMSE  33.17   27.99   25.84   24.21   23.78   25.06   39.01
                                       ^ bracketed on both sides
```

That is a clean interior minimum with steep curvature either side.

**f is not a free knob.** It encodes a maintenance record — installed 1969,
last repaired 1996, storm damage 2003, fillet peaks 2004. Making the module a
static trapping rate was tested and rejected on measurement: at M = 60, setting
f = 1 degrades period 2 from 14.90 to 17.87 m (+20%), and moves the modelled
D5−D6 differential the *wrong way*.

---

## What this pair does NOT claim

**It is not uniquely determined.** 13 cells lie within 0.5 m of the best,
spanning M = 40–95 and f = 0.4–1.0.

**The rig does not confirm M — corrected 2026-08-30.** An earlier note said the
rig bracketed both parameters. It does not. Its M improves *monotonically* to
the last value that runs:

```
M      20      40      60      70          80          >=100
RMSE  53.9    33.3    23.8    320-378     556-766     (crashes)
```

The jump at M = 70 is a 13× discontinuity — the instability `GROIN_PLAN.md`
records, not a fit degradation. **M = 60 is the largest M the rig can hold, not
its optimum.** The rig corroborates f; it is only *consistent with* M. (The
stability wall is rig-specific: all 36 production cells including M = 70 and
M = 80 ran clean. Do not quote that ceiling for production runs.)

"Independent" is also generous. Since the 2026-08-30 repoint both routes use
the same topography product (`1984-start/v1`), the same wave climate, the same
BRIE physics and the same wet/dry shoreline family. They are independent
**windows**, not independent **evidence**.

**The window does real work.** Across defensible windows M ranges 40–95:

```
window    best M    f     RMSE   no-groin   gain
D4-D8         60   1.0    10.09     12.19    2.10   <- production
D3-D9         60   1.0    13.93     15.98    2.05   <- supports it
D4-D7         70   1.0    10.54     13.63    3.09
D3-D8         95   0.8    12.91     16.80    3.89
D4-D9         40   1.0    12.23     12.65    0.41
D5-D7          0   0.0     6.18      6.18    0.00   <- see below
D5-D8          0   0.0     6.18      6.18    0.00
```

**And the sharpest result: the domains closest to the groin are the ones the
groin explains least.** D5–D7 and D5–D8, the windows tightest on the structure,
return M = 0 at a gain of exactly 0.00. The 2.10 m gain at D4–D8 comes from D4
and D8, the window's *outer* edges. D6–D7 carries an erosion trough peaking one
domain **north** of the structure that a groin actively worsens by pushing D6
seaward. This should be said plainly rather than buried.

**It is conditional on be1, not independent of it.** The D4–D8 optimum trades
off monotonically with the edge rate: M = 50 at be1 = −46, 60 at −42.6, 80 at
−34, 125 at −10. What makes this acceptable is that the surface is smooth and
monotonic, and the global minimum over the whole (be1, M, f) space sits at
be1 = −42.6 — the independently calibrated value. The groin fit prefers the
same edge rate the source/sink fit arrived at.

**It is void without edgeBE.** Under zeroBE the fit improves monotonically to
the grid edge, because with no edge forcing the groin simply absorbs the
missing source/sink term. A zeroBE groin fit is meaningless.

**It does not reproduce the shape.** Inside the fit window the observed profile
peaks at D6, dips at D7 and peaks again at D8; every cell draws a smooth
monotonic rise. The gain comes from matching the overall D4–D8 slope. Read the
residual as the split between what the groin explains and what the source/sink
calibration absorbs.

**What the groin explains:** period 1, +17.2 m of an observed +52 m (33%).
Period 2, ≈ 0 of an observed −76 m. It runs in period 2 for consistency of the
structure's calendar timeline, not because it explains that period.

---

## The trap

`joint_fit.json` — here and in the parent directory — is **pinned by hand**.
Its own ranking returned **M = 160, f = 0.8**, railed at a grid bound, because
it scores both periods jointly and period 2 records a release the module cannot
produce at any (M, f). `HAT_run_all.py` stage 6 reads that file and passes
whatever it holds to every groin run in the matrix. On 2026-08-30 it was found
holding the railed values.

**Re-running stage 5 overwrites it. Re-pin afterwards.** The superseded ranking
is kept at `../archive/joint_fit_RAILED_ranking_20260830.json` and each entry
in the pinned file carries its own `superseded_ranking` block.

Two other places the pair is written, both corrected 2026-08-30:
`scripts/hatteras_ms/hat_run.yaml` (was M = 50 / f = 0.9 placeholders) and
`scripts/hatteras_ms/HAT_hindcast_config.py` (the same values as fallback
defaults). A matrix run reads `joint_fit.json`; a single run reads the yaml.
**Check both before quoting a groin run's parameters.**

---

## What the sediment budget shows — figure 10

`groin_diagnostics.csv` has logged `cumulative_updrift_m` and
`cumulative_downdrift_m` every year of every groin run since the module was
written, and nothing plotted them until 2026-08-30. Two things follow.

**The transfer is exactly volume-neutral, and the real structure is not.** Over
the rig's 50 years the module moves **±28.7 million m³** — 718,500 m³/yr while
the structure is sound, falling to 431,100 m³/yr after 2003 — and the net across
the pair is exactly zero by construction. The real Buxton structure accretes
updrift with **no measurable downdrift deficit**: observed downdrift extent is
0 m against the model's 2,500 m, because the sand comes from the cape, not from
D5. This is the module's largest structural assumption. It has been tested —
removing the sink halves the fillet and quadruples reach bias — so deleting it
is not the fix.

**Almost none of what the dipole injects is retained, and this changes how the
affordability number should be read.** The module applies **2,400 m** of
cumulative one-sided displacement over the run and holds a fillet of **69 m** at
the end (peak 129 m). That is **~3% retained**; BRIE's alongshore diffusion
removes the rest.

So M is **not the rate at which sand is impounded** — it is the rate needed to
*sustain* a fillet against diffusion, which is a much larger number. The
719,000 m³/yr affordability figure quoted in `GROIN_PLAN.md` and in the run
reports is therefore a **gross restoring rate set against a net transport
budget**, and the two are not like for like. That does not make the comparison
wrong; it is a deliberate, documented diagnostic, and it remains a reasonable
argument for preferring M = 60 over M = 95. But **"marginally above the drift
band" must not be read as "impounds more sand than the coast carries."**

This is the quantitative version of `GROIN_PLAN.md`'s existing warning that M is
an effective, grid-specific, field-aggregate rate and "not a sediment flux."

---

## Gaps worth filling — figures that do not exist yet

Ranked by what they would actually settle, not by how easy they are. Nothing
here is required for the decision; the pair is locked. These are for making the
module's behaviour legible to someone who was not in the room.

**1. ~~The sediment budget over time.~~ BUILT 2026-08-30** — figure 10, from
`HAT_groin_sediment_budget_figure.py`. It turned out to carry more than
expected: see the section above. The ~3% retention result was not anticipated
and is worth propagating into how M is described elsewhere.

**2. The residual decomposition.** Period 1 observes +52 m; the groin supplies
≈ +17.2 m (33%) and the source/sink calibration absorbs the rest. A stacked bar
or waterfall per period — observed, groin's share, source/sink's share, Cape
Point's unrepresented share — would state the division of labour that the whole
calibration rests on. Right now it exists only as numbers in prose.

**3. Updrift and downdrift plotted separately.** Every figure here plots the
fillet, which is a *difference*. But `GROIN_PLAN.md`'s central correction is
that both sides eroded and the sheltered side eroded less — and that ~85% of
the post-2004 collapse is the *updrift* side failing, not sand draining
downdrift. A two-panel D5 and D6 absolute-position plot would make that
visible instead of asserted, and would show whether the model gets the right
answer for the right reason or cancels two errors in the difference.

**4. A be1 × M surface.** The pair is conditional on be1 = −42.6, and the D4–D8
optimum slides 50 → 125 as be1 goes −46 → −10. The cells exist in
`1984_2004_edgeBE/sweep_results.jsonl`. One contour plot with the production
be1 marked, and the global minimum marked, would show the reader the single
largest caveat on M — and the modest cross-check that the groin fit prefers the
same edge rate the source/sink fit independently arrived at.

**5. The full window-sensitivity panel.** `08_window_sensitivity_D4D7.png`
covers one alternative window. The table in this README covers seven. A small
multiple — one axis per window, best M annotated — would show honestly that
M ranges 40–95 and that **D5–D7 and D5–D8, the windows tightest on the
structure, return M = 0 at a gain of exactly 0.00.** That last result deserves
a figure rather than a footnote.

**6. A side-by-side gif.** The comparison gifs here are separate files, so
comparing them means flipping between windows. One animation with the chosen
run and two comparison runs in stacked panels on a shared clock would do the
job they currently do badly. *(Partly addressed 2026-08-30: the full-life gif
now animates the observations, so the chosen run at least has a moving target
to be judged against. The M-to-M comparison is still file-flipping.)*

Scripts for 1, 2 and 3 would each read data already on disk. 4 and 5 read the
existing sweep jsonl. None needs a new model run.
