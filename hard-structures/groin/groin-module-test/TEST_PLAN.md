# Groin module test: idealized rig, design framework

Draft, 2026-09-11. Four design choices were settled on 2026-09-11 and are recorded
in section 4: **sweep wave height rather than pick one**, **40 domains -- 20
working plus 10 buffer a side**, **200 years maximum**, and **sweep M / r_ipl
rather than an M × f grid**. A second round on 2026-09-11 settled the test
matrix in section 5 and the handling of the blocker in 5.1. Three choices in
section 4 are still open.

This tests the MODULE, not Buxton. Buxton is in `../GROIN_PLAN.md`, and the two
documents answer different questions: that one asks whether a groin at Buxton
explains the observed fillet, this one asks whether `cascade.groin` does what it
was designed to do, on a coastline with nothing else going on.

---

## 1. What is being tested

`GroinCallback` adds `-M` to one domain and `+M·f` to its neighbour each year,
before BRIE's alongshore solve. Everything else about a groin in CASCADE is
emergent: the fillet's size, its alongshore extent, its taper, how long it
persists after the structure fails. "Emergent" is a claim that needs testing, and
three separate layers can each be wrong:

| layer | what could be wrong | instrument |
|---|---|---|
| **bookkeeping** | signs, indices, the install and deterioration schedule, the diagnostics | 5-domain rig, short runs |
| **the solver** | BRIE's diffusion does something other than spread a dipole | offline audit, no CASCADE |
| **the coupling** | Barrier3D's cross-shore response changes the answer | wide rig, long runs |

A single rig cannot separate those, which is why section 3 proposes three.

---

## 2. Stage 0 is already done, and it changes the design

`0-solver-audit/HAT_groin_solver_audit.py` reproduces BRIE's alongshore solve
exactly -- same diffusivity table, taken from a live `Brie` instance, same matrix
assembly, same periodic wrap, same clip -- with Barrier3D removed. A 200-year,
six-axis experiment runs in about a minute. Every number below is printed by that
script, so re-run it rather than trusting this table.

**Seven findings, in order of how much they matter.**

### 2.1 The groin has a stability limit, and the default wave climate sits on it

BRIE's diffusion number is clipped at zero, and its wave-climate-averaged
diffusivity turns negative once the shoreline angle across one cell passes a
critical value. Past that angle a cell stops exchanging sand with its neighbours,
keeps receiving its share of the dipole, and the whole reach translates landward
at roughly M metres a year forever.

| wave height | diffusion number | M = 40 | M = 60 | M = 80 |
|---|---|---|---|---|
| Hs 1.0 (CASCADE default) | 0.196 | 101 m | **runaway, year 18** | **runaway, year 9** |
| Hs 2.5 (Hatteras) | 0.793 | 20 m | 30 m | 42 m |

A 5-domain rig at default settings, run far into the future at the calibrated
M = 60, **does not show you a groin**. It shows you a numerical runaway that looks
like a groin of implausible strength. After roughly 1,000 years it stops with
`IndexError: index 180 is out of bounds` from `brie.py` line 1297, because the
diffusivity lookup clips its index to 180 on a length-180 table. That takes a
57 km offset across one 500 m cell, which a runaway reaches and a groin never
does.

The critical angle is about -23 degrees, a 212 m offset across one domain. **The
observed Buxton fillet peaks at 155 m.** The site sits at roughly three quarters
of the shutdown threshold, which is close enough that it belongs in the paper.

### 2.2 M is not a free parameter, it is M / diffusion number

Every stable cell collapses onto one line:

```
fillet  =  0.42 · M / r_ipl          (32 cells, spread ±0.04, at 200 yr)
```

and `r_ipl` scales as Hs^2.4. So a fitted M is a statement about a wave climate
and nothing else. At Hs 2.5 the fitted M = 60 gives a 30 m fillet; the same 30 m
costs M = 15 at Hs 1.0 and M = 95 at Hs 3.5. The Hs = 3.5 experiment recorded in
`GROIN_PLAN.md` section 4 was therefore also an M experiment, whether or not it
was meant to be.

This is the sweep axis to use. One dimensionless number replaces two.

### 2.3 The shutdown threshold is set by the wave angles, not the wave height

The critical angle does not move at all with Hs. It moves a lot with the angular
distribution:

| `wave_angle_high_fraction` | shutdown angle | offset across one cell |
|---|---|---|
| 0.1 | -25 deg | 233 m |
| **0.2 (used)** | **-23 deg** | **212 m** |
| 0.3 | -20 deg | 182 m |
| 0.4 | -14 deg | 125 m |
| 0.5 | unstable at zero angle | the straight coast does not diffuse |

At 0.4 the threshold is **below the observed Buxton fillet**, and a groin at
M = 60 runs away in year 8. Raising wave asymmetry does the same thing more
gently: 0.5 gives a 664 m threshold, 0.9 gives 162 m.

So the module's stability depends on two wave parameters nobody has treated as
groin parameters. They belong in the sensitivity design next to M and f.

### 2.4 The solve is not volume conserving, and the error is the size of the signal

The matrix is scaled row by row by a per-domain diffusion number, so its column
sums are not unity and the scheme does not conserve sand once the shoreline stops
being straight. The cleanest demonstration is `f = 1`, where the dipole injects
exactly nothing and the reach mean must therefore not move. At 40 domains it moves
10.2 m in 40 years and 63 m in 200, all of it invented.

At the production `f = 0.6` the reach mean should move, and the invented part of
that movement is the difference between what moved and what the injection implies:

| reach | invented drift, 40 years | 200 years |
|---|---|---|
| 5 domains | 50.0 m | 252 m |
| 40 domains (chosen) | 6.2 m | 31 m |
| 121 domains (Hatteras production) | 2.0 m | 10.2 m |

At production width over a hindcast window the invented drift is about 2 m,
against a groin signal of about 17 m. That is 12%, it is a **uniform offset**
across the reach, and the production score is demeaned -- so the fit is protected,
for a reason beyond the one `GROIN_PLAN.md` records for demeaning. Worth saying
out loud, because it means demeaning is load-bearing rather than tidy.

It also reframes "only 3% of what the dipole injects is retained". Diffusion in a
closed periodic reach cannot remove sand, only spread it. The missing 97% went
into the reach mean, into Barrier3D's cross-shore response, or into this closure
error, and those are three different statements. Stage B below separates them.

### 2.5 Five domains costs 18% of the fillet and inflates the drift 24-fold

The solve is periodic, so a short reach wraps its own fillet into its own
downdrift notch. At M = 60, f = 0.6, Hs 2.5, 200 years:

| domains | reach | fillet | vs 121 domains | spurious mean drift |
|---|---|---|---|---|
| 5 | 2.5 km | 24.9 m | **-18.3%** | -3.54 m/yr |
| 11 | 5.5 km | 28.2 m | -7.5% | -1.61 m/yr |
| 21 | 10.5 km | 29.6 m | -3.2% | -0.85 m/yr |
| **40 (chosen)** | **20.0 km** | **30.2 m** | **-0.9%** | **-0.45 m/yr** |
| 81 | 40.5 km | 30.5 m | 0.0% | -0.22 m/yr |
| 121 | 60.5 km | 30.5 m | 0% | -0.15 m/yr |

The fillet converges fast, which is good news for a small rig. The reach mean does
not, because the net source from `f < 1` is spread over fewer cells. **Five
domains is fine for testing mechanics and wrong for anything involving the mean
shoreline, the barrier's survival, or the alongshore extent.** This table, and not
any Hatteras precedent, is the reason the morphology rig is 40 domains wide.

### 2.6 The buffers separate, they do not absorb

BRIE's solve is **periodic**, so a buffer domain is not an open boundary and does
not soak anything up. Its only job is to keep the groin away from the wrap. The
dipole's diffusive reach grows as `dy·sqrt(2·r_ipl·t)` and contains no M at all,
so whether a 10-domain buffer is wide enough is decided by the wave climate and
the run length, never by the structure. At 200 years:

| Hs | diffusion number | reach at 200 yr | inside a 10-domain buffer? |
|---|---|---|---|
| 1.0 | 0.196 | 8.8 domains | yes |
| 1.5 | 0.368 | 12.1 | no |
| 2.0 | 0.569 | 15.1 | no |
| 2.5 | 0.793 | 17.8 | no |
| 3.0 | 1.038 | 20.4 | no |

So above Hs 1.0 the fillet's tail comes round the back of a 40-domain rig before
the run ends. Section 2.5 says that costs under 1% of the fillet, so the **fillet**
measurements stand. The **reach mean** and any attempt to measure the fillet's
**alongshore extent** do not: test C2 needs either a wider grid or an earlier
readout. Worth knowing before the extent result gets reported as a failure of the
module.

### 2.7 A groin field behaves like a groin field, emergently

Dipoles every `spacing` domains, 200 years, with the fillet at the most updrift
structure against the one in the middle:

| structures | spacing | apart | updrift fillet | interior fillet | interior / updrift |
|---|---|---|---|---|---|
| 5 | 1 | 0.5 km | 61.7 m | 28.9 m | 0.47 |
| 5 | 2 | 1.0 km | 50.2 m | 29.8 m | 0.59 |
| 5 | 8 | 4.0 km | 41.1 m | 29.0 m | 0.71 |
| 5 | 16 | 8.0 km | 35.2 m | 29.4 m | 0.83 |

The most updrift structure holds the largest fillet and the interior ones are
starved -- the textbook behaviour of a real field, and nobody put it in. The
structures stop interacting at about 8 km apart, which is the diffusion length
over the run. That is a genuine emergent result and a good figure.

Two cautions on the same table. Co-located dipoles sum **exactly** -- four
structures at M = 60 in one cell are bit-for-bit one structure at M = 240 --
which is not a test of anything, it is a demonstration that a field inside one
cell is **unidentifiable** from a single stronger groin. And M = 240 runs away at
Hs 2.5, so "M = 60 is a field aggregate, not divisible by four" is not a
convenience, it is required.

A mirrored pair facing each other gives 37.2 m on one side and 32.5 m on the
other. The geometry is symmetric and the answer is not, by 14%, because BRIE takes
its shoreline angle as a forward difference while its Laplacian is centred. A
groin's strength in this model depends slightly on which way it faces.

One limit the 40-domain rig imposes on this: five structures 16 domains apart span
64 domains and do not fit. Spacings up to 4 domains fit inside the 20 working
cells with three structures, and the wider spacings need the audit's 121-domain
grid or a wider CASCADE rig. The spacing study therefore lives mostly in tier 0.

---

## 3. The rig, in three tiers

Tier 0 exists and is cheap. Tiers 1 and 2 are CASCADE runs and need your
decisions first.

**Tier 0 -- solver audit.** `0-solver-audit/`. No Barrier3D, no storms, no sea
level. Seconds per experiment. Answers: what does the diffusion do with a dipole.

**Tier 1 -- unit rig, 5 domains.** `1-unit-rig/`. Default Barrier3D domains, one
groin, short runs, every diagnostic written out. Answers: does the module keep its
own books. Deliberately small and fast enough to run on every change to
`cascade/groin.py`. Its fillet is 18% low and its reach mean is meaningless, and
that is acceptable for what it tests.

**Tier 2 -- morphology rig, 40 domains: 20 working, 10 buffer a side.**
`2-morphology-rig/`. The fillet is within 0.9% of converged and the spurious drift
is 0.45 m/yr. This is where the long runs, the sensitivity sweeps and the groin
fields go. The buffers separate the groin from the periodic wrap, they do not
absorb anything, and above Hs 1.0 the dipole's tail exceeds them within the run --
see 2.6, which decides which measurements that does and does not spoil.

Both CASCADE tiers get **their own copy of the Barrier3D parameter YAML inside
their own directory**. The shared `data/default_cascade_variables/barrier3d-default-parameters.yaml`
has already been overwritten by Hatteras runs -- it currently points at a Hatteras
storm file and carries `TMAX: 19` -- and `set_yaml` rewrites it on every CASCADE
construction. Per-rig copies are what make the sweeps safe to run in parallel,
which lifts the one-orchestrator-at-a-time rule for this tree.

---

## 4. Assumptions

Four decided on 2026-09-11, three still open. All of them change what the model
ingests, which is why they are here and not in the code.

### Decided

**4.1 Wave height is a SWEEP AXIS, not a setting.** Hs 1.0 to 3.0. This is the
better answer than picking Hatteras or the default, because section 2.2 shows the
fillet is bought by `M / r_ipl` and `r_ipl` scales as Hs^2.4 -- so sweeping wave
height *is* sweeping the cost of M, and the runaway boundary gets mapped instead of
assumed. The rig is a module test, not a Hatteras test, so no single site's wave
climate is privileged. Hatteras values stay in the table as one row among five.

**4.2 Two rigs: 5 domains for mechanics, 40 for morphology.** 40 = 20 working plus
10 buffer a side. Justified by the convergence table in 2.5 -- at 200 years a
40-domain reach is within 0.9% of a 121-domain one -- and not by any Hatteras
precedent. The buffers separate rather than absorb, per 2.6.

**4.3 200 years is the ceiling.** Per advisor. This costs less than it sounds: at
Hs 2.5 the fillet is within 2% of its steady state by year 40, and every runaway in
the audit declares itself between years 4 and 18. What 200 years cannot show is the
slow post-failure relaxation, so test D4 reports a fitted e-folding time rather
than a completed decay. At 0.004 m/yr the barrier takes 0.8 m of sea level over the
window and survives it, which a 1,000-year run would not have.

**4.4 Primary axis is `M / r_ipl`, with f and the wave-angle parameters
orthogonal.** Replaces the 6 × 6 M × f factorial, which spent 36 runs on a slice of
a three-dimensional space and never touched the axis that decides stability.

**4.8 The progradation ceiling is a RESULT, not a bug to be patched.** Barrier3D is
left alone. Stage P measures where the ceiling falls, the wider barrier is used to
push it up, and the forced stages run below it. The alternative, fixing
`route_overwash`, means editing a shared dependency and proving every existing run
is unchanged, and it would hide a limit that is genuinely part of what this module
can do.

**4.9 Two barriers, one per tier.** `b3d_pt45_802yrs_high` for the morphology rig,
at 410 m of interior, which is the widest the repo holds and therefore the highest
ceiling available. `b3d_pt75_3284yrs_low` for the unit rig, so stages A to C
exercise the same initial condition as `tests/test_coupling.py`.

**4.10 The Hs 3.0 attribution gets re-tested**, one run, inside stage P.

### Still open

**4.5 Does "perfectly straight" mean exactly straight?** BRIE seeds its initial
shoreline with a sub-metre random jitter from a fixed seed (1973), so the rig's
coast is reproducible but not flat. Recommend **zeroing it explicitly** in the rig
so the angle-dependent diffusivity starts uniform and every metre of structure in
the answer is attributable to the groin. That needs a two-line override after
construction and a note that Barrier3D's per-domain shoreface length was set from
the jittered values, leaving a sub-metre inconsistency.

**4.6 Do human modules stay off?** Recommend **off for tiers 0 and 1, and a
dedicated arm in tier 2.** A groin that delays a road relocation is the actual
management question and nobody has run it, but it has no business in a test of
whether the dipole's signs are right.

**4.7 Does the test target the sandbox `Cascade` only?** The hook lives in
`cascade/cascade_groin.py`; `cascade/cascade.py` does not have it, and the two
files differ by that hook plus cosmetic renames. Recommend **adding a tier-1 test
that the two produce bit-identical output with no groin attached**, which is
asserted in a code comment and has never been checked.

---

## 5. The test matrix

Settled in the interview of 2026-09-11. The deliverable is a **module behaviour
section in its own right**, three of the four candidate experiments get compute,
forcing is **off for mechanics and on for science with a ten-member ensemble for
island response**, and twins are paired for a **selected subset** only.

| stage | rig | forcing | runs | domain-years | status |
|---|---|---|---|---|---|
| **A** bookkeeping | 5 domains | off | 11 | 1,400 | ready |
| **B** budget closure | 5 and 40 | off | 4 | 18,000 | ready |
| **C** prediction check | 40 | off | 1 | 8,000 | ready |
| **P** progradation ceiling | 40 and 5 | on | 10 | 8,000 | ready, and it gates D to G |
| **D** release | 40 | on | 3 | 24,000 | **gated by P** |
| **E** sensitivity | 40 | on | 36 | 288,000 | **gated by P** |
| **F** island response | 40 | on | 20 | 160,000 | **gated by P** |
| **G** groin fields | 40 | on | 6 | 48,000 | **gated by P** |
| | | | **91** | **555,000** | |

**Measured cost.** A 5-domain storm-forced run times at 0.145 seconds per
domain-year on one core, so the whole matrix is about 22 core-hours. On eight
cores that is an evening, not the multi-night exposure that has killed sweeps
here before. `num_cores` parallelises across domains inside a run, so a 40-domain
run uses them all.

**Declined, and why it is recorded anyway.** The fourth candidate experiment, a
runaway bracket in full CASCADE, was not funded. The sensitivity grid in stage E
already contains cells past the shutdown boundary at the low wave height, and
those are reported as shutdown rather than trimmed, so the boundary's location in
CASCADE arrives as a by-product. The management arm stays out until decision 4.6
is made, and is written down at the end of this section so it is not forgotten.

Each line below says what would falsify it, because a test that cannot fail is a
demonstration.

### 5.1 A blocker found while validating the recipe

Barrier3D's overwash router crashes when the groin progrades a domain, and the
three forced stages cannot run at realistic M until that is settled.

**What happens.** `route_overwash` in `barrier3d.py` reads out of bounds. It is
numba-jitted with bounds checking off, so the failure arrives as a **silent
segmentation fault**, not an exception. Set `NUMBA_BOUNDSCHECK=1` and it becomes
an ordinary `IndexError`, which is how it was localized.

**What triggers it.** Measured on the 5-domain default rig at Hs 1.0:

| configuration | outcome |
|---|---|
| no groin | runs clean |
| M = 5 | runs clean through year 8 |
| M = 8, M = 10, M = 20 | crashes entering year 3 |
| M = 20 with the sink removed | crashes entering year 3 |
| M = 20 at Hs 2.5 | crashes entering year 3 |
| M = 60, storms suppressed | runs clean through year 7 |

So it needs the groin and it needs a storm, and it is not the downdrift sink:
a pure-source groin crashes on its own. The domain that fails is the **prograding**
one, which is the only domain whose interior row count has *grown*. The other
domains are unremarkable at the moment of failure.

**Why it matters beyond this rig.** `GROIN_PLAN.md` records an unexplained
failure in the same family: Hs = 3.0 crashed at `beach_dune_manager.py:254
filter_overwash` on an interior grid of 136 rows against 137. That is a row-count
mismatch on a changing interior, which is what this is. The Hatteras note
attributes it to the wave height. It may instead be the groin's progradation,
which would mean the Hs = 3.0 exclusion was misattributed.

**A second ceiling, this one physical.** Storm-free at M = 60 the run survives,
but the downdrift domain narrows from 26 interior rows to 15 in seven years and
its dunes fall to 9 cm. The default pathways barrier is about 260 m of interior,
and a 60 m/yr sink erases the downdrift cell inside a decade. **M does not
transfer between barriers**, and the rig's initial barrier sets the usable range
before Barrier3D is even involved.

**What is not blocked.** Stages A, B and C are specified storm-free, and
`route_overwash` is only reached when a storm occurs. They can run now, at any M.

**How this was settled (decision 4.8).** Barrier3D is not patched. The ceiling
becomes stage P, a measured result, and the forced stages run beneath it on the
widest barrier available. A silent segmentation fault is a poor way for a model to
report a limit, so stage P runs with `NUMBA_BOUNDSCHECK=1` and records the year and
the state at every failure.

**Which barrier hosts the groin.** The interior width sets the physical ceiling, so
the choice is part of the experiment rather than a detail. What the repo has, in
interior rows above water:

| initial barrier | rows | interior | dunes |
|---|---|---|---|
| `b3d_pt75_3284yrs_low` | 26 | 260 m | 0.40 to 0.60 m, uniform |
| `b3d_pt75_829yrs_high` | 31 | 310 m | 0.07 m, uniform |
| `b3d_pt45_8757yrs_low` | 39 | 390 m | 0.40 to 0.60 m |
| `b3d_pt45_802yrs_high` | 41 | 410 m | 0.07 to 0.43 m |

The first is what `tests/test_coupling.py` uses and what the measurements above
were taken on. The last is 58% wider and should carry a proportionately higher M
before the downdrift cell is erased. Decision 4.9 takes the last for the morphology
rig and the first for the unit rig.

**A second, smaller blocker.** `sandbag_management_on` defaults to `False`, is
stored unwrapped at `cascade_groin.py:312`, and is indexed per domain inside
`update()`. Any caller that leaves the default raises `TypeError: 'bool' object is
not subscriptable`. `cascade.py` has it at the same two lines, so it is not a
sandbox-only defect; the Hatteras runner only escapes it by passing a list. The rig
passes `[False] * ny` and touches no library code.

### A. Bookkeeping -- 5 domains, 20 years, storms off

| test | fails if |
|---|---|
| **A1 sign** | the updrift domain does not move seaward and the downdrift domain landward |
| **A2 schedule** | the applied rate is nonzero before `install_year`, or does not hold at `M·f` after the ramp |
| **A3 ledger** | `diagnostics_frame` cumulative totals disagree with the year-by-year sum of what reached `x_s_dt` |
| **A4 inert hook** | a run with no callback differs from `cascade.cascade.Cascade` in any digit |
| **A5 determinism** | two identical runs differ. BRIE's wave-angle generator is **unseeded**, and its output feeds only the inlet model, which is off -- so runs should be reproducible. This test is what makes that claim safe to rely on |
| **A6 position invariance** | moving the groin to a different domain pair changes the fillet by more than the closure error |
| **A7 mirror** | swapping updrift and downdrift gives something other than a mirrored fillet. Expect a ~14% asymmetry from the forward-difference angle, per 2.7. Quantify it, do not fix it |
| **A8 guard rails** | non-adjacent domains, out-of-range indices, `f` outside [0,1], or `instant` with a nonzero ramp fail to raise |

### B. Budget closure -- 5 and 40 domains, 200 years, storms and sea level off

The one thing this rig exists for. Per year, account for every metre:

```
injected by the dipole  =  held in the fillet
                        +  spread into the reach mean
                        +  taken by Barrier3D cross-shore
                        +  closure error of the solve
```

Tier 0 gives the last two terms with Barrier3D absent. Running the same
configuration with Barrier3D present and differencing isolates the cross-shore
term. **B1** fails if the four terms do not sum to the injection. **B2** settles
"only 3% is retained" by naming which term holds the other 97%.

### C. Analytical prediction -- 40 domains

`predict_fillet` is documented as the module's independent scientific test:
amplitude from M, extent from diffusion and time alone. Nobody has checked it.

| test | fails if |
|---|---|
| **C1 amplitude** | the measured fillet is not linear in M below the runaway boundary. Tier 0 says it is to about 10%, with a mild upward curvature: at Hs 2.5 the fillet per unit M is 0.50 up to M = 60 and 0.55 at M = 120 |
| **C2 extent** | the measured alongshore extent departs from `dy·sqrt(2·r·t)`, which contains no M. **Read this at Hs 1.0, or early in the run, or on a wider grid**: per 2.6 the reach exceeds a 10-domain buffer by year 63 at Hs 2.5, and a wrapped tail has no extent to measure |
| **C3 the factor of two** | `predict_fillet` returns 18.9 m where the fillet is 30.2 m. The docstring does not say whether the return is the one-sided offset or the gap, and those differ by two. **This is a documentation bug to fix either way** |

### P. The progradation ceiling -- forced, short runs. 10 runs

Runs first, because it sets the M range every forced stage is allowed. Short runs,
20 years or until failure, with `NUMBA_BOUNDSCHECK=1` so a failure is an exception
with a year attached rather than a segmentation fault.

| test | what it measures |
|---|---|
| **P1 ceiling on the wide barrier** | M at 5, 10, 20, 40, 60, 80 on `b3d_pt45_802yrs_high`. The largest M that survives 20 forced years is the ceiling the forced stages must respect. 6 runs |
| **P2 ceiling against interior width** | the same bisection on the two narrower barriers, enough cells to place the ceiling. Turns a number into a relationship, and says whether M must be rescaled per barrier. 3 runs |
| **P3 the Hs 3.0 attribution** | Hs 3.0 with the groin detached, per decision 4.10. If it survives, the crash `GROIN_PLAN.md` blames on wave height was the groin's progradation, and that note needs rewording. 1 run |

**What a failure is allowed to mean.** A crash is reported as a crash, with M, the
barrier, the year and the interior row counts. It is never reported as drowning,
and the fillet from a crashed cell is never quoted. The distinction matters because
`GROIN_PLAN.md` already records M of 100 or more as having "drowned the barrier" on
the old rig, and that may have been this instead.

### D. Release -- 40 domains, 200 years, forced. 3 runs

Install at **year 1** and fail at **year 100**, not year 150. The fillet is steady
by year 40, so failing at the midpoint buys a hundred years of decay instead of
fifty at no extra cost. Both failure modes, because only one of them is ever
exercised at Buxton.

| test | what it measures |
|---|---|
| **D1 instant failure** | `deterioration_mode="instant"` to zero at year 100. The cleanest decay curve, and the one to fit a time constant to |
| **D2 ramped failure** | `linear_ramp` over the Buxton 7 years. The difference from D1 is the part of the post-2004 record the ramp is meant to explain |
| **D3 reference** | the same rig with no groin, same storms and same seed, so the decay is measured against a drifting baseline rather than against zero |

`GROIN_PLAN.md` says diffusion closes the observed gap at 1.5% of the observed
post-2004 rate. This replaces that with an e-folding time. At a 200-year ceiling
the decay will not finish, so report the fitted constant, not a recovery.

Tier 0 supplies the alongshore-only twin for both failure runs at no compute cost,
which is one of the selected pairings.

### E. Sensitivity -- 40 domains, 200 years

Three axes, per decision 4.4:

| axis | values | why |
|---|---|---|
| `M / r_ipl` | Hs 1.0 to 3.0 crossed with M 10 to 80, reported as the ratio | the fillet is bought by this ratio and nothing else (2.2), and the grid brackets the runaway boundary on the low-Hs side |
| `f` | 0, 0.2, 0.4, 0.6, 0.8, 1.0 | scales the sink only, so it sets both the fillet and whether the groin is a net source (2.4) |
| `wave_angle_high_fraction`, `wave_asymmetry` | 0.1 to 0.4, and 0.5 to 0.9 | decide the shutdown threshold outright (2.3), and move the fillet 63% from `ahf` 0.2 to 0.3 |

Report the fillet, the extent, the reach mean, the closure error and the shutdown
year for every cell. Cells flagged as shutdown are reported, never dropped.

**E1** fails if the fillet does not collapse onto `M / r_ipl`. **E2** asks whether
the 0.10 m RMSE differences that chose M = 60 over M = 70 in `GROIN_PLAN.md` are
larger than this rig's noise floor, which A5 will have measured.

36 runs: 25 for the 5 by 5 wave-height-by-M grid, 5 for the f axis at one wave
height, and 6 for the two angular parameters. The paired alongshore-only twin is
taken from tier 0 for the calibrated cell only.

### F. Island response -- 40 domains, 200 years, forced, 10-member ensemble. 20 runs

The experiment that justifies doing this in CASCADE rather than in BRIE alone.
Does updrift accretion widen the barrier, raise the dunes and cut the overwash
flux, and does the downdrift notch do the reverse? **A groin that moves the
shoreline and nothing else is a shoreline model wearing a barrier model's
clothes**, and that is worth reporting either way.

Ten storm realizations from `data/pathways_init_data/`, because overwash is noisy
and the groin's effect on it has to beat that noise. Ten more for the no-groin
twin, on the same realizations and the same seed.

| measure | where it comes from |
|---|---|
| barrier width | `InteriorWidth_AvgTS` |
| dune crest height | `DuneDomain` maximum per year |
| overwash flux | `QowTS` |
| proximity to drowning | shoreface slope and interior row count |

**Two controls, both reported.** The no-groin twin is the clean control. The rig's
own buffer domains are an in-run control that feels no groin, and above Hs 1.0 they
are contaminated from year 63 by the dipole's own tail, per 2.6. Where the two
controls disagree, the disagreement measures that contamination.

### G. Several groins -- 40 domains, 200 years, forced. 6 runs

Tier 0 already shows the field behaving like a field, so these runs test whether
Barrier3D changes that, not whether it happens.

| test | what it measures |
|---|---|
| **G1 spacing** | 3 structures at 1, 2 and 4 domains apart, which is what fits the 20 working cells. The 8 and 16-domain spacings stay in tier 0 on its 121-domain grid. 3 runs |
| **G2 starvation** | the updrift-to-interior gradient out of G1, against real groin-field surveys. Tier 0 gives 0.47 at half a kilometre rising to 0.83 at eight. No extra runs |
| **G3 staged construction** | build the three one per decade, updrift-first then downdrift-first. The order should matter and the module has no way to know about it. 2 runs |
| **G4 sequential failure** | let them deteriorate one at a time. Does the field collapse progressively or all at once. 1 run |
| **G5 aggregate equivalence** | what single M reproduces the three-structure field. Tier 0 says co-located dipoles sum exactly, so any departure here is Barrier3D and not the solve. No extra runs |
| **G6 opposed pair** | two structures facing each other, the Cape Point geometry the Buxton work keeps running into. Tier 0 only, for now |

### Deferred: the management arm

Not funded in the interview, and held against decision 4.6. It has no business in
a test of whether the dipole's signs are right, and it is the only question a
coastal manager would ask, so it is written down rather than forgotten.

| test | what it measures |
|---|---|
| road | does a groin delay the first road relocation, and by how many years |
| nourishment | does a groin reduce the fill volume needed, updrift and downdrift separately. The downdrift answer is the one the volume-neutral dipole is most likely to get wrong |
| transfer | whether a groin moves the management cost alongshore rather than reducing it. Nearly a tautology at `f = 1`, an open question at `f < 1`, which is a clean way to show what `f` means |

---

## 6. What this rig cannot test

Stated up front so nobody reports them as results.

- **Sub-grid structure.** One domain is 500 m and a real fillet is about 190 m.
  Nothing here resolves an individual groin, at any tier.
- **Grid resolution.** `r_ipl` scales as 1/dy², so the grid sets the fillet. BRIE
  fixes `dy = 500` with a "do not change" comment in the coupler, so M is
  grid-specific and the rig cannot show how specific.
- **Real drift rates.** The dipole is a shoreline displacement, not a flux. The
  rig cannot convert M into a sediment budget, which is the same limitation
  `GROIN_PLAN.md` records about the 719,000 m³/yr figure.
- **Bypassing.** There is no leakage term, so the rig cannot test one.

---

## 7. Layout

```
groin-module-test/
  TEST_PLAN.md              this file
  0-solver-audit/           BRIE's solve, no CASCADE -- exists and runs
    HAT_groin_solver_audit.py
    solver_audit_*.csv
  1-unit-rig/               5 domains, stages A to C
  2-morphology-rig/         40 domains (20 working + 10 buffer a side), P and D to G
  figures/                  house style, captions in CAPTIONS.md
```

Figures come under `scripts/site_layer/hat_figure_style.py` and put their prose in
`CAPTIONS.md` beside them, with numbers computed at draw time, as the rest of the
groin tree now does.

Each CASCADE tier owns its parameter YAML, so the runs do not fight over the
shared file and sweeps can go in parallel.
