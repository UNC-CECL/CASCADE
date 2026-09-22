# Decision: keep Hs = 2.5

**Decided 2026-09-01.** Evidence in this directory, plus the sensitivity sweep
under `output/calibration/sensitivity/`.

---

## The question

Does raising the wave height make the model reproduce the observed shoreline
change well enough on its own that the calibrated source/sink field needs to do
less work? If so, Hs is compensating for a wave-transport deficit and should be
raised; if not, the field is correcting something else and raising Hs only moves
the error around.

## What was run

Two pass-0 calibrations of the same model, differing only in wave height. Same
script (`be_zone_residual_fit.py`), same base scenario
(`edgeBE road_bdm`, groin on), **groin frozen at the production M = 60, f = 0.6**
so the comparison isolates Hs. Both arms derived fresh; neither read from the
stored calibration, and `HAT_BE_OUTPUT_DIR` kept both out of it.

## Result: a ~6% smaller field, over the same number of zones

| | 1984–2004 | 2004–2024 |
|---|---|---|
| total \|BE\| required | 64.86 → **61.09** (−5.8%) | 118.26 → **111.24** (−5.9%) |
| RMS residual | 1.216 → 1.155 (−5.0%) | 1.760 → 1.752 (−0.5%) |
| domains needing correction | 43 → **46** (+3) | 63 → **61** (−2) |

The magnitude of correction falls about 6% in both periods — the most robust
number here, being an integral rather than a threshold count. **The number of
zones does not fall.** Period 1 needs three more, period 2 two fewer.

## And it redistributes rather than reduces

| zone | 1984–2004 | 2004–2024 |
|---|---|---|
| Pea Island NWR (D84–90) | **−24%** | **−37%** |
| Buxton–Avon Transition (D11–20) | **−14%** | **−12%** |
| Avon (D21–31) | +4% | +2% |
| Mid-island (D32–59) | **+8%** | **+6%** |

Four cells, two independent 20-year windows, same direction in all four: the far
north and the Buxton–Avon transition improve, the mid-island and Avon get worse.
That is a reproducible spatial signal, not noise. Wimble Shoals (+19% / −4%) and
Cape Point (+4% / +29%) disagree between periods and should not be read.

Reading: more alongshore diffusivity helps where alongshore gradients are steep
— the north end, near Oregon Inlet — and hurts through the middle reach, where
the residual is not transport-shaped and the extra smoothing works against the
fit. So the mid-island correction is **not** standing in for wave transport.

## Why that is not enough to change Hs

**The shoreline skill difference is inside the noise.** Interior RMSE gains at
Hs 3.0 are 0.84 standard errors (period 1) and 0.11 (period 2); 2.5 sits within
one SE of the minimum in both. Those SEs are optimistic — the LOESS-10 target
correlates neighbouring domains, so effective *n* is far below 88.

**It was already tested higher and rejected.** `GROIN_PLAN.md`: at Hs 3.5 the
fillet decay improves from −15.3 m to −20.6 m against −24.4 m observed —
confirming the mechanism — **but the reach RMSE goes 15.97 → 36.58.**

**3.0 is further from the observed wave climate, not closer.** USGS hindcast
mean is ~1.2–1.3 m, morphologically effective ~1.4–1.5 m. The model preferring
more wave energy than the coast receives is a finding about the model, most
likely a missing erosion mechanism. Raising Hs would hide that signal.

**The cost is a full re-derivation.** Groin sweep → joint fit → base run → BE
zone analysis → config → matrix. Six steps, and calibBE RMSE would land back
near 0.54 whatever Hs is chosen, because a 48-domain residual fit absorbs
whatever the physics does not supply.

## What would overturn this

- The **mid-island degradation reversing** at some Hs between 2.5 and 3.0. Only
  two points were tested; the managed-arm sweep put the middle reach's own
  optimum near Hs 1.2, so the interior may want *less* energy, not more.
- An **independent** constraint on Hs — a wave hindcast or buoy record for this
  reach — putting the true value near 3 m. That would make this a model-error
  question rather than a tuning one.
- The **roadway module's +0.19 m/yr bias** turning out to be a defect. Fixing it
  moves the managed arm's residual, which is what the calibration rests on, and
  the Hs question should then be asked again from scratch.

## Recorded limits

- Two wave heights, not a curve.
- The zone rule has a hard 0.5 m/yr threshold, so counts move for domains
  sitting either side of it without much meaning. The magnitude totals do not
  have this problem and should be preferred.
- Period 2's residual is dominated by a near-uniform +0.99 m/yr offset that Hs
  barely moves (→ +0.93). Its insensitivity was predicted before the run.
- The groin was held fixed. This says nothing about what M should be at Hs 3.0,
  and M was itself fitted at 2.5.

## Note

Hs = 3.0 could not be tested before today: it crashed on a pre-existing CASCADE
bug in `beach_dune_manager.resize_interior_domain` (recorded in `GROIN_PLAN.md`
as "interior grid 136 vs 137 rows"). That bug was fixed on 2026-09-01 and
verified to reproduce prior runs byte-for-byte, which is what made this
measurement possible at all.

---

# ADDENDUM: THE EDGE RE-SOLVE — 2026-09-01

**The decision above stands, and period 1 now supports it instead of merely
failing to overturn it.**

## What the measurement above left open

Both arms were scored with the two locked end domains held at values solved at
**Hs = 2.5**. `Recorded limits` names the frozen groin but not the frozen
edges, and for `edgeBE` that omission is not a detail — the preset imposes
background erosion at GIS 1 and GIS 90 and **nowhere else**, so those two
numbers are the entire forcing it applies. The test arm was therefore wearing
the control arm's boundary condition, which is the one confound the ~6%
result could not speak to.

## What was run

Newton-on-a-secant at both end domains, per period, at Hs = 3.0. The gain was
**re-measured** with a +3.0 m/yr probe rather than read off the Hs 2.5 table,
because BRIE's diffusivity scales as `Hs^2.4 / (h_b_crit + d_sf)` with
`d_sf = 8.9*Hs` (`brie.py:383-385`, `:270`) — about +31% from 2.5 to 3.0 — so
more of an imposed edge rate diffuses away and the gain must fall. It did:

| | gain @ Hs 2.5 | gain @ Hs 3.0 | BE solved at Hs 3.0 |
|---|---|---|---|
| P1 GIS 1  | 0.109 | 0.0989 (−9%)  | −42.6 → **−51.99** |
| P1 GIS 90 | 0.104 | 0.1048 (+1%)  | +13.0 → **+27.66** |
| P2 GIS 1  | 0.096 | 0.0896 (−7%)  | +50.3 → **+42.42** |
| P2 GIS 90 | 0.099 | 0.0888 (−10%) | +46.7 → **+52.38** |

One Newton step was enough in both periods. Edge residuals against the CoastSat
target, before and after:

| | P1 GIS 1 | P1 GIS 90 | P2 GIS 1 | P2 GIS 90 |
|---|---|---|---|---|
| Hs 3.0, edges from the 2.5 solve | −0.928 | +1.537 | −0.706 | +0.504 |
| Hs 3.0, edges re-solved | **+0.153** | **+0.211** | **−0.007** | **−0.005** |
| Hs 2.5, the control arm's own | +0.559 | +0.572 | +0.234 | −0.510 |

Note the third row. **The re-solved Hs 3.0 edges are converged TIGHTER than the
control arm's**, so what follows cannot be explained by the test arm being the
less-converged of the two.

## Result: the period-1 gain was the stale edges, and it reverses

Interior RMSE, m/yr:

| | 1984–2004 | 2004–2024 |
|---|---|---|
| Hs 2.5, edges solved | **1.2283** | **1.7798** |
| Hs 3.0, edges from the 2.5 solve | 1.1592 (−5.6%) | 1.7666 (−0.7%) |
| Hs 3.0, **edges re-solved** | **1.3223 (+7.7%)** | **1.7651 (−0.8%)** |

**Period 1 changes sign.** The 5.6% improvement that the sensitivity sweep
reported, and that this directory's first measurement inherited, was the
Hs = 3.0 run wearing a boundary condition fitted for a different wave climate.
Solve the boundary properly and Hs 3.0 is 7.7% **worse** than Hs 2.5.

**Period 2 remains a non-result**: 0.8%, against the estimate above that its
gains run ~0.11 standard errors. Its mean interior bias is **−0.939 m/yr** with
the edges solved — still the near-uniform offset that sits there at every wave
setting. Wave forcing is not the period-2 skill lever, as predicted before the
run.

## Why the interior gets worse when the edges get better

The end domains are not isolated from the interior: their BE diffuses inward.
At Hs 3.0 the edges needed far larger corrections — GIS 90 moves +13.0 → +27.66
in period 1 — and the interior pays for them. The note above measures this
trade at Hs 2.5 as "0.03 RMSE at the interior for 0.61 at the edge", about
20:1. Here it is ~13:1 (0.163 interior for ~2.1 of edge residual): the same
order, still contracting rather than oscillating, but the interior bill is
roughly **five times larger** because the steps themselves were so much bigger.

That is the mechanism, and it is consistent with the redistribution finding
above rather than a separate effect. More alongshore diffusivity moves the
model's error around; it does not remove it.

## What this adds, and what it still does not

**Adds.** The edge confound is closed. The verdict no longer rests on a ~6%
parsimony gain that could have been an artefact of the boundary — period 1 now
argues against Hs 3.0 on skill directly, under the non-circular `edgeBE` arm.

**Does not.** The groin is STILL frozen at M = 60, f = 0.6. The refit was
planned and then dropped as moot: at +31% diffusivity, holding the fitted
fillet size would drive M from 60 to about 78, further past the reach sediment
budget, so refitting could only compound a loss it cannot reverse. If Hs is
ever revisited, that refit is still owed.

**Still only two wave heights**, and only one Newton iteration. Period 1's
residuals stop at +0.15 / +0.21 rather than zero; a second pass would move the
edges by ~1.5–2 m/yr against first steps of 9–15, so it is very unlikely to
close a 0.094 m/yr interior gap, but it was not run.

## Provenance

Runs: `output/raw_runs/waveHs3_probe/` (gain probes) and
`output/raw_runs/waveHs3_edge1/` (solved edges), both periods, `edgeBE`,
`full_management`, groin on at M = 60 / f = 0.6, topography `1984-start/v1`
and `2004-start/v1`. `be_values_digest` 46ef9e025244 (P1) and 5496589cb749 (P2)
— distinct from the base runs' 0fd4a97a4912 / ac6a7aca07f6, so an edge-solved
run cannot be mistaken for the preset it started from.

**The production config was never edited.** Two variables were added to the
runner for this: `HAT_BE_OVERRIDE` forces per-domain BE for one run, and
`HAT_ARM_TAG` scopes `output/raw_runs/<arm>/` and the `arm` column in
`run_index.csv` from one string, so a probe cannot land on the run it probes.
With both unset the runner reproduces its pre-change output **bit for bit** in
BOTH periods -- rate CSV and shoreline matrix, sha256 identical, max |delta|
0.0 m/yr -- so nothing here changes any existing run. Worth noting that period
2 matched EXACTLY, where `hatteras_site_config.py` records it as normally
reproducible only to ~1e-4; that tolerance was measured on the sweep-worker
path, and this re-run went through the runner.

**One bug was introduced and fixed in the course of this.** `append_run_index`
keyed on `("run_name", "Hs_m")`; a probe shares both with the run it probes and
differs only in `arm`, so filing the first probes DELETED the two base rows.
The key now includes `arm`, the rule being that it must name every component of
`OUTPUT_BASE_DIR` that the run name does not. The two rows were restored by
re-running the base runs, which is what produced the bit-for-bit check above.
