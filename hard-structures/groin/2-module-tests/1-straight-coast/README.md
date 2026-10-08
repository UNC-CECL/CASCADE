# 2026-10-08-straight-coast

**Question.** Before fitting any groin to Buxton, does each version of the groin module behave like a groin on a coast with nothing else going on? It is tested on straight coasts at orientations from −40° to +40°.

**Setup.** `straight_coast_test.py` is BRIE's alongshore step alone: its own `coast_diff` table, sparse indices and clip at zero, under the option-A waves (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle 0.5). No storms, no Barrier3D, no sea level. The coast is a straight line at angle θ0 plus a periodic perturbation, so the tilt never meets BRIE's periodic wrap. There are 61 domains and one groin at the middle face, run for 50 years. One step of the emulator matches a real `Brie.update()` exactly (max difference 0.0 m on an arbitrary shoreline).

| folder | module | what it does each year |
|---|---|---|
| `1-source-sink-dipole/` | `GroinCallback` | −M updrift, +M downdrift (M 12), whatever the coast does |
| `2-trapping-pinned/` | `BlockingGroinCallback`, the 10-05 pin | cancels b (0.6) of each cell's own diffusive coupling to the face |
| `3-trapping-conserving/` | the 10-08 variant | the same, one face coefficient on both sides |
| `4-trapping-drift/` | new, here only | cancels b of the wave-climate net drift Q_net arriving at the groin |

Each folder holds `tables/` (by year, and profiles with the updrift side at positive offsets) and `figures/`. `comparison/` holds the cross-module figures, the year-50 summary and a Pelnard-Considère check.

**The drift.** BRIE moves the shoreline with a diffusivity D(θ). It never uses the net drift itself, which it computes only for inlets (`_coast_qs`). Q_net here is the wave pdf convolved with `_coast_qs`, aligned as BRIE aligns `coast_diff`. That reproduces BRIE's diffusivity: D = −(1/depth) dQ_net/dθ, to within one 1° bin (depth = h_b_crit + d_sf = 19.8 m). On a coast aligned with the grid, 0.91 M m³/yr moves toward lower domain numbers (southward at Hatteras). The drift falls to zero at about +32° and reverses beyond it. Module 4 reads Q_net at the angle of the updrift approach face, the Pelnard-Considère boundary condition. A first version read it across the groin's own step, which rotates the wrong way as the fillet grows; that was corrected.

## Findings

1. **High-angle waves freeze a grid-aligned coast in BRIE.** At 0°, D is −119 m²/yr: high-angle waves (fraction 0.5) cancel the low-angle smoothing. BRIE clips it to 0, which also discards the instability high-angle waves would drive. The drift is still 0.91 M m³/yr. The pinned and conserving groins block only the diffusive part of the transport, so at 0° they do exactly nothing. The drift groin traps there. The dipole adds sand regardless.
2. **The pinned and conserving groins trap on the side set by the coast's tilt, not by the drift.** They block D·slope, which is Q_net(θ) − Q_net(0). For θ0 > 0 that difference runs against the real drift, so they build the "fillet" on the downdrift side (panels d-f of `comparison/figures/straight_coast_profiles_year50.png`). For θ0 < 0 it has the right sign, but it is a small fraction of the drift.
3. **BRIE's alongshore solve does not conserve sand, and a groin is exactly what exposes it.** The row-scaled scheme takes each cell's diffusivity from its own forward face. Where neighbouring faces differ by many degrees, as they do at a groin, it creates sand (θ0 < 0) or destroys it (θ0 > 0). After 50 years this amounts to thousands of metres of summed shoreline for the dipole, conserving and drift modules. On a coast at +10° to +30°, the solve erodes away everything the groin traps.
4. **The pinned groin is the only one that conserves sand at system level.** Its own imbalance offsets the solve's: at every orientation, the total change across the reach stays within 30 m of summed shoreline in 50 years. It cancels a share of BRIE's own row-scaled exchange, so it inherits BRIE's asymmetry with the opposite sign. The "conserving" variant balances the module on paper, but it leaves the solve's error in place. The same holds in the Buxton runs, 1996-2009, all 120 domains, total change against no groin: pinned −22 m, conserving +186 m, dipole +138 m. So the 2026-10-08 finding that the pinned groin "loses 186 m" measured only the module's own books, and the conserving-groin grid was stopped before it ran.
5. **Pelnard-Considère does not apply here.** The drift is so large relative to D·depth that tan α is above 0.2 at every orientation below +25°, so linear theory is invalid. Where it is valid (+25° to +40°, near zero drift), the trapped sand is too small to compare against.

## What this means for the module

None of the four is right yet. The pinned groin conserves sand but traps the wrong quantity. The drift groin traps the right quantity but sits on a solve that does not conserve sand. A groin that both blocks the real drift and conserves sand would need the shoreline step in flux form, (Q_net(θ_i) − Q_net(θ_{i−1})) / (dy · depth), or the drift groin's injection corrected for the solve's error at its face. Either is a change to BRIE or to how CASCADE drives it, so it is a decision for later, not part of this test.

Run: `python straight_coast_test.py [--years 50]` (about a minute).

`straight_coast_gif.py` animates the same runs year by year: six orientations, the four modules overlaid, fixed y ranges, and the total sand change across the reach per module in each panel. Output: `comparison/figures/straight_coast_over_time.gif` (51 frames, 5 per second).
