# 5-fine-grid-reference — a real-sized groin field as the yardstick

**Question.** A real groin field is smaller than one 500 m domain. Which version of the groin module, with one groin on BRIE's grid, best reproduces what a real-sized field does to the domains either side of it?

**Why.** The scorecard in `../4-groin-field/` checks behaviour (right side, sand conserved, a sawtooth, a cut below the field). It cannot say whether the amounts and timing are right, and its groins are 1 to 10 km apart, while real groins are a few hundred metres apart. This builds the comparison the groin literature implies: a fine-grid shoreline model with real-sized groins that bypass sand once full, as in GENESIS (Kraus, Hanson & Blomgren 1994; see `../4-groin-field/README.md`), averaged to 500 m.

**The reference.** A one-line shoreline model in flux form on 25 m cells, 61 domains long, explicit in time. It is driven by the same wave-climate physics as BRIE, option A waves:
- BRIE's diffusivity, clipped at zero as BRIE clips it, is integrated into a drift curve, anchored at the real net drift of the undisturbed coast. Clipping keeps the high-angle-wave instability from growing at the 25 m grid scale, as BRIE's clip does at 500 m.
- **Check:** with no groins, a 40 m shoreline bump (σ 2 km) run 10 years on both grids differs by 0.03 to 0.19 m RMS at −20°, 0° and +20°. Sand is conserved exactly (net change below 10⁻¹¹ m in every reference run).

**The field.** Three groins 200 m apart, the downdrift one on the coarse groin face and the other two inside the updrift domain, as at Buxton (three groin positions about 200 m apart inside GIS 6; digitized lengths 90 to 155 m, `../../1-observations/gis_data/groins_hatteras.geojson`). Each groin passes (1 − b)·Q of the drift Q arriving from updrift. It blocks b = 0.6 × (1 − a/L)^(2/3), where a is how far the cell just updrift has built out and L is the groin length. This is GENESIS's bypassing with the depth at the tip from a Dean profile, scaled so the groin starts at the modules' blocking fraction (0.6), and passes everything once the fillet reaches the tip. L is 110 m, with 55 and 220 m as a sensitivity.

**The modules.** Dipole, pinned, drift and corrected, one groin on the coarse face, run with `../4-groin-field/groin_field_test.py`'s `run_field()` unchanged. They are compared on the gap between the updrift and downdrift domains over years 1 to 50, the quantity measured at Buxton, and on the shoreline over 22 domains at year 50.

```
fine_grid_reference.py           the reference runs, the module runs and the tables (about 7 minutes)
fine_grid_reference_figures.py   both figures, from the tables (seconds)
tables/                          fine_reference_by_year.csv (gap, both sides and net change, every source and year),
                                 fine_reference_profiles.csv (year-50 shorelines), fine_reference_scores.csv (errors)
figures/                         fine_reference_1_gap_and_profiles, fine_reference_2_errors; PDFs and CAPTIONS.md in supporting/
run_log.txt                      the 2026-10-09 run's printed tables
```

    python fine_grid_reference.py
    python fine_grid_reference_figures.py

## Findings

| module | gap error, mean over orientations | shoreline error |
|---|---|---|
| corrected | **71 m** | **25 m** |
| drift | 103 m | 121 m |
| dipole | 116 m | 51 m |
| pinned | 241 m | 67 m |

1. **The real-sized field is limited by its groin length, not by the waves.** It fills to its tips within 5 to 20 years at −30° to +20°, then holds a gap of about twice the groin length at every one of those orientations: 108 to 113 m for 55 m groins, 213 to 227 m for 110 m groins, 441 to 461 m for 220 m groins. Only at +30°, where the drift is a tenth of its grid-aligned value, is it still filling at year 50 (136 m). This is the literature's picture: a compartment fills to the tip and then bypasses.
2. **The corrected version is closest, by a wide margin on coasts turned against the drift.** At −30° and −20° its gap ends at 215 and 214 m against the reference's 225 and 227 m (gap error 23 and 20 m, shoreline error 8 and 9 m). That is partly coincidence: its fillet levels off because the shoreline turns toward the waves, and here that limit happens to land near what 110 m groins hold.
3. **Without a groin length, every module misses at other orientations.** At 0° the corrected version overshoots, 328 m against 223 m, and is still growing. At +20° it undershoots, 44 m against 222 m: BRIE's diffusivity there is five times larger than at −20° and spreads the fillet. The reference still fills to its tips. The drift version overshoots at most orientations (343 m at −20°) and shows a numerical zigzag. The pinned version builds its gap on the wrong side from +10° on (−230 m at +20°).
4. **What this suggests for Buxton.** The observed gap reached 117 m by 1978, about ten years after installation, and peaked at 150 to 155 m in 1995-2004. That is the speed and the kind of levelling-off the reference shows. A gap of about twice the groin length would put Buxton's effective length near 75 m, shorter than the digitized 90 to 155 m, which include the part of each groin on the beach. This is a hypothesis for the hindcast stage, not a result: the record also carries storms, the 2003 damage and the post-1995 deterioration.

**What this means for the module.** Of the existing versions, the corrected one is the best stand-in for a real-sized field, and it is the right base. To match the reference across orientations it needs a groin length: blocking that falls to zero as the updrift cell builds out toward the tip, as in the reference. That is a change to the module. The next step would be to add it as a seventh straight-coast variant and score it against this reference before any Buxton run.
