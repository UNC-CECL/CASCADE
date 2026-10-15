# 7-mirror — does the groin handle drift from either direction the same way?

**Question.** Buxton's drift runs south, so every run so far has the drift arriving from one side. A sign or
side error that only shows when the drift comes from the other direction would go unseen. The standard check:
run the mirror image of a coast and see whether the groin gives the mirror image of the result. Asked for on
2026-10-15, with the split-sample validation (`../../4-hindcast/2-blocking-1996-2025/2026-10-15-split-sample/`).

**Setup.** The straight-coast rig of `../6-stress/` with the runner's `CorrectedGroinCallback`. Mirroring
reverses the domain order: the coast at θ0 becomes −θ0, the groin face LO|LO+1 becomes NY−2−LO, and the waves
come from the other side. To keep BRIE's table discretization out of it, the mirror gets BRIE's tables
reflected exactly (part C says why).

```
mirror_check.py   parts A-C (about 2 minutes)
tables/           mirror_one_year.csv (A), mirror_runs.csv (B), brie_table_symmetry.csv (C)
```

    python mirror_check.py

## Findings

1. **The groin's own direction handling is exact (A).** One call of the callback on a shoreline and on its
   reversal returns exactly reversed changes: 42 of 42 cases (seven orientations, strengths 0.25 and 1,
   no step or a ±150 m step) agree to 0.0 m, while the groin moves up to 110 m. Its choice of side, its sign
   and the face it reads the arriving drift on all reflect correctly.
2. **No whole run mirrors exactly, because BRIE does not (B).** BRIE's alongshore step takes each row's
   diffusion number from that cell's forward face, so reversing the coast changes the arithmetic. With no
   groin at all, a ±150 m step and its mirror end 1-109 m apart after 100 years. On a straight coast, with
   nothing for BRIE's asymmetry to act on, groin runs and their mirrors end within 0-6.5 m of each other,
   against changes of 2-231 m. With a starting step the groin adds at most 24 m to BRIE's own asymmetry
   (−10°, strength 1), and at 0° it removes 36-69 m of it, since its correction pulls BRIE's step toward the
   face-centred (symmetric) form. At Buxton-like settings (strength 0.25, −30° to −10°) it adds at most 10 m.
3. **BRIE's own tables are not exactly mirror-symmetric (C).** Its tables for asymmetry 1 − a differ from an
   exact reflection of those for a by up to 4.5% (diffusivity) and 3.9% (drift), because its 180 angle bins
   are 1.006° wide but are looked up as if 1° wide. A coast with the drift arriving from the north and one with
   it from the south are therefore not quite the same model in BRIE, groin or not.

**What this means.** The groin treats both drift directions identically. What does not mirror in a whole run
comes from BRIE's alongshore step and its angle tables, not from the groin. For reporting: the module is
direction-independent; BRIE's alongshore solve carries a small directional bias of its own.
