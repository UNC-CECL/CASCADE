# photos_vs_coastsat — the two records of the gap, side by side

**The comparison.** The gap across the Buxton groin is measured twice: by the aerial photos
(`../1-aerial_photo_gap/`, 24 dated surveys 1967-2023, change since 1967 against a fixed offshore datum)
and by CoastSat (`../2-coastsat_gap/`, annual means 1984-2025, each transect about its own 1984-2016
mean). This folder puts them on one scale, says where they agree and where they don't, and reads the
two together for the groin's state in 1996, the start of the calibration period. Observations only; no
model runs. Made 2026-10-14 from figures 0, 5, 7 and 8 of what was `2-gap_across_groins/` (7 now folded into figure 3); the numbers are
unchanged.

```
photos_vs_coastsat.py           reads both folders' tables -> tables/
photos_vs_coastsat_figures.py   figures -> figures/photos_vs_coastsat_<n>_<what>.png; borrows the
                                plotting helpers from ../2-coastsat_gap/coastsat_gap_figures.py
tables/                         photo_vs_coastsat.csv, coastsat_sides_on_1967_datum.csv,
                                photo_minus_coastsat_by_side.csv, coastsat_one_day_spread.csv, era_rates_both.csv
figures/                        five figures in reading order; captions in figures/supporting/CAPTIONS.md
```

```
python ../1-aerial_photo_gap/photo_gap.py               # first, if its tables are missing
python ../2-coastsat_gap/coastsat_gap_across_groins.py  # likewise
python photos_vs_coastsat.py
python photos_vs_coastsat_figures.py
```

## Method

- **Putting them on one scale.** The two datums differ, so the photo gap is shifted onto the CoastSat
  whole-domain gap (GIS 6 minus GIS 5) by the mean difference over the photo years inside 1984-2016,
  the pre-fill window both cover. Only the shapes compare. For figures 2, 3 and 5 the shift is run the other way,
  CoastSat onto the photos' 1967 datum: each side onto the same side's photo positions over those years
  (figures 2-3), or the gap onto the photo gap (+117 m, figure 5; the two give the same gap). CoastSat starts
  in 1984, 15 years into the groin's life, so only the photos can say how big the gap was; zeroing CoastSat
  on its own first year would draw the gap as starting from nothing in 1984.
- **Era rates.** Least-squares slopes over the same eras for both records: to the 1995 last repair, 1996
  to Isabel (2003), and 2004 to the 2017 fill. CoastSat on annual means, the photos on the surveys
  inside each era.

## Main findings

1. **They agree except at one date.** The yardstick is how far one CoastSat day sits from its year's
   mean (a photo is one day): a typical SD of 11 m, so ±23 m at 2 SD. Before the fill all 13 north and
   13 south photo-minus-CoastSat differences, and 12 of 13 gap differences, fall inside it (RMS 10, 10
   and 14 m). The exception is the 2004 gap, +42 m (north +22, south −19). 2004 is also CoastSat's most
   scattered year, after Isabel, so against that year's own scatter it is 1.8 SD. That single survey is
   still the only support for "the groin held until 2004", and the date that pinned b in the 2026-10-05
   calibration. *Figure 3.*
2. **After the fill both sides of the photos sit landward of CoastSat, but the gap still agrees.** North
   −18 m and south −24 m on average (5 of 12 side differences outside the band), gap +6 m (all inside):
   a shift common to both sides, not a disagreement about the groin. A wet/dry line and a satellite
   waterline could part that way on the wider, re-graded fill beach; nothing here confirms it. *Figure 3.*
3. **Their rates agree before the repair and split after 2004.** 1984-95: +4.5 ± 1.1 m/yr (CoastSat) and
   +4.9 ± 1.5 m/yr (photos). 1996-2003: −2.2 ± 1.8 m/yr from CoastSat; the photos have only two surveys.
   2004-16: CoastSat flat (+0.2 ± 0.6) while the photos fall −3.3 ± 0.8 m/yr, a fall that starts from the
   high 2004 survey. *Figure 4.*
4. **The groin in 1996: a full-sized gap, no longer widening, not failed.** 135 m since 1967 in the
   photos (87% of the 155 m peak in 1995), a year after the last repair; CoastSat breaks in 1995 and
   changes −2.2 ± 1.8 m/yr over 1996-2003, and nothing marks a collapse before the 2015-2021 losses.
   The 1995 peak is partly Hurricane Gordon's (1994) loss on the south side. *Figure 5.*

**Caveats.** A photo is one day, while a CoastSat value is a year's mean. CoastSat before 1999 rests on
5-13 Landsat images a year per transect. Values here read the two folders' 2-decimal tables, so they can
differ from a full-precision run by 0.01 m.

## Where should the model groin sit: GIS 5/6 or 6/7? (2026-10-14)

The groin field fills almost all of GIS 6: the southern groin is about 75 m inside the 5/6 boundary, the
northern groins end 25-35 m short of 6/7. So the two candidate boundaries are the two ends of one field.

On one scale (figure 6) the 5/6 gap is the groin's signal: it built to 155 m by 1995 and widened +4.5 m/yr
in 1984-95. The 6/7 gap peaked at 44 m (1993), was 31 m in 1996 and never widened while the groins worked.
Over the calibration (1996-2009) both fell by a similar amount (5/6 −19 m, 6/7 −16 m); over the test
(2009-2025) 5/6 fell 60 m and 6/7 only 2 m. For the runs:
- a groin at 6/7 should be scored on the 6/7 gap (`tables/placement_targets.csv`, column `6/7`), not the
  5/6 one every fit has used so far;
- whichever boundary it sits on, the other boundary's observed gap is a check: a 6/7 groin should still
  leave GIS 6 holding against GIS 5 (5/6 gap ~140 m in 1996), and a 5/6 groin should leave a small
  6/7 gap.
The photo and CoastSat profiles behind this are `../1-aerial_photo_gap/` figure 6 and `../2-coastsat_gap/`
figure 10.

## The figures, in reading order

| # | figure | question | answer |
|---|---|---|---|
| 1 | `photos_vs_coastsat_1_what_the_gap_is` | what is the gap, and how does each record measure it? | the north-side position minus the south-side position; both sides retreated (by 1996 the north 49 m, the south 183 m since 1967), and the 135 m between them is the groin's effect. Was `gap_across_groins_0`; its animation is now `../1-aerial_photo_gap/figures/photo_gap_4_two_sides.gif` |
| 2 | `photos_vs_coastsat_2_one_scale` | how are the two records put on one scale? | before: CoastSat sits about its own 1984-2016 mean, far from the photos; after: each CoastSat side moved by one constant (north −45 m, south −162 m) onto the same side's photos over their 13 shared pre-fill years, and the lines run through the squares. Added 2026-10-14 |
| 3 | `photos_vs_coastsat_3_gap_and_agreement` | how wide was the gap in both records, and how far is each photo from CoastSat? | 142 m in 1996 (CoastSat), a slow slip with no step at Isabel, 64 m by 2026. Before the fill every difference is inside the ±23 m one-day scatter but the 2004 gap (+42 m); after it both sides of the photos sit 18-24 m landward while the gap agrees. Added 2026-10-14; it replaces `photos_vs_coastsat_3_do_they_agree` (was `gap_across_groins_7`), whose differences are panel (b) here |
| 4 | `photos_vs_coastsat_4_era_rates` | do the two records give the same rates? | yes before the repair; after 2004 CoastSat is flat while the photos fall. Was `gap_across_groins_5`, without the near-field gap (in `../2-coastsat_gap/` figure 8) |
| 5 | `photos_vs_coastsat_5_state_in_1996` | what state was the groin in at 1996? | a full-sized gap a year after the last repair, no longer widening, no collapse before 2015-21. Was `gap_across_groins_8` |
| 6 | `photos_vs_coastsat_6_placement` | what is the observed gap across GIS 5/6 and across 6/7, on one scale, at the calibration and test keys? | 5/6: 142 m in 1996, 123 m in 2009, 63 m in 2025. 6/7: 31, 15 and 13 m; it never widened while the groins worked (−0.2 ± 0.6 m/yr over 1984-95, against +4.5 ± 1.1 at 5/6). The two records agree to 14 m RMS at 5/6 and 5 m at 6/7. Targets in `tables/placement_targets.csv`. Added 2026-10-14 |

The figures follow the analysis: what the gap is in each record (1), putting them on one scale (2),
comparing them (3), their rates (4), what they say together about 1996 (5).

## Tables

| file | contents |
|---|---|
| `photo_vs_coastsat.csv` | each photo year: the photo gap since 1967, shifted onto CoastSat, the CoastSat whole-domain gap, and photo minus CoastSat. Was `../2-coastsat_gap/tables/photo_gap.csv` |
| `coastsat_sides_on_1967_datum.csv` | each CoastSat year: north and south moved onto the photos' 1967 datum, side by side, and the gap |
| `photo_minus_coastsat_by_side.csv` | each photo year CoastSat covers: photo minus CoastSat for north, south and the gap, on the 1967 datum |
| `coastsat_one_day_spread.csv` | each CoastSat year: the SD of one image day's side position (and gap) about the year's mean, and the image count; the yardstick for agreement |
| `boundary_gaps_on_1967_datum.csv` | the CoastSat gap across GIS 5/6 and 6/7, each moved onto the photos' gap at that boundary (figure 6) |
| `placement_targets.csv` | for each boundary: the shift, the photo-CoastSat RMS, the gap at 1996, 2009 and 2025 (CoastSat, and photos where a survey exists), its change over the calibration and the test, the era rates, and the photo peak. The numbers a placement run is scored on |
| `era_rates_both.csv` | gap rate ± SE per era: CoastSat whole domains (copied from `../2-coastsat_gap/tables/era_rates.csv`) and the photos |
