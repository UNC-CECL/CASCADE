# coastsat_gap — when the groin stopped holding the shoreline, from CoastSat

**The narrow view.** Only the CoastSat transects either side of the groins, reduced to one number per year: the north (updrift) shoreline minus the south (downdrift) one. Its wide counterpart is `../4-shoreline_rates_by_era/`: change rates on every transect for ±60 km, by era and decade, and how far along the coast the groin's effect reaches. CoastSat only: the aerial photos on their own are `../1-aerial_photo_gap/`, and the two records side by side (including what the gap is, where they disagree, and the groin's state in 1996) are `../3-photos_vs_coastsat/`. There is no two-sides-and-gap figure like the photos' here, on purpose: CoastSat starts in 1984, 15 years into the groin's life, so it can say how the gap changed but not how big it was. Any zero taken from inside the record (1984, or the 1984-2016 mean the figures use) would draw the gap as starting from nothing partway through. That figure, with CoastSat moved onto the photos' 1967 datum, is `../3-photos_vs_coastsat/` figures 2-3. Until 2026-10-08 this folder was `coastsat_groin_condition/`; until 2026-10-14 it was `2-gap_across_groins/` and also held the photo and comparison figures (its figures 0, 7 and 8, now in `../3-photos_vs_coastsat/`).

What condition was the Buxton groin in when the calibration period starts in 1996? This folder answers that from CoastSat, independently of the wet/dry photo record behind the current model schedule. It's an observational study with no model runs. It started 2026-10-08 (moved here from the DEM-to-DEM fits the same day, because it is an observation, not a model test).

**Why it matters.** The pinned blocking groin (b 0.6, f 0.6) runs at full strength until an instant failure at the 2004 step. That timing comes from the photos, where the gap seemed to hold until 2004. The maintenance record has no repair after 1995.

**Decision rule (agreed 2026-10-08).** Fit a piecewise breakpoint to the gap. If the turn falls at 2003–04, keep full strength until 2004. If it falls around 1996 or earlier, move to an earlier or ramped failure and refit.

## Main findings (2026-10-08)

1. **The gap stopped widening in 1995, at the last repair, not at Hurricane Isabel (2003).** Both gap definitions break at 1995 (90% ranges 1992–1998 and 1991–2000). The 2002–2005 timing the model assumes holds only 1–2% of the bootstrap fits, and the break stays at 1994–96 under every robustness check. *Figures 5, 6, 8.*
2. **After 1995 the groin's effect fades slowly; it does not collapse.** The gap widened about +4 m/yr before the break and slipped about −1.4 to −2.4 m/yr after it. The 1996–2003 rate is not distinguishable from zero, and a separate 2004 drop (about −14 m) is within two standard errors of zero in every fit. *Figures 7, 8.*
3. **The 1995 peak is partly storm loss on the south side, not trapping.** The 1991–95 widening is mostly GIS 5 retreating (about −98 m against −24 m on the north side), with Hurricane Gordon in 1994. *Figure 4.*
4. **The photos disagree at one date.** The 2004 survey, which pinned b in the 2026-10-05 calibration, sits 42 m above CoastSat; see `../3-photos_vs_coastsat/`.
5. **For the model:** by the decision rule above, this points to an earlier failure. The schedule refit that followed (`../../../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`) found that, scored on CoastSat, the best groin is weak (b × f ≈ 0.32–0.36), and that failing from 1996 fits slightly best and matches the repair record. Its post-failure strength equals the current pin's, so test and forward runs would not change. The choice of schedule is still open.

**Caveats.** CoastSat before 1999 rests on 5–13 Landsat images a year per transect.

```
coastsat_gap_across_groins.py   analysis -> tables/
coastsat_gap_figures.py         figures from the tables -> figures/coastsat_gap_<n>_<what>.png; its
                                plotting helpers are imported by ../3-photos_vs_coastsat/
tables/                         every number quoted below
figures/                        eight figures in the order of the analysis; captions in figures/supporting/CAPTIONS.md
```

```
python coastsat_gap_across_groins.py   # ~30 s (bootstraps)
python coastsat_gap_figures.py
```

## Method

- **Annual means.** A calendar-year mean per transect, from at least 3 images. Each transect is taken about its own 1984–2016 mean.
- **Sides.** Each side is the mean of its transects, in years where at least half of them have data.
- **Gap.** Updrift minus downdrift, seaward positive, so a rising gap means the north side is holding. It is measured two ways:
  - **GIS 6 minus GIS 5** (6 and 13 transects): the quantity the model is fitted on;
  - **500 m either side of the field** (4 and 12): the condition check.
- **Fit window.** 1984–2016, ending before the Buxton 2017 fill lands on GIS 6.
- **Break model.** A continuous one-break hinge, with the break-year interval from 2000 residual bootstraps.
- **Alternatives.** A two-break hinge, and a separate level step at 2004, both compared by BIC.
- **Robustness checks.**
  - starting in 1988, which drops the sparse early Landsat years (5–10 images a year);
  - leaving out the 1995 peak;
  - ending in 2013.

## The figures, in reading order

Each panel title is the question the panel answers. Each caption in `figures/supporting/CAPTIONS.md` says what the figure **tests**, **how to read** it, and what it **shows**.

| # | figure | question | answer |
|---|---|---|---|
| 1 | `coastsat_gap_1_transect_map` | which transects make each gap, and where is Buxton? | 6 + 13 transects for the domain gap, 4 + 12 within 500 m. The southernmost GIS 6 transect lies just south of the southern groin; the near-field gap has no such transect |
| 2 | `coastsat_gap_2_one_transect` | how does one transect's images become a value per year? | the calendar-year mean of the images (at least 3), taken about the transect's own 1984-2016 mean; one image scatters tens of metres about its year's mean, and the early years have the fewest images. Added 2026-10-14 |
| 3 | `coastsat_gap_3_transects_to_sides` | how do the transects become one position per side? | the mean of a side's transect anomalies in years where at least half have data; before the fill they move together, after it transect 0077 (just south of the southern groin) rises with GIS 5. Added 2026-10-14 |
| 4 | `coastsat_gap_4_side_positions` | which side of the groins moved? | the 1991–95 widening is mostly the south side retreating (GIS 5 about −98 m against GIS 6 about −24 m, with Gordon in 1994). After 2017 the south side gains fill sand |
| 5 | `coastsat_gap_5_gap_and_break` | when did the gap stop widening? | 1995 in both gaps (90% ranges 1992–98 and 1991–2000) |
| 6 | `coastsat_gap_6_break_year` | which break year fits best, and how certain is it? | 1995; 2002–05, the model's failure timing, holds 1–2% of the bootstraps |
| 7 | `coastsat_gap_7_era_rates` | how are the era rates fitted, and was the gap widening or narrowing? | a least-squares line through each era's annual gap; widening before the repair, slightly narrowing 1996–2003 (not distinguishable from zero), flat 2004–16 |
| 8 | `coastsat_gap_8_robustness` | does the break hold under other fits, and is there a separate 2004 drop? | the break stays at 1994–96; no 2004 step is clearly different from zero |
| 9 | `coastsat_gap_9_placement` | which boundary does the CoastSat rate pattern point to, GIS 5/6 or 6/7? | over 1984-95 the rate steps at the southern groin, −4.7 to −7.9 m/yr south of it against −1.3 to −3.4 north, with no second step at 6/7 (GIS 6 −2.4, GIS 7 −3.3, GIS 8 −2.4 m/yr). The 5/6 gap widened +4.5 ± 1.1 m/yr, the 6/7 gap −0.2 ± 0.6. Later eras reverse the step at the southern groin and still show none at 6/7. Added 2026-10-14 |

The figures follow the analysis: transects (1), one transect to annual means (2), transects to sides (3),
the two sides (4), the gap and its break (5-6), era rates (7), robustness (8). Figures 2-3 read
`tables/transect_annual.csv`. Until 2026-10-14 the figures were `gap_across_groins_<n>_*` (numbered 0-8),
then for a day `coastsat_gap_1`-`_6`; 2-6 of those are now 4-8.


An earlier figure 8, the failure schedules, is a model setting, not an observation. Since 2026-10-08 it is `schedule_refit_0_failure_schedules.png` in `../../../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/figures/`.

## Result

| | GIS 6 minus GIS 5 | 500 m either side |
|---|---|---|
| break year (90% bootstrap) | **1995 (1992–1998)** | **1995 (1991–2000)** |
| bootstrap share at 2002–2005 / ≤ 1997 | 1% / 91% | 2% / 86% |
| rate before / after the break, m/yr | +4.2 / −1.4 | +3.3 / −2.4 |
| gap LRR 1984–95 / 1996–2003 / 2004–16, m/yr | +4.5 ± 1.1 / −2.2 ± 1.8 / +0.2 ± 0.6 | +4.2 ± 1.5 / −2.8 ± 2.5 / −0.5 ± 0.6 |
| separate level step at 2004 | −14 ± 9 m | −14 ± 12 m |
| BIC: break vs break + 2004 step vs two breaks | 170.6 / 171.4 / 172.7 | 187.7 / 189.6 / 188.8 |

**Reading.**
1. **The gap stopped widening at the 1995 last repair, not at Isabel.** Both gap definitions put the break at 1995. It stays at 1994–96 under every robustness check, and the 2003–04 timing gets 1–2% of the bootstrap mass. On the decision rule, this points to an earlier failure.
2. **What follows is a slow decline, not a collapse.** The 1996–2003 rate is not distinguishable from zero. The 2004 step is less than two standard errors in every fit. Adding the step lowers the BIC in only one of the eight fits (GIS 6 minus GIS 5, 1995 left out), and then by 0.1.
3. **The 1995 peak is partly the south side being eroded.** The gap's rise in 1991–95 is mainly GIS 5 retreating, so the peak may be partly storm damage downdrift (Gordon) rather than trapping. The photos show it too (`../1-aerial_photo_gap/`).
4. **The photos and CoastSat agree, except at 2004:** see `../3-photos_vs_coastsat/` (figure 3).
5. **For the model:** the 2026-10-05 calibration fit says b is "well pinned down by the 2004 date", which is the date CoastSat contradicts. The refit should score b and f against the annual CoastSat gap and try schedules that start weakening in 1996 (done: `../../../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`).

**Caveats.**
- CoastSat before 1999 rests on 5–13 Landsat images a year per transect.

## Last repair date (checked 2026-10-08)

CSE 2013 (Coastal Science & Engineering, *Shoreline Erosion Assessment & Plan for Beach Restoration, Rodanthe & Buxton (NC)*, report 2403-PHASE1-FR, Nov 2013, p. 20) gives the last repair as **1995**: 184 ft of steel sheet piling on the south groin, after Hurricane Gordon in 1994. This comes from the page-cited event list in the deleted `HAT_groin_zone_investigation.py` (in git at a229c8ff). The report PDF itself is not on this machine, so it wasn't re-read. No source in the repository gives 1996. "1996" appears to be the old ramp's onset (1969 + 27 years delay), read as a repair date.

Fixed in `hard-structures/groin/GROIN_PLAN.md` (now `../../data/structure_history.md`), `scripts/hatteras_ms/HAT_hindcast_methods.md` and `scripts/site_layer/README.md`.

Left at 1996 on purpose, because these are the old ramp's onset and the pinned groin's failure is the 2004 step:
- `GROIN_LAST_REPAIR_YEAR` in `scripts/hatteras_ms/groin-sweep/HAT_groin_sweep_config.py`;
- the `last_repair_year` demo in `cascade/groin.py`;
- the event labels in the `HAT-groin-figures` (now the three `figures/` folders) and `groin-sweep` figure scripts;
- the notebook's ramp comment.

These labels say "last repair" at 1996. Relabel or move them if the figures are redrawn.

## Next

Done 2026-10-08: `3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`. All three schedules reach about 10 m on the annual CoastSat gap (no groin 41.6 m). The best is failure from 1996 at b × f 0.36, a weak, constant groin. The current pin fits the photos (4.0 m) but scores 29.2 m on CoastSat. From 2009 the two are identical, because both have b × f 0.36. The choice of schedule is still open; see that README.

## Tables

| file | contents |
|---|---|
| `gap_annual.csv` | annual gap and side anomalies, with the transect count behind each value |
| `gap_breakpoint.csv` | hinge fits, bootstrap interval, two-break and step tests, BIC |
| `gap_robustness.csv` | break year, its interval, the step and BIC under each check |
| `break_profile.csv`, `break_bootstrap.csv`, `hinge_fit.csv` | the curves behind figures 5–6 |
| `era_rates.csv` | gap rate ± SE per era, both CoastSat gaps |
| `transects_used.csv` | transect ids, side and UTM ends |
| `boundary_gaps.csv` | the GIS 5, 6 and 7 side means (each about its 1984-2016 mean) and the gap across 5/6 and 6/7 (figure 9) |
| `alongshore_rates.csv` | every GIS 4-8 transect's least-squares rate in each era and after the fill, with its northing from the southern groin (figure 9) |
| `transect_annual.csv` | each transect's annual mean, image count, whether the year is used (at least 3 images), its 1984-2016 mean and the anomaly (figures 2-3) |
