# gap_across_groins — when the groin stopped holding the shoreline

**The narrow view.** Only the CoastSat transects either side of the groins, reduced to one number per year: the north (updrift) shoreline minus the south (downdrift) one. Its wide counterpart is `../shoreline_rates_by_era/`: change rates on every transect for ±60 km, by era and decade, and how far along the coast the groin's effect reaches. Until 2026-10-08 this folder was `coastsat_groin_condition/`.

What condition was the Buxton groin in when the calibration period starts in 1996? This folder answers that from CoastSat, independently of the wet/dry photo record behind the current model schedule. It's an observational study with no model runs. It started 2026-10-08 (moved here from the DEM-to-DEM fits the same day, because it is an observation, not a model test).

**Why it matters.** The pinned blocking groin (b 0.6, f 0.6) runs at full strength until an instant failure at the 2004 step. That timing comes from the photos, where the gap seemed to hold until 2004. The maintenance record has no repair after 1995.

**Decision rule (agreed 2026-10-08).** Fit a piecewise breakpoint to the gap. If the turn falls at 2003–04, keep full strength until 2004. If it falls around 1996 or earlier, move to an earlier or ramped failure and refit.

## Main findings (2026-10-08)

1. **The gap stopped widening in 1995, at the last repair, not at Hurricane Isabel (2003).** Both gap definitions break at 1995 (90% ranges 1992–1998 and 1991–2000). The 2002–2005 timing the model assumes holds only 1–2% of the bootstrap fits, and the break stays at 1994–96 under every robustness check. *Figures 3, 4, 6.*
2. **After 1995 the groin's effect fades slowly; it does not collapse.** The gap widened about +4 m/yr before the break and slipped about −1.4 to −2.4 m/yr after it. The 1996–2003 rate is not distinguishable from zero, and a separate 2004 drop (about −14 m) is within two standard errors of zero in every fit. *Figures 5, 6.*
3. **The 1995 peak is partly storm loss on the south side, not trapping.** The 1991–95 widening is mostly GIS 5 retreating (about −98 m against −24 m on the north side), with Hurricane Gordon in 1994. *Figure 2.*
4. **The aerial photos agree with CoastSat except at one date.** Before the fill they agree to within 13 m in every photo year except 2004, which sits 42 m above CoastSat. That single survey is the only support for "the groin held until 2004", and the date that pinned b in the 2026-10-05 calibration. *Figure 7.*
5. **For the model:** by the decision rule above, this points to an earlier failure. The schedule refit that followed (`../../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`) found that, scored on CoastSat, the best groin is weak (b × f ≈ 0.32–0.36), and that failing from 1996 fits slightly best and matches the repair record. Its post-failure strength equals the current pin's, so test and forward runs would not change. The choice of schedule is still open.

**Caveats.** CoastSat before 1999 rests on 5–13 Landsat images a year per transect; a photo is one day while a CoastSat value is a year's mean; the photo series is aligned to CoastSat because their datums differ, so only its shape is compared.

```
coastsat_gap_across_groins.py   analysis -> tables/
gap_across_groins_figures.py    figures from the tables -> figures/
tables/                         every number quoted below
figures/                        seven figures in reading order; captions in figures/supporting/CAPTIONS.md
```

```
python coastsat_gap_across_groins.py   # ~5 min (bootstraps)
python gap_across_groins_figures.py
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
| 1 | `gap_across_groins_1_transect_map` | which transects make each gap, and where is Buxton? | 6 + 13 transects for the domain gap, 4 + 12 within 500 m. The southernmost GIS 6 transect lies just south of the southern groin; the near-field gap has no such transect |
| 2 | `gap_across_groins_2_side_positions` | which side of the groins moved? | the 1991–95 widening is mostly the south side retreating (GIS 5 about −98 m against GIS 6 about −24 m, with Gordon in 1994). After 2017 the south side gains fill sand |
| 3 | `gap_across_groins_3_gap_and_break` | when did the gap stop widening? | 1995 in both gaps (90% ranges 1992–98 and 1991–2000); the photos follow the same shape |
| 4 | `gap_across_groins_4_break_year` | which break year fits best, and how certain is it? | 1995; 2002–05, the model's failure timing, holds 1–2% of the bootstraps |
| 5 | `gap_across_groins_5_era_rates` | was the gap widening or narrowing in each era? | widening before the repair in both sources; after 2004 CoastSat is flat while the photos fall |
| 6 | `gap_across_groins_6_robustness` | does the break hold under other fits, and is there a separate 2004 drop? | the break stays at 1994–96; no 2004 step is clearly different from zero |
| 7 | `gap_across_groins_7_photos_vs_coastsat` | do the photos and CoastSat agree? | within 13 m before the fill, except the 2004 survey at +42 m |

The failure-schedule figure that was figure 8 is a model setting, not an observation. Since 2026-10-08 it is `schedule_refit_0_failure_schedules.png` in `../../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/figures/`.

## Result

| | GIS 6 minus GIS 5 | 500 m either side |
|---|---|---|
| break year (90% bootstrap) | **1995 (1992–1998)** | **1995 (1991–2000)** |
| bootstrap share at 2002–2005 / ≤ 1997 | 1% / 91% | 2% / 86% |
| rate before / after the break, m/yr | +4.2 / −1.4 | +3.3 / −2.4 |
| gap LRR 1984–95 / 1996–2003 / 2004–16, m/yr | +4.5 ± 1.1 / −2.2 ± 1.8 / +0.2 ± 0.6 | +4.2 ± 1.5 / −2.8 ± 2.5 / −0.5 ± 0.6 |
| separate level step at 2004 | −14 ± 9 m | −14 ± 12 m |
| BIC: break vs break + 2004 step vs two breaks | 170.6 / 171.4 / 172.7 | 187.7 / 189.6 / 188.8 |

Photo-gap rates on the same eras: +4.9 ± 1.5 m/yr (1985–95), no rate for 1996–2003 (only two surveys), and −3.3 ± 0.8 m/yr (2004–16).

**Reading.**
1. **The gap stopped widening at the 1995 last repair, not at Isabel.** Both gap definitions put the break at 1995. It stays at 1994–96 under every robustness check, and the 2003–04 timing gets 1–2% of the bootstrap mass. On the decision rule, this points to an earlier failure.
2. **What follows is a slow decline, not a collapse.** The 1996–2003 rate is not distinguishable from zero. The 2004 step is less than two standard errors in every fit. Adding the step lowers the BIC in only one of the eight fits (GIS 6 minus GIS 5, 1995 left out), and then by 0.1.
3. **The 1995 peak is partly the south side being eroded.** The gap's rise in 1991–95 is mainly GIS 5 retreating, so the peak may be partly storm damage downdrift (Gordon) rather than trapping. Both CoastSat and the photos show it.
4. **The photos and CoastSat agree, except at 2004.** The photos' "held to 2004" rests on the single 2004 survey, which sits 42 m above CoastSat's 2004 mean. That is the largest pre-fill difference (figure 7). After 2004 the photos fall 3.3 m/yr while CoastSat is flat, and that photo decline starts from the high 2004 point.
5. **For the model:** the 2026-10-05 calibration fit says b is "well pinned down by the 2004 date", which is the date CoastSat contradicts. The refit should score b and f against the annual CoastSat gap and try schedules that start weakening in 1996 (done: `../../3-hindcast/2-blocking-1996-2025/2026-10-08-schedule-refit/`).

**Caveats.**
- CoastSat before 1999 rests on 5–13 Landsat images a year per transect.
- A photo is one day, while a CoastSat value is a year's mean.
- The photo series is shifted onto CoastSat over the shared years, because their datums differ, so only its shape is compared.

## Last repair date (checked 2026-10-08)

CSE 2013 (Coastal Science & Engineering, *Shoreline Erosion Assessment & Plan for Beach Restoration, Rodanthe & Buxton (NC)*, report 2403-PHASE1-FR, Nov 2013, p. 20) gives the last repair as **1995**: 184 ft of steel sheet piling on the south groin, after Hurricane Gordon in 1994. This comes from the page-cited event list in the deleted `HAT_groin_zone_investigation.py` (in git at a229c8ff). The report PDF itself is not on this machine, so it wasn't re-read. No source in the repository gives 1996. "1996" appears to be the old ramp's onset (1969 + 27 years delay), read as a repair date.

Fixed in `hard-structures/groin/GROIN_PLAN.md` (now `../structure_history.md`), `scripts/hatteras_ms/HAT_hindcast_methods.md` and `scripts/site_layer/README.md`.

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
| `break_profile.csv`, `break_bootstrap.csv`, `hinge_fit.csv` | the curves behind figures 3–4 |
| `era_rates.csv` | gap rate ± SE per era, CoastSat and photos |
| `photo_gap.csv` | the photo gap, shifted, and its difference from CoastSat |
| `transects_used.csv` | transect ids, side and UTM ends |
