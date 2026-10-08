# 1-observations/coastsat_groin_condition

What condition was the Buxton groin in when the calibration period starts in 1996? This folder answers that from CoastSat, independently of the wet/dry photo record behind the current model schedule. It's an observational study with no model runs. It started 2026-10-08 (moved here from the DEM-to-DEM fits the same day, because it is an observation, not a model test).

**Why it matters.** The pinned blocking groin (b 0.6, f 0.6) runs at full strength until an instant failure at the 2004 step. That timing comes from the photos, where the gap seemed to hold until 2004. The maintenance record has no repair after 1995.

**Decision rule (agreed 2026-10-08).** Fit a piecewise breakpoint to the gap. If the turn falls at 2003–04, keep full strength until 2004. If it falls around 1996 or earlier, move to an earlier or ramped failure and refit.

```
coastsat_groin_gap.py        analysis -> tables/
groin_condition_figures.py   figures from the tables -> figures/
tables/                      every number quoted below
figures/                     eight figures in reading order; captions in figures/supporting/CAPTIONS.md
```

```
python coastsat_groin_gap.py        # ~5 min (bootstraps)
python groin_condition_figures.py
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

| # | figure | what it shows |
|---|---|---|
| 1 | `groin_condition_1_transect_map` | which transects make each gap. The southernmost GIS 6 transect lies just south of the southern groin; the near-field gap has no such transect |
| 2 | `groin_condition_2_side_positions` | each side on its own. The 1991–95 widening is mostly the south side retreating (GIS 5 about −98 m against GIS 6 about −24 m, with Gordon in 1994). After 2017 the south side gains fill sand |
| 3 | `groin_condition_3_gap_and_break` | the gap with its hinge and break year; the photos overlaid on the domain gap |
| 4 | `groin_condition_4_break_year` | misfit by break year, and the bootstrap distribution of the break; 2002–05 holds 1–2% |
| 5 | `groin_condition_5_era_rates` | the gap rate up to the last repair, from 1996 to Isabel, and from 2004 to the fill, for CoastSat and the photos |
| 6 | `groin_condition_6_robustness` | the break year and the 2004 step under the four checks |
| 7 | `groin_condition_7_photos_vs_coastsat` | photo minus CoastSat in each photo year; the 2004 survey stands out at +42 m |
| 8 | `groin_condition_8_candidate_schedules` | the three failure schedules for the refit |

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
5. **For the model:** the 2026-10-05 calibration fit says b is "well pinned down by the 2004 date", which is the date CoastSat contradicts. The refit should score b and f against the annual CoastSat gap and try schedules that start weakening in 1996 (figure 8).

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
