# 2026-10-05-be-domain-solve-1996_2009

**Question.** What source/sink rate does each of the 90 domains need so that the calibration run (1996 → 2009) ends where CoastSat says the shoreline went? This is BE set 1 of the DEM-to-DEM plan.

**Setup.**
- Driver: `scripts/hatteras_ms/HAT_be_domain_solve_net_change.py solve --period 1996`.
- Start (step 0): the pinned calibration run, `experiments/groin/2026-10-05-blocking-fit-dem-to-dem/b0.60_f0.6/`. It has the solved end rates (+1.4981 / +10.5659), the blocking groin (b 0.6, f 0.6) and zero elsewhere.
- Every pass: full management, groin on, relocations off, current (option A) waves. The whole field is passed through `HAT_BE_OVERRIDE`.
- Residual: the net-change target (`coastsat/net_change/1996_2009`, 7-domain LOWESS, GIS 1–10 raw) minus the model's raw end-minus-start shoreline.
- Update: rate += 0.8 × residual / gain. The gain is 13 m per m/yr (the run length), or 4.4 at GIS 90.
- Stop: interior RMSE < 0.5 m and every domain within 1 m, or 10 passes.

**Result.** It stopped at the 10-pass limit with the interior RMSE converged (0.47–0.49 m from step 8 on, bias ±0.1 m). Five domains still miss by more than 1 m: GIS 10–12 (−1.4 to −2.0 m), GIS 64 and GIS 81 (+1.9 m). GIS 10–12 barely respond to their rate in the last passes. The per-pass record is `solve_log.csv`; fields and residuals are in `fields/` (`be_field_step10.csv` is the answer).

| step | interior bias | interior RMSE | worst miss |
|---|---|---|---|
| 0 | −5.85 m | 18.58 m | 47.8 m (GIS 82) |
| 1 | +1.53 | 3.54 | 12.7 (GIS 81) |
| 4 | −0.10 | 0.97 | 5.0 (GIS 82) |
| 8 | −0.09 | 0.48 | 2.0 (GIS 11) |
| 10 | −0.06 | 0.49 | 2.0 (GIS 11) |

**The field (step 10).** Interior mean −0.27 m/yr, sd 1.11. The range is −4.85 m/yr (GIS 82) to +12.43 m/yr (GIS 90, which rose from the end solve's +10.57 as the interior changed). It follows the observed erosion and accretion bands: −1 to −2 at GIS 7–14 and 21–25, +1.4 to +1.8 at 17–18 and 28–31, −2 to −2.4 at 77–79. **GIS 81–85 alternate sharply** (+1.13, −4.85, −1.41, −0.19, −3.44): the field is cancelling domain-scale structure in the model's own response near Rodanthe, not a smooth observed signal. Summed over the reach (500 m domains × the 18.77 m active profile) it is a net sink of about 90,000 m³/yr.

**Status.** BE set 1, not yet written to the site config.
