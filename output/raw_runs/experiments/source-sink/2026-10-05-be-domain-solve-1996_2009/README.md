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

**Status.** BE set 1, stored 2026-10-06 as the `domainBE` preset (`HATTERAS_BE_RATES_DOMAIN` in the site config, 1996 and 2009). The matrix run `matrix/1996_2009/domainBE/HAT_1996_2009_domainBE_offsetmetres_road_bdm_groinblock` reproduces step 10 bit for bit.

## Option b: the step-10 field smoothed once (2026-10-05)

`HAT_be_domain_solve_net_change.py smooth --step 10`: 7-domain LOWESS over GIS 2–89 with **no robust passes** (`it=0`); GIS 2–10 and both end rates are kept. The first try used statsmodels' default robust reweighting (`it=3`). That treated GIS 80–81 as outliers and passed GIS 82's −4.85 through untouched, making the field more extreme, not less (RMSE 8.2 m). It is kept in `superseded_robust_lowess/` as the record.

| field | sd (GIS 2–89) | range | interior bias | interior RMSE | worst miss |
|---|---|---|---|---|---|
| step 10, solved | 1.11 m/yr | −4.85 to +12.43 | −0.06 m | 0.49 m | 2.0 m (GIS 11) |
| step 10, LOWESS-7 | 0.85 m/yr | −1.93 to +12.43 | −0.61 m | **4.33 m** | 19.6 m (GIS 82) |

**Reading.** Smoothing costs about 9× the RMSE, and only partly at Rodanthe: without GIS 76–88 the RMSE is still 3.4 m. Away from Rodanthe, the misses are the observed erosion and accretion bands themselves (GIS 12–13, 17–18, 28–30, 35–36, 64–68). BRIE's alongshore diffusion spreads the response to a rate, so the field has to be sharper than the target it produces. At Rodanthe, GIS 82 (−19.6 m) and GIS 85 (−11.5 m) need their extra sink. Those cells respond differently from their neighbours inside the model, so the sawtooth cancels real model-internal structure there rather than noise.
