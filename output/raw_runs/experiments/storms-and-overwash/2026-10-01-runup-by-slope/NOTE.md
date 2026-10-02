# Run-up from each domain's beach slope (2026-10-01)

**Question (Hannah: "test the spatially varying runup with beach slope").** Every domain sees the same storm: one gauge, one wave point, one slope (0.06). Does run-up from each domain's own foreshore slope improve the overwash skill?

**Setup** (12 runs, all complete, on `hatteras/adopted@d6546c7`; driver `scripts/hatteras_ms/experiments/HAT_storm_runup_by_slope.py`; details in the experiments README):

- **Slopes:** measured on the 2009 lidar for both windows (`tables/beach_slope_by_domain.csv`), MHW to berm along each row, corrected for the shoreline angle. They are clipped to Stockdon's 0.01–0.16 and smoothed with a 5-domain running median, giving 0.05–0.13 with median 0.075. The 1996 graft DEM was unusable: 19 domains could not be measured, and slopes reached 0.72.
- **Storm series:** the same events and hours as the adopted series, with Rhigh and Rlow recomputed hourly from each domain's slope (`tables/series.csv`).
  - `control`: 0.06 everywhere, identical to the adopted file in all 90 domains.
  - `measured`: each domain's own slope.
  - `pattern`: the measured slopes rescaled to median 0.06.
- **Runs:** managed and natural, both windows, edgeBE, site-config ends (not re-solved).

**The control reproduces the adopted matrix exactly** (skill 0.572 / 0.262).

**Results** (`tables/scores.csv`; outcomes per domain in `tables/outcomes_by_domain.csv`):

| scenario | window | slopes | hit rate | false-alarm rate | skill | hits / misses / false alarms | alongshore r | RMSE (m/yr) | bias (m/yr) |
|---|---|---|---|---|---|---|---|---|---|
| managed | 1996–2010 | **0.06** | 0.74 | 0.17 | **0.57** | 51 / 18 / 123 | 0.65 | 1.17 | +0.06 |
| | | measured | 0.78 | 0.27 | 0.51 | 54 / 15 / 198 | 0.50 | 1.40 | −0.28 |
| | | pattern | 0.71 | 0.20 | 0.51 | 49 / 20 / 146 | 0.51 | 1.23 | −0.01 |
| managed | 2010–2024 | **0.06** | 0.51 | 0.25 | **0.26** | 59 / 56 / 132 | 0.19 | 2.07 | −1.37 |
| | | measured | 0.67 | 0.33 | 0.34 | 77 / 38 / 174 | 0.36 | 2.54 | −1.77 |
| | | pattern | 0.57 | 0.26 | 0.31 | 66 / 49 / 138 | 0.31 | 2.20 | −1.43 |
| natural | 1996–2010 | **0.06** | 0.77 | 0.21 | **0.56** | 53 / 16 / 155 | 0.60 | 1.15 | −0.10 |
| | | measured | 0.81 | 0.29 | 0.52 | 56 / 13 / 214 | 0.55 | 1.36 | −0.42 |
| | | pattern | 0.75 | 0.25 | 0.51 | 52 / 17 / 180 | 0.55 | 1.23 | −0.11 |
| natural | 2010–2024 | **0.06** | 0.57 | 0.27 | **0.29** | 65 / 50 / 143 | 0.24 | 2.75 | −2.22 |
| | | measured | 0.73 | 0.34 | 0.39 | 84 / 31 / 177 | 0.42 | 3.09 | −2.49 |
| | | pattern | 0.63 | 0.28 | 0.34 | 72 / 43 / 149 | 0.37 | 2.67 | −2.13 |

**Findings:**

- **2010–2024: varying run-up helps.** The pattern alone raises skill by +0.05 in both scenarios. The alongshore r of where overwash happens rises from 0.19 to 0.31 (managed) and from 0.24 to 0.37 (natural), with only 6 more false alarms. Most of the gain is misses recovered on the Irene image: in the village zone they fall from 12 to 8 (pattern) or 4 (measured). The measured slopes, which are steeper on average, gain more skill (+0.08 to +0.10), but they add 34–42 false alarms and worsen RMSE by 0.3–0.5 m/yr.
- **1996–2010: it hurts.** Skill falls by 0.05–0.06 and the alongshore r falls from 0.65 to 0.51. The pattern adds 34 false alarms on the steep beaches (slope > 0.07) and recovers no hits. The slopes come from the 2009 lidar, 13 years after this window starts; beach slope is a property of a particular time, so the 1996 beaches may not have had this pattern.
- **Shoreline:** neutral to worse. The pattern changes RMSE by −0.08 to +0.13 m/yr. More overwash means more retreat, and the ends were not re-solved.

**Reading.** Beach slope from the same era as the window carries real spatial signal for overwash (2010–2024). The same pattern 13 years out of date does harm (1996–2010). Whether 1996 has a usable beach slope source is the open question: the ALACE graft is not one.

**Not adopted.** It would mean a CASCADE change (one storm series per domain; CASCADE takes one storm file) and re-solved ends. That decision is Hannah's.
