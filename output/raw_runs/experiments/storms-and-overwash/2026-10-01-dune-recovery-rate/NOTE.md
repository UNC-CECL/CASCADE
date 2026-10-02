# Faster dune recovery from flat (2026-10-01)

**Question.** The dune-recovery diagnosis (`../2026-10-01-dune-recovery-diagnosis/`) found that a flattened dune cell restarts at 7.5 cm and needs 8–9 years to clear a typical year's storm. Does a dune that recovers a fixed amount each year improve the overwash skill?

The neighbour-floored ceiling, meant to be tested alongside, was dropped before running (Hannah: "recovery only"). The 2010 north-end lows are whole low stretches of the 2009 lidar, not narrow breaches.

**Setup** (16 runs, all complete, on `hatteras/adopted@d6546c7`; driver `scripts/hatteras_ms/experiments/HAT_dune_recovery_rate.py`):

- **The change:** Barrier3D's dune growth plus A·max(0, 1 − D/ceiling) per year. A goes to the front dune row and A/3 to the back row. It is applied in-process only.
- **Recovery rates:** A = 0 (control), 0.15, 0.30 and 0.50 m/yr.
- **Runs:** managed and natural, both windows, edgeBE, site-config ends, the adopted storm series.
- **Scoring:** storms_vs_overwash's rule.

**The control reproduces the adopted matrix exactly**: skill 0.572 / 0.262, false alarms 123 / 132.

**Results** (`tables/scores.csv`):

| scenario | window | A (m/yr) | hit rate | false-alarm rate | skill | false alarms (Cape / villages / road) | dune cell-years flat | RMSE (m/yr) | bias (m/yr) |
|---|---|---|---|---|---|---|---|---|---|
| managed | 1996–2010 | **0** | 0.74 | 0.17 | 0.57 | 123 (25 / 44 / 54) | 3.8% | 1.17 | +0.06 |
| | | 0.15 | 0.74 | 0.17 | 0.57 | 126 | 3.4% | 1.17 | +0.06 |
| | | 0.30 | 0.77 | 0.17 | 0.59 | 128 | 3.1% | 1.18 | +0.06 |
| | | 0.50 | 0.75 | 0.17 | 0.59 | 121 (22 / 44 / 55) | 2.4% | 1.18 | +0.06 |
| managed | 2010–2024 | **0** | 0.51 | 0.25 | 0.26 | 132 (16 / 37 / 79) | 8.0% | 2.07 | −1.37 |
| | | 0.15 | 0.50 | 0.25 | 0.26 | 129 | 7.7% | 2.06 | −1.36 |
| | | 0.30 | 0.50 | 0.25 | 0.26 | 130 | 7.3% | 2.05 | −1.35 |
| | | 0.50 | 0.50 | 0.24 | 0.26 | 126 (16 / 35 / 75) | 6.7% | 2.02 | −1.32 |
| natural | 1996–2010 | **0** | 0.77 | 0.21 | 0.56 | 155 | 8.1% | 1.15 | −0.10 |
| | | 0.50 | 0.77 | 0.20 | 0.57 | 146 | 4.2% | 1.14 | −0.03 |
| natural | 2010–2024 | **0** | 0.57 | 0.27 | 0.29 | 143 | 19.7% | 2.75 | −2.22 |
| | | 0.50 | 0.56 | 0.27 | 0.29 | 139 | 13.7% | 2.44 | −1.96 |

**Findings:**

- **Recovery barely changes the overwash skill.** Skill moves by at most ±0.02 and false alarms by at most 9, even at 0.5 m/yr. Fewer dune cells sit flat (2010 natural 20% → 14%), but the hotspot domains still overwash. Their recovered crests, capped at ceilings of about 2.5–4.4 m MHW, are still overtopped by the larger storms in most image windows. Being flattened explains why the cells sit at the berm. It is not what makes those domains false alarms.
- **The false alarms at Cape Point don't change.** There, progradation wipes the dune rows, and this term doesn't touch that.
- **The shoreline improves in the natural 2010–2024 run** (RMSE 2.75 → 2.44, bias −2.22 → −1.96): less overwash means less retreat. It is unchanged in the managed runs.

**Not adopted.** Recovery is not a lever for overwash skill at these rates.
