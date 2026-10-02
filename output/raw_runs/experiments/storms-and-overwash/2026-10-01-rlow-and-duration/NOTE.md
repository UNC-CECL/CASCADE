# Rlow and duration as Barrier3D expects them (2026-10-01)

**Question (Hannah: "set up the Rlow + duration experiment").** Two storm columns don't mean what Barrier3D assumes:

- **Rlow.** The builder takes min(TWL) over the event. Upstream MSSM (`UNC-CECL/MSSM_RENCI`, `multivariateSeaStorm_NCB.m`) uses the event maximum of hourly TWL − S/2. Barrier3D uses Rlow only to choose inundation routing.
- **Duration.** Barrier3D holds the peak Rexcess for every hour of `duration`. The adopted series trims each event to 24 h around its peak.

**Setup** (16 runs, all complete, on `hatteras/adopted@d6546c7`; driver `scripts/hatteras_ms/experiments/HAT_storm_rlow_duration.py`; details in the experiments README):

- **Storm series, all on split12:**
  - `split12`: the adopted series, rebuilt and checked identical.
  - `rlow`: MSSM Rlow, with the 24 h trim kept.
  - `dureq`: impulse-equivalent duration, sum((TWL − berm)^1.5) / (Rhigh − berm)^1.5, with no cap.
  - `rlow_dureq`: both changes.
- **Runs:** managed and natural, both windows, edgeBE, site-config ends (not re-solved).

**The series** (`tables/series.csv`): dureq cuts storm-hours about in half (1996: 2,626 → 1,278; 2010: 3,542 → 1,869). Median duration drops from 20–22 h to 7–8 h. The longest event is now 45 h (Mar 2018). MSSM Rlow raises Isabel (1.55 → 2.91 m MHW) but lowers the long nor'easters (Sandy 2.66 → 2.10, Mar 2018 2.84 → 2.25).

**Results** (`tables/scores.csv`):

| scenario | window | storms | PSS | POD | POFD | RMSE (m/yr) | bias (m/yr) | overwash (m³/m) | inundation share |
|---|---|---|---|---|---|---|---|---|---|
| managed | 1996–2010 | **split12** | 0.60 | 0.78 | 0.19 | 1.17 | +0.06 | 1,614 | 0.00 |
| | | rlow | 0.60 | 0.78 | 0.19 | 1.18 | +0.07 | 1,511 | 0.00 |
| | | dureq | 0.58 | 0.77 | 0.19 | 1.18 | +0.10 | 986 | 0.00 |
| | | rlow_dureq | 0.58 | 0.77 | 0.19 | 1.19 | +0.11 | 918 | 0.00 |
| managed | 2010–2024 | **split12** | 0.17 | 0.57 | 0.40 | 2.07 | −1.37 | 3,648 | 0.01 |
| | | rlow | 0.18 | 0.57 | 0.39 | 2.06 | −1.40 | 4,028 | 0.01 |
| | | dureq | 0.18 | 0.57 | 0.40 | 2.05 | −1.34 | 3,052 | 0.01 |
| | | rlow_dureq | 0.18 | 0.57 | 0.39 | 2.05 | −1.38 | 3,508 | 0.01 |
| natural | 1996–2010 | **split12** | 0.56 | 0.81 | 0.25 | 1.15 | −0.10 | 3,867 | 0.01 |
| | | rlow | 0.56 | 0.81 | 0.25 | 1.15 | −0.07 | 3,484 | 0.01 |
| | | dureq | 0.56 | 0.81 | 0.25 | 1.14 | +0.01 | 2,213 | 0.01 |
| | | rlow_dureq | 0.56 | 0.81 | 0.25 | 1.15 | +0.02 | 2,095 | 0.01 |
| natural | 2010–2024 | **split12** | 0.12 | 0.61 | 0.49 | 2.75 | −2.22 | 11,851 | 0.03 |
| | | rlow | 0.12 | 0.61 | 0.49 | 2.62 | −2.12 | 10,350 | 0.04 |
| | | dureq | 0.12 | 0.61 | 0.49 | 2.47 | −1.98 | 8,155 | 0.03 |
| | | rlow_dureq | 0.12 | 0.61 | 0.49 | 2.42 | −1.94 | 7,492 | 0.04 |

**Findings:**

- **Where and when overwash happens is unchanged by every variant.** PSS, POD and POFD move by at most 0.02. This is the same lesson as the trim-length and splitting checks: storm peaks against the per-cell crests decide overwash presence.
- **Rlow is not a lever.** Barrier3D routes 0–4% of storm-gaps as inundation under either definition. With crests near 5 m, Rlow rarely clears a gap's mean crest whichever way it is defined. The definition mismatch is real but inert.
- **dureq cuts overwash volume by 25–45%.**
  - Managed shoreline skill is unchanged: RMSE within 0.02 m/yr, and the 1996 bias moves +0.06 → +0.10.
  - Natural 2010–2024 improves: RMSE 2.75 → 2.47 (2.42 with rlow), bias −2.22 → −1.94.
- **Caveat:** the ends were solved on split12 and not re-solved here. The volume change would shift them.

**Not adopted.** Adopting means a builder option, new files in all 4 windows, re-solved ends and a matrix re-run. That decision is Hannah's.
