# Storm trim length on the adopted setup (2026-09-28)

**Question (Hannah: "why do we limit to 24 hr?").** trim24 was chosen on the old dunes (held near 3 m MHW), and 24 h was the shortest length tried. Barrier3D applies a storm's peak Rhigh, and the gap discharge it drives, for every hour of the storm, so the length controls how much sand overwash moves. Does the choice hold on the setup being adopted?

**Setup** (10 runs, all ran to the end on `hatteras/adopted@2b8f167`; driver `scripts/hatteras_ms/experiments/HAT_trim_length_adopted.py`):

- **Barrier3D:** the local branch `hatteras/adopted`, which has the three overwash fixes plus per-cell dune ceilings. The ceilings are switched on in-process (`DuneCeilingFromStart`, floor 0.5 m).
- **Storms:** every event kept, trimmed to 12, 24, 48 or 72 h, or full length. The 12 h files come from the builder (`--long-events trim --max-duration 12`) and are stored here; 24 h is the adopted `hindcast_storms/*_v3_trim24`; 48 h, 72 h and full are the storm-length selection's files.
- **Runs:** managed, both windows, at the site config's LOESS-7 end rates.
- **Scoring:** as before. Overwash against the imagery is scored under the same Barrier3D, and RMSE is against LOESS-7 via `run_registry.skill_vs_target`.
- **Interruption:** the 2010–2024 runs were stopped and re-run from scratch at Hannah's request. The partial folders were removed first.

| window | trim | storm-hours | PSS | POD | POFD | timing r | space r | RMSE | bias | overwash total (m³/m) | Irene low / rest |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1996–2010 | 12 h | 1,717 | 0.58 | 0.77 | 0.19 | 0.83 | 0.62 | 1.18 | +0.03 | 756 | |
| | **24 h** | 2,839 | **0.60** | 0.78 | 0.19 | **0.84** | 0.63 | **1.16** | −0.03 | 1,695 | |
| | 48 h | 3,832 | 0.59 | 0.78 | 0.20 | 0.82 | 0.62 | **1.16** | −0.14 | 3,127 | |
| | 72 h | 4,186 | 0.58 | 0.78 | 0.20 | 0.82 | 0.65 | 1.21 | −0.22 | 4,289 | |
| | full | 4,613 | 0.58 | 0.78 | 0.20 | 0.82 | 0.64 | 1.30 | −0.39 | 6,138 | |
| 2010–2024 | 12 h | 2,049 | **0.17** | 0.57 | 0.40 | 0.40 | 0.19 | **2.03** | −1.32 | 1,840 | 29/25 · 16/47 |
| | **24 h** | 3,345 | **0.17** | 0.57 | 0.41 | 0.39 | 0.21 | 2.11 | −1.44 | 3,614 | 29/25 · 16/47 |
| | 48 h | 4,730 | 0.15 | 0.57 | 0.42 | 0.36 | 0.21 | 2.27 | −1.66 | 6,789 | 28/25 · 16/47 |
| | 72 h | 5,406 | 0.14 | 0.56 | 0.42 | 0.35 | 0.21 | 2.38 | −1.78 | 8,420 | 28/25 · 16/47 |
| | full | 5,992 | 0.11 | 0.53 | 0.42 | 0.34 | 0.19 | 2.54 | −1.95 | 10,192 | 28/25 · 16/47 |

The Irene columns are modelled / observed domains: the low-dune third (lidar crest < 4.33 m), then the rest.

**What it shows.**

1. **Trim length does not change where or when overwash happens.** PSS, POD, POFD and the timing and space correlations hardly move between 12 h and full length, and Irene is identical at every length. Whether a domain overwashes depends on storm peaks and dune heights.
2. **It changes how much sand moves, and so the shoreline.** Total overwash grows 6–8× from 12 h to full length. The interior bias moves landward steadily (1996: +0.03 to −0.39 m/yr; 2010: −1.32 to −1.95), and RMSE rises with it.
3. **12 h and 24 h are effectively tied.**
   - 1996–2010 favours 24 h: PSS 0.60 against 0.58, RMSE 1.16 against 1.18.
   - 2010–2024 favours 12 h on RMSE (2.03 against 2.11), with PSS equal.
   - The 2010–2024 differences are smaller than what re-solving the end rates will shift, and that window's RMSE is dominated by its −1.3 to −1.4 m/yr bias (see the 2021 CoastSat step).
   - Longer trims are worse in both windows.
4. **Against the per-cell ceiling alone** (trim24, `../2026-09-28-dune-ceiling-per-domain/`), adding the three overwash fixes lifts 2010–2024 PSS from 0.14 to 0.17, POD from 0.49 to 0.57 and Irene's higher-dune domains from 10 to 16 of 47, at POFD 0.34 → 0.41 and RMSE 2.09 → 2.11. In 1996–2010, PSS goes from 0.57 to 0.60 and RMSE from 1.17 to 1.16.

**Reading.** Keep the trim short. 24 h is supported: best in 1996–2010, tied in 2010–2024. Physically it is about two tidal cycles around the storm peak. The imagery records whether washover occurred, not how much, so between 12 h and 24 h the evidence is the shoreline fit, and that is a draw.

**Status: record.** It confirms `v3_trim24` for the adoption; nothing else changed.
