# Buxton groin — the observed history

What the structure is and what the shoreline did around it. None of this depends on the model. Split from `GROIN_PLAN.md` (reference note of 2026-08-24) on 2026-10-08. The dipole fit that note also held is now `../3-hindcast/1-dipole-1967-2017/dipole_fit_notes.md`.

**Structure history:** installed **1969** · last repaired **1995** (south groin, 184 ft of steel sheet piling after Hurricane Gordon 1994; CSE 2013, report 2403-PHASE1-FR, p. 20) · storm damage **2003**. Until 2026-10-08 this line said 1996; no source gave 1996. The old ramp's onset of 1996 is the first year after the repair. The field spans northing 3901373–3901789, entirely inside GIS 6, so the groin is the face between GIS 5 (downdrift, south) and GIS 6 (updrift, north).

**The one thing to keep straight:** the "fillet" is the GAP between the updrift and downdrift shorelines, not a volume of new beach. Both sides eroded; the sheltered side eroded less. By 2004 the downdrift side had retreated 175 m since 1967 and the updrift side only 25 m. That 150 m difference is the fillet. The post-2004 collapse is ~85% the updrift side eroding once the structure failed, not impounded sand draining downdrift.

## What the wet/dry surveys show

Fillet = shoreline offset `GIS 5 − GIS 6` against a fixed **1967 datum**, from 24 dated wet/dry surveys (`wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv`). Landward-positive: a rising curve means the updrift side is holding while the downdrift side retreats, which is what a groin builds.

| phase | window | fillet | rate |
|---|---|---|---|
| **BUILD** | 1967 → 1978 | 0 → 117 m | **+10.7 m/yr** |
| slow growth | 1978 → 2004 | 117 → 150 m | +1.3 m/yr |
| **RELEASE** | 2004 → 2023 | 150 → 74 m | **−4.0 m/yr** |

The fillet peaks in **2004** and declines from there. The turning point matches the storm, not anything in the model. The widest gap in the record is the 1995 survey at 155 m, a single-survey spike driven by the downdrift side; the 2004 value is 150 m.

### The 1984–2004 and 2004–2023 windows carry opposite signals

| window | fillet change | rate | groin is |
|---|---|---|---|
| 1985 → 2004 | **+52.0 m** | +2.74 m/yr | **still trapping** |
| 2004 → 2023 | **−76.4 m** | −4.02 m/yr | **releasing** |

When the gap stopped widening, from CoastSat rather than the photos, is in `gap_across_groins/`.

`figures/groin_two_shorelines.png` is the one-picture version: the two shorelines and the gap between them.
