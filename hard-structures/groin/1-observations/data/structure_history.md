# Buxton groin — the observed history

What the structure is and what the shoreline did around it. None of this depends on the model. Split from `GROIN_PLAN.md` (reference note of 2026-08-24) on 2026-10-08. The dipole fit that note also held is now `../../3-hindcast/1-dipole-1967-2017/dipole_fit_notes.md`.

**Structure history:** installed **1969** · last repaired **1995** (south groin, 184 ft of steel sheet piling after Hurricane Gordon 1994; CSE 2013, report 2403-PHASE1-FR, p. 20) · storm damage **2003**. Until 2026-10-08 this line said 1996; no source gave 1996. The old ramp's onset of 1996 is the first year after the repair. The field spans northing 3901373–3901789, entirely inside GIS 6, so the groin is the face between GIS 5 (downdrift, south) and GIS 6 (updrift, north).

**After 1995, from press reports (checked 2026-10-09).**
- The U.S. Navy built three groins in 1969 and 1970. The northernmost failed soon after; a December 1972 photo shows its centre section gone. The southernmost was repaired in 1975 and 1994 and extended landward in the early 1980s. (The GIS layer has four digitized lines, two of them at one site.)
- Maintenance largely stopped after the lighthouse was moved in 1999.
- Beach users reported "increasingly significant holes and structural losses beginning around 2015 and 2016". By 2021 "large portions of the structures had been lost" (Buxton Civic Association, quoted by the Island Free Press).
- Today the northern and middle groins are submerged most days. Only the southernmost is more than 50% intact, the threshold for repair; the other two cannot be rebuilt.
- Dare County's manager attributed the faster sand loss after the 2022 fill partly to the southernmost groin's deterioration ("there's no way to prove that one hundred percent").
- The southernmost groin is being rebuilt in 2026 on its old footprint, about 640 ft (195 m) into the ocean, with substantial completion due by 15 December 2026.
- Both the 2017-18 and 2022 fills ran from the Haulover south to the southernmost groin, so both were placed across the whole groin field. No source found says whether they buried or overtopped the groins.

Sources:
- Island Free Press, "Decades of shoreline engineering: the long history of a changing Buxton beach": https://islandfreepress.org/blog/decades-of-shoreline-engineering-the-long-history-of-a-changing-buxton-beach/
- Outer Banks Voice, 2026-09-17: https://www.outerbanksvoice.com/2026/09/17/latest-update-on-buxton-beach-nourishment-and-groin-repair/
- Coastal Review, 2026-05: https://coastalreview.org/2026/05/final-federal-permit-clears-way-for-buxton-groin-repair-project/
- Outer Banks Voice, 2018-03-01: https://outerbanksvoice.com/2018/03/01/delayed-buxton-beach-nourishment-project-is-finally-done/

**The one thing to keep straight:** the "fillet" is the GAP between the updrift and downdrift shorelines, not a volume of new beach. Both sides eroded; the sheltered side eroded less. By 2004 the downdrift side had retreated 175 m since 1967 and the updrift side only 25 m. That 150 m difference is the fillet. The post-2004 collapse is ~85% the updrift side eroding once the structure failed, not impounded sand draining downdrift.

## What the wet/dry surveys show

Fillet = shoreline offset `GIS 5 − GIS 6` against a fixed **1967 datum**, from 24 dated wet/dry surveys (`../analysis/1-wetdry_photo_positions/Change_from_wetdry_1967_D2_D12.csv`). Landward-positive: a rising curve means the updrift side is holding while the downdrift side retreats, which is what a groin builds.

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

When the gap stopped widening, from CoastSat rather than the photos, is in `../analysis/2-gap_across_groins/`.

`../analysis/1-wetdry_photo_positions/figures/groin_two_shorelines.png` is the one-picture version: the two shorelines and the gap between them.
