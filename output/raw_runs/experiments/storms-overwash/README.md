# storms-overwash

The storm series and how the modelled overwash is scored against the imagery. Split out of `storms-overwash-and-dunes/` on 2026-10-20 (`../RENAMES.md`); the dune studies are in [`../dunes/`](../dunes/README.md). Not to be confused with the pre-10-05 topic `storms-and-overwash/`, which is the old name of `storms-overwash-and-dunes/` (on D:).

## Studies, oldest first

| study | question | answer | status |
|---|---|---|---|
| `2026-10-13-overwash-scoring-tolerance` | How much do the overwash scores depend on how strictly model and imagery are matched? | Allowing one image either side lifts calibration skill 0.59 → 0.71; one domain either side to 0.81. The test period stays low (0.12-0.26). | record |
| `2026-10-13-soundside-overwash` | How much of the observed overwash could not have come from the ocean, and so is likely sound-side? | 24% (1996-2009) and 45% (2009-2025) of observed overwash had ocean run-up below the lowest surveyed crest; the model hits none of these. | record |
| `2026-10-19-false-alarm-locations` | Where does slower dune growth add false alarms? | All 86 added false alarms are in 2005-2009, after Isabel. 54 of them are at real low spots that stayed dry. | record |
| `2026-10-19-isabel-pattern-scores` | Does Isabel's overwash pattern pick the growth rate? | No: every ceiling ≤ 5.5 m gives the same Isabel (45-47 of 48 hits), so the dune record sets the growth rate. | record |
| `2026-10-20-storm-notebook-rule` | Does running the storms by the notebook's own rule (`v3_240`: no split, no trim) change the skill, against `v3_split12_trim24`? | `v3_240` adopted 2026-10-20; scores in `tables/scores.csv`. | **current**; the runner's storms |

**Status**: **current** means its answer is in use now. **record** means a finished check, kept so the number can be traced.

The 10-13 and 10-19 studies ran on the old storms and dune settings. Their model states were deleted on 2026-10-20; the tables and logs are kept. Studies from before 2026-10-05 are on `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\storms-overwash-and-dunes\`.

Back to [the map](../README.md).
