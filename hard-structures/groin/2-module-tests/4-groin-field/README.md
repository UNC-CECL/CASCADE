# 4-groin-field — several groins on one coast

**Question.** With several groins along one coast, which version of the groin module behaves most like a groin field?

**Setup.** The straight-coast emulator from `../1-straight-coast/` (imported, unchanged): BRIE's alongshore step alone under the option-A waves, no storms, no Barrier3D, no sea level. Here there are 161 domains (80 km) instead of 61, with four groins 2, 5, 10 or 20 domains apart (1 to 10 km), and a single groin on the same coast as the reference. Each layout runs for 50 years at seven orientations, −30° to +30°. With the single groin on the 61-domain rig, `run_field()` reproduces `straight_coast_test.run()` to 0.0 m for all four modules.

Four of the six straight-coast modules, chosen 2026-10-09: the dipole as the reference, and the three trapping versions still in contention (pinned, drift, corrected). The conserving and hybrid versions were dropped after the straight-coast test. Settings are those of the straight-coast test: dipole M 12 m/yr, blocking fraction 0.6. Where two groins' correction zones overlap (spacing under 10 domains), the corrected version applies its leak correction once per cell, because the correction is a property of the cell, not of the groin.

```
groin_field_test.py      the runs, the scores and both figures (about 5 s)
tables/                  groin_field_runs.csv (every run), groin_field_scores.csv (the field runs with each test),
                         groin_field_scorecard.csv (passes per module), groin_field_profiles.csv (final shorelines)
figures/                 groin_field_profiles (each spacing row zoomed to its field), groin_field_scorecard;
                         PDFs and CAPTIONS.md in supporting/
```

    python groin_field_test.py [--years 50]

## The four tests

The tests were set on 2026-10-09, before the runs. Each is scored at year 50 for every orientation and spacing (28 cases per module). Every orientation run has a drift of at least a tenth of the grid-aligned drift, so all are scored.

| test | passes when |
|---|---|
| conserves sand | the net change summed over the coast is at most 10% of the gross change |
| holds sand at every groin | at each of the four groins, the cell just updrift is more than 1 m seaward of the cell just downdrift: a sawtooth |
| cuts the coast below the field | some cell among the five just downdrift of the last groin has retreated more than 1 m |
| close groins act as one | each downdrift groin's offset, as a share of an isolated groin's, is smaller at the closest spacing than at the widest, and at least 0.5 at the widest; counted only where the isolated groin traps updrift; one result per orientation |

**Two changes after the first run, both on 2026-10-09.**
- *A logic fix.* The "act as one" ratio was first counted wherever the isolated groin's offset was larger than 1 m in either direction. That let a groin trapping on the wrong side pass, a negative over a negative: the pinned version passed at +10° to +30°. The ratio now counts only where the isolated groin holds sand updrift. No threshold changed.
- *Wider spacing, asked for to make the effects easier to see.* The first run had the groins 1, 2, 5 and 10 domains apart on 121 domains. At 1 domain apart there is only one cell between groins, so the drift and corrected versions read the approach angle across the next groin's step. That pushed the ratio back toward 1, and the corrected version failed "act as one" at −20°, −10°, +20° and +30°. The first run's scorecard (conserves / every groin / lee cut / act as one, of 28) was: corrected 25 / 27 / 28 / 12, dipole 1 / 26 / 22 / 28, drift 1 / 24 / 22 / 24, pinned 25 / 12 / 12 / 12. With the spacings at 2 to 20 the test compares 2 apart with 20 apart.

## Findings

| module | conserves | every groin | lee cut | act as one | total |
|---|---|---|---|---|---|
| corrected | 25 | 28 | 28 | 24 | **105 / 112** |
| dipole | 0 | 27 | 20 | 28 | 75 |
| drift | 0 | 25 | 20 | 28 | 73 |
| pinned | 24 | 12 | 12 | 12 | 60 |

1. **The corrected version behaves most like a groin field.** It builds a sawtooth and cuts the coast below the field in every case. It conserves sand everywhere except +20° at 2 to 10 apart, where so little sand moves that a net change of 200 m of summed shoreline is 29% of the gross. Elsewhere its largest net change is 135 m against a gross change in the thousands. It fails "act as one" only at +30°, where the drift is a tenth of its grid-aligned value.
2. **The pinned version still traps on the side set by the tilt.** It conserves sand (largest net change 125 m), but four of its 24 passes are at 0°, where it moves nothing. It builds a sawtooth only at −30° to −10°. From +10° on it builds out below each groin instead of above it, the straight-coast finding 2 repeated at every groin.
3. **The dipole and drift versions do not conserve sand in a field.** The drift version creates sand on coasts turned against the drift and destroys it on coasts turned with it: up to 18,800 m of summed shoreline at +20° with the groins 10 apart, where the whole field erodes (panels i and l of `groin_field_profiles.png`). The dipole adds or removes up to 7,000 m. Both build a sawtooth because each groin's own transfer forces one, not because the field holds sand.
4. **Close groins share their sand only weakly in every version.** In a real field the compartments between close groins are short, so each groin holds a small offset and the field acts as one wider structure. Here, at 1 km apart, each of the corrected version's downdrift groins still holds 71% to 92% of an isolated groin's offset at −30° to +20° (45% at +10°). By 5 km apart it is 80% to 101%. None of the modules has a groin length, so nothing caps the offset, and at 0.5 km apart (the first run) each groin held as much as an isolated one at −20°, −10° and +20°. The test passes, but the sharing is much weaker than a real field's.

**What this means for the module.** If a model run has more than one groin, the corrected version is the one to use, with the same reservation as for one groin: its leak correction also acts on any natural bend near a groin (`../3-real-planform/corrected_groin/`). A groin length, a cap on the offset a groin can hold, would be needed before close groins share sand as a real field does. That is a module change, so it is not part of this test.
