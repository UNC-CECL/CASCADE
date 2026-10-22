# groin-wave-driven

The wave-driven groin (`cascade.groin.GroinCallback`, once "corrected"): the runner's only groin since 2026-10-14. Its strength per window, its whole-life strength, where its face sits and the shape of its fillet. Split out of `groin/` on 2026-10-20 (`../RENAMES.md`). Write-ups are in `hard-structures/groin/4-hindcast/3-wave-driven-1967-2025/`; the column on the right names the write-up folder each study belongs to.

## Studies, oldest first

| study | question | answer | status | write-up |
|---|---|---|---|---|
| `2026-10-13-pregroin-planform` | How strong is the groin when working, from the 1967 shoreline (before the groins)? | With the next two studies: 0.25 from installation to 2009, 0.20–0.25 to 2016, off from 2017. | **current**; the whole-life strength (0.25) for runs that start before the fillet exists | `2026-10-13-pregroin-planform/` |
| `2026-10-13-pregroin-full-life` | The same strengths run over the groin's full life. | See the row above; also the runs `2026-10-14-face-6-7` and `-fillet-shape` read. | **current** | `2026-10-13-pregroin-planform/` |
| `2026-10-13-continuous-1967-2025` (+ `-seed1`, `-seed2`) | One continuous run 1967–2025, no reset at 1996 or 2009; the seed variants resample the 1968–1983 storms instead of replaying them. | Follows the record to 17 m over 1967–2025. | **current** | `2026-10-13-pregroin-planform/` |
| `2026-10-14-repin-check-experiment`, `-runner` | The in-process experiment at 0.25 against the runner at its defaults, run when the groin moved into the runner. | No write-up; the run folders and logs are the record. | record | none |
| `2026-10-14-window-refit` | What strength does each window need, from its observed shoreline on today's runner? | **0.16** for 1996–2009, **0.10** for 2009–2025 (0.25 overshoots). | **current**; `hat_run.yaml` `groin.strength_1996` / `_2009` | `2026-10-14-window-refit/` |
| `2026-10-14-window-refit-check` | Does the runner at its defaults reproduce the 0.16 and 0.10 refit runs? | Yes, to 0.0 m. | record | `2026-10-14-window-refit/` |
| `2026-10-14-face-6-7` | Where should the face be? The groins all sit in GIS 6. | A field across GIS 6 fits the build-up best (21.7 vs 25.7 m), but the windows prefer the 5\|6 face, which is kept. | record | `2026-10-14-face-6-7/` |
| `2026-10-14-field-windows` | The GIS 6 field (half the sand kept in it) run on the hindcast windows. | The windows prefer the 5\|6 face (see the row above). | record | `2026-10-14-face-6-7/` |
| `2026-10-15-split-sample` | Run unchanged on the other window, does each strength still beat no groin? | Yes: skill 0.59 (0.16 on 2009–2025) and 0.84 (0.10 on 1996–2009). | **current** | `2026-10-15-split-sample/` |
| `2026-10-20-tapered-fillet-mixed-settings` | First attempt at the tapered fillet; the runner's settings changed mid-study. | Its "any groin hurts 2009–2025" result was the settings change, not the groin. | record; don't compare across its runs | `2026-10-20-tapered-fillet/` |
| `2026-10-20-groin-2009-diagnosis` | Why did every strength end 2009–2025 ~45 m landward of no groin at GIS 3–4? | The no-groin run predated commit ecd6f0ff; on one setting the groin's changes stay in GIS 1–9. | record | `2026-10-20-tapered-fillet/` |
| `2026-10-20-tapered-fillet` | Should the fillet taper north (1/2, 1/3, 1/6 in GIS 6, 7, 8)? | Ties on the gap, worse on net change near the groin (34.1 vs 23.1 m in 1996–2009). Runner's layout kept. | record | `2026-10-20-tapered-fillet/` |

**Status**: **current** means its answer is in use now. **record** means a finished check, kept so the number can be traced. **superseded** means a later study replaced the answer.

Every study before 10-20 ran on the old storms (`v3_split12_trim24`), beach slope 0.04 and Dmaxel 5.5 m. The strengths in use (0.16 / 0.10 / 0.25) still come from them and need refitting on the 10-20 settings.

Studies here from before 2026-10-13 were moved on 2026-10-20 to `D:\hatteras_ch1\model_runs\output\raw_runs\experiments\groin-wave-driven\`, with this topic's table as it stood in the README there.

Back to [the map](../README.md).
