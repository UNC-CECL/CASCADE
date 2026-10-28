# 4-hindcast — the groin module fitted to Buxton

**Code these studies ran against.** The fixed-shift and blocking (`BlockingGroinManager`) groins that the scripts in `1-dipole-1967-2017/` and `2-blocking-1996-2025/` use or patch were removed from `cascade/groin.py` on 2026-10-14, when the wave-driven groin became the runner's only groin. The scripts ran against the code at commit `4d937b79`; check that commit out to rerun one.

The module inside real model runs of Hatteras, scored against `../1-observations/`. Numbered in the order they were done; the first two are superseded.

```
1-dipole-1967-2017/     SUPERSEDED. The fixed-shift groin (the M, f dipole) fitted on a 1967-2017 rig of GIS 2-12
                        and on the 1984-2004 / 2004-2024 hindcast periods. Pinned M 60,
                        f 0.6 on 2026-08-30; replaced 2026-10-05
2-blocking-1996-2025/   SUPERSEDED. The blocking groin on the DEM-to-DEM periods: calibration
                        1996-2009, test 2009-2025. Pinned b 0.6, f 0.6 on 2026-10-05; replaced
                        2026-10-14. Also the drift and conserving groins tried on the way
3-wave-driven-1967-2025/ The wave-driven groin, the runner's groin: the window strengths
                        (0.16 / 0.10) and the whole-life 0.25 from 1967
```

Run output for all three is under `output/raw_runs/experiments/groin-blocking/` and `groin-wave-driven/` (`groin/` on D: and before 2026-10-20) and `output/calibration/`; only scripts, tables, figures and logs are here.
