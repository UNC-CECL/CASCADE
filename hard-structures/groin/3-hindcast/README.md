# 3-hindcast — the groin module fitted to Buxton

The module inside real model runs of Hatteras, scored against `../1-observations/`. Numbered in the order they were done; the first is superseded.

```
1-dipole-1967-2017/     SUPERSEDED. The dipole (M, f) fitted on a 1967-2017 rig of GIS 2-12
                        and on the 1984-2004 / 2004-2024 hindcast periods. Pinned M 60,
                        f 0.6 on 2026-08-30; replaced 2026-10-05
2-blocking-1996-2025/   The blocking groin on the DEM-to-DEM periods: calibration 1996-2009,
                        test 2009-2025. Pinned b 0.6, f 0.6 on 2026-10-05
```

Run output for both is under `output/raw_runs/experiments/groin/` and `output/calibration/`; only scripts, tables, figures and logs are here.
