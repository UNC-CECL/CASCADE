# hard-structures — engineered structures in CASCADE

One folder per structure. Each is laid out like `groin/`: numbered by the kind
of claim it holds (1 observations, 2 module tests, 3 hindcast, 4 forecast), with
code beside its data. Model output goes under
`output/raw_runs/experiments/<structure>/`.

| folder | structure | status |
|---|---|---|
| `groin/` | the Buxton groin field | fitted; the wave-driven groin (`GroinManager`) is in the runner |
| `revetment/` | a shore-parallel wall | simple version under test on a straight coast (2026-10-12) |
| `breakwater/` | a detached offshore breakwater | not started; notes only |

All three are, or will be, CASCADE modules (like `cascade.groin.GroinManager`)
or schedules: BRIE and Barrier3D are never edited.
