# hard-structures — engineered structures in CASCADE

One folder per structure. Each is laid out like `groin/`: numbered by the kind
of claim it holds (1 observations, 2 module tests, 3 hindcast, 4 forecast), with
code beside its data. Model output goes under
`output/raw_runs/experiments/<structure>/`.

| folder | structure | status |
|---|---|---|
| `groin/` | the Buxton groin field | fitted; the blocking groin is pinned in the runner |
| `revetment/` | a shore-parallel wall | simple version under test on a straight coast (2026-10-12) |
| `breakwater/` | a detached offshore breakwater | not started; notes only |

All three are, or will be, CASCADE callbacks or schedules: BRIE and Barrier3D
are never edited.
