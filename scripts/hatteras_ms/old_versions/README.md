# old_versions — superseded, kept for the record

| File | What it is |
|---|---|
| `HAT_hindcast_1984_2024.py` | an earlier copy of the runner |
| `HAT_hindcast_1984_2024_newdomains.py` | the same, during the domain rebuild |
| `HAT_hindcast_1984_2024_newplot.py` | the same, during a plotting change |
| `benton_script.py` | **a different lineage** -- the original driver this project inherited, predating the runner. See its header. |

The first three are versions of the current runner. The fourth is not: it runs
105 real domains against today's 90 and builds Cascade directly, so it shares
no geometry and no configuration with anything else here.

Nothing in this folder is maintained or expected to run.
