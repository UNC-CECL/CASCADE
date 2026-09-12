# parameter_backups — NOT the file the model reads

The live parameter file is one level up:

```
data/hatteras_init/Hatteras-CASCADE-parameters.yaml     <- the model reads THIS
data/hatteras_init/reference_yaml_hatteras.yaml         <- annotated reference,
                                                          read by nobody
```

`scripts/cascade_pipeline/run_registry.py` names that path, and it is tracked.
Everything in this folder is a backup and nothing reads it.

| File | What it is |
|---|---|
| `Hatteras-CASCADE-parameters.yaml.older` | an earlier, longer version of the live file |
| `Hatteras-CASCADE-parameters.yaml.bak` | the same vintage, kept under a second name |
| `Hatteras-CASCADE-parameters.yaml.corrupt.bak` | kept deliberately, recording a file that failed to parse |

## Why this folder exists

Until 2026-09-12 these sat in two sibling folders named `yaml/` and
`yaml_backup/`. Between them they held no live file: `yaml/` contained only a
backup, and `yaml_backup/` held an older copy plus the corrupt one. Two folders
whose names suggested authority, neither of which the model reads, directly
beside the file it does.

The backups are **longer** than the live file (40 lines against 24), so the
difference is not a truncation of these — the live file is the shorter one by
intent. Compare before assuming either is stale.
