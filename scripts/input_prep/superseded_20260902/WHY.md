# old_source_sink_search — retired

The scripts for the parametric source/sink search; its results are in the data tree.

**Dated 2026-09-02**, from git: the last change before the 2026-09 reorganisation, which is the best available evidence of when work on it stopped. It is not necessarily the day a decision was taken.

Nothing here is maintained or expected to run. Rule 4 of `ORGANIZATION.md`: retirement is a dated folder with a note, and deletes nothing.

## What is in here

```
parametric_ss/
    HAT_source_sink_bothperiods.py            30.9 kB
    HAT_source_sink_bothperiods_fulldomain.py
    HAT_source_sink_datadriven.py
    old_versions/
        HAT_source_sink_bothperiods.py        25.5 kB  -- NOT the same file
        HAT_source_sink_bothperiods_v4.py
```

Added 2026-09-22, because rule 7 asks a parent note to actually name its
children, and this one did not.

**Mind the duplicate.** `HAT_source_sink_bothperiods.py` exists at both depths
with different contents (30.9 kB against 25.5 kB, different md5). The one in
`old_versions/` is the earlier of the two. Nothing imports either, so this is
a reading hazard rather than a running one, but it is the reason to check
which depth you opened before quoting anything from it.

The `old_versions/` name is left as it was found. It sits inside this dated
folder, so it inherits this note, and renaming a file or folder inside a
retirement edits a record for no gain.
