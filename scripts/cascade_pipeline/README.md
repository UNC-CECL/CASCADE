# cascade_pipeline - the library the runner imports

Not a stage and not a script: this is the shared code behind
`scripts/hatteras_ms/HAT_hindcast_1984_2024.py` and everything that reads its
output.

| Module | Holds |
|---|---|
| `hindcast.py` | building and running CASCADE, the run name, the shoreline target |
| `run_registry.py` | where a run lives, and the run index. **Address runs through this**, never by joining paths |
| `domains.py` | the padded and GIS domain geometry, and the conversions between them |
| `roadway.py`, `nourishment.py` | the management forcing a period carries |
| `coastsat_loess.py` | the observed rate series and its smoothing |
| `annotations.py`, `plotting/` | the geography layer and the figure types |
| `reports.py` | the blocks the runner prints, so a run log states how it was driven |

Changing anything here changes every run made afterwards. The run index records
enough per run - preset values, digest, topography version, git commit - to
tell runs made before a change from runs made after.
