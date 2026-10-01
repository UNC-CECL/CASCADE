# HAT-groin-buxton-input - what the groin study is fed

Inputs for the Buxton groin study (`../README.md`): the 1967 island offset for
domains D2-D12, the storm files for the 1967-1997 test and the 1967-2017
hindcast, the observed dune-line target the groin must reproduce, and a
geometric check on the shoreline positions. Products sit beside the code, as
everywhere in this study (`../README.md`, "Where it departs").

```
groin_init/                       products, plus one repair tool
    fix_cascade_yaml.py               repair a parameter YAML with numpy scalars in it
    island_offset/                    1967 D2-D12 offsets (input geojsons, raw, padded)
    storms/                           input_storms/, 1967_1997/, 1967_2017/, 1967_2024/
    target/                           HAT_target_1967_1997_* (the calibration target)
input_prep/                       the scripts that make them
    HAT_target_shoreline_change.py    observed dune-line change and rate -> groin_init/target/
    island_offset/
        island_offset_hybrid_1967.py  1967 offset, padded to 41 -> see its README
    shoreline_position/
        HAT_geometric_distance_sanity_check.py  distance to datum + by-eye figures
                                      -> ../HAT-groin-buxton-output/shoreline_position_output/
    storms/
        HAT_resample_grointest_storms.py  1984-2004 storms resampled to 30 yr -> storms/1967_1997/
        HAT_build_1967_2017_storms.py     stitched 1967-2017 storms -> storms/1967_2017/
```

Run order, where it matters: `HAT_resample_grointest_storms.py` before
`HAT_build_1967_2017_storms.py` (the build reads the resample from
`storms/input_storms/`, a copy of the `1967_1997/` product).

Every folder that holds a script has its own README with the details.
