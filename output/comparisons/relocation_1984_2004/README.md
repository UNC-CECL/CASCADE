# relocation_1984_2004 — does CASCADE relocate NC-12 where and when history did?

Written by `scripts/hatteras_ms/HAT_relocation_comparison.py`: arm A (the
roadway module deciding on its own) against arm B (the recorded 1989 and 1999
relocations prescribed), one set per pair of runs. This folder is gitignored
except the READMEs; the numbers that matter are quoted in
`scripts/hatteras_ms/RELOCATION_COMPARISON_RESULTS.md`, and every set can be
regenerated from the runs it names.

## Layout (2026-09-09): the dune-topo version first

```
<version>/                    the dune-topo version BOTH arms were run on, read from the
                              runs' metadata by the script, never typed
    <preset>/                 groin off (the arm names carry `nogroin`)
    <preset>_groin/           groin on
    dune_position_check/      HAT_relocation_dune_position_check.py, one PNG per preset
```

Each set holds the same files: `report.txt` (the console output, headed by the
identity of the two runs — name, date, topo version, commit), `confusion.csv`
(recall and false positives at each tolerance), `first_relocation_year.csv`,
`near_miss_margin.csv`, `setback_by_year.csv`, `setback_summary.csv`,
`indexing_check.csv`, `road_outcomes.csv`, and six GIFs (road relocation and
road-on-topography, for the whole island and for each event block).

| version | what it is | sets here | runs |
|---|---|---|---|
| `v1/` | the original pick set (2026-08-27); **superseded** as CURRENT by v2 on 2026-09-04 | zeroBE, edgeBE, calibBE, each with and without the groin; the dune-position check | the calibration tree, 2026-09-01 |
| `v2/` | the re-pick base, CURRENT | calibBE_groin | `output/raw_runs/version-pair/v2/`, 2026-09-09 |
| `v3/` | the 1984 reconstruction built from v2 (rows behind the road, copy fill, 1984 setbacks) | calibBE_groin | `output/raw_runs/version-pair/v3/`, 2026-09-09 |

The v2 and v3 sets are the like-for-like pair: same scenario, same code, same
day, nothing but the topography and its setback CSV different. The broader
v2-against-v3 figure (setbacks, relocations, island geometry through time) is
in the reconstruction's result step,
`data/hatteras_init/1-barrier3d-domains/1984-start/2-domain-reconstruction-1984/6-result/`.

## Regenerate

```
python scripts/hatteras_ms/HAT_relocation_comparison.py --preset <preset>          # nogroin pair from the calibration tree
python scripts/hatteras_ms/HAT_relocation_comparison.py --preset calibBE \
    --arm-a output/raw_runs/version-pair/v3/1984_2004/calibBE/HAT_1984_2004_calibBE_road_bdm_groin \
    --arm-b output/raw_runs/version-pair/v3/1984_2004/calibBE/HAT_1984_2004_calibBE_road_reloc_bdm_groin
```

The output folder follows from the arms: `<version>/<preset>[_groin]/`. Pass
`--out` only to put a set somewhere else on purpose. A pair of arms on
different versions is refused.
