# 1984-start `dune-topo` — version index

Four versions. **Two are full extractions of all 90 domains; two are small
additive layers on top of one of them.** That is the whole structure:

```
v1  extraction, OLD pick set (2026-08-27)          67 MB   ── superseded, kept
                                                            for the pea1989 arms

v3  extraction, RE-PICK with the road visible      67 MB   ── THE BASE. What
    (2026-09-02). 39 domains got a taller crest;              topo_dirs() resolves
    GIS 85 +1.52 m. No inserted rows.                         to today.
     │
     ├── v4  v3 + rows at the two relocation        4.9 MB  ── superseded by v5;
     │       blocks only (GIS 9-14, 84-87)                     kept for the
     │       8 domains get rows                                v4-vs-v5 scope figure
     │
     └── v5  v3 + rows at EVERY domain whose        4.9 MB  ── THE ONE TO TAKE
             measured 1984→1997 shift rounds                   FORWARD
             to ≥ 1 cell. 38 domains, 98 rows.
```

`v4` and `v5` have **identical N at all ten block domains** (verified
2026-09-03). The only difference between them is scope: `v4` scoped the insert
by *where NC-12 was relocated*, `v5` by *where the measurement says 1984 land is
missing from the 1996 survey*. Under `v4`, 30 domains had N ≥ 1 and were skipped
by a scope decision rather than by a measurement, so "unchanged" meant two
different things along the island.

## Telling them apart at a glance

| | full extraction? | pick set | rows added | domains with rows | own setback CSV |
|---|---|---|---|---|---|
| `v1` | yes, 90 domains | old (08-27) | no | — | no |
| `v3` | yes, 90 domains | **re-pick (09-02)** | no | — | no |
| `v4` | no, layer on v3 | re-pick | yes | 8 | yes |
| `v5` | no, layer on v3 | re-pick | yes | **38** | yes |

A version with inserted rows carries `README.md`, `RUN_MANIFEST.txt`, a
per-domain `HAT_seaward_row_insert_audit.csv`, and a
`RoadSetback_1984_dunestart.csv` matched to *those* arrays. **Mixing a setback
CSV with another version's topography reintroduces exactly the off-by-N error
the inserts remove.**

## Which one loads — and why `CURRENT` currently does nothing

`scripts/hat_topo_version.py` resolves in this order:

1. an explicit `override=` argument
2. the `HAT_TOPO_VERSION_1984_START` environment variable
3. **the extractor's own `VERSION` literal**
4. the `CURRENT` file in this folder
5. the only version present, if there is exactly one

> **`CURRENT` says `v5`, but `resolve_version("1984-start")` returns `v3`.**
> Verified 2026-09-03. Rule 3 fires first: `HAT_dune_topo_extractor.py` has
> `VERSION = "v3"`, and that outranks the file. `CURRENT` was set to `v5` to
> record the intent, but **it is inert and will stay inert** until one of the
> two things below happens.

That double duty is the known weak point of this design: the extractor's
`VERSION` says both *what it writes* and *what everyone reads*. Do **not**
"fix" it by editing the literal to `v5` — the extractor writes to that name, so
a re-run would overwrite v5's arrays with a fresh extraction.

### To actually run on v5

Both steps, or you get a mixed state:

```
set HAT_TOPO_VERSION_1984_START=v5
copy dune-topo\v5\RoadSetback_1984_dunestart.csv ^
     ..\..\4-mgmt-forcing\road_offset\dunestart_offset\1984\
```

`hatteras_site_config.py:142` hardcodes that setback path, which is why the
copy is not optional.

## What was removed, and when

**`v2` — deleted 2026-09-03.** It was `v1 + rows, block scope`: the same insert
as `v4`, but built on the pre-re-pick `v1` extraction. Superseded twice over.
Its `blocksdate` arm was dropped from `HAT_run_crest_experiment.py` and
`HAT_plot_seaward_insert_compare.py`'s default pair moved from `v1` vs `v2` to
`v3` vs `v5` — the same comparison on the current pick set.

**The completed run output survives** at `output/raw_runs/blocksdate` and
`blocksdatenoreloc`. Only re-running it from its inputs is no longer possible.

Sizes and reasons are in `../../archive_purge_20260903.csv`; the lineage entry
is in `../../LINEAGE.md`.

Retired experimental variants live in `../dune-topo-experiments/` and are
deliberately outside this folder — `hat_topo_version.versions()` lists the
directories here, so keeping them in would make every resolution message
unreadable.

## No figures in this folder any more

The fifteen loose PNGs that used to sit here were all insert figures, and all
eight plotters defaulted their output here. On 2026-09-03 five stale ones were
deleted and the rest moved to
`../row-insert-scope/figures/` (two unregenerable ones into its `frozen/`
subfolder). The plotter defaults were repointed via
`hat_topo_version.insert_figures_dir`, so new runs land there too.

This folder now holds versions and this index. Each version keeps its own
`figures/` subfolder, which was never part of the loose pile.

Written 2026-09-03.
