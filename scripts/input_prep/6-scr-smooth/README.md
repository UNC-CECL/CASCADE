# 6-scr-smooth - LOESS smoothing of the CoastSat shoreline change rates

Takes the per-transect and per-domain LRR tables that `5-scr` produces and
smooths them along the coast, so the noisy transect signal becomes a curve the
model can be scored against.

    5-scr/CoastSat/<period>/transect_lrr_full.csv     906 rows, one per transect (~50 m spacing)
    5-scr/CoastSat/<period>/domain_lrr_summary.csv    90 rows, one per GIS domain

**Nothing in this stage writes a model input.** It produced a *decision* -
smooth transects first, then average to domains, at a 7-10 domain window - and
the hindcast implements that decision itself in
`scripts/cascade_pipeline/coastsat_loess.py`. No CSV crosses from here into a
run. See "What the hindcast actually runs" below.

## The stage

```
HAT_loess_method_comparison.py   CURRENT. Transect-first vs domain-first smoothing.
HAT_loess_dsas_vs_coastsat.py    Side question: DSAS vs CoastSat as sources.
hindcast_loess_snapshot/         Read-only copy of the module the runs use.
old_smoothing/                   The lineage, v1-v4, provenance only.
```

Both scripts run from anywhere - every path is anchored on the `pyproject.toml`
at the repo root, not typed as an absolute literal.

Every script here is `HAT_loess_<what it compares>`, matching the `HAT_` prefix
used across `0-elevation`, `3-env-forcings`, `4-mgmt-forcings` and
`7-source-sink`. Each writes its products to
`data/hatteras_init/6-scr-smooth/<its own name>_output/` - scripts live under
`scripts/`, products under `data/`, the same split every other stage uses, and
those folders are gitignored because one run regenerates them. The one exception
is the snapshot, which keeps the live module's filename `coastsat_loess.py` on
purpose - it is a mirror, and the name is what makes it one.

Both reconfigure stdout to UTF-8 at import. Their progress lines use arrows and
en-dashes, and a cp1252 Windows console cannot encode those: the DSAS script
used to die in its closing summary *after* writing every figure, which looks
like a failed run that had in fact finished.

### HAT_loess_method_comparison.py

The one that settled the method. Runs **both** smoothers on every call:

- **transect-first** - LOESS over 906 transects using along-coast metres as
  x, then average the smoothed signal down to 90 domains;
- **domain-first** - LOESS over the 90 pre-averaged domain means directly.

`LOESS_WINDOW_KM = 3.5` (7 domains) is the primary window;
`COMPARE_WINDOWS_KM = [2.5, 3.5, 5.0]` (5, 7, 10 domains) is the sensitivity set.
Writes to `data/hatteras_init/6-scr-smooth/HAT_loess_method_comparison_output/`
in four numbered subfolders:
`01_transect_based`, `02_domain_averaged`, `03_cascade_inputs`,
`04_method_comparison`.

The two methods agree across most of the island and separate only where the
gradient is steepest - Rodanthe around domains 76-80 in period 1, Buxton around
7-8 in period 2, and the first few domains. That is the argument for
transect-first: domain-first pre-averages away the structure exactly where it
matters.

Verified against the module the runs use: at a matched 7-domain window this
script and `cascade_pipeline/coastsat_loess.py` agree to 8.9e-16 m/yr. The
decision recorded here and the code implementing it are the same method, not
two lookalikes.

Its transect-space figures carry the same geographic marks as the domain-space
ones. The `T_*` positions are in metres, derived from the domain-space values
using this script's own conventions - domain *d* occupies `[(d-1)*500, d*500)`
and plots at its band centre `(d-0.5)*500` - so an inclusive span `lo..hi`
becomes `((lo-1)*500, hi*500)` and a point feature becomes `(d-0.5)*500`.
Re-derive them if `DOMAIN_SPACING_M` changes. One mismatch is deliberately left
and commented: `GROINS` puts the Buxton groin at domain 6, band centre 2750 m,
while the hindcast places it at GIS 5.5 - the domain 5/6 boundary, 2500 m.

Note the folder name `03_cascade_inputs` overstates it: nothing reads that CSV,
and it is smoothed at this script's primary window (7 domains), whereas the
hindcast scores against the 10-domain curve. Treat it as a table to read, not
as an input to a run.

### HAT_loess_dsas_vs_coastsat.py

A different question - not method, but **source**: does DSAS agree with CoastSat?
It runs on the older **1978-1997 / 1997-2019** period pair (from
`5-scr/CoastSat/old_time_periods/`), domain-level only, at a fixed
`LOESS_FRAC = 0.15`. Not part of the 1984-2024 lineage; keep the period
difference in mind before reading its figures next to the others. The legends
name the periods and nothing else - they used to say "Calibration Period" and
"Validation Period", which name the hindcast's split, not this pair - and both
two-panel figures carry a footnote saying which pair they are.

**What it shows.** DSAS reads systematically more accretional than CoastSat in
1978-1997: +1.74 m/yr averaged over all 90 domains, r = 0.64. The offset
largely closes in 1997-2019 (-0.69 m/yr, r = 0.90). The disagreement is
structured, not noise.

Three ways these figures used to mislead, each now marked:

- **DSAS has no 1997-2019 data for domains 1-7** (83 of 90). Its line simply
  started late beside a full-span CoastSat line, which reads as agreement
  rather than absence. `shade_missing_dsas()` hatches and labels the gap on
  every panel that draws DSAS. It finds the runs from the data, so it stays
  correct if coverage changes, and draws nothing when a period is complete. It
  matters most on the difference panels, where an absent bar reads as *zero
  difference* - a meaningful value on that axis.

- **The smoothed scatter panel is not an agreement statistic.** Both series are
  LOESS-smoothed along the same domain order, so neighbouring smoothed values
  come from overlapping windows and are not independent samples: r goes
  0.64 -> 0.87 in 1978-1997 without either source having moved. The bias does
  *not* shrink under smoothing (+1.74 -> +1.95 m/yr), which is the tell - the
  offset is real and only the apparent scatter is manufactured. Both panels now
  carry their r; the smoothed one also carries the raw r and a count of
  effectively independent spans (a window holds `frac*n` domains, so there are
  about `1/frac` of them, ~7 here, not the 90 points plotted).

- **LOESS extrapolates at the southern edge.** It is a local *linear* fit, so
  at the very edge it runs past the data: CoastSat 1978-1997 smoothed to
  -6.21 m/yr at domain 1, where the raw value is -0.59 and the whole local
  spread over domains 1-6 is -7.01..-0.59. `SKIP_SOUTHERN_DOMAINS = 10`
  withholds the smoothed curves there - the same guard, at the same width, as
  the hindcast's `LoessConfig(skip_southern_domains=10)`, and for the same
  reason: Oregon Inlet dominates that zone. Display only, as it is there. The
  fit still runs over all 90 domains, so the southern data still pulls the
  values just north of the cut; only the result is withheld, and the raw series
  still plots across the whole island. Set it to `0` to smooth everywhere.

`comparison_table_smoothed.csv` follows the figures rather than the fit: its
three `*_smooth` columns are blank across domains 1-10, so the table cannot be
read as evidence the plots withheld. Raw columns are complete at 90 of 90 -
except DSAS 1997-2019 at 83, which is that source's real coverage gap and not
the guard.

## What the hindcast actually runs

`scripts/cascade_pipeline/coastsat_loess.py`, imported by
`scripts/hatteras_ms/HAT_hindcast_1984_2024.py` (section 8) and by the groin
sweep, the sensitivity plotter and `7-source-sink/loess_smooth/`. It reads the
same raw `transect_lrr_full.csv` and applies the same transect-first method
this stage chose, configured as:

    LoessConfig(window_domains=(7, 10), skip_southern_domains=10)
    TARGET_WINDOW = 10          # rate_comparison uses max(window_domains)

`skip_southern_domains=10` is **display-only**: LOESS still fits over all
transects and only the result is truncated across GIS 1-10, so the southern
transects still pull the values just north of the cut.

`hindcast_loess_snapshot/coastsat_loess.py` is a read-only copy of that module,
kept here so this folder shows the version the method ended on. It is never
imported. After editing the live module, refresh it:

    python - <<'PY'
    import pathlib
    live = pathlib.Path("scripts/cascade_pipeline/coastsat_loess.py")
    snap = pathlib.Path("scripts/input_prep/6-scr-smooth/hindcast_loess_snapshot/coastsat_loess.py")
    banner = snap.read_bytes().split(b"\n")[:14]
    snap.write_bytes(b"\n".join(banner) + b"\n" + live.read_bytes())
    PY

and update the commit named in its banner. To check whether it has drifted:

    diff <(tail -n +15 hindcast_loess_snapshot/coastsat_loess.py) \
         ../../cascade_pipeline/coastsat_loess.py

## old_smoothing/ - how the method got here

Provenance only; none of these run against the current tree (their paths still
name `input_preperation`). Oldest first:

The `v1`-`v4` order is reconstructed from file dates (2026-04-07, 05-14, 07-20,
08-06), not from anything the files themselves declare - so read the ordering as
likely, and the description as the reliable part.

| file | was | what it added |
|---|---|---|
| `HAT_loess_v1_domain_only.py` | `coastsat_smoothed_single LOESS.py` | The start. Domain-averaged only, one hardcoded `frac = 0.111`, no sensitivity. |
| `HAT_loess_v2_transect_mode_switch.py` | `coastsat_smoothed_transect_domain_dottedline.py` | First transect-first implementation, behind a hand-set `SMOOTHING_MODE` switch. |
| `HAT_loess_v3_window_sensitivity.py` | `coastsat_smoothed_final.py` | Window sizes expressed in *domains* rather than a raw frac, plus `window_comparison.png`. Still domain-only. Named "final" - it was not. |
| `HAT_loess_v4_both_modes.py` | `coastsat_smoothed_transect_domain.py` | Drops the switch - both modes every run, windows in km. Tested 5/6/7/8 domains. No method-comparison figure yet. |

The current script is that last one plus the method-comparison plots and the
5/7/10 domain window set that the hindcast's `(7, 10)` came from.
