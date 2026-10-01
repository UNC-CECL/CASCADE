# 7-source-sink — scripts

The background-erosion (BE) calibration: the per-domain source/sink field that
CASCADE carries as `DOMAIN_BE_RATES`. It is derived from what the modules could
NOT explain — the residual between the LOWESS-smoothed CoastSat rate and the
model's own LRR — so every value has a named physical zone behind it.

`scripts/site_layer/hatteras_site_config.py` is the source of truth for the field; the
copies under `data/hatteras_init/7-source-sink/` are exported FROM it by stage
4, never maintained alongside it.

## The stages

| stage | file | what it does |
|---|---|---|
| `1-prepare/` | `backfill_run_lrr.py` | Adds `lrr_m_yr` / `lrr_r2` to run rate CSVs written before the LRR estimator existed, recovering them exactly from each run's own `*_shoreline_matrix.npy`. A precondition, not a fit: stage 2 reads the LRR column and stops if a run lacks it. All 191 current runs already carry it — this is for a restored archive. |
| `2-calibrate/` | `be_zone_residual_fit.py` | The calibration. LOWESS-smooths the observed rate, differences it against the base-run LRR, identifies spatially coherent zones with a named mechanism, and writes `DOMAIN_BE_RATES*.txt` plus the metrics tables. |
| | `be_apply_fit_to_config.py` | Writes those rates into `hatteras_site_config.py`, preserving the two locked ends and the zone labels that a bare paste would destroy. `--add` accumulates instead of replacing. |
| | `be_edge_domain_solve.py` | The two locked end domains (GIS 1, 90), which are solved separately by Newton steps on a secant rather than fitted from a residual. Reads finished runs and prints the next probe; it does not run the model. |
| `3-figures/` | `plot_be_convergence.py` | Did the fixed-point solve converge, and was the zone set fixed before it ran. |
| | `plot_be_zones.py` | Which domains were eligible at all, and how much of each rate came from the one-shot solve versus the iteration. |
| | `plot_groin_reserved_residual.py` | Why the largest residual in the hindcast (D6) is deliberately left uncorrected. |
| `4-export/` | `export_be_calibration.py` | Publishes the converged field to `data/hatteras_init/7-source-sink/` — dicts, per-domain CSV, figures, provenance README. Refuses to run on figures older than the newest calibBE run. |

## Naming

Scripts here are `be_<what it solves>` for the calibration, `plot_<what>` for a
figure and `export_<what>` for the publish step. The `HAT_` prefix they all
carried until 2026-09-22 was dropped to match `5-scr` and `6-scr-smooth`; see
`../README.md` for which stages are bare and which are not.

The environment variables kept theirs: `HAT_BE_BASE_PRESET` and
`HAT_BE_OUTPUT_DIR` are unchanged, and are not script names.

Two files here are loaded BY PATH from elsewhere in the repo, by a filename
typed as a literal:

    be_zone_residual_fit.py    5 loaders, in 3-figures/, 4-export/ and
                               figure_making/model_output/
    be_edge_domain_solve.py    3 loaders, in 2-calibrate/ and
                               hatteras_ms/experiments/

Grep the name before renaming or moving either. A stale literal here does not
raise where it is written -- `spec_from_file_location` builds a spec from a
path that no longer exists, and the failure arrives at `exec_module`.

## Stage 2 is a loop, not three steps in a row

The numbering says stage 2 comes after stage 1 and before stage 3. It does not
say the three files inside it run once each in listed order. The documented
method (`hatteras_site_config.py`, METHOD) interleaves the interior fit with the
edge solve:

```
pass 0     interior from the edgeBE base runs          replace
pass 1-3   interior from the calibBE base runs         --add
GIS 90     re-solved after the interior settled        Newton, +3.0 probe
pass 4     final interior pass at the final edge values --add
```

Each pass is fit → apply → re-run the model → fit again. The additive form is
the point: imposing X m/yr of background erosion at a domain moves that domain's
rate by well under X once BRIE has diffused it alongshore, so each pass closes a
fraction of whatever misfit remains and no estimate of that fraction is ever
needed.

Zone MEMBERSHIP is identified once and frozen (`FROZEN_ZONE_DOMAINS`). Only
magnitude iterates. Re-deriving zones each pass would let the arithmetic rewrite
the science — later passes would start correcting the spillover of earlier ones,
which never terminates.

## Run

```
cd scripts/input_prep/7-source-sink

python 1-prepare/backfill_run_lrr.py --check       # only after restoring an archive

python 2-calibrate/be_zone_residual_fit.py                      # pass 0
python 2-calibrate/be_apply_fit_to_config.py --check
python 2-calibrate/be_apply_fit_to_config.py
#   re-run the model at calibBE, then for each further pass:
HAT_BE_BASE_PRESET=calibBE python 2-calibrate/be_zone_residual_fit.py
python 2-calibrate/be_apply_fit_to_config.py --add
#   and for the ends:
python 2-calibrate/be_edge_domain_solve.py --period 1996 --run <run_name>

python 3-figures/plot_be_convergence.py
python 3-figures/plot_be_zones.py
python 3-figures/plot_groin_reserved_residual.py

python 4-export/export_be_calibration.py --check
```

`--check` writes nothing anywhere it is offered. Use it first: stage 2 edits the
site config in place and stage 4 moves the previous export to `superseded_<date>/`.

## HAT_BE_OUTPUT_DIR

An exploratory pass — a different Hs, a trial base run — MUST set
`HAT_BE_OUTPUT_DIR`. The production directory holds a converged calibration
whose stopping point is a recorded scientific claim, and a what-if pass
overwriting it would destroy the provenance silently.

Both `be_zone_residual_fit.py` and `be_apply_fit_to_config.py` honour it,
so a redirected pass is applied from the directory it actually wrote. The apply
step prints the path it read before it reads it, so which calibration is being
applied is never left to be inferred from the environment.

## Where things land

Tables go to `data/hatteras_init/7-source-sink/2-calibrate/<pair>/` (mirroring the
code folder that writes them), figures to `.../7-source-sink/3-figures/<pair>/`,
and the stage 4 export to `.../7-source-sink/4-export/`, with its README at the
top of `7-source-sink/`. A pair is `<p1start>_<p1end>__<p2start>_<p2end>`, and
every pair has one, the default `1984_2004__2004_2024` included (2026-09-18;
before that the default wrote to the unlabelled root of both folders). Config
backups go to `2-calibrate/prebe/`, shared by every pair. Every script resolves
these through `scripts/site_layer/hat_source_sink.py`; do not type them.

Style: `scripts/site_layer/hat_figure_style.py`. No in-image titles or footnotes; the words
are in `CAPTIONS.md`.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### 1-prepare/backfill_run_lrr.py

Add lrr_m_yr and lrr_r2 to run rate CSVs written before the LRR estimator existed.

From the script's original header:

```text
Adds lrr_m_yr / lrr_r2 to run rate CSVs written before the LRR estimator.

WHY THIS EXISTS RATHER THAN A RE-RUN. Every finished run already carries the
whole trajectory it was scored from: `*_shoreline_matrix.npy` is the
[state x padded domain] array that `compute_change_rate` reduced to a single
column. The LRR is a second reduction of that same array, so it can be
recovered exactly -- bit for bit identical to what the run would have written
had the estimator existed at the time -- without re-simulating anything.

It is a backfill and not a fix: `change_rate_m_yr` is left exactly as the run
wrote it, and this script refuses to touch a CSV whose endpoint column does
not reproduce from the matrix. If those two disagree, the CSV and the matrix
came from different runs, and guessing which one is authoritative is not this
script's job.

Usage:
    python scripts/input_prep/7-source-sink/1-prepare/backfill_run_lrr.py [--check]

    --check  report what would change and write nothing.
```

Notes that were in the code:

```text
The endpoint column must reproduce from the matrix to this tolerance for the
pair to count as the same run. Both sides are float64 reductions of the same
array, so the only expected difference is CSV round-tripping.
```

```text
The matrix stays at the run root in both layouts, so its parent is
the run folder; the CSV's parent is tables/ once migrated.
```

<details><summary>Function notes (the original docstrings)</summary>

**`find_pairs()`**

```text
Every (rate_csv, shoreline_matrix) pair under a raw-runs tree.

Runs are DISCOVERED by their metadata file, not by the rate CSV. The rate
CSV moved into the run folder's tables/ subfolder and dropped the run-name
prefix, so a `*_shoreline_change_rate.csv` glob no longer finds a migrated
run; the metadata file stays at the run root with its prefix, and
run_layout.resolve then finds the CSV and the matrix in whichever layout
the folder is in.

Args:
    root: Directory to walk. Period and preset nesting is not assumed --
        the files are matched by belonging to the same run folder.

Returns:
    A list of (csv_path, npy_path) tuples, sorted by run directory.
```

**`backfill_one()`**

```text
Recomputes one run's LRR columns from its shoreline matrix.

Args:
    csv_path: The run's `*_shoreline_change_rate.csv`.
    npy_path: The run's `*_shoreline_matrix.npy`.
    geometry: DomainGeometry, for the real-domain slice.
    check: Report only; do not write.

Returns:
    A status string: "written", "would write", "already", or a reason the
    run was left alone.
```

</details>

### 2-calibrate/be_apply_fit_to_config.py

Write be_zone_residual_fit.py's fitted rates into hatteras_site_config.py.

From the script's original header:

```text
Writes be_zone_residual_fit.py's fitted rates into the site config.

The calibration prints a DOMAIN_BE_RATES block "ready to paste". Pasting it by
hand loses two things every time, so this does it instead:

  - THE LOCKED END DOMAINS. The generator writes GIS 1 and 90 as 0.0 with a
    "use your solved value" comment, because they are solved separately by
    buffer-cell reproduction rather than fit from a residual. Pasting the block
    verbatim silently zeroes both, and HATTERAS_BE_RATES_EDGE is SLICED from
    this table -- so it would zero the edge preset too.
  - THE ZONE LABEL ON A ZERO-VALUED DOMAIN. The generator omits the trailing
    comment where the correction solved to 0.0, so a domain's physical zone
    disappears from the config purely because its rate came out zero.

Both are preserved here: locked lines are copied verbatim from the config, and
a missing zone label falls back to the one already there.

Usage (from scripts/input_prep/7-source-sink):
    python 2-calibrate/be_apply_fit_to_config.py [--check] [--add]
```

Notes that were in the code:

```text
Resolved the way the generator resolves its OUTPUT_DIR, not typed again: the
products moved from scripts/ to the data tree on 2026-09-12 and this path was
missed, so the whole iteration step raised FileNotFoundError until 09-14. A
path spelled out in two places is a pair that can disagree, and did.

HAT_BE_OUTPUT_DIR is honoured for the same reason the generator honours it.
A what-if pass redirects its rates elsewhere; reading production anyway would
quietly apply the PREVIOUS calibration and report success. Either way the
resolved path is printed before anything is written, so which calibration is
being applied is never left to be inferred from the environment.
```

```text
THE PAIR BEING APPLIED -- read from the SAME variable the generator reads, so
the fit and its application cannot be pointed at different periods.
```

```text
THE ONE WAY TO SEED A NEW PERIOD'S GIS 90. It is not derivable from the
edge preset, so it has to be stated: HAT_BE_GIS90="1996=12.3,2010=45.6".
Naming it explicitly is the point -- the value is a separate solve, and
typing it is the acknowledgement that it was done.
```

```text
Same rule as the generator: every pair keeps its products in its own
folder (the default one too, since 2026-09-18), so applying one never reads
the other's file.
```

```text
The generator writes this file through the Windows console encoding, not
UTF-8, so its em-dashes arrive as cp1252 bytes.
```

```text
EVERY PERIOD THE CONFIG ALREADY HOLDS. The block used to be rebuilt from
the two periods this script knew, so applying a different pair would have
silently deleted the solved ones (2026-09-18).
```

```text
HARD BLOCKER, CHECKED FIRST: HATTERAS_BE_RATES_EDGE is SLICED from
the calibrated table and injects HATTERAS_BE_EDGE_D90[period]. Adding
a period the D90 table does not have produces a config that raises
KeyError on import -- which is exactly what the first real write of a
1996/2010 field did, taking every script that imports the config with
it (found and reverted 2026-09-18).
```

```text
A period with no block yet needs its two END domains, which are
solved by buffer-cell reproduction, not fitted from a residual.
```

```text
NOT "does edge have a GIS 90" -- it does, for every period. The point
is that edge's GIS 90 is a DIFFERENT QUANTITY from the calibrated
one, so its presence is no help. Only an explicitly supplied
calibBE value gets past here (corrected 2026-09-18: the first cut
tested presence and duly seeded 10.0 and 31.3 out of the edge table,
which is precisely the conflation this guard exists to prevent).
```

```text
The snapshot is this step's OUTPUT -- the BE field as it stood going in
-- so it files with the rest of 2-calibrate, not beside the config it
copies. Until 2026-09-18 these landed at the root of scripts/, where they
read as stray copies of the config: that is how the 08-24 pass-0 backup
came to be discarded, and plot_be_zones.py could not be drawn.

PRODUCTION, not HAT_BE_OUTPUT_DIR. A what-if pass redirects its rates but
still overwrites the real config, so its backup is a real backup and
belongs with the others; sending it to the what-if folder would leave the
only copy of the overwritten field outside the lineage.
```

```text
The union, in period order: a period the config already had and this
pass did not select is re-emitted exactly as it was.
```

<details><summary>Function notes (the original docstrings)</summary>

**`parse_generated()`**

```text
Reads the P1 and P2 blocks out of the calibration's output.

Args:
    path: DOMAIN_BE_RATES.txt.

Returns:
    {period: {gis: (rate, zone_label)}} for the configured pair. The file also
    holds three forecast blocks; those are ignored -- this writes the
    hindcast presets only.
```

**`existing_block()`**

```text
Locates the calibrated table in the config source.

Returns:
    (start, end, block_text) character offsets into `source`.
```

**`render()`**

```text
Renders one period's dict, keeping locked lines and zone labels.

Args:
    add: If True, the generated value is a DELTA measured against a base
        run that already carries `old_rates`, so it is added rather than
        replacing. See --add in main().
```

</details>

### 2-calibrate/be_dune_edgesolve_loop.py

Drive the dune-line end-domain solve to convergence: ask for the next probe, run it, repeat.

From the script's original header:

```text
Drive the dune-line end-domain solve to convergence: for each (window,
reading) chain, ask be_edge_domain_solve.py for the next probe, run it,
and repeat until both ends sit within --tol of the dune-line target. Written
2026-09-18 for the re-solve after the 1997, 2009 and 2023 dune lines were
re-digitized; the 09-16 solve was stepped by hand.

WHAT IT DOES NOT CHANGE
    The arithmetic is be_edge_domain_solve.py's (target from the stored
    5-scr/3-rates/duneline/endpoint product, the local secant through the last
    two runs, --estimator endpoint). This only runs the probes it prints.

LOCKSTEP
    Every live chain runs its next probe AT THE SAME TIME (one runner process
    each), then all are waited for, then every solver is read. The solver
    reads each run's imposed rates from run_index.csv, so reading only after a
    whole step has finished keeps it from reading an index a still-finishing
    run is rewriting.

BRACKETS (step 0, not re-run)
    1996, 2010   the current matrix zeroBE and edgeBE full-management runs
    2004         the 09-16 brackets (experiments/end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline/
                 brackets): the 2004-start inputs did not change on 09-18

OUTPUT   output/raw_runs/experiments/<exp>/<reading>/step<k>/<window>/edgeBE/<run>/
         output/raw_runs/experiments/<exp>/logs/<reading>_step<k>_<start>.log
         output/raw_runs/experiments/<exp>/loop_log.csv   one row per step
                                                          per chain

COASTSAT TARGET (2026-09-19)
    --target coastsat runs the same loop against the CoastSat target, as the
    matrix end values were solved (model lrr_m_yr against target_lrr_m_yr,
    GIS 1 raw, GIS 90 LOWESS-10). One chain per window, filed under the
    reading name "coastsat"; --smooth is ignored. E.g. --exp
    end-domain-boundaries/2026-09-19-end-domains-2010-recheck --windows 2010 --target coastsat.

USAGE
    python be_dune_edgesolve_loop.py --exp end-domain-boundaries/2026-09-18-end-domains-solved-on-redigitized-duneline \
        --windows 1996 2004 2010 --smooth raw mean3
```

Notes that were in the code:

```text
The matrix run names carry the offset token since the metres offset became
the default (2026-09-24; the option A matrix, 2026-09-27). The /10 brackets
the 09-16 to 09-19 solves started from are in raw_runs/archive/2026-09-24-pre-metres/.
```

<details><summary>Function notes (the original docstrings)</summary>

**`solve()`**

```text
Run the solver over `runs`; return (residual GIS 1, residual GIS 90,
override string, full text). smooth == "coastsat" is the CoastSat target.
```

**`imposed()`**

```text
{gis: rate} a finished run imposed at the two ends, from its log line
'GIS <n>  <preset> -> <imposed> m/yr' or, when not overridden, the preset.
```

**`merge_override()`**

```text
The solver prints only the ends it wants MOVED; an end it leaves out
keeps the value the chain's last run imposed. Passing the solver's string
through as-is (the first version of this driver, 2026-09-18) reset such an
end to the edgeBE PRESET -- +32.2 at GIS 1 for 1996 -- for one probe.
```

**`existing_steps()`**

```text
[(run_name, 'experiment', tag, run_dir), ...] for every FINISHED step
already on disk, oldest first, stopping at the first gap.
```

</details>

### 2-calibrate/be_dune_edgesolve_results.py

Close the books on a dune-line end-domain solve: the standing step, the pair, and how the interior scores.

From the script's original header:

```text
Close the books on experiments/end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline: which step stands
as each solve's answer, what the pair is, and how the run scores against BOTH
observations.

WHAT IT WRITES  (under output/raw_runs/experiments/end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline/)
    solved.csv      one row per (window, smooth): the solved step, its run
                    name, the pair at GIS 1 / 90, and the CoastSat-solved pair
                    it replaces. rate_windows.py reads this to find the
                    dune-solved runs (model sets dune-mean3, dune-raw).
    skill.csv       every solve AND its CoastSat-solved counterpart scored the
                    same two ways over GIS 2-89: against the CoastSat LRR
                    target (model lrr_m_yr, as run_index.csv scores) and
                    against the dune-line endpoint rate (model
                    change_rate_m_yr, as rate_windows.py vs_duneline/endpoint_net_change
                    scores).
    RESULTS.md      the two tables, rendered.

USAGE
    python be_dune_edgesolve_results.py --solved 1984:raw:3 1984:mean3:2 ...
        # window start year : smoothing : the step that converged
    python be_dune_edgesolve_results.py --exp end-domain-boundaries/2026-09-18-end-domains-solved-on-redigitized-duneline         --solved 1984:raw:3@end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline 1996:raw:2 ...
        # --exp is where the files are written and where a bare spec's runs
        # sit; "@<experiment>" carries a solve over from an earlier one (the
        # 09-18 re-solve kept 1984-2004, whose lines did not change)

    The dune target is read from 5-scr/3-rates/duneline/endpoint/ (2026-09-18),
    the same stored product rate_windows.py draws.
```

Notes that were in the code:

```text
The 1984 and 2004 brackets were run once, under the 09-16 experiment; their
inputs did not change with the 09-18 re-digitization (the 1984 and 2004
lines, and the run itself never reads a dune line), so every re-solve
reuses them from there.
```

```text
The CoastSat-solved counterpart of each window: the matrix run for 1996
and 2010, the fresh bracket for 1984 and 2004 (same code, same versions).
1996 and 2010 carry the offset token since the metres offset (the option A
matrix, 2026-09-27); 1984 and 2004 are the /10-era brackets and have no
metres run.
```

```text
full per-domain dune rate, raw (the interior score does not smooth),
from the stored product
```

### 2-calibrate/be_edge_domain_solve.py

Solve the background-erosion rate the two locked end domains (GIS 1 and 90) carry, for one period.

Notes that were in the code:

```text
be_edge_domain_solve.py

What background-erosion rate should the two LOCKED END DOMAINS carry, for one
hindcast period?

WHAT THE END VALUES ARE FOR
GIS 1 and GIS 90 sit on the open boundaries of the modelled reach and
absorb the artefact there. They are boundary-artefact absorbers, not a
sediment budget -- see the end-domain note in hatteras_site_config.py,
which this script does not restate and must not contradict.

WHY A SCRIPT
The solve is a Newton iteration: run, read the residual at the two ends,
step, run again. It was done by hand for 1984 and 2004, and the arithmetic
between runs -- which target column, which estimator, which secant -- was
carried in a person's head. Two periods were added on 2026-09-11 and both
need the same solve, so the arithmetic is written down here.

It does NOT run the model. It reads runs that have already happened and
prints the next probe, so every step stays a deliberate act.

THE TWO TARGETS ARE DIFFERENT ESTIMATORS, DELIBERATELY
GIS 1  the raw per-domain transect mean. LowessConfig.skip_southern_domains
is 10, so D1-D10 are drawn raw rather than smoothed.
GIS 90 the LOWESS value (TARGET_WINDOW, 7 since 2026-09-28; 10 before),
which is what is drawn everywhere north of D10.
That splice is what the rate-comparison figure draws, so fitting against
the same table means fit and figure cannot disagree. Both come out of
build_target_table, so neither is computed here.

THE GAIN IS SMALL AND NOT CONSTANT
d(LRR)/d(BE) ran 0.092 to 0.123 across the four solved cases: only about a
tenth of an imposed edge rate survives in that domain's own shoreline, the
rest being diffused alongshore by BRIE. So each value is roughly ten times
the misfit it closes, and a single global slope should not be assumed --
which is why the second step uses the LOCAL secant through two real runs
rather than the nominal gain again.

AN EXTENDED GEOMETRY (2026-09-16, the Pea Island extension experiment)
The ends are wherever HATTERAS_DOMAINS puts them -- GIS 1 and 115 under
HAT_GEOMETRY=n115 -- and the target for a
domain beyond GIS 90 comes from the window's extension rate table
(coastsat_lrr/<window>/ext/transect_lrr_with_base.csv, the surveyed
transects plus the extension's, one LOWESS over the whole reach), which is
exactly the table the runner grades that geometry's ends against. Run
this script with the SAME HAT_GEOMETRY as the runs it reads. A run that
imposed nothing at an end (a zeroBE probe) reads as 0.0 there.

THE DUNE LINE AS THE TARGET (2026-09-16, Hannah: "what if we used the dune
line change" to set the ends)
--target duneline reads the observation from the two digitised dune lines
of the window instead of CoastSat: the end vintage's line minus the start
vintage's, per domain, over the survey interval (coastsat_vs_duneline
.KNOWN_SURVEY_DATES; a missing date is mid-year), seaward positive --
exactly what rate_windows.py draws under vs_duneline/endpoint_net_change. Two readings of
it at an end domain, --dune-smooth raw (the domain's own value) and mean3
(the mean of it and its two inward neighbours, GIS 1-3 / 88-90). A dune
line is two surveys, so --estimator endpoint reads change_rate_m_yr on
the model side; the default lrr keeps the CoastSat protocol. The runs of
that solve are under output/raw_runs/experiments/end-domain-boundaries/2026-09-16-end-domains-solved-on-duneline/.

USAGE
One run -- report the residual and a first step at the nominal gain:
python 2-calibrate/be_edge_domain_solve.py --period 1996 --run <run_name>

Two or more -- local secant through the last two, and the next step:
python 2-calibrate/be_edge_domain_solve.py --period 1996 --run <first> --run <second>

Runs are named as they appear in run_index.csv and located through
run_registry.find_run_dir: --kind and --tag say where (the matrix by
default; an experiment's runs need --kind experiment --tag <set>/<member>,
one --tag per --run when the members differ). The rate each was run
under is read from the index, not retyped.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-30
```

```text
Section 8 of the runner builds the target this way. Kept identical rather
than imported from it, because importing that file RUNS a hindcast.
```

```text
The model side of the residual. Must be the OLS slope, matching the
observed side; change_rate_m_yr is the endpoint difference and every preset
fitted against it before 2026-08-22 is not reproducible from this pipeline.
```

```text
d(LRR)/d(BE), used ONLY for the first step, when there is nothing to take a
secant through. Mid-range of the four solved cases.
```

```text
The index column holding the rate each run imposed, per end domain. The
runner writes be_rate_gis<N>_m_yr for each of HATTERAS_BE_EDGE_DOMAINS.
```

```text
raises, listing windows, if absent; the extension table if the
geometry reaches beyond GIS 90, as the runner's section 8 does
```

```text
The preset folder the run sits under is the preset it ran, which
the index knows: a zeroBE stage-0 run and the edgeBE probes after
it belong to one solve and are read together.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_target()`**

```text
The target table the runner grades against, as {gis: rate}. `window`
(start, end) grades against another CoastSat LRR window instead, e.g.
the full 1996-2024 rate for a 1996-2010 run (2026-09-19, Hannah).
```

**`load_dune_target()`**

```text
The dune-line endpoint rate per domain, seaward positive, read at each
end domain raw or as a three-domain mean. Returns ({gis: rate}, note).

READ FROM the stored product 5-scr/3-rates/duneline/endpoint/<window>/
(2026-09-18), the same numbers rate_windows.py draws, rather than
recomputed here from the raw offsets.
```

**`imposed_rates()`**

```text
What this run imposed at each end domain, from its index row.

A run from a preset that names no rate at an end (zeroBE, or an edgeBE
whose column predates this end) reads as 0.0 there, which is what it
imposed.
```

</details>

### 2-calibrate/be_zone_residual_fit.py

Derive the source/sink (background erosion) field from the residual between smoothed CoastSat and a base run.

From the script's original header:

```text
Derives defensible background erosion (BE) source/sink corrections for CASCADE
from the residual between a LOWESS-smoothed CoastSat observed shoreline change
rate and the CASCADE base-run LRR.

Philosophy
Corrections are applied only where THREE conditions are simultaneously met:
  1. The residual (smoothed observed rate minus CASCADE LRR) exceeds a
     significance threshold (not just noise)
  2. The signal is spatially coherent across multiple adjacent domains
  3. A physical mechanism can be named for that zone

This minimises the number of free parameters and maximises scientific
defensibility — every correction has a name and a reason.

Workflow
  1. Load CoastSat domain-averaged LRR (P1 and P2) — the observed shoreline
     change rate
  2. Load CASCADE base-run LRR from NPZ (P1 and P2, management included, BE=0)
  3. LOWESS-smooth the observed shoreline change rate (7-domain window; 10 until 2026-09-28),
     excluding domains 1-GROIN_EXCLUDE_THROUGH_DOMAIN (Buxton groin influence
     zone) from the fit entirely — those domains pass through with their raw,
     unsmoothed rate
  4. Compute residual = smoothed CoastSat rate - CASCADE LRR, per domain per
     period (the fully raw, unsmoothed residual is also retained, for
     diagnostic comparison only — it no longer drives any decision)
  5. Identify significant zones (|residual| > SIGNIFICANCE_THRESHOLD,
     spatially coherent over >= MIN_ZONE_WIDTH adjacent domains)
  6. Classify each domain: stable correction vs shifting between periods
  7. Output:
       - Diagnostic figures showing raw/smoothed rate, residual, and zone
         identification
       - DOMAIN_BE_RATES dicts for P1, P2, and three forecast scenarios
       - Summary CSV with all metrics and physical zone assignments

Forecast scenarios (for domains where P1 ≠ P2 correction)
  "continue" : use P2 correction (current trajectory continues into future)
  "revert"   : use P1 correction (system returns to pre-2004 state)
  "neutral"  : use mean(P1, P2) (no prior on future state)

Usage
  python 2-calibrate/be_zone_residual_fit.py

Dependencies
  pip install pandas numpy matplotlib scipy statsmodels tqdm
```

Notes that were in the code:

```text
The pipeline's own code builds both inputs below. Imported, never copied:
this file used to reimplement the observed smoothing and the modelled LRR
extraction, and both had drifted from what the model is scored against.
```

```text
CoastSat TRANSECT-level LRR. Not the domain-averaged summary: the target
is smoothed at transect resolution and only then averaged, and doing it in
the other order gives a measurably different curve.
Resolved through hat_observed_rates.py (2026-09-18), not typed.
```

```text
THE PAIR BEING FITTED. Set HAT_BE_PERIODS to two comma-separated period
STARTS to fit a different pair; each end comes from HATTERAS_PERIODS, and
the CoastSat product for a window lives under <start>_<end>/, so naming
the starts names the observations too. Default is the pair this field was
originally solved on (generalised 2026-09-18).

P1/P2 mean EARLIER and LATER throughout this file -- including the _p1 /
_p2 columns in the metrics CSV -- not 1984 and 2004.
```

```text
The base run, per period. edgeBE is the preset that is ZERO at every domain
being solved (2-89) while carrying the independently solved values at the
locked ends -- so an interior residual is the whole correction rather than
an increment on an existing one, and the open-boundary artifact at GIS 1/90
is absorbed there instead of diffusing alongshore into interior terms.

full_management because the observed CoastSat rate is from a real island
that WAS managed. nogroin because domains 1-10 are already excluded from
the LOWESS fit for groin influence, and running the base with the groin on
would push its signal into the residual and double-count it against the
separate M/f sweep.
ITERATION. The calibration is a ONE-SHOT solve: it measures the residual of a
base run and imposes it as the BE field. That is exact only if imposing X m/yr
moves the domain's LRR by X m/yr, and it does not -- BRIE diffuses an imposed
rate alongshore, so a domain keeps only a fraction g of what it is given.
Measured on 2026-08-24, g is not a constant: a contiguous same-signed block of
corrections passes at g ~ 0.8-1.2, while a pattern that alternates sign at the
grid scale is damped to g ~ 0.1. One pass therefore closes 42% (P1) and 57%
(P2) of the misfit rather than all of it, and the shortfall is uneven.

The fix is to iterate rather than to guess g: point this at the CURRENT
calibBE runs, measure what residual is left, and ADD it to the field already
in place (be_apply_fit_to_config.py --add). Each pass closes fraction g of whatever
remains, so it converges whatever g turns out to be, and it needs no estimate
of g at all. Amplifying by 1/g instead was rejected: g varies by an order of
magnitude with wavelength, so narrow features would be amplified ~10x into
rates that are indefensible read as sediment fluxes.

pass 0   HAT_BE_BASE_PRESET unset -> edgeBE   ->  be_apply_fit_to_config.py
pass 1+  HAT_BE_BASE_PRESET=calibBE           ->  be_apply_fit_to_config.py --add

Stop when no zone clears SIGNIFICANCE_THRESHOLD, which is then the tolerance
the field is converged to.
```

```text
BASE_SCENARIO was here until 2026-09-18: defined once, read nowhere, and
wrong for a period with fills (those runs carry a nourish token). The
base run is resolved by globbing the run stem instead -- see _base_run.
```

```text
Where a CONCLUDED experiment's forcing arms are kept. A live calibration
probe is still written into raw_runs/ under HAT_ARM_TAG -- that is what stops
it overwriting the base run it probes -- and is moved here when the question
it was answering is settled, so raw_runs/ stays the production matrix.
```

```text
The section 8 settings, matching the runner. TARGET_WINDOW is the widest
window; `rate_comparison` resolves the reference the same way.
```

```text
HAT_BE_OUTPUT_DIR redirects every output -- be_zone_metrics.csv,
DOMAIN_BE_RATES*.txt, convergence_history.json, the figures. A what-if pass
(a different Hs, a trial base run) MUST set it: the production directory
holds a converged calibration whose stopping point is a recorded scientific
claim, and an exploratory pass silently overwriting it would destroy the
provenance without anyone noticing.
Products moved to the data tree 2026-09-12; only the script lives under
scripts/. HAT_BE_OUTPUT_DIR still redirects a what-if pass anywhere.
EVERY PAIR WRITES TO ITS OWN FOLDER, the default one included (2026-09-18;
hat_source_sink.py). HAT_BE_OUTPUT_DIR always wins. Before, only a
non-default pair got a folder, added after a run on another pair wrote over
the committed 1984/2004 field; the default wrote to the unlabelled root.
```

```text
Figures are read out of the data tree, not out of scripts/. The tables above
stay with the calibration that produced them; the two PNGs go where the rest
of the section 7 figures live, so a re-run refreshes the copies people open.
A what-if pass with HAT_BE_OUTPUT_DIR set keeps its figures with its tables.
```

```text
── Correction thresholds ─────────────────────────────────────────────────────
Minimum smoothed residual magnitude to warrant any correction at all.
Below this = within model noise/uncertainty, leave at zero.
```

```text
Minimum number of adjacent domains with |residual| > threshold to form a zone.
Prevents correcting isolated noisy domains.
WHICH MODEL COLUMN THE RESIDUAL IS BUILT FROM. Must be the same
estimator as the observed target, which is a per-transect OLS slope;
see load_model_lrr. "change_rate_m_yr" reproduces a pre-2026-08-22
calibration and nothing else.
```

```text
WHICH BASE RUN THE RESIDUAL IS DERIVED FROM.
True   prefer the groin-ON base run carrying the joint fit's (M, f),
so the source/sink carries only what the MODULES could not
explain. Falls back to the no-groin run, loudly, when the
joint fit has not run or no matching run exists.
False  always the no-groin run. Reproduces every calibration before
2026-08-22, and hands the source/sink the groin's D5/D6
signal to absorb.

ORDERING. Fit the groin against zeroBE/edgeBE, which impose nothing at
D5/D6, then recalibrate the source/sink against a run carrying the
fitted groin. Running this with the switch on BEFORE the joint fit is
harmless -- it warns and falls back -- but the result is a pre-groin
calibration and should not be quoted as a post-groin one.
```

```text
If |P1_correction - P2_correction| exceeds this, the zone is "shifting"
and needs period-specific values + forecast scenarios.
```

```text
── LOWESS smoothing ───────────────────────────────────────────────────────────
Fraction of data used for each local regression (larger = smoother).
7-domain window matches the CoastSat LOWESS calibration window used
throughout the dissertation (10 until 2026-09-28). As of this version, smoothing is applied to
the OBSERVED SHORELINE CHANGE RATE itself (before differencing against
CASCADE) — not to the residual. LOWESS_FRAC is calibrated for the full
90-domain array; see smooth_shoreline_rate() for how it's re-derived once
the groin zone is excluded from the fit.
```

```text
Domains 1 through this value are excluded ENTIRELY from the shoreline-rate
LOWESS fit (Buxton groin influence zone) — the groin's localized signal
would otherwise bleed into the smoothed estimate at neighbouring domains.
These domains always keep their raw, unsmoothed CoastSat rate.
```

```text
── Manual overrides ─────────────────────────────────────────────────────────
Domain-level corrections that override the LOWESS-derived value.
Use sparingly — only where the smoothing window demonstrably under/over-corrects
and you have a clear physical justification for the different value.
Format: domain → (p1_override, p2_override, reason)
Set either value to None to keep the LOWESS-derived value for that period.
Forecast scenarios for overridden domains use the same logic as normal:
continue = p2, revert = p1, neutral = mean(p1, p2)
```

```text
No overrides — pure LOWESS-derived values against the current baseline.
This script is for DISCOVERY: see what the data says before any
manual intervention. Once you have chosen final values (informed by
this comparison plus your own judgement), enter them in
HAT_be_apply_final_rates.py to generate the final figures and
forecast scenarios.
```

```text
── Locked domains ────────────────────────────────────────────────────────────
Domains where the BE rate has already been solved independently (e.g. D1 and
D90, found by reproducing the buffer-cell shoreline change rate directly).
These are forced to 0.0 in ALL DOMAIN_BE_RATES outputs and excluded from the
significance/strategy calculation entirely — the script will not suggest a
correction for these domains, since you will always supply your own value.
Add a short note here for your own records of what each locked value is.
── Frozen zone set ───────────────────────────────────────────────────────────
THE ZONES ARE THE SCIENCE; THE MAGNITUDE IS THE ARITHMETIC. Zone membership
says "this stretch of coast has a real sediment-budget deficit, and here is
the process". Magnitude says "deliver the amount you diagnosed, given that
BRIE diffuses roughly half of it away". Iterating BOTH lets the second
quietly rewrite the first: each pass re-derives zones from a NEW residual, so
as the coherent features are satisfied, progressively less coherent ones
cross the 0.5 m/yr threshold and get corrected. Worse, adding BE at a domain
pushes sediment into its neighbours and changes THEIR residuals -- so later
passes partly correct the alongshore spillover of earlier passes. That is
bookkeeping, not geomorphology, and it never terminates.

Measured 2026-08-24: two unmasked passes corrected 19 domains in 1984-2004
and 12 in 2004-2024 that the pass-0 zone identification never selected --
including D5-D7 in period 2, the groin's own footprint.

So zones are identified ONCE, from the edgeBE residual under the ordinary
significance and width rules, and then held fixed. Everything outside stays
at 0.0 no matter what its residual does; that residual remains in the results
as honest unexplained variance, which is the defensible number anyway.

Regenerate ONLY by re-running pass 0 against edgeBE and re-deriving the set;
do not edit a domain in here to chase a residual.

CORRECTED 2026-09-14 -- SEVEN DOMAINS WERE IN THE WRONG PERIOD'S TUPLE.
D8, D48, D57 sat in 2004 and D9, D22, D44, D62 in 1984, and every one of the
seven is warranted by the ordinary rules in the OTHER period -- where each was
already a legitimate member (D22 inside period 2's D8-D22 run, D48 inside
period 1's D48-D57). They had been written into both tuples. Seven for seven
is a transposition when the sets were assembled, not drift in the base runs.

It showed up as five corrections ONE DOMAIN WIDE -- D22, D44, D62 here in
1984, D48 and D57 in 2004 -- which MIN_ZONE_WIDTH = 3 exists to make
impossible. D8 and D9 are the same error, invisible because they land beside
legitimate members.

The sets below are now exactly what identify_correction_zones returns for the
pass-0 edgeBE residual, plus D5-D7 in 1984 where the width rule puts them
anyway. Neither has an isolated member. Re-derived, not hand-pruned: the width
rule was re-applied to the residuals directly, because compute_be_rates
overwrites the verdict at D5-D7 after the rule has run, which makes D8 read as
isolated in the metrics CSV when it is not.

The superseded field and the full account are in
data/hatteras_init/7-source-sink/superseded_20260914/.
```

```text
── Domains reserved for the groin module ─────────────────────────────────────
D5-D7 is the Buxton groin's own footprint, and the residual left there is the
GROIN's residual, not a source/sink term. Measured 2026-08-24, after two
iteration passes:

period 1   D6 = +1.72 m/yr   observed is more seaward than modelled --
M = 60 does not build enough fillet
period 2   D6 = -2.21 m/yr   modelled is more seaward than observed --
the fillet does not release, which a module
with trapping bounded at >= 0 CANNOT do

Letting the BE field absorb those is exactly the double-count that
GROIN_AWARE_BASE_RUN exists to prevent: the source/sink term would quietly do
the groin's job and the groin would look better calibrated than it is.

THIS WAS ALREADY HAPPENING, BY ACCIDENT. D5-D6 (period 1) and D6-D7 (period 2)
each form a significant run only 2 domains wide, and MIN_ZONE_WIDTH = 3 means
no zone forms, so no correction is ever applied. That is the right outcome
reached by a rule that knows nothing about the groin -- and it would silently
reverse if MIN_ZONE_WIDTH were ever retuned. Naming the domains here makes the
behaviour intentional and survives that.

FREEZE, NOT ZERO. These emit 0.0, which under `be_apply_fit_to_config.py --add` means
"add nothing" -- the value already in the config is kept. That matters: at D5
only ~16% of the standing correction is the groin, the other ~84% being Cape
Point background measured with the groin switched OFF. Zeroing would discard
it. Under the pass-0 replace path 0.0 would zero them, so reserve domains only
once the field they should keep is already in place.
```

```text
── Physical zone definitions ─────────────────────────────────────────────────
These are your prior hypotheses about where physical mechanisms operate.
The script will test whether the residual data supports them.
Format: zone_name → (domain_start, domain_end, mechanism_description)
```

```text
Display-only shortenings for the in-place zone strip labels (fig_be_rates).
The canonical PHYSICAL_ZONES names above are kept everywhere else — CSV
comparison, DOMAIN_BE_RATES comments, mechanism lookups — since the longer
names carry more information there. Add more entries here if you want
other zones shortened on the chart too.
```

```text
── Alongshore annotation ────────────────────────────────────────────────────
The village spans come from the site config through `town_bands()`, so this
file cannot disagree with it. What is left here is the structures: the two
piers and the Buxton groin, drawn as rulers in the muted ink.
```

```text
Type comes from hat_figure_style.apply_style(); only the two sizes that are
deliberately smaller than the 8 pt tick default are named here.
```

```text
The floor is not stored as its own field. The runner writes the
deterioration as a HUMAN-READABLE STRING -- "linear_ramp, floor 0.6"
-- so `deterioration_fraction` is always None and, before this,
every candidate was skipped and the analysis fell back to the
no-groin base while reporting "no groin base run at the fitted
M, f". The runs were correct; only the lookup was wrong, and the
fallback is a warning rather than an error, so it produced a
complete pre-groin calibration that LOOKED like a post-groin one.
```

```text
hatteras_ms is not on the path for this script the way SCRIPTS_DIR is.
Added inside the function, mirroring HAT_groin_sweep_config._wave_scope,
so the module does not gain an import-time dependency on the runner.
```

```text
Runs are filed [<forcing arm>/]<period>/<preset>/, the arm being absent
at the calibration climate. HAT_BE_HS points this at the tree for a
different wave height, which is what lets the SAME calibration method be
run against a differently-forced model and the two compared.

BOTH THE ARM'S SPELLING AND THE LAYOUT COME FROM SHARED CODE. This block
previously spelled "waveHs3" itself and joined the path itself, so it
held a second copy of two rules the runner also implements -- and a
second copy is how the two drift into reading different directories.

WHICH ROOT, THOUGH. The calibration arm lives in raw_runs/; an
off-calibration arm does NOT. The Hs 3.0 arms were moved to
hs_experiment/runs/ on 2026-09-02, beside the DECISION.md they are the
evidence for, so that raw_runs/ holds only the production matrix and its
sensitivity cells and no run name appears there twice. The layout INSIDE
each root is identical, which is why one preset_dir_for call serves both.
```

```text
The matrix, wherever the registry files it (raw_runs/matrix/ since
2026-09-16, with the older layouts still readable).
```

```text
hs_experiment/runs/ keeps its 2026-09-02 shape, <arm>/<period>/
<preset>/, and is a closed experiment; spelled here because the
registry no longer knows that layout.
```

```text
GROIN-ON BASE RUN, WHEN ONE EXISTS AT THE FITTED (M, f).

The source/sink field is meant to carry what the MODULES could not
explain. Deriving it from a groin-OFF run hands it the groin's whole
signal at D5/D6, and a later calibBE-plus-groin run then applies both
-- the double count section 7 of HAT_hindcast_methods.md warns about.
It is not hypothetical: the 2026-08-22 calibration put -1.4 m/yr at
D6 in 2004-2024, absorbing exactly the fillet relaxation the groin is
now known to be able to produce.

THE (M, f) MATCH IS THE POINT, not merely finding a groin run. A seed
run exists at PROVISIONAL values (M = 50, f = 0.9) purely to give the
sweep a drift-guard reference; calibrating against that would fit the
source/sink to a groin nobody has fitted yet. So the run must carry
the values in joint_fit.json, and until the joint fit has run there
is nothing to match and this falls back -- loudly.
```

```text
RESOLVED, NOT GLOBBED. The rate CSV is tables/shoreline_change_rate.csv
in the new run layout and {run}_shoreline_change_rate.csv in the old,
so a glob on the old name silently finds nothing in a migrated run.
```

```text
Groin-reserved domains: emit 0.0 and skip the strategy logic, so an
iteration pass adds nothing here and the standing value is kept.
See GROIN_RESERVED_DOMAINS for why this residual is not ours.
```

```text
Locked domains: force to 0.0, skip all significance/strategy logic.
These domains already have an independently-solved BE rate that you
will always supply yourself — the script should not suggest a value.
```

```text
Correction value to apply: use smoothed residual (spatially coherent)
Only apply where warranted
```

```text
EXPLORATORY PASS. With HAT_BE_FREEZE=off the zone set is not applied, so
a period that has no frozen set yet can still be RUN and its candidate
zones read off the metrics CSV and the diagnostic figure. The output is a
diagnosis, not a field: every domain that cleared the significance and
coherence tests keeps its correction, including the grid-scale wiggles
the frozen set exists to withhold. It must not be applied to the config
(2026-09-18, added with the period generalisation so the message that
points here is true).
```

```text
Hold the zone set fixed -- see FROZEN_ZONE_DOMAINS. Applied here rather
than inside the loop so the metrics CSV still records the residual and
the significance verdict for every domain: the diagnosis stays visible,
only the correction is withheld.
```

```text
Shade warranted zones, distinguishing APPLIED from WITHHELD.
correction_warranted_* is deliberately left unmasked so the metrics
CSV keeps the diagnosis for every domain -- but a domain outside
FROZEN_ZONE_DOMAINS is diagnosed and NOT corrected, and shading the
two alike would tell the reader a correction was applied where none
was. The withheld class is hatched rather than given a third
saturated colour: it is secondary to both.
```

```text
Wimble Shoals is named in panel (c); a third grey band here would be one
too many. The significance band belongs to the residual, not to these
rates, so it is off as well.
```

```text
Domains whose correction differs between the periods, as a strip at the
foot: as a full-height wash it covered half the panel and swamped the
bars it was meant to qualify.
Which domains hold DIFFERENT values in the two periods -- read off the
field itself rather than from results["strategy"], which records how this
pass classified the residual and so changes pass to pass. SHIFT_THRESHOLD
is the same rule the strategy column applies.
```

```text
The three scenarios are DEFINITIONS over the two hindcast fields, not
separate fits: continue = carry period 2 forward, revert = restore period
1, neutral = the mean. Derived here from the config field for the same
reason panel (a) is -- the results[] columns carry this pass's increment,
so a forecast drawn from them forecasts the leftover.
```

```text
-- Panel 3: the physical zones, and where a correction was applied -----
The zone strip is not a categorical colour scale: each zone is named in
place, so the fill is free to carry the one thing the names cannot, which
is whether the calibration actually corrected that domain. A seven-colour
palette here would also collide with the vintage pair above.
```

```text
"Correction applied" means the FIELD carries a nonzero value in either
period, which is what the legend claims. It used to read
results["strategy"] != "zero" -- this pass's significance verdict, which
at a late pass marks the domains still moving rather than the domains
corrected, and disagreed visibly with the bars above it.
```

```text
Direct in-place zone labels, fitted to each zone's own width -- drawn
after a draw() so the pixel-width measurements used to size each label
reflect the figure's final axes geometry, not a stale pre-layout one.
```

```text
The observed curve is already smoothed -- `build_target_table` did it at
transect resolution. `smooth_shoreline_rate` is deliberately NOT called
here any more: running it would LOWESS an already-LOWESSed curve, and the
second pass would flatten exactly the coherent zones this script exists
to detect. The function is kept for reference by the diagnostic figure.
```

```text
── Raw residual — fully unsmoothed on both sides, kept for diagnostic
comparison only. It no longer drives any decision below. ────────────────
```

```text
── Residual from the smoothed observed rate — THIS is what drives zone
identification, significance testing, and BE corrections from here on ──
```

<details><summary>Function notes (the original docstrings)</summary>

**`_never_die_on_a_print()`**

```text
Stop a console encoding from killing a finished computation.

This script prints arrows, ellipses and plus-minus signs. A Windows
console is cp1252, which cannot encode any of them, so `print` raises
UnicodeEncodeError -- and on 2026-08-28 that happened at the very last
status line, AFTER the whole calibration had been computed and BEFORE
DOMAIN_BE_RATES.txt was written. The run looked like a crash, the numbers
were gone, and the stale file left behind still carried the previous
pass's date, so nothing downstream would have noticed it was old.

Reconfiguring is preferred to ASCII-ifying every print: the next arrow
someone types would reintroduce the bug, and the failure mode is silent
data loss rather than a wrong character. If UTF-8 cannot be set, fall
back to errors="replace" so an unencodable character degrades to "?"
instead of raising.
```

**`frozen_zones()`**

```text
The domains that may receive a correction in this period.

THIS TABLE IS A SCIENTIFIC JUDGEMENT, NOT A COMPUTATION. A correction is
applied only where the residual is significant, spatially coherent AND a
physical mechanism can be named for the zone; the first two this script
measures, the third a person decides. So a period with no entry raises
rather than defaulting to "everywhere" (which would apply every
grid-scale wiggle) or to "nowhere" (which would silently fit nothing).

To add a period: run the fit once to see which zones clear
SIGNIFICANCE_THRESHOLD and MIN_ZONE_WIDTH, name a mechanism for each one
you accept, and list its domains here (2026-09-18).
```

**`_fitted_groin()`**

```text
The (M, fraction) the joint fit settled on, or None.

Args:
    preset: Which preset's fit to read. The base run is BASE_PRESET,
        so its own fit is the one that describes it.

Returns:
    An (M, fraction) tuple, or None if the joint fit has not run or
    holds no entry for this preset.
```

**`_groin_run_at()`**

```text
The groin base run carrying exactly the fitted (M, f), if present.

Matched on the run's own metadata rather than on its directory name:
the name records that a groin ran, not which one.

Args:
    period_dir: <period>/<preset> directory to search.
    stem: Invariant leading part of the run name.
    fitted: (M, fraction) to match.
    tolerance: Absolute tolerance on both values.

Returns:
    The matching Path, or None.
```

**`_wave_arm()`**

```text
The forcing arm HAT_BE_HS selects, as run_registry spells arms.

CALIBRATION_ARM at the calibration wave climate, so every run made before
arms existed resolves exactly where it always did. The token comes from
wave_climate_token -- the same function the runner derives its own arm tag
from -- rather than being spelled here, so the two cannot disagree about
which directory a run was filed in.

Returns:
    An arm name for preset_dir_for.
```

**`base_run_dir()`**

```text
The base run directory for one period, resolved from what is on disk.

The scenario tokens are NOT the same in both periods: period 2 has
nourishment scheduled, so its full_management run carries a `nourish`
token that period 1 has no reason to. Globbing the invariant part and
excluding the arms that must not be picked is what keeps this from
silently resolving to the wrong run when the token set changes again.

Raises:
    FileNotFoundError: If no base run is present. Loud rather than
        falling back to another preset -- a calibration silently derived
        against the wrong base would look entirely normal downstream.
    RuntimeError: If more than one candidate matches, rather than
        guessing which run the calibration should rest on.
```

**`load_model_lrr()`**

```text
Per-GIS-domain modelled LRR, m/yr, (+) seaward.

Read from the run's own shoreline change rate CSV rather than
re-derived from the .npz. The pipeline writes that file from the same
array section 12 scores, so this cannot disagree with the model about
sign, units, or padded-index alignment.

Takes `lrr_m_yr`, the OLS slope through the run's annual states --
NOT `change_rate_m_yr`, which is (x[-1] - x[0]) / span. The residual
this calibration turns into a background-erosion rate is model minus
observed, and the observed side is a per-transect OLS slope, so the
model side has to be the same estimator or the residual carries the
difference between two estimators as if it were a sediment budget.
Every BE preset before 2026-08-22 was fit on the endpoint column;
this function is the reason those values are not reproducible from
the current pipeline without setting RATE_COLUMN back.
```

**`load_observed()`**

```text
(raw_per_domain_mean, target) for one period, m/yr, (+) seaward.

`target` is the curve the runner's section 8 builds and section 12 grades
against: LOWESS at transect resolution over along-coast distance, averaged
to domains, with GIS 1..skip_southern_domains spliced in as raw means.
`raw` is the unsmoothed per-domain mean, kept for the diagnostic residual
only -- it drives nothing.
```

**`smooth_shoreline_rate()`**

```text
LOWESS-smooth a domain-indexed OBSERVED shoreline change rate, excluding
domains 1..exclude_through entirely from the fit (Buxton groin influence
zone) and passing those domains through unchanged with their raw rate.

Domains > exclude_through are smoothed using ONLY data from domains
> exclude_through — the groin-zone values never enter the regression at
all, so they cannot bleed into the smoothed estimate near the zone
boundary (this is stronger than merely overwriting the comparison for
domains 1..exclude_through after smoothing over the full array).

frac is re-derived here rather than reusing LOWESS_FRAC directly: LOWESS_FRAC
(7/90) is calibrated to give a 7-domain window when fit over all 90
domains. Once the groin zone is excluded, only 80 domains remain in the
fit, so reusing 7/90 unchanged would narrow the window to ~6.2 domains.
Recomputing frac = window_domains / n_valid preserves the true ~7-domain
(~3.5 km) window this dissertation uses everywhere else.
```

**`identify_correction_zones()`**

```text
Find contiguous runs of domains where |smoothed_residual| > threshold
and the run is at least min_width domains wide.

Returns array of booleans: True = correction warranted.
```

**`compute_be_rates()`**

```text
For each domain, determine the appropriate BE correction strategy.

Rules:
  - If smoothed residual not significant in either period → BE = 0
  - If significant in one or both periods:
      - If |correction_P1 - correction_P2| < SHIFT_THRESHOLD → stable → use mean
      - Otherwise → shifting → flag for scenario treatment, provide P1/P2/mean
```

**`annotate_ax()`**

```text
The alongshore furniture every panel shares.

Villages come from the site config through `town_bands()` as a strip along
the top edge -- a full-height wash cannot be told apart from the zone
shading these panels already carry. Wimble Shoals gets the matching strip
along the bottom, on the panels that do not name it some other way. The
piers and the groin are rulers, so they are drawn in the muted ink rather
than a colour of their own, and their names sit at the foot of the panel
under a white halo so they never have to fight the data for a place.

`thresholds` draws the significance band. It belongs on a residual panel
and nowhere else: on a panel of background-erosion rates the same pair of
lines would imply a test that was never applied to those numbers.
```

**`find_zone_runs()`**

```text
Collapse a per-domain physical_zone Series into contiguous
(start_domain, end_domain, zone_name) runs, in domain order.
```

**`_contiguous()`**

```text
Collapse a sorted list of domains into (first, last) runs.

A local convenience so a set of domains can be handed to `town_bands()`
as spans; `find_zone_runs` above needs a per-domain Series instead.
Worth lifting into hat_figure_style if another script wants it.
```

**`label_zone_runs()`**

```text
Place one centered text label per zone run directly on the strip,
shrinking that label's own font (and only that one) until it fits
within its own zone's width — so a narrow zone with a long name
(e.g. "Cape Point / Shoal Dynamics") never bleeds into its neighbour.
```

**`plot_diagnostic()`**

```text
Five stacked panels: the two rates, the two residuals, the strategy.

A working diagnostic rather than a manuscript figure -- it is drawn to the
house size and type so it can sit beside the others, and no further.
```

**`_field_from_config()`**

```text
The CALIBRATED FIELD as the model actually carries it, from the config.

NOT results["be_hindcast_p*"]. That column is what THIS pass proposes, and
the two are the same thing only at pass 0, where the apply step replaces
rather than adds. At every later pass it is an increment, so a figure drawn
from it is a picture of the last step rather than of the field -- and at
convergence, when the increment is near zero by definition, it is a picture
of almost nothing under a title that says "hindcast field". That is what
this figure showed on 2026-09-14: 5 nonzero domains at max 0.9 m/yr against
a real field of 43 at max 5.3.

Reading the config instead makes the figure say the same thing whenever it
is drawn, and independent of which pass happened to run last. It is also
where plot_be_zones.py already reads from, so the two agree by
construction.

The locked ends are zeroed here for DRAWING ONLY -- D1 and D90 carry rates
about ten times the interior because they are boundary absorbers rather
than sediment budgets, and leaving them in flattens every real feature into
the axis. They are not part of this figure's subject.
```

**`print_be_dicts()`**

```text
Print ready-to-paste DOMAIN_BE_RATES dicts for all scenarios
and optionally write to a txt file.
```

</details>

### 3-figures/plot_be_convergence.py

Did the source/sink calibration converge, and was its zone set fixed in advance?

From the script's original header:

```text
Did the source/sink calibration converge, and was its zone set fixed in advance?

Those are the two questions a reader has to be able to answer, because the
calibration is a FIXED-POINT SOLVE rather than a closed-form one. Where it
stopped is a scientific claim -- "this is the model's limit" -- and a stopping
point is only meaningful if the sequence was contracting and the target was not
moving while it ran.

WHY ITERATE AT ALL
    The ordinary calibration measures the residual of a base run and imposes it
    as the background-erosion field, which assumes that giving a domain X m/yr
    moves that domain's shoreline rate by X m/yr. It does not: BRIE diffuses an
    imposed rate alongshore and the domain keeps only a fraction of it. Measured
    here, one pass closes 42% of the misfit in period 1 and 57% in period 2 --
    so a one-shot residual mixes "the model cannot reproduce this" with "the
    correction was only half applied", and no reader can separate them.
    Iterating removes the second, leaving a residual that means one thing.

    Iterating also needs no estimate of the surviving fraction g, which matters
    because g is not a constant: a contiguous same-signed block of corrections
    passes at ~0.8-1.2 while a pattern alternating at the grid scale is damped
    to ~0.1. Dividing by g instead would amplify narrow features roughly tenfold
    into rates that are indefensible read as sediment fluxes.

WHY THE ZONES ARE FROZEN, AND WHY THE LEFT PANEL SHOWS THE RUN THAT SCORED
BETTER
    Zone membership is the scientific step: it says this stretch of coast has a
    real sediment-budget deficit and here is the process. Magnitude is
    arithmetic. Iterating both lets the arithmetic rewrite the science, because
    each pass re-derives zones from a NEW residual -- so as coherent features
    are satisfied, less coherent ones cross the threshold. And since adding BE
    at a domain pushes sediment into its neighbours, later passes partly correct
    the spillover of earlier ones, which never terminates.

    The unmasked run is drawn because it scored BETTER (dashed, right of the
    converged points). Hiding it would be the wrong kind of tidy: the gap is the
    fit available only by correcting outside justifiable zones, and the argument
    for this calibration is that the gap was declined deliberately, which the
    reader can only weigh by seeing its size.

WHAT THE RIGHT PANEL IS FOR
    To show the zone set was fixed BEFORE the iteration ran, not grown to fit.
    D5-D7 are marked separately: they are the groin's own footprint, reserved so
    the source/sink field cannot absorb the groin's shortfall and double-count
    against the M/f fit. D6 carries the largest residual in both periods and is
    deliberately never corrected.

Usage:
    python 3-figures/plot_be_convergence.py

Reads  2-calibrate/1984_2004__2004_2024/convergence_history.json, and the live FROZEN_ZONE_DOMAINS /
       GROIN_RESERVED_DOMAINS / HATTERAS_BE_RATES_CALIBRATED, so the figure
       cannot drift from the calibration it documents.
Writes data/hatteras_init/7-source-sink/3-figures/1984_2004__2004_2024/2-method/fig_be_convergence.png (and the
       PDF beside it); the caption is written to CAPTIONS.md in that folder.
```

Notes that were in the code:

```text
Reads and writes the data tree, where the fit now puts its products: the
default pair's folder (resolved by hat_source_sink.py since 2026-09-18).
```

```text
The figure belongs with the rest of the section 7 figures, in the data tree;
the iteration's own record stays beside the calibration that wrote it.
```

```text
Labels sit at SEGMENT MIDPOINTS, not on the markers. A gain belongs to
the step, not the endpoint, and at the markers the two periods' labels
collided with each other and with the lines.
```

```text
The earlier period runs BELOW its line and the later one above,
so the two sets of gains cannot meet in the middle.
```

```text
Scaled to the SEQUENCE. edgeBE and zeroBE are 2-4x these values and drawing
them as lines squashed the whole iteration into the bottom fifth of the
panel, which defeats the point of the figure; they are in the caption.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_calibration()`**

```text
The live constants, imported rather than copied.

A figure that hardcodes the zone set would keep rendering happily after
someone edited the calibration, which is the failure mode this exists to
guard against.
```

</details>

### 3-figures/plot_be_zones.py

Which domains were eligible for correction, and how much each one got.

From the script's original header:

```text
Which domains were eligible for correction, and how much each one got.

Two questions a reader of the calibrated field has to be able to answer, and
neither is visible in a table of 90 numbers:

    WHICH DOMAINS QUALIFIED
        A 0.0 in the field is ambiguous on its face -- it can mean "no residual
        here" or "this domain was never eligible". Those are opposite claims.
        The top panel separates them: coloured by physical zone where the
        domain was inside the frozen zone set, grey where it was withheld
        however large its residual, orange at D5-D7 where the groin owns the
        misfit, purple at the two locked ends.

    HOW MUCH CORRECTION IT GOT, AND FROM WHERE
        The lower panels split each final rate into the ONE-SHOT solve and what
        the ITERATION added on top. That split is the case for iterating at
        all: if the light bars were most of the height, the ordinary
        calibration was already converged and the extra passes bought nothing.
        They are not -- the iteration roughly doubles the field, because
        imposing X m/yr of background erosion moves a domain's rate by far less
        than X once BRIE has diffused it alongshore.

WHY ZONE MEMBERSHIP IS FIXED
    Zone identification is the scientific step -- this stretch has a real
    sediment-budget deficit, here is the process. Magnitude is arithmetic.
    Iterating both lets the arithmetic rewrite the science: each pass
    re-derives zones from a new residual, so as coherent features are satisfied
    less coherent ones cross the threshold, and since adding background erosion
    at one domain changes its neighbours' residuals, later passes start
    correcting the spillover of earlier ones. So zones are identified once and
    held (`FROZEN_ZONE_DOMAINS`).

Usage:
    python 3-figures/plot_be_zones.py

Reads  the live FROZEN_ZONE_DOMAINS / GROIN_RESERVED_DOMAINS / PHYSICAL_ZONES,
       the calibrated field from hatteras_site_config.py, and the pass-0 field
       from the masked iteration's first backup -- so the figure cannot drift
       from the calibration it documents.
Writes data/hatteras_init/7-source-sink/3-figures/1984_2004__2004_2024/2-method/fig_be_zones_and_corrections.png
       (and the PDF beside it); the caption goes to CAPTIONS.md in that folder.

REGENERABLE AGAIN SINCE 2026-09-14. It was not, for three weeks: the pass-0
field came from `scripts/hatteras_site_config_prebe_20260824_223143.py`, which
was never committed and is not in the tree, and no later backup could stand in
because each belongs to a different lineage. `convergence_history.json` records
only the per-pass RMSE, not the per-domain field, so it cannot substitute
either.

The frozen-zone correction on 2026-09-14 re-ran the masked iteration from pass
0 and kept every backup it wrote, so the split is recoverable from this
lineage's own pass-0 file (see PASS0_BACKUP below). The lesson stands: the
pass-0 backup is the only record of the one-shot half, nothing reconstructs it
after the fact, and it must be kept with the field it produced.
```

Notes that were in the code:

```text
THE PASS-0 FIELD, AND WHY IT IS THIS FILE AND NOT THE ONE BESIDE IT.
The apply step writes its backup BEFORE it writes, so a `prebe` file holds the
field as it stood going INTO that pass, not coming out. The one-shot solve is
therefore the backup taken before PASS 1, not the one before pass 0 -- that
earlier file holds the superseded field the re-derivation replaced.

..._175853  the 2026-08-24 field, retired  (46 nonzero P1, 65 P2)
..._180700  pass 0, the one-shot solve     (43 nonzero P1, 63 P2)  <- this
..._181336  pass 1
..._182007  pass 2

Re-pointed 2026-09-14. The previous target, ..._20260824_223143.py, was the
pass-0 backup of the retired lineage and was never committed, which is why
this figure could not be drawn and why be_pass0_* / iteration_added_* were
empty in the exported CSV. The corrected iteration kept its backups, so the
split is recoverable again. Do not delete these four files.

Moved 2026-09-18 from scripts/ into the data tree. They are the calibrate
step's output -- a snapshot of the BE field, not code -- and sitting loose at
the root of scripts/ they read as stray config copies, which is how the
08-24 one came to be discarded. They now file beside the rest of 2-calibrate.
One definition, in hat_source_sink.py, which the export reads too.
```

```text
---- TOP: eligibility, one row per period ----------------------------
The fill says what happened to the domain, not which zone it is in: the
zones are named on the ruler beneath, and seven categorical colours here
would collide with the two period colours the bars below depend on.
```

```text
Zone extents named once on a ruler under the strip, in the short forms
the analysis module keeps for exactly this, and on two staggered rows:
set on one row the neighbouring names overlap wherever a zone is narrow.
```

```text
Stacked in the SAME direction as pass 0, so the bar reads as a total
rather than a difference; where the iteration reversed a sign the
segment simply crosses zero, which is itself worth seeing.
```

```text
No village bands here: panel (a) sits directly above on the same
x-scale and carries them, and a second grey could not be told from
the wash that marks the withheld domains.
```

### 3-figures/plot_groin_reserved_residual.py

Why the largest residual in the hindcast, at D6, is deliberately left uncorrected.

From the script's original header:

```text
Why the largest residual in the hindcast is deliberately left uncorrected.

At convergence the source/sink calibration leaves its biggest misfit at D6 --
2.00 m/yr in period 1 and 2.59 m/yr in period 2, roughly twice the next worst
domain in either. That looks like a calibration failure and is not one. D5-D7
are the Buxton groin's footprint, held in GROIN_RESERVED_DOMAINS, and the
residual there is the GROIN's shortfall rather than a background-erosion term.

WHY IT WOULD BE WRONG TO CORRECT IT
    The groin's trapping rate M and deterioration floor f were fitted against
    the observed shoreline, and the source/sink field is then derived from what
    the modules could NOT explain -- which is why the calibration runs against a
    groin-ON base run in the first place (GROIN_AWARE_BASE_RUN). Letting BE
    absorb the residual at D5-D7 would close the same gap twice: the groin would
    score as well-calibrated because a source term was quietly doing its work,
    and the M/f fit could never be falsified by the hindcast.

    So the number stays visible. It is the honest statement of what the groin
    module cannot do.

WHAT THE TWO SIGNS MEAN, AND WHY THEY ARE OPPOSITE
    period 1   residual POSITIVE -- observed is more seaward than modelled.
               The model does not build enough fillet. M = 60 m/yr is the most
               the sediment budget will support (719,000 m3/yr against a
               5-7e5 littoral drift), so this is a bound, not a missed fit.

    period 2   residual NEGATIVE -- modelled is more seaward than observed.
               The real fillet RELEASED after the 2003 storm damage; the module
               cannot, because trapping is bounded at >= 0, so it can stop
               adding sand but never remove it. This is outside the
               parameterisation at any (M, f), not a badly chosen one.

    The opposite signs are the point. A source/sink term fitted to close both
    would have to change sign between periods at the same domain, which is a
    fitted constant standing in for a structure that was built, damaged and
    left -- exactly the kind of thing the zone rules exist to keep out.

Usage:
    python 3-figures/plot_groin_reserved_residual.py

Reads  the converged calibBE full_management runs, groin on and off, plus the
       live GROIN_RESERVED_DOMAINS.
Writes data/hatteras_init/7-source-sink/3-figures/1984_2004__2004_2024/3-limits/fig_groin_reserved_residual.png
       (and the PDF beside it); the caption goes to CAPTIONS.md in that folder.
```

Notes that were in the code:

```text
RESOLVED, NOT JOINED: the rate CSV is tables/shoreline_change_rate.csv in
the new run layout and {run}_shoreline_change_rate.csv in the old one.
```

```text
The gap itself is drawn; its size is a number, and numbers belong in
the caption and on the bars at the right, not floating over a curve.
```

```text
The reserved reach is the groin module's own footprint, so it takes
the accent tint rather than a grey: the villages already occupy the
grey strip along the top, and two greys on one panel cannot be told
apart. Named once, on the upper panel only.
```

```text
groin-off as an outline behind: the gap between the two IS the groin.
One legend entry only: the two periods' outlines are visually
identical, so labelling both just doubles the legend.
```

```text
Room under the bars for their value labels, and a clear band above them
for the key: at "lower right" the key landed on the D6 label.
```

### 4-export/export_be_calibration.py

Export the converged source/sink field to data/hatteras_init/7-source-sink, with a README.

From the script's original header:

```text
Export the converged source/sink field to data/hatteras_init/7-source-sink.

WHY THIS EXISTS
    `hatteras_site_config.py` is the source of truth -- it is what the runs
    import -- but it is a Python module in the scripts tree, which is a poor
    place to look for "what were the calibrated values". The data directory
    already held per-period copies, and they had drifted: written 2026-06-15,
    they carry GIS 1 = -40 against the current -41.8, zeros across D2-D11 that
    are now +1.4 to +2.6, and the 2004 file was truncated mid-dict. The config
    itself carries a comment warning readers not to trust them.

    So this regenerates them FROM the config rather than alongside it, and adds
    the context a bare dict cannot carry: which domains were eligible to be
    corrected at all, which were withheld and why, and how much of each final
    value came from the one-shot solve versus the iteration.

WHAT IT WRITES (paths from site_layer/hat_source_sink.py since 2026-09-18)
    4-export/be_rates_<period>.py        the dict, same shape as the files
                                         it replaces
    4-export/be_calibration_domains.csv  one row per domain: zone,
                                         eligibility, the pass-0 and final
                                         rate, what the iteration added,
                                         and the residual still standing
    README.md                            at the top of 7-source-sink/:
                                         provenance, and the caveats that
                                         matter

    It READS the default pair's 2-calibrate/1984_2004__2004_2024/ and
    3-figures/1984_2004__2004_2024/. The figures are written straight there
    by stage 3, so the data directory is self-contained -- someone handed
    just this folder can see what was done, not only what came out.

    convergence_history.json is NOT copied up. It lives once, in the pair's
    2-calibrate/ folder, beside the calibration that wrote it -- a second
    byte-identical copy at the top level gave one fact two owners, and the
    two would drift the first time a pass was re-run without re-exporting.

    The superseded 2026-06-15 files are MOVED to archive/superseded_<date>/
    rather than deleted -- they are what earlier runs were built against, so
    they are history, not clutter.

Usage:
    python scripts/input_prep/7-source-sink/4-export/export_be_calibration.py [--check]
```

Notes that were in the code:

```text
Resolved by hat_source_sink.py since 2026-09-18. The README stays at the
top of 7-source-sink/; the exported field goes to 4-export/, the step that
writes it; the calibration and figures read are the DEFAULT pair's folders.
```

```text
The calibration products moved into the data tree 2026-09-12, so this
reads them from there rather than from beside the script that made them.
```

```text
WHERE THE FIGURES COME FROM (corrected 2026-09-12). This used to read them
out of the calibration OUTPUT directory and copy them into 3-figures/ --
but be_zone_residual_fit.py writes its figures straight there, and
the copies in the output directory were older. So the export overwrote the
CURRENT figures with SUPERSEDED ones, quietly, every time it ran.

3-figures/ is the record now, and the staleness check below reads the same
place. The superseded copies were moved under superseded_20260825/ (in
data/hatteras_init/7-source-sink/archive/, deleted 2026-10-01).
```

```text
Paths RELATIVE TO the pair's figure folder (FIGURE_DIR,
3-figures/1984_2004__2004_2024/), which is why each carries its subfolder.
The figures were grouped on 2026-09-14 by the question each answers:
1-field is what the calibration produced and what it was fitted against,
2-method is the evidence that the way it was produced holds up, and
3-limits is what it deliberately does not do. A flat folder gave those
three equal weight, and the limits figure is the one most often read as a
calibration failure rather than a stated boundary.
```

```text
The masked-iteration lineage. The unmasked attempt (backups 214732, 220039)
is deliberately NOT here: it was abandoned for correcting outside the
geomorphological zone set, and mixing its values into the record would
misrepresent what was shipped.
The one-shot solve, for the be_pass0_* / iteration_added_* split. The apply
step backs up BEFORE it writes, so the pass-0 field is the backup taken before
PASS 1 -- the file before pass 0 holds the superseded field, not this lineage.
Re-pointed 2026-09-14 with the corrected iteration, whose backups were kept;
the previous target was the retired lineage's pass-0 file and was never
committed, which is why these two columns exported empty.
This pointed at scripts/, where the backup has not been since it was filed
under 2-calibrate/prebe/, so the next export would have written the pass-0
columns empty. One definition now, shared with plot_be_zones.py
(2026-09-18).
```

```text
Depth-agnostic on purpose. This was a fixed-depth
"*/calibBE/*/*_shoreline_change_rate.csv" glob, which reaches
<period>/calibBE/<run>/ but NOT the <arm>/<period>/calibBE/<run>/ a run
forced off the calibration wave climate is filed under. A missed run can
only ever make `newest` OLDER, so the failure is the staleness check
passing a figure it should have caught -- silently, and in the direction
that publishes the stale file.

Runs are found by their metadata file, which stays at the run folder's
root under its full name; the rate CSV itself has moved into tables/ and
dropped the prefix, so globbing for it would miss a migrated run -- again
in the direction that publishes the stale file. run_layout.resolve reads
either layout.
```

```text
Computed, not typed: these read 42/57 for the lineage retired on
2026-09-14 and were still saying so after the field had been re-derived.
A number stated in prose beside the table it contradicts is worse than no
number at all.
```

```text
THE PASS-0 FIELD IS OPTIONAL NOW (2026-09-12). It used to be required,
which made this exporter unrunnable: the backup it names was never
committed and is not on disk, so the two columns it feeds cannot be
reconstructed. Refusing to export at all meant the VALUES stayed stale
for the sake of a provenance column -- and stale values are the failure
this file exists to prevent. Absent, the pass-0 and iteration-added
columns are written empty and the README says why.
```

```text
BEFORE anything is written, and before --check returns, so a dry run
reports the same refusal a real one would. Checking it later left a
half-finished export on disk: values written, figures not.
```

```text
3-figures/ IS the record, so there is nothing to copy into it -- the
analysis scripts write here directly. This block now only reports what
is present, which is what the README claims.
```

```text
Named rather than skipped silently: a figure absent from the record
is indistinguishable from one that was never made.
```

```text
convergence_history.json is NOT copied up (dropped 2026-09-14). It lives
once, in 2-calibrate/, beside the calibration that wrote it. The copy
that used to sit here was byte-identical the day it was made and would
have diverged the first time a pass ran without a re-export -- at which
point a reader has two files with one name and no way to tell which the
figures were drawn from. The README points at the one copy instead.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_never_die_on_a_print()`**

```text
Stop a console encoding from killing a finished export.

This file prints en- and em-dashes, which a Windows cp1252 console cannot
encode, so `print` raises UnicodeEncodeError. Its sibling
be_zone_residual_fit.py lost a completed calibration to exactly that
on 2026-08-28 -- the crash landed between computing the numbers and
writing them out. Reconfigure rather than ASCII-ify, so the next dash
someone types cannot reintroduce it.
```

**`stale_against_runs()`**

```text
Which of `paths` predate the newest calibBE run output.

The exporter COPIES whatever is on disk; it does not regenerate. So a
figure left over from before the last hindcast would be published into the
data directory looking exactly as authoritative as a current one. This is
the check that catches that -- it is the failure mode that actually
happened on 2026-08-25, when fig_be_convergence.png was three minutes older
than the run that had just been rebuilt.
```

**`caveats_block()`**

```text
What this export could NOT carry, and why.

Generated rather than hand-written into the README, because the README is
rewritten on every export and a hand-added note would vanish on the next
run -- which is how the file came to disagree with the config in the first
place.
```

</details>
