# tools

Read or repair the run record; none of them run the model. Start with
`HAT_period_input_check.py` before a run and `HAT_list_runs.py` after one.

| Script | What it answers |
|---|---|
| `HAT_period_input_check.py` | is this period runnable, and what is missing |
| `HAT_list_runs.py` | what is on disk, grouped |
| `HAT_run_supersession_report.py` | which runs a later one has superseded |

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_list_runs.py

What is under output/raw_runs, grouped by whatever makes runs comparable.

From the script's original header:

```text
What is under output/raw_runs, grouped by whatever makes runs comparable.

WHY A LISTING RATHER THAN A FOLDER
    98 of 164 runs sit on a topography that is no longer CURRENT, and opening
    a preset folder does not say which. The obvious fix is to nest by
    topography version, and it is the wrong one: 2004-start has only one
    version, so 63 runs would gain a level that says nothing; a sweep run is
    already five levels down; and topography is not the only axis that decides
    whether two runs are comparable, so nesting one of them just moves the
    question.

    Every run records all of it -- topo_product, topo_dune_version,
    be_values_digest -- in its metadata and in run_index.csv. So the question
    "which runs are on v1" is a query, and this is the query.

    `--stamp` writes the same facts as a one-line BUILT_ON.txt inside each run
    folder, so a folder opened on its own also answers it.

USAGE
    python HAT_list_runs.py                    # by topography (the default)
    python HAT_list_runs.py --by erosion
    python HAT_list_runs.py --by preset --detail
    python HAT_list_runs.py --only-stale       # just what is not CURRENT
    python HAT_list_runs.py --stamp            # write BUILT_ON.txt per run
```

<details><summary>Function notes (the original docstrings)</summary>

**`newest_digest_per_preset()`**

```text
The digest the most recent run of each (preset, period) carries.

PER PERIOD, not per preset: the background-erosion field is calibrated for
one period at a time, so one preset legitimately has a different digest in
1984-2004 than in 2004-2024. Keying on the preset alone called that
by-design difference "superseded".
```

**`parent_bucket()`**

```text
The folder a run sits in, collapsed to what says its purpose.

Purpose layout (2026-09-16): matrix/<period>/<preset>, sensitivity/<axis>,
experiments/<tag>, versions/<tag>. The 09-10 layout's sweeps/<family>
collapses to 'sweeps' as before.
```

</details>

### HAT_period_input_check.py

Does every hindcast period resolve everything it needs, and does each file exist on disk?

Notes that were in the code:

```text
HAT_period_input_check.py

Does every hindcast period resolve everything it needs, and does each of
those things exist on disk?

WHY THIS EXISTS
HATTERAS_PERIODS names six forcing inputs per period as RELATIVE PATHS.
Nothing checks them until a run is most of the way through section 5, and
two of the six are read later still -- so a period wired against a file
that was never built fails partway into a run that has already spent
minutes building an island. Four periods since 2026-09-11, two of them
wired ahead of inputs still being digitised, made that a certainty rather
than a risk.

It is also the answer to "is this period runnable yet", which is otherwise
answered by starting a run and waiting.

WHAT IT DOES NOT DO
It does not validate VALUES. A setback of the right shape measured against
the wrong topography is a real failure this cannot see; that is what
HAT_road_setback_audit.py and the extractor audits are for. This checks
resolution and existence, the class of failure that costs a run rather
than a result.

EXIT CODE
0 when every period is runnable, 1 when any is blocked, so a batch driver
can gate on it.

python HAT_period_input_check.py
python HAT_period_input_check.py --period 2010

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-18
```

```text
Anchored by SEARCHING UPWARD for the project root rather than by
counting parent directories (2026-09-13). A counted depth is correct
only while the file stays where it was written, and these moved into
subfolders of hatteras_ms. Six files here already did it this way.
```

```text
Resolved, not built: hat_observed_rates owns where the rate fits live, so
this check cannot look somewhere the runner does not.
```

```text
Not a failure, but worth seeing: a zero setback puts the road on the dune
line, where one cell of retreat re-fires a relocation.
```

<details><summary>Function notes (the original docstrings)</summary>

**`check_presets()`**

```text
Which source/sink presets are solved for this period.

Not blocking: zeroBE alone is enough to run. A period missing edgeBE simply
cannot be run under it, which is a fact about the calibration rather than
about the inputs.
```

**`events_in_window()`**

```text
What the record fires inside this window.

The window is start..end-1, matching run_cascade_simulation: an event dated
exactly on the end year belongs to the next period and never fires here.
```

</details>

### HAT_run_supersession_report.py

Which runs under output/raw_runs are candidates for retirement, and why (read-only).

From the script's original header:

```text
Which runs under output/raw_runs are candidates for retirement, and why.

READ-ONLY. It moves nothing and deletes nothing. Hannah's instruction on
2026-09-10 was "flag candidates, retire nothing yet", and a script that can
only report cannot be misread as one that acts.

WHAT "SUPERSEDED" MEANS HERE, AND WHAT IT DOES NOT
    A run is flagged when something it was built on has since moved, which
    makes it incomparable with a run made today. That is not the same as
    wrong: a run is a faithful record of the inputs it had, and the two
    1984-2004 archives were kept for exactly that reason. The judgement about
    whether an incomparable run is still worth keeping is the reader's.

THE CHECKS
    topography   the run's dune-topo version against the product's CURRENT.
                 The 2026-09-03 re-pick moved 1984-start from v1 to v2, so a
                 v1 run measures a different island from a run made today.
    calibration  the source/sink `values_digest` within one preset AND ONE
                 PERIOD. The field is calibrated per period, so a preset having
                 a different digest in each period is by design; two digests
                 inside one period would mean it was recalibrated and only some
                 runs remade. (Grouping on the preset alone reported the
                 by-design difference as 55 superseded runs on 2026-09-10.)
    duplicates   one run name in more than one place. Usually legitimate --
                 an arm is a different forcing of the same scenario -- so
                 these are listed for orientation, not flagged.

WHERE   output/raw_runs/SUPERSEDED_CANDIDATES.md

USAGE
    python HAT_run_supersession_report.py
    python HAT_run_supersession_report.py --print
```

Notes that were in the code:

```text
An arm component that names a dune-topo version, e.g. the v3 of
arms/version-pair/v3.
```

```text
ARMS THAT NAME A VERSION are judged against that version, not CURRENT: they
exist to compare two islands, so their state is kept whatever CURRENT says.
```

```text
--- calibration ------------------------------------------------------
WITHIN A PERIOD. The background-erosion field is calibrated per period,
so one preset legitimately carries a different digest for 1984-2004 than
for 2004-2024. Grouping on the preset alone reported that by-design
difference as 55 superseded runs on 2026-09-10, which was wrong; the
split that would matter is two digests for one preset in ONE period.
```

```text
4
---- model state ----------------------------------------------------
```

<details><summary>Function notes (the original docstrings)</summary>

**`expected_version()`**

```text
The version this run SHOULD be on, which is not always CURRENT.

An arm can name a version -- `version-pair/v3` holds v3 against v2, and
`behindroad-copy` was built on the v3 footprint layer. Those runs are on
that version deliberately, so measuring them against CURRENT reports a
deliberate choice as drift, and acting on it would destroy the comparison
the arm exists for (flagged 2026-09-11).
```

**`_state_bytes()`**

```text
Size of this run's model state, or 0 if it was not saved.

The .npz is the ONLY artifact that lets a deep-dive figure be re-derived
without re-running. Everything else -- the rate table, the shoreline
matrix, the metadata -- is written regardless and is small.
```

**`state_verdict()`**

```text
Does the keep rule keep this run's model state?

THE RULE (Hannah, 2026-09-14): keep the state of runs on the CURRENT
topography, and of arms that name a version deliberately. Drop the rest.

Why it is defensible: a run on an island a re-pick behind cannot be
compared against a run made today, so re-plotting it deeply is answering a
question nobody can ask. Its tables and metadata survive either way.

Why it is not automatic: a run is one to three minutes of compute, but
re-running reproduces it at TODAY'S code, not the code it was made with.
That is what the state is really insuring against, and it is why nothing
here deletes anything.
```

</details>

## Retired 2026-10-23

Removed when the hindcast moved to CASCADE-obx: `HAT_migrate_run_layout.py` (the one-off 2026-09-16
move into the current run layout, done), `HAT_rerun_arm.py` (re-ran the stored matrix runs with
relocations on, for one study; the matrix it read was archived) and `HAT_compare_rerun.py` (compared
those re-runs). Recover any of them with

    git show a7c4358c:scripts/real-world-cases/Hatteras/hatteras_ms/tools/<script>.py

Their sections, as they were:

#### HAT_compare_rerun.py

What changed between a stored run and the same run made under today's code?

Notes that were in the code:

```text
HAT_compare_rerun.py

What changed between a stored run and the same run made under today's code?

Pairs each run in an ARM against the stored run of the same name, and reports
the differences that matter: the skill metrics, whether a road drowned, how
many relocations fired, and the largest per-domain rate change.

WHY IT IS NOT ENOUGH TO COMPARE SKILL. An island-wide RMSE can sit still
while a road drowns: drowning stops roadway management for that domain from
that year on, which changes what happens to the interior without necessarily
moving the shoreline much. So the road table is differenced too.

python HAT_compare_rerun.py --tag code-checks/2026-09-14-relocation-arm-rerun-new-code

Since 2026-09-16 a re-run is an EXPERIMENT (raw_runs/experiments/<tag>/) and
the stored run a MATRIX row; the index is keyed on (run_name, kind, tag).

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
```

#### HAT_migrate_run_layout.py

Move output/raw_runs into the purpose layout of 2026-09-16 (moves only), then rebuild the index.

From the script's original header:

```text
Move output/raw_runs into the PURPOSE layout of 2026-09-16. Moves only;
nothing is deleted, rewritten or renamed in content, and run_index.csv is
rebuilt from the runs afterwards rather than edited.

THE MOVES
    <period>/<preset>/<run>                       -> matrix/<period>/<preset>/<run>
    <period>/<preset>/sweeps/<axis>/<run>_<tok>   -> sensitivity/<axis>/<period>/<preset>/<run>_<tok>
    arms/<arm>/<period>/<preset>/<run>            -> versions/<tag>/... or experiments/<tag>/...
                                                     as run_registry.LEGACY_ARMS says
    arms/waveHs<x>/1996_2010/<preset>/<run>       LEFT IN PLACE, and listed. Those are
                                                     the twelve 1996 wave cells filed by
                                                     the 09-01 rule; they are re-run as
                                                     sensitivity cells (the token back in
                                                     the name) and then deleted, because
                                                     renaming a run's files and metadata
                                                     is exactly the in-content edit this
                                                     tool refuses to make.

    Why (Hannah, 2026-09-16): a wave sweep had fanned out into twelve
    top-level arms; nineteen of thirty arms were finished one-off experiments
    nothing marked as finished; and the layout could not tell those from the
    version comparisons that are kept on purpose. A run's folder now says what
    it is FOR. The whole design is in output/raw_runs/README.md.

WHAT ELSE IT WRITES
    experiments/<tag>/NOTE.md   one per experiment set, seeded from the index
                                and the write-ups that already name these
                                runs. Says what was asked, where the answer
                                is, and whether the runs may be deleted. Not
                                overwritten if present.

WHY IT IS SAFE TO RUN, AND TO INTERRUPT
    Every destination is checked to be free before anything moves, so a name
    collision is reported and nothing happens. All three layouts are readable:
    run_registry.find_run_dir tries the purpose path, then the 09-10 one, then
    the flat one, so a tree that is half moved -- or interrupted here -- still
    reads. Running it twice is a no-op.

    The 2026-09-10 pass (files inside each run into figures/, animations/,
    tables/) is retained as --files, and is a no-op on a migrated run.

USAGE
    python HAT_migrate_run_layout.py --dry-run       # show every move
    python HAT_migrate_run_layout.py
    python HAT_migrate_run_layout.py --files         # the 09-10 in-run file pass
```

Notes that were in the code:

```text
What each experiment set was for, from the run index and the write-ups
that cite it. Seeded here so the NOTE.md a folder gets on migration is not
blank; edit the file, not this table, once the folder exists.
```

<details><summary>Function notes (the original docstrings)</summary>

**`tree_moves()`**

```text
(source, destination) for each run folder that changes place, and the
folders that stay with why.
```

</details>

#### HAT_rerun_arm.py

Re-run a set of existing runs under today's code, into an arm, so the two can be compared.

Notes that were in the code:

```text
HAT_rerun_arm.py

Re-run a set of existing runs under TODAY'S code, into an arm, so the old
results survive and the two can be compared.

WHY
A run records the git commit it was made at, and the index shows those
commits drifting apart. The relocation arm of the 1984-2004 matrix was made
on 2026-09-01 from a dirty tree; spot-checking one of its cells on current
code drowned NC-12 at GIS 11, where the stored run reports none. That is
either a real change in the model or a change in the inputs, and the only
way to tell which runs is to re-run and difference.

PINNING THE TOPOGRAPHY IS HALF THE JOB, and the half that is easy to miss.
A road setback is metres landward of interior row 0, so it belongs to the
extraction it was measured on. `--topo-version v1` pins the ARRAYS but the
run still reads whatever setback file the period table names, which is the
v2-era measurement. Those two differ at 27 of 82 domains, up to 205 m, and
the mismatch moves every domain's rate by up to 0.005 m/yr.

That is exactly what happened on the first use of this script: a twelve-cell
comparison that looked like code drift was mostly a mismatched pair. Pinning
a version means pinning BOTH halves; the archived setbacks are under
road_offset/superseded_<date>/.

IT WRITES INTO AN ARM, NEVER OVER THE ORIGINAL. An arm scopes the output
directory, so the stored runs and their index rows are untouched and the
comparison is reversible. Promoting a re-run to the production path is a
separate, deliberate act.

python HAT_rerun_arm.py --list
python HAT_rerun_arm.py --arm recode-20260914 --topo-version v1
python HAT_rerun_arm.py --arm recode-20260914 --topo-version v1 --limit 2

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-27
```

```text
The switches a run name encodes, recovered from the index columns rather than
parsed out of the name: the name is DERIVED from the switches, so reading it
back would be inverting a lossy function.
```
