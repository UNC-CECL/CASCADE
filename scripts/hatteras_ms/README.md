# hatteras_ms — the hindcast, and everything that reads its output

Four kinds of thing live here, and until 2026-09-13 they were one flat list of
31 items sharing a `HAT_` prefix that sorted them without separating them.

```
HAT_hindcast_1984_2024.ipynb   THE RUN. The notebook is authoritative;
HAT_hindcast_1984_2024.py      the .py is its headless mirror, with no
                               features of its own. A change goes in BOTH.
HAT_hindcast_config.py         which run happens, and where the value came
hat_run.yaml                   from: env > yaml > the default in the module
HAT_run_all.py                 the batch driver: the matrix, without and with the groin
HAT_hindcast_methods.md        the written method
DEM_TO_DEM_CALIBRATION.md      the calibration/test plan (1996-2009, 2009-2025):
                               what is done, the open decisions, where to resume
HINDCAST_PLAN.md               planning notes: the order the notebook builds
                               the run in. Was `HAT_hindcast_plan`, a file
                               with no extension that this line described as
                               a folder (2026-09-22)

tools/        read or repair the run record; none of them run the model
experiments/  one-off studies, each a run-it then plot-it pair
```

## Which Barrier3D the run uses

Barrier3D is a separate repository, installed editable, so the branch checked out at `../Barrier3D` is the model. Every run records the branch and commit it imported (`run_registry.barrier3d_provenance`). None of the fixes below has been pushed upstream (Hannah, 2026-09-28).

**To do: push the Barrier3D fixes eventually** (Hannah, 2026-09-29). They stay local for now. Until they are pushed, this branch runs only on a machine that has `hatteras/adopted` checked out at `../Barrier3D`: the runner refuses a Barrier3D without per-cell ceilings. The full record, with evidence paths, is `HATTERAS_FIXES.md` on `hatteras/adopted`. Changes to CASCADE's own `cascade/` package are recorded in `HATTERAS_CASCADE_CHANGES.md` at the repository root.

| branch / tag | what it has | in use? |
|---|---|---|
| **`hatteras/adopted`** (2b8f167 code; 8a588ea records it) | everything below merged: the route_overwash fix, the three overwash fixes, and per-cell dune ceilings (`DuneCeilingFromStart`, switched on in `data/hatteras_init/Hatteras-CASCADE-parameters.yaml`) | **yes, since 2026-09-28**: `../Barrier3D` is on it. The runner refuses a Barrier3D without per-cell ceilings. The matrix before this is in `output/raw_runs/archive/2026-09-28-pre-ceiling/`. |
| `feature/per-cell-dune-ceiling` (d343461) | per-cell dune ceilings alone, on 49fd069; worktree `../Barrier3D-dune-ceiling` | merged into `hatteras/adopted` |
| `fix/overwash-gaps-momentum` (db0ba30, tag `hat-fix-overwash-gaps-momentum`; worktree `../Barrier3D-overwashfix`) | `DuneGaps` dropping cells; the gap discharge slice; the inundation momentum constant reset to 0 | merged into `hatteras/adopted` |
| `fix/route-overwash-axis-swap` (49fd069, tag `hat-fix-route-overwash`) | the `route_overwash` index swap: wrong cells read, and a silent crash on long storms | merged; every run from 2026-09-24 to 09-28 used it alone |

The storm series moved to `v3_trim24` on the same day (`hat_env_forcings.DEFAULT_STORM_VARIANT`). Every run records its storm file and dune-ceiling mode in its metadata. To reproduce a run made before 2026-09-28, check out `fix/route-overwash-axis-swap` and use that date's parameter template and `v3_72` storms.

`../Barrier3D-prefix-ce36866` is a detached worktree of the code **before** the route_overwash fix. It is kept only for the storm-duration cause test (`output/raw_runs/experiments/storms-overwash-and-dunes/2026-09-28-storm-max-duration/`) and can be removed with `git -C ../Barrier3D worktree remove ../Barrier3D-prefix-ce36866`.

## tools/

| Script | What it answers |
|---|---|
| `HAT_period_input_check.py` | is this period runnable, and what is missing |
| `HAT_index_runs.py` | rebuild `run_index.csv` from the runs on disk |
| `HAT_list_runs.py` | what is on disk, grouped |
| `HAT_run_supersession_report.py` | which runs a later one has superseded |
| `HAT_migrate_run_layout.py` | move runs to the current directory layout |

Start with the period check before a run and `HAT_list_runs.py` after one.

## experiments/

Two live studies, plus one retired.

* **the crest experiment** -- `HAT_run_crest_experiment.py` and its plotter.
* **the relocation set** -- `HAT_relocation_comparison.py` (takes `--period`
  since 2026-09-15: 1984 or 1996, one output root per window),
  `HAT_relocation_period_compare.py` (the 1999 event read across both
  windows, from the per-period tables),
  `HAT_relocation_dune_position_check.py`, `HAT_score_relocation_timing.py`
  and `HAT_score_road_position.py`. `RELOCATION_COMPARISON_RESULTS.md` is
  what it concluded.
* `superseded_20260907/` -- the 1984 seaward row-insert set, which could not
  be re-run (its topography layers, outputs and arms were gone). Deleted
  2026-10-01 with the other retirement folders below; git has the scripts.

A driver spawns the runner as a subprocess with `HAT_IGNORE_SETTINGS=1`, so
whatever is sitting in `hat_run.yaml` cannot reach an experiment.

## figures/ — moved

The five figure scripts that were here (`hindcast_final_figure_lowess`,
`scenario_grid`, `rerender_run_figures`, `planview_evolution_gif`,
`gis11_relocation_drown_figure`) and their `superseded_20260914/` (since
deleted, 2026-10-01) moved to
`scripts/figure_making/model_output/` on 2026-09-18, so every figure script is
in one tree.

## Deleted 2026-10-14

One groin model is kept, the corrected groin (`cascade.groin.CorrectedGroinCallback`), set by
its strength per window and its end year in `hat_run.yaml` (`groin.strength_1996` 0.16,
`groin.strength_2009` 0.10, `groin.end_year` 2017; `groin.strength` overrides both). Removed with
the dipole (`GroinCallback`, M) and blocking (`BlockingGroinCallback`, b) groins and their
deterioration floor f:

* `groin-sweep/` -- the M and f sweeps, the joint fit and their figures
* `HAT_run_all.py` stages 3-5 (seed runs, sweeps, joint fit); the groin matrix, stage 6, is now
  stage 3 and reads the groin from the config, not `joint_fit.json`
* `tests/test_groin_blocking.py` (now `tests/test_groin.py`, for the corrected groin)

The last commit with all of them is `4d937b79`. The dated groin experiments under
`hard-structures/groin/` that build on the old classes ran against that code; check it out to
rerun one.

## Deleted 2026-10-01

Hannah chose to delete the retirement folders rather than keep them (git keeps
them):

* `superseded_20260819/` -- the old 1984-2024 runners
  (`HAT_hindcast_1984_2024.py`, `_newdomains`, `_newplot`) and the inherited
  driver `benton_script.py` (that lineage lives on `origin/benton-cascade`)
* `superseded_20260911/` -- `HAT_hindcast_1984_2024_versioncontrol.py` and
  `HAT_hindcast_QC_skeleton.ipynb`
* `experiments/superseded_20260907/` -- the 1984 seaward row-insert set

To recover one, by the path it was last committed under:

```
git log --diff-filter=D --oneline -- scripts/hatteras_ms/<folder>/<file>
git show <commit>^:scripts/hatteras_ms/<folder>/<file>
```

## Paths: search upward, do not count

Every script here finds the project root by walking up until it sees a marker,
never by counting parent directories:

```python
REPO = next(p for p in Path(__file__).resolve().parents
            if (p / "pyproject.toml").exists())
```

Six files already did this; the other fourteen counted depth and were
converted on 2026-09-13, before the move, so that the move itself could not
break them. Keep it that way — a counted depth is correct only until the file
is filed somewhere better.

The runner stays at the top level because the config module, the settings yaml
and the batch driver all reach it as a sibling. A script that needs it names
`REPO / "scripts" / "hatteras_ms"` rather than its own folder.

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_hindcast_config.py

Run-selecting settings for the Hatteras hindcast, in one place: hat_run.yaml, overridable by HAT_* variables.

From the script's original header:

```text
Run-selecting settings for the Hatteras hindcast, in one place.

WHERE A SETTING IS TYPED
    `hat_run.yaml`, beside this file. Edit it, save it, run the hindcast.
    That file is the interface; this module is the machinery that reads it,
    and its own values are only the fallbacks.

WHY THIS MODULE EXISTS
    `HAT_hindcast_1984_2024.ipynb` and its headless mirror
    `HAT_hindcast_1984_2024.py` used to carry these values as literals in
    sections 1, 3, 7, 9 and 11. Running the scenario matrix therefore meant
    hand-editing a source file between every run, which is how the two
    published edgeBE runs in the retired `run_index.csv` ended up disagreeing
    with each other on the background-erosion values they were fit under.

    Both files now read the values from here, so a run can be selected by
    editing one untracked-by-the-model settings file, or driven through the
    environment without editing anything at all.

PRECEDENCE
    environment variable  >  hat_run.yaml  >  the default in this file

    `HAT_IGNORE_SETTINGS=1` drops the middle term. `HAT_run_all.py` sets it on
    every run it launches, so a batch run is described entirely by the driver
    plus this file's defaults, and a half-finished experiment left in
    hat_run.yaml can never reach the comparison matrix or the sweep.

    `describe()` reports which of the three each value came from, so a run log
    states how it was driven rather than leaving it to be inferred.

RE-READING, FOR THE NOTEBOOK
    `load_run_config()` re-reads hat_run.yaml on every call. Section 3 of both
    files calls it rather than using the module-level RUN_CONFIG, so editing
    the yaml and re-running the cell picks the change up without restarting
    the kernel -- a module-level constant would be cached by the import system
    and the edit would silently not apply.

WHAT IS *NOT* VALIDATED HERE
    Value legality: scenario names, periods and offset modes each have exactly
    one home in the code (`SCENARIOS` in section 3, `HATTERAS_PERIODS`,
    `ISLAND_OFFSET_MODES`), and a second copy here is a second thing to drift.
    A bad value raises in section 3 with the valid list attached. What IS
    validated here is what this module owns: that every key in the yaml is a
    key it knows, and that each value casts to the right type.

THE SYNC RULE STILL APPLIES
    This module is imported by BOTH the notebook and the .py. Adding a field
    here is only part of the change -- the yaml needs a documented key and the
    matching section of both files has to read it, or the two drift again.
    See the module docstring of the .py.
```

Notes that were in the code:

```text
The model .npz, for the preflight cost line. Measured, not guessed: the
saved states dominate and the file is ~160 MB for either period.
```

```text
Each raises ValueError on malformed input. Deliberately fatal rather than
falling back to the default: a driver or a yaml that misspells a value would
otherwise run the default configuration under the name of the one it asked
for.
```

```text
(attribute, yaml path, cast, default)

The yaml path is a tuple so a nested block reads naturally in the file
(`groin.enabled`) while staying a flat attribute in code (`groin_enabled`).
The environment name is always ENV_PREFIX + the attribute upper-cased, and
is NOT derived from the yaml path -- the ten names that existed before this
file did are load-bearing in `HAT_run_all.py` and in shell history, and
renaming them to match a yaml layout would break both.
```

```text
1996 SINCE 2026-09-17 (Hannah: "moving forward the main time period,
for figures and future runs, is 1996-2010-2024"). 1984 was the pair
the project ran first. Every other start is still reachable by naming
it; this is only what a run gets when nothing does.
```

```text
METRES SINCE 2026-09-24 (Hannah): the offset file is metres and so is
BRIE's x_s, so the file goes in as it is. "asrun" (offset / 10, the
units error every run before this date carries) is still reachable by
naming it, so those runs stay reproducible -- with the offset build they
were made from, HAT_OFFSET_VERSION_<year>=superseded_20260924_pre-metres/v1,
since the current v1 closes the buffer differently (the numbering
restarted at v1 that day). The name rule is unchanged: every mode but asrun
earns an `offset<mode>` token, so a metres run can never take the name
of the /10 run it replaces. Study:
output/raw_runs/experiments/island-offset/2026-09-24-div10-vs-metres-wave-sweep/.
```

```text
THE REACH (2026-09-16): a name from hat_extension_domains.GEOMETRIES.
"base" is GIS 1-90. hatteras_site_config reads the same HAT_GEOMETRY
from the environment to build HATTERAS_DOMAINS, and the runner refuses
a run where the two disagree (a yaml `geometry:` with no env var).
```

```text
WHERE A RUN IS FILED (2026-09-16). run_kind is one of run_registry.KINDS
-- matrix (the default), sensitivity, experiment, version -- and run_tag
names the experiment or version. A matrix run has no tag; a sensitivity
cell's tag is derived from its name's sweep token. HAT_RUN_KIND and
HAT_RUN_TAG in the environment; HAT_ARM_TAG is the old spelling and the
runner reads it as an experiment tag.
```

```text
THE GROIN, 2026-10-14. The corrected groin, by two numbers: the share of
the drift arriving at the GIS 5|6 face it holds back (strength, 0.25), and the
first year it no longer works (end_year, 2017, the Buxton fill). Both fitted on
one continuous 1967-2025 run against the wet/dry photo dates and annual
CoastSat (hard-structures/groin/3-hindcast/2-blocking-1996-2025/
2026-10-13-pregroin-planform/). The dipole M, the blocking b and the floor f
went with their groins.
```

```text
WAVE CLIMATE. All four are forcing, not management: they change what is
simulated, so a run that moves any of them off the value below earns a
`wave...` token in its name and lands in its own directory. Without that
token a sensitivity cell would derive the SAME name as the matrix run it
is being compared against, and the last one to finish would be left
wearing the production name -- the failure output/calibration/groin/README.md
documents for the rig sweep.

The three below `hs` were literals in section 11 of the .py until
2026-09-01 (FIXED_WAVE_PERIOD and friends). They are fields now for one
reason: scripts/sensitivity_analysis sweeps them, and a sweep that has to
edit the model source between cells is the hand-editing failure this
module exists to remove.

OPTION A, ADOPTED 2026-09-27 (Hannah): Hs 2.0 m, Tp 7.5 s, asymmetry 0.6,
high-angle 0.5, the same in both windows. Chosen on the metres offset
from the fixed-ends wave grid, on the raw score; the edge ends in
hatteras_site_config.HATTERAS_BE_EDGE_ONLY were solved at exactly these
four values and are not valid at any other. Record:
output/raw_runs/experiments/wave-climate/2026-09-27-wave-recommendation/.
Option B (Hs 2.5 in 2010-2024 only) is recorded, not wired:
hatteras_site_config.HATTERAS_WAVE_OPTION_B.
Until 2026-09-27 these read 2.5 / 8.0 / 0.7 / 0.1, the /10-offset
calibration; every run named before then is named against those.
```

```text
WHERE A RELOCATED ROAD GOES, in metres behind the dune line. This is the
RELOCATION TARGET ONLY -- the road's position at t = 0 always comes from
the period's measured RoadSetback_<year>_dunestart.csv and is untouched
by this.

CASCADE has no separate parameter for the two: cascade_groin.py:689
re-assigns `road_relocation_setback = road_setback` every year, so the
target is whatever the road's MEASURED 1984/2004 offset happened to be.
That is observed geometry, not a design standard, and it ranges 0-430 m
across the 55 road domains. At GIS 85 and 86 it is 0 m, so a relocation
puts the road back on the dune line with no clearance and the next 10 m
of retreat re-fires it: 7 relocations for 7.3 cells of retreat at GIS 85,
6 for 6.0 at GIS 86. 13 of the 18 events in the 1984-2004 calibBE groin
run are that ratchet.

20.0, decided 2026-09-01: "when the road is rebuilt, it is rebuilt to a
standard clearance". CASCADE's own default is 30 (cascade_groin.py:135)
and the matrix was first run at that, but 30 drowned NC-12 at GIS 11 in
all eight 1984-2004 reloc arms. The cause is NOT the clearance itself --
a prescribed historical relocation is stored as a DISPLACEMENT and
`_apply_relocation` adds it to the model's CURRENT setback, so raising
the emergent target raised where the 1999 event landed too: 0 + 77 = 77 m
became 20 + 77 = 97 m, two cells further back, past the point where 24%
of the bordering row is at or below MHW and `bulldoze` gives the road up.
At 20 m the same event lands at 87 m and all eight drownings go away,
with relocation counts across the twelve arms moving only 26 -> 28. See
output/comparisons/relocation/standard_setback/.

That coupling is a real weakness and 87 m clears the threshold by ONE
CELL, so this value is not robust to different forcing. Anchoring
`_apply_relocation` to an absolute setback would remove the coupling and
let this be chosen on its merits alone; it has not been done.

Set to `measured` for the pre-2026-08-31 behaviour, where every domain
relocated to its own measured offset. Setbacks quantise to whole 10 m
cells (`road_start = int(setback / 10)`), so use multiples of 10 -- 20
and 29 are the same model.
```

```text
Management, not physics: a defence someone decides to build. Top-level
in the yaml for that reason, and not in the `scenario` table because no
historical sandbag campaign is reconstructed for either period.
```

```text
Effectively fixed, and deliberately absent from hat_run.yaml -- see the
"not settable here" block at the foot of that file. It stays a field so
HAT_USE_SANDBOX_CASCADE can still force the installed package for a
one-off A/B, and so describe() records which model a run actually built.
It must NOT be derived from groin_enabled: section 12.3's paired
baseline has to be the same model as the groin run in every respect.
```

```text
Aliases kept so an environment variable that predates the yaml still works.
HAT_SOURCE_SINK_PRESET is the name HAT_run_all.py sets; the attribute is the
same, so this is only about the env spelling being longer than the yaml key.
```

```text
A yaml `key:` with nothing after it parses to None. For a
nullable field that is a deliberate "unset"; for any other it is
an unfinished edit, and taking the default silently would run
something the file does not say.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_as_bool()`**

```text
Casts a boolean, strictly.

Accepts real booleans (which is what the yaml parser produces) and the
unambiguous string spellings a shell driver emits. Anything else raises
rather than being silently truthy, because `bool("False")` is True and
that failure mode is invisible in a run log.
```

**`_as_opt_bool()`**

```text
Casts a boolean that may also be explicitly unset.

`null` in the yaml and "" or "none" in the environment mean "leave the
decision to whoever reads this" -- for `relocations` that is the named
scenario, for `show_figures` it is the notebook/.py split.
```

**`_as_opt_float()`**

```text
Casts a float that may also be explicitly unset.

`null` in the yaml and "", "none" or "measured" in the environment mean
"no standard -- use the per-domain measured setbacks". Spelling it
"measured" is allowed because that is what the alternative IS, and a run
log reading `relocation_setback_m: measured` says what happened where a
bare `none` would only say what did not.
```

**`field_default()`**

```text
The code default for one field, ignoring the yaml and the environment.

This is the CALIBRATION value of a setting, not the value the current run
is using. `RUN_CONFIG.hs` answers "what is this run doing"; this answers
"what is it being varied away from", which is the question a sensitivity
sweep has to ask before it can tell a cell from the baseline.

Kept here rather than being re-typed in the sweep script because _FIELDS is
already the one home for these numbers. A second copy in a sweep would go
stale silently and every cell would then be measured against a value the
model no longer uses.

Args:
    name: A RunConfig attribute name, e.g. "hs".

Returns:
    The default for that field.

Raises:
    KeyError: If no such field exists. A typo must not return None and be
        mistaken for a field whose default is genuinely unset.
```

**`_flatten()`**

```text
Flattens a nested yaml mapping to {path tuple: value}.

A block that is itself a known field's parent (e.g. `groin`) flattens into
its leaves; a block that is not is reported as an unknown key by the
caller, path and all, rather than being silently skipped.
```

**`_load_settings_file()`**

```text
Reads hat_run.yaml, or returns nothing if it is absent or suppressed.

Returns:
    (flat mapping of yaml path -> value, the path actually read or None).

Raises:
    RuntimeError: If the file exists but pyyaml is not installed -- the
        settings would otherwise be silently ignored and the run would use
        defaults under the name of whatever the file asked for.
    ValueError: If the file holds a key this module does not know. A
        misspelled key is the failure this catches: ignoring it produces a
        run that used the default while its settings file says otherwise.
```

**`RunConfig()`**

```text
The values that select which run the hindcast performs.

Attributes:
    start_year: A key of HATTERAS_PERIODS -- 1984, 1996, 2004 or 2010.
        Selects the period, and every forcing that follows from it.
    source_sink_preset: "zeroBE", "edgeBE" or "calibBE".
    scenario: A key of the SCENARIOS table in section 3.
    relocations: Overrides the scenario's historical-relocation switch.
        None leaves the scenario preset in charge.
    offset_mode: Which shoreline_offset variant to build: "metres"
        (the default since 2026-09-24), "asrun" (offset / 10, the old
        units error) or "detrended".
    groin_enabled: Whether the groin callback is attached.
    groin_strength_1996, groin_strength_2009: The share of the arriving
        drift the groin holds back, for a run starting in that year.
    groin_strength: Overrides both when set (None by default); needed
        for any other start year with the groin on.
    groin_end_year: The first year the groin no longer works.
    hs: Significant wave height, m.
    wave_period_s: Peak wave period, s.
    wave_asymmetry: Fraction of waves from the left of shore-normal.
    wave_angle_high_fraction: Fraction approaching at more than 45
        degrees. With hs, the four make up the wave climate; a value
        off its default earns the run a `wave...` name token.
    relocation_setback_m: Where a relocated road is rebuilt, m behind
        the dune line. None means each domain uses its own measured
        offset. Off its default it earns an `rset` name token.
    sandbags: Whether sandbag placement is enabled.
    show_figures: True renders figures inline. None means the reader
        decides -- the .py uses False, the notebook True.
    make_gifs: Whether section 9's shoreline animations are built.
    save_model_state: Whether the ~160 MB model .npz is written.
    overwrite: True empties an existing run directory and reuses it.
    use_sandbox_cascade: True imports cascade.cascade_groin, which carries
        the hook the groin callback needs.
    origins: attribute -> "default" | "file" | "env", for describe().
    settings_path: The yaml actually read, or None.
```

**`load_run_config()`**

```text
Reads the settings afresh.

Call this rather than using the module-level RUN_CONFIG when the file may
have changed since import -- which in a notebook is every time, since the
kernel caches the module and an edit to the yaml would otherwise not
apply until a restart.
```

**`describe()`**

```text
Renders the settings and their provenance as a printable block.

Section 3 of both the notebook and the .py prints this, so every run log
records not just what was run but whether each value was typed in the
settings file, driven from the environment, or left at this module's
default -- the distinction that matters when a matrix run and a
hand-iterated run land in the same index.

Returns:
    A multi-line string, no trailing newline.
```

**`_runtime_estimate()`**

```text
Median wall-clock of prior runs of this period, from run_index.csv.

Measured rather than assumed: the index records `runtime_min` for every
run that has completed, so the estimate is this machine's own history for
this period and no constant has to be maintained here.

Returns:
    (median minutes, number of runs it was taken over), or (None, 0) when
    the index is absent or holds no run of this period.
```

**`preflight()`**

```text
Renders what this run will produce, before it produces it.

Answers the three questions worth asking before a long run starts: what
will it be called, where will it land, and does something already live
there. The name is the section 3 preview, not a name built here -- section
7.5 derives the authoritative one from what the modules actually built and
raises if the two disagree, so this can never quietly become the thing
that names the directory.

The collision line WARNS rather than raises. `guard_run_dir` in section 11
is the authority on that, and a second gate here would be a second place
for the rule to live.

Args:
    run_name: RUN_NAME_PREVIEW from section 3.
    run_dir: Directory the run will write to.
    config: Settings to report. Defaults to the module-level RUN_CONFIG.
    index_path: run_index.csv, for the runtime estimate. Defaults to the
        repo's own.

Returns:
    A multi-line string, no trailing newline.
```

</details>

### HAT_run_all.py

End-to-end unattended driver: the comparison matrix, without and with the groin.

Since 2026-10-14 three stages: 1 archive, 2 the no-groin matrix, 3 the groin matrix, at the
corrected groin's strength and end year from the config (a groin run's manifest key carries
both). The seed runs, the M/f sweeps and the joint fit were stages 3-5 and went with the
dipole and blocking groins; the header below is kept as it was written.

From the script's original header:

```text
End-to-end unattended driver: comparison matrix, groin sweeps, joint fit.

Runs the whole thing in one command, in the order the dependencies require,
and can be re-invoked at any point to pick up where it stopped.

THE ORDER, AND WHY IT IS THIS ORDER
    1  archive        The retired run_index.csv is moved aside. Its rows
                      describe runs whose directories no longer exist, and two
                      of them disagree with each other about the background
                      erosion they were fit under, so merging new runs into it
                      would produce an index of mixed provenance.

    2  matrix nogroin The scenarios that are DISTINCT in each period x
                      preset x 2 periods x the relocations arm, groin off.
                      These are the comparison baselines every groin result is
                      read against. They need no M or f, so they can run
                      before anything is fitted -- and running them first also
                      satisfies the rule that a scenario's no-groin run must
                      exist before its groin run, which is how section 12.3
                      resolves its paired baseline.

                      Two guards decide which cells exist, and both skip
                      rather than force a name onto what would be a duplicate
                      run:

                      `scenario_applies` drops `full_no_fill` in 1984-2004.
                      All three nourishment projects are 2014 or later, so
                      there it differs from `full_management` in a switch with
                      nothing to act on: the two build the same modules,
                      derive the same RUN_NAME, and the collision guard
                      refuses the second. The fills axis exists only in
                      period 2.

                      `relocation_applies` drops the reloc arm wherever it
                      would be inert. Relocations need roadway management to
                      have a setback to move, so `natural` and
                      `beachdune_only` carry no arm; and both historical
                      events are 1989 and 1999, so 2004-2024 carries none
                      either -- the 2022 bridge event fires whatever the
                      switch says, which is exactly why turning the switch on
                      there would produce an identical run under a name
                      claiming a contrast.

                      Both presets, both periods, groin off: 22 runs.
                      Per preset, 1984-2004 contributes 4 scenarios + 2 reloc
                      arms and 2004-2024 contributes 5 scenarios + 0 arms, so
                      11 each. `--presets zeroBE` covers half of it.

    3  seed           2 runs: edgeBE / full_management / groin on, one per
                      period, at provisional M and f. These exist only so the
                      sweep has something to validate its duplicated model
                      code against -- the drift guard differences a rate curve
                      against a published matrix run, and there is none on
                      disk. Stage 6 re-runs both at the fitted values.

    4  sweeps         344 cells: 4 sweeps (2 periods x 2 presets). Period 1
                      edgeBE carries the be1 axis and is 215 of them.

    5  joint fit      Intersects the two periods' surfaces per preset. M and f
                      are not separable within one period, so this is where
                      the parameters are actually chosen.

    6  matrix groin   18 runs: the same 18 scenario cells with the groin on,
                      at the fitted (M, f) for their preset.

RESUME
    Every finished job is appended to `driver_manifest.jsonl` the moment it
    lands, keyed by (stage, period, preset, scenario, groin, reloc) -- plus
    (M, fraction) when the groin is on, and plus the SOURCE/SINK DIGEST when
    the preset imposes any, so a revised fit invalidates the runs made at the
    old values instead of silently reusing them. Re-invoking skips anything
    already recorded complete and retries anything recorded failed. The
    manifest -- not the run directory -- is the source of truth, so this file
    never has to reimplement the runner's RUN_NAME construction, which is
    derived in the runner's section 7.5 from what its modules actually built
    and would drift if copied.

    THE KEY MUST NAME EVERY INPUT THAT CAN CHANGE UNDER A FIXED NAME. That
    rule has now been learned three times: the reloc axis, then (M, fraction)
    for groin runs on 2026-08-24, then the source/sink digest on 2026-08-28.
    Each time the symptom was identical -- a run made at superseded values
    matched an unchanged key and was reported "SKIP (done)" -- and each time
    the skip happens BEFORE the collision guard, so --overwrite cannot reach
    it. If a fourth quantity is ever edited between runs while its preset
    keeps its name, put it in the key rather than remembering to clear the
    manifest by hand.

    A manifest written before the M/f fields existed carries six-field groin
    keys that can never match the eight-field ones. Migrate it rather than
    letting every groin run repeat: the rows already RECORD M and fraction, so
    the new key can be rebuilt from each row exactly.

    The digest is NOT retrofittable the same way: rows written before
    2026-08-28 do not carry `be_digest`, and it cannot be recovered from the
    row because the whole point is that the values have since changed. Those
    rows simply stop matching, which re-runs them -- the safe direction. Only
    zeroBE is unaffected, because an "empty" digest is never appended and
    zeroBE imposes {} in every period.

CONCURRENCY
    Matrix runs are SERIAL here for simplicity; since 2026-09-16 nothing
    forces it (run_index.csv is rebuilt from disk, not appended, and each
    run has its own parameter yaml), so two drivers may overlap. The sweep is parallel -- its
    cells write only to their own directories and are collected through a
    single thread -- and it is where nearly all the wall clock goes.

Usage:
    python HAT_run_all.py [--workers N] [--dry-run] [--no-model-state]
                          [--stages 2,4,5] [--presets zeroBE]
```

Notes that were in the code:

```text
The matrix's preset axis is WIDER than the default pair. PRESETS (edgeBE,
zeroBE) is what the matrix covers unless told otherwise, but a hindcast run is
valid under any preset the site config defines, so the rest are reachable
through --presets. Order follows PRESETS first so the default schedule is
unchanged.
```

```text
A 20-year run is a few minutes; an hour is a wide margin that still stops a
wedged run from holding an overnight job open until morning.
```

```text
Hs JOINED THE KEY ON 2026-09-01, the fourth instance of the failure this
docstring already records three of. Without it an Hs = 3.0 matrix matches
every key the 2.5 m matrix wrote and is reported "SKIP (done)" -- and the
skip happens before the collision guard, so --overwrite cannot reach it
either. Appended only when off the default, so every key already on disk
stays byte-identical.
```

```text
The line the runner prints once its section 7.5 has derived the authoritative
name from what sections 5 and 6 actually built.
```

```text
Hs is one of the six settings this driver deliberately leaves to the code
default, so --hs is the only way to move it. Passed explicitly here, which
means a stray HAT_HS in the calling shell still cannot reach a run.
```

```text
reloc-off first, for the same reason the no-groin matrix
runs before the groin one: a reloc arm is read against the
arm that is identical to it in every other token, and the
baseline should be on disk before the arm that cites it.
```

```text
Recorded as ok so a re-invocation does not retry it, with
skipped=True so the summary can tell "did not need running"
from "ran and succeeded".
```

```text
Off by default: the guard is what stops one matrix run silently
replacing another. --overwrite is for the case the guard cannot
distinguish -- a directory that exists but is STALE, because the
forcing or the reported quantity changed underneath it.
```

```text
A `reloc` token only where the arm carries one, matching how the
runner derives RUN_NAME in 7.5. Spelling it "noreloc" everywhere
else would rename every existing log for a switch that is off.
```

```text
NOT recorded on a dry run. The manifest means "this job is done",
and a dry run did not do it -- a row written here would be read by
`already_done` on the next real invocation and skip the very cell
the preview was checking. Previewing a plan must never cancel it.
```

```text
Read back from the runner's output, not derived --
see run_name_from_log. This is what lets a manifest
row be joined to its run_index.csv row.
```

```text
Order follows PRESETS, not the order they were typed, so two
invocations naming the same presets schedule them the same way.
```

```text
Same treatment as --presets: validated against the runner's own table
and reordered into its order, so two invocations naming the same
scenarios schedule them identically.
```

```text
Held for the whole invocation rather than per stage: the parameter file
is shared by every stage, not just the matrix.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_scenario_table()`**

```text
Reads the runner's SCENARIOS table out of its source, without running it.

`HAT_hindcast_1984_2024.py` is a script, not a module: importing it runs a
whole hindcast, so the table cannot simply be imported. It is parsed
instead, because the literal dict in the runner's section 3 is the one
thing that decides which management switches a scenario name stands for,
and a copy typed out here would be free to drift from it.

Drift in the `roadway` flag would be the damaging kind. `relocation_applies`
below reads it to decide whether a reloc arm is a distinct run; if it said
True where the runner says False, the runner would force the switch off,
derive its NON-reloc name in 7.5, and the arm would collide with the very
baseline it was meant to be read against.

Returns:
    {scenario_name: {switch_name: bool}}, in the runner's own order.

Raises:
    RuntimeError: If the table is absent or is no longer written as
        literal dict(...) calls. Loud on purpose: the alternative is a
        silently empty or partial scenario list, which reads downstream
        as "there was nothing to run".
```

**`period_has_fill()`**

```text
Whether any nourishment project falls inside a period.

Resolved with the same `build_schedule` call the runner makes, rather than
by re-implementing the date filter, so this cannot disagree with what
section 6 actually builds.
```

**`period_relocation_years()`**

```text
Relocation-event years that fall inside a period.

Derived from HATTERAS_ROAD_EVENTS rather than listed here, so adding an
event adds its runs. The window matches the loop's in
`run_cascade_simulation`: `current_year` runs from `start_year` to
`start_year + run_years - 1`, so an event dated exactly on a period's end
year belongs to the NEXT period and never fires in this one.

BridgeEvents are deliberately not counted. `apply_historical_event`
applies them whatever `relocations_enabled` says, so they are not part of
what the switch turns on.
```

**`relocation_applies()`**

```text
Whether relocations-on is a DISTINCT run for this cell.

Two ways it is not, and either would put two byte-identical runs on disk
under names claiming to be a contrast:

ROADWAY MANAGEMENT OFF
    The runner's `_RELOCATIONS_FORCED_OFF` turns the switch off wherever
    nothing manages the road -- a relocation cannot move a setback that no
    RoadwayManager reads. The forced-off run then derives its NON-reloc
    name in 7.5 and the collision guard refuses it as a duplicate.

NO RELOCATION EVENT IN THE PERIOD
    Both historical relocations are 1989 and 1999, so 2004-2024 has none.
    The 2022 Jug Handle bridge event is not gated by the switch, so a
    period-2 reloc run simulates exactly what its non-reloc twin does.

Same shape and same motive as `scenario_applies`: skipping is the honest
handling, because running it under a forced name files two identical runs
as if they were a contrast.

Returns:
    (applies, reason). reason is None when it applies.
```

**`scenario_applies()`**

```text
Whether a scenario is a DISTINCT run in this period.

`full_no_fill` exists to isolate the nourishment fill against
`full_management`. In a period with no fill scheduled, the two differ in a
switch that has nothing to act on: they build the same modules, produce
the same result, and -- correctly -- derive the same RUN_NAME, so the
runner's collision guard refuses the second one.

1984-2004 has no fill (all three projects are 2014 or later), so the
fills axis simply does not exist there. Skipping is the honest handling:
running it under a forced name would file two identical runs as if they
were a contrast.

Returns:
    (applies, reason). reason is None when it applies.
```

**`DriverLock()`**

```text
Refuses to start while another driver is running.

Every run rewrites PARAMETER_FILE and reads it back, so two drivers
interleave writes to one file and tear each other's reads. The lock is an
O_EXCL create, which is atomic on Windows and POSIX alike; the PID inside
is for the human reading the error, not for liveness checking.

A crashed driver leaves the file behind. That is deliberate -- a stale
lock reports what to delete rather than being cleared automatically, since
"the lock looks old" is exactly the reasoning that reintroduces the race.
```

**`be_digest_for()`**

```text
Fingerprint of the source/sink VALUES a (period, preset) pair imposes.

Not the preset name -- the numbers behind it. `values_digest` is imported
from the run registry rather than reimplemented so the manifest key and
the `be_values_digest` column in run_index.csv can never disagree about
what a run carried.

Returns "empty" for a preset that imposes nothing, and for the stages
(archive, sweep) whose keys carry no period/preset pair.
```

**`job_key()`**

```text
Stable identity for one unit of work.

The relocations flag joined the key when the reloc axis was added. A key
written before that has five fields and can never match a six-field one,
so an older manifest does not mark the new cells done -- which is the
correct outcome, since those runs were made without the axis and the
reloc arms among them do not exist.

M AND f JOINED THE KEY FOR GROIN RUNS ON 2026-08-24, for the same reason
and after the same failure. The key identified a groin cell by scenario
alone, so when the fitted values changed the old run still matched and was
SKIPPED -- and the skip happens before the collision guard, so --overwrite
could not reach it either. A run made at M = 110, f = 0.4 survived into a
matrix that was supposed to be M = 50, f = 0.6, reported as "SKIP (done)".
Including the values means a revised fit invalidates exactly the runs it
should and leaves the rest alone. Since 2026-10-14 a groin run carries the
corrected groin's strength and end year in their place.

Only groin runs carry them: a no-groin cell has no M or f, and appending
"None|None" to its key would invalidate every groin-off run already
recorded for no reason.

THE SOURCE/SINK DIGEST JOINED THE KEY ON 2026-08-28, third instance of
the same failure. The key named the PRESET but not its VALUES, and the
source/sink table is edited between runs -- that is the whole shape of an
edge solve, where `edgeBE` keeps its name while GIS 1 and 90 move at every
Newton step. So a run made at the previous edge values still matched, and
was reported "SKIP (done)" into a matrix that was supposed to carry the
new ones. Caught before it corrupted anything only because the edgeBE leg
was killed and restarted by hand.

Derived here from (period, preset) rather than passed in, so a call site
cannot forget it. `values_digest` returns "empty" for a preset that
imposes nothing, and an "empty" digest is NOT appended -- that keeps every
zeroBE key byte-identical to what is already on disk, which is correct:
zeroBE imposes {} in every period and no edit to the calibrated table can
change what it ran.
```

**`run_name_from_log()`**

```text
The run name the runner REPORTED, read back out of its own output.

NOT DERIVED HERE. Everything else in this file is deliberately built from
the job key rather than from a run name, because reimplementing the
runner's section 7.5 derivation is how the two drift apart -- the module
docstring says so and that rule stands. This reads back what the runner
itself printed, which is the opposite of reimplementing it: if 7.5 changes,
this follows automatically, and if the line is ever absent the manifest
simply carries no name rather than a wrong one.

It exists so a manifest row can be joined to run_index.csv. The manifest is
keyed on (stage, period, preset, scenario, groin, reloc, M, f, digest) --
the scenario VOCABULARY -- while the index is keyed on the derived name, so
the two could previously only be matched by hand.

Args:
    log_path: Path to the run's captured stdout, or None on a dry run.

Returns:
    The run name, or None if the log is missing or does not carry the line
    (a failed run that died before 7.5, for instance).
```

**`run_hindcast()`**

```text
Runs the hindcast once, driven entirely through the environment.

The runner reads these through HAT_hindcast_config, so no tracked source
is edited and the notebook stays byte-identical to what a driven run uses.

HAT_RELOCATIONS is set on EVERY run, including the reloc-off ones, rather
than only when it departs from the scenario preset. Every named scenario
sets relocations=False, so an explicit "false" is not a departure and the
runner reports none -- but `describe()` then records the value as having
come from the environment, so a run log states that the driver chose it
instead of leaving the reader to infer a default.

THE SETTINGS FILE IS SHUT OUT, TWICE OVER
    HAT_IGNORE_SETTINGS=1 stops the runner reading `hat_run.yaml`, so a
    batch run is described entirely by the values set here plus
    HAT_hindcast_config's defaults. The driver sets only nine of the
    fifteen settings; without this, the other six -- offset_mode,
    make_gifs, Hs, sandbags, use_sandbox_cascade, show_figures -- would be
    taken from whatever experiment was last left in the file, and every
    cell of the matrix would silently inherit it.

    Stray HAT_* variables in the calling shell are dropped for the same
    reason, and reported when they are, since `os.environ.copy()` would
    otherwise carry them into every child.

Returns:
    (ok, seconds, detail).
```

**`matrix_stage()`**

```text
Runs the matrix cells for one groin state, serially.

Three axes: period x preset x scenario, plus a relocations arm on every
cell where relocations-on is a distinct run. `scenario_applies` and
`relocation_applies` between them decide which cells exist; a cell they
refuse is recorded as skipped WITH its reason rather than dropped
silently, so the manifest distinguishes "did not need running" from "was
never scheduled".

Args:
    stage: "matrix_nogroin" or "matrix_groin".
    groin: Whether the groin is attached.
    manifest: Loaded manifest, for resume.
    args: Parsed CLI arguments. `args.presets` and `args.scenarios`
        select which slice of the matrix this invocation covers.

Returns:
    (completed, failed) counts.
```

</details>

### tools/HAT_compare_rerun.py

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

### tools/HAT_index_runs.py

Rebuild output/raw_runs/run_index.csv from the runs themselves, with one ledger of retired runs.

From the script's original header:

```text
Rebuild output/raw_runs/run_index.csv from the runs themselves, and keep a
single ledger of runs that have been retired.

WHY THE INDEX IS DERIVED (Hannah, 2026-09-11; runs stopped appending 09-16)
    Every run writes its own metadata JSON, and that is the source of truth:
    it is produced by the run, beside the run, from the values the run used.
    The index restates those facts in one table so that a question across
    runs -- which topography, which preset, what skill -- is one read instead
    of two hundred.

    A restatement can drift from what it restates. Until 2026-09-16 every run
    also APPENDED its row to the file, which is how a smoke test on 09-10
    truncated floats in five columns of an unrelated row, and why two runs
    could never be in flight at once. Since 09-16 a run writes its row INTO
    its metadata (the "index row" section) and calls the same rebuild this
    tool runs, so the file is regenerated from disk every time and never
    edited in place. `--check` turns "is it right" into a question with an
    answer.

WHAT THE REBUILD DOES (run_registry.rebuild_run_index)
    One row per *_run_metadata.json under raw_runs. A run made since 09-16
    supplies its own row; an older run keeps the row the existing file holds
    for it, found by the old (run_name, Hs_m, arm) key. Every row gets `kind`
    and `tag` from where the run sits (the purpose layout, or the two older
    layouts translated) and a `status`: current, superseded (a matrix or
    sensitivity run on a topography that is no longer its product's CURRENT),
    or archived. The old `arm` column is dropped.

WHAT IS NOT DERIVED
    A run deleted from disk leaves no metadata to rebuild from, so a rebuild
    would drop its row silently. `retired_runs.csv` is the append-only record
    of those: a row that was in the index, whose run is gone. It is written
    here and never rewritten, so the history of what was removed survives a
    rebuild. `--adopt-archives` seeds it from the pre-2026-09-11 archived
    copies of the index.

USAGE
    python HAT_index_runs.py --check          # compare, change nothing
    python HAT_index_runs.py --dry-run
    python HAT_index_runs.py                  # rebuild, retiring vanished rows
    python HAT_index_runs.py --adopt-archives # seed the ledger from the copies
```

<details><summary>Function notes (the original docstrings)</summary>

**`append_ledger()`**

```text
Append retired rows. Never rewrites a row: the ledger is the one
record of what was removed, and a rebuild must not be able to erase it.
The header gained kind/tag on 2026-09-16; an older ledger is widened
once, keeping every row.
```

**`rebuild()`**

```text
Compare the index with disk and, unless asked not to, rewrite it.

Returns the exit code: 1 under --check when they differ, else 0.
```

</details>

### tools/HAT_list_runs.py

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

### tools/HAT_migrate_run_layout.py

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
    run_index.csv               rebuilt (HAT_index_runs), every row with
                                kind, tag and status.

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

### tools/HAT_period_input_check.py

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

### tools/HAT_rename_experiments_20260925.py

Rename the older experiments and group every experiment by theme (the 2026-09-25 one-off).

From the script's original header:

```text
Rename the older experiments and group every experiment by theme (2026-09-25).

Hannah, 2026-09-25: "rename the old experiments so their names are clearer
about what they tested", then "organize the experiment folders ... by what
the theme or investigation is". Names and themes approved the same day:

    experiments/<theme>/<date>-<what it tested>/

What this does, while nothing is running:
  1. removes the parameters-only folders a paused sweep leaves behind (a run
     killed at start-up), which would block its resume
  2. moves each study to <theme>/<new name>
  3. rewrites the old path/name inside every text file of every study
     (run-metadata tags, NOTE/README/FINDINGS, logs and CSVs naming runs),
     fixing the relative links between studies (one level deeper now) and
     notes the old name at the top of the renamed studies' NOTE.md
  4. rewrites the old names in every tracked text file outside the archive
     (scripts -- a study's folder is also its run tag -- comparison READMEs
     and captions) and in the memory notes
  5. writes experiments/README.md (the map) and a README per theme
  6. rebuilds the run index and checks every run is under its new name
The archive and retired_runs.csv are history and are left as they were.
cascade_pipeline.run_registry.check_tag allows 4 tag levels since this change.

    python scripts/hatteras_ms/tools/HAT_rename_experiments_20260925.py [--dry-run]
```

Notes that were in the code:

```text
2010_2024 only: the pause stops the sweep as its 2010-2024 cells start.
A DROWNED run also leaves only its parameters file, and those (1996-2010,
Hs 0.75 / Tp 10) are results that must stay.
```

```text
Resumable (2026-09-25: the first pass stopped at a folder held open by
a shell's working directory): a study already at its new place is
skipped; the text rewrite cannot double a prefix, and the rename note
is written once.
```

<details><summary>Function notes (the original docstrings)</summary>

**`replacements()`**

```text
Ordered (old, new) string pairs. Relative links first: from inside a
study, ../<other> becomes ../../<theme>/<other> and the chain INDEX moves
one level up; then every remaining bare name gets its theme.
```

</details>

### tools/HAT_rerun_arm.py

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

### tools/HAT_run_supersession_report.py

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

### HAT_hindcast_1984_2024.py

Hatteras Island CASCADE hindcast runner: the headless twin of HAT_hindcast_1984_2024.ipynb.

From the script's original header:

```text
HATTERAS ISLAND -- CASCADE hindcast runner (1984-2004 / 2004-2024).

Headless twin of `HAT_hindcast_1984_2024.ipynb`. The notebook is the source of
truth: this file mirrors it section for section, in the same order, calling the
same `cascade_pipeline` / `hatteras_site_config` code with the same values, so
the two produce the same run.

WHAT THAT MEANS IN PRACTICE

  - Nothing is duplicated that the packages already own. The background-erosion
    presets, road events, nourishment projects, community zones and annotations
    all come from `hatteras_site_config`; the loaders, audits, verifiers and
    plotting come from `cascade_pipeline`. An earlier version of this file kept
    its own copies of all of it, and they drifted.
  - The functions the notebook defines in its own cells (`build_cascade`,
    `run_cascade_simulation`, `build_target_table`, `brie_r_ipl`,
    `measure_groin_extent`, ...) are copied here verbatim. Edit one, edit both.
  - Executed top to bottom at module level, which is what a notebook does.
    Nothing imports this file.

THE TWO DELIBERATE DIFFERENCES FROM THE NOTEBOOK

  - The notebook's QC plots are not here: island orientation, initialization
    plan view, RSLR, storms, source/sink, roadway plan view, beach/dune
    footprints, groin setting, CoastSat target. They display, they do not feed
    the run. Every *check* that sits beside them is kept -- the Barrier3D units
    contract, the setback audit, the nourishment audit, the annotation guard.
  - `SHOW_FIGURES` (below) defaults to False. The notebook renders its final
    figures inline with `show=True`; a headless run cannot. The files written
    to disk are identical either way.

WHERE THE SWITCHES LIVE

  `hat_run.yaml`, beside this file. Edit it, save it, run. That file is the
  interface; the sections below only read it, through HAT_hindcast_config:

    section 1   output.show_figures -- read at import, because it selects a
                matplotlib backend
    section 3   start_year, source_sink, scenario, relocations, offset_mode
                (scenario expands to ENABLE_ROADWAY_MANAGEMENT,
                ENABLE_BEACH_DUNE_MANAGEMENT, ENABLE_NOURISHMENT_FILLS and
                ENABLE_HISTORICAL_ROAD_RELOCATIONS; commented override lines
                sit directly below the table)
    section 7   groin.enabled, groin.trapping_M, groin.deterioration_f
    section 9   output.make_gifs
    section 11  physics.wave_height_Hs, sandbags, output.overwrite,
                output.save_model_state

  An environment variable beats the file, so `HAT_run_all.py` can drive a
  matrix without editing anything; it also sets HAT_IGNORE_SETTINGS=1, so a
  half-finished experiment left in the yaml cannot reach a batch run.

  STILL TYPED IN THIS FILE, because they are properties of the study rather
  than of a run: dune thresholds and datums (section 11), FLIP_SIGN_MODEL and
  PLOT_REAL_DOMAINS_ONLY (section 9), the groin's geometry and deterioration
  schedule (section 7), NUM_CORES (section 11).

RUN_NAME is derived from those switches in section 7.5 -- it is not typed by
hand, so the output directory cannot disagree with what was simulated.
```

Notes that were in the code:

```text
cascade_pipeline and hatteras_site_config live in scripts/, which isn't
installed. The notebook walks up from cwd to find pyproject.toml; a script
knows where it is.
```

```text
HAT_hindcast_config sits beside this file. Running the .py puts that
directory on sys.path automatically; importing it from the notebook does
not, so it is added explicitly and both files reach the module the same way.
```

```text
--- the settings file, read before anything it selects ----------------------
`hat_run.yaml` is where a run is chosen; see HAT_hindcast_config for the
precedence rules. It is read HERE rather than in section 3 because two of
its values are settled at import time and cannot be changed afterwards:

use_sandbox_cascade  decides which Cascade class cascade_pipeline.hindcast
binds, which happens the moment it is imported below.
It is pinned True and has no line in the yaml -- see
the "not settable here" block at the foot of that
file for why it must not follow groin.enabled.
show_figures         decides the matplotlib backend, which is fixed by the
first figure created.

Everything else is re-read in section 3, so an interactive edit to the yaml
applies on a re-run of that cell without restarting the kernel. Section 3
also checks these two against the file as it stands then, and raises if they
have changed -- a stale sandbox flag is the failure where the groin silently
does nothing.
```

```text
The notebook draws its final figures inline. A headless run cannot, so the
figures are saved and not shown. This is the only behavioural difference in
the output path; the files on disk are identical.

`output.show_figures: null` in hat_run.yaml means "whichever file is being
run decides", which here is False. Setting it true in the yaml forces
figures open in a headless run, which is why it is not the default.
```

```text
NO LONGER period-independent (changed 2026-08-25). Until then both periods
read one topography and this section said so. They now start from different
DEMs, so the product is selected from HATTERAS_PERIODS:

1984  ->  1984-start   from DEM 2009-2014-1996
1996  ->  1984-start   from DEM 2009-2014-1996
2004  ->  2004-start   from DEM 2009-2014
2010  ->  2004-start   from DEM 2009-2014

Two products, four periods (2026-09-11): each new period reads the surface
whose survey is nearest its own start, which is a product that already
exists. The mapping is YEAR_PRODUCT's, not this comment's -- ask it.

The VERSION within a product is still resolved, never pinned.
```

```text
There is no TOPO_DUNE_INIT_YEAR any more. The arrays carry no year - the
period is the PRODUCT DIRECTORY - and this runner no longer builds their
names at all; build_domain_file_paths() delegates to the resolver. See the
note at the top of hat_topo_version.py for why a per-period tag was tried
and reverted.
```

```text
WHICH EXTRACTION -- resolved, not pinned. topo_dirs() reads VERSION out of
HAT_dune_topo_extractor.py, so the runner, the dune-start road setbacks and
the audits always describe the same extraction, and a version missing from
disk is a loud error rather than a silently stale read.

Pinning it by hand is what let this runner sit on 2009_v3 while the setbacks
were rebuilt on v4 (2026-08-19) and then v5 (2026-08-20). The setback is
measured from interior row 0, so that combination places the road against a
row that does not exist on the grid being run.

To reproduce an older run deliberately:
topo_dirs("2004-start", override="v3").
```

```text
_BOOT_CONFIG, not RUN_CONFIG: the period must be known HERE, and RUN_CONFIG is
not loaded until section 3. The two are compared a few lines below section 3's
reload, and section 3.0 re-asserts that this product matches the period that
actually ran, so a boot/run divergence cannot silently pick the wrong barrier.
```

```text
The rate fits are DATA and the model reads them; where they live is
hat_observed_rates.py's to say (2026-09-18), not this file's.
```

```text
Taken from what topo_dirs() RETURNED rather than re-joined from parts. The
old line rebuilt the path independently, which is how a resolver gets bypassed
without anyone noticing - the same failure mode as HAT_road_elevation.py.
```

```text
Runs are filed per period: section 3 builds OUTPUT_BASE_DIR =
OUTPUT_ROOT / "<start>_<end>" once START_YEAR is expanded, so a
1984-2004 run and a 2004-2024 run of one scenario cannot land beside
each other. run_index.csv stays at the root, covering both periods.
```

```text
The PRODUCT is passed, not the directories: build_domain_file_paths now
delegates to hat_topo_version.domain_arrays(), which resolves the directory
and the filename together. Nothing here spells an array name.
```

```text
--- 2.3 verify every file exists --------------------------------------------
A stale TOPO_DUNE_VERSION or a moved data folder fails here with a count and
the first offender, rather than as an opaque traceback inside Barrier3D.
```

```text
--- 2.4 units check against Barrier3D's input contract ----------------------
load_input.py converts the SCALAR yaml parameters but loads the elevation and
dune .npy files verbatim, so the arrays on disk must already be in decameters
relative to MHW. A file written in metres would run without error and model an
island 10x too tall. Checked on the RAW arrays, before any DAM_TO_M scaling.
```

```text
START_YEAR is the one flip. It selects a period from HATTERAS_PERIODS and
everything in section 4 follows from it: run length, RSLR rate, storm series,
background-erosion preset, road-setback file, nourishment settings.

The run-selecting values in this section come from `hat_run.yaml`, read
through HAT_hindcast_config: edit that file, save it, run. A driver can also
set them through the environment without editing anything, and
interactively you can still type over any assignment below -- the loaded
value is only the default, and the assignment is still the last word.

load_run_config() RE-READS the yaml on every call rather than reusing the
object section 1 built. In a notebook the module is cached by the import
system, so a fresh read is the only way an edit to the yaml applies without
restarting the kernel.
```

```text
The two values section 1 already spent. They select an import and a
matplotlib backend, so a yaml edited after section 1 ran would be reported
in the block below while the process is still running the old choice --
and for the sandbox flag that means a groin-on run whose groin silently
does nothing. Raised, not warned.
```

```text
The source/sink axis of the run matrix. Each name states a hypothesis about
where the alongshore sediment budget is unresolved:
"zeroBE"   nothing imposed anywhere
"edgeBE"   only the two end domains, absorbing the open-boundary artifact
"calibBE"  the full per-domain fit against the CoastSat LRR
```

```text
The management axis of the run matrix. Every switch that decides whether a
CASCADE management module runs at all lives here, in one block, because a
switch kept beside the section that uses it is how a run ends up half-managed
and named for a scenario it did not simulate. Sections 5 and 6 read these;
they do not define their own.

ENABLE_ROADWAY_MANAGEMENT     roadway_manager: bulldozing, dune rebuild to
the design elevation, setback tracking, road
drowning. Off leaves the road as forcing that
nothing acts on.
ENABLE_BEACH_DUNE_MANAGEMENT  beach_dune_manager, module and all: overwash
filtering, the fixed dune line
(dune_migration_on = False), the 50 m
community-width drowning check, and fills.
Off is the only way to get natural shoreline
behaviour in the village domains.
ENABLE_NOURISHMENT_FILLS      Within an enabled beach_dune_manager, whether
historical fill is actually spent. False
leaves the module and its footprint exactly
as they were, so a fills-on / fills-off pair
differs in the fill and nothing else.

Both modules off is the natural-dynamics run. The groin is deliberately not
in this block: it is a structure, not a management module, and keeps its own
switch in 7.1 so groin-only scenarios stay reachable.
```

```text
The management axis of the run matrix, named rather than typed as four
booleans, so a scenario is chosen in one word and cannot be assembled wrong
by accident. What each switch controls:

roadway      roadway_manager: bulldozing, dune rebuild to the design
elevation, setback tracking, road drowning. Off leaves NC-12
as forcing that nothing acts on -- still loaded, still
audited in 5.1, never handed to a RoadwayManager.
beach_dune   beach_dune_manager, module and all: overwash filtering, the
fixed dune line (dune_migration_on = False), the 50 m
community-width drowning check, and fills. Off is the ONLY
way to get natural shoreline behaviour in the village
domains -- see the section 6 markdown on what is always-on.
fills        Within an enabled beach_dune_manager, whether the historical
fill is actually spent. False leaves the module and its
footprint exactly as they are, so a fills-on / fills-off pair
differs in the fill and nothing else.
relocations  Historical NC-12 relocation events. False in every named
scenario; read 5.1 before overriding it on.
```

```text
Read here rather than beside the offset load in section 4: the run-name
preview below needs it, and a value the name depends on must be settled
before the name is predicted.
```

```text
THE REACH (2026-09-16, the Pea Island extension experiment).
hatteras_site_config built HATTERAS_DOMAINS from HAT_GEOMETRY in the
environment at import, before this config was read; the config carries
the same name so a yaml `geometry:` cannot silently disagree with the
arrays already in memory. "base" is GIS 1-90, every matrix run.
```

```text
The groin is NOT part of the scenario, deliberately. 12.3 measures its effect
against a paired no-groin baseline identical in every other token, so every
scenario is run twice -- False first to create the baseline, then True.
Folding it into the table would double the table to say the same thing.
```

```text
ONE GROIN since 2026-10-14, the corrected groin, so its runs carry the plain
`groin` token. The dipole's runs used `groin` too and the blocking groin's
`groinblock`; both are in the dated groin experiments and on the offload drive,
not in the matrix.
```

```text
Where a RELOCATED roadway is rebuilt, in metres behind the dune line. This is
the relocation target only: the road's position at t = 0 always comes from the
period's measured RoadSetback_<year>_dunestart.csv, and build_cascade applies
this after construction so that measured geometry is untouched.

CASCADE itself has no separate parameter -- cascade_groin.py:689 re-assigns the
relocation target from the initial setback every year -- so without this every
domain relocates to wherever its road happened to sit in 1984/2004. That is
observed geometry, not a design standard, and at GIS 85 and 86 it is 0 m, which
returns a relocated road to the dune line with no clearance and re-fires on the
next 10 m of retreat. `measured` in hat_run.yaml restores that behaviour.
```

```text
--- sensitivity tokens ------------------------------------------------------
The wave climate and the relocation target are FORCING, not management: they
change what is simulated without changing which modules get built, so nothing
in SCENARIO_SWITCHES sees them. A sensitivity cell that moves one of them
would therefore derive the matrix run's exact name and overwrite it. These
two tokens are what keep a cell in its own directory, and they are None at
the calibration values, so no existing run is renamed.

Read from RUN_CONFIG here rather than from section 11's constants because 7.5
needs them and 7.5 runs first. Unlike the management switches there is no
built module for these to disagree with -- they are numbers handed straight
to build_cascade -- so the 7.5 preview check is an identity for them, and
section 11 asserts it is still handing over the values named here.
```

```text
Emitted into the run NAME again since 2026-09-16 (it scoped the directory
between 09-01 and 09-16): a sensitivity cell is its baseline's name plus
this token, filed under sensitivity/<axis>/ by run_registry.
```

```text
HAT_RELOCATIONS overrides the scenario when set; None leaves the preset
in charge. Either way the departure is detected just below and the
`reloc` token still comes from the switch, not from SCENARIO.
```

```text
--- one-off overrides -------------------------------------------------------
Uncomment to depart from the named scenario for a single run. The departure
is detected and printed below, and the run name is still derived from the
switches rather than from SCENARIO, so an overridden run cannot be filed
under the scenario label it departed from.
ENABLE_ROADWAY_MANAGEMENT = False
ENABLE_BEACH_DUNE_MANAGEMENT = False
ENABLE_NOURISHMENT_FILLS = False
ENABLE_HISTORICAL_ROAD_RELOCATIONS = True
```

```text
Normalised to the canonical key. The deprecated aliases ("base",
"calibrated") still run, but it is the canonical name that reaches RUN_NAME
in 7.5 -- an alias can never put a stale token in a directory name.
```

```text
Section 2 picked the topography from _BOOT_CONFIG.start_year, before
RUN_CONFIG existed. If those disagree the run would model the wrong barrier
entirely, so it is checked rather than assumed.
```

```text
--- resolve the combinations that cannot both be true -----------------------
A fill cannot land in a domain with no BeachDuneManager: 6.3 would report it
as dropped and the run would nourish nothing. A relocation event cannot move
a setback nothing reads. Both are resolved here and announced below, rather
than left as a contradiction for a later section to trip over.
```

```text
Run name: the period stem only. The scenario suffix is derived from the active
management switches in 7.5 -- a hand-typed label is how a groin-off run ends
up in a directory named for a groin-on one.
```

```text
Run directories are filed by period. Resolved here, not in section 1,
because the period is not known until START_YEAR is expanded above.
Every later section derives its paths from OUTPUT_BASE_DIR -- RUN_DIR
in 9, the paired groin baseline in 12.3 -- so scoping it here scopes
all of them, and the baseline lookup can no longer resolve to a run
from the other period.
```

```text
Filed by period, then by source/sink preset. The preset directory is
redundant with the preset token in RUN_NAME on purpose: the token is what
run_index.csv, the logs and the figure captions all key on, and the
directory is only there so the three presets of one period can be read
side by side instead of interleaved in one listing of thirty-odd runs.
Filed by PURPOSE (2026-09-16): matrix/, sensitivity/<axis>/,
experiments/<tag>/, versions/<tag>/ -- see run_registry. The KIND says what
the run is for and the TAG names the experiment or version; a matrix run has
no tag, and a sensitivity cell's tag is the axis its name token belongs to.
Before this a forcing value, an input version and an experiment label were
all "arms", and nothing said which was which.
```

```text
HAT_ARM_TAG is the pre-09-16 spelling. Read as an experiment tag so an old
driver still files somewhere sensible, and say so in the log.
```

```text
A run forced off the calibration wave climate is NOT a matrix run, whatever
the environment says: filing it in matrix/ under a token-bearing name would
put an experiment beside the production runs. Refused, not re-filed.
```

```text
Same rule for the island offset SOURCE (2026-09-22). A run built from the
CoastSat shoreline instead of the dune line derives exactly the same name as
the matrix run it is being compared against -- nothing in the name carries
the source -- so filing it in matrix/ would overwrite that run's outputs with
a different model input behind them. Refused, not re-filed.
```

```text
The model state (~99% of a run's size) is kept for MATRIX runs only unless
output.save_model_state was set explicitly: a sweep cell or an experiment
is re-runnable in minutes and its product is a figure, not the pickle.
```

```text
Built by run_registry, not joined here. The reader that resolves a finished
run and the writer that files one have to agree about the layout, and the
only way they cannot drift is to be the same code.
```

```text
WHICH BARRIER3D (2026-09-24). Barrier3D is installed editable, so the branch
checked out in its repository is the model. Since 2026-09-24 that must be
fix/route-overwash-axis-swap, which corrects the route_overwash indexing
bug (experiments/code-checks/2026-09-24-metres-3-barrier3d-overwash-fix/NOTE.md); `git checkout master`
there would silently put runs back on it. Recorded in the metadata and the
index; warned about here, not refused, so a deliberate unfixed run (to
reproduce an old one) is still possible.
```

```text
THE ADOPTED MODEL (2026-09-28, Hannah). The parameter template switches on
per-cell dune ceilings (DuneCeilingFromStart). A Barrier3D without them
ignores that setting and grows every dune toward the 3.4 m NAVD88 default,
the fault that made the model overwash far more than the imagery -- so a
run on one is refused, not warned about. To reproduce a run made before
2026-09-28, check out fix/route-overwash-axis-swap in the Barrier3D
repository and use the parameter template and storm variant (v3_72) of
that date. experiments/storms-overwash-and-dunes/ holds the record.
```

```text
--- the name this scenario will produce, predicted from the switches --------
Advisory only. 7.5 derives the authoritative RUN_NAME_BASE from what sections
5 and 6 actually built and raises if the two disagree, so this preview can
never quietly become the thing that names the directory. Token order matches
SCENARIO_SWITCHES in 7.5 exactly -- that is what makes the comparison valid.
The period's fill is resolved with build_schedule, the same call 6 makes,
rather than by re-implementing the date filter here.
```

```text
What this run will be called, where it will land, whether something already
lives there, and roughly how long it will take -- reported HERE rather than
at the section 11 guard so a name collision is visible before sections 4-10
do their work. Advisory: `guard_run_dir` in 11 is still the authority on
what a collision does, and this does not duplicate that decision.
```

```text
--- 3.1 island offsets ------------------------------------------------------
One value per padded domain giving that domain's cross-shore starting
position. Becomes `shoreline_offset` on the Cascade() call.
```

```text
OFFSET_MODE selects which shoreline_offset variant is built; see
cascade_pipeline.hindcast.build_island_offset. The default since
2026-09-24 is "metres": the file is metres and Cascade wants metres, so it
goes in as it is. "asrun" reproduces the old units error (offset / 10) for
the runs made before that date.
```

```text
island_offset is what Cascade is handed: metres, except in asrun, where
it is the file / 10. Multiplying by DAM_TO_M unconditionally printed the
asrun span right and every metres span ten times too large ("0-62186 m").
```

```text
--- 4.2 storm series --------------------------------------------------------
One row per storm: time (1-based model step), Rhigh, Rlow (decameters),
period, duration (hours).
```

```text
--- 4.3 source/sink (background erosion) ------------------------------------
Per-domain rate in m/yr, passed to Barrier3D as `Rat`. Sign convention from
cascade/brie_coupler.py: (-) = erosion, (+) = accretion. Presets are sparse --
a domain absent from a preset gets 0.0 m/yr.
```

```text
be_rates() rather than a dict lookup: since 2026-09-11 not every wired
period is calibrated, and it says which fit is missing instead of
raising a bare KeyError on the year.

THE FIRST EDGE PROBE OF AN UNSOLVED PERIOD (2026-09-16, the 2010 solve).
An edgeBE run on a period with no edge entry yet has nothing to look up,
but that is exactly the run a Newton step is: HAT_BE_OVERRIDE below
supplies both ends. So an unsolved period starts from an empty mapping
WHEN an override is present, and the guard after the override block
still refuses the run if either end is left without a nonzero rate.
Without an override the refusal stands, naming the missing fit.
```

```text
HAT_BE_OVERRIDE -- per-domain rates for THIS run only, "gis=rate" pairs, e.g.
HAT_BE_OVERRIDE="1=-45.2,90=11.8". Unset, nothing below runs and the preset
is used exactly as it comes out of the config.

WHY THIS EXISTS. Solving the two locked end domains is a Newton iteration:
run, read the residual at GIS 1 and 90, step, run again. Without an override
every step means editing hatteras_site_config.py, and the config is GLOBAL --
so a solve for one forcing arm silently redefines what every other run on the
machine means, and an interrupted solve leaves the production preset holding
a probe value with nothing on disk saying so. That is the exact failure this
file has hit before with a stale topography pin.

It does not weaken provenance. `values_digest(DOMAIN_BE_RATES)` fingerprints
the mapping AFTER this block, so an overridden run gets a different
be_values_digest from the preset it started from, and section 12.4 writes
that digest into run_index.csv. A run forced this way cannot be mistaken for
a config run; it is separated by the same column that separates two edits of
the config from each other.
```

```text
Copied before mutating: HATTERAS_BE_PRESETS hands back the config's own
dict, so writing into it would edit the preset in memory for anything
else importing it in this process.
```

```text
An edgeBE run imposes a value at BOTH ends. The config guarantees that
for the base geometry; in an extended geometry (2026-09-16) the new
end has no solved value until the Newton steps are done, and each step
supplies its probe through HAT_BE_OVERRIDE. Refused here rather than
run as a zeroBE in edgeBE clothing.
```

```text
A derived fact, not a switch: the preset decides it. Checked because the two
used to be separate knobs that could disagree, and the run name carried both.
```

```text
NC-12's forcing is three things: where the road sits (setback, by period), how
high it is (elevation, period-independent), and which domains are managed.

Two conventions, both Barrier3D's:
- road_ele is metres MHW-relative, not NAVD88 -- bulldoze writes it straight
into the interior grid, which the extractor stores MHW-relative.
- relocation events carry a DISPLACEMENT, not an absolute setback. CASCADE
already decrements the setback by dune migration each year, so adding the
measured displacement counts the retreat once; an absolute setback
referenced to an older dune line counts it twice.
```

```text
ENABLE_ROADWAY_MANAGEMENT and ENABLE_HISTORICAL_ROAD_RELOCATIONS are set in
section 3, with the other management switches. The forcing below is loaded
either way: with management off the setbacks and elevations are still read
and audited in 5.1, they simply never reach a RoadwayManager.
```

```text
--- 5.1 pre-flight audit: which road_offset will not survive year one -------------
bulldoze tests the two rows FLANKING the road -- never the road's own cells --
and drowns it when either flank is more than 20% water. A drowned road is not
a warning: roadway_manager sets _drown_break and returns immediately on every
later year, so the domain becomes an unmanaged barrier wearing a road label.
```

```text
Two different things arrive bundled in one CASCADE module.

ALWAYS-ON wherever the module is enabled, every year: a percentage of overwash
deposition is removed from the interior and returned to the shoreface, the
dune line is held fixed (dune_migration_on = False), and the community is
abandoned if the average interior width falls below 50 m.

EVENT-DRIVEN, only where and when nourish_now says so: sand is added to the
shoreface. There is no way to get the second without the first.

Three conventions, all places this pipeline has been wrong before:
- overwash_filter is a PERCENT, not a fraction. filter_overwash divides by
100. A value of 0.4 filters 0.4% of overwash, which is indistinguishable
from no filtering. BeachDuneConfig now refuses the fraction scale.
- The per-year volume goes to cascade.nourishment_volume, which is where
CASCADE reads it from. A volume written onto the BeachDuneManager instance
lands on the attribute CASCADE overwrites one line before the manager
reads it, and the fill quietly spends the Cascade() init default.
- The manager's time series are offset by one; NourishmentSchedule.time_index
is the single place that conversion lives.

The module footprint is the UNION of the permanent community zones and every
domain that receives fill. Where that overlaps the roadway footprint, both
managers run on the same domain -- reported here, verified in section 12.
```

```text
--- schedule: one project list, period-filtered ------------------------------
Every Hatteras project falls in 2004-2024, so a 1984 run builds an empty
schedule from this same list rather than needing a period-keyed one.
```

```text
--- the schedule the model is actually driven by -----------------------------
ENABLE_NOURISHMENT_FILLS (section 3). Suppressing fills builds an EMPTY
schedule for the same period rather than skipping apply_to_cascade: the loop
still rewrites nourish_now and nourishment_volume to zero every year, and the
audit below and the verification in 12.2 both check against what was driven
rather than what was intended. BN_SCHEDULE itself is left intact and still
defines the module footprint below, so turning fills off changes the fill and
not the footprint.
```

```text
--- what CASCADE is handed ---------------------------------------------------
The filter is inert with the module off -- only BeachDuneManager applies it --
so it is built as zeros there rather than left at its community values. The
array handed to CASCADE should state what the run does, not what it would do.
```

```text
Placeholder for the Cascade() call. Every value is rewritten each model year
by BN_SCHEDULE.apply_to_cascade(), so this is 0 rather than
PERIOD["nourishment_volume"]: if the schedule ever fails to reach the model, a
year should nourish nothing rather than quietly nourish the default.
```

```text
cascade/groin.py attaches through cascade._groin_callback and is called once
per model year from inside Cascade.update(), immediately before the alongshore
transport solve. Each active year it adds -M to the updrift domain and +M to
the downdrift domain of x_s_dt. BRIE's implicit diffusion solve spreads that
dipole in the same step, so the fillet's taper and extent are EMERGENT -- only
its amplitude is imposed.

It is a forcing, not a barrier: nothing here blocks alongshore transport. The
pair is volume-neutral by construction, trapping never saturates on state, and
install_year is inert (1969 precedes both periods).

Two independent estimates of M disagree by an order of magnitude -- the
shoreline-position fit says 50 m/yr, the sediment budget says that is
unaffordable. Configured at the sweep's value; the breach is REPORTED, not
corrected. Section 12 reports the resulting misfit.
```

```text
GROIN_ENABLED is set in section 3, beside the scenario, so one cell decides
what a run simulates. The guard below stays here: this is where a wrongly
configured groin would silently do nothing.
```

```text
The pre-AST hook exists ONLY in cascade/cascade_groin.py. Attaching a callback
to a Cascade built from the real package is a silent no-op: the run succeeds,
the groin does nothing, and the output looks like a valid groin-on run.
```

```text
Net alongshore transport on Hatteras is southward, so updrift = north. The two
domains must be adjacent -- they share the blocked boundary at GIS 5.5.
```

```text
--- the groin's two numbers ------------------------------------------------
The share of the drift arriving at the face it holds back, and the first year
it no longer works, both from hat_run.yaml; the callback is built with them as
a strength history, [(1969, strength), (end_year, 0)].
```

```text
--- 7.2 build the callback ---------------------------------------------------
GROIN_CB is built unconditionally so 7.4's report renders either way.
GROIN_CALLBACK is the one attached to the model in section 11.
```

```text
--- 7.5 scenario summary, and the run name derived from it ------------------
RUN_NAME_SUFFIX used to be typed by hand: run with the groin off, forget to
retype the label, and the output lands in a directory named for a different
experiment. It is now derived from the switches themselves.
```

```text
No token when "asrun": every run predating the shoreline_offset
unit finding used it, and adding a token would rename them all.
```

```text
No token of its own: it is implied by the preset, and checked against it
in 4.3. Emitting both produced names like "..._base_noBE_..." that said
the same thing twice without saying which zero it was.
```

```text
"nonourish" only when there was fill to withhold and a module to
withhold it in: with beach_dune_manager off, "nobdm" already says no
fill, and 1984 has no project to suppress. The project count is
deliberately not in the name -- it is a property of the period, not of
the scenario, and run_index.csv carries it as nourishment_projects.
```

```text
Forcing, tokened only when off the calibration value -- see section 3.
LAST, so a sensitivity cell's trailing token is the wave one and
run_registry.sweep_family reads the axis off it.
```

```text
The section 3 preview was predicted from the switches; this name is derived
from what sections 5 and 6 actually built. A difference means a switch did
not reach the module it names -- the failure that produces a run filed under
a scenario it did not simulate. Raised, not warned.
```

```text
The observational target. The "calibBE" source/sink preset was fit against
the curve produced here, and section 12 draws the model against it.

LOWESS runs at TRANSECT resolution using along-coast distance as x, then
averages to domain resolution. The widest window is the target, implicitly:
rate_comparison picks max(window_domains). TARGET_WINDOW names that choice.

skip_southern_domains = 10 suppresses LOWESS across GIS 1-10. That is
DISPLAY-ONLY -- LOWESS still fits over all transects and only the result is
truncated, so the southern transects still pull the values just north of the
cut. Left as built, because the calibrated preset was fit against this curve.

The Buxton groin sits at GIS 5.5, so both its flanking domains fall in the
unsmoothed zone: the groin's observational target is a raw mean over a handful
of transects, not a point on the reference curve.
```

```text
--- 8.1 load both periods ----------------------------------------------------
Both periods always load. This is observational data, not simulation forcing,
and section 12 needs the non-active period for reference styling -- so it is
deliberately NOT folded into HATTERAS_PERIODS.
```

```text
The two periods added 2026-09-11. Every dataset loads on every run --
section 12 styles the non-active ones as reference curves -- so all four
windows are drawn whichever one is being simulated.
```

```text
The 7-domain smoothing ALONE (Hannah, 2026-09-28: "the smoothing rate for
everything in this project should be 7 domains, that is what our group chose
as the range"). Every skill target is built at TARGET_WINDOW = 7, asserted
just below. From 2026-09-10 to 09-28 this was 10 alone; runs made before
09-28 were scored against the LOWESS-10 target, and the end-domain rates in
HATTERAS_BE_EDGE_ONLY were re-solved for LOWESS-7 the same day.
```

```text
Named here rather than left implicit. rate_comparison resolves the reference
window as max(window_domains); this makes that choice visible, and the
assertion below catches a window list whose maximum is not what was intended.
```

```text
AN EXTENDED GEOMETRY (2026-09-16) has CoastSat rows beyond GIS 90 of its
own (coastsat_extension_lrr.py). Its active window loads the surveyed
table WITH those rows appended, so one LOWESS runs over the whole reach
and the new end domain has a target. The surveyed table stays the
interior score's target (COASTSAT_TARGET_BASE, below) in every geometry.
```

```text
The surveyed reach's own target: what the interior score is graded
against in every geometry (2026-09-16). The same object as
COASTSAT_TARGET in the base geometry; in an extended one, the surveyed
table smoothed over GIS 1-90 alone, exactly as the baseline run saw it.
```

```text
Not the plotting functions -- the configuration they need, gathered before
section 12 calls them.

Every plotting entry point defaults to the PACKAGE defaults, and
DEFAULT_ANNOTATIONS is an AnnotationConfig() with every field empty. Omit
`annotations=HATTERAS_ANNOTATIONS` at any one call site and that figure comes
out with no villages, no piers, no groin line and no shoal zones -- nothing
raises. Hence the keyword bundles below, splatted into every call.
```

```text
WHICH ESTIMATOR THE FIGURES DRAW. "lrr" is the OLS slope through every
annual state; "endpoint" is (x[-1] - x[0]) / RUN_YEARS. The observational
target is an LRR -- CoastSat's transect_lrr_full.csv is a per-transect OLS
slope with r_squared and unc_m_yr beside it -- so "lrr" is the only setting
that puts the same quantity on both sides of the comparison. "endpoint"
exists to redraw a pre-2026-08-22 figure, not as an alternative.
```

```text
Below this, a domain's LRR is reported as a poor summary of its trajectory
rather than passed over. 0.50 is a reporting threshold and nothing else --
no value is dropped, filtered, or flagged in the CSV on account of it.
```

```text
--- 9.2 animation jobs -------------------------------------------------------
range: "real" | "all" | "groin" | "groin_span" | (gis_lo, gis_hi)
mode:  "position" | "displacement" | "difference"
pad:   half-width in domains, read only by "groin" / "groin_span"

"groin" fans out into one GIF per structure in annotations.groins, so new
structures are picked up without editing this list.
```

```text
`output.make_gifs: false` in hat_run.yaml empties this list rather than
skipping the section 12 call. The shoreline matrix .npy is written by that
same call, OUTSIDE the job loop, and section 12.3's paired groin baseline
and scenario_grid.py both read it -- so short-circuiting the call would
cost the run its matrix, while an empty job list costs only the animations.
```

```text
--- 9.4 validation target -- the surveyed island position in the end year ----
Read from the RAW transect CSVs, not the padded offset files: the padded files
are each zeroed on their own most-seaward domain, so differencing two years
subtracts a constant and flips the sign of the mean, which would put the
target on the wrong side of the model. The raw CSVs share a fixed offshore
datum, so differencing them is a real shoreline change.
```

```text
--- 9.5 report, and the annotation guard ------------------------------------
A swapped-in empty AnnotationConfig would strip every figure's geography
silently; fail here instead, where the cause is one line away.
```

```text
Both live in cascade_pipeline/hindcast.py, imported in section 1, and are
shared with HAT_groin_sweep_worker.py -- which used to carry its own copy.

The split is what lets section 11 hold a built-but-unstepped Cascade: BRIE's
diffusivity and the groin's fillet prediction are only meaningful as initial
conditions, and a prediction printed after the run is not one.

`run_years` is TRANSITIONS, not states. See the module comment there for the
off-by-one this replaced, which ran 19 updates for a 20-year period while
dividing by 20.
```

```text
--- dune rebuild thresholds, m MHW ------------------------------------------
Both are floored by roadway_manager on the first step:
dune_design_elevation  = max(passed, BermEl * 10 + 1.0)
dune_minimum_elevation = max(passed, BermEl * 10 + 0.3)
Stated in the documented unit rather than the original's 0.01 "# dam", which
was the wrong unit for the argument and far below the floor anyway.
```

```text
--- wave climate: one configuration, no sweep -------------------------------
A sweep over these is a separate script (scripts/sensitivity_analysis), which
drives this file once per cell through the environment. All four were
literals here until 2026-09-01; they read hat_run.yaml's `physics` block now
so that sweep never has to edit this source between cells.

FIXED_ is kept in the three names: they are still fixed FOR A RUN. Hs reaches
run_index.csv as Hs_m and the other three reach the run metadata, so a
changed value is recoverable from the index rather than only from the file.
```

```text
The run is named for these in 7.5, before they are read here. If the two ever
came apart, the directory would describe a wave climate the model was not
forced with -- the one failure the name token exists to prevent.
```

```text
--- sandbags: off for the hindcast ------------------------------------------
From hat_run.yaml (top-level `sandbags`: a management decision, not a
property of the coast); reaches run_index.csv as sandbags_on.
```

```text
--- run identity, from section 7.5's derived name ----------------------------
OVERWRITE is checked here, before the model is built, so a name collision
costs nothing either way.

False  the matrix default. A directory that already holds a result
stops the run. 12.3 resolves the groin baseline by directory
name, so replacing one run's output silently redefines the
paired run's answer -- which is why this refuses rather than
asks.
True   iterating on one scenario: tweak a value, re-run, read the
figures, tweak again. The directory is EMPTIED and reused, so
it never mixes two trials' files, and run_index.csv's row for
this name is replaced. The previous trial is gone. Set it back
to False before running the matrix.
```

```text
Read before the guard runs: with OVERWRITE=True the guard empties the
directory, and a silent wipe is how a trial gets mistaken for the run it
replaced.
```

```text
CASCADE REWRITES THE PARAMETER FILE while it constructs (brie_coupler's
set_yaml writes the shoreface, RSLR and file paths into it), so the tracked
copy in data/hatteras_init used to be a shared scratch file: two runs at
once corrupted each other's, and a run killed mid-write left it empty
(2026-09-16). Each run now gets its own copy in its directory and hands
CASCADE that path; the tracked file is a read-only template and the tree no
longer goes dirty on every run. Absolute on purpose: brie_coupler joins it
onto datadir with pathlib, which yields the absolute path unchanged, and
Barrier3D opens "<prefix>-parameters.yaml" by that prefix.
```

```text
--- 11.3 pre-run diagnostics ------------------------------------------------
Reads the constructed model rather than the yaml: the active profile height,
the dune thresholds as roadway_manager will apply them, and r_ipl at the groin
cell. The corrected groin's trapped volume depends on the drift arriving, so
there is no a-priori fillet to write down; the extent is measured after the
run against the paired no-groin run.
```

```text
r_ipl is read at the groin cell's own starting angle when a groin is attached,
not shore-normal. Under option A BRIE's diffusivity at 0 deg is negative
(-119 m2/yr), which the r_ipl report cannot use, while the real coast at
GIS 6 sits near -12 deg where it is positive. x_s is a uniform base plus
island_offset, so the offset alone gives the angle BRIE's first solve reads
(forward difference, as brie.py). No groin: shore-normal, as before.
Found 2026-09-29: hard-structures/groin/2-module-tests/2-solver-audit/
2026-09-29-option-a-real-planform/.
```

```text
Verification comes before figures, deliberately, and the nourishment check
ASSERTS. A fill that never reached the model invalidates the run, and that
exact failure went unnoticed for a long time precisely because it produced
plausible-looking output. Nothing is lost when it fires --
run_cascade_simulation saves the model and its logs before returning.
```

```text
--- 12.1 run ----------------------------------------------------------------
Barrier3D seeds x_s_TS with one entry and appends one per update, so a model
that has already been stepped has more than one. Stepping it again would run
past TMAX and raise from deep inside Barrier3D.
```

```text
The denominator that produced a 5% bias in every earlier run. Checked, not
trusted -- see the run_years comment in cascade_pipeline/hindcast.py.
```

```text
TWO ESTIMATORS OF THE SAME THING, and they are not interchangeable.

change_rate  (x[-1] - x[0]) / RUN_YEARS. A net displacement over a span.
Reads two of the RUN_YEARS+1 states, so a single-year
excursion still present in the final state arrives at full
amplitude. Kept because it is the only column that exactly
conserves the period's net shoreline movement, and because
every run before 2026-08-22 was scored on it.
model_lrr    OLS slope through every state. This is what section 8's
target IS -- CoastSat's transect_lrr_full.csv holds a
per-transect OLS slope -- so this is the column the figures
and the skill metrics use.

The distinction is not cosmetic here. A nourishment fill enters as an
instantaneous step in x_s, and BRIE's Crank-Nicolson alongshore solve
answers a step with a grid-scale mode that alternates sign along the coast
and decays only ~39%/yr (see compute_lrr). The 2022 Avon and Buxton fills
are two years from the end of the 2004-2024 run, so an endpoint rate
reports that ringing as a +-1.7 m/yr sawtooth through GIS 6-15 and 22-28.
The LRR reports about a quarter of it, and the residue is the real fill
edge rather than the solver.
```

```text
One name for whichever estimator section 12.5 draws, resolved once here so
no figure call picks its own.
```

```text
Both estimators ship, so a consumer states which one it wants rather than
inheriting whichever the pipeline happened to write. lrr_r2 rides along
because a slope through a domain that stepped rather than trended is a
summary worth flagging at the point of use -- the nourished domains come
out near 0.00.
```

```text
--- skill against the section 8 target --------------------------------------
Both spans are reported. The end domains carry the locked source/sink values
(tens of m/yr), so an island-wide RMSE for a calibBE or edgeBE run is
dominated by two domains that were pinned rather than predicted -- and a
zeroBE run has no such term. Ranking the presets on the island-wide number
alone would mostly rank their boundary treatment.
SKILL is the LRR one: the target is an LRR, so this is the only pairing
that compares like with like. SKILL_ENDPOINT is kept beside it because
calibBE and the groin M were fit before this distinction was drawn, and a
preset's provenance is unreadable once the metric it was fit on stops being
recorded. Expect SKILL to be slightly WORSE than SKILL_ENDPOINT on most
runs -- modelled erosion decelerates across a period, so an all-years slope
is more erosive than an endpoint difference. That is the estimator changing,
not the model.
```

```text
The interior number is GIS 2-89 against the SURVEYED target in every
geometry (2026-09-16), so an extended run and its 90-domain baseline are
graded alike; in the base geometry that is exactly the margin-1 interior
and nothing changes. The extended reach's own interior, its two ends
dropped, is kept beside it as the "reach" interior.
```

```text
Where a straight line is a poor summary of the modelled trajectory, so the
LRR is read with that in mind rather than silently. A nourished domain sits
flat for most of the period and then steps, which no slope describes well.
r2 is a variance ratio, so a domain that barely moved lands here too -- the
two cases are told apart by whether the domain took a fill.
```

```text
--- run metadata: the scenario, plus what distinguishes this run -----------
One structure renders both files: the .txt to read, the .json to parse.
Built together so they cannot disagree -- the earlier version emitted only
prose, so anything downstream had to re-parse it.
```

```text
--- what the source/sink preset actually was -------------------------
The preset's name and its numbers are separate facts.
HATTERAS_BE_RATES_EDGE is a slice of HATTERAS_BE_RATES_CALIBRATED at
the end domains, so editing an end domain to test edgeBE changes
calibBE too -- and both keep their names. Recording the values is what
makes two trials of one preset tell themselves apart afterwards.
```

```text
The three that drift silently across a multi-day batch: which
Cascade class ran, which extractor version built the init surface,
and which commit produced both.
```

```text
BOTH are recorded, and the product is not optional. Version numbers
restart at v1 per product (2026-08-26), so "v1" alone no longer
identifies a topography - 1984-start/v1 and 2004-start/v1 are
different surfaces. A run that logs only the version cannot be traced
back to the arrays it read.
```

```text
The island offset is versioned too (2026-09-15) and the run name
does not say which one was read, so it is recorded here.
```

```text
The reach (2026-09-16): "base" is GIS 1-90; an extended geometry
is not a matrix run and its rows must say so.
```

```text
Barrier3D is a separate repository, installed editable: its branch
is part of the model (2026-09-24, the route_overwash fix).
```

```text
What the model was actually built with, read off the constructed
model rather than the template (2026-09-28).
```

```text
--- cross-run index: one row per run, for comparing the matrix --------------
One row per run; a re-run replaces its row because the rebuild reads disk.
This is what makes 12 runs comparable without opening 12 metadata files.
At OUTPUT_ROOT rather than the period directory: run_name already
carries the period stem, so one file covers the whole matrix and the
two periods stay comparable in a single table.
```

```text
f was absent from this index until 2026-08-31, so no groin run before
that date records which deterioration fraction it used -- and the pair
is quoted as (M, f), not as M alone. The seed runs turned out to be
f = 0.9, recoverable only from a comment in
be_zone_residual_fit.py. be1 needs no column: it is already here
as be_rate_gis1_m_yr.
```

```text
Where the run is filed, spelled out rather than inferred from the
path, so a row and a directory match without reconstructing either.
```

```text
The row goes INTO the metadata, under "index row", and run_index.csv is
REBUILT from every run's metadata rather than appended to (2026-09-16).
Appending was the one shared write that forced runs to be serial, and a
rebuild after a move produced skeleton rows; a derived file has neither
problem, and two runs finishing together both write the same table.
```

```text
--- 12.5 figures -------------------------------------------------------------
Site config arrives via section 9's bundles; omitting them would silently
strip the geographic layer (DEFAULT_ANNOTATIONS is empty).
```

```text
Section 9.4's validation target, resolved now that the run has a year 0.
Buffers are NaN, so the comparison below is over the real domains only.
```

```text
Roadway relocations, so the shoreline GIFs mark the year the road moved.
Assembled here rather than inside the plotting module: the module takes an
array and should not know that a relocation lives on a RoadwayManager.
None when the run has no roadways, which draws no markers and leaves the
GIFs byte-for-byte as they were.
```

```text
NOT `getattr(...) or []`: the attribute is a numpy array and
`array or []` raises "truth value of an array is ambiguous".
```

```text
DUNE_MINIMUM_ELEVATION_M = 0.0   # rebuild trigger: let CASCADE's berm floor
                                 # govern, explicitly rather than by accident
```
