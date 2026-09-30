# hatteras_ms — the hindcast, and everything that reads its output

Four kinds of thing live here, and until 2026-09-13 they were one flat list of
31 items sharing a `HAT_` prefix that sorted them without separating them.

```
HAT_hindcast_1984_2024.ipynb   THE RUN. The notebook is authoritative;
HAT_hindcast_1984_2024.py      the .py is its headless mirror, with no
                               features of its own. A change goes in BOTH.
HAT_hindcast_config.py         which run happens, and where the value came
hat_run.yaml                   from: env > yaml > the default in the module
HAT_run_all.py                 the batch driver: the matrix, then the sweep
HAT_hindcast_methods.md        the written method
HINDCAST_PLAN.md               planning notes: the order the notebook builds
                               the run in. Was `HAT_hindcast_plan`, a file
                               with no extension that this line described as
                               a folder (2026-09-22)

tools/        read or repair the run record; none of them run the model
experiments/  one-off studies, each a run-it then plot-it pair
figures/      figures built from finished runs
groin-sweep/  the M and f fit, its own config and worker
old_drafts/   superseded
old_versions/ superseded runners, and the inherited driver that predates them
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

`../Barrier3D-prefix-ce36866` is a detached worktree of the code **before** the route_overwash fix. It is kept only for the storm-duration cause test (`output/raw_runs/experiments/storms-and-overwash/2026-09-28-storm-max-duration/`) and can be removed with `git -C ../Barrier3D worktree remove ../Barrier3D-prefix-ce36866`.

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
* `superseded_20260907/` -- the 1984 seaward row-insert set, which cannot be
  re-run: the topography layers it studied were deleted, its output folders
  are empty, and its arms are not in the run tree. Kept because the four
  scripts are the only record of how it was driven. Its `WHY.md` has the
  evidence.

A driver spawns the runner as a subprocess with `HAT_IGNORE_SETTINGS=1`, so
whatever is sitting in `hat_run.yaml` cannot reach an experiment.

## figures/ — moved

The five figure scripts that were here (`hindcast_final_figure_lowess`,
`scenario_grid`, `rerender_run_figures`, `planview_evolution_gif`,
`gis11_relocation_drown_figure`) and their `superseded_20260914/` moved to
`scripts/figure_making/model_output/` on 2026-09-18, so every figure script is
in one tree.

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
The decided pair, 2026-08-30. M from period 1 (D4-D8 demeaned), f from
the 1967-2018 rig -- NOT jointly fitted. These defaults matter only when
hat_run.yaml is absent or omits the block; they used to read 50.0 / 0.9,
which were placeholders and f = 0.9 was never fitted at all.
```

```text
WHICH GROIN, 2026-09-29: "dipole" (GroinCallback, +/-M a year) or
"blocking" (BlockingGroinCallback, intercepts a fraction b of the
alongshore transport at the face). b = 0.6 is the option A emulator's
best joint fit with f = 0.3 under the instant 2004 failure, not yet
confirmed in the full model.
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
    groin_trapping_rate_m_yr: M, the groin amplitude knob.
    groin_deterioration_fraction: f, the post-deterioration floor as a
        fraction of M.
    groin_kind: "dipole" or "blocking" -- which groin form is attached.
    groin_blocking_fraction: b, the blocking groin's intercepted fraction
        of alongshore transport; f applies to it as it does to M.
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

End-to-end unattended driver: the comparison matrix, the groin sweeps, the joint fit.

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
The sweep, its config and the joint fit live in groin-sweep/. That is a
plain directory and not a package -- the hyphen in the name makes it
unimportable as one -- so its path goes on sys.path and the config is
imported by module name, exactly as it was when the four files sat beside
this one. Without this the driver dies on ModuleNotFoundError before it
reaches stage 1.
```

```text
The matrix's preset axis is WIDER than the sweep's. PRESETS is the pair the
groin sweep fits M and f against, and it stays the default because the
comparison matrix exists to feed that fit. But a hindcast run is valid under
any preset the site config defines, and calibBE is the only one that puts a
source/sink term on the interior domains -- so it is reachable through
--presets even though no sweep is fitted for it. Order follows PRESETS first
so the default schedule is unchanged.
```

```text
The seed runs exist only to give the drift guard a reference. Their M and f
are provisional and are NOT a prior: stage 6 re-runs both cells at the
fitted values, so nothing downstream ever reads these numbers.
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
The two seed cells already hold a run at provisional values.
Only those are overwritten; the other 18 keep the collision
guard that stops a matrix run from silently replacing another.
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
reloc off: the seed exists only to give the sweep's drift guard a
rate curve to difference against, and the sweep's own cells are
built without relocations.
```

```text
Period 1 edgeBE is 215 of the 344 cells; running it first means the
longest leg starts earliest and a morning check finds it either done or
visibly progressing, rather than not yet begun.
```

```text
Read by HAT_groin_sweep_worker, and by sweep_output_dir, which
sends the results to their own directory.
```

```text
Streamed to a log rather than captured: a 215-cell sweep prints a
progress line per cell, and holding four hours of that in memory to
write it at the end would lose all of it if the driver were killed.
```

```text
Order follows PRESETS, not the order they were typed, so two
invocations naming the same presets schedule them the same way.
```

```text
Said once, here, rather than left for stage 4 to discover: a preset with
no sweep grid still runs every matrix cell, and only the stages that
need a fitted (M, f) skip it.
```

```text
Same treatment as --presets: validated against the runner's own table
and reordered into its order, so two invocations naming the same
scenarios schedule them identically.
```

```text
Stage 6 needs fitted (M, f). Stage 5 returns them, but running stage 6
ALONE used to hand it {} and it skipped every groin run with "no fitted
(M, f)" -- 21 of 22 on 2026-08-24. Falling back to the file makes stage 6
runnable on its own, which matters whenever the values in joint_fit.json
were NOT produced by the joint fit's own ranking: re-running stage 5 to
satisfy stage 6 would overwrite them with the ranking's answer.
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
should and leaves the rest alone.

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
    fits: {preset: {"M": ..., "fraction": ...}} for the groin stage.

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

### groin-sweep/HAT_calibration_summary_figures.py

Summary figures documenting the 2026-08-30 calibration, and how it was tested.

From the script's original header:

```text
Summary figures documenting the 2026-08-30 calibration, and how it was tested.

WHY THESE EXIST
    The reasoning behind the source/sink and groin calibrations lives in prose,
    in comments in scripts/site_layer/hatteras_site_config.py. That is the right home for
    the conclusions, but two things were not recorded anywhere at all:

      * the BE convergence sequence. `--overwrite` REPLACES a run's row in
        run_index.csv, so only the last pass survives there, and
        convergence_history.json still carries 2026-08-24 baselines from the
        pre-restructure topography. The pass-by-pass numbers existed only as
        text typed into a comment.
      * the three-target comparison. That the groin's fitted M is set by the
        CHOICE OF TARGET rather than by the data is the methodological result
        of the exercise, and nothing on disk showed it.

    These figures are the durable record of both.

WHAT EACH ONE SHOWS
    fig_three_targets.png    the same 61 sweep cells scored three ways. The
                             fillet says M = 95, the D1-D12 profile says M = 0,
                             D4-D8 demeaned says M = 60. Identical model runs.
    fig_Mf_identifiability.png  D4-D8 demeaned RMSE over the (M, f) grid, with
                             iso-M*f contours. Built to TEST the claim that
                             "only the product M*f is identified" -- and it
                             refuted it: corr(RMSE, M*f) = -0.07 against
                             +0.61 for M and -0.49 for f, and equal-product
                             cells score 10.4 to 12.5 m. But the REPLACEMENT
                             claim ("M and f each weakly constrained") was
                             also wrong: per GROIN_PLAN.md the invariant is
                             period-1 cumulative trapping, M(15.5 + 4.5f).
    fig_top_profiles.png     the top cells and the no-groin baseline against
                             the observed change profile, fit window marked.
    fig_be_convergence.png   interior RMSE per calibration pass, both periods,
                             with the GIS 90 re-solve marked.
    fig_period2_and_bug.png  why period 2 is not fitted, and what the
                             topography-product bug was worth.

A CAVEAT THAT NO LONGER APPLIES, WITHDRAWN 2026-08-31
    This said the D1-D12 panel came from the 40-year window on the PRE-FIX
    topography while the other two were period 1 on the corrected one, so
    the three were not a controlled comparison. True when written; not true
    now. The fullperiod sweep was re-run 2026-08-30 18:20, five hours AFTER
    the worker topography fix (562c75c, 13:01), and this figure was rebuilt
    at 23:52 from those cells.

    VERIFIED, not assumed: re-running cell M60_f0.50 through the worker
    reproduces its stored result to 8.3e-05 on rates of ~2.9 m/yr -- the
    SAME noise floor a period-2 cell shows against itself (1.0e-04 between
    two re-runs), which the 1984-2024 window inherits because it contains
    period 2. A wrong-island run would differ in the first or second
    decimal, not the fifth.

    All three panels are on the same corrected topography.

STYLE, 2026-09-11
    Under the project house style (`scripts/site_layer/hat_figure_style.py`), which
    replaced this file's own INK/MUTED/ACCENT/FOIL palette and its local
    rcParams block. Three consequences worth knowing before reading an older
    copy of these images side by side with a new one:

      * The canvases were 8.6-15 in wide and are now a 190 mm printed column,
        so the type is the size it claims to be on a page.
      * The two calibration periods are drawn in the house VINTAGE pair -- the
        earlier period red, the later blue -- in every panel that shows both.
        They were pink/teal here and pink meant "the answer" elsewhere in the
        same figure set.
      * The suptitles, the italic per-panel verdicts and the footnote
        paragraphs are off the canvas and in CAPTIONS.md beside the images.
        The five captions carry every sentence they used to, so nothing in the
        argument was dropped to make room.

Usage:
    python HAT_calibration_summary_figures.py

Writes output/calibration/groin/figures/ (untracked; the PDFs beside the PNGs are).
Reasoning and results: CALIBRATION_FIGURES.md, beside this file.
```

Notes that were in the code:

```text
The LIVE full-period sweep, not the 2026-08-28 archive. This pointed into
`superseded_20260828/` (now `output/archive/2026-08-28_full-tree/`), which is a "do not use for analysis" tree, and
that one line was the only thing keeping 12 GB of superseded output
undeletable. The live file carries the same 16 columns and the same 43 rows.
```

```text
Semantic colours, from the house palette. ACCENT is the target that works
and the answer it gives; BASE is a target that fails or a baseline; REF is a
reference construction laid over the data (the iso-product curves, a
tolerance line); the VINTAGE pair is the two calibration periods.
```

```text
Five top cells as one family: a ramp of the accent, so they read as
variations of the same thing rather than five unrelated series. viridis was
used here until 2026-09-11 and put the best cell in the same green as the
reference curves in the neighbouring figure.
```

```text
Rates are m/yr; the observed target is CHANGE over the period, so the
model side is scaled by the period length. Demeaned because a uniform
level offset in the groin's neighbourhood belongs to the source/sink
term, not to the groin -- what the groin must get right is the shape.
```

```text
LANDWARD-POSITIVE, so erosion is UP and this panel reads as a plan view,
matching the gifs and the other profile figures. Both sources are
SEAWARD-positive, so both are negated -- at the PLOTTING layer only.
score_d48 is computed upstream and is unaffected.
fig_three_targets needs no flip: it plots error against M, not a profile.
```

```text
Darkest first: `top` is sorted best to worst, and the caption says the
best cell is the darkest.
```

```text
Outside the axes: seven entries inside the panel covered the observed
line across D1-D3, which is the part of the profile the caption is about.
```

```text
Transcribed from the pass-by-pass table in hatteras_site_config.py. NOT
derivable from run_index.csv: --overwrite replaces a run's row, so only the
final pass survives there. This figure is the durable record.
```

### groin-sweep/HAT_d4d7_window_figure.py

The D4-D7 window: the largest RMSE gain of any window tested, and why that is misleading.

From the script's original header:

```text
The D4-D7 window -- largest RMSE gain of any window, and why that is misleading.

D4-D7 returns M = 70 with a 3.09 m gain, the largest of the eight windows
tested on 2026-08-30. This figure exists to show that a larger gain is not a
better fit: it decomposes the gain per domain, and the decomposition says the
groin improves the domain OUTSIDE the dipole and degrades the downdrift domain
the structure actually acts on.

STYLE, 2026-09-11
    Under the house style (`scripts/site_layer/hat_figure_style.py`). The footnote
    paragraph is now the caption in output/calibration/groin/figures/CAPTIONS.md, the
    canvas is a 190 mm printed column instead of 14.5 in, and the bars are the
    house BASE grey for the baseline against ACCENT purple for the run under
    test. The per-domain deltas keep a green/red split, which is the one place
    in this figure set where those two mean better and worse rather than the
    1984/1997 vintages -- nothing on this figure is a vintage, and the whole
    point of panel (b) is which domains moved the wrong way.

Writes output/calibration/groin/figures/fig_d4d7_window.png (and .pdf)
```

Notes that were in the code:

```text
BETTER / WORSE on panel (b) only; ACC is the run under test, FOIL the
baseline it is scored against.
```

```text
LANDWARD-POSITIVE from here on, so erosion is UP and panel (a) reads as a
plan view, matching the gifs. rate_D* and observed_change_profile are
both SEAWARD-positive at source, so both are negated -- at the PLOTTING
layer only. The scoring pipeline (FLIP_SIGN_MODEL, the sweep worker) is
untouched, and panel (b) is unchanged either way because it plots
|residual|.
```

```text
Panel (a)'s handles only: panel (b)'s bars are the same two things in the
same two colours, and a figure-wide legend listed each of them twice.
```

### groin-sweep/HAT_fullperiod_figures.py

The four sweep outputs for the continuous 1984-2024 groin calibration, the set the 1967 rig produced.

From the script's original header:

```text
The four sweep outputs for the continuous 1984-2024 groin calibration.

Deliberately the same set the 1967 rig produced, because that set answered the
question the scalar-fillet figures could not: a heatmap showing whether the
optimum is INTERIOR, and profile plots showing whether the winning cell has the
right SHAPE and not merely the right magnitude at one point.

    heatmap.png            profile RMSE over the M-f grid, best cell marked,
                           cells the model refused drawn as gaps rather than
                           silently dropped
    best_fit_profile.png   the winning cell against the observed change profile
    top_n_profiles.png     the best N cells together, so the spread shows how
                           sharply the metric discriminates

SIGN CONVENTION, AND WHY IT DIFFERS FROM THE RIG FIGURE
    These plot SEAWARD-POSITIVE change, because that is what the CoastSat
    chainage target is measured in and converting for display would put two
    conventions in one workflow. The 1967 rig's figures are landward-positive
    ("+ = landward"), so a curve that rises here falls there. The axis label
    states it on every panel rather than relying on the reader to remember.

WHAT AN INTERIOR OPTIMUM WOULD MEAN
    Every earlier attempt at this calibration railed: the best cell sat on a
    grid edge, which means the search wanted to keep going and ran out of grid,
    not that it found a minimum. A best cell with neighbours on all four sides
    is the evidence that this window and this metric can actually identify the
    pair. The heatmap flags the outcome either way.

Usage:
    python HAT_fullperiod_figures.py [--top-n 5]

Reads  output/calibration/groin/fullperiod_1984_2024/results.csv
Writes output/calibration/groin/fullperiod_1984_2024/figures/
```

Notes that were in the code:

```text
House colours (2026-09-11): observations in INK; the structure is a
muted guide line, not the vintage red it used to be drawn in.
```

```text
Two greys, not two hues: these bands say WHERE the structure is, and
the colour on this figure belongs to what is plotted.
```

```text
The caveat belongs in the TITLE, not a footnote. This window nets period
1's fillet build against period 2's collapse, so a module whose trapping
is bounded at >= 0 can never win here regardless of how it is scored. The
sweep bounds M from above; it does not test whether a groin operated.
Fitting is done on period 1 -- see why_M60_f06.png.
```

### groin-sweep/HAT_fullperiod_sweep.py

Groin sweep over one continuous 1984-2024 window, scored on the change profile.

From the script's original header:

```text
Groin sweep over ONE continuous 1984-2024 window, scored on the change profile.

This is the 1967-rig method applied to the production geometry and the full
hindcast span. It exists because the two-period, scalar-fillet approach failed
in a specific and repeatable way, and this one did not.

WHAT WENT WRONG BEFORE, AND WHAT IS DIFFERENT HERE

    scalar target -> profile target
        Matching one number (the fillet) left a RIDGE of equally good (M, f)
        pairs; the ranking then returned whichever cell sat at a grid edge.
        Scoring the per-domain CHANGE PROFILE constrains shape as well as
        magnitude. That is what gave the rig an interior optimum in both knobs.

    two 20-year windows -> one 40-year window
        The fillet builds through 1984-2004 (+52 m) and declines through
        2004-2024 (-76 m). Each 20-year window sees only one of those, so M and
        f trade off within it. One continuous run sees both, and they constrain
        different combinations of the pair.

    coarse then fine
        A wide coarse pass, then a fine pass centred on its best cell. M is
        capped at 80: on the rig every cell at M >= 100 DROWNED the barrier
        partway through, and M >= 70 produced RMSE an order of magnitude above
        its neighbours. Sweeping into that region wastes runs on a model that
        refuses the parameter.

ASSUMPTIONS, STATED

    static background erosion
        `background_erosion` is written once at construction, so a 40-year run
        carries ONE field. The two periods' calibrated edge values differ in
        sign (be1 -41.8 vs +50.3), but over the whole window they largely
        cancel: the surveyed change at D1 is -3.2 m in 40 years. be1 is
        therefore SOLVED against that, not inherited from either period.

    spliced storms
        1984-2003 from the 1984-2004 series, 2004-2024 from the 2004-2024 one,
        with the model-year index remapped. Checked at the join: wave height,
        runup and period agree to within 2%.

    this is a rig, not a hindcast
        One static background field cannot reproduce a field that reverses
        mid-window. The run is fit for calibrating a LOCAL structure, where a
        smooth regional field largely cancels in the profile's shape. It is not
        a substitute for the two-period hindcast.

Usage:
    python HAT_fullperiod_sweep.py [--workers 6] [--be1 -5.0] [--stage coarse|fine|both]

Writes output/calibration/groin/fullperiod_1984_2024/:
    results.csv          one row per cell, ranked by profile RMSE
    <combo>/             shoreline matrix per cell
    figures/             heatmap, best-fit profile, top-N profiles
```

Notes that were in the code:

```text
Scoped by wave climate: the worker reads HAT_SWEEP_HS, and without this a
run at another Hs would write its cells over the 2.5 m results that the
recorded -15.3 m decay and the M upper bound both came from.
```

```text
Capped at 80: M >= 100 drowned the barrier on every rig cell, M >= 70 went
unstable there. The production grid is better buffered so the threshold may
differ, but sweeping far past a known failure mode buys nothing.
```

```text
An empty list gives a DataFrame with no columns, and sorting on a column
that does not exist raises KeyError -- which is what a first invocation,
or a --collate-only before anything has run, would hit.
```

### groin-sweep/HAT_fullperiod_target.py

The observed shoreline-change profile for a continuous 1984-2024 window.

From the script's original header:

```text
The observed shoreline-change profile for a continuous 1984-2024 window.

WHY A PROFILE AND NOT A SCALAR
    Every scalar target tried on this calibration produced a RIDGE: many (M, f)
    pairs matched one number equally well, and the ranking then picked whichever
    cell sat at a grid edge. A per-domain CHANGE PROFILE constrains shape as
    well as magnitude, which is what separated the knobs on the 1967 rig -- the
    one configuration that produced an interior optimum in both.

WHY A CONTINUOUS WINDOW
    The fillet BUILDS through the first half of 1984-2024 (+52 m by 2004, from
    the fixed-datum wet/dry surveys) and DECLINES through the second (-76 m by
    2023). A 20-year window sees only one of those, so M and f trade off inside
    it. One 40-year run sees both, and they constrain different combinations.

HOW THE TARGET IS BUILT
    From `groin_analysis_chainage_all.csv` -- 793,259 CoastSat shoreline
    POSITIONS in the window, every one of the 90 domains carrying at least 20.
    For each domain an OLS line is fitted to chainage against decimal year and
    evaluated at both ends; the difference is that domain's change.

    FITTED, NOT DIFFERENCED. Two year-bins differenced give a noisy endpoint
    estimate -- at D1 that route disagreed with the published LRR by a factor
    of two. The OLS fit is also the estimator the rest of the project uses, so
    the target and the model rates are the same kind of number.

SIGN
    Chainage is SEAWARD-POSITIVE, verified by correlating its 1984->2004 rate
    against the published CoastSat LRR (+0.74). Barrier3D's x_s is
    landward-positive, so the model side is negated before comparison.
```

<details><summary>Function notes (the original docstrings)</summary>

**`observed_change_profile()`**

```text
Observed shoreline change per domain over the window, in metres.

Args:
    start, end: Window bounds, inclusive.
    domains: GIS domain ids to return.

Returns:
    A dict of domain -> change in metres, SEAWARD-POSITIVE.

Raises:
    FileNotFoundError: If the chainage table is absent.
    ValueError: If any requested domain has too few observations.
```

**`model_change_profile()`**

```text
Modelled shoreline change per domain, in metres, seaward-positive.

Args:
    shoreline_m: [state x padded domain] array, metres, landward-positive.
    geometry: DomainGeometry for GIS -> pad translation.
    domains: GIS domain ids to return.

Returns:
    A dict of domain -> change in metres.
```

</details>

### groin-sweep/HAT_fullperiod_windows.py

Does a groin signal survive in a narrow window? The full-period sweep rescored several ways.

From the script's original header:

```text
Does a groin signal survive in a NARROW window? Rescores the sweep several ways.

THE QUESTION
    On the full D1-D12 window the no-groin cell fits best and RMSE rises
    monotonically with trapping. But that window is dominated by two features
    the groin module cannot produce:

        D2-D4   observed ACCRETES (+3 to +7 m over 40 years); the model erodes
                there. Cape Point, which the parameterisation does not
                represent.
        D6-D7   observed erodes 48 and 63 m, a trough peaking one domain NORTH
                of the structure. A groin pushes D6 seaward, i.e. the wrong
                way, so raising M makes the single largest misfit worse.

    Neither is a groin signal, and together they swamp one. The question this
    file answers is whether a groin signal exists underneath, in the few
    domains the structure actually reaches.

TWO SCORES PER WINDOW, AND WHY BOTH ARE NEEDED
    raw         RMSE of modelled against observed change. As the window
                narrows this is increasingly dominated by a constant OFFSET:
                if the model sits 20 m low across D4-D8, raw RMSE is ~20 m
                whatever M does, and the groin's contribution is invisible.

    detrended   A straight line is fitted across the window and removed from
                BOTH profiles first, leaving only SHAPE. A groin's dipole is a
                shape -- a step across the structure -- so this is the score
                that can actually see it. A cell that gets the local shape
                right while sitting on a biased background shows up here and
                nowhere else.

    Reported together on purpose: a groin that improves the detrended score
    while leaving the raw one unchanged has explained the local pattern
    without fixing the level, which is a precise and reportable statement
    rather than a pass or a fail.

WHAT WOULD COUNT AS A SIGNAL
    An INTERIOR optimum at M > 0 that beats the M = 0 baseline. If every
    window at both scores still prefers M = 0, that is a much stronger
    negative result than one window alone -- it says the groin's effect is not
    merely swamped at reach scale but absent at the scale it acts on.

Usage:
    python HAT_fullperiod_windows.py

Reads  output/calibration/groin/fullperiod_1984_2024/results.csv
Writes fit_windows.csv beside it, and prints the comparison.
```

Notes that were in the code:

```text
The groin field occupies D6; the pair is D5/D6. Measured influence is 2.25 km
updrift (~4.5 domains) and ZERO downdrift, so windows widen mostly northward
in spirit even though they are drawn symmetrically for simplicity.
```

```text
Two points define the line exactly, so detrending would zero them
and every cell would score identically. Return as-is and let the
caller see the raw score instead of a meaningless zero.
```

<details><summary>Function notes (the original docstrings)</summary>

**`demean()`**

```text
Removes a CONSTANT offset, preserving every gradient and step.

This is the score that matches how the model is actually applied: a
per-domain source/sink correction is fitted afterwards, so a uniform level
error in the groin's neighbourhood is absorbed downstream and is not the
groin's job to fix. What the groin must get right is the SHAPE.

Deliberately weaker than `detrend`. Removing a linear trend across a short
window also removes the gradient that a dipole PRODUCES -- fit a line
through D4-D8 and the step across D5/D6 is partly absorbed into it, so a
working groin can be scored as no better than none. Subtracting the mean
cannot do that.
```

</details>

### groin-sweep/HAT_groin_choice_figure.py

Why M = 60, f = 0.6: the period-1 fit that supports it.

From the script's original header:

```text
Why M = 60, f = 0.6 -- the period-1 fit that supports it.

REPLACES AN EARLIER VERSION OF THIS FIGURE, and the reason matters.

    The previous version argued that the fillet metric was unsatisfiable, so M
    had to be bounded by physics (sediment budget, barrier stability) with the
    reach fit refining inside those bounds. Two things were wrong with it:

      1. The premise is superseded. M = 60 is now supported by a DIRECT FIT on
         period 1 over D4-D8 -- a 25% improvement on no groin -- rather than by
         an argument from constraints.

      2. It shaded "unstable M >= 70" and "barrier drowns M >= 100" on a
         PRODUCTION-geometry plot. Those thresholds were measured on the
         41-domain rig and do not transfer: all 36 production cells, including
         the full M = 70 and M = 80 rows, ran clean. Shading them was
         misleading.

WHY PERIOD 1, AND WHY D4-D8
    Period 1 is the only window in the hindcast where the observed gap between
    the groin's two flanks WIDENS (+52 m). That is the only behaviour the
    module can produce -- trapping is bounded at >= 0, so it can widen the gap
    or stop widening it, never close it. Period 2 closes, and the continuous
    1984-2024 window nets the two against each other, so neither can fit a
    groin however it is scored.

    D4-D8 excludes D1, where the cape's shoreline change over period 1 is
    81-104 m -- roughly five times the groin's ~17 m signal. On the full window
    with a raw score the cape swamps the groin and no-groin wins by 0.18 m; on
    D4-D8 the groin wins by 5.06 m. The signal was always there; the window was
    hiding it.

WHY THE SCORE IS DEMEANED
    A uniform level offset in the groin's neighbourhood is absorbed by the
    source/sink calibration that runs afterwards, so correcting it is not the
    groin's job. What the groin must get right is the SHAPE. Demeaning removes
    a constant and keeps every gradient -- unlike a linear detrend, which would
    partly absorb the dipole's own gradient and hide a working groin.

WHAT THE FIGURE DOES NOT CLAIM
    Not that (60, 0.6) is uniquely determined. Fourteen cells lie within 0.5 m
    of the best, spanning M = 40-95 and f = 0.4-1.0. The ridge is in
    PERIOD-1 CUMULATIVE TRAPPING, M(15.5 + 4.5f) -- not in M*f, which
    fig_Mf_identifiability.png tested and refuted (corr(RMSE, M*f) = -0.07).
    See CALIBRATION_FIGURES.md. The chosen pair is not the top-scoring one
    (M = 50, f = 1.0 scores 11.59 m, 0.10 m better) because f = 1.0 asserts the
    groin never deteriorated, which the GIS record contradicts outright. Within
    the tied band, (60, 0.6) is the cell that carries a real deterioration floor
    and still stays inside the affordable drift. The right panel draws the band rather
    than a single star, because a star would assert precision the data does not
    support.

Usage:
    python HAT_groin_choice_figure.py

Writes output/calibration/groin/figures/why_M60_f06.png
```

Notes that were in the code:

```text
LANDWARD-POSITIVE, so erosion is UP in the left panel and it reads as a
plan view, matching the gifs and the other profile figures.
observed_change_profile is SEAWARD-positive, so it is negated; x_s is
landward-positive already, so the negation on `change` below is gone.
Both flip together, so every RMSE here is unchanged.
```

```text
The two flanks of the structure, named rather than colour-coded: the
blue/red washes here were the vintage pair doing a third job.
```

### groin-sweep/HAT_groin_full_life_figure.py

The groin module over the structure's whole life, against every survey.

From the script's original header:

```text
The groin module over the STRUCTURE'S WHOLE LIFE, against every survey.

WHY THIS FIGURE EXISTS
    `HAT_groin_timeseries_check.py` plots the two hindcast windows, and both
    of them start 15 years after the groin went in. They therefore show the
    fillet's DECAY and never its CREATION -- which is why they cannot see f,
    and why they make the module look worse than it is. The 1967-2018 rig is
    the only window that contains the build phase, the 1996 repair, the 2003
    storm damage and the decline that follows. This figure is that window.

    It is also the figure that justifies f = 0.6: the deterioration ramp is
    visible here and nowhere else.

WHAT IS PLOTTED
    (a) FILLET THROUGH TIME. The surveyed fillet (D5 - D6 offset against the
        fixed 1967 datum, from 24 dated wet/dry surveys) as markers -- markers
        only, because the record is irregular and a joining line would imply
        samples that do not exist -- against the rig's modelled fillet. The
        structure's documented timeline is marked on top.

    (b) WHAT THE MODULE WAS DOING. The trapping rate the module actually
        applied each year, read from the run's own groin_diagnostics.csv
        rather than recomputed. This is the schedule f parameterises: zero
        before install, M while the structure is sound, a linear ramp down
        from the 1996 repair to the 2003 damage, then a hold at M*f.

    Reading the two together is the point. Panel (b) explains the shape of
    the model curve in (a), and the observed peak in (a) at 2004 is what
    fixes the end of the ramp in (b).

HOW TO READ THE RESIDUAL
    The module reproduces the SIGN and the TIMING of the build, and it
    undershoots the AMPLITUDE. That is expected and documented: no admissible
    M matches the fillet on this grid, because the real fillet is ~190 m wide
    against a 500 m domain and the dipole is volume-neutral where the real
    structure is not (observed downdrift extent 0 m, the model's 2,500 m).
    The gap is the part the source/sink calibration and the Cape Point
    dynamics absorb -- not a failed fit.

    Note also that the rig runs a 1967 window off 1984 topography
    (RIG_TOPO_PRODUCT = "1984-start"), a deliberate anachronism accepted
    because the target is a shoreline OFFSET rather than an elevation.

Usage:
    python HAT_groin_full_life_figure.py

Writes output/calibration/groin/figures/full_life_1967_2017.png
```

Notes that were in the code:

```text
The rig pads 11 real domains (D2-D12) with 15 buffer domains either side, so
D2 -> 15 and D5 -> 18, D6 -> 19. This is _gis_to_pad() in
HAT_groin_hindcast_1967_2017.py:76, restated rather than imported because
importing that module builds a CASCADE run.
```

```text
The rig lives in output/calibration/groin_rig/, not output/raw_runs/ (moved 2026-08-31;
under calibration/ since 2026-09-18).
It is a DIFFERENT GRID -- 41 domains against production's 120 -- and M is
grid-specific, so mixing the two invited quoting a rig number as a
production one. raw_runs is production only, and run_index.csv covers it.
```

```text
A rig run directory does NOT name its own parameters -- the sweep writes every
cell into one run name, so whatever survives is the last cell that finished.
On 2026-08-30 this directory was found holding an UNSTABLE M = 70 cell whose
fillet ran away to 444 m, while being named as though it were the calibrated
run. So the applied rate is read back from the diagnostics and checked here
rather than trusted from the label.
```

```text
House colours (2026-09-11): the surveys are INK, the run under test the
ACCENT, the groin-off run BASE grey, and the structure's dated events are
guide lines in muted ink rather than a fourth hue. These were near-black,
an orange, a blue, a pink and a grey chosen in this file.
```

```text
Inside the axes: above the top edge these ran through the panel
title at the printed width.
```

<details><summary>Function notes (the original docstrings)</summary>

**`observed_fillet_by_year()`**

```text
Surveyed fillet against the fixed 1967 datum, {year: metres}.

Same convention as HAT_groin_timeseries_check.py: downdrift minus
updrift, so a rising curve means the updrift side is holding while the
downdrift side retreats -- what a groin builds.
```

</details>

### groin-sweep/HAT_groin_full_life_gif.py

Animated D2-D12 shoreline over the groin's whole life, 1967-2017.

From the script's original header:

```text
Animated D2-D12 shoreline over the groin's whole life, 1967-2017.

WHAT THIS ADDS OVER `HAT_groin_zoom_gifs.py`
    That script animates the two hindcast windows and draws the observed change
    as a FIXED endpoint target, because inside 1984-2004 the observation is a
    single endpoint. Over the full life it is not: the wet/dry record carries
    25 dated surveys between 1967 and 2023 across D2-D12, so THE OBSERVATIONS
    CAN ANIMATE TOO. Each frame shows the most recent survey at or before that
    model year, with the surveys already passed left behind as fading ghosts.

    That is the comparison the hindcast windows cannot show -- the module
    building a fillet from a standing start at installation, holding it, and
    then losing it through the deterioration ramp, against a survey record that
    is doing the same thing on its own clock.

ORIENTATION: EROSION IS UP, SO THE PANEL READS AS A PLAN VIEW
    The y axis is LANDWARD-POSITIVE. A retreating shoreline moves UP the panel
    and an accreting one moves down, so the reader is looking down on the
    island with the ocean below the axis and the island above it. This is the
    OPPOSITE of `HAT_groin_zoom_gifs.py`, which is seaward-positive -- the two
    must not be read side by side without noticing.

    Both source arrays are already landward-positive and are therefore NOT
    negated here:
      `Change_from_wetdry_1967_*.csv`  landward-positive (a rising value is
          retreat) -- which is why the fillet is built as downdrift minus
          updrift throughout this directory.
      `shoreline_matrix.npy`           Barrier3D's x_s, landward-positive and
          already in METRES -- no dam conversion. (A dam->m rescale added to
          the zoom gifs on 2026-08-30 made every curve ten times too large;
          anything that rescales this must be re-checked against the cell's own
          shoreline_change_rate.csv.)

    The fillet is therefore reported as D5 - D6, matching GROIN_PLAN.md, rather
    than the D6 - D5 the earlier seaward-positive revision of this script used.

WHAT IS PLOTTED
      model, M = 60 f = 0.6   the rig's groin run
      model, groin OFF        the paired baseline at the same edge calibration
      observed                the most recent survey <= this year, with its own
                              year and coverage; the previous five as ghosts

    NOT DEMEANED, unlike the zoom gifs. The question here is how well the
    module reproduces the MEASURED POSITION, so the level is left in. A uniform
    alongshore offset is owned by the source/sink calibration, not the groin --
    so read a whole-curve offset as that term's business, and read the SHAPE
    around D5/D6 as the groin's.

Usage:
    python HAT_groin_full_life_gif.py

Writes output/calibration/groin/figures/full_life_gif/shoreline_D2-D12_1967_2017.gif
```

Notes that were in the code:

```text
The rig pads 11 real domains (D2-D12) with 15 buffer either side: D2 -> 15,
D5 -> 18, D6 -> 19, D12 -> 25. The RIG's convention, which differs from
production's -- see HAT_groin_hindcast_1967_2017.py:76.
```

```text
The rig lives in output/calibration/groin_rig/, not output/raw_runs/ (moved 2026-08-31;
under calibration/ since 2026-09-18).
It is a DIFFERENT GRID -- 41 domains against production's 120 -- and M is
grid-specific, so mixing the two invited quoting a rig number as a
production one. raw_runs is production only, and run_index.csv covers it.
```

```text
House palette (2026-09-11), replacing an Okabe-Ito set chosen in this file:
the run under test is the ACCENT, the groin-off baseline BASE grey, the
surveys INK, and the structure a muted guide line rather than a fourth hue.
```

```text
Documented structure history -- GROIN_PLAN.md section 1. The rig's own
install year is 1970 (it keeps 1967-69 as a free control window); the
documented installation is 1969. Both are shown rather than reconciled.
The strip is a ruler, not data: four greys, light to dark as the structure
degrades, so the timeline reads in one glance without spending four hues on
it. It was a grey/blue/orange/red set until 2026-09-11.
```

```text
An all-NaN column for a domain is a real gap in the survey, not an
error -- keep the NaN and let the line break there.
```

```text
---- static furniture, drawn once ------------------------------------
Title and subtitle stay ON the canvas: a GIF is watched standalone, with
no caption file beside it in a viewer. House sizes, not a 16.5 pt banner.
```

```text
The footnote paragraph goes to CAPTIONS.md beside the GIF, written after
the animation is saved.
```

```text
labelpad clears the orientation cues at x = -0.088, which the
24 pt pad of the pre-2026-09-11 layout ran through.
```

```text
Bottom-RIGHT: the structure label occupies the bottom-left, and the
accretion half of the panel is empty in every frame.
```

```text
A three-year segment cannot hold a label; the "installed"
event marker already names that stretch.
```

<details><summary>Function notes (the original docstrings)</summary>

**`observed_by_year()`**

```text
{survey year: change since 1967 over DOMAINS, LANDWARD-positive}.

NaNs are preserved: several surveys cover only 8-10 of the 11 domains, and
interpolating across a gap would invent a shoreline.
```

</details>

### groin-sweep/HAT_groin_joint_fit.py

Fit M and f together from the two periods' sweep surfaces.

From the script's original header:

```text
Intersects the two periods' sweep surfaces to fit M and f together.

WHY THIS IS A SEPARATE STEP
    Neither period identifies both knobs on its own. Cumulative trapping per
    unit M, measured from `cascade.groin.GroinCallback` with the documented
    1969 install / 1996 onset / 7-year ramp schedule:

        1984-2004    16 + 4f      f moves this only 16.0 -> 20.0, so holding
                                  cumulative trapping fixed while sliding f
                                  across its whole range needs just a 25%
                                  change in M. f is nearly free here.
        2004-2024    20f          the run sits entirely past the 2003 ramp,
                                  so M and f enter only as their product.

    One period gives one constraint on two unknowns. Two periods intersect.
    This script does that intersection, per source/sink preset.

HOW be1 IS TREATED
    be1 is a NUISANCE parameter, not a fitted result: it is swept in period 1
    under edgeBE and absent everywhere else. For each (M, f) the period-1
    contribution is MINIMISED over be1 -- a profile likelihood -- and the
    winning be1 is reported alongside. Summing over be1, or fixing it at one
    value, would both charge the groin for background erosion the model was
    free to place elsewhere.

WHAT THE ANSWER WILL LOOK LIKE, AND WHY
    Period 2's observed D6 - D5 is negative (-2.47 m/yr). The source/sink pair
    adds -M updrift and +M downdrift, so the modelled differential is
    non-negative at any M >= 0 and no cell can reach that target. Period 2
    therefore contributes a monotone penalty in M*f, pushing the joint
    solution toward f = 0, and period 1 then sets M through M*(16 + 4f).

    Algebraically, with A the period-1 constraint and B the period-2 one:

        f = 4*B / (5*A - B)

    B pinned near zero gives f near zero. That is a RESULT -- the structure
    stopped trapping after the 2003 storm -- but it is reached by railing to
    a grid edge, so this script flags every fitted value that lands on a
    bound rather than in the interior. A railed value is a bound, not an
    optimum, and must not be quoted as a fitted parameter.

Usage:
    python HAT_groin_joint_fit.py [--preset edgeBE] [--no-figures]

Reads   output/calibration/groin/<period>_<preset>/sweep_results.csv  (all four)
Writes  output/calibration/groin/joint_fit.json    fitted (M, f, be1) per preset
        output/calibration/groin/joint_fit.csv     the full joint surface
        output/calibration/groin/figures/joint_<preset>_surface.png
        output/calibration/groin/joint_<preset>_constraints.png
```

Notes that were in the code:

```text
parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
The guard below is what makes a future move fail here, loudly, instead of
resolving to scripts/scripts and surfacing as a missing data file several
imports deeper.
```

```text
Figures live in a subdirectory; joint_fit.json does NOT. That file is read by
HAT_run_all.py (stage 6 takes its fitted M and f from it) and by
be_zone_residual_fit.py (which uses it to find a groin-aware base run),
both of which pin the top-level path. Moving it would break stage 6 silently.
```

```text
House colours (2026-09-11): the two periods are the vintage pair, the fitted
point is the ACCENT, and the error surface is greyscale so the marks on it
stay findable. MODEL_COLOR/GROIN_COLOR/OBSERVED_COLOR were an orange, a dark
red and near-black chosen here, and the dark red was the 1984 vintage colour
doing a second job.
```

```text
RANKED ON PERIOD 1 ALONE. Cumulative trapping is M*(16 + 4f) in
1984-2004 but 20*M*f in 2004-2024, so period 2 sees only the PRODUCT
M*f and cannot separate the two parameters. Every bit of information
distinguishing M from f lives in period 1, the one window straddling
the 1996-2003 ramp.

Summing both legs let a window that cannot identify the parameters pull
them anyway, and in the wrong direction: holding period-1 trapping
fixed, driving f from 0.9 to 0 forces M from 50 to 61.3 m/yr -- away
from what the sediment budget can afford, not toward it.

joint_err is still computed and reported, so the change is visible in
the output rather than only in this comment.
```

```text
A value sitting on the edge of its grid is a bound, not an optimum: the
search wanted to keep going and ran out of grid. Reported explicitly so
a railed result cannot be quoted as a fitted parameter.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_period()`**

```text
Loads one sweep's scored results.

Returns:
    A DataFrame, or None if that sweep has not produced a CSV yet.
```

**`period_surface()`**

```text
Reduces one period's rows to a value per (M, f) cell.

`err` is the FILLET-SIZE error where the sweep recorded one, the
differential error otherwise. The fillet saturates, so its level
carries the information about M and its slope -- which is what the
differential measures -- carries almost none.

be1 is profiled out: for each (M, f) the best-scoring be1 is kept, and
which one it was is carried along so the fitted background erosion can be
reported with the fitted groin.

The M = 0 baseline carries no f, so it is broadcast across every f -- with
no groin attached, all f give the identical run, and leaving the row at a
single f would punch a hole in the surface at M = 0.

Returns:
    A DataFrame indexed by (M, fraction) with columns err, be1,
    differential.
```

**`joint_fit()`**

```text
Builds the joint surface for one preset and picks its best cell.

Returns:
    (surface, fit) where surface is a DataFrame over (M, f) and fit is a
    dict describing the winning cell, or (None, reason) if a period's
    sweep is missing.
```

**`plot_constraints()`**

```text
Draws each period's own best-fit valley in (M, f), and where they cross.

This is the figure that shows WHY the joint fit lands where it does: each
period contributes a valley, period 1's running along M*(16 + 4f) = const
and period 2's along M*f = const, and the fitted point is their crossing.
```

**`_pinned_presets()`**

```text
Presets in an existing joint_fit.json that were set by hand, not fitted.

A pinned entry carries `provenance` naming how it got there. This ranking
does not write that key, so its presence is the marker.

Args:
    path: joint_fit.json, which may not exist.

Returns:
    Sorted preset names that are pinned. Empty if the file is absent,
    unreadable, or holds only ranking output -- an unreadable file must not
    be allowed to block a legitimate write.
```

**`_write_fits()`**

```text
Writes the ranking's answer, unless it would clobber a hand pin.

WHY THIS GUARD EXISTS. This ranking scores BOTH periods jointly, and
period 2 records a fillet RELEASE the module cannot produce at any (M, f).
So it rails: on 2026-08-30 it returned edgeBE M = 160 / f = 0.8 and zeroBE
M = 140 / f = 1.0, both at a grid bound. Fitting period 2 is the wrong
thing to attempt -- see hard-structures/groin/GROIN_PLAN.md -- so the file
is pinned by hand to M = 60, f = 0.6.

HAT_run_all.py stage 6 passes whatever this file holds to every groin run
in the matrix. A stage-5 re-run would therefore silently rebuild the whole
matrix on the railed values, and nothing downstream would notice. That was
found on 2026-08-31 with the railed pair sitting in the file.

Refusing rather than warning is deliberate: stage 5 exits 0 and stage 6
then reads the PRESERVED pin, so the pipeline does the right thing
unattended. The ranking is not lost -- it goes to a sidecar.

Args:
    fits: {preset: fit dict} this run computed.
    force: Overwrite a pinned file anyway.
```

</details>

### groin-sweep/HAT_groin_position_figure.py

Observed against modelled shoreline position across the groin, per period.

From the script's original header:

```text
Observed vs modelled SHORELINE POSITION across the groin, per period.

Every other figure in this directory plots a RATE. This one plots position:
where the shoreline started, where it ended, and where the model put it at the
fitted (M, f). A rate figure can hide a run that gets the trend right from the
wrong place; a position figure cannot.

WHERE THE OBSERVATIONS COME FROM
    `HAT-groin-gis-analysis/.../groin_analysis_chainage_all.csv` -- 904k
    shoreline observations carrying `chainage_m`, the cross-shore distance from
    the project's offshore datum line, already mapped to CASCADE domains.
    Chainage is SEAWARD-POSITIVE: verified on 2026-08-23 by differencing the
    1984 and 2004 domain means and correlating against the published CoastSat
    LRR (+0.74). The model's x_s is landward-positive, so it is negated.

    POSITIONS ARE FITTED, NOT BINNED. Two year-bins differenced give a noisy
    endpoint rate -- at D1 that route gave -1.6 m/yr against a published LRR of
    -4.2. Instead an OLS line is fitted to chainage against decimal year across
    the whole period and evaluated at both ends, so the plotted start and end
    are consistent with the LRR that the sweep is scored on. Same estimator,
    same answer.

WHY THE START LINES COINCIDE -- read this before reading the figure
    The model's cross-shore origin is Barrier3D's own, not a real datum, so the
    model cannot be placed on a surveyed axis independently. Following
    `build_shoreline_target` in cascade_pipeline.hindcast, the model is plotted
    as OBSERVED START + MODEL CHANGE. The consequence is unavoidable and worth
    stating plainly: model and observed start at the same place BY
    CONSTRUCTION. Only the separation at the END year carries information. The
    figure annotates this rather than letting a reader mistake a shared start
    for a validated initial condition.

THE ZOOM PANEL, AND WHY IT IS NOT A FAILURE
    The groin field -- four structures, northing 3901373-3901789 -- sits
    entirely inside D6, which spans 3901298-3901798. Its fillet is ~190 m
    wide while a model domain is 500 m, so the model CANNOT resolve it: its
    dipole is 500 m wide by construction. The right panel plots the
    observations at transect resolution against the model's 500 m steps so the
    mismatch is visible as a resolution limit rather than being averaged away.
    A reader who sees the model miss a 190 m notch should know the model was
    never able to draw one.

Usage:
    python HAT_groin_position_figure.py
    python HAT_groin_position_figure.py --preset zeroBE

Writes to output/calibration/groin/figures/:
    position_<preset>.png
```

Notes that were in the code:

```text
The groin field's real footprint, from HAT_groin_shoreline_analysis_v2.py's
metadata. Drawn on the zoom panel so the structure's extent is visible
rather than implied by a single line.
```

```text
CoastSat only: nc_state and wet_dry are different features measured to
different definitions, and mixing them into one fitted line would put
the source's offset into the trend.
```

```text
DOMAIN-COORDINATE shading belongs to the reach panel ONLY. The zoom
panel's x axis is alongshore METRES, so a span drawn at 5.5-6.5
lands six metres from the origin and drags the axis back to zero --
which is exactly what it did before this was split.
```

```text
--- zoom panel: transect resolution ------------------------------
X IS REAL ALONGSHORE DISTANCE, NOT THE DOMAIN ID. Plotting transects
against their integer domain puts every transect in a cell on the
same x, so a domain's whole spread draws as one vertical spike and
the fillet -- the entire point of this panel -- is unreadable.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_chainage()`**

```text
The shoreline-position observations, CoastSat only.

Returns:
    A DataFrame with domain, decimal_year, chainage_m, alongshore_m.

Raises:
    FileNotFoundError: If the GIS analysis output is absent.
```

**`fitted_positions()`**

```text
Start and end shoreline position per domain, by OLS over the period.

The same estimator the sweep is scored on (see `compute_lrr`), so the
endpoints of this line and the LRR are the same statement about the same
data. Differencing two year-bins instead gives an endpoint rate that
disagrees with the published LRR by a factor of two at the noisiest
domains.

Args:
    chainage: Observation frame from `load_chainage`.
    period: 1984 or 2004.
    domains: Iterable of GIS domain ids.

Returns:
    (start, end, n_obs) dicts keyed by domain, in metres seaward-positive.
    Domains with too few observations are absent from all three.
```

**`transect_positions()`**

```text
Per-transect start and end positions, for the zoom panel.

Same OLS treatment as `fitted_positions` but without the domain averaging,
so a fillet narrower than a domain survives.

Returns:
    A DataFrame with alongshore_m, start_m, end_m, one row per transect.
```

**`domain_alongshore_bounds()`**

```text
Alongshore extent of each domain, measured from the observations.

Derived from the transects rather than from the domain polygons because
the zoom panel plots observations on an alongshore axis and the model's
cells must be drawn on the SAME axis. Taking the extent from the data that
is being plotted keeps the two aligned even where a domain's transect
coverage is partial.

Args:
    chainage: Observation frame from `load_chainage`.
    domains: (first, last) inclusive domain range.

Returns:
    {domain: (min_alongshore_m, max_alongshore_m)}.
```

**`model_change()`**

```text
Modelled start->end shoreline change per domain, seaward-positive.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".
    combo: Combination directory name.
    domains: GIS domain ids wanted.

Returns:
    A dict of domain -> change in metres, or None if the run is absent.
```

</details>

### groin-sweep/HAT_groin_sediment_budget_figure.py

The groin module's sediment budget: what it moves, and what it keeps.

From the script's original header:

```text
The groin module's sediment budget: what it moves, and what it keeps.

WHY THIS FIGURE EXISTS
    `groin_diagnostics.csv` has recorded `cumulative_updrift_m` and
    `cumulative_downdrift_m` every year of every groin run since the module was
    written, and nothing has ever plotted them. Two things that matter for
    reading M were therefore visible nowhere:

      1. THE DIPOLE IS VOLUME-NEUTRAL AND THE REAL STRUCTURE IS NOT. The module
         takes exactly as much from the downdrift cell as it gives the updrift
         one. The real Buxton structure accretes updrift with NO measurable
         downdrift deficit -- observed downdrift extent 0 m against the model's
         2,500 m -- because the sand comes from the cape, not from D5. This is
         the module's largest structural assumption and it had no figure.

      2. ALMOST NONE OF WHAT THE DIPOLE INJECTS IS RETAINED. Over the rig's 50
         years at M = 60, f = 0.6 the module applies ~2,400 m of cumulative
         one-sided displacement and holds a fillet of ~69 m. BRIE's alongshore
         diffusion removes the rest. So M is NOT the rate at which sand is
         impounded -- it is the rate needed to SUSTAIN a fillet against
         diffusion, which is a much larger number.

    Panel (c) is the reason this figure is worth having. It reframes the
    affordability comparison that GROIN_PLAN.md and the run reports both make:
    719,000 m3/yr at M = 60 against a 5-7e5 m3/yr littoral drift is a GROSS
    restoring rate set against a NET transport budget, and they are not like
    for like. That does not make the comparison wrong -- it is a deliberate,
    documented diagnostic -- but it does mean "marginally above the drift band"
    should not be read as "impounds more sand than the coast carries."

WHAT IS PLOTTED
    (a) ANNUAL INTERCEPTION. M_eff each year converted to a volume by the
        repo's own `implied_interception_m3_yr` (M * dy * profile height),
        against the shaded 5-7e5 m3/yr littoral drift band. Shows the
        deterioration schedule carrying the module from marginally above the
        band down into it.
    (b) CUMULATIVE VOLUME, updrift and downdrift as exact mirror images. The
        symmetry IS the assumption; the annotation is where it departs from the
        field evidence.
    (c) GROSS AGAINST NET. Cumulative applied displacement against the fillet
        actually realised, same axis, same units.

    Volumes use the repo's conversion and the same profile height as
    HAT_groin_choice_figure.py, so the numbers reconcile with the affordability
    figures quoted elsewhere.

Usage:
    python HAT_groin_sediment_budget_figure.py

Writes output/calibration/groin/figures/sediment_budget.png
```

Notes that were in the code:

```text
The rig pads 11 real domains (D2-D12) with 15 buffer either side, so D5 -> 18
and D6 -> 19. This is the RIG's convention and differs from production's
(D5 -> 19, D6 -> 20) -- see HAT_groin_hindcast_1967_2017.py:76.
```

```text
The rig lives in output/calibration/groin_rig/, not output/raw_runs/ (moved 2026-08-31;
under calibration/ since 2026-09-18).
It is a DIFFERENT GRID -- 41 domains against production's 120 -- and M is
grid-specific, so mixing the two invited quoting a rig number as a
production one. raw_runs is production only, and run_index.csv covers it.
```

```text
A rig run directory does not name its own parameters -- the sweep writes every
cell into one run name. Checked against the diagnostics rather than trusted.
```

```text
The two sides of the structure, in the house semantics: the sheltered
updrift cell is the ACCENT (the thing the module does), the downdrift
cell it takes from is BASE grey. INK/UP_C/DOWN_C/BAD/MUTED/GRID were a
navy, an orange, a blue, a pink and two greys chosen in this file.
```

```text
Sign convention in the CSV: updrift is negative (seaward in the model's
landward-positive frame), downdrift positive. Magnitudes are identical by
construction -- that identity is the point of panel (b).
```

### groin-sweep/HAT_groin_sweep.py

Groin / background-erosion grid search for one period and preset.

From the script's original header:

```text
Groin / background-erosion grid search for ONE period and preset.

Sweeps the groin trapping rate M against the deterioration floor f -- and, in
period 1 under edgeBE, against the GIS 1 background-erosion rate be1 -- scoring
each cell on that period's CoastSat LRR over D1-D12. Every reported number
comes from a real CASCADE run; no surrogate, no interpolation.

This script fits ONE period. M and f are not separable within a single period
(see HAT_groin_sweep_config's JOINT IDENTIFIABILITY note), so the actual
parameter choice is made by `HAT_groin_joint_fit.py`, which intersects the two
periods' surfaces. What this script produces is one period's surface.

WHY M, f AND be1 TOGETHER
    The groin sits at GIS 5.5, four domains from the modelled reach's southern
    boundary. Its downdrift sink does not stay local: at M = 70 it imposes
    roughly -0.5 m/yr at D1 and -0.7 m/yr at D5. The be1 edge source pushes on
    exactly the same domains, so fitting either with the other held fixed
    charges the same erosion to whichever knob happens to be free. f enters
    through the cumulative trapping the schedule actually delivers. All three
    are swept together because none of them is separable from the others.

WHAT IS RANKED, AND WHAT IS ONLY REPORTED
    Ranked:   |differential - observed|, where differential is the modelled
              D6 - D5 rate. It is the only metric that identifies M. Over
              D1-D12 the profile RMSE moves 7% while M moves 4x; the
              differential moves 5x over the same span.
    Reported: RMSE over the D5/D6 pair, RMSE and bias over D1-D12, the
              per-domain modelled rates, and the fillet extent.
    Never fit: the extent. `cascade.groin`'s design argument is that M sets
              amplitude only and the alongshore extent falls out for free, so
              the emergent extent is the module's one independent test. It is
              measured against the paired M = 0 run at the same be1 and
              written to the CSV, and nothing ranks on it.

WHAT THIS SWEEP CANNOT DO, STATED UP FRONT
    Period 1's residual over D1-D12 does not have the shape of a single-domain
    edge source. At the joint least-squares optimum it still runs -1.7 m/yr at
    D6 and +1.7 m/yr at D11: the model under-erodes the Cape Point reach by
    2-3 m/yr, be1's influence decays from 0.121 to 0.013 m/yr per unit across
    the window, and no value of be1 has the reach to flatten that. Expect a
    floor around 1.0 m/yr RMSE. That floor is the finding -- it is what
    `calibBE`'s Cape Point entries were invented to absorb, and holding to
    strict edgeBE leaves it visible instead of hiding it in the interior.

    Period 2's observed D6 - D5 is NEGATIVE (-2.47 m/yr): the updrift domain
    eroded faster than the downdrift one. The source/sink pair cannot produce
    that at any M >= 0, so the period-2 leg reports a BOUND, not an optimum,
    and says so in its summary. Do not read its best cell as a fitted value.

    The model restarts its dipole fresh at the start of each period, from an
    observed surface that already carries 15 years (1984) or 35 years (2004)
    of accumulated fillet. Treating M as one structure-level parameter across
    both windows -- which the joint fit does -- is a modelling decision, not
    something these runs establish.

DESIGN
    One subprocess per combination. An earlier in-process sweep died partway
    through with a Windows access violation (0xC0000005) from state
    accumulating across many Cascade constructions; process isolation makes
    the OS reclaim everything between combinations.

    Resumable. Combinations with a scored result are skipped; rows recorded
    as failed are retried automatically on the next invocation. Each result
    is appended as a JSON line the moment it lands, so an interrupted sweep
    costs at most one run.

    Parallel. NUM_CORES is 1 in the worker (>1 has crashed on this
    configuration) and CASCADE's internal joblib.Parallel therefore does not
    fan out, so concurrent workers do not oversubscribe. The pool is capped if
    the requested width will not fit in available RAM.

    Self-validating. Before reporting anything, the worker's duplicated copy
    of build_cascade / run_cascade_simulation is run against the period's
    published matrix run and the two rate curves differenced. The validation
    cell lives in a per-period `_validation/` directory and is shared by both
    presets, because what it checks is the CODE, not the preset.

Usage:
    python HAT_groin_sweep.py --period 1984 --preset edgeBE [--workers N]
                              [--dry-run] [--skip-validation]

Writes to output/calibration/groin/<start>_<end>_<preset>/:
    sweep_results.jsonl      one JSON line per combination, written as it lands
    sweep_results.csv        the same rows as a table, plus extent, at the end
    <combo>/                 shoreline matrix + rate curve per combination
Nothing is written to output/raw_runs/ or run_index.csv: a sweep combination
is not a run of the scenario matrix and must not be filed as one.
```

Notes that were in the code:

```text
parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
The guard below is what makes a future move fail here, loudly, instead of
resolving to scripts/scripts and surfacing as a missing data file several
imports deeper.
```

```text
A 20-year run is ~1.5 min; 10 minutes is a wide margin for a loaded machine
and stops one wedged worker from holding the pool open indefinitely.
```

```text
A sweep with only failures has no differential_err column at all, and
every caller filters on it.
```

```text
One thread per worker. numpy/BLAS defaults to one thread per core, so N
concurrent workers each spawn N threads and the pool spends its time
context-switching instead of running: measured, four workers took 4-5x
longer per cell than one, for almost no net throughput. CASCADE's own
NUM_CORES is already 1, so nothing here wants the extra threads.
```

```text
No result line: an access violation exits non-zero with nothing on
stdout, so the stderr tail is the only diagnostic there is.
```

```text
1.8 GB/worker, measured from a 120-domain 20-year run, plus 2 GB left
for the OS and this process.
```

```text
The worker writes nothing to stdout when it dies, so
its stderr tail is the only record of WHY. Printed on
the first failure of a sweep rather than only into the
results file: a sweep where every cell fails otherwise
produces a log of 215 identical "exit 1" lines and no
diagnosis, which is how four sweeps came to be broken
without anyone being able to say what broke them.
```

```text
The metadata stores the deterioration as prose, so the floor has to be
parsed back out. An unreadable string is reported rather than defaulted:
guessing the floor would validate the worker against a groin schedule
the reference did not run.
```

```text
A retried combination appears twice in the JSONL. The last line wins:
it is the one whose files are on disk.
```

```text
RANKED ON FILLET SIZE. differential_err is retained and reported,
but it scores the fillet's SLOPE, which is near zero once the
fillet saturates and therefore nearly uninformative about M. Falls
back to the differential only if no fillet could be measured --
which means the M = 0 baselines are missing, and the sweep should
be re-run rather than ranked on the weaker metric silently.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load_results()`**

```text
Loads every recorded result, or an empty frame if none exist.

Returns:
    A DataFrame of previous results with guaranteed `combo` and
    `differential_err` columns. Failed rows carry NaN in
    `differential_err` and a message in `error`.
```

**`append_result()`**

```text
Records one result immediately, as a JSON line.

JSONL rather than appending to the CSV: a failure record has far fewer
keys than a success record, and `to_csv(mode="a")` writes values
positionally under the existing header -- so one failure after a run of
successes silently shifts every later row's columns. JSON lines carry
their own keys, so a mixed-schema log stays readable and resumable.
```

**`run_worker()`**

```text
Runs one combination in its own subprocess.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".
    combo: An (M, be1, fraction) tuple.
    out_root: Directory the combination's subdirectory goes under.

Returns:
    A result dict on success, or a failure record carrying `error` and
    the worker's last stderr lines.
```

**`safe_worker_count()`**

```text
Caps the pool width at what available RAM will hold.

Each worker builds 120 Barrier3D domains and holds the run's full state.
Oversubscribing memory does not fail cleanly -- it swaps, and a swapping
sweep is slower than a serial one.

Returns:
    The width to actually use. Falls back to `requested` if psutil is not
    installed, since a missing optional dependency should not block a
    sweep the user explicitly sized.
```

**`read_reference_config()`**

```text
Reads the reference matrix run's own M, f and be1 from its metadata.

Not pinned in the config module on purpose: the reference run is re-run
with new values as the fit is refined, and a pinned pair would report
drift the first time it changed.

Returns:
    ((M, be1, fraction), None) on success, or (None, message) if the
    reference is missing, mid-run, or not comparable to what the worker
    builds.
```

**`validate_against_matrix_run()`**

```text
Checks the worker reproduces the period's reference matrix run exactly.

The worker holds a third copy of `build_cascade` and
`run_cascade_simulation` -- the notebook has one, the hindcast runner has
another. Copies drift, and drift here is silent: the sweep would keep
producing plausible numbers from a model that no longer matches the
hindcast. Running the reference's own configuration and differencing the
rate curves catches that numerically.

The validation cell lives in a per-period `_validation/` directory shared
by both presets, and always runs under edgeBE regardless of which preset
is being swept: what it checks is the duplicated CODE, not the preset.

Returns:
    (ok, message). ok is False if the sweep must not proceed.
```

**`attach_fillets()`**

```text
Adds the fillet size and its error to every M > 0 row.

THIS IS THE RANKING METRIC, and it is computed here rather than in the
worker for the same reason the extent is: it needs the paired M = 0 run,
and pairing is a property of the grid, not of one combination. Each cell
pairs with the baseline at the SAME be1 -- pairing across be1 would report
the background-erosion difference as fillet.

Args:
    frame: Scored results, one row per combination.
    out_root: Directory holding the per-combination subdirectories.
    period: 1984 or 2004, selecting the observed fillet to score against.

Returns:
    A copy of `frame` with `fillet_m` and `fillet_err` columns. Both are
    NaN on the M = 0 baseline rows, which have no fillet by definition.
```

**`attach_extents()`**

```text
Adds the emergent fillet extent to every M > 0 row.

Each groin cell is paired with the M = 0 cell at the SAME be1. Pairing
across be1 would report the background-erosion difference as fillet.
Nothing here feeds the ranking -- the extent is the pre-registered
independent check from `cascade.groin`.
```

</details>

### groin-sweep/HAT_groin_sweep_comparison.py

Cross-reference figures spanning every groin sweep at once.

From the script's original header:

```text
Cross-reference figures spanning every groin sweep at once.

`HAT_groin_sweep_figures.py` draws one sweep at a time, into that sweep's own
directory. That is the right place for a diagnostic and the wrong place for a
comparison: answering "does edgeBE put the optimum where zeroBE does?" or
"do the two periods agree about f?" currently means opening four folders and
holding four colour scales in your head.

This file draws the four sweeps together.

WHAT THE COMPARISON IS FOR
    The sweep grid is the same in every cell -- same M values, same f values,
    same observed target per period -- so the four surfaces are directly
    comparable and the interesting content is where they DIFFER:

    across presets (zeroBE vs edgeBE)
        Same period, same groin, different background erosion. If the optimum
        moves, the fitted groin is absorbing background erosion rather than
        describing the structure.

    across periods (1984-2004 vs 2004-2024)
        Same structure, different window. Period 1 straddles the 1996 repair
        and the 2003 storm and is the only window that can separate M from f;
        period 2 sits entirely past the ramp and sees only the product M*f.
        Disagreement here is the scientific result, not a defect.

WHY THE PANELS DO NOT SHARE A COLOUR SCALE
    Period 1's fillet error spans roughly 0-40 m and period 2's roughly 43-95
    m, because period 2's observed fillet is NEGATIVE (-43.2 m: the fillet
    relaxed) and no M >= 0 can build a negative fillet, so every cell carries
    at least that much error. Forcing one scale would flatten period 1 --
    where the actual optimum lives -- into a single colour. Each panel is
    scaled to its own range and the numbers are given on the panel, so the
    comparison is read from the annotations rather than from the hue.

TIES ARE DRAWN, NOT RESOLVED
    In both period-2 sweeps the whole f = 0 row scores identically: a fully
    deteriorated groin traps nothing, so M has no effect and seven cells tie
    to within 5e-5 m. Marking one of them as "best" would report a fitted M
    that is really just whichever cell sorted first. Tied sets are drawn as
    open circles and labelled as unconstrained.

Usage:
    python HAT_groin_sweep_comparison.py
    python HAT_groin_sweep_comparison.py --top-n 3

Writes to output/calibration/groin/figures/:
    comparison_surfaces.png   the four M-f error surfaces side by side
    comparison_optima.png     every sweep's optimum in one (M, f) plane
    comparison_profiles.png   each sweep's best LRR curve against CoastSat

Sweeps that have not run yet are drawn as labelled placeholders rather than
skipped, so a missing panel reads as "not swept" instead of silently
shrinking the figure.
```

Notes that were in the code:

```text
Imported, never copied. The single-sweep figures and these comparisons must
reduce a sweep the SAME way -- same be1 profiling, same tie rule, same
ranking metric -- or the two figure sets would disagree about which cell won.
```

```text
One colour per sweep, stable across all three figures so a reader who learns
"orange is 1984 edgeBE" on one figure keeps it on the next.
Colour is the PERIOD, marker is the preset (2026-09-11). Four invented hues
were in use here -- a blue, an orange, a green and a dark red -- which spent
two colours on a distinction the marker already carries, and put the vintage
red on one arbitrary sweep. The house pair means the same thing on every
figure in the project: red is the earlier period, blue the later.
```

```text
The valley floor shows the SHAPE of each sweep's constraint, which is
what makes two sweeps comparable even when their optima coincide: a
flat floor means the axis is unconstrained, a steep one means it bites.
```

<details><summary>Function notes (the original docstrings)</summary>

**`collect()`**

```text
Loads every sweep that has results.

Returns:
    {(period, preset): surface}, where surface is one row per (M, f) with
    be1 profiled out. Sweeps with no results are absent from the dict.
```

</details>

### groin-sweep/HAT_groin_sweep_config.py

Shared constants for the groin / background-erosion sweep: target, fit window and cell naming.

From the script's original header:

```text
Shared constants for the groin / background-erosion sweep, both periods.

The orchestrator and the worker must agree on the target, the fit window and
the combination naming, or the orchestrator will rank rows the worker scored
against something else. Both import from here so there is one copy of each.

WHAT CHANGED FROM THE 1984-2004-ONLY VERSION
    Three things, all forced by fitting M and f jointly across both periods
    rather than fitting M alone in one period:

    1. Everything period-specific is now keyed by start year. There is no
       module-level START_YEAR; a caller states which period it means.
    2. The deterioration fraction f is a SWEPT AXIS, not a two-point bracket.
       The old version's argument for bracketing it -- that over 1984-2004 the
       cumulative trapping is M*(16 + 4f), so f barely separates from M -- is
       still true and is why f cannot be fit from period 1 alone. It is fit
       from the two periods together (see JOINT IDENTIFIABILITY below).
    3. The observed targets are COMPUTED from the CoastSat transect files
       rather than typed in. The 1984-2004 values are asserted against the
       numbers the old version carried, so the change is provably a
       refactor and not a new target.

JOINT IDENTIFIABILITY -- why two periods are needed for two knobs
    Cumulative trapping per unit M, measured from `cascade.groin.GroinCallback`
    with the documented 1969 install / 1996 onset / 7-year ramp schedule:

        1984-2004    16 + 4f      (f moves this only 16.0 -> 20.0)
        2004-2024    20f          (the run sits entirely past the 2003 ramp,
                                   so only the product M*f is identifiable)

    Period 1 alone leaves f nearly free: sliding f across its whole range needs
    only a 25% change in M to hold cumulative trapping fixed. Period 2 alone
    cannot separate M from f at all. Together they intersect:

        f = 4*B / (5*A - B),   A = period-1 constraint, B = period-2

    THIS BLOCK USED TO SAY that period 2's negative differential pinned B near
    zero and drove f to 0 -- "the structure stopped trapping after the 2003
    storm" -- and called that a RESULT rather than a fitting artifact. That was
    wrong, and it was wrong because of the AGGREGATION, not the physics.
    Corrected 2026-08-23 against the CoastSat transects:

        the fillet is ~190 m wide; one model domain is 500 m, so the whole
        dipole fits inside a single cell. Domain averaging blends the fillet's
        peak with its own far field, and in 2004-2024 that does not merely
        weaken the signal -- it INVERTS it. Transect-scale the dipole is
        +0.48 m/yr (updrift accreting, groin-like); the domain-mean D6-D5
        difference is -2.47 m/yr (anti-groin).

    Scored against the domain-mean target the sweep drove f to 0, because the
    honest way to match a negative fillet is to trap nothing. Scored against
    the transects, the SAME runs show the groin holding ~44% of its period-1
    strength, which is what the imagery shows: deteriorated, still working.

    So `observed_fillet_m` is now built from transects (see below), and the
    period-2 leg constrains the PRODUCT M*f rather than reporting a bound at
    zero. `differential`, being a domain-mean quantity, carries the same
    aggregation defect and is retained for reporting only -- never ranked on.
```

Notes that were in the code:

```text
parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
The guard below is what makes a future move fail here, loudly, instead of
resolving to scripts/scripts and surfacing as a missing data file several
imports deeper.
```

```text
Every groin-sweep product lives under here: the sweeps, joint_fit.json (an
INPUT that HAT_run_all.py stage 6 reads), the figures and SELECTED_. Moved
from output/calibration/groin/ on 2026-09-18; build paths from this, never by hand.
```

```text
THE CANONICAL CHAIN IS 1996 -> 2010 -> 2024 (Hannah, 2026-09-17), so the
matrix and the sweep iterate these two starts. Both are in
HATTERAS_PERIODS and END_YEAR derives from it, so the windows cannot
drift from the site config.

TWO THINGS THIS DOES NOT CARRY OVER, and neither is silent:
* calibBE is solved for 1984 and 2004 ONLY. PRESETS below is
(edgeBE, zeroBE), both of which are solved for all four starts, so
the matrix is unaffected; a calibBE run on 1996 or 2010 raises the
explicit "not solved for that period" error from be_rates().
* the pinned groin fit (M = 60, f = 0.6 in output/calibration/groin/
joint_fit.json) was fitted on PERIOD 1 = 1984-2004 and records
fit_period 1984. Stage 5 run under this pair intersects the 1996 and
2010 surfaces instead and would produce a different pair -- and that
file is pinned, with its own warning that re-running stage 5
overwrites it. The existing fit stays valid for the window it was
fitted on; it is simply not re-derived by this matrix.
```

```text
FOUR STRUCTURES, ONE DIPOLE -- BY DESIGN, NOT BY OVERSIGHT.
The Buxton groin field is four groins spanning northing 3901373-3901789
(HAT_groin_shoreline_analysis_v2.py metadata, n_groin_features = 4). D6 spans
3901298-3901798, so the ENTIRE field falls inside a single model domain.

The module therefore represents the field's CUMULATIVE effect as one
source/sink pair rather than four. There is nothing to gain from splitting
them: four dipoles inside one 500 m cell would sum to exactly the one dipole
the cell can express, and `GroinCallback` requires adjacent domains anyway,
so a per-structure representation is not expressible on this grid.

WHAT THIS MAKES M. M is the trapping rate of the FIELD AS A WHOLE, not of a
groin. It cannot be divided by four to get a per-structure rate, and it
cannot be compared against a published single-groin trapping rate. Combined
with the resolution note in `observed_fillet_m` -- the field's fillet is
~190 m wide inside a 500 m cell -- M is an effective, grid-specific,
field-aggregate quantity. The defensible structure-level statement is the
RATIO of trapping between periods, which is what f carries.
```

```text
The SAME schedule is used in both periods, deliberately. M and f are one
structure-level pair, so period 2 must not be given a special "no
deterioration, M already reduced" configuration -- that would hard-code the
collapse to M*f instead of letting it emerge, and would make the period-2
fit a different parameter from the period-1 one.
```

```text
Fraction of the peak effect that defines the fillet's edge. Matches the
hindcast runner so a sweep extent and a matrix extent mean the same thing.
```

```text
M = 0 is the paired-baseline column, not a weak groin: measure_groin_extent
needs a no-groin run at the SAME be1, and a baseline at a different be1
would report the background-erosion difference as fillet. f is meaningless
when M = 0, so the grid builder collapses those cells to one per be1.

The ceiling is 110, not the 80 the M-only sweep used. Freeing f forces it:
period-1 cumulative trapping is M*(16 + 4f), so reaching at f = 0 what
M = 80 delivers at f = 0.9 needs M ~= 98. An 80 ceiling would clip the
ridge exactly where the joint solution is heading.
EXTENDED ABOVE 110 ON 2026-08-24, to document the ridge rather than argue it.
The 1984 edgeBE fit landed on M = 110, f = 0.4 -- the grid maximum -- with a
fillet error of 0.0026 m. That is NOT "the grid ran out before the optimum":
fillet is monotonically increasing in M at every f, so M > 110 at f = 0.4
OVERSHOOTS the 22.11 m target. The at_grid_bound flag is mechanical.

What the grid edge hides is a RIDGE. Reading the be1 = -34 surface, the target
is matched at roughly (M 46, f 1.0), (M 57, f 0.8), (M 75, f 0.6) and
(M 110, f 0.4) -- every one an equally good endpoint fit. Extending to 160
adds the f = 0.2 arm of that ridge (which needs M ~ 200 to reach the target),
so the figure shows a ridge running off the grid instead of a peak sitting on
its edge. It resolves nothing on its own; the trajectory metric is what
separates these cells, because they reach the same fillet along very
different paths.
```

```text
be1 is swept in period 1 only. The 2004 edgeBE values are taken as given --
there is no prior fit behind a 2004 bracket, and inventing one would add an
axis that resolves nothing. See the orchestrator docstring.

THE SHALLOW END EXISTS TO BRACKET, NOT TO SEARCH. The site config's own
edgeBE value at GIS 1 for 1984 is -24.0, which sat between the old
bracket's first two rungs and within one step of its shallow edge. A fit
landing on -22.0 could then mean either 'the optimum is -22' or 'the
optimum is shallower than -22 and the grid ran out', and joint_fit's
at_grid_bound check cannot tell those apart -- it can only say the value
railed. -16 and -10 were added on 2026-08-22 so the configured value is
bracketed on BOTH sides and a shallow optimum is a result rather than a
bound. They cost 2 x 43 = 86 extra cells in the 1984 edgeBE sweep.
-42.6 IS THE CALIBRATED VALUE ITSELF, added 2026-08-30. The 2026-08-29
sweep bracketed it (-46, -40) but never evaluated it, and the two
brackets disagreed badly about M (95 vs 160), so interpolating across
them would have invented a number the grid never measured. Fixing be1
at the calibrated value turns M into a one-parameter question.
```

```text
Raw per-domain transect means over D1-D12, not the LOWESS-smoothed table:
LowessConfig's skip_southern_domains is 10, so D1-D10 are raw in
COASTSAT_TARGET anyway, and taking D11-D12 raw as well makes the whole
window one construction instead of two. Sign: + is seaward/accreting.
```

```text
The values the M-only sweep carried as a literal table. Kept ONLY as an
assertion target: if computing them from the transect file no longer
reproduces these, the target moved and every published period-1 number
needs re-checking, so that must fail loudly rather than quietly re-fit.
```

```text
Provable refactor: the computed 1984 target must reproduce the table the
M-only sweep was scored against, to the 2 dp it was written at.
Loaded on its own if the sweep no longer iterates 1984: the guard is about
the PUBLISHED period-1 numbers, not about which periods this run happens to
sweep, so it must not lapse when the canonical chain moves (2026-09-17).
```

```text
The ranking metric's target: observed updrift minus downdrift. Derived from
OBSERVED_LRR rather than retyped, so the two cannot disagree.
```

```text
WHAT THE DIFFERENTIAL ACTUALLY MEASURES. LRR[D6] - LRR[D5] is algebraically
the OLS slope of the fillet, x_s[D5] - x_s[D6] -- its TREND across the
window, not its size. Three consequences, all counter-intuitive, all
measured on the 2026-08-22 1984 zeroBE sweep:

- the fillet SATURATES (~18 m at M = 40, reached by 1990, where BRIE's
diffusion balances the trapping), so a healthy CONSTANT groin scores
near ZERO;
- a DEGRADING groin scores NEGATIVE, because its fillet is relaxing;
- sign therefore reports building-vs-failing, not strong-vs-weak, and two
very different groins can score identically.

THIS BLOCK USED TO CLAIM a negative differential was unreachable at M >= 0,
on the grounds that the callback adds -M updrift and +M downdrift. That
holds for the fillet's SIZE, not its slope. Falsified twice:

1. cell M40_beNA_f0.00 scores -0.713, MORE negative than the no-groin
baseline's -0.229 -- a groin at M = 40 driving the differential below
no groin at all;
2. the period-2 seed run inherits a 25.3 m fillet at t = 0 (the real
Buxton fillet is in the 2004 initial shoreline) and relaxes to 21.9 m,
scoring -0.168 m/yr at M = 50, f = 0.9. Negative, groin attached.

So the observed 2004-2024 differential of -2.47 m/yr IS reachable, and is
what a fillet relaxing after the 2003 storm damage looks like. Period 2
failing to match it is informative rather than guaranteed, which is why it
is now used out-of-sample rather than summed into the objective.

The flag is kept, renamed to say what it really tests: whether the OBSERVED
differential carries the sign of a groin still BUILDING its fillet across
the window. Useful to report; not a statement about what the model can do.
```

```text
Old name kept so nothing importing it breaks. It no longer means 'the model
cannot reach this'.
```

```text
The worker duplicates build_cascade and run_cascade_simulation from the
hindcast runner, which is itself a mirror of the notebook -- three copies
that can drift. The guard runs the reference matrix run's own configuration
through the worker and differences the two rate curves, so drift is caught
numerically rather than by hashing source text.

The reference's PARAMETERS ARE NOT PINNED HERE. They are read from that
run's metadata JSON at sweep time, because the matrix run gets re-run with
new values as the fit is refined -- and a pinned M or be1 here would report
"drift" the first time it did, when nothing had drifted.
SET FROM THE MODEL'S MEASURED NOISE FLOOR, NOT FROM ZERO.

This was 1e-6, which assumes the model is deterministic. It is not. Running
the SAME worker combination (2004 edgeBE, M = 50, f = 0.9, be1 = 50.3) twice
on 2026-08-23 gave two rate curves differing by 7.7e-4 m/yr -- the same
magnitude, and at the same domain, as the "drift" the guard was reporting
against the published matrix run. Reference-vs-run2 (1.6e-4) came out
SMALLER than run1-vs-run2, which is run-to-run scatter, not a code offset.

WHERE THE SCATTER COMES FROM. The difference is float noise (median 1.8e-8
across the reach) everywhere except a smooth bump centred on D16, which is
the first domain OUTSIDE the beach-dune manager's footprint -- the 2022
Buxton fill covers D6-D15, so D15/D16 is the sharpest discontinuity in the
model. A threshold decision there turns last-bit differences into something
a diffusion solver then spreads over ~10 domains.

So a 1e-6 guard can never pass however well-synced the code is, and every
sweep aborts on physics rather than on drift. 5e-3 sits ~6x above the
measured floor and still catches real drift by orders of magnitude: a
genuine divergence in build_cascade or the update loop moves rates by
tenths of a m/yr, not thousandths.

7.4e-4 m/yr inside the fit window is ~1.5 cm of shoreline over 20 years,
far below the resolution of a fit whose observed target cannot even
reproduce the SHAPE of the groin pair (see observed_fillet_m).
```

```text
Scenario tokens that NEGATE a module. The validation reference is the full
management run, so a name carrying any of these is a different scenario and
must not match -- `..._road_bdm_nonourish_groin` (the full_no_fill arm) sat
beside `..._road_bdm_nourish_groin` in 2004-2024 and both matched, so
validation_run_dir returned "2 candidates" and period-2 sweep validation could
not run at all.

MEMBERSHIP, NOT PREFIX. Testing `token.startswith("no")` looks equivalent and
is wrong in the one place it matters: `nourish` itself begins with "no", so a
prefix test rejects the very run being looked for. The tokens are enumerated
for that reason. `nogroin` is in the set, which subsumes the endswith check
this replaced.
```

```text
hatteras_ms is not on the path for this package the way SCRIPTS_DIR is;
added here rather than at import time so nothing else in this module
gains a dependency on the runner's config.
```

```text
Runs are filed <wave scope>/<period>/<preset>/. The stem below already
pins the preset to edgeBE, so the directory pins it to the same thing.

THE WAVE SCOPE IS NOT OPTIONAL. sweep_output_dir and joint_fit_paths are
already scoped by HAT_SWEEP_HS; this lookup was not, so a sweep at
Hs = 3.0 validated its duplicated model code against the Hs = 2.5 matrix
run -- the one arm it must never be compared with. That failure is
SILENT: the reference exists, it is readable, and it is the wrong island.

The component is taken from wave_climate_token, the same function the
runner derives HS_TAG from, rather than being spelled here. A second copy
of the "3.0 -> waveHs3" rule is how the two would drift apart.
```

```text
The offset token sits after the preset in a run name
(HAT_1996_2010_edgeBE_offsetmetres_road_bdm_...), and without it an
asrun sweep would validate against a metres run or the reverse.
```

```text
A sweep run at a non-default Hs gets its OWN directory. The worker reads
HAT_SWEEP_HS and would otherwise write its results over the 2.5 m sweep
the fitted M = 60 came from -- the fit would be silently replaced by one
made under different forcing, with nothing on disk recording that it had
happened. Same reasoning as the run-name tokens on the hindcast side.
```

```text
The offset mode likewise (2026-09-24): a metres sweep is a different
island from the /10 sweeps already on disk, not a re-run of them.
```

```text
WHY SIZE AND NOT SLOPE. `differential` (LRR[D6] - LRR[D5]) is algebraically
the fillet's OLS SLOPE. The fillet saturates -- ~18 m at M = 40, by 1990 --
so its slope over a 20-year window is near zero and carries almost no
information about M, while its LEVEL carries all of it. Scored on slope the
1984 zeroBE pilot railed at M = 110, f = 1.0, both grid edges, at 115% of
the reach sediment budget and still 35% short. Rescored on size, the SAME 43
runs put the optimum at M = 50, f = 0.6 -- interior on both axes, 0.1 m
error, and an f consistent with the 1996 repair and 2003 storm damage.

THE TOPOGRAPHY BUG, AND WHY THE VALUE SURVIVED IT. Until 2026-08-30 the
worker resolved its topography without naming a product, so DEFAULT_PRODUCT
("2004-start") answered and every period-1 cell was built on the 2004
barrier. Fixed; the drift guard now reproduces its reference to 0 m/yr.
Re-scored on the corrected island at the calibrated be1 = -42.6, the D4-D8
demeaned profile puts the best cell at 11.44 m and M = 60 / f = 0.6 at
11.58 m, against 15.20 m with no groin -- so the production pair is within
0.14 m of optimal and the bug moved it less than the ridge it sits on.

DO NOT RE-FIT ON THE FILLET. An attempt on 2026-08-30 returned M = 95 with f
railed at the grid bound and the best M swinging 95 -> 160 between adjacent
be1 values. No admissible M can match the fillet on this grid
(HAT_groin_timeseries_check.py:29); fitting it anyway produces a rail, not an
optimum. D4-D8 DEMEANED is the target that works, and only the PRODUCT M*f
is identified by it -- not M and f separately. See the groin note in
scripts/site_layer/hatteras_site_config.py.
```

```text
THE GROIN IS A SUB-GRID FEATURE. Measured from the CoastSat transects on
2026-08-23, the fillet's observable width is ~190 m while one model domain
is 500 m: the whole dipole fits inside a single cell.

That is why this target is built from TRANSECTS and not from OBSERVED_LRR.
The domain means average the fillet's peak against its own far field and
what survives is the regional gradient, not the groin. In 2004-2024 that
does not merely weaken the signal, it INVERTS it: the transect-scale dipole
is +0.48 m/yr (groin-like, updrift accreting) while the domain-mean
D6-D5 difference is -2.47 m/yr (anti-groin). Scored against the domain-mean
version the sweep drove f to 0 -- "the groin stopped trapping" -- because
the honest way to match a negative fillet is to trap nothing. Scored against
the transects the same runs show the groin retaining ~44% of its period-1
strength, which is what the imagery shows.
```

```text
A zero peak means the two runs are identical, so the threshold is
also zero, every "abs(value) < 0" test is False, and the walk would
report the whole island as fillet.
```

<details><summary>Function notes (the original docstrings)</summary>

**`be1_values()`**

```text
Background-erosion values to sweep at GIS 1 for one period/preset.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".

Returns:
    A list of be1 values in m/yr, or [None] when the preset has no free
    background-erosion knob. None means "use the preset as it stands",
    which the worker turns into an all-zero field for zeroBE and into the
    site-config table values for a period whose be1 is not swept.
```

**`be_gis90()`**

```text
The fixed north-end rate for one period, from the site config.

Read from HATTERAS_BE_RATES_EDGE rather than pinned here. The old version
pinned 15.0, the table has since been re-solved to 10.0 for 1984, and the
orchestrator's drift guard aborts on exactly that mismatch -- a pinned
copy is a stale copy waiting to happen.

Args:
    period: 1984 or 2004.

Returns:
    The rate at GIS 90 in m/yr.
```

**`_load_observed_lrr()`**

```text
Computes the per-domain observed LRR over the fit window.

Args:
    period: 1984 or 2004.

Returns:
    A dict of GIS domain ID to mean LRR in m/yr.

Raises:
    FileNotFoundError: If the period's transect file is absent.
    ValueError: If any domain in the fit window has no transects.
```

**`_wave_scope()`**

```text
The raw_runs path component for this sweep's wave climate.

"" at the calibration climate, so every pre-existing run resolves exactly
where it always did, and e.g. "waveHs3" otherwise. Mirrors HS_TAG in
section 3 of the runner by calling the same token function, so the two
cannot disagree about where a run was filed.

Returns:
    A string, empty at the calibration wave climate.
```

**`sweep_offset_mode()`**

```text
The island-offset mode a sweep runs under.

HAT_SWEEP_OFFSET_MODE, else the runner's own default (field_default), so
a sweep and the matrix run it validates against build the offset the same
way. Every sweep before 2026-09-24 ran at "asrun" (offset / 10): the worker
had its own loader that divided by ten and no mode at all.

Returns:
    One of cascade_pipeline.hindcast.ISLAND_OFFSET_MODES.
```

**`_offset_token()`**

```text
"offset<mode>" for any mode but asrun, else "": the runner's run-name
rule, so a sweep folder, its joint fit and the matrix run it validates
against are all named the same way. asrun keeps no token so the sweeps
made before 2026-09-24 keep their names and are never written over.
```

**`validation_run_dir()`**

```text
Directory of the matrix run this period's sweep validates against.

RESOLVED BY LOOKUP, NOT CONSTRUCTED. The runner's RUN_NAME is derived in
its section 7.5 from what sections 5-6 actually built, and it carries a
`nourish` token only when the period has fill scheduled -- so the
reference is `..._edgeBE_road_bdm_groin` in 1984-2004 but
`..._edgeBE_road_bdm_nourish_groin` in 2004-2024. Building the name here
means reimplementing that derivation, and getting it wrong is silent: the
sweep reports "reference missing" for a run that is sitting on disk under
a name one token different.

Args:
    period: 1984 or 2004.

Returns:
    (path, None) for exactly one match, else (None, message).
```

**`combo_dir_name()`**

```text
Directory / row name for one combination.

Sweep combinations must not use the matrix's RUN_NAME scheme: that name
carries a preset token ("edgeBE") but no be1 or f value, so every
combination would resolve to the same directory and silently overwrite
the last.

M = 0 collapses the f axis -- with no groin attached there is nothing to
deteriorate, so every f would produce a byte-identical run. Those cells
are named without an f token so the grid builder cannot enqueue six
duplicates of the same baseline.

Args:
    M: Groin trapping rate, m/yr.
    be1: Background erosion at GIS 1 in m/yr, or None when the period /
        preset has no swept be1.
    fraction: Groin deterioration floor.

Returns:
    A filesystem-safe name, e.g. "M60_be-34_f0.40" or "M0_beNA".
```

**`sweep_output_dir()`**

```text
Where one period/preset sweep writes its results.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".

Returns:
    The output directory Path.
```

**`joint_fit_paths()`**

```text
Where the joint fit writes, scoped by wave climate.

The default file holds HAND-PINNED values (M = 60, f = 0.6) whose note says
re-running the fit will overwrite them. A fit made at a different Hs must
therefore not land there: it is not a better answer to the same question,
it is the answer to a different one, and the pinned record of how the
production pair was chosen would be gone.

Returns:
    (json_path, csv_path).
```

**`build_grid()`**

```text
Every combination for one period/preset sweep.

The M = 0 baseline is emitted once per be1 rather than once per (be1, f):
with no groin attached f cannot apply, so the extra cells would be exact
duplicates and would inflate the run count by five per be1 for nothing.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".

Returns:
    A list of (M, be1, fraction) tuples.
```

**`measure_fillet()`**

```text
Fillet size at the end of a run, relative to a paired baseline.

The fillet is `x_s[downdrift] - x_s[updrift]`: positive when the updrift
domain sits seaward of the downdrift one, which is what a groin builds.
Differenced against the M = 0 run at the SAME be1, so what is reported is
the groin's contribution and not the coast's natural shape -- D5/D6 sit on
a curving shoreline at Cape Point and carry an offset with no groin at all.

Args:
    shoreline_m: [time x padded domain] array for the groin run.
    baseline_m: The same for its paired M = 0 run.
    geometry: DomainGeometry describing the padded array.
    updrift_gis, downdrift_gis: The groin pair; module defaults if None.

Returns:
    Fillet size in metres at the final year, groin minus baseline.
```

**`observed_fillet_m()`**

```text
Fillet change over one period, from the fixed 1967 wet/dry datum.

THIRD TARGET IN TWO DAYS, and the reasons the first two failed are the
reason this one is built from a different dataset entirely:

  1. Domain-mean LRR anomaly (-43.2 m for 2004-2024). The fillet is ~190 m
     wide inside a 500 m cell, so domain averaging blends its peak with its
     own far field.
  2. Transect LRR peak-to-trough (+22.1 / +9.7 m). `max(one side) -
     min(other side)` is biased POSITIVE: simulated at the observed scatter
     the bias alone is +27 m and +42 m, exceeding the signal it claimed to
     measure. It inverted period 2's sign.

`Change_from_wetdry_1967_D2_D12.csv` avoids both. It is a purpose-built
table of shoreline change against a FIXED 1967 datum, per domain, at 24
dates from 1967 to 2023 -- so the fillet at D5/D6 is differenced directly,
with no smoothing window to blend it and no order statistic to bias it.

A TREND, NOT TWO ENDPOINTS. The year-to-year scatter is 11-14 m, so
differencing the two end years throws away eight intermediate observations
and inherits the noise of both. An OLS fit across the period gives
+3.46 +/- 0.70 m/yr for 1984-2004 and -3.85 +/- 0.76 m/yr for 2004-2024 --
4.9 and 5.1 sigma. Both periods carry a real, well-determined signal, which
neither earlier target could show.

SIGN matches `measure_fillet`: the table is landward-positive (+ = erosion),
so D5 minus D6 is positive when the downdrift domain has retreated further
than the updrift one -- which is what a groin builds.

Args:
    period: 1984 or 2004.
    trend_order: Polynomial order for the fit across the period.

Returns:
    Fillet CHANGE in metres over the period.

Raises:
    FileNotFoundError: If the wet/dry change table is absent.
    ValueError: If a period has too few dated observations to fit.
```

**`measure_groin_extent()`**

```text
Alongshore extent of the groin's effect, from a paired baseline run.

Copied from the hindcast runner's section 12.3. It lives here rather than
in the worker because the orchestrator needs it too -- pairing a groin run
with its baseline is a property of the grid, not of one combination -- and
importing the worker would pull CASCADE into the orchestrator process for
what is a few lines of numpy.

The baseline MUST be the M = 0 run at the same be1. Pairing across
different be1 would report the background-erosion difference as fillet.

Args:
    shoreline_m: [time x domain] matrix from the groin run.
    baseline_m: [time x domain] matrix from the paired M = 0 run.
    geometry: DomainGeometry describing the padded array.
    updrift_gis, downdrift_gis: The groin's flanking domains.
    threshold_frac: Fraction of the peak effect defining the edge.

Returns:
    A dict with peak_m, threshold_m, and the updrift/downdrift extents in
    domains and meters.
```

</details>

### groin-sweep/HAT_groin_sweep_figures.py

Per-period diagnostic figures for one groin sweep: whether a winning cell is worth believing.

From the script's original header:

```text
Per-period diagnostic figures for one groin sweep.

The sweep orchestrator writes numbers; the joint fit draws the two-period
surface. Neither draws the per-period diagnostics that say WHETHER a winning
cell is any good: what the M-f error surface looks like on its own, and what
the winning shoreline-change curve looks like against CoastSat. This file
draws those, reading only what the sweep already wrote.

WHAT IS PLOTTED, AND WHY THAT AND NOT SOMETHING ELSE
    The heatmap carries TWO panels, not one, because the sweep's ranking
    metric and the reach-scale fit are different questions and can disagree:

        fillet_err   |modelled fillet - observed fillet| at the groin pair.
                     This is what the sweep RANKS on. The fillet saturates
                     (~18 m by 1990 at M = 40, where BRIE's diffusion balances
                     the trapping), so its LEVEL carries the information about
                     M while its slope -- which is what `differential` measures
                     -- carries almost none.
        rmse_window  RMSE of the modelled LRR against CoastSat across D1-D12.
                     A whole-reach number, dominated by domains the groin
                     never touches.

    If the two panels pick different cells that is a RESULT worth seeing, not
    a defect to average away: it means the groin parameters that reproduce the
    local fillet are not the ones that reproduce the reach.

    The profile figures plot the per-domain LRR, because that is the quantity
    the sweep is scored against -- the figure and the ranking then agree. The
    alternative (end-of-run shoreline position) shows the fillet's shape more
    directly but has no observed counterpart to overlay.

THE NOTCH CAVEAT, DRAWN RATHER THAN FOOTNOTED
    `observed_fillet_m` in the sweep config derives the observed fillet by
    fitting a regional trend across the fit window EXCLUDING the groin pair
    and taking the pair's departure from it. That anomaly is ONE-SIGNED --
    both D5 and D6 sit ABOVE the trend in 1984-2004 (+0.79 and +1.48 m/yr).
    A volume-neutral source/sink dipole must put a NOTCH at the downdrift
    domain, and the model does at every cell on the grid.

    So the observed shape is not a groin dipole, and the scalar fillet is the
    best available summary of a feature whose SHAPE the model cannot
    reproduce. Every profile figure that hits this case says so on its face,
    because a reader looking at a mismatched D5 should be told it is a
    property of the target and not a bad fit.

M = 0 IS NOT A COLUMN
    The M = 0 cells are the paired baselines the fillet is differenced
    against, so they have no fillet by definition and would be a blank column
    on the ranking panel. They are reported as a reference line in the panel
    subtitle instead, which is also the honest way to read them: the number to
    beat, not a candidate.

Usage:
    python HAT_groin_sweep_figures.py                     # every swept cell
    python HAT_groin_sweep_figures.py --period 2004 --preset zeroBE
    python HAT_groin_sweep_figures.py --top-n 8

Writes to output/calibration/groin/<start>_<end>_<preset>/figures/:
    heatmap.png              fillet error and reach RMSE over the M-f grid
    best_fit_profile.png     winning cell's LRR against CoastSat
    top_n_profiles.png       the best N cells on the same axes
    period2_surface.png      2004-2024 only: the M*f ridge, with contours
```

Notes that were in the code:

```text
parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
The guard below is what makes a future move fail here, loudly, instead of
resolving to scripts/scripts and surfacing as a missing data file several
imports deeper.
```

```text
Shared with HAT_groin_joint_fit.py so a groin is the same colour in every
figure the sweep produces.
House colours (2026-09-11): the observations are INK, the modelled cell
under test the ACCENT, and the structure's own marks a muted guide. An
orange, a dark red and near-black were chosen in this file, and the dark
red was the 1984 vintage colour doing a second job.
```

```text
The whole island, not just the fit window. Loaded lazily and cached: this
reads the transect file, and the module is imported by the comparison and
position scripts too.
```

```text
Two greys, not two hues: these bands say WHERE, and the colour in this
figure is reserved for WHAT is plotted.
```

```text
Both conditions matter: a strong correlation with a NEGLIGIBLE bias is
just a well-fit reach responding mildly to M.
```

```text
--- full reach ------------------------------------------------------
The fit window is 12 of 90 domains. A cell that matches the fillet says
nothing on its own about the other 78, and the groin's own influence is
measured (not fitted) out to a few km -- so this panel is where an
emergent extent can be checked against the observations that were never
part of the objective.
```

```text
One colour family, darkest is best: these are variations of one thing,
and `autumn` put the best cell in a red that means 1984 elsewhere.
```

```text
Constant-M*f hyperbolae, anchored on the products the GRID spans rather
than on the best cell. Anchoring on the best cell silently drew nothing
here: period 2's optimum sits at f = 0, so every product was 0 and each
contour was skipped as non-positive -- the one figure whose entire point
is to show the M*f ridge came out with no ridge on it.
```

```text
Draw the whole tied set, not one arbitrary member of it. A single star
on a tie reads as a fitted value; a row of them reads as what it is.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_matplotlib()`**

```text
Imports pyplot with a headless backend.

Deferred rather than imported at module scope so `--list` and the argument
parsing stay usable on a machine where matplotlib is missing.
```

**`load_scored()`**

```text
Loads one sweep's scored results.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".

Returns:
    A DataFrame of scored rows, or None if that sweep has no CSV yet.

Raises:
    ValueError: If the CSV predates the fillet-size rescore. Ranking such
        a frame would silently fall back to the differential -- the
        fillet's SLOPE -- which is the metric that railed at both grid
        edges. Re-running the orchestrator over a finished sweep is free
        (every cell resumes from disk) and adds the column, so this is a
        fixable error rather than a reason to plot the weaker metric.
```

**`profile_be1()`**

```text
Reduces rows to one per (M, f) by keeping the best-scoring be1.

be1 is a nuisance axis here: it exists in the 1984 edgeBE sweep only, and
the figure's subject is (M, f). Profiling it out -- keeping, for each
(M, f), the be1 that scored best -- is what the joint fit does, so the
surface drawn here and the surface fitted there are the same reduction.
The winning be1 is carried along so it can be annotated rather than lost.

Args:
    frame: Scored rows for one period/preset.

Returns:
    A DataFrame with one row per (M, f), sorted by the ranking metric.
```

**`rate_curve()`**

```text
The modelled per-domain LRR for one row, over the fit window.

Args:
    row: One row of a scored frame.

Returns:
    (gis, rate) arrays over FIT_DOMAINS_GIS.
```

**`observed_anomaly_is_one_signed()`**

```text
Whether both groin domains depart the regional trend the SAME way.

Mirrors the trend removal in `observed_fillet_m`. When true, the observed
pair is not a dipole and no volume-neutral source/sink can reproduce its
shape -- which the profile figures state on their face.

Args:
    period: 1984 or 2004.
    trend_order: Polynomial order for the regional trend. Must match the
        default in `observed_fillet_m` or the two disagree about the same
        data.

Returns:
    (is_one_signed, updrift_anomaly, downdrift_anomaly) in m/yr.
```

**`observed_curve_full()`**

```text
CoastSat per-domain LRR across ALL 90 domains.

`observed_curve` covers D1-D12 because that is the window the sweep is
SCORED on. Nothing about the data stops at D12 -- both periods have
transects in all 90 domains -- so the full reach is available for showing
where a groin fitted on 12 domains leaves the other 78.

Args:
    period: 1984 or 2004.

Returns:
    (gis, rate) arrays over the domains that have transects.
```

**`model_curve_full()`**

```text
One cell's modelled LRR across all 90 domains.

Read from the cell's own `shoreline_change_rate.csv` rather than from the
rate_D1..rate_D12 columns of sweep_results.csv, which carry the fit window
only. `lrr_m_yr` is the column, not `change_rate_m_yr`: the sweep is scored
on the OLS slope, and rate_D* is that same quantity.

Returns:
    (gis, rate) arrays, or (None, None) if the cell has no rate file.
```

**`_footnote()`**

```text
Registers `text` as the figure's caption.

Until 2026-09-11 this drew the text on the canvas, wrapped with `textwrap`
at a fixed column. The house rule is that nothing on the image belongs in
a caption, so it now lands in a CAPTIONS.md beside the PNG, written when
the figure is saved. Every caller is unchanged, and `width` is accepted and
ignored -- a caption file has no columns. Calling it twice on one figure
replaces the text, so the second call must carry everything: see
`_notch_note`, which appends rather than overwriting.
```

**`reach_panel_is_bias_driven()`**

```text
Whether the reach RMSE panel is tracking bias rather than groin skill.

THE FAILURE THIS CATCHES, measured on the 1984 zeroBE sweep. With
background erosion switched off the whole reach is far too accretional
(bias +2.50 m/yr against CoastSat), and `rmse_window` is dominated by that
offset rather than by shape. Trapping sediment adds net erosion, so
raising M shaves the bias and the reach RMSE falls MONOTONICALLY with M --
corr(M, bias) = -0.96 -- until it rails at the largest M on the grid.

A reader seeing that panel rail at M = 110 would reasonably conclude the
reach fit wants an enormous groin. It does not: it wants background
erosion, and M is the only knob on the grid that can supply any. Flagged
on the figure so the railing cannot be quoted as a groin result.

Args:
    groin: Scored rows with M > 0 for one period/preset.
    threshold: Correlation at or below which bias is called dominant.

Returns:
    (is_bias_driven, correlation, mean_absolute_bias).
```

**`tied_best()`**

```text
The best cell and every cell statistically tied with it.

WHY THIS IS NOT `idxmin`. In 2004-2024 the observed fillet is NEGATIVE
(-43.2 m: the fillet RELAXED across the window), and no M >= 0 can build a
negative fillet. The whole f = 0 row therefore scores identically -- a
fully deteriorated groin traps nothing, so M has no effect at all -- and on
the 2026-08-23 sweep seven cells from M = 40 to M = 110 tied to within
5e-5 m. `idxmin` picks one of them arbitrarily and the figure then
announces "best: M=95", which is not a fitted value: it is whichever tied
cell numpy happened to reach first.

Reporting the tie is the honest version, because the tie IS the result --
it says this period cannot constrain M.

Args:
    groin: Scored rows with M > 0.
    column: Metric to rank on.
    rel_tol: Fraction of the best score within which a cell counts as tied.

Returns:
    (best_row, tied_frame). `tied_frame` always contains at least the best
    row, so `len(tied) > 1` is the test for an unidentified parameter.
```

**`fig_heatmap()`**

```text
Two-panel M-f error surface: the ranking metric and the reach fit.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".
    surface: One row per (M, f), from `profile_be1`.
    out_dir: Directory to write into.

Returns:
    The written path.
```

**`fig_top_n_profiles()`**

```text
The best N cells on one set of axes, against CoastSat.

Shows how tightly the ranking discriminates: curves that lie on top of
each other mean the metric cannot tell those cells apart, which is a
statement about identifiability rather than about the fit.
```

**`fig_period2_surface()`**

```text
The second period's error surface with constant-M*f contours drawn on.

WHY THIS FIGURE EXISTS. 2004-2024 sits entirely past the 2003 end of the
deterioration ramp, so its cumulative trapping is 20*M*f: only the PRODUCT
is identifiable and the surface is a valley running along a hyperbola, not
a bowl with a minimum. Drawing the constant-M*f contours makes that
visible -- if the valley floor follows a contour, the non-identifiability
is shown rather than asserted in a caption.

Args:
    period: Must be the second period; the figure is meaningless for the
        first, which straddles the ramp.
    preset: "edgeBE" or "zeroBE".
    surface: One row per (M, f), from `profile_be1`.
    out_dir: Directory to write into.

Returns:
    The written path.
```

**`figures_for()`**

```text
Draws every figure for one swept period/preset.

Args:
    period: 1984 or 2004.
    preset: "edgeBE" or "zeroBE".
    top_n: How many cells the top-N profile figure carries.

Returns:
    A list of written paths, empty if that sweep has no results yet.
```

</details>

### groin-sweep/HAT_groin_sweep_worker.py

One cell of the groin / background-erosion sweep, in its own process.

From the script's original header:

```text
One combination of the groin / background-erosion sweep.

Runs a SINGLE (period, preset, M, deterioration_fraction, be1) combination in
its own fresh process, scores it against that period's CoastSat target, prints
one machine-readable line, and exits. Launched as a subprocess by
`HAT_groin_sweep.py` -- never in a loop inside one interpreter.

WHY A SUBPROCESS PER COMBINATION
    An earlier in-process sweep (all combinations back-to-back in one Python
    process) died partway through with a Windows access violation
    (0xC0000005): the signature of state accumulating across many repeated
    Cascade/Barrier3D constructions, not a bug in any single run. Process
    isolation makes the OS reclaim everything between combinations.

WHAT IS SWEPT, AND WHAT IS HELD FIXED
    Swept:  M (trapping rate) and f (deterioration floor), jointly across both
            periods, plus be1 in period 1 under edgeBE only.
    Fixed:  scenario    full_management (roadway + beach_dune + fills)
            Hs          2.5 m
            groin       GIS 6 updrift / GIS 5 downdrift, installed 1969,
                        linear_ramp deterioration, onset 1996, ramp 7 yr --
                        the SAME schedule in both periods, so period 2's
                        collapse to M*f emerges rather than being hard-coded
            be90        from HATTERAS_BE_RATES_EDGE, not pinned here
    See the orchestrator's docstring for why M and f need both periods.

WHAT IS SHARED, AND WHAT GUARDS THE REST
    `build_cascade` is imported from `cascade_pipeline.hindcast` -- the same
    definition the notebook and `HAT_hindcast_1984_2024.py` use, so a model
    built here is built exactly as a published run is. Only
    `run_cascade_simulation` below is this file's own, and it differs on
    purpose (see below).

    That difference is still guarded numerically, not by a source hash: the
    orchestrator reads the period's published matrix run's own M, f and be1
    out of its metadata, runs that combination through this worker, and
    refuses to report anything if the two rate curves differ. The reference's
    parameters are read rather than pinned, because the reference gets re-run
    as the fit is refined.

    Everything else that CAN come from the shared packages does: the presets,
    road events, nourishment projects and community zones from
    `hatteras_site_config`; the loaders, audits and shoreline maths from
    `cascade_pipeline`.

TWO DELIBERATE DIFFERENCES FROM THE HINDCAST RUNNER
    - `cascade.save()` is NOT called. It writes a ~140 MB npz per run; the
      344 cells of the full sweep would be ~48 GB of state nothing downstream
      reads. The shoreline matrix and the rate curve are written instead
      (~50 KB).
    - No figures, no GIFs, no run_index.csv row. A sweep combination is not a
      run of the matrix and must not be filed as one.

Usage:  HAT_groin_sweep_worker.py <period> <preset> <M> <fraction> <be1|none>
                                  <out_dir>
Prints: RESULT_JSON={...} on success. Non-zero exit means the combination
        failed; the orchestrator records it and moves on.
```

Notes that were in the code:

```text
parents[3], not [2]: this file lives in scripts/hatteras_ms/groin-sweep/.
The guard below is what makes a future move fail here, loudly, instead of
resolving to scripts/scripts and surfacing as a missing data file several
imports deeper.
```

```text
The pre-AST groin hook lives only in the sandbox copy. Attaching a callback
to stock Cascade silently does nothing, so which Cascade is imported is not
optional -- cascade_pipeline.hindcast makes that choice once, for the
notebook, the runner and this worker alike, via USE_SANDBOX_CASCADE.
```

```text
The period cannot come from a module constant any more: one worker serves
both. It is parsed here, at the top, because every path and forcing file in
section 1 is derived from it -- reading it later would mean two sources of
truth for the same value for the length of the module.

Parsed at import rather than inside main() for the same reason, and made
fatal rather than defaulted: a worker that silently fell back to 1984 would
fill a 2004 sweep directory with period-1 results that look entirely valid.
```

```text
Resolved from the extractor, not pinned -- the same source the hindcast
runner and the road setbacks use, so the sweep cannot drift out of step with
them. See scripts/site_layer/hat_topo_version.py.
```

```text
THE PRODUCT MUST BE NAMED. topo_dirs() with no argument resolves
DEFAULT_PRODUCT, which is "2004-start" -- so a 1984 sweep built its model on
the 2004 island while its reference matrix run used 1984-start. All 90
domains differ between the two products and 65 differ in interior SHAPE, so
the sweep was fitting M and f against a different barrier. The orchestrator's
drift guard caught it on 2026-08-29 at 0.126 m/yr, 25x its tolerance, and
blamed "section 3" -- which was in sync all along.

This line used to carry a comment saying the bare call "still resolves
2004-start, which is what this sweep read before the 2026-08-25
restructure". That was true and it was the bug: preserving pre-restructure
behaviour is exactly wrong once the periods stopped sharing a topography.
```

```text
Taken from what topo_dirs() RETURNED rather than re-joined from parts -
re-joining is how a resolver gets bypassed without anyone noticing.
```

```text
Cascade resolves its parameter file relative to cwd, exactly as the runner
does. Every path above is already absolute, so this only affects CASCADE.
```

```text
CONTINUOUS-WINDOW OVERRIDES. Unset, these change nothing and the worker runs
exactly the period the site config defines. Set, they let ONE run span both
hindcast periods -- 1984-2024 rather than 1984-2004 -- which is what makes M
and f separable: the fillet builds in the first half and declines in the
second, so the two knobs are constrained by different stretches of one curve
instead of trading off inside a single 20-year window.

They are environment variables and not a new HATTERAS_PERIODS entry on
purpose. The site config is keyed by START YEAR, so a 1984-2024 entry would
collide with 1984-2004, and every consumer that loops over HATTERAS_PERIODS
-- the run matrix above all -- would silently gain a period it was never
meant to run.
```

```text
SCENARIO = "full_management", expanded. Written out rather than looked up so
a sweep combination cannot inherit a scenario-table edit made for the matrix.
```

```text
ASYMMETRIC SOURCE/SINK. 1.0 is the volume-neutral pair the module has always
used; 0.0 makes the groin a pure updrift source, the physical statement that
the trapped sand comes from outside the pair rather than from the downdrift
beach. Read from the environment so a variant is a separate sweep in its own
directory, never mixed into the symmetric results.
```

```text
Hs IS A CALIBRATION VALUE, NOT AN OBSERVATION. The runner records it as
"calibration value 2.5" and notes that a sweep over it belongs in a separate
script -- so overriding it here is using the knob as intended, not bending a
measured constant.

It matters for the groin because BRIE's alongshore diffusivity scales roughly
as Hs^2.5, and diffusivity sets two things that pull in opposite directions:

fillet DECAY rate  ~ diffusivity        (too slow at Hs = 2.5: the
modelled fillet relaxes -15.3 m
over 1984-2024 against -24.4 m
observed, even with no groin)
fillet EXTENT      ~ sqrt(diffusivity)  (currently 2,000 m modelled
against 2,250 m observed, i.e.
already slightly narrow)

So raising Hs speeds the decay toward the observations while widening the
fillet away from them. The diagnostic exists to measure that trade rather
than argue it. Unset, this changes nothing.
```

```text
build_domain_file_paths() used to be COPIED here from the hindcast runner,
names and all. Both copies spelled f"domain_{gis_id}_topography_{init_year}"
and both had to be edited in step. It now comes from cascade_pipeline, which
itself delegates to hat_topo_version.domain_arrays() - one definition of the
directory and one of the filename, for every reader in the repo.
```

```text
The island offset comes from cascade_pipeline.hindcast.build_island_offset,
the runner's own builder, in the sweep's offset mode. Until 2026-09-24 this
file carried a copy of the runner's old loader, which divided the metre file
by ten and had no mode: every groin sweep ran on offset / 10 whatever the
runner did.
```

```text
No product argument: this worker sweeps the DEFAULT_PRODUCT, the same one
topo_dirs() resolved above for DUNE_TOPO_DIR.
The PRODUCT is passed, exactly as the runner does at its
section 3. Omitting it silently selects DEFAULT_PRODUCT.
```

```text
--- source/sink --------------------------------------------------------
Built directly rather than through resolve_be_preset(): where be1 is
swept the sweep varies a value the presets hold fixed, so calling it
"edgeBE" would file a swept run under a preset name whose numbers it
does not use. The SHAPE is the preset's and nothing else --
zeroBE  nothing imposed anywhere, so there is no knob to sweep and
the groin has to account for the southern deficit alone.
edgeBE  two end domains, interior untouched. be1 is swept in period 1
and taken from the site config in period 2, where no prior fit
stands behind a bracket.
```

```text
`build_cascade` is NOT here: it is imported from cascade_pipeline.hindcast,
the same definition the runner and the notebook build from, so the three can
no longer disagree about how a model is constructed.

`run_cascade_simulation` below is a deliberately different function, not a
copy -- it collects the per-year events instead of printing them, writes no
artifacts, and RAISES on a missing groin hook where the runner warns. The
module docstring states why.
```

```text
`row` already carries `year` (apply_to_cascade puts it
there), so passing year=current_year as well is a duplicate
keyword and raises. The runner spreads the row without
re-supplying year; this now matches it.
```

```text
`row` carries its own `kind` ("bridge", "relocation", ...) --
the runner switches on it. A dict literal lets the row win
rather than colliding, which `dict(kind=..., **row)` did.
```

```text
The failure mode this guards is silent: with the hook missing, a groin
run produces a valid-looking no-groin result and the sweep would fit M
against a model M never touched.
```

```text
IN A CONTINUOUS WINDOW THESE TARGETS DO NOT APPLY. OBSERVED_LRR and
OBSERVED_DIFFERENTIAL are keyed by START_YEAR and describe that period's
own 20-year window. Scoring a 40-year run against them would silently
compare a 1984-2024 rate curve to a 1984-2004 target and write the result
into the CSV as if it meant something. The modelled rates are still
returned -- they are a property of the run, not of any target -- and the
continuous sweep scores the per-domain CHANGE PROFILE itself.

NOTE FOR ANY FUTURE CALLER: HAT_groin_sweep.py treats a NaN
differential_err as "this cell failed", so the production orchestrator
must NOT be pointed at a continuous window. Use the continuous driver.
```

```text
Extent is measured by the orchestrator, not here: pairing a groin run with
its M = 0 baseline is a property of the grid, and the baseline may not have
finished when this combination does. See measure_groin_extent in
HAT_groin_sweep_config.
```

```text
CASCADE's Cascade.__init__ reaches brie_coupler.initialize_equal, which calls
set_yaml("TMAX", ...) -- and set_yaml opens the SHARED parameter file for
reading, then reopens it for writing:

with open(file_name) as f: doc = full_load(f)
doc[var_name] = new_vals
with open(file_name, "w") as f: dump(doc, f)

Every worker points at the same data/hatteras_init/Hatteras-CASCADE-
parameters.yaml. Run two at once and one truncates the file while the other
is reading it: full_load returns None and the next line raises
"TypeError: 'NoneType' object does not support item assignment". It is
invisible in a serial run and fires on essentially every cell in a parallel
one -- 215 of 215 in the first attempt at this sweep.

Serialising CONSTRUCTION alone fixes it. Construction is seconds; the
20-year update loop it precedes is minutes and stays fully parallel, so the
pool keeps almost all of its speedup.

The alternative -- giving each worker a private copy of the datadir -- was
rejected as the bigger change: the yaml sits beside the init surfaces the
runs read, and duplicating that per worker trades a lock for a much larger
surface of things that can go stale.
```

```text
Windows refuses to unlink a file another process still has open, and both
the release and the stale-steal below race exactly that. Unguarded, a single
WinError 32 killed the worker AND left the lock behind, after which every
later cell blocked on it -- 21 of 43 cells lost to one transient collision,
and 17 hours of wall clock, on 2026-08-22. Retry briefly, then give up
quietly: an orphan is self-healing via CONSTRUCT_LOCK_STALE_S, whereas an
exception here is not.
```

```text
A failed steal is not fatal -- another worker may be
stealing the same lock this instant. Fall through and let
the loop retry rather than dying on the collision.
```

```text
M = 0 means no groin, not a groin of zero strength: GroinCallback would
still be attached and still append diagnostics, and a "0 m/yr groin run"
in the record is a run someone will later mistake for a failed one.
```

```text
BRIE's diffusion number at the initial state, for the analytic fillet
prediction. Recorded per combination because it is what makes the
emergent extent checkable against something written down beforehand.
```

```text
Both estimators, and the sweep is SCORED on the LRR. OBSERVED_LRR is a
per-transect OLS slope, so fitting M against an endpoint difference
would fit the groin to the gap between two estimators as well as to
the groin. The endpoint column is still written: the drift guard in
HAT_groin_sweep.py compares it against the published matrix run.
```

```text
Road drowning happens inside CASCADE's roadway_manager, not through the
historical-event list, so it has to be read back off the finished model
-- the same source the hindcast runner's section 12.2 uses. Counting
`events` here would report 0 for a run that drowned eight road_offset.
```

```text
WHAT THIS CELL WAS BUILT ON. Added 2026-08-31, and the reason is
specific: result.json recorded M, be1, f, period, preset and scores
and NOTHING about its inputs, so a sweep cell could not be told apart
from a cell of the same combination built on different topography.

That is not hypothetical. The worker resolved topo_dirs() without
naming a product until 2026-08-30, so DEFAULT_PRODUCT ("2004-start")
answered and a 1984 sweep built on the 2004 barrier. When the
question "are these cells stale?" was asked on 2026-08-31, the
production runs answered in seconds from their own topo_product and
be_values_digest columns; the sweep could only be answered by
RE-RUNNING CELLS and diffing, at ~6 minutes each.

These four fields are named to match run_index.csv exactly, so a
sweep cell and a matrix run can be compared without a translation
step.
```

```text
The offset mode, since 2026-09-24: before that every cell was offset
/ 10 and nothing said so. A result.json without this field is asrun.
```

```text
The effective rate the groin actually applied, summed over the run
and divided by its length. This is the quantity the two periods have
in common: period 1 delivers M*(16 + 4f)/20 and period 2 delivers
M*f. Recorded so the joint fit does not have to re-derive it from
the schedule, and so a period-2 row states plainly that M and f
reached it only as a product.
```

```text
Nested dict flattened so the orchestrator's DataFrame gets one column
per domain rather than a column of dicts.
```

```text
argv was parsed and validated at import (section 0), because the period
decides every path in section 1. Nothing is re-read here.
```

<details><summary>Function notes (the original docstrings)</summary>

**`_parse_cli()`**

```text
Reads the six positional arguments, or exits with usage.

Returns:
    (period, preset, M, fraction, be1, out_dir). be1 is None when the
    caller passed 'none', meaning "use the preset as it stands".
```

**`build_background_erosion()`**

```text
Expands sparse per-GIS background-erosion rates onto the padded array.

Copied from the hindcast runner. See its section 4.3.

Args:
    be_rates: Mapping of GIS domain ID to rate in m/yr. Domains absent
        from the mapping get 0.0.
    geometry: DomainGeometry describing the padded array.

Returns:
    A list of geometry.total_domains rates, ready to pass to Cascade().

Raises:
    ValueError: If a GIS ID falls outside the padded array.
```

**`assemble_forcing()`**

```text
Builds every per-domain forcing array the model needs.

Everything here is independent of the groin, so it is identical across
combinations that share `be1`. It is rebuilt per process anyway -- the
whole point of process isolation is that nothing is shared.

Args:
    be1: Background-erosion rate at GIS 1, m/yr, or None to take the
        site-config value for this period. Ignored under zeroBE.

Returns:
    A dict of the arrays build_cascade needs, plus the schedule and
    management masks the run loop and the audits use.
```

**`run_cascade_simulation()`**

```text
Steps a built Cascade through its period.

Differs from the runner's copy only in what it writes: no artifacts are
saved here (the caller writes the two small files it needs) and the
per-year event log is collected rather than printed.

Args:
    cascade: A Cascade from build_cascade, not yet stepped.
    run_years: Annual transitions to simulate.
    name: Run name, recorded on nourishment log rows.
    start_year: Calendar year of the run's first state.
    geometry: DomainGeometry, for GIS <-> pad translation.
    alongshore_section_count: Padded domain count.
    historical_road_events: RelocationEvent / BridgeEvent sequence.
    relocations_enabled: Global toggle for relocation events.
    setback_check: {gis: measured_setback_m} reported beside relocations.
    nourishment_schedule: A NourishmentSchedule, or None for no fills.
    groin_callback: The attached GroinCallback, for the drive check.

Returns:
    A (cascade, events) tuple. `events` is a list of dicts recording the
    nourishment and roadway actions the loop took.

Raises:
    RuntimeError: If a groin was attached but never called -- the
        pre-AST hook is missing and the run is silently a no-groin run.
```

**`score_combo()`**

```text
Scores one finished run against the CoastSat 1984-2004 target.

Four metrics, all computed, one ranked on. The differential is the ranking
metric because it is the only one that identifies M: over D1-D12 the
profile RMSE moves 7% while M moves 4x, whereas the differential moves 5x.
The profile numbers are diagnostics, and the extent is neither -- it is
the pre-registered check from `cascade.groin`, reported so it stays a test
rather than a target.

Args:
    model_lrr: Padded model LRR array, m/yr, seaward-positive -- the
        OLS slope through every annual state, the same estimator
        OBSERVED_LRR is. NOT the endpoint difference.
    geometry: DomainGeometry describing the padded array.

Returns:
    A dict of the four metrics plus the modelled per-domain rates over
    FIT_DOMAINS_GIS.
```

**`_release_construct_lock()`**

```text
Removes the construction lock, tolerating Windows sharing errors.

Returns:
    True if the lock is gone, False if it could not be removed. False is
    not an error: the stale-steal path reclaims it after
    CONSTRUCT_LOCK_STALE_S, which is strictly better than raising and
    taking the worker down with the lock still held.
```

**`cascade_construction_lock()`**

```text
Serialises Cascade construction across worker processes.

A stale lock is stolen rather than waited on: a worker killed mid-
construction (the 0xC0000005 this design already expects) would otherwise
block every later cell until the sweep's own timeout.

Raises:
    RuntimeError: If the lock cannot be taken within the timeout.
```

**`ensure_parameter_file_intact()`**

```text
Repairs the shared parameter file if a previous run left it broken.

The lock stops two workers interleaving their writes, but it cannot help
if a worker is KILLED between set_yaml's truncating open and its dump --
and 0xC0000005 does exactly that on this configuration. The file is left
half-written, and because it is shared, every remaining cell of the sweep
then dies on a YAML scanner error. That is how a single transient crash
turns into an overnight run with zero results.

Every value CASCADE writes here (TMAX, the shoreface terms, the three
file paths) is recomputed at the next construction, so restoring a
snapshot loses nothing: only the hand-set parameters matter and those are
identical in every version of this file.

Must be called while holding the construction lock.
```

**`build_cascade_locked()`**

```text
`build_cascade` under the construction lock.

A wrapper rather than a lock inside `build_cascade`, so the shared
definition in cascade_pipeline.hindcast stays free of sweep-only
concurrency machinery.
```

**`run_combo()`**

```text
Builds, runs and scores one combination.

Args:
    M: Groin trapping rate, m/yr. 0.0 attaches no groin at all -- the
        paired baseline column.
    be1: Background erosion at GIS 1, m/yr.
    fraction: Groin deterioration floor, in [0, 1].
    out_dir: Directory for this combination's two output files.

Returns:
    A result dict, written to disk as result.json and printed as
    RESULT_JSON= for the orchestrator to parse.
```

</details>

### groin-sweep/HAT_groin_timeseries_check.py

Does the chosen groin hold up through time, not just at the end year?

From the script's original header:

```text
Does the chosen groin hold up THROUGH TIME, not just at the end year?

Every other figure here compares a single end state. That cannot distinguish a
groin that tracks the observations year by year from one that wanders and
happens to arrive in the right place -- and it cannot show WHEN a fit starts to
fail. This one plots the fillet's whole trajectory against the surveys.

WHAT IS PLOTTED
    Three curves per panel:

      observed   the surveyed fillet from `Change_from_wetdry_1967_D2_D12.csv`,
                 re-referenced to the period's start year so it begins at zero
                 like the model does. Markers only on the years actually
                 surveyed -- the record is irregular and joining it with a line
                 would imply samples that do not exist.
      groin ON   the chosen cell, M = 60, f = 0.6.
      groin OFF  the paired M = 0 run at the same be1. The gap between the two
                 model curves is the groin's contribution; the gap from ON to
                 observed is what the source/sink field is left to absorb.

WHY THE MODEL IS RE-REFERENCED TOO
    A run starting in 1984 or 2004 inherits the real fillet in its initial
    shoreline, so its absolute D5-D6 offset is not comparable to a survey
    measured from a 1967 datum. Differencing both sides against their own start
    year removes the inherited part and leaves the CHANGE, which is the only
    quantity the two share.

WHAT THIS FIGURE IS EXPECTED TO SHOW
    A shortfall, and a documented one. M = 60 / f = 0.6 was not chosen to match
    the fillet -- no admissible M can, on this grid -- but by a DIRECT FIT to
    the period-1 D4-D8 change profile (demeaned RMSE 11.69 m against 15.58 m
    with no groin), bounded above by affordability (719,000 m3/yr against a
    5-7e5 littoral drift). The stability bounds an earlier version of this text
    cited -- "M >= 70 unstable, M >= 100 drowns" -- were measured on the
    41-domain rig and do NOT transfer: every production cell through M = 160
    ran clean.
    The residual this figure shows IS the quantity that calibration absorbs,
    together with the Cape Point dynamics the dipole does not represent. Read
    it as the split between the two, not as a failed fit.

Usage:
    python HAT_groin_timeseries_check.py
    python HAT_groin_timeseries_check.py --M 50 --fraction 0.6

Writes output/calibration/groin/figures/timeseries_check.png
```

Notes that were in the code:

```text
House colours (2026-09-11): surveys in INK, the run under test the
ACCENT, the groin-off run BASE grey.
```

```text
be1 is swept only in the 1984 edgeBE sweep. -40 is the grid value nearest
production's -41.8, and is the be1 the period-1 D4-D8 fit was pinned at, so
this figure and the M = 60 choice rest on the same cell. (An earlier version
used -34, the reach-RMSE minimiser, which is no longer what M is fitted on.)
```

```text
Re-reference the survey to this period's start: the model begins at
zero, so the observations must too.
```

```text
Whole years: the default locator put 1987.5 on a year axis, and
five of those labels collide at the printed width.
```

### groin-sweep/HAT_groin_zoom_gifs.py

Animated D2-D12 shoreline for a selection of (M, f) cells, for eyeballing the fit.

From the script's original header:

```text
Animated D2-D12 shoreline for a selection of (M, f) cells, for eyeballing the fit.

The existing run GIFs cover D1-D15 and exist only for the pair that was run as a
hindcast. Every SWEEP cell carries its own shoreline_matrix.npy, so any (M, f)
can be animated -- which is what is needed to judge whether M = 60 / f = 0.6 is
actually the best pair rather than just the best-scoring one.

Each frame shows shoreline CHANGE since 1984, demeaned over D4-D8 so the
alongshore shape is what the eye compares (a uniform level offset belongs to the
source/sink term, not the groin). The no-groin cell at the same be1 is drawn on
every frame as the reference, and the observed 1984->2004 change is drawn as a
fixed target so the endpoint can be judged.

    THE OBSERVED TARGET IS FIXED, AND THAT IS A LIMIT OF THIS WINDOW, NOT A
    CHOICE. Inside 1984-2004 the observation IS a single endpoint (an OLS fit
    to CoastSat chainage evaluated at both ends). The full-life companion,
    HAT_groin_full_life_gif.py, animates the observations too, because the
    1967 wet/dry record carries 19 dated surveys inside its window.

ORIENTATION: EROSION IS UP, SO THE PANEL READS AS A PLAN VIEW
    Changed 2026-08-30. These gifs were SEAWARD-positive and the full-life gif
    was LANDWARD-positive, so two animations in the same folder had opposite y
    axes. They are all landward-positive now: a retreating shoreline moves UP
    and the reader is looking down on the island with the ocean below the axis.

    Sign handling, which differs per source and is the easy thing to get wrong:
      `shoreline_matrix.npy`        Barrier3D's x_s, ALREADY landward-positive
          and already in METRES. Used as-is -- the negation that used to be
          here is gone.
      `observed_change_profile()`   chainage, SEAWARD-positive (verified in
          HAT_fullperiod_target.py against the published CoastSat LRR). It is
          NEGATED here.
    So the flip did not simply move a minus sign; it moved it from one series
    to the other.

DOMAINS D2-D12, matching the extent of the 1967 wet/dry survey and centred on
the groin at D5/D6. D1 is excluded: the cape's change over period 1 is 81-104 m,
about five times the groin's signal, and it swamps the axis.

Writes output/calibration/groin/figures/zoom_gifs_D2_D12/
```

Notes that were in the code:

```text
NO dam->m CONVERSION. shoreline_matrix.npy is written in METRES already.
Multiplying by 10 on 2026-08-30 made every curve ten times too large and put
it off a +/-70 m axis. Checked against the cell's own shoreline_change_rate.csv:
(m[-1] - m[0]) gives 82.9 m at D4 where the CSV reports 84.1, agreeing to
the endpoint-vs-LRR estimator difference. Anything that rescales this must be
re-checked against that CSV.
```

```text
Okabe-Ito, colour-vision-safe and muted enough to print. Shared with
HAT_groin_full_life_gif.py so the two animations read as one set.
House palette (2026-09-11), replacing an Okabe-Ito set chosen in this file.
```

```text
The pairs worth comparing: the chosen value, its f-neighbours, and M values
either side -- plus a high-M cell where the module actually draws a dipole.
```

```text
What each cell is in the folder for. Shown as the subtitle so a gif opened on
its own still says why it exists.
```

```text
The house style is the whole of it now; the local rcParams block that used
to sit here set its own type stack, ink and tick sizes.
```

```text
observed_change_profile is SEAWARD-positive; negate it onto the landward-
positive axis these panels now use.
```

```text
Title, subtitle and legend stay ON the canvas: a GIF is watched
standalone, with no caption file beside it in a viewer.
```

```text
-98 not -70: M = 160 dives past -75 at D12 and ran into the
fillet readout. All six cells share one limit so they stay
directly comparable.
```

```text
26, not 40: the axes start at 0.135 of the figure width, so a
bigger pad pushes the label off the canvas. The orientation
cues moved further out instead.
```

<details><summary>Function notes (the original docstrings)</summary>

**`load()`**

```text
Change since 1984 per domain, LANDWARD-positive, metres.

Two directory spellings, because the no-groin cell has no deterioration
floor to name: the swept cells are `M<M>_be<be>_f<f>` while the baseline on
disk is `M0_be<be>`. Only the suffixed name was tried until 2026-09-11, so
the baseline never loaded and the legend advertised a dotted line these
gifs never drew.
```

</details>

### groin-sweep/HAT_period1_top_n_figure.py

Top-N sweep results against observed change, period 1, production geometry.

From the script's original header:

```text
Top-N sweep results against observed change -- PERIOD 1, production geometry.

The direct counterpart to the 1967 rig's `HAT_groin_sweep_top_n_profiles.png`,
so the two calibrations can be read side by side. Same question: do the
best-scoring cells actually reproduce the observed alongshore SHAPE, or do they
only match a summary number?

THREE DIFFERENCES FROM THE RIG FIGURE, ALL DELIBERATE

    window      1984-2004, not 1967-2018. Period 1 is the only window in the
                hindcast where the observed gap between the groin's flanks
                WIDENS, which is the only behaviour a module with trapping
                >= 0 can produce.

    geometry    120-domain production grid, not the rig's 41. M is
                grid-specific -- a confined array preserves dipole amplitude
                that an open one diffuses away -- so a value fitted here
                transfers to the hindcast and one fitted on the rig does not.

    score       DEMEANED, and ranked on D4-D8 only. A uniform level offset in
                the groin's neighbourhood is absorbed by the source/sink
                calibration downstream, so it is not the groin's job; what the
                groin must get right is the shape. D1 is excluded from the
                ranking because the cape's change over period 1 is 81-104 m,
                roughly five times the groin's signal, and it swamps it.

WHAT TO LOOK FOR
    The no-groin baseline is drawn alongside the top cells. If the groin is
    doing real work the top cells should sit closer to the observations than
    that grey line does, in the shaded fit window. They do: 15.58 -> 11.69 m
    for the chosen cell (M = 60, f = 0.6), a 25% reduction.

    Watch also how tightly the top five bundle together. They span M = 40-95
    and f = 0.4-1.0 yet differ by under 0.5 m, which is the visual statement of
    the ridge in period-1 cumulative trapping, M(15.5 + 4.5f). An earlier
    version of this caption said "the metric identifies a product, not a
    pair"; fig_Mf_identifiability.png tested that and refuted it
    (corr(RMSE, M*f) = -0.07). See CALIBRATION_FIGURES.md and GROIN_PLAN.md.

Usage:
    python HAT_period1_top_n_figure.py [--top-n 5]

Writes output/calibration/groin/figures/period1_top_n_profiles.png
```

Notes that were in the code:

```text
The CALIBRATED value, re-solved 2026-08-28. This was -40.0, which was the
nearest grid point to the old -41.8 rather than the value itself; -42.6
is now on the grid and is what production spends.
```

```text
House semantics: the cell taken forward is the ACCENT, the no-groin baseline
is BASE, the structure's position is a guide line in muted ink. GROIN_COLOR
was a dark red here, which is the 1984 vintage colour elsewhere.
```

```text
LANDWARD-POSITIVE, so erosion is UP and the panel reads as a plan view,
matching the gifs. observed_change_profile is SEAWARD-positive, so it is
negated; x_s is landward-positive already, so the negation that used to
sit on `change` below is gone. Both series flip together, so `score` is
unchanged, and the scoring pipeline itself is untouched.
```

```text
Demeaned over the FIT window, so every curve is centred the same way the
score centres it. Centring on the plotted window instead would show a
different quantity from the one that was ranked.
```

```text
One colour family for the top cells, darkest is best, so they read as
variations of one thing. viridis was used here until 2026-09-11.
```

### groin-sweep/HAT_profiles_by_M_figure.py

Every M value as its own panel, with all six f curves on it.

From the script's original header:

```text
Every M value as its own panel, with all six f curves drawn on it.

The top-N overlay (fig_top_profiles.png) shows only the best cells, and they sit
so close together that nothing about the parameter response is visible. This
draws the whole grid instead: one panel per M, six f curves inside it, the same
observed and no-groin reference on every panel, and a shared y axis so panels
can be read against each other.

WHAT TO LOOK FOR
    * within a panel: how much f moves the profile at fixed M. Period 1 mostly
      PRECEDES the 1996-2003 deterioration ramp, so f should move it little --
      period-1 cumulative trapping is M(15.5 + 4.5f), which f changes by only
      29% across its whole range.
    * across panels: M lifts the whole curve rather than building a local
      fillet at D5/D6. That is the finding the per-domain decomposition in
      fig_d4d7_window.png makes numerically -- the groin's gain comes from D4,
      outside the dipole, while D5 (downdrift) gets worse.

Writes output/calibration/groin/figures/profiles_by_M/
    fig_all_M_profiles.png     the grid, for comparison across M
    fig_M<value>.png           one file per M, for detail
```

Notes that were in the code:

```text
The f family as one colour family, so a panel reads as "one M, six f" and
not as six unrelated series; darkest is the largest f. viridis was used here
until 2026-09-11 and shared its green with the reference marks elsewhere.
```

```text
Top widened on the 2026-08-30 flip to landward-positive: the high-M cells
(M >= 110) reach past +70 at D12 and were clipping.
```

```text
The per-cell error is panel-specific: in the grid these handles feed
ONE shared legend, so quoting a number there would attach panel (a)'s
errors to every panel.
```

```text
Four to a row: the letter and a full sentence collide, so the panel
carries the M and the caption carries the best f per panel.
```

### groin-sweep/HAT_rig_f_bracket_figure.py

What the 1967 rig can and cannot say: f is bracketed, M is railed.

From the script's original header:

```text
What the 1967 rig can and cannot say: f is bracketed, M is railed.

WHY THIS FIGURE EXISTS
    The rig is quoted as corroborating BOTH parameters. It does not, and the
    difference matters enough to draw. `scripts/site_layer/hatteras_site_config.py` said
    until 2026-08-30 that in the rig "f = 0.6 is a clean INTERIOR minimum ...
    So is M ... Neither railed." The rig's own sweep CSV refuses the second
    half: RMSE improves monotonically to M = 60 and then the model blows up.

    So the rig resolves f and is only CONSISTENT with M. M is set by the
    production period-1 fit instead. This figure is the evidence for both
    halves of that sentence, in one place, so the claim cannot drift back.

WHAT IS PLOTTED
    (a) f AT M = 60 -- a clean interior minimum at f = 0.6, bracketed on both
        sides with steep curvature. This is the parameter the rig owns: it is
        the only window containing the 1996-2003 deterioration ramp, because
        both hindcast windows begin 15 years after the structure went in.

    (b) M ACROSS THE WHOLE GRID, on a log axis because the failure spans three
        orders of magnitude. Every f-series improves monotonically to M = 60
        and then jumps ~13x at M = 70. Cells at M >= 100 do not complete at
        all and are drawn on the "crashed" rule at the top. M = 60 is the LAST
        VALUE THAT RUNS, not the value where the fit stops improving.

    The stability wall is RIG-SPECIFIC and does not transfer: all 36 production
    cells, including the full M = 70 and M = 80 rows, ran clean on the
    120-domain grid. Do not quote this ceiling for production runs.

Usage:
    python HAT_rig_f_bracket_figure.py

Writes output/calibration/groin/figures/rig_f_bracket.png
```

Notes that were in the code:

```text
House colours (2026-09-11). The series under test is the ACCENT; a cell
that did not complete is drawn in the vintage red, the one warm colour in
the palette, because "this run does not exist" has to be findable at a
glance; the f families in panel (b) are a ramp of the accent so they read
as one family. GOOD/BAD/INK/GRID were an orange, a pink and a navy chosen
in this file.
```

```text
A blank rmse is a cell that did not complete -- the barrier drowned or
the solver diverged. Those are information, not missing data.
```

```text
Outside: eight entries inside panel (b) sat on top of the low-M cells,
which are the ones the panel is about.
```

