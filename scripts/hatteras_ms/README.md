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
