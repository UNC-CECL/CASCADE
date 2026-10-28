# End domains and buffers

Kept on 2026-10-23: HAT_buffer_closure_test.py, HAT_cape_source_shape.py. The other scripts this README describes were retired (see ../RETIRED_STUDIES.md; recover from commit a7c4358c).

Solving the GIS 1 and GIS 90 end rates, the buffer ring beyond them, the Cape Point source shape and the Pea Island extension.

Each study drives the unchanged runner (`../../HAT_hindcast.py`); see `../README.md` for how
the experiments are run and filed, and how they import each other across folders.

| Script | What it asks |
|---|---|
| `HAT_resolve_ends_metres.py` | Re-solve the two end domains under the metres offset. |
| `HAT_resolve_ends_on_position_change.py` | Solve the two end domains on position change, not LRR. |
| `HAT_position_change_ends_figure.py` | Figure 3 of the wave recommendation, redrawn on the ends solved on position change. |
| `HAT_buffer_closure_test.py` | Do smoother or longer buffers stop the ring closure eroding the island's end domains? |
| `HAT_buffer_orientation.py` | Does the northern buffer's shape matter? The first 10 buffer cells follow the real coast north of GIS 90. |
| `HAT_cape_source_shape.py` | Which shape should the Cape Point sediment source take: all at GIS 1, or tapered over GIS 1-5? |
| `HAT_peaisland_extension.py` | Does the buffer's orientation matter? The domain set extended north over Pea Island to GIS 115. |

## The scripts in detail

Each script's header says what it does and how to run it. Below, for each,
is its original header and any notes that were in its code, kept word for
word when the scripts were brought in line with `scripts/STYLE.md`
(2026-09-30).

### HAT_resolve_ends_metres.py

Re-solve the two end domains under the metres offset.

From the script's original header:

```text
Re-solve the two end domains under the metres offset (2026-09-27).

Hannah, 2026-09-26/27: sweep the waves with the end source/sink terms FIXED,
"but we might need to resolve the ends now that we fixed the offset". The
stored values (HATTERAS_BE_EDGE_ONLY: 1996 +32.2 / +10.0, 2010 +72.6 / +31.3
m/yr) were solved at Hs 2.5 on the /10 offset. Chosen with Hannah:

    reference  the step-2 baseline wave climate, Hs 1.0 m, Tp 8 s,
               asymmetry 0.8, high-angle 0.45: a neutral setting, not the
               winner of either search, so the fixed ends do not pre-favour
               the sweep that uses them
    scenario   full management (road, beach and dune management, fills; no
               relocations, no groin), as the matrix end values always were;
               the pair is then used for natural and managed runs alike
    target     each window's CoastSat LRR: GIS 1 against the raw domain mean,
               GIS 90 against the LOWESS value (10 domains until 2026-09-28,
               7 since: common.SMOOTH_DOMAINS)
    solve      step 0 a fresh zeroBE run; then HAT_wave_shortlist_ends_solved's
               safeguarded step (secant capped at +-30 m/yr until a probe lies
               on each side of the target, then interpolation held inside the
               bracket), first step from the metres response measured on
               2026-09-26 (about 0.24 m/yr of residual per m/yr imposed at
               GIS 1, 0.20 at GIS 90); converged at |residual| <= 0.02 m/yr,
               at most MAX_STEPS probes
    output     tables/ends.json -- {period: {"1": rate, "90": rate}} -- read by
               HAT_wave_grid_fixed_ends.py. The config (HATTERAS_BE_EDGE_ONLY)
               is NOT changed.

WHERE: output/raw_runs/experiments/end-domain-boundaries/2026-09-27-ends-resolved-metres-offset/

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/ends/HAT_resolve_ends_metres.py

    --tag   file the solve under another study (2026-09-28: the LOWESS-7 re-solve)
    --seed  "1996=4.8394,17.545;2010=18.8,24.535": step 1 probes these ends
            instead of the first-gain guess from zeroBE, so a re-solve near a
            known answer starts there
```

Notes that were in the code:

```text
Matched on the wave settings (fixed 2026-09-27): the folder holds the
step-0 run of every reference solved so far, and taking the first one
found gave the Hs-2 and Hs-2.5 solves the Hs-1 run's residuals as their
step 0. The values those solves adopted were each confirmed by a probe
run at the right waves, so only the solver's path was affected.
```

```text
Re-solve one window at other waves (2026-09-27, Hannah: "re-solve the
2010 ends at Hs 2"): --periods and the wave settings; the result is
merged into ends.json, the replaced value kept under "history". The
closest probe (smallest worst-end residual) is the answer, and --accept
is the residual the caller will act on (the 0.02 target is not reachable
at GIS 1 in 2010-2024, whose response is not monotonic below ~0.1).
```

### HAT_resolve_ends_on_position_change.py

Solve the two end domains on position change, not LRR.

From the script's original header:

```text
Solve the two end domains on position change, not LRR (2026-09-28).

Hannah, 2026-09-28: "calculate the edge source sink values we should be
matching based on the position change instead" -- option 1 of three: an
experiment, the solver and the config unchanged. The LRR-solved ends match
the rates to 0.15 m/yr but overshoot the observed end-minus-start change at
GIS 1 1996-2010 (+38 m) and GIS 90 2010-2024 (+21 m), because the observed
ends do not move linearly and the model does.

    target     the CoastSat total change, mean of the end year's images minus
               mean of the start year's (5-scr/3-rates/coastsat/total_change/
               <window>/smoothed/tables/domain_smoothed.csv): GIS 1 raw, GIS 90
               LOWESS over common.SMOOTH_DOMAINS (7) -- the same treatment the
               LRR target gives each end (Hannah: "use the smoothed value at
               GIS 90"; built by HAT_wave_recommendation_figures.
               observed_change_smoothed)
    model      end-minus-start position, the run's change_rate_m_yr x 14
    residual   (model - observed) / 14 years, in m/yr, so HAT_resolve_ends_
               metres' solver, gains and 0.02 m/yr tolerance (0.28 m of
               change) apply as they are
    waves      option A (Hs 2.0, Tp 7.5, asymmetry 0.6, high-angle 0.5),
               full management, as the LRR solve
    seed       the LOWESS-7 LRR-solved ends (2026-09-28-ends-resolved-lowess7)

WHERE: output/raw_runs/experiments/end-domain-boundaries/2026-09-28-ends-solved-on-position-change/

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/ends/HAT_resolve_ends_on_position_change.py
```

<details><summary>Function notes (the original docstrings)</summary>

**`residuals()`**

```text
(model - observed) end-minus-start change, per year. `targets` (the
LRR target R.main builds) is not used.
```

</details>

### HAT_position_change_ends_figure.py

Figure 3 of the wave recommendation, redrawn on the ends solved on position change.

From the script's original header:

```text
Figure 3 redrawn on the position-change ends (2026-09-28).

Hannah: "redraw the figure with the position-change ends when done". The solve
(HAT_resolve_ends_on_position_change.py) runs full management only; this runs
the natural scenario at the same ends (option A waves), then draws
HAT_wave_recommendation_figures.fig3 on the two runs per window. The LRR-ends
figure is left as it is, for comparison.

    python scripts/real-world-cases/Hatteras/hatteras_ms/experiments/ends/HAT_position_change_ends_figure.py

Natural runs land under <study>/runs/natural_final/; the figure goes to the
study's figures/ folder.
```

### HAT_peaisland_extension.py

Does the buffer's orientation matter? The domain set extended north over Pea Island to GIS 115.

From the script's original header:

```text
HAT_peaisland_extension.py -- does the buffer's orientation matter?
THE QUESTION (Hannah, 2026-09-16). The hindcast models GIS 1-90 and pads
each end with 15 invented domains that extrapolate the local shoreline
slope, then bridge back to close BRIE's periodic ring; edgeBE pins GIS 1
and 90 to their observed rates with boundary source/sink terms that take up
whatever the buffer gets wrong (+32.2 and +10.0 m/yr in 1996-2010). What if the
buffer carried the REAL coast instead -- Pea Island north to GIS 115 --
with the ends re-solved there? (A one-domain southern extension was run and
removed the same evening: no domain polygon lies south of GIS 1.)

THE DESIGN
    geometries   n115 (GIS 1-115) against base
    offset modes asrun (the compressed planform every calibrated run uses)
                 and detrended (the planform at full strength; no calibrated
                 baseline exists, so base-detrended is solved here too)
    period       1996-2010, full_management, no groin, relocations off, Hs 2.5
    topography   1984-start CURRENT for GIS 1-90; the buffer profile beyond
    management   none on the extension (no road, fills, relocation, BE)
    stage 0      zeroBE on every member: the orientation effect with no
                 boundary term anywhere, against the zeroBE matrix run
    stage 1      the end domains solved by Newton steps (edgeBE, one probe
                 per step through HAT_BE_OVERRIDE), against the edgeBE
                 matrix run
    score        interior RMSE on GIS 2-89 against the SURVEYED target (the
                 runner's rmse_interior_m_yr, identical in meaning for every
                 geometry), the solved end values, and the rates on GIS 80-90

WHERE THINGS ARE
    inputs   2-brie-offset/1996/ext/<geometry>/          the offsets
             5-scr/3-rates/coastsat/lrr/1996_2010/ext/            the targets
    runs     output/raw_runs/experiments/topography-and-domains/2026-09-16-pea-island-domain-extension/<member>/
             one member per <geometry>-<mode>: its 1996_2010/zeroBE/ run is
             stage 0, its step<k>/ folders are the Newton probes, and SOLVED
             names the step that stands as the solved run
    logs     output/raw_runs/experiments/topography-and-domains/2026-09-16-pea-island-domain-extension/logs/<member>/
    answer   RESULTS.md and figures/ beside NOTE.md in that folder

USAGE
    python HAT_peaisland_extension.py stage0                  # the 5 zeroBE runs
    python HAT_peaisland_extension.py check                   # base geometry, edgeBE:
                                                              # must reproduce the matrix row
    python HAT_peaisland_extension.py probe --member n115-asrun --step 1 \
        --override "1=32.2,115=12.0"                          # one Newton probe
    python HAT_peaisland_extension.py next --member n115-asrun  # the next probe, from
                                                              # the runs so far
    python HAT_peaisland_extension.py score                   # RESULTS.md + figures

Each run is the ordinary hindcast runner driven through the environment,
exactly as HAT_run_matrix.py drives the matrix; nothing here reimplements a
run. Runs are never overwritten (pass --overwrite to redo one).
```

Notes that were in the code:

```text
the detrended planform has no calibrated 90-domain run to compare
against, so its baseline is solved here alongside the extensions
```

```text
the regression check: base geometry, compressed, edgeBE -- must land
on the matrix row's numbers exactly
```

```text
the zeroBE stage-0 run sits under zeroBE/, the probes under edgeBE/;
the solve script reads each run's preset off the index, so every run
is passed with its own tag and nothing else
```

```text
Per end: the solve script prints CONVERGED once the residual is inside
tolerance, and "no secant" once two probes imposed the same value --
which only happens after a converged end's step rounded to 0.0. Either
is done. Converged means every end is.
```

```text
An end the script gave no step for (converged, or no secant)
keeps its last imposed value: the runner needs a rate at every end.
```

```text
a probe row is "<run name> <imposed> <model> <residual>"; the
"next probe" line also starts with HAT_ and has no numbers
```

```text
A probe whose |value| exceeds this is not a boundary term any more; the
matrix values are tens of m/yr and Barrier3D's overwash router has been seen
to die silently under runaway progradation. Stop and say so instead.
```

<details><summary>Function notes (the original docstrings)</summary>

**`run_env()`**

```text
The environment one run reads, built the way HAT_run_matrix builds it:
every HAT_* variable named here, none inherited from the shell.
```

**`launch()`**

```text
One hindcast run, filed under experiments/<TAG>/<member>/ (stage 0)
or experiments/<TAG>/<member>/step<k>/ (a Newton probe).
```

**`solve_history()`**

```text
(run_name, preset, tag) of the member's stage-0 run and every probe
on disk, oldest first, read off the experiment folder.
```

**`_member_step()`**

```text
(member, step) from a run's tag: <TAG>/<member> is stage 0,
<TAG>/<member>/step<k> the k-th Newton probe. (None, None) otherwise.
```

**`mark_solved()`**

```text
A SOLVED file in the member folder naming the step that is the
solved run, for a reader who is not going to parse the index.
```

**`next_probe()`**

```text
The solve script over the member's runs so far, with the member's
geometry in the environment. Returns (override or None, converged,
text).
```

**`cmd_solve()`**

```text
Newton steps for one member until both ends converge or --max-steps
is reached: next probe, run it, repeat.
```

**`_target()`**

```text
The CoastSat target as the runner builds it: the surveyed GIS 1-90
table (what the interior score uses), or the extension's table over
GIS 0-115 (what an extended run's ends are solved against).
```

**`collect()`**

```text
Every run of the experiment plus the two matrix baselines, as rows:
member, step, preset, index row, and the per-domain LRR.
```

**`draw_compare()`**

```text
The figure to read first: the main approach (base geometry, edgeBE,
compressed planform, the matrix run) against the solved n115 extension
on GIS 1-90 only. One panel (Hannah, 2026-09-16: no difference or
residual panels, no southern extension).
```

**`draw()`**

```text
The whole extended reach: the alongshore LRR of the solved runs
against the target, a panel per offset mode, in house style.
```

</details>
