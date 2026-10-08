# 2-solver-audit — the dipole in BRIE's alongshore solve alone

What BRIE's alongshore solve does with the groin dipole, measured without
running CASCADE. The emulator reproduces BRIE's implicit shoreline-diffusion
step (its own `coast_diff` table and sparse indices) plus the +/-M dipole that
`cascade.groin.GroinCallback` injects into `x_s_dt`, and nothing else, so a
1,000-year sweep costs seconds and anything the full model does differently is
Barrier3D's cross-shore feedback.

| file | what it is |
|---|---|
| `HAT_groin_solver_audit.py` | the six audits; also the emulator the later studies import |
| `solver_audit_<audit>.csv` | its tables (2026-09-11): `diffusivity`, `amplitude_runaway`, `sink_fraction`, `domain_count`, `groin_field`, `chosen_rig` |
| `../3-real-planform/` | the follow-up under option A waves, on the real planform (its own README) |

Run it from anywhere; the tables go beside the script unless `--outdir` says
otherwise:

    python HAT_groin_solver_audit.py [--years 200] [--outdir DIR]

The scripts in `../3-real-planform/` and `../1-straight-coast/` import
`brie_diffusivity`, `DY_M` and `DT_YR` from this script, so a change to those
reaches them.

## The scripts in detail

The header of each script says what it does and how to run it; the reasoning
behind it is here. Moved out of the script on 2026-10-01, when the groin study
was brought in line with `scripts/STYLE.md`.

### HAT_groin_solver_audit.py

**What it is.** An emulator, not a reimplementation for production use. It
mirrors the `if self._ast_model_on:` block of `brie.brie.Brie.update()`: the
same `coast_diff` table (taken from a real `Brie` instance, never recomputed),
the same row-scaled matrix assembly from `_di`/`_dj`, the same periodic wrap,
and the same `np.maximum(0, ...)` clip on the diffusion number. If BRIE
changes, it has to be re-checked against it.

**The two fixed constants.** `DY_M = 500` and `DT_YR = 1` are fixed by BRIE in
the CASCADE coupler, whose comment says "do not change", so they are constants
here rather than arguments.

**The two wave climates.** `DEFAULT_CLIMATE` is the CASCADE / Barrier3D
default; `HATTERAS_CLIMATE` is the Hatteras hindcast setting of the time. They
are not interchangeable for this test: the diffusivity table shows they differ
by a factor of four in the restoring rate, and so in the M a given fillet
costs.

**The solve.** `solve_reach` starts from a perfectly straight shoreline at
x_s = 0, so every metre of structure in the answer came from a dipole. The
reach is periodic, so the domain count is a circumference. The net source a
field injects per year is seaward-negative; a field with every f = 1 is volume
neutral and must leave the reach mean alone. Each year the diffusion number is
taken from BRIE's forward-difference shoreline angle and clipped at zero, as
BRIE does; the clip is what lets the scheme stop diffusing rather than go
unstable.

**The audits**, in the order `main()` runs them (the header said "five"; the
sixth, the chosen rig, was added after):

1. diffusivity and shutdown angle for each climate;
2. fillet against M, and the runaway boundary;
3. fillet and volume closure against f;
4. domain-count convergence;
5. groin fields: stacked, spaced and opposed structures;
6. the chosen rig: buffer size and the runaway boundary across wave heights.

From the script's original header, as it stood before 2026-10-01, word for word:

```text
Stage 0 of the groin module test: audit BRIE's solve without running CASCADE.

WHY THIS EXISTS
    Every claim about what a groin "does" in CASCADE is a claim about what
    BRIE's implicit alongshore-diffusion solve does with the +/-M dipole that
    `cascade.groin.GroinCallback` injects into `x_s_dt`. That solve is cheap:
    one sparse tridiagonal-plus-corners system per year. Reproducing it here
    WITHOUT Barrier3D means a 1,000-year, many-axis experiment costs seconds
    instead of days, and -- more importantly -- it separates the groin's
    alongshore behaviour from the cross-shore behaviour Barrier3D adds on top.
    Anything the full rig shows that this does not is a Barrier3D feedback.

    This is an EMULATOR, not a reimplementation for production use. It mirrors
    `brie.brie.Brie.update()` (the `if self._ast_model_on:` block, lines
    ~1290-1330) exactly: the same `coast_diff` table, taken from a real Brie
    instance rather than recomputed, the same row-scaled matrix assembly from
    `_di`/`_dj`, the same periodic wrap, the same `np.maximum(0, ...)` clip on
    the diffusion number. If BRIE changes, this must be re-checked against it.

WHAT IT MEASURES
    1. DIFFUSIVITY AND SHUTDOWN ANGLE -- BRIE's diffusion number is clipped at
                                         zero, and the wave-climate-averaged
                                         diffusivity goes NEGATIVE past a
                                         critical shoreline angle. Past that
                                         angle a cell stops exchanging sand
                                         with its neighbours and the dipole
                                         accumulates without limit.
    2. FILLET AMPLITUDE vs M          -- is the response linear in M, as
                                         `groin.predict_fillet` asserts, and
                                         where is the runaway boundary?
    3. VOLUME CLOSURE vs f            -- the matrix is row-scaled by a
                                         per-domain diffusion number, so its
                                         column sums are not unity and the
                                         scheme is NOT exactly conservative
                                         once the shoreline is not straight.
                                         This reports the spurious mean drift
                                         against the drift the injected volume
                                         actually implies.
    4. DOMAIN-COUNT CONVERGENCE       -- the solve is PERIODIC in the
                                         alongshore, so a short reach wraps the
                                         fillet into its own downdrift notch.
                                         This reports how few domains the
                                         fillet tolerates and how badly the
                                         mean drift is inflated by a short one.
    5. GROIN FIELDS                   -- several dipoles at a given spacing.
                                         Two questions: do they superpose (a
                                         field of N reads as one groin of N*M),
                                         and at what spacing does each
                                         structure still hold its own fillet
                                         rather than the field holding one?

HOW TO READ "RUNAWAY"
    A cell whose diffusion number has been clipped to zero keeps receiving its
    share of the dipole and has no way to pass it on, so the reach translates at
    roughly M metres per year indefinitely. In the full model this presents as a
    barrier that migrates absurdly or drowns; in BRIE alone it eventually
    presents as `IndexError: index 180 is out of bounds` from brie.py, because
    the diffusivity lookup indexes `coast_diff` (length 180) with `90 - theta`
    clipped to 180 rather than 179. That needs a 57 km offset across one 500 m
    cell, which a runaway reaches in roughly 1,000 years at M = 60. It is the
    runaway surfacing, not a separate coding mistake.

Author: Hannah A. Henry, UNC CECL
```

<details><summary>Function notes (the original docstrings)</summary>

**`brie_diffusivity()`**

```text
Return (coast_diff, di, dj) from a real Brie instance.

The diffusivity table and the sparse index arrays are taken from BRIE
rather than recomputed, so this audit cannot quietly drift away from the
model it is auditing.
```

**`shutdown_angle_deg()`**

```text
Shoreline angles bounding the band where BRIE's diffusivity is positive.

Outside that band the diffusion number is clipped to zero and the cell
decouples from its neighbours. Returns (negative_cutoff, positive_cutoff)
in degrees.
```

**`solve_reach()`**

```text
Integrate BRIE's shoreline diffusion with groin dipoles, nothing else.

Starts from a perfectly straight shoreline at x_s = 0, so every metre of
structure in the answer came from a dipole.

Parameters
----------
ny : int
    Domain count. The reach is PERIODIC, so this is a circumference.
groins : list of (updrift, downdrift, M, f)
    One tuple per structure. `updrift` and `downdrift` are domain indices
    and must be adjacent, matching `GroinCallback`'s own check. Sign
    follows `cascade.groin`: the updrift cell gets -M (seaward advance),
    the downdrift cell gets +M*f (landward retreat).
years : int
    Run length.
climate : dict
    Keys Hs, Tp, asym, ahf.
record_years : iterable of int
    Years whose full shoreline profile to keep.

Returns
-------
frame : DataFrame
    One row per year: fillet at each structure, reach mean, the mean a
    conservative scheme would give, the closure error, and the minimum
    diffusion number anywhere.
shutdown_year : int or None
    First year the diffusion number hit zero anywhere.
profiles : dict
    Shoreline profiles for `record_years`.
```

**`one_groin()`**

```text
A single structure at the middle of the reach, drift from high index.
```

**`audit_diffusivity()`**

```text
Diffusivity and shutdown angle for each wave climate.
```

**`audit_amplitude_and_runaway()`**

```text
Fillet against M across wave heights; flag the runaway boundary.
```

**`audit_sink_fraction()`**

```text
Fillet and volume closure against f, at fixed M.
```

**`audit_domain_count()`**

```text
How the fillet and the mean drift depend on reach length.
```

**`audit_groin_field()`**

```text
Several structures: do they superpose, and at what spacing do they merge?

Two separate questions, deliberately kept apart.

STACKED. Several dipoles on the SAME pair of domains is the case
the dipole fit notes (`../../3-hindcast/1-dipole-1967-2017/dipole_fit_notes.md`) call "four groins, one dipole, deliberately" -- the real Buxton
field fits inside one 500 m cell. In a LINEAR solve N stacked dipoles of
amplitude M are exactly one dipole of amplitude N*M. The solve is not
linear, because the diffusion number depends on the shoreline angle, so
this measures how far from N*M the answer actually lands.

SPACED. Dipoles every `spacing` domains. A field whose structures are far
apart holds one fillet each; a field whose structures are close holds one
fillet for the whole field, with the interior ones doing nothing. The
interior-to-end fillet ratio says which regime a spacing is in.
```

**`audit_chosen_rig()`**

```text
The rig as designed: 20 working domains, `n_buffer` buffers per side.

Two things this has to establish, because the rig's width rests on them.

WHETHER THE BUFFER IS BIG ENOUGH. It is NOT an absorbing boundary -- BRIE's
solve is periodic, so a buffer separates the groin from the wrap, it does
not soak anything up. The diffusive reach of the dipole grows as
dy*sqrt(2*r_ipl*t) with no dependence on M, so the buffer is exceeded after
a time that depends only on the wave climate. Past that point the fillet's
own tail has come round the back. Section 4 of this audit says that costs
about 3% of the fillet, and rather more of the reach mean and of any attempt
to measure the alongshore extent.

WHERE THE RUNAWAY BOUNDARY FALLS ACROSS THE WAVE-HEIGHT RANGE. The fillet is
bought by M / r_ipl, so sweeping wave height IS sweeping the cost of M. This
prints the grid the CASCADE sweep should cover and flags which cells are
past the shutdown.
```

</details>

Notes that were comments in the code:

- BRIE fixes `DY_M` and `DT_YR` in the CASCADE coupler, whose comment says "do
  not change", so they are constants here rather than arguments.
- Wave climate: the first is the CASCADE / Barrier3D default, the second the
  Hatteras hindcast setting. They are NOT interchangeable for this test: see
  the diffusivity table, where they differ by a factor of four in the restoring
  rate and therefore in the M a given fillet costs.
- The net source the field injects per year is seaward-negative. A field with
  every f = 1 is volume neutral and must leave the reach mean alone.
- BRIE's forward-difference shoreline angle, and the diffusion number it
  selects. The clip at zero is BRIE's, and it is what lets the scheme stop
  diffusing rather than go unstable.
- The section banner read "The five audits"; there are six since the chosen
  rig was added.
