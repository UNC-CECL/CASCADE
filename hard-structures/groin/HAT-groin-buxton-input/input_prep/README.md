# input_prep - the scripts that make the groin study's inputs

Map of the whole input tree: `../README.md`.

```
HAT_target_shoreline_change.py    observed dune-line change and rate, D2-D12 -> ../groin_init/target/
island_offset/                    the 1967 island offset (its own README)
shoreline_position/               distance-to-datum check (its own README)
storms/                           the two storm builders (their own README)
```

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning and choices behind it are here. Moved out of the scripts on 2026-10-01, when they were brought in line with `scripts/STYLE.md` (code unchanged, proven with `style_equivalence_check.py`; path fixes listed per script).

### HAT_target_shoreline_change.py

Builds the target the groin module must reproduce: the observed dune-line
position change and rate per domain, D2-D12, between `CHANGE_FROM` (1967) and
`CHANGE_TO` (1997), from the raw ArcGIS offset files in
`data/hatteras_init/2-brie-offset/raw_offsets/`.

**Shared baseline.** `island_offset_hybrid_1967.py` references each year to its
own minimum, which is right for a single CASCADE input but destroys
cross-year comparability. This script keeps every year on the shared raw
`ORIG_LEN` reference, so positions can be differenced; `DISPLAY_REFERENCE`
re-references only the plots (`"none"`, `("year", Y)` or `("domain", D)`).
Differences and rates do not depend on it.

**Sign.** `ORIG_LEN` increases landward, like CASCADE's x_s, so a positive
`position_change_m` is erosion. With `FLIP_FOR_RATE` the reported
`rate_m_per_yr` is + seaward, matching the hindcast plots.

**Quantification.** `quantify()` reduces the curve to zone metrics (updrift and
downdrift mean rates, the downdrift peak, the erosion differential) and a
fold-and-sum split across the groin at `GROIN_BOUNDARY` = 5.5: each updrift
domain u is paired with its mirror p = 2*5.5 - u; A = (rate_u - rate_p)/2 is the
groin's antisymmetric signal, B = (rate_u + rate_p)/2 the background with the
groin removed. A is the background-free target M reproduces; a flat B (range
under 1.5 m/yr) supports the uniform-background dipole model.

**Figures.** The target curve (positions by year over the rate curve), the
fold-sum figure, a quantification dashboard, a trajectory panel (position
against year, updrift validation domains D6-D10 bold, the rest muted) and a GIF
holding each snapshot for a time proportional to the real gap to the next
(`GIF_SECONDS_PER_YEAR`, floor `GIF_MIN_HOLD_S`; skipped if pillow is missing).

**Path fix (2026-10-01).** `OUT_DIR` was the dead literal
`C:\Users\hanna\PycharmProjects\CASCADE\scripts\groin_module\hindcast_groin_test\groin_init\target`,
so the script failed on its first save. It now resolves from the repo root to
`hard-structures/groin/HAT-groin-buxton-input/groin_init/target/`, where its
`HAT_target_1967_1997_*` products already are.

From the script's original header:

```text
HAT_target_shoreline_change.py
==============================
Build the TARGET the groin module must reproduce: the observed dune-line
POSITION CHANGE (and rate) across years, per domain, from the raw ArcGIS offset
files. This is what you compare the CASCADE groin runs against.

Why this exists
---------------
The groin's job is to bend the modeled shoreline so it matches the historical
1967->1997 differential (holding updrift, eroding downdrift). To calibrate M you
need that historical signal as a per-domain curve. This script produces it.

Key difference from island_offset_hybrid_1967.py
------------------------------------------------
That pipeline references EACH year to its OWN minimum -- correct for a single
CASCADE input file, but it destroys cross-year comparability (every year gets a
different zero). To compare POSITIONS across years, all years must share ONE
baseline. This script:
  1. extracts per-domain mean ORIG_LEN (raw cross-shore position, m) per year,
  2. keeps them on a SHARED raw reference (no per-year re-zeroing),
  3. computes position change between any two years, and the rate (m/yr),
  4. optionally re-references the whole set to one year (e.g. 1967) or one
     domain, purely for display -- differences/rates are reference-invariant.

Sign convention
---------------
ORIG_LEN increases LANDWARD (larger = more retreat), same as CASCADE x_s. So a
POSITIVE position change = landward = EROSION. With FLIP_FOR_RATE=True the
reported RATE is flipped to + = seaward/accretion, matching your hindcast plots.

Output
------
  <OUT>_positions_by_year.csv   per-domain position each year (shared raw ref)
  <OUT>_change_<A>_<B>.csv      position change + rate between year A and B
  <OUT>_target_curve.png        the target rate curve (what the groin must match)
```

<details><summary>Function notes (the original docstrings)</summary>

**`per_domain_positions()`**

```text
Mean POSITION_COL per domain (raw, shared reference). Series indexed by domain.
```

**`build_position_table()`**

```text
DataFrame: index = domain (DOMAINS), columns = years, values = raw position (m).
```

**`apply_display_reference()`**

```text
Return a copy shifted per DISPLAY_REFERENCE (for the positions plot only).
```

**`compute_change()`**

```text
Per-domain position change and rate between two years. Reference-invariant.
```

**`plot_trajectories()`**

```text
Position vs YEAR, one line per domain -- the honest 'over time' view.
Years sit at their true spacing on the x-axis, so uneven gaps show correctly.
Updrift validation domains (D6-D10) are bold/saturated; downdrift (D2-D5,
not validated) and far-updrift edge (D11-D12) are muted.
```

**`make_gif()`**

```text
Proportional-duration GIF: each snapshot is held for a time proportional
to the real gap to the NEXT snapshot, so uneven year spacing is respected.
Falls back gracefully (prints a note) if pillow isn't available.
```

**`quantify()`**

```text
Reduce the per-domain target to defensible scalar metrics + a fold-and-sum
decomposition (A = groin signal, B = background) across the groin.

Returns (zone_summary_dict, foldsum_dataframe).
```

**`plot_foldsum()`**

```text
Two-line figure: A(y) = groin signal (what M reproduces), B(y) = background.
```

**`plot_quantification_dashboard()`**

```text
One-glance summary of the quantified target: per-domain rate bars colored
by zone, a downdrift-vs-updrift mean comparison with the differential ratio,
and the A/B fold-sum decomposition.
```

</details>

<details><summary>Notes that were comments in the code</summary>

- Anchored 2026-09-14: absolute into a home directory, or into a tree renamed since. Rule 5 of ORGANIZATION.md.
- Display reference (does NOT affect change/rate, only the positions plot): "none"        -> keep shared raw ORIG_LEN ("year", Y)   -> subtract year Y's positions (so Y becomes the zero line) ("domain", D) -> subtract each year's value at domain D
- Distinct saturated colors for the updrift validation zone (D6-D10); muted grays for downdrift (not validated) and the D11-D12 edge.
- Build the frame schedule: repeat each year's frame in proportion to the gap to the next year (last year gets the same hold as the previous gap).

</details>
