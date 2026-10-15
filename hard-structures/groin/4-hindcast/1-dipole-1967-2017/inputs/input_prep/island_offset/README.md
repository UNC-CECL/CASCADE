# island_offset - the 1967 island offset for the groin domains

One script: `island_offset_hybrid_1967.py`, the subset version of the main
`scripts/input_prep/2-brie-offset/1-produce/island_offset_hybrid.py`, run for
1967 and domains D2-D12 (Cape Point to Buxton). Map: `../../README.md`.

## The scripts in detail

Each script's header says what it does and how to run it; the reasoning and choices behind it are here. Moved out of the scripts on 2026-10-01, when they were brought in line with `scripts/STYLE.md` (code unchanged, proven with `style_equivalence_check.py`; path fixes listed per script).

### island_offset_hybrid_1967.py

Reads one raw dune-baseline intersection CSV, computes the relative dune offset
per domain (metres, baseline = the minimum of the domains in the run), and pads
it for CASCADE: 11 real domains + 15 buffers per side = 41.

**Per-domain value.** Every transect present in a domain is used, the first
record within each transect; the domain value is the mean, so thin edge domains
simply average fewer points (flagged `<- thin` under three transects). A
missing domain stops the script: a short array would misalign every CASCADE
domain downstream.

**Hybrid buffers.** The innermost `SLOPE_BUFFER_DOMAINS` (5) on each side follow
the local coastline slope, fitted on `EXTRAP_FIT_DOMAINS` (5) real domains at
that edge and anchored exactly on the edge domain, so the step into the buffer
is one slope step |m|, not zero: the buffer continues the trend rather than
repeating the edge value. The remaining outer buffers are one linear bridge
between the two slope tails (a single linspace including both tails, split
into the two outer blocks), which keeps the wrap-around array continuous.
What happens in the outer buffer does not affect the real domains.

Index map (n_pad per side, n_slope slope, n_bridge bridge): idx 0 is the
outermost left buffer; [0, n_bridge) left bridge; [n_bridge, n_pad) left slope
(n_pad-1 touches D_first); [n_pad, n_pad+n_real) real; then the right slope
(touching D_last) and the right bridge to idx target_length-1.

**`CLIP_BUFFERS_AT_ZERO = False`.** Real offsets are >= 0 by construction, but
extrapolating outward from the edge that holds the run minimum goes negative at
once and would be flattened into an artificial zero-gradient shelf. For the
1967 D2-D12 subset the minimum is D12, the northern edge, so the north buffer
is allowed to go negative and continue the real D8-D12 trend.

**`RE_REFERENCE_PADDED = True`.** The padded array is shifted so its minimum is
0. The offsets are relative, so this uniform translation of real and buffer
values changes no alongshore gradient; it only acts when the padded array goes
negative. Off, CASCADE gets raw negatives and the real values match the
unpadded file exactly.

**Output location.** Until 2026-10-01 `OUTPUT_DIR` named
`hard-structures/groin/HAT-buxton-hindcast-groin-test/groin_init`, a folder that
no longer exists (the script would have created it). It now names
`../../groin_init/island_offset/`, where the products in use already were.

From the script's original header:

```text
Hatteras CASCADE Dune Offset Pipeline — subset-capable
======================================================

This script:
1. Reads a single raw dune–baseline intersection CSV.
2. Calculates the relative dune raw_offset per domain (meters, baseline = minimum
   of the domains actually included in the run).
3. Pads the result for CASCADE using a hybrid buffer strategy:
     - inner buffer domains follow the local coastline slope extrapolated
       outward from each real edge (anchored exactly at the edge domain),
     - outer buffer domains are a linear bridge connecting the two slope
       tails, keeping the wrap-around array continuous.
4. Saves a diagnostic figure of the full padded raw_offset profile.

Generalized for ANY contiguous domain subset via START_DOMAIN / END_DOMAIN.
Current configuration: 1967, domains 2–12 (Cape Point → Buxton).

Author: Hannah A. Henry (extrapolation buffer version, subset-capable)
```

<details><summary>Function notes (the original docstrings)</summary>

**`calculate_relative_offset()`**

```text
Compute mean relative dune raw_offset per domain from raw CSV.

For each domain in `grids`, every transect present is used; within a
transect the first record is taken. The domain value is the mean across
whatever transects exist — domains with fewer transects (e.g. partial
edge domains) are handled naturally and simply average fewer points.
```

**`pad_for_cascade()`**

```text
Pad raw_offset array for CASCADE using a hybrid buffer strategy.

Let D_first / D_last be the first and last REAL domains of the run
(south → north order), and n_pad = padding_zeros.

  Left buffer (n_pad domains, outermost → D_first):
    - Innermost `slope_buffer_domains` (closest to D_first): follow the
      local coastline slope, anchored exactly at D_first so there is no
      gap at the boundary.
    - Remaining outer domains: linear bridge toward the right buffer.

  Right buffer (n_pad domains, D_last → outermost):
    - Innermost `slope_buffer_domains` (closest to D_last): follow the
      local coastline slope, anchored exactly at D_last.
    - Remaining outer domains: the same bridge, approached from the
      other end.

This guarantees:
  - No discontinuity at either real-domain boundary.
  - The buffer interior is connected (no jumps anywhere).
  - Behaviour in the outer buffer is irrelevant to the real simulation.

Index map (n_real real domains, n_slope slope, n_bridge bridge per side):
    idx 0                     = outermost left buffer
    idx [0, n_bridge)         = left bridge
    idx [n_bridge, n_pad)     = left slope   (idx n_pad-1 touches D_first)
    idx [n_pad, n_pad+n_real) = real domains
    idx [n_pad+n_real, +n_slope)          = right slope (touches D_last)
    idx [n_pad+n_real+n_slope, +n_bridge) = right bridge
    idx target_length-1       = outermost right buffer

Returns
-------
padded : pd.DataFrame
diag   : dict  (arrays for diagnostic plot)
```

**`clip_zones_to_run()`**

```text
Filter/clip full-island zone definitions to the run's domain window.
```

**`plot_buffer_diagnostic()`**

```text
Save a diagnostic figure of the full padded raw_offset profile.

Layout:
  - Full profile across all padded indices
  - Buffer zones shaded; slope vs bridge segments distinguished
  - Real-domain zone annotations along the top
  - Inset zoom panels for left and right buffer transitions
```

</details>

<details><summary>Notes that were comments in the code</summary>

- Anchored 2026-09-14: this named a home directory, or a tree renamed since. Rule 5 of ORGANIZATION.md.
- Number of real domains from each edge used to fit the local extrapolation trend. Must be <= N_REAL (clamped automatically with a warning).
- Number of buffer domains on each side that follow the local coastline slope. The remaining (PADDING_ZEROS - SLOPE_BUFFER_DOMAINS) buffer domains on each side are filled by a linear bridge connecting the two slope tails.
- Clip buffer offsets at zero. Relative offsets are >= 0 by construction for real domains (baseline = minimum), but an outward extrapolation from the edge that holds the run minimum will immediately go negative and be flattened into an artificial zero-gradient shelf. For the D2-D12/1967 subset the minimum sits on D12 (the northern edge), so this is set False: the north buffer is allowed to go negative and continue the real D8-D12 trend.
- Re-reference the PADDED array to its own minimum so the file CASCADE reads is non-negative. The offsets are relative, so this is a uniform translation of the whole array (real + buffers) and does not change any alongshore gradient. Only takes effect if the padded array actually goes negative. Set False to ship raw negatives if CASCADE tolerates them — in that case the real-domain values stay exactly as they appear in the unpadded file.
- Community zone annotations for the diagnostic figure. Defined on the FULL island (real domain numbers, 1-indexed). Zones are automatically clipped/filtered to the START_DOMAIN–END_DOMAIN window.
- Use whatever transects are actually present (robust to gaps and to partial edge domains with only one or two transects).
- 2. Left slope segment (innermost, closest to D_first) Steps: -1 (adjacent to D_first) … -slope_buffer_domains
- 4. Linear bridge ONE linspace from the left slope TAIL to the right slope TAIL, including both tails as endpoints, then split into the two outer buffer blocks. Left/right slope arrays are ordered inner-first, so [-1] of each is the tail the bridge anchors to.
- 5. Assemble full left and right buffer arrays left_full  : idx 0 = outermost, idx -1 = adjacent to D_first right_full : idx 0 = adjacent to D_last, idx -1 = outermost
- 6b. Re-reference to the padded minimum The offsets are relative, so adding a constant to every element (real AND buffer) is a uniform translation that leaves every alongshore gradient untouched. This is what makes it safe to let the edge extrapolation go negative: the trend at the boundary is preserved, and the file CASCADE reads is still >=0.
- Boundary step = the jump from the edge real domain into its adjacent buffer. Because the slope segment is anchored ON the edge domain, this should equal exactly one slope step (|m|), not zero — i.e. the buffer continues the local trend rather than repeating the edge value.
- Hard stop if any requested domain is missing — a silently short array would misalign every downstream CASCADE domain.

</details>
