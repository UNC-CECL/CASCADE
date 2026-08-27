"""
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
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

YEAR = 1967
RAW_FILE = r"/data/hatteras_init/2-brie-offset/raw_offsets/1967_duneline_offset_raw.csv"

OUTPUT_DIR      = r"/scripts/groin_module_noBE/HAT-hindcast-groin-test/groin_init"
OUTPUT_BASENAME = "Island_Dune_Offsets_1967_D2_D12"

# -------------------------------------------------------------------------
# Domain range for this run (contiguous subset of the 1–90 island grid)
# -------------------------------------------------------------------------
START_DOMAIN = 2
END_DOMAIN   = 12
B3D_GRIDS    = list(range(START_DOMAIN, END_DOMAIN + 1))

N_REAL = END_DOMAIN - START_DOMAIN + 1                     # 11

PADDING_ZEROS = 15                                          # buffers per side
TARGET_LENGTH = N_REAL + 2 * PADDING_ZEROS                  # 11 + 30 = 41

# Number of real domains from each edge used to fit the local extrapolation
# trend. Must be <= N_REAL (clamped automatically with a warning).
EXTRAP_FIT_DOMAINS = 5

# Number of buffer domains on each side that follow the local coastline slope.
# The remaining (PADDING_ZEROS - SLOPE_BUFFER_DOMAINS) buffer domains on each
# side are filled by a linear bridge connecting the two slope tails.
SLOPE_BUFFER_DOMAINS = 5

# Clip buffer offsets at zero. Relative offsets are >= 0 by construction for
# real domains (baseline = minimum), but an outward extrapolation from the
# edge that holds the run minimum will immediately go negative and be
# flattened into an artificial zero-gradient shelf. For the D2-D12/1967 subset
# the minimum sits on D12 (the northern edge), so this is set False: the north
# buffer is allowed to go negative and continue the real D8-D12 trend.
CLIP_BUFFERS_AT_ZERO = False

# Re-reference the PADDED array to its own minimum so the file CASCADE reads
# is non-negative. The offsets are relative, so this is a uniform translation
# of the whole array (real + buffers) and does not change any alongshore
# gradient. Only takes effect if the padded array actually goes negative.
# Set False to ship raw negatives if CASCADE tolerates them — in that case the
# real-domain values stay exactly as they appear in the unpadded file.
RE_REFERENCE_PADDED = True

COL_MAP = {
    "Domain_ID": "domain_id",
    "Distance":  "ORIG_LEN",
    "Transect":  "LineID",
}

# Community zone annotations for the diagnostic figure.
# Defined on the FULL island (real domain numbers, 1-indexed). Zones are
# automatically clipped/filtered to the START_DOMAIN–END_DOMAIN window.
COMMUNITY_ZONES = [
    (1,   6,  "Cape Point"),
    (7,   8,  "Buxton"),
    (9,  20,  "Buxton–Avon"),
    (21, 31,  "Avon"),
    (32, 67,  "Avon–Tri-Village"),
    (68, 83,  "Tri-Village"),
    (84, 90,  "Pea Island NWR"),
]

# =============================================================================
# 2. FUNCTIONS
# =============================================================================

def calculate_relative_offset(file_path, year, col_map, grids):
    """
    Compute mean relative dune raw_offset per domain from raw CSV.

    For each domain in `grids`, every transect present is used; within a
    transect the first record is taken. The domain value is the mean across
    whatever transects exist — domains with fewer transects (e.g. partial
    edge domains) are handled naturally and simply average fewer points.
    """
    print(f"\n--- Processing {year} ---")
    print(f"Input file: {file_path}")

    try:
        # utf-8-sig strips the BOM ArcGIS writes on the first column header.
        raw_df = pd.read_csv(file_path, encoding="utf-8-sig")
    except FileNotFoundError:
        print(f"ERROR: File not found: {file_path}")
        return None

    for key, col in col_map.items():
        if col not in raw_df.columns:
            print(f"ERROR: Column '{col}' not found.")
            return None

    data_df = pd.DataFrame({
        "B3D_Grid": raw_df[col_map["Domain_ID"]],
        "Distance": raw_df[col_map["Distance"]],
        "Transect": raw_df[col_map["Transect"]],
    })

    mean_distances = []
    seen_domains   = []
    n_transects    = []

    for grid_id in grids:
        subset = data_df[data_df["B3D_Grid"] == int(grid_id)]
        if subset.empty:
            print(f"  Warning: No data for domain {grid_id}.")
            continue

        # Use whatever transects are actually present (robust to gaps and to
        # partial edge domains with only one or two transects).
        transect_ids = sorted(subset["Transect"].dropna().unique())
        if not transect_ids:
            print(f"  Warning: No valid transect IDs for domain {grid_id}. Skipping.")
            continue

        distances = []
        for t_id in transect_ids:
            t_vals = subset[subset["Transect"] == t_id]
            if not t_vals.empty:
                distances.append(t_vals.iloc[0]["Distance"])

        if not distances:
            continue

        mean_distances.append(float(np.mean(distances)))
        seen_domains.append(grid_id)
        n_transects.append(len(distances))

    if not mean_distances:
        print(f"FATAL: No valid domains processed for year {year}.")
        return None

    baseline         = min(mean_distances)
    relative_offsets = np.subtract(mean_distances, baseline)

    print(f"  {len(mean_distances)} domains processed "
          f"(requested {len(grids)}: D{grids[0]}–D{grids[-1]}).")
    print(f"  Baseline distance = {baseline:.3f} m "
          f"(min mean, at domain {seen_domains[int(np.argmin(mean_distances))]}).")

    # Per-domain transect coverage — flags thin edge domains explicitly.
    print("\n  Domain  n_transects   mean_dist (m)   rel_offset (m)")
    for d, n, md, ro in zip(seen_domains, n_transects, mean_distances, relative_offsets):
        flag = "  <- thin" if n < 3 else ""
        print(f"  {d:>6}  {n:>11}   {md:>13.2f}   {ro:>14.2f}{flag}")

    return pd.DataFrame({"Domain_ID": seen_domains, str(year): relative_offsets})


def pad_for_cascade(df, padding_zeros, target_length,
                    extrap_fit_domains=5, slope_buffer_domains=5,
                    clip_at_zero=True, re_reference=True, start_domain=1):
    """
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
    """
    data_columns = list(df.columns)
    col          = data_columns[0]
    values       = df[col].to_numpy(dtype=float)   # shape (n_real,)
    n_real       = len(values)

    # ------------------------------------------------------------------ #
    # 0. Guards — matter for short subset runs                            #
    # ------------------------------------------------------------------ #
    if slope_buffer_domains >= padding_zeros:
        raise ValueError(
            f"slope_buffer_domains ({slope_buffer_domains}) must be < "
            f"padding_zeros ({padding_zeros}); no room left for the bridge."
        )

    if extrap_fit_domains > n_real:
        print(f"  Warning: EXTRAP_FIT_DOMAINS ({extrap_fit_domains}) > number of "
              f"real domains ({n_real}). Clamping to {n_real}.")
        extrap_fit_domains = n_real

    if extrap_fit_domains < 2:
        raise ValueError("extrap_fit_domains must be >= 2 to fit a slope.")

    if extrap_fit_domains > n_real // 2:
        print(f"  Note: EXTRAP_FIT_DOMAINS ({extrap_fit_domains}) covers more than "
              f"half of the {n_real} real domains — the two edge slope fits "
              f"overlap and are not independent.")

    d_first = start_domain
    d_last  = start_domain + n_real - 1

    # ------------------------------------------------------------------ #
    # 1. Estimate local slopes at each edge                               #
    # ------------------------------------------------------------------ #
    x_fit = np.arange(extrap_fit_domains, dtype=float)

    # Left edge slope (fit the first `extrap_fit_domains`, anchor at values[0])
    m_left, _ = np.polyfit(x_fit, values[:extrap_fit_domains], 1)

    # Right edge slope (fit the last `extrap_fit_domains`, anchor at values[-1])
    m_right, _ = np.polyfit(x_fit, values[-extrap_fit_domains:], 1)

    # ------------------------------------------------------------------ #
    # 2. Left slope segment (innermost, closest to D_first)               #
    #    Steps: -1 (adjacent to D_first) … -slope_buffer_domains          #
    # ------------------------------------------------------------------ #
    left_slope_steps  = np.arange(-1, -(slope_buffer_domains + 1), -1, dtype=float)
    left_slope_raw    = values[0] + m_left * left_slope_steps
    left_slope_values = np.clip(left_slope_raw, 0.0, None) if clip_at_zero else left_slope_raw

    # ------------------------------------------------------------------ #
    # 3. Right slope segment (innermost, closest to D_last)               #
    # ------------------------------------------------------------------ #
    right_slope_steps  = np.arange(1, slope_buffer_domains + 1, dtype=float)
    right_slope_raw    = values[-1] + m_right * right_slope_steps
    right_slope_values = np.clip(right_slope_raw, 0.0, None) if clip_at_zero else right_slope_raw

    n_clip_left  = int(np.sum(left_slope_raw  < 0.0)) if clip_at_zero else 0
    n_clip_right = int(np.sum(right_slope_raw < 0.0)) if clip_at_zero else 0

    # ------------------------------------------------------------------ #
    # 4. Linear bridge                                                    #
    #    ONE linspace from the left slope TAIL to the right slope TAIL,   #
    #    including both tails as endpoints, then split into the two outer #
    #    buffer blocks. Left/right slope arrays are ordered inner-first,  #
    #    so [-1] of each is the tail the bridge anchors to.               #
    # ------------------------------------------------------------------ #
    n_bridge_each = padding_zeros - slope_buffer_domains

    full_bridge = np.linspace(
        left_slope_values[-1],      # L: left slope tail
        right_slope_values[-1],     # R: right slope tail
        n_bridge_each * 2 + 1,      # includes both endpoints
    )

    # Left bridge:  [fb[n-1], ..., fb[0]]  — outermost to innermost
    left_bridge_values  = full_bridge[n_bridge_each - 1::-1]

    # Right bridge: [fb[2n], ..., fb[n+1]] — innermost to outermost
    right_bridge_values = full_bridge[2 * n_bridge_each:n_bridge_each:-1]

    # ------------------------------------------------------------------ #
    # 5. Assemble full left and right buffer arrays                       #
    #  left_full  : idx 0 = outermost, idx -1 = adjacent to D_first       #
    #  right_full : idx 0 = adjacent to D_last, idx -1 = outermost        #
    # ------------------------------------------------------------------ #
    left_full  = np.concatenate([left_bridge_values, left_slope_values[::-1]])
    right_full = np.concatenate([right_slope_values, right_bridge_values])

    left_block  = pd.DataFrame({col: left_full})
    right_block = pd.DataFrame({col: right_full})

    # ------------------------------------------------------------------ #
    # 6. Assemble and validate                                            #
    # ------------------------------------------------------------------ #
    padded = pd.concat(
        [left_block, df[[col]], right_block],
        ignore_index=True,
    )

    # ------------------------------------------------------------------ #
    # 6b. Re-reference to the padded minimum                              #
    #     The offsets are relative, so adding a constant to every element #
    #     (real AND buffer) is a uniform translation that leaves every    #
    #     alongshore gradient untouched. This is what makes it safe to    #
    #     let the edge extrapolation go negative: the trend at the        #
    #     boundary is preserved, and the file CASCADE reads is still >=0. #
    # ------------------------------------------------------------------ #
    shift = 0.0
    if re_reference and padded[col].min() < 0.0:
        shift = float(-padded[col].min())
        padded[col] = padded[col] + shift
        left_full  = left_full  + shift
        right_full = right_full + shift
        values     = values     + shift

    print(f"\nPadding summary (hybrid slope + bridge):")
    print(f"  Real domains           : D{d_first}–D{d_last} ({n_real} domains)")
    print(f"  Slope domains per side : {slope_buffer_domains}  "
          f"(slope fit on {extrap_fit_domains} real domains)")
    print(f"  Bridge domains per side: {n_bridge_each}")
    print(f"  Edge slopes            : m_left={m_left:+.2f} m/domain  "
          f"m_right={m_right:+.2f} m/domain")
    # Boundary step = the jump from the edge real domain into its adjacent
    # buffer. Because the slope segment is anchored ON the edge domain, this
    # should equal exactly one slope step (|m|), not zero — i.e. the buffer
    # continues the local trend rather than repeating the edge value.
    step_left  = abs(left_full[-1] - values[0])
    step_right = abs(right_full[0] - values[-1])
    print(f"  D{d_first} boundary step       : buffer={left_full[-1]:.2f}  "
          f"D{d_first}={values[0]:.2f}  step={step_left:.2f} m "
          f"(expected |m_left|={abs(m_left):.2f})")
    print(f"  D{d_last} boundary step      : buffer={right_full[0]:.2f}  "
          f"D{d_last}={values[-1]:.2f}  step={step_right:.2f} m "
          f"(expected |m_right|={abs(m_right):.2f}{', reduced by 0-clip' if n_clip_right else ''})")
    print(f"  Left  buffer range : {left_full.min():.2f} – {left_full.max():.2f} m")
    print(f"  Right buffer range : {right_full.min():.2f} – {right_full.max():.2f} m")
    print(f"  Padded length      : {len(padded)} (target: {target_length})")

    if n_clip_left or n_clip_right:
        print(f"  NOTE: clipped at 0 — left slope: {n_clip_left}/{slope_buffer_domains} "
              f"domains, right slope: {n_clip_right}/{slope_buffer_domains} domains. "
              f"The outward extrapolation went negative there (expected when an "
              f"edge domain holds the run minimum); those buffers are flat at 0.")

    if shift > 0.0:
        print(f"  RE-REFERENCED: padded array shifted by +{shift:.2f} m so the "
              f"minimum is 0 (uniform translation; no gradient changed).")
        print(f"    -> D{d_first}={values[0]:.2f} m, D{d_last}={values[-1]:.2f} m in the PADDED file.")
        print(f"    -> The unpadded files keep the original D{d_last}=0 reference.")
    elif padded[col].min() < 0.0:
        print(f"  NOTE: padded array contains negative offsets "
              f"(min={padded[col].min():.2f} m) and RE_REFERENCE_PADDED is off.")

    if len(padded) != target_length:
        raise ValueError(
            f"Length mismatch: got {len(padded)}, expected {target_length}. "
            f"Check START_DOMAIN/END_DOMAIN vs PADDING_ZEROS."
        )

    diag = {
        "left_values":          left_full,             # outermost→D_first
        "real_values":          values,                # D_first→D_last
        "right_values":         right_full,            # D_last→outermost
        "slope_buffer_domains": slope_buffer_domains,
        "extrap_fit_domains":   extrap_fit_domains,
        "m_left":               m_left,
        "m_right":              m_right,
        "start_domain":         d_first,
        "end_domain":           d_last,
        "shift":                shift,
    }

    return padded, diag


def clip_zones_to_run(community_zones, start_domain, end_domain):
    """Filter/clip full-island zone definitions to the run's domain window."""
    out = []
    for d_start, d_end, label in community_zones:
        s = max(d_start, start_domain)
        e = min(d_end, end_domain)
        if s <= e:
            out.append((s, e, label))
    return out


def plot_buffer_diagnostic(diag, year, padding_zeros, community_zones,
                           output_dir, output_basename):
    """
    Save a diagnostic figure of the full padded raw_offset profile.

    Layout:
      - Full profile across all padded indices
      - Buffer zones shaded; slope vs bridge segments distinguished
      - Real-domain zone annotations along the top
      - Inset zoom panels for left and right buffer transitions
    """
    lv = diag["left_values"]     # outermost first
    rv = diag["real_values"]     # D_first → D_last
    rb = diag["right_values"]    # D_last-adjacent first

    n_slope  = diag["slope_buffer_domains"]
    n_bridge = padding_zeros - n_slope
    n_real   = len(rv)
    d_first  = diag["start_domain"]
    d_last   = diag["end_domain"]
    shift    = diag.get("shift", 0.0)

    # Full padded array in plot order (index 0 = leftmost buffer)
    full = np.concatenate([lv, rv, rb])

    left_buf_x  = np.arange(padding_zeros)
    real_x      = np.arange(padding_zeros, padding_zeros + n_real)
    right_buf_x = np.arange(padding_zeros + n_real, len(full))

    r_start = padding_zeros + n_real

    left_bridge_x  = np.arange(0, n_bridge)
    left_slope_x   = np.arange(n_bridge, padding_zeros)
    right_slope_x  = np.arange(r_start, r_start + n_slope)
    right_bridge_x = np.arange(r_start + n_slope, r_start + padding_zeros)

    left_bridge_v  = full[left_bridge_x]
    left_slope_v   = full[left_slope_x]
    right_slope_v  = full[right_slope_x]
    right_bridge_v = full[right_bridge_x]

    # ------------------------------------------------------------------ #
    # Figure layout                                                        #
    # ------------------------------------------------------------------ #
    fig = plt.figure(figsize=(18, 7), facecolor="white")
    ref_note = (f"  |  re-referenced +{shift:.0f} m" if shift > 0 else "")
    fig.suptitle(
        f"Buffer Diagnostic — Dune Offset Profile  |  {year}  |  "
        f"Real: D{d_first}–D{d_last} ({n_real})  |  "
        f"Slope: ±{n_slope} domains  |  Bridge: {n_bridge} domains each side"
        f"{ref_note}",
        fontsize=13, fontweight="bold", color="#1a1a2e", y=0.98,
    )

    ax_main  = fig.add_axes([0.05, 0.18, 0.68, 0.68])
    ax_left  = fig.add_axes([0.76, 0.55, 0.21, 0.32])
    ax_right = fig.add_axes([0.76, 0.13, 0.21, 0.32])

    BUF_COLOR    = "#f4a582"   # salmon  — slope segment
    BRIDGE_COLOR = "#b2abd2"   # purple  — linear bridge
    REAL_COLOR   = "#2166ac"   # blue    — real domains
    ZOOM_COLOR   = "#d6604d"   # red     — inset

    # ---- Main plot ---------------------------------------------------- #
    ax_main.axvspan(left_buf_x[0]  - 0.5, left_buf_x[-1]  + 0.5,
                    color=BUF_COLOR, alpha=0.10, zorder=0)
    ax_main.axvspan(right_buf_x[0] - 0.5, right_buf_x[-1] + 0.5,
                    color=BUF_COLOR, alpha=0.10, zorder=0)

    for xv in [padding_zeros - 0.5, padding_zeros + n_real - 0.5]:
        ax_main.axvline(xv, color="#444444", lw=1.0, ls="--", zorder=2)

    # Bridge segments (outer, purple dashed)
    ax_main.plot(left_bridge_x,  left_bridge_v,  "s--", color=BRIDGE_COLOR,
                 ms=4, lw=1.5, zorder=3,
                 label=f"Bridge ({n_bridge} domains/side, linear connect)")
    ax_main.plot(right_bridge_x, right_bridge_v, "s--", color=BRIDGE_COLOR,
                 ms=4, lw=1.5, zorder=3)

    # Slope segments (inner, salmon solid)
    ax_main.plot(left_slope_x,  left_slope_v,  "o-", color=BUF_COLOR,
                 ms=6, lw=2.0, zorder=4,
                 label=f"Slope (±{n_slope} domains, local coastline trend)")
    ax_main.plot(right_slope_x, right_slope_v, "o-", color=BUF_COLOR,
                 ms=6, lw=2.0, zorder=4)

    # Real domains
    ax_main.plot(real_x, rv, "o-", color=REAL_COLOR,
                 ms=4, lw=1.8, zorder=5,
                 label=f"Real domains (D{d_first}–D{d_last})")

    # Community zone annotations along the top
    span_v = full.max() - full.min()
    ax_main.set_ylim(bottom=min(-0.5, full.min() - 0.04 * span_v))
    if full.min() < 0:
        ax_main.axhline(0.0, color="#999999", lw=0.8, ls=":", zorder=1)
    _, ymax = ax_main.get_ylim()
    zone_colors = ["#e8e8f0", "#d8d8ec"] * 10

    zones = clip_zones_to_run(community_zones, d_first, d_last)
    for zi, (d_start, d_end, label) in enumerate(zones):
        px_start = (d_start - d_first) + padding_zeros
        px_end   = (d_end   - d_first) + padding_zeros
        ax_main.axvspan(px_start - 0.5, px_end + 0.5,
                        color=zone_colors[zi], alpha=0.35, zorder=1)
        span = px_end - px_start + 1
        ax_main.text((px_start + px_end) / 2, ymax * 0.97, label,
                     ha="center", va="top", fontsize=7.5,
                     color="#333355", rotation=90 if span < len(label) * 0.6 else 0,
                     clip_on=True)

    ax_main.set_xlabel(
        f"Padded domain index  (0–{len(full) - 1};  "
        f"real = {padding_zeros}–{padding_zeros + n_real - 1})", fontsize=11)
    ax_main.set_ylabel("Dune raw_offset (m)", fontsize=11)
    ax_main.set_xlim(-0.5, len(full) - 0.5)
    ax_main.legend(fontsize=9, loc="upper center", framealpha=0.9)

    # Top axis: real domain numbers, tick spacing scaled to run length
    ax2 = ax_main.twiny()
    ax2.set_xlim(ax_main.get_xlim())
    step = 1 if n_real <= 15 else (2 if n_real <= 30 else 10)
    real_tick_padded = [padding_zeros + i for i in range(0, n_real, step)]
    real_tick_labels = [str(d_first + i) for i in range(0, n_real, step)]
    ax2.set_xticks(real_tick_padded)
    ax2.set_xticklabels(real_tick_labels, fontsize=8)
    ax2.set_xlabel(f"Real domain (D{d_first} south → D{d_last} north)",
                   fontsize=9, labelpad=4)

    for spine in ["top", "right"]:
        ax_main.spines[spine].set_visible(False)

    # ---- Left buffer zoom --------------------------------------------- #
    n_zoom_real = min(n_slope + 3, n_real)
    zoom_real_x = real_x[:n_zoom_real]
    zoom_real_v = rv[:n_zoom_real]

    ax_left.plot(zoom_real_x,   zoom_real_v,   "o-",  color=REAL_COLOR,   ms=5, lw=2,   label="Real")
    ax_left.plot(left_slope_x,  left_slope_v,  "o-",  color=ZOOM_COLOR,   ms=6, lw=2,   label=f"Slope (±{n_slope})")
    ax_left.plot(left_bridge_x, left_bridge_v, "s--", color=BRIDGE_COLOR, ms=5, lw=1.5, label="Bridge")
    ax_left.axvline(padding_zeros - 0.5, color="#444444", lw=1.0, ls="--")
    ax_left.axvline(padding_zeros - n_slope - 0.5, color=BRIDGE_COLOR, lw=0.8, ls=":", alpha=0.7)
    ax_left.axvspan(left_buf_x[0] - 0.5, left_buf_x[-1] + 0.5, color=BUF_COLOR, alpha=0.12)
    ax_left.set_title(f"Left buffer (S of D{d_first})", fontsize=9, fontweight="bold")
    ax_left.set_xlabel("Padded index", fontsize=8)
    ax_left.set_ylabel("Offset (m)", fontsize=8)
    ax_left.tick_params(labelsize=8)
    ax_left.legend(fontsize=7)
    for spine in ["top", "right"]:
        ax_left.spines[spine].set_visible(False)

    # ---- Right buffer zoom -------------------------------------------- #
    zoom_real_x2 = real_x[-n_zoom_real:]
    zoom_real_v2 = rv[-n_zoom_real:]

    ax_right.plot(zoom_real_x2,   zoom_real_v2,   "o-",  color=REAL_COLOR,   ms=5, lw=2,   label="Real")
    ax_right.plot(right_slope_x,  right_slope_v,  "o-",  color=ZOOM_COLOR,   ms=6, lw=2,   label=f"Slope (±{n_slope})")
    ax_right.plot(right_bridge_x, right_bridge_v, "s--", color=BRIDGE_COLOR, ms=5, lw=1.5, label="Bridge")
    ax_right.axvline(padding_zeros + n_real - 0.5, color="#444444", lw=1.0, ls="--")
    ax_right.axvline(padding_zeros + n_real + n_slope - 0.5, color=BRIDGE_COLOR, lw=0.8, ls=":", alpha=0.7)
    ax_right.axvspan(right_buf_x[0] - 0.5, right_buf_x[-1] + 0.5, color=BUF_COLOR, alpha=0.12)
    ax_right.set_title(f"Right buffer (N of D{d_last})", fontsize=9, fontweight="bold")
    ax_right.set_xlabel("Padded index", fontsize=8)
    ax_right.set_ylabel("Offset (m)", fontsize=8)
    ax_right.tick_params(labelsize=8)
    ax_right.legend(fontsize=7)
    for spine in ["top", "right"]:
        ax_right.spines[spine].set_visible(False)

    # ---- Save ---------------------------------------------------------- #
    fig_path = os.path.join(output_dir, f"{output_basename}_buffer_diagnostic.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\nDiagnostic figure saved to:\n  {fig_path}")
    return fig_path


# =============================================================================
# 3. MAIN
# =============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # --- 3.1. Compute relative offsets ---
    result = calculate_relative_offset(
        file_path=RAW_FILE,
        year=YEAR,
        col_map=COL_MAP,
        grids=B3D_GRIDS,
    )

    if result is None:
        print("Offset calculation failed. Exiting.")
        return

    # Hard stop if any requested domain is missing — a silently short array
    # would misalign every downstream CASCADE domain.
    if len(result) != N_REAL:
        missing = sorted(set(B3D_GRIDS) - set(result["Domain_ID"]))
        print(f"\nERROR: expected {N_REAL} domains (D{START_DOMAIN}–D{END_DOMAIN}), "
              f"got {len(result)}. Missing: {missing}")
        print("Exiting rather than writing a misaligned array.")
        return

    # Save unpadded file with Domain_ID
    unpadded_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_CASCADE_Input_unpadded.csv")
    result.to_csv(unpadded_path, index=False)
    print(f"\nUnpadded offsets saved to:\n  {unpadded_path}")

    # CASCADE-format: raw_offset column only
    cascade_df = result[[str(YEAR)]].copy()
    cascade_unpadded_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_CASCADE_Input.csv")
    cascade_df.to_csv(cascade_unpadded_path, index=False)
    print(f"Unpadded CASCADE-format file saved to:\n  {cascade_unpadded_path}")

    # --- 3.2. Pad with hybrid slope + linear bridge buffers ---
    padded_df, diag = pad_for_cascade(
        df=cascade_df,
        padding_zeros=PADDING_ZEROS,
        target_length=TARGET_LENGTH,
        extrap_fit_domains=EXTRAP_FIT_DOMAINS,
        slope_buffer_domains=SLOPE_BUFFER_DOMAINS,
        clip_at_zero=CLIP_BUFFERS_AT_ZERO,
        re_reference=RE_REFERENCE_PADDED,
        start_domain=START_DOMAIN,
    )

    padded_path = os.path.join(OUTPUT_DIR, f"{OUTPUT_BASENAME}_PADDED_{TARGET_LENGTH}.csv")
    padded_df.to_csv(padded_path, index=False)
    print(f"\nSUCCESS: Padded CASCADE input saved to:\n  {padded_path}")

    # --- 3.3. Diagnostic figure ---
    plot_buffer_diagnostic(
        diag=diag,
        year=YEAR,
        padding_zeros=PADDING_ZEROS,
        community_zones=COMMUNITY_ZONES,
        output_dir=OUTPUT_DIR,
        output_basename=OUTPUT_BASENAME,
    )


if __name__ == "__main__":
    main()
