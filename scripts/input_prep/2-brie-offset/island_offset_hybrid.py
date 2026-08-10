"""
Hatteras CASCADE Dune Offset Pipeline
====================================

This script:
1. Reads a single raw dune–baseline intersection CSV (1984).
2. Calculates the relative dune raw_offset per domain (meters, baseline = minimum).
3. Pads the result for CASCADE using outward linear extrapolation from each
   real edge, so that buffer domains closest to the real coast follow the
   local shoreline trend rather than bridging to the opposite end.
4. Saves a diagnostic figure showing the full 120-domain raw_offset profile.

Author: Hannah A. Henry (extrapolation buffer version)
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# 1. USER CONFIGURATION
# =============================================================================

YEAR = 2004
RAW_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\2004\2004_duneline_offset_raw.csv"

OUTPUT_DIR    = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\hindcast_2004_v3"
OUTPUT_BASENAME = "Island_Dune_Offsets_2004"

START_DOMAIN = 1
END_DOMAIN   = 90
B3D_GRIDS    = list(range(START_DOMAIN, END_DOMAIN + 1))

PADDING_ZEROS = 15
TARGET_LENGTH = (END_DOMAIN - START_DOMAIN + 1) + 2 * PADDING_ZEROS  # 120

# Number of real domains from each edge used to fit the local extrapolation trend.
EXTRAP_FIT_DOMAINS = 10

# Number of buffer domains on each side that follow the local coastline slope.
# The remaining (PADDING_ZEROS - SLOPE_BUFFER_DOMAINS) buffer domains on each
# side are filled by a linear bridge connecting the two slope tails.
SLOPE_BUFFER_DOMAINS = 10

COL_MAP = {
    "Domain_ID": "domain_id",
    "Distance":  "ORIG_LEN",
    "Transect":  "LineID",
}

# Community zone annotations for the diagnostic figure
# (real domain numbers, 1-indexed)
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
    """Compute mean relative dune raw_offset per domain from raw CSV."""
    print(f"\n--- Processing {year} ---")
    print(f"Input file: {file_path}")

    try:
        raw_df = pd.read_csv(file_path)
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

    for grid_id in grids:
        subset = data_df[data_df["B3D_Grid"] == int(grid_id)]
        if subset.empty:
            print(f"  Warning: No data for domain {grid_id}.")
            continue

        try:
            min_t = int(subset["Transect"].min())
            max_t = int(subset["Transect"].max())
        except ValueError:
            print(f"  Warning: Invalid transect data for domain {grid_id}. Skipping.")
            continue

        distances = []
        for t_id in range(min_t, max_t + 1):
            t_vals = subset[subset["Transect"] == t_id]
            if not t_vals.empty:
                distances.append(t_vals.iloc[0]["Distance"])

        if not distances:
            continue

        mean_distances.append(float(np.mean(distances)))
        seen_domains.append(grid_id)

    if not mean_distances:
        print(f"FATAL: No valid domains processed for year {year}.")
        return None

    baseline         = min(mean_distances)
    relative_offsets = np.subtract(mean_distances, baseline)

    print(f"  {len(mean_distances)} domains processed.")
    print(f"  Baseline distance = {baseline:.3f} m (min mean).")

    return pd.DataFrame({"Domain_ID": seen_domains, str(year): relative_offsets})


def pad_for_cascade(df, padding_zeros, target_length,
                    extrap_fit_domains=5, slope_buffer_domains=5):
    """
    Pad raw_offset array for CASCADE using a hybrid buffer strategy:

      Left buffer  (15 domains, outermost → D1):
        - Innermost `slope_buffer_domains` (closest to D1): follow local
          coastline slope, anchored exactly at D1 so no gap at boundary.
        - Remaining outer domains: linear bridge from the slope tail to the
          right buffer's outer end, keeping the full array connected.

      Right buffer (15 domains, D90 → outermost):
        - Innermost `slope_buffer_domains` (closest to D90): follow local
          coastline slope, anchored exactly at D90.
        - Remaining outer domains: linear bridge (same bridge as above,
          approached from the other end).

    This guarantees:
      - No discontinuity at either real-domain boundary.
      - The buffer interior is connected (no jumps anywhere).
      - Behaviour in the outer buffer is irrelevant to the real simulation.

    Parameters
    ----------
    df : pd.DataFrame
        Single raw_offset column, real domains D1–D90, south → north order.
    padding_zeros : int
        Buffer domains on each side (15).
    target_length : int
        Total padded length (120).
    extrap_fit_domains : int
        Real domains used to estimate local edge slope (default 5).
    slope_buffer_domains : int
        Buffer domains on each side that follow the local slope (default 5).
        Must be < padding_zeros.

    Returns
    -------
    padded : pd.DataFrame
    diag   : dict  (arrays for diagnostic plot)
    """
    assert slope_buffer_domains < padding_zeros, \
        "slope_buffer_domains must be less than padding_zeros"

    data_columns = list(df.columns)
    col          = data_columns[0]
    values       = df[col].to_numpy(dtype=float)   # shape (90,)

    # ------------------------------------------------------------------ #
    # 1. Estimate local slopes at each edge                               #
    # ------------------------------------------------------------------ #
    x_fit   = np.arange(extrap_fit_domains, dtype=float)

    # Left edge slope (fit D1..D5, anchor at D1 = values[0])
    m_left, _ = np.polyfit(x_fit, values[:extrap_fit_domains], 1)

    # Right edge slope (fit D86..D90, anchor at D90 = values[-1])
    m_right, _ = np.polyfit(x_fit, values[-extrap_fit_domains:], 1)

    # ------------------------------------------------------------------ #
    # 2. Left slope segment (innermost, closest to D1)                    #
    #    Steps: -1 (adjacent to D1) … -slope_buffer_domains (outermost   #
    #    of slope segment).  Anchor: values[0].                          #
    # ------------------------------------------------------------------ #
    left_slope_steps  = np.arange(-1, -(slope_buffer_domains + 1), -1, dtype=float)
    left_slope_values = values[0] + m_left * left_slope_steps   # anchored at D1
    left_slope_values = np.clip(left_slope_values, 0.0, None)
    # left_slope_values[0]  → adjacent to D1  (step = -1)
    # left_slope_values[-1] → outermost of slope segment

    # ------------------------------------------------------------------ #
    # 3. Right slope segment (innermost, closest to D90)                  #
    # ------------------------------------------------------------------ #
    right_slope_steps  = np.arange(1, slope_buffer_domains + 1, dtype=float)
    right_slope_values = values[-1] + m_right * right_slope_steps
    right_slope_values = np.clip(right_slope_values, 0.0, None)
    # right_slope_values[0]  → adjacent to D90
    # right_slope_values[-1] → outermost of slope segment

    # ------------------------------------------------------------------ #
    # 4. Linear bridge + final assembly                                   #
    #                                                                     #
    #  Think of the full padded array as one continuous sequence:         #
    #    idx 0            = outermost left buffer                         #
    #    idx 14           = adjacent to D1                                #
    #    idx 15–104       = real domains D1–D90                           #
    #    idx 105          = adjacent to D90                               #
    #    idx 119          = outermost right buffer                        #
    #                                                                     #
    #  Slope segments (length = slope_buffer_domains = 5):               #
    #    Left slope  : idx 10–14  (inner left buffer, touching D1)       #
    #    Right slope : idx 105–109 (inner right buffer, touching D90)    #
    #                                                                     #
    #  Bridge (length = n_bridge_each per side = 10):                    #
    #    Left bridge : idx 0–9   (outer left buffer)                     #
    #    Right bridge: idx 110–119 (outer right buffer)                  #
    #                                                                     #
    #  The bridge is ONE linspace from left slope TAIL (idx 10 value)    #
    #  to right slope TAIL (idx 109 value), with 20 interior points.     #
    #  Left bridge  = first 10 points (running left→right, ascending)    #
    #  Right bridge = last  10 points (running left→right, continuing)   #
    #                                                                     #
    #  Key: left_slope_values is ordered adjacent-to-D1 first,           #
    #  so left_slope_values[-1] is the OUTERMOST left slope point        #
    #  (i.e. the value at padded idx 10, the bridge's left anchor).      #
    # ------------------------------------------------------------------ #
    n_bridge_each = padding_zeros - slope_buffer_domains   # 10

    # ------------------------------------------------------------------ #
    # Build bridge as 21 points INCLUDING both slope tails as endpoints.  #
    # full_bridge[0]               = L = left  slope tail                 #
    # full_bridge[n_bridge_each]   = M = midpoint                         #
    # full_bridge[2*n_bridge_each] = R = right slope tail                 #
    #                                                                      #
    # Left bridge stored outermost→innermost in padded space (idx 0→9):   #
    #   idx 9 (innermost) = full_bridge[0] = L  ← zero-gap connect to    #
    #   idx 10 = left_slope_values[-1] = L                                #
    #   idx 0 (outermost) = full_bridge[9] (near midpoint)                #
    #                                                                      #
    # Right bridge stored innermost→outermost in padded space (idx 110→119):
    #   idx 110 (innermost) = full_bridge[20] = R ← zero-gap connect to  #
    #   idx 109 = right_slope_values[-1] = R                              #
    #   idx 119 (outermost) = full_bridge[11] (near midpoint)             #
    # ------------------------------------------------------------------ #
    full_bridge = np.linspace(
        left_slope_values[-1],      # L: left slope tail
        right_slope_values[-1],     # R: right slope tail
        n_bridge_each * 2 + 1,      # 21 points including both endpoints
    )

    # Left bridge:  [fb[9], fb[8], ..., fb[0]]  — outermost to innermost
    left_bridge_values  = full_bridge[n_bridge_each - 1::-1]

    # Right bridge: [fb[20], fb[19], ..., fb[11]] — innermost to outermost
    right_bridge_values = full_bridge[2 * n_bridge_each:n_bridge_each:-1]

    # ------------------------------------------------------------------ #
    # 5. Assemble full left and right buffer arrays                       #
    #  left_full  (idx 0=outermost, idx 14=adjacent to D1):              #
    #    [left_bridge] + [left_slope reversed]                            #
    #    left_full[9]  = fb[0] = L                                       #
    #    left_full[10] = left_slope_values[-1] = L  ← ZERO GAP           #
    #  right_full (idx 0=adjacent to D90, idx 14=outermost):             #
    #    [right_slope] + [right_bridge]                                   #
    #    right_full[4] = right_slope_values[-1] = R                      #
    #    right_full[5] = fb[20] = R              ← ZERO GAP              #
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

    print(f"\nPadding summary (hybrid slope + bridge):")
    print(f"  Slope domains per side : {slope_buffer_domains}")
    print(f"  Bridge domains per side: {padding_zeros - slope_buffer_domains}")
    print(f"  D1  boundary check     : left_buf[-1]={left_full[-1]:.2f}  D1={values[0]:.2f}  "
          f"diff={abs(left_full[-1]-values[0]):.4f} m")
    print(f"  D90 boundary check     : right_buf[0]={right_full[0]:.2f}  D90={values[-1]:.2f}  "
          f"diff={abs(right_full[0]-values[-1]):.4f} m")
    print(f"  Left  buffer range : {left_full.min():.2f} – {left_full.max():.2f} m")
    print(f"  Right buffer range : {right_full.min():.2f} – {right_full.max():.2f} m")
    print(f"  Padded length      : {len(padded)} (target: {target_length})")

    if len(padded) != target_length:
        print(f"  ERROR: length mismatch!")

    diag = {
        "left_values":         left_full,           # outermost→D1, length=15
        "real_values":         values,               # D1→D90, length=90
        "right_values":        right_full,           # D90→outermost, length=15
        "slope_buffer_domains": slope_buffer_domains,
        "extrap_fit_domains":  extrap_fit_domains,
        "m_left":              m_left,
        "m_right":             m_right,
    }

    return padded, diag


def plot_buffer_diagnostic(diag, year, padding_zeros, community_zones,
                            output_dir, output_basename):
    """
    Save a diagnostic figure of the full padded raw_offset profile (120 domains).

    Layout:
      - Full profile across all 120 padded indices
      - Buffer zones shaded in salmon
      - Real-domain zone annotations along the top
      - Inset zoom panels for left and right buffer transitions
    """
    lv = diag["left_values"]    # shape (15,) outermost first
    rv = diag["real_values"]    # shape (90,)
    rb = diag["right_values"]   # shape (15,) D90-adjacent first

    n_slope  = diag["slope_buffer_domains"]
    n_bridge = padding_zeros - n_slope

    # Full padded array in plot order (index 0 = leftmost buffer)
    full = np.concatenate([lv, rv, rb])   # 120 values

    # Index ranges (padded space)
    left_buf_x   = np.arange(padding_zeros)
    real_x       = np.arange(padding_zeros, padding_zeros + len(rv))
    right_buf_x  = np.arange(padding_zeros + len(rv), len(full))

    # Index directly into `full` — no re-slicing of diag sub-arrays.
    # Left buffer  [0 .. padding_zeros-1]:
    #   [0 .. n_bridge-1]      = bridge (outermost, far from D1)
    #   [n_bridge .. 14]       = slope  (adjacent to D1)
    # Right buffer [padding_zeros+90 .. 119]:
    #   [+0 .. +n_slope-1]     = slope  (adjacent to D90)
    #   [+n_slope .. +14]      = bridge (outermost, far from D90)
    r_start = padding_zeros + len(rv)

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
    fig.suptitle(
        f"Buffer Diagnostic — Dune Offset Profile  |  {year}  |  "
        f"Slope: ±{n_slope} domains  |  Bridge: {n_bridge} domains each side",
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

    for xv in [padding_zeros - 0.5, padding_zeros + len(rv) - 0.5]:
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
                 ms=3, lw=1.8, zorder=5, label="Real domains (D1–D90)")

    # Community zone annotations along the top
    ax_main.set_ylim(bottom=-0.5)
    _, ymax = ax_main.get_ylim()
    zone_colors = ["#e8e8f0", "#d8d8ec"] * 10

    for zi, (d_start, d_end, label) in enumerate(community_zones):
        px_start = (d_start - 1) + padding_zeros
        px_end   = (d_end   - 1) + padding_zeros
        ax_main.axvspan(px_start - 0.5, px_end + 0.5,
                        color=zone_colors[zi], alpha=0.35, zorder=1)
        ax_main.text((px_start + px_end) / 2, ymax * 0.97, label,
                     ha="center", va="top", fontsize=7.5,
                     color="#333355", rotation=90 if len(label) > 8 else 0,
                     clip_on=True)

    ax_main.set_xlabel("Padded domain index  (0–119;  real = 15–104)", fontsize=11)
    ax_main.set_ylabel("Dune raw_offset (m)", fontsize=11)
    ax_main.set_xlim(-0.5, len(full) - 0.5)
    ax_main.legend(fontsize=9, loc="upper center", framealpha=0.9)

    ax2 = ax_main.twiny()
    ax2.set_xlim(ax_main.get_xlim())
    real_tick_padded = [padding_zeros + i for i in range(0, 90, 10)]
    real_tick_labels = [str(i + 1) for i in range(0, 90, 10)]
    ax2.set_xticks(real_tick_padded)
    ax2.set_xticklabels(real_tick_labels, fontsize=8)
    ax2.set_xlabel("Real domain (D1 = Cape Point → D90 = Rodanthe)", fontsize=9, labelpad=4)

    for spine in ["top", "right"]:
        ax_main.spines[spine].set_visible(False)

    # ---- Left buffer zoom (show bridge + slope + first few real domains) #
    n_show = n_slope + n_bridge + 3
    zoom_real_x  = real_x[:n_slope + 3]
    zoom_real_v  = rv[:n_slope + 3]

    ax_left.plot(zoom_real_x,  zoom_real_v,   "o-", color=REAL_COLOR,   ms=5, lw=2, label="Real")
    ax_left.plot(left_slope_x, left_slope_v,  "o-", color=ZOOM_COLOR,   ms=6, lw=2, label=f"Slope (±{n_slope})")
    ax_left.plot(left_bridge_x,left_bridge_v, "s--",color=BRIDGE_COLOR, ms=5, lw=1.5,label="Bridge")
    ax_left.axvline(padding_zeros - 0.5, color="#444444", lw=1.0, ls="--")
    ax_left.axvline(padding_zeros - n_slope - 0.5, color=BRIDGE_COLOR, lw=0.8, ls=":", alpha=0.7)
    ax_left.axvspan(left_buf_x[0] - 0.5, left_buf_x[-1] + 0.5, color=BUF_COLOR, alpha=0.12)
    ax_left.set_title("Left buffer (S cape)", fontsize=9, fontweight="bold")
    ax_left.set_xlabel("Padded index", fontsize=8)
    ax_left.set_ylabel("Offset (m)", fontsize=8)
    ax_left.tick_params(labelsize=8)
    ax_left.legend(fontsize=7)
    for spine in ["top", "right"]:
        ax_left.spines[spine].set_visible(False)

    # ---- Right buffer zoom -------------------------------------------- #
    zoom_real_x2 = real_x[-(n_slope + 3):]
    zoom_real_v2 = rv[-(n_slope + 3):]

    ax_right.plot(zoom_real_x2,  zoom_real_v2,   "o-", color=REAL_COLOR,   ms=5, lw=2, label="Real")
    ax_right.plot(right_slope_x, right_slope_v,  "o-", color=ZOOM_COLOR,   ms=6, lw=2, label=f"Slope (±{n_slope})")
    ax_right.plot(right_bridge_x,right_bridge_v, "s--",color=BRIDGE_COLOR, ms=5, lw=1.5,label="Bridge")
    ax_right.axvline(padding_zeros + len(rv) - 0.5, color="#444444", lw=1.0, ls="--")
    ax_right.axvline(padding_zeros + len(rv) + n_slope - 0.5, color=BRIDGE_COLOR, lw=0.8, ls=":", alpha=0.7)
    ax_right.axvspan(right_buf_x[0] - 0.5, right_buf_x[-1] + 0.5, color=BUF_COLOR, alpha=0.12)
    ax_right.set_title("Right buffer (N Rodanthe)", fontsize=9, fontweight="bold")
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
    print(f"\n✓ Diagnostic figure saved to:\n  {fig_path}")
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
    )

    if padded_df is None:
        print("Padded comparison not created due to errors.")
        return

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
