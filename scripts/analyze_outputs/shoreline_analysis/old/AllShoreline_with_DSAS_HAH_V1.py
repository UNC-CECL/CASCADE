import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# --- SETTINGS: Update these values only ---
# ============================================================
NPZ_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_Natural_State_ModStorms.csv\HAT_1978_1997_Natural_State_ModStorms.csv.npz"

# Root directory where all shoreline snapshot outputs will be stored.
OUTPUT_ROOT_DIR = "../HAT_78_97_Mod_Storm"

START_YEAR = 1978  # Label for first time step
END_YEAR = 1997  # Label for last time step
RELATIVE_CHANGE = False  # True = plot shoreline change relative to t=0, False = absolute position
TO_METERS = True  # True = convert dam → m by multiplying by 10

# Domain indices for buffers
LEFT_BUFFER_END = 14  # last left-buffer CASCADE domain index (0–14 = 15 left buffers)
RIGHT_BUFFER_START = 105  # first right-buffer CASCADE domain index

# Optional: cleaned DSAS per-domain CSV to overlay (set to None to skip)
DSAS_DOMAIN_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_CLEAN.csv"
DSAS_DOMAIN_COL = "domain_id"  # GIS domain ID field in CLEAN file (1–90)
DSAS_CHANGE_COL = "obs_total_change_m"  # total shoreline change column in CLEAN file (m)

# NEW: Diagnostic settings
FLIP_CASCADE_ALONGSHORE = False  # Set to True if CASCADE and GIS are indexed in opposite directions


# ============================================================


def extract_run_name(npz_path):
    """Return filename without extension."""
    base = os.path.basename(npz_path)
    return os.path.splitext(base)[0]


def load_cascade(npz_path):
    """Load the CASCADE object from a .npz file."""
    data = np.load(npz_path, allow_pickle=True)
    return data["cascade"][0]


def get_x_s_TS(b3d):
    """Get shoreline time series from a Barrier3D object."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True, flip_alongshore=False):
    """
    Build absolute shoreline[t, domain] matrix.

    Parameters:
    -----------
    flip_alongshore : bool
        If True, flip the alongshore direction
    """
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)
    nt = len(get_x_s_TS(b3d_list[0]))

    shoreline = np.zeros((nt, ndom))

    for j in range(ndom):
        xs = get_x_s_TS(b3d_list[j])
        shoreline[:, j] = xs

    if to_meters:
        shoreline = shoreline * 10.0  # dam → m

    if flip_alongshore:
        shoreline = np.flip(shoreline, axis=1)
        print("*** CASCADE alongshore direction flipped ***")

    return shoreline


def verify_data_alignment(shoreline_abs, obs_domains, obs_change,
                          left_buffer_end, right_buffer_start, start_year, end_year):
    """Print comprehensive diagnostic info to verify units, alignment, and orientation."""
    print("\n" + "=" * 70)
    print("DATA VERIFICATION DIAGNOSTICS")
    print("=" * 70)

    nt, ndom = shoreline_abs.shape

    # Calculate time span
    time_span_years = end_year - start_year
    if time_span_years == 0:
        time_span_years = nt - 1

    # Basic info
    print(f"\n--- BASIC INFO ---")
    print(f"Time steps: {nt} ({start_year} to {end_year})")
    print(f"Time span: {time_span_years} years")
    print(f"Total CASCADE domains: {ndom}")
    print(f"Real island domains: {left_buffer_end + 1} to {right_buffer_start - 1} "
          f"({right_buffer_start - left_buffer_end - 1} domains)")

    # Check CASCADE data range (absolute positions)
    print(f"\n--- CASCADE SHORELINE POSITIONS (ABSOLUTE) ---")
    print(f"Initial (year {start_year}) range: {shoreline_abs[0].min():.1f} to {shoreline_abs[0].max():.1f} m")
    print(f"Final (year {end_year}) range: {shoreline_abs[-1].min():.1f} to {shoreline_abs[-1].max():.1f} m")

    # Calculate total change and rate
    cascade_total_change = shoreline_abs[-1] - shoreline_abs[0]
    cascade_rate = cascade_total_change / time_span_years

    print(f"\n--- CASCADE TOTAL CHANGE (FINAL - INITIAL) ---")
    print(f"Full domain change range: {cascade_total_change.min():.1f} to {cascade_total_change.max():.1f} m")
    print(f"Full domain mean change: {cascade_total_change.mean():.2f} m")

    print(f"\n--- CASCADE CHANGE RATE ---")
    print(f"Full domain rate range: {cascade_rate.min():.2f} to {cascade_rate.max():.2f} m/yr")
    print(f"Full domain mean rate: {cascade_rate.mean():.2f} m/yr")

    # Real island portion only (excluding buffers)
    real_island_change = cascade_total_change[left_buffer_end + 1:right_buffer_start]
    real_island_rate = cascade_rate[left_buffer_end + 1:right_buffer_start]

    print(f"\nReal island change range: {real_island_change.min():.1f} to {real_island_change.max():.1f} m")
    print(f"Real island mean change: {real_island_change.mean():.2f} m")
    print(f"Real island std dev: {real_island_change.std():.2f} m")

    print(f"\nReal island rate range: {real_island_rate.min():.2f} to {real_island_rate.max():.2f} m/yr")
    print(f"Real island mean rate: {real_island_rate.mean():.2f} m/yr")
    print(f"Real island rate std dev: {real_island_rate.std():.2f} m/yr")

    # Buffer regions
    left_buffer_change = cascade_total_change[0:left_buffer_end + 1]
    right_buffer_change = cascade_total_change[right_buffer_start:]
    print(f"\nLeft buffer mean change: {left_buffer_change.mean():.2f} m")
    print(f"Right buffer mean change: {right_buffer_change.mean():.2f} m")

    # Check for suspicious patterns
    print(f"\n--- PATTERN CHECKS ---")
    if abs(real_island_change.mean()) > 100:
        print("⚠️  WARNING: Mean change > 100m seems large. Check units conversion.")
    if abs(real_island_change.mean()) < 1:
        print("⚠️  WARNING: Mean change < 1m seems small. Check units conversion.")

    # DSAS comparison
    if obs_change is not None and obs_domains is not None:
        # Convert DSAS total change to rate
        obs_rate = obs_change / time_span_years

        print(f"\n--- OBSERVED DSAS DATA ---")
        print(f"Number of DSAS observations: {len(obs_change)}")
        print(f"DSAS change range: {obs_change.min():.1f} to {obs_change.max():.1f} m")
        print(f"DSAS mean change: {obs_change.mean():.2f} m")
        print(f"DSAS std dev: {obs_change.std():.2f} m")
        print(f"\nDSAS rate range: {obs_rate.min():.2f} to {obs_rate.max():.2f} m/yr")
        print(f"DSAS mean rate: {obs_rate.mean():.2f} m/yr")
        print(f"DSAS rate std dev: {obs_rate.std():.2f} m/yr")
        print(f"DSAS CASCADE domain range: {obs_domains.min()} to {obs_domains.max()}")

        # Get corresponding CASCADE values at observed locations
        cascade_at_obs = cascade_total_change[obs_domains.astype(int)]
        cascade_rate_at_obs = cascade_rate[obs_domains.astype(int)]

        print(f"\n--- MODEL vs OBSERVED COMPARISON ---")
        print(f"CASCADE at observed locations - mean change: {cascade_at_obs.mean():.2f} m")
        print(f"DSAS observed - mean change: {obs_change.mean():.2f} m")
        print(f"CASCADE at observed locations - mean rate: {cascade_rate_at_obs.mean():.2f} m/yr")
        print(f"DSAS observed - mean rate: {obs_rate.mean():.2f} m/yr")

        ratio = cascade_at_obs.mean() / obs_change.mean() if obs_change.mean() != 0 else np.inf
        print(f"\nRatio (CASCADE/DSAS total change): {ratio:.3f}")

        ratio_rate = cascade_rate_at_obs.mean() / obs_rate.mean() if obs_rate.mean() != 0 else np.inf
        print(f"Ratio (CASCADE/DSAS rate): {ratio_rate:.3f}")

        if 0.9 <= ratio <= 1.1:
            print("✓ GOOD: Ratio near 1.0 - units likely correct")
        elif 9 <= ratio <= 11:
            print("⚠️  WARNING: Ratio ~10 - CASCADE might already be in meters (set TO_METERS=False)")
        elif 0.09 <= ratio <= 0.11:
            print("⚠️  WARNING: Ratio ~0.1 - CASCADE conversion might be wrong (check if dam→m is correct)")
        else:
            print("⚠️  WARNING: Unexpected ratio - check units and alignment")

        # Correlation check
        correlation = np.corrcoef(cascade_at_obs, obs_change)[0, 1]
        print(f"\nCorrelation (CASCADE vs DSAS): {correlation:.3f}")

        if correlation > 0.7:
            print("✓ GOOD: Strong positive correlation - spatial patterns align")
        elif correlation < -0.7:
            print("⚠️  WARNING: Strong negative correlation - might need FLIP_CASCADE_ALONGSHORE=True")
        else:
            print("⚠️  WARNING: Weak correlation - check alignment or data quality")

        # Sample comparison at a few points
        print(f"\n--- SAMPLE POINT COMPARISON ---")
        print(f"{'GIS Domain':<12} {'CASCADE Idx':<14} {'CASCADE Δ (m)':<16} {'DSAS Δ (m)':<14} {'Difference':<12}")
        print("-" * 70)
        sample_indices = [0, len(obs_domains) // 4, len(obs_domains) // 2, 3 * len(obs_domains) // 4, -1]
        for idx in sample_indices:
            gis_dom = obs_domains[idx] - left_buffer_end
            print(f"{gis_dom:<12} {obs_domains[idx]:<14} {cascade_at_obs[idx]:<16.2f} "
                  f"{obs_change[idx]:<14.2f} {cascade_at_obs[idx] - obs_change[idx]:<12.2f}")

        print("\n" + "=" * 70)
        print("RECOMMENDATIONS:")
        print("-" * 70)

        if ratio > 5:
            print("• Consider setting TO_METERS = False")
        elif ratio < 0.2:
            print("• Check if CASCADE units are actually in meters, not decameters")

        if correlation < -0.5:
            print("• Consider setting FLIP_CASCADE_ALONGSHORE = True")
        elif 0.5 < correlation < 0.7:
            print("• Moderate correlation - verify domain index mapping")

    print("• Visually inspect the 'shoreline_change_rate_over_run.png' plot")
    print("• Check if spatial patterns match known Hatteras geography")
    print("=" * 70 + "\n")


def plot_shoreline_snapshots_by_domain(
        shoreline_abs,
        run_name,
        outdir,
        start_year=0,
        relative=False,
        units="m",
):
    """For each time step (year), plot shoreline vs domain number."""
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline_abs.shape
    domain_numbers = np.arange(ndom)

    # For plotting, optionally convert to relative change
    if relative:
        shoreline = shoreline_abs - shoreline_abs[0, :]
    else:
        shoreline = shoreline_abs

    for t in range(nt):
        year = start_year + t

        plt.figure(figsize=(12, 4))
        plt.plot(domain_numbers, shoreline[t, :], marker="o", linewidth=1.2)

        plt.xticks(domain_numbers[::10])
        plt.xlabel("Domain number (alongshore index)")

        if relative:
            plt.ylabel(f"Shoreline change ({units})")
            plt.title(f"{run_name} – Shoreline Change Across Domains – Year {year}")
        else:
            plt.ylabel(f"Shoreline position ({units})")
            plt.title(f"{run_name} – Shoreline Position Across Domains – Year {year}")

        plt.grid(alpha=0.3)
        plt.tight_layout()

        fname = os.path.join(outdir, f"snapshot_year_{year}.png")
        plt.savefig(fname, dpi=200)
        plt.close()

        print(f"Saved: {fname}")


def plot_mean_shoreline_timeseries(
        shoreline_abs,
        run_name,
        outdir,
        start_year=0,
        units="m",
):
    """Plot the mean absolute shoreline (across all domains) as a time series."""
    os.makedirs(outdir, exist_ok=True)

    nt, _ = shoreline_abs.shape
    years = np.arange(nt) + start_year

    mean_shoreline = shoreline_abs.mean(axis=1)

    plt.figure(figsize=(8, 4))
    plt.plot(years, mean_shoreline, linewidth=1.8)

    plt.xlabel("Year")
    plt.ylabel(f"Mean shoreline position ({units})")
    plt.title(f"{run_name} – Mean Shoreline Position Over Time")

    plt.grid(alpha=0.3)
    plt.tight_layout()

    fname = os.path.join(outdir, "shoreline_mean_timeseries.png")
    plt.savefig(fname, dpi=200)
    plt.close()

    print(f"Saved mean shoreline time series: {fname}")


def plot_total_shoreline_change_over_run(
        shoreline_abs,
        run_name,
        outdir,
        units="m",
        start_year=1978,
        end_year=1997,
        left_buffer_end=14,
        right_buffer_start=105,
        obs_domains=None,
        obs_change=None,
):
    """
    Plot shoreline change RATE over the whole run for each domain.

    Bottom x-axis: CASCADE domain index.
    Top x-axis:    GIS domain number (1–90) over the real-island portion.
    """
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline_abs.shape
    domain_numbers = np.arange(ndom)

    # Calculate time span
    time_span_years = end_year - start_year
    if time_span_years == 0:
        time_span_years = nt - 1

    # Total change = final - initial
    total_change = shoreline_abs[-1, :] - shoreline_abs[0, :]

    # Convert to rate (m/yr)
    change_rate = total_change / time_span_years

    # Save both total change and rate to CSV
    model_csv = os.path.join(outdir, "total_shoreline_change_model.csv")
    pd.DataFrame({
        "domain": domain_numbers,
        "model_total_change_m": total_change,
        "model_rate_m_per_yr": change_rate
    }).to_csv(model_csv, index=False)
    print(f"Saved model total shoreline change CSV: {model_csv}")

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(domain_numbers, change_rate, marker="o", linewidth=1.5, label="Model", markersize=4)

    # Domain ticks every 10
    ax.set_xticks(domain_numbers[::10])
    ax.set_xlabel("Domain number (alongshore index)")

    # Zero-change reference line
    ax.axhline(0, linestyle="--", linewidth=1, alpha=0.7, color='gray')

    # Shade buffer regions
    if left_buffer_end is not None and left_buffer_end >= 0:
        ax.axvspan(domain_numbers[0] - 0.5,
                   left_buffer_end + 0.5,
                   alpha=0.15, color='red', label='Buffers')
    if right_buffer_start is not None and right_buffer_start < ndom:
        ax.axvspan(right_buffer_start - 0.5,
                   domain_numbers[-1] + 0.5,
                   alpha=0.15, color='red')

    # Overlay observed DSAS change rate if provided
    if obs_domains is not None and obs_change is not None:
        obs_rate = obs_change / time_span_years
        ax.plot(obs_domains, obs_rate, marker="x", linewidth=1.2,
                label="Observed", markersize=6, linestyle='--')

    # Secondary x-axis for GIS domain numbers (top)
    def cas_to_gis(x):
        return x - left_buffer_end

    def gis_to_cas(x):
        return x + left_buffer_end

    top_ax = ax.secondary_xaxis('top', functions=(cas_to_gis, gis_to_cas))
    gis_ticks = np.arange(10, 91, 10)
    top_ax.set_xticks(gis_ticks)
    top_ax.set_xlabel("GIS domain number (1–90)")
    top_ax.tick_params(labelsize=8)

    # Labels and title
    ax.set_ylabel(f"Shoreline change rate (m/yr)")
    ax.set_title(f"{run_name} – Shoreline Change Rate ({start_year}–{end_year})")

    # Optional text labels for regions
    y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
    y_top = ax.get_ylim()[1] - 0.05 * y_range
    if left_buffer_end is not None and right_buffer_start is not None:
        mid_real = (left_buffer_end + right_buffer_start) / 2
        ax.text(mid_real, y_top, "Real island",
                ha="center", va="top", fontsize=9, fontweight='bold')
        ax.text((domain_numbers[0] + left_buffer_end) / 2,
                y_top, "Left buffer",
                ha="center", va="top", fontsize=8, style='italic')
        ax.text((right_buffer_start + domain_numbers[-1]) / 2,
                y_top, "Right buffer",
                ha="center", va="top", fontsize=8, style='italic')

    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    fig.tight_layout()

    fname = os.path.join(outdir, "shoreline_change_rate_over_run.png")
    fig.savefig(fname, dpi=200)
    plt.close(fig)

    print(f"Saved shoreline change rate plot: {fname}")


def main():
    run_name = extract_run_name(NPZ_FILE)
    run_outdir = os.path.join(OUTPUT_ROOT_DIR, run_name)

    print(f"\n{'=' * 70}")
    print(f"CASCADE SHORELINE ANALYSIS: {run_name}")
    print(f"{'=' * 70}\n")

    cascade = load_cascade(NPZ_FILE)

    # Build absolute shoreline matrix
    shoreline_abs = build_shoreline_matrix(
        cascade,
        to_meters=TO_METERS,
        flip_alongshore=FLIP_CASCADE_ALONGSHORE,
    )

    units = "m" if TO_METERS else "dam"

    # Load observed per-domain DSAS change
    obs_domains = None
    obs_change = None
    if DSAS_DOMAIN_CSV is not None:
        print(f"Loading DSAS observations from: {DSAS_DOMAIN_CSV}")
        dsas = pd.read_csv(DSAS_DOMAIN_CSV)

        # Raw GIS domain IDs (1–90)
        gis_domains = dsas[DSAS_DOMAIN_COL].to_numpy()

        # Map GIS domains (1–90) → CASCADE real domains (15–104)
        obs_domains_all = gis_domains + LEFT_BUFFER_END
        obs_change_all = dsas[DSAS_CHANGE_COL].to_numpy()

        # Keep only observed points that fall on real-island CASCADE domains
        mask = (obs_domains_all > LEFT_BUFFER_END) & (obs_domains_all < RIGHT_BUFFER_START)
        obs_domains = obs_domains_all[mask]
        obs_change = obs_change_all[mask]

        print(f"Loaded {len(obs_change)} DSAS observations for real island domains\n")

    # DIAGNOSTIC VERIFICATION
    verify_data_alignment(
        shoreline_abs,
        obs_domains,
        obs_change,
        LEFT_BUFFER_END,
        RIGHT_BUFFER_START,
        START_YEAR,
        END_YEAR
    )

    # GENERATE PLOTS
    print(f"\n{'=' * 70}")
    print("GENERATING PLOTS")
    print(f"{'=' * 70}\n")

    # 1) Snapshots across domains for each year
    print("Creating yearly snapshots...")
    plot_shoreline_snapshots_by_domain(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        relative=RELATIVE_CHANGE,
        units=units,
    )

    # 2) Mean shoreline time series
    print("\nCreating mean shoreline time series...")
    plot_mean_shoreline_timeseries(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        units=units,
    )

    # 3) Shoreline change rate comparison
    print("\nCreating shoreline change rate comparison...")
    plot_total_shoreline_change_over_run(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        units=units,
        start_year=START_YEAR,
        end_year=END_YEAR,
        left_buffer_end=LEFT_BUFFER_END,
        right_buffer_start=RIGHT_BUFFER_START,
        obs_domains=obs_domains,
        obs_change=obs_change,
    )

    print(f"\n{'=' * 70}")
    print(f"ANALYSIS COMPLETE")
    print(f"All outputs saved to: {run_outdir}")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()