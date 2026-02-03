import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# --- SETTINGS: Update these values only ---
# ============================================================
NPZ_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_wh2.5_wa0.1\HAT_1978_1997_wh2.5_wa0.1.npz"

# Root directory where all shoreline snapshot outputs will be stored.
OUTPUT_ROOT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\output\shoreline_analysis\HAT_1978_1997_wh2.5_wa0.1"

START_YEAR = 1978  # Label for first time step
END_YEAR = 1997  # Label for last time step
RELATIVE_CHANGE = False  # True = plot shoreline change relative to t=0, False = absolute position
TO_METERS = True  # True = convert dam → m by multiplying by 10

# Domain indices for buffers
LEFT_BUFFER_END = 14  # last left-buffer CASCADE domain index (0–14 = 15 left buffers)
RIGHT_BUFFER_START = 105  # first right-buffer CASCADE domain index

# *** UPDATED: DSAS file now contains annual rate (LRR), not total change ***
DSAS_DOMAIN_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_domain_means_SIMPLE.csv"
DSAS_DOMAIN_COL = "domain_id"  # GIS domain ID field in CLEAN file (1–90)
DSAS_RATE_COL = "annual_rate_m_per_yr"  # *** CHANGED: Now using annual rate directly from LRR ***

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


def verify_data_alignment(shoreline_abs, obs_domains, obs_rate,
                          left_buffer_end, right_buffer_start, start_year, end_year):
    """
    Print comprehensive diagnostic info to verify units, alignment, and orientation.
    
    *** UPDATED: Now expects obs_rate (annual rate in m/yr) instead of obs_change (total change) ***
    """
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

    # *** UPDATED: DSAS comparison using rate (LRR) ***
    if obs_rate is not None and obs_domains is not None:
        # Calculate total change from observed rate for comparison
        obs_total_change = obs_rate * time_span_years

        print(f"\n--- OBSERVED DSAS DATA (from LRR) ---")
        print(f"Number of DSAS observations: {len(obs_rate)}")
        print(f"DSAS rate range: {obs_rate.min():.2f} to {obs_rate.max():.2f} m/yr")
        print(f"DSAS mean rate: {obs_rate.mean():.2f} m/yr")
        print(f"DSAS rate std dev: {obs_rate.std():.2f} m/yr")
        print(f"\nDSAS implied total change range: {obs_total_change.min():.1f} to {obs_total_change.max():.1f} m")
        print(f"DSAS implied mean total change: {obs_total_change.mean():.2f} m")
        print(f"DSAS CASCADE domain range: {obs_domains.min()} to {obs_domains.max()}")

        # Get corresponding CASCADE values at observed locations
        cascade_at_obs = cascade_total_change[obs_domains.astype(int)]
        cascade_rate_at_obs = cascade_rate[obs_domains.astype(int)]

        print(f"\n--- MODEL vs OBSERVED COMPARISON ---")
        print(f"CASCADE at observed locations - mean change: {cascade_at_obs.mean():.2f} m")
        print(f"DSAS observed (implied) - mean change: {obs_total_change.mean():.2f} m")
        print(f"CASCADE at observed locations - mean rate: {cascade_rate_at_obs.mean():.2f} m/yr")
        print(f"DSAS observed - mean rate: {obs_rate.mean():.2f} m/yr")

        ratio = cascade_at_obs.mean() / obs_total_change.mean() if obs_total_change.mean() != 0 else np.inf
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

        # Correlation check (using rates for better comparison)
        correlation = np.corrcoef(cascade_rate_at_obs, obs_rate)[0, 1]
        print(f"\nCorrelation (CASCADE vs DSAS rates): {correlation:.3f}")

        if correlation > 0.7:
            print("✓ GOOD: Strong positive correlation - spatial patterns align")
        elif correlation < -0.7:
            print("⚠️  WARNING: Strong negative correlation - might need FLIP_CASCADE_ALONGSHORE=True")
        else:
            print("⚠️  WARNING: Weak correlation - check alignment or data quality")

        # Sample comparison at a few points
        print(f"\n--- SAMPLE POINT COMPARISON ---")
        print(f"{'GIS Domain':<12} {'CASCADE Idx':<14} {'CASCADE rate':<16} {'DSAS rate':<14} {'Difference':<12}")
        print("-" * 70)
        sample_indices = np.linspace(0, len(obs_domains) - 1, min(10, len(obs_domains)), dtype=int)
        for i in sample_indices:
            gis_dom = obs_domains[i] - left_buffer_end
            casc_idx = obs_domains[i]
            casc_r = cascade_rate_at_obs[i]
            obs_r = obs_rate[i]
            diff = casc_r - obs_r
            print(f"{gis_dom:<12.0f} {casc_idx:<14.0f} {casc_r:<16.2f} {obs_r:<14.2f} {diff:<12.2f}")

    print("=" * 70)


def plot_shoreline_snapshots_by_domain(
        shoreline_abs,
        run_name="run",
        outdir="output",
        start_year=1978,
        relative=False,
        units="m",
):
    """Plot each year's shoreline vs domain number."""
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline_abs.shape
    domain_numbers = np.arange(ndom)

    for t in range(nt):
        year = start_year + t
        if relative:
            y_vals = shoreline_abs[t, :] - shoreline_abs[0, :]
            ylabel = f"Shoreline change since t=0 [{units}]"
            fname_base = f"shoreline_snapshot_relative_year{year}.png"
        else:
            y_vals = shoreline_abs[t, :]
            ylabel = f"Shoreline position [{units}]"
            fname_base = f"shoreline_snapshot_abs_year{year}.png"

        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(domain_numbers, y_vals, marker="o", linewidth=1.5)
        ax.set_xlabel("Domain number (alongshore index)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"{run_name} – Year {year}")
        ax.grid(alpha=0.3)
        fig.tight_layout()

        fname = os.path.join(outdir, fname_base)
        fig.savefig(fname, dpi=150)
        plt.close(fig)

    print(f"Created {nt} shoreline snapshots.")


def plot_mean_shoreline_timeseries(
        shoreline_abs,
        run_name="run",
        outdir="output",
        start_year=1978,
        units="m",
):
    """Plot mean shoreline vs time."""
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline_abs.shape
    mean_shoreline = shoreline_abs.mean(axis=1)
    years = np.arange(start_year, start_year + nt)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(years, mean_shoreline, marker="o", linewidth=1.5)
    ax.set_xlabel("Year")
    ax.set_ylabel(f"Mean shoreline position [{units}]")
    ax.set_title(f"{run_name} – Mean Shoreline Over Time")
    ax.grid(alpha=0.3)
    fig.tight_layout()

    fname = os.path.join(outdir, "mean_shoreline_timeseries.png")
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"Saved mean shoreline time series: {fname}")


def plot_total_shoreline_change_over_run(
        shoreline_abs,
        run_name="run",
        outdir="output",
        units="m",
        start_year=1978,
        end_year=1997,
        left_buffer_end=14,
        right_buffer_start=105,
        obs_domains=None,
        obs_rate=None,  # *** CHANGED: Now expects obs_rate instead of obs_change ***
):
    """
    Plot shoreline change RATE over the whole run for each domain.

    Bottom x-axis: CASCADE domain index.
    Top x-axis:    GIS domain number (1–90) over the real-island portion.
    
    *** UPDATED: Now expects obs_rate (annual rate) instead of obs_change (total change) ***
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

    # *** UPDATED: Overlay observed DSAS rate directly (no conversion needed) ***
    if obs_domains is not None and obs_rate is not None:
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
    ax.set_title("Observed vs. Modeled Shoreline Change Rate (1978-1997)\nHatteras Island, North Carolina") #(f"{run_name} – Shoreline Change Rate ({start_year}–{end_year})")

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

    # *** UPDATED: Load observed LRR (annual rate) instead of total change ***
    obs_domains = None
    obs_rate = None
    if DSAS_DOMAIN_CSV is not None:
        print(f"Loading DSAS observations from: {DSAS_DOMAIN_CSV}")
        dsas = pd.read_csv(DSAS_DOMAIN_CSV)

        # Raw GIS domain IDs (1–90)
        gis_domains = dsas[DSAS_DOMAIN_COL].to_numpy()

        # Map GIS domains (1–90) → CASCADE real domains (15–104)
        obs_domains_all = gis_domains + LEFT_BUFFER_END
        obs_rate_all = dsas[DSAS_RATE_COL].to_numpy()  # *** CHANGED: Now loading rate directly ***

        # Keep only observed points that fall on real-island CASCADE domains
        mask = (obs_domains_all > LEFT_BUFFER_END) & (obs_domains_all < RIGHT_BUFFER_START)
        obs_domains = obs_domains_all[mask]
        obs_rate = obs_rate_all[mask]  # *** CHANGED: Now using rate ***

        print(f"Loaded {len(obs_rate)} DSAS observations (annual rates) for real island domains\n")

    # DIAGNOSTIC VERIFICATION
    verify_data_alignment(
        shoreline_abs,
        obs_domains,
        obs_rate,  # *** CHANGED: Passing rate instead of change ***
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
        obs_rate=obs_rate,  # *** CHANGED: Passing rate instead of change ***
    )

    print(f"\n{'=' * 70}")
    print(f"ANALYSIS COMPLETE")
    print(f"All outputs saved to: {run_outdir}")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
