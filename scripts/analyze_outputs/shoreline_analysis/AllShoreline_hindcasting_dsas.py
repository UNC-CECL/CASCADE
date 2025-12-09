import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# --- SETTINGS: Update these values only ---
# ============================================================
NPZ_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"

# Root directory where all shoreline snapshot outputs will be stored.
# A subfolder will be created inside this for each run, based on the .npz file name.
OUTPUT_ROOT_DIR = "shoreline_snapshots_dsas"

START_YEAR = 0           # Label for first time step (e.g., 1978)
RELATIVE_CHANGE = False  # True = plot shoreline change relative to t=0, False = absolute position
TO_METERS = True         # True = convert dam → m by multiplying by 10

# Domain indices for buffers
LEFT_BUFFER_END = 14      # last left-buffer CASCADE domain index (0–14 = 15 left buffers)
RIGHT_BUFFER_START = 105  # first right-buffer CASCADE domain index

# Optional: cleaned DSAS per-domain CSV to overlay (set to None to skip)
# DSAS file should have one row per GIS domain (1–90) with a total-change column.
DSAS_DOMAIN_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_CLEAN.csv"
DSAS_DOMAIN_COL = "domain_id"             # GIS domain ID field in CLEAN file (1–90)
DSAS_CHANGE_COL = "obs_total_change_m" # total shoreline change column in CLEAN file (m)
# ============================================================


def extract_run_name(npz_path):
    """Return filename without extension, e.g., 'HAT_1978_1998_Natural_State'."""
    base = os.path.basename(npz_path)
    return os.path.splitext(base)[0]


def load_cascade(npz_path):
    """Load the CASCADE object from a .npz file."""
    data = np.load(npz_path, allow_pickle=True)
    return data["cascade"][0]


def get_x_s_TS(b3d):
    """
    Get shoreline time series from a Barrier3D object.

    CASCADE sometimes stores shoreline as x_s_TS or _x_s_TS.
    """
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    elif hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    else:
        raise AttributeError("No shoreline time series found on Barrier3D object.")


def build_shoreline_matrix(cascade, to_meters=True):
    """
    Build absolute shoreline[t, domain] matrix (no relative change).
    """
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)

    # Time dimension from domain 0
    nt = len(get_x_s_TS(b3d_list[0]))

    shoreline = np.zeros((nt, ndom))

    for j in range(ndom):
        xs = get_x_s_TS(b3d_list[j])  # dam
        shoreline[:, j] = xs

    if to_meters:
        shoreline = shoreline * 10.0  # dam → m

    return shoreline


def plot_shoreline_snapshots_by_domain(
    shoreline_abs,
    run_name,
    outdir,
    start_year=0,
    relative=False,
    units="m",
):
    """
    For each time step (year), plot shoreline vs domain number.
    """
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
    """
    Plot the mean absolute shoreline (across all domains) as a time series.
    """
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
    left_buffer_end=14,
    right_buffer_start=105,
    obs_domains=None,
    obs_change=None,
):
    """
    Plot total shoreline change over the whole run for each domain.

    Bottom x-axis: CASCADE domain index.
    Top x-axis:    GIS domain number (1–90) over the real-island portion.
    """
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline_abs.shape
    domain_numbers = np.arange(ndom)

    # Total change = final - initial (using absolute shoreline)
    total_change = shoreline_abs[-1, :] - shoreline_abs[0, :]

    # Also save to CSV for record / later analysis
    model_csv = os.path.join(outdir, "total_shoreline_change_model.csv")
    pd.DataFrame(
        {"domain": domain_numbers, "model_total_change_m": total_change}
    ).to_csv(model_csv, index=False)
    print(f"Saved model total shoreline change CSV: {model_csv}")

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(domain_numbers, total_change, marker="o", linewidth=1.5, label="Model")

    # Domain ticks every 10
    ax.set_xticks(domain_numbers[::10])
    ax.set_xlabel("Domain number (alongshore index)")

    # Zero-change reference line
    ax.axhline(0, linestyle="--", linewidth=1, alpha=0.7)

    # Shade buffer regions (if indices provided and in range)
    if left_buffer_end is not None and left_buffer_end >= 0:
        ax.axvspan(domain_numbers[0] - 0.5,
                   left_buffer_end + 0.5,
                   alpha=0.1)
    if right_buffer_start is not None and right_buffer_start < ndom:
        ax.axvspan(right_buffer_start - 0.5,
                   domain_numbers[-1] + 0.5,
                   alpha=0.1)

    # Overlay observed DSAS change if provided
    if obs_domains is not None and obs_change is not None:
        ax.plot(obs_domains, obs_change, marker="x", linewidth=1.2, label="Observed")

    # ---------- Secondary x-axis for GIS domain numbers (top) ----------
    # Mapping: GIS = CASCADE - left_buffer_end ; CASCADE = GIS + left_buffer_end
    def cas_to_gis(x):
        return x - left_buffer_end

    def gis_to_cas(x):
        return x + left_buffer_end

    top_ax = ax.secondary_xaxis('top', functions=(cas_to_gis, gis_to_cas))
    # Ticks every 10 GIS domains; adjust if you want every 5 or every 1
    gis_ticks = np.arange(10, 91, 10)
    top_ax.set_xticks(gis_ticks)
    top_ax.set_xlabel("GIS domain number (1–90)")
    top_ax.tick_params(labelsize=8)

    # Labels and title
    ax.set_ylabel(f"Total shoreline change over run ({units})")
    ax.set_title(f"{run_name} – Total Shoreline Change Over Simulation")

    # Optional text labels for regions
    y_top = ax.get_ylim()[1]
    if left_buffer_end is not None and right_buffer_start is not None:
        mid_real = (left_buffer_end + right_buffer_start) / 2
        ax.text(mid_real, y_top * 0.9, "Real island",
                ha="center", va="top", fontsize=9)
        ax.text((domain_numbers[0] + left_buffer_end) / 2,
                y_top * 0.9, "Left buffer",
                ha="center", va="top", fontsize=8)
        ax.text((right_buffer_start + domain_numbers[-1]) / 2,
                y_top * 0.9, "Right buffer",
                ha="center", va="top", fontsize=8)

    if obs_domains is not None and obs_change is not None:
        ax.legend()

    ax.grid(alpha=0.3)
    fig.tight_layout()

    fname = os.path.join(outdir, "total_shoreline_change_over_run.png")
    fig.savefig(fname, dpi=200)
    plt.close(fig)

    print(f"Saved total shoreline change plot: {fname}")


def main():
    run_name = extract_run_name(NPZ_FILE)
    run_outdir = os.path.join(OUTPUT_ROOT_DIR, run_name)

    cascade = load_cascade(NPZ_FILE)

    # Absolute shoreline matrix in desired units
    shoreline_abs = build_shoreline_matrix(
        cascade,
        to_meters=TO_METERS,
    )

    units = "m" if TO_METERS else "dam"

    # Optional: load observed per-domain DSAS change
    obs_domains = None
    obs_change = None
    if DSAS_DOMAIN_CSV is not None:
        dsas = pd.read_csv(DSAS_DOMAIN_CSV)

        # Raw GIS domain IDs (1–90)
        gis_domains = dsas[DSAS_DOMAIN_COL].to_numpy()

        # Map GIS domains (1–90) → CASCADE real domains (15–104)
        # CASCADE index = GIS domain + LEFT_BUFFER_END
        obs_domains_all = gis_domains + LEFT_BUFFER_END
        obs_change_all = dsas[DSAS_CHANGE_COL].to_numpy()

        # Keep only observed points that fall on real-island CASCADE domains
        mask = (obs_domains_all > LEFT_BUFFER_END) & (obs_domains_all < RIGHT_BUFFER_START)
        obs_domains = obs_domains_all[mask]
        obs_change = obs_change_all[mask]

    # 1) Snapshots across domains for each year
    plot_shoreline_snapshots_by_domain(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        relative=RELATIVE_CHANGE,
        units=units,
    )

    # 2) Mean shoreline time series
    plot_mean_shoreline_timeseries(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        units=units,
    )

    # 3) Total shoreline change over run per domain, with optional observed overlay
    plot_total_shoreline_change_over_run(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        units=units,
        left_buffer_end=LEFT_BUFFER_END,
        right_buffer_start=RIGHT_BUFFER_START,
        obs_domains=obs_domains,
        obs_change=obs_change,
    )


if __name__ == "__main__":
    main()
