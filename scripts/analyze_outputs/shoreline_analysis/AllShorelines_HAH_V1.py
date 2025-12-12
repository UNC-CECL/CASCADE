import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# --- SETTINGS: Update these values only ---
# ============================================================
NPZ_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"

# Root directory where all shoreline snapshot outputs will be stored.
# A subfolder will be created inside this for each run, based on the .npz file name.
OUTPUT_ROOT_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\analyze_outputs\shoreline_analysis\AllShorelines_HAH_V1"

START_YEAR = 0           # Label for first time step (e.g., 1978)
RELATIVE_CHANGE = False  # True = plot shoreline change relative to t=0, False = absolute position
TO_METERS = True         # True = convert dam → m by multiplying by 10

# Domain indices for buffers (edit if needed)
LEFT_BUFFER_END = 14      # last left-buffer domain index
RIGHT_BUFFER_START = 105  # first right-buffer domain index
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

    Parameters
    ----------
    cascade : Cascade object
        Loaded from the .npz file.
    to_meters : bool
        If True, convert dam → m by multiplying by 10.
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

    Parameters
    ----------
    shoreline_abs : ndarray
        Absolute shoreline array of shape (nt, ndom), shoreline[t, domain].
    run_name : str
        Name derived from the .npz file for labeling and organization.
    outdir : str
        Output directory for PNGs.
    start_year : int
        Label for starting year.
    relative : bool
        If True, plot shoreline change relative to t=0 (per domain).
    units : str
        Units string for y-axis label (e.g., 'm').
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

        # Domain number ticks every 10 domains
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

    Parameters
    ----------
    shoreline_abs : ndarray
        Absolute shoreline array of shape (nt, ndom).
    run_name : str
        Name derived from the .npz file.
    outdir : str
        Output directory for the PNG.
    start_year : int
        Label for starting year.
    units : str
        Units string for y-axis label.
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
):
    """
    Plot total shoreline change over the whole run for each domain.

    X-axis: domain number
    Y-axis: x_final - x_initial (total change over time period)

    Includes:
    - Zero-change reference line
    - Shaded buffer regions
    - Text labels for buffers and real island
    """
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline_abs.shape
    domain_numbers = np.arange(ndom)

    # Total change = final - initial (using absolute shoreline)
    total_change = shoreline_abs[-1, :] - shoreline_abs[0, :]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(domain_numbers, total_change, marker="o", linewidth=1.5)

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

    # Labels and title
    ax.set_ylabel(f"Total shoreline change over run ({units})")
    ax.set_title(f"{run_name} – Total Shoreline Change Over Simulation")

    # Optional text labels for regions (placed near top of plot)
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

    # 1) Snapshots across domains for each year (absolute or relative depending on RELATIVE_CHANGE)
    plot_shoreline_snapshots_by_domain(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        relative=RELATIVE_CHANGE,
        units=units,
    )

    # 2) Mean shoreline time series over whole island (absolute)
    plot_mean_shoreline_timeseries(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        units=units,
    )

    # 3) Total shoreline change over run per domain (final - initial),
    #    with buffers shaded and labeled
    plot_total_shoreline_change_over_run(
        shoreline_abs,
        run_name=run_name,
        outdir=run_outdir,
        units=units,
        left_buffer_end=LEFT_BUFFER_END,
        right_buffer_start=RIGHT_BUFFER_START,
    )


if __name__ == "__main__":
    main()
