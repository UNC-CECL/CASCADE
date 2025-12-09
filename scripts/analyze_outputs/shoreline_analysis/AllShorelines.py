import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# --- SETTINGS: Update these values only ---
# ============================================================
NPZ_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1998_Natural_State\HAT_1978_1998_Natural_State.npz"

# Root directory where all shoreline snapshot outputs will be stored.
# A subfolder will be created inside this for each run, based on the .npz file name.
OUTPUT_ROOT_DIR = "shoreline_snapshots"

START_YEAR = 0           # Label for first time step (e.g., 1978)
RELATIVE_CHANGE = False  # True = shoreline change relative to t=0, False = absolute position
TO_METERS = True         # True = convert dam → m by multiplying by 10
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


def build_shoreline_matrix(cascade, to_meters=True, relative=False):
    """
    Build shoreline[t, domain] matrix.

    Parameters
    ----------
    cascade : Cascade object
        Loaded from the .npz file.
    to_meters : bool
        If True, convert dam → m by multiplying by 10.
    relative : bool
        If True, subtract t=0 shoreline from all time steps (per domain).
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

    if relative:
        # Subtract the initial shoreline position at t=0 (per domain)
        shoreline = shoreline - shoreline[0, :]

    return shoreline


def plot_shoreline_snapshots_by_domain(
    shoreline,
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
    shoreline : ndarray
        Array of shape (nt, ndom), shoreline[t, domain].
    run_name : str
        Name derived from the .npz file for labeling and organization.
    outdir : str
        Output directory for PNGs.
    start_year : int
        Label for starting year.
    relative : bool
        If True, y-axis is 'change' instead of absolute position.
    units : str
        Units string for y-axis label (e.g., 'm').
    """
    os.makedirs(outdir, exist_ok=True)

    nt, ndom = shoreline.shape
    domain_numbers = np.arange(ndom)

    for t in range(nt):
        year = start_year + t

        plt.figure(figsize=(12, 4))
        plt.plot(domain_numbers, shoreline[t, :], marker="o", linewidth=1.2)

        # --- Domain number ticks every 10 domains ---
        plt.xticks(domain_numbers[::10])
        plt.xlabel("Domain number (alongshore index)")

        # Y-axis label + title
        if relative:
            plt.ylabel(f"Shoreline change ({units})")
            plt.title(f"{run_name} – Shoreline Change Across Domains – Year {year}")
        else:
            plt.ylabel(f"Shoreline position ({units})")
            plt.title(f"{run_name} – Shoreline Position Across Domains – Year {year}")

        plt.grid(alpha=0.3)
        plt.tight_layout()

        # Output filename
        fname = os.path.join(outdir, f"snapshot_year_{year}.png")
        plt.savefig(fname, dpi=200)
        plt.close()

        print(f"Saved: {fname}")


def main():
    run_name = extract_run_name(NPZ_FILE)
    run_outdir = os.path.join(OUTPUT_ROOT_DIR, run_name)

    cascade = load_cascade(NPZ_FILE)

    shoreline = build_shoreline_matrix(
        cascade,
        to_meters=TO_METERS,
        relative=RELATIVE_CHANGE,
    )

    units = "m" if TO_METERS else "dam"

    plot_shoreline_snapshots_by_domain(
        shoreline,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        relative=RELATIVE_CHANGE,
        units=units,
    )


if __name__ == "__main__":
    main()
