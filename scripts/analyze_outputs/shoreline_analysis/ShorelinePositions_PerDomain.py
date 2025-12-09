import os
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# --- SETTINGS: Update these two lines only ---
# ============================================================
NPZ_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1998_Natural_State\HAT_1978_1998_Natural_State.npz"

# Root directory where all shoreline outputs will be stored.
# A subfolder will be created inside this for each run, based on the .npz file name.
OUTPUT_ROOT_DIR = "shoreline_outputs"

START_YEAR = 0
RELATIVE_CHANGE = False   # Set True to plot x_s relative to first value
# ============================================================


def extract_run_name(npz_path):
    """Return filename without extension, e.g., 'HAT_1978_1998_Natural_State'."""
    base = os.path.basename(npz_path)
    return os.path.splitext(base)[0]


def load_cascade(npz_path):
    data = np.load(npz_path, allow_pickle=True)
    return data["cascade"][0]


def get_x_s_TS(b3d):
    """CASCADE sometimes stores shoreline as x_s_TS or _x_s_TS."""
    if hasattr(b3d, "x_s_TS"):
        return np.array(b3d.x_s_TS)
    if hasattr(b3d, "_x_s_TS"):
        return np.array(b3d._x_s_TS)
    raise AttributeError("No shoreline time series found on Barrier3D object.")


def plot_shoreline_each_domain(cascade, run_name, outdir, start_year=0, relative=False):
    """
    Plot shoreline position (or change) for each domain and save PNGs.

    Parameters
    ----------
    cascade : Cascade object
        Loaded from the .npz file.
    run_name : str
        Name of the run, typically derived from the .npz filename.
    outdir : str
        Directory where output PNGs will be saved.
    start_year : int
        Starting year for the time axis.
    relative : bool
        If True, plot shoreline change relative to the first time step.
    """
    b3d_list = cascade.barrier3d
    ndom = len(b3d_list)

    os.makedirs(outdir, exist_ok=True)

    for i in range(ndom):
        b3d = b3d_list[i]
        x_s = get_x_s_TS(b3d)

        if relative:
            x_s = x_s - x_s[0]

        years = np.arange(len(x_s)) + start_year

        plt.figure(figsize=(10, 4))
        plt.plot(years, x_s, linewidth=1.5)
        plt.grid(alpha=0.3)
        plt.title(f"{run_name} – Shoreline Position – Domain {i}")
        plt.xlabel("Year")
        plt.ylabel("Shoreline (dam)" if not relative else "Shoreline Change (dam)")
        plt.tight_layout()

        # Inside a run-specific folder, filenames can just be by domain
        fname = os.path.join(outdir, f"domain_{i}.png")
        plt.savefig(fname, dpi=200)
        plt.close()

        print(f"Saved: {fname}")


def main():
    # Derive run name and run-specific output directory
    run_name = extract_run_name(NPZ_FILE)
    run_outdir = os.path.join(OUTPUT_ROOT_DIR, run_name)

    cascade = load_cascade(NPZ_FILE)

    plot_shoreline_each_domain(
        cascade,
        run_name=run_name,
        outdir=run_outdir,
        start_year=START_YEAR,
        relative=RELATIVE_CHANGE
    )


if __name__ == "__main__":
    main()
