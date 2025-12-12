"""
Export and visualize overwash flux (Qow) by domain and year from a CASCADE run.

Functionality:
- Load a saved CASCADE object from disk
- Export overwash flux (Qow) per domain and year to CSV (tidy format)
- Build a 2D matrix (domains × years) of overwash flux
- Plot:
    (1) Heatmap: domain vs. year
    (2) Alongshore profile for a single year (optional)
    (3) Time series for a single domain (optional)

Organization:
- Automatically creates an output subfolder based on the .npz file name
  e.g., if input = HAT_1978_1997_Natural_State.npz
        output folder = <OUTPUT_ROOT_DIR>/HAT_1978_1997_Natural_State/
- All CSVs and figures are saved into that folder.
"""

import csv
import pickle
from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# CONFIG
# =============================================================================

# Path to your saved CASCADE object (.npz for CASCADE)
PATH_TO_CASCADE_FILE = (
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
    r"\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"
)

# Root directory where you want all overwash analyses to live.
# A subfolder named after the npz file will be created inside this.
OUTPUT_ROOT_DIR = (
    r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\analyze_outputs\overwash_analysis"
)

# Number of years (time steps) to export, starting from year index 0
NUM_YEARS = 19

# Domain indices to include (Python indices, inclusive)
FIRST_DOMAIN_INDEX = 15
LAST_DOMAIN_INDEX = 105

# X-axis tick spacing for year indices on plots (set to 2 or 5, etc.)
X_TICK_INTERVAL = 2


# =============================================================================
# CORE EXPORT FUNCTION
# =============================================================================

def export_overwash_flux_by_domain_year(
    cascade,
    output_file: Path,
    num_years: int,
    first_domain: int,
    last_domain: int,
) -> None:
    """
    Export overwash flux (Qow) for each domain and year to a CSV.

    Parameters
    ----------
    cascade : Cascade-like
        Loaded Cascade object with a `barrier3d` attribute (list of domain models).
    output_file : Path
        Path to the CSV file to write.
    num_years : int
        Number of years (time steps) to export, starting at year index 0.
    first_domain : int
        First domain index (Python index) to include, e.g., 1.
    last_domain : int
        Last domain index (Python index) to include, inclusive, e.g., 71.
    """
    barrier_models = cascade.barrier3d

    output_rows = []

    # Loop over domains first (so CSV is grouped by domain)
    for domain_idx, domain_segment in enumerate(barrier_models):

        # Skip domains outside the chosen alongshore range
        if not (first_domain <= domain_idx <= last_domain):
            continue

        # Loop over years for this domain
        for year in range(num_years):

            # Safeguard: only read Qow if this time index exists
            if year < len(domain_segment._QowTS):
                overwash_flux = domain_segment._QowTS[year]
            else:
                overwash_flux = None

            output_rows.append(
                {
                    "domain_index": domain_idx,  # model index (0-based)
                    "year_index": year,          # time step index
                    "overwash_flux": overwash_flux,
                }
            )

    # Write CSV
    output_file = Path(output_file)
    with output_file.open("w", newline="") as csvfile:
        writer = csv.DictWriter(
            csvfile,
            fieldnames=["domain_index", "year_index", "overwash_flux"],
        )
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"\nCSV successfully written → {output_file.resolve()}")
    print(f"Total rows: {len(output_rows)}")


# =============================================================================
# MATRIX BUILDING FOR VISUALIZATION
# =============================================================================

def build_overwash_matrix(
    cascade,
    num_years: int,
    first_domain: int,
    last_domain: int,
):
    """
    Build a 2D matrix of overwash flux (Qow) for visualization.

    Returns
    -------
    overwash_matrix : np.ndarray
        2D array with shape (n_domains, num_years).
        Rows correspond to alongshore domains; columns correspond to years.
        Missing values are filled with np.nan.
    domain_indices : list of int
        List of domain indices (Python indices) used in the rows.
    """
    barrier_models = cascade.barrier3d

    domain_indices = list(range(first_domain, last_domain + 1))
    n_domains = len(domain_indices)

    # Initialize with NaN so missing years appear as blank in the heatmap
    overwash_matrix = np.full((n_domains, num_years), np.nan, dtype=float)

    for row_i, dom_idx in enumerate(domain_indices):
        segment = barrier_models[dom_idx]
        ts = segment._QowTS  # time series of overwash flux for this domain

        max_year = min(num_years, len(ts))
        for year in range(max_year):
            overwash_matrix[row_i, year] = ts[year]

    return overwash_matrix, domain_indices


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_overwash_heatmap(
    overwash_matrix: np.ndarray,
    domain_indices,
    title: str,
    save_path: Optional[Path] = None,
):
    """
    Plot a heatmap of overwash flux.

    Axes:
        x-axis: year index (0, 1, 2, ...)
        y-axis: domain index (alongshore position)
        color: overwash flux magnitude
    """
    n_domains, num_years = overwash_matrix.shape

    plt.figure(figsize=(10, 6))
    im = plt.imshow(
        overwash_matrix,
        aspect="auto",
        origin="lower",
        extent=[0, num_years - 1, domain_indices[0], domain_indices[-1]],
    )
    cbar = plt.colorbar(im)
    cbar.set_label("Overwash flux (Qow)")

    # Set x-axis ticks at regular intervals of X_TICK_INTERVAL years
    xticks = np.arange(0, num_years, X_TICK_INTERVAL)
    plt.xticks(xticks)

    plt.xlabel("Year index")
    plt.ylabel("Domain index")
    plt.title(title)
    plt.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Heatmap figure saved → {save_path.resolve()}")

    plt.show()


def plot_overwash_alongshore_for_year(
    overwash_matrix: np.ndarray,
    domain_indices,
    year_index: int,
    save_path: Optional[Path] = None,
):
    """
    Plot overwash flux alongshore (by domain) for a single year.
    """
    if year_index < 0 or year_index >= overwash_matrix.shape[1]:
        raise ValueError(f"year_index must be between 0 and {overwash_matrix.shape[1] - 1}")

    plt.figure(figsize=(8, 4))
    plt.plot(domain_indices, overwash_matrix[:, year_index])
    plt.xlabel("Domain index")
    plt.ylabel("Overwash flux (Qow)")
    plt.title(f"Overwash flux alongshore (year index = {year_index})")
    plt.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Alongshore profile figure saved → {save_path.resolve()}")

    plt.show()


def plot_overwash_timeseries_for_domain(
    overwash_matrix: np.ndarray,
    domain_indices,
    domain_index: int,
    save_path: Optional[Path] = None,
):
    """
    Plot overwash flux through time (by year) for a single domain.
    """
    if domain_index not in domain_indices:
        raise ValueError(f"domain_index {domain_index} not in domain_indices list.")

    row_i = domain_indices.index(domain_index)
    series = overwash_matrix[row_i, :]

    plt.figure(figsize=(8, 4))
    plt.plot(range(series.size), series)

    # Set x-axis ticks for the time series as well
    xticks = np.arange(0, series.size, X_TICK_INTERVAL)
    plt.xticks(xticks)

    plt.xlabel("Year index")
    plt.ylabel("Overwash flux (Qow)")
    plt.title(f"Overwash timeseries (domain index = {domain_index})")
    plt.tight_layout()

    if save_path is not None:
        save_path = Path(save_path)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Timeseries figure saved → {save_path.resolve()}")

    plt.show()


# =============================================================================
# HELPER: LOAD CASCADE FROM FILE
# =============================================================================

def load_cascade_from_file(cascade_path: Path):
    """
    Load a CASCADE object from a .npz (preferred) or fallback to pickle.

    For .npz, this assumes the archive contains an array under the key
    'cascade' that holds the Cascade object (common CASCADE pattern).
    """
    suffix = cascade_path.suffix.lower()

    if suffix == ".npz":
        data = np.load(cascade_path, allow_pickle=True)
        # Typical pattern: np.savez(..., cascade=self)
        if "cascade" in data:
            cascade = data["cascade"].item()
            print("Loaded CASCADE object from npz key 'cascade'.")
            return cascade
        else:
            # Help you debug if the key is different
            available_keys = list(data.keys())
            raise KeyError(
                f"'cascade' key not found in npz file.\n"
                f"Available keys: {available_keys}\n"
                f"Update load_cascade_from_file() to use the correct key."
            )

    # Fallback: try pickle for .pkl or other extensions
    with cascade_path.open("rb") as f:
        cascade = pickle.load(f)
    print("Loaded CASCADE object using pickle.load(...).")
    return cascade


# =============================================================================
# SCRIPT ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    # --- 1. Paths & run name ---
    cascade_path = Path(PATH_TO_CASCADE_FILE)

    if not cascade_path.exists():
        raise FileNotFoundError(
            f"Could not find CASCADE file at:\n{cascade_path}\n"
            f"Update PATH_TO_CASCADE_FILE in this script."
        )

    # Run name from the file, e.g., "HAT_1978_1997_Natural_State"
    run_name = cascade_path.stem

    # Root output directory and run-specific subfolder
    output_root = Path(OUTPUT_ROOT_DIR)
    run_output_dir = output_root / run_name
    run_output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Outputs will be saved in:\n  {run_output_dir.resolve()}")

    # Predefine key output paths
    csv_path = run_output_dir / f"{run_name}_overwash_flux_by_domain_year.csv"
    heatmap_path = run_output_dir / f"{run_name}_overwash_heatmap.png"
    alongshore_profile_path = run_output_dir / f"{run_name}_overwash_alongshore_last_year.png"
    timeseries_path = run_output_dir / f"{run_name}_overwash_timeseries_domain30.png"

    # --- 2. Load the Cascade object ---
    cascade = load_cascade_from_file(cascade_path)

    # --- 3. Export overwash flux to CSV ---
    export_overwash_flux_by_domain_year(
        cascade=cascade,
        output_file=csv_path,
        num_years=NUM_YEARS,
        first_domain=FIRST_DOMAIN_INDEX,
        last_domain=LAST_DOMAIN_INDEX,
    )

    # --- 4. Build matrix for visualization ---
    overwash_matrix, domain_indices = build_overwash_matrix(
        cascade=cascade,
        num_years=NUM_YEARS,
        first_domain=FIRST_DOMAIN_INDEX,
        last_domain=LAST_DOMAIN_INDEX,
    )

    # --- 5. Visualizations ---

    # (a) Heatmap: domains × years
    plot_overwash_heatmap(
        overwash_matrix,
        domain_indices,
        title=f"Overwash flux by domain and year ({run_name}, subset)",
        save_path=heatmap_path,
    )

    # (b) Optional: alongshore profile for a specific year (here, last year)
    # Uncomment if you want this saved automatically:
    # plot_overwash_alongshore_for_year(
    #     overwash_matrix,
    #     domain_indices,
    #     year_index=NUM_YEARS - 1,  # last simulated year
    #     save_path=alongshore_profile_path,
    # )

    # (c) Optional: time series for a specific domain (e.g., domain 30)
    # Uncomment and adjust domain_index if desired:
    # plot_overwash_timeseries_for_domain(
    #     overwash_matrix,
    #     domain_indices,
    #     domain_index=30,
    #     save_path=timeseries_path,
    # )
