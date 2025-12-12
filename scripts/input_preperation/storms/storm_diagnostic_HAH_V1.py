"""
CASCADE Storm–Dune–Overwash Diagnostic Dashboard (HAH v1)

This script helps you understand why overwash is (or is not) happening
by visualizing, for selected domains:

  - Total Water Level (TWL), Water Level (WL), and Hs through time
  - Dune crest elevation
  - Overwash flux (Qow) through time
  - Optional: per-year max Rhigh from the storm file

Assumptions:
- CASCADE object saved in a .npz with key 'cascade'
- Each barrier3d segment may contain attributes like:
    TWL_TS, WL_TS, Hs_TS, _QowTS, z_dune (or similar)
- Storm file is a NumPy array with columns:
    [Year_Index, Rhigh, Rlow, Wave Period, Duration]
"""

from pathlib import Path
from typing import Optional, List, Dict

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd

# =============================================================================
# CONFIG: UPDATE THESE FOR YOUR RUN
# =============================================================================

# 1. Paths
PATH_TO_CASCADE_FILE = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
    r"\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"
)

# OPTIONAL: path to your storm .npy file (set to None if you don't want to use it)
PATH_TO_STORM_FILE: Optional[Path] = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms"
    r"\storms_1978_1997.npy"
)

# 2. Hindcast period
HINDCAST_START_YEAR = 1978   # first model year
NUM_YEARS = 20               # number of time steps in this run (1978–1997 ⇒ 20)

# 3. Domains to inspect in detail
DOMAINS_TO_PLOT: List[int] = [30, 40, 50, 60]

# 4. Storm file column names (must match your storm npy structure)
STORM_COLS = ["Year_Index", "Rhigh", "Rlow", "Wave Period", "Duration"]

# 5. Whether to overlay per-year max Rhigh from storm file on the Qow panel
PLOT_RHIGH_ON_QOW = True


# =============================================================================
# HELPERS
# =============================================================================

def load_cascade_from_file(cascade_path: Path):
    """Load a CASCADE object from a .npz (preferred) or pickle."""
    if not cascade_path.exists():
        raise FileNotFoundError(f"CASCADE file not found at: {cascade_path}")

    suffix = cascade_path.suffix.lower()

    if suffix == ".npz":
        data = np.load(cascade_path, allow_pickle=True)
        if "cascade" in data:
            cascade = data["cascade"].item()
            print("Loaded CASCADE object from npz key 'cascade'.")
            return cascade
        else:
            keys = list(data.keys())
            raise KeyError(
                f"'cascade' key not found in npz file.\n"
                f"Available keys: {keys}"
            )

    # Fallback to pickle
    with cascade_path.open("rb") as f:
        cascade = pickle.load(f)
    print("Loaded CASCADE object using pickle.load(...).")
    return cascade


def load_storm_dataframe(storm_path: Optional[Path]) -> Optional[pd.DataFrame]:
    """Load the storm .npy file into a DataFrame with Calendar_Year."""
    if storm_path is None:
        return None
    if not storm_path.exists():
        print(f"Storm file not found at: {storm_path}")
        return None

    arr = np.load(storm_path)
    df = pd.DataFrame(arr, columns=STORM_COLS)
    df["Calendar_Year"] = df["Year_Index"].astype(int) + HINDCAST_START_YEAR

    print(f"\nLoaded storm file: {storm_path}")
    print(f"Storm array shape: {arr.shape}")
    print(df[["Calendar_Year", "Rhigh", "Rlow"]].head())

    return df


def get_dune_crest(segment) -> Optional[float]:
    """
    Try common attribute names for dune crest elevation.

    Adjust this list if your CASCADE object uses a different name.
    """
    candidate_attrs = [
        "z_dune", "dune_crest_elev", "z_crest", "z_dune_crest"
    ]
    for attr in candidate_attrs:
        if hasattr(segment, attr):
            return float(getattr(segment, attr))
    return None


def get_time_series(segment) -> Dict[str, np.ndarray]:
    """
    Safely extract TWL, WL, Hs, and Qow time series (if they exist) from a segment.
    Returns a dict like {"TWL": ..., "WL": ..., "Hs": ..., "Qow": ...}
    """
    ts_dict: Dict[str, np.ndarray] = {}

    def grab(attr_name: str, key_name: str):
        if hasattr(segment, attr_name):
            try:
                ts = np.array(getattr(segment, attr_name), dtype=float)
                ts_dict[key_name] = ts
            except Exception as e:
                print(f"  Could not convert {attr_name} to float array: {e}")

    grab("TWL_TS", "TWL")
    grab("WL_TS", "WL")
    grab("Hs_TS", "Hs")
    grab("_QowTS", "Qow")

    return ts_dict


def compute_year_axis(n_steps: int) -> np.ndarray:
    """Return calendar year array for n_steps based on HINDCAST_START_YEAR."""
    return HINDCAST_START_YEAR + np.arange(n_steps)


def compute_yearly_rhigh(storm_df: Optional[pd.DataFrame]) -> Optional[pd.Series]:
    """
    Compute per-year max Rhigh for the hindcast period.
    Returns a Series indexed by Calendar_Year.
    """
    if storm_df is None:
        return None

    # Restrict to hindcast years
    end_year = HINDCAST_START_YEAR + NUM_YEARS - 1
    mask = (storm_df["Calendar_Year"] >= HINDCAST_START_YEAR) & \
           (storm_df["Calendar_Year"] <= end_year)
    subset = storm_df.loc[mask]

    if subset.empty:
        print("No storms within hindcast years in the storm file.")
        return None

    yearly_max = subset.groupby("Calendar_Year")["Rhigh"].max()
    return yearly_max


# =============================================================================
# DASHBOARD PLOTTING
# =============================================================================

def plot_domain_dashboard(
    dom_index: int,
    segment,
    storm_df: Optional[pd.DataFrame],
):
    """
    Make a 2-panel figure for one domain:
      Panel 1: TWL / WL / Hs + dune crest
      Panel 2: Qow + optional per-year max Rhigh
    """
    print(f"\n=== DOMAIN {dom_index} DIAGNOSTICS ===")
    dune_crest = get_dune_crest(segment)
    ts_dict = get_time_series(segment)

    # Figure out how many time steps we have from whichever series is longest
    if not ts_dict:
        print("  No time series (TWL/WL/Hs/Qow) found on this segment.")
        return

    max_len = max(len(ts) for ts in ts_dict.values())
    cal_years = compute_year_axis(max_len)

    # Print summary info
    if dune_crest is not None:
        print(f"  Dune crest elevation: {dune_crest:.2f} m")
    else:
        print("  Dune crest elevation: NOT FOUND (check attribute names).")

    for name, ts in ts_dict.items():
        nonzero = np.count_nonzero(np.isfinite(ts) & (ts != 0))
        print(f"  {name}: min={np.nanmin(ts):.3f}, max={np.nanmax(ts):.3f}, "
              f"nonzero={nonzero}")

    # --- Figure layout ---
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(11, 8), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]}
    )

    # -----------------------
    # Panel 1: Forcing vs dune
    # -----------------------
    if "TWL" in ts_dict:
        twl = ts_dict["TWL"]
        ax1.plot(
            cal_years[:len(twl)],
            twl,
            label="TWL",
            linewidth=1.5,
        )
    if "WL" in ts_dict:
        wl = ts_dict["WL"]
        ax1.plot(
            cal_years[:len(wl)],
            wl,
            label="WL",
            linestyle="--",
            linewidth=1.2,
        )
    if "Hs" in ts_dict:
        hs = ts_dict["Hs"]
        # Optional: rescale Hs so it fits visually
        hs_scaled = hs * 0.5
        ax1.plot(
            cal_years[:len(hs)],
            hs_scaled,
            label="0.5×Hs (scaled)",
            linestyle=":",
            linewidth=1.2,
        )

    if dune_crest is not None:
        ax1.axhline(
            dune_crest,
            color="gray",
            linestyle="--",
            label=f"Dune crest ({dune_crest:.2f} m)",
        )

    ax1.set_ylabel("Water Level / Hs (m)")
    ax1.set_title(f"Domain {dom_index}: Storm Forcing vs Dune Crest")
    ax1.grid(alpha=0.3)
    ax1.legend(loc="upper right", fontsize=9)

    # -----------------------
    # Panel 2: Overwash vs Rhigh
    # -----------------------
    if "Qow" in ts_dict:
        qow = ts_dict["Qow"]
        ax2.plot(
            cal_years[:len(qow)],
            qow,
            label="Qow",
            color="tab:purple",
            linewidth=1.5,
        )
    else:
        ax2.text(
            0.5, 0.5,
            "No Qow time series found",
            transform=ax2.transAxes,
            ha="center",
            va="center",
        )
        ax2.set_ylabel("Qow")
        ax2.set_xlabel("Calendar Year")
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        plt.show()
        return

    ax2.set_ylabel("Overwash flux (Qow)")
    ax2.grid(alpha=0.3)

    # Optional: overlay per-year max Rhigh on a secondary y-axis
    if PLOT_RHIGH_ON_QOW and storm_df is not None:
        yearly_max_rhigh = compute_yearly_rhigh(storm_df)
        if yearly_max_rhigh is not None and not yearly_max_rhigh.empty:
            ax3 = ax2.twinx()
            years_rhigh = yearly_max_rhigh.index.values
            ax3.plot(
                years_rhigh,
                yearly_max_rhigh.values,
                color="tab:orange",
                marker="o",
                linestyle="-",
                label="Max Rhigh per year",
            )
            ax3.set_ylabel("Rhigh (storm intensity)")
            # Combine legends
            lines2, labels2 = ax2.get_legend_handles_labels()
            lines3, labels3 = ax3.get_legend_handles_labels()
            ax3.legend(lines2 + lines3, labels2 + labels3,
                       loc="upper right", fontsize=9)
        else:
            ax2.legend(loc="upper right", fontsize=9)
    else:
        ax2.legend(loc="upper right", fontsize=9)

    ax2.set_xlabel("Calendar Year")
    fig.suptitle(
        f"CASCADE Storm–Dune–Overwash Diagnostics (Domain {dom_index})",
        fontsize=14,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    # Load CASCADE
    cascade = load_cascade_from_file(PATH_TO_CASCADE_FILE)
    barrier_models = cascade.barrier3d

    # Load storms (optional)
    storm_df = load_storm_dataframe(PATH_TO_STORM_FILE)

    # Safety check
    max_dom_index = len(barrier_models) - 1
    print(f"\nCASCADE has {len(barrier_models)} alongshore domains "
          f"(0–{max_dom_index}).")

    for dom in DOMAINS_TO_PLOT:
        if dom < 0 or dom > max_dom_index:
            print(f"\nDomain {dom} is out of range, skipping.")
            continue
        seg = barrier_models[dom]
        plot_domain_dashboard(dom, seg, storm_df)
