"""
Domain Dashboard for a Single Barrier3D Segment

What this does:
    - Loads a CASCADE object from a .npz (or pickle).
    - Selects one alongshore domain: cascade.barrier3d[DOMAIN_INDEX].
    - Computes:
        * Dune crest elevation time series (from DuneDomain)
        * Overwash flux time series (QowTS / _QowTS)
        * Shoreline position time series (x_s_TS / _x_s_TS)
        * Storm intensity (Rhigh) per storm from StormSeries / _StormSeries
        * Freeboard per year = max(Rhigh) - dune crest
    - Plots a "dashboard" with:
        (1) Dune crest elevation vs year
        (2) Overwash flux vs year + storm intensity + freeboard
        (3) Shoreline position vs year

Assumptions:
    - DuneDomain has shape (time, cross_shore, 2) and field 0 is dune elevation.
    - StormSeries columns are [Year_Index, Rhigh, Rlow, Wave Period, Duration]
    - All vertical & cross-shore quantities in CASCADE are in decameters (dam)
    - Time index 0 → START_YEAR (e.g., 1978)
"""

from pathlib import Path
import numpy as np
import pickle
import matplotlib.pyplot as plt

# ==========================
# CONFIG: UPDATE THESE
# ==========================

# Path to your saved CASCADE object (.npz preferred)
PATH_TO_CASCADE_FILE = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
    r"\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"
)

# Domain index to inspect (0-based)
DOMAIN_INDEX = 40

# Start calendar year for the hindcast (used to convert year index → year)
START_YEAR = 1978

# Whether to label x-axis with calendar years instead of year indices
USE_CALENDAR_YEARS = True

# CASCADE units: decameters → meters
DAM_TO_M = 10.0


# ==========================
# HELPERS
# ==========================

def load_cascade_from_file(cascade_path: Path):
    """Load a CASCADE object from a .npz with key 'cascade', or pickle as fallback."""
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

    # Fallback: try pickle
    with cascade_path.open("rb") as f:
        cascade = pickle.load(f)
    print("Loaded CASCADE object using pickle.load(...).")
    return cascade


def get_dune_crest_ts(segment):
    """
    Compute dune crest elevation time series for a single Barrier3D segment.

    Uses segment.DuneDomain, assumed shape (time, cross_shore, 2),
    and treats field 0 as dune elevation (field 1 is identical in your run).

    Returns
    -------
    crest_ts_m : np.ndarray of shape (T,)
        Max dune elevation at each time step, in meters.
    """
    if not hasattr(segment, "DuneDomain"):
        raise AttributeError("Segment has no DuneDomain attribute.")

    dune = np.array(segment.DuneDomain, dtype=float)  # (time, x, 2)

    if dune.ndim != 3 or dune.shape[2] < 1:
        raise ValueError(f"Unexpected DuneDomain shape: {dune.shape}")

    dune_elev_dam = dune[:, :, 0]               # decameters
    crest_ts_dam = np.nanmax(dune_elev_dam, axis=1)
    crest_ts_m = crest_ts_dam * DAM_TO_M        # convert to meters

    return crest_ts_m


def get_qow_ts(segment):
    """
    Get overwash flux time series for a segment, trying both QowTS and _QowTS.
    Units depend on CASCADE config; we leave them as-is.
    """
    if hasattr(segment, "QowTS"):
        qow = np.array(segment.QowTS, dtype=float)
    elif hasattr(segment, "_QowTS"):
        qow = np.array(segment._QowTS, dtype=float)
    else:
        raise AttributeError("Segment has neither QowTS nor _QowTS.")
    return qow


def get_shoreline_ts(segment):
    """
    Get shoreline position time series (x_s_TS or _x_s_TS), in meters.
    """
    if hasattr(segment, "x_s_TS"):
        xs_dam = np.array(segment.x_s_TS, dtype=float)
    elif hasattr(segment, "_x_s_TS"):
        xs_dam = np.array(segment._x_s_TS, dtype=float)
    else:
        raise AttributeError("Segment has neither x_s_TS nor _x_s_TS.")

    xs_m = xs_dam * DAM_TO_M
    return xs_m


def get_storm_series(segment):
    """
    Get StormSeries array for this segment (StormSeries or _StormSeries).

    Returns
    -------
    storms : np.ndarray of shape (n_storms, 5)
        Columns assumed to be [Year_Index, Rhigh, Rlow, Wave Period, Duration],
        still in decameters for the vertical quantities.
    """
    if hasattr(segment, "StormSeries"):
        storms = np.array(segment.StormSeries, dtype=float)
    elif hasattr(segment, "_StormSeries"):
        storms = np.array(segment._StormSeries, dtype=float)
    else:
        raise AttributeError("Segment has neither StormSeries nor _StormSeries.")
    return storms


def compute_year_axis(n_t, start_year, use_calendar):
    """
    Build x-axis for plotting.
    """
    year_idx = np.arange(n_t)
    if use_calendar:
        return start_year + year_idx
    else:
        return year_idx


# ==========================
# MAIN
# ==========================

if __name__ == "__main__":
    # --- Load CASCADE object ---
    cascade = load_cascade_from_file(PATH_TO_CASCADE_FILE)

    barrier_models = cascade.barrier3d
    n_domains = len(barrier_models)
    print(f"CASCADE has {n_domains} alongshore domains (0–{n_domains - 1}).")

    if not (0 <= DOMAIN_INDEX < n_domains):
        raise IndexError(
            f"DOMAIN_INDEX {DOMAIN_INDEX} is out of range. "
            f"Must be between 0 and {n_domains - 1}."
        )

    seg = barrier_models[DOMAIN_INDEX]
    print(f"Building dashboard for barrier3d[{DOMAIN_INDEX}]")

    # --- Core time series ---
    crest_ts_m = get_dune_crest_ts(seg)  # meters
    qow_ts = get_qow_ts(seg)             # model units
    xs_ts_m = get_shoreline_ts(seg)      # meters

    # Make sure lengths line up (trim to common length if needed)
    n_t = min(len(crest_ts_m), len(qow_ts), len(xs_ts_m))
    crest_ts_m = crest_ts_m[:n_t]
    qow_ts = qow_ts[:n_t]
    xs_ts_m = xs_ts_m[:n_t]

    x_axis = compute_year_axis(n_t, START_YEAR, USE_CALENDAR_YEARS)

    # --- Storm series (per-event Rhigh, grouped by year index) ---
    storms = get_storm_series(seg)  # shape (n_storms, 5) expected
    if storms.size > 0:
        storm_year_idx = storms[:, 0].astype(int)
        storm_rhigh_dam = storms[:, 1]          # decameters
        storm_rhigh_m = storm_rhigh_dam * DAM_TO_M

        # VALIDATE: keep only storms within [0, n_t-1]
        valid_event_mask = (storm_year_idx >= 0) & (storm_year_idx < n_t)
        storm_year_idx_valid = storm_year_idx[valid_event_mask]
        storm_rhigh_valid_m = storm_rhigh_m[valid_event_mask]

        # Aggregate: max Rhigh per year index (meters)
        max_rhigh_by_year_m = {}
        for y_idx, rh_m in zip(storm_year_idx_valid, storm_rhigh_valid_m):
            if y_idx not in max_rhigh_by_year_m:
                max_rhigh_by_year_m[y_idx] = rh_m
            else:
                max_rhigh_by_year_m[y_idx] = max(max_rhigh_by_year_m[y_idx], rh_m)

        # Build arrays aligned with model years
        yearly_rhigh_m = np.full(n_t, np.nan)
        for yi, rh_m in max_rhigh_by_year_m.items():
            yearly_rhigh_m[yi] = rh_m

        storm_x = compute_year_axis(n_t, START_YEAR, USE_CALENDAR_YEARS)

        # --- Freeboard: Rhigh_max - dune crest (m) ---
        # Only meaningful where we actually have storms
        freeboard_m = yearly_rhigh_m - crest_ts_m
    else:
        storm_year_idx_valid = np.array([], dtype=int)
        storm_rhigh_valid_m = np.array([], dtype=float)
        yearly_rhigh_m = np.full(n_t, np.nan)
        freeboard_m = np.full(n_t, np.nan)
        storm_x = x_axis

    # --- Print a quick text summary ---
    print("\n=== Domain Dashboard Summary ===")
    print(f"  Domain index: {DOMAIN_INDEX}")
    print(f"  Time steps (TMAX): {n_t}")
    print(f"  Dune crest elev (m): min={np.nanmin(crest_ts_m):.3f}, max={np.nanmax(crest_ts_m):.3f}")
    print(f"  Shoreline x_s (m):   min={np.nanmin(xs_ts_m):.3f}, max={np.nanmax(xs_ts_m):.3f}")
    print(f"  QowTS:               min={np.nanmin(qow_ts):.3e}, max={np.nanmax(qow_ts):.3e}")
    if storms.size > 0 and storm_year_idx_valid.size > 0:
        print(f"  Storm events (raw):  {storms.shape[0]} total")
        print(f"  Storm events (used): {storm_year_idx_valid.size} within model time")
        print(f"  Rhigh (used, m):     min={np.nanmin(storm_rhigh_valid_m):.3f}, "
              f"max={np.nanmax(storm_rhigh_valid_m):.3f}")
        print(f"  Freeboard (m):       min={np.nanmin(freeboard_m):.3f}, "
              f"max={np.nanmax(freeboard_m):.3f}")
    else:
        print("  Storm events:        none found or none within model time")

    # --- Plot dashboard ---
    fig, axs = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    fig.suptitle(f"Domain {DOMAIN_INDEX} Dashboard", fontsize=16)

    # Panel 1: Dune crest elevation
    axs[0].plot(x_axis, crest_ts_m, marker="o")
    axs[0].set_ylabel("Dune crest (m)")
    axs[0].grid(True, linestyle="--", alpha=0.5)
    axs[0].set_title("Dune crest elevation vs time")

    # Panel 2: Qow + storms + freeboard
    axs[1].bar(x_axis, qow_ts, width=0.6, alpha=0.7, label="Qow")
    axs[1].set_ylabel("Overwash flux (Qow)")
    axs[1].grid(True, linestyle="--", alpha=0.5)

    # Twin axis for storm intensity + freeboard
    ax2b = axs[1].twinx()
    if storm_year_idx_valid.size > 0:
        # Individual storms as scatter (gray)
        ax2b.scatter(
            storm_x[storm_year_idx_valid],
            storm_rhigh_valid_m,
            s=25,
            alpha=0.4,
            color="gray",
            label="Rhigh (events)"
        )
        # Yearly max Rhigh as line
        ax2b.plot(
            storm_x,
            yearly_rhigh_m,
            "-o",
            linewidth=2,
            alpha=0.9,
            color="tab:blue",
            label="Max Rhigh per year"
        )
        # Freeboard as line
        ax2b.plot(
            storm_x,
            freeboard_m,
            "-s",
            linewidth=2,
            alpha=0.9,
            color="tab:orange",
            label="Freeboard (Rhigh - crest)"
        )
    ax2b.set_ylabel("Storm intensity / freeboard (m)")
    axs[1].set_title("Overwash flux, storm intensity, and freeboard")

    # Handle legend for twin axes
    lines1, labels1 = axs[1].get_legend_handles_labels()
    lines2, labels2 = ax2b.get_legend_handles_labels()
    if lines1 or lines2:
        axs[1].legend(lines1 + lines2, labels1 + labels2, loc="upper right")

    # Panel 3: shoreline position
    axs[2].plot(x_axis, xs_ts_m, marker="o")
    axs[2].set_ylabel("Shoreline x_s (m)")
    axs[2].grid(True, linestyle="--", alpha=0.5)
    axs[2].set_title("Shoreline position vs time")

    # X-axis label
    if USE_CALENDAR_YEARS:
        axs[2].set_xlabel("Calendar year")
    else:
        axs[2].set_xlabel("Year index")

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()
