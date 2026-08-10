#!/usr/bin/env python3
"""
CASCADE Overwash Flux (Qow) — Standalone Plotter
=================================================
Reads a saved CASCADE run NPZ file and produces overwash figures WITHOUT
re-running the model or importing CASCADE.

How it works
------------
cascade.save() writes a single NPZ containing a pickled Cascade object under
the key 'cascade'.  That object holds a list of Barrier3D objects at
cascade._barrier3d[i], each of which stores its Qow time series at
b3d._QowTS (shape = [nt]).  This script uses a flexible unpickler to read
that data without needing the cascade or barrier3d packages installed.

Figures produced
----------------
  1. Heatmap      — domain (x) × year (y), colour = Qow magnitude
  2. Lines by year — Qow profile along the island, one line per year
  3. Time series  — island-mean and island-total Qow over time

All figures are saved into the same folder as the NPZ file.

Usage
-----
1. Point NPZ_PATH to your cascade run file.
2. Set START_YEAR, END_YEAR, and domain parameters to match the run.
3. python plot_overwash.py
"""

import os
import io
import pickle
import types
import zipfile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable

# =============================================================================
# SECTION 1: RUN CONFIGURATION  <- edit this every time
# =============================================================================

# Full path to the single .npz file produced by cascade.save()
NPZ_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1984_2004_basestorms_Hs2p0\HAT_1984_2004_basestorms_Hs2p0.npz"

# Simulation period -- used for year-axis labels only
START_YEAR = 1984
END_YEAR   = 2004

# Domain configuration -- must match the run that produced this NPZ
NUM_REAL_DOMAINS    = 90
NUM_BUFFER_DOMAINS  = 15
FIRST_GIS_DOMAIN_ID = 1   # GIS ID of the southernmost real domain

# =============================================================================
# SECTION 2: PLOT & OUTPUT CONFIGURATION
# =============================================================================

# Which figures to produce
PLOT_HEATMAP    = True   # domain x year heatmap
PLOT_LINES      = True   # per-year Qow profile lines
PLOT_TIMESERIES = True   # island-mean and total Qow over time

# Save figures alongside the NPZ (True) or just display (False)
SAVE_FIGURES = True

# Export the Qow matrix as a CSV
SAVE_CSV = True

# Lines plot: only draw every Nth year to reduce clutter (1 = every year)
PLOT_EVERY_N_YEARS = 1

# X-axis tick spacing (GIS domain IDs)
DOMAIN_TICK_STEP = 5

# Heatmap colormap -- any sequential matplotlib colormap name
HEATMAP_CMAP = "YlOrRd"

# =============================================================================
# SECTION 3: DERIVED DOMAIN CONSTANTS  (do not edit)
# =============================================================================

TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS
START_REAL_INDEX   = NUM_BUFFER_DOMAINS
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS
LAST_GIS_DOMAIN_ID = FIRST_GIS_DOMAIN_ID + NUM_REAL_DOMAINS - 1
RUN_YEARS          = END_YEAR - START_YEAR

# =============================================================================
# NPZ LOADER
# =============================================================================

# Class registry for the flexible unpickler.
# Maps (module, name) to a dummy class so the pickled Cascade/Barrier3D
# objects can be loaded without those packages being installed.
_class_registry: dict = {}


def _make_dummy_class(module: str, name: str) -> type:
    """Return a cached dummy class that captures pickle state as attributes."""
    key = (module, name)
    if key not in _class_registry:
        def _setstate(self, state):
            if isinstance(state, dict):
                self.__dict__.update(state)
            else:
                self._unpickled_state = state
        cls = type(
            f"{module}.{name}",
            (),
            {"__init__":    lambda self, *a, **k: None,
             "__setstate__": _setstate},
        )
        _class_registry[key] = cls
    return _class_registry[key]


class _FlexUnpickler(pickle.Unpickler):
    """Unpickler that substitutes dummy classes for any missing module."""
    def find_class(self, module, name):
        try:
            return super().find_class(module, name)
        except (ModuleNotFoundError, AttributeError):
            return _make_dummy_class(module, name)


def load_cascade_npz(npz_path: str):
    """
    Load a CASCADE .npz save file and return the unpickled Cascade object.

    The NPZ contains a single member 'cascade.npy', which is a numpy object
    array (shape [1]) wrapping the pickled Cascade instance.  The .npy v1
    header is 128 bytes; the pickle payload follows immediately after.

    Returns
    -------
    cascade_obj : dummy Cascade instance with all attributes accessible
    """
    if not os.path.isfile(npz_path):
        raise FileNotFoundError(f"NPZ file not found:\n  {npz_path}")

    size_mb = os.path.getsize(npz_path) / 1e6
    print(f"Loading: {os.path.basename(npz_path)}  ({size_mb:.1f} MB) ...")

    with zipfile.ZipFile(npz_path) as zf:
        members    = zf.namelist()
        npy_member = next((m for m in members if m.endswith(".npy")), None)
        if npy_member is None:
            raise ValueError(
                f"No .npy member found in {npz_path}.\n"
                f"Members present: {members}"
            )
        raw = zf.read(npy_member)

    # Skip the 128-byte NPY v1.0 header to reach the pickle payload
    NPY_HEADER_LEN = 128
    buf        = io.BytesIO(raw[NPY_HEADER_LEN:])
    wrapper    = _FlexUnpickler(buf).load()   # numpy object array, shape (1,)

    cascade_obj = wrapper.flat[0] if isinstance(wrapper, np.ndarray) else wrapper

    if not hasattr(cascade_obj, "_barrier3d"):
        raise AttributeError(
            "Loaded object has no '_barrier3d' attribute.\n"
            "This may not be a CASCADE save file, or the internal format has changed."
        )

    print(f"  done.")
    return cascade_obj


# =============================================================================
# QOW EXTRACTION
# =============================================================================

# Candidate attribute names for the Qow time series on a Barrier3D object.
_QOW_ATTR_CANDIDATES = ["_QowTS", "QowTS", "_Qow_TS", "Qow_TS", "_Qow", "Qow"]


def _extract_qow_ts(b3d, domain_idx: int):
    """
    Extract the 1-D Qow time series from one Barrier3D object.

    Returns
    -------
    arr  : 1-D float array, shape (nt,)
    attr : str -- the attribute name that was used
    """
    for attr in _QOW_ATTR_CANDIDATES:
        if hasattr(b3d, attr):
            arr = np.asarray(getattr(b3d, attr), dtype=float).squeeze()
            if arr.ndim == 1 and arr.size > 1:
                return arr, attr

    ow_attrs = [k for k in vars(b3d) if "qow" in k.lower() or "ow" in k.lower()]
    raise AttributeError(
        f"Domain {domain_idx}: could not find a Qow time series.\n"
        f"  Tried: {_QOW_ATTR_CANDIDATES}\n"
        f"  OW-related attrs found: {ow_attrs}"
    )


def build_qow_matrix(cascade_obj):
    """
    Build the [time x real_domain] Qow matrix from the Cascade object.

    Returns
    -------
    years          : int array (nt,)
    Q              : float array (nt, ndom) -- Qow values
    gis_domain_ids : int array (ndom,)      -- GIS IDs of real domains
    qow_attr       : str                    -- attribute name used
    """
    b3d_list = cascade_obj._barrier3d
    n_loaded = len(b3d_list)

    print(f"Barrier3D objects in file: {n_loaded}  (TOTAL_DOMAINS={TOTAL_DOMAINS})")

    # Determine which indices are real domains
    if n_loaded == TOTAL_DOMAINS:
        real_indices = list(range(START_REAL_INDEX, END_REAL_INDEX))
    elif n_loaded == NUM_REAL_DOMAINS:
        print("  Note: file contains real domains only (no buffers).")
        real_indices = list(range(NUM_REAL_DOMAINS))
    else:
        print(
            f"  Warning: expected {TOTAL_DOMAINS} (with buffers) or "
            f"{NUM_REAL_DOMAINS} (without), found {n_loaded}.\n"
            "  Check NUM_REAL_DOMAINS / NUM_BUFFER_DOMAINS."
        )
        real_indices = list(range(min(n_loaded, NUM_REAL_DOMAINS)))

    qow_attr = None
    rows     = []

    for padded_idx in real_indices:
        ts, attr = _extract_qow_ts(b3d_list[padded_idx], padded_idx)
        if qow_attr is None:
            qow_attr = attr
        rows.append(ts)

    lengths = {r.size for r in rows}
    if len(lengths) > 1:
        raise ValueError(
            f"Qow time series have inconsistent lengths: {lengths}. "
            "The run may have stopped early in some domains."
        )

    Q              = np.column_stack(rows)
    nt             = Q.shape[0]
    years          = np.arange(START_YEAR, START_YEAR + nt, dtype=int)
    gis_domain_ids = np.arange(
        FIRST_GIS_DOMAIN_ID,
        FIRST_GIS_DOMAIN_ID + len(real_indices),
        dtype=int,
    )

    print(f"Qow matrix: {nt} years x {Q.shape[1]} domains  |  attr='{qow_attr}'")
    print(f"Qow range : {np.nanmin(Q):.4f} .. {np.nanmax(Q):.4f} dam3/yr")
    return years, Q, gis_domain_ids, qow_attr


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def _run_label() -> str:
    return os.path.splitext(os.path.basename(NPZ_PATH))[0]


def _out_dir() -> str:
    return os.path.dirname(NPZ_PATH)


def _save_or_show(fig, suffix: str):
    if SAVE_FIGURES:
        out = os.path.join(_out_dir(), f"{_run_label()}_{suffix}.png")
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved: {out}")
    plt.show()


# --------------------------------------------------------------------------
# Plot 1 -- Heatmap
# --------------------------------------------------------------------------

def plot_heatmap(years, Q, gis_domain_ids):
    """
    2-D heatmap: x = GIS domain ID, y = year, colour = Qow.
    Reveals persistent spatial hotspots and high-overwash years at a glance.
    """
    nt    = len(years)
    fig_h = max(5.0, nt * 0.30)
    fig, ax = plt.subplots(figsize=(14, fig_h), constrained_layout=True)
    fig.patch.set_facecolor("white")

    extent = [
        gis_domain_ids[0] - 0.5, gis_domain_ids[-1] + 0.5,
        years[-1] + 0.5,          years[0]  - 0.5,
    ]
    im = ax.imshow(
        Q, aspect="auto", cmap=HEATMAP_CMAP,
        extent=extent, origin="upper", interpolation="nearest",
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("Overwash Flux  Qow  (dam\u00b3/yr)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    xtick_vals = np.arange(gis_domain_ids[0], gis_domain_ids[-1] + 1, DOMAIN_TICK_STEP)
    ax.set_xticks(xtick_vals)
    ax.set_xticklabels(xtick_vals, rotation=45, ha="right", fontsize=9)

    y_step = 1 if nt <= 25 else 2
    ax.set_yticks(years[::y_step])
    ax.set_yticklabels(years[::y_step], fontsize=9)

    ax.set_xlabel("GIS Domain ID", fontsize=11, fontweight="bold", labelpad=6)
    ax.set_ylabel("Year", fontsize=11, fontweight="bold", labelpad=6)
    ax.set_title(
        f"Overwash Flux (Qow)  |  {_run_label()}\n"
        f"{START_YEAR}\u2013{END_YEAR}  |  {NUM_REAL_DOMAINS} real domains",
        fontsize=12, fontweight="bold", color="#1a2a3a", pad=8,
    )
    _save_or_show(fig, "Qow_heatmap")


# --------------------------------------------------------------------------
# Plot 2 -- Lines by year
# --------------------------------------------------------------------------

def plot_lines_by_year(years, Q, gis_domain_ids):
    """
    One line per model year, coloured by year (viridis).
    Shows how the spatial pattern of overwash shifts over time.
    """
    year_indices = list(range(0, len(years), PLOT_EVERY_N_YEARS))
    n_lines      = len(year_indices)
    cmap_lines   = plt.get_cmap("viridis", n_lines)
    norm_lines   = mcolors.Normalize(vmin=years[0], vmax=years[-1])

    fig, ax = plt.subplots(figsize=(13, 5.5), constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    for k, i in enumerate(year_indices):
        col = cmap_lines(k / max(n_lines - 1, 1))
        ax.plot(gis_domain_ids, Q[i, :], lw=1.4, alpha=0.85, color=col)

    sm = ScalarMappable(cmap="viridis", norm=norm_lines)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01)
    cbar.set_label("Year", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    ax.axhline(0, color="#444444", lw=0.8, ls="--", alpha=0.5)

    xtick_vals = np.arange(gis_domain_ids[0], gis_domain_ids[-1] + 1, DOMAIN_TICK_STEP)
    ax.set_xticks(xtick_vals)
    ax.set_xticklabels(xtick_vals, rotation=45, ha="right", fontsize=9)
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis="both", which="major", direction="in", length=5, labelsize=10)
    ax.tick_params(axis="both", which="minor", direction="in", length=3)
    ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_linewidth(1.1)

    ax.set_xlabel("GIS Domain ID", fontsize=11, fontweight="bold", labelpad=6)
    ax.set_ylabel("Overwash Flux  Qow  (dam\u00b3/yr)", fontsize=11,
                  fontweight="bold", labelpad=6)
    every_str = f"  |  every {PLOT_EVERY_N_YEARS} yr(s) shown" if PLOT_EVERY_N_YEARS > 1 else ""
    ax.set_title(
        f"Overwash Flux Profile by Year  |  {_run_label()}\n"
        f"{START_YEAR}\u2013{END_YEAR}{every_str}",
        fontsize=12, fontweight="bold", color="#1a2a3a", pad=8,
    )
    _save_or_show(fig, "Qow_lines")


# --------------------------------------------------------------------------
# Plot 3 -- Time series
# --------------------------------------------------------------------------

def plot_timeseries(years, Q):
    """
    Top panel: island-mean Qow per year.
    Bottom panel: island-total (sum) Qow per year.
    """
    q_mean  = np.nanmean(Q, axis=1)
    q_total = np.nansum(Q,  axis=1)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7),
                                    sharex=True, constrained_layout=True)
    fig.patch.set_facecolor("white")

    panels = [
        (ax1, q_mean,  "#1F4E79", "Mean Qow  (dam\u00b3/yr)",  "Mean"),
        (ax2, q_total, "#C0392B", "Total Qow  (dam\u00b3/yr)", "Total"),
    ]
    for ax, vals, col, ylabel, prefix in panels:
        ax.set_facecolor("white")
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_linewidth(1.1)
        ax.tick_params(axis="both", which="major", direction="in",
                       length=5, labelsize=10)
        ax.grid(True, which="major", ls=":", lw=0.6, alpha=0.4, color="gray")

        ax.fill_between(years, vals, alpha=0.18, color=col)
        ax.plot(years, vals, lw=2.2, color=col, marker="o", markersize=4,
                label=f"{prefix} Qow")
        pmean = np.nanmean(vals)
        ax.axhline(pmean, ls="--", lw=1.1, color=col, alpha=0.55,
                   label=f"Period mean = {pmean:.3f}")
        ax.set_ylabel(ylabel, fontsize=11, fontweight="bold", labelpad=6)
        ax.legend(fontsize=9, framealpha=0.9)

    y_step = 1 if RUN_YEARS <= 20 else 2
    ax2.set_xticks(years[::y_step])
    ax2.set_xticklabels(years[::y_step], rotation=45, ha="right", fontsize=9)
    ax2.set_xlabel("Year", fontsize=11, fontweight="bold", labelpad=6)

    fig.suptitle(
        f"Island-wide Overwash Flux  |  {_run_label()}\n{START_YEAR}\u2013{END_YEAR}",
        fontsize=12, fontweight="bold", color="#1a2a3a",
    )
    _save_or_show(fig, "Qow_timeseries")


# --------------------------------------------------------------------------
# CSV export
# --------------------------------------------------------------------------

def export_csv(years, Q, gis_domain_ids):
    """Write the Qow matrix as CSV: rows = years, columns = GIS domain IDs."""
    out    = os.path.join(_out_dir(), f"{_run_label()}_Qow_by_year_domain.csv")
    header = ["Year"] + [f"GIS_{d}" for d in gis_domain_ids]
    with open(out, "w") as f:
        f.write(",".join(header) + "\n")
        for i, yr in enumerate(years):
            row = [str(int(yr))] + [f"{v:.6f}" for v in Q[i, :]]
            f.write(",".join(row) + "\n")
    print(f"  Saved CSV: {out}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("CASCADE Overwash Plotter -- Standalone")
    print("=" * 70)
    print(f"NPZ    : {NPZ_PATH}")
    print(f"Period : {START_YEAR}-{END_YEAR}  ({RUN_YEARS} years)")
    print(f"Domains: {NUM_REAL_DOMAINS} real  |  {NUM_BUFFER_DOMAINS} buffers each side  "
          f"|  GIS {FIRST_GIS_DOMAIN_ID}-{LAST_GIS_DOMAIN_ID}")
    print()

    cascade_obj = load_cascade_npz(NPZ_PATH)
    print()

    years, Q, gis_domain_ids, qow_attr = build_qow_matrix(cascade_obj)
    print()

    # Trim if model ran longer than END_YEAR - START_YEAR
    if len(years) > RUN_YEARS:
        years = years[:RUN_YEARS]
        Q     = Q[:RUN_YEARS, :]
        print(f"Note: trimmed to {RUN_YEARS} years ({START_YEAR}-{END_YEAR}).\n")

    print("Generating figures...")
    if PLOT_HEATMAP:
        plot_heatmap(years, Q, gis_domain_ids)
    if PLOT_LINES:
        plot_lines_by_year(years, Q, gis_domain_ids)
    if PLOT_TIMESERIES:
        plot_timeseries(years, Q)
    if SAVE_CSV:
        export_csv(years, Q, gis_domain_ids)

    print("\nDone.")


if __name__ == "__main__":
    main()
