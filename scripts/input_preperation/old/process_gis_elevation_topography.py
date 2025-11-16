# ==============================================================================
# Hatteras CASCADE Dune and Topography Extractor (MHW-relative, fixed)
# ------------------------------------------------------------------------------
# Assumptions per your check:
#   - rows = alongshore
#   - columns = cross-shore
#   - ocean is on the RIGHT (land on the left)
# Outputs (decameters):
#   - topography: (N_ROWS x alongshore)
#   - dune:       (alongshore,)
# ==============================================================================

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- Pick a single test domain to figure out numpy orientation ---
f = Path(r"/data/hatteras_init/elevations/old_study_area_2009/domain_1_resampled.npy")
a = np.load(f)

# --- Plot it raw ---
plt.figure(figsize=(8,6))
plt.imshow(a, cmap="terrain", origin="upper", aspect="auto")
plt.colorbar(label="Elevation (m NAVD88)")
plt.title(f.name)
plt.show()

# --- Edge diagnostics ---
edge = 5  # width of edge window in pixels
left_low   = np.nanpercentile(a[:, :edge], 25)
right_low  = np.nanpercentile(a[:, -edge:], 25)
top_low    = np.nanpercentile(a[:edge, :], 25)
bottom_low = np.nanpercentile(a[-edge:, :], 25)

print("25th percentile near edges (lower ~ ocean side):")
print(f"  Left   = {left_low:.2f}")
print(f"  Right  = {right_low:.2f}")
print(f"  Top    = {top_low:.2f}")
print(f"  Bottom = {bottom_low:.2f}")

# === INPUT/OUTPUT PATHS =====================================================

LOAD_PATH = r'/data/hatteras_init/elevations/old_study_area_2009'
BASE_OUT  = r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init'
TOPO_SAVE_PATH = os.path.join(BASE_OUT, 'topography')
DUNE_SAVE_PATH = os.path.join(BASE_OUT, 'dunes')
os.makedirs(TOPO_SAVE_PATH, exist_ok=True)
os.makedirs(DUNE_SAVE_PATH, exist_ok=True)

# === PARAMETERS =============================================================

# Hydro / units
MHW = 0.26               # meters NAVD88 (use your chosen value)
BERM_ELEV_NAVD = 1.9     # meters NAVD88 (CASCADE NC standard)
BERM_ELEV_REL  = BERM_ELEV_NAVD - MHW   # berm relative to MHW (meters); e.g., 1.64 if MHW=0.26

DX_M = 10.0              # cross-shore pixel size [m] (set correctly for your grids)
NODATA = -9999           # sentinel in input arrays (if present)

# Orientation enforcement (optional safety)
# If any arrays export with ocean on the LEFT, set this True to flip them.
FORCE_OCEAN_LEFT = False

# Dune finding / beach logic
BEACH_THRESH_MHW = 0.0   # meters above MHW considered "dry beach" (0.0–0.5 common)
DUNE_WINDOW_M = 80.0     # seaward search width from first dry-beach sample (meters)
DUNE_MIN_DAM = 0.1       # minimum dune height in decameters for CASCADE

# Output geometry
N_ROWS = 100             # samples landward of crest (dune -> back-barrier), per profile

# === HELPERS ================================================================

def natural_key(name: str) -> int:
    m = re.search(r'domain_(\d+)', name)
    return int(m.group(1)) if m else 10**9

def to_float(a: np.ndarray) -> np.ndarray:
    a = a.astype(float, copy=False)
    if NODATA is not None:
        a[a == NODATA] = np.nan
    return a

def ensure_ocean_right(b: np.ndarray) -> np.ndarray:
    """
    Ensure columns = cross-shore with ocean on the RIGHT.
    Uses low-percentile near left/right edges to infer ocean side.
    """
    if FORCE_OCEAN_LEFT:
        return np.fliplr(b)

    edge = min(5, max(1, b.shape[1] // 20))  # small window near edges
    left_low  = np.nanpercentile(b[:, :edge],  25)
    right_low = np.nanpercentile(b[:, -edge:], 25)
    ocean_is_left = left_low < right_low
    return np.fliplr(b) if ocean_is_left else b

def resample_to_length(vals: np.ndarray, out_len: int) -> np.ndarray:
    """Linearly resample a 1D array to exactly out_len (ignoring trailing NaNs)."""
    if vals.size == 0:
        return np.full(out_len, np.nan)
    finite_mask = np.isfinite(vals)
    if not finite_mask.any():
        return np.full(out_len, np.nan)
    last = np.where(finite_mask)[0][-1] + 1
    vals = vals[:last]
    x = np.linspace(0, 1, num=vals.size)
    xi = np.linspace(0, 1, num=out_len)
    return np.interp(xi, x, vals)

# === CORE LOGIC =============================================================

def process_topo_arrays(input_path, topo_out_path, dune_out_path,
                        mhw=MHW, berm_rel=BERM_ELEV_REL):
    files = sorted(
        [f for f in os.listdir(input_path) if f.endswith('.npy')],
        key=natural_key
    )

    dune_win_cols = max(1, int(round(DUNE_WINDOW_M / DX_M)))

    for fname in files:
        fpath = os.path.join(input_path, fname)
        base  = Path(fname).stem  # e.g., domain_1_resampled

        # Load & prepare
        a = np.load(fpath)
        a = to_float(a)
        a = ensure_ocean_right(a)  # land-left, ocean-right
        a = a - mhw                # now elevations are MHW-relative (0 ≈ shoreline)

        n_rows_along, n_cols_cross = a.shape  # rows=alongshore, cols=cross-shore

        # Outputs: (distance samples) x (alongshore rows)  — meters (convert later)
        topo_m = np.full((N_ROWS, n_rows_along), np.nan)
        dune_m = np.full((n_rows_along,), np.nan)

        for i in range(n_rows_along):
            # Cross-shore profile for alongshore row i (land -> ocean)
            x = a[i, :]

            # 1) First dry-beach index (>= threshold relative to MHW)
            idxs = np.where(np.isfinite(x) & (x >= BEACH_THRESH_MHW))[0]
            if idxs.size == 0:
                continue
            beach_idx = int(idxs[0])

            # 2) Dune crest in fixed window seaward of beach
            win_end = min(x.size, beach_idx + dune_win_cols)
            window = x[beach_idx:win_end]
            if window.size == 0 or not np.isfinite(window).any():
                continue
            crest_off = int(np.nanargmax(window))
            crest_idx = beach_idx + crest_off
            crest_elev = float(x[crest_idx])   # MHW-relative meters

            # 3) Landward segment = from land edge up to crest (exclusive),
            #    reversed so it runs dune -> back-barrier
            landward = x[:crest_idx][::-1]

            # 4) Resample to fixed N_ROWS and store (meters, MHW-relative)
            topo_m[:, i] = resample_to_length(landward, N_ROWS)

            # 5) Dune height above berm (both MHW-relative)
            dune_m[i] = max(crest_elev - berm_rel, 0.0)

        # Convert to decameters for CASCADE + min threshold in dam
        topo_dm = topo_m * 0.1
        dune_dm = dune_m * 0.1
        dune_dm = np.where(np.isfinite(dune_dm), np.maximum(dune_dm, DUNE_MIN_DAM), np.nan)

        # Save
        np.save(os.path.join(topo_out_path, f"{base}_topography.npy"), topo_dm)
        np.save(os.path.join(dune_out_path,  f"{base}_dune.npy"),        dune_dm)

        # Console summary
        topo_valid = np.isfinite(topo_dm).mean()
        dune_valid = np.isfinite(dune_dm).mean()
        print(f"{base}: topo valid={topo_valid:6.2%}, dune valid={dune_valid:6.2%}, "
              f"shape={topo_dm.shape[0]}x{topo_dm.shape[1]} (rows x alongshore)")

# === EXECUTE ================================================================

if __name__ == "__main__":
    process_topo_arrays(LOAD_PATH, TOPO_SAVE_PATH, DUNE_SAVE_PATH)
