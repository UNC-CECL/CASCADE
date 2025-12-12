"""
Profile Visualization Tool

Plots actual cross-shore elevation profiles to see what the dune extraction
algorithm is doing. Helps diagnose if dune heights are correct.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==========================
# CONFIG
# ==========================
NPZ_FILE = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\elevations\2009_pea_hatteras\domain_90.npy")

MHW_M = 0.26
BEACH_START_THR_M = 0.50
DUNE_WINDOW_PX = 8
BERM_SEARCH_START = 1
BERM_SEARCH_END = 4

# Which alongshore profiles to plot (e.g., [0, 10, 25, 40])
PROFILES_TO_PLOT = [0, 12, 25, 37, 49]


def load_and_prep_array(filepath):
    """Load array and convert to MHW-relative."""
    arr = np.load(filepath).astype(float)
    arr = arr - MHW_M
    arr[arr < -1.0] = -3.0
    return arr


def process_single_profile(prof, profile_idx):
    """
    Process one profile and return all the key points for visualization.

    Returns dict with:
    - profile: the elevation data
    - start_beach: index where beach starts
    - dune_loc: index of dune peak
    - dune_elev: elevation of dune peak
    - berm_elev: elevation of berm
    - dune_height: calculated height
    - berm_window_indices: indices used to find berm
    """
    # Flip so ocean is at index 0
    prof = np.flip(prof)

    result = {'profile': prof, 'profile_idx': profile_idx}

    # Find beach start
    idx = np.where(prof > BEACH_START_THR_M)[0]
    if idx.size == 0:
        return None
    start_beach = int(idx[0])
    result['start_beach'] = start_beach

    # Find dune peak
    end_dune_search = min(start_beach + DUNE_WINDOW_PX, prof.size)
    dune_search_window = prof[start_beach:end_dune_search]
    if dune_search_window.size == 0:
        return None

    dune_elev = float(np.max(dune_search_window))
    dune_loc = start_beach + np.argmax(dune_search_window)
    result['dune_elev'] = dune_elev
    result['dune_loc'] = dune_loc
    result['dune_search_window'] = (start_beach, end_dune_search)

    # Find berm
    berm_search_start = start_beach + BERM_SEARCH_START
    berm_search_end = min(start_beach + BERM_SEARCH_END, end_dune_search)

    if berm_search_end <= berm_search_start:
        berm_elev = float(prof[start_beach])
        berm_window_indices = [start_beach]
    else:
        berm_window = prof[berm_search_start:berm_search_end]
        berm_elev = float(np.min(berm_window))
        berm_window_indices = list(range(berm_search_start, berm_search_end))

    result['berm_elev'] = berm_elev
    result['berm_window'] = berm_window_indices

    # Calculate dune height
    dune_height = dune_elev - berm_elev
    result['dune_height'] = dune_height

    return result


def plot_profiles(arr, profiles_to_plot, domain_name):
    """Plot multiple profiles with annotations."""
    n_profiles = len(profiles_to_plot)
    fig, axes = plt.subplots(n_profiles, 1, figsize=(14, 4 * n_profiles))

    if n_profiles == 1:
        axes = [axes]

    for ax, prof_idx in zip(axes, profiles_to_plot):
        if prof_idx >= arr.shape[0]:
            print(f"Profile {prof_idx} out of range (max: {arr.shape[0] - 1})")
            continue

        prof_data = arr[prof_idx, :]
        result = process_single_profile(prof_data, prof_idx)

        if result is None:
            ax.text(0.5, 0.5, f"Profile {prof_idx}: No beach found",
                    ha='center', va='center', transform=ax.transAxes)
            continue

        prof = result['profile']
        x = np.arange(len(prof))

        # Plot the profile
        ax.plot(x, prof, 'k-', linewidth=2, label='Elevation profile')
        ax.axhline(0, color='blue', linestyle='--', alpha=0.5, label='MHW (0m)')
        ax.axhline(BEACH_START_THR_M, color='cyan', linestyle=':', alpha=0.5,
                   label=f'Beach start threshold ({BEACH_START_THR_M}m)')

        # Mark key points
        start_beach = result['start_beach']
        dune_loc = result['dune_loc']
        dune_elev = result['dune_elev']
        berm_elev = result['berm_elev']
        dune_height = result['dune_height']

        # Beach start
        ax.plot(start_beach, prof[start_beach], 'go', markersize=10,
                label=f'Beach start (x={start_beach})')

        # Dune search window
        ds_start, ds_end = result['dune_search_window']
        ax.axvspan(ds_start, ds_end, alpha=0.2, color='orange',
                   label='Dune search window')

        # Berm window
        berm_indices = result['berm_window']
        if len(berm_indices) > 0:
            ax.axvspan(berm_indices[0], berm_indices[-1], alpha=0.3, color='yellow',
                       label='Berm search window')

        # Dune peak
        ax.plot(dune_loc, dune_elev, 'r^', markersize=15,
                label=f'Dune peak ({dune_elev:.2f}m)')

        # Berm elevation line
        ax.axhline(berm_elev, color='green', linestyle='-.', linewidth=2,
                   label=f'Berm elevation ({berm_elev:.2f}m)')

        # Dune height annotation
        mid_x = (start_beach + dune_loc) / 2
        ax.annotate('', xy=(mid_x, dune_elev), xytext=(mid_x, berm_elev),
                    arrowprops=dict(arrowstyle='<->', color='red', lw=2))
        ax.text(mid_x + 2, (dune_elev + berm_elev) / 2,
                f'Dune height\n{dune_height:.2f}m',
                color='red', fontweight='bold', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_xlabel('Cross-shore position (pixels, ocean at left)')
        ax.set_ylabel('Elevation (m, MHW-relative)')
        ax.set_title(f'{domain_name} - Profile {prof_idx}')
        ax.legend(loc='upper right', fontsize=8)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, min(100, len(prof)))

        # Print stats
        print(f"\nProfile {prof_idx}:")
        print(f"  Beach start: x={start_beach}, z={prof[start_beach]:.2f}m")
        print(f"  Berm elevation: {berm_elev:.2f}m")
        print(f"  Dune peak: x={dune_loc}, z={dune_elev:.2f}m")
        print(f"  Dune HEIGHT: {dune_height:.2f}m")

    plt.tight_layout()
    plt.show()


def main():
    print(f"Loading: {NPZ_FILE}")
    arr = load_and_prep_array(NPZ_FILE)

    print(f"Array shape: {arr.shape} (alongshore x cross-shore)")
    print(f"Plotting profiles: {PROFILES_TO_PLOT}")

    domain_name = NPZ_FILE.stem
    plot_profiles(arr, PROFILES_TO_PLOT, domain_name)


if __name__ == "__main__":
    main()