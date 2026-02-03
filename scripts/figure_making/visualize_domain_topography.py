"""
Visualize CASCADE Domain Topography
====================================
Creates a cross-shore vs alongshore elevation heatmap for a single CASCADE domain,
with labeled dune domain and interior domain regions.

Usage:
    python visualize_domain_topography.py

Author: Hannah (UNC Chapel Hill)
Date: January 2025
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# CONFIGURATION
# ============================================================================

# Domain to visualize
DOMAIN_NUMBER = 45  # Change this to visualize different domains

# Time step to visualize
TIME_STEP = -1  # -1 for final time step, or specify a year (0, 1, 2, etc.)

# File paths
CASCADE_OUTPUT_FILE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_natural\HAT_1978_1997_natural.npz"

# Domain parameters (adjust to match your CASCADE setup)
DY = 10.0  # Cross-shore resolution (m) - CASCADE cells are 10m × 10m
DX = 10.0  # Alongshore resolution (m) - CASCADE cells are 10m × 10m
DUNE_WIDTH = 2  # Number of cross-shore cells in dune domain (first 2 rows)

# Cross-shore extent options
CROP_CROSSSHORE = True   # If True, crop to MAX_CROSSSHORE_M; if False, use CASCADE's actual width
MAX_CROSSSHORE_M = 500    # Maximum cross-shore distance to display when CROP_CROSSSHORE = True (m)
                          # Common options: 300 (nearshore focus), 500 (moderate), 1000 (wide view)
                          # Set CROP_CROSSSHORE = False to use CASCADE's evolved island width

# Figure parameters
FIGURE_WIDTH = 10
FIGURE_HEIGHT = 6
DPI = 150
AUTO_SCALE_Z = True   # If True, automatically set Z_LIM based on data range
Z_LIM = 8.0           # Maximum elevation for colorbar (m MSL) - used only if AUTO_SCALE_Z = False
                      # Your data ranges from -1 to ~7.3 m
SAVE_FIGURE = True
OUTPUT_FILE = f"domain_{DOMAIN_NUMBER}_topography.png"

# ============================================================================
# FUNCTIONS
# ============================================================================

def load_elevation_data(filepath, domain_idx, time_step=-1, max_crossshore_m=None, dy=10.0):
    """
    Load elevation data for a specific domain from CASCADE output.
    
    Parameters:
    -----------
    filepath : str
        Path to NPZ file containing CASCADE output
    domain_idx : int
        Domain index to extract (0-indexed)
    time_step : int
        Time step to extract. Default -1 gets the final time step.
    max_crossshore_m : float or None
        Maximum cross-shore distance to include (m). If None, includes full domain.
    dy : float
        Cross-shore resolution (m) for calculating crop index
        
    Returns:
    --------
    elevation : ndarray
        2D array of elevation values (cross-shore x alongshore) in meters MSL
    """
    # Load CASCADE output
    data = np.load(filepath, allow_pickle=True)
    
    # Extract CASCADE object
    cascade = data["cascade"][0]
    barrier3d = cascade.barrier3d
    
    # Get the domain
    domain = barrier3d[domain_idx]
    
    # Get domain topography at specified time step
    # DomainTS is in units of dam (decameters), so multiply by 10 to get meters
    interior_elevation = domain.DomainTS[time_step] * 10
    
    # Get dune domain at specified time step
    # DuneDomain is also in dam, convert to meters and add berm elevation
    dune_elevation = (domain.DuneDomain[time_step, :, :] + domain.BermEl) * 10
    
    # Rotate and flip dune domain to match orientation
    dune_elevation = np.rot90(dune_elevation)
    dune_elevation = np.flipud(dune_elevation)
    
    # Combine dune and interior domains
    # Stack dune domain (first few rows) with interior domain
    elevation = np.vstack([dune_elevation, interior_elevation])
    
    # Crop cross-shore extent if requested
    if max_crossshore_m is not None:
        max_crossshore_cells = int(max_crossshore_m / dy)
        if max_crossshore_cells < elevation.shape[0]:
            elevation = elevation[:max_crossshore_cells, :]
            print(f"Cropped to {max_crossshore_m} m cross-shore ({max_crossshore_cells} cells)")
    
    # Set any negative values to -1 (below sea level indicator)
    elevation[elevation < 0] = -1
    
    return elevation


def create_topography_plot(elevation, dy, dx, dune_width, domain_num, z_lim):
    """
    Create a heatmap visualization of domain topography.
    
    Parameters:
    -----------
    elevation : ndarray
        2D array of elevation values (cross-shore x alongshore)
    dy : float
        Cross-shore resolution (m)
    dx : float
        Alongshore resolution (m)
    dune_width : int
        Number of cross-shore cells in dune domain
    domain_num : int
        Domain number for title
    z_lim : float
        Maximum elevation for colorbar
        
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    # Create coordinate arrays
    num_crossshore, num_alongshore = elevation.shape
    crossshore_width = np.arange(num_crossshore) * dy
    alongshore_length = np.arange(num_alongshore) * dx
    
    # Calculate actual physical dimensions
    total_crossshore_m = num_crossshore * dy
    total_alongshore_m = num_alongshore * dx
    
    # Calculate figure dimensions to maintain proper aspect ratio
    # aspect_ratio = physical_width / physical_height
    aspect_ratio = total_alongshore_m / total_crossshore_m
    
    # Set figure size based on desired height and calculated width
    fig_height = FIGURE_HEIGHT
    fig_width = fig_height * aspect_ratio * 1.2  # 1.2 factor adds space for colorbar
    
    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Plot elevation heatmap using terrain colormap (matches GIF script)
    im = ax.pcolormesh(alongshore_length, crossshore_width, elevation,
                       cmap='terrain', vmin=-1.1, vmax=z_lim, shading='auto')
    
    # Set equal aspect ratio so cells display with true proportions
    ax.set_aspect('equal', adjustable='box')
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, label='Elevation (m MSL)')
    
    # Mark dune domain boundary
    dune_boundary = dune_width * dy
    ax.axhline(y=dune_boundary, color='red', linewidth=2, 
               linestyle='-', alpha=0.8)
    
    # Add region labels
    # Dune domain label (near bottom)
    dune_center_y = dune_boundary / 2
    dune_center_x = num_alongshore * dx / 2
    ax.text(dune_center_x, dune_center_y, 'Dune Domain',
            ha='center', va='center', fontsize=12, fontweight='bold',
            color='white', bbox=dict(boxstyle='round', facecolor='red', 
                                    alpha=0.7, edgecolor='none'))
    
    # Interior domain label (middle of interior)
    interior_center_y = dune_boundary + (num_crossshore * dy - dune_boundary) / 2
    ax.text(dune_center_x, interior_center_y, 'Interior Domain',
            ha='center', va='center', fontsize=12, fontweight='bold',
            color='white', bbox=dict(boxstyle='round', facecolor='black', 
                                    alpha=0.5, edgecolor='none'))
    
    # Labels and title
    ax.set_xlabel('Alongshore Length (m)', fontsize=12)
    ax.set_ylabel('Cross-shore Width (m)', fontsize=12)
    ax.set_title(f'CASCADE Domain {domain_num} Topography', 
                 fontsize=14, fontweight='bold')
    
    # Adjust layout
    plt.tight_layout()
    
    return fig, ax


def print_domain_statistics(elevation, dy, dx, dune_width):
    """
    Print useful statistics about the domain topography.
    
    Parameters:
    -----------
    elevation : ndarray
        2D array of elevation values
    dy : float
        Cross-shore resolution (m)
    dx : float
        Alongshore resolution (m)
    dune_width : int
        Number of cross-shore cells in dune domain
    """
    num_crossshore, num_alongshore = elevation.shape
    
    print("\n" + "="*60)
    print("DOMAIN TOPOGRAPHY STATISTICS")
    print("="*60)
    print(f"Domain dimensions: {num_crossshore} x {num_alongshore} cells")
    print(f"Cross-shore extent: {num_crossshore * dy:.1f} m")
    print(f"Alongshore extent: {num_alongshore * dx:.1f} m")
    print(f"Dune domain: {dune_width} cells ({dune_width * dy:.1f} m)")
    print(f"Interior domain: {num_crossshore - dune_width} cells "
          f"({(num_crossshore - dune_width) * dy:.1f} m)")
    print()
    
    # Elevation statistics
    print("Elevation statistics (entire domain):")
    print(f"  Min elevation: {np.nanmin(elevation):.2f} m MSL")
    print(f"  Max elevation: {np.nanmax(elevation):.2f} m MSL")
    print(f"  Mean elevation: {np.nanmean(elevation):.2f} m MSL")
    print(f"  Std elevation: {np.nanstd(elevation):.2f} m")
    print()
    
    # Dune domain statistics
    dune_elevation = elevation[:dune_width, :]
    print("Dune domain statistics:")
    print(f"  Mean elevation: {np.nanmean(dune_elevation):.2f} m MSL")
    print(f"  Max elevation: {np.nanmax(dune_elevation):.2f} m MSL")
    print()
    
    # Interior domain statistics
    interior_elevation = elevation[dune_width:, :]
    print("Interior domain statistics:")
    print(f"  Mean elevation: {np.nanmean(interior_elevation):.2f} m MSL")
    print(f"  Max elevation: {np.nanmax(interior_elevation):.2f} m MSL")
    print("="*60 + "\n")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print(f"\nVisualizing CASCADE Domain {DOMAIN_NUMBER} Topography")
    print(f"Time step: {TIME_STEP if TIME_STEP >= 0 else 'Final'}")
    if CROP_CROSSSHORE:
        print(f"Cross-shore extent: Cropped to {MAX_CROSSSHORE_M} m")
    else:
        print(f"Cross-shore extent: CASCADE's actual evolved width")
    print(f"Loading data from: {CASCADE_OUTPUT_FILE}\n")
    
    try:
        # Determine crop parameters
        max_crossshore = MAX_CROSSSHORE_M if CROP_CROSSSHORE else None
        
        # Load elevation data
        elevation = load_elevation_data(CASCADE_OUTPUT_FILE, DOMAIN_NUMBER, TIME_STEP, 
                                       max_crossshore_m=max_crossshore, dy=DY)
        
        # Print statistics
        print_domain_statistics(elevation, DY, DX, DUNE_WIDTH)
        
        # Determine Z_LIM based on auto-scaling setting
        if AUTO_SCALE_Z:
            max_elev = np.nanmax(elevation)
            z_lim_to_use = np.ceil(max_elev)  # Round up to next integer
            print(f"AUTO_SCALE_Z enabled: Setting Z_LIM to {z_lim_to_use} m\n")
        else:
            z_lim_to_use = Z_LIM
            # Check if Z_LIM is appropriate for the data
            max_elev = np.nanmax(elevation)
            if max_elev > Z_LIM:
                print(f"WARNING: Maximum elevation ({max_elev:.2f} m) exceeds Z_LIM ({Z_LIM} m)")
                print(f"         Consider increasing Z_LIM to {np.ceil(max_elev)} or higher")
                print(f"         for accurate color representation.\n")
        
        # Create visualization
        fig, ax = create_topography_plot(elevation, DY, DX, 
                                         DUNE_WIDTH, DOMAIN_NUMBER, z_lim_to_use)
        
        # Save figure if requested
        if SAVE_FIGURE:
            plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight')
            print(f"Figure saved to: {OUTPUT_FILE}")
        
        # Display figure
        plt.show()
        
        print("\nVisualization complete!")
        
    except FileNotFoundError:
        print(f"ERROR: Could not find file: {CASCADE_OUTPUT_FILE}")
        print("Please update the CASCADE_OUTPUT_FILE path in the configuration section.")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
