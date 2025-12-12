"""
Dune Growth Diagnostic

Investigates why dunes are so tall and how they evolved over time.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ==========================
# CONFIG
# ==========================
NPZ_FILE = Path(
    r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs"
    r"\HAT_1978_1997_Natural_State\HAT_1978_1997_Natural_State.npz"
)

START_YEAR = 1978
DAM_TO_M = 10.0
DOMAIN_START = 15
DOMAIN_END = 104


# ==========================
# FUNCTIONS
# ==========================

def load_cascade(npz_path):
    """Load CASCADE object from .npz file."""
    data = np.load(npz_path, allow_pickle=True)
    return data["cascade"][0]


def get_dune_domain_timeseries(segment):
    """Get full DuneDomain time series."""
    if hasattr(segment, "DuneDomain"):
        return np.array(segment.DuneDomain, dtype=float)
    elif hasattr(segment, "_DuneDomain"):
        return np.array(segment._DuneDomain, dtype=float)
    else:
        raise AttributeError("No DuneDomain found")


def analyze_dune_evolution(cascade, domain_idx):
    """Analyze how dunes evolved over time for a specific domain."""
    seg = cascade.barrier3d[domain_idx]

    # Get dune domain time series (time, cross_shore, 2)
    dune_ts = get_dune_domain_timeseries(seg)
    n_t, n_x, _ = dune_ts.shape

    # Extract elevations (field 0)
    dune_elev_dam = dune_ts[:, :, 0]
    dune_elev_m = dune_elev_dam * DAM_TO_M

    # Calculate statistics at each time step
    max_height = np.nanmax(dune_elev_m, axis=1)
    mean_height = np.nanmean(dune_elev_m, axis=1)
    width = np.sum(~np.isnan(dune_elev_m), axis=1)  # number of non-nan cells

    years = np.arange(n_t) + START_YEAR

    return years, max_height, mean_height, width, dune_elev_m


def check_dune_parameters(seg):
    """Extract dune-related parameters from a segment."""
    params = {}

    dune_params = [
        'Dmax', '_Dmax',  # Maximum dune height
        'DuneWidth', '_DuneWidth',  # Dune width
        'HdDiffu', '_HdDiffu',  # Dune diffusion coefficient
        'rmin', '_rmin', 'rmax', '_rmax',  # Growth rates
        'DuneRestart', '_DuneRestart',  # Restart height
        'Dmaxel', '_Dmaxel',  # Max dune elevation
        'BermEl', '_BermEl',  # Berm elevation
    ]

    for param in dune_params:
        if hasattr(seg, param):
            val = getattr(seg, param)
            if isinstance(val, (int, float, np.integer, np.floating)):
                params[param] = val

    return params


def plot_dune_evolution_dashboard(cascade, sample_domains):
    """Create dashboard showing dune evolution across multiple domains."""
    n_domains = len(sample_domains)
    fig = plt.figure(figsize=(16, 4 * n_domains))

    for idx, dom_idx in enumerate(sample_domains):
        # Analyze this domain
        years, max_h, mean_h, width, dune_profile = analyze_dune_evolution(cascade, dom_idx)

        # Get parameters
        params = check_dune_parameters(cascade.barrier3d[dom_idx])

        # Create subplots for this domain
        ax1 = plt.subplot(n_domains, 3, idx * 3 + 1)
        ax2 = plt.subplot(n_domains, 3, idx * 3 + 2)
        ax3 = plt.subplot(n_domains, 3, idx * 3 + 3)

        # Plot 1: Dune height evolution
        ax1.plot(years, max_h, 'o-', linewidth=2, label='Max height')
        ax1.plot(years, mean_h, 's--', alpha=0.7, label='Mean height')

        # Add horizontal lines for parameters
        if 'Dmax' in params or '_Dmax' in params:
            dmax = params.get('Dmax', params.get('_Dmax', None))
            if dmax is not None:
                dmax_m = dmax * DAM_TO_M
                ax1.axhline(dmax_m, color='red', linestyle='--',
                            label=f'Dmax = {dmax_m:.2f}m', alpha=0.7)

        if 'Dmaxel' in params or '_Dmaxel' in params:
            dmaxel = params.get('Dmaxel', params.get('_Dmaxel', None))
            if dmaxel is not None:
                dmaxel_m = dmaxel * DAM_TO_M
                ax1.axhline(dmaxel_m, color='orange', linestyle='--',
                            label=f'Dmaxel = {dmaxel_m:.2f}m', alpha=0.7)

        ax1.set_xlabel('Year')
        ax1.set_ylabel('Dune height (m)')
        ax1.set_title(f'Domain {dom_idx}: Dune Height Evolution')
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3)

        # Plot 2: Growth rate
        growth_rate = np.diff(max_h)
        ax2.bar(years[1:], growth_rate, alpha=0.7, width=0.8)
        ax2.axhline(0, color='black', linewidth=1)
        ax2.set_xlabel('Year')
        ax2.set_ylabel('Annual height change (m/yr)')
        ax2.set_title(f'Domain {dom_idx}: Growth Rate')
        ax2.grid(alpha=0.3, axis='y')

        # Plot 3: Dune profile evolution (heatmap)
        im = ax3.imshow(dune_profile.T, aspect='auto', origin='lower',
                        extent=[years[0], years[-1], 0, dune_profile.shape[1]],
                        cmap='terrain', vmin=0, vmax=5)
        ax3.set_xlabel('Year')
        ax3.set_ylabel('Cross-shore position')
        ax3.set_title(f'Domain {dom_idx}: Dune Profile Evolution')
        plt.colorbar(im, ax=ax3, label='Elevation (m)')

    plt.tight_layout()
    plt.show()


def print_parameter_summary(cascade, domain_indices):
    """Print summary of dune parameters across domains."""
    print("\n" + "=" * 80)
    print("DUNE PARAMETER SUMMARY")
    print("=" * 80)

    all_params = {}
    for dom_idx in domain_indices:
        params = check_dune_parameters(cascade.barrier3d[dom_idx])
        for key, val in params.items():
            if key not in all_params:
                all_params[key] = []
            all_params[key].append(val)

    print("\nParameter values across analyzed domains:")
    print("-" * 80)
    print(f"{'Parameter':<20} {'Min':<12} {'Max':<12} {'Mean':<12} {'Unit':<10}")
    print("-" * 80)

    for param, values in all_params.items():
        vals_array = np.array(values)
        unit = "dam" if param.endswith(('max', 'maxel', 'El', 'Width')) else ""
        unit_display = f"({unit})" if unit else ""

        print(f"{param:<20} {np.min(vals_array):<12.4f} {np.max(vals_array):<12.4f} "
              f"{np.mean(vals_array):<12.4f} {unit_display:<10}")

        # Also show in meters if applicable
        if unit == "dam":
            print(f"  → in meters:       {np.min(vals_array) * DAM_TO_M:<12.2f} "
                  f"{np.max(vals_array) * DAM_TO_M:<12.2f} {np.mean(vals_array) * DAM_TO_M:<12.2f} (m)")

    print("=" * 80)


def identify_problematic_growth(cascade, domain_indices):
    """Identify domains with excessive dune growth."""
    print("\n" + "=" * 80)
    print("DUNE GROWTH ANALYSIS")
    print("=" * 80)

    results = []

    for dom_idx in domain_indices:
        years, max_h, mean_h, width, _ = analyze_dune_evolution(cascade, dom_idx)

        initial_height = max_h[0]
        final_height = max_h[-1]
        total_growth = final_height - initial_height
        growth_rate = total_growth / (years[-1] - years[0])

        results.append({
            'domain': dom_idx,
            'initial_m': initial_height,
            'final_m': final_height,
            'growth_m': total_growth,
            'rate_m_per_yr': growth_rate
        })

    # Sort by total growth
    results.sort(key=lambda x: x['growth_m'], reverse=True)

    print("\nTop 10 domains by total dune growth:")
    print("-" * 80)
    print(f"{'Domain':<10} {'Initial (m)':<15} {'Final (m)':<15} {'Growth (m)':<15} {'Rate (m/yr)':<15}")
    print("-" * 80)

    for r in results[:10]:
        print(f"{r['domain']:<10} {r['initial_m']:<15.2f} {r['final_m']:<15.2f} "
              f"{r['growth_m']:<15.2f} {r['rate_m_per_yr']:<15.3f}")

    # Check for unrealistic growth
    excessive = [r for r in results if r['rate_m_per_yr'] > 0.5]
    if excessive:
        print(f"\n⚠️  WARNING: {len(excessive)} domains show excessive growth (>0.5 m/yr)")
        print("    This suggests dune growth parameters may be unrealistic.")
        print("    Typical dune vertical growth rates: 0.1-0.3 m/yr")

    # Check initial conditions
    avg_initial = np.mean([r['initial_m'] for r in results])
    print(f"\nAverage initial dune height: {avg_initial:.2f} m")
    if avg_initial > 4:
        print("⚠️  WARNING: Initial dune heights are very high (>4m)")
        print("    Check your initial DuneDomain configuration.")

    print("=" * 80)


def main():
    print(f"Loading CASCADE from: {NPZ_FILE}")
    cascade = load_cascade(NPZ_FILE)

    # Analyze subset of domains
    domain_indices = list(range(DOMAIN_START, min(DOMAIN_END + 1, len(cascade.barrier3d))))

    # Select representative domains to plot (avoid too many plots)
    sample_domains = [
        domain_indices[0],  # First
        domain_indices[len(domain_indices) // 4],  # Quarter
        domain_indices[len(domain_indices) // 2],  # Middle
        domain_indices[3 * len(domain_indices) // 4],  # Three-quarters
        domain_indices[-1]  # Last
    ]

    # Print parameter summary
    print_parameter_summary(cascade, domain_indices)

    # Analyze growth patterns
    identify_problematic_growth(cascade, domain_indices)

    # Plot evolution for sample domains
    print(f"\nPlotting dune evolution for domains: {sample_domains}")
    plot_dune_evolution_dashboard(cascade, sample_domains)


if __name__ == "__main__":
    main()