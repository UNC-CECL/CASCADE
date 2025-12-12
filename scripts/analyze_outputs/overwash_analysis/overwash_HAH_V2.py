"""
Overwash Diagnostic Tool

Investigates why overwash is not occurring in CASCADE model runs.
Checks storm intensity vs dune heights across all domains.
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

# Domain range to analyze (set to None to analyze all)
DOMAIN_START = 15  # First real island domain (after left buffer)
DOMAIN_END = 104  # Last real island domain (before right buffer)


# ==========================
# FUNCTIONS
# ==========================

def load_cascade(npz_path):
    """Load CASCADE object from .npz file."""
    data = np.load(npz_path, allow_pickle=True)
    return data["cascade"][0]


def get_dune_crest_ts(segment):
    """Get dune crest elevation time series in meters."""
    if hasattr(segment, "DuneDomain"):
        dune = np.array(segment.DuneDomain, dtype=float)
    elif hasattr(segment, "_DuneDomain"):
        dune = np.array(segment._DuneDomain, dtype=float)
    else:
        raise AttributeError("No DuneDomain found")

    dune_elev_dam = dune[:, :, 0]
    crest_ts_dam = np.nanmax(dune_elev_dam, axis=1)
    return crest_ts_dam * DAM_TO_M


def get_storm_series(segment):
    """Get storm series array."""
    if hasattr(segment, "StormSeries"):
        return np.array(segment.StormSeries, dtype=float)
    elif hasattr(segment, "_StormSeries"):
        return np.array(segment._StormSeries, dtype=float)
    else:
        raise AttributeError("No StormSeries found")


def get_qow_ts(segment):
    """Get overwash flux time series."""
    if hasattr(segment, "QowTS"):
        return np.array(segment.QowTS, dtype=float)
    elif hasattr(segment, "_QowTS"):
        return np.array(segment._QowTS, dtype=float)
    else:
        raise AttributeError("No QowTS found")


def analyze_overwash_potential(cascade, domain_start=None, domain_end=None):
    """
    Analyze overwash potential across domains.

    Returns statistics on:
    - Dune heights
    - Storm intensities
    - Freeboard (storm - dune)
    - Actual overwash occurrence
    """
    barrier_models = cascade.barrier3d
    n_domains = len(barrier_models)

    if domain_start is None:
        domain_start = 0
    if domain_end is None:
        domain_end = n_domains

    domain_range = range(domain_start, min(domain_end + 1, n_domains))

    # Storage for results
    results = {
        'domain_ids': [],
        'mean_dune_height': [],
        'max_dune_height': [],
        'min_dune_height': [],
        'max_storm_rhigh': [],
        'mean_storm_rhigh': [],
        'n_storms': [],
        'max_freeboard': [],
        'mean_freeboard': [],
        'total_overwash': [],
        'n_overwash_events': [],
        'overwash_years': []
    }

    print(f"\nAnalyzing domains {domain_start} to {domain_end}...")
    print("=" * 80)

    for dom_idx in domain_range:
        seg = barrier_models[dom_idx]

        # Get data
        try:
            crest_ts_m = get_dune_crest_ts(seg)
            storms = get_storm_series(seg)
            qow_ts = get_qow_ts(seg)

            # Storm data
            if storms.size > 0:
                storm_rhigh_m = storms[:, 1] * DAM_TO_M
                storm_year_idx = storms[:, 0].astype(int)

                # Filter valid storms
                n_t = len(crest_ts_m)
                valid_mask = (storm_year_idx >= 0) & (storm_year_idx < n_t)
                storm_year_idx_valid = storm_year_idx[valid_mask]
                storm_rhigh_valid_m = storm_rhigh_m[valid_mask]

                # Aggregate max Rhigh per year
                yearly_rhigh_m = np.full(n_t, np.nan)
                for y_idx, rh_m in zip(storm_year_idx_valid, storm_rhigh_valid_m):
                    if np.isnan(yearly_rhigh_m[y_idx]):
                        yearly_rhigh_m[y_idx] = rh_m
                    else:
                        yearly_rhigh_m[y_idx] = max(yearly_rhigh_m[y_idx], rh_m)

                # Freeboard (where we have storms)
                freeboard_m = yearly_rhigh_m - crest_ts_m
                valid_freeboard = freeboard_m[~np.isnan(freeboard_m)]

                # Overwash analysis
                overwash_events = np.sum(qow_ts > 0)
                total_overwash = np.sum(qow_ts)
                overwash_year_indices = np.where(qow_ts > 0)[0]

                # Store results
                results['domain_ids'].append(dom_idx)
                results['mean_dune_height'].append(np.mean(crest_ts_m))
                results['max_dune_height'].append(np.max(crest_ts_m))
                results['min_dune_height'].append(np.min(crest_ts_m))
                results['max_storm_rhigh'].append(np.max(storm_rhigh_valid_m))
                results['mean_storm_rhigh'].append(np.mean(storm_rhigh_valid_m))
                results['n_storms'].append(len(storm_rhigh_valid_m))
                results['max_freeboard'].append(np.max(valid_freeboard) if len(valid_freeboard) > 0 else np.nan)
                results['mean_freeboard'].append(np.mean(valid_freeboard) if len(valid_freeboard) > 0 else np.nan)
                results['total_overwash'].append(total_overwash)
                results['n_overwash_events'].append(overwash_events)
                results['overwash_years'].append(overwash_year_indices.tolist())

        except Exception as e:
            print(f"Error processing domain {dom_idx}: {e}")
            continue

    # Convert to arrays
    for key in results:
        if key != 'overwash_years':
            results[key] = np.array(results[key])

    return results


def print_summary_statistics(results):
    """Print summary statistics."""
    print("\n" + "=" * 80)
    print("OVERWASH DIAGNOSTIC SUMMARY")
    print("=" * 80)

    print(f"\n--- DUNE HEIGHTS (meters) ---")
    print(f"Mean dune height across domains: {np.mean(results['mean_dune_height']):.3f} m")
    print(
        f"Range of mean dune heights: {np.min(results['mean_dune_height']):.3f} - {np.max(results['mean_dune_height']):.3f} m")
    print(f"Max dune height observed: {np.max(results['max_dune_height']):.3f} m")

    print(f"\n--- STORM INTENSITIES (meters) ---")
    print(f"Mean Rhigh across all storms: {np.mean(results['mean_storm_rhigh']):.3f} m")
    print(f"Max Rhigh across all storms: {np.max(results['max_storm_rhigh']):.3f} m")
    print(f"Total storms per domain: {np.mean(results['n_storms']):.1f} (avg)")

    print(f"\n--- FREEBOARD (Rhigh - Dune Height, meters) ---")
    valid_freeboard = results['max_freeboard'][~np.isnan(results['max_freeboard'])]
    print(f"Max freeboard across domains: {np.max(valid_freeboard):.3f} m")
    print(f"Mean max freeboard: {np.mean(valid_freeboard):.3f} m")

    if np.max(valid_freeboard) < 0:
        print("\n⚠️  CRITICAL: Maximum freeboard is NEGATIVE!")
        print("    This means NO storms are tall enough to overtop any dunes.")
        print("    Overwash requires positive freeboard (Rhigh > dune height).")
    elif np.max(valid_freeboard) < 0.1:
        print("\n⚠️  WARNING: Maximum freeboard is very small (<10cm).")
        print("    Storms are barely overtopping dunes. Limited overwash expected.")
    else:
        print("\n✓ Freeboard values suggest overwash SHOULD be occurring.")

    print(f"\n--- ACTUAL OVERWASH ---")
    domains_with_overwash = np.sum(results['n_overwash_events'] > 0)
    print(f"Domains with overwash: {domains_with_overwash} / {len(results['domain_ids'])}")
    print(f"Total overwash events: {np.sum(results['n_overwash_events'])}")
    print(f"Total Qow magnitude: {np.sum(results['total_overwash']):.6f}")

    if domains_with_overwash == 0:
        print("\n⚠️  NO OVERWASH OCCURRING IN ANY DOMAIN!")

    # Find domains with positive freeboard but no overwash
    positive_freeboard_mask = results['max_freeboard'] > 0
    no_overwash_mask = results['n_overwash_events'] == 0
    problematic = positive_freeboard_mask & no_overwash_mask

    if np.any(problematic):
        prob_domains = results['domain_ids'][problematic]
        print(f"\n⚠️  {len(prob_domains)} domains have positive freeboard but NO overwash:")
        print(f"    Domains: {prob_domains[:10]}{'...' if len(prob_domains) > 10 else ''}")
        print("    This suggests a model configuration issue, not just storm/dune mismatch.")

    print("=" * 80)


def plot_overwash_diagnostics(results):
    """Create diagnostic plots."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("Overwash Diagnostic Analysis", fontsize=16, fontweight='bold')

    domains = results['domain_ids']

    # Plot 1: Dune heights vs Storm intensities
    ax1 = axes[0, 0]
    ax1.plot(domains, results['mean_dune_height'], 'o-', label='Mean dune height', linewidth=2)
    ax1.plot(domains, results['max_dune_height'], 's--', alpha=0.6, label='Max dune height')
    ax1.plot(domains, results['max_storm_rhigh'], '^-', label='Max storm Rhigh', linewidth=2, color='red')
    ax1.set_xlabel("Domain index")
    ax1.set_ylabel("Elevation (m)")
    ax1.set_title("Dune Heights vs Storm Intensity")
    ax1.legend()
    ax1.grid(alpha=0.3)

    # Plot 2: Freeboard
    ax2 = axes[0, 1]
    ax2.plot(domains, results['max_freeboard'], 'o-', linewidth=2)
    ax2.axhline(0, color='red', linestyle='--', linewidth=2, label='Zero freeboard')
    ax2.fill_between(domains, 0, results['max_freeboard'],
                     where=results['max_freeboard'] > 0, alpha=0.3, color='green',
                     label='Overwash possible')
    ax2.fill_between(domains, results['max_freeboard'], 0,
                     where=results['max_freeboard'] < 0, alpha=0.3, color='red',
                     label='Dunes too high')
    ax2.set_xlabel("Domain index")
    ax2.set_ylabel("Max freeboard (m)")
    ax2.set_title("Freeboard (Storm Rhigh - Dune Height)")
    ax2.legend()
    ax2.grid(alpha=0.3)

    # Plot 3: Overwash occurrence
    ax3 = axes[1, 0]
    ax3.bar(domains, results['n_overwash_events'], alpha=0.7)
    ax3.set_xlabel("Domain index")
    ax3.set_ylabel("Number of overwash events")
    ax3.set_title("Overwash Event Count by Domain")
    ax3.grid(alpha=0.3, axis='y')

    # Plot 4: Scatter - Freeboard vs Overwash
    ax4 = axes[1, 1]
    scatter = ax4.scatter(results['max_freeboard'], results['n_overwash_events'],
                          c=results['max_storm_rhigh'], s=50, alpha=0.6, cmap='viridis')
    ax4.axvline(0, color='red', linestyle='--', linewidth=2, alpha=0.5)
    ax4.set_xlabel("Max freeboard (m)")
    ax4.set_ylabel("Overwash events")
    ax4.set_title("Freeboard vs Overwash Occurrence")
    ax4.grid(alpha=0.3)
    cbar = plt.colorbar(scatter, ax=ax4)
    cbar.set_label("Max Rhigh (m)")

    plt.tight_layout()
    plt.show()


def investigate_specific_domain(cascade, domain_idx):
    """Deep dive into a specific domain."""
    print(f"\n{'=' * 80}")
    print(f"DETAILED INVESTIGATION: Domain {domain_idx}")
    print(f"{'=' * 80}")

    seg = cascade.barrier3d[domain_idx]

    # Get data
    crest_ts_m = get_dune_crest_ts(seg)
    storms = get_storm_series(seg)
    qow_ts = get_qow_ts(seg)

    # Check for threshold parameters
    print("\n--- CASCADE Overwash Parameters ---")
    threshold_params = ['_Rin_r', 'Rin_r', '_threshold_in', 'threshold_in',
                        '_Rin_i', 'Rin_i', '_Ki', 'Ki', '_Kr', 'Kr']
    for param in threshold_params:
        if hasattr(seg, param):
            print(f"{param}: {getattr(seg, param)}")

    # Analyze storms
    if storms.size > 0:
        storm_rhigh_m = storms[:, 1] * DAM_TO_M
        storm_year_idx = storms[:, 0].astype(int)
        n_t = len(crest_ts_m)

        print(f"\n--- Storm Analysis ---")
        print(f"Total storms in series: {len(storms)}")
        print(f"Storm Rhigh range: {storm_rhigh_m.min():.3f} - {storm_rhigh_m.max():.3f} m")

        # Year-by-year analysis
        print(f"\n--- Year-by-Year Comparison ---")
        print(f"{'Year':<8} {'Dune (m)':<12} {'Max Rhigh (m)':<15} {'Freeboard (m)':<15} {'Qow':<12}")
        print("-" * 70)

        for t in range(n_t):
            year_storms = storms[storm_year_idx == t]
            if len(year_storms) > 0:
                max_rhigh = np.max(year_storms[:, 1]) * DAM_TO_M
                freeboard = max_rhigh - crest_ts_m[t]
                print(f"{START_YEAR + t:<8} {crest_ts_m[t]:<12.3f} {max_rhigh:<15.3f} "
                      f"{freeboard:<15.3f} {qow_ts[t]:<12.6f}")


def main():
    # Load CASCADE
    print(f"Loading CASCADE from: {NPZ_FILE}")
    cascade = load_cascade(NPZ_FILE)

    n_domains = len(cascade.barrier3d)
    print(f"Loaded CASCADE with {n_domains} domains")

    # Analyze overwash potential
    results = analyze_overwash_potential(cascade, DOMAIN_START, DOMAIN_END)

    # Print summary
    print_summary_statistics(results)

    # Create plots
    plot_overwash_diagnostics(results)

    # Deep dive into a specific domain (choose one with positive freeboard but no overwash)
    positive_freeboard_mask = results['max_freeboard'] > 0
    no_overwash_mask = results['n_overwash_events'] == 0
    problematic = positive_freeboard_mask & no_overwash_mask

    if np.any(problematic):
        investigate_domain = results['domain_ids'][problematic][0]
        investigate_specific_domain(cascade, investigate_domain)
    elif len(results['domain_ids']) > 0:
        # Just investigate the middle domain
        mid_idx = len(results['domain_ids']) // 2
        investigate_specific_domain(cascade, results['domain_ids'][mid_idx])


if __name__ == "__main__":
    main()