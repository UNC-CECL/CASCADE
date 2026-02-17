"""
CASCADE Shoreline Sign Convention Diagnostic
============================================
This script loads CASCADE NPZ output files and compares all shoreline-related
variables at known control domains to definitively determine sign convention.

Run this against your HAT_1978_1997 output to resolve the Rodanthe sign question.

Usage:
    python cascade_sign_diagnostic.py

Set the paths in SECTION 1 below, then run.

Author: Diagnostic tool for Hannah Henry / UNC Chapel Hill
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# =============================================================================
# SECTION 1: CONFIGURE THESE PATHS AND DOMAINS
# =============================================================================

# Path to your CASCADE run output directory (the folder with the .npz files)
OUTPUT_DIR = r'C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1978_1997_wh2.5_wa0.1'

# Number of buffer domains on each side
NUM_BUFFER = 15

# --- Control domains (real domain numbers, 1-indexed) ---
# Pick domains where you are CERTAIN of the sign from DSAS observations.
# These are the "anchor points" to ground-truth the sign convention.

CONTROL_DOMAINS = {
    # Domain number : (label, expected_sign, dsas_rate_m_per_yr)
    # Expected sign: 'accretion' or 'erosion'
    # Buxton: known accretion zone (Cape Hatteras headland effect)
    65: ('South of Salvo (Cape Hatteras accretion)', 'accretion', +4),

    # Rodanthe: known erosion zone (the problem area)
    85: ('Rodanthe (known erosion)', 'erosion', -1.5),

    # Oregon Inlet vicinity: known strong erosion
    10: ('Buxton', 'erosion', -1.5),
}

# DSAS observed total shoreline change over the period (m), for comparison
# Positive = accretion (seaward movement), Negative = erosion (landward)
# Update these with your actual DSAS values for 1978-1997
DSAS_OBSERVED_CHANGE_M = {
    65: +76,   # Example: ~+x m/yr * 19 yr
    85: -28.0,   # Example: ~-x m/yr * 19 yr
    10: -28.0,   # Example: ~-x m/yr * 19 yr
}

# =============================================================================
# SECTION 2: LOAD OUTPUT AND EXTRACT SHORELINE VARIABLES
# =============================================================================

def load_domain(output_dir, domain_number, num_buffer):
    """Load NPZ file for a real domain (1-indexed domain number)."""
    # Array index = domain_number - 1 + num_buffer
    array_index = (domain_number - 1) + num_buffer
    
    # Try common file naming conventions
    patterns = [
        f'barrier3d-domain-{array_index}.npz',
        f'cascade-domain-{array_index}.npz',
        f'domain_{array_index:03d}.npz',
    ]
    
    for pattern in patterns:
        fpath = os.path.join(output_dir, pattern)
        if os.path.exists(fpath):
            print(f"  ✓ Found: {pattern}")
            return np.load(fpath, allow_pickle=True)
    
    # If nothing found, list what's actually there
    files = [f for f in os.listdir(output_dir) if f.endswith('.npz')]
    print(f"  ✗ No NPZ found for domain {domain_number} (array index {array_index})")
    print(f"    Available NPZ files (first 5): {files[:5]}")
    return None


def extract_shoreline_variables(npz_data, domain_label):
    """
    Extract all shoreline-related variables and compute change with BOTH sign conventions.
    Returns a dict with findings.
    """
    results = {'label': domain_label, 'variables': {}}
    
    print(f"\n  --- Available keys in NPZ ---")
    available_keys = list(npz_data.keys())
    shoreline_keys = [k for k in available_keys if any(
        term in k.lower() for term in ['shore', 'x_s', 'scr', 'change']
    )]
    print(f"  Shoreline-related keys: {shoreline_keys}")
    print(f"  All keys: {available_keys[:20]}{'...' if len(available_keys) > 20 else ''}")
    
    # --- Extract x_s_TS (raw shoreline position) ---
    if 'x_s_TS' in npz_data:
        x_s_ts = np.array(npz_data['x_s_TS'])
        # Remove initial condition (index 0)
        x_s_change = x_s_ts[-1] - x_s_ts[0]
        
        print(f"\n  x_s_TS (raw shoreline position, in DAM):")
        print(f"    Initial x_s: {x_s_ts[0]:.4f} dam")
        print(f"    Final x_s:   {x_s_ts[-1]:.4f} dam")
        print(f"    Change (x_s[-1] - x_s[0]): {x_s_change:.4f} dam = {x_s_change*10:.2f} m")
        print(f"    >> x_s INCREASES landward (erosion direction)")
        print(f"    >> Direct change is POSITIVE for erosion (OPPOSITE of coastal convention)")
        
        results['variables']['x_s_TS'] = {
            'initial': float(x_s_ts[0]),
            'final': float(x_s_ts[-1]),
            'raw_change_dam': float(x_s_change),
            'raw_change_m': float(x_s_change * 10),
            'corrected_change_m': float(-x_s_change * 10),  # flip sign for coastal convention
            'sign_as_stored': 'POSITIVE=erosion (must negate for coastal convention)',
        }
    else:
        print("  ✗ x_s_TS not found in NPZ")
    
    # --- Extract ShorelineChangeTS (per-year change, should be conventional sign) ---
    if 'ShorelineChangeTS' in npz_data:
        sc_ts = np.array(npz_data['ShorelineChangeTS'])
        cumulative = np.sum(sc_ts)
        
        print(f"\n  ShorelineChangeTS (per-year change, in DAM cells):")
        print(f"    Length: {len(sc_ts)} time steps")
        print(f"    Cumulative sum: {cumulative:.4f} dam = {cumulative*10:.2f} m")
        print(f"    >> Positive = ACCRETION (conventional coastal sign)")
        
        results['variables']['ShorelineChangeTS'] = {
            'cumulative_dam': float(cumulative),
            'cumulative_m': float(cumulative * 10),
            'sign_convention': 'POSITIVE=accretion (conventional)',
        }
    else:
        print("  ✗ ShorelineChangeTS not found in NPZ")

    # --- Extract ShorelineChange (total cumulative) ---
    if 'ShorelineChange' in npz_data:
        sc_total = float(npz_data['ShorelineChange'])
        print(f"\n  ShorelineChange (total cumulative, DAM):")
        print(f"    Value: {sc_total:.4f} dam = {sc_total*10:.2f} m")
        print(f"    >> Positive = ACCRETION (conventional coastal sign)")
        
        results['variables']['ShorelineChange'] = {
            'value_dam': sc_total,
            'value_m': sc_total * 10,
            'sign_convention': 'POSITIVE=accretion (conventional)',
        }
    else:
        print("  ✗ ShorelineChange not found in NPZ")

    return results


# =============================================================================
# SECTION 3: RUN COMPARISON AND PLOT
# =============================================================================

def run_diagnostic():
    print("=" * 70)
    print("CASCADE SHORELINE SIGN CONVENTION DIAGNOSTIC")
    print("=" * 70)
    print(f"Output directory: {OUTPUT_DIR}\n")
    
    if not os.path.exists(OUTPUT_DIR):
        print(f"ERROR: Output directory not found: {OUTPUT_DIR}")
        print("Update OUTPUT_DIR in SECTION 1.")
        return
    
    all_results = {}
    
    for domain_num, (label, expected, dsas_rate) in CONTROL_DOMAINS.items():
        print(f"\n{'='*60}")
        print(f"DOMAIN {domain_num}: {label}")
        print(f"  Expected: {expected} | DSAS rate: {dsas_rate:+.1f} m/yr")
        print(f"{'='*60}")
        
        data = load_domain(OUTPUT_DIR, domain_num, NUM_BUFFER)
        if data is None:
            continue
        
        results = extract_shoreline_variables(data, label)
        results['expected'] = expected
        results['dsas_total_m'] = DSAS_OBSERVED_CHANGE_M.get(domain_num, None)
        all_results[domain_num] = results
    
    # --- Summary comparison table ---
    print("\n" + "=" * 70)
    print("SIGN CONVENTION SUMMARY")
    print("=" * 70)
    print(f"{'Domain':<8} {'Label':<35} {'Expected':<12} {'x_s change':>12} {'SC change':>12} {'DSAS':>10}")
    print(f"{'':8} {'':35} {'':12} {'(m, raw)':>12} {'(m)':>12} {'(m)':>10}")
    print("-" * 95)
    
    for domain_num, res in all_results.items():
        label = CONTROL_DOMAINS[domain_num][0][:34]
        expected = res['expected']
        dsas = res.get('dsas_total_m', 'N/A')
        
        xs_change = res['variables'].get('x_s_TS', {}).get('raw_change_m', 'N/A')
        sc_change = res['variables'].get('ShorelineChange', {}).get('value_m',
                    res['variables'].get('ShorelineChangeTS', {}).get('cumulative_m', 'N/A'))
        
        if isinstance(xs_change, float):
            xs_str = f"{xs_change:+.1f}"
        else:
            xs_str = str(xs_change)
            
        if isinstance(sc_change, float):
            sc_str = f"{sc_change:+.1f}"
        else:
            sc_str = str(sc_change)
            
        dsas_str = f"{dsas:+.1f}" if isinstance(dsas, (int, float)) else str(dsas)
        
        print(f"{domain_num:<8} {label:<35} {expected:<12} {xs_str:>12} {sc_str:>12} {dsas_str:>10}")
    
    print("\n>> INTERPRETATION:")
    print("   If x_s raw change sign MATCHES DSAS → plotting script must NEGATE x_s for coastal convention")
    print("   If x_s raw change sign OPPOSES DSAS  → x_s sign is already wrong; use ShorelineChangeTS instead")
    print("   If ShorelineChange sign MATCHES DSAS → use ShorelineChange/ShorelineChangeTS for plotting")
    
    # --- Create comparison figure ---
    _make_figure(all_results)


def _make_figure(all_results):
    """Plot per-domain time series comparison."""
    n_domains = len(all_results)
    if n_domains == 0:
        return
    
    fig = plt.figure(figsize=(14, 4 * n_domains))
    gs = gridspec.GridSpec(n_domains, 2, figure=fig, wspace=0.35, hspace=0.45)
    
    for i, (domain_num, res) in enumerate(all_results.items()):
        label = CONTROL_DOMAINS[domain_num][0]
        expected = res['expected']
        dsas_total = res.get('dsas_total_m', None)
        
        # Left panel: x_s_TS time series
        ax1 = fig.add_subplot(gs[i, 0])
        if 'x_s_TS' in res['variables']:
            # We stored only summary stats; note - in production you'd store the full array
            xs_change = res['variables']['x_s_TS']['raw_change_m']
            xs_corrected = res['variables']['x_s_TS']['corrected_change_m']
            ax1.bar(['x_s raw\nchange (m)', 'x_s corrected\n(-1 × raw)'],
                    [xs_change, xs_corrected],
                    color=['steelblue', 'darkorange'])
            if dsas_total is not None:
                ax1.axhline(dsas_total, color='red', linewidth=2, linestyle='--',
                           label=f'DSAS observed: {dsas_total:+.1f} m')
                ax1.legend(fontsize=9)
            ax1.axhline(0, color='k', linewidth=0.5)
            ax1.set_ylabel('Total change (m)')
            ax1.set_title(f'Domain {domain_num}: x_s_TS\n{label}', fontsize=9)
            color = 'green' if expected == 'accretion' else 'red'
            ax1.text(0.02, 0.95, f'Expected: {expected}', transform=ax1.transAxes,
                    color=color, fontsize=9, va='top', fontweight='bold')
        
        # Right panel: ShorelineChange
        ax2 = fig.add_subplot(gs[i, 1])
        sc_val = res['variables'].get('ShorelineChange', {}).get('value_m',
                 res['variables'].get('ShorelineChangeTS', {}).get('cumulative_m', None))
        
        if sc_val is not None:
            ax2.bar(['ShorelineChange\n(m)'], [sc_val], color='mediumpurple')
            if dsas_total is not None:
                ax2.axhline(dsas_total, color='red', linewidth=2, linestyle='--',
                           label=f'DSAS observed: {dsas_total:+.1f} m')
                ax2.legend(fontsize=9)
            ax2.axhline(0, color='k', linewidth=0.5)
            ax2.set_ylabel('Total change (m)')
            ax2.set_title(f'Domain {domain_num}: ShorelineChangeTS\n{label}', fontsize=9)
            color = 'green' if expected == 'accretion' else 'red'
            ax2.text(0.02, 0.95, f'Expected: {expected}', transform=ax2.transAxes,
                    color=color, fontsize=9, va='top', fontweight='bold')
    
    fig.suptitle('CASCADE Sign Convention Diagnostic\n'
                 'Red dashed = DSAS observed | Bars = CASCADE output\n'
                 'Which bar matches DSAS? → That is the correct variable to use.',
                 fontsize=11, y=1.02)
    
    out_fig = os.path.join(os.path.dirname(OUTPUT_DIR), 'sign_diagnostic.png')
    plt.savefig(out_fig, dpi=150, bbox_inches='tight')
    print(f"\n✓ Figure saved to: {out_fig}")
    plt.show()


# =============================================================================
# SECTION 4: QUICK TIME-SERIES PLOTTER
# (For a single domain - detailed x_s_TS trajectory)
# =============================================================================

def plot_xsTS_timeseries(domain_number, run_years=19):
    """
    Load a single domain's x_s_TS and plot both the raw and sign-corrected version.
    This shows you the full evolution trajectory, not just start/end.
    
    Useful for: visually confirming which sign matches reality.
    """
    print(f"\nLoading domain {domain_number} for time-series plot...")
    data = load_domain(OUTPUT_DIR, domain_number, NUM_BUFFER)
    if data is None:
        return
    
    if 'x_s_TS' not in data:
        print("x_s_TS not found in NPZ")
        return
    
    x_s = np.array(data['x_s_TS'])
    years = np.arange(len(x_s))
    
    # Compute change relative to initial (in meters)
    x_s_change_raw = (x_s - x_s[0]) * 10          # raw: positive = landward = erosion
    x_s_change_corrected = -x_s_change_raw          # corrected: positive = seaward = accretion
    
    sc_ts = np.array(data.get('ShorelineChangeTS', np.zeros(len(x_s))))
    sc_cumulative = np.cumsum(sc_ts) * 10            # in meters, positive = accretion
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    
    axes[0].plot(years, x_s * 10, 'b-o', markersize=3, label='x_s (m)')
    axes[0].set_ylabel('x_s position (m)\nIncreases LANDWARD')
    axes[0].set_title(f'Domain {domain_number}: Raw x_s_TS\n'
                     '(Slope UP = erosion; slope DOWN = accretion)')
    axes[0].legend(); axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(years, x_s_change_raw, 'r-o', markersize=3, label='x_s change (raw)')
    axes[1].plot(years, x_s_change_corrected, 'g-o', markersize=3, label='x_s change (×-1 = coastal convention)')
    axes[1].axhline(0, color='k', linewidth=0.5)
    axes[1].set_ylabel('Cumulative change (m)')
    axes[1].set_title('x_s change: Raw vs. Corrected\n'
                     'Green (corrected) should match DSAS sign: + = accretion')
    axes[1].legend(); axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(years, sc_cumulative, 'm-o', markersize=3, label='ShorelineChangeTS cumulative')
    axes[2].axhline(0, color='k', linewidth=0.5)
    axes[2].set_ylabel('Cumulative change (m)\n+ = accretion')
    axes[2].set_xlabel('Year')
    axes[2].set_title('ShorelineChangeTS (cumulative sum)\n'
                     'Already conventional: + = accretion')
    axes[2].legend(); axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    label = CONTROL_DOMAINS.get(domain_number, (f'Domain {domain_number}',))[0]
    plt.suptitle(f'Sign Diagnostic Time Series: {label}', y=1.01, fontsize=11)
    
    out_fig = os.path.join(os.path.dirname(OUTPUT_DIR), f'sign_diagnostic_domain_{domain_number}.png')
    plt.savefig(out_fig, dpi=150, bbox_inches='tight')
    print(f"✓ Time series figure saved to: {out_fig}")
    plt.show()


# =============================================================================
# RUN
# =============================================================================

if __name__ == '__main__':
    # Main diagnostic: compare control domains
    run_diagnostic()
    
    # Optional: plot full time series for one domain
    # Uncomment and set domain number to inspect a single domain in detail
    # plot_xsTS_timeseries(domain_number=70)  # Rodanthe
    # plot_xsTS_timeseries(domain_number=50)  # Buxton accretion
