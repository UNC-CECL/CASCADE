"""
HAT_source_sink_bothperiods_fulldomain.py
==========================================
Joint two-period parametric sweep for source/sink correction.

Finds parametric correction values that perform well across BOTH time periods:
  Period 1 (calibration): 1984-2004
  Period 2 (validation):  2004-2024

Evaluation strategy:
  - RMSE computed over ALL 90 domains (no zone exclusions).
  - The optimizer sees the full-island residual pattern and finds the
    best-fit correction without the analyst pre-specifying where it matters.
  - Pareto-optimal solutions show the period-1 vs period-2 trade-off
    across the complete domain.
  - A diagnostic plot of raw CoastSat vs CASCADE LRR (before any correction)
    is included so the data can reveal the spatial signal.

Parametric correction shapes are unchanged from the zone-focused version:
  - Wimble Shoals triangular shadow (domains WS_X_START → WS_X_PEAK → WS_X_END)
  - Rodanthe Gaussian bell (centered on pier at ROD_CENTER)
  - Avon linear ramp (domains AVON_D_SOUTH → AVON_D_NORTH)
These shapes are the hypotheses; the full-domain RMSE tests whether they
actually describe the observed residual.

Outputs:
  - sweep_results_joint_fulldomain.csv    : all parameter combos with P1 RMSE, P2 RMSE, combined
  - pareto_front_fulldomain.csv           : Pareto-optimal solutions
  - joint_best_fit_fulldomain.png         : diagnostic + period comparison + Pareto plot
  - DOMAIN_BE_RATES dict                  : printed, ready to paste into CASCADE run script

Author: Hannah Henry, UNC Chapel Hill
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import product
import os

# =============================================================================
# CONFIGURATION — edit paths and weights here
# =============================================================================

# Period 1: Calibration (1984-2004)
P1_COASTSAT_CSV     = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\1984_2004\domain_lrr_1984_2004_summary.csv"
P1_CASCADE_BASELINE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_1984_2004_BASE_Hs2.5\HAT_1984_2004_BASE_Hs2.5.npz"
P1_START, P1_END    = 1984, 2004
P1_LABEL            = "1984–2004 (calibration; 2.5)"

# Period 2: Validation (2004-2024)
P2_COASTSAT_CSV     = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_2004_2024_summary.csv"
P2_CASCADE_BASELINE = r"C:\Users\hanna\PycharmProjects\CASCADE\output\raw_runs\HAT_2004_2024_BASE_BN_Hs2.5\HAT_2004_2024_BASE_BN_Hs2.5.npz"
P2_START, P2_END    = 2004, 2024
P2_LABEL            = "2004–2024 (test; 2.5)"

# Output directory
OUTPUT_DIR = r"/scripts/input_prep/original/parametric_ss/results_fulldomain_bothperiods_Hs2.5"

# CASCADE npz comparison variable name for shoreline position
CASCADE_SHORELINE_KEY = 'x_s_TS'

# CASCADE x_s_TS sign convention:
#   In Barrier3D, x_s_TS is measured from the LANDWARD side.
#   Positive slope = shoreline retreating LANDWARD = EROSION
#   This is the OPPOSITE of standard convention (where erosion is negative).
#   Set CASCADE_X_S_SIGN = -1 to flip to standard convention before computing residual.
#   Set CASCADE_X_S_SIGN = +1 if your CASCADE measures x_s seaward (accretion = positive).
CASCADE_X_S_SIGN = -1

# CASCADE x_s_TS units:
#   'dam' → x_s_TS in decameters, multiply by 10 to get m/yr  (standard for CASCADE/Barrier3D)
#   'm'   → x_s_TS already in meters, no conversion needed
CASCADE_X_S_UNITS = 'dam'

# CoastSat CSV column names
CS_DOMAIN_COL = 'domain_number'
CS_LRR_COL    = 'mean_lrr'

# CASCADE domain setup
N_DOMAINS = 90
N_PAD     = 15

# Weighting between periods (0.5 = equal; increase P1_WEIGHT to prioritise calibration period)
P1_WEIGHT = 0.5
P2_WEIGHT = 1.0 - P1_WEIGHT

# =============================================================================
# EVALUATION DOMAINS
# Full-domain evaluation: all 90 domains, no exclusions.
# The optimizer grades itself on the complete island — the data drives
# where the correction needs to be large or small.
# =============================================================================

EVAL_DOMAINS   = list(range(1, N_DOMAINS + 1))   # all 90 domains
EXCLUSION_ZONES = set()                            # no exclusions

# =============================================================================
# PHYSICALLY ANCHORED BOUNDARIES — fixed from bathymetric evidence
# These define the SHAPE of the correction; the sweep finds the magnitude.
# =============================================================================

WS_X_START = 38    # Wimble Shoals smooth shelf begins
WS_X_PEAK  = 68    # Wimble Shoals platform directly offshore (6-12 fathom depth)
WS_X_END   = 80    # Northern taper end
ROD_CENTER = 79    # Rodanthe Fishing Pier (fixed — structure location)
AVON_D_SOUTH, AVON_D_NORTH = 25, 30
AVON_VAL_SOUTH = -0.10  # Fixed: slight downdrift erosion south of Avon Pier (dam/yr)

# =============================================================================
# SWEEP RANGES — free parameters to optimize
# Scaled to match observed residual magnitude (~0.1–0.5 dam/yr).
# Run once at coarse resolution, then narrow around best-fit values.
# =============================================================================

WS_PEAK_VALUES    = np.round(np.arange(0.05, 0.70, 0.05), 3)   # Wimble Shoals peak (dam/yr)
ROD_AMPLITUDES    = np.round(np.arange(-0.60, -0.05, 0.05), 3)  # Rodanthe amplitude (dam/yr, negative)
ROD_WIDTHS        = np.arange(3, 8, 1)                           # Rodanthe half-width (domains)
AVON_NORTH_VALUES = np.round(np.arange(0.05, 0.60, 0.05), 3)   # Avon north value (dam/yr)

# Unit multiplier: sweep runs in dam/yr internally; DOMAIN_BE_RATES expects m/yr.
BE_UNIT_MULTIPLIER = 10.0   # dam/yr × 10 = m/yr

# Edge boundary conditions — period- and Hs-specific, set manually.
# 1984-2004: domain 1 = -65 (Hs2.5), -130 (Hs5.0)  | domain 90 = 12 (Hs2.5), 40 (Hs5.0)
# 2004-2024: domain 1 =   6 (Hs2.5),  -40 (Hs5.0)  | domain 90 = 40 (Hs2.5), 75 (Hs5.0)
EDGE_D1  =   6    # <-- update per run
EDGE_D90 =  40    # <-- update per run

# =============================================================================
# PARAMETRIC CORRECTION FUNCTIONS
# =============================================================================

def wimble_shadow(domain, x_start, x_peak, x_end, peak_value):
    """
    Triangular wave shadow (dam/yr).
    Linear ramp: 0 at x_start → peak at x_peak → 0 at x_end.
    Physical basis: diffraction shadow attenuates monotonically from shoal.
    All adjacent values determined by same 4 parameters.
    """
    if domain < x_start or domain > x_end:
        return 0.0
    elif domain <= x_peak:
        return peak_value * (domain - x_start) / (x_peak - x_start)
    else:
        return peak_value * (x_end - domain) / (x_end - x_peak)


def rodanthe_bell(domain, center, width, amplitude):
    """
    Gaussian erosion bell (dam/yr, amplitude negative).
    Physical basis: shoreline curvature → divergent transport;
    Rodanthe Pier concentrates erosion at domain 79.
    """
    return amplitude * np.exp(-0.5 * ((domain - center) / width) ** 2)


def avon_correction(domain, d_south, d_north, val_south, val_north):
    """
    Linear ramp for Avon zone (dam/yr).
    Physical basis: offshore shoal features + Avon Pier groin effect.
    """
    if domain < d_south or domain > d_north:
        return 0.0
    t = (domain - d_south) / (d_north - d_south)
    return val_south + t * (val_north - val_south)


def build_correction(domains, ws_peak, rod_amp, rod_width, avon_north):
    """Assemble total correction array from parameters."""
    return {
        d: (wimble_shadow(d, WS_X_START, WS_X_PEAK, WS_X_END, ws_peak)
            + rodanthe_bell(d, ROD_CENTER, rod_width, rod_amp)
            + avon_correction(d, AVON_D_SOUTH, AVON_D_NORTH,
                              AVON_VAL_SOUTH, avon_north))
        for d in domains
    }

# =============================================================================
# DATA LOADING
# =============================================================================

def load_coastsat(csv_path, domain_col=CS_DOMAIN_COL, lrr_col=CS_LRR_COL):
    """Load CoastSat LRR CSV → dict {domain: lrr_m_per_yr}."""
    df = pd.read_csv(csv_path)
    print(f"  Columns found: {df.columns.tolist()}")
    print(f"  Rows: {len(df)}  |  Looking for '{domain_col}' and '{lrr_col}'")
    if domain_col not in df.columns or lrr_col not in df.columns:
        raise KeyError(
            f"\nColumn mismatch in {csv_path}\n"
            f"  Expected: '{domain_col}' and '{lrr_col}'\n"
            f"  Found:    {df.columns.tolist()}\n"
            f"  → Update CS_DOMAIN_COL and CS_LRR_COL at the top of the script."
        )
    return dict(zip(df[domain_col].astype(int), df[lrr_col].astype(float)))


def extract_cascade_lrr(npz_path, start_year, end_year,
                         n_domains=N_DOMAINS, n_pad=N_PAD,
                         shoreline_key=CASCADE_SHORELINE_KEY):
    """
    Extract per-domain LRR (m/yr) from CASCADE comparison npz.

    CASCADE npz structure:
      comparison['cascade'][0]              → cascade object
      cascade.barrier3d                 → list of Barrier3D objects (all domains inc. padding)
      cascade.barrier3d[i].x_s_TS      → shoreline position time series for domain i
      Real domains: indices n_pad to n_pad + n_domains (1-indexed: 1 to n_domains)
    """
    data = np.load(npz_path, allow_pickle=True)
    print(f"  NPZ keys found: {list(data.keys())}")

    # Unpack cascade object — stored as array, take index [0]
    cascade_obj = data['cascade'][0]
    barrier3d   = cascade_obj.barrier3d
    n_total     = len(barrier3d)
    print(f"  Total barrier3d domains (inc. padding): {n_total}")
    print(f"  Real domains: indices {n_pad} to {n_pad + n_domains - 1} "
          f"(CASCADE domains 1–{n_domains})")

    # Check shoreline attribute exists
    test_b3d = barrier3d[n_pad]
    available = [a for a in dir(test_b3d) if not a.startswith('_')]
    if not hasattr(test_b3d, shoreline_key):
        raise AttributeError(
            f"\nAttribute '{shoreline_key}' not found on Barrier3D object.\n"
            f"  Available attributes (sample): "
            f"{[a for a in available if 'x_s' in a.lower() or 'shore' in a.lower()]}\n"
            f"  → Update CASCADE_SHORELINE_KEY at the top of the script."
        )

    # Build time axis (1 timestep = 1 year in CASCADE)
    sample_ts = np.array(getattr(barrier3d[n_pad], shoreline_key))
    n_steps   = len(sample_ts)
    time      = np.linspace(start_year, end_year, n_steps)
    print(f"  '{shoreline_key}' length: {n_steps} timesteps  "
          f"({start_year}–{end_year})")

    # Extract LRR per real domain via linear regression
    lrr = {}
    for i in range(n_domains):
        b3d_idx = i + n_pad                          # index into full barrier3d list
        domain  = i + 1                              # 1-indexed CASCADE domain number
        y = np.array(getattr(barrier3d[b3d_idx], shoreline_key), dtype=float)
        valid = ~np.isnan(y)
        if valid.sum() > 2:
            slope, _ = np.polyfit(time[valid], y[valid], 1)
            lrr[domain] = slope   # units match x_s_TS units per year
        else:
            lrr[domain] = np.nan

    # Sanity check — print a few values
    sample_vals = {d: round(lrr[d], 3) for d in [1, 10, 30, 50, 70, 90]
                   if d in lrr and not np.isnan(lrr[d])}
    print(f"  Sample LRR values (check units — should be m/yr): {sample_vals}")
    return lrr


def compute_residual(coastsat_lrr, cascade_lrr, domains,
                      cascade_units=CASCADE_X_S_UNITS,
                      cascade_sign=CASCADE_X_S_SIGN):
    """
    Observed residual = CoastSat LRR - CASCADE baseline LRR, expressed in dam/yr.

    Sign convention:
      CoastSat LRR:  negative = erosion  (standard convention, from your CSV in m/yr)
      CASCADE x_s_TS: positive = landward retreat = EROSION (Barrier3D convention)
      CASCADE_X_S_SIGN = -1 flips CASCADE to standard convention before subtraction.

    Unit conversion:
      If CASCADE_X_UNITS = 'dam', multiply cascade LRR by 10 to get m/yr.
      Both are then in m/yr, residual divided by 10 to express in dam/yr.

    Positive residual → CoastSat more accretional than CASCADE → need source (+)
    Negative residual → CoastSat more erosional than CASCADE  → need sink (-)
    """
    m_per_unit = 10.0 if cascade_units == 'dam' else 1.0
    residual = {}
    for d in domains:
        cs  = coastsat_lrr.get(d, np.nan)
        cas = cascade_lrr.get(d, np.nan)
        if not (np.isnan(cs) or np.isnan(cas)):
            # Convert CASCADE to standard m/yr (apply sign + unit conversion)
            cas_standard_m = cas * m_per_unit * cascade_sign
            residual[d] = (cs - cas_standard_m) / 10.0   # result in dam/yr
    return residual

# =============================================================================
# RMSE + COMBINED OBJECTIVE
# =============================================================================

def rmse(residual, correction, eval_domains, exclusions):
    """RMSE between observed residual and parametric correction (dam/yr)."""
    errors = [
        residual[d] - correction[d]
        for d in eval_domains
        if d not in exclusions
        and d in residual
        and not np.isnan(residual[d])
    ]
    return np.sqrt(np.mean(np.array(errors) ** 2)) if errors else np.inf


def combined_score(rmse_p1, rmse_p2, w1=P1_WEIGHT, w2=P2_WEIGHT):
    """Weighted combined RMSE across both periods."""
    return w1 * rmse_p1 + w2 * rmse_p2

# =============================================================================
# PARETO FRONT
# =============================================================================

def find_pareto_front(df, col1='rmse_p1', col2='rmse_p2'):
    """
    Find Pareto-optimal solutions: no other point is better on BOTH metrics.
    Returns boolean mask.
    """
    is_pareto = np.ones(len(df), dtype=bool)
    vals = df[[col1, col2]].values
    for i, v in enumerate(vals):
        if is_pareto[i]:
            dominated = np.all(vals <= v, axis=1) & np.any(vals < v, axis=1)
            is_pareto[dominated] = False
    return is_pareto

# =============================================================================
# MAIN SWEEP
# =============================================================================

def run_joint_sweep(res_p1, res_p2):
    """
    Sweep all parameter combinations, computing RMSE for both periods
    and a weighted combined score.

    Evaluated over ALL 90 domains with no exclusions.
    """
    all_domains = list(range(1, N_DOMAINS + 1))

    n_total = (len(WS_PEAK_VALUES) * len(ROD_AMPLITUDES)
               * len(ROD_WIDTHS)   * len(AVON_NORTH_VALUES))
    print(f"Running {n_total} parameter combinations over {len(EVAL_DOMAINS)} domains...")

    results = []
    for ws_peak, rod_amp, rod_width, avon_north in product(
            WS_PEAK_VALUES, ROD_AMPLITUDES, ROD_WIDTHS, AVON_NORTH_VALUES):

        corr = build_correction(all_domains, ws_peak, rod_amp, rod_width, avon_north)

        r1 = rmse(res_p1, corr, EVAL_DOMAINS, EXCLUSION_ZONES)
        r2 = rmse(res_p2, corr, EVAL_DOMAINS, EXCLUSION_ZONES)
        sc = combined_score(r1, r2)

        results.append({
            'ws_peak_dam_yr':    ws_peak,
            'rod_amplitude':     rod_amp,
            'rod_width_domains': rod_width,
            'avon_north_dam_yr': avon_north,
            'rmse_p1':           round(r1, 4),
            'rmse_p2':           round(r2, 4),
            'combined_score':    round(sc, 4),
        })

    df = (pd.DataFrame(results)
            .sort_values('combined_score')
            .reset_index(drop=True))

    # Flag Pareto-optimal solutions
    df['pareto_optimal'] = find_pareto_front(df)
    return df

# =============================================================================
# DOMAIN_BE_RATES DICT GENERATOR
# =============================================================================

def print_be_rates(params, base_be=0.0,
                   multiplier=BE_UNIT_MULTIPLIER,
                   edge_d1=EDGE_D1, edge_d90=EDGE_D90):
    """
    Print DOMAIN_BE_RATES dict ready to paste into CASCADE run script.

    Automatically converts sweep comparison (dam/yr) → CASCADE units (m/yr)
    by applying BE_UNIT_MULTIPLIER (default ×10).

    Values are grouped by physical zone with inline comments, and
    6 domains per line for compact readability.
    """
    all_domains = list(range(1, N_DOMAINS + 1))
    corr = build_correction(all_domains, **params)

    # Build value dict: base_be + parametric correction, converted to m/yr
    vals = {}
    for d in all_domains:
        raw = base_be + corr[d]
        vals[d] = round(raw * multiplier, 3)

    def fmt(d):
        v = vals[d]
        s = f"{v:6.3f}" if v != 0 else f"{'0.000':>6}"
        return f"{d}: {s}"

    def row(domains):
        return "    " + ",   ".join(fmt(d) for d in domains) + ","

    ws_peak_m  = round(params['ws_peak']    * multiplier, 2)
    rod_amp_m  = round(params['rod_amp']    * multiplier, 2)
    avon_n_m   = round(params['avon_north'] * multiplier, 2)

    print(f"\n# ── DOMAIN_BE_RATES (full-domain joint fit, converted to m/yr) ────────────")
    print(f"# WS peak = {ws_peak_m} m/yr  |  "
          f"Rod amp = {rod_amp_m} m/yr  |  "
          f"Rod width = {params['rod_width']} domains  |  "
          f"Avon north = {avon_n_m} m/yr")
    print(f"# P1 weight = {P1_WEIGHT}  |  P2 weight = {P2_WEIGHT}  |  "
          f"Eval domains = 1–{N_DOMAINS} (no exclusions)  |  "
          f"Converted from dam/yr × {multiplier}")
    print()
    print("DOMAIN_BE_RATES = {")

    # ── Avon zone ────────────────────────────────────────────────────────────
    avon_doms = list(range(AVON_D_SOUTH, AVON_D_NORTH + 1))
    avon_nonzero = [d for d in avon_doms if vals[d] != 0.0]
    if avon_nonzero:
        print(f"    # Avon zone — linear ramp (offshore shoals + pier)")
        for i in range(0, len(avon_doms), 6):
            chunk = avon_doms[i:i+6]
            print(row(chunk))
        print()

    # ── Wimble Shoals southern taper ─────────────────────────────────────────
    ws_south = list(range(WS_X_START, WS_X_PEAK + 1))
    print(f"    # Wimble Shoals southern taper — 0 at domain {WS_X_START}, "
          f"peaks at domain {WS_X_PEAK} ({ws_peak_m} m/yr)")
    for i in range(0, len(ws_south), 6):
        chunk = ws_south[i:i+6]
        print(row(chunk))
    print()

    # ── Wimble Shoals northern taper ─────────────────────────────────────────
    ws_north = list(range(WS_X_PEAK + 1, WS_X_END + 1))
    print(f"    # Wimble Shoals northern taper — returns to 0 at domain {WS_X_END}")
    for i in range(0, len(ws_north), 6):
        chunk = ws_north[i:i+6]
        print(row(chunk))
    print()

    # ── Rodanthe bell ─────────────────────────────────────────────────────────
    rod_half  = int(round(params['rod_width'] * 2))
    rod_start = max(1, ROD_CENTER - rod_half)
    rod_end   = min(N_DOMAINS, ROD_CENTER + rod_half)
    rod_doms  = list(range(rod_start, rod_end + 1))
    rod_nonzero = [d for d in rod_doms if vals[d] != 0.0]
    if rod_nonzero:
        print(f"    # Rodanthe bell — Gaussian erosion centered on pier at domain {ROD_CENTER}")
        print(f"    # (values ADD to Wimble Shoals taper where zones overlap)")
        for i in range(0, len(rod_doms), 6):
            chunk = rod_doms[i:i+6]
            print(row(chunk))
        print()

    # ── Edge boundary conditions ──────────────────────────────────────────────
    print(f"    # Edge boundary conditions — period/Hs specific, set manually")
    print(f"    # 1984-2004: d1=-65(Hs2.5)/-130(Hs5.0)  d90=12(Hs2.5)/40(Hs5.0)")
    print(f"    # 2004-2024: d1=6(Hs2.5)/-40(Hs5.0)     d90=40(Hs2.5)/75(Hs5.0)")
    print(f"    1:  {edge_d1},")
    print(f"    90: {edge_d90},")
    print("}")
    print()

# =============================================================================
# PLOTTING
# =============================================================================

def plot_joint_results(best_params, res_p1, res_p2,
                        cs_p1, cs_p2, cas_p1, cas_p2, df):
    """
    4-row figure:
      Row 0: Raw CoastSat vs CASCADE LRR for both periods (diagnostic — data before correction)
      Row 1: Parametric correction shape (stacked components + total)
      Row 2: Period residual vs best-fit correction (P1 left, P2 right)
      Row 3: Pareto front (left) + WS peak sensitivity (right)
    """
    all_d   = list(range(1, N_DOMAINS + 1))
    corr    = build_correction(all_d, **best_params)
    corr_v  = [corr[d] for d in all_d]

    # Component curves
    ws_v  = [wimble_shadow(d, WS_X_START, WS_X_PEAK, WS_X_END,
                           best_params['ws_peak']) for d in all_d]
    rod_v = [rodanthe_bell(d, ROD_CENTER, best_params['rod_width'],
                           best_params['rod_amp']) for d in all_d]
    av_v  = [avon_correction(d, AVON_D_SOUTH, AVON_D_NORTH,
                             AVON_VAL_SOUTH, best_params['avon_north'])
             for d in all_d]

    fig = plt.figure(figsize=(18, 18))
    gs  = gridspec.GridSpec(4, 2, figure=fig, hspace=0.45, wspace=0.28)

    colors = {'p1': '#1f77b4', 'p2': '#ff7f0e', 'corr': '#2ca02c',
              'cs': '#9467bd', 'cas': '#8c564b'}

    # ── Row 0: Raw data diagnostic — CoastSat vs CASCADE LRR ─────────────────
    # Convert CASCADE LRR to standard m/yr for direct comparison with CoastSat
    m_per_unit = 10.0 if CASCADE_X_S_UNITS == 'dam' else 1.0

    ax_raw1 = fig.add_subplot(gs[0, 0])
    cs1_domains = sorted([d for d in all_d if d in cs_p1 and not np.isnan(cs_p1[d])])
    cas1_domains = sorted([d for d in all_d if d in cas_p1 and not np.isnan(cas_p1[d])])
    ax_raw1.plot(cs1_domains,
                 [cs_p1[d] for d in cs1_domains],
                 color=colors['cs'], lw=2, label='CoastSat LRR (m/yr)')
    ax_raw1.plot(cas1_domains,
                 [cas_p1[d] * m_per_unit * CASCADE_X_S_SIGN for d in cas1_domains],
                 color=colors['cas'], lw=2, ls='--', label='CASCADE baseline (m/yr)')
    ax_raw1.axhline(0, color='k', lw=0.8, ls='--')
    ax_raw1.axvspan(WS_X_START, WS_X_END, alpha=0.07, color='C0', label='Wimble Shoals zone')
    ax_raw1.set_ylabel('LRR (m/yr)', fontsize=10)
    ax_raw1.set_title(f'Raw data — {P1_LABEL}', fontsize=10)
    ax_raw1.legend(fontsize=8)
    ax_raw1.set_xlim(1, 90)
    ax_raw1.grid(True, alpha=0.3)

    ax_raw2 = fig.add_subplot(gs[0, 1])
    cs2_domains = sorted([d for d in all_d if d in cs_p2 and not np.isnan(cs_p2[d])])
    cas2_domains = sorted([d for d in all_d if d in cas_p2 and not np.isnan(cas_p2[d])])
    ax_raw2.plot(cs2_domains,
                 [cs_p2[d] for d in cs2_domains],
                 color=colors['cs'], lw=2, label='CoastSat LRR (m/yr)')
    ax_raw2.plot(cas2_domains,
                 [cas_p2[d] * m_per_unit * CASCADE_X_S_SIGN for d in cas2_domains],
                 color=colors['cas'], lw=2, ls='--', label='CASCADE baseline (m/yr)')
    ax_raw2.axhline(0, color='k', lw=0.8, ls='--')
    ax_raw2.axvspan(WS_X_START, WS_X_END, alpha=0.07, color='C0', label='Wimble Shoals zone')
    ax_raw2.set_ylabel('LRR (m/yr)', fontsize=10)
    ax_raw2.set_title(f'Raw data — {P2_LABEL}', fontsize=10)
    ax_raw2.legend(fontsize=8)
    ax_raw2.set_xlim(1, 90)
    ax_raw2.grid(True, alpha=0.3)

    # ── Row 1: Parametric correction shape ────────────────────────────────────
    ax1 = fig.add_subplot(gs[1, :])
    ax1.stackplot(all_d, ws_v, rod_v, av_v,
                  labels=['Wimble Shoals shadow (triangular)',
                          'Rodanthe bell (Gaussian)',
                          'Avon zone (linear)'],
                  alpha=0.5, colors=['C0', 'C3', 'C2'])
    ax1.plot(all_d, corr_v, 'k', lw=2.5, label='Total correction')
    ax1.axhline(0, color='k', lw=0.8, ls='--')
    ax1.axvspan(WS_X_START, WS_X_END, alpha=0.07, color='C0')
    ax1.set_ylabel('Correction (dam/yr)', fontsize=11)
    ax1.set_title(
        f'Joint best-fit parametric correction  |  '
        f'WS peak = {best_params["ws_peak"]} dam/yr  |  '
        f'P1 RMSE = {df.iloc[0]["rmse_p1"]:.3f}  |  '
        f'P2 RMSE = {df.iloc[0]["rmse_p2"]:.3f} dam/yr  |  '
        f'Full domain (1–90, no exclusions)',
        fontsize=10)
    ax1.legend(loc='upper left', fontsize=9)
    ax1.set_xlim(1, 90)
    ax1.grid(True, alpha=0.3)

    # ── Row 2 left: Period 1 residual vs correction ───────────────────────────
    ax2 = fig.add_subplot(gs[2, 0])
    rd1 = sorted(res_p1.keys())
    ax2.plot(rd1, [res_p1[d] for d in rd1],
             color=colors['p1'], lw=2, label=f'Residual {P1_LABEL}')
    ax2.plot(all_d, corr_v, color=colors['corr'], lw=2, ls='--',
             label='Parametric correction')
    ax2.axhline(0, color='k', lw=0.8, ls='--')
    ax2.axvspan(WS_X_START, WS_X_END, alpha=0.07, color='C0')
    ax2.set_ylabel('(dam/yr)', fontsize=10)
    ax2.set_title(f'Residual vs correction — {P1_LABEL}', fontsize=10)
    ax2.legend(fontsize=8)
    ax2.set_xlim(1, 90)
    ax2.grid(True, alpha=0.3)

    # ── Row 2 right: Period 2 residual vs correction ──────────────────────────
    ax3 = fig.add_subplot(gs[2, 1])
    rd2 = sorted(res_p2.keys())
    ax3.plot(rd2, [res_p2[d] for d in rd2],
             color=colors['p2'], lw=2, label=f'Residual {P2_LABEL}')
    ax3.plot(all_d, corr_v, color=colors['corr'], lw=2, ls='--',
             label='Parametric correction')
    ax3.axhline(0, color='k', lw=0.8, ls='--')
    ax3.axvspan(WS_X_START, WS_X_END, alpha=0.07, color='C0')
    ax3.set_ylabel('(dam/yr)', fontsize=10)
    ax3.set_title(f'Residual vs correction — {P2_LABEL}', fontsize=10)
    ax3.legend(fontsize=8)
    ax3.set_xlim(1, 90)
    ax3.grid(True, alpha=0.3)

    # ── Row 3 left: Pareto front ──────────────────────────────────────────────
    ax4 = fig.add_subplot(gs[3, 0])
    pareto     = df[df['pareto_optimal']].sort_values('rmse_p1')
    non_pareto = df[~df['pareto_optimal']]
    ax4.scatter(non_pareto['rmse_p1'], non_pareto['rmse_p2'],
                c='lightgrey', s=15, zorder=1, label='All combinations')
    ax4.scatter(pareto['rmse_p1'], pareto['rmse_p2'],
                c='C3', s=40, zorder=3, label='Pareto-optimal')
    best = df.iloc[0]
    ax4.scatter(best['rmse_p1'], best['rmse_p2'],
                c='k', s=100, marker='*', zorder=4,
                label=f'Best combined\n(P1={best["rmse_p1"]:.3f}, P2={best["rmse_p2"]:.3f})')
    ax4.set_xlabel(f'RMSE {P1_LABEL} (dam/yr)', fontsize=10)
    ax4.set_ylabel(f'RMSE {P2_LABEL} (dam/yr)', fontsize=10)
    ax4.set_title('Pareto front — period 1 vs period 2 trade-off', fontsize=10)
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    # ── Row 3 right: WS peak sensitivity ─────────────────────────────────────
    ax5 = fig.add_subplot(gs[3, 1])
    sens = df.groupby('ws_peak_dam_yr').agg(
        rmse_p1=('rmse_p1', 'min'),
        rmse_p2=('rmse_p2', 'min'),
        combined=('combined_score', 'min')
    )
    ax5.plot(sens.index, sens['rmse_p1'], 'C0o-', lw=2, label=P1_LABEL)
    ax5.plot(sens.index, sens['rmse_p2'], 'C1s-', lw=2, label=P2_LABEL)
    ax5.plot(sens.index, sens['combined'], 'ks--', lw=1.5, label='Combined (weighted)')
    ax5.axvline(best['ws_peak_dam_yr'], color='C3', ls='--', lw=1.5,
                label=f'Best = {best["ws_peak_dam_yr"]} dam/yr')
    ax5.set_xlabel('WS peak value (dam/yr)', fontsize=10)
    ax5.set_ylabel('Min RMSE (dam/yr)', fontsize=10)
    ax5.set_title('RMSE sensitivity to Wimble Shoals peak', fontsize=10)
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    fig.suptitle(
        f'Joint two-period parametric sweep — Hatteras Island\n'
        f'Full-domain evaluation (domains 1–90, no exclusions)  |  '
        f'P1 weight={P1_WEIGHT}  P2 weight={P2_WEIGHT}',
        fontsize=12, y=0.99)

    return fig

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load data for both periods ────────────────────────────────────────────
    print("Loading Period 1 data...")
    cs_p1  = load_coastsat(P1_COASTSAT_CSV)
    cas_p1 = extract_cascade_lrr(P1_CASCADE_BASELINE, P1_START, P1_END)

    print("Loading Period 2 data...")
    cs_p2  = load_coastsat(P2_COASTSAT_CSV)
    cas_p2 = extract_cascade_lrr(P2_CASCADE_BASELINE, P2_START, P2_END)

    # ── Compute residuals ─────────────────────────────────────────────────────
    all_domains = list(range(1, N_DOMAINS + 1))
    res_p1 = compute_residual(cs_p1, cas_p1, all_domains)
    res_p2 = compute_residual(cs_p2, cas_p2, all_domains)

    # ── Residual diagnostic — verify sign and units before sweeping ───────────
    print(f"\nResidual diagnostic (CASCADE_X_S_SIGN={CASCADE_X_S_SIGN}, "
          f"CASCADE_X_S_UNITS='{CASCADE_X_S_UNITS}'):")
    print(f"  Positive residual = CoastSat more accretional than CASCADE → correction should be POSITIVE (source)")
    print(f"  Negative residual = CoastSat more erosional than CASCADE  → correction should be NEGATIVE (sink)")
    print(f"  Evaluation: all 90 domains, no exclusions")
    for label, dom_range in [
        ("South island (1-20)",    list(range(1, 21))),
        ("Avon zone   (25-30)",    list(range(25, 31))),
        ("Mid-island  (31-37)",    list(range(31, 38))),
        ("Wimble Shoals (38-80)",  list(range(38, 81))),
        ("Rodanthe    (74-84)",    list(range(74, 85))),
        ("North island (81-90)",   list(range(81, 91))),
    ]:
        for period, res in [("P1", res_p1), ("P2", res_p2)]:
            vals = [res[d] for d in dom_range if d in res and not np.isnan(res[d])]
            if vals:
                print(f"  {label} {period}: "
                      f"mean={np.mean(vals):+.3f}  "
                      f"range=[{min(vals):+.3f}, {max(vals):+.3f}] dam/yr")
    print(f"\n  If Wimble Shoals residuals are negative but you expect a source,")
    print(f"  flip CASCADE_X_S_SIGN to +1 and rerun.")

    # ── Run sweep ─────────────────────────────────────────────────────────────
    df = run_joint_sweep(res_p1, res_p2)

    # Save results
    df.to_csv(os.path.join(OUTPUT_DIR, "sweep_results_joint_fulldomain.csv"), index=False)
    df[df['pareto_optimal']].to_csv(
        os.path.join(OUTPUT_DIR, "pareto_front_fulldomain.csv"), index=False)

    print(f"\nTop 10 combined-score parameter combinations:")
    print(df.head(10)[['ws_peak_dam_yr', 'rod_amplitude', 'rod_width_domains',
                        'avon_north_dam_yr', 'rmse_p1', 'rmse_p2',
                        'combined_score', 'pareto_optimal']].to_string(index=False))

    print(f"\nPareto-optimal solutions ({df['pareto_optimal'].sum()} found):")
    print(df[df['pareto_optimal']].sort_values('rmse_p1')[
        ['ws_peak_dam_yr', 'rod_amplitude', 'rod_width_domains',
         'avon_north_dam_yr', 'rmse_p1', 'rmse_p2']].to_string(index=False))

    # ── Best parameters ───────────────────────────────────────────────────────
    best_row = df.iloc[0]
    best_params = {
        'ws_peak':    best_row['ws_peak_dam_yr'],
        'rod_amp':    best_row['rod_amplitude'],
        'rod_width':  best_row['rod_width_domains'],
        'avon_north': best_row['avon_north_dam_yr'],
    }
    print(f"\nBest joint parameters (full-domain evaluation): {best_params}")
    print(f"  Period 1 RMSE: {best_row['rmse_p1']:.4f} dam/yr "
          f"({best_row['rmse_p1']*10:.2f} m/yr)")
    print(f"  Period 2 RMSE: {best_row['rmse_p2']:.4f} dam/yr "
          f"({best_row['rmse_p2']*10:.2f} m/yr)")
    print(f"  Combined score: {best_row['combined_score']:.4f} dam/yr")

    # ── Generate DOMAIN_BE_RATES ──────────────────────────────────────────────
    print_be_rates(best_params, base_be=0.0)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig = plot_joint_results(best_params, res_p1, res_p2,
                              cs_p1, cs_p2, cas_p1, cas_p2, df)
    out_plot = os.path.join(OUTPUT_DIR, "joint_best_fit_fulldomain.png")
    fig.savefig(out_plot, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {out_plot}")
    plt.show()
