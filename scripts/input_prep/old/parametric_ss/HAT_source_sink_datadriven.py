"""
HAT_source_sink_datadriven.py
==============================
Data-driven joint two-period parametric sweep for source/sink correction.

Evaluates over ALL 90 domains with NO exclusions.
Correction shapes were selected from the residual pattern observed in
HAT_source_sink_bothperiods_fulldomain.py:

  Shape 1 — South island linear ramp   (domains  1–20)
    Physical basis: Oregon Inlet boundary influence + sparse early Landsat
    in the south; large negative P1 residual, mostly absent in P2.
    ⚠ Flag: P1/P2 asymmetry means this may reflect data quality, not
    persistent morphodynamics. Inspect the printed per-period RMSE
    contribution carefully before applying to CASCADE runs.

  Shape 2 — Mid-island linear fill      (domains 31–37)
    Physical basis: gap between Avon Pier groin effect and Wimble Shoals
    shadow onset; large positive P2 residual consistent with post-Isabel
    morphodynamic recovery in this zone.

  Shape 3 — Wimble Shoals triangle      (domains 38–80)
    Physical basis: diffraction shadow from fixed bathymetric feature.
    Boundaries fixed to bathymetric evidence; peak magnitude is free.

  Shape 4 — Rodanthe Gaussian bell      (center domain 79)
    Physical basis: shoreline curvature → divergent transport at pier.
    Center fixed to structure location; amplitude and width are free.

  Shape 5 — Avon linear ramp            (domains 25–30)
    Physical basis: offshore shoal + Avon Pier groin effect.
    South anchor fixed; north value is free.

Optimization:
  - Combined RMSE = 0.5 × P1_RMSE + 0.5 × P2_RMSE (equal weighting)
  - Graded over all 90 domains, no exclusions
  - Pareto front reported for P1 vs P2 trade-off

Outputs:
  - sweep_results_datadriven.csv     : all combinations with P1, P2, combined RMSE
  - pareto_front_datadriven.csv      : Pareto-optimal solutions
  - datadriven_best_fit.png          : 5-row diagnostic figure
  - DOMAIN_BE_RATES dict             : printed, ready to paste into CASCADE run script

Author: Hannah Henry, UNC Chapel Hill
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from itertools import product
import os

# =============================================================================
# CONFIGURATION
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
OUTPUT_DIR = r"/scripts/input_prep/original/parametric_ss/results_datadriven_Hs2.5"

# CASCADE npz comparison variable name for shoreline position
CASCADE_SHORELINE_KEY = 'x_s_TS'

# CASCADE sign convention:
#   Positive x_s_TS slope = landward retreat = EROSION (Barrier3D convention)
#   CASCADE_X_S_SIGN = -1 flips to standard (erosion negative) before residual
CASCADE_X_S_SIGN  = -1
CASCADE_X_S_UNITS = 'dam'   # 'dam' → multiply by 10 for m/yr

# CoastSat CSV column names
CS_DOMAIN_COL = 'domain_number'
CS_LRR_COL    = 'mean_lrr'

# CASCADE domain setup
N_DOMAINS = 90
N_PAD     = 15

# Equal weighting between periods
P1_WEIGHT = 0.5
P2_WEIGHT = 0.5

# Evaluate over full domain, no exclusions
EVAL_DOMAINS    = list(range(1, N_DOMAINS + 1))
EXCLUSION_ZONES = set()

# Unit multiplier: dam/yr → m/yr for DOMAIN_BE_RATES comparison
BE_UNIT_MULTIPLIER = 10.0

# Edge boundary conditions — period/Hs specific, set manually
EDGE_D1  =   6
EDGE_D90 =  40

# =============================================================================
# FIXED SHAPE BOUNDARIES — from bathymetric/structural evidence
# =============================================================================

# Shape 1: South island ramp
SI_X_START =  1    # Oregon Inlet / domain 1
SI_X_END   = 20    # Southern taper end

# Shape 2: Mid-island fill
MI_X_START = 31    # Just north of Avon zone
MI_X_END   = 37    # Just south of Wimble Shoals onset

# Shape 3: Wimble Shoals triangle
WS_X_START = 38
WS_X_PEAK  = 68
WS_X_END   = 80

# Shape 4: Rodanthe Gaussian
ROD_CENTER = 79

# Shape 5: Avon ramp
AVON_D_SOUTH   = 25
AVON_D_NORTH   = 30
AVON_VAL_SOUTH = -0.10   # Fixed south anchor (dam/yr)

# =============================================================================
# SWEEP RANGES
# All values in dam/yr. Ranges informed by observed residual magnitudes.
#
# South island: P1 residual mean ~ -0.1 to -0.3 dam/yr → sweep negative
# Mid-island:   P2 residual mean ~ +0.3 dam/yr          → sweep positive
# Wimble Shoals: small positive signal, as before        → narrow positive range
# Rodanthe:     small negative signal                    → narrow negative range
# Avon north:   positive, as before
#
# Coarse sweep first — narrow around best-fit for a finer pass.
# =============================================================================

SI_AMP_VALUES     = np.round(np.arange(-0.50,  0.05, 0.05), 3)  # South island amplitude (dam/yr)
MI_AMP_VALUES     = np.round(np.arange( 0.05,  0.55, 0.05), 3)  # Mid-island fill amplitude (dam/yr)
WS_PEAK_VALUES    = np.round(np.arange( 0.05,  0.35, 0.05), 3)  # Wimble Shoals peak (dam/yr)
ROD_AMPLITUDES    = np.round(np.arange(-0.30, -0.05, 0.05), 3)  # Rodanthe amplitude (dam/yr, negative)
ROD_WIDTHS        = np.arange(3, 7, 1)                           # Rodanthe half-width (domains)
AVON_NORTH_VALUES = np.round(np.arange( 0.10,  0.55, 0.05), 3)  # Avon north (dam/yr)

# =============================================================================
# PARAMETRIC CORRECTION FUNCTIONS
# =============================================================================

def south_island_ramp(domain, x_start, x_end, amplitude):
    """
    Linear ramp across the south island (dam/yr).
    Peak at x_start (domain 1, Oregon Inlet end), tapers to 0 at x_end.
    Physical basis: Oregon Inlet boundary influence attenuates northward.
    ⚠ P1/P2 asymmetry — inspect per-period contributions before applying.
    """
    if domain < x_start or domain > x_end:
        return 0.0
    # Ramp: full amplitude at domain 1, zero at domain 20
    t = (domain - x_start) / (x_end - x_start)
    return amplitude * (1.0 - t)


def mid_island_fill(domain, x_start, x_end, amplitude):
    """
    Flat fill across the mid-island gap (dam/yr).
    Constant value across domains 31–37.
    Physical basis: post-Isabel recovery zone between Avon Pier
    groin effect and Wimble Shoals shadow onset.
    """
    if domain < x_start or domain > x_end:
        return 0.0
    return amplitude


def wimble_shadow(domain, x_start, x_peak, x_end, peak_value):
    """
    Triangular wave shadow (dam/yr).
    Linear ramp: 0 at x_start → peak at x_peak → 0 at x_end.
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
    """
    return amplitude * np.exp(-0.5 * ((domain - center) / width) ** 2)


def avon_correction(domain, d_south, d_north, val_south, val_north):
    """
    Linear ramp for Avon zone (dam/yr).
    """
    if domain < d_south or domain > d_north:
        return 0.0
    t = (domain - d_south) / (d_north - d_south)
    return val_south + t * (val_north - val_south)


def build_correction(domains, si_amp, mi_amp, ws_peak, rod_amp, rod_width, avon_north):
    """Assemble total correction from all five shapes."""
    return {
        d: (south_island_ramp(d,  SI_X_START, SI_X_END,  si_amp)
          + mid_island_fill(d,    MI_X_START, MI_X_END,  mi_amp)
          + wimble_shadow(d,      WS_X_START, WS_X_PEAK, WS_X_END, ws_peak)
          + rodanthe_bell(d,      ROD_CENTER, rod_width, rod_amp)
          + avon_correction(d,    AVON_D_SOUTH, AVON_D_NORTH, AVON_VAL_SOUTH, avon_north))
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
    """Extract per-domain LRR from CASCADE comparison npz."""
    data = np.load(npz_path, allow_pickle=True)
    print(f"  NPZ keys found: {list(data.keys())}")
    cascade_obj = data['cascade'][0]
    barrier3d   = cascade_obj.barrier3d
    print(f"  Total barrier3d domains (inc. padding): {len(barrier3d)}")

    test_b3d = barrier3d[n_pad]
    if not hasattr(test_b3d, shoreline_key):
        available = [a for a in dir(test_b3d) if 'x_s' in a.lower() or 'shore' in a.lower()]
        raise AttributeError(
            f"\nAttribute '{shoreline_key}' not found. "
            f"Shoreline-related attributes: {available}\n"
            f"  → Update CASCADE_SHORELINE_KEY."
        )

    sample_ts = np.array(getattr(barrier3d[n_pad], shoreline_key))
    n_steps   = len(sample_ts)
    time      = np.linspace(start_year, end_year, n_steps)
    print(f"  '{shoreline_key}' length: {n_steps} timesteps ({start_year}–{end_year})")

    lrr = {}
    for i in range(n_domains):
        b3d_idx = i + n_pad
        domain  = i + 1
        y = np.array(getattr(barrier3d[b3d_idx], shoreline_key), dtype=float)
        valid = ~np.isnan(y)
        if valid.sum() > 2:
            slope, _ = np.polyfit(time[valid], y[valid], 1)
            lrr[domain] = slope
        else:
            lrr[domain] = np.nan

    sample_vals = {d: round(lrr[d], 3) for d in [1, 10, 30, 50, 70, 90]
                   if d in lrr and not np.isnan(lrr[d])}
    print(f"  Sample LRR (dam/yr raw): {sample_vals}")
    return lrr


def compute_residual(coastsat_lrr, cascade_lrr, domains,
                      cascade_units=CASCADE_X_S_UNITS,
                      cascade_sign=CASCADE_X_S_SIGN):
    """
    Residual = CoastSat LRR - CASCADE baseline LRR, in dam/yr.
    Positive → CoastSat more accretional → need source (+)
    Negative → CoastSat more erosional   → need sink (-)
    """
    m_per_unit = 10.0 if cascade_units == 'dam' else 1.0
    residual = {}
    for d in domains:
        cs  = coastsat_lrr.get(d, np.nan)
        cas = cascade_lrr.get(d, np.nan)
        if not (np.isnan(cs) or np.isnan(cas)):
            cas_m      = cas * m_per_unit * cascade_sign
            residual[d] = (cs - cas_m) / 10.0   # dam/yr
    return residual

# =============================================================================
# RMSE + COMBINED OBJECTIVE
# =============================================================================

def rmse(residual, correction, eval_domains, exclusions):
    errors = [
        residual[d] - correction[d]
        for d in eval_domains
        if d not in exclusions
        and d in residual
        and not np.isnan(residual[d])
    ]
    return np.sqrt(np.mean(np.array(errors) ** 2)) if errors else np.inf


def combined_score(r1, r2, w1=P1_WEIGHT, w2=P2_WEIGHT):
    return w1 * r1 + w2 * r2

# =============================================================================
# PARETO FRONT
# =============================================================================

def find_pareto_front(df, col1='rmse_p1', col2='rmse_p2'):
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

def run_sweep(res_p1, res_p2):
    all_domains = list(range(1, N_DOMAINS + 1))

    n_total = (len(SI_AMP_VALUES) * len(MI_AMP_VALUES) * len(WS_PEAK_VALUES)
               * len(ROD_AMPLITUDES) * len(ROD_WIDTHS) * len(AVON_NORTH_VALUES))
    print(f"\nRunning {n_total:,} parameter combinations over {len(EVAL_DOMAINS)} domains...")
    print(f"(This may take a minute — 6 free parameters)\n")

    results = []
    for si_amp, mi_amp, ws_peak, rod_amp, rod_width, avon_north in product(
            SI_AMP_VALUES, MI_AMP_VALUES, WS_PEAK_VALUES,
            ROD_AMPLITUDES, ROD_WIDTHS, AVON_NORTH_VALUES):

        corr = build_correction(all_domains, si_amp, mi_amp, ws_peak,
                                rod_amp, rod_width, avon_north)
        r1 = rmse(res_p1, corr, EVAL_DOMAINS, EXCLUSION_ZONES)
        r2 = rmse(res_p2, corr, EVAL_DOMAINS, EXCLUSION_ZONES)
        sc = combined_score(r1, r2)

        results.append({
            'si_amp_dam_yr':     si_amp,
            'mi_amp_dam_yr':     mi_amp,
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
    df['pareto_optimal'] = find_pareto_front(df)
    return df

# =============================================================================
# DOMAIN_BE_RATES DICT GENERATOR
# =============================================================================

def print_be_rates(params, base_be=0.0,
                   multiplier=BE_UNIT_MULTIPLIER,
                   edge_d1=EDGE_D1, edge_d90=EDGE_D90):
    """Print DOMAIN_BE_RATES dict ready to paste into CASCADE run script."""
    all_domains = list(range(1, N_DOMAINS + 1))
    corr = build_correction(all_domains, **params)

    vals = {d: round((base_be + corr[d]) * multiplier, 3) for d in all_domains}

    def fmt(d):
        v = vals[d]
        return f"{d}: {v:6.3f}"

    def row(doms):
        return "    " + ",   ".join(fmt(d) for d in doms) + ","

    print(f"\n# ── DOMAIN_BE_RATES (data-driven joint fit, converted to m/yr) ─────────────")
    print(f"# SI amp    = {round(params['si_amp']*multiplier,2)} m/yr  |  "
          f"MI amp = {round(params['mi_amp']*multiplier,2)} m/yr  |  "
          f"WS peak = {round(params['ws_peak']*multiplier,2)} m/yr")
    print(f"# Rod amp   = {round(params['rod_amp']*multiplier,2)} m/yr  |  "
          f"Rod width = {params['rod_width']} domains  |  "
          f"Avon north = {round(params['avon_north']*multiplier,2)} m/yr")
    print(f"# P1 weight = {P1_WEIGHT}  |  P2 weight = {P2_WEIGHT}  |  "
          f"Eval: all 90 domains  |  dam/yr × {multiplier}")
    print()
    print("DOMAIN_BE_RATES = {")

    # South island ramp
    si_doms = list(range(SI_X_START, SI_X_END + 1))
    if any(vals[d] != 0.0 for d in si_doms):
        print(f"    # South island ramp — peaks at domain 1 (Oregon Inlet end), "
              f"tapers to 0 at domain {SI_X_END}")
        print(f"    # ⚠ P1/P2 asymmetric signal — flag for advisor review")
        for i in range(0, len(si_doms), 6):
            print(row(si_doms[i:i+6]))
        print()

    # Mid-island fill
    mi_doms = list(range(MI_X_START, MI_X_END + 1))
    if any(vals[d] != 0.0 for d in mi_doms):
        print(f"    # Mid-island fill — flat correction domains {MI_X_START}–{MI_X_END}")
        print(f"    # (post-Isabel recovery gap between Avon Pier and Wimble Shoals onset)")
        for i in range(0, len(mi_doms), 6):
            print(row(mi_doms[i:i+6]))
        print()

    # Avon zone
    avon_doms = list(range(AVON_D_SOUTH, AVON_D_NORTH + 1))
    if any(vals[d] != 0.0 for d in avon_doms):
        print(f"    # Avon zone — linear ramp (offshore shoals + pier)")
        for i in range(0, len(avon_doms), 6):
            print(row(avon_doms[i:i+6]))
        print()

    # Wimble Shoals
    ws_doms = list(range(WS_X_START, WS_X_END + 1))
    print(f"    # Wimble Shoals — triangle peak at domain {WS_X_PEAK} "
          f"({round(params['ws_peak']*multiplier,2)} m/yr)")
    for i in range(0, len(ws_doms), 6):
        print(row(ws_doms[i:i+6]))
    print()

    # Rodanthe bell (just the non-WS portion for clarity)
    rod_half  = int(round(params['rod_width'] * 2))
    rod_start = max(1, ROD_CENTER - rod_half)
    rod_end   = min(N_DOMAINS, ROD_CENTER + rod_half)
    rod_doms  = [d for d in range(rod_start, rod_end + 1) if d not in ws_doms]
    if rod_doms and any(vals[d] != 0.0 for d in rod_doms):
        print(f"    # Rodanthe bell — Gaussian beyond Wimble Shoals zone "
              f"(center domain {ROD_CENTER})")
        for i in range(0, len(rod_doms), 6):
            print(row(rod_doms[i:i+6]))
        print()

    print(f"    # Edge boundary conditions — period/Hs specific, set manually")
    print(f"    # 1984-2004: d1=-65(Hs2.5)  d90=12(Hs2.5)")
    print(f"    # 2004-2024: d1=6(Hs2.5)    d90=40(Hs2.5)")
    print(f"    1:  {edge_d1},")
    print(f"    90: {edge_d90},")
    print("}")

# =============================================================================
# PLOTTING
# =============================================================================

def plot_results(best_params, res_p1, res_p2,
                  cs_p1, cs_p2, cas_p1, cas_p2, df):
    """
    5-row figure:
      Row 0: Raw CoastSat vs CASCADE LRR (P1 left, P2 right)
      Row 1: Full correction shape (stacked components)
      Row 2: P1 residual vs correction | P2 residual vs correction
      Row 3: Corrected residual (residual − correction) for both periods
      Row 4: Pareto front | Sensitivity to south island amplitude
    """
    all_d  = list(range(1, N_DOMAINS + 1))
    corr   = build_correction(all_d, **best_params)
    corr_v = [corr[d] for d in all_d]

    # Component arrays
    si_v  = [south_island_ramp(d, SI_X_START, SI_X_END, best_params['si_amp'])
              for d in all_d]
    mi_v  = [mid_island_fill(d, MI_X_START, MI_X_END, best_params['mi_amp'])
              for d in all_d]
    ws_v  = [wimble_shadow(d, WS_X_START, WS_X_PEAK, WS_X_END, best_params['ws_peak'])
              for d in all_d]
    rod_v = [rodanthe_bell(d, ROD_CENTER, best_params['rod_width'], best_params['rod_amp'])
              for d in all_d]
    av_v  = [avon_correction(d, AVON_D_SOUTH, AVON_D_NORTH,
                             AVON_VAL_SOUTH, best_params['avon_north'])
              for d in all_d]

    m_per_unit = 10.0 if CASCADE_X_S_UNITS == 'dam' else 1.0
    colors = {'p1': '#1f77b4', 'p2': '#ff7f0e', 'corr': '#2ca02c',
              'cs': '#9467bd', 'cas': '#8c564b'}

    fig = plt.figure(figsize=(18, 22))
    gs  = gridspec.GridSpec(5, 2, figure=fig, hspace=0.48, wspace=0.28)

    # ── Row 0: Raw data ───────────────────────────────────────────────────────
    for col, (label, cs_lrr, cas_lrr) in enumerate([
            (P1_LABEL, cs_p1, cas_p1), (P2_LABEL, cs_p2, cas_p2)]):
        ax = fig.add_subplot(gs[0, col])
        cs_d  = sorted([d for d in all_d if d in cs_lrr  and not np.isnan(cs_lrr[d])])
        cas_d = sorted([d for d in all_d if d in cas_lrr and not np.isnan(cas_lrr[d])])
        ax.plot(cs_d,  [cs_lrr[d] for d in cs_d],
                color=colors['cs'],  lw=2, label='CoastSat LRR (m/yr)')
        ax.plot(cas_d, [cas_lrr[d] * m_per_unit * CASCADE_X_S_SIGN for d in cas_d],
                color=colors['cas'], lw=2, ls='--', label='CASCADE baseline (m/yr)')
        ax.axhline(0, color='k', lw=0.8, ls='--')
        for span, color, alpha in [
                ((SI_X_START, SI_X_END),   'C3',  0.10),
                ((MI_X_START, MI_X_END),   'C4',  0.10),
                ((WS_X_START, WS_X_END),   'C0',  0.07),
                ((AVON_D_SOUTH, AVON_D_NORTH), 'C2', 0.10)]:
            ax.axvspan(span[0], span[1], alpha=alpha, color=color)
        ax.set_ylabel('LRR (m/yr)', fontsize=10)
        ax.set_title(f'Raw data — {label}', fontsize=10)
        ax.legend(fontsize=8)
        ax.set_xlim(1, 90)
        ax.grid(True, alpha=0.3)

    # ── Row 1: Correction shape ────────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[1, :])
    ax1.stackplot(all_d, si_v, mi_v, ws_v, rod_v, av_v,
                  labels=['South island ramp', 'Mid-island fill',
                          'Wimble Shoals triangle', 'Rodanthe bell',
                          'Avon ramp'],
                  alpha=0.5, colors=['C3', 'C4', 'C0', 'C1', 'C2'])
    ax1.plot(all_d, corr_v, 'k', lw=2.5, label='Total correction')
    ax1.axhline(0, color='k', lw=0.8, ls='--')
    ax1.set_ylabel('Correction (dam/yr)', fontsize=11)
    ax1.set_title(
        f'Best-fit correction  |  '
        f'SI={best_params["si_amp"]} dam/yr  |  '
        f'MI={best_params["mi_amp"]} dam/yr  |  '
        f'WS={best_params["ws_peak"]} dam/yr  |  '
        f'Rod={best_params["rod_amp"]} dam/yr  |  '
        f'Avon N={best_params["avon_north"]} dam/yr',
        fontsize=10)
    ax1.legend(loc='upper right', fontsize=8, ncol=3)
    ax1.set_xlim(1, 90)
    ax1.grid(True, alpha=0.3)

    # ── Row 2: Residual vs correction ──────────────────────────────────────────
    for col, (label, res, color_key) in enumerate([
            (P1_LABEL, res_p1, 'p1'), (P2_LABEL, res_p2, 'p2')]):
        ax = fig.add_subplot(gs[2, col])
        rd = sorted(res.keys())
        ax.plot(rd, [res[d] for d in rd],
                color=colors[color_key], lw=2, label=f'Residual {label}')
        ax.plot(all_d, corr_v, color=colors['corr'], lw=2, ls='--',
                label='Parametric correction')
        ax.axhline(0, color='k', lw=0.8, ls='--')
        ax.set_ylabel('(dam/yr)', fontsize=10)
        ax.set_title(f'Residual vs correction — {label}', fontsize=10)
        ax.legend(fontsize=8)
        ax.set_xlim(1, 90)
        ax.grid(True, alpha=0.3)

    # ── Row 3: Corrected residual (what CASCADE still needs to explain) ────────
    ax_r1 = fig.add_subplot(gs[3, 0])
    ax_r2 = fig.add_subplot(gs[3, 1])
    for ax, res, label, color_key in [
            (ax_r1, res_p1, P1_LABEL, 'p1'),
            (ax_r2, res_p2, P2_LABEL, 'p2')]:
        rd = sorted(res.keys())
        remaining = [res[d] - corr[d] for d in rd if d in corr]
        ax.bar(rd, remaining, color=colors[color_key], alpha=0.6, width=0.8)
        ax.axhline(0, color='k', lw=1)
        rmse_val = np.sqrt(np.mean(np.array(remaining)**2))
        ax.set_ylabel('Remaining residual (dam/yr)', fontsize=10)
        ax.set_title(f'Corrected residual — {label}  |  RMSE={rmse_val:.3f} dam/yr',
                     fontsize=10)
        ax.set_xlim(1, 90)
        ax.grid(True, alpha=0.3)

    # ── Row 4: Pareto front + SI amplitude sensitivity ────────────────────────
    ax4 = fig.add_subplot(gs[4, 0])
    pareto     = df[df['pareto_optimal']].sort_values('rmse_p1')
    non_pareto = df[~df['pareto_optimal']]
    ax4.scatter(non_pareto['rmse_p1'], non_pareto['rmse_p2'],
                c='lightgrey', s=10, zorder=1, label='All combinations')
    ax4.scatter(pareto['rmse_p1'], pareto['rmse_p2'],
                c='C3', s=40, zorder=3, label='Pareto-optimal')
    best = df.iloc[0]
    ax4.scatter(best['rmse_p1'], best['rmse_p2'],
                c='k', s=120, marker='*', zorder=4,
                label=f'Best combined\n(P1={best["rmse_p1"]:.3f}, P2={best["rmse_p2"]:.3f})')
    ax4.set_xlabel(f'RMSE {P1_LABEL} (dam/yr)', fontsize=10)
    ax4.set_ylabel(f'RMSE {P2_LABEL} (dam/yr)', fontsize=10)
    ax4.set_title('Pareto front — P1 vs P2 trade-off', fontsize=10)
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    ax5 = fig.add_subplot(gs[4, 1])
    sens = df.groupby('si_amp_dam_yr').agg(
        rmse_p1=('rmse_p1', 'min'),
        rmse_p2=('rmse_p2', 'min'),
        combined=('combined_score', 'min')
    )
    ax5.plot(sens.index * 10, sens['rmse_p1'], 'C0o-', lw=2, label=P1_LABEL)
    ax5.plot(sens.index * 10, sens['rmse_p2'], 'C1s-', lw=2, label=P2_LABEL)
    ax5.plot(sens.index * 10, sens['combined'], 'ks--', lw=1.5, label='Combined')
    ax5.axvline(best['si_amp_dam_yr'] * 10, color='C3', ls='--', lw=1.5,
                label=f'Best = {round(best["si_amp_dam_yr"]*10,2)} m/yr')
    ax5.set_xlabel('South island amplitude (m/yr)', fontsize=10)
    ax5.set_ylabel('Min RMSE (dam/yr)', fontsize=10)
    ax5.set_title('RMSE sensitivity to south island correction', fontsize=10)
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    fig.suptitle(
        f'Data-driven joint parametric sweep — Hatteras Island\n'
        f'5-shape library  |  Full domain (1–90, no exclusions)  |  '
        f'P1={P1_WEIGHT}  P2={P2_WEIGHT}  equal weighting',
        fontsize=12, y=0.995)

    return fig

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # ── Load data ─────────────────────────────────────────────────────────────
    print("Loading Period 1 data...")
    cs_p1  = load_coastsat(P1_COASTSAT_CSV)
    cas_p1 = extract_cascade_lrr(P1_CASCADE_BASELINE, P1_START, P1_END)

    print("\nLoading Period 2 data...")
    cs_p2  = load_coastsat(P2_COASTSAT_CSV)
    cas_p2 = extract_cascade_lrr(P2_CASCADE_BASELINE, P2_START, P2_END)

    # ── Compute residuals ─────────────────────────────────────────────────────
    all_domains = list(range(1, N_DOMAINS + 1))
    res_p1 = compute_residual(cs_p1, cas_p1, all_domains)
    res_p2 = compute_residual(cs_p2, cas_p2, all_domains)

    # ── Residual diagnostic ───────────────────────────────────────────────────
    print(f"\nResidual diagnostic (sign={CASCADE_X_S_SIGN}, units='{CASCADE_X_S_UNITS}'):")
    print(f"  Positive residual → CoastSat more accretional → need source (+)")
    print(f"  Negative residual → CoastSat more erosional   → need sink (−)")
    zones = [
        ("South island  ( 1–20)", list(range(1,  21))),
        ("Mid-island    (21–30)", list(range(21, 31))),
        ("Mid-island    (31–37)", list(range(31, 38))),
        ("Wimble Shoals (38–80)", list(range(38, 81))),
        ("Rodanthe      (74–84)", list(range(74, 85))),
        ("North island  (81–90)", list(range(81, 91))),
    ]
    for label, doms in zones:
        for period, res in [("P1", res_p1), ("P2", res_p2)]:
            vals = [res[d] for d in doms if d in res and not np.isnan(res[d])]
            if vals:
                print(f"  {label} {period}: "
                      f"mean={np.mean(vals):+.3f}  "
                      f"range=[{min(vals):+.3f}, {max(vals):+.3f}] dam/yr")

    # ── South island asymmetry warning ────────────────────────────────────────
    si_p1 = np.mean([res_p1[d] for d in range(1, 21) if d in res_p1])
    si_p2 = np.mean([res_p2[d] for d in range(1, 21) if d in res_p2])
    print(f"\n⚠ South island residual asymmetry:")
    print(f"  P1 mean = {si_p1:+.3f} dam/yr  |  P2 mean = {si_p2:+.3f} dam/yr")
    print(f"  Large P1/P2 difference may reflect sparse early Landsat coverage")
    print(f"  rather than persistent morphodynamics. Review south island correction")
    print(f"  carefully before applying to CASCADE runs.")

    # ── Run sweep ─────────────────────────────────────────────────────────────
    df = run_sweep(res_p1, res_p2)

    df.to_csv(os.path.join(OUTPUT_DIR, "sweep_results_datadriven.csv"), index=False)
    df[df['pareto_optimal']].to_csv(
        os.path.join(OUTPUT_DIR, "pareto_front_datadriven.csv"), index=False)

    print(f"\nTop 10 combinations:")
    print(df.head(10)[['si_amp_dam_yr', 'mi_amp_dam_yr', 'ws_peak_dam_yr',
                        'rod_amplitude', 'rod_width_domains', 'avon_north_dam_yr',
                        'rmse_p1', 'rmse_p2', 'combined_score']].to_string(index=False))

    print(f"\nPareto-optimal solutions ({df['pareto_optimal'].sum()} found):")
    print(df[df['pareto_optimal']].sort_values('rmse_p1')[
        ['si_amp_dam_yr', 'mi_amp_dam_yr', 'ws_peak_dam_yr',
         'rod_amplitude', 'rod_width_domains', 'avon_north_dam_yr',
         'rmse_p1', 'rmse_p2']].to_string(index=False))

    # ── Best parameters ───────────────────────────────────────────────────────
    best_row = df.iloc[0]
    best_params = {
        'si_amp':    best_row['si_amp_dam_yr'],
        'mi_amp':    best_row['mi_amp_dam_yr'],
        'ws_peak':   best_row['ws_peak_dam_yr'],
        'rod_amp':   best_row['rod_amplitude'],
        'rod_width': best_row['rod_width_domains'],
        'avon_north':best_row['avon_north_dam_yr'],
    }
    print(f"\nBest joint parameters: {best_params}")
    print(f"  P1 RMSE: {best_row['rmse_p1']:.4f} dam/yr  "
          f"({best_row['rmse_p1']*10:.2f} m/yr)")
    print(f"  P2 RMSE: {best_row['rmse_p2']:.4f} dam/yr  "
          f"({best_row['rmse_p2']*10:.2f} m/yr)")
    print(f"  Combined: {best_row['combined_score']:.4f} dam/yr")

    # ── DOMAIN_BE_RATES ───────────────────────────────────────────────────────
    print_be_rates(best_params)

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig = plot_results(best_params, res_p1, res_p2,
                        cs_p1, cs_p2, cas_p1, cas_p2, df)
    out_plot = os.path.join(OUTPUT_DIR, "datadriven_best_fit.png")
    fig.savefig(out_plot, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to {out_plot}")
    plt.show()
