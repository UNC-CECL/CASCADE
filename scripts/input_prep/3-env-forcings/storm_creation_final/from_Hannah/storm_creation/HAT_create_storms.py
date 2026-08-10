# =============================================================================
# HAT_create_storm_file.py
# CASCADE-Formatted Storm File Generator — Hatteras Island
# -----------------------------------------------------------------------------
# Description:
#   Generates a CASCADE-compatible storm file (.npy) from NOAA tide gauge data
#   and WIS hindcast wave data. Storm events are identified as periods where
#   the computed Total Water Level (TWL) exceeds a threshold. For each event,
#   five parameters are extracted and stored in CASCADE's native unit system:
#
#       [Year_Index, Rhigh, Rlow, Wave_Period, Duration]
#
#   Physics (Stockdon et al., 2006):
#       Rhigh = eta_obs + R2%          (2% exceedance total runup, in dam)
#       Rlow  = still water level at storm peak, floored at MHW (in dam)
#       Wave Period = mean peak Tp during storm (s)
#       Duration = storm length in hours (h), capped at MAX_STORM_DURATION_HR
#       Year_Index = calendar year of storm peak minus START_YEAR
#
#   Units note: CASCADE uses decameters (dam) internally. All elevations are
#   divided by 10 before saving. Wave period and duration remain in their
#   native units (seconds and hours respectively).
#
# -----------------------------------------------------------------------------
# KEY PARAMETER DECISIONS (informed by two working Outer Banks reference files):
#
#   MAX_STORM_DURATION_HR = 36
#       CASCADE ran successfully with max durations of exactly
#       36 hours. Hannah's file crashed with durations up to 189 hours.
#       triggering a C-level access violation (0xC0000005).
#
#   MIN_INTER_STORM_GAP_HR = 48
#       Increasing to 48 hours keeps separate
#       storm events distinct, producing realistic individual durations.
#       MAX_STORM_DURATION_HR acts as a backstop for any events that still
#       exceed the ceiling after the gap change.
#
#   STORM_THRESHOLD = 1.7 m
#       This produces a realistic distribution of forcing intensities.
#
#   Rlow = max(water_level_at_peak_TWL, MHW_M) / 10
#       Rlow is the still water level (tide + surge, no runup or setup) at the
#       storm's peak intensity moment, floored at MHW. This prevents near-zero
#       or negative Rlow values. Rlow_m and Rlow_dam in the readable CSV are
#       always consistent: Rlow_dam = Rlow_m / 10.
#
#   MAX_STORMS_PER_YEAR = 5
#       Highest-Rhigh events are retained when trimming is needed.
#
# Adapted from: Storm_Creation.ipynb
# Author: Hannah Henry
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
import utide
from noaa_coops import Station

warnings.filterwarnings('ignore')

# =============================================================================
# USER CONFIGURATION — edit everything in this block
# =============================================================================

# --- Time period ---
# Run once per hindcast period, changing the four variables below each time.
# Period 1: BEGIN_DATE="19840101", END_DATE="20041231", START_YEAR=1984, OUTPUT_NAME="storms_1984_2004.npy"
# Period 2: BEGIN_DATE="20040101", END_DATE="20241231", START_YEAR=2004, OUTPUT_NAME="storms_2004_2024_base.npy"
BEGIN_DATE  = "20040101"
END_DATE    = "20241231"
START_YEAR  = 2004
OUTPUT_NAME = "storms_2004_2024.npy"

# --- NOAA Tide Gauge ---
# Duck, NC (Station 8651370) — closest long-record gauge on the Outer Banks.
# Downloaded automatically in annual chunks.
NOAA_STATION_ID = "8651370"   # Duck, NC
NOAA_DATUM      = "NAVD"      # NAVD88 — matches WIS and dune elevation data
NOAA_LAT        = 36.183      # Duck, NC latitude (for utide nodal corrections)

# --- WIS Wave Data ---
# CSV from WIS Generic Export (https://wisportal.erdc.dren.mil/), station ST63228.
# One file covering the full 1984–2024 period works for both hindcast periods.
WIS_PATH = Path(
    r"/scripts/input_prep/storm_creation_final/from_Hannah/storm_creation\WIS_raw_data\ST63228-Generic_Export-20260427T11T11_48.csv")

# --- Hatteras site parameters ---
# MHW_M must match MHW_ELEVATION in HAT_hindcast_1984_2024_old version.py.
# Set at berm crest so only events that reach or overtop the berm are included.
MHW_M           = 0.36        # [m NAVD88] — Duck NC gauge; must match hindcast script
BEACH_SLOPE     = 0.06        # foreshore slope (dimensionless) #0.2
BERM_CREST_M    = 1.7         # berm crest elevation [m NAVD88] (set to 1.5 before?)
STORM_THRESHOLD = 1.7   # berm crest [m NAVD88] — aligns with CASCADE's collision threshold

# --- Storm identification parameters ---
MIN_STORM_DURATION_HR  = 6    # discard events shorter than this [hours]
MIN_INTER_STORM_GAP_HR = 48   # merge events separated by less than this [hours]
                              # FIX: increased from 24 → 48 to prevent separate
                              # nor'easters from merging into 100–189 h mega-events

# --- Duration cap ---
# Hard ceiling on storm duration in the saved .npy file.
# Both working reference files (colleague and Benton/Ocracoke) have max
# duration of exactly 36 hours. Set here as the cap; any event longer is
# truncated. Acts as a backstop for the gap merging fix above.
MAX_STORM_DURATION_HR = 36    # [hours] — must not exceed ~36 for Barrier3D stability

# --- Per-year storm cap ---
# Retains highest-Rhigh events when a year exceeds MAX_STORMS_PER_YEAR.
# Barrier3D's internal limit is ~5 storms/year.
MAX_STORMS_PER_YEAR = 5

# --- Surge multiplier ---
# Scales non-tidal residual only — tidal signal unchanged:
#   water_level = eta_A  +  (SURGE_MULTIPLIER * eta_NTR)
SURGE_MULTIPLIER = 1.0

# --- Output ---
OUTPUT_DIR = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\fixed_storms")

# --- Plotting ---
SHOW_PLOTS = True
SAVE_PLOTS = True

# =============================================================================
# CONSTANTS
# =============================================================================

G = 9.81   # gravitational acceleration [m/s²]

# =============================================================================
# STEP 1 — NOAA TIDE GAUGE: DOWNLOAD, TIDAL DECOMPOSITION, SURGE SEPARATION
# =============================================================================

def load_noaa_water_levels(station_id: str, begin: str, end: str,
                            datum: str, lat: float) -> pd.DataFrame:
    """
    Download hourly NOAA water levels and decompose into tidal prediction
    (eta_A) and non-tidal residual (eta_NTR = storm surge proxy).

    Returns DataFrame with columns:
        observed_wl : raw observed water level (m NAVD88)
        eta_A       : tidal prediction (m)
        eta_NTR     : non-tidal residual / storm surge (m)
    """
    print(f"\n{'='*70}")
    print(f"STEP 1: Downloading NOAA Station {station_id} ({begin}–{end})")
    print(f"{'='*70}")

    station  = Station(id=station_id)
    begin_dt = pd.Timestamp(begin)
    end_dt   = pd.Timestamp(end)
    chunks   = []

    year = begin_dt.year
    while year <= end_dt.year:
        chunk_begin = max(begin_dt, pd.Timestamp(f"{year}-01-01"))
        chunk_end   = min(end_dt,   pd.Timestamp(f"{year}-12-31"))
        try:
            chunk = station.get_data(
                begin_date=chunk_begin.strftime("%Y%m%d"),
                end_date=chunk_end.strftime("%Y%m%d"),
                product="hourly_height",
                datum=datum,
                units="metric",
                time_zone="gmt"
            )
            chunks.append(chunk)
            print(f"  Downloaded {year}: {len(chunk)} records")
        except Exception as e:
            print(f"  WARNING: Failed to download {year} — {e}")
        year += 1

    df = pd.concat(chunks)
    df['water_level'] = pd.to_numeric(df['v'], errors='coerce')
    df['water_level'] = df['water_level'].interpolate(method='linear', limit=3)

    valid_mask  = ~np.isnan(df['water_level'])
    valid_time  = df.index[valid_mask]
    valid_water = df['water_level'][valid_mask].values

    print(f"\nTotal valid records: {valid_mask.sum():,} / {len(df):,}")
    print(f"Date range: {valid_time.min()} → {valid_time.max()}")

    print("Solving tidal constituents with utide (this may take ~1–2 min)...")
    coef = utide.solve(
        valid_time, valid_water,
        lat=lat, method='ols', conf_int='none', nodal=True, verbose=False
    )
    tide = utide.reconstruct(valid_time, coef, verbose=False)

    eta_A   = pd.Series(tide['h'],               index=valid_time, name='eta_A')
    eta_NTR = pd.Series(valid_water - tide['h'], index=valid_time, name='eta_NTR')

    result = pd.DataFrame({
        'observed_wl': valid_water,
        'eta_A':       eta_A.values,
        'eta_NTR':     eta_NTR.values
    }, index=valid_time)

    print(f"NTR range: {eta_NTR.min():.3f} to {eta_NTR.max():.3f} m")
    return result


# =============================================================================
# STEP 2 — WIS WAVE DATA: LOAD FROM CSV
# =============================================================================

def load_wis_data(filepath: Path) -> pd.DataFrame:
    """
    Load WIS hindcast CSV from wisportal.erdc.dren.mil Generic Export.
    Expected columns: time, waveHs, waveTp, waveMeanDirection.

    Returns DataFrame with columns: Hs (m), Tp (s), WAVD (deg).
    """
    print(f"\n{'='*70}")
    print(f"STEP 2: Loading WIS data from {filepath.name}")
    print(f"{'='*70}")

    df = pd.read_csv(filepath, parse_dates=['time'], index_col='time')

    out = pd.DataFrame({
        'Hs':   pd.to_numeric(df['waveHs'],           errors='coerce'),
        'Tp':   pd.to_numeric(df['waveTp'],            errors='coerce'),
        'WAVD': pd.to_numeric(df['waveMeanDirection'], errors='coerce')
    }, index=df.index)

    print(f"  Loaded {len(out):,} rows")
    print(f"  Date range: {out.index.min()} → {out.index.max()}")
    print(f"  Hs: {out['Hs'].min():.2f}–{out['Hs'].max():.2f} m  "
          f"| Tp: {out['Tp'].min():.2f}–{out['Tp'].max():.2f} s")
    return out


# =============================================================================
# STEP 3 — MERGE AND COMPUTE RUNUP (STOCKDON ET AL., 2006)
# =============================================================================

def compute_twl(tide_df: pd.DataFrame, wis_df: pd.DataFrame,
                slope: float, surge_mult: float) -> pd.DataFrame:
    """
    Merge tide and wave data; compute Stockdon runup components and TWL.

    Key columns produced:
        water_level : tide + scaled surge, no runup (m NAVD88) — used for Rlow
        Rhigh_m     : full 2% exceedance TWL (m NAVD88)
        TWL         : same as Rhigh_m (used for storm identification threshold)
    """
    print(f"\n{'='*70}")
    print(f"STEP 3: Merging data and computing runup (slope={slope})")
    print(f"{'='*70}")

    df = pd.DataFrame(index=wis_df.index)
    df['Hs']   = wis_df['Hs']
    df['Tp']   = wis_df['Tp']
    df['WAVD'] = wis_df['WAVD']

    df['eta_A']   = tide_df['eta_A'].reindex(wis_df.index).interpolate(
        method='linear', limit=6, limit_area='inside')
    df['eta_NTR'] = tide_df['eta_NTR'].reindex(wis_df.index).interpolate(
        method='linear', limit=6, limit_area='inside')

    if surge_mult != 1.0:
        print(f"  Applying surge multiplier: {surge_mult}x to eta_NTR")

    # Still water level: tide + scaled surge only — no wave component
    df['eta_NTR_scaled'] = df['eta_NTR'] * surge_mult
    df['water_level']    = df['eta_A'] + df['eta_NTR_scaled']

    df = df.dropna(subset=['Hs', 'Tp', 'water_level'])

    # Stockdon et al. (2006) runup components
    L0      = (G * df['Tp']**2) / (2 * np.pi)
    setup   = 0.35 * slope * np.sqrt(df['Hs'] * L0)
    S_inc   = 0.75 * slope * np.sqrt(df['Hs'] * L0)
    S_infra = 0.06 * np.sqrt(df['Hs'] * L0)
    S_total = np.sqrt(S_inc**2 + S_infra**2)
    R2      = 1.1 * (setup + S_total / 2)

    df['setup']   = setup
    df['R2']      = R2
    df['TWL']     = df['water_level'] + R2
    df['Rhigh_m'] = df['TWL']

    print(f"  Records after merge: {len(df):,}")
    print(f"  TWL range: {df['TWL'].min():.3f} to {df['TWL'].max():.3f} m")
    return df


# =============================================================================
# STEP 4 — STORM IDENTIFICATION
# =============================================================================

def identify_storms(df: pd.DataFrame, threshold: float,
                    min_duration_hr: int, min_gap_hr: int) -> list:
    """
    Identify discrete storm events as contiguous periods where TWL > threshold.

    1. Find all hours where TWL exceeds threshold.
    2. Merge events separated by fewer than min_gap_hr hours.
       min_gap_hr = 48 prevents separate nor'easters from merging into
       unrealistically long events that crash Barrier3D.
    3. Discard events shorter than min_duration_hr hours.

    Returns list of (start_timestamp, end_timestamp) pairs.
    """
    print(f"\n{'='*70}")
    print(f"STEP 4: Storm identification (threshold={threshold} m, "
          f"min_dur={min_duration_hr} h, min_gap={min_gap_hr} h)")
    print(f"{'='*70}")

    exceed = (df['TWL'] > threshold).astype(int)
    print(f"  Exceedance rate: {exceed.mean()*100:.2f}% of hours above threshold")

    diff       = exceed.diff().fillna(exceed)
    starts_raw = df.index[diff == 1].tolist()
    ends_raw   = df.index[diff == -1].tolist()

    if exceed.iloc[0] == 1:
        starts_raw.insert(0, df.index[0])
    if exceed.iloc[-1] == 1:
        ends_raw.append(df.index[-1])

    if len(starts_raw) == 0:
        print("  WARNING: No storm events identified. Lower STORM_THRESHOLD.")
        return []

    merged_starts = [starts_raw[0]]
    merged_ends   = [ends_raw[0]]
    for s, e in zip(starts_raw[1:], ends_raw[1:]):
        gap_hours = (s - merged_ends[-1]).total_seconds() / 3600
        if gap_hours <= min_gap_hr:
            merged_ends[-1] = e
        else:
            merged_starts.append(s)
            merged_ends.append(e)

    events = []
    for s, e in zip(merged_starts, merged_ends):
        dur_hr = (e - s).total_seconds() / 3600
        if dur_hr >= min_duration_hr:
            events.append((s, e))

    print(f"  Raw events:            {len(starts_raw)}")
    print(f"  After gap merging:     {len(merged_starts)}")
    print(f"  After duration filter: {len(events)}")
    return events


# =============================================================================
# STEP 5 — EXTRACT PER-STORM PARAMETERS
# =============================================================================

def extract_storm_params(df: pd.DataFrame, events: list,
                          start_year: int, mhw_m: float,
                          max_duration_hr: float) -> tuple:
    """
    For each storm event, extract CASCADE storm parameters.

    Rlow (KEY FIX):
        Rlow = max(water_level_at_peak_TWL, mhw_m) / 10
        where water_level = tide + surge only (no runup, no setup).
        Flooring at MHW ensures Rlow is never below the tidal datum.
        Rlow_m and Rlow_dam are always consistent: Rlow_dam = Rlow_m / 10.

    Duration (KEY FIX):
        Duration = min(actual exceedance hours, max_duration_hr)
        Both working reference files cap at 36 h. Raw duration before capping
        is preserved in the readable CSV as Raw_Duration_h for QA.

    Returns:
        arr      : (N, 5) float64 array: [Year_Index, Rhigh_dam, Rlow_dam, Tp_s, Duration_h]
        readable : DataFrame with human-readable QA columns
    """
    print(f"\n{'='*70}")
    print(f"STEP 5: Extracting per-storm parameters (duration cap={max_duration_hr} h)")
    print(f"{'='*70}")

    rows              = []
    readable_rows     = []
    n_duration_capped = 0

    for s, e in events:
        storm    = df.loc[s:e]
        peak_idx = storm['TWL'].idxmax()

        year_index = peak_idx.year - start_year

        # Duration: actual exceedance length, capped at max_duration_hr
        raw_duration = (e - s).total_seconds() / 3600
        duration     = min(raw_duration, max_duration_hr)
        if raw_duration > max_duration_hr:
            n_duration_capped += 1

        # Rhigh: maximum 2% exceedance TWL during storm (m → dam)
        rhigh_m   = storm['Rhigh_m'].max()
        rhigh_dam = rhigh_m / 10.0

        # Rlow: still water level at peak TWL, floored at MHW
        rlow_m   = max(storm.loc[peak_idx, 'water_level'], mhw_m)
        rlow_dam = rlow_m / 10.0

        # Sanity: Rlow must be < Rhigh
        if rlow_dam >= rhigh_dam:
            print(f"  WARNING: Rlow ({rlow_dam:.4f}) >= Rhigh ({rhigh_dam:.4f}) "
                  f"at {peak_idx} — clamping Rlow to 0.95 × Rhigh")
            rlow_dam = 0.95 * rhigh_dam
            rlow_m   = rlow_dam * 10.0

        tp_mean = storm['Tp'].mean()
        hs_peak = storm['Hs'].max()

        rows.append([year_index, rhigh_dam, rlow_dam, tp_mean, duration])
        readable_rows.append({
            'Storm_Start'    : s.strftime('%Y-%m-%d %H:%M'),
            'Storm_End'      : e.strftime('%Y-%m-%d %H:%M'),
            'Peak_TWL_Time'  : peak_idx.strftime('%Y-%m-%d %H:%M'),
            'Calendar_Year'  : peak_idx.year,
            'Year_Index'     : int(year_index),
            'Rhigh_m'        : round(rhigh_m, 3),
            'Rlow_m'         : round(rlow_m, 3),       # still WL at peak TWL, >= MHW
            'Rhigh_dam'      : round(rhigh_dam, 4),
            'Rlow_dam'       : round(rlow_dam, 4),      # always = Rlow_m / 10
            'Wave_Period_s'  : round(tp_mean, 2),
            'Duration_h'     : duration,                # capped at max_duration_hr
            'Raw_Duration_h' : round(raw_duration, 1),  # actual exceedance (QA only)
            'Peak_Hs_m'      : round(hs_peak, 3),
        })

    arr      = np.array(rows, dtype=np.float64)
    readable = pd.DataFrame(readable_rows)

    print(f"  Total events identified: {len(arr)}")
    print(f"  Duration-capped events:  {n_duration_capped} "
          f"(raw > {max_duration_hr} h, truncated to cap)")
    if len(arr) > 0:
        print(f"  Year_Index:  {arr[:,0].min():.0f}–{arr[:,0].max():.0f} "
              f"({int(arr[:,0].min())+start_year}–{int(arr[:,0].max())+start_year})")
        print(f"  Rhigh (dam): {arr[:,1].min():.4f}–{arr[:,1].max():.4f}")
        print(f"  Rlow  (dam): {arr[:,2].min():.4f}–{arr[:,2].max():.4f}")
        print(f"  Duration(h): {arr[:,4].min():.1f}–{arr[:,4].max():.1f}  "
              f"(before cap: {readable['Raw_Duration_h'].max():.1f} h max)")

    return arr, readable


# =============================================================================
# STEP 5b — PER-YEAR STORM CAP
# =============================================================================

def cap_storms_per_year(readable: pd.DataFrame, max_storms: int,
                         start_year: int) -> pd.DataFrame:
    """
    Limit storms to max_storms per Year_Index, retaining highest-Rhigh events.

    Barrier3D's internal per-year array limit is ~5 storms. When trimming,
    the most morphologically significant storms (highest Rhigh) are kept.
    """
    before = len(readable)

    capped = (
        readable
        .sort_values('Rhigh_m', ascending=False)
        .groupby('Year_Index', group_keys=False)
        .head(max_storms)
        .sort_values(['Year_Index', 'Rhigh_m'], ascending=[True, False])
        .reset_index(drop=True)
    )

    after         = len(capped)
    trimmed_years = readable.groupby('Year_Index').size()
    trimmed_years = trimmed_years[trimmed_years > max_storms]

    print(f"\n{'='*70}")
    print(f"STEP 5b: Per-year storm cap (max={max_storms})")
    print(f"{'='*70}")
    print(f"  Storms before cap: {before}")
    print(f"  Storms after cap:  {after}  (removed {before - after})")
    if len(trimmed_years) > 0:
        print(f"  Years trimmed:")
        for yi, n in trimmed_years.items():
            print(f"    Year_Index {int(yi):2d} ({int(yi)+start_year}): "
                  f"{n} → {min(n, max_storms)} storms kept")
    else:
        print("  No years exceeded the cap — all storms retained.")
    return capped


# =============================================================================
# STEP 6 — DIAGNOSTICS AND VISUALIZATION
# =============================================================================

def plot_twl_timeseries(df: pd.DataFrame, events: list,
                         threshold: float, berm_crest: float,
                         save_dir: Path = None):
    """TWL time series with identified storm events highlighted."""
    fig, ax = plt.subplots(figsize=(18, 5))
    ax.plot(df.index, df['TWL'], color='steelblue', lw=0.5, alpha=0.8, label='TWL')
    ax.axhline(threshold,  color='orange', ls='--', lw=1.5,
               label=f'Threshold ({threshold} m)')
    ax.axhline(berm_crest, color='red',    ls='--', lw=1.5,
               label=f'Berm crest ({berm_crest} m)')
    ax.axhline(MHW_M,      color='gray',   ls=':',  lw=1.0,
               label=f'MHW ({MHW_M} m)')
    for s, e in events:
        ax.axvspan(s, e, alpha=0.25, color='red')
    ax.set_ylabel('Elevation (m NAVD88)')
    ax.set_xlabel('Date')
    ax.set_title('Total Water Level and Identified Storm Events')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_dir and SAVE_PLOTS:
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / "storm_twl_timeseries.png", dpi=200, bbox_inches='tight')
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


def plot_storm_distributions(arr: np.ndarray, start_year: int,
                              save_dir: Path = None):
    """Four-panel distribution plot of storm parameters."""
    df = pd.DataFrame(arr, columns=['Year_Index','Rhigh','Rlow','Wave Period','Duration'])
    df['Calendar_Year'] = df['Year_Index'].astype(int) + start_year
    yr_min = int(df['Calendar_Year'].min())
    yr_max = int(df['Calendar_Year'].max())

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f"Storm Parameter Distributions — Hatteras "
                 f"(N={len(df)}, {yr_min}–{yr_max})", fontsize=14)

    df['Rhigh'].hist(ax=axs[0,0], bins=20, color='#1f77b4', edgecolor='k', alpha=0.8)
    axs[0,0].set(title='Rhigh Distribution (dam)', xlabel='Rhigh (dam)')

    df['Wave Period'].hist(ax=axs[0,1], bins=20, color='#ff7f0e', edgecolor='k', alpha=0.8)
    axs[0,1].set(title='Wave Period Distribution', xlabel='Period (s)')

    df['Duration'].hist(ax=axs[1,0], bins=20, color='#2ca02c', edgecolor='k', alpha=0.8)
    axs[1,0].set(title=f'Duration Distribution (cap={MAX_STORM_DURATION_HR} h)',
                 xlabel='Duration (h)')
    axs[1,0].axvline(MAX_STORM_DURATION_HR, color='red', ls='--', lw=1.5,
                     label=f'Cap ({MAX_STORM_DURATION_HR} h)')
    axs[1,0].legend(fontsize=8)

    axs[1,1].scatter(df['Rlow'], df['Rhigh'], alpha=0.7, color='#d62728', s=30)
    lim = max(df['Rhigh'].max(), df['Rlow'].max()) * 1.05
    axs[1,1].plot([0, lim], [0, lim], 'k--', lw=0.8, alpha=0.4, label='1:1')
    axs[1,1].set(title='Rlow vs Rhigh', xlabel='Rlow (dam)', ylabel='Rhigh (dam)')
    axs[1,1].legend(fontsize=8)
    axs[1,1].grid(True, ls='--', alpha=0.5)

    for ax in axs.flat:
        ax.grid(alpha=0.3)
        ax.set_ylabel('Count')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    if save_dir and SAVE_PLOTS:
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / "storm_distributions_hatteras.png", dpi=300, bbox_inches='tight')
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


def plot_storms_by_year(readable: pd.DataFrame, start_year: int,
                        berm_crest: float, surge_mult: float,
                        save_dir: Path = None):
    """Two-panel annual summary for historical cross-checking."""
    years     = sorted(readable['Calendar_Year'].unique())
    all_years = list(range(int(min(years)), int(max(years)) + 1))
    counts    = readable.groupby('Calendar_Year').size().reindex(all_years, fill_value=0)
    yr_min, yr_max = min(all_years), max(all_years)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9),
                                    gridspec_kw={'height_ratios': [1, 2.5]})
    fig.suptitle(
        f"Storm record — Hatteras Island ({yr_min}–{yr_max}, "
        f"surge mult={surge_mult}×, N={len(readable)}, "
        f"dur cap={MAX_STORM_DURATION_HR} h, gap={MIN_INTER_STORM_GAP_HR} h)",
        fontsize=12
    )

    bar_colors = ['#d62728' if c > 0 else '#cccccc' for c in counts.values]
    ax1.bar(all_years, counts.values, color=bar_colors, edgecolor='white', lw=0.5)
    ax1.set_ylabel('Storms / year', fontsize=10)
    ax1.set_ylim(0, max(counts.max() + 1, MAX_STORMS_PER_YEAR + 2))
    if counts[counts > 0].any():
        ax1.axhline(counts[counts > 0].mean(), color='gray', ls='--', lw=1,
                    label=f'Mean ({counts[counts>0].mean():.1f}/yr, active years)')
    ax1.axhline(MAX_STORMS_PER_YEAR, color='steelblue', ls=':', lw=1.5,
                label=f'Cap ({MAX_STORMS_PER_YEAR}/yr)')
    ax1.legend(fontsize=8, loc='upper right')
    ax1.grid(axis='y', alpha=0.3)
    ax1.set_xticks(all_years)
    ax1.set_xticklabels(all_years, rotation=45, ha='right', fontsize=8)

    sc = ax2.scatter(
        readable['Calendar_Year'], readable['Rhigh_m'],
        c=readable['Peak_Hs_m'], s=readable['Duration_h'] * 2,
        cmap='YlOrRd', alpha=0.85, edgecolors='k', linewidths=0.4, zorder=3
    )
    ax2.axhline(berm_crest,       color='steelblue', ls='--', lw=1.5,
                label=f'Berm crest ({berm_crest} m)')
    ax2.axhline(STORM_THRESHOLD,  color='orange',    ls=':',  lw=1.5,
                label=f'Threshold ({STORM_THRESHOLD} m)')
    ax2.set_ylabel('Rhigh (m NAVD88)', fontsize=10)
    ax2.set_xlabel('Year', fontsize=10)
    ax2.set_xlim(yr_min - 0.8, yr_max + 0.8)
    ax2.set_xticks(all_years)
    ax2.set_xticklabels(all_years, rotation=45, ha='right', fontsize=8)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=8, loc='upper left')

    for yr, grp in readable.groupby('Calendar_Year'):
        peak = grp.loc[grp['Rhigh_m'].idxmax()]
        ax2.annotate(f"{peak['Rhigh_m']:.2f}", xy=(yr, peak['Rhigh_m']),
                     xytext=(0, 6), textcoords='raw_offset points',
                     ha='center', fontsize=7, color='#333333')

    cbar = plt.colorbar(sc, ax=ax2, pad=0.01, shrink=0.8)
    cbar.set_label('Peak Hs (m)', fontsize=9)
    for dur, label in [(12, '12 h'), (24, '24 h'), (36, '36 h')]:
        ax2.scatter([], [], s=dur*2, c='gray', alpha=0.6,
                    edgecolors='k', linewidths=0.4, label=f'Duration {label}')
    ax2.legend(fontsize=8, loc='upper left', ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    if save_dir and SAVE_PLOTS:
        save_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_dir / "storm_annual_summary.png", dpi=300, bbox_inches='tight')
    if SHOW_PLOTS:
        plt.show()
    else:
        plt.close()


def print_annual_summary(readable: pd.DataFrame, start_year: int):
    """Print per-year summary for cross-checking against historical record."""
    yr_min = int(readable['Calendar_Year'].min())
    yr_max = int(readable['Calendar_Year'].max())
    print(f"\n{'='*60}")
    print("STORMS PER YEAR (post-cap)")
    print(f"{'='*60}")
    for yr in range(yr_min, yr_max + 1):
        sub = readable[readable['Calendar_Year'] == yr]
        n   = len(sub)
        if n > 0:
            print(f"  {yr}: {n:2d} storm(s)  "
                  f"Rhigh_max={sub['Rhigh_m'].max():.3f} m  "
                  f"dur_max={sub['Duration_h'].max():.0f} h  "
                  f"(raw_max={sub['Raw_Duration_h'].max():.0f} h)")
        else:
            print(f"  {yr}:  0 storms")
    total = len(readable)
    span  = yr_max - yr_min + 1
    print(f"\nTotal: {total} events over {span} years ({total/span:.1f}/year)")


def validate_storm_array(arr: np.ndarray, readable: pd.DataFrame,
                          mhw_m: float, berm_m: float,
                          max_dur: float) -> bool:
    """Pre-flight checks before saving. Must all pass before file is written."""
    print(f"\n{'='*70}")
    print("VALIDATION CHECKS")
    print(f"{'='*70}")

    mhw_dam = mhw_m / 10.0
    passed  = True

    # Rlow < Rhigh
    bad = arr[:, 2] >= arr[:, 1]
    if bad.any():
        print(f"  ❌ {bad.sum()} storm(s) with Rlow >= Rhigh")
        passed = False
    else:
        print(f"  ✓ All Rlow < Rhigh")

    # Rlow >= MHW
    below_mhw = arr[:, 2] < mhw_dam
    if below_mhw.any():
        print(f"  ❌ {below_mhw.sum()} storm(s) with Rlow < MHW ({mhw_dam:.4f} dam)")
        passed = False
    else:
        print(f"  ✓ All Rlow >= MHW ({mhw_dam:.4f} dam)")

    # No negatives
    if (arr[:, 1:3] < 0).any():
        print(f"  ❌ Negative Rhigh or Rlow found")
        passed = False
    else:
        print(f"  ✓ No negative Rhigh/Rlow")

    # Duration cap
    over_dur = arr[:, 4] > max_dur
    if over_dur.any():
        print(f"  ❌ {over_dur.sum()} storm(s) exceed duration cap ({max_dur} h)")
        passed = False
    else:
        print(f"  ✓ All durations <= {max_dur} h")

    # Per-year count
    by_year  = readable.groupby('Year_Index').size()
    over_lim = by_year[by_year > 5]
    if len(over_lim) > 0:
        print(f"  ❌ {len(over_lim)} year(s) exceed Barrier3D ~5 storm limit")
        passed = False
    else:
        print(f"  ✓ All years within Barrier3D limit (max={by_year.max()} storms/year)")

    # Rlow_m / Rlow_dam consistency
    expected = (readable['Rlow_m'] / 10).round(4)
    mismatch = ~np.isclose(readable['Rlow_dam'], expected, rtol=0.01)
    if mismatch.any():
        print(f"  ❌ {mismatch.sum()} row(s) where Rlow_dam != Rlow_m/10")
        passed = False
    else:
        print(f"  ✓ Rlow_dam = Rlow_m/10 for all {len(readable)} rows")

    # Isabel 2003 check
    isabel = readable[readable['Calendar_Year'] == 2003]
    if len(isabel) > 0:
        max_yr  = readable.loc[readable['Rhigh_m'].idxmax(), 'Calendar_Year']
        max_val = readable['Rhigh_m'].max()
        if max_yr == 2003:
            print(f"  ✓ Isabel 2003 is max Rhigh ({max_val:.3f} m) — as expected")
        else:
            print(f"  ⚠ Max Rhigh is {max_val:.3f} m in {max_yr}, not 2003 — verify")
    else:
        print(f"  ⚠ No 2003 storms found — Isabel should appear in 1984–2004 record")

    # Rlow/Rhigh ratio
    ratio = arr[:, 2] / arr[:, 1]
    print(f"  Rlow/Rhigh ratio: mean={ratio.mean():.3f}  "
          f"min={ratio.min():.3f}  max={ratio.max():.3f}  "
          f"(reference files: ~0.29–0.35)")

    print()
    if passed:
        print("  ✅ All checks passed — storm file ready for CASCADE.")
    else:
        print("  ❌ Validation issues found — resolve before running CASCADE.")
    return passed


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\n" + "="*70)
    print("HAT_create_storm_file.py — CASCADE Storm File Generator")
    print(f"Period:           {BEGIN_DATE} → {END_DATE}  |  Start year: {START_YEAR}")
    print(f"Surge multiplier: {SURGE_MULTIPLIER}x")
    print(f"Threshold:        {STORM_THRESHOLD} m  |  Min gap: {MIN_INTER_STORM_GAP_HR} h")
    print(f"Duration cap:     {MAX_STORM_DURATION_HR} h  |  Year cap: {MAX_STORMS_PER_YEAR} storms/yr")
    print("="*70)

    # 1. NOAA water levels
    tide_df = load_noaa_water_levels(
        NOAA_STATION_ID, BEGIN_DATE, END_DATE, NOAA_DATUM, NOAA_LAT
    )

    # 2. WIS wave data
    wis_df = load_wis_data(WIS_PATH)
    wis_df = wis_df.loc[BEGIN_DATE:END_DATE]

    # 3. Merge + compute runup
    df = compute_twl(tide_df, wis_df, BEACH_SLOPE, SURGE_MULTIPLIER)
    df = df.loc[BEGIN_DATE:END_DATE]

    # 4. Storm identification (48 h gap prevents nor'easter merging)
    events = identify_storms(df, STORM_THRESHOLD,
                              MIN_STORM_DURATION_HR, MIN_INTER_STORM_GAP_HR)
    if not events:
        print("\nNo storm events found. Try lowering STORM_THRESHOLD.")
        return None

    # 5. Extract parameters — Rlow fixed, duration capped at 36 h
    storm_arr, storm_readable = extract_storm_params(
        df, events, START_YEAR, MHW_M, MAX_STORM_DURATION_HR
    )

    # 5b. Per-year cap
    storm_readable_capped = cap_storms_per_year(
        storm_readable, MAX_STORMS_PER_YEAR, START_YEAR
    )

    # Rebuild npy from capped readable to keep them in sync
    storm_arr_capped = storm_readable_capped[
        ['Year_Index', 'Rhigh_dam', 'Rlow_dam', 'Wave_Period_s', 'Duration_h']
    ].to_numpy(dtype=np.float64)

    # 6. Diagnostics
    print_annual_summary(storm_readable_capped, START_YEAR)
    passes = validate_storm_array(
        storm_arr_capped, storm_readable_capped, MHW_M, BERM_CREST_M, MAX_STORM_DURATION_HR
    )

    if SHOW_PLOTS or SAVE_PLOTS:
        plot_twl_timeseries(df, events, STORM_THRESHOLD, BERM_CREST_M,
                            OUTPUT_DIR if SAVE_PLOTS else None)
        plot_storm_distributions(storm_arr_capped, START_YEAR,
                                  OUTPUT_DIR if SAVE_PLOTS else None)
        plot_storms_by_year(storm_readable_capped, START_YEAR, BERM_CREST_M,
                            SURGE_MULTIPLIER, OUTPUT_DIR if SAVE_PLOTS else None)

    # 7. Save — only if all validation checks pass
    if not passes:
        print("\n⚠️  Save skipped — resolve validation issues first.")
        return None

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    out_path = OUTPUT_DIR / OUTPUT_NAME
    np.save(out_path, storm_arr_capped)
    print(f"\n✓ Saved NPY: {out_path}")
    print(f"  Shape: {storm_arr_capped.shape}  |  dtype: {storm_arr_capped.dtype}")
    print(f"  Columns: [Year_Index, Rhigh_dam, Rlow_dam, Wave_Period_s, Duration_h]")

    csv_name = OUTPUT_NAME.replace('.npy', '_readable.csv')
    csv_path = OUTPUT_DIR / csv_name
    storm_readable_capped.to_csv(csv_path, index=False)
    print(f"\n✓ Saved CSV: {csv_path}")
    print(f"  Raw_Duration_h = actual exceedance length before duration cap (QA only)")
    print(f"  Rlow_m  = still water level at peak TWL, floored at MHW ({MHW_M} m)")
    print(f"  Rlow_dam = Rlow_m / 10 — consistent with stored .npy values")

    return storm_arr_capped


if __name__ == "__main__":
    storm_arr = main()