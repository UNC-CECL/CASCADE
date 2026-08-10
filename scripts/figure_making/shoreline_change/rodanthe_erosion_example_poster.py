"""
CoastSat Shoreline Erosion Trends — Rodanthe Area
Hatteras Island, NC

Produces a single square publication-quality figure for poster use showing:
  - Annual mean shoreline position per domain (thin lines) — shows the raw signal
  - LRR trend line per domain (bold, full-period) — shows the erosion rate

Two messages conveyed:
  1. Rodanthe has been eroding consistently over the full record
  2. Different parts of Rodanthe erode at different rates

Data source: CoastSat satellite-derived shorelines

Usage
-----
  1. Edit the CONFIG section below.
  2. Run: python plot_coastsat_rodanthe_erosion.py
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# ============================================================================
# CONFIGURATION
# ============================================================================

ROOT_DATA_DIR = r"/scripts/input_prep/CoastSat/coastsat_timeseries"
LOOKUP_CSV    = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"

DOMAIN_MIN  = 77
DOMAIN_MAX  = 83
SITE_FILTER = "usa_NC"
START_DATE  = "1984-01-01"
END_DATE    = "2024-12-31"
MIN_OBS     = 5

OUTPUT_DIR  = r"/scripts/input_prep/CoastSat/rodanthe_plots"
OUTPUT_FILE = "rodanthe_erosion_trends.png"

FIG_SIZE = (8, 8)

DOMAIN_COLORS = {
    77: '#1a6b8a',
    78: '#2196a6',
    79: '#43b89c',
    80: '#8fbe6a',
    81: '#d4a030',
    82: '#d4622a',
    83: '#c0392b',
}

# ============================================================================
# DATA LOADING
# ============================================================================

def collect_csv_map(root_dir, site_filter=""):
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"WARNING: ROOT_DATA_DIR not found: {root_dir}")
        return csv_map
    for d in sorted(os.listdir(root_dir)):
        subfolder = os.path.join(root_dir, d)
        if os.path.isdir(subfolder) and (site_filter == "" or site_filter in d):
            for fpath in glob.glob(os.path.join(subfolder, "*.csv")):
                stem = os.path.splitext(os.path.basename(fpath))[0]
                csv_map[stem] = fpath
    print(f"Found {len(csv_map)} total time-series CSVs.")
    return csv_map


def load_timeseries(filepath):
    df = pd.read_csv(filepath, header=0)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "chainage_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True)
    df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
    return df.dropna(subset=["date", "chainage_m"]).sort_values("date").reset_index(drop=True)


def filter_dates(df, start, end):
    if start:
        df = df[df["date"] >= pd.Timestamp(start, tz="UTC")]
    if end:
        df = df[df["date"] <= pd.Timestamp(end, tz="UTC")]
    return df.reset_index(drop=True)


def load_data(root_dir, lookup_csv, domain_min, domain_max,
              site_filter, start_date, end_date, min_obs):
    lookup = pd.read_csv(lookup_csv)
    lookup.columns = [c.split('csv')[-1] if 'csv' in c else c for c in lookup.columns]
    mask   = (lookup["domain_number"] >= domain_min) & (lookup["domain_number"] <= domain_max)
    lookup = lookup[mask].copy()
    print(f"Domains {domain_min}–{domain_max}: {len(lookup)} transects in lookup.")

    csv_map = collect_csv_map(root_dir, site_filter)
    records = []
    for _, row in lookup.iterrows():
        tid    = str(row["transect_id"])
        domain = int(row["domain_number"])
        if tid not in csv_map:
            continue
        try:
            df_raw  = load_timeseries(csv_map[tid])
            df_filt = filter_dates(df_raw, start_date, end_date)
            if len(df_filt) < min_obs:
                continue
            df_filt = df_filt.copy()
            df_filt["transect_id"]   = tid
            df_filt["domain_number"] = domain
            records.append(df_filt)
        except Exception as e:
            print(f"  ERROR {tid}: {e}")

    if not records:
        raise RuntimeError("No data loaded. Check ROOT_DATA_DIR and LOOKUP_CSV.")

    df_all = pd.concat(records, ignore_index=True)
    df_all["year"] = df_all["date"].dt.year
    print(f"Loaded {len(df_all)} obs from {df_all['transect_id'].nunique()} transects, "
          f"{df_all['domain_number'].nunique()} domains.")
    return df_all


# ============================================================================
# COMPUTE DOMAIN ANNUAL MEANS AND LRR TREND LINES
# ============================================================================

def compute_domain_stats(df_all):
    domain_year = (
        df_all.groupby(["domain_number", "year"])["chainage_m"]
        .mean()
        .reset_index()
        .rename(columns={"chainage_m": "mean_chainage"})
    )

    # Offset each domain to start at 0 so slope differences are readable
    first_vals  = domain_year.groupby("domain_number")["mean_chainage"].first()
    domain_year = domain_year.merge(first_vals.rename("first_val"), on="domain_number")
    domain_year["position"] = domain_year["mean_chainage"] - domain_year["first_val"]

    lrr_stats = {}
    for domain, grp in domain_year.groupby("domain_number"):
        grp = grp.sort_values("year")
        x   = grp["year"].values.astype(float)
        y   = grp["position"].values
        if len(x) < 3:
            continue
        slope, intercept, r, p, se = stats.linregress(x, y)
        dof = len(x) - 2
        from scipy.stats import t as t_dist
        unc = t_dist.ppf(0.975, dof) * se if dof > 0 else np.nan
        lrr_stats[domain] = {
            "lrr_m_yr": round(slope, 2),
            "unc"      : round(unc, 2),
            "r2"       : round(r**2, 3),
            "x_line"   : np.array([x.min(), x.max()]),
            "y_line"   : slope * np.array([x.min(), x.max()]) + intercept,
        }

    return domain_year, lrr_stats


# ============================================================================
# PLOT
# ============================================================================

def make_figure(domain_year, lrr_stats, domain_colors, domain_min, domain_max,
                start_date, end_date, output_path, fig_size):

    domains = sorted(domain_year["domain_number"].unique())

    fig, ax = plt.subplots(figsize=fig_size, facecolor='white')
    ax.set_facecolor('white')

    for domain in domains:
        color = domain_colors.get(domain, '#888888')
        grp   = domain_year[domain_year["domain_number"] == domain].sort_values("year")

        # Thin, transparent annual mean line — no markers
        ax.plot(grp["year"], grp["position"],
                color=color, linewidth=0.9, alpha=0.25, zorder=2)

        # Bold LRR trend line
        if domain in lrr_stats:
            st = lrr_stats[domain]
            ax.plot(st["x_line"], st["y_line"],
                    color=color, linewidth=2.8, alpha=0.95, zorder=4,
                    solid_capstyle='round')

    # Zero reference
    ax.axhline(0, color='#2c2c2c', linewidth=1.0, linestyle='--', alpha=0.5, zorder=3)

    # Light erosion zone shading
    ax.fill_between(
        [pd.Timestamp(start_date).year, pd.Timestamp(end_date).year],
        [ax.get_ylim()[0]] * 2, [0, 0],
        color='#d6e8f0', alpha=0.15, zorder=0
    )

    # Axes formatting
    ax.set_xlim(pd.Timestamp(start_date).year - 0.5,
                pd.Timestamp(end_date).year + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='both', which='major', labelsize=10, direction='in', length=5)
    ax.tick_params(axis='both', which='minor', direction='in', length=3)
    ax.grid(True, which='major', linestyle=':', linewidth=0.5, alpha=0.4, color='gray')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(1.2)

    ax.set_xlabel("Year", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel("Shoreline position change (m)", fontsize=12, fontweight='bold', labelpad=8)

    # Erosion / Accretion labels
    ybot, ytop = ax.get_ylim()
    zero_frac  = (0 - ybot) / (ytop - ybot)
    ax.text(1.01, zero_frac + (1 - zero_frac) * 0.45, 'Accretion ▲',
            transform=ax.transAxes, fontsize=9, color='#555555',
            ha='left', va='center', style='italic')
    ax.text(1.01, zero_frac * 0.45, 'Erosion ▼',
            transform=ax.transAxes, fontsize=9, color='#555555',
            ha='left', va='center', style='italic')

    # Legend — one entry per domain with rate, plus a style note
    domain_handles = [
        Line2D([0], [0], color=domain_colors.get(d, '#888888'), linewidth=2.5,
               label=f'Domain {d}  ({lrr_stats[d]["lrr_m_yr"]:+.1f} m/yr)'
               if d in lrr_stats else f'Domain {d}')
        for d in domains
    ]
    style_note = Line2D([0], [0], color='#888888', linewidth=0.9, alpha=0.4,
                        label='Annual domain mean')
    ax.legend(handles=[style_note] + domain_handles,
              fontsize=9, framealpha=0.95, edgecolor='#cccccc',
              loc='lower left', ncol=1)

    ax.set_title(
        f'Shoreline Change — Rodanthe, Hatteras Island\n'
        f'CoastSat satellite-derived shorelines  '
        f'({pd.Timestamp(start_date).year}–{pd.Timestamp(end_date).year})',
        fontsize=12, fontweight='bold', color='#1a2a3a', pad=10
    )

    fig.text(
        0.01, 0.005,
        'Shoreline position relative to each domain\'s first observation year. '
        'Thin lines = annual domain mean; bold lines = LRR trend. '
        'Shoreline data: CoastSat (Vos et al., 2019).',
        fontsize=7.5, color='#666666', style='italic', ha='left', va='bottom'
    )

    plt.tight_layout(rect=[0, 0.03, 0.97, 1])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"\nSaved: {output_path}")
    plt.show()
    plt.close()


# ============================================================================
# MAIN
# ============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print("=" * 60)
    print(f"Rodanthe Erosion Trends  (Domains {DOMAIN_MIN}–{DOMAIN_MAX})")
    print(f"Period: {START_DATE} → {END_DATE}")
    print("=" * 60)

    df_all = load_data(
        ROOT_DATA_DIR, LOOKUP_CSV, DOMAIN_MIN, DOMAIN_MAX,
        SITE_FILTER, START_DATE, END_DATE, MIN_OBS
    )

    domain_year, lrr_stats = compute_domain_stats(df_all)

    print("\nDomain LRR summary:")
    print(f"  {'Domain':>8}  {'LRR (m/yr)':>10}  {'±95% CI':>8}  {'R²':>6}")
    print(f"  {'-'*8}  {'-'*10}  {'-'*8}  {'-'*6}")
    for d, st in sorted(lrr_stats.items()):
        print(f"  {d:>8}  {st['lrr_m_yr']:>+10.2f}  {st['unc']:>8.2f}  {st['r2']:>6.3f}")

    make_figure(domain_year, lrr_stats, DOMAIN_COLORS,
                DOMAIN_MIN, DOMAIN_MAX, START_DATE, END_DATE,
                os.path.join(OUTPUT_DIR, OUTPUT_FILE), FIG_SIZE)


if __name__ == "__main__":
    main()
