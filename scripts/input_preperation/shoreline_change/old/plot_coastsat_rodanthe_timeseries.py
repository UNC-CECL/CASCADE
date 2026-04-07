"""
CoastSat Shoreline Change Visualization — Rodanthe Area
Hatteras Island, NC

Produces three square publication-quality figures for poster use,
all focused on the Rodanthe community (CASCADE domains 77–83):

  Figure 1: Space-time heatmap
            X = time, Y = domain/transect, color = shoreline position anomaly
            Shows when and where shoreline change occurred.

  Figure 2: Spaghetti plot
            Cross-shore position vs. time for each transect, one line per transect,
            colored by domain. Shows the raw time series signal.

  Figure 3: Mean shoreline position per domain per year
            Aggregated to annual means per domain, shown as colored lines.
            Cleaner version of Figure 2 for poster use.

Data source: CoastSat satellite-derived shorelines

Usage
-----
  1. Edit the CONFIG section below (ROOT_DATA_DIR, LOOKUP_CSV, date range).
  2. Run: python plot_coastsat_rodanthe_timeseries.py
"""

import os
import glob
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# ============================================================================
# CONFIGURATION — edit these paths and parameters
# ============================================================================

# Root folder containing site subfolders of CoastSat CSVs
# e.g. C:/…/coastsat_timeseries/usa_NC_0032_timeseries/usa_NC_0032_0011.csv
ROOT_DATA_DIR = r"/scripts/input_preperation/CoastSat/coastsat_timeseries"

# Transect → domain lookup table (produced by coastsat_domain_mapping.py)
LOOKUP_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\transect_domain_lookup.csv"

# Domains to include (Rodanthe area)
DOMAIN_MIN = 77
DOMAIN_MAX = 83

# Full time range to display (all available CoastSat data)
START_DATE = "1984-01-01"
END_DATE   = "2024-12-31"

# Optional: filter to just one period for the plots — set both to None to show all
# e.g. START_DATE = "2004-01-01", END_DATE = "2024-12-31"

# Site filter for subfolder discovery
SITE_FILTER = "usa_NC"

# Minimum observations per transect to include in plots
MIN_OBS = 5

# Output directory
OUTPUT_DIR = r"/scripts/input_preperation/CoastSat/rodanthe_plots"

# Figure size — square for poster
FIG_SIZE = (8, 8)

# Color palette for domain lines (one color per domain 77–83)
DOMAIN_COLORS = {
    77: '#1a6b8a',  # dark teal
    78: '#2196a6',  # mid teal
    79: '#43b89c',  # teal-green
    80: '#8fbe6a',  # green
    81: '#d4a030',  # amber
    82: '#d4622a',  # orange-red
    83: '#c0392b',  # deep red
}

# ============================================================================
# DATA LOADING — reusing logic from coastsat_lrr_analysis.py
# ============================================================================

def collect_csv_map(root_dir, site_filter=""):
    """Auto-discover all time-series CSVs under root_dir."""
    csv_map = {}
    if not os.path.isdir(root_dir):
        print(f"WARNING: ROOT_DATA_DIR not found: {root_dir}")
        return csv_map
    subfolders = [
        os.path.join(root_dir, d)
        for d in sorted(os.listdir(root_dir))
        if os.path.isdir(os.path.join(root_dir, d))
        and (site_filter == "" or site_filter in d)
    ]
    for sf in subfolders:
        for fpath in glob.glob(os.path.join(sf, "*.csv")):
            stem = os.path.splitext(os.path.basename(fpath))[0]
            csv_map[stem] = fpath
    print(f"Found {len(csv_map)} total time-series CSVs.")
    return csv_map


def load_timeseries(filepath):
    """Load a CoastSat transect CSV → DataFrame with [date, chainage_m]."""
    df = pd.read_csv(filepath, header=0)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(columns={df.columns[0]: "date", df.columns[1]: "chainage_m"})
    df["date"] = pd.to_datetime(df["date"], utc=True)
    df["chainage_m"] = pd.to_numeric(df["chainage_m"], errors="coerce")
    df = df.dropna(subset=["date", "chainage_m"]).sort_values("date").reset_index(drop=True)
    return df


def filter_dates(df, start, end):
    """Clip to [start, end] date range."""
    if start:
        df = df[df["date"] >= pd.Timestamp(start, tz="UTC")]
    if end:
        df = df[df["date"] <= pd.Timestamp(end, tz="UTC")]
    return df.reset_index(drop=True)


def load_rodanthe_data(root_dir, lookup_csv, domain_min, domain_max,
                       site_filter, start_date, end_date, min_obs):
    """
    Load all CoastSat transect time series for the specified domain range.
    Returns a long-format DataFrame: [date, chainage_m, transect_id, domain_number]
    """
    # Load lookup table
    lookup = pd.read_csv(lookup_csv)
    # Handle corrupted header (filename accidentally prepended)
    lookup.columns = [c.split('csv')[-1] if 'csv' in c else c for c in lookup.columns]

    # Filter to target domains
    mask = (lookup["domain_number"] >= domain_min) & (lookup["domain_number"] <= domain_max)
    lookup_sub = lookup[mask].copy()
    print(f"Domains {domain_min}–{domain_max}: {len(lookup_sub)} transects in lookup.")

    # Collect CSV map
    csv_map = collect_csv_map(root_dir, site_filter)

    # Load each transect
    records = []
    for _, row in lookup_sub.iterrows():
        tid = str(row["transect_id"])
        domain = int(row["domain_number"])
        if tid not in csv_map:
            continue
        try:
            df_raw = load_timeseries(csv_map[tid])
            df_filt = filter_dates(df_raw, start_date, end_date)
            if len(df_filt) < min_obs:
                continue
            df_filt = df_filt.copy()
            df_filt["transect_id"] = tid
            df_filt["domain_number"] = domain
            records.append(df_filt)
        except Exception as e:
            print(f"  ERROR loading {tid}: {e}")

    if not records:
        raise RuntimeError("No transect data loaded. Check your ROOT_DATA_DIR and LOOKUP_CSV paths.")

    df_all = pd.concat(records, ignore_index=True)
    df_all["year"] = df_all["date"].dt.year
    df_all["decimal_year"] = df_all["year"] + df_all["date"].dt.dayofyear / 365.25
    print(f"Loaded {len(df_all)} observations from {df_all['transect_id'].nunique()} transects "
          f"across {df_all['domain_number'].nunique()} domains.")
    return df_all


# ============================================================================
# FIGURE 1: SPACE-TIME HEATMAP
# ============================================================================

def plot_heatmap(df_all, domain_min, domain_max, output_dir, fig_size):
    """
    Space-time heatmap: X = time, Y = domain, color = mean annual
    shoreline position anomaly (relative to each domain's long-term mean).
    Each cell is the mean chainage across all transects in that domain for that year.
    """
    print("\nBuilding Figure 1: Space-time heatmap...")

    # Aggregate to domain × year means
    domain_year = (
        df_all.groupby(["domain_number", "year"])["chainage_m"]
        .mean()
        .reset_index()
    )

    # Remove anomaly relative to each domain's long-term mean
    # so spatial differences in absolute position don't dominate
    domain_means = domain_year.groupby("domain_number")["chainage_m"].mean().rename("domain_mean")
    domain_year = domain_year.merge(domain_means, on="domain_number")
    domain_year["anomaly"] = domain_year["chainage_m"] - domain_year["domain_mean"]

    # Pivot to matrix: rows = domains (sorted S→N), cols = years
    pivot = domain_year.pivot(index="domain_number", columns="year", values="anomaly")
    pivot = pivot.sort_index(ascending=False)  # south at top → north at bottom (Hatteras convention)

    years = pivot.columns.values
    domains = pivot.index.values

    fig, ax = plt.subplots(figsize=fig_size)
    fig.patch.set_facecolor('white')

    # Diverging colormap centered on zero
    vmax = np.nanpercentile(np.abs(pivot.values), 95)
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu_r

    im = ax.pcolormesh(years, np.arange(len(domains)), pivot.values,
                       cmap=cmap, norm=norm, shading='auto')

    # Y-axis: domain labels
    ax.set_yticks(np.arange(len(domains)) + 0.5)
    ax.set_yticklabels([f"Domain {d}" for d in domains], fontsize=9)

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Shoreline position anomaly (m)", fontsize=10)
    cbar.ax.tick_params(labelsize=9)

    # Period shading if both dates set
    ax.set_xlim(years.min(), years.max())
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='x', which='major', labelsize=10)
    ax.tick_params(axis='both', which='minor', length=3)
    ax.grid(True, axis='x', which='major', linestyle=':', linewidth=0.5, alpha=0.4)

    ax.set_xlabel("Year", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel("CASCADE Domain", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_title(f"Shoreline Position Anomaly — Rodanthe (Domains {domain_min}–{domain_max})\n"
                 f"CoastSat satellite-derived shorelines",
                 fontsize=12, fontweight='bold', color='#1a2a3a', pad=10)

    ax.spines[['top', 'right']].set_visible(False)

    fig.text(0.01, 0.005,
             'Annual mean cross-shore position anomaly relative to each domain\'s long-term mean. '
             'Red = seaward (accretion), Blue = landward (erosion).',
             fontsize=7.5, color='#666666', style='italic', ha='left', va='bottom')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = os.path.join(output_dir, "rodanthe_heatmap.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out_path}")
    plt.show()
    plt.close()


# ============================================================================
# FIGURE 2: SPAGHETTI PLOT (raw transect time series)
# ============================================================================

def plot_spaghetti(df_all, domain_min, domain_max, domain_colors, output_dir, fig_size):
    """
    Raw cross-shore position vs. time for every transect, colored by domain.
    Each transect is normalized to its own mean so all lines are on the same scale.
    """
    print("\nBuilding Figure 2: Spaghetti plot...")

    # Normalize each transect to its own mean (removes spatial offset between transects)
    df_plot = df_all.copy()
    transect_means = df_plot.groupby("transect_id")["chainage_m"].transform("mean")
    df_plot["chainage_norm"] = df_plot["chainage_m"] - transect_means

    fig, ax = plt.subplots(figsize=fig_size)
    fig.patch.set_facecolor('white')

    domains_present = sorted(df_plot["domain_number"].unique())

    for domain in domains_present:
        color = domain_colors.get(domain, '#888888')
        transects = df_plot[df_plot["domain_number"] == domain]["transect_id"].unique()
        for tid in transects:
            t_df = df_plot[df_plot["transect_id"] == tid].sort_values("decimal_year")
            ax.plot(t_df["decimal_year"], t_df["chainage_norm"],
                    color=color, linewidth=0.7, alpha=0.55, zorder=2)

    # Zero line
    ax.axhline(0, color='#2c2c2c', linewidth=1.0, linestyle='--', alpha=0.6, zorder=3)

    # Legend — one entry per domain
    legend_elements = [
        Line2D([0], [0], color=domain_colors.get(d, '#888888'), linewidth=2.0,
               label=f'Domain {d}')
        for d in domains_present
    ]
    ax.legend(handles=legend_elements, fontsize=9, framealpha=0.95,
              edgecolor='#cccccc', loc='upper left')

    ax.set_xlim(df_plot["decimal_year"].min() - 0.5,
                df_plot["decimal_year"].max() + 0.5)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='both', which='major', labelsize=10, direction='in', length=5)
    ax.tick_params(axis='both', which='minor', direction='in', length=3)
    ax.grid(True, which='major', linestyle=':', linewidth=0.5, alpha=0.4, color='gray')
    ax.spines[['top', 'right']].set_visible(False)

    ax.set_xlabel("Year", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel("Cross-shore position anomaly (m)", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_title(f"CoastSat Shoreline Time Series — Rodanthe (Domains {domain_min}–{domain_max})",
                 fontsize=12, fontweight='bold', color='#1a2a3a', pad=10)

    fig.text(0.01, 0.005,
             'Each line is one CoastSat transect, normalized to its own long-term mean. '
             'Positive = seaward (accretion); negative = landward (erosion).',
             fontsize=7.5, color='#666666', style='italic', ha='left', va='bottom')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = os.path.join(output_dir, "rodanthe_spaghetti.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out_path}")
    plt.show()
    plt.close()


# ============================================================================
# FIGURE 3: ANNUAL DOMAIN MEAN SHORELINE POSITION
# ============================================================================

def plot_domain_annual_means(df_all, domain_min, domain_max, domain_colors, output_dir, fig_size):
    """
    Annual mean shoreline position anomaly per domain — cleaner than the
    spaghetti plot, one bold line per domain. Good for poster.
    """
    print("\nBuilding Figure 3: Annual domain mean shoreline positions...")

    # Annual mean per domain, normalized to each domain's long-term mean
    domain_year = (
        df_all.groupby(["domain_number", "year"])["chainage_m"]
        .mean()
        .reset_index()
    )
    domain_means = domain_year.groupby("domain_number")["chainage_m"].mean().rename("domain_mean")
    domain_year = domain_year.merge(domain_means, on="domain_number")
    domain_year["anomaly"] = domain_year["chainage_m"] - domain_year["domain_mean"]

    fig, ax = plt.subplots(figsize=fig_size)
    fig.patch.set_facecolor('white')

    domains_present = sorted(domain_year["domain_number"].unique())

    for domain in domains_present:
        color = domain_colors.get(domain, '#888888')
        d_df = domain_year[domain_year["domain_number"] == domain].sort_values("year")
        ax.plot(d_df["year"], d_df["anomaly"],
                color=color, linewidth=2.2, zorder=3, label=f'Domain {domain}')
        ax.scatter(d_df["year"], d_df["anomaly"],
                   color=color, s=20, zorder=4, alpha=0.8)

    # Zero line
    ax.axhline(0, color='#2c2c2c', linewidth=1.0, linestyle='--', alpha=0.6, zorder=2)

    # Shaded fill for each domain
    for domain in domains_present:
        color = domain_colors.get(domain, '#888888')
        d_df = domain_year[domain_year["domain_number"] == domain].sort_values("year")
        ax.fill_between(d_df["year"], d_df["anomaly"], 0,
                        where=(d_df["anomaly"] >= 0), alpha=0.04, color=color, interpolate=True)
        ax.fill_between(d_df["year"], d_df["anomaly"], 0,
                        where=(d_df["anomaly"] < 0), alpha=0.07, color=color, interpolate=True)

    ax.legend(fontsize=9, framealpha=0.95, edgecolor='#cccccc',
              loc='upper left', ncol=2)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(axis='both', which='major', labelsize=10, direction='in', length=5)
    ax.tick_params(axis='both', which='minor', direction='in', length=3)
    ax.grid(True, which='major', linestyle=':', linewidth=0.5, alpha=0.4, color='gray')
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_linewidth(1.2)

    # Erosion / accretion labels
    ybot, ytop = ax.get_ylim()
    zero_frac = (0 - ybot) / (ytop - ybot)
    ax.text(1.0, zero_frac + (1 - zero_frac) / 2, 'Accretion ▲',
            transform=ax.transAxes, fontsize=9, color='#555555',
            ha='right', va='center', style='italic')
    ax.text(1.0, zero_frac / 2, 'Erosion ▼',
            transform=ax.transAxes, fontsize=9, color='#555555',
            ha='right', va='center', style='italic')

    ax.set_xlabel("Year", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_ylabel("Shoreline position anomaly (m)", fontsize=12, fontweight='bold', labelpad=8)
    ax.set_title(f"Annual Mean Shoreline Position — Rodanthe (Domains {domain_min}–{domain_max})\n"
                 f"CoastSat satellite-derived shorelines",
                 fontsize=12, fontweight='bold', color='#1a2a3a', pad=10)

    fig.text(0.01, 0.005,
             'Annual mean cross-shore position per domain, normalized to each domain\'s long-term mean. '
             'Shoreline data: CoastSat (Vos et al., 2019).',
             fontsize=7.5, color='#666666', style='italic', ha='left', va='bottom')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    out_path = os.path.join(output_dir, "rodanthe_annual_domain_means.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"  Saved: {out_path}")
    plt.show()
    plt.close()


# ============================================================================
# MAIN
# ============================================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print(f"CoastSat Rodanthe Visualization (Domains {DOMAIN_MIN}–{DOMAIN_MAX})")
    print(f"Period: {START_DATE} → {END_DATE}")
    print("=" * 60)

    # Load all transect data for the Rodanthe domains
    df_all = load_rodanthe_data(
        ROOT_DATA_DIR, LOOKUP_CSV, DOMAIN_MIN, DOMAIN_MAX,
        SITE_FILTER, START_DATE, END_DATE, MIN_OBS
    )

    # Produce all three figures
    plot_heatmap(df_all, DOMAIN_MIN, DOMAIN_MAX, OUTPUT_DIR, FIG_SIZE)
    plot_spaghetti(df_all, DOMAIN_MIN, DOMAIN_MAX, DOMAIN_COLORS, OUTPUT_DIR, FIG_SIZE)
    plot_domain_annual_means(df_all, DOMAIN_MIN, DOMAIN_MAX, DOMAIN_COLORS, OUTPUT_DIR, FIG_SIZE)

    print("\nAll figures saved to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
