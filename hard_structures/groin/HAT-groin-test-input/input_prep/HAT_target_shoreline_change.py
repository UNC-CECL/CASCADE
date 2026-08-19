"""
HAT_target_shoreline_change.py
==============================
Build the TARGET the groin module must reproduce: the observed dune-line
POSITION CHANGE (and rate) across years, per domain, from the raw ArcGIS offset
files. This is what you compare the CASCADE groin runs against.

Why this exists
---------------
The groin's job is to bend the modeled shoreline so it matches the historical
1967->1997 differential (holding updrift, eroding downdrift). To calibrate M you
need that historical signal as a per-domain curve. This script produces it.

Key difference from island_offset_hybrid_1967.py
------------------------------------------------
That pipeline references EACH year to its OWN minimum -- correct for a single
CASCADE input file, but it destroys cross-year comparability (every year gets a
different zero). To compare POSITIONS across years, all years must share ONE
baseline. This script:
  1. extracts per-domain mean ORIG_LEN (raw cross-shore position, m) per year,
  2. keeps them on a SHARED raw reference (no per-year re-zeroing),
  3. computes position change between any two years, and the rate (m/yr),
  4. optionally re-references the whole set to one year (e.g. 1967) or one
     domain, purely for display -- differences/rates are reference-invariant.

Sign convention
---------------
ORIG_LEN increases LANDWARD (larger = more retreat), same as CASCADE x_s. So a
POSITIVE position change = landward = EROSION. With FLIP_FOR_RATE=True the
reported RATE is flipped to + = seaward/accretion, matching your hindcast plots.

Output
------
  <OUT>_positions_by_year.csv   per-domain position each year (shared raw ref)
  <OUT>_change_<A>_<B>.csv      position change + rate between year A and B
  <OUT>_target_curve.png        the target rate curve (what the groin must match)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================== CONFIG ==============================

# Raw ArcGIS dune-line offset files, keyed by year.
RAW_DIR = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\raw_offsets"
RAW_FILES = {
    1967: os.path.join(RAW_DIR, "1967_duneline_offset_raw.csv"),
    1978: os.path.join(RAW_DIR, "1978_duneline_offset_raw.csv"),
    1984: os.path.join(RAW_DIR, "1984_duneline_offset_raw.csv"),
    1997: os.path.join(RAW_DIR, "1997_duneline_offset_raw.csv"),
}

# Column names in the raw files (match island_offset_hybrid_1967.py COL_MAP).
DOMAIN_COL   = "domain_id"
POSITION_COL = "ORIG_LEN"     # cross-shore position (m); increases landward

# Domain window to compare. D2-D12 is the range with 1967 coverage.
START_DOMAIN = 2
END_DOMAIN   = 12
DOMAINS = list(range(START_DOMAIN, END_DOMAIN + 1))

# The two years whose change is the primary groin TARGET.
CHANGE_FROM = 1967
CHANGE_TO   = 1997

# Display reference (does NOT affect change/rate, only the positions plot):
#   "none"        -> keep shared raw ORIG_LEN
#   ("year", Y)   -> subtract year Y's positions (so Y becomes the zero line)
#   ("domain", D) -> subtract each year's value at domain D
DISPLAY_REFERENCE = ("year", 1967)

# Report rate as + = seaward/accretion (matches your hindcast plots).
FLIP_FOR_RATE = True

# --- Time-axis trajectory panel + GIF ---
MAKE_TRAJECTORY_PANEL = True    # position-vs-year, one line per domain
MAKE_GIF              = True    # proportional-duration animation of the snapshots
GIF_SECONDS_PER_YEAR  = 0.18    # each real year of gap -> this many seconds of hold
GIF_MIN_HOLD_S        = 1.2     # floor so even short gaps are readable
GIF_FPS               = 10      # smoothness of the hold (frames are static repeats)

# Groin annotation.
GROIN_BOUNDARY = 5.5
GROIN_COLOR    = "#B71C1C"

OUT_DIR      = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\groin_module\hindcast_groin_test\groin_init\target"
OUT_BASENAME = "HAT_target_1967_1997"


# ============================ EXTRACTION ============================

def per_domain_positions(path):
    """Mean POSITION_COL per domain (raw, shared reference). Series indexed by domain."""
    df = pd.read_csv(path, encoding="utf-8-sig")
    for c in (DOMAIN_COL, POSITION_COL):
        if c not in df.columns:
            raise KeyError(f"'{c}' not in {os.path.basename(path)}")
    df = df[df[DOMAIN_COL].isin(DOMAINS)]
    return df.groupby(DOMAIN_COL)[POSITION_COL].mean()


def build_position_table():
    """DataFrame: index = domain (DOMAINS), columns = years, values = raw position (m)."""
    cols = {}
    for yr, path in RAW_FILES.items():
        if not os.path.isfile(path):
            print(f"  [WARN] missing raw file for {yr}: {path}")
            continue
        cols[yr] = per_domain_positions(path)
    tbl = pd.DataFrame(cols).reindex(DOMAINS)
    tbl.index.name = "domain_id"
    return tbl.sort_index(axis=1)


def apply_display_reference(tbl):
    """Return a copy shifted per DISPLAY_REFERENCE (for the positions plot only)."""
    ref = DISPLAY_REFERENCE
    if ref == "none" or ref is None:
        return tbl.copy()
    kind, val = ref
    if kind == "year":
        if val not in tbl.columns:
            print(f"  [WARN] display ref year {val} absent; showing raw positions.")
            return tbl.copy()
        return tbl.sub(tbl[val], axis=0)
    if kind == "domain":
        if val not in tbl.index:
            print(f"  [WARN] display ref domain {val} absent; showing raw positions.")
            return tbl.copy()
        return tbl.sub(tbl.loc[val], axis=1)
    raise ValueError(f"bad DISPLAY_REFERENCE: {ref!r}")


def compute_change(tbl, yr_from, yr_to):
    """Per-domain position change and rate between two years. Reference-invariant."""
    if yr_from not in tbl.columns or yr_to not in tbl.columns:
        raise ValueError(f"need both {yr_from} and {yr_to} in the position table.")
    span = yr_to - yr_from
    change_m = tbl[yr_to] - tbl[yr_from]              # + = landward = erosion
    rate = change_m / span
    if FLIP_FOR_RATE:
        rate = -rate                                  # + = seaward/accretion
    out = pd.DataFrame({
        "domain_id": tbl.index,
        f"pos_{yr_from}_m": tbl[yr_from].values,
        f"pos_{yr_to}_m": tbl[yr_to].values,
        "position_change_m": change_m.values,         # raw sign: + = erosion
        "rate_m_per_yr": rate.values,                 # flipped: + = accretion
    })
    return out


# ============================== PLOT ==============================

def plot_target(tbl_display, change_df):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 9), constrained_layout=True)

    # -- top: positions over time (display-referenced) --
    for yr in tbl_display.columns:
        ax1.plot(tbl_display.index, tbl_display[yr], marker="o", ms=4, lw=1.6,
                 label=str(yr))
    ax1.axvline(GROIN_BOUNDARY, color=GROIN_COLOR, ls="--", lw=1.5, alpha=0.9)
    ax1.text(GROIN_BOUNDARY, ax1.get_ylim()[1], " Buxton groin", color=GROIN_COLOR,
             fontsize=8, rotation=90, va="top", ha="left")
    ref = DISPLAY_REFERENCE
    ref_txt = (f"referenced to {ref[1]}" if isinstance(ref, tuple) else "raw ORIG_LEN")
    ax1.set_ylabel(f"Dune-line position (m)\n[{ref_txt}; + = landward]")
    ax1.set_title(f"Dune-line position by year  D{START_DOMAIN}-D{END_DOMAIN}")
    ax1.grid(alpha=0.3); ax1.legend(title="Year", fontsize=8)
    ax1.invert_yaxis()  # landward down, so the plot reads like a map (sea at top)

    # -- bottom: the TARGET rate curve --
    ax2.plot(change_df["domain_id"], change_df["rate_m_per_yr"],
             marker="o", ms=5, lw=2, color="#08519C",
             label=f"Observed {CHANGE_FROM}-{CHANGE_TO}")
    ax2.axhline(0, color="gray", ls="--", lw=1, alpha=0.7)
    ax2.axvline(GROIN_BOUNDARY, color=GROIN_COLOR, ls="--", lw=1.5, alpha=0.9)
    ax2.axvspan(START_DOMAIN - 0.5, GROIN_BOUNDARY, alpha=0.06, color="firebrick")
    ax2.axvspan(GROIN_BOUNDARY, END_DOMAIN + 0.5, alpha=0.06, color="seagreen")
    ax2.text(GROIN_BOUNDARY, ax2.get_ylim()[1], " groin", color=GROIN_COLOR,
             fontsize=8, rotation=90, va="top", ha="left")
    ax2.set_xlabel(f"GIS Domain ID (D{START_DOMAIN}-D{END_DOMAIN})")
    ax2.set_ylabel("Target shoreline change rate (m/yr)\n[+ = seaward]")
    ax2.set_title(f"TARGET for groin calibration  ({CHANGE_FROM}-{CHANGE_TO})  "
                  f"| downdrift amplified (red), updrift suppressed (green)")
    ax2.grid(alpha=0.3); ax2.legend(fontsize=8)
    return fig


# ========================= TRAJECTORY + GIF =========================

def plot_trajectories(tbl_display):
    """Position vs YEAR, one line per domain -- the honest 'over time' view.
    Years sit at their true spacing on the x-axis, so uneven gaps show correctly.
    Updrift validation domains (D6-D10) are bold/saturated; downdrift (D2-D5,
    not validated) and far-updrift edge (D11-D12) are muted."""
    years = list(tbl_display.columns)
    fig, ax = plt.subplots(figsize=(11, 6), constrained_layout=True)

    # Distinct saturated colors for the updrift validation zone (D6-D10);
    # muted grays for downdrift (not validated) and the D11-D12 edge.
    updrift_focus = [6, 7, 8, 9, 10]
    focus_colors = plt.cm.autumn(np.linspace(0.0, 0.75, len(updrift_focus)))
    focus_map = dict(zip(updrift_focus, focus_colors))

    for dom in tbl_display.index:
        traj = tbl_display.loc[dom].values
        if dom in focus_map:
            ax.plot(years, traj, marker="o", ms=5, lw=2.4, color=focus_map[dom],
                    zorder=5, label=f"D{dom} (updrift*)")
        else:
            side = "downdrift" if dom <= GROIN_BOUNDARY else "updrift edge"
            ax.plot(years, traj, marker="o", ms=3, lw=1.2, color="0.65",
                    alpha=0.7, zorder=2,
                    label=f"D{dom} ({side})" if dom in (2, 11) else None)

    ax.axhline(0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.set_xlabel("Year")
    ref = DISPLAY_REFERENCE
    ref_txt = (f"referenced to {ref[1]}" if isinstance(ref, tuple) else "raw ORIG_LEN")
    ax.set_ylabel(f"Dune-line position (m)\n[{ref_txt}; + = landward]")
    ax.set_title("Dune-line position over time, by domain   "
                 "(downdrift amplified, updrift suppressed = groin signature)")
    ax.invert_yaxis()   # landward down -> lines fall as the shore retreats
    ax.set_xticks(years)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2, title="* = updrift validation zone")
    return fig


def make_gif(tbl_display, out_path):
    """Proportional-duration GIF: each snapshot is held for a time proportional
    to the real gap to the NEXT snapshot, so uneven year spacing is respected.
    Falls back gracefully (prints a note) if pillow isn't available."""
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except Exception as e:
        print(f"  [GIF skipped] matplotlib animation unavailable: {e}")
        return False

    years = list(tbl_display.columns)
    doms  = list(tbl_display.index)

    # Build the frame schedule: repeat each year's frame in proportion to the
    # gap to the next year (last year gets the same hold as the previous gap).
    gaps = [years[i + 1] - years[i] for i in range(len(years) - 1)]
    gaps.append(gaps[-1] if gaps else 1)
    holds_s = [max(GIF_MIN_HOLD_S, g * GIF_SECONDS_PER_YEAR) for g in gaps]
    frames_per = [max(1, int(round(h * GIF_FPS))) for h in holds_s]
    frame_years = []
    for yr, n in zip(years, frames_per):
        frame_years.extend([yr] * n)

    # Fixed axes so the animation doesn't jump.
    all_vals = tbl_display.values
    ymin, ymax = np.nanmin(all_vals), np.nanmax(all_vals)
    pad = 0.08 * (ymax - ymin if ymax > ymin else 1)

    fig, ax = plt.subplots(figsize=(9, 5.5))

    def draw(frame_idx):
        ax.clear()
        yr = frame_years[frame_idx]
        # faded trails for all earlier snapshots
        for prev in years:
            if prev < yr:
                ax.plot(doms, tbl_display[prev].values, marker="o", ms=3,
                        lw=1, color="0.75", alpha=0.5, zorder=1)
        ax.plot(doms, tbl_display[yr].values, marker="o", ms=6, lw=2.4,
                color="#08519C", zorder=3)
        ax.axvline(GROIN_BOUNDARY, color=GROIN_COLOR, ls="--", lw=1.5, alpha=0.9)
        ax.axvspan(min(doms) - 0.5, GROIN_BOUNDARY, alpha=0.06, color="firebrick")
        ax.axvspan(GROIN_BOUNDARY, max(doms) + 0.5, alpha=0.06, color="seagreen")
        ax.set_ylim(ymax + pad, ymin - pad)   # inverted: landward down
        ax.set_xlim(min(doms) - 0.5, max(doms) + 0.5)
        ax.set_xlabel(f"GIS Domain ID (D{min(doms)}-D{max(doms)})")
        ref = DISPLAY_REFERENCE
        ref_txt = (f"ref {ref[1]}" if isinstance(ref, tuple) else "raw")
        ax.set_ylabel(f"Dune-line position (m) [{ref_txt}; + landward]")
        ax.set_title(f"Dune-line position  \u2014  {yr}")
        ax.text(0.02, 0.06, str(yr), transform=ax.transAxes, fontsize=28,
                fontweight="bold", color="#08519C", alpha=0.85, va="bottom")
        ax.grid(alpha=0.3)

    anim = FuncAnimation(fig, draw, frames=len(frame_years), interval=1000 / GIF_FPS)
    try:
        anim.save(out_path, writer=PillowWriter(fps=GIF_FPS))
        plt.close(fig)
        print(f"  Saved GIF ({len(years)} snapshots, "
              f"holds {[round(h,1) for h in holds_s]} s): {out_path}")
        return True
    except Exception as e:
        plt.close(fig)
        print(f"  [GIF skipped] could not write: {e}")
        return False


# ========================= QUANTIFICATION =========================

def quantify(change_df):
    """Reduce the per-domain target to defensible scalar metrics + a fold-and-sum
    decomposition (A = groin signal, B = background) across the groin.

    Returns (zone_summary_dict, foldsum_dataframe).
    """
    rate = dict(zip(change_df["domain_id"], change_df["rate_m_per_yr"]))  # + = seaward
    pos  = dict(zip(change_df["domain_id"], change_df["position_change_m"]))  # + = landward

    up_doms = [d for d in DOMAINS if d > GROIN_BOUNDARY]
    dn_doms = [d for d in DOMAINS if d < GROIN_BOUNDARY]
    up_rate = [rate[d] for d in up_doms]
    dn_rate = [rate[d] for d in dn_doms]

    zone = {
        "updrift_domains":   f"D{up_doms[0]}-D{up_doms[-1]}",
        "downdrift_domains": f"D{dn_doms[0]}-D{dn_doms[-1]}",
        "updrift_mean_rate_m_yr":   float(np.mean(up_rate)),
        "downdrift_mean_rate_m_yr": float(np.mean(dn_rate)),
        "updrift_rate_at_groin_D%d" % up_doms[0]: float(rate[up_doms[0]]),
        "downdrift_peak_erosion_m_yr": float(min(dn_rate)),
        "downdrift_peak_domain": f"D{dn_doms[int(np.argmin(dn_rate))]}",
        "erosion_differential_ratio": float(np.mean(dn_rate) / np.mean(up_rate))
                                      if np.mean(up_rate) != 0 else float("nan"),
        "updrift_total_change_m":   float(np.sum([pos[d] for d in up_doms])),
        "downdrift_total_change_m": float(np.sum([pos[d] for d in dn_doms])),
    }

    # Fold-and-sum: pair updrift domain u with downdrift partner (mirror across groin).
    rows = []
    for u in up_doms:
        p = int(round(2 * GROIN_BOUNDARY - u))   # mirror domain
        if p in rate:
            A = (rate[u] - rate[p]) / 2.0        # groin's antisymmetric signal
            B = (rate[u] + rate[p]) / 2.0        # background (groin removed)
            rows.append({
                "dist_from_groin": u - GROIN_BOUNDARY,
                "updrift_domain": u, "downdrift_domain": p,
                "updrift_rate": rate[u], "downdrift_rate": rate[p],
                "A_groin_signal_m_yr": A,
                "B_background_m_yr": B,
            })
    foldsum = pd.DataFrame(rows)

    # Integrated groin signal (area under A vs distance) -- the number M reproduces.
    if not foldsum.empty:
        zone["integrated_A_groin_signal_m_yr_per_domain"] = float(foldsum["A_groin_signal_m_yr"].sum())
        zone["A_peak_m_yr"] = float(foldsum["A_groin_signal_m_yr"].max())
        zone["B_mean_m_yr"] = float(foldsum["B_background_m_yr"].mean())
        zone["B_range_m_yr"] = float(foldsum["B_background_m_yr"].max()
                                     - foldsum["B_background_m_yr"].min())
    return zone, foldsum


def plot_foldsum(foldsum):
    """Two-line figure: A(y) = groin signal (what M reproduces), B(y) = background."""
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    d = foldsum["dist_from_groin"]
    ax.plot(d, foldsum["A_groin_signal_m_yr"], marker="o", ms=6, lw=2.2,
            color="#08519C", label="A = groin signal (antisymmetric)")
    ax.plot(d, foldsum["B_background_m_yr"], marker="s", ms=6, lw=2.2,
            color="#B71C1C", ls="--", label="B = background (symmetric)")
    ax.axhline(0, color="gray", ls=":", lw=1, alpha=0.7)
    ax.set_xlabel("Distance from groin (domains)")
    ax.set_ylabel("Rate (m/yr)")
    ax.set_title("Fold-and-sum decomposition   |   A is the target your M reproduces\n"
                 "(flat B would confirm uniform-background dipole model)")
    # annotate the pairing
    for _, r in foldsum.iterrows():
        ax.annotate(f"D{int(r['updrift_domain'])}/D{int(r['downdrift_domain'])}",
                    (r["dist_from_groin"], r["A_groin_signal_m_yr"]),
                    textcoords="offset points", xytext=(0, 8),
                    fontsize=7, ha="center", color="#08519C")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    return fig


def plot_quantification_dashboard(change_df, zone, foldsum):
    """One-glance summary of the quantified target: per-domain rate bars colored
    by zone, a downdrift-vs-updrift mean comparison with the differential ratio,
    and the A/B fold-sum decomposition."""
    rate = dict(zip(change_df["domain_id"], change_df["rate_m_per_yr"]))
    doms = list(change_df["domain_id"])

    fig = plt.figure(figsize=(13, 8), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0])
    ax_bars = fig.add_subplot(gs[0, :])     # top: per-domain rate bars
    ax_mean = fig.add_subplot(gs[1, 0])     # bottom-left: zone means
    ax_fold = fig.add_subplot(gs[1, 1])     # bottom-right: A/B decomposition

    # ── TOP: per-domain rate bars, colored by zone ──
    colors = ["#B71C1C" if d < GROIN_BOUNDARY else "#08519C" for d in doms]
    bars = ax_bars.bar(doms, [rate[d] for d in doms], color=colors,
                       edgecolor="black", linewidth=0.5, width=0.7)
    ax_bars.axhline(0, color="gray", lw=1)
    ax_bars.axvline(GROIN_BOUNDARY, color=GROIN_COLOR, ls="--", lw=1.6, zorder=5)
    ax_bars.text(GROIN_BOUNDARY, ax_bars.get_ylim()[1] * 0.9, " groin",
                 color=GROIN_COLOR, fontsize=9, rotation=90, va="top")
    for d in doms:
        ax_bars.annotate(f"{rate[d]:+.1f}", (d, rate[d]),
                         textcoords="offset points",
                         xytext=(0, 4 if rate[d] >= 0 else -12),
                         ha="center", fontsize=7)
    ax_bars.set_xticks(doms)
    ax_bars.set_xlabel(f"GIS Domain ID (D{DOMAINS[0]}-D{DOMAINS[-1]})")
    ax_bars.set_ylabel("Rate (m/yr)  [+ seaward]")
    ax_bars.set_title(f"Per-domain shoreline change rate  {CHANGE_FROM}-{CHANGE_TO}   "
                      f"(red = downdrift / eroding, blue = updrift)")
    ax_bars.grid(alpha=0.3, axis="y")

    # ── BOTTOM-LEFT: zone means + differential ──
    means = [zone["downdrift_mean_rate_m_yr"], zone["updrift_mean_rate_m_yr"]]
    labels = [f"Downdrift\n{zone['downdrift_domains']}", f"Updrift\n{zone['updrift_domains']}"]
    b = ax_mean.bar(labels, means, color=["#B71C1C", "#08519C"],
                    edgecolor="black", linewidth=0.5, width=0.6)
    ax_mean.axhline(0, color="gray", lw=1)
    for rect, v in zip(b, means):
        ax_mean.annotate(f"{v:+.2f}", (rect.get_x() + rect.get_width() / 2, v),
                         textcoords="offset points",
                         xytext=(0, 6 if v >= 0 else -14), ha="center",
                         fontsize=10, fontweight="bold")
    ax_mean.set_ylabel("Mean rate (m/yr)")
    ax_mean.set_title(f"Zone means  |  downdrift erodes "
                      f"{zone['erosion_differential_ratio']:.1f}\u00d7 the updrift")
    ax_mean.grid(alpha=0.3, axis="y")

    # ── BOTTOM-RIGHT: A/B decomposition ──
    if not foldsum.empty:
        d = foldsum["dist_from_groin"]
        ax_fold.plot(d, foldsum["A_groin_signal_m_yr"], marker="o", ms=6, lw=2.2,
                     color="#08519C", label="A = groin signal")
        ax_fold.plot(d, foldsum["B_background_m_yr"], marker="s", ms=6, lw=2.2,
                     color="#B71C1C", ls="--", label="B = background")
        ax_fold.axhline(0, color="gray", ls=":", lw=1)
        for _, r in foldsum.iterrows():
            ax_fold.annotate(f"D{int(r['updrift_domain'])}/D{int(r['downdrift_domain'])}",
                             (r["dist_from_groin"], r["A_groin_signal_m_yr"]),
                             textcoords="offset points", xytext=(0, 8),
                             fontsize=7, ha="center", color="#08519C")
        b_flat = zone.get("B_range_m_yr", 0) < 1.5
        ax_fold.set_title(f"Fold-sum: A=groin, B=background\n"
                          f"(B {'flat -> uniform model OK' if b_flat else 'trends -> varies alongshore'})",
                          fontsize=10)
        ax_fold.set_xlabel("Distance from groin (domains)")
        ax_fold.set_ylabel("Rate (m/yr)")
        ax_fold.grid(alpha=0.3); ax_fold.legend(fontsize=8)
    else:
        ax_fold.text(0.5, 0.5, "no paired domains", ha="center", va="center")
        ax_fold.axis("off")

    return fig


# ============================== MAIN ==============================

def main():
    print("=" * 70)
    print("Building groin calibration target from dune-line positions")
    print("=" * 70)

    tbl = build_position_table()
    if tbl.shape[1] < 2:
        print("Need at least two years of data. Check RAW_FILES paths.")
        return
    print(f"\nYears loaded: {list(tbl.columns)}")
    print(f"Domains: D{START_DOMAIN}-D{END_DOMAIN}\n")

    change_df = compute_change(tbl, CHANGE_FROM, CHANGE_TO)
    tbl_display = apply_display_reference(tbl)

    print(f"Target change {CHANGE_FROM}->{CHANGE_TO}:")
    print(change_df[["domain_id", "position_change_m", "rate_m_per_yr"]]
          .to_string(index=False))

    # ── Quantify ─────────────────────────────────────────────────────────────
    zone, foldsum = quantify(change_df)
    print("\n" + "-" * 60)
    print("QUANTIFIED TARGET")
    print("-" * 60)
    print(f"  downdrift {zone['downdrift_domains']:<8} mean rate "
          f"{zone['downdrift_mean_rate_m_yr']:+.2f} m/yr  "
          f"(peak {zone['downdrift_peak_erosion_m_yr']:+.2f} at {zone['downdrift_peak_domain']})")
    print(f"  updrift   {zone['updrift_domains']:<8} mean rate "
          f"{zone['updrift_mean_rate_m_yr']:+.2f} m/yr")
    print(f"  erosion differential: downdrift erodes "
          f"{zone['erosion_differential_ratio']:.1f}x the updrift mean")
    if not foldsum.empty:
        print(f"\n  Fold-and-sum (A = groin signal, B = background):")
        print(foldsum[["dist_from_groin", "updrift_domain", "downdrift_domain",
                       "A_groin_signal_m_yr", "B_background_m_yr"]].to_string(index=False))
        print(f"\n  A peak = {zone['A_peak_m_yr']:+.2f} m/yr   "
              f"integrated A = {zone['integrated_A_groin_signal_m_yr_per_domain']:+.2f} m/yr/domain")
        print(f"  B mean = {zone['B_mean_m_yr']:+.2f} m/yr   "
              f"B range = {zone['B_range_m_yr']:.2f} m/yr "
              f"({'flat -> uniform-background model OK' if zone['B_range_m_yr'] < 1.5 else 'trends -> background varies alongshore'})")

    os.makedirs(OUT_DIR, exist_ok=True)
    pos_csv = os.path.join(OUT_DIR, f"{OUT_BASENAME}_positions_by_year.csv")
    chg_csv = os.path.join(OUT_DIR, f"{OUT_BASENAME}_change_{CHANGE_FROM}_{CHANGE_TO}.csv")
    fig_png = os.path.join(OUT_DIR, f"{OUT_BASENAME}_target_curve.png")
    zone_csv = os.path.join(OUT_DIR, f"{OUT_BASENAME}_quantified_summary.csv")
    fold_csv = os.path.join(OUT_DIR, f"{OUT_BASENAME}_foldsum.csv")

    tbl.to_csv(pos_csv)
    change_df.to_csv(chg_csv, index=False)
    pd.DataFrame([zone]).to_csv(zone_csv, index=False)
    if not foldsum.empty:
        foldsum.to_csv(fold_csv, index=False)
    plot_target(tbl_display, change_df).savefig(fig_png, dpi=200, bbox_inches="tight")

    print(f"\nSaved:")
    print(f"  {pos_csv}")
    print(f"  {chg_csv}")
    print(f"  {fig_png}")
    print(f"  {zone_csv}")
    if not foldsum.empty:
        print(f"  {fold_csv}")

    # -- fold-and-sum figure --
    if not foldsum.empty:
        fold_png = os.path.join(OUT_DIR, f"{OUT_BASENAME}_foldsum.png")
        plot_foldsum(foldsum).savefig(fold_png, dpi=200, bbox_inches="tight")
        print(f"  {fold_png}")

    # -- quantification dashboard (the one-glance summary figure) --
    dash_png = os.path.join(OUT_DIR, f"{OUT_BASENAME}_quantification_dashboard.png")
    plot_quantification_dashboard(change_df, zone, foldsum).savefig(
        dash_png, dpi=200, bbox_inches="tight")
    print(f"  {dash_png}")

    # -- time-axis trajectory panel --
    if MAKE_TRAJECTORY_PANEL:
        traj_png = os.path.join(OUT_DIR, f"{OUT_BASENAME}_trajectories.png")
        plot_trajectories(tbl_display).savefig(traj_png, dpi=200, bbox_inches="tight")
        print(f"  {traj_png}")

    # -- proportional-duration GIF --
    if MAKE_GIF:
        gif_path = os.path.join(OUT_DIR, f"{OUT_BASENAME}_position_over_time.gif")
        make_gif(tbl_display, gif_path)

    print("\nUse the rate_m_per_yr column as the per-domain target when tuning M.")
    print("Use A_groin_signal_m_yr (fold-sum) as the background-free groin target.")


if __name__ == "__main__":
    main()
