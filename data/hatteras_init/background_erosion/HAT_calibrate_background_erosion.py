#!/usr/bin/env python3
"""
HATTERAS ISLAND: BACKGROUND EROSION RATE CALIBRATION SCRIPT

Purpose:
    Derives per-domain background erosion (BE) rates by computing the spatial
    residual between DSAS-observed and CASCADE-modeled shoreline change rates.
    The residual represents processes not captured by wave forcing alone (e.g.,
    inlet effects, sediment bypassing, structural erosion/accretion trends).

Method:
    1. Load DSAS observed rates and CASCADE modeled rates (from a no-BE run)
    2. Compute residual: DSAS - Model  (positive = needs more accretion in model)
    3. Mask domains where mismatch has a known non-BE cause (piers, shoals,
       inlet boundary effects) — these should NOT be corrected with BE
    4. Smooth the usable residuals (rolling mean) to avoid chasing noise
    5. Convert m/yr → dam/yr for CASCADE input
    6. Output a ready-to-paste BACKGROUND_EROSION_RATES list + diagnostic plots

Outputs:
    - Diagnostic plot (residual, smoothed BE, exclusion zones)
    - Printed BACKGROUND_EROSION_RATES list ready to paste into the hindcast script
    - CSV of per-domain BE values for records

Usage:
    Run AFTER completing a no-BE CASCADE hindcast run.
    Paste the printed BACKGROUND_EROSION_RATES into HAT_hindcast_1978_1997.py
    and re-run to assess fit improvement.

Validation note:
    BE rates calibrated on 1978-1997 should be tested against the 1997-2019
    period. If fit degrades substantially, BE may be absorbing wave calibration
    error rather than representing real persistent processes.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.ndimage import uniform_filter1d

# =============================================================================
# SECTION 1: CONFIGURATION
# =============================================================================

PROJECT_BASE_DIR   = r"C:\Users\hanna\PycharmProjects\CASCADE"
HATTERAS_DATA_BASE = os.path.join(PROJECT_BASE_DIR, "data", "hatteras_init")
OUTPUT_BASE_DIR    = os.path.join(PROJECT_BASE_DIR, "output", "raw_runs")
CALIB_OUTPUT_DIR   = os.path.join(PROJECT_BASE_DIR, "output", "background_rates")

# -----------------------------------------------------------------------------
# CASCADE model rate CSV (output from a no-BE hindcast run)
# Set to the run you want to calibrate against — typically your best-fit
# wave parameter run with BACKGROUND_EROSION_RATES all set to 0.
# -----------------------------------------------------------------------------
MODEL_RATE_CSV = os.path.join(
    OUTPUT_BASE_DIR,
    "HAT_1978_1997_SQ_filter_Hs2p0",
    "HAT_1978_1997_SQ_filter_Hs2p0_shoreline_change_rate.csv",
)
MODEL_RATE_COL = "model_rate_m_per_yr"

# -----------------------------------------------------------------------------
# DSAS observed rate CSV (same period as the model run above)
# -----------------------------------------------------------------------------
DSAS_CSV = os.path.join(
    HATTERAS_DATA_BASE, "shoreline_change",
    "dsas_1978_1997_domain_means_SIMPLE.csv",
)
DSAS_DOMAIN_COL = "domain_id"          # GIS domain IDs (1-90)
DSAS_RATE_COL   = "annual_rate_m_per_yr"

# -----------------------------------------------------------------------------
# Domain configuration (must match hindcast script)
# -----------------------------------------------------------------------------
NUM_REAL_DOMAINS   = 90
NUM_BUFFER_DOMAINS = 15
TOTAL_DOMAINS      = NUM_BUFFER_DOMAINS + NUM_REAL_DOMAINS + NUM_BUFFER_DOMAINS  # 120
START_REAL_INDEX   = NUM_BUFFER_DOMAINS   # = 15
END_REAL_INDEX     = START_REAL_INDEX + NUM_REAL_DOMAINS  # = 105
FIRST_GIS_ID       = 1
LAST_GIS_ID        = 90

# -----------------------------------------------------------------------------
# SMOOTHING
# Window size in domains (odd number recommended).
# 5 = ~2.5 km either side; 7 = ~3.5 km.
# Larger = smoother BE, less risk of overfitting spatial noise.
# -----------------------------------------------------------------------------
SMOOTH_WINDOW = 5

# -----------------------------------------------------------------------------
# EXCLUSION ZONES
# Domains where DSAS-model mismatch has a known non-BE cause.
# These are interpolated across rather than used as BE targets.
# Add or remove ranges as your understanding of the system evolves.
# -----------------------------------------------------------------------------
#
# Key:  (gis_start, gis_end, reason)
#
EXCLUSION_ZONES = [
    (1,   6,  "Oregon Inlet boundary / sediment shadow"),
    (24,  28, "Avon Fishing Pier (1963) — updrift accumulation"),
    (29,  46, "Wimble Shoals influence — offshore sediment sink"),
    (77,  81, "Rodanthe Fishing Pier (1960) — updrift accumulation"),
]

# -----------------------------------------------------------------------------
# CLAMP: maximum absolute BE magnitude to allow (dam/yr)
# Residuals beyond this threshold are likely wave calibration error,
# not genuine background processes. Flag but do not apply.
# Reference: site-wide calibrated value is -1.091 dam/yr = -0.1091 dam/yr
# -----------------------------------------------------------------------------
BE_CLAMP_DAM = 0.30   # dam/yr  (~3 m/yr)

# =============================================================================
# SECTION 2: LOAD DATA
# =============================================================================

os.makedirs(CALIB_OUTPUT_DIR, exist_ok=True)

print("=" * 70)
print("HATTERAS ISLAND — BACKGROUND EROSION CALIBRATION")
print("=" * 70)

# Load model rates
if not os.path.exists(MODEL_RATE_CSV):
    print(f"ERROR: Model rate CSV not found:\n  {MODEL_RATE_CSV}")
    sys.exit(1)
model_df = pd.read_csv(MODEL_RATE_CSV)
print(f"Loaded model rates: {len(model_df)} rows from {os.path.basename(MODEL_RATE_CSV)}")

# Load DSAS rates
if not os.path.exists(DSAS_CSV):
    print(f"ERROR: DSAS CSV not found:\n  {DSAS_CSV}")
    sys.exit(1)
dsas_df = pd.read_csv(DSAS_CSV)
print(f"Loaded DSAS rates:  {len(dsas_df)} rows from {os.path.basename(DSAS_CSV)}")

# =============================================================================
# SECTION 3: MERGE AND COMPUTE RESIDUALS
# =============================================================================

# Keep only real domains from model output
model_real = model_df.dropna(subset=["gis_domain_id"]).copy()
model_real["gis_domain_id"] = model_real["gis_domain_id"].astype(int)

# Merge on GIS domain ID
merged = model_real.merge(
    dsas_df[[DSAS_DOMAIN_COL, DSAS_RATE_COL]],
    left_on="gis_domain_id",
    right_on=DSAS_DOMAIN_COL,
    how="left",
)
merged = merged.rename(columns={DSAS_RATE_COL: "dsas_rate_m_per_yr"})
merged = merged.sort_values("gis_domain_id").reset_index(drop=True)

n_missing = merged["dsas_rate_m_per_yr"].isna().sum()
if n_missing > 0:
    print(f"  Warning: {n_missing} GIS domains have no DSAS observation (will interpolate)")

# Residual: positive = DSAS more accretional than model → BE should be positive
# Negative = DSAS more erosional than model → BE should be negative
merged["residual_m_per_yr"] = merged["dsas_rate_m_per_yr"] - merged[MODEL_RATE_COL]

print(f"\nResidual stats (m/yr, all real domains):")
print(f"  Mean : {merged['residual_m_per_yr'].mean():.3f}")
print(f"  Std  : {merged['residual_m_per_yr'].std():.3f}")
print(f"  Min  : {merged['residual_m_per_yr'].min():.3f}")
print(f"  Max  : {merged['residual_m_per_yr'].max():.3f}")

# =============================================================================
# SECTION 4: MASK EXCLUSION ZONES
# =============================================================================

# Build a boolean mask: True = use this domain for BE estimation
merged["use_for_BE"] = True
for gis_lo, gis_hi, reason in EXCLUSION_ZONES:
    mask = (merged["gis_domain_id"] >= gis_lo) & (merged["gis_domain_id"] <= gis_hi)
    merged.loc[mask, "use_for_BE"] = False

n_excluded = (~merged["use_for_BE"]).sum()
print(f"\nExclusion zones:")
for gis_lo, gis_hi, reason in EXCLUSION_ZONES:
    n = ((merged["gis_domain_id"] >= gis_lo) & (merged["gis_domain_id"] <= gis_hi)).sum()
    print(f"  GIS {gis_lo:2d}-{gis_hi:2d}  ({n} domains)  {reason}")
print(f"  Total excluded: {n_excluded} | Usable: {len(merged) - n_excluded}")

# =============================================================================
# SECTION 5: SMOOTH AND INTERPOLATE
# =============================================================================

# Working series indexed on GIS domain ID (1-90)
gis_ids = merged["gis_domain_id"].values
residual_raw = merged["residual_m_per_yr"].values.copy()

# NaN out excluded zones before smoothing
residual_masked = residual_raw.copy().astype(float)
residual_masked[~merged["use_for_BE"].values] = np.nan

# Linear interpolation across excluded/missing gaps before smoothing
# (uniform_filter1d cannot handle NaNs)
s = pd.Series(residual_masked, index=gis_ids)
s_interp = s.interpolate(method="linear", limit_direction="both")

# Rolling mean smooth
be_smooth_m = (
    s_interp
    .rolling(window=SMOOTH_WINDOW, center=True, min_periods=2)
    .mean()
    .values
)

# Convert m/yr → dam/yr
be_smooth_dam = be_smooth_m / 10.0

# =============================================================================
# SECTION 6: CLAMP AND FLAG OUTLIERS
# =============================================================================

be_final_dam = be_smooth_dam.copy()
clamped_domains = []
for i, (gis_id, val) in enumerate(zip(gis_ids, be_final_dam)):
    if np.isnan(val):
        be_final_dam[i] = 0.0
        continue
    if abs(val) > BE_CLAMP_DAM:
        clamped_domains.append((gis_id, val))
        be_final_dam[i] = np.sign(val) * BE_CLAMP_DAM

if clamped_domains:
    print(f"\n  WARNING: {len(clamped_domains)} domains clamped to ±{BE_CLAMP_DAM} dam/yr:")
    for gid, val in clamped_domains:
        print(f"    GIS {gid:2d}  raw BE = {val:.4f} dam/yr  → clamped")
    print("  Consider reviewing wave calibration for these domains.")
else:
    print(f"\n  No domains exceeded clamp threshold (±{BE_CLAMP_DAM} dam/yr). Good.")

# =============================================================================
# SECTION 7: BUILD FULL PADDED ARRAY (120 domains)
# =============================================================================

# Pad domain: buffers = 0, real domains = calibrated BE
be_padded = np.zeros(TOTAL_DOMAINS)
for i, gis_id in enumerate(gis_ids):
    pad_idx = START_REAL_INDEX + (gis_id - FIRST_GIS_ID)
    be_padded[pad_idx] = be_final_dam[i]

# =============================================================================
# SECTION 8: PRINT READY-TO-PASTE OUTPUT
# =============================================================================

# Community/segment labels for comments
def _segment_label(gis_id):
    if gis_id < 1:
        return "buf"
    if gis_id <= 6:
        return f"N.undev GIS{gis_id}"
    if gis_id <= 8:
        return f"Buxton  GIS{gis_id}"
    if gis_id <= 20:
        return f"Mid.undev GIS{gis_id}"
    if gis_id <= 31:
        return f"Avon    GIS{gis_id}"
    if gis_id <= 67:
        return f"C.undev GIS{gis_id}"
    if gis_id <= 83:
        return f"SWR     GIS{gis_id}"
    if gis_id <= 90:
        return f"S.undev GIS{gis_id}"
    return "buf"

print("\n" + "=" * 70)
print("READY-TO-PASTE BACKGROUND_EROSION_RATES")
print("(Copy into Section 5 of HAT_hindcast_1978_1997.py)")
print("=" * 70)
print()
print("# Background erosion rates calibrated from DSAS residuals")
print(f"# Source run : {os.path.basename(MODEL_RATE_CSV)}")
print(f"# DSAS period: {os.path.basename(DSAS_CSV)}")
print(f"# Smooth window: {SMOOTH_WINDOW} domains | Clamp: ±{BE_CLAMP_DAM} dam/yr")
print("#")
print("# NOTE: Buffer labels (_B) are always 0 — do not change.")
print("# Replace any value with a literal to override a specific domain.")
print()

# Print the named shorthand labels (informational — actual values are in the list)
print("_B   = 0.0   # Buffer (always 0 - do not change)")
print()
print("#                         pad:  0     1     2     3     4     5     6     7     8     9")
print("BACKGROUND_EROSION_RATES = [")

rows_per_line = 10
for row_start in range(0, TOTAL_DOMAINS, rows_per_line):
    row_end = min(row_start + rows_per_line, TOTAL_DOMAINS)
    vals = be_padded[row_start:row_end]

    # Determine GIS range for this row
    gis_range_parts = []
    for pad_idx in range(row_start, row_end):
        gis_id = pad_idx - START_REAL_INDEX + FIRST_GIS_ID
        if START_REAL_INDEX <= pad_idx < END_REAL_INDEX:
            gis_range_parts.append(gis_id)

    # Build value string
    val_strs = []
    for pad_idx in range(row_start, row_end):
        v = be_padded[pad_idx]
        is_buf = pad_idx < START_REAL_INDEX or pad_idx >= END_REAL_INDEX
        if is_buf:
            val_strs.append(" _B, ")
        else:
            val_strs.append(f"{v:5.3f},")

    val_line = " ".join(val_strs)

    # Comment
    if gis_range_parts:
        g0, g1 = gis_range_parts[0], gis_range_parts[-1]
        # Check if any domain in range was excluded
        exc_in_row = any(
            (merged.loc[merged["gis_domain_id"] == g, "use_for_BE"] == False).any()
            for g in gis_range_parts
            if g in merged["gis_domain_id"].values
        )
        exc_note = " [interp]" if exc_in_row else ""
        comment = f"# GIS {g0:2d}-{g1:2d}{exc_note}"
    else:
        buf_end = min(row_end - 1, START_REAL_INDEX - 1) if row_start < START_REAL_INDEX else row_start
        comment = f"# buf {row_start:3d}-{row_end - 1:3d}"

    print(f"    {val_line}  {comment}")

print("]")
print()
print("=" * 70)
print(f"Summary: {np.sum(be_padded != 0)} real domains with non-zero BE")
print(f"  Mean non-zero BE : {np.mean(be_padded[be_padded != 0]):.4f} dam/yr")
print(f"  Range            : {be_padded.min():.4f} to {be_padded.max():.4f} dam/yr")
print("=" * 70)

# =============================================================================
# SECTION 9: SAVE CSV
# =============================================================================

out_df = pd.DataFrame({
    "cascade_padded_index": np.arange(TOTAL_DOMAINS),
    "gis_domain_id": (
        [None] * START_REAL_INDEX
        + list(range(FIRST_GIS_ID, LAST_GIS_ID + 1))
        + [None] * (TOTAL_DOMAINS - END_REAL_INDEX)
    ),
    "be_dam_per_yr": be_padded,
})
csv_out = os.path.join(CALIB_OUTPUT_DIR, "HAT_calibrated_BE_rates.csv")
out_df.to_csv(csv_out, index=False)
print(f"\nSaved BE CSV: {csv_out}")

# =============================================================================
# SECTION 10: DIAGNOSTIC PLOTS
# =============================================================================

fig, axes = plt.subplots(3, 1, figsize=(16, 13), constrained_layout=True,
                          sharex=True)

# ── Panel 1: Raw rates (DSAS vs Model) ───────────────────────────────────────
ax = axes[0]
ax.plot(gis_ids, merged["dsas_rate_m_per_yr"], "o--", color="darkorange",
        linewidth=1.4, markersize=4, label="DSAS observed (m/yr)", zorder=3)
ax.plot(gis_ids, merged[MODEL_RATE_COL], "-", color="steelblue",
        linewidth=2, label="CASCADE model — no BE (m/yr)", zorder=3)
ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.7)

for gis_lo, gis_hi, reason in EXCLUSION_ZONES:
    ax.axvspan(gis_lo - 0.5, gis_hi + 0.5, alpha=0.12, color="red", zorder=0)

ax.set_ylabel("Shoreline change rate (m/yr)")
ax.set_title("Panel 1 — Observed vs Modeled rates (no background erosion)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# ── Panel 2: Raw residual + smoothed BE ──────────────────────────────────────
ax = axes[1]
ax.bar(gis_ids, residual_raw, color="gray", alpha=0.4, width=0.8,
       label="Raw residual DSAS - Model (m/yr)")
ax.plot(gis_ids, be_smooth_m, "-", color="firebrick", linewidth=2,
        label=f"Smoothed residual (window={SMOOTH_WINDOW}) (m/yr)")
ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.7)
ax.axhline( BE_CLAMP_DAM * 10, color="firebrick", linewidth=0.8,
            linestyle=":", alpha=0.6, label=f"Clamp ±{BE_CLAMP_DAM * 10:.1f} m/yr")
ax.axhline(-BE_CLAMP_DAM * 10, color="firebrick", linewidth=0.8,
            linestyle=":", alpha=0.6)

for gis_lo, gis_hi, reason in EXCLUSION_ZONES:
    ax.axvspan(gis_lo - 0.5, gis_hi + 0.5, alpha=0.12, color="red", zorder=0)
    mid = (gis_lo + gis_hi) / 2
    ax.text(mid, ax.get_ylim()[0] if ax.get_ylim()[0] < 0 else -0.5,
            reason.split("—")[0].strip(), ha="center", va="top",
            fontsize=6.5, color="red", style="italic")

ax.set_ylabel("Residual (m/yr)")
ax.set_title("Panel 2 — Residual and smoothed BE estimate (exclusion zones in red)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)

# ── Panel 3: Final calibrated BE (dam/yr) ────────────────────────────────────
ax = axes[2]
colors = ["firebrick" if v < 0 else "steelblue" for v in be_final_dam]
ax.bar(gis_ids, be_final_dam, color=colors, alpha=0.75, width=0.8,
       label="Calibrated BE (dam/yr)")
ax.axhline(0, color="gray", linewidth=0.8, linestyle="--", alpha=0.7)
ax.axhline( BE_CLAMP_DAM, color="black", linewidth=0.8, linestyle=":",
            alpha=0.5, label=f"Clamp ±{BE_CLAMP_DAM} dam/yr")
ax.axhline(-BE_CLAMP_DAM, color="black", linewidth=0.8, linestyle=":", alpha=0.5)

# Community span labels
community_spans = [
    (7,  8,  "Buxton"),
    (21, 31, "Avon"),
    (68, 83, "SWR"),
]
for gis_lo, gis_hi, label in community_spans:
    ax.axvspan(gis_lo - 0.5, gis_hi + 0.5, alpha=0.08, color="steelblue", zorder=0)
    ax.text((gis_lo + gis_hi) / 2, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 0.05,
            label, ha="center", va="bottom", fontsize=7,
            color="steelblue", style="italic")

for gis_lo, gis_hi, reason in EXCLUSION_ZONES:
    ax.axvspan(gis_lo - 0.5, gis_hi + 0.5, alpha=0.10, color="red", zorder=0)

erosion_patch   = mpatches.Patch(color="firebrick", alpha=0.75, label="Erosional BE (negative)")
accretion_patch = mpatches.Patch(color="steelblue", alpha=0.75, label="Accretional BE (positive)")
ax.legend(handles=[erosion_patch, accretion_patch], fontsize=9)

ax.set_ylabel("BE rate (dam/yr)")
ax.set_xlabel("GIS Domain ID (1-90)")
ax.set_title("Panel 3 — Final calibrated background erosion rates (CASCADE input)")
ax.grid(alpha=0.3)

# Shared x-axis ticks
xticks = np.arange(FIRST_GIS_ID, LAST_GIS_ID + 1, 5)
axes[2].set_xticks(xticks)
axes[2].set_xticklabels([str(i) for i in xticks], rotation=45, ha="right", fontsize=9)

# Exclusion zone legend (shared)
exc_patch = mpatches.Patch(color="red", alpha=0.12,
                            label="Excluded zone (interpolated)")
comm_patch = mpatches.Patch(color="steelblue", alpha=0.08,
                             label="Developed community")
fig.legend(handles=[exc_patch, comm_patch], loc="lower center",
           ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.02))

fig.suptitle(
    "Hatteras Island — Background Erosion Rate Calibration\n"
    f"Source: {os.path.basename(MODEL_RATE_CSV)}  |  "
    f"Smooth window: {SMOOTH_WINDOW} domains  |  Clamp: ±{BE_CLAMP_DAM} dam/yr",
    fontsize=11, fontweight="bold",
)

fig_out = os.path.join(CALIB_OUTPUT_DIR, "HAT_BE_calibration_diagnostic.png")
fig.savefig(fig_out, dpi=300, bbox_inches="tight")
print(f"Saved diagnostic plot: {fig_out}")
plt.show()
