"""
CoastSat Pipeline Data Verification
=====================================
Runs targeted checks on the outputs of the CoastSat → CASCADE pipeline
to confirm that transect LRR values, domain averages, and transect counts
are internally consistent before these values enter CASCADE.

Checks performed
----------------
  1. NaN audit        — how many transects have missing LRR, and why
  2. Domain mean match — do the pre-computed domain means equal the manual mean
                         of transect LRRs in the same CSV?
  3. Transect count   — does n_transects in the summary match the actual count
                         in the transect file?
  4. Per-domain detail — prints a side-by-side table for every domain so you
                         can spot any domain that looks off

Usage
-----
Edit the CONFIG section to point at your files, then run:
    python coastsat_verify_pipeline.py

Outputs
-------
  Console report (always)
  verification_report.csv  — detailed per-domain comparison table
"""

# ============================================================
# CONFIG — edit these paths to match your actual file locations
# ============================================================

# Full transect-level LRR results (comparison of coastsat_domain_lrr_fixed.py)
TRANSECT_LRR_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\transect_lrr_full.csv"

# Domain-level summary (comparison of coastsat_domain_lrr_fixed.py)
DOMAIN_SUMMARY_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat\2004_2024\domain_lrr_summary.csv"

# Where to save the per-domain comparison table
OUTPUT_REPORT_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\scripts\input_preperation\CoastSat_verification\verification_report_2004_2024.csv"

# Tolerance for floating-point comparison in Check 2 (m/yr)
# Differences smaller than this are treated as matching
TOLERANCE = 0.001

# ============================================================
# IMPORTS
# ============================================================
import os
import numpy as np
import pandas as pd

# ============================================================
# HELPERS
# ============================================================

PASS = "  ✓ PASS"
WARN = "  ⚠ WARN"
FAIL = "  ✗ FAIL"

def section(title):
    print("\n" + "=" * 60)
    print(f"  {title}")
    print("=" * 60)

def load_data():
    """Load both CSVs with clear error messages if files are missing."""
    missing = []
    for path, name in [(TRANSECT_LRR_CSV, "transect_lrr_full.csv"),
                       (DOMAIN_SUMMARY_CSV, "domain_lrr_summary.csv")]:
        if not os.path.exists(path):
            missing.append(f"  NOT FOUND: {path}  ({name})")
    if missing:
        print("\nERROR — Cannot run validation. Missing files:")
        for m in missing:
            print(m)
        raise SystemExit(1)

    t = pd.read_csv(TRANSECT_LRR_CSV)
    d = pd.read_csv(DOMAIN_SUMMARY_CSV)

    # Standardise domain column name (handles "domain_lrr_summary.csvdomain_number" corruption)
    d.columns = [c.split(".csv")[-1] if ".csv" in c else c for c in d.columns]
    for col in ["domain_number", "mean_lrr", "std_lrr", "n_valid", "n_transects"]:
        if col in d.columns:
            d[col] = pd.to_numeric(d[col], errors="coerce")

    for col in ["domain_number", "lrr_m_yr", "n_obs"]:
        if col in t.columns:
            t[col] = pd.to_numeric(t[col], errors="coerce")

    return t, d


# ============================================================
# CHECK 1 — NaN audit
# ============================================================

def check_nan_lrr(t):
    """
    How many transects are missing LRR values?
    A large number here usually means transect IDs didn't match CSV filenames.
    """
    section("CHECK 1 — Missing LRR values (NaN audit)")

    total      = len(t)
    n_valid    = t["lrr_m_yr"].notna().sum()
    n_nan      = t["lrr_m_yr"].isna().sum()
    pct_nan    = n_nan / total * 100

    print(f"  Total transects       : {total:,}")
    print(f"  With valid LRR        : {n_valid:,}  ({100 - pct_nan:.1f}%)")
    print(f"  Missing LRR (NaN)     : {n_nan:,}  ({pct_nan:.1f}%)")

    if n_nan == 0:
        print(f"\n{PASS}  All transects have valid LRR values.")
        nan_by_domain = pd.Series(dtype=int)
    else:
        # Show which domains are most affected
        nan_by_domain = (
            t[t["lrr_m_yr"].isna()]
            .groupby("domain_number")["transect_id"]
            .count()
            .sort_values(ascending=False)
        )
        print(f"\n  Domains with most NaN transects:")
        print(nan_by_domain.head(10).to_string())

        # Check n_obs for NaN rows — small n_obs = too few observations
        if "n_obs" in t.columns:
            nan_rows = t[t["lrr_m_yr"].isna()]
            zero_obs = (nan_rows["n_obs"] == 0).sum()
            few_obs  = (nan_rows["n_obs"].between(1, 2)).sum()
            print(f"\n  Of {n_nan} NaN transects:")
            print(f"    {zero_obs:>4} had 0 observations (CSV not found or empty)")
            print(f"    {few_obs:>4} had 1–2 observations (below MIN_OBS threshold)")
            print(f"    {n_nan - zero_obs - few_obs:>4} had other/unknown reason")

        if pct_nan > 20:
            print(f"\n{FAIL}  More than 20% of transects are missing LRR.")
            print(f"       Likely cause: transect IDs in the lookup table do not")
            print(f"       match the filenames of the time-series CSVs.")
            print(f"       Fix: verify TRANSECT_ID_COL in coastsat_domain_mapping.py")
            print(f"            matches the CSV filename format (e.g. usa_NC_0032_0021).")
        else:
            print(f"\n{WARN}  Some NaN values present. Review the domains listed above.")

    return nan_by_domain


# ============================================================
# CHECK 2 — Domain mean consistency
# ============================================================

def check_domain_means(t, d):
    """
    Does the mean_lrr in domain_lrr_summary.csv equal the manual mean
    of lrr_m_yr values for each domain in transect_lrr_full.csv?

    A mismatch means something changed between when the summary was computed
    and the current transect file, or a different set of transects was used.
    """
    section("CHECK 2 — Domain mean consistency")

    # Manual domain means from transect file (valid LRR only)
    manual = (
        t[t["lrr_m_yr"].notna()]
        .groupby("domain_number")["lrr_m_yr"]
        .mean()
        .rename("manual_mean_lrr")
        .reset_index()
    )
    manual["domain_number"] = manual["domain_number"].astype(float)

    # Summary means
    summary = d[["domain_number", "mean_lrr"]].copy()
    summary["domain_number"] = summary["domain_number"].astype(float)

    merged = summary.merge(manual, on="domain_number", how="outer")
    merged["difference"] = (merged["mean_lrr"] - merged["manual_mean_lrr"]).abs()
    merged["match"] = merged["difference"] <= TOLERANCE

    n_match    = merged["match"].sum()
    n_mismatch = (~merged["match"]).sum()
    n_domains  = len(merged)

    print(f"  Domains compared      : {n_domains}")
    print(f"  Matching (diff < {TOLERANCE} m/yr) : {n_match}")
    print(f"  Mismatching           : {n_mismatch}")

    if n_mismatch == 0:
        print(f"\n{PASS}  All domain means in the summary match manual calculation.")
        print(f"       Max difference : {merged['difference'].max():.6f} m/yr")
    else:
        print(f"\n{FAIL}  {n_mismatch} domain(s) have inconsistent means.")
        print(f"       These domains may have been computed from a different transect set.")
        bad = merged[~merged["match"]].sort_values("difference", ascending=False)
        print("\n  Mismatching domains:")
        print(bad[["domain_number", "mean_lrr", "manual_mean_lrr", "difference"]]
              .to_string(index=False))

    return merged


# ============================================================
# CHECK 3 — Transect count consistency
# ============================================================

def check_transect_counts(t, d):
    """
    Does n_transects in domain_lrr_summary.csv match the actual count
    of rows per domain in transect_lrr_full.csv?

    Includes both valid and NaN-LRR transects (n_transects = total in domain).
    """
    section("CHECK 3 — Transect count per domain")

    # Actual counts from transect file
    actual = (
        t.groupby("domain_number")["transect_id"]
        .count()
        .rename("actual_count")
        .reset_index()
    )
    actual["domain_number"] = actual["domain_number"].astype(float)

    # Summary counts
    if "n_transects" not in d.columns:
        print(f"\n{WARN}  'n_transects' column not found in domain summary. Skipping.")
        return actual

    summary_counts = d[["domain_number", "n_transects"]].copy()
    summary_counts["domain_number"] = summary_counts["domain_number"].astype(float)

    merged = summary_counts.merge(actual, on="domain_number", how="outer")
    merged["count_diff"] = (merged["n_transects"] - merged["actual_count"]).abs()
    merged["match"] = merged["count_diff"] == 0

    n_match    = merged["match"].sum()
    n_mismatch = (~merged["match"]).sum()

    print(f"  Domains compared      : {len(merged)}")
    print(f"  Matching counts       : {n_match}")
    print(f"  Mismatching counts    : {n_mismatch}")

    # Overall distribution
    print(f"\n  Transects per domain (from transect file):")
    desc = merged["actual_count"].describe()
    print(f"    Min    : {desc['min']:.0f}")
    print(f"    Median : {desc['50%']:.0f}")
    print(f"    Mean   : {desc['mean']:.1f}")
    print(f"    Max    : {desc['max']:.0f}")
    print(f"    Total  : {merged['actual_count'].sum():.0f}")

    if n_mismatch == 0:
        print(f"\n{PASS}  Transect counts match for all domains.")
    else:
        print(f"\n{FAIL}  Count mismatch in {n_mismatch} domains.")
        bad = merged[~merged["match"]].sort_values("count_diff", ascending=False)
        print(bad[["domain_number", "n_transects", "actual_count", "count_diff"]]
              .to_string(index=False))

    # Flag any domain with unusually high count (potential data issue)
    median_count = merged["actual_count"].median()
    outliers = merged[merged["actual_count"] > median_count * 5]
    if len(outliers) > 0:
        print(f"\n{WARN}  Domains with unusually high transect count (>5× median of {median_count:.0f}):")
        print(outliers[["domain_number", "actual_count"]].to_string(index=False))
        print(f"       This may indicate correct geography (e.g. Wimble Shoals area)")
        print(f"       or a data issue — check these domains in the validation map.")

    return merged


# ============================================================
# CHECK 4 — Per-domain detail table
# ============================================================

def build_detail_table(t, d, count_check):
    """
    Build a side-by-side per-domain table combining:
      - Values from domain_lrr_summary.csv
      - Manually computed values from transect_lrr_full.csv
    """
    section("CHECK 4 — Per-domain detail table")

    t_valid = t[t["lrr_m_yr"].notna()].copy()

    manual = (
        t_valid.groupby("domain_number")["lrr_m_yr"]
        .agg(
            manual_n_valid  = "count",
            manual_mean_lrr = "mean",
            manual_std_lrr  = "std",
            manual_min_lrr  = "min",
            manual_max_lrr  = "max",
        )
        .reset_index()
    )
    manual["domain_number"] = manual["domain_number"].astype(float)

    # Total transect count (including NaN LRR)
    total_count = (
        t.groupby("domain_number")["transect_id"]
        .count()
        .rename("actual_total_transects")
        .reset_index()
    )
    total_count["domain_number"] = total_count["domain_number"].astype(float)

    d2 = d.copy()
    d2["domain_number"] = d2["domain_number"].astype(float)

    detail = (
        d2
        .merge(manual, on="domain_number", how="outer")
        .merge(total_count, on="domain_number", how="outer")
        .sort_values("domain_number")
    )

    detail["mean_diff"]  = (detail["mean_lrr"] - detail["manual_mean_lrr"]).abs().round(4)
    detail["count_ok"]   = detail["n_transects"] == detail["actual_total_transects"]
    detail["mean_ok"]    = detail["mean_diff"] <= TOLERANCE

    # Flag any domain where both match → True, otherwise → False
    detail["all_ok"] = detail["count_ok"] & detail["mean_ok"]

    n_problem = (~detail["all_ok"]).sum()
    if n_problem == 0:
        print(f"\n{PASS}  All {len(detail)} domains are internally consistent.")
    else:
        print(f"\n{FAIL}  {n_problem} domain(s) have inconsistencies — review the table below.")

    # Print full table (truncated columns for readability)
    display_cols = [
        "domain_number",
        "n_transects", "actual_total_transects",
        "n_valid",     "manual_n_valid",
        "mean_lrr",    "manual_mean_lrr", "mean_diff",
        "std_lrr",
        "all_ok",
    ]
    display_cols = [c for c in display_cols if c in detail.columns]

    pd.set_option("display.max_rows", 100)
    pd.set_option("display.width", 120)
    pd.set_option("display.float_format", "{:+.3f}".format)
    print()
    print(detail[display_cols].to_string(index=False))

    return detail


# ============================================================
# MAIN
# ============================================================

def main():
    print("\n" + "=" * 60)
    print("  CoastSat Pipeline Verification")
    print("=" * 60)

    # ── Load ──────────────────────────────────────────────────
    print("\nLoading files...")
    print(f"  Transect CSV : {os.path.basename(TRANSECT_LRR_CSV)}")
    print(f"  Summary CSV  : {os.path.basename(DOMAIN_SUMMARY_CSV)}")
    t, d = load_data()
    print(f"  Transect rows: {len(t):,}")
    print(f"  Domain rows  : {len(d):,}")

    # ── Checks ────────────────────────────────────────────────
    nan_by_domain = check_nan_lrr(t)
    mean_check    = check_domain_means(t, d)
    count_check   = check_transect_counts(t, d)
    detail        = build_detail_table(t, d, count_check)

    # ── Summary ───────────────────────────────────────────────
    section("VERIFICATION SUMMARY")
    checks = {
        "Check 1 — NaN audit"         : t["lrr_m_yr"].isna().sum() == 0,
        "Check 2 — Domain means match": mean_check["match"].all(),
        "Check 3 — Transect counts"   : ("count_diff" in count_check.columns and
                                         count_check["count_diff"].sum() == 0),
        "Check 4 — All domains OK"    : detail["all_ok"].all(),
    }
    all_pass = all(checks.values())
    for name, passed in checks.items():
        print(f"  {'✓' if passed else '✗'}  {name}")

    if all_pass:
        print(f"\n{PASS}  ALL CHECKS PASSED — data appears internally consistent.")
    else:
        print(f"\n{FAIL}  ONE OR MORE CHECKS FAILED — review the output above.")
        print(f"       Focus on checks marked ✗ before using this data in CASCADE.")

    # ── Save report ───────────────────────────────────────────
    detail.to_csv(OUTPUT_REPORT_CSV, index=False)
    print(f"\nDetailed report saved to:")
    print(f"  {OUTPUT_REPORT_CSV}")

    print("\n" + "=" * 60 + "\n")


if __name__ == "__main__":
    main()
