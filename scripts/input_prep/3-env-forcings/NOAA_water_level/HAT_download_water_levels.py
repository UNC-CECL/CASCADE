# =============================================================================
# HAT_download_water_levels.py
# NOAA CO-OPS Water Level Downloader — Duck, NC (Station 8651370)
# -----------------------------------------------------------------------------
# Purpose:
#   Produce a gap-checked, completeness-verified hourly water level CSV for use
#   as the `water_levels_file` input to historical_storm_creation_v3_HAT.py.
#
#   Output columns match that script's config exactly:
#       t = datetime (GMT)          -> t_name_water = "t"
#       v = water level [m NAVD88]  -> water_name   = "v"
#
# -----------------------------------------------------------------------------
# WHY THIS EXISTS:
#   The previous Duck CSV was missing the record's three largest events
#   (Gloria 1985, Halloween 1991, Isabel 2003). Root cause was not the API —
#   it was that download_noaa_water_levels_resumable.py accepted ANY non-empty
#   response as a complete month, cached it permanently, wrote the comparison CSV
#   even when months were missing, and checked completeness at the YEAR level
#   (< 6000 records) — too coarse to ever see a missing month (8760 - 720 =
#   8040, well above the threshold). v3's load_data() then dropna()'d the
#   absent hours in silence, so the storm file lost Isabel with no error.
#
# WHAT CHANGED, AND WHY:
#   Monthly chunks         Kept from the resumable script — one month (~720
#                          records) is far lighter on CO-OPS than a large
#                          window and times out much less. This was the right
#                          call and is preserved.
#
#   Cache reuse            Same folder and same {station}_{datum}_{YYYYMM}.csv
#                          naming, so months already downloaded are NOT re-
#                          fetched. Existing cache carries over.
#
#   Cache VALIDATION       NEW. Every month — cached or freshly fetched — is
#                          checked against its expected hour count. A truncated
#                          month is re-fetched rather than trusted forever.
#                          This is the fix for the actual bug.
#
#   Short-month handling   A month that stays short across retries is treated
#                          as a genuine gauge outage (not API truncation) and
#                          logged, not retried forever.
#
#   Fail loudly            Output CSV is NOT written if any month is empty or
#                          if any checked storm window is uncovered.
#
#   Conservative fill      Only gaps <= INTERP_LIMIT_HR are interpolated.
#                          Linear fill across a multi-hour outage at a storm
#                          peak flattens the peak and biases Rhigh low — the
#                          quiet version of the bug we just found.
#
# Author: Hannah Henry
# =============================================================================

import json
import os
import time
from calendar import monthrange
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from noaa_coops import Station
except ImportError:
    raise SystemExit("noaa_coops is not installed.  Install it:  pip install noaa_coops")

# =============================================================================
# USER CONFIGURATION
# =============================================================================

STATION_ID = "8651370"                  # Duck, NC
DATUM      = "NAVD"                     # NAVD88 — matches WIS and dune elevations
UNITS      = "metric"
TIME_ZONE  = "gmt"                      # WIS is UTC; do not change without checking WIS

BEGIN = "1984-01-01 00:00:00"
END   = "2024-12-31 23:00:00"           # note the 23:00 — the original script stopped at 00:00
                                        # and silently lost the last 23 hours

# --- Paths ---
# CACHE_DIR intentionally matches download_noaa_water_levels_resumable.py so the
# existing cache is reused. Point it at your real cache folder.
CACHE_DIR   = Path(r"/scripts/input_prep/3-env-forcings/NOAA_water_level/noaa_cache_8651370")
OUTPUT_DIR  = Path(r"/scripts/input_prep/3-env-forcings/NOAA_water_level")
# Derived from BEGIN/END so the filename always states its own span. This is a
# guard, not cosmetics: the previous broken CSV was named for the full span, and
# writing a different span under that name is how a stale file gets read as fresh.
OUTPUT_NAME = None   # None -> auto: {station}_DUCK_{begin}_{end}_{datum}.csv

# --- Request tuning ---
MAX_RETRIES   = 5
RETRY_WAIT    = 10                      # base seconds; grows each attempt
REQUEST_PAUSE = 1.0                     # seconds between API calls

# --- Completeness ---
MIN_MONTH_FRAC = 0.90                   # a month needs >= this fraction of its hours
REVALIDATE_CACHE = True                 # re-check cached months against expected hours.
                                        # Leave True — this is the bug fix. Set False only
                                        # to skip the (cheap, local) recount.

# --- Gap handling ---
INTERP_LIMIT_HR     = 3                 # fill gaps <= this; leave longer gaps NaN
REPORT_GAPS_OVER_HR = 3

# --- Storm windows that MUST be covered ---
# If the gauge is missing here, Rhigh is wrong in a way nothing downstream sees.
STORM_CHECKS = {
    "Gloria 1985":    ("1985-09-26 00:00", "1985-09-28 00:00"),
    "Halloween 1991": ("1991-10-30 00:00", "1991-11-02 00:00"),
    "Emily 1993":     ("1993-08-31 00:00", "1993-09-01 12:00"),
    "Floyd 1999":     ("1999-09-15 00:00", "1999-09-17 00:00"),
    "Isabel 2003":    ("2003-09-17 00:00", "2003-09-19 12:00"),
    "Irene 2011":     ("2011-08-26 00:00", "2011-08-28 12:00"),
    "Sandy 2012":     ("2012-10-28 00:00", "2012-10-30 12:00"),
}

# Storms known to be uncoverable from this gauge. The gate still fires for anything
# NOT listed here — so a new gap can never pass silently — but a documented outage
# does not block the run. The reason string is your methods-section note.
ACKNOWLEDGED_GAPS = {
    "Sandy 2012": "Duck gauge offline 2012-10-29 -> 2012-11-29 (751 h). CO-OPS has "
                  "no record; five retries returned an identical 34/61 hours. "
                  "Options if Sandy matters for Period 2: splice Oregon Inlet "
                  "(8652587) with an raw_offset/amplitude correction fit on overlapping "
                  "Duck data, reconstruct from predicted tide only (loses surge, "
                  "biases Rhigh low), or document as a known omission. Ask Laura.",
}

ABORT_IF_STORM_MISSING = True           # applies only to gaps NOT in ACKNOWLEDGED_GAPS


# =============================================================================
# STEP 1 — MONTH PLANNING
# =============================================================================

def month_chunks(begin, end):
    """Yield (label 'YYYYMM', start_ts, end_ts, expected_hours) per calendar month."""
    b, e = pd.Timestamp(begin), pd.Timestamp(end)
    cur = b.replace(day=1, hour=0, minute=0, second=0)
    out = []
    while cur <= e:
        nxt     = cur + pd.offsets.MonthBegin(1)
        m_start = max(b, cur)
        m_end   = min(e, nxt - pd.Timedelta(hours=1))
        expected = int((m_end - m_start).total_seconds() // 3600) + 1
        out.append((cur.strftime("%Y%m"), m_start, m_end, expected))
        cur = nxt
    return out


def cache_path(cache_dir, label):
    return cache_dir / f"{STATION_ID}_{DATUM}_{label}.csv"


# --- known-short registry -----------------------------------------------------
# A month that stays short across all retries is a real gauge outage, not API
# truncation. Without a record of that, REVALIDATE_CACHE re-fetches it (5 attempts
# with backoff, ~100 s) on EVERY run, forever. This registry remembers the
# confirmed count so the month is accepted from cache next time — but re-fetches
# if the count ever changes, in case CO-OPS backfills the record.

def short_registry_path(cache_dir):
    return cache_dir / "_known_short.json"


def load_known_short(cache_dir):
    p = short_registry_path(cache_dir)
    if not p.exists():
        return {}
    try:
        with open(p) as f:
            return json.load(f)
    except Exception:
        return {}


def save_known_short(cache_dir, reg):
    try:
        with open(short_registry_path(cache_dir), "w") as f:
            json.dump(reg, f, indent=2, sort_keys=True)
    except Exception as e:
        print(f"    (could not write known-short registry: {e})")


def read_cached_month(path):
    """Return DataFrame with columns t, v — or None if unreadable/empty."""
    try:
        df = pd.read_csv(path, parse_dates=["t"])
        if "t" not in df.columns or "v" not in df.columns or df.empty:
            return None
        df["v"] = pd.to_numeric(df["v"], errors="coerce")
        return df
    except Exception:
        return None


# =============================================================================
# STEP 2 — FETCH WITH RETRY, KEEPING THE BEST RESPONSE
# =============================================================================

def fetch_month(station, m_start, m_end, expected):
    """
    Download one month with retry/backoff, keeping the LARGEST response seen.

    Retrying matters because CO-OPS under load returns a truncated month rather
    than an error. Keeping the largest response distinguishes a transient
    truncation (a later attempt returns more) from a genuine gauge outage
    (every attempt returns the same short count).

    Returns (DataFrame or None, n_valid).
    """
    best, best_n = None, -1
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            raw = station.get_data(
                begin_date=m_start.strftime("%Y%m%d %H:%M"),
                end_date=m_end.strftime("%Y%m%d %H:%M"),
                product="hourly_height",
                datum=DATUM,
                units=UNITS,
                time_zone=TIME_ZONE,
            )
            if raw is not None and len(raw):
                vcol = "v" if "v" in raw.columns else raw.columns[0]
                # NOTE: raw.index is a DatetimeIndex NAMED 't'. Passing the Index
                # and a Series into DataFrame({...}) makes pandas adopt the Series'
                # index as the frame index — giving a frame with BOTH a column 't'
                # and an index named 't', which makes sort_values('t') ambiguous.
                # .values strips both, leaving a clean RangeIndex.
                df = pd.DataFrame({
                    "t": pd.to_datetime(raw.index).values,
                    "v": pd.to_numeric(raw[vcol], errors="coerce").values,
                })
                n = int(df["v"].notna().sum())
                if n > best_n:
                    best, best_n = df, n
                if n >= expected * MIN_MONTH_FRAC:
                    return best, best_n           # good enough, stop early
                print(f"    attempt {attempt}: only {n}/{expected} hours — retrying")
            else:
                print(f"    attempt {attempt}: empty response")
        except Exception as e:
            print(f"    attempt {attempt} failed ({str(e)[:45]})")

        if attempt < MAX_RETRIES:
            time.sleep(RETRY_WAIT * attempt)

    return best, max(best_n, 0)


# =============================================================================
# STEP 3 — DOWNLOAD LOOP
# =============================================================================

def download_all(cache_dir):
    """
    Walk every month: reuse valid cache, re-fetch invalid or missing.
    Returns (DataFrame[t, v], status DataFrame).
    """
    print("=" * 74)
    print(f"STEP 1: Monthly download — station {STATION_ID} ({BEGIN[:10]} -> {END[:10]})")
    print("=" * 74)

    station = Station(id=STATION_ID)
    existed = cache_dir.exists()
    cache_dir.mkdir(parents=True, exist_ok=True)
    chunks = month_chunks(BEGIN, END)

    # --- pre-flight: what is already on disk, and how long will this take? ---
    n_have = sum(1 for label, _, _, _ in chunks if cache_path(cache_dir, label).exists())
    n_todo = len(chunks) - n_have
    print(f"  cache: {cache_dir}")
    if not existed:
        print("         (created — no cache found, so this is a full download)")
    print(f"  months to cover:   {len(chunks)}")
    print(f"  already cached:    {n_have}")
    print(f"  to download:       {n_todo}")
    if n_todo:
        # ~1-3 s per call in practice, plus REQUEST_PAUSE between them
        lo = n_todo * (1.0 + REQUEST_PAUSE) / 60
        hi = n_todo * (3.0 + REQUEST_PAUSE) / 60
        print(f"  rough runtime:     {lo:.0f}-{hi:.0f} min (longer if CO-OPS is retrying)")
        print(f"  safe to interrupt: each month is cached as it lands; re-run resumes")
    print()

    parts, status = [], []
    known_short = load_known_short(cache_dir)
    if known_short:
        print(f"  known short months (confirmed gauge outages, not re-fetched): "
              f"{len(known_short)}\n")
    n_cache = n_fetch = n_refetch = n_short = n_empty = 0

    for label, m_start, m_end, expected in chunks:
        cf  = cache_path(cache_dir, label)
        df  = read_cached_month(cf) if cf.exists() else None
        src = "cache"

        n_valid = int(df["v"].notna().sum()) if df is not None else 0

        # Re-fetch if the cached month is short — UNLESS we already confirmed this
        # exact count is all CO-OPS has. Then it's a real outage, so accept it.
        confirmed_short = known_short.get(label) == n_valid
        stale = (REVALIDATE_CACHE
                 and n_valid < expected * MIN_MONTH_FRAC
                 and not confirmed_short)

        if df is None or stale:
            if df is not None and stale:
                print(f"  {label}: cached but only {n_valid}/{expected} hours — re-fetching")
                n_refetch += 1
            else:
                print(f"  {label}: downloading...")
            new_df, new_n = fetch_month(station, m_start, m_end, expected)
            time.sleep(REQUEST_PAUSE)

            # keep whichever is more complete: the cache or the new pull
            if new_df is not None and new_n >= n_valid:
                df, n_valid, src = new_df, new_n, "api"
                df.to_csv(cf, index=False)
                n_fetch += 1
            elif df is not None:
                src = "cache (re-fetch no better — likely a real gauge outage)"
        else:
            n_cache += 1

        if df is None or n_valid == 0:
            print(f"  {label}: NO DATA")
            n_empty += 1
            status.append({"month": label, "valid": 0, "expected": expected,
                           "frac": 0.0, "source": "empty"})
            continue

        if n_valid < expected * MIN_MONTH_FRAC:
            tag = "confirmed outage, from cache" if confirmed_short else src
            print(f"  {label}: SHORT — {n_valid}/{expected} hours "
                  f"({100*n_valid/expected:.0f}%)  [{tag}]")
            n_short += 1
            # remember it so the next run doesn't burn 5 retries on it again
            if known_short.get(label) != n_valid:
                known_short[label] = n_valid
                save_known_short(cache_dir, known_short)

        parts.append(df)
        status.append({"month": label, "valid": n_valid, "expected": expected,
                       "frac": n_valid / expected, "source": src})

    if not parts:
        raise SystemExit("Nothing downloaded. If a browser test URL also fails, CO-OPS is down.")

    # ignore_index=True discards every part's index at the stitch, so no part can
    # reintroduce an index level named 't' and make sort_values('t') ambiguous.
    full = (pd.concat(parts, ignore_index=True)
              .drop_duplicates(subset="t")
              .sort_values("t")
              .reset_index(drop=True)
              .set_index("t"))

    print(f"\n  months from cache: {n_cache}  |  fetched: {n_fetch}  "
          f"(of which re-fetched as stale: {n_refetch})")
    print(f"  short months: {n_short}  |  empty months: {n_empty}")
    return full, pd.DataFrame(status)


# =============================================================================
# STEP 4 — HOURLY GRID, GAP REPORT, CONSERVATIVE FILL
# =============================================================================

def build_hourly(df, begin, end, interp_limit_hr, report_over_hr):
    print("\n" + "=" * 74)
    print("STEP 2: Hourly grid, gap report, conservative fill")
    print("=" * 74)

    grid = pd.date_range(pd.Timestamp(begin), pd.Timestamp(end), freq="h")
    out  = df.reindex(grid)
    out.index.name = "t"

    n_total   = len(out)
    n_missing = int(out["v"].isna().sum())
    print(f"  Grid hours:       {n_total:,}")
    print(f"  Missing pre-fill: {n_missing:,}  ({100*n_missing/n_total:.2f}%)")

    # contiguous NaN runs
    isna, gaps, i = out["v"].isna().values, [], 0
    while i < len(isna):
        if isna[i]:
            j = i
            while j + 1 < len(isna) and isna[j + 1]:
                j += 1
            gaps.append({"start": out.index[i], "end": out.index[j], "hours": j - i + 1})
            i = j + 1
        else:
            i += 1
    gaps_df = pd.DataFrame(gaps, columns=["start", "end", "hours"])

    if len(gaps_df):
        big = gaps_df[gaps_df["hours"] > report_over_hr]
        print(f"  Gap runs: {len(gaps_df)}  |  longer than {report_over_hr} h: {len(big)}")
        if len(big):
            print(f"\n  Longest gaps (NOT filled — left as NaN):")
            for _, g in big.sort_values("hours", ascending=False).head(15).iterrows():
                print(f"    {g['start']:%Y-%m-%d %H:%M} -> {g['end']:%Y-%m-%d %H:%M}  "
                      f"({int(g['hours'])} h)")
            if len(big) > 15:
                print(f"    ... and {len(big) - 15} more (see gap log)")
    else:
        print("  No gaps.")

    # Fill short gaps only.
    #
    # CAREFUL: pandas interpolate(limit=N) fills the first N NaNs of EVERY run,
    # including long ones — a 37 h outage would get 3 fabricated hours at its
    # leading edge. So interpolate, then restore NaN across every run longer
    # than the limit. Only wholly-short runs survive as filled.
    out["v"] = out["v"].interpolate(method="linear",
                                    limit=interp_limit_hr,
                                    limit_area="inside")
    if len(gaps_df):
        long_gaps = gaps_df[gaps_df["hours"] > interp_limit_hr]
        for _, g in long_gaps.iterrows():
            out.loc[g["start"]:g["end"], "v"] = np.nan
        print(f"\n  Restored NaN across {len(long_gaps)} gap(s) > {interp_limit_hr} h "
              f"(no partial shoulder fill)")

    n_after = int(out["v"].isna().sum())
    print(f"  Filled gaps <= {interp_limit_hr} h  |  remaining NaN: {n_after:,} "
          f"({100*n_after/n_total:.2f}%)")
    print("  Remaining NaNs stay NaN by design — v3's load_data() drops those hours,")
    print("  which is correct provided none sit on a storm peak (checked next).")
    return out, gaps_df


# =============================================================================
# STEP 5 — STORM COVERAGE
# =============================================================================

def check_storm_coverage(df, checks, acknowledged):
    """
    Returns (unexpected_missing, acknowledged_missing) — lists of storm names.
    Only unexpected_missing should ever block a run.
    """
    print("\n" + "=" * 74)
    print("STEP 3: Coverage during known storms")
    print("=" * 74)

    unexpected, ack_hit = [], []
    for name, (s, e) in checks.items():
        s, e = pd.Timestamp(s), pd.Timestamp(e)
        if s < df.index.min() or e > df.index.max():
            print(f"  --  {name:<16} outside requested period, skipped")
            continue
        win   = df.loc[s:e, "v"]
        n_exp = len(win)
        n_ok  = int(win.notna().sum())
        pct   = 100 * n_ok / n_exp if n_exp else 0
        peak  = win.max()
        known = name in acknowledged
        if pct == 100:
            print(f"  OK  {name:<16} {n_ok}/{n_exp} hours   peak WL {peak:.3f} m")
        elif pct >= 80:
            tag = "KNOWN" if known else "PARTIAL, peak may be clipped"
            print(f"  !!  {name:<16} {n_ok}/{n_exp} hours ({pct:.0f}%)  "
                  f"peak WL {peak:.3f} m  — {tag}")
            (ack_hit if known else unexpected).append(name)
        else:
            tag = "MISSING (acknowledged)" if known else "MISSING"
            print(f"  XX  {name:<16} {n_ok}/{n_exp} hours ({pct:.0f}%)  — {tag}")
            (ack_hit if known else unexpected).append(name)

    print()
    if not unexpected and not ack_hit:
        print("  All checked storms covered.")
    if ack_hit:
        print("  Acknowledged gaps (documented, not blocking):")
        for n in ack_hit:
            print(f"    {n}: {acknowledged[n]}")
    if unexpected:
        print("  UNEXPECTED gaps — Rhigh will be biased low there:")
        for n in unexpected:
            print(f"    {n}")
    return unexpected, ack_hit


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("\nHAT_download_water_levels.py")
    print(f"Station {STATION_ID} | {DATUM} | {TIME_ZONE.upper()} | "
          f"{BEGIN[:10]} -> {END[:10]}\n")

    out_name = OUTPUT_NAME or (
        f"{STATION_ID}_DUCK_{pd.Timestamp(BEGIN):%Y%m%d}_"
        f"{pd.Timestamp(END):%Y%m%d}_{DATUM}.csv")

    raw, status = download_all(CACHE_DIR)
    hourly, gaps_df = build_hourly(raw, BEGIN, END, INTERP_LIMIT_HR, REPORT_GAPS_OVER_HR)
    unexpected, ack_hit = check_storm_coverage(hourly, STORM_CHECKS, ACKNOWLEDGED_GAPS)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    empty_months = status[status["frac"] == 0]
    short_months = status[(status["frac"] > 0) & (status["frac"] < MIN_MONTH_FRAC)]

    # --- logs always written, even on abort ---
    status.to_csv(OUTPUT_DIR / out_name.replace(".csv", "_month_status.csv"), index=False)
    if len(gaps_df):
        gaps_df.to_csv(OUTPUT_DIR / out_name.replace(".csv", "_gaps.csv"), index=False)

    print("\n" + "=" * 74)
    if len(short_months):
        print(f"SHORT MONTHS ({len(short_months)}) — likely real gauge outages, "
              f"re-fetch did not improve them:")
        for r in short_months.itertuples():
            print(f"  {r.month}: {int(r.valid)}/{int(r.expected)} ({100*r.frac:.0f}%)")

    abort = False
    if len(empty_months):
        print(f"\nABORT: {len(empty_months)} month(s) returned no data at all:")
        for r in empty_months.itertuples():
            print(f"  {r.month}")
        print("  Re-run — completed months are cached and will be skipped.")
        abort = True

    if unexpected and ABORT_IF_STORM_MISSING:
        print(f"\nABORT: {len(unexpected)} storm window(s) uncovered and NOT acknowledged:")
        for n in unexpected:
            print(f"  {n}")
        print("  Writing this CSV is how Isabel went missing last time. Either fix the")
        print("  gap, or add the storm to ACKNOWLEDGED_GAPS with a written reason so")
        print("  the omission is documented rather than silent.")
        abort = True

    # --- stale-comparison guard ---------------------------------------------------
    # If we are not going to write out_name this run, any pre-existing file with
    # that name is stale and will be silently read by v3_HAT.py. Quarantine it.
    out_path = OUTPUT_DIR / out_name
    if abort and out_path.exists():
        stamp  = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
        stale  = OUTPUT_DIR / out_name.replace(".csv", f"_STALE_{stamp}.csv")
        out_path.rename(stale)
        print(f"\n  QUARANTINED the pre-existing {out_name}")
        print(f"    -> {stale.name}")
        print(f"    It predates this run and cannot be trusted. Renaming it stops")
        print(f"    v3_HAT.py from reading it and reporting a stale result as fresh.")

    if abort:
        partial = OUTPUT_DIR / out_name.replace(".csv", "_PARTIAL.csv")
        hourly.to_csv(partial)
        print(f"\n  Wrote {partial.name} instead — the _PARTIAL suffix keeps it from")
        print(f"  being consumed downstream by mistake.")
        print("=" * 74)
        return None

    hourly.to_csv(out_path)
    print(f"Saved: {out_path}")
    print(f"  Rows: {len(hourly):,}  |  valid: {int(hourly['v'].notna().sum()):,}")
    print(f"  WL range: {hourly['v'].min():.3f} to {hourly['v'].max():.3f} m NAVD88")
    print(f"  Columns: t, v  ->  in v3_HAT.py set:")
    print(f"      water_levels_file = r\"{out_path}\"")
    print(f"      t_name_water = \"t\"")
    print(f"      water_name   = \"v\"")
    print("=" * 74)
    return hourly


if __name__ == "__main__":
    water_levels = main()
