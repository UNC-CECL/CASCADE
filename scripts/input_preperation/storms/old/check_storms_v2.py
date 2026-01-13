# =============================================================================
# CASCADE Storm File Viewer (no assumed column names)
# -----------------------------------------------------------------------------
# - Loads .npy / .npz storm files
# - Infers likely fields (year, occurrence/probability, intensity, duration)
# - Prints tables & summary; plots distributions & time-like series
# Author: Hannah Henry
# Last modified: 2025-09-25
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ======= USER INPUT =======
STORM_PATH = Path(r"/data/hatteras_init/storms/hindcast_storms/OCR_1974_2022_Final_Hindcast_Storms.npy")
YEARS_KNOWN = (1978, 2022)   # only used for per-year averages if no explicit year column

# ======= LOADERS =======
def load_numpy(path: Path):
    if path.suffix.lower() == ".npz":
        z = np.load(path, allow_pickle=True)
        # pick largest array as main
        key = max(z.files, key=lambda k: z[k].size)
        return np.asarray(z[key]), {"npz_key": key, "all_keys": list(z.files)}
    arr = np.load(path, allow_pickle=True)
    return np.asarray(arr), {}

def is_structured(arr):
    return isinstance(arr, np.ndarray) and arr.dtype.names is not None

# ======= INFERENCE =======
def looks_like_year(x):
    s = pd.Series(x).dropna()
    if s.empty: return False
    ints = (np.abs(s - np.round(s)) < 1e-6).mean() > 0.95
    return ints and 1800 < s.median() < 2100

def looks_like_occurrence(x):
    s = pd.Series(x).dropna()
    if s.empty: return False
    uniq = np.unique(s.astype(float))
    return set(np.round(uniq).astype(int)).issubset({0,1})

def looks_like_probability(x):
    s = pd.Series(x).dropna()
    return (not s.empty) and s.min() >= -1e-6 and s.max() <= 1+1e-6 and s.var() > 0

def build_dataframe(obj):
    """
    Returns (df, notes) where df is a tidy DataFrame with whatever columns we can infer.
    If unknown, keeps generic c0..cN names.
    """
    notes = []
    if is_structured(obj):
        df = pd.DataFrame({k: obj[k] for k in obj.dtype.names})
        notes.append(f"Structured array with fields: {list(df.columns)}")
        return df, notes

    if isinstance(obj, np.ndarray) and obj.dtype == object and len(obj) and hasattr(obj[0], "keys"):
        df = pd.DataFrame([{k: d.get(k, np.nan) for k in obj[0].keys()} for d in obj])
        notes.append(f"Array of dict-like rows with keys: {list(df.columns)}")
        return df, notes

    arr = np.asarray(obj)
    if arr.ndim == 1:
        df = pd.DataFrame({"c0": arr})
        notes.append("1D numeric array; stored in column 'c0'.")
        return df, notes

    # 2D numeric
    df = pd.DataFrame(arr, columns=[f"c{i}" for i in range(arr.shape[1])])

    # Try to assign semantic names
    cols = list(df.columns)
    used = set()

    # year
    for c in cols:
        if looks_like_year(df[c].values):
            df.rename(columns={c: "year"}, inplace=True)
            used.add(c)
            break

    # occurrence (0/1) then probability (0..1)
    for c in cols:
        if c in used: continue
        if looks_like_occurrence(df[c].values):
            df.rename(columns={c: "occurrence"}, inplace=True)
            used.add(c); break

    for c in cols:
        if c in used: continue
        if looks_like_probability(df[c].values):
            df.rename(columns={c: "probability"}, inplace=True)
            used.add(c); break

    # intensity: pick remaining column with highest variance
    rem = [c for c in df.columns if c not in {"year", "occurrence", "probability"}]
    if rem:
        var_col = max(rem, key=lambda c: pd.Series(df[c]).var())
        if var_col not in {"year", "occurrence", "probability"}:
            df.rename(columns={var_col: "intensity"}, inplace=True)
            used.add(var_col)

    # duration: next highest variance
    rem = [c for c in df.columns if c not in {"year", "occurrence", "probability", "intensity"}]
    if rem:
        var_col2 = max(rem, key=lambda c: pd.Series(df[c]).var())
        if var_col2 not in {"year", "occurrence", "probability", "intensity"}:
            df.rename(columns={var_col2: "duration"}, inplace=True)
            used.add(var_col2)

    notes.append(f"Inferred columns: {list(df.columns)}")
    return df, notes

# ======= MAIN =======
arr, meta = load_numpy(STORM_PATH)
print(f"Loaded: {STORM_PATH}")
if meta: print(f"(npz info) selected key: {meta.get('npz_key')}, all keys: {meta.get('all_keys')}")

df, notes = build_dataframe(arr)
print("\n".join(notes))

# Coerce numerics
for c in df.columns:
    if df[c].dtype == object:
        df[c] = pd.to_numeric(df[c], errors="coerce")

print("\n=== HEAD ===")
print(df.head(12).to_string(index=False))

print("\n=== SUMMARY ===")
print(df.describe(percentiles=[.01,.25,.5,.75,.99]).T)

# Basic warnings
if "occurrence" in df:
    vals = df["occurrence"].dropna().unique()
    if not set(np.round(vals).astype(int)).issubset({0,1}):
        print("NOTE: 'occurrence' doesn’t look like 0/1; treating as probability-per-event instead.")
if "probability" in df:
    p = df["probability"].dropna()
    if not p.empty and (p.min() < 0 or p.max() > 1):
        print("WARNING: probability outside [0,1].")

# ======= PLOTS =======
plots = []

# Intensity
if "intensity" in df:
    plt.figure(figsize=(10,4))
    df["intensity"].dropna().hist(bins=20)
    plt.title("Storm intensity distribution")
    plt.xlabel("Intensity (units from file)")
    plt.ylabel("Count"); plt.tight_layout(); plt.show(); plots.append("intensity")

# Duration
if "duration" in df:
    plt.figure(figsize=(10,4))
    df["duration"].dropna().hist(bins=20, color="orange")
    plt.title("Storm duration distribution")
    plt.xlabel("Duration (hours?)"); plt.ylabel("Count")
    plt.tight_layout(); plt.show(); plots.append("duration")

# Event-wise probability/occurrence
for name in ("occurrence","probability"):
    if name in df:
        plt.figure(figsize=(10,3.5))
        plt.plot(df[name].values, lw=1)
        plt.title(f"{name.capitalize()} by event index")
        plt.xlabel("Event index"); plt.ylabel(name)
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show(); plots.append(name)

# Yearly aggregation if year exists
if "year" in df:
    grp = df.groupby("year", as_index=True)
    per_year = pd.DataFrame({
        "events": grp.size(),
        "occ_sum": grp["occurrence"].sum() if "occurrence" in df else None,
        "prob_sum": grp["probability"].sum() if "probability" in df else None,
        "mean_intensity": grp["intensity"].mean() if "intensity" in df else None,
        "max_intensity": grp["intensity"].max() if "intensity" in df else None,
    })
    print("\n=== PER-YEAR SUMMARY ===")
    print(per_year.head(15).to_string())

    if "events" in per_year:
        plt.figure(figsize=(10,3.5))
        per_year["events"].plot(marker="o")
        plt.title("Storm events per year")
        plt.xlabel("Year"); plt.ylabel("Count")
        plt.grid(alpha=0.3); plt.tight_layout(); plt.show(); plots.append("per_year")

# If no explicit year, compute average storms/year from occurrence/probability
if "year" not in df:
    n_years = YEARS_KNOWN[1] - YEARS_KNOWN[0] + 1
    occ_like = None
    if "occurrence" in df:
        occ_like = df["occurrence"].astype(float).clip(lower=0).sum()
    elif "probability" in df:
        occ_like = df["probability"].astype(float).clip(lower=0, upper=1).sum()
    if occ_like is not None:
        print(f"\nYears assumed: {YEARS_KNOWN[0]}–{YEARS_KNOWN[1]} ({n_years} years)")
        print(f"Total occurrence-like sum: {occ_like:.2f}")
        print(f"Average expected storms per year: {occ_like / n_years:.2f}")

# Fallback matrix preview if nothing else plotted
if not plots and isinstance(arr, np.ndarray) and arr.ndim == 2:
    plt.figure(figsize=(10,4))
    im = plt.imshow(arr, aspect="auto", cmap="viridis")
    plt.colorbar(im, label="value")
    plt.title("Storm matrix preview")
    plt.tight_layout(); plt.show()
