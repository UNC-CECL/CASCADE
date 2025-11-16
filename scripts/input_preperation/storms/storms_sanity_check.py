# =============================================================================
# CASCADE Storm QC: tables + plots for .npy / .npz storms
# =============================================================================
import os, math
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- INPUT STORM DATA ----------
STORM_PATH = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\HAT_1978_2022_Final_Hindcast_Storms.npy")
OUT_CSV    = STORM_PATH.with_suffix("").as_posix() + "_qc.csv"
# --------------------------------

def load_numpy(path: Path):
    if path.suffix.lower() == ".npz":
        z = np.load(path, allow_pickle=True)
        if len(z.files) == 1:
            return z[z.files[0]]
        return {k: z[k] for k in z.files}
    return np.load(path, allow_pickle=True)

def looks_like_year(col):
    c = pd.Series(col).dropna()
    if c.empty: return False
    ints = (np.abs(c - np.round(c)) < 1e-6).mean() > 0.95
    rng  = c.max() - c.min()
    return ints and (r  := c.median()) > 1800 and r < 2100 and rng >= 5

def looks_like_occurrence(col):
    c = pd.Series(col).dropna()
    if c.empty: return False
    vals = np.unique(np.clip(c.astype(int), 0, 1))
    return set(vals.tolist()).issubset({0,1}) and c.min() >= -1e-6 and c.max() <= 1+1e-6

def looks_like_probability(col):
    c = pd.Series(col).dropna()
    if c.empty: return False
    return c.min() >= -1e-6 and c.max() <= 1+1e-6 and c.var() > 0

def infer_table(obj):
    """
    Return (df, notes). df columns subset of: year, occurrence, probability, intensity, duration, c0..cN
    """
    notes = []
    if isinstance(obj, dict):
        # pick largest 2D array
        key = max(obj.keys(), key=lambda k: obj[k].size)
        arr = np.asarray(obj[key])
        notes.append(f"NPZ detected; using key '{key}' (shape={arr.shape}).")
    else:
        arr = np.asarray(obj)

    if arr.ndim == 1:
        # maybe array of dicts/records
        if arr.dtype == object and hasattr(arr[0], "keys"):
            df = pd.DataFrame([{k: d.get(k, np.nan) for k in arr[0].keys()} for d in arr])
            return df, notes + ["Array of dict-like rows detected."]
        # otherwise, single column numeric
        df = pd.DataFrame({"value": arr})
        return df, notes + ["1D numeric array; stored under 'value'."]

    # 2D numeric: try to label columns
    ncol = arr.shape[1]
    colnames = [f"c{i}" for i in range(ncol)]
    df = pd.DataFrame(arr, columns=colnames)

    # guesses
    used = set()
    # year
    for i in range(ncol):
        if looks_like_year(df.iloc[:, i].values):
            df = df.rename(columns={f"c{i}": "year"}); used.add(i); break
    # occurrence/probability
    for i in range(ncol):
        if i in used: continue
        if looks_like_occurrence(df.iloc[:, i].values):
            df = df.rename(columns={f"c{i}": "occurrence"}); used.add(i); break
    for i in range(ncol):
        if i in used: continue
        if looks_like_probability(df.iloc[:, i].values):
            df = df.rename(columns={f"c{i}": "probability"}); used.add(i); break
    # intensity (pick the remaining column with largest variance)
    cand = [i for i in range(ncol) if i not in used]
    if cand:
        var_idx = max(cand, key=lambda j: pd.Series(df.iloc[:, j]).var())
        df = df.rename(columns={f"c{var_idx}": "intensity"}); used.add(var_idx)
    # duration (next best variance)
    cand = [i for i in range(ncol) if i not in used]
    if cand:
        var2 = max(cand, key=lambda j: pd.Series(df.iloc[:, j]).var())
        df = df.rename(columns={f"c{var2}": "duration"}); used.add(var2)

    notes.append(f"Inferred columns: {list(df.columns)}")
    return df, notes

# ---- Run QC ----
if not STORM_PATH.exists():
    raise FileNotFoundError(f"Storm file not found: {STORM_PATH}")

raw = load_numpy(STORM_PATH)
df, notes = infer_table(raw)
print("\n".join(notes))

# Basic cleaning
for col in df.columns:
    if df[col].dtype == object:
        with np.errstate(all='ignore'):
            df[col] = pd.to_numeric(df[col], errors="coerce")

# Tables
print("\n=== HEAD ===")
print(df.head(10).to_string(index=False))

print("\n=== SUMMARY (numeric) ===")
print(df.describe(percentiles=[.01,.25,.5,.75,.99]).T)

# Per-year aggregation (if year present)
if "year" in df.columns:
    yr = (
        df
        .assign(storm=(df.get("occurrence", pd.Series([np.nan]*len(df))).fillna(0) > 0).astype(int))
        .groupby("year", as_index=False)
        .agg(storms=("storm","sum"),
             mean_intensity=("intensity","mean"),
             max_intensity=("intensity","max"),
             mean_duration=("duration","mean"))
    )
    print("\n=== PER-YEAR SUMMARY ===")
    print(yr.head(15).to_string(index=False))
    # export
    try:
        yr.to_csv(OUT_CSV, index=False)
        print(f"\nSaved per-year summary → {OUT_CSV}")
    except Exception as e:
        print(f"(could not save CSV: {e})")

# Flags / sanity checks
def flag(msg): print("WARNING:", msg)

if "occurrence" in df.columns:
    vals = df["occurrence"].dropna().unique()
    if not set(np.unique(vals)).issubset({0,1}):
        flag("occurrence contains values other than 0/1.")
if "probability" in df.columns:
    p = df["probability"].dropna()
    if not p.empty and (p.min() < 0 or p.max() > 1):
        flag("probability outside [0,1].")

# ---- Plots ----
plots = 0
if "year" in df.columns and "occurrence" in df.columns:
    plt.figure(figsize=(10,4))
    yearly = df.groupby("year")["occurrence"].sum()
    plt.plot(yearly.index.values, yearly.values)
    plt.title("Storm count per year")
    plt.xlabel("Year")
    plt.ylabel("Count")
    plt.grid(alpha=0.3)
    plt.tight_layout(); plt.show(); plots += 1

if "intensity" in df.columns:
    plt.figure(figsize=(10,4))
    v = df["intensity"].dropna().values
    bins = max(10, int(np.sqrt(len(v)))) if len(v) > 0 else 10
    plt.hist(v, bins=bins)
    plt.title("Storm intensity distribution")
    plt.xlabel("Intensity (units as in file)")
    plt.ylabel("Frequency")
    plt.tight_layout(); plt.show(); plots += 1

if "probability" in df.columns:
    plt.figure(figsize=(10,4))
    plt.plot(df.index.values, df["probability"].values)
    plt.title("Storm probability (timeseries)")
    plt.xlabel("Row index")
    plt.ylabel("Probability (0–1)")
    plt.grid(alpha=0.3)
    plt.tight_layout(); plt.show(); plots += 1

# If nothing inferred, still show a clean heatmap of the matrix
if plots == 0 and isinstance(raw, np.ndarray) and raw.ndim == 2:
    plt.figure(figsize=(10,4))
    im = plt.imshow(raw, aspect='auto')
    plt.colorbar(im, label="value")
    plt.title("Storm matrix preview")
    plt.tight_layout(); plt.show()
