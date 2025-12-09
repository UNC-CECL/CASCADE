import pandas as pd

# Input dune-line CSVs (raw, per-transect)
dune1978 = pd.read_csv(
    r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1978\1978_duneline_offset_raw.csv"
)
dune1997 = pd.read_csv(
    r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\island_offset\1997\1997_duneline_offset_raw.csv"
)

# ---- 1. Aggregate each year to one value per domain (mean ORIG_LEN) ----
d78 = (
    dune1978
    .groupby("domain_id")["ORIG_LEN"]
    .mean()
    .reset_index()
    .rename(columns={"ORIG_LEN": "dune_1978_m"})
)

d97 = (
    dune1997
    .groupby("domain_id")["ORIG_LEN"]
    .mean()
    .reset_index()
    .rename(columns={"ORIG_LEN": "dune_1997_m"})
)

# ---- 2. Merge per-domain means on domain_id ----
df = d78.merge(d97, on="domain_id", how="inner")

# ---- 3. Compute dune-line change (1997 - 1978) ----
# NEW: positive = seaward / accretion, negative = landward / retreat
df['dune_change_m'] = df['dune_1978_m'] - df['dune_1997_m']


# ---- 4. Sort and save ----
df = df.sort_values("domain_id")

outpath = (
    r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\duneline_change"
    r"\dune_change_1978_1997_by_domain.csv"
)
df.to_csv(outpath, index=False)

print("Saved dune-line change table to:", outpath)
print(df.head())