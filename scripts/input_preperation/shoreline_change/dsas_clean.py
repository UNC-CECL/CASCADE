import pandas as pd

# Paths – update these
RAW_DSAS_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_rates.csv"
OUTPUT_DOMAIN_CSV = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1978_1997_CLEAN.csv"

# Load
df = pd.read_csv(RAW_DSAS_CSV)

# Keep only relevant columns
cols_keep = ["domain_id", "TransectId", "TCD", "LRR", "LR2", "NSM"]
df = df[cols_keep].copy()

# Optional quality filter: keep only good regression fits
df = df[df["LR2"] >= 0.5]   # you can tighten to 0.7 or 0.8 if you want

# Sort (not strictly required, but nice)
df = df.sort_values(["domain_id", "TCD"])

# Aggregate to domain: mean NSM (m) over the period
domain_means = (
    df.groupby("domain_id")["NSM"]
      .mean()
      .reset_index()
      .rename(columns={"NSM": "obs_total_change_m"})
)

# Save for overlay with model
domain_means.to_csv(OUTPUT_DOMAIN_CSV, index=False)

print("Saved domain-level DSAS shoreline change to:")
print(OUTPUT_DOMAIN_CSV)
