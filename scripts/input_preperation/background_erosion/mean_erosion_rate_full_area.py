import pandas as pd

csv_path = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\shoreline_change\dsas_1997_2019_domain_means_SIMPLE.csv"

df = pd.read_csv(csv_path)

col = "annual_rate_m_per_yr"  # change if needed

print(f"Number of domains: {len(df)}")
print(f"Average erosion rate: {df[col].mean():.3f} m/yr")
print(f"Std dev: {df[col].std():.3f}")

# Results
# Number of domains: 90
# Average erosion rate: -1.091 m/yr
# Std dev: 1.555