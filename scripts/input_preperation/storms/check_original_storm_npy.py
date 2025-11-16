# =============================================================================
# CASCADE Storm File Viewer (v3)
# -----------------------------------------------------------------------------
# Description:
#   Loads and visualizes a CASCADE-format storm numpy file, converting the
#   relative year index (0, 1, 2, ...) to the absolute calendar year (1978, ...).
# =============================================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load required library for plotting statistics
from scipy.stats import lognorm

# === USER INPUTS =============================================================
# The start year of the hindcast period, inferred from the file name.
START_YEAR = 1978
STORM_PATH = r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\HAT_1978_2022_Final_Hindcast_Storms.npy"

# Column names confirmed by your labmate:
COLS = ["Year_Index", "Rhigh", "Rlow", "Wave Period", "Duration"]

# === LOAD AND PREPARE DATA ===================================================
storm_arr = np.load(STORM_PATH)

print(f"Loaded storm file: {STORM_PATH}")
print(f"Array shape: {storm_arr.shape}")

# Convert to DataFrame
df = pd.DataFrame(storm_arr, columns=COLS)

# --- CRITICAL ADJUSTMENT ---
# Convert the relative Year_Index (0, 1, 2, ...) to the absolute Calendar Year (1978, 1979, 1980, ...)
df['Calendar_Year'] = df['Year_Index'].astype(int) + START_YEAR
# Drop the original index column
df = df.drop(columns=['Year_Index'])

# Clean up the remaining Year column (if you want to keep the index, rename it)
# df['Year_Index'] = df['Year_Index'].astype(int) # This column is now gone

print("\n=== HEAD (with Absolute Year) ===")
# Display the absolute year in the head
print(df[['Calendar_Year', 'Rhigh', 'Rlow', 'Wave Period', 'Duration']].head(10))

print("\n=== SUMMARY (numeric) ===")
print(df[["Rhigh", "Rlow", "Wave Period", "Duration"]].describe(percentiles=[0.25, 0.5, 0.75, 0.95, 0.99]))

# === VISUALIZATION (Updated for Absolute Years) ==============================

# 1. Histograms for key parameters
fig, axs = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle(f"Storm Parameter Distributions (N={len(df)} Events, {START_YEAR}-2022)", fontsize=16)

# Rhigh Distribution
df["Rhigh"].hist(ax=axs[0,0], bins=15, color="#1f77b4", edgecolor='black', alpha=0.8)
axs[0,0].set_title("1. Storm Intensity (Rhigh) Distribution")
axs[0,0].set_xlabel("Rhigh Value")
axs[0,0].set_ylabel("Count")

# Wave Period Distribution
df["Wave Period"].hist(ax=axs[0,1], bins=15, color="#ff7f0e", edgecolor='black', alpha=0.8)
axs[0,1].set_title("2. Wave Period Distribution")
axs[0,1].set_xlabel("Wave Period (s)")
axs[0,1].set_ylabel("Count")

# Duration Distribution
df["Duration"].hist(ax=axs[1,0], bins=15, color="#2ca02c", edgecolor='black', alpha=0.8)
axs[1,0].set_title("3. Storm Duration Distribution")
axs[1,0].set_xlabel("Duration (Hours)")
axs[1,0].set_ylabel("Count")

# Rlow vs Rhigh Scatter Plot
axs[1,1].scatter(df["Rlow"], df["Rhigh"], color="#d62728", alpha=0.7)
axs[1,1].set_title("4. Rlow vs. Rhigh Intensity Comparison")
axs[1,1].set_xlabel("Rlow Value")
axs[1,1].set_ylabel("Rhigh Value")
axs[1,1].grid(True, linestyle='--', alpha=0.6)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# 2. Time Series Plot (Now accurate with Calendar Year)
plt.figure(figsize=(15, 6))

# Group Rhigh by year and plot the average Rhigh for all storms in that year
yearly_avg_rhigh = df.groupby('Calendar_Year')['Rhigh'].mean().reset_index()

plt.plot(df['Calendar_Year'], df['Rhigh'], 'o', label='Individual Storm Rhigh', alpha=0.5, color='gray')
plt.plot(yearly_avg_rhigh['Calendar_Year'], yearly_avg_rhigh['Rhigh'], '-', linewidth=3, label='Average Rhigh per Year', color='#9467bd')

plt.title("Storm Intensity (Rhigh) Over Time (1978-2022)")
plt.xlabel("Calendar Year")
plt.ylabel("Rhigh Value (Intensity Index)")
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()

# 3. Simple Table Summary (using Calendar Year)
print("\n=== TOP 5 MOST INTENSE STORMS (by Rhigh) ===")
print(df.sort_values(by='Rhigh', ascending=False).head(5))

# 4. Quick check: storm count per year
storm_counts = df['Calendar_Year'].value_counts().sort_index()
print(f"\nYears in dataset: {df['Calendar_Year'].min()} to {df['Calendar_Year'].max()}")
print(f"Total Unique Years with recorded storms: {len(storm_counts)}")
print(f"Total Storms (Events): {len(df)}")
print("\nStorms per Year (Top 5):\n", storm_counts.head())