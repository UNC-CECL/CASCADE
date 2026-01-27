import numpy as np
import pandas as pd
from pathlib import Path

# Input and output paths
STORM_PATH = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\storms\hindcast_storms\storms_1978_1997.npy")
EXCEL_PATH = STORM_PATH.with_suffix('.xlsx')  # Same location, .xlsx extension

# Load and convert
arr = np.load(STORM_PATH)
df = pd.DataFrame(arr, columns=["Year_Index", "Rhigh", "Rlow", "Wave Period", "Duration"])
df['Calendar_Year'] = df['Year_Index'].astype(int) + 1978
df = df[['Calendar_Year', 'Rhigh', 'Rlow', 'Wave Period', 'Duration']]

# Export
df.to_excel(EXCEL_PATH, index=False)
print(f"✓ Exported {len(df)} storms to: {EXCEL_PATH}")