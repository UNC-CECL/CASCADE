import os
from collections import Counter

folder = "/Users/rsahrae/PycharmProjects/Flood_Analysis/FLOOD_ANALYSIS/100_years_run_SQ/Input_Data/Runs"

files = [f for f in os.listdir(folder) if os.path.isfile(os.path.join(folder, f))]
counts = Counter(files)

duplicates = [f for f, c in counts.items() if c > 1]

if duplicates:
    print("Duplicate filenames found:")
    for f in duplicates:
        print(f)
else:
    print("No duplicate filenames.")
