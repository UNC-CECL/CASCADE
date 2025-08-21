import os
import shutil
from glob import glob

# source parent directory (that has many subfolders)
source_dir = "/Users/rsahrae/PycharmProjects/Flood_Analysis/FLOOD_ANALYSIS/100_years_run_SQ/Input_Data"
# destination folder where you want to collect files
dest_dir = "/Users/rsahrae/PycharmProjects/Flood_Analysis/FLOOD_ANALYSIS/100_years_run_SQ/Input_Data/Runs"

os.makedirs(dest_dir, exist_ok=True)

# loop over all subfolders and move files
for folder in glob(os.path.join(source_dir, "*")):
    if os.path.isdir(folder):
        for file in glob(os.path.join(folder, "*")):
            if os.path.isfile(file):
                shutil.move(file, dest_dir)

print("✅ All files moved successfully!")
