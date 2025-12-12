import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ---- 1. Load a resampled .npy raster ----
f = Path(r"C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\elevations\2009_pea_hatteras\domain_60_resampled.npy")
arr = np.load(f)

print(f"Array shape: {arr.shape}")
print(f"Min={np.nanmin(arr):.2f} m, Max={np.nanmax(arr):.2f} m")

# ---- 2. Visualize elevation grid ----
plt.figure(figsize=(8,5))
plt.imshow(arr, cmap="terrain", origin="upper", aspect="auto")
plt.colorbar(label="Elevation (m NAVD88)")
plt.title(f"{f.name}: elevation grid")
plt.show()

# ---- 3. Quick profile checks ----
edge = min(5, arr.shape[1] // 20)
left_q  = np.nanpercentile(arr[:, :edge], 25)
right_q = np.nanpercentile(arr[:, -edge:], 25)
print(f"25th percentile near edges: left={left_q:.2f}, right={right_q:.2f}")
if right_q < left_q:
    print("✅ Ocean is likely on the RIGHT side (lower elevations).")
else:
    print("⚠️ Ocean might be on the LEFT side (check orientation).")

# ---- 4. Plot a few cross-shore profiles ----
for i in [0, arr.shape[0]//2, arr.shape[0]-1]:
    prof_lr = arr[i, :]
    prof = np.flip(prof_lr)  # flip if ocean should be at index 0
    plt.figure(figsize=(7,3))
    plt.plot(prof, lw=1)
    plt.axhline(0.26, ls='--', color='k', label='MHW (0.26 m)')
    plt.title(f"Row {i}: Ocean→Land profile")
    plt.legend()
    plt.tight_layout()
    plt.show()
