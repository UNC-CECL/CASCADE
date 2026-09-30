"""
Show one saved .npy elevation array as an image.

    python scripts/figure_making/tools/npy_view.py

The path is typed at the top and points into data/hatteras_init/topography/2009_FIXED/,
a folder that no longer exists: edit it before use. Details: scripts/figure_making/tools/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-14
"""

import numpy as np
import matplotlib.pyplot as plt

# Load a sample array
array = np.load(r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\topography\2009_FIXED\domain_1_topography_2009.npy')

# Visualize it
plt.imshow(array, cmap='terrain')
plt.title("Domain 1 Elevation (raw)")
plt.colorbar(label="Elevation (m NAVD88)")
plt.show()
