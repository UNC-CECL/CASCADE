import numpy as np
import matplotlib.pyplot as plt

# Load a sample array
array = np.load(r'C:\Users\hanna\PycharmProjects\CASCADE\data\hatteras_init\topography\2009_FIXED\domain_1_topography_2009.npy')

# Visualize it
plt.imshow(array, cmap='terrain')
plt.title("Domain 1 Elevation (raw)")
plt.colorbar(label="Elevation (m NAVD88)")
plt.show()
