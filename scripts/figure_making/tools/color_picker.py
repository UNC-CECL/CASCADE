"""
Preview a linear colour gradient between two hex colours, with each step's hex code.

    python scripts/figure_making/tools/color_picker.py

Edit the two colours and the step count at the top; shows the swatches in a
window. Needs the rgb_gradient package. Details: scripts/figure_making/tools/README.md.

Author:  Hannah A. Henry, Coastal Environmental Change Lab,
         University of North Carolina at Chapel Hill
Contact: hahenry@unc.edu
Version: 2026-09-14
"""

from rgb_gradient import get_linear_gradient
import matplotlib.pyplot as plt

# Define start and end colors
intermediate_colors = [
    '#1B9E77',  #
    '#5E3C99'   #
]

# Generate gradient
gradient = get_linear_gradient(
    colors=intermediate_colors,
    nb_colors=9,
    return_format='hex'
)

# Print gradient colors
print("Generated Gradient:")
print(gradient)

# Display gradient
fig, ax = plt.subplots(figsize=(8, 2))

for i, color in enumerate(gradient):
    ax.add_patch(
        plt.Rectangle((i, 0), 1, 1, color=color)
    )

    ax.text(
        i + 0.5,
        -0.15,
        color,
        ha='center',
        va='top',
        fontsize=8,
        rotation=45
    )

# Formatting
ax.set_xlim(0, len(gradient))
ax.set_ylim(-0.4, 1)
ax.set_aspect('equal')
ax.axis('off')

plt.tight_layout()
plt.show()