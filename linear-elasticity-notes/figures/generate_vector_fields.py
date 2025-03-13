#!/usr/bin/env python3
"""
This script generates SVG plots of 2D linear vector fields for several matrices.
The matrices represent:
  - a rotational transformation,
  - a pure shear transformation,
  - and a combined shear transformation (the sum of the above).

The vector field at point (x, y) is given by:
    v(x, y) = M * [x, y]^T.
Each plot is saved as an SVG file in the figures/ subdirectory.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

# Ensure the output directory exists.
output_dir = "figures"
os.makedirs(output_dir, exist_ok=True)

# Define the matrices.

# Compression component: represents the scalar part.
compression = np.array([[-1, 0],
                       [0, -1]])

# Rotational component: represents the antisymmetric part.
rotation = np.array([[0, 1],
                       [-1, 0]])

# Pure shear component: represents the symmetric, traceless part.
pure_shear = np.array([[0, 1],
                       [1, 0]])

# Combined shear: the sum of the above (e.g. corresponds to the shear matrix in the text)
traditional_shear = rotation + pure_shear

# You can add more matrices in this dictionary if desired.
matrices = {
    "compression": compression,
    "rotation": rotation,
    "pure_shear": pure_shear,
    "traditional_shear": traditional_shear
}

# Define a grid for the vector field.
x = np.linspace(-1, 1, 8)
y = np.linspace(-1, 1, 8)
X, Y = np.meshgrid(x, y)

# Generate plots for each matrix.
for name, M in matrices.items():
    # Compute the vector field: v = M * [x, y]^T at each grid point.
    U = M[0, 0]*X + M[0, 1]*Y
    V = M[1, 0]*X + M[1, 1]*Y

    # Create a figure.
    plt.figure(figsize=(5, 5))
    plt.quiver(X, Y, U, V, color="blue")
    
    # Set title with larger font
    plt.title(f"{name.replace('_', ' ').title()}", fontsize=14, pad=10)
    
    # Set axis labels
    plt.xlabel("x", fontsize=12)
    plt.ylabel("y", fontsize=12)
    
    # Set integer ticks
    plt.xticks([-1, 0, 1])
    plt.yticks([-1, 0, 1])
    
    plt.axis("equal")
    plt.grid(True, linestyle="--", alpha=0.5)

    # Save the figure as an SVG file.
    filename = os.path.join(output_dir, f"vector_field_{name}.svg")
    plt.savefig(filename, format="svg", bbox_inches='tight')
    plt.close()
    print(f"Saved vector field for {name} to {filename}")
