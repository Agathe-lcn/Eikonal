import numpy as np
import matplotlib.pyplot as plt 
import os

def compute_error():
    # Settings
    n = 200
    m = 200
    length = 1.0
    h = length / n

    # Creating the figure with 2 or 3 subplots
    #fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig, axes = plt.subplots(1, 2, figsize=(12,5))

    # Loading the matrix calculated by FIM
    filename = "matrix_fim.txt"
    matrix = np.loadtxt(filename, dtype=str)

    # Conversion to float
    T_fim = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            if matrix[i,j] == 'inf' or matrix[i,j] == '1e+99':
                T_fim[i,j] = np.nan
            else:
                T_fim[i,j] = float(matrix[i,j])

    # Loading source coordinates
    coord_file = "coords_mesh.txt"
    if os.path.exists(coord_file):
        coords = np.loadtxt(coord_file)
        
        x_src, y_src = coords[0], coords[1]

    # Euclidean distance
    T_euclid = np.zeros((n,m))
    for i in range(n):
        y = i * h 
        dy = y - y_src
        dy2 = dy * dy
        for j in range(m):
            x = j * h 
            dx = x - x_src
            dx2 = dx * dx
            T_euclid[i,j] = np.sqrt(dx2 + dy2)

    # Error = | T_euclid - T_fim |
    error = np.abs(T_euclid - T_fim)
    
    # Statistics
    valid_error = error[~np.isnan(error)]
    print("\nError statistics |T_euclid - T_fim|\n")
    print(f"Maximum error: {np.max(valid_error):.6e}")
    print(f"Minimal error: {np.min(valid_error):.6e}")
    print(f"Average error: {np.mean(valid_error):.6e}")
    print(f"Median error: {np.median(valid_error):.6e}")

    # Visualization
    x = np.linspace(0, length, m)
    y = np.linspace(0, length, n)
    X, Y = np.meshgrid(x,y)

    # Figure 1: FIM solution
    im1 = axes[0].pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    axes[0].set_title("FIM-distance map")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    axes[0].scatter(x_src, y_src, color='red', s=10, marker='.', label="Source")
    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0])

    # Figure 2: Euclidean distance
    #im2 = axes[1].pcolormesh(X, Y, T_euclid, shading='nearest', cmap='viridis')
    #axes[1].set_title("Euclidean distance")
    #axes[1].set_xlabel("X")
    #axes[1].set_ylabel("Y")
    #axes[1].scatter(x_src, y_src, color='red', s=10, marker='.', label="Source")
    #axes[1].legend()
    #axes[1].axis('equal')
    #plt.colorbar(im2, ax=axes[1], label="Distance")


    # Figure 3: Error | T_euclid - T_fim |
    im3 = axes[1].pcolormesh(X, Y, error, shading='nearest', cmap='hot')
    axes[1].set_title("FIM-error")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    axes[1].scatter(x_src, y_src, color='cyan', s=10, marker='.', label="Source")
    axes[1].legend()
    axes[1].axis('equal')
    plt.colorbar(im3, ax=axes[1])

    plt.tight_layout()
    plt.savefig("error_fim.png", dpi=300, bbox_inches="tight")
    plt.show()



if __name__ == "__main__":
    compute_error()