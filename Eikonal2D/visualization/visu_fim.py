import numpy as np 
import matplotlib.pyplot as plt 

def visualize_fim():
    n = 200
    m = 200
    length = 1.0

    # Loading the time matrix
    filename = "matrix_fim.txt"
    matrix = np.loadtxt(filename, dtype=str)

    # Converting 'inf' to np.inf
    matrix_float = np.zeros(matrix.shape, dtype=float)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i,j] == 'inf':
                matrix_float[i,j] = np.inf
            else:
                matrix_float[i,j] = float(matrix[i,j])

    # Validation
    if matrix_float.shape != (n,m):
        raise ValueError(f"Invalid shape: {matrix_float.shape}, expected: ({n}, {m})")

    # Creating the mesh grid
    x = np.linspace(0, length, m)
    y = np.linspace(0, length, n)
    X, Y = np.meshgrid(x, y)

    print(f"Matrix loaded: {matrix_float.shape}")

    # Loading source coordinates
    points_file = "coords_mesh.txt"
    try:
        coords = np.loadtxt(points_file)
        if coords.ndim == 1:
            x_points = [coords[0]]
            y_points = [coords[1]]
        else:
            x_points, y_points = coords[:, 0], coords[:, 1]
        print(f"{len(x_points)} loaded sources")
    except:
        print("Coordinates file not found.")
        x_points, y_points = [], []

    # Replace 'inf' with a high value for the visualization
    matrix_display = np.where(np.isinf(matrix_float), np.nan, matrix_float)

    # Visualization
    plt.figure(figsize=(10, 8))

    plt.pcolormesh(X, Y, matrix_display, shading='nearest', cmap='viridis')
    plt.colorbar(label="Arrival Time")
    plt.scatter(x_points, y_points, color='red', s=10, marker='*', label="Sources", edgecolors='black', linewidth=1)
    plt.title("FIM - Vertical wall in the centrer \n Arrival time")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend()
    plt.axis('equal')

    # Save
    plt.savefig("visualization_fim.png", dpi=400, bbox_inches="tight")
    print("Figure saved: visualization_fim.png")

    plt.show()



if __name__ == "__main__":
    visualize_fim()