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
    T_fim = np.zeros(matrix.shape, dtype=float)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i,j] == 'inf' or matrix[i,j] == '1e+99':
                T_fim[i,j] = np.nan
            else:
                T_fim[i,j] = float(matrix[i,j])

    # chargement coords
    coords = np.loadtxt("coords_mesh.txt")
    if coords.ndim == 1:
        sources = [coords]
    else:
        sources = coords

    # Creating the mesh grid
    x = np.linspace(0, length, m)
    y = np.linspace(0, length, n)
    X, Y = np.meshgrid(x, y)

    # Visualization
    plt.figure(figsize=(10, 8))

    plt.pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    plt.colorbar()

    for (xs, ys) in sources:
        plt.scatter(xs, ys, color = 'red', s=80, marker='.', edgecolors='white', linewidth=1)
    plt.title("FIM-distance map")
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