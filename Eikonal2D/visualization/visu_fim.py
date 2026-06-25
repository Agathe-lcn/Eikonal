import numpy as np 
import matplotlib.pyplot as plt 
import os

def visualize_fim():
    # Loading the FIM matrix
    fim_file = "matrix_fim.txt"
    if not os.path.exists(fim_file):
        print(f"Error: File {fim_file} not found")
        return

    try:
        matrix = np.loadtxt(fim_file, dtype=str)
    except Exception as e:
        print(f"Error: Unable to load {fim_file}: {e}")
        return

    #retrieving n, m, and h
    n,m = matrix.shape

    if os.path.exists("config.txt"):
        with open("config.txt", "r") as f:
            for line in f:
                line = line.strip()
                if "h =" in line:
                    h = float(line.split("=")[1].strip())

    # Convert to float
    T_fim = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            if matrix[i,j] == 'inf':
                T_fim[i,j] = np.inf
            else:
                T_fim[i,j] = float(matrix[i,j])


    # Loading source coordinates
    coord_file = "coords_source.txt"
    sources=[]
    if os.path.exists(coord_file):
        try:
            coords = np.loadtxt(coord_file)
            if coords.ndim == 1:
                sources = [coords]
            else:
                sources = coords
        except:
            print(f"Error: Unable to load the sources: {e}")
    else:
        print(f"Error: File {coord_file} not found")


    # Creating the mesh grid
    x = np.arange(m)
    y = np.arange(n)
    X, Y = np.meshgrid(x,y)

    # Visualization
    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    plt.title("FIM-distance map")
    plt.xlabel("X")
    plt.ylabel("Y")

    # Sources
    for (xs,ys) in sources:
        # Convert real coordinates to indices
        if m > 1:
            j_src = int(round(xs/h))
        else:
            j_src = 0

        if n > 1:
            i_src = int(round(ys/h))
        else:
            i_src = 0

        # Check the limits
        j_src = max(0, min(j_src, m-1))
        i_src = max(0, min(i_src, n-1))
        plt.scatter(j_src, i_src, color='red', s=10, marker='.', edgecolors='red', linewidth=1)

    plt.legend()
    plt.axis('equal')
    plt.colorbar()

    # Save
    plt.savefig("visualization_fim.png", dpi=400, bbox_inches="tight")
    print("Figure saved: visualization_fim.png")

    plt.show()



if __name__ == "__main__":
    visualize_fim()