import numpy as np
import matplotlib.pyplot as plt 
import os

def compute_error():
    # Paramètres
    n = 200
    m = 200
    length = 1.0
    h = length / n

    # Creation de la figure avec 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    # Chargement de la matrice calculée par la FIM
    filename = "matrix_fim.txt"
    matrix = np.loadtxt(filename, dtype=str)

    # Conversion en flottant
    T_fim = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            if matrix[i,j] == 'inf' or matrix[i,j] == '1e+99':
                T_fim[i,j] = np.nan
            else:
                T_fim[i,j] = float(matrix[i,j])

    # Chargement des coordonnées des sources
    coord_file = "coords_source.txt"
    if not os.path.exists(coord_file):
        raise FileNotFoundError(f"{coord_file} pas trouvé")

    # Chargement de toutes les sources
    coords = np.loadtxt(coord_file)
    if coords.ndim == 1:
        # Source unique
        sources = [coords]
    else:
        # Plusieurs sources
        sources = coords

    # Distance euclidienne
    T_euclid = np.zeros((n,m))
    for i in range(n):
        y = i * h 
        for j in range(m):
            x = j * h
            # Calcule la distance par rapport à chaque source et garde la plus petite
            min_dist = np.inf
            for (x_src, y_src) in sources:
                dx = x - x_src
                dy = y - y_src
                dist = np.sqrt(dx*dx + dy*dy)
                if dist < min_dist:
                    min_dist = dist
            T_euclid[i,j] = min_dist

    # Erreur = | T_euclid - T_fim |
    error = np.abs(T_euclid - T_fim)
    
    # Statistiques
    valid_error = error[~np.isnan(error)]
    print("\nError statistics |T_euclid - T_fim|\n")
    print(f"Maximum error: {np.max(valid_error):.6e}")
    print(f"Minimal error: {np.min(valid_error):.6e}")
    print(f"Average error: {np.mean(valid_error):.6e}")
    print(f"Median error: {np.median(valid_error):.6e}")

    # Visualisation
    x = np.linspace(0, length, m)
    y = np.linspace(0, length, n)
    X, Y = np.meshgrid(x,y)

    # Figure 1: solution de la FIM
    im1 = axes[0].pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    axes[0].set_title("FIM-distance map")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    # Affiche toutes les sources
    for (x_src, y_src) in sources:
        axes[0].scatter(x_src, y_src, color='red', s=10, marker='.', label="Source" if len(sources) == 1 else "")
    if len(sources) > 1:
        axes[0].scatter([], [], color='red', s=10, marker='.', label="Sources")
    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0])

    # Figure 2: Erreur | T_euclid - T_fim |
    im2 = axes[1].pcolormesh(X, Y, error, shading='nearest', cmap='hot')
    axes[1].set_title("FIM-error")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    # Affiche toutes les sources
    for (x_src, y_src) in sources:
        axes[1].scatter(x_src, y_src, color='cyan', s=10, marker='.', label="Source" if len(sources) == 1 else "")
    if len(sources) > 1:
        axes[1].scatter([], [], color='cyan', s=10, marker='.', label="Sources")
    axes[1].legend()
    axes[1].axis('equal')
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    plt.savefig("error_fim.png", dpi=300, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    compute_error()