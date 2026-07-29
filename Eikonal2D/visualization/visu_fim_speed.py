import numpy as np 
import matplotlib.pyplot as plt 
import os

def visualize_fim_F():
    # Charge la matrice de la FIM
    fim_file = "matrix_fim.txt"
    if not os.path.exists(fim_file):
        print(f"Error: File {fim_file} not found")
        return

    try:
        matrix = np.loadtxt(fim_file, dtype=str)
    except Exception as e:
        print(f"Error: Unable to load {fim_file}: {e}")
        return

    # Retrouve n, m, and h
    n,m = matrix.shape

    h = 1   # Valeur par défaut
    if os.path.exists("config.txt"):
        with open("config.txt", "r") as f:
            for line in f:
                line = line.strip()
                if "h =" in line:
                    h = float(line.split("=")[1].strip())
    
    # Conversion en flottant
    T_fim = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            if matrix[i,j] == 'inf':
                T_fim[i,j] = np.nan
            else:
                T_fim[i,j] = float(matrix[i,j])


    # Chargement de la matrice de vitesse
    speed_file = "speed.txt"
    if not os.path.exists(speed_file):
        print(f"Erreur: Fichier {speed_file} pas trouvé")
        return

    try:
        F = np.loadtxt(speed_file)
    except Exception as e:
        print(f"Erreur: Impossible de charger {speed_file}: {e}")
        return

    # Chargement des coordonnées des sources
    coord_file = "coords_source.txt"
    sources=[]
    if os.path.exists(coord_file):
        try:
            coords = np.loadtxt(coord_file)
            if coords.ndim == 1:
                sources = [coords]
            else:
                sources = coords
        except Exception as e:
            print(f"Erreur: UImpossible de charger les sources: {e}")
    else:
        print(f"Erreur: Fichier {coord_file} pas trouvé")


    # Création de la grille
    x = np.arange(m) * h
    y = np.arange(n) * h
    X, Y = np.meshgrid(x,y)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Figure 1: solution de la FIM
    im1 = axes[0].pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    axes[0].set_title("FIM-carte de distance")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    for (xs, ys) in sources:
        axes[0].scatter(xs, ys, color='red', s=10, marker='.', edgecolors='red', linewidth=1, label='Source')
    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0])

    # Ajout des isocontours
    T_min = np.nanmin(T_fim)
    T_max = np.nanmax(T_fim)
    nb_contours = 15
    levels = np.linspace(T_min, T_max, nb_contours)
    axes[0].contour(X, Y, T_fim, levels, colors = 'white', linewidths=0.8, alpha=0.7)

    # Figure 2: Vitesse F
    im2 = axes[1].pcolormesh(X, Y, F, shading='nearest', cmap='viridis')
    axes[1].set_title("Vitesse F")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    for (xs, ys) in sources:
        axes[1].scatter(xs, ys, color='red', s=10, marker='.', edgecolors='red', linewidth=1, label='Source')
    axes[1].legend()
    axes[1].axis('equal')
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    plt.savefig("visu_fim_speed.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    visualize_fim_F()