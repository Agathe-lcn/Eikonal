import numpy as np 
import matplotlib.pyplot as plt 
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors

def visualize_fim():
    # Chargement de la matrice FIM
    fim_file = "matrix_fim.txt"
    if not os.path.exists(fim_file):
        print(f"Erreur: Fichier {fim_file} pas trouvé")
        return

    try:
        matrix = np.loadtxt(fim_file, dtype=str)
    except Exception as e:
        print(f"Erreur: Impossible de charger {fim_file}: {e}")
        return

    # Chargement de n, m, et h
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
            print(f"Erreur: Impossible de charger les sources {e}")
    else:
        print(f"Erreur: Fichier {coord_file} pas trouvé")


    # Création de la grille
    x = np.arange(m) * h
    y = np.arange(n) * h
    X, Y = np.meshgrid(x,y)

    # Visualisation
    plt.figure(figsize=(10,8))
    plt.pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    plt.title("FIM-carte de distance")
    plt.xlabel("X")
    plt.ylabel("Y")

    # Ajout des isocontours
    T_min = np.nanmin(T_fim)
    T_max = np.nanmax(T_fim)
    nb_contours = 15
    levels = np.linspace(T_min, T_max, nb_contours)
    plt.contour(X, Y, T_fim, levels, colors = 'white', linewidths=0.8, alpha=0.7)

    # Sources
    for (xs,ys) in sources:
        plt.scatter(xs, ys, color='red', s=10, marker='.', label="Source" if len(sources) == 1 else "", edgecolors='red', linewidth=1)
    if len(sources) > 1:
        plt.scatter([], [], color='red', s=10, marker='.', label="Sources")
    
    plt.legend()
    plt.axis('equal')
    plt.colorbar()

    # Enregistrement
    plt.savefig("visualization_fim.png", dpi=400, bbox_inches="tight")
    print("Figure enregistrée: visualization_fim.png")

    plt.show()



if __name__ == "__main__":
    visualize_fim()