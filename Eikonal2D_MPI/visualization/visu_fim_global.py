import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import glob


def visualize_checkpoint():
    # Recherche du fichier checkpoint
    checkpoint_files = glob.glob("*.chkpt")
    if not checkpoint_files:
        print("Erreur: aucun fichier .chkpt trouvé")
        return

    fim_file = checkpoint_files[0]

    # Chargement de n, m, h et des résultats de la fim
    try:
        with open(fim_file, "rb") as f:
            n = struct.unpack('i', f.read(4))[0]
            m = struct.unpack('i', f.read(4))[0]
            h = struct.unpack('d', f.read(8))[0]
            data = np.frombuffer(f.read(), dtype=np.float64, count=n * m)
    except Exception as e:
        print(f"Erreur: Impossible de charger {fim_file}: {e}")
        return

    T_fim = data.reshape(n, m).copy()
    T_fim[T_fim >= 1e300] = np.nan  # masque les cellules jamais atteintes

    # Chargement des coordonnées des sources
    coord_file = "coords_source.txt"
    sources = []
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
    X, Y = np.meshgrid(x, y)

    # Visualisation
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    plt.title("FIM MPI - carte de distance globale")
    plt.xlabel("X")
    plt.ylabel("Y")

    # Ajout des isocontours
    T_min = np.nanmin(T_fim)
    T_max = np.nanmax(T_fim)
    nb_contours = 15
    levels = np.linspace(T_min, T_max, nb_contours)
    plt.contour(X, Y, T_fim, levels, colors='white', linewidths=0.8, alpha=0.7)

    # Sources
    for (xs, ys) in sources:
        plt.scatter(xs, ys, color='red', s=10, marker='.', label="Source" if len(sources) == 1 else "", edgecolors='red', linewidth=1)
    if len(sources) > 1:
        plt.scatter([], [], color='red', s=10, marker='.', label="Sources")
    plt.legend()
    plt.axis('equal')
    plt.colorbar()

    # Enregistrement
    plt.savefig("visualization_fim_mpi.png", dpi=400, bbox_inches="tight")
    print("Figure enregistrée: visualization_fim_mpi.png")

    plt.show()


if __name__ == "__main__":
    visualize_checkpoint()