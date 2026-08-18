import numpy as np 
import matplotlib.pyplot as plt 
import os
from datetime import datetime

# Création d'un dossier avec la date et l'heure actuelle
def create_output_directory():
    now = datetime.now()
    dir_name = now.strftime("%Y-%m-%d_%H-%M-%S")
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def read_config(config_file):
    n = None
    m = None
    h = 1.0
    sources_ij = []
    section = None

    if not os.path.exists(config_file):
        print(f"Erreur: Fichier {config_file} pas trouvé, h=1 utilisé par défaut")
        return n, m, h, sources_ij

    with open(config_file, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line == "sources:":
                section = "sources"
                continue
            if line == "walls:":
                section = "walls"
                continue

            if section == "sources":
                parts = line.split()
                if len(parts) == 2:
                    try:
                        sources_ij.append((int(parts[0]), int(parts[1])))
                    except ValueError:
                        pass
                continue

            if "=" not in line:
                continue

            key, value = (part.strip() for part in line.split("=", 1))
            if key == "n":
                n = int(value)
            elif key == "m":
                m = int(value)
            elif key == "h":
                h = float(value)

    return n, m, h, sources_ij


def sources_physical_coordinates(sources_xy, h):
    return [(x * h, y * h) for (x, y) in sources_xy]

def visualize_fim_F(config_file="config.txt", output_dir=None):
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

    # Chargement de n, m, h et des sources dans l'ordre (x, y) du fichier de config
    n_matrix, m_matrix = matrix.shape
    n_cfg, m_cfg, h, sources_xy_indices = read_config(config_file)
    n = n_cfg if n_cfg is not None else n_matrix
    m = m_cfg if m_cfg is not None else m_matrix

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
 
    sources_xy = sources_physical_coordinates(sources_xy_indices, h)


    # Création de la grille
    x = np.arange(n) * h
    y = np.arange(m) * h
    X, Y = np.meshgrid(x,y)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    # Figure 1: solution de la FIM
    im1 = axes[0].pcolormesh(X, Y, T_fim.T, shading='nearest', cmap='viridis')
    axes[0].set_title("FIM - carte de distance globale")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    for (xs, ys) in sources_xy:
        axes[0].scatter(xs, ys, color='red', s=24, marker='o', label="Source" if len(sources_xy) == 1 else "", edgecolors='white', linewidth=0.6, zorder=3)
    if len(sources_xy) > 1:
        axes[0].scatter([], [], color='red', s=24, marker='o', label="Sources")

    if sources_xy:
        axes[0].legend(loc="upper right")
    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0])

    # Ajout des isocontours
    T_min = np.nanmin(T_fim)
    T_max = np.nanmax(T_fim)
    nb_contours = 15
    levels = np.linspace(T_min, T_max, nb_contours)
    axes[0].contour(X, Y, T_fim.T, levels, colors = 'white', linewidths=0.8, alpha=0.7)

    # Figure 2: Vitesse F
    im2 = axes[1].pcolormesh(X, Y, F.T, shading='nearest', cmap='viridis')
    axes[1].set_title("Vitesse F")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    for (xs, ys) in sources_xy:
        axes[1].scatter(xs, ys, color='red', s=24, marker='o', label="Source" if len(sources_xy) == 1 else "", edgecolors='white', linewidth=0.6, zorder=3)
    if len(sources_xy) > 1:
        axes[1].scatter([], [], color='red', s=24, marker='o', label="Sources")

    if sources_xy:
        axes[1].legend(loc="upper right")
    axes[1].legend()
    axes[1].axis('equal')
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    output_filename = os.path.join(output_dir, "visu_fim_speed.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    output_dir = create_output_directory()
    visualize_fim_F(output_dir=output_dir)