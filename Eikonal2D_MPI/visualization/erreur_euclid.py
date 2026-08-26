import numpy as np
import matplotlib.pyplot as plt 
import os
from datetime import datetime

FLOAT64_MAX = 1.7976931348623157e+308

# Création d'un dossier avec la date et l'heure actuelle
def create_output_directory():
    now = datetime.now()
    dir_name = now.strftime("%Y-%m-%d_%H-%M-%S")
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def load_chkpt(path):
    with open(path, "rb") as f:
        header = np.fromfile(f, dtype=np.int32, count=2)
        n, m = int(header[0]), int(header[1])
        np.fromfile(f, dtype=np.float64, count=1)
        data = np.fromfile(f, dtype=np.float64, count=n * m)
    data = data.reshape((n, m))
    data = np.where(data == FLOAT64_MAX, 0, data)
    return data


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


def compute_error(config_file="config.txt", output_dir=None):
    # Chargement de la matrice FIM
    fim_file = "matrix_fim.txt.chkpt"
    if not os.path.exists(fim_file):
        print(f"Erreur: Fichier {fim_file} pas trouvé")
        return

    try:
        matrix = load_chkpt(fim_file)
    except Exception as e:
        print(f"Erreur: Impossible de charger {fim_file}: {e}")
        return

    # Chargement de n, m, h et des sources dans l'ordre (x, y) du fichier de config
    n_matrix, m_matrix = matrix.shape
    n_cfg, m_cfg, h, sources_xy_indices = read_config(config_file)
    n = n_cfg if n_cfg is not None else n_matrix
    m = m_cfg if m_cfg is not None else m_matrix

    # Conversion en flottant
    T_fim = np.zeros((n_matrix, m_matrix))
    for i in range(n_matrix):
        for j in range(m_matrix):
            if matrix[i, j] == 'inf':
                T_fim[i, j] = np.nan
            else:
                T_fim[i, j] = float(matrix[i, j])

    sources_xy = sources_physical_coordinates(sources_xy_indices, h)

    # Création de la grille avec x associe a la premiere coordonnee du fichier de config
    x = np.arange(n) * h
    y = np.arange(m) * h
    X, Y = np.meshgrid(x, y)

    # Distance euclidienne
    T_euclid = np.zeros((n,m))
    for i in range(n):
        y = i * h 
        for j in range(m):
            x = j * h
            # Calcule la distance par rapport à chaque source et garde la plus petite
            min_dist = np.inf
            for (xs, ys) in sources_xy:
                dx = x - xs
                dy = y - ys
                dist = np.sqrt(dx*dx + dy*dy)
                if dist < min_dist:
                    min_dist = dist
            T_euclid[i,j] = min_dist

    # Erreur = | T_euclid - T_fim |
    error = np.abs(T_euclid - T_fim)
    
    # Statistiques
    valid_error = error[~np.isnan(error)]
    print("\nStatistiques sur l'erreur |T_euclid - T_fim|\n")
    print(f"Erreur maximale: {np.max(valid_error):.6e}")
    print(f"Erreur minimale: {np.min(valid_error):.6e}")
    print(f"Erreur moyenne: {np.mean(valid_error):.6e}")
    print(f"Erreur médiane: {np.median(valid_error):.6e}")

    # Création de la figure avec 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Figure 1: solution de la FIM
    im1 = axes[0].pcolormesh(X, Y, T_fim.T, shading='nearest', cmap='viridis')
    axes[0].set_title("FIM - carte de distance globale")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    # Ajout des isocontours
    T_min = np.nanmin(T_fim)
    T_max = np.nanmax(T_fim)
    nb_contours = 15
    levels = np.linspace(T_min, T_max, nb_contours)
    axes[0].contour(X, Y, T_fim.T, levels, colors = 'white', linewidths=0.8, alpha=0.7)
    # Affiche toutes les sources
    for (xs, ys) in sources_xy:
        axes[0].scatter(xs, ys, color='red', s=24, marker='o', label="Source" if len(sources_xy) == 1 else "", edgecolors='white', linewidth=0.6, zorder=3)
    if len(sources_xy) > 1:
        axes[0].scatter([], [], color='red', s=24, marker='o', label="Sources")
    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0])

    # Figure 2: Erreur | T_euclid - T_fim |
    im2 = axes[1].pcolormesh(X, Y, error.T, shading='nearest', cmap='hot')
    axes[1].set_title("FIM - erreur")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    # Affiche toutes les sources
    for (xs, ys) in sources_xy:
        axes[1].scatter(xs, ys, color='cyan', s=15, marker='o', label="Source" if len(sources_xy) == 1 else "")
    if len(sources_xy) > 1:
        axes[1].scatter([], [], color='cyan', s=24, marker='o', label="Sources")
    axes[1].legend()
    axes[1].axis('equal')
    plt.colorbar(im2, ax=axes[1])

    plt.tight_layout()
    output_filename = os.path.join(output_dir, "error_fim.png")
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    output_dir = create_output_directory()
    compute_error(output_dir=output_dir)