import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import glob
from datetime import datetime

# Création d'un dossier avec la date et l'heure actuelle
def create_output_directory():
    now = datetime.now()
    dir_name = now.strftime("%Y-%m-%d_%H-%M-%S")
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def read_config(config_file="config.txt"):
    n = None
    m = None
    h = 1.0
    sources_xy = []
    section = None

    if not os.path.exists(config_file):
        return n, m, h, sources_xy

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
                        sources_xy.append((int(parts[0]), int(parts[1])))
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

    return n, m, h, sources_xy


def sources_physical_coordinates(sources_xy, h):
    return [(x * h, y * h) for (x, y) in sources_xy]


def visualize_checkpoint(output_dir=None):
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

    n_cfg, m_cfg, h_cfg, sources_xy_indices = read_config("config.txt")
    if n_cfg is not None:
        n = n_cfg
    if m_cfg is not None:
        m = m_cfg
    h = h_cfg
    sources = sources_physical_coordinates(sources_xy_indices, h)

    # Création de la grille avec x associe a la premiere coordonnee du fichier de config
    x = np.arange(n) * h
    y = np.arange(m) * h
    X, Y = np.meshgrid(x, y)

    # Visualisation
    fig, ax = plt.subplots(figsize=(10, 8))
    mesh = ax.pcolormesh(X, Y, T_fim.T, shading='nearest', cmap='viridis')
    ax.set_title("FIM MPI - carte de distance globale")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    # Ajout des isocontours
    T_min = np.nanmin(T_fim)
    T_max = np.nanmax(T_fim)
    nb_contours = 15
    levels = np.linspace(T_min, T_max, nb_contours)
    ax.contour(X, Y, T_fim.T, levels, colors='white', linewidths=0.8, alpha=0.7)

    # Sources
    for (xs, ys) in sources:
        ax.scatter(xs, ys, color='red', s=24, marker='o', label="Source" if len(sources) == 1 else "", edgecolors='red', linewidth=1)
    if len(sources) > 1:
        ax.scatter([], [], color='red', s=24, marker='o', label="Sources")
    ax.legend()
    ax.set_xlim(0.0, n * h)
    ax.set_ylim(0.0, m * h)
    ax.set_aspect('equal')
    fig.colorbar(mesh, ax=ax)

    output_filename = os.path.join(output_dir,"visualization_fim_mpi.png")

    # Enregistrement
    fig.savefig(output_filename, dpi=400, bbox_inches="tight")
    print(f"Figure enregistrée : {output_filename}")

    plt.show()


if __name__ == "__main__":
    output_dir = create_output_directory()
    visualize_checkpoint(output_dir=output_dir)