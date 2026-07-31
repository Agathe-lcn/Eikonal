import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys


def read_config_values():
    # Chargement de n, m et h
    n, m, h = None, None, 1.0
    if not os.path.exists("config.txt"):
        print(f"Erreur: Fichier config.txt pas trouvé, h=1 par défaut")
        return n, m, h

    with open("config.txt", "r") as f:
        for line in f:
            line = line.strip()
            if "=" not in line:
                continue
            key, value = (part.strip() for part in line.split("=", 1))
            if key == "n":
                n = int(value)
            elif key == "m":
                m = int(value)
            elif key == "h":
                h = float(value)

    return n, m, h


# Calcul de la distribution (Q,R) de la topologie pour retrouver la position globale de chaque rang
def compute_row_distribution(n, nproc):
    Q, R = divmod(n, nproc)
    starts = []

    for rank in range(nproc):
        if rank < R:
            n_owned = Q + 1
            i_start = rank * (Q + 1)
        else:
            n_owned = Q
            i_start = R * (Q + 1) + (rank - R) * Q
        starts.append((i_start, n_owned))

    return starts


# Chargement des coordonnées des sources
def load_sources():
    coord_file = "coords_source.txt"
    sources = []
    if os.path.exists(coord_file):
        try:
            coords = np.loadtxt(coord_file)
            if coords.ndim == 1:
                sources = [coords]
            elif coords.ndim == 2:
                sources = coords
        except Exception as e:
            print(f"Erreur: impossible de charger les sources: {e}")
    else:
        print(f"Attention: fichier {coord_file} non trouvé")

    return sources


# Taille de la figure correspond à la "taille" des données (m*h en largeur et n_local*h en hauteur)
def compute_figsize(n_local, m, h, max_dim=10.0, min_dim=3.0):
    width = m * h
    height = n_local * h
    if width <= 0 or height <= 0:
        return (max_dim, max_dim)

    if width >= height:
        width = max_dim
        height = max_dim * height / width
    else:
        height = max_dim
        width = max_dim * width / height

    width = max(width, min_dim)
    height = max(height, min_dim)
    
    return (width, height)


def visualize_fim_rank(rank, h, i_start_overlap, all_sources):
    # Chargement de la matrice FIM
    fim_file = f"result_test_rank{rank}.txt"
    if not os.path.exists(fim_file):
        print(f"Erreur: fichier {fim_file} non trouvé")
        return

    try:
        matrix = np.loadtxt(fim_file)
    except Exception as e:
        print(f"Erreur: impossible de charger {fim_file}: {e}")
        return

    n_local, m = matrix.shape

    T_fim = matrix.astype(float).copy()
    T_fim[T_fim >= 1e300] = np.nan   # masque les cellules jamais atteintes

    # Grille locale : x et y démarrent à 0
    x = np.arange(m) * h
    y = np.arange(n_local) * h
    X, Y = np.meshgrid(x, y)

    figsize = compute_figsize(n_local, m, h)
    plt.figure(figsize=figsize)
    plt.pcolormesh(X, Y, T_fim, shading='nearest', cmap='viridis')
    plt.title(f"FIM MPI — Processus {rank} (grille locale {n_local}x{m})")
    plt.xlabel("X")
    plt.ylabel("Y (local, recouvrement inclus)")

    # Ajout des isocontours
    if np.isfinite(T_fim).any():
        T_min = np.nanmin(T_fim)
        T_max = np.nanmax(T_fim)
        nb_contours = 15
        levels = np.linspace(T_min, T_max, nb_contours)
        plt.contour(X, Y, T_fim, levels, colors = 'white', linewidths=0.8, alpha=0.7)

    # Sources
    local_sources = []
    if i_start_overlap is not None:
        for (xs, ys) in all_sources:
            i_source = ys / h
            if i_start_overlap <= i_source < i_start_overlap + n_local:
                local_sources.append((xs, ys - i_start_overlap * h))
    else:
        local_sources = all_sources

    for (xs, ys) in local_sources:
        plt.scatter(xs, ys, color='red', s=10, marker='.', label="Source" if len(local_sources) == 1 else "", edgecolors='red', linewidth=1)
    if len(local_sources) > 1:
        plt.scatter([], [], color='red', s=10, marker='.', label="Sources")
    if local_sources:
        plt.legend()

    plt.axis('equal')
    plt.colorbar()

    outname = f"visualization_fim_rank{rank}.png"
    plt.savefig(outname, dpi=400, bbox_inches="tight")
    print(f"Figure enregistrée: {outname}")


def visualize_fim_mpi():
    rank_files = glob.glob("result_test_rank*.txt")
    if not rank_files:
        print("Erreur: aucun fichier result_test_rank*.txt trouvé")
        return

    ranks = sorted(int(f.replace("result_test_rank", "").replace(".txt", "")) for f in rank_files)
    nproc = max(ranks) + 1

    # Chargement des coordonnées des sources
    n, m, h = read_config_values()

    # Overlap passé en argument
    overlap = None
    if len(sys.argv) > 1:
        try:
            overlap = int(sys.argv[1])
        except ValueError:
            print("Attention: overlap invalide, positionnement des sources désactivé")

    if (overlap is not None and n):
        starts = compute_row_distribution(n, nproc) 
    else: 
        starts = None

    all_sources = load_sources()

    for rank in ranks:
        i_start_overlap = None
        if starts is not None:
            i_start_owned, _ = starts[rank]
            top_ghost = overlap if rank > 0 else 0
            top_ghost = min(top_ghost, i_start_owned)
            i_start_overlap = i_start_owned - top_ghost
            
        visualize_fim_rank(rank, h, i_start_overlap, all_sources)

    plt.show()


if __name__ == "__main__":
    visualize_fim_mpi()