import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os


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

def visualize_tags(config_file="config.txt"):
    # Chargement des tags
    tags_file = "source_tags.txt"
    if not os.path.exists(tags_file):
        print(f"Erreur: Fichier {tags_file} pas trouvé")
        return
    fim_tags = np.loadtxt(tags_file, dtype=int)

    # Chargement de n, m, h et des sources dans l'ordre (x, y) du fichier de config
    n_tag, m_tag = fim_tags.shape
    n_cfg, m_cfg, h, sources_xy_indices = read_config(config_file)
    n = n_cfg if n_cfg is not None else n_tag
    m = m_cfg if m_cfg is not None else m_tag

    # Chargement des coordonnées des sources
    sources_xy = sources_physical_coordinates(sources_xy_indices, h)

    # Création de la grille
    x = np.arange(n) * h
    y = np.arange(m) * h
    X, Y = np.meshgrid(x, y)

    # Tag distance euclidienne: source la plus proche pour chaque point de la grille
    euclid_tags = np.zeros((n, m), dtype=int)
    for i in range(n):
        xi = i * h
        for j in range(m):
            yj = j * h
            min_dist = None
            min_s = 0
            for s in range(len(sources_xy)):
                xs, ys = sources_xy[s][0], sources_xy[s][1]
                dx = xi - xs
                dy = yj - ys
                dist = np.sqrt(dx * dx + dy * dy)
                if min_dist is None or dist < min_dist:
                    min_dist = dist
                    min_s = s
            euclid_tags[i, j] = min_s

    # Carte de l'erreur
    mismatch = (fim_tags != euclid_tags).astype(int)

    # Palette de coulours pour les tags
    colors = ['#00FF00', '#0000FF', '#FFA500', '#800080', '#00FFFF', '#FF1493', '#008000', '#FFD700', '#4B0082', '#FF69B4', '#00FF7F', '#FF4500', '#1E90FF', '#FF00FF']
    cmap_colors = ['black']
    for s in range(len(sources_xy)):
        cmap_colors.append(colors[s % len(colors)])
    cmap_tags = ListedColormap(cmap_colors)
    norm_tags = plt.Normalize(vmin=-1, vmax=len(sources_xy) - 0.5)

    cmap_compare = ListedColormap(['black', 'red'])

    # Création de la figure avec 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Figure 1: carte des tags de la FIM
    im1 = axes[0].pcolormesh(X, Y, fim_tags.T, shading='nearest', cmap=cmap_tags, norm=norm_tags)
    axes[0].set_title("Tag de la source la plus proche (FIM)")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")

    for (xs, ys) in sources_xy:
        axes[0].scatter(xs, ys, color='red', s=24, marker='o', label="Source" if len(sources_xy) == 1 else "", edgecolors='white', linewidth=0.6, zorder=3)
    if len(sources_xy) > 1:
        axes[0].scatter([], [], color='red', s=24, marker='o', label="Sources")

    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0], label="Source ID")

    # Figure 2: Comparison FIM vs Euclidean
    im2 = axes[1].pcolormesh(X, Y, mismatch.T, shading='nearest', cmap=cmap_compare, vmin=0, vmax=1)
    axes[1].set_title("Comparaison: Source la plus proche, FIM vs euclidienne")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")

    for (xs, ys) in sources_xy:
        axes[1].scatter(xs, ys, color='cyan', s=24, marker='o', label="Source" if len(sources_xy) == 1 else "", edgecolors='white', linewidth=0.6, zorder=3)
    if len(sources_xy) > 1:
        axes[1].scatter([], [], color='cyan', s=24, marker='o', label="Sources")
    
    axes[1].legend()
    axes[1].axis('equal')
    cbar2 = plt.colorbar(im2, ax=axes[1], ticks=[0, 1])
    cbar2.ax.set_yticklabels(['Match', 'Mismatch'])

    plt.tight_layout()
    plt.savefig("tags_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    visualize_tags()