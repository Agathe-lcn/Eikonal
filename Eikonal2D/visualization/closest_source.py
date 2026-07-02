import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os

def visualize_tags():
    # Loading tags
    filename = "source_tags.txt"
    fim_tags = np.loadtxt(filename, dtype=int)
    n, m = fim_tags.shape

    # Reading h from config.txt
    h = 1.0
    if os.path.exists("config.txt"):
        with open("config.txt", "r") as f:
            for line in f:
                line = line.strip()
                if "h =" in line:
                    h = float(line.split("=")[1].strip())

    # Loading source coordinates
    coord_file = "coords_source.txt"
    if not os.path.exists(coord_file):
        raise FileNotFoundError(f"{coord_file} not found")
    sources = np.loadtxt(coord_file)
    if sources.ndim == 1:
        sources = sources.reshape(1, -1)
    n_sources = sources.shape[0]

    # Grid for visualization
    x = np.arange(m) * h
    y = np.arange(n) * h
    X, Y = np.meshgrid(x, y)

    # Euclidean tag: nearest source for each grid point
    euclid_tags = np.zeros((n, m), dtype=int)
    for i in range(n):
        yi = i * h
        for j in range(m):
            xj = j * h
            min_dist = None
            min_s = 0
            for s in range(n_sources):
                xs, ys = sources[s, 0], sources[s, 1]
                dx = xj - xs
                dy = yi - ys
                dist = np.sqrt(dx * dx + dy * dy)
                if min_dist is None or dist < min_dist:
                    min_dist = dist
                    min_s = s
            euclid_tags[i, j] = min_s

    # Mismatch map
    mismatch = (fim_tags != euclid_tags).astype(int)

    # Color palette for tags
    colors = ['#00FF00', '#0000FF', '#FFA500', '#800080', '#00FFFF', '#FF1493', '#008000', '#FFD700', '#4B0082', '#FF69B4', '#00FF7F', '#FF4500', '#1E90FF', '#FF00FF']
    cmap_colors = ['black']
    for s in range(n_sources):
        cmap_colors.append(colors[s % len(colors)])
    cmap_tags = ListedColormap(cmap_colors)
    norm_tags = plt.Normalize(vmin=-1, vmax=n_sources - 0.5)

    cmap_compare = ListedColormap(['black', 'red'])

    # Creating the figure with 2 subplots
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Figure 1: FIM tag map
    im1 = axes[0].pcolormesh(X, Y, fim_tags, shading='nearest', cmap=cmap_tags, norm=norm_tags)
    axes[0].set_title("Tag of the closest source (FIM)")
    axes[0].set_xlabel("X")
    axes[0].set_ylabel("Y")
    axes[0].scatter(sources[:, 0], sources[:, 1], color='red', s=10, marker='.', label="Sources")
    axes[0].legend()
    axes[0].axis('equal')
    plt.colorbar(im1, ax=axes[0], label="Source ID")

    # Figure 2: Comparison FIM vs Euclidean
    im2 = axes[1].pcolormesh(X, Y, mismatch, shading='nearest', cmap=cmap_compare, vmin=0, vmax=1)
    axes[1].set_title("Comparison: FIM vs Euclidean nearest source")
    axes[1].set_xlabel("X")
    axes[1].set_ylabel("Y")
    axes[1].scatter(sources[:, 0], sources[:, 1], color='cyan', s=10, marker='.', label="Sources")
    axes[1].legend()
    axes[1].axis('equal')
    cbar2 = plt.colorbar(im2, ax=axes[1], ticks=[0, 1])
    cbar2.ax.set_yticklabels(['Match', 'Mismatch'])

    plt.tight_layout()
    plt.savefig("tags_comparison.png", dpi=300, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    visualize_tags()