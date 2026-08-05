import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import re


def read_config_values():
    # Chargement de n, m, h et des sources dans l'ordre (x, y) du fichier de config
    n, m, h = None, None, 1.0
    sources_xy = []
    section = None
    if not os.path.exists("config.txt"):
        print("Erreur: Fichier config.txt pas trouvé, h=1 par défaut")
        return n, m, h, sources_xy

    with open("config.txt", "r") as f:
        for line in f:
            line = line.strip()
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


# Taille de la figure correspond à la "taille" des données (n_local*h en largeur et m*h en hauteur)
def compute_figsize(n_local, m, h, max_dim=10.0, min_dim=3.0):
    width = n_local * h
    height = m * h
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


def load_rank_metadata(rank):
    meta_file = f"result_rank{rank}_meta.txt"
    if not os.path.exists(meta_file):
        print(f"Erreur: fichier {meta_file} non trouvé")
        return None

    metadata = {}
    with open(meta_file, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or "=" not in line:
                continue
            key, value = line.split("=", 1)
            metadata[key.strip()] = value.strip()

    int_keys = {
        "rank", "n_global", "m_global", "overlap", "n_local", "top_ghost",
        "bottom_ghost", "i_start_overlap", "i_end_overlap", "i_owned_start",
        "i_owned_end", "n_owned"
    }
    for key in int_keys:
        if key in metadata:
            metadata[key] = int(metadata[key])
    if "h" in metadata:
        metadata["h"] = float(metadata["h"])

    return metadata


def load_rank_payloads():
    rank_files = sorted(glob.glob("result_rank*.txt"))
    payloads = []

    for rank_file in rank_files:
        match = re.fullmatch(r"result_rank(\d+)\.txt", os.path.basename(rank_file))
        if not match:
            continue

        rank = int(match.group(1))
        meta = load_rank_metadata(rank)
        if meta is None:
            continue

        try:
            matrix = np.loadtxt(rank_file)
        except Exception as e:
            print(f"Erreur: impossible de charger {rank_file}: {e}")
            continue

        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)

        matrix = matrix.astype(float).copy()
        matrix[matrix >= 1e300] = np.nan
        payloads.append({"rank": rank, "matrix": matrix, "meta": meta})

    return payloads


def collect_color_limits(payloads):
    finite_values = []
    for payload in payloads:
        matrix = payload["matrix"]
        if np.isfinite(matrix).any():
            finite_values.append(matrix[np.isfinite(matrix)])

    if not finite_values:
        return None, None

    merged = np.concatenate(finite_values)
    return float(np.nanmin(merged)), float(np.nanmax(merged))


def draw_overlap_bands(ax, meta):
    h = meta["h"]
    y_max = meta["m_global"] * h

    x_owned_min = meta["i_owned_start"] * h
    x_owned_max = (meta["i_owned_end"] + 1) * h
    x_overlap_min = meta["i_start_overlap"] * h
    x_overlap_max = (meta["i_end_overlap"] + 1) * h

    if meta["top_ghost"] > 0:
        ax.axvspan(x_overlap_min, x_owned_min, color="red", alpha=0.14, lw=0)
    if meta["bottom_ghost"] > 0:
        ax.axvspan(x_owned_max, x_overlap_max, color="red", alpha=0.14, lw=0)

    ax.axvline(x_owned_min, color="red", linestyle="--", linewidth=1.2)
    ax.axvline(x_owned_max, color="red", linestyle="--", linewidth=1.2)

    label_y = min(0.06 * y_max, y_max - 0.1 * h) if y_max > 0 else 0.0
    if meta["top_ghost"] > 0:
        ax.text(0.5 * (x_overlap_min + x_owned_min), label_y, "overlap", color="red", fontsize=8, ha="center")
    if meta["bottom_ghost"] > 0:
        ax.text(0.5 * (x_owned_max + x_overlap_max), label_y, "overlap", color="red", fontsize=8, ha="center")
    ax.text(x_owned_min + 0.02 * max(h, x_overlap_max - x_overlap_min), label_y, "owned", color="red", fontsize=8, ha="left")


def visualize_rank_payload(payload, all_sources, color_limits):
    rank = payload["rank"]
    meta = payload["meta"]
    matrix = payload["matrix"]

    h = meta["h"]
    n_local, m = matrix.shape
    x = (meta["i_start_overlap"] + np.arange(n_local)) * h
    y = np.arange(m) * h
    X, Y = np.meshgrid(x, y)

    figsize = compute_figsize(n_local, m, h)
    fig, ax = plt.subplots(figsize=figsize)
    vmin, vmax = color_limits
    mesh = ax.pcolormesh(X, Y, matrix.T, shading="nearest", cmap="viridis", vmin=vmin, vmax=vmax)
    ax.set_title(f"FIM MPI - sous-processus {rank}")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    if np.isfinite(matrix).any() and vmin is not None and vmax is not None and vmax > vmin:
        levels = np.linspace(vmin, vmax, 15)
        ax.contour(X, Y, matrix.T, levels, colors="white", linewidths=0.8, alpha=0.7)

    visible_sources = []
    x_min = meta["i_start_overlap"] * h
    x_max = (meta["i_end_overlap"] + 1) * h
    for (xs, ys) in all_sources:
        if x_min <= xs <= x_max:
            visible_sources.append((xs, ys))

    for (xs, ys) in visible_sources:
        ax.scatter(xs, ys, color="red", s=24, marker="o", edgecolors="white", linewidth=0.6, zorder=3, label="Source" if len(visible_sources) == 1 else "")
    if len(visible_sources) > 1:
        ax.scatter([], [], color="red", s=24, marker="o", label="Sources")
    if visible_sources:
        ax.legend(loc="upper right")

    draw_overlap_bands(ax, meta)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0.0, meta["m_global"] * h)
    ax.set_aspect("equal")
    fig.colorbar(mesh, ax=ax)

    outname = f"visualization_fim_rank{rank}.png"
    fig.savefig(outname, dpi=400, bbox_inches="tight")
    print(f"Figure enregistrée: {outname}")


def visualize_fim_mpi():
    payloads = load_rank_payloads()
    if not payloads:
        print("Erreur: aucun fichier result_rank*.txt trouvé")
        return

    _, _, h, sources_xy_indices = read_config_values()
    all_sources = sources_physical_coordinates(sources_xy_indices, h)
    color_limits = collect_color_limits(payloads)

    for payload in payloads:
        visualize_rank_payload(payload, all_sources, color_limits)

    plt.show()


if __name__ == "__main__":
    visualize_fim_mpi()