import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt


def read_config_values():
    n, m, h = None, None, 1.0
    if not os.path.exists("config.txt"):
        print("Erreur: Fichier config.txt pas trouvé")
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


def load_rank_matrix(rank):
    fname = f"result_test_rank{rank}.txt"
    if not os.path.exists(fname):
        return None
    try:
        matrix = np.loadtxt(fname)
    except Exception as e:
        print(f"Erreur: impossible de lire {fname}: {e}")
        return None
    return np.array(matrix, dtype=float)


def build_global_from_local(rank_files, starts, overlap, n, m):
    global_grid = np.full((n, m), np.nan, dtype=float)
    coverage = np.zeros((n, m), dtype=int)

    for rank, fname in enumerate(rank_files):
        matrix = load_rank_matrix(rank)
        if matrix is None:
            continue

        i_start_owned, n_owned = starts[rank]
        top_ghost = min(overlap, i_start_owned) if rank > 0 else 0
        i_start_overlap = i_start_owned - top_ghost
        n_local = matrix.shape[0]

        for local_row in range(n_local):
            global_row = i_start_overlap + local_row
            if 0 <= global_row < n:
                coverage[global_row, :] += 1
                if i_start_owned <= global_row < i_start_owned + n_owned:
                    global_grid[global_row, :] = matrix[local_row, :]

    return global_grid, coverage


def write_summary(starts, overlap, nproc):
    lines = []
    for rank, (i_start_owned, n_owned) in enumerate(starts):
        top_ghost = min(overlap, i_start_owned) if rank > 0 else 0
        i_start_overlap = i_start_owned - top_ghost
        lines.append(
            f"rank {rank}: owned[{i_start_owned}, {i_start_owned + n_owned - 1}] overlap[{i_start_overlap}, {i_start_overlap + n_owned + top_ghost - 1}]"
        )
    with open("visualization_fim_mpi_debug_summary.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print("Résumé enregistré: visualization_fim_mpi_debug_summary.txt")


def main():
    n, m, h = read_config_values()
    if n is None or m is None:
        print("Impossible de lire la configuration, arrêt.")
        return

    rank_files = sorted(glob.glob("result_test_rank*.txt"))
    if not rank_files:
        print("Erreur: aucun fichier result_test_rank*.txt trouvé")
        return

    nproc = len(rank_files)
    overlap = 1
    if len(sys.argv) > 1:
        try:
            overlap = int(sys.argv[1])
        except ValueError:
            print("Attention: overlap invalide, valeur par défaut = 1")

    starts = compute_row_distribution(n, nproc)
    global_grid, coverage = build_global_from_local(rank_files, starts, overlap, n, m)
    write_summary(starts, overlap, nproc)

    x = np.arange(m) * h
    y = np.arange(n) * h
    X, Y = np.meshgrid(x, y)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)

    ax = axes[0]
    masked = np.ma.masked_invalid(global_grid)
    c0 = ax.pcolormesh(X, Y, masked, shading='nearest', cmap='viridis')
    ax.set_title("Vue globale (cellules owned only)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal')
    fig.colorbar(c0, ax=ax, label='T')

    for rank, (i_start_owned, n_owned) in enumerate(starts):
        ax.axhline(i_start_owned * h, color='white', linewidth=0.6, alpha=0.7)
        ax.axhline((i_start_owned + n_owned) * h, color='white', linewidth=0.6, alpha=0.7)

    ax = axes[1]
    c1 = ax.pcolormesh(X, Y, coverage, shading='nearest', cmap='magma')
    ax.set_title("Couverture par sous-domaines")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect('equal')
    fig.colorbar(c1, ax=ax, label='Nombre de sous-domaines couvrant la cellule')

    fig.savefig("visualization_fim_mpi_debug.png", dpi=300, bbox_inches='tight')
    print("Figure enregistrée: visualization_fim_mpi_debug.png")

    # Vue synthétique par rang
    fig2, axes2 = plt.subplots(1, nproc, figsize=(2.5 * nproc, 3.5), squeeze=False)
    for rank in range(nproc):
        matrix = load_rank_matrix(rank)
        if matrix is None:
            continue
        ax = axes2[0, rank]
        if matrix.ndim == 1:
            matrix = matrix.reshape(1, -1)
        img = ax.pcolormesh(matrix, shading='nearest', cmap='viridis')
        ax.set_title(f"Rang {rank}")
        ax.set_xlabel("col")
        ax.set_ylabel("ligne locale")
        ax.set_aspect('auto')
        fig2.colorbar(img, ax=ax, shrink=0.9)

    fig2.savefig("visualization_fim_mpi_ranks.png", dpi=300, bbox_inches='tight')
    print("Figure enregistrée: visualization_fim_mpi_ranks.png")

    plt.show()


if __name__ == "__main__":
    main()
