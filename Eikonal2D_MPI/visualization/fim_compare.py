import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

OUTPUT = "comparaison_fim.png"
FLOAT64_MAX = 1.7976931348623157e+308

def load_chkpt(path):
    with open(path, "rb") as f:
        header = np.fromfile(f, dtype=np.int32, count=2)
        n, m = int(header[0]), int(header[1])
        np.fromfile(f, dtype=np.float64, count=1)
        data = np.fromfile(f, dtype=np.float64, count=n * m)
    data = data.reshape((n, m))
    data = np.where(data == FLOAT64_MAX, 0, data)
    return data


def load_matrix(path):
    if path.endswith(".chkpt"):
        return load_chkpt(path)
    data = np.loadtxt(path)
    data = np.where(data == FLOAT64_MAX, 0, data)
    return data


def read_config(config_file="config.txt"):
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


# Création d'un dossier avec la date et l'heure actuelle
def create_output_directory():
    now = datetime.now()
    dir_name = now.strftime("%Y-%m-%d_%H-%M-%S")
    os.makedirs(dir_name, exist_ok=True)
    return dir_name


if len(sys.argv) != 3:
    print("Usage: fim_compare.py fichier1 fichier2")
    sys.exit(1)

file1, file2 = sys.argv[1], sys.argv[2]
label1, label2 = os.path.basename(file1), os.path.basename(file2)

a = load_matrix(file1)
b = load_matrix(file2)

if a.shape != b.shape:
    print(f"ERREUR: dimensions différentes -> {label1}={a.shape} {label2}={b.shape}")
    sys.exit(1)

# Lecture du fichier de configuration pour obtenir h, n, m
n_cfg, m_cfg, h, sources_xy_indices = read_config("config.txt")
n, m = a.shape

# Création de la grille avec les bonnes échelles
x = np.arange(n) * h
y = np.arange(m) * h
X, Y = np.meshgrid(x, y)

# Calcul de la différence et de l'erreur absolue
diff = b - a
abs_diff = np.abs(diff)
mask = np.isfinite(a) & np.isfinite(b)
max_abs_diff = np.abs(b[mask] - a[mask]).max(initial=0.0)

# Création du dossier de sortie
output_dir = create_output_directory()
output_path = os.path.join(output_dir, OUTPUT)

# Création de la figure avec seulement l'erreur absolue
fig, ax = plt.subplots(figsize=(7, 6))

# Affichage de l'erreur absolue avec les bonnes échelles X et Y
bound = max(max_abs_diff, 1e-15)
im = ax.pcolormesh(X, Y, abs_diff.T, shading='nearest', cmap='plasma', vmin=0, vmax=bound)
ax.set_title("Comparaison entre le séquentiel et le parallèle\n")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.axis('equal')
fig.colorbar(im, ax=ax)

fig.tight_layout()
fig.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"\nPlot enregistré: {output_path}")

plt.show()