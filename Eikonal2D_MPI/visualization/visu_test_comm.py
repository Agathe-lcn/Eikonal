import numpy as np
import matplotlib.pyplot as plt
import struct
import os

def load_chkpt_global(filename):
    """Charge le checkpoint binaire de la vue globale."""
    if not os.path.exists(filename):
        return None, 0, 0, 0
    with open(filename, 'rb') as f:
        n, m = struct.unpack('ii', f.read(8))
        h = struct.unpack('d', f.read(8))[0]
        data = struct.unpack(f'{n*m}d', f.read(n*m*8))
        matrix = np.array(data).reshape((n, m))
        matrix[matrix > 1e10] = np.nan
        return matrix, n, m, h

def load_local_txt(filename):
    """Charge la mémoire brute d'un processus."""
    if not os.path.exists(filename):
        return None
    with open(filename, 'r') as f:
        lines = f.readlines()
        if not lines:
            return None
        n, m = map(int, lines[0].split())
        data = []
        for line in lines[1:]:
            row = [np.nan if val == "INF" else float(val) for val in line.split()]
            if row:
                data.append(row)
        return np.array(data)

steps = [
    ("Étape 1 : Avant Comm 1", "step1_avant_comm1"),
    ("Étape 2 : Après Comm 1", "step2_apres_comm1"),
    ("Étape 3 : Synchro Cohérente", "step3_synchro_coherente"),
    ("Étape 4 : Cycle 2 (Avant Comm 2)", "step4_avant_comm2")
]

OVERLAP = 2  # Bande de recouvrement de 2 pixels

fig, axes = plt.subplots(4, 3, figsize=(18, 19))

for i, (title_step, prefix) in enumerate(steps):
    glob_matrix, n, m, h = load_chkpt_global(f"{prefix}.chkpt")
    p0_matrix = load_local_txt(f"{prefix}_proc0.txt")
    p1_matrix = load_local_txt(f"{prefix}_proc1.txt")

    ## Vue globale
    ax_g = axes[i, 0]
    if glob_matrix is not None:
        im_g = ax_g.imshow(glob_matrix, cmap='viridis', origin='lower', extent=[0, m*h, 0, n*h])
        ax_g.set_title(f"{title_step}", fontsize=9, fontweight='bold', color='navy')
        ax_g.set_ylabel("Y (Global)")
        ax_g.axhline(y=(n//2)*h, color='red', linestyle='--', linewidth=1.5, label="Interface MPI")
        plt.colorbar(im_g, ax=ax_g, fraction=0.046, pad=0.04)
    else:
        ax_g.text(0.5, 0.5, "Global indisponible", ha='center', va='center')

    ## Vue locale proc 0
    ax_p0 = axes[i, 1]
    if p0_matrix is not None:
        im_p0 = ax_p0.imshow(p0_matrix, cmap='viridis', origin='lower')
        ax_p0.set_title(f"{title_step} [Proc 0 - BAS]", fontsize=9, fontweight='bold', color='darkgreen')
        
        n_p0 = p0_matrix.shape[0]
        # Frontière 1 : Séparation Domaine Possédé / Ghost Cells (Haut)
        b1_p0 = n_p0 - OVERLAP - 0.5
        # Frontière 2 : Limite de la bande de 2px possédée et accessible par Proc 1
        b2_p0 = n_p0 - 2 * OVERLAP - 0.5
        
        ax_p0.axhline(y=b1_p0, color='red', linestyle='-', linewidth=1.5)
        ax_p0.axhline(y=b2_p0, color='orange', linestyle='--', linewidth=1.5)
        
        # Annotations textuelles
        ax_p0.text(0.5, b1_p0 + 0.3, "▲ Ghost Cells (Reçues du Proc 1)", color='red', fontsize=7, fontweight='bold')
        ax_p0.text(0.5, b1_p0 - 0.7, "▼ Overlap Possédé (Envoyé au Proc 1)", color='orange', fontsize=7, fontweight='bold')
        ax_p0.text(0.5, b2_p0 - 0.7, "▼ Domaine Intérieur Propre", color='gray', fontsize=7)
        
        plt.colorbar(im_p0, ax=ax_p0, fraction=0.046, pad=0.04)
    else:
        ax_p0.text(0.5, 0.5, "Proc 0 TXT non trouvé", ha='center', va='center')

    # Vue locale proc 1
    ax_p1 = axes[i, 2]
    if p1_matrix is not None:
        im_p1 = ax_p1.imshow(p1_matrix, cmap='viridis', origin='lower')
        ax_p1.set_title(f"{title_step} [Proc 1 - HAUT]", fontsize=9, fontweight='bold', color='darkred')
        
        # Frontière 1 : Séparation Ghost Cells (Bas) / Domaine Possédé
        b1_p1 = OVERLAP - 0.5
        # Frontière 2 : Limite de la bande de 2px possédée et accessible par Proc 0
        b2_p1 = 2 * OVERLAP - 0.5
        
        ax_p1.axhline(y=b1_p1, color='red', linestyle='-', linewidth=1.5)
        ax_p1.axhline(y=b2_p1, color='orange', linestyle='--', linewidth=1.5)
        
        # Annotations textuelles
        ax_p1.text(0.5, b1_p1 - 0.7, "▼ Ghost Cells (Reçues du Proc 0)", color='red', fontsize=7, fontweight='bold')
        ax_p1.text(0.5, b1_p1 + 0.3, "▲ Overlap Possédé (Envoyé au Proc 0)", color='orange', fontsize=7, fontweight='bold')
        ax_p1.text(0.5, b2_p1 + 0.3, "▲ Domaine Intérieur Propre", color='gray', fontsize=7)
        
        plt.colorbar(im_p1, ax=ax_p1, fraction=0.046, pad=0.04)
    else:
        ax_p1.text(0.5, 0.5, "Proc 1 TXT non trouvé", ha='center', va='center')

    if i == 3:
        ax_g.set_xlabel("X")
        ax_p0.set_xlabel("X (Local)")
        ax_p1.set_xlabel("X (Local)")

plt.tight_layout()
plt.suptitle("Analyse FIM MPI : Décomposition Complète des Zones de Recouvrement (Possédées & Fantômes)", fontsize=12, fontweight='bold', y=1.01)

output_file = "visu_4steps_complete_overlap.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"Graphique enregistré sous : {output_file}")
plt.show()