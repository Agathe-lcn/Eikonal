import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import glob

def read_checkpoint(filename):
    """Lit un fichier checkpoint généré par FIMIO_Checkpoint"""
    try:
        with open(filename, 'rb') as f:
            # Lire l'en-tête : n (int), m (int), h (double)
            n = struct.unpack('i', f.read(4))[0]
            m = struct.unpack('i', f.read(4))[0]
            h = struct.unpack('d', f.read(8))[0]
            
            print(f"En-tête: n={n}, m={m}, h={h}")
            
            # Lire toutes les données restantes
            data_bytes = f.read()
            data = np.frombuffer(data_bytes, dtype=np.float64)
            
            expected_size = n * m
            if len(data) != expected_size:
                print(f"Attention: attendu {expected_size}, reçu {len(data)}")
                if len(data) >= expected_size:
                    data = data[:expected_size]
                else:
                    data = np.concatenate([data, np.full(expected_size - len(data), np.nan)])
            
            data = data.reshape(n, m)
            return n, m, h, data
            
    except Exception as e:
        print(f"Erreur: {e}")
        return None, None, None, None

def load_sources(coord_file="coords_source.txt", h=1.0):
    """Charge les sources depuis le fichier coords_source.txt"""
    sources = []
    if os.path.exists(coord_file):
        try:
            coords = np.loadtxt(coord_file)
            if coords.ndim == 1 and len(coords) == 2:
                sources = [coords]
            elif coords.ndim == 2:
                sources = coords
            print(f"Sources chargées: {len(sources)} sources")
            for s in sources:
                i = int(round(s[1] / h))
                j = int(round(s[0] / h))
                print(f"  Source: ({i}, {j})")
        except Exception as e:
            print(f"Erreur lors du chargement des sources: {e}")
    else:
        print("Fichier coords_source.txt non trouvé")
    return sources

def visualize_checkpoint():
    # ============================================================
    # 1. Chercher le fichier checkpoint
    # ============================================================
    checkpoint_files = glob.glob("*.chkpt")
    print(f"Fichiers .chkpt trouvés: {checkpoint_files}")
    
    if not checkpoint_files:
        print("Aucun fichier .chkpt trouvé")
        return
    
    # Prendre le premier fichier trouvé
    filename = checkpoint_files[0]
    print(f"\nLecture de: {filename}")
    
    # ============================================================
    # 2. Lire le checkpoint
    # ============================================================
    n, m, h, data = read_checkpoint(filename)
    
    if data is None:
        return
    
    # ============================================================
    # 3. Charger les sources depuis config.txt
    # ============================================================
    # Lire h depuis config.txt pour les sources
    h_source = h
    if os.path.exists("config.txt"):
        with open("config.txt", "r") as f:
            for line in f:
                line = line.strip()
                if "h =" in line:
                    try:
                        h_source = float(line.split("=")[1].strip())
                    except:
                        pass
    
    sources = load_sources("coords_source.txt", h_source)
    
    # ============================================================
    # 4. Nettoyer les données (inf -> NaN)
    # ============================================================
    data_clean = np.where(np.isinf(data) | (data > 1e30), np.nan, data)
    
    # Statistiques
    valid = data_clean[~np.isnan(data_clean)]
    print(f"\nStatistiques:")
    print(f"  Matrice: {data.shape}")
    print(f"  Cellules valides: {len(valid)} / {data.size}")
    if len(valid) > 0:
        print(f"  Min: {np.nanmin(data_clean):.4f}")
        print(f"  Max: {np.nanmax(data_clean):.4f}")
    else:
        print("  ⚠️ Aucune valeur valide")
    
    # ============================================================
    # 5. Visualisation de la grille complète
    # ============================================================
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Création de la grille
    x = np.arange(m) * h
    y = np.arange(n) * h
    X, Y = np.meshgrid(x, y)
    
    # Afficher la matrice
    if len(valid) > 0:
        vmin = np.nanmin(data_clean)
        vmax = np.nanmax(data_clean)
        if vmin == vmax:
            vmin = vmin - 1
            vmax = vmax + 1
        im = ax.pcolormesh(X, Y, data_clean, shading='nearest', cmap='viridis',
                           vmin=vmin, vmax=vmax)
    else:
        im = ax.pcolormesh(X, Y, data_clean, shading='nearest', cmap='viridis')
    
    ax.set_title(f"Grille complète - {n} x {m}", fontsize=14)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    
    # ============================================================
    # 6. Ajout des isocontours
    # ============================================================
    if len(valid) > 10:
        T_min = np.nanmin(data_clean)
        T_max = np.nanmax(data_clean)
        if T_min < T_max:
            nb_contours = 15
            levels = np.linspace(T_min, T_max, nb_contours)
            contour = ax.contour(X, Y, data_clean, levels, colors='white', 
                                linewidths=0.8, alpha=0.7)

    
    # ============================================================
    # 7. Ajout des sources
    # ============================================================
    for (xs, ys) in sources:
        # Convertir les coordonnées en indices
        i_global = int(round(ys / h))
        j_global = int(round(xs / h))
        
        # Vérifier que la source est dans la grille
        if 0 <= i_global < n and 0 <= j_global < m:
            # Afficher la source
            ax.scatter(xs, ys, color='red', s=80, marker='.', 
                      edgecolors='white', linewidth=1.5, zorder=5, 
                      label='Source' if len(sources) > 0 else '')
            print(f"Source affichée à ({xs:.2f}, {ys:.2f})")
        else:
            print(f"Source hors de la grille: ({xs:.2f}, {ys:.2f})")
    
    # Ajouter une légende pour les sources (si au moins une source est affichée)
    if len(sources) > 1:
        ax.scatter([], [], color='red', s=80, marker='.', 
                  edgecolors='white', linewidth=1.5, label='Sources')
    
    # Ajouter la légende
    if len(sources) > 0:
        ax.legend(loc='upper right')
    
    ax.axis('equal')
    
    # Ajouter la colorbar
    cbar = plt.colorbar(im, ax=ax, label='T')
    cbar.ax.tick_params(labelsize=10)
    
    # ============================================================
    # 8. Enregistrement et affichage
    # ============================================================
    output_file = "grid_complete.png"
    plt.savefig(output_file, dpi=150, bbox_inches="tight")
    print(f"\nFigure enregistrée: {output_file}")
    plt.show()

if __name__ == "__main__":
    visualize_checkpoint()