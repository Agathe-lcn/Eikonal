import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os 
from datetime import datetime
import io

# Création d'un dossier avec la date et l'heure actuelle
def create_output_directory():
    now = datetime.now()
    dir_name = now.strftime("%Y-%m-%d_%H-%M-%S")
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def plot_benchmark_time(data, output_dir):
    # Calcul des moyennes par overlap
    grouped = data.groupby('overlap').agg({'temps_sec': 'mean'}).reset_index()
    
    overlaps = grouped['overlap'].values
    avg_temps = grouped['temps_sec'].values
    
    # Création de la figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Courbe : Temps de calcul
    ax.set_xlabel('Overlap', fontsize=12)
    ax.set_ylabel('Temps moyen (secondes)', fontsize=12)
    ax.plot(overlaps, avg_temps, 'o-', color='blue', linewidth=2, markersize=8, label='Temps de calcul')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)
    
    # Ajout des valeurs sur les points
    for i, (x, y) in enumerate(zip(overlaps, avg_temps)):
        ax.annotate(f'{y:.1f}s', (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=9)
    
    plt.title('Temps de calcul en fonction de l\'overlap', fontsize=14)
    
    # Ajustement et sauvegarde
    fig.tight_layout()
    output_filename = os.path.join(output_dir, 'benchmark_time.png')
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.show()
    
    return output_filename

def plot_benchmark_communications(data, output_dir):
    # Calcul des moyennes par overlap
    grouped = data.groupby('overlap').agg({'global_halo_exchange_rounds': 'mean'}).reset_index()
    
    overlaps = grouped['overlap'].values
    avg_rounds = grouped['global_halo_exchange_rounds'].values
    
    # Création de la figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Courbe : Communications
    ax.set_xlabel('Overlap', fontsize=12)
    ax.set_ylabel("Nombre d'échanges", fontsize=12)
    ax.plot(overlaps, avg_rounds, 's-', color='red', linewidth=2, markersize=8, label='Communications')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right', fontsize=11)
    
    # Ajout des valeurs sur les points
    for i, (x, y) in enumerate(zip(overlaps, avg_rounds)):
        ax.annotate(f'{y:,}', (x, y), textcoords="offset points", xytext=(0,10), ha='center', fontsize=9)
    
    plt.title('Nombre d\'échanges en fonction de l\'overlap', fontsize=14)
    
    # Ajustement et sauvegarde
    fig.tight_layout()
    output_filename = os.path.join(output_dir, 'benchmark_communications.png')
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.show()
    
    return output_filename

def main():
    # Création du dossier avec la date
    output_dir = create_output_directory()
    
    # Lecture des données (solution 2 avec nettoyage)
    with open('benchmark_overlap_results.txt', 'r') as f:
        lines = f.readlines()
    
    # Supprimer les commentaires et les espaces inutiles
    clean_lines = []
    for line in lines:
        if not line.startswith('#'):
            # Remplacer les espaces multiples par un seul espace
            clean_line = ' '.join(line.split())
            if clean_line.strip():  # Ignorer les lignes vides
                clean_lines.append(clean_line)
    
    # Lecture des données nettoyées
    data = pd.read_csv(io.StringIO('\n'.join(clean_lines)), sep=' ')

    # Génération des 2 graphiques séparés
    time_file = plot_benchmark_time(data, output_dir)
    
    comm_file = plot_benchmark_communications(data, output_dir)

if __name__ == "__main__":
    main()