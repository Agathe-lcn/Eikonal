import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np  # <-- AJOUT: import de numpy

def plot_benchmark():
    # Lecture des données
    data = pd.read_csv('results_fim.txt', comment='#', sep='\t')
    
    # Séparer avec et sans seuil
    without_threshold = data[data['max_radius'] < 0]
    with_threshold = data[data['max_radius'] > 0]
    
    # Trier par taille
    without_threshold = without_threshold.sort_values('n')
    with_threshold = with_threshold.sort_values('n')
    
    # Convertir en tableaux numpy pour matplotlib
    n_without = without_threshold['n'].values
    time_without = without_threshold['avg_time_seconds'].values
    
    n_with = with_threshold['n'].values
    time_with = with_threshold['avg_time_seconds'].values
    
    # Création de la figure
    plt.figure(figsize=(10,7))
    
    # Tracer les courbes avec les tableaux numpy
    plt.plot(n_without, time_without, 'o-', color='blue', label='Sans seuil', linewidth=2, markersize=8)
    plt.plot(n_with, time_with, 's-', color='red', label='Seuil 30 pixels', linewidth=2, markersize=8)
    
    # Mise en forme
    plt.xlabel('Dimension de la grille (n x n)', fontsize=12)
    plt.ylabel('Temps moyen (secondes)', fontsize=12)
    plt.title('Temps de calcul de la FIM en fonction de la taille de la grille', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)
    
    # Sauvegarde et affichage
    plt.tight_layout()
    plt.savefig("fim_benchmark_results.png", dpi=400, bbox_inches="tight")
    plt.show()

if __name__ == "__main__":
    plot_benchmark()