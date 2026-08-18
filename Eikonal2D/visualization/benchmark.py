import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def plot_dataset(data, dataset_name, output_filename):
    # On ne garde que les lignes correspondant à ce dataset
    subset = data[data['dataset'] == dataset_name]

    # Séparer avec et sans seuil
    without_threshold = subset[subset['max_depth'] < 0]
    with_threshold = subset[subset['max_depth'] > 0]

    # Trier par taille
    without_threshold = without_threshold.sort_values('size')
    with_threshold = with_threshold.sort_values('size')

    # Convertir en tableaux numpy pour matplotlib
    n_without = without_threshold['size'].values
    time_without = without_threshold['avg_time_seconds'].values
    n_with = with_threshold['size'].values
    time_with = with_threshold['avg_time_seconds'].values

    # Récupère la valeur de max_depth utilisée (pour la légende), si elle existe
    depth_value = with_threshold['max_depth'].iloc[0] if len(with_threshold) > 0 else None

    # Création de la figure
    plt.figure(figsize=(10, 7))

    # Tracer les courbes avec les tableaux numpy
    plt.plot(n_without, time_without, 'o-', color='blue', label='Sans seuil', linewidth=2, markersize=8)
    label_with = f'Seuil {int(depth_value)} pixels' if depth_value is not None else 'Avec seuil'
    plt.plot(n_with, time_with, 's-', color='red', label=label_with, linewidth=2, markersize=8)

    # Mise en forme
    plt.xlabel('Dimension de la grille (n x n)', fontsize=12)
    plt.ylabel('Temps moyen (secondes)', fontsize=12)
    plt.title(f'Temps de calcul de la FIM en fonction de la taille de la grille ({dataset_name})', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)

    # Sauvegarde et affichage
    plt.tight_layout()
    plt.savefig(output_filename, dpi=400, bbox_inches="tight")
    plt.show()


def plot_benchmark():
    # Lecture des données
    data = pd.read_csv('results_fim.txt', comment='#', sep='\t')

    # Un graphique pour le scénario en cercle
    plot_dataset(data, 'circle', 'fim_benchmark_circle.png')

    # Un graphique pour le scénario en ligne
    plot_dataset(data, 'line', 'fim_benchmark_line.png')


if __name__ == "__main__":
    plot_benchmark()