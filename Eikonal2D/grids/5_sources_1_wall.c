// Génère une grille 2D contenant 5 sources aléatoires et un mur.
// Enregistre la matrice des temps dans un fichier .txt

#include "../include/FIM2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 1e-12

int main(){
    // Paramètres de la grille
    int n = 200;
    int m = 200;
    double h = 1.0/n;
    int nsources = 5;

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(n, m, h);
    if (!g){
        printf("Erreur: Impossible de créer la grille.\n");
        return 1;
    }

    // Vitesse constante égale à 1
    eikonal_grid_set_speed_constant(g, 1.0);

    // Mur vertical au centre
    int wall_x = m/2;
    for (int i=0; i < n; i++){
        eikonal_grid_set_obstacle(g, i, wall_x);
    }

    // Fichier pour stocker les informations sur la source
    FILE* coord_file = fopen("coords_source.txt", "w");
    if (!coord_file){
        printf("Erreur: Impossible de créer le fichier coord_source.txt.\n");
        eikonal_grid_free(g);
        return 1;
    }

    // Génération des sources de façon aléatoire
    srand(time(NULL));
    int* src_i = (int*)malloc(nsources * sizeof(int));
    int* src_j = (int*)malloc(nsources * sizeof(int));

    for (int s = 0; s < nsources; s++){
        src_i[s] = rand() % n;
        src_j[s] = rand() % m;

        // On évite de placer les sources dans le mur
        while (src_j[s] == wall_x)
            src_j[s] = rand() % m;

        double x = src_j[s] * h;
        double y = src_i[s] * h;
        fprintf(coord_file, "%.6f %.6f\n", x, y);
    }

    fclose(coord_file);

    // Exécution de la FIM
    fim_solve(g, src_i, src_j, nsources, EPSILON, -1.0, -1);

    eikonal_save_matrix(g, "matrix_fim.txt");

    // Nettoyage
    free(src_i);
    free(src_j);
    eikonal_grid_free(g);

    return 0;
}