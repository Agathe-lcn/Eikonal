// Génère une grille 2D avec une source aléatoire
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
    int nsources = 1;

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(n, m, h);
    if (!g){
        printf("Erreur: Impossible de créer la grille. \n");
        return 1;
    }

    // Vitesse constante égale à 1
    eikonal_grid_set_speed_constant(g, 1.0);
    
    // Génération de la source
    int* src_i = (int *)malloc(nsources * sizeof(int));
    int* src_j = (int *)malloc(nsources * sizeof(int));

    src_i[0] = n/2;
    src_j[0] = m/2;

    // Fichier pour stocker les informations de la source
    FILE* coord_file = fopen("coords_source.txt", "w");
    if (!coord_file){
        printf("Erreur: Impossible de créer le fichier coord_source.txt. \n");
        free(src_i);
        free(src_j);
        eikonal_grid_free(g);
        return 1;
    }

    double x = src_j[0] * h;
    double y = src_i[0] * h;
    fprintf(coord_file, "%.6f %.6f\n", x, y);
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