#include "../include/FIM2D.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define EPSILON 1e-12

int main(){
    // Paramètres
    int n = 200;
    int m = 200;
    double length = 1.0;
    double h = length/n;
    int nsources = 5;

    // Création de la grille
    EikonalGrid* g = eikonal_grid_create(n,m,h);
    if (!g){
        printf("Erreur: impossible de créer la grille\n");
        return 1;
    }

    // vitesse cste = 1
    eikonal_grid_set_speed_constant(g, 1.0);

    // mur vertical au centre
    int wall_x = m/2;
    for (int i=0; i < n; i++){
        eikonal_grid_set_obstacle(g, i, wall_x);
    }
    printf("mur vertical à la colonne %d\n", wall_x);

    // sources aléatoires
    srand(time(NULL));
    int* src_i = (int*)malloc(nsources * sizeof(int));
    int* src_j = (int*)malloc(nsources * sizeof(int));

    FILE* coord_file = fopen("coords_mesh.txt", "w");
    if (!coord_file){
        printf("erreur: impossible de créer coord_mesh.txt\n");
        free(src_i);
        free(src_j);
        eikonal_grid_free(g);
        return 1;
    }

    for (int s = 0; s < nsources; s++){
        src_i[s] = rand() % n;
        src_j[s] = rand() % m;

        // éviter mettre source dans le mur
        while (src_j[s] == wall_x)
            src_j[s] = rand() % m;

        double x = src_j[s] * h;
        double y = src_i[s] * h;
        fprintf(coord_file, "%.6f %.6f\n", x, y);
    }

    fclose(coord_file);

    // exécution FIM
    fim_solve(g, src_i, src_j, nsources, EPSILON);

    // sauvegarde
    FILE* mat_file = fopen("matrix_fim.txt", "w");
    // AJOUTER ERREUR SI OUVRE PAS
    for (int i=0; i<n; i++){
        for (int j=0; j<m; j++){
            int index = i*m+j;
            double T = g->T[index];

            if (g->F[index] <= EIKONAL_EPS || T >= EIKONAL_INF /2.0)
                fprintf(mat_file, "%.12f ", "inf");
            else 
                fprintf(mat_file, "%.12f ", T);
        }
        fprintf(mat_file, "\n");
    }
    fclose(mat_file);

    //nettoyage
    free(src_i);
    free(src_j);
    eikonal_grid_free(g);

    return 0;
}