// Generates a 2D grid with five random sources with one wall
// Saves the time matrix to a .txt file

#include "../include/FIM2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 1e-12

int main(){
    // Grid settings
    int n = 200;
    int m = 200;
    double h = 1.0/n;
    int nsources = 5;

    // Creating the grid
    EikonalGrid* g = eikonal_grid_create(n, m, h);
    if (!g){
        printf("Error: Unable to create the grid.\n");
        return 1;
    }

    // Constant speed of 1
    eikonal_grid_set_speed_constant(g, 1.0);

    // Vertical wall in the center
    int wall_x = m/2;
    for (int i=0; i < n; i++){
        eikonal_grid_set_obstacle(g, i, wall_x);
    }

    // File for storing source information
    FILE* coord_file = fopen("coords_source.txt", "w");
    if (!coord_file){
        printf("Error: Unable to create the coord_source.txt file.\n");
        eikonal_grid_free(g);
        return 1;
    }

    // Generation of the random sources
    srand(time(NULL));
    int* src_i = (int*)malloc(nsources * sizeof(int));
    int* src_j = (int*)malloc(nsources * sizeof(int));

    for (int s = 0; s < nsources; s++){
        src_i[s] = rand() % n;
        src_j[s] = rand() % m;

        // Avoid placing the sources in the wall
        while (src_j[s] == wall_x)
            src_j[s] = rand() % m;

        double x = src_j[s] * h;
        double y = src_i[s] * h;
        fprintf(coord_file, "%.6f %.6f\n", x, y);
    }

    fclose(coord_file);

    // Execution of FIM
    fim_solve(g, src_i, src_j, nsources, EPSILON, -1.0);

    eikonal_save_matrix(g, "matrix_fim.txt");

    // Cleaning
    free(src_i);
    free(src_j);
    eikonal_grid_free(g);

    return 0;
}