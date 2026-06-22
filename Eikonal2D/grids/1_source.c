// Generates a 2D grid with one random source
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
    double length = 1.0;
    double h = length/n;

    // Number of random sources
    int nsources = 1;

    // Creating the grid
    EikonalGrid* g = eikonal_grid_create(n, m, h);
    if (!g){
        printf("Error: Unable to create the grid. \n");
        return 1;
    }

    // Constant speed of 1
    eikonal_grid_set_speed_constant(g, 1.0);
    
    // Generation of random sources
    int* src_i = (int *)malloc(nsources * sizeof(int));
    int* src_j = (int *)malloc(nsources * sizeof(int));

    src_i[0] = n/2;
    src_j[0] = m/2;

    // File for storing source information
    FILE* coord_file = fopen("coords_mesh.txt", "w");
    if (!coord_file){
        printf("Error: Unable to create the coord_mesh.txt file. \n");
        free(src_i);
        free(src_j);
        eikonal_grid_free(g);
        return 1;
    }

    double x = src_j[0] * h;
    double y = src_i[0] * h;
    fprintf(coord_file, "%.6f %.6f\n", x, y);
    fclose(coord_file);

    // Execution of FIM
    fim_solve(g, src_i, src_j, nsources, EPSILON);

    eikonal_save_matrix(g, "matrix_fim.txt");

    // Cleaning
    free(src_i);
    free(src_j);
    eikonal_grid_free(g);
    
    return 0;
}