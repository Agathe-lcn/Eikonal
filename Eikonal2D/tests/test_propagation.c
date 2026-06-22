#include "../include/FIM2D.h"
#include <stdio.h>

int main(){
    int n = 3;
    int m = 3;
    double h = 0.2;
    EikonalGrid* g = eikonal_grid_create(n,m,h);
    eikonal_grid_set_speed_constant(g,1.0);

    // Source in the center
    int src_i[] = {1};
    int src_j[] = {1};

    // Manual initialization of the source
    int index = src_i[0] * m + src_j[0];
    g->T[index] = 0.0;

    printf("Before FIM:\n");
    for (int i=0; i < n; i++){
        for (int j=0; j<m; j++)
            printf("%8.4f ", g->T[i * m + j]);
    }
    printf("\n");

    fim_solve(g, src_i, src_j, 1, 1e-12);

    printf("\nAfter FIM:\n");
    for (int i=0; i < n; i++){
        for (int j=0; j < m; j++)
            printf("%8.4f ", g->T[i * m + j]);
        printf("\n");
    }

    eikonal_grid_free(g);
    return 0;
}