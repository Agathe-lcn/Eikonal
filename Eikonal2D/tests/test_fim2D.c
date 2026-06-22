#include "../include/FIM2D.h"

#include <stdio.h>

int main(){
    printf("\n Test FIM2D\n");

    // Simple grid
    int n = 5;
    int m = 5;
    EikonalGrid* g = eikonal_grid_create(n, m, 0.2);
    if (!g)
        return 1;

    // Constant speed
    eikonal_grid_set_speed_constant(g, 1.0);

    // A single source at the center
    int src_i[] = {2};
    int src_j[] = {2};

    // FIM Execution
    fim_solve(g, src_i, src_j, 1, 1e-12);

    for (int i=0; i < n; i++){
        for (int j=0; j < m; j++)
            printf("%8.4f ", g->T[i * m + j]);
    }
    printf("\n");

    eikonal_grid_free(g);
    return 0;
}